import argparse
import json
import os
import string
from pathlib import Path
from typing import Any

import yaml
from oneq.api import Gigabytes, OneQ, OneQConfig

DEFAULT_VERSION = os.environ.get("ALPHAFOLD_RUNNER_VERSION") or "latest"
CACHE_DIR = "/mnt/msa-server"


def parse_fasta(fasta: Path, remove_prefix: bool = False) -> dict[str, str]:
    # read fasta lines
    lines = fasta.read_text().splitlines()

    # define indices of lines that contain sequence identifiers
    border_lines = [i for i, line in enumerate(lines) if line.startswith(">")]
    border_lines.append(len(lines))

    # collect lines between identifiers into sequences
    sequences = {}
    for start, end in zip(border_lines[:-1], border_lines[1:]):
        sequence_id = lines[start][1:]
        if remove_prefix:
            sequence_id = sequence_id.split(":")[-1]
        sequence = "".join(map(lambda line: line.strip(), lines[start + 1 : end]))
        sequences[sequence_id] = sequence

    # reassign new ids if any of current chain ids has length != 1 (AlphaFold requirement)
    if any((len(chain_id) != 1 for chain_id in sequences.keys())):
        sequences = {
            string.ascii_uppercase[i]: sequence
            for i, sequence in enumerate(sequences.values())
        }

    return sequences


def json_data_from_fasta(fasta: Path, model_seeds: list[int] = [1]) -> dict[str, Any]:
    # parse fasta and remove prefixes if present
    # (older versions required UID:CHAIN_ID format)
    sequences = parse_fasta(fasta, remove_prefix=True)

    # form dictionary for AlphaFold 3
    json_dict = {
        "name": fasta.stem,
        "sequences": [
            {"protein": {"id": id, "sequence": sequence}}
            for id, sequence in sequences.items()
        ],
        "modelSeeds": model_seeds,
        "dialect": "alphafold3",
        "version": 1,
    }
    return json_dict

def af3_input_sanitised_name(name: str) -> str:
    # Copied from Input.sanitized_name from alphafold3/src/alphafold3/common/folding_input.py::919
    lower_spaceless_name = name.lower().replace(' ', '_')
    allowed_chars = set(string.ascii_lowercase + string.digits + '_-.')
    return ''.join(l for l in lower_spaceless_name if l in allowed_chars)


def create_command(
    json_path: Path, use_xla: bool = True, max_template_date: str = "2042-01-01", ignore_msa_cache: bool = False,
) -> str:
    # JSON evolution: input -> precomputed msa enrichment -> AF3 pipeline enrichment
    input_json_path = f"$DATADIR/{json_path.name}"
    enriched_json_path = f"$DATADIR/{json_path.stem}_enriched.json"
    task_name = af3_input_sanitised_name(json.loads(json_path.read_text())["name"])
    af3_output_json_path = f"$DATADIR/{task_name}/{task_name}_data.json"

    manifest_path = f"$DATADIR/{json_path.stem}_enriched.json.manifest.json"

    commands = [f"export MSA_SERVER_CACHE_DIR={CACHE_DIR}"]
    # Step 1: enrich input json with precomputed MSA and templates
    if not ignore_msa_cache:
        enrich_cmd = [
            "msa-server enrich",
            f"--input {input_json_path}",
            f"--output {enriched_json_path}",
            "--identity-threshold 1",
        ]
        commands.append(' '.join(enrich_cmd))

    # Step 2: run AF3 pipeline (use enriched json if available, otherwise original input)
    af3_input = enriched_json_path if not ignore_msa_cache else input_json_path
    af3_cmd = [
        "cd /app2/alphafold &&",
        "python run_alphafold.py",
        f"--json_path={af3_input}",
        "--output_dir=$DATADIR/",
        "--model_dir=/mnt/models/",
        "--db_dir=/mnt/db/af3_data/",
        # TODO: default database name changed to `mmcif_files`
        "--pdb_database_path=/mnt/db/af3_data/pdb_2022_09_28_mmcif_files.tar",
        f"--max_template_date={max_template_date}",
    ]
    if use_xla:
        af3_cmd.append("--flash_attention_implementation=xla")

    commands.append(' '.join(af3_cmd))

    # Step 3: ingest final json file into MSA cache
    ingest_cmd = [
        "msa-server ingest",
        f"--af3-output {af3_output_json_path}",
    ]
    if not ignore_msa_cache:
        # Use manifest from enrichment step to skip already-cached sequences
        ingest_cmd.append(f"--manifest {manifest_path}")
    else:
        # No enrichment was done, so ingest all sequences
        ingest_cmd.append("--include-all")
    commands.append(' '.join(ingest_cmd))

    # Step 4: print current cache stats
    commands.append("msa-server stats")

    return ' && '.join(commands)


def get_oneq_config(config_path: Path) -> OneQConfig:
    config_yaml = yaml.full_load(config_path.open("r"))
    oneq_config = OneQConfig(**config_yaml)
    return oneq_config


def start_oneq_task(
    oneq_config_path: Path,
    json_path: Path,
    max_template_date: str = "2042-01-01",
    ignore_msa_cache: bool = False,
    cpu: int = 8,
    gpu_model: str = "v100",
    priority: int = 0,
    wait: bool = True,
    image: str = "dock.biocad.ru/alphafold-oneq:latest",
) -> None:
    # special parameters for v100 GPU, not needed for a100, a100i3, a10
    use_xla = False
    env = None
    if gpu_model == "v100":
        use_xla = True
        env = {
            "XLA_FLAGS": "--xla_disable_hlo_passes=custom-kernel-fusion-rewriter",
        }

    command = create_command(json_path, use_xla, max_template_date, ignore_msa_cache)

    # start OneQ task
    oneq_config = get_oneq_config(oneq_config_path)
    oneq = OneQ(oneq_config, service_name=Path(image).name)

    task_id = oneq.start_task(
        command=command,
        cpu=cpu,
        gpu=1,
        gpu_model=gpu_model,
        priority=priority,
        memory=Gigabytes(192),
        files=[json_path.absolute()],
        image=image,
        env=env,
        volumes=[
            ("453b8b39-1863-4505-a734-27cab8ac21a8", "/mnt/db"),
            ("228eda68-ad8c-422e-9308-d079c6b9a5b2", "/mnt/models"),
            ("54fd7dfa-9847-4874-adb8-65552088a1da", CACHE_DIR),
        ],
    )
    print(f"TASK_ID: {task_id}")
    if wait:
        oneq.handle_tasks([task_id], json_path.parent, unarchive=True)
        print(f"Done! Results are saved to {json_path.parent}")
    else:
        print(f"Use oneq get-results {task_id} to download folding results")


def main() -> None:
    parser = argparse.ArgumentParser(description="Alphafold runner")
    parser.add_argument(
        "-i",
        "--input-data",
        type=Path,
        help="Fasta file with chains to fold or JSON file suitable for AlphaFold3",
        required=True,
    )
    parser.add_argument(
        "--ignore-msa-cache",
        action="store_true",
        help="Flag to ignore cached MSA and templates and recompute them",
    )
    parser.add_argument(
        "--max-template-date",
        type=str,
        default="2042-01-01",
        help="Max template date in the YYYY-MM-DD format",
    )
    parser.add_argument(
        "--image-tag",
        type=str,
        default=DEFAULT_VERSION,
        help="Tag of the alphafold-oneq image",
    )
    parser.add_argument(
        "--oneq-config-path",
        type=Path,
        required=False,
        default=Path.home() / ".config/oneq/config.yaml",
    )
    parser.add_argument(
        "--cpu",
        type=int,
        required=False,
        default=8,
        help="Number of CPUs",
    )
    parser.add_argument(
        "--gpu-model",
        choices=["v100", "a100", "a100i3", "a10"],
        required=False,
        default="v100",
        help="GPU model used for inference",
    )
    parser.add_argument(
        "--priority",
        type=int,
        required=False,
        default=0,
        help="OneQ task proirity, a number from 0 to 9.",
    )
    parser.add_argument(
        "--memory",
        type=int,
        required=False,
        default=96,
        help="Memory in GB for the OneQ task (default: 128)",
    )
    parser.add_argument(
        "--wait",
        action="store_true",
        help="Wait for the alphafold task to finish and download folding results alongside fasta file",
    )

    args = parser.parse_args()

    input_data: Path = args.input_data
    if input_data.suffix == ".fasta":
        json_data = json_data_from_fasta(input_data)
        json_path = input_data.parent / f"{input_data.stem}.json"
        if json_path.is_file():
            raise FileExistsError(
                f"File {json_path} already exists. Move it somewhere to avoid overwriting."
            )
        json_path.write_text(json.dumps(json_data))

    elif input_data.suffix == ".json":
        json_path = input_data

    else:
        raise ValueError(
            f"Extension `{input_data.suffix}` not supported for input data."
        )

    start_oneq_task(
        oneq_config_path=args.oneq_config_path,
        json_path=json_path,
        max_template_date=args.max_template_date,
        ignore_msa_cache=args.ignore_msa_cache,
        cpu=args.cpu,
        gpu_model=args.gpu_model,
        priority=args.priority,
        wait=args.wait,
        image=f"dock.biocad.ru/alphafold-oneq:{args.image_tag}",
    )


if __name__ == "__main__":
    main()
