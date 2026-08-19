from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

import gemmi

from antibody_dataset.annotation import annotate_fasta
from antibody_dataset.contacts import calculate_contacts
from antibody_dataset.filtering import classify_complex, find_antibody_pairs
from antibody_dataset.repeats import detect_repeated_units
from antibody_dataset.chains import get_model_chains, read_structure
from antibody_dataset.repeats import detect_repeated_units


def find_fasta(
    fasta_dir: Path,
    structure_path: Path,
) -> Path | None:
    """Looks for FASTA with the same stem as CIF.

    If structure has _1 or _2 suffix, also tries base name.
    """
    stem = structure_path.stem

    # Try exact stem first
    candidates = [
        fasta_dir / f"{stem}.fasta",
    ]

    # If stem has _1 or _2 suffix, also try base name
    if "_" in stem:
        base_name = stem.rsplit("_", 1)[0]
        candidates.append(fasta_dir / f"{base_name}.fasta")

    for candidate in candidates:
        if candidate.is_file():
            return candidate

    return None


def process_structure(
    structure_path: Path,
    fasta_path: Path,
    cutoff: float,
) -> dict:
    """Processes one CIF/FASTA pair."""
    structure = read_structure(structure_path)
    structure_chains = get_model_chains(structure)
    polymer_chain_ids = {chain.chain_id for chain in structure_chains}
    annotations = annotate_fasta(fasta_path)

    contacts = calculate_contacts(
        structure=structure,
        polymer_chain_ids=polymer_chain_ids,
        annotations=annotations,
        cutoff=cutoff,
    )

    antibody_pairs = find_antibody_pairs(
        contacts=contacts,
    )

    repeat_analysis = detect_repeated_units(
        annotations=annotations,
        contacts=contacts,
        antibody_pairs=antibody_pairs,
        polymer_chain_ids=polymer_chain_ids,
    )

    filter_result = classify_complex(
        polymer_chain_ids=polymer_chain_ids,
        annotations=annotations,
        contacts=contacts,
        antibody_pairs=antibody_pairs,
        repeat_analysis=repeat_analysis,
    )

    row = filter_result.as_dict()
    row["structure"] = structure_path.stem
    row["fasta"] = fasta_path.name

    return row


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Mass classification of CIF complexes"
    )
    parser.add_argument(
        "--cif-dir",
        type=Path,
        default=Path("CIFs-filtered-new"),
    )
    parser.add_argument(
        "--fasta-dir",
        type=Path,
        default=Path("fasta-filtered-new"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "results/reports/filtering.tsv"
        ),
    )
    parser.add_argument(
        "--cutoff",
        type=float,
        default=4.0,
    )

    args = parser.parse_args()

    structure_paths = sorted(
        args.cif_dir.glob("*.cif")
    )

    rows: list[dict] = []
    errors: list[dict] = []

    for structure_path in structure_paths:
        fasta_path = find_fasta(
            args.fasta_dir,
            structure_path,
        )

        if fasta_path is None:
            errors.append(
                {
                    "structure": structure_path.stem,
                    "error": "fasta_not_found",
                }
            )
            continue

        try:
            row = process_structure(
                structure_path=structure_path,
                fasta_path=fasta_path,
                cutoff=args.cutoff,
            )
            rows.append(row)

            print(
                f"{structure_path.stem}: "
                f"{row['classification']}"
            )

        except Exception as error:
            errors.append(
                {
                    "structure": structure_path.stem,
                    "error": repr(error),
                }
            )
            print(
                f"{structure_path.stem}: ERROR: {error}"
            )

    args.output.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    if rows:
        fieldnames = list(rows[0].keys())

        with args.output.open(
            "w",
            newline="",
            encoding="utf-8",
        ) as file:
            writer = csv.DictWriter(
                file,
                fieldnames=fieldnames,
                delimiter="\t",
            )
            writer.writeheader()
            writer.writerows(rows)

    error_path = args.output.with_name(
        f"{args.output.stem}_errors.tsv"
    )

    if errors:
        with error_path.open(
            "w",
            newline="",
            encoding="utf-8",
        ) as file:
            writer = csv.DictWriter(
                file,
                fieldnames=["structure", "error"],
                delimiter="\t",
            )
            writer.writeheader()
            writer.writerows(errors)

    print()
    print(f"Processed: {len(rows)}")
    print(f"Errors: {len(errors)}")
    print(f"Report: {args.output}")


if __name__ == "__main__":
    main()
