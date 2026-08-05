from pathlib import Path

from antibody_dataset.chains import get_structure_chains, read_structure
from antibody_dataset.io import (
    build_fasta_index,
    find_cif_files,
)
from antibody_dataset.reporting import write_tsv


PROJECT_DIR = Path(__file__).resolve().parent.parent
CIF_DIR = PROJECT_DIR / "CIFs-filtered-new"
FASTA_DIR = PROJECT_DIR / "fasta-filtered-new"
OUTPUT_PATH = PROJECT_DIR / "results" / "reports" / "chain_inventory.tsv"


def get_chain_sequences(cif_path: Path) -> dict[str, str]:
    """Возвращает последовательности цепей по остаткам с координатами."""
    structure = read_structure(cif_path)
    model = structure[0]

    sequences = {}

    for chain in model:
        residues = [
            residue
            for residue in chain.get_polymer()
            if len(residue) > 0
        ]

        sequence = []

        for residue in residues:
            sequence_info = residue.name
            sequence.append(
                {
                    "ALA": "A",
                    "ARG": "R",
                    "ASN": "N",
                    "ASP": "D",
                    "CYS": "C",
                    "GLN": "Q",
                    "GLU": "E",
                    "GLY": "G",
                    "HIS": "H",
                    "ILE": "I",
                    "LEU": "L",
                    "LYS": "K",
                    "MET": "M",
                    "PHE": "F",
                    "PRO": "P",
                    "SER": "S",
                    "THR": "T",
                    "TRP": "W",
                    "TYR": "Y",
                    "VAL": "V",
                }.get(sequence_info, "X")
            )

        sequences[chain.name] = "".join(sequence)

    return sequences


def main() -> None:
    fasta_index = build_fasta_index(FASTA_DIR)
    rows = []

    cif_files = find_cif_files(CIF_DIR)

    if not cif_files:
        raise FileNotFoundError(
            f"В папке не найдены CIF-файлы: {CIF_DIR}"
        )

    for cif_path in cif_files:
        fasta_path = fasta_index.get(cif_path.stem)
        fasta_name = fasta_path.name if fasta_path else ""

        try:
            chain_infos = get_structure_chains(cif_path)
            sequences = get_chain_sequences(cif_path)

            for chain_info in chain_infos:
                rows.append(
                    {
                        "structure": cif_path.name,
                        "fasta": fasta_name,
                        "chain": chain_info.chain_id,
                        "residue_count": chain_info.residue_count,
                        "first_residue": (
                            f"{chain_info.first_residue.number}"
                            f"{chain_info.first_residue.insertion_code}"
                        ),
                        "last_residue": (
                            f"{chain_info.last_residue.number}"
                            f"{chain_info.last_residue.insertion_code}"
                        ),
                        "sequence": sequences.get(
                            chain_info.chain_id,
                            "",
                        ),
                    }
                )

        except Exception as error:
            rows.append(
                {
                    "structure": cif_path.name,
                    "fasta": fasta_name,
                    "chain": "ERROR",
                    "sequence": str(error),
                }
            )

    write_tsv(rows, OUTPUT_PATH)

    print(f"CIF-файлов обработано: {len(cif_files)}")
    print(f"Цепей записано: {len(rows)}")
    print(f"Отчёт: {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
