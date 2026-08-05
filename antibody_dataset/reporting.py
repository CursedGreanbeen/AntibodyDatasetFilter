import csv
from pathlib import Path
from typing import Iterable


FIELDNAMES = [
    "structure",
    "fasta",
    "chain",
    "residue_count",
    "first_residue",
    "last_residue",
    "sequence",
]


def write_tsv(
    rows: Iterable[dict],
    output_path: Path,
) -> None:
    """Записывает строки в TSV-файл."""
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(
            file,
            fieldnames=FIELDNAMES,
            delimiter="\t",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)
