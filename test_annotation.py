from __future__ import annotations

import argparse
from pathlib import Path

from antibody_dataset.annotation import annotate_fasta


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Тестовая аннотация FASTA с помощью ANARCI"
    )
    parser.add_argument(
        "fasta",
        type=Path,
        help="Путь к FASTA-файлу",
    )

    args = parser.parse_args()

    if not args.fasta.is_file():
        raise FileNotFoundError(
            f"FASTA-файл не найден: {args.fasta}"
        )

    annotations = annotate_fasta(args.fasta)

    print(
        "chain_id\tchain_type\trole\t"
        "status\tnanobody\theader"
    )

    for annotation in annotations:
        for chain_id in annotation["chain_ids"]:
            print(
                f"{chain_id}\t"
                f"{annotation['chain_type'] or '-'}\t"
                f"{annotation['role']}\t"
                f"{annotation['status']}\t"
                f"{annotation['is_nanobody']}\t"
                f"{annotation['header']}"
            )


if __name__ == "__main__":
    main()