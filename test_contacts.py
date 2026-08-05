from __future__ import annotations

import argparse
import sys
from pathlib import Path

import gemmi

# Add scripts directory to path for imports
scripts_dir = Path(__file__).parent.resolve()
if str(scripts_dir) not in sys.path:
    sys.path.insert(0, str(scripts_dir))

from antibody_dataset.annotation import annotate_fasta
from antibody_dataset.contacts import calculate_contacts
from antibody_dataset.filtering import classify_complex


def display(value: object) -> str:
    """Заменяет None на прочерк."""
    return "-" if value is None else str(value)


def get_polymer_chain_ids(structure) -> set[str]:
    """Возвращает chain ID цепей с полимерным содержимым."""
    if len(structure) == 0:
        return set()

    model = structure[0]
    polymer_chain_ids: set[str] = set()

    for chain in model:
        polymer = chain.get_polymer()

        if len(polymer) > 0:
            polymer_chain_ids.add(chain.name)

    return polymer_chain_ids


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Тестовый расчёт межцепочечных контактов"
    )
    parser.add_argument(
        "cif",
        type=Path,
        help="Путь к CIF/mmCIF-файлу",
    )
    parser.add_argument(
        "fasta",
        type=Path,
        help="Путь к FASTA-файлу",
    )
    parser.add_argument(
        "--cutoff",
        type=float,
        default=5.0,
        help="Порог расстояния контакта в A",
    )

    args = parser.parse_args()

    if not args.cif.is_file():
        raise FileNotFoundError(
            f"CIF-файл не найден: {args.cif}"
        )

    if not args.fasta.is_file():
        raise FileNotFoundError(
            f"FASTA-файл не найден: {args.fasta}"
        )

    structure = gemmi.read_structure(str(args.cif))
    annotations = annotate_fasta(args.fasta)

    print("Аннотация цепей:")
    print("chain_id\tchain_type\trole\tstatus")

    for annotation in annotations:
        for chain_id in annotation["chain_ids"]:
            print(
                f"{chain_id}\t"
                f"{display(annotation['chain_type'])}\t"
                f"{annotation['role']}\t"
                f"{annotation['status']}"
            )

    contacts = calculate_contacts(
        structure=structure,
        annotations=annotations,
        cutoff=args.cutoff,
    )

    polymer_chain_ids = get_polymer_chain_ids(structure)
    filter_result = classify_complex(
        annotations=annotations,
        contacts=contacts,
        other_polymer_chain_ids=polymer_chain_ids,
    )

    print("\nКлассификация:")
    for key, value in filter_result.as_dict().items():
        print(f"  {key}: {value}")

    print(f"\nКонтакты, cutoff = {args.cutoff:.1f} A:")
    print(
        "chain_a\trole_a\tchain_b\trole_b\t"
        "atom_contacts\tresidue_contacts_a\t"
        "residue_contacts_b\tresidue_contacts\t"
        "minimum_distance"
    )

    for contact in sorted(
        contacts,
        key=lambda item: (
            item.chain_a,
            item.chain_b,
        ),
    ):
        data = contact.as_dict()

        print(
            f"{data['chain_a']}\t"
            f"{data['role_a']}\t"
            f"{data['chain_b']}\t"
            f"{data['role_b']}\t"
            f"{data['atom_contacts']}\t"
            f"{data['residue_contacts_a']}\t"
            f"{data['residue_contacts_b']}\t"
            f"{data['residue_contacts']}\t"
            f"{display(data['minimum_distance'])}"
        )


if __name__ == "__main__":
    main()
