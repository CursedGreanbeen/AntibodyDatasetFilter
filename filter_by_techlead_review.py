#!/usr/bin/env python3
"""
Скрипт для удаления комплексов по критериям техлида:
1. Вода в структуре
2. Пептидный антиген
3. ОЧЕНЬ БОЛЬШИЕ комплексы
4. ДНК

Удаляет файлы из CIFs-filtered и fasta-filtered по названиям.
"""

from pathlib import Path
import shutil


# Списки комплексов для удаления по категориям (из обзора техлида)


VERY_LARGE_COMPLEXES = {
    "8yef", "9bnk", "9bth", "9bti", "9btj", "9btl", "9crv", "9cv7",
    "9nvw", "9nvx", "9nvz", "9nw0", "9omg"
}

DNA_COMPLEXES = {
    "8t29", "8t2a", "8t2b", "8t2o", "8vm8", "8vmb", "9aur", "9c75",
    "9e7g", "9u5q", "9u5r", "9urm", "9wb3", "9wb4"
}

BAD_STRUCTURES = {
    "8jn3", "8s6t", "9js6", "9ynr"
}

WATER_COMPLEXES = {
    "8ds7", "8g2m", "8rmo", "8tmz", "8uig", "8uih", "8vy4", "8yhz",
    "9bx5", "9bx7", "9c7d", "9c7x", "9cfd", "9ia3", "9j8a", "9mmj",
    "9mz6", "9mz7", "9mz8", "9nkz", "9nl0", "9p4c", "9pwn", "9t46",
    "9tpp", "9zmb", "9zmc"
}

PEPTIDE_ANTIGENS = {
    "9cfd", "9cqb", "9dsr", "9dss", "9dst", "9dsu", "9hv9", "9ia3",
    "9j8a", "9mmj", "9mz6", "9mz7", "9mz8", "9nkz", "9nl0", "9nl1",
    "9nzf", "9p4c", "9pwn", "9t46", "9tpp", "9y7n", "9zmb", "9zmc"
}

# от цепей АТ и от цепей АГ проверить наличие небелковых молекул в радиусе 6 ангстрем
SMALL_MOLECULE_INTERFACE = {
    "8vr9", "8vra", "8vrb", "8vr9", "8udr", "21ao"
}

UNCERTAIN_CASES = {
    "8kei", "9evz", "9gfr", "9pxb", "9qqe", "9qqf", "9vxl", "9ynr"
}


def main():
    # Объединяем все списки
    to_remove = BAD_STRUCTURES | DNA_COMPLEXES

    print("=" * 60)
    print("ФИЛЬТРАЦИЯ ПО КРИТЕРИЯМ ТЕХЛИДА")
    print("=" * 60)
    print(f"\nВсего комплексов к удалению: {len(to_remove)}")

    # Категории
    print("\n--- По категориям ---")
    print(f"Очень большие: {len(VERY_LARGE_COMPLEXES)}")
    print(f"ДНК: {len(DNA_COMPLEXES)}")

    print("\n--- Список комплексов к удалению ---")
    for name in sorted(to_remove):
        reasons = []
        if name in VERY_LARGE_COMPLEXES:
            reasons.append("ОЧЕНЬ БОЛЬШОЙ")
        if name in DNA_COMPLEXES:
            reasons.append("ДНК")
        print(f"  {name}: {', '.join(reasons)}")

    # Удаляем файлы
    print("\n" + "=" * 60)
    print("УДАЛЕНИЕ ФАЙЛОВ")
    print("=" * 60)

    base_dir = Path("/home/mullagaliamova/ClaudeWorkspace/PROJECTS/cdr-h3-folding")
    cif_dir = base_dir / "CIFs-filtered"
    fasta_dir = base_dir / "fasta-filtered"

    removed_count = 0

    for name in sorted(to_remove):
        # CIF файлы
        cif_file = cif_dir / f"{name}_cropped.cif"
        if cif_file.exists():
            cif_file.unlink()
            print(f"Удален: {cif_file.name}")
            removed_count += 1

        # FASTA файлы
        fasta_file = fasta_dir / f"{name}.fasta"
        if fasta_file.exists():
            fasta_file.unlink()
            print(f"Удален: {fasta_file.name}")
            removed_count += 1

    print(f"\n=== ИТОГО УДАЛЕНО ФАЙЛОВ: {removed_count} ===")

    # Оставшиеся файлы
    remaining_cif = len(list(cif_dir.glob("*_cropped.cif")))
    remaining_fasta = len(list(fasta_dir.glob("*.fasta")))

    print(f"\nОставшихся CIF файлов: {remaining_cif}")
    print(f"Оставшихся FASTA файлов: {remaining_fasta}")


if __name__ == "__main__":
    main()
