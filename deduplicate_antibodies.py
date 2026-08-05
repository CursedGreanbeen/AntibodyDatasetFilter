#!/usr/bin/env python3
"""
Скрипт для дедупликации антителых комплексов.

Если H+L пара последовательностей встречается в нескольких комплексах,
оставляется только первый комплекс (по алфавиту), остальные удаляются.
"""

from pathlib import Path
from collections import defaultdict
import csv


def main():
    table_file = Path("/home/mullagaliamova/ClaudeWorkspace/PROJECTS/cdr-h3-folding/FOLDING_DATA/sequence_occurrence_table.tsv")
    base_dir = Path("/home/mullagaliamova/ClaudeWorkspace/PROJECTS/cdr-h3-folding")
    cif_dir = base_dir / "CIFs-filtered-new"
    fasta_dir = base_dir / "fasta-filtered-new"

    # Структура: (h_seq, l_seq) -> list of (complex_name, h_seq, l_seq)
    ab_pairs = defaultdict(list)

    # Читаем таблицу
    with open(table_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            seq_type = row['Chain_Type']
            if seq_type not in ('H', 'L'):
                continue

            # Парсим колонку Complexes
            complexes_str = row['Complexes']
            complexes = [c.strip() for c in complexes_str.split(';')]

            # Для каждого комплекса добавляем информацию о цепи
            for complex_entry in complexes:
                # Формат: "9k0c(A)" или "9k0c(unknown)"
                parts = complex_entry.split('(')
                if len(parts) != 2:
                    continue
                complex_name = parts[0].strip()
                chain_id = parts[1].rstrip(')').strip()

                ab_pairs[complex_name].append({
                    'chain_type': seq_type,
                    'sequence': row['Sequence'],
                    'chain_id': chain_id
                })

    # Теперь группируем по парам (h_seq, l_seq)
    pair_to_complexes = defaultdict(list)

    for complex_name, chains in ab_pairs.items():
        h_seq = None
        l_seq = None
        for chain in chains:
            if chain['chain_type'] == 'H':
                h_seq = chain['sequence']
            elif chain['chain_type'] == 'L':
                l_seq = chain['sequence']

        if h_seq and l_seq:
            pair_key = (h_seq, l_seq)
            pair_to_complexes[pair_key].append(complex_name)

    # Находим дубликаты (пары, встречающиеся в >1 комплексах)
    duplicates = {pair: comps for pair, comps in pair_to_complexes.items() if len(comps) > 1}

    print("=" * 70)
    print("ДЕДУПЛИКАЦИЯ АНТИTEЛЬНЫХ КОМПЛЕКСОВ")
    print("=" * 70)
    print(f"\nВсего уникальных H+L пар: {len(pair_to_complexes)}")
    print(f"Пар с дубликатами (>1 комплекс): {len(duplicates)}")

    # Определяем, какие комплексы удалить
    to_remove = set()
    kept = set()

    for pair, complexes in duplicates.items():
        sorted_comps = sorted(complexes)
        keep = sorted_comps[0]  # первый по алфавиту
        remove = sorted_comps[1:]  # остальные

        kept.add(keep)
        to_remove.update(remove)

        print(f"\nH+L пара в {len(complexes)} комплексах:")
        print(f"  ОСТАВЛЕН: {keep}")
        if remove:
            print(f"  УДАЛЕНЫ: {', '.join(remove)}")

    print("\n" + "=" * 70)
    print("УДАЛЕНИЕ ФАЙЛОВ")
    print("=" * 70)
    print(f"\nКомплексов к удалению: {len(to_remove)}")
    print(f"Комплексов оставлено: {len(kept)}")

    removed_count = 0

    for name in sorted(to_remove):
        # CIF файлы
        cif_file = cif_dir / f"{name}.cif"
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
    remaining_cif = len(list(cif_dir.glob("*.cif")))
    remaining_fasta = len(list(fasta_dir.glob("*.fasta")))

    print(f"\nОставшихся CIF файлов: {remaining_cif}")
    print(f"Оставшихся FASTA файлов: {remaining_fasta}")


if __name__ == "__main__":
    main()
