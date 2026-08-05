#!/usr/bin/env python3
"""
Скрипт для создания таблицы последовательностей с информацией о дублировании.
Особенно отслеживает полные дублирования антител (H + L цепи вместе).
"""

from pathlib import Path
from collections import defaultdict
import re


def parse_fasta(fasta_path: Path) -> list[dict]:
    """
    Парсит FASTA файл и возвращает список записей с информацией.
    """
    entries = []
    current_header = None
    current_sequence = []

    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.rstrip()
            if line.startswith('>'):
                if current_header:
                    entries.append({
                        'header': current_header,
                        'sequence': ''.join(current_sequence)
                    })
                current_header = line
                current_sequence = []
            else:
                current_sequence.append(line)

        if current_header:
            entries.append({
                'header': current_header,
                'sequence': ''.join(current_sequence)
            })

    return entries


def extract_chain_info(header: str) -> tuple[str, str, str]:
    """
    Извлекает из заголовка: chain_id, chain_type (H/L/other), protein_name.
    """
    # Извлекаем chain ID
    chain_match = re.search(r'Chain\s+([A-Z])', header)
    chain_id = chain_match.group(1) if chain_match else "unknown"

    # Определяем тип цепи (H = Heavy, L = Light)
    header_lower = header.lower()
    if 'heavy' in header_lower or chain_id in 'HL':
        if 'light' not in header_lower:
            chain_type = 'H'
        else:
            chain_type = 'other'
    elif 'light' in header_lower:
        chain_type = 'L'
    else:
        chain_type = 'other'

    # Извлекаем название белка (берём часть после последнего | до |Homo/Mus/synthetic)
    parts = header.split('|')
    protein_name = parts[-1].strip() if len(parts) > 1 else header[1:]
    # Убираем организмы
    protein_name = re.sub(r'\s*\([^(]+\)\s*$', '', protein_name)

    return chain_id, chain_type, protein_name


def main():
    source_dir = Path("/home/mullagaliamova/ClaudeWorkspace/PROJECTS/cdr-h3-folding/fasta-filtered")
    output_file = Path("/home/mullagaliamova/ClaudeWorkspace/PROJECTS/cdr-h3-folding/sequence_occurrence_table.tsv")

    # Структура: sequence -> список (filename, header, chain_id, chain_type)
    seq_data = defaultdict(list)

    # Читаем все FASTA файлы
    for fasta_file in sorted(source_dir.glob("*.fasta")):
        entries = parse_fasta(fasta_file)
        for entry in entries:
            chain_id, chain_type, protein_name = extract_chain_info(entry['header'])
            seq_data[entry['sequence']].append({
                'file': fasta_file.stem,
                'header': entry['header'],
                'chain_id': chain_id,
                'chain_type': chain_type
            })

    # Сортируем по количеству вхождений (убывание)
    sorted_seqs = sorted(seq_data.items(), key=lambda x: -len(x[1]))

    # Находим полные дубликаты антител (H + L вместе в тех же комплексах)
    # Группируем по парам (heavy_seq, light_seq)
    ab_pairs = defaultdict(set)
    for seq, occurrences in seq_data.items():
        types_in_seq = set(o['chain_type'] for o in occurrences)
        files_in_seq = set(o['file'] for o in occurrences)
        if 'H' in types_in_seq or 'L' in types_in_seq:
            for occ in occurrences:
                if occ['chain_type'] == 'H':
                    # Ищем light chain в том же файле
                    for other_seq, other_occ in seq_data.items():
                        for other in other_occ:
                            if other['file'] == occ['file'] and other['chain_type'] == 'L':
                                pair_key = (seq, other_seq)
                                ab_pairs[pair_key].add(occ['file'])

    # Пишем таблицу
    with open(output_file, 'w') as f:
        # Заголовок
        f.write("Chain_Header\tOccurrence_Count\tChain_ID\tChain_Type\tProtein_Name\tComplexes\tSequence\n")

        for seq, occurrences in sorted_seqs:
            # Уникальные заголовки (первый)
            first_occ = occurrences[0]
            header = first_occ['header']
            chain_id = first_occ['chain_id']
            chain_type = first_occ['chain_type']

            # Название белка
            parts = header.split('|')
            protein_name = parts[-1].strip() if len(parts) > 1 else header[1:]
            protein_name = re.sub(r'\s*\([^(]+\)\s*$', '', protein_name)

            # Список комплексов с указанием цепи
            complexes = []
            for occ in occurrences:
                complexes.append(f"{occ['file']}({occ['chain_id']})")
            complexes_str = '; '.join(complexes)

            # Количество вхождений
            count = len(occurrences)

            # Экранируем табы в последовательности (их там быть не должно, но на всякий случай)
            seq_clean = seq.replace('\t', ' ')

            f.write(f"{header}\t{count}\t{chain_id}\t{chain_type}\t{protein_name}\t{complexes_str}\t{seq_clean}\n")

    print(f"Таблица сохранена: {output_file}")
    print(f"Всего уникальных последовательностей: {len(seq_data)}")

    # Дополнительная статистика по антителам
    print("\n=== Статистика по антителам ===")

    # Считаем H и L цепи
    h_seqs = sum(1 for seq, occs in seq_data.items() if any(o['chain_type'] == 'H' for o in occs))
    l_seqs = sum(1 for seq, occs in seq_data.items() if any(o['chain_type'] == 'L' for o in occs))
    print(f"Уникальных Heavy chain последовательностей: {h_seqs}")
    print(f"Уникальных Light chain последовательностей: {l_seqs}")

    # Находим полные дубликаты антител
    full_ab_dups = [(pair, files) for pair, files in ab_pairs.items() if len(files) > 1]
    print(f"\nПолных дубликатов антител (H+L в одних и тех же комплексах): {len(full_ab_dups)}")

    if full_ab_dups:
        print("\nПримеры полных дубликатов антител:")
        for (h_seq, l_seq), files in sorted(full_ab_dups, key=lambda x: -len(x[1]))[:10]:
            print(f"  {len(files)} комплексов: {', '.join(sorted(files)[:5])}...")


if __name__ == "__main__":
    main()
