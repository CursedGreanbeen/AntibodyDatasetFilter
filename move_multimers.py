#!/usr/bin/env python3
"""
Скрипт для перемещения FASTA файлов с более чем 4 цепями из fasta-filtered в multimers.
"""

import re
import shutil
from pathlib import Path


def extract_chain_ids(header: str) -> list[str]:
    """
    Извлекает auth chain IDs из заголовка FASTA.
    Обрабатывает форматы:
    - Chain A
    - Chain A [auth H]
    - Chains A, B, C
    - Chains A [auth H], B [auth L], C
    """
    multi_pattern = r'Chains?\s+([A-Z](?:\[[^\]]*\])?(?:\s*,\s*[A-Z](?:\[[^\]]*\])?)*)'
    single_pattern = r'Chain\s+([A-Z](?:\[[^\]]*\])?)'

    auth_chains = []

    multi_match = re.search(multi_pattern, header)
    if multi_match:
        chain_str = multi_match.group(1)
        entries = [c.strip() for c in chain_str.split(',')]
        for entry in entries:
            auth_match = re.search(r'\[auth\s*([A-Z])\]', entry)
            if auth_match:
                auth_chains.append(auth_match.group(1))
            else:
                match = re.match(r'^([A-Z])', entry)
                if match:
                    auth_chains.append(match.group(1))
    else:
        single_match = re.search(single_pattern, header)
        if single_match:
            entry = single_match.group(1)
            auth_match = re.search(r'\[auth\s*([A-Z])\]', entry)
            if auth_match:
                auth_chains.append(auth_match.group(1))
            else:
                match = re.match(r'^([A-Z])', entry)
                if match:
                    auth_chains.append(match.group(1))

    return auth_chains


def count_unique_chains(fasta_path: Path) -> set[str]:
    """
    Читает FASTA файл и возвращает множество уникальных chain IDs.
    """
    all_chains = set()

    with open(fasta_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                chains = extract_chain_ids(line)
                all_chains.update(chains)

    return all_chains


def main():
    source_dir = Path("/home/mullagaliamova/ClaudeWorkspace/PROJECTS/cdr-h3-folding/fasta-filtered")
    target_dir = Path("/home/mullagaliamova/ClaudeWorkspace/PROJECTS/cdr-h3-folding/multimers")
    threshold = 4  # Перемещать файлы с БОЛЕЕ чем 4 цепями

    # Создаем целевую папку если не существует
    target_dir.mkdir(exist_ok=True)

    moved_count = 0

    for fasta_file in source_dir.glob("*.fasta"):
        unique_chains = count_unique_chains(fasta_file)
        num_chains = len(unique_chains)

        if num_chains > threshold:
            target_path = target_dir / fasta_file.name
            shutil.move(str(fasta_file), str(target_path))
            print(f"Перемещен: {fasta_file.name} ({num_chains} цепей: {sorted(unique_chains)})")
            moved_count += 1

    print(f"\nИтого перемещено файлов: {moved_count}")


if __name__ == "__main__":
    main()
