from __future__ import annotations
from dataclasses import dataclass, field
from pathlib import Path
import re


# ---------------------------------------------------------------------------
# Датакласс
# ---------------------------------------------------------------------------

@dataclass
class FastaRecord:
    """Одна запись FASTA-файла.

    Несколько ключей в словаре могут указывать на один и тот же экземпляр,
    если запись охватывает несколько цепей (Chains A, B).
    """
    header: str
    sequence: str
    chain_ids: list[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Вспомогательные функции для заголовков (перенесены из старого fasta_handler)
# ---------------------------------------------------------------------------

def extract_chain_ids(header: str) -> list[str]:
    """Извлекает auth chain IDs из заголовка FASTA."""
    multi_pattern = (
        r'Chains?\s+([A-Z](?:\[[^\]]*\])?(?:\s*,\s*[A-Z](?:\[[^\]]*\])?)*)'
    )
    single_pattern = r'Chain\s+([A-Z](?:\[[^\]]*\])?)'

    auth_chains = []

    def _resolve_auth(entry: str) -> str | None:
        auth_match = re.search(r'\[auth\s*([A-Z])\]', entry)
        if auth_match:
            return auth_match.group(1)
        m = re.match(r'^([A-Z])', entry.strip())
        return m.group(1) if m else None

    multi_match = re.search(multi_pattern, header)
    if multi_match:
        for entry in multi_match.group(1).split(','):
            cid = _resolve_auth(entry)
            if cid:
                auth_chains.append(cid)
    else:
        single_match = re.search(single_pattern, header)
        if single_match:
            cid = _resolve_auth(single_match.group(1))
            if cid:
                auth_chains.append(cid)

    return auth_chains


def edit_fasta_header(header: str, chains_to_remove: set[str]) -> str:
    """Удаляет цепи из заголовка.

    Возвращает отредактированный заголовок или пустую строку,
    если все цепи были удалены (запись подлежит удалению).
    """
    _DELETED = "__DELETE_RECORD__"

    def _get_check_id(entry: str) -> str:
        auth_match = re.search(r'\[auth\s*([A-Z])\]', entry)
        if auth_match:
            return auth_match.group(1)
        m = re.match(r'^([A-Z])', entry.strip())
        return m.group(1) if m else entry

    def replace_multi(match):
        entries = [e.strip() for e in match.group(2).split(',')]
        remaining = [e for e in entries if _get_check_id(e) not in chains_to_remove]
        if not remaining:
            return _DELETED
        if len(remaining) == 1:
            return f"Chain {remaining[0]}"
        return f"Chains {', '.join(remaining)}"

    def replace_single(match):
        entry = match.group(2)
        auth_match = re.search(r'\[auth\s*([A-Z])\]', entry)
        cid = auth_match.group(1) if auth_match else re.match(r'^([A-Z])', entry).group(1)
        return _DELETED if cid in chains_to_remove else match.group(0)

    result = re.sub(
        r'(Chains)\s+([A-Z](?:\[[^\]]*\])?(?:,\s*[A-Z](?:\[[^\]]*\])?)*)',
        replace_multi,
        header,
    )
    if _DELETED in result:
        return ""

    result = re.sub(r'(Chain)\s+([A-Z](?:\[[^\]]*\])?)', replace_single, result)
    return "" if _DELETED in result else result


# ---------------------------------------------------------------------------
# Парсер и сборщик
# ---------------------------------------------------------------------------

def parse_fasta(content: str) -> dict[str, FastaRecord]:
    """Разбирает FASTA-контент в словарь {chain_id: FastaRecord}.

    Несколько ключей могут указывать на один объект FastaRecord,
    если запись охватывает несколько цепей.
    """
    records: dict[str, FastaRecord] = {}
    current_header: str | None = None
    current_lines: list[str] = []

    def _flush():
        if current_header is None:
            return
        sequence = "".join(current_lines)
        chain_ids = extract_chain_ids(current_header)
        record = FastaRecord(
            header=current_header,
            sequence=sequence,
            chain_ids=chain_ids,
        )
        if chain_ids:
            for cid in chain_ids:
                records[cid] = record
        if not chain_ids:
            print(f"[fasta_handler] Предупреждение: не удалось извлечь chain ids из заголовка: {current_header!r}")

    for line in content.splitlines():
        if line.startswith('>'):
            _flush()
            current_header = line
            current_lines = []
        else:
            stripped = line.strip()
            if stripped:
                current_lines.append(stripped)

    _flush()
    return records


def serialize_fasta(records: dict[str, FastaRecord]) -> str:
    """Сериализует словарь записей обратно в FASTA-строку.

    Дедуплицирует записи по идентичности объектов (несколько ключей →
    один объект выводится один раз).
    """
    seen: set[int] = set()
    parts: list[str] = []

    for record in records.values():
        if id(record) in seen:
            continue
        seen.add(id(record))
        parts.append(record.header)
        parts.append(record.sequence)

    return "\n".join(parts)


# ---------------------------------------------------------------------------
# Операции над словарём записей
# ---------------------------------------------------------------------------

def remove_chain(records: dict[str, FastaRecord], chain_id: str) -> None:
    """Удаляет цепь из словаря записей in-place.

    Обновляет заголовок и chain_ids записи; если цепей в записи не осталось —
    удаляет все связанные ключи.
    """
    record = records.get(chain_id)
    if record is None:
        return

    new_header = edit_fasta_header(record.header, {chain_id})

    if not new_header:
        # Удаляем все ключи, указывавшие на эту запись
        keys_to_delete = [k for k, v in records.items() if v is record]
        for k in keys_to_delete:
            del records[k]
    else:
        record.header = new_header
        record.chain_ids = [c for c in record.chain_ids if c != chain_id]
        del records[chain_id]


def update_sequence(records: dict[str, FastaRecord], chain_id: str, new_sequence: str) -> None:
    """Обновляет последовательность для цепи (только sequence, заголовок не трогаем)."""
    record = records.get(chain_id)
    if record is not None:
        record.sequence = new_sequence


# ---------------------------------------------------------------------------
# Чтение / запись файла
# ---------------------------------------------------------------------------

def read_fasta_file(path: Path) -> dict[str, FastaRecord]:
    """Читает FASTA-файл и возвращает словарь записей."""
    with open(path, 'r') as f:
        return parse_fasta(f.read())


def write_fasta_file(path: Path, records: dict[str, FastaRecord]) -> None:
    """Записывает словарь записей в FASTA-файл."""
    with open(path, 'w') as f:
        f.write(serialize_fasta(records))
