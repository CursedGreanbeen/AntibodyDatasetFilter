from __future__ import annotations

from pathlib import Path
from typing import Any

from .anarci_annotator import chain_role, run_anarci
from .fasta_manager import FastaRecord, read_fasta_file


NANOBODY_WORDS = (
    "nanobody",
    "nanoantibody",
    "vhh",
    "camelid single-domain",
    "single domain antibody",
)


def is_nanobody_header(header: str) -> bool:
    """Определяет нанотело по описанию FASTA."""
    text = header.lower()
    return any(word in text for word in NANOBODY_WORDS)


def annotate_record(
    record: FastaRecord,
) -> dict[str, Any]:
    """Запускает ANARCI для одной FASTA-записи."""
    chain_type, start, end, status = run_anarci(
        record.sequence
    )

    nanobody = is_nanobody_header(record.header)

    if nanobody:
        role = "nanobody"
    else:
        role = chain_role(chain_type)

    return {
        "header": record.header,
        "sequence": record.sequence,
        "chain_ids": record.chain_ids,
        "chain_type": chain_type,
        "role": role,
        "start": start,
        "end": end,
        "status": status,
        "is_nanobody": nanobody,
    }


def annotate_fasta(
    fasta_path: Path,
) -> list[dict[str, Any]]:
    """
    Аннотирует все уникальные FASTA-записи файла.

    Если одна FASTA-запись относится к нескольким цепям,
    она анализируется ANARCI только один раз.
    """
    records = read_fasta_file(fasta_path)

    unique_records: list[FastaRecord] = []
    seen_records: set[int] = set()

    for record in records.values():
        record_id = id(record)

        if record_id in seen_records:
            continue

        seen_records.add(record_id)
        unique_records.append(record)

    annotations = []

    for record in unique_records:
        annotations.append(
            annotate_record(record)
        )

    return annotations


def expand_annotations_by_chain(
    annotations: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """
    Делает отдельную запись для каждой цепи.

    Это удобно для дальнейшего объединения с chains.py
    и расчёта контактов.
    """
    expanded: list[dict[str, Any]] = []

    for annotation in annotations:
        for chain_id in annotation["chain_ids"]:
            row = dict(annotation)
            row["chain_id"] = chain_id
            expanded.append(row)

    return expanded


