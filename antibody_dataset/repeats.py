from __future__ import annotations
from dataclasses import dataclass

from .contacts import ContactMetrics, roles_by_chain


@dataclass(frozen=True)
class SequenceGroup:
    """
    Группа структурных цепей с одинаковой ролью
    и одинаковой аминокислотной последовательностью.
    """

    role: str
    sequence: str
    chain_ids: tuple[str, ...]


@dataclass(frozen=True)
class RepeatUnit:
    """
    Одна потенциальная копия комплекса.

    heavy_chain и light_chain — цепи одного антитела.
    antigen_chains — цепи полимерного партнёра,
    контактирующие с этой парой.
    """

    heavy_chain: str
    light_chain: str
    antigen_chains: tuple[str, ...] = ()


@dataclass(frozen=True)
class RepeatAnalysis:
    """
    Результат поиска повторяющихся комплексов.

    units:
        все найденные полноценные H–L–Ag-юниты.

    repeated_units:
        группа наиболее похожих друг на друга юнитов.

    copy_count:
        число копий в repeated_units.

    is_repeated_copy:
        True, если найдено минимум две одинаковые копии.

    representative_unit:
        один юнит, выбранный как представитель группы.

    reasons:
        пояснения для отчёта.
    """

    units: tuple[RepeatUnit, ...]
    repeated_units: tuple[RepeatUnit, ...]
    copy_count: int
    is_repeated_copy: bool
    representative_unit: RepeatUnit | None
    reasons: tuple[str, ...] = ()


def normalize_sequence(sequence: str) -> str:
    """
    Удаляет пробелы и переводы строк,
    чтобы последовательности можно было корректно сравнивать.
    """
    return "".join(sequence.split()).upper()


REPEATABLE_ROLES = {
    "heavy",
    "light",
    "other",
}


def group_chains_by_sequence(
    annotations: list[dict],
) -> list[SequenceGroup]:
    """Группирует цепи по роли и аминокислотной последовательности."""
    grouped: dict[tuple[str, str], set[str]] = {}

    for annotation in annotations:
        role = annotation["role"]

        if role not in REPEATABLE_ROLES:
            continue

        sequence = normalize_sequence(
            annotation["sequence"]
        )

        if not sequence:
            continue

        key = (role, sequence)

        grouped.setdefault(key, set()).update(
            annotation["chain_ids"]
        )

    return [
        SequenceGroup(
            role=role,
            sequence=sequence,
            chain_ids=tuple(sorted(chain_ids)),
        )
        for (role, sequence), chain_ids in sorted(
            grouped.items()
        )
    ]


def index_groups_by_chain(
    groups: list[SequenceGroup],
) -> dict[str, SequenceGroup]:
    """
    Создаёт быстрый поиск группы по ID цепи.

    Было:
        SequenceGroup → chain_ids

    Станет:
        chain_id → SequenceGroup
    """
    index: dict[str, SequenceGroup] = {}

    for group in groups:
        for chain_id in group.chain_ids:
            index[chain_id] = group

    return index


def build_repeat_unit(
    heavy_chain: str,
    light_chain: str,
    contacts: list[ContactMetrics],
    annotations: list[dict],
    polymer_chain_ids: set[str],
) -> RepeatUnit | None:
    """
    Создаёт RepeatUnit для одной пары heavy–light.

    В antigen_chains включаются только полимерные цепи,
    которые имеют роль "other" и контактируют с heavy
    или light данной пары.
    """
    # Получаем роль каждой цепи из результатов ANARCI.
    roles_by_chain: dict[str, str] = {}

    for annotation in annotations:
        for chain_id in annotation["chain_ids"]:
            roles_by_chain[chain_id] = annotation["role"]

    antibody_chains = {
        heavy_chain,
        light_chain,
    }

    antigen_chains: set[str] = set()

    for contact in contacts:
        contact_chain_ids = {
            contact.chain_a,
            contact.chain_b,
        }

        # Контакт не относится к данной паре антитела.
        if not contact_chain_ids & antibody_chains:
            continue

        # Убираем heavy/light этой пары.
        partner_chain_ids = (
            contact_chain_ids - antibody_chains
        )

        for chain_id in partner_chain_ids:
            # В RepeatUnit включаем только полимерные цепи.
            if chain_id not in polymer_chain_ids:
                continue

            # Heavy/light других антител и nanobody
            # не считаются антигеном.
            if roles_by_chain.get(chain_id) != "other":
                continue

            antigen_chains.add(chain_id)

    # Без антигенного партнёра полноценный RepeatUnit
    # построить нельзя.
    if not antigen_chains:
        return None

    return RepeatUnit(
        heavy_chain=heavy_chain,
        light_chain=light_chain,
        antigen_chains=tuple(
            sorted(antigen_chains)
        ),
    )


def build_repeat_units(
    antibody_pairs: list[tuple[str, str, ContactMetrics]],
    contacts: list[ContactMetrics],
    annotations: list[dict],
    polymer_chain_ids: set[str],
) -> list[RepeatUnit]:
    """
    Строит RepeatUnit для каждой выбранной пары heavy–light.

    antibody_pairs имеет формат:
        (heavy_chain, light_chain, ContactMetrics)

    Контактный объект в самой паре здесь не используется:
    он нужен был на этапе выбора пары. Для построения юнита
    используются все контакты, чтобы найти антигенные цепи.
    """
    units: list[RepeatUnit] = []

    for heavy_chain, light_chain, _ in antibody_pairs:
        unit = build_repeat_unit(
            heavy_chain=heavy_chain,
            light_chain=light_chain,
            contacts=contacts,
            annotations=annotations,
            polymer_chain_ids=polymer_chain_ids,
        )

        # Пара без определённого антигенного партнёра
        # не считается полноценным юнитом.
        if unit is None:
            continue

        units.append(unit)

    return units


def _unit_signature(
    unit: RepeatUnit,
    groups_by_chain: dict[str, SequenceGroup],
) -> tuple[str, str, tuple[str, ...]] | None:
    """
    Создаёт последовательностную сигнатуру одного юнита.

    Сигнатура содержит:
        - последовательность heavy chain;
        - последовательность light chain;
        - отсортированный набор последовательностей
          антигенных цепей.

    None означает, что для одной из цепей
    не удалось найти SequenceGroup.
    """
    heavy_group = groups_by_chain.get(unit.heavy_chain)
    light_group = groups_by_chain.get(unit.light_chain)

    if heavy_group is None or light_group is None:
        return None

    antigen_sequences: list[str] = []

    for chain_id in unit.antigen_chains:
        antigen_group = groups_by_chain.get(chain_id)

        if antigen_group is None:
            return None

        antigen_sequences.append(
            antigen_group.sequence
        )

    return (
        heavy_group.sequence,
        light_group.sequence,
        tuple(sorted(antigen_sequences)),
    )


def _group_units_by_signature(
    units: list[RepeatUnit],
    groups_by_chain: dict[str, SequenceGroup],
) -> list[tuple[RepeatUnit, ...]]:
    """
    Объединяет RepeatUnit с одинаковым составом последовательностей.

    Юниты с одинаковой сигнатурой считаются копиями
    одного и того же комплекса.
    """
    grouped: dict[
        tuple[str, str, tuple[str, ...]],
        list[RepeatUnit],
    ] = {}

    for unit in units:
        signature = _unit_signature(
            unit=unit,
            groups_by_chain=groups_by_chain,
        )

        # Если для цепи нет SequenceGroup, юнит нельзя
        # надёжно сравнить с другими.
        if signature is None:
            continue

        grouped.setdefault(signature, []).append(unit)

    return [
        tuple(unit_group)
        for unit_group in grouped.values()
    ]


def detect_repeated_units(
    annotations: list[dict],
    contacts: list[ContactMetrics],
    antibody_pairs: list[
        tuple[str, str, ContactMetrics]
    ],
    polymer_chain_ids: set[str],
) -> RepeatAnalysis:
    """
    Находит повторяющиеся копии одного H–L–Ag-комплекса.
    """
    sequence_groups = group_chains_by_sequence(
        annotations
    )

    groups_by_chain = index_groups_by_chain(
        sequence_groups
    )

    units = build_repeat_units(
        antibody_pairs=antibody_pairs,
        contacts=contacts,
        annotations=annotations,
        polymer_chain_ids=polymer_chain_ids,
    )

    if not units:
        return RepeatAnalysis(
            units=(),
            repeated_units=(),
            copy_count=0,
            is_repeated_copy=False,
            representative_unit=None,
            reasons=("no_complete_repeat_unit",),
        )

    unit_groups = _group_units_by_signature(
        units=units,
        groups_by_chain=groups_by_chain,
    )

    if not unit_groups:
        return RepeatAnalysis(
            units=tuple(units),
            repeated_units=(),
            copy_count=0,
            is_repeated_copy=False,
            representative_unit=None,
            reasons=("units_not_comparable",),
        )

    # Выбираем наиболее крупную группу одинаковых юнитов.
    largest_group = max(
        unit_groups,
        key=len,
    )

    copy_count = len(largest_group)
    is_repeated_copy = copy_count >= 2

    reasons: list[str] = []

    if is_repeated_copy:
        reasons.append("repeated_identical_copies")

    if len(unit_groups) > 1:
        reasons.append("multiple_unit_types")

    return RepeatAnalysis(
        units=tuple(units),
        repeated_units=tuple(largest_group),
        copy_count=copy_count,
        is_repeated_copy=is_repeated_copy,
        representative_unit=largest_group[0],
        reasons=tuple(reasons),
    )
