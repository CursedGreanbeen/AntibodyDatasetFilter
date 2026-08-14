from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable

from .contacts import ContactMetrics


ANTIBODY_ROLES = {"heavy", "light"}
MIN_ANTIBODY_ATOM_CONTACTS = 10
MIN_INTERFACE_ATOM_CONTACTS = 20
MIN_ANTIBODY_RESIDUES_PER_CHAIN = 5


@dataclass
class FilterResult:
    """Результат классификации одной структуры."""

    classification: str
    reasons: list[str] = field(default_factory=list)

    heavy_chains: list[str] = field(default_factory=list)
    light_chains: list[str] = field(default_factory=list)
    antibody_pairs: list[tuple[str, str]] = field(
        default_factory=list
    )

    antigen_chains: list[str] = field(default_factory=list)
    nanobody_chains: list[str] = field(default_factory=list)
    unknown_contact_chains: list[str] = field(
        default_factory=list
    )

    antibody_contacts: list[ContactMetrics] = field(
        default_factory=list
    )
    antigen_contacts: list[ContactMetrics] = field(
        default_factory=list
    )

    def as_dict(self) -> dict:
        """Преобразует результат в строковый словарь для отчёта."""
        return {
            "classification": self.classification,
            "reasons": ";".join(self.reasons),
            "heavy_chains": ",".join(self.heavy_chains),
            "light_chains": ",".join(self.light_chains),
            "antibody_pairs": ";".join(
                f"{heavy}-{light}"
                for heavy, light in self.antibody_pairs
            ),
            "antigen_chains": ",".join(self.antigen_chains),
            "nanobody_chains": ",".join(
                self.nanobody_chains
            ),
            "unknown_contact_chains": ",".join(
                self.unknown_contact_chains
            ),
        }


def _roles_by_chain(
    annotations: list[dict],
) -> dict[str, str]:
    """Создаёт отображение chain_id -> role."""
    roles: dict[str, str] = {}

    for annotation in annotations:
        for chain_id in annotation["chain_ids"]:
            roles[chain_id] = annotation["role"]

    return roles


def _annotation_by_chain(
    annotations: list[dict],
) -> dict[str, dict]:
    """Создаёт отображение chain_id -> аннотация."""
    result: dict[str, dict] = {}

    for annotation in annotations:
        for chain_id in annotation["chain_ids"]:
            result[chain_id] = annotation

    return result


def _contact_chain_ids(
    contact: ContactMetrics,
) -> tuple[str, str]:
    return contact.chain_a, contact.chain_b


def _other_chain(
    contact: ContactMetrics,
    antibody_chain_ids: set[str],
) -> str | None:
    """Возвращает вторую цепь интерфейса относительно антитела."""
    chain_a, chain_b = _contact_chain_ids(contact)

    if chain_a in antibody_chain_ids:
        return chain_b

    if chain_b in antibody_chain_ids:
        return chain_a

    return None


def _has_sufficient_contacts(
    contact: ContactMetrics,
    minimum: int,
) -> bool:
    return contact.atom_contacts >= minimum


def classify_complex(
    annotations: list[dict],
    contacts: Iterable[ContactMetrics],
    polymer_chain_ids: set[str] | None = None,
    minimum_pair_atom_contacts: int = (
        MIN_ANTIBODY_ATOM_CONTACTS
    ),
    minimum_interface_atom_contacts: int = (
        MIN_INTERFACE_ATOM_CONTACTS
    ),
) -> FilterResult:
    """
    Классифицирует один комплекс.

    `annotations` должны быть результатом annotate_fasta().
    `contacts` должны быть результатом calculate_contacts().
    """
    if polymer_chain_ids is None:
        polymer_chain_ids = set()
        
    contacts = list(contacts)
    roles = _roles_by_chain(annotations)
    annotation_by_chain = _annotation_by_chain(annotations)

    heavy_chains = sorted(
        chain_id
        for chain_id, role in roles.items()
        if role == "heavy"
    )

    light_chains = sorted(
        chain_id
        for chain_id, role in roles.items()
        if role == "light"
    )

    nanobody_chains = sorted(
        chain_id
        for chain_id, role in roles.items()
        if role == "nanobody"
    )

    result = FilterResult(
        classification="reject",
        heavy_chains=heavy_chains,
        light_chains=light_chains,
        nanobody_chains=nanobody_chains,
    )

    if not heavy_chains:
        result.reasons.append("no_heavy_chain")

    if not light_chains:
        result.reasons.append("no_light_chain")

    if not heavy_chains or not light_chains:
        return result

    antibody_pairs: list[tuple[str, str]] = []
    antibody_contacts: list[ContactMetrics] = []

    selected_pairs = find_antibody_pairs(
    contacts=contacts,
    minimum_atom_contacts=minimum_pair_atom_contacts,
    )

    antibody_pairs = [
        (heavy_id, light_id)
        for heavy_id, light_id, _ in selected_pairs
    ]

    antibody_contacts = [
        contact
        for _, _, contact in selected_pairs
    ]

    result.antibody_pairs = antibody_pairs
    result.antibody_contacts = antibody_contacts

    # Убираем возможные повторы пар.
    antibody_pairs = sorted(set(antibody_pairs))

    result.antibody_pairs = antibody_pairs
    result.antibody_contacts = antibody_contacts

    if not antibody_pairs:
        result.reasons.append("no_heavy_light_interface")
        return result

    if len(antibody_pairs) > 1:
        result.reasons.append("multiple_antibody_pairs")
        return result

    if nanobody_chains:
        result.reasons.append("nanobody_present")

    heavy, light = antibody_pairs[0]
    antibody_chain_ids = {heavy, light}

    antigen_chains: set[str] = set()
    unknown_contact_chains: set[str] = set()
    antigen_contacts: list[ContactMetrics] = []

    for contact in contacts:
        if not _has_sufficient_contacts(
            contact,
            minimum_pair_atom_contacts,
        ):
            continue

        other_chain = _other_chain(
            contact,
            antibody_chain_ids,
        )

        if other_chain is None:
            continue

        if other_chain in antibody_chain_ids:
            continue

        other_role = roles.get(other_chain)

        if other_role == "other":
            antigen_chains.add(other_chain)
            antigen_contacts.append(contact)

        elif other_role == "nanobody":
            # Уже отражено в nanobody_present.
            continue

        else:
            # Включая цепи, отсутствующие в FASTA.
            unknown_contact_chains.add(other_chain)

    result.antigen_chains = sorted(antigen_chains)
    result.unknown_contact_chains = sorted(
        unknown_contact_chains
    )
    result.antigen_contacts = antigen_contacts

    other_polymer_chain_ids = (
        polymer_chain_ids - antibody_chain_ids
    )
    
    if not antigen_chains:
        if unknown_contact_chains:
            pass
        elif not other_polymer_chain_ids:
            result.reasons.append("no_other_polymer_chain")
        else:
            result.reasons.append("no_antigen_interface")

    if len(antigen_chains) > 1:
        result.reasons.append("multiple_antigen_chains")

    if unknown_contact_chains:
        result.reasons.append("unknown_contact_chain")

    if (
        len(antigen_chains) == 1
        and not nanobody_chains
        and not unknown_contact_chains
        and not result.reasons
    ):
        result.classification = "strict"
        return result

    if (
        "no_antigen_interface" in result.reasons
        or "no_other_polymer_chain" in result.reasons
    ):
        result.classification = "reject"
    else:
        result.classification = "candidate"

    return result


def _is_heavy_light_contact(
    contact: ContactMetrics,
) -> bool:
    """Проверяет, является ли контакт парой heavy–light."""
    return {
        contact.role_a,
        contact.role_b,
    } == {"heavy", "light"}


def _contact_chain_ids(
    contact: ContactMetrics,
) -> tuple[str, str]:
    """Возвращает heavy и light в фиксированном порядке."""
    if contact.role_a == "heavy":
        return contact.chain_a, contact.chain_b

    return contact.chain_b, contact.chain_a


def _pair_score(
    contact: ContactMetrics,
) -> tuple[int, int, int]:
    """
    Ключ качества интерфейса.

    Сначала учитывается число атомных контактов,
    затем число контактирующих остатков каждой цепи
    и общее число контактирующих остатков.
    """
    return (
        contact.atom_contacts,
        min(
            len(contact.residue_contacts_a),
            len(contact.residue_contacts_b),
        ),
        contact.residue_contacts,
    )


def find_antibody_pairs(
    contacts: list[ContactMetrics],
    minimum_atom_contacts: int = MIN_ANTIBODY_ATOM_CONTACTS,
    minimum_residues_per_chain: int = (
        MIN_ANTIBODY_RESIDUES_PER_CHAIN
    ),
) -> list[tuple[str, str, ContactMetrics]]:
    """
    Находит взаимно лучшие пары heavy–light.

    Возвращает список кортежей:
        (heavy_chain_id, light_chain_id, contact_metrics)
    """
    candidates: list[
        tuple[str, str, ContactMetrics]
    ] = []

    for contact in contacts:
        if not _is_heavy_light_contact(contact):
            continue

        if contact.atom_contacts < minimum_atom_contacts:
            continue

        residues_a = len(contact.residue_contacts_a)
        residues_b = len(contact.residue_contacts_b)

        if min(residues_a, residues_b) < minimum_residues_per_chain:
            continue

        heavy_id, light_id = _contact_chain_ids(contact)

        candidates.append(
            (heavy_id, light_id, contact)
        )

    if not candidates:
        return []

    by_heavy: dict[
        str,
        list[tuple[str, str, ContactMetrics]],
    ] = {}

    by_light: dict[
        str,
        list[tuple[str, str, ContactMetrics]],
    ] = {}

    for candidate in candidates:
        heavy_id, light_id, _ = candidate

        by_heavy.setdefault(
            heavy_id,
            [],
        ).append(candidate)

        by_light.setdefault(
            light_id,
            [],
        ).append(candidate)

    best_light_for_heavy: dict[str, str] = {}

    for heavy_id, options in by_heavy.items():
        best = max(
            options,
            key=lambda item: _pair_score(item[2]),
        )
        best_light_for_heavy[heavy_id] = best[1]

    best_heavy_for_light: dict[str, str] = {}

    for light_id, options in by_light.items():
        best = max(
            options,
            key=lambda item: _pair_score(item[2]),
        )
        best_heavy_for_light[light_id] = best[0]

    pairs: list[tuple[str, str, ContactMetrics]] = []

    for heavy_id, light_id, contact in candidates:
        if best_light_for_heavy.get(heavy_id) != light_id:
            continue

        if best_heavy_for_light.get(light_id) != heavy_id:
            continue

        pairs.append(
            (heavy_id, light_id, contact)
        )

    return sorted(
        pairs,
        key=lambda item: (
            item[0],
            item[1],
        ),
    )
