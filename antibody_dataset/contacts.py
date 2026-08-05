from __future__ import annotations

from dataclasses import dataclass, field

import gemmi


ANTIBODY_ROLES = {"heavy", "light"}
DEFAULT_CUTOFF = 5.0


@dataclass
class ContactMetrics:
    """Метрики интерфейса между двумя цепями."""

    chain_a: str
    role_a: str
    chain_b: str
    role_b: str

    atom_contacts: int = 0
    residue_contacts_a: set[str] = field(default_factory=set)
    residue_contacts_b: set[str] = field(default_factory=set)
    minimum_distance: float | None = None

    _atom_pairs: set[tuple] = field(
        default_factory=set,
        repr=False,
    )

    @property
    def residue_contacts(self) -> int:
        """Общее число уникальных контактирующих остатков."""
        return len(
            self.residue_contacts_a | self.residue_contacts_b
        )

    def add_contact(
        self,
        atom_pair: tuple,
        residue_a: str,
        residue_b: str,
        distance: float,
    ) -> None:
        """Добавляет уникальную пару контактирующих атомов."""
        if atom_pair in self._atom_pairs:
            return

        self._atom_pairs.add(atom_pair)
        self.atom_contacts += 1

        self.residue_contacts_a.add(residue_a)
        self.residue_contacts_b.add(residue_b)

        if (
            self.minimum_distance is None
            or distance < self.minimum_distance
        ):
            self.minimum_distance = distance

    def as_dict(self) -> dict:
        """Возвращает метрики в виде словаря."""
        return {
            "chain_a": self.chain_a,
            "role_a": self.role_a,
            "chain_b": self.chain_b,
            "role_b": self.role_b,
            "atom_contacts": self.atom_contacts,
            "residue_contacts_a": len(self.residue_contacts_a),
            "residue_contacts_b": len(self.residue_contacts_b),
            "residue_contacts": self.residue_contacts,
            "minimum_distance": self.minimum_distance,
        }


def residue_id(residue: gemmi.Residue) -> str:
    """Формирует идентификатор остатка внутри цепи."""
    number = residue.seqid.num
    insertion_code = residue.seqid.icode.strip()

    return f"{number}{insertion_code}"


def roles_by_chain(
    annotations: list[dict],
) -> dict[str, str]:
    """
    Создаёт отображение auth chain ID → роль.

    Одна FASTA-запись может соответствовать нескольким цепям.
    """
    result: dict[str, str] = {}

    for annotation in annotations:
        role = annotation["role"]

        for chain_id in annotation["chain_ids"]:
            result[chain_id] = role

    return result


def canonical_pair(
    chain_a: str,
    chain_b: str,
) -> tuple[str, str]:
    """Возвращает пару цепей в стабильном порядке."""
    return tuple(sorted((chain_a, chain_b)))


def calculate_contacts(
    structure: gemmi.Structure,
    annotations: list[dict],
    cutoff: float = DEFAULT_CUTOFF,
    min_distance: float = 0.1,
) -> list[ContactMetrics]:
    """
    Считает контакты антительных цепей с другими цепями.

    Учитываются:

    - heavy–light;
    - heavy/light–antigen;
    - heavy/light–other.

    Контакты между двумя неантительными цепями не считаются.
    Nanobody не считается обычной антительной цепью.
    """
    if len(structure) == 0:
        raise ValueError("Структура не содержит моделей")

    model = structure[0]
    chain_roles = roles_by_chain(annotations)

    antibody_chain_ids = {
        chain_id
        for chain_id, role in chain_roles.items()
        if role in ANTIBODY_ROLES
    }

    if not antibody_chain_ids:
        return []

    neighbor_search = gemmi.NeighborSearch(
        model,
        structure.cell,
        cutoff,
    )
    neighbor_search.populate(include_h=False)

    metrics: dict[tuple[str, str], ContactMetrics] = {}

    for source_chain in model:
        source_id = source_chain.name

        if source_id not in antibody_chain_ids:
            continue

        source_role = chain_roles[source_id]

        for source_residue_index, source_residue in enumerate(source_chain):
            for source_atom_index, source_atom in enumerate(source_residue):
                neighbors = neighbor_search.find_atoms(
                    source_atom.pos,
                    min_dist=min_distance,
                    radius=cutoff,
                )

                for neighbor in neighbors:
                    target_chain = model[neighbor.chain_idx]
                    target_id = target_chain.name

                    if target_id == source_id:
                        continue

                    if target_id not in chain_roles:
                        target_role = "unknown"
                    else:
                        target_role = chain_roles[target_id]

                    # Для пары heavy-light считаем только один раз.
                    if target_role in ANTIBODY_ROLES:
                        if source_id > target_id:
                            continue

                    target_residue = target_chain[neighbor.residue_idx]
                    target_atom = target_residue[neighbor.atom_idx]

                    pair = canonical_pair(source_id, target_id,)

                    if pair not in metrics:
                        chain_a, chain_b = pair

                        metrics[pair] = ContactMetrics(
                            chain_a=chain_a,
                            role_a=chain_roles.get(
                                chain_a,
                                "unknown",
                            ),
                            chain_b=chain_b,
                            role_b=chain_roles.get(
                                chain_b,
                                "unknown",
                            ),
                        )

                    source_atom_key = (
                        source_id,
                        source_residue_index,
                        source_atom_index,
                    )

                    target_atom_key = (
                        target_id,
                        neighbor.residue_idx,
                        neighbor.atom_idx,
                    )

                    atom_pair = tuple(
                        sorted(
                            (
                                source_atom_key,
                                target_atom_key,
                            )
                        )
                    )

                    source_residue_id = residue_id(source_residue)
                    target_residue_id = residue_id(target_residue)

                    if source_id == pair[0]:
                        residue_a = source_residue_id
                        residue_b = target_residue_id
                    else:
                        residue_a = target_residue_id
                        residue_b = source_residue_id

                    metrics[pair].add_contact(
                        atom_pair=atom_pair,
                        residue_a=residue_a,
                        residue_b=residue_b,
                        distance=source_atom.pos.dist(neighbor.pos),
                    )

    return list(metrics.values())

