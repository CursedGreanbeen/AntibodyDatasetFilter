from dataclasses import dataclass
from pathlib import Path

import gemmi


@dataclass
class ResidueInfo:
    name: str
    number: int
    insertion_code: str
    one_letter_code: str


@dataclass
class ChainInfo:
    chain_id: str
    residue_count: int
    first_residue: ResidueInfo
    last_residue: ResidueInfo


def read_structure(cif_path: Path) -> gemmi.Structure:
    """Загружает CIF/mmCIF-файл."""
    structure = gemmi.read_structure(str(cif_path))

    if len(structure) == 0:
        raise ValueError(f"В файле нет моделей: {cif_path}")

    return structure


def get_model_chains(
    structure: gemmi.Structure,
) -> list[ChainInfo]:
    """Возвращает информацию обо всех непустых полимерных цепях."""
    if len(structure) == 0:
        raise ValueError("Структура не содержит моделей")

    chains: list[ChainInfo] = []

    for chain in structure[0]:
        chain_info = get_chain_info(chain)

        if chain_info is not None:
            chains.append(chain_info)

    return chains


def get_polymer_residues(
    chain: gemmi.Chain,
) -> list[gemmi.Residue]:
    """
    Возвращает содержащие атомы аминокислотные остатки.

    Сначала используется Gemmi get_polymer(). Если CIF не содержит
    entity-информации и get_polymer() возвращает пустой результат,
    выполняется fallback по всем остаткам цепи.
    """
    polymer_residues = [
        residue
        for residue in chain.get_polymer()
        if len(residue) > 0
    ]

    if polymer_residues:
        return polymer_residues

    fallback_residues = []

    for residue in chain:
        if len(residue) == 0:
            continue

        residue_info = gemmi.find_tabulated_residue(
            residue.name
        )

        if residue_info.is_amino_acid():
            fallback_residues.append(residue)

    return fallback_residues


def residue_to_info(residue: gemmi.Residue) -> ResidueInfo:
    """Преобразует объект Gemmi Residue в удобную структуру данных."""
    residue_info = gemmi.find_tabulated_residue(residue.name)

    if residue_info is None:
        one_letter_code = "X"
    else:
        one_letter_code = residue_info.one_letter_code

    return ResidueInfo(
        name=residue.name,
        number=residue.seqid.num,
        insertion_code=residue.seqid.icode.strip(),
        one_letter_code=one_letter_code,
    )


def get_chain_info(chain: gemmi.Chain) -> ChainInfo | None:
    """Возвращает информацию о цепи или None для пустой цепи."""
    residues = get_polymer_residues(chain)

    if not residues:
        return None

    return ChainInfo(
        chain_id=chain.name,
        residue_count=len(residues),
        first_residue=residue_to_info(residues[0]),
        last_residue=residue_to_info(residues[-1]),
    )


def get_structure_chains(
    cif_path: Path,
) -> list[ChainInfo]:
    """Загружает CIF и возвращает информацию о полимерных цепях."""
    structure = read_structure(cif_path)

    return get_model_chains(structure)
