#!/usr/bin/env python3
"""
Analyze non-protein components at antibody-antigen interface.

Examines CIF files for ligands, ions, carbohydrates, nucleic acids,
and other hetero-components within 6 Angstroms of protein chains.

Usage:
    python scripts/analyze_interface_components.py

Output:
    TSV report with columns:
    pdb_code, label_chain_id, resname, resnum, dist_to_antibody, dist_to_antigen,
    nearest_target, category, decision, decision_reason
"""

import csv
import math
import gemmi
import re
from pathlib import Path
from dataclasses import dataclass
from typing import Set, Tuple, Optional


# Ручные маппинги для специфичных структур (label_asym_id из CIF)
CHAIN_GROUPS = {
    "8x2l": {
        "antibody": {"C", "D"},
        "antigen": {"B"},
    },
}


def parse_fasta_for_chains(fasta_dir: Path, pdb_code: str) -> Tuple[Set[str], Set[str]]:
    """
    Parse FASTA file to identify antibody (H/L) and antigen chains.

    FASTA headers have multiple formats:

    Format 1 (single chain with auth):
    >8GAT_2|Chain B[auth M]|Fab 1G01, heavy chain|Homo sapiens (9606)
    >8CBF_2|Chain B|Beta-49 light chain|Homo sapiens (9606)

    Format 2 (multiple chains, no auth):
    >10OP_8|Chains T, V, X|294S Light Chain|Homo sapiens (9606)
    >10OP_2|Chains B, D, F|Envelope glycoprotein 2|Sudan ebolavirus

    Format 3 (multiple chains with shared auth):
    >8G8A_2|Chains B, D[auth L]|DH1317.8 light chain|Homo sapiens (9606)
    >8FHY_1|Chains A, D, G[auth I]|Spike protein S1|SARS-CoV-2

    Key: auth chain ID in brackets [auth H] or [auth L] = antibody
         Or description contains "heavy chain" or "light chain"
    Returns: (antibody_chains, antigen_chains) as label_asym_id
    """
    antibody_chains: Set[str] = set()
    antigen_chains: Set[str] = set()

    # Try multiple possible FASTA directories
    fasta_patterns = [
        fasta_dir / f"{pdb_code}.fasta",
        fasta_dir / f"{pdb_code}_cropped.fasta",
    ]

    for fasta_path in fasta_patterns:
        if not fasta_path.exists():
            continue

        with open(fasta_path) as f:
            for line in f:
                line = line.strip()
                if not line.startswith('>'):
                    continue

                description = ""

                # Format 1: Chain X[auth Y] or Chain X
                match1 = re.match(r'>\S+\|Chain\s+(\S+)(?:\[auth\s+(\w+)\])?\|(.+)', line)
                # Format 2/3: Chains X, Y, Z or Chains X, Y[auth Z]
                match2 = re.match(r'>\S+\|Chains?\s+(.+?)\|(.+)', line)

                if match1:
                    label_chain = match1.group(1)
                    auth_chain = match1.group(2)
                    description = match1.group(3).lower()

                    # Remove any trailing bracket from label (e.g., "B[auth M]" -> "B")
                    label_chain = re.sub(r'\[.*\]', '', label_chain).strip()

                    # Antibody: auth H (heavy) or auth L (light) in brackets
                    if auth_chain and auth_chain.upper() in ('H', 'L'):
                        antibody_chains.add(label_chain)
                    elif 'heavy chain' in description or 'light chain' in description:
                        antibody_chains.add(label_chain)
                    else:
                        antigen_chains.add(label_chain)

                elif match2:
                    chains_part = match2.group(1)
                    description = match2.group(2).lower()

                    # Check for shared auth: "B, D[auth L]" -> chains=['B', 'D'], auth='L'
                    auth_match = re.search(r'\[auth\s+(\w+)\]', chains_part)
                    shared_auth = auth_match.group(1).upper() if auth_match else None

                    # Extract chain IDs (remove auth brackets)
                    chains_str = re.sub(r'\[auth\s+\w+\]', '', chains_part)
                    chain_ids = [c.strip() for c in chains_str.split(',')]

                    # Antibody if:
                    # - shared auth is H or L, OR
                    # - description says "heavy chain" or "light chain"
                    is_antibody = (shared_auth in ('H', 'L') or
                                   'heavy chain' in description or
                                   'light chain' in description)

                    for chain_id in chain_ids:
                        if is_antibody:
                            antibody_chains.add(chain_id)
                        else:
                            antigen_chains.add(chain_id)

        if antibody_chains or antigen_chains:
            break

    return antibody_chains, antigen_chains


STANDARD_AA = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
    "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
    "THR", "TRP", "TYR", "VAL",
}

# Остатки, которые являются модификациями аминокислот
MODIFIED_AA = {
    "MSE", "SEC", "PYL", "CSO", "TPO", "PTR", "CME",
    "HYP", "HID", "HIE", "HIP", "ASH", "GLH", "LYN",
}

WATER = {"HOH", "WAT", "DOD"}

METAL_IONS = {
    "LI", "NA", "K", "RB", "CS",
    "MG", "CA", "SR", "BA",
    "MN", "FE", "FE2", "FE3", "CO", "NI",
    "CU", "CU1", "CU2", "ZN", "CD", "HG",
}

SIMPLE_IONS = {
    "CL", "BR", "IOD", "I",
    "F", "NO3", "SCN",
    "SO4", "PO4", "CO3",
}

CARBOHYDRATES = {
    "NAG", "NDG", "BMA", "MAN", "GAL", "GLC",
    "FUC", "SIA", "NEU", "NAN", "API",
    "XYS", "BGC", "A2G", "G6D",
}

NUCLEIC_ACIDS = {
    "A", "C", "G", "U",
    "DA", "DC", "DG", "DT",
    "DI", "DU",
}

COFACTORS = {
    "ATP", "ADP", "AMP",
    "GTP", "GDP", "GMP",
    "NAD", "NADH", "NAP",
    "FAD", "FMN", "COA",
    "PLP", "TPP", "HEM",
    "SAM", "SAH", "GSH",
}

CRYSTALLIZATION_ADDITIVES = {
    "EDO", "GOL", "MPD", "DMS",
    "PEG", "PG4", "PGE",
    "ACT", "FMT", "MES",
    "HEP", "TRS", "TRIS",
    "BIS", "CAP", "PIP",
    "CAC", "IMD",
    "SO4", "PO4",
}


@dataclass
class Component:
    pdb_code: str
    label_chain_id: str
    resname: str
    resnum: str
    min_dist_to_antibody: float | None
    min_dist_to_antigen: float | None
    nearest_target: str
    category: str
    decision: str
    decision_reason: str


def atom_xyz(atom):
    """Return atom coordinates as a tuple."""
    return atom.pos.x, atom.pos.y, atom.pos.z


def squared_distance(a, b):
    return (
        (a[0] - b[0]) ** 2
        + (a[1] - b[1]) ** 2
        + (a[2] - b[2]) ** 2
    )


def residue_number(residue) -> str:
    """
    Preserve residue number and insertion code.
    Examples: 25, 25A.
    """
    number = str(residue.seqid.num)
    insertion_code = str(residue.seqid.icode).strip()

    if insertion_code and insertion_code not in {".", "?"}:
        return f"{number}{insertion_code}"

    return number


def is_protein_like(residue) -> bool:
    """
    Identify standard and modified amino-acid residues.
    """
    name = residue.name.upper()

    if name in STANDARD_AA:
        return True

    if name in MODIFIED_AA:
        return True

    return False


def classify_component(resname: str) -> str:
    """
    Classification is heuristic and should be manually reviewed
    for unknown or chemically ambiguous components.
    """
    name = resname.upper()

    if name in WATER:
        return "water"

    if name in MODIFIED_AA:
        return "modified_amino_acid"

    if name in METAL_IONS:
        return "metal_ion"

    if name in CARBOHYDRATES:
        return "carbohydrate"

    if name in NUCLEIC_ACIDS:
        return "nucleic_acid"

    if name in COFACTORS:
        return "cofactor"

    if name in CRYSTALLIZATION_ADDITIVES:
        return "buffer_or_crystallization_additive"

    if name in SIMPLE_IONS:
        return "simple_ion"

    return "unknown_ligand"


def decision_for_component(category: str):
    """
    Conservative AF3 decision.
    'review' means that the script cannot decide reliably.
    """
    if category in {
        "water",
        "buffer_or_crystallization_additive",
        "simple_ion",
    }:
        return "remove", "usually crystallization/solvent component"

    if category in {
        "carbohydrate",
        "nucleic_acid",
        "cofactor",
        "unknown_ligand",
    }:
        return "review", "possible biological relevance"

    if category == "metal_ion":
        return "review", "may be structurally or biologically important"

    if category == "modified_amino_acid":
        return "review", "check whether it is part of the protein chain"

    return "review", "manual inspection required"


def get_atoms_from_chains(model, chain_ids):
    """
    Collect atoms from protein-like residues in selected chains only.
    """
    atoms = []

    for chain in model:
        if chain.name not in chain_ids:
            continue

        for residue in chain:
            if not is_protein_like(residue):
                continue

            for atom in residue:
                atoms.append(atom_xyz(atom))

    return atoms


def minimum_distance(component_atoms, target_atoms, cutoff):
    """
    Return minimum atom-atom distance, or None if no atom is within cutoff.
    """
    if not component_atoms or not target_atoms:
        return None

    cutoff_sq = cutoff ** 2
    min_dist_sq = float("inf")

    for component_atom in component_atoms:
        for target_atom in target_atoms:
            dist_sq = squared_distance(component_atom, target_atom)

            if dist_sq < min_dist_sq:
                min_dist_sq = dist_sq

            # This is only an optimization; exact minimum is still retained.
            if dist_sq == 0:
                return 0.0

    if min_dist_sq <= cutoff_sq:
        return math.sqrt(min_dist_sq)

    return None


def analyze_cif_file(cif_path: Path, cutoff: float = 6.0,
                     fasta_dir: Optional[Path] = None) -> Optional[list]:
    """
    Analyze a CIF file for interface components.

    Returns:
        - List of Components if processed (even if empty)
        - None if skipped due to missing chain info
    """
    pdb_code = cif_path.stem.lower().replace("_cropped", "")

    # Try manual mapping first, then FASTA
    if pdb_code in CHAIN_GROUPS:
        groups = CHAIN_GROUPS[pdb_code]
        antibody_chains = groups["antibody"]
        antigen_chains = groups["antigen"]
    elif fasta_dir:
        antibody_chains, antigen_chains = parse_fasta_for_chains(fasta_dir, pdb_code)
        if not antibody_chains and not antigen_chains:
            return None  # Skipped - no chain info
    else:
        return None  # Skipped - no chain info

    structure = gemmi.read_structure(str(cif_path))
    model = structure[0]

    antibody_atoms = get_atoms_from_chains(model, antibody_chains)
    antigen_atoms = get_atoms_from_chains(model, antigen_chains)

    existing_chains = {chain.name for chain in model}

    missing_ab = antibody_chains - existing_chains
    missing_ag = antigen_chains - existing_chains

    if missing_ab:
        print(f"WARNING {pdb_code}: missing antibody chains {missing_ab}")

    if missing_ag:
        print(f"WARNING {pdb_code}: missing antigen chains {missing_ag}")

    results = []

    for chain in model:
        for residue in chain:
            resname = residue.name.upper()

            # Белковые остатки не являются отдельными небелковыми компонентами.
            if is_protein_like(residue):
                continue

            # Воду обычно не передают в AF3.
            if resname in WATER:
                continue

            component_atoms = [
                atom_xyz(atom)
                for atom in residue
            ]

            dist_ab = minimum_distance(
                component_atoms,
                antibody_atoms,
                cutoff,
            )

            dist_ag = minimum_distance(
                component_atoms,
                antigen_atoms,
                cutoff,
            )

            # Компонент интересен, если близок хотя бы к одной целевой группе.
            if dist_ab is None and dist_ag is None:
                continue

            if dist_ab is not None and dist_ag is not None:
                if dist_ab <= dist_ag:
                    nearest_target = "antibody"
                else:
                    nearest_target = "antigen"
            elif dist_ab is not None:
                nearest_target = "antibody"
            else:
                nearest_target = "antigen"

            category = classify_component(resname)
            decision, reason = decision_for_component(category)

            results.append(
                Component(
                    pdb_code=pdb_code,
                    label_chain_id=chain.name,
                    resname=resname,
                    resnum=residue_number(residue),
                    min_dist_to_antibody=dist_ab,
                    min_dist_to_antigen=dist_ag,
                    nearest_target=nearest_target,
                    category=category,
                    decision=decision,
                    decision_reason=reason,
                )
            )

    return results


def format_distance(value):
    return "" if value is None else f"{value:.2f}"


def main():
    cif_dir = Path(
        "/home/mullagaliamova/ClaudeWorkspace/"
        "PROJECTS/cdr-h3-folding/CIFs-filtered-new"
    )

    fasta_dir = Path(
        "/home/mullagaliamova/ClaudeWorkspace/PROJECTS/cdr-h3-folding/fasta-new-cropped"
    )

    output_file = cif_dir.parent / "interface_components_report.tsv"
    issues_file = cif_dir.parent / "processing_issues.tsv"

    cif_files = sorted(cif_dir.glob("*.cif"))

    print("=" * 70)
    print("INTERFACE COMPONENT ANALYSIS")
    print("=" * 70)
    print(f"Input directory: {cif_dir}")
    print(f"FASTA directory: {fasta_dir}")
    print(f"Found CIF files: {len(cif_files)}")
    print("Cutoff: 6.0 A")
    print()

    all_components = []
    issues = []  # Track skipped/problematic PDBs

    for cif_path in cif_files:
        pdb_code = cif_path.stem.lower().replace("_cropped", "")
        print(f"Processing {cif_path.name}...", end=" ")

        try:
            components = analyze_cif_file(cif_path, cutoff=6.0, fasta_dir=fasta_dir)
            if components is None:
                # Skipped due to missing chain info
                issues.append((pdb_code, "no_chain_info", "No chain mapping found in FASTA"))
                print("skipped (no chain info)")
            elif not components:
                # Processed but no components found
                issues.append((pdb_code, "no_components", "No interface components within cutoff"))
                print("found 0 components")
            else:
                all_components.extend(components)
                print(f"found {len(components)} components")

        except Exception as error:
            issues.append((pdb_code, "error", str(error)))
            print(f"ERROR: {error}")

    fields = [
        "pdb_code",
        "label_chain_id",
        "resname",
        "resnum",
        "min_dist_to_antibody",
        "min_dist_to_antigen",
        "nearest_target",
        "category",
        "decision",
        "decision_reason",
    ]

    with open(output_file, "w", newline="") as file:
        writer = csv.writer(file, delimiter="\t")
        writer.writerow(fields)

        for component in all_components:
            writer.writerow([
                component.pdb_code,
                component.label_chain_id,
                component.resname,
                component.resnum,
                format_distance(component.min_dist_to_antibody),
                format_distance(component.min_dist_to_antigen),
                component.nearest_target,
                component.category,
                component.decision,
                component.decision_reason,
            ])

    # Write issues report
    with open(issues_file, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["pdb_code", "issue_type", "details"])
        for pdb_code, issue_type, details in issues:
            writer.writerow([pdb_code, issue_type, details])

    # Print summary
    print()
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Report written to: {output_file}")
    print(f"Total components: {len(all_components)}")
    print(f"Processing issues written to: {issues_file}")
    print()
    print(f"Total PDBs processed: {len(cif_files)}")
    print(f"PDBs with components: {len(set(c.pdb_code for c in all_components))}")
    print(f"PDBs with issues: {len(issues)}")

    # Count by issue type
    issue_counts = {}
    for _, issue_type, _ in issues:
        issue_counts[issue_type] = issue_counts.get(issue_type, 0) + 1

    print()
    print("Issues by type:")
    for issue_type, count in sorted(issue_counts.items()):
        print(f"  {issue_type}: {count}")


if __name__ == "__main__":
    main()
