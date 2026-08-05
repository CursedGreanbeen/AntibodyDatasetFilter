#!/usr/bin/env python3
"""Compare two PDB/cif structures and calculate RMSD per chain."""

from pathlib import Path
import numpy as np
from Bio.PDB import MMCIFParser, Superimposer

def get_chain_info(structure):
    """Extract CA coordinates and chain info."""
    chains = {}
    for model in structure:
        for chain in model:
            chain_id = chain.id
            atoms = []
            residue_ids = []
            for residue in chain:
                if 'CA' in residue:
                    atoms.append(residue['CA'].get_coord())
                    residue_ids.append(residue.resname)
            if atoms:
                chains[chain_id] = {
                    'coords': np.array(atoms),
                    'length': len(atoms),
                    'resnames': residue_ids
                }
    return chains

def calculate_rmsd_with_superposition(coords1, coords2):
    """Calculate RMSD after optimal superposition using Kabsch algorithm."""
    min_len = min(len(coords1), len(coords2))
    c1 = coords1[:min_len].copy()
    c2 = coords2[:min_len].copy()

    # Center both structures at origin
    c1_center = np.mean(c1, axis=0)
    c2_center = np.mean(c2, axis=0)
    c1 -= c1_center
    c2 -= c2_center

    # Kabsch algorithm
    H = np.dot(c1.T, c2)
    U, S, Vt = np.linalg.svd(H)
    R = np.dot(Vt.T, U.T)

    # Ensure right-handed coordinate system
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = np.dot(Vt.T, U.T)

    # Apply rotation
    c1_rotated = np.dot(c1, R)

    # Calculate RMSD
    diff = c1_rotated - c2
    rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
    return rmsd, min_len

def match_chains(ref_chains, pred_chains):
    """Match chains based on user-verified correspondence."""
    # User-verified mapping: Ref H->Pred A, Ref L->Pred C, Ref K->Pred B
    # Ref A (C2b fragment) may not have a match in prediction
    user_mapping = {
        'H': 'A',  # Heavy chain
        'L': 'C',  # Light chain
        'K': 'B',  # Nanobody
    }

    matches = []
    for ref_id, pred_id in user_mapping.items():
        if ref_id in ref_chains and pred_id in pred_chains:
            matches.append((ref_id, pred_id, ref_chains[ref_id], pred_chains[pred_id]))

    return matches

def main():
    parser = MMCIFParser(QUIET=True)

    ref_path = Path("/home/mullagaliamova/cdr-h3-folding/sequences/s3data/8acf/8ACF.cif")
    pred_path = Path("/home/mullagaliamova/cdr-h3-folding/sequences/s3data/8acf/8acf_model.cif")

    ref_structure = parser.get_structure("reference", str(ref_path))
    pred_structure = parser.get_structure("predicted", str(pred_path))

    ref_chains = get_chain_info(ref_structure)
    pred_chains = get_chain_info(pred_structure)

    print("=" * 70)
    print("STRUCTURE COMPARISON REPORT")
    print("=" * 70)
    print(f"\nReference: {ref_path.name}")
    print(f"Predicted: {pred_path.name}")
    print(f"\nReference chains: {list(ref_chains.keys())} (lengths: {[(k, v['length']) for k, v in ref_chains.items()]})")
    print(f"Predicted chains: {list(pred_chains.keys())} (lengths: {[(k, v['length']) for k, v in pred_chains.items()]})")
    print("\nMapping (based on visual verification):")
    print("  H (Heavy chain) -> A")
    print("  L (Light chain) -> C")
    print("  K (Nanobody) -> B")
    print("\n" + "-" * 70)
    print(f"{'Ref Chain':<12} {'Pred Chain':<12} {'RMSD (Å)':<15} {'Ref Len':<10} {'Pred Len':<10}")
    print("-" * 70)

    all_rmsds = []
    matches = match_chains(ref_chains, pred_chains)

    for ref_id, pred_id, ref_data, pred_data in matches:
        rmsd, n_residues = calculate_rmsd_with_superposition(ref_data['coords'], pred_data['coords'])
        all_rmsds.append(rmsd)
        print(f"{ref_id:<12} -> {pred_id:<12} {rmsd:<15.3f} {ref_data['length']:<10} {pred_data['length']:<10}")

    print("-" * 70)
    if all_rmsds:
        print(f"{'OVERALL (mean)':<24} {np.mean(all_rmsds):<15.3f}")
        print(f"{'Min RMSD':<24} {np.min(all_rmsds):<15.3f}")
        print(f"{'Max RMSD':<24} {np.max(all_rmsds):<15.3f}")
    print("=" * 70)

if __name__ == "__main__":
    main()
