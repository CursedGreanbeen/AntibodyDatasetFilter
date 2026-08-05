#!/usr/bin/env python3
"""Count unique PDB IDs from SABDab TSV file."""

import argparse


def count_unique_pdbs(tsv_path: str, column_index: int = 0) -> None:
    """Count unique PDB IDs in a TSV file."""
    pdbs = set()
    with open(tsv_path, 'r') as f:
        header = f.readline().strip().split('\t')
        print(f"Columns: {header}")
        print(f"Using column {column_index}: '{header[column_index]}'")
        for line in f:
            parts = line.strip().split('\t')
            if parts and parts[column_index]:
                pdbs.add(parts[column_index])
    print(f"Unique PDB IDs: {len(pdbs)}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("tsv_path", help="Path to TSV file")
    parser.add_argument("--column", type=int, default=0, help="Column index (0-based)")
    args = parser.parse_args()

    count_unique_pdbs(args.tsv_path, args.column)
