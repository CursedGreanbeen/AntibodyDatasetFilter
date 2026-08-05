#!/usr/bin/env python3
"""
Fetch antibody sequences from RCSB PDB API for a list of PDB codes.
Uses the PDB REST API to retrieve protein sequences in FASTA format.
"""

import argparse
import csv
import os
import time
from pathlib import Path

import requests
from urllib3.exceptions import InsecureRequestWarning
from urllib3 import disable_warnings
from io import StringIO

# Disable SSL warnings
disable_warnings(InsecureRequestWarning)


PDB_API_URL = "https://data.rcsb.org/rest/v1/core/entry/{}"
PDB_SEQ_API_URL = "https://data.rcsb.org/rest/v1/sequence/{}"


def fetch_pdb_entry(pdb_code: str) -> dict | None:
    """Fetch PDB entry metadata from RCSB API."""
    url = PDB_API_URL.format(pdb_code.lower())
    try:
        response = requests.get(url, timeout=30, verify=False)
        if response.status_code == 200:
            return response.json()
        else:
            print(f"  Warning: HTTP {response.status_code} for {pdb_code}")
            return None
    except requests.RequestException as e:
        print(f"  Error fetching {pdb_code}: {e}")
        return None


def fetch_polymer_sequences(pdb_code: str) -> list[dict] | None:
    """
    Fetch polymer (protein/DNA/RNA) sequences from PDB.
    Returns list of dicts with entity_id, sequence, description.
    """
    url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_code.lower()}/1"
    try:
        response = requests.get(url, timeout=30)
        if response.status_code == 200:
            return response.json()
        return None
    except requests.RequestException as e:
        print(f"  Error fetching sequences for {pdb_code}: {e}")
        return None


def get_antibody_chains(pdb_code: str, output_dir: Path) -> int:
    """
    Fetch and save antibody sequences for a PDB code.
    Returns number of chains fetched.
    Uses RCSB FASTA download endpoint.
    """
    # Use RCSB FASTA endpoint
    url = f"https://www.rcsb.org/fasta/entry/{pdb_code.lower()}"

    try:
        response = requests.get(url, timeout=30, verify=False)
        if response.status_code != 200:
            return 0

        # Create FASTA file for this PDB
        fasta_path = output_dir / f"{pdb_code.lower()}.fasta"

        with open(fasta_path, 'w') as f:
            f.write(response.text)

        # Count chains from FASTA headers
        chain_count = response.text.count('>')

        return chain_count

    except requests.RequestException as e:
        print(f"  Error: {e}")
        return 0


def download_cif_file(pdb_code: str, output_dir: Path) -> bool:
    """
    Download .cif (mmCIF) file for a PDB code from RCSB PDB.
    Returns True if successful, False otherwise.
    """
    # RCSB PDB provides CIF files at this endpoint
    url = f"https://files.rcsb.org/download/{pdb_code.lower()}.cif"

    try:
        response = requests.get(url, timeout=300, verify=False)  # 5min timeout for large files
        if response.status_code != 200:
            print(f"  Warning: Could not fetch CIF for {pdb_code} (HTTP {response.status_code})")
            return False

        # Create CIF file for this PDB
        cif_path = output_dir / f"{pdb_code.lower()}.cif"

        with open(cif_path, 'w') as f:
            f.write(response.text)

        return True

    except requests.RequestException as e:
        print(f"  Error downloading CIF for {pdb_code}: {e}")
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Fetch antibody sequences from PDB API for a list of PDB codes"
    )
    parser.add_argument(
        "input_tsv",
        help="Input TSV file with PDB codes (first column)"
    )
    parser.add_argument(
        "-ofasta", "--output-fasta",
        default="sequences",
        help="Output directory for FASTA files (default: fasta-sequences/)"
    )
    parser.add_argument(
        "-ocif", "--output-cif",
        default="sequences",
        help="Output directory for CIF files (default: CIFs/)"
    )
    parser.add_argument(
        "--rate-limit",
        type=float,
        default=0.2,
        help="Delay between requests in seconds (default: 0.2s)"
    )
    parser.add_argument(
        "--fasta",
        action="store_true",
        help="Download FASTA files"
    )
    parser.add_argument(
        "--cif",
        action="store_true",
        help="Download .cif (mmCIF) files"
    )

    args = parser.parse_args()

    # At least one format must be specified
    if not args.fasta and not args.cif:
        parser.error("At least one of --fasta or --cif is required")

    download_fasta = args.fasta
    download_cif = args.cif

    output_fasta = Path(args.output_fasta)
    output_fasta.mkdir(parents=True, exist_ok=True)
    output_cif = Path(args.output_cif)
    output_cif.mkdir(parents=True, exist_ok=True)

    # Read PDB codes from TSV
    pdb_codes = []
    with open(args.input_tsv, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            pdb_codes.append(row['pdb'])

    # Get unique codes
    unique_codes = list(dict.fromkeys(pdb_codes))  # Preserve order
    print(f"Found {len(unique_codes)} unique PDB codes")

    # Fetch sequences
    success_count = 0
    fail_count = 0

    for i, pdb_code in enumerate(unique_codes, 1):
        pdb_code_lower = pdb_code.lower()

        # Download FASTA if requested
        if download_fasta:
            fasta_path = output_fasta / f"{pdb_code_lower}.fasta"

            if fasta_path.exists():
                print(f"[{i}/{len(unique_codes)}] {pdb_code}: FASTA already exists (skipping)")
                success_count += 1
            else:
                chains = get_antibody_chains(pdb_code, output_fasta)

                if chains > 0:
                    print(f"[{i}/{len(unique_codes)}] {pdb_code}: {chains} chains fetched")
                    success_count += 1
                else:
                    print(f"[{i}/{len(unique_codes)}] {pdb_code}: FASTA fetch FAILED")
                    fail_count += 1

        # Download CIF if requested
        if download_cif:
            cif_path = output_cif / f"{pdb_code_lower}.cif"
            if cif_path.exists():
                print(f"[{i}/{len(unique_codes)}] {pdb_code}: CIF already exists (skipping)")
                success_count += 1
            else:
                if download_cif_file(pdb_code, output_cif):
                    print(f"[{i}/{len(unique_codes)}] {pdb_code}: CIF downloaded")
                else:
                    print(f"[{i}/{len(unique_codes)}] {pdb_code}: CIF fetch FAILED")
                    fail_count += 1

        # Rate limiting
        if i < len(unique_codes):
            time.sleep(args.rate_limit)

    print(f"\nDone! Success: {success_count}, Failed: {fail_count}")
    print(f"FASTA files saved to: {output_fasta}")
    print(f"CIF files saved to: {output_cif}")


if __name__ == "__main__":
    main()
