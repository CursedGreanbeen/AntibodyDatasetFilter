#!/usr/bin/env python3
"""
Filter FASTA files to keep only heavy, light chains, and antigen chains.
Uses SabDab TSV table to identify chain types (Hchain, Lchain, antigen_chain columns).
Input: FASTA files from sequences/ directory
Output: Filtered FASTA files in fasta-filtered/ directory, multimers in multimers/
"""

import argparse
import csv
import re
import shutil
from pathlib import Path


def parse_fasta(fasta_path: Path) -> list[tuple[str, str, str]]:
    """
    Parse FASTA file and return list of (header, sequence, full_header).
    """
    entries = []
    current_header = None
    current_sequence = []

    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_header is not None:
                    entries.append((current_header, ''.join(current_sequence), full_header))
                full_header = line[1:]  # Remove '>'
                current_header = full_header.lower()
                current_sequence = []
            else:
                current_sequence.append(line)

        # Don't forget the last entry
        if current_header is not None:
            entries.append((current_header, ''.join(current_sequence), full_header))

    return entries


def extract_chain_ids(full_header: str) -> list[str]:
    """
    Extract chain ID(s) from FASTA header.
    Expected formats:
      - >...|Chain X[auth Y]|... → returns ['Y']
      - >...|Chains X, Y, Z|... → returns ['X', 'Y', 'Z']
    Returns list of chain IDs (e.g., ['A'], ['T', 'V', 'X'])
    """
    chain_ids = []

    # Try format: Chain X[auth Y] - extract the auth letter
    match = re.search(r'Chain\s\w+\[auth\s+([A-Za-z])\]', full_header, re.IGNORECASE)
    if match:
        return [match.group(1).upper()]

    # Try format: Chains X, Y, Z - extract all letters after "Chains"
    match = re.search(r'Chains?\s+([A-Za-z,\s]+)', full_header, re.IGNORECASE)
    if match:
        chains_str = match.group(1)
        chain_ids = [c.strip().upper() for c in re.split(r'[,\s]+', chains_str) if c.strip()]

    return chain_ids


def load_chain_data(tsv_path: Path) -> tuple[dict[str, dict], dict[str, int]]:
    """
    Load chain data from SabDab TSV file.
    Returns:
        - dict: pdb_code -> {
            'heavy_chains': set(),
            'light_chains': set(),
            'antigen_chains': set()
        }
        - dict: pdb_code -> row_count (number of rows in TSV for each PDB)
    """
    chain_data = {}
    row_counts = {}

    with open(tsv_path, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            pdb_code = row['pdb'].lower()
            hchain = row.get('Hchain', 'NA').strip()
            lchain = row.get('Lchain', 'NA').strip()
            antigen_chain = row.get('antigen_chain', 'NA')
            antigen_type = row.get('antigen_type', 'NA')

            # Initialize entry for this PDB if not exists
            if pdb_code not in chain_data:
                chain_data[pdb_code] = {
                    'heavy_chains': set(),
                    'light_chains': set(),
                    'antigen_chains': set()
                }
                row_counts[pdb_code] = 0
            row_counts[pdb_code] += 1

            # Add heavy chain if present
            if hchain != 'NA' and hchain:
                chain_data[pdb_code]['heavy_chains'].add(hchain.upper())

            # Add light chain if present
            if lchain != 'NA' and lchain:
                chain_data[pdb_code]['light_chains'].add(lchain.upper())

            # Add antigen chains if present
            if antigen_chain != 'NA' and antigen_type != 'NA':
                # Parse chains and types (they are separated by ' | ')
                chains = [c.strip() for c in antigen_chain.split('|')]
                types = [t.strip() for t in antigen_type.split('|')]

                for chain, chain_type in zip(chains, types):
                    if chain_type.lower() == 'protein':
                        chain_data[pdb_code]['antigen_chains'].add(chain.upper())

    return chain_data, row_counts


def filter_fasta_file(input_path: Path, chain_data: dict[str, dict], row_counts: dict[str, int]) -> dict:
    """
    Filter a single FASTA file.
    Returns dict with: pdb_code, has_heavy, has_light, antigen_count, chain_count, success, is_multimer
    """
    pdb_code = input_path.stem.lower()  # Filename without extension

    # Check if PDB exists in chain_data
    if pdb_code not in chain_data:
        return {
            'pdb_code': pdb_code,
            'has_heavy': False,
            'has_light': False,
            'heavy_count': 0,
            'light_count': 0,
            'antigen_count': 0,
            'chain_count': 0,
            'success': False,
            'is_multimer': False,
            'not_in_table': True
        }

    entries = parse_fasta(input_path)
    pdb_data = chain_data[pdb_code]
    heavy_chains = pdb_data['heavy_chains']
    light_chains = pdb_data['light_chains']
    antigen_chains = pdb_data['antigen_chains']

    heavy_entries = []
    light_entries = []
    antigen_entries = []

    for header, sequence, full_header in entries:
        chain_ids = extract_chain_ids(full_header)

        # Determine chain type based on SabDab table
        for chain_id in chain_ids:
            if chain_id in heavy_chains:
                heavy_entries.append((full_header, sequence))
                break
            elif chain_id in light_chains:
                light_entries.append((full_header, sequence))
                break
            elif chain_id in antigen_chains:
                antigen_entries.append((full_header, sequence))

    has_heavy = len(heavy_entries) > 0
    has_light = len(light_entries) > 0
    antigen_count = len(antigen_entries)

    # Check if multimer (multiple rows in TSV for this PDB)
    is_multimer = row_counts.get(pdb_code, 1) > 1

    result = {
        'pdb_code': pdb_code,
        'has_heavy': has_heavy,
        'has_light': has_light,
        'heavy_count': len(heavy_entries),
        'light_count': len(light_entries),
        'antigen_count': antigen_count,
        'chain_count': len(heavy_entries) + len(light_entries) + antigen_count,
        'success': has_heavy and has_light,
        'is_multimer': is_multimer,
        'not_in_table': False
    }

    return result


def main():
    parser = argparse.ArgumentParser(
        description="Filter FASTA files to keep heavy, light, and antigen chains based on SabDab table"
    )
    parser.add_argument(
        "input_dir",
        help="Input directory with FASTA files (e.g., sequences/)"
    )
    parser.add_argument(
        "-o", "--output-dir",
        default="fasta-filtered",
        help="Output directory for filtered FASTA files (default: fasta-filtered/)"
    )
    parser.add_argument(
        "-t", "--tsv",
        required=True,
        help="TSV file with Hchain, Lchain, antigen_chain columns (e.g., sabdab_summary_2024_plus.tsv)"
    )

    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Create multimers directory
    multimers_dir = output_dir.parent / "multimers"
    multimers_dir.mkdir(parents=True, exist_ok=True)

    # Load chain data from TSV
    chain_data, row_counts = load_chain_data(Path(args.tsv))
    print(f"Loaded chain data for {len(chain_data)} PDB entries")

    # Find all FASTA files
    fasta_files = list(input_dir.glob("*.fasta"))
    print(f"Found {len(fasta_files)} FASTA files")

    # Process each file
    results = []
    success_count = 0
    skip_count = 0
    multimer_count = 0
    antigen_included_count = 0

    for i, fasta_file in enumerate(fasta_files, 1):
        result = filter_fasta_file(fasta_file, chain_data, row_counts)
        results.append(result)

        pdb_code = result['pdb_code']

        if result['not_in_table']:
            print(f"[{i}/{len(fasta_files)}] {pdb_code}: not found in TSV -> skipped")
            skip_count += 1
        elif result['success']:
            if result['is_multimer']:
                # Copy entire file to multimers/
                output_path = multimers_dir / f"{pdb_code}.fasta"
                import shutil
                shutil.copy2(fasta_file, output_path)
                msg = f"[{i}/{len(fasta_files)}] {pdb_code}: "
                msg += f"{result['heavy_count']} heavy + {result['light_count']} light"
                if result['antigen_count'] > 0:
                    msg += f" + {result['antigen_count']} antigen"
                msg += " chains -> multimer (saved to multimers/)"
                print(msg)
                multimer_count += 1
            else:
                # Write filtered file to output_dir
                output_path = output_dir / f"{pdb_code}.fasta"
                heavy_entries = []
                light_entries = []
                antigen_entries = []

                # Re-parse to get entries (simplified - in production would return from filter_fasta_file)
                entries = parse_fasta(fasta_file)
                pdb_data = chain_data[pdb_code]

                for header, sequence, full_header in entries:
                    chain_ids = extract_chain_ids(full_header)
                    for chain_id in chain_ids:
                        if chain_id in pdb_data['heavy_chains']:
                            heavy_entries.append((full_header, sequence))
                            break
                    for chain_id in chain_ids:
                        if chain_id in pdb_data['light_chains']:
                            light_entries.append((full_header, sequence))
                            break
                    for chain_id in chain_ids:
                        if chain_id in pdb_data['antigen_chains']:
                            antigen_entries.append((full_header, sequence))
                            break

                with open(output_path, 'w') as f:
                    for full_header, sequence in heavy_entries + light_entries + antigen_entries:
                        f.write(f">{full_header}\n{sequence}\n")

                msg = f"[{i}/{len(fasta_files)}] {pdb_code}: "
                msg += f"{result['heavy_count']} heavy + {result['light_count']} light"
                if result['antigen_count'] > 0:
                    msg += f" + {result['antigen_count']} antigen"
                msg += " chains -> saved"
                print(msg)
                success_count += 1
                if result['antigen_count'] > 0:
                    antigen_included_count += 1
        else:
            print(f"[{i}/{len(fasta_files)}] {pdb_code}: "
                  f"missing {'heavy' if not result['has_heavy'] else 'light'} chain -> skipped")
            skip_count += 1

    # Write summary TSV
    summary_path = output_dir / "filtered_summary.tsv"
    with open(summary_path, 'w') as f:
        f.write("pdb_code\thas_heavy\thas_light\theavy_count\tlight_count\tantigen_count\tchain_count\tis_multimer\n")
        for r in results:
            f.write(f"{r['pdb_code']}\t{r['has_heavy']}\t{r['has_light']}\t"
                    f"{r['heavy_count']}\t{r['light_count']}\t{r['antigen_count']}\t{r['chain_count']}\t{r['is_multimer']}\n")

    print(f"\nDone!")
    print(f"Success: {success_count} files with both heavy and light chains")
    print(f"Multimers: {multimer_count} files saved to multimers/")
    print(f"Antigen chains included: {antigen_included_count} files")
    print(f"Skipped: {skip_count} files")
    print(f"Filtered FASTA files saved to: {output_dir}")
    print(f"Summary saved to: {summary_path}")


if __name__ == "__main__":
    main()
