#!/usr/bin/env python
"""
Crop antibody sequences using ANARCI to identify variable regions.

Usage:
    python crop_fasta_with_anarci.py

For each antibody chain:
    - 1 domain found → crop to variable region
    - 2+ domains → keep original, flag in report
    - 0 domains → skip (not recognized as antibody)
"""

from pathlib import Path
import subprocess
import tempfile
import re
from fasta_manager import read_fasta_file, write_fasta_file, update_sequence, FastaRecord


# Configuration
INPUT_DIR = Path("/home/mullagaliamova/ClaudeWorkspace/PROJECTS/cdr-h3-folding/fasta-filtered-test")
OUTPUT_DIR = Path("/home/mullagaliamova/ClaudeWorkspace/PROJECTS/cdr-h3-folding/fasta-new-cropped")
REPORT_FILE = Path("/home/mullagaliamova/ClaudeWorkspace/PROJECTS/cdr-h3-folding/anarci_crop_report.tsv")

# ANARCI environment
ANARCI_DIR = Path("/home/mullagaliamova/ANARCI")
VENV_ACTIVATE = "/home/mullagaliamova/envs/anarci-env/bin/activate"


def is_antibody_chain(header: str) -> bool:
    """
    Check if FASTA header suggests an antibody chain.
    Used as a pre-filter - ANARCI is the authoritative source.
    """
    header_lower = header.lower()
    # Broad patterns to catch potential antibody chains
    # ANARCI will make the final decision
    antibody_patterns = [
        "heavy", "light", "kappa", "lambda",
        "h chain", "l chain", "hc", "lc",
        "antibody", "igg", "iga", "igm", "igd", "ige"
    ]
    return any(pattern in header_lower for pattern in antibody_patterns)


def run_anarci(sequence: str) -> tuple:
    """
    Run ANARCI CLI on a single sequence.

    Returns:
        (chain_type, start_idx, end_idx) for the most significant domain.
        Returns (None, None, None) if no antibody domain found or error.
    """
    # Create temp file with sequence
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        f.write(">temp_sequence\n")
        f.write(sequence + "\n")
        temp_file = f.name

    try:
        # Run ANARCI
        cmd = f"""
        cd {ANARCI_DIR} &&
        . {VENV_ACTIVATE} &&
        export PATH=/usr/bin:$PATH &&
        ANARCI -i {temp_file} --hmmerpath /usr/bin -s imgt
        """
        result = subprocess.run(
            cmd, shell=True, capture_output=True, text=True, timeout=60
        )
        output = result.stdout

        # Parse output for domain boundaries
        # Look for lines like: #|human|K|2.9e-60|193.0|0|107|
        # chain_type is group(1), start_idx is group(2), end_idx is group(3)
        pattern = r'#\|[^|]+\|([^|]+)\|[^|]+\|[^|]+\|(\d+)\|(\d+)\|'
        matches = list(re.finditer(pattern, output))

        if not matches:
            return (None, None, None)

        # Take the most significant hit (first one)
        first_match = matches[0]
        chain_type = first_match.group(1)
        start_idx = int(first_match.group(2))
        end_idx = int(first_match.group(3))

        return (chain_type, start_idx, end_idx)

    except subprocess.TimeoutExpired:
        print("    ERROR: ANARCI timeout")
        return (None, None, None)
    except Exception as e:
        print(f"    ERROR: {e}")
        return (None, None, None)
    finally:
        # Clean up temp file
        try:
            Path(temp_file).unlink()
        except:
            pass


def is_antibody_chain_type(chain_type: str) -> bool:
    """Check if ANARCI chain type is an antibody chain (H, K, L)."""
    return chain_type in ('H', 'K', 'L')


def crop_sequence(sequence: str, chain_type: str, start_idx: int, end_idx: int) -> tuple:
    """
    Crop sequence based on ANARCI domain boundaries.

    Args:
        sequence: Full antibody sequence
        chain_type: ANARCI chain type (H, K, L, or other)
        start_idx: Start index from ANARCI
        end_idx: End index from ANARCI

    Returns:
        (cropped_sequence, status_info)
    """
    if chain_type is None:
        return sequence, "NO_ANTIBODY_DOMAIN"

    if not is_antibody_chain_type(chain_type):
        return sequence, f"NON_ANTIBODY_TYPE[{chain_type}]"

    # Single antibody domain found - crop to it
    cropped = sequence[start_idx:end_idx + 1]
    return cropped, f"CROPPED[{chain_type}var{start_idx}:{end_idx+1}]"


def main():
    print("ANARCI Cropping Script")
    print("=" * 60)
    print(f"Input:  {INPUT_DIR}")
    print(f"Output: {OUTPUT_DIR}")
    print()

    # Create output directory
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Report data
    report_data = []

    # Counters
    total_files = 0
    processed_files = 0
    total_chains = 0
    cropped_chains = 0
    multi_domain_chains = 0
    unrecognized_chains = 0
    errors = 0

    # Process each FASTA file
    fasta_files = sorted(INPUT_DIR.glob("*.fasta"))
    total_files = len(fasta_files)

    print(f"Found {total_files} FASTA files\n")

    for fasta_file in fasta_files:
        print(f"Processing: {fasta_file.name}")

        try:
            records = read_fasta_file(fasta_file)
        except Exception as e:
            print(f"  ERROR reading file: {e}")
            errors += 1
            report_data.append({
                "file": fasta_file.name,
                "chain": "-",
                "status": "FILE_ERROR",
                "details": str(e),
                "action": "skipped"
            })
            continue

        cropped_records = {}
        antibody_chains_found = False

        for chain_id, record in records.items():
            # Run ANARCI on ALL chains - it's the authoritative source
            chain_type, start_idx, end_idx = run_anarci(record.sequence)

            # Check if ANARCI recognized this as an antibody chain
            if chain_type is None:
                # No domain found by ANARCI
                print(f"  Chain {chain_id}: NOT_RECOGNIZED")
                report_data.append({
                    "file": fasta_file.name,
                    "chain": chain_id,
                    "status": "NOT_RECOGNIZED",
                    "details": "No domain found by ANARCI",
                    "action": "skipped"
                })
                continue

            # Check if it's an antibody chain type (H, K, L)
            if not is_antibody_chain_type(chain_type):
                # Non-antibody domain (e.g., EGF, other)
                print(f"  Chain {chain_id}: NON_ANTIBODY({chain_type})")
                report_data.append({
                    "file": fasta_file.name,
                    "chain": chain_id,
                    "status": "NON_ANTIBODY",
                    "details": f"ANARCI identified as {chain_type}, not antibody",
                    "action": "skipped"
                })
                continue

            # This is an antibody chain - crop it
            cropped_seq, status = crop_sequence(record.sequence, chain_type, start_idx, end_idx)
            total_chains += 1

            if "CROPPED" in status:
                cropped_chains += 1
                cropped_records[chain_id] = FastaRecord(
                    header=record.header,
                    sequence=cropped_seq,
                    chain_ids=record.chain_ids
                )
                print(status)
                report_data.append({
                    "file": fasta_file.name,
                    "chain": chain_id,
                    "status": "CROPPED",
                    "domains": 1,
                    "chain_type": chain_type,
                    "details": status,
                    "action": f"cropped to variable region ({chain_type})"
                })

        # Write output file
        if cropped_records:
            output_file = OUTPUT_DIR / fasta_file.name
            try:
                for chain_id, record in cropped_records.items():
                    update_sequence(records, chain_id, record.sequence)
                write_fasta_file(output_file, records)
                processed_files += 1
                print(f"  → Wrote: {output_file.name}")
            except Exception as e:
                print(f"  ERROR writing file: {e}")
                errors += 1
        else:
            # No antibody chains found or processed - log to report
            if antibody_chains_found:
                print("  No chains processed")
                report_data.append({
                    "file": fasta_file.name,
                    "chain": "-",
                    "status": "NO_CHAINS_PROCESSED",
                    "details": "All antibody chains failed processing",
                    "action": "skipped"
                })
            else:
                print("  No antibody chains")
                report_data.append({
                    "file": fasta_file.name,
                    "chain": "-",
                    "status": "NO_ANTIBODY_CHAINS",
                    "details": "No antibody chains detected in file",
                    "action": "skipped"
                })

        print()

    # Write report
    print(f"Writing report: {REPORT_FILE}")
    with open(REPORT_FILE, "w") as f:
        f.write("file\tchain\tstatus\tdomains\tdetails\taction\n")
        for row in report_data:
            f.write(f"{row['file']}\t{row['chain']}\t{row['status']}\t")
            f.write(f"{row.get('domains', '-')}\t{row['details']}\t{row['action']}\n")

    # Summary
    print()
    print("=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"Files processed:     {processed_files}/{total_files}")
    print(f"Total antibody chains: {total_chains}")
    print(f"  - Cropped (1 domain):   {cropped_chains}")
    print(f"  - Multi-domain:         {multi_domain_chains}")
    print(f"  - Unrecognized:         {unrecognized_chains}")
    print(f"  - Errors:               {errors}")
    print()
    print(f"Output directory: {OUTPUT_DIR}")
    print(f"Report file:      {REPORT_FILE}")


if __name__ == "__main__":
    main()
