#!/usr/bin/env python3
"""
Script to replace CIF files in CIFs-filtered/ with original files from CIFs/.

Usage:
    python scripts/replace_cif_files.py

The script will prompt for a comma-separated list of PDB codes.
For each code, the original file from CIFs/{pdb_code}.cif is copied to
CIFs-filtered/{pdb_code}_cropped.cif.

A log file is created at CIFs-filtered/replacement_log.txt.
"""

import shutil
from datetime import datetime
from pathlib import Path

# Paths relative to script location
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
CIFS_DIR = PROJECT_ROOT / "CIFs"
CIFS_FILTERED_DIR = PROJECT_ROOT / "CIFs-filtered"
LOG_FILE = CIFS_FILTERED_DIR / "replacement_log.txt"


def parse_pdb_codes(input_string: str) -> list[str]:
    """Parse comma-separated PDB codes, stripping whitespace."""
    if not input_string.strip():
        return []
    return [code.strip().lower() for code in input_string.split(",") if code.strip()]


def main() -> None:
    # Ensure output directory exists
    CIFS_FILTERED_DIR.mkdir(parents=True, exist_ok=True)

    # Get PDB codes from user
    print("Enter PDB codes (comma-separated, e.g., 4ghg, 5lgb, 3l9g):")
    user_input = input("> ")

    pdb_codes = parse_pdb_codes(user_input)

    if not pdb_codes:
        print("No PDB codes provided. Exiting.")
        return

    print(f"\nProcessing {len(pdb_codes)} PDB code(s): {', '.join(pdb_codes)}\n")

    # Track results for logging
    successful: list[str] = []
    failed: list[str] = []

    for pdb_code in pdb_codes:
        source_file = CIFS_DIR / f"{pdb_code}.cif"
        dest_file = CIFS_FILTERED_DIR / f"{pdb_code}_cropped.cif"

        if not source_file.exists():
            print(f"  [SKIP] {pdb_code}.cif not found in CIFs/")
            failed.append(pdb_code)
            continue

        try:
            shutil.copy2(source_file, dest_file)
            print(f"  [OK] Copied {pdb_code}.cif -> {pdb_code}_cropped.cif")
            successful.append(pdb_code)
        except Exception as e:
            print(f"  [ERROR] Failed to copy {pdb_code}.cif: {e}")
            failed.append(pdb_code)

    # Write log file
    log_content = f"""CIF File Replacement Log
========================
Timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
Input PDB codes: {', '.join(pdb_codes)}

Summary
-------
Total codes: {len(pdb_codes)}
Successful: {len(successful)}
Failed/Skipped: {len(failed)}

Successful replacements:
{chr(10).join(f'  - {code} -> {code}_cropped.cif' for code in successful)}

Failed/Skipped:
{chr(10).join(f'  - {code}: file not found or copy error' for code in failed)}
"""

    LOG_FILE.write_text(log_content)

    print(f"\nLog written to: {LOG_FILE}")
    print(f"Completed: {len(successful)} successful, {len(failed)} failed/skipped")


if __name__ == "__main__":
    main()
