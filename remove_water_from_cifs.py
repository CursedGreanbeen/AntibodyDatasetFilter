#!/usr/bin/env python3
"""
Script to remove water (HOH) from CIF files in CIFs-filtered/ using PyMOL.

Usage:
    python scripts/remove_water_from_cifs.py
"""

from pathlib import Path
import subprocess
import sys

WATER_COMPLEXES = {
    "8ds7", "8g2m", "8rmo", "8tmz", "8uig", "8uih", "8vy4", "8yhz",
    "9bx5", "9bx7", "9c7d", "9c7x", "9cfd", "9ia3", "9j8a", "9mmj",
    "9mz6", "9mz7", "9mz8", "9nkz", "9nl0", "9p4c", "9pwn", "9t46",
    "9tpp", "9zmb", "9zmc"
}


def remove_water_from_cif(pdb_code: str, cif_dir: Path, temp_dir: Path) -> bool:
    """Remove water (HOH) from a CIF file using PyMOL command line."""
    cif_file = cif_dir / f"{pdb_code}_cropped.cif"

    if not cif_file.exists():
        print(f"  [SKIP] {cif_file.name} not found")
        return False

    # Create a temporary PyMOL script
    pymol_script = temp_dir / f"{pdb_code}_remove_hoh.pml"
    pymol_script.write_text(f"""
load {cif_file}, structure
remove resn HOH
save {cif_file}, structure
quit
""")

    try:
        result = subprocess.run(
            ['/usr/bin/pymol', '-cq', str(pymol_script)],
            capture_output=True,
            text=True,
            timeout=120
        )
        if result.returncode == 0:
            print(f"  [OK] Water removed from {cif_file.name}")
            return True
        else:
            print(f"  [ERROR] PyMOL error for {pdb_code}: {result.stderr}")
            return False
    except subprocess.TimeoutExpired:
        print(f"  [ERROR] Timeout for {pdb_code}")
        return False
    except Exception as e:
        print(f"  [ERROR] {e}")
        return False


def main():
    cif_dir = Path("/home/mullagaliamova/ClaudeWorkspace/PROJECTS/cdr-h3-folding/CIFs-filtered")
    temp_dir = Path("/tmp/pymol_scripts")
    temp_dir.mkdir(exist_ok=True)

    print("=" * 60)
    print("УДАЛЕНИЕ ВОДЫ (HOH) ИЗ CIF ФАЙЛОВ")
    print("=" * 60)
    print(f"\nВсего комплексов: {len(WATER_COMPLEXES)}")

    successful = []
    failed = []

    for pdb_code in sorted(WATER_COMPLEXES):
        if remove_water_from_cif(pdb_code, cif_dir, temp_dir):
            successful.append(pdb_code)
        else:
            failed.append(pdb_code)

    print("\n" + "=" * 60)
    print(f"ИТОГИ: {len(successful)} успешно, {len(failed)} ошибка")
    if failed:
        print(f"\nНеудачные: {', '.join(failed)}")


if __name__ == "__main__":
    main()
