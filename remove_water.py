from Bio.PDB import MMCIFParser, PDBIO, Select
import os

WATER_COMPLEXES = {
    "8ds7", "8g2m", "8rmo", "8tmz", "8uig", "8uih", "8vy4", "8yhz",
    "9bx5", "9bx7", "9c7d", "9c7x", "9cfd", "9ia3", "9j8a", "9mmj",
    "9mz6", "9mz7", "9mz8", "9nkz", "9nl0", "9p4c", "9pwn", "9t46",
    "9tpp", "9zmb", "9zmc"
}

class NoWaterSelect(Select):
    def accept_residue(self, residue):
        return residue.resname not in ("HOH", "WAT", "H2O", "DOD")

parser = MMCIFParser(QUIET=True)
io = PDBIO()

input_dir = "./CIFs-filtered"  # папка с исходными файлами
output_dir = "./CIFs-filtered"  # та же папка для записи

for pdb_id in WATER_COMPLEXES:
    input_path = os.path.join(input_dir, f"{pdb_id}_cropped.cif")
    output_path = os.path.join(output_dir, f"{pdb_id}_cropped.cif")

    if not os.path.exists(input_path):
        print(f"[SKIP] {pdb_id} — файл не найден")
        continue

    structure = parser.get_structure(pdb_id, input_path)
    io.set_structure(structure)
    io.save(output_path, NoWaterSelect())
    print(f"[OK] {pdb_id} → {output_path}")
