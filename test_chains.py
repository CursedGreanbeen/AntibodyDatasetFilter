from pathlib import Path

from antibody_dataset.chains import get_structure_chains


cif_path = next(Path("CIFs-filtered-new").glob("*.cif"))

for chain in get_structure_chains(cif_path):
    print(
        chain.chain_id,
        chain.residue_count,
        chain.first_residue.number,
        chain.last_residue.number,
    )
