import gemmi

struct = gemmi.read_structure("CIFs-filtered-new/10op.cif")
model = struct[0]
cell = struct.cell

# Check distances between chains
chain_y = model.get_chain('Y')
chain_a = model.get_chain('A')

if chain_y and chain_a:
    # Get first atom from each
    res_y = list(chain_y)[0]
    atom_y = list(res_y)[0]

    res_a = list(chain_a)[0]
    atom_a = list(res_a)[0]

    dist = atom_y.pos.distance(atom_a.pos)
    print(f"Distance between Chain Y residue 1 and Chain A residue 1: {dist:.2f} Å")

    # Check all Y atoms against all A atoms
    min_dist = float('inf')
    for res_y in chain_y:
        for atom_y in res_y:
            for res_a in chain_a:
                for atom_a in res_a:
                    d = atom_y.pos.distance(atom_a.pos)
                    if d < min_dist:
                        min_dist = d
    print(f"Minimum distance between Chain Y and Chain A: {min_dist:.2f} Å")
