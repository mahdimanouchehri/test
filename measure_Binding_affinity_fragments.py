# -*- coding: utf-8 -*-

from pdb_manip_py import pdb_manip
from docking_py import docking
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdmolfiles
import pandas as pd


df = pd.read_csv("/content/test/geom_multifrag_test_table.csv")

# Create a Coor object
coor_1hsg = pdb_manip.Coor()
coor_1hsg.get_PDB('8X5Y', 'data/8X5Y.pdb')


# Select res_name MK1
lig_coor = coor_1hsg.select_part_dict(selec_dict={'res_name': ['MK1']})
# Save the ligand coordinates
lig_coor.write_pdb('data/lig.pdb')

# Keep only the amino acids
rec_coor = coor_1hsg.select_part_dict(selec_dict={'res_name': pdb_manip.PROTEIN_RES})
rec_coor.write_pdb('data/rec.pdb')




for i in range(len(df)):
  fragment_list = df["fragments"][0].split(".")
  j = 1
  for fragment in fragment_list:
    mol = Chem.MolFromSmiles(fragment)

    # Add hydrogens to the molecule
    mol = Chem.AddHs(mol)

    # Generate 3D conformers
    AllChem.EmbedMolecule(mol)

    # Save the molecule to a PDB file
    rdmolfiles.MolToPDBFile(mol, f"output{i}{j}.pdb")



#    df[f"binding_affinity_frag{j}"][i] = binding_affinity1

    test_dock = docking.Docking('test', lig_pdb=f"output{i}{j}.pdb", rec_pdb='data/rec.pdb')
    test_dock.prepare_ligand()
    test_dock.prepare_receptor()

    test_dock.run_docking(out_pdb=f'test_dock{i}{j}.pdb',
                      num_modes=10,
                      energy_range=10,
                      exhaustiveness=16,
                      dock_bin='smina')
    print(test_dock.affinity)
    df.loc[i, f"binding_affinity_frag{j}"] = test_dock.affinity[1]["affinity"]

    j += 1


df.to_csv("geom_multifrag_test_table_with_binding_affinity.csv")




