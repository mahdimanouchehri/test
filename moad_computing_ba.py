# -*- coding: utf-8 -*-
"""moad_Computing_BA.ipynb

Automatically generated by Colab.

Original file is located at
    https://colab.research.google.com/drive/1qDhv0x2g3QcJp0rq6lPz3f7QOvUfvWin
"""

from pdb_manip_py import pdb_manip
from docking_py import docking
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdmolfiles
import pandas as pd
import subprocess
import glob
import os
def create_dataframe(*args):
  # input column name
  # dictionary of lists
  dict = {}
  for i in args :
    dict[i] = []

  df = pd.DataFrame(dict)
  return df

df = create_dataframe("idx", "BA protein1", "BA protein2")
pr_list = ["8X5Y", "1JQE"]

import os

dir_path = "/content/drive/MyDrive/These/GNN_raw_data"
for index_path in os.listdir(dir_path):
  full_path = os.path.join(dir_path, index_path)
  lock_file = f"{full_path}/{index_path}.lock"
  # Check if the file already exists
  if os.path.exists(f"{full_path}/{index_path}.csv"):
    continue

  # If the lock file exists (another system is using the folder), skip to the next iteration
  if os.path.exists(lock_file):
    continue

  # Create a lock file to signal other systems that this folder is in use
  open(lock_file, 'a').close()
  print(f"working on {index_path}")
  for ligand_path in os.listdir(full_path):
    if ligand_path.endswith('.pdb'):
      save_inf= []
      ligand_name = os.path.splitext(os.path.basename(ligand_path))[0]
      save_inf.append(ligand_name)
      for protein in pr_list :
        coor_1hsg= pdb_manip.Coor()
        coor_1hsg.get_PDB(f'{protein}', f'data/{protein}.pdb')

        # Keep only the amino acids
        rec_coor= coor_1hsg.select_part_dict(selec_dict={'res_name': pdb_manip.PROTEIN_RES})
        rec_coor.write_pdb(f'data/rec_{protein}.pdb')

        test_dock = docking.Docking('test', lig_pdb= f"{full_path}/{ligand_path}", rec_pdb= f'data/rec_{protein}.pdb')
        test_dock.prepare_ligand()
        test_dock.prepare_receptor()
        try :
          test_dock.run_docking(out_pdb=f'results/{ligand_name}_{protein}.pdb',
                              num_modes=3,
                              energy_range=5,
                              exhaustiveness=8,
                              dock_bin='smina')
          #print(test_dock.affinity)
          save_inf.append(test_dock.affinity[1]["affinity"])
        except:
          save_inf.append(0)

      df.loc[len(df.index)] = save_inf
  # Remove the lock file when done
  os.remove(lock_file)
  df.to_csv(f"{full_path}/{index_path}.csv")

