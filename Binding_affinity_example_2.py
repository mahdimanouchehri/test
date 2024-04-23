# -*- coding: utf-8 -*-
import time


from pdb_manip_py import pdb_manip
from docking_py import docking

# Create a Coor object
coor_1hsg = pdb_manip.Coor()
coor_1hsg.get_PDB('8x5y', 'data/8x5y.pdb')

# Select res_name MK1
lig_coor = coor_1hsg.select_part_dict(selec_dict={'res_name': ['MK1']})
print(lig_coor)
# Save the ligand coordinates
lig_coor.write_pdb('data/lig2.pdb')

# Keep only the amino acids
rec_coor = coor_1hsg.select_part_dict(selec_dict={'res_name': pdb_manip.PROTEIN_RES})
rec_coor.write_pdb('data/rec2.pdb')



test_dock = docking.Docking('test', lig_pdb='data/lig.pdb', rec_pdb='data/rec2.pdb')
test_dock.prepare_ligand()
test_dock.prepare_receptor()
# test_dock.prepare_grid(out_folder = "test" , 
#                       grid_npts =[11 , 15 ,16],
#                       center=[178 , 173 , 138])

print("start...")
start = time.time()

test_dock.run_docking(out_pdb='test_dock122.pdb',
                      num_modes=3,
                      energy_range=10,
                      exhaustiveness=3,
                      center = [178.5959, 173.27124, 138.28264],
                      grid_size = [11.634003 * 2, 15.052002 *2 , 16.612 *2],
                      dock_bin='qvinaw')

rmsd_list = test_dock.compute_dock_rmsd(test_dock.lig_pdbqt)

print(rmsd_list)
print(test_dock.affinity)

end = time.time()
print("measured time:", end - start)
