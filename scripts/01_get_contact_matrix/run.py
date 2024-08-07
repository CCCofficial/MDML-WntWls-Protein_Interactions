import sys
from scripts.util.distances import *
import MDAnalysis as mda
import numpy as np
import os
import argparse

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description='''
-------------------------------------------------------------------------------------------------
Generate distance contact matrices for the wnt system specified. Matrices are output in chunks, same format as 
input trajectories
---------------------------------------------------------------------------------------------------------------------
Example Command: python3 run.py -sys 5a -dir . -dt 20
---------------------------------------------------------------------------------------------------------------------''')
parser.add_argument('-sys', dest = 'SYSTEM', help='Wnt system directory (from box)', type=str, required=True)
parser.add_argument('-dir', dest = 'workingDir', help='Working directory - used to store results', type=str, required=True)
parser.add_argument('-dt', dest = 'timestep', help='saving frequency (in ps)', type=int, required=True)
args = parser.parse_args()

proteinnames = [args.SYSTEM]
store_freq = args.timestep
workingDir = args.workingDir

iter = ["1"]
n = [0]
nProteins = len(proteinnames)

dcd_arr = []
psf_arr = []
# Get the input paths

for option in range(len(proteinnames)):
    data_dir = f"/Users/masauer2/Library/CloudStorage/Box-Box/Summerinternship_2024/WNT-WLS-project/Data_Backup/wnt{proteinnames[option]}/copy0{iter[option]}"
    psf_arr.append(f"{data_dir}/Wnt{proteinnames[option]}WlsPc_copy_0{iter[option]}.psf")
    dcd_dir = f"{data_dir}/dcd_files"
    dcd_names = sorted(os.listdir(dcd_dir))
    n.append(len(dcd_names))
    dcd_arr.append([f"{dcd_dir}/{name}" for name in dcd_names])

dcd_arr = [x for xs in dcd_arr for x in xs]
print(dcd_arr)

for filenum in range(len(dcd_arr)):
    print(f"Reading trajectory at iteration {filenum}.")
    # Selection choices
    WNT_CALPHA = "(segid PROA and name CA)"
    WNTLESS_CALPHA = "(segid PROB and name CA)"
    SYSTEM_CALPHA = "(segid PROA and name CA) or (segid PROB and name CA)"

    # A list that holds all dcd filenames
    dcd_list = [dcd_arr[filenum]]

    # Current universe containing all dcd files
    u = [mda.Universe(psf_arr[0], dcd_list)]

    WNT_atoms = [c.select_atoms(WNT_CALPHA) for c in u]  # AtomGroup associated w/ WNT
    WNTLESS_atoms = [c.select_atoms(WNTLESS_CALPHA) for c in u]  # AtomGroup associated w/ WNTLESS

    nFrames = len(u[0].trajectory)

    # Step 1 - get the distance matrix with a high stride (10)
    stack, res_list = gather_matrix(u, WNT_atoms, WNTLESS_atoms, nFrames, nProteins, 1, pair_indeces=None, equilibration=False)

    print(f"From {stack.shape[0]} frames ({store_freq * stack.shape[0] / 1000} ns) we have {stack.shape[1]} input features!")


    # Output the distance matrix
    # Reshape such that the array shape is (number_of_frames x number_of_contact_pairs)
    np.save(f"{workingDir}/01_get_contact_matrix/output/WNT{proteinnames[0]}_distances_iter{filenum}.npy", stack)

# Output the residue list
# Reshape such that the array shape is (number_of_frames x number_of_contact_pairs)
np.savetxt(f"{workingDir}/01_get_contact_matrix/output/WNT{proteinnames[0]}_resids.txt", res_list, fmt="%s")
