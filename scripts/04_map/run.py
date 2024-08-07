
import numpy as np
import pandas as  pd
import sys
import argparse
from scripts.util.map import *

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description='''
-------------------------------------------------------------------------------------------------
Map the residue indeces of a wnt system onto wnt8a system. Default RMSD cutoff of 4.5 Angstroms is used.
If a best mapping is above the RMSD cutoff, a residue value of -1 is output - this represents no fit.
---------------------------------------------------------------------------------------------------------------------
Example Command: python3 run_sequences.py -sys 5a -cutoff 12 -dir .
---------------------------------------------------------------------------------------------------------------------''')
parser.add_argument('-sys', dest = 'SYSTEM', help='Wnt system directory (from box)', type=str, required=True)
parser.add_argument('-cutoff', dest = 'distance_threshhold', help='Cutoff for contacts (in angstroms)', type=int, required=True)
parser.add_argument('-dir', dest = 'workingDir', help='Working directory - used to store results', type=str, required=True)
args = parser.parse_args()

SYSTEM= args.SYSTEM
distance_threshhold = args.distance_threshhold
currDir = args.workingDir

# Import in the labels
fitwnt_resids = np.loadtxt(f'{currDir}/03_finalize_dataset/output/WNT{SYSTEM}_threshhold{distance_threshhold}_labels.txt',dtype=str)

# Res idx hashmap
hashmap_df = pd.read_csv(f"{currDir}/00_map/output/{SYSTEM}_to_8a_map.csv")
print(f"There are {fitwnt_resids.shape[0]} contact pairs in Wnt{SYSTEM}.")

# Go through the data files - figure out what the unique wnt residues are
parsed_fitwnt = delete_resnames(fitwnt_resids)

converted_fit_into_8a = map_onto_ref(hashmap_df, parsed_fitwnt, SYSTEM, "8a")
converted_fit_into_8a_str = arr_to_str(converted_fit_into_8a)

print(f"After conversion, there are {converted_fit_into_8a_str.shape[0]} contact pairs in Wnt{SYSTEM}.")
print(f"{np.unique(converted_fit_into_8a[:,0]).shape[0]} unique contact pairs in Wnt{SYSTEM} and {np.unique(converted_fit_into_8a[:,1]).shape[0]} unique pairs on WntLess.")

np.savetxt(f"{currDir}/04_map/output/WNT{SYSTEM}_threshhold{distance_threshhold}_labels.txt", converted_fit_into_8a_str, fmt='%s')

