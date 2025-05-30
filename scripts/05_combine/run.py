import numpy as np
from util.map import *
import pandas as  pd
import argparse
import os

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description='''
-------------------------------------------------------------------------------------------------
Map the residue indeces of a wnt system onto wnt8a system. Default RMSD cutoff of 4.5 Angstroms is used.
If a best mapping is above the RMSD cutoff, a residue value of -1 is output - this represents no fit.
---------------------------------------------------------------------------------------------------------------------
Example Command: python3 run_sequences.py -cutoff 12 
---------------------------------------------------------------------------------------------------------------------''')

parser.add_argument('-cutoff', dest = 'distance_threshhold', help='Cutoff for contacts (in angstroms)', type=int, required=True)
parser.add_argument('-dir', dest = 'workingDir', help='Working directory - used to store results', type=str, required=True)
args = parser.parse_args()

workingDir = args.workingDir
distance_threshhold = args.distance_threshhold


# The labels that we want to merge
wnt1_labels = np.loadtxt(f"{workingDir}/04_map/output/WNT1_threshhold{distance_threshhold}_labels.txt", dtype=str)
wnt3a_labels = np.loadtxt(f"{workingDir}/04_map/output/WNT3a_threshhold{distance_threshhold}_labels.txt", dtype=str)
wnt5a_labels = np.loadtxt(f"{workingDir}/04_map/output/WNT5a_threshhold{distance_threshhold}_labels.txt", dtype=str)
ref_labels = np.loadtxt(f"{workingDir}/04_map/output/WNT8a_threshhold{distance_threshhold}_labels.txt", dtype=str) # ref = 8a

# Labels that we are going to keep
keep_these_wnt8a = []
keep_these_wnt5a = []
keep_these_wnt3a = []
keep_these_wnt1 = []

# For each feature in wnt8a (ref)
for j in range(len(ref_labels)):

    # If that feature exists in all 4 of the systems
    if str(ref_labels[j]) in wnt1_labels and str(ref_labels[j]) in wnt3a_labels and str(ref_labels[j]) in wnt5a_labels:

        # We keep the feature
        keep_these_wnt8a.append(j)
        keep_these_wnt1.append(np.where(wnt1_labels == ref_labels[j])[0][0])
        keep_these_wnt3a.append(np.where(wnt3a_labels == ref_labels[j])[0][0])
        keep_these_wnt5a.append(np.where(wnt5a_labels == ref_labels[j])[0][0])

# The data we are going to keep
wntfit_contact_matrix = np.load(f"{workingDir}/03_finalize_dataset/output/WNT8a_threshhold{distance_threshhold}_matrix.npy")
wnt1_contact_matrix = np.load(f"{workingDir}/03_finalize_dataset/output/WNT1_threshhold{distance_threshhold}_matrix.npy")
wnt3a_contact_matrix = np.load(f"{workingDir}/03_finalize_dataset/output/WNT3a_threshhold{distance_threshhold}_matrix.npy")
wnt5a_contact_matrix = np.load(f"{workingDir}/03_finalize_dataset/output/WNT5a_threshhold{distance_threshhold}_matrix.npy")

# The parsed dataframes
wnt8a_contact_matrix_parsed = pd.DataFrame(wntfit_contact_matrix[:,keep_these_wnt8a], columns=ref_labels[keep_these_wnt8a])
wnt5a_contact_matrix_parsed = pd.DataFrame(wnt5a_contact_matrix[:,keep_these_wnt5a], columns=wnt5a_labels[keep_these_wnt5a])
wnt3a_contact_matrix_parsed = pd.DataFrame(wnt3a_contact_matrix[:,keep_these_wnt3a], columns=wnt3a_labels[keep_these_wnt3a])
wnt1_contact_matrix_parsed = pd.DataFrame(wnt1_contact_matrix[:,keep_these_wnt1], columns=wnt1_labels[keep_these_wnt1])

# The final data set - the order is fixed
# 1 . 75000 frames from Wnt1
# 2. 75000 frames from Wnt3a
# 3. 75000 frames from Wnt5a
# 4. 75000 frames from Wnt8a

final = pd.concat([wnt1_contact_matrix_parsed,wnt3a_contact_matrix_parsed,wnt5a_contact_matrix_parsed,wnt8a_contact_matrix_parsed])
print(f"final input shape: {final.shape[0]} by {final.shape[1]}")

output_dir = f"{workingDir}/05_combine/output"
os.makedirs(output_dir,exist_ok=True)

final.to_csv(f"{output_dir}/ML_input.csv", index=False)
