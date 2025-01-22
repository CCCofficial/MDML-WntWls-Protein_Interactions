import sys
from scripts.util.distances import *
import os
import argparse
import re

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description='''
-------------------------------------------------------------------------------------------------
Use the distance threshhold to compute all atom pairs within the threshhold.
---------------------------------------------------------------------------------------------------------------------
Example Command: python3 run.py -sys 5a -cutoff 12 -dir . -dt 20 -tf 1500000
---------------------------------------------------------------------------------------------------------------------''')
parser.add_argument('-sys', dest = 'SYSTEM', help='Wnt system directory (from box)', type=str, required=True)
parser.add_argument('-cutoff', dest = 'distance_threshhold', help='Cutoff for contacts (in angstroms)', type=int, required=True)
parser.add_argument('-dir', dest = 'workingDir', help='Working directory - used to store results', type=str, required=True)
parser.add_argument('-dt', dest = 'timestep', help='saving frequency (in ps)', type=int, required=True)
parser.add_argument('-tf', dest = 'tFinal', help='Trajectory time to read in (in ps)', type=int, required=True)

args = parser.parse_args()

SYSTEM = args.SYSTEM
distance_threshhold = args.distance_threshhold
workingDir = args.workingDir
timestep = args.timestep
tFinal =args.tFinal
nFrames = int(tFinal/timestep)

idx = np.loadtxt(f"{workingDir}/01_get_contact_matrix/output/WNT{SYSTEM}_resids.txt", dtype=str)

dcd_names = [i for i in sorted(os.listdir(f"{workingDir}/01_get_contact_matrix/output")) if f"WNT{SYSTEM}" in str(i) and ".npy" in str(i)]
idx_arr = np.array([re.findall(rf'\d+', i) for i in dcd_names])
dcd_names = list(np.array(dcd_names)[np.argsort(idx_arr[:,1].astype(int))])
print(dcd_names)

index_keeper = []
currFrame = 0

for i in range(len(dcd_names)):
    curr = np.load(f"{workingDir}/01_get_contact_matrix/output/{dcd_names[i]}")

    if currFrame < nFrames/2 and currFrame + int(curr.shape[0]) >= nFrames / 2:
        print(f"Number of frames currently read in: {currFrame} ({currFrame*20/1000} ns) -> Start at frame {currFrame + int((nFrames/2) - currFrame)} ({(currFrame + int((nFrames/2) - currFrame))*20/1000} ns)")
        pair_indeces, pair_names = identify_contacts(curr[int((nFrames/2) - currFrame):], idx, distance_threshhold);
        for j in pair_indeces:
            if j not in index_keeper:
                index_keeper.append(j)

    if currFrame >= nFrames/2 and currFrame < nFrames and currFrame + int(curr.shape[0]) < nFrames:
        print(f"Number of frames read in: {currFrame} ({currFrame*20/1000} ns)")
        pair_indeces, pair_names = identify_contacts(curr, idx, distance_threshhold);
        for j in pair_indeces:
            if j not in index_keeper:
                index_keeper.append(j)

    if currFrame >= nFrames/2 and currFrame < nFrames and currFrame + int(curr.shape[0]) >= nFrames:
        print(f"Number of frames read in: {currFrame} ({currFrame*20/1000} ns) -> End at frame {currFrame + int(nFrames - currFrame)} ({(currFrame + int(nFrames - currFrame))*20/1000} ns)")
        pair_indeces, pair_names = identify_contacts(curr[:int(nFrames - currFrame)], idx, distance_threshhold);
        for j in pair_indeces:
            if j not in index_keeper:
                index_keeper.append(j)

    currFrame += int(curr.shape[0])

np.savetxt(f"{workingDir}/02_apply_threshhold/output/WNT{SYSTEM}_idx_thresh{distance_threshhold}.txt", index_keeper, fmt="%s")

