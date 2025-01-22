from MDAnalysis.analysis import distances
import numpy as np


def get_res_list(WNT_atoms, WNTLESS_atoms, idx):
    res_list = np.empty((len(WNT_atoms.resids), len(WNTLESS_atoms.resids)),
                        dtype=object)  # store the index names of the distance matrix
    for wnt_name in range(len(WNT_atoms.resids)):
        for wntless_name in range(len(WNTLESS_atoms.resids)):
            res_list[
                wnt_name, wntless_name] = f"{WNT_atoms.resnames[wnt_name]}{WNT_atoms.resids[wnt_name]}_{WNTLESS_atoms.resnames[wntless_name]}{WNTLESS_atoms.resids[wntless_name]}"
    res_list = res_list.reshape(len(WNT_atoms.resids) * len(WNTLESS_atoms.resids))
    if idx is None:
        return res_list
    else:
        return res_list[idx]


def get_frames(traj_length, num_proteins, iter, equi):
    inv_length = int(traj_length)
    if equi is True:
        start = int((inv_length / 2) + (
                    inv_length * iter))  # Only consider second half of the trajectory for determining whether atom pair is within threshhold
        end = int((inv_length) + (inv_length * iter))
    if equi is False:
        start = int(
            inv_length * iter)  # Only consider second half of the trajectory for determining whether atom pair is within threshhold
        end = int((inv_length) + (inv_length * iter))
    return start, end


def init_distance_matrix(WNT_atoms, WNTLESS_atoms, num_frames, num_proteins, stride, idx, equi):
    if idx is None:
        if equi is True:
            return np.empty((int(np.ceil((num_frames / 2) * num_proteins / stride)),
                             len(WNT_atoms) * len(WNTLESS_atoms)))  # distance matrix
        else:
            return np.empty((int(np.ceil(num_frames * num_proteins / stride)),
                             len(WNT_atoms) * len(WNTLESS_atoms)))  # distance matrix
    else:
        return np.empty((int(np.ceil(num_frames * num_proteins / stride)), len(idx)))


def get_distance_matrix(u, WNT_atoms, WNTLESS_atoms, idx):
    # NOTE: MINIMUM IMAGE CONVENTION NEEDS TO BE APPLIED!!!!!!!
    if idx is None:
        return distances.distance_array(WNT_atoms.positions, WNTLESS_atoms.positions, box=u.dimensions).reshape(
            len(WNT_atoms) * len(WNTLESS_atoms))
    else:
        return distances.distance_array(WNT_atoms.positions, WNTLESS_atoms.positions, box=u.dimensions).reshape(
            len(WNT_atoms) * len(WNTLESS_atoms))[idx]  # compute the distance matrix at the current timestep


def check_input(traj, num_frames, num_proteins):
    for iteration in range(num_proteins):
        if len(traj[iteration].trajectory) < num_frames:
            raise Exception("Trying to read more frames than exist in trajectory. Please try again.")


def gather_matrix(u, WNT_atoms, WNTLESS_atoms, nFrames, nProteins, stride, pair_indeces=None, equilibration=False):
    """
    Function used to calculate the distance matrix between two AtomGroups
    u -> MDAnalysis Trajectory
    WNT_atoms -> AtomGroup 1
    WNTLESS_atoms -> AtomGroup 2
    nFrames -> The number of frames that we want to read from EACH trajectory
    nProteins -> The number of trajectories that we want ot read from
    stride -> How often (in frames) we want to calculate the distance matrix
    """

    # Edge case, if the number of frames is longer than the length of the trajectory, throw exception and redo
    check_input(u, nFrames, nProteins)


    # How we want to store the data
    res_list = get_res_list(WNT_atoms[0], WNTLESS_atoms[0], pair_indeces)
    frames = init_distance_matrix(WNT_atoms[0], WNTLESS_atoms[0], nFrames, nProteins, stride, pair_indeces,
                                  equilibration)

    i = -1  # keep a counter
    for iteration in range(nProteins):
        start, end = get_frames(nFrames, nProteins, 0, equilibration)

        # We subset the trajectory based on tFinal
        for _ in u[iteration].trajectory[start:end:stride]:
            i = i + 1  # update the counter
            frames[i] = get_distance_matrix(u[iteration], WNT_atoms[iteration], WNTLESS_atoms[iteration], pair_indeces)

    return frames, res_list


def identify_contacts(frames, res_list, distance_threshhold):
    """
    Function used to figure out what contact distances are within the threshhold

    nFrames
    """
    contacts_within_threshhold = ""  # string storing contact pairs we are keeping
    index_within_threshhold = ""  # string storing the indeces that we are keeping

    for contact_pair in range(frames.shape[1]):
        # For a given pair of atoms, if there is at least one frame where the contact distance is less
        # than the short range LJ cutoff (12 Angstroms), consider this in the input set
        #
        # Equilibration Time: Last 750 ns of simulation.

        if np.min(frames[:, contact_pair]) < distance_threshhold:
            contacts_within_threshhold += res_list[contact_pair] + "\n"
            index_within_threshhold += str(contact_pair) + "\n"
    pair_indeces = np.array(index_within_threshhold.split("\n")[:-1], dtype=int)
    pair_names = np.array([pair for pair in contacts_within_threshhold.split("\n")[:-1]])
    return pair_indeces, pair_names
