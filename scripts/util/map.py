import numpy as np

def delete_resnames(arr):
    parsed_arr = np.zeros((len(arr),2))
    for i in range(len(arr)):
        splitstr = arr[i].split('_')
        parsed_arr[i,0] = int(splitstr[0][3:])
        parsed_arr[i,1] = int(splitstr[1][3:])
    return parsed_arr

def arr_to_str(arr):
    parsed_arr = np.empty(len(arr), dtype=object)
    for i in range(len(arr)):
        parsed_arr[i] = f"{str(int(arr[i, 0]))}_{str(int(arr[i, 1]))}"
    return parsed_arr

def idxlist_to_vmd(arr):
    selstr = "segname PROA and ("
    for i in arr:
        selstr += f"resid {i} or "
    selstr = selstr[:-4] + ")"
    return selstr

def map_onto_ref(df, fit_from_arr, fitstr, refstr):
    map_arr = np.zeros((len(fit_from_arr), 2))
    for i in range(len(fit_from_arr)):
        curridx = fit_from_arr[i,0]
        map_arr[i,0] = df[df[fitstr] == curridx][refstr].item()
        map_arr[i,1] = fit_from_arr[i,1]
    return map_arr