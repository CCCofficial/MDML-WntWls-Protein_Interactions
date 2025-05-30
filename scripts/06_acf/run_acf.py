from statsmodels.tsa.stattools import acovf
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from joblib import Parallel, delayed


def corr_final(x, tau):
    xo = x - np.mean(x, axis=0)
    output = []
    for dtau in tau:
        if dtau == 0:
            numer = np.dot(xo, xo)
            tot = np.sum(numer) / len(xo)
        else:
            numer = np.dot(xo[dtau:], xo[:-(dtau)])
            tot = np.sum(numer) / len(xo)
        output.append(tot)
    return np.array(output)


def corr_final_norm(x, tau):
    xo = x - np.mean(x, axis=0)
    output = []
    for dtau in tau:
        if dtau == 0:
            tot = 1
        else:
            numer = np.dot(xo[(dtau):], xo[:-(dtau)])
            denom = np.dot(xo, xo)
            tot = np.sum(numer) / np.sum(denom)
        output.append(tot)
    return np.array(output)


def corr_final_matt(x, tau):
    xo = x - np.mean(x, axis=0)
    output = []
    taumax = tau[-1]
    for dtau in tau:
        if dtau == 0:
            tot = 1
        else:
            start = dtau
            end = taumax + dtau
            numer = np.dot(xo[start:end], xo[:taumax])
            denom = np.dot(xo, xo)
            tot = np.sum(numer) / np.sum(denom)
        output.append(tot)
    return np.array(output)


def corr_final_matt_slide(x, tau):
    ## use sliding window method
    xo = x - np.mean(x, axis=0)
    n_sample = xo.shape[0]

    taumax = tau[-1]
    win_size = taumax +1 # size of sliding window, also the number of taus to store
    iter_end = n_sample - win_size # the first idx of the last window
    acf_array = np.zeros((iter_end+1,win_size)) # row for each window, column for each tau

    for i in range(iter_end+1):
        xo_window = xo[i:i+win_size]
        dot_product = np.dot(xo_window,xo_window[0])
        acf_array[i] = dot_product

    acf_final = np.sum(acf_array,axis=0) # sum over all windows
    acf_final = acf_final / acf_final[0]
    return acf_final

# Define function to process one feature
def process_feature_numpy(feature, data):
    numpy_result = acovf(data[:, feature], adjusted=False, demean=True, fft=True, missing="none", nlag=total_tau-1)
    numpy_result_norm = numpy_result / numpy_result[0]
    #acf_final_sliding = corr_final_matt_slide(data[:, feature], window_arr)
    
    return numpy_result_norm  #, acf_final_sliding# Define function to process one feature

# Define function to process one feature
def process_feature_window(feature, data, window_arr):
#    numpy_result = acovf(data[:, feature], adjusted=False, demean=True, fft=True, missing="none", nlag=total_tau-1)
#    numpy_result_norm = numpy_result / numpy_result[0]
    acf_final_sliding = corr_final_matt_slide(data[:, feature], window_arr)
    
    return acf_final_sliding  # Define function to process one feature
    




contact = np.array(pd.read_csv(f"../05_combine/output/ML_input.csv"))

protein_names = ["1", "3a", "5a", "8a"]
nFrames = 75000
total_tau = 35000
tau_arr = np.arange(0, total_tau, 1)
window_arr = np.arange(0, total_tau, 1)
for i in range(len(protein_names)):
    print(f"--------------\nLooking at protien {protein_names[i]}\n-------------------")
    data = contact[nFrames * i:nFrames * (i + 1)]
    total_numpy = np.zeros((len(tau_arr), data.shape[1]))
#    total_matt = np.zeros((len(window_arr), data.shape[1]))
    total_slidewindow = np.zeros((len(window_arr), data.shape[1]))

    results_numpy = Parallel(n_jobs=50)(
        delayed(process_feature_numpy)(feature, data) for feature in range(data.shape[1])
    )
    results_slidewindow = Parallel(n_jobs=50)(
        delayed(process_feature_window)(feature, data, window_arr) for feature in range(data.shape[1])
    )
    
    # Convert results to NumPy array and transpose
    total_numpy = np.array(results_numpy).T


    total_slidewindow = np.array(results_slidewindow).T

#    for feature in range(data.shape[1]):
#        print(f"Iteration {feature}")
#        numpy_result = acovf(data[:, feature], adjusted=False, demean=True, fft=True, missing="none", nlag=total_tau-1)
#        numpy_result_norm = numpy_result / numpy_result[0]
##        my_result_2 = corr_final_matt(data[:, feature], window_arr)
##        my_result_2_norm = my_result_2 / my_result_2[0]
#        acf_final_sliding = corr_final_matt_slide(data[:,feature],window_arr)
#
#        print(f"acf at t=0 for acovf: {numpy_result[0]} \n For sliding window: {acf_final_sliding[0]}")
#        total_numpy[:, feature] = numpy_result_norm
##        total_matt[:, feature] = my_result_2_norm
#        total_slidewindow[:,feature] = acf_final_sliding
    np.save(f"output/wnt{protein_names[i]}_acf_t12_numpy.npy", total_numpy)
#    np.save(f"output/wnt{protein_names[i]}_acf_t12_matt.npy", total_matt)
    np.save(f"output/wnt{protein_names[i]}_acf_t12_slidewindow.npy", total_slidewindow)


