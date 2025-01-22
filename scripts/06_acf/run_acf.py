from statsmodels.tsa.stattools import acovf
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def corr_final(x, tau):
    xo = x - np.mean(x, axis=0)
    output = []
    for dtau in tau:
        if dtau == 0:
            numer = np.dot(xo, xo)
            tot = np.sum(numer) / len(xo)
        else:
            numer = np.dot(xo[dtau + 1:], xo[:-(dtau + 1)])
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
            numer = np.dot(xo[(dtau + 1):], xo[:-(dtau + 1)])
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
            start = dtau + 1
            end = taumax + dtau + 1
            numer = np.dot(xo[start:end], xo[:taumax])
            denom = np.dot(xo, xo)
            tot = np.sum(numer) / np.sum(denom)
        output.append(tot)
    return np.array(output)


contact = np.array(pd.read_csv(f"../05_combine/output/ML_input.csv"))

protein_names = ["1", "3a", "5a", "8a"]
nFrames = 75000
tau_arr = np.arange(0, 75000, 1)
window_arr = np.arange(0, 35000, 1)
for i in range(len(protein_names)):
    print(f"--------------\nLooking at protien {protein_names[i]}\n-------------------")
    data = contact[nFrames * i:nFrames * (i + 1)]
    total_numpy = np.zeros((len(tau_arr), data.shape[1]))
    total_matt = np.zeros((len(window_arr), data.shape[1]))
    for feature in range(10):
        print(f"Iteration {feature}")
        numpy_result = acovf(data[:, feature], adjusted=False, demean=True, fft=False, missing="none", nlag=74999)
        numpy_result_norm = numpy_result / numpy_result[0]
        my_result_2 = corr_final_matt(data[:, feature], window_arr)
        my_result_2_norm = my_result_2 / my_result_2[0]
        total_numpy[:, feature] = numpy_result_norm
        total_matt[:, feature] = my_result_2_norm
    np.save(f"output/wnt{protein_names[i]}_acf_t12_numpy.npy", total_numpy)
    np.save(f"output/wnt{protein_names[i]}_acf_t12_matt.npy", total_matt)