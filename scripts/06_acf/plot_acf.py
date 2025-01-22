import numpy as np
import matplotlib.pyplot as plt

protein = ["1","3a","5a","8a"]
colors = ['red', 'orange', 'cyan', 'purple']
tau_arr = np.arange(0, 75000, 1)
c = -1
fig, ax=plt.subplots(figsize=(7,5), layout='tight')
for p in protein:
    c = c + 1
    total = np.load(f"output/wnt{p}_acf_t12_numpy.npy")
    total2 = []
    for column in range(total.shape[1]):
        if total[0,column] > 0.9 and total[1,column] <1.1:
            total2.append(total[:,column])

    total2 = np.array(total2)

    avg_total2 = np.average(total2,axis=0)
    std_total2 = np.std(total2,axis=0)
    bool = True
    inflection_pt = 0
    for val in range(len(avg_total2) - 1):
        if avg_total2[val] > 0 and avg_total2[val+1] < 0 and bool is True:
            bool = False
            inflection_pt = val
    print(f" Protein {p} is closest to zero at {inflection_pt}")

    plt.plot(avg_total2, label = f'Wnt{p} - {inflection_pt} frames', color=colors[c])
    plt.fill_between(tau_arr, avg_total2-std_total2, avg_total2+std_total2, alpha=0.1, color=colors[c])
    plt.axvline(inflection_pt, linestyle='dashed', color=colors[c])
    plt.xlabel(f"$\\tau$ (frames)", fontsize=16)
    plt.ylabel(f"C ($\\tau$) ", fontsize=16)
    plt.title("Average Contact Distance ACF (over 985 features)", fontsize=16)
    plt.legend(bbox_to_anchor=(0.6, 0.5), loc='center left')

plt.savefig(f"corr_acf_t12_numpy.jpg", dpi=400)