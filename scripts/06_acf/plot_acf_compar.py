import numpy as np
import matplotlib.pyplot as plt

protein = ["1","3a","5a","8a"]
colors = ['red', 'orange', 'cyan', 'purple']
tau_arr = np.arange(0, 35000, 1)


npy_files = ['numpy','matt','slidewindow']

fig, axes=plt.subplots(3,1,figsize=(8,16), layout='tight')

for i,npy in enumerate(npy_files):
    ax = axes[i]
    c = -1
    for p in protein:
        c = c + 1
        total = np.load(f"output/wnt{p}_acf_t12_{npy}.npy")

        if npy =="matt": # only calculate the first 10 features
            total = total[:,:10]
        avg_total = np.average(total,axis=1)
        std_total = np.std(total,axis=1)

        inflection_pt = 0
        for val in range(len(avg_total) - 1):
            if avg_total[val] > 0 and avg_total[val+1] < 0:
                inflection_pt = val
                break
        print(f" Protein {p} is closest to zero at {inflection_pt}")

        ax.plot(avg_total, label = f'Wnt{p} - {inflection_pt} frames', color=colors[c])
        ax.fill_between(tau_arr, avg_total-std_total, avg_total+std_total, alpha=0.1, color=colors[c])
        ax.axvline(inflection_pt, linestyle='dashed', color=colors[c])
        ax.set_xlabel(f"$\\tau$ (frames)", fontsize=14)
        ax.set_ylabel(f"C ($\\tau$) ", fontsize=14)
        ax.set_title(f"Method {npy}: Average Contact Distance ACF (over 1153 features)", fontsize=14)
        ax.legend(bbox_to_anchor=(0.6, 0.5), loc='center left')


plt.savefig(f"corr_acf_t12_compar.jpg", dpi=400)


# Now randomly chose one feature and plot and compare
# less than 10
feature_idx = 8

fig, axes=plt.subplots(3,1,figsize=(7,15), layout='tight')

for i,npy in enumerate(npy_files):
    ax = axes[i]
    c = -1
    for p in protein:
        c = c + 1
        total = np.load(f"output/wnt{p}_acf_t12_{npy}.npy")

        avg_total=total[:,feature_idx]


#        std_total = np.std(total,axis=1)

        inflection_pt = 0
        for val in range(len(avg_total) - 1):
            if avg_total[val] > 0 and avg_total[val+1] < 0:
                inflection_pt = val
                break
        print(f" Protein {p} feature {feature_idx} is closest to zero at {inflection_pt}")

        ax.plot(avg_total, label = f'Wnt{p} - {inflection_pt} frames', color=colors[c])
#        ax.fill_between(tau_arr, avg_total-std_total, avg_total+std_total, alpha=0.1, color=colors[c])
        ax.axvline(inflection_pt, linestyle='dashed', color=colors[c])
        ax.set_xlabel(f"$\\tau$ (frames)", fontsize=14)
        ax.set_ylabel(f"C ($\\tau$) ", fontsize=14)
        ax.set_title(f"Method {npy}: Distance ACF for contact {feature_idx}", fontsize=14)
        ax.legend(bbox_to_anchor=(0.6, 0.5), loc='center left')


plt.savefig(f"corr_acf_t12_feature{feature_idx}_compar.jpg", dpi=400)
