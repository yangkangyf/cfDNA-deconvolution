from matplotlib import pyplot as plt
import matplotlib.cm as cm
import numpy as np

pregn_cfdna_list = ["preg_cfdna_" + str(i+1) + ".txt" for i in range(5)]

meth_vals = []
for file in pregn_cfdna_list:
    with open(file) as f:
        meth_vals.append([float(next(f).split()[4]) for x in range(100000)])

meth_vals_rounded = []
for sublist in meth_vals:
    meth_vals_rounded.append([round(elem, 3) for elem in sublist])


colors = cm.rainbow(np.linspace(0, 1, len(meth_vals_rounded)))
for i, c in zip(range(len(meth_vals_rounded)), colors):
    plt.xlim(0, 1000)
    plt.scatter(range(len(meth_vals_rounded[i])), meth_vals_rounded[i], color=c, label="pregnant cfDNA: " + str(i))
    plt.xlabel("CpG site")
    plt.ylabel("Proportion methylated")
# plt.legend(loc=0)
plt.savefig("ex.ps", format="ps", dpi=1000)
