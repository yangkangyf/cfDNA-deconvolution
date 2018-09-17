from astropy.stats import jackknife_stats
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

jackknife_estimate = np.loadtxt("/Users/Christa.Caggiano/Documents/UCSF_year2/research/cfDNA_deconvolution/output/jackknife_estimate_test.txt")

all_estimate = np.loadtxt("/Users/Christa.Caggiano/Documents/UCSF_year2/research/cfDNA_deconvolution/output/ctrl2_merged_new_ref.txt_qp")

error = np.zeros(all_estimate.shape[0])

for row in jackknife_estimate:
	error += np.square(row - all_estimate)

overall_error = np.sqrt(error*((error.shape[0]-1)/error.shape[0]))

x_pos = np.arange(len(overall_error))
x_labels = pd.read_table("/Users/Christa.Caggiano/Documents/UCSF_year2/research/cfDNA_deconvolution/data/tissues-new-ref-2.txt", header=None)

fig, ax = plt.subplots()
ax.bar(x_pos, all_estimate, yerr=overall_error, align='center', alpha=0.5, ecolor='black', capsize=10)
ax.set_ylabel("Proportion of tissue")
ax.set_xticks(x_pos)
ax.set_xticklabels(x_labels.ix[:, 0])
plt.xticks(rotation=90)
ax.yaxis.grid(True)

# Save the figure and show
plt.tight_layout()
plt.savefig('bar_plot_with_error_bars.png')
plt.show()

x_labels.ix[:, 0]
x_pos