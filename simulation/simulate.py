import matplotlib.pyplot as plt

from optimization.naive_optimization import perform_optimization as naive
from optimization.quadratic_programming_optimization import perform_optimization as qp
from simulation.simulation_utils import *
from utils.error import corr, mse, scaled_mse
from utils.run_utils import *


def generate_simulated_optimization(individuals, sites, tissues, read_depth, method, noise, fixed_proportion=None):

    proportions = np.asarray(generate_proportion_fixed(individuals, tissues, fixed_proportion))
    # proportions = generate_proportion(individuals, tissues).as_matrix()  # initialized proportions of tissue
    # reference = generate_reference(tissues, sites).as_matrix()  # cpg methylation fraction
    reference = read_roadmap_reference("/Users/Christa.Caggiano/Desktop/zaitlen_lab_desktop/roadmap_top10000_var_DMRS.txt")

    observed = np.matmul(proportions, reference)  # observed estimated is just reference times the proportions
    observed = observed + np.random.normal(0, noise, observed.shape)  # add small amounts of noise to observed
    observed = round_to_one(observed)
    depth = generate_depth(sites, individuals, read_depth)
    methylated = generate_counts(depth, observed, sites, individuals)
    unmethylated = depth - methylated

    proportions_est = np.zeros((individuals, tissues)) + 0.5  # start with random estimation of proportions
    proportions_est = proportions_est/(np.sum(proportions_est))

    return proportions, method(proportions_est, reference, methylated, unmethylated)


# global simulation parameters
individuals = 50  # number of people optimizing for
sites = 9999  # number of cpg sites
tissues = 40  # number of tissues
read_depth = 100  # read depth (methylated/unmethylated counts)
noise = 0.1

# generate_simulated_optimization(1, sites, tissues, read_depth, qp, noise)

naive_error = []
qp_error = []
read_depth_range = [10, 100, 1000, 10000]
fixed_proportion_range = [10, 1, 0.1, 0.01]

for proportion in fixed_proportion_range:
    read_depth_error = []
    for read_depth in read_depth_range:
        error = 0
        for individual in range(individuals):
            truth, guess = generate_simulated_optimization(1, sites, tissues, read_depth, qp, noise, proportion)
            error += (corr(truth[0], guess))
        read_depth_error.append(error/individuals)
        # print(read_depth_error)
    qp_error.append(read_depth_error)

plt.rcParams.update({'font.size': 12})
for proportion_error, proportion in zip(qp_error, fixed_proportion_range):
    # # # @TODO box and whisker
    plt.plot(read_depth_range, proportion_error, '-o', label=str(proportion) + "%")
    plt.xscale("log")

plt.legend()
plt.show()