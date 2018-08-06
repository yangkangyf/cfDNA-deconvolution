from optimization.quadratic_programming_optimization import perform_optimization as qp
from optimization.naive_optimization import perform_optimization as naive

from utils.run_utils import *


def generate_optimization(reference, methylated, unmethylated, method):

    proportions_est = np.zeros((1, reference.shape[0])) + 0.5  # start with random estimation of proportions
    proportions_est = proportions_est / (np.sum(proportions_est))
    # print(proportions_est.shape)
    return method(proportions_est, reference, methylated, unmethylated)


if __name__ == "__main__":


    for patient in ["deconv_exp1_mix3.txt"]:
        reference, methylated, unmethylated = generate_matrices("/Users/Christa.Caggiano/Desktop/zaitlen_lab_desktop/" + patient, "data/exp_tissues.txt")
        np.savetxt("reference.txt", reference)
        x = generate_optimization(reference, methylated, unmethylated, qp)
        print(x)
        np.savetxt("output/" + patient + "_qp", x)