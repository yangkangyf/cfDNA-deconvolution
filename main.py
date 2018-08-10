from optimization.quadratic_programming_optimization import perform_optimization as qp
from optimization.naive_optimization import perform_optimization as naive

from utils.run_utils import *


def generate_optimization(reference, methylated, unmethylated, method):

    proportions_est = np.zeros((1, reference.shape[0])) + 0.5  # start with random estimation of proportions
    proportions_est = proportions_est / (np.sum(proportions_est))
    # print(proportions_est.shape)
    return method(proportions_est, reference, methylated, unmethylated)


if __name__ == "__main__":


    for patient in ["ctrl2_merged_new_ref.txt", "ctrl3_merged_new_ref.txt", "ctrl4_merged_new_ref.txt",
                    "als2_merged_new_ref.txt", "als3_merged_new_ref.txt", "als4_merged_new_ref.txt"]:
        reference, methylated, unmethylated = generate_matrices("/Users/Christa.Caggiano/Desktop/zaitlen_lab_desktop/" + patient, "data/tissues-new-ref-2.txt")
        x = generate_optimization(reference, methylated, unmethylated, qp)
        print(x)
        np.savetxt("output/" + patient + "_qp", x)
