# naive optimization that constrains by summing row and divide by sum
# May 2018
# author <christacaggiano@ucsf.edu>

# imports
import numpy as np
from scipy.optimize import minimize
from scipy.stats import binom
from sklearn.metrics import mean_squared_error


def log_likelihood(proportions, observed, reference):
    """
    calculate log likelihood for optimization

    :param proportions: estimation of proportions
    :param observed: observed methylation states
    :param reference: reference methylation by tissue
    :return: negative log likelihood
    """

    proportions = proportions/np.sum(proportions)  # manually dividing by sum of rows to constrain
    b = np.transpose(np.matmul(proportions, reference))
    sigma = np.var(observed-b)
    N = len(proportions)

    # negative log likelihood function, copied from n. zaitlen's r code
    neg_ll = np.sum(-N/2*np.log(2*np.pi) - N/2 * np.log(sigma) - 1/(2*sigma**2) * (np.sum((observed - b)**2)))

    return neg_ll


def counts_log_likelihood(proportions, methylated, unmethylated, reference):

    proportions = proportions/(np.sum(proportions))

    b = np.matmul(proportions, reference)

    ll = np.sum(binom.logpmf(methylated, methylated+unmethylated, b, loc=0))

    return -ll


def perform_optimization(proportions_est, reference, methylated, unmethylated):
    """
    use scipy optimize to perform minimization using BFGS

    :param individuals: number of individuals (always pretty much 1)
    :param tissues: number of tissues
    :param proportions_est:
    :return: truth, guess
    """

    bounds = tuple((0, 1) for x in range(np.shape(proportions_est)[1]))  # constrains that values in array must be prop

    # perform minimization using scipy optimization, BFGS technique. Max iterations==10,000
    prop_guess = minimize(counts_log_likelihood, proportions_est, args=(methylated, unmethylated, reference),
                          bounds=bounds, method="L-BFGS-B", options={'maxiter': 10000, 'ftol': 1e-08})
    print(prop_guess)
    return prop_guess["x"]


