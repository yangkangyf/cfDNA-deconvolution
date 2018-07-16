# naive optimization that constrains by summing row and divide by sum
# May 2018
# author <christacaggiano@ucsf.edu>

# imports

import numpy as np
from scipy.optimize import minimize
from sklearn.metrics import mean_squared_error
from scipy.stats import binom


def log_likelihood(proportions, observed, reference):
    """
    calculate log likelihood for optimization

    :param proportions: current estimation of proportions
    :param observed: observed methylation states
    :param reference: reference methylation by tissue
    :return: negative log likelihood
    """

    b = np.transpose(np.matmul(proportions, reference))
    sigma = np.var(observed-b)
    N = len(proportions)


    # negative log likelihood function, copied from n. zaitlen's r code
    ll = np.sum(-N/2*np.log(2*np.pi) - N/2 * np.log(sigma) - 1/(2*sigma**2) * (np.sum((observed - b)**2)))
    return ll


def counts_log_likelihood(proportions, methylated, unmethylated, reference):
    b = np.matmul(proportions, reference)

    ll = np.sum(binom.logpmf(methylated, methylated+unmethylated, b, loc=0))
    # print(ll)
    return -ll/1000


def perform_optimization(proportions_est, reference, methylated, unmethylated):
    """
    performs optimization using quadratic programming

    :param proportions_est: estimation of proportions
    :param proportions: 'true' proportions
    :param observed: methylation patterns for cpgs
    :param reference: reference methylation patterns for all tissues
    :return: truth, guess
    """
    bounds = tuple([0, 1] for x in range(np.shape(proportions_est)[1]))
    cons = ({'type': 'eq', 'fun': lambda x: 1 - sum(x)})

    prop_guess = minimize(counts_log_likelihood, proportions_est, args=(methylated, unmethylated, reference),
                            bounds=bounds, constraints=cons, method="SLSQP", options={'maxiter': 10000, 'ftol': 1e-010})

    # print(prop_guess)
    return prop_guess["x"]

