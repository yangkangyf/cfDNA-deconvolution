import pandas as pd
import random
import numpy as np
from scipy.optimize import minimize
from sklearn.metrics import mean_squared_error
from scipy.stats import binom
import theano


def log_likelihood(proportions, observed, reference):
    """
    calculate log likelihood for optimization

    :param proportions: estimation of proportions
    :param observed: observed methylation states
    :param reference: reference methylation by tissue
    :return: log likelihood
    """

    b = np.transpose(np.matmul(proportions, reference))
    sigma = np.var(observed-b)
    N = len(proportions)

    # log likelihood function, copied from n. zaitlen's r code
    ll = np.sum(-N/2*np.log(2*np.pi) - N/2 * np.log(sigma) - 1/(2*sigma**2) * (np.sum((observed - b)**2)))

    return ll


def counts_log_likelihood(proportions, methylated, unmethylated, reference):

    b = np.matmul(proportions, reference)

    ll = np.sum(binom.logpmf(methylated, methylated+unmethylated, b, loc=0))

    return -ll


def compute_projection(proportions, z=1):

    proportions.sort()
    prop_sums = 0
    p = 0

    for j in range(1, proportions.shape[1]):
        prop_sums += sum(proportions[0][:j])
        p = proportions[0][j-1] - (1 / j)*(prop_sums - z)

        if p < 0:
            break

    for item in range(0, proportions.shape[1]):
        proportions[0][item] = max(0.0, item - p)

    return proportions


def batch_gradient_descent(likelihood, projection, learning_rate=0.01):
    gradient = theano.grad(cost=likelihood, wrt=projection)
    print(gradient)


if __name__=="__main__":
    batch_gradient_descent(28, [1, 2, 3], 0.01)









