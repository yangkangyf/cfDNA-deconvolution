import numpy as np
from sklearn.metrics import mean_squared_error


def mse(truth, guess):

    return mean_squared_error(truth, guess)


def scaled_mse(truth, guess):

    return mean_squared_error(truth, guess)/(np.sum(np.square(truth)))


def corr(truth, guess):
    return np.corrcoef(truth, guess)[0][1]