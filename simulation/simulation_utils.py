import pandas as pd
import random
import numpy as np


def generate_proportion(individuals, tissue):
    """
    generates the contribution of different tissues for
    each individual

    :param tissue: number of tissues to be simulated
    :param individuals: number of individuals
    :return: pandas dataframe of proportion of tissue, of shape individuals by tissue
    """
    rows = []
    for i in range(individuals):
        vals = np.random.choice(10, tissue)  # pick tissue number of values from 1 to 10 to be proportions
        rows.append([x/sum(vals) for x in vals])  # proportions must sum to 1

    return pd.DataFrame.from_records(rows)


def generate_proportion_fixed(individuals, tissue, percentage):

    rows = []
    for i in range(individuals):
        vals = np.random.multinomial(100-percentage, np.ones(tissue-1)/(tissue-1), size=1)[0].tolist()
        vals.append(percentage)
        rows.append([x/100 for x in vals])  # proportions must sum to 1

    return rows


def generate_reference(tissue, sites):
    """
    generate reference methylation values for each tissue

    :param tissue: number of tissues
    :param sites: number of CpG sites
    :return: pandas dataframe of proportion of methylated for each site
    """
    rows = []
    for t in range(tissue):  # generate for each tissue
        vals = []
        for s in range(sites):
            vals.append(random.uniform(0, 1))
        rows.append(vals)
    return pd.DataFrame(rows)


def generate_depth(sites, individuals, read_depth):

    return np.random.poisson(lam=read_depth, size=(individuals, sites))


def generate_counts(depths, observed, sites, individuals):

    counts = np.zeros((individuals, sites))

    for i in range(individuals):
        counts[i, ] = np.random.binomial(depths[i, ], observed[i, ], size=sites)

    return counts


def round_to_one(array):
    for i in np.nditer(array, op_flags=['readwrite']):
        if i > 1:
            i[...] = 1
        if i < 0:
            i[...] = 0
    return array
