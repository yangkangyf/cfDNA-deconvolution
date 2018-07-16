import pandas as pd
import numpy as np


def generate_matrices(intersect_file, tissue_file):
    data = pd.read_table(intersect_file, low_memory=False)

    reference = data

    reference = reference.replace(to_replace="NA", value=np.nan)
    reference = reference.fillna(data.mean())

    # reference = reference.dropna()

    meth_array = np.round(reference["meth"].values, 0)
    unmeth_array = np.round(reference["unmeth"].values, 0)

    reference = reference.drop(["meth", "unmeth"], axis=1)

    tissues = read_tissues(tissue_file)
    average_reference = calculate_tissue_means(reference, tissues)
    return average_reference.T, np.reshape(meth_array, (1, len( meth_array))), np.reshape(unmeth_array, (1, len(unmeth_array)))


def fill_mode(df):
    for column in df:
        # print(df[column].mode()[0])
        mode = df[column].mode()[0]
        df[column].fillna(mode, inplace=True)

    return df

def read_tissues(tissue_file):
    tissues = {}
    with open(tissue_file) as f:
        for line in f:
            line = line.strip("\n").split("\t")
            tissues[line[0]] = [x for x in line[1:] if not x == ""]
    return tissues


def calculate_tissue_means(reference,  tissue_dict):
    tissue_average = pd.DataFrame()
    for tissue in tissue_dict:
        df = reference[tissue_dict[tissue]]
        df = df.astype("float64")
        tissue_average = tissue_average.assign(**{tissue:df.mean(axis=1, skipna=True)})
    print(list(tissue_average))
    return tissue_average.values


def read_roadmap_reference(file):
    roadmap = pd.read_table(file)

    return (roadmap.as_matrix()).T
