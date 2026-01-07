import numpy as np
import pandas as pd
import os

file_path = os.path.dirname(os.path.realpath(__file__))


def get_A(filename):
    """Read adjacency matrix of a graph and cluster size

    Args:
        filename (_type_): _description_

    Returns:
        raw_data, adjacency_matrix, k
    """
    A_raw = np.loadtxt(filename, dtype=int, usecols=(0, 1))
    n, k = A_raw[0, :]
    
    A = np.zeros(shape=(n,n))
    A[A_raw[1:,0]-1, A_raw[1:,1]-1] = 1
    A = A+A.transpose()

    return A_raw, A, k

def get_solutions():
    """Read solutions.

    Returns:
        pandas dataframe with solutions
    """
    with open(f'{file_path}/../solutions.html') as f:
        solutions = pd.read_html(f, header=0)[0]
    
    old_columns = ['Instance', 'optimum', 'Instance.1', 'optimum.1', 'Instance.2',
           'optimum.2', 'Instance.3', 'optimum.3', 'Instance.4', 'optimum.4',
           'Instance.5', 'optimum.5']
    
    new_columns = ['instance_40', 'optimum_40', 'instance_80', 'optimum_80', 'instance_100',
           'optimum_100', 'instance_120', 'optimum_120', 'instance_140', 'optimum_140',
           'instance_160', 'optimum_160']

    return solutions.rename(columns=dict(zip(old_columns, new_columns)))