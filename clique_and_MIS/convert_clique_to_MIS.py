import os
from scipy.io import mmread, mmwrite
from scipy.sparse import coo_matrix
import numpy as np
import glob
from pathlib import Path


def read_mtx_instance(f):
    """Read a .mtx instance using a built-in numpy function mmread"""
    return mmread(f)

def write_mtx_instance(f, M, **kwargs):
    """Write a .mtx instance using a built-in numpy function mmread"""
    return mmwrite(f, M, **kwargs)


def build_complement_instances(input_directory, output_directory, input_mask, comment=''):
    """
    Reads .mtx instances available at https://networkrepository.com/networks.php
    and creates a .txt file for the complementary graph. 

    The output file is of the form
    n m     # num. of vertices, num. of edges
    i j     # edge i, j in the complement
    ...
    """
    files = glob.glob(f'{input_directory}/{input_mask}')
    
    for i in files:
        M = read_mtx_instance(i).toarray()
        M_complement = 1-M
        np.fill_diagonal(M_complement, 0)
        result_file = i.replace(input_directory, output_directory)
        Path(os.path.dirname(result_file)).mkdir(parents=True, exist_ok=True)
        write_mtx_instance(result_file, coo_matrix(M_complement), field='pattern', symmetry='symmetric', comment=comment)

        print(f"Converted the max clique grpah {i} to: {result_file}")


# Example use
build_complement_instances('Clique', 'MIS', '*.mtx', comment='\n This is a complement of the original instance \n')
 