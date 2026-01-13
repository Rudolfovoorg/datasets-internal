import os
from scipy.io import mmread
import numpy as np

def extract_nonzero_pairs(A):
    """
    Extracts nonzero pairs of a matrix A. These pairs are edges in the corresponding graph.
    Since A is symmetric, we check nonzero elements only in the upper triangle of A.

    The numpy/python indices (i, j) are converted to indices (int(i) + 1, int(j) + 1) so that
    they are aligned with the classical vertex numbering system {1, 2, ..., n}. 
    """
    rows, cols = np.nonzero(np.triu(A, 1)) # Since A is symmetric, we want to have each edge counter only once
    pairs = list(zip(rows, cols))
    corrected_pairs = [(int(i) + 1, int(j) + 1) for (i, j) in pairs]
    return corrected_pairs


def read_mtx_instance(f):
    """Read a .mtx instance using a built-in numpy function mmread"""
    M = mmread(f).toarray()
    return M


def build_complement_instances(input_directory, output_directory):
    """
    Reads .mtx instances available at https://networkrepository.com/networks.php
    and creates a .txt file for the complementary graph. 

    The output file is of the form
    n m     # num. of vertices, num. of edges
    i j     # edge i, j in the complement
    ...
    """
    all_files = sorted(os.listdir(input_directory))
    for instance_name in all_files:
        instance_file = os.path.join(input_directory, instance_name)
        print(f"Converting file: {instance_name}")

        result_file = os.path.join(output_directory, instance_name.replace(".mtx", ".txt"))

        M = read_mtx_instance(instance_file)

        M_complement = 1-M
        np.fill_diagonal(M_complement, 0)
        #print(f"Original matrix: {M} \n")
        #print(f"Complement matrix: {M_complement} \n")

        original_edges = extract_nonzero_pairs(M)
        complement_edges = extract_nonzero_pairs(M_complement)

        num_original_edges = len(original_edges)
        num_complement_edges = len(complement_edges)
        
        n, _ = M.shape
        print(f"Number of vertices: {n}")

        all_edges = n*(n-1)/2
        print(f"Number of all edges: {all_edges}")

        print(f"Number of edges in the max clique graph: {num_original_edges}, \nNumber of edges in the complement graph: {num_complement_edges}")
        print(num_original_edges + num_complement_edges == all_edges) # Test that the sum is exactly the num. of edges of a complete graph

        with open(result_file, 'w') as output:
            output.write(f"{n} {num_complement_edges}\n")
            for (i, j) in complement_edges:
                output.write(f"{i} {j}\n")

        print(f"\Converted the max clique grpah to: {result_file}")


# Example use
build_complement_instances('Clique', 'MIS')
