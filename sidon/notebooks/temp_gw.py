import glob
import os
import json
import pandas as pd
import glob
import os
import json
import pandas as pd
import glob
import os
from neal import SimulatedAnnealingSampler
import json
import networkx as nx
import numpy as np
import time
import dimod
import glob
import os
import json


from sidion_qubo import load_data, get_G_ising, solve_ising, export_results_to_json, Ising

import cvxpy as cvx
from qbo_sdp import qbo_sdp_01, goemans_williamson_loop



def solve_GW(Q_offset):
    Q, offset = Q_offset
    Q /= 28
    
    (value, (x, X, Y, status, runtime)) = qbo_sdp_01(
        Q, solver=cvx.MOSEK, 
        linear_constraint=None, 
        quadratic_constraints=None, verbose=True)
       
    gw_x, gw_info = goemans_williamson_loop(
        Y, Q,
        linear_constraint=None,
        quadratic_constraints=None,
        tol_lin=1e-6, tol_quad=1e-6,
        max_trials=10, seed=42,
        local_improve=None,               # or "hill-climber", function for local improvement, TBD
        stop_at_first_feasible=False)

    return {'value': 28*gw_info["ub_value"]+offset, 'sdp': 28*value+offset}


def ising2qubo(ising):
    Q, nodes = nx.attr_matrix(ising, edge_attr='J')
    Q = np.triu(Q, 1)
    node_h = nx.get_node_attributes(ising, 'h')
    L = np.array([node_h[node] for node in nodes])
    
    QUBO = 4*Q
    
    offset = -L.sum() + Q.sum()
    np.fill_diagonal(QUBO, 2*L - 2*Q.sum(axis=0) - 2*Q.sum(axis=1))

    return QUBO, offset

def solve_qubo(Q_offset, sampler, **kwargs):
    Q, offset = Q_offset
    
    time_b = time.time()
    bqm = dimod.BQM.from_qubo(Q, offset=offset)
    sampleset = sampler.sample(bqm, **kwargs)
    running_time = time.time() - time_b
    
    return {'value': sampleset.first.energy, 'X': sampleset.first.sample, 'sampler_response': sampleset, 'ising': ising, 
            'time': running_time,
            'sampler': {'class': type(sampler).__name__, 'kwargs': kwargs}}


solver_type = 'Advantage2_system1.8'
case_path = '../data/{}/send_{}/m{}/case{}/beta{}/'
list_of_cases = sorted(glob.glob(case_path.format(solver_type, '8', '*', '*', 0.0)))
list_of_cases = sorted(glob.glob(case_path.format(solver_type, '8', '6', '0', 0.0)))
#isings = ['J_log', 'J_qac']
isings = ['J_log']
num_reads = [10]


for case in list_of_cases:
    for ising_type in isings:
        with open(f'{case}{ising_type}.json') as f:
            data = json.load(f)

        ising = Ising().deserialize(data)

        for n in num_reads:
            print(f'Solving {case} with GW ...')

            
            result = solve_GW(ising2qubo(ising))
            result['case'] = case
            print(' ... Done.')

            
            output_file_folder = case.replace(f'data/{solver_type}', f'results/{solver_type}/GW')
            os.makedirs(os.path.dirname(output_file_folder), exist_ok=True)
            output_file_name = f"{output_file_folder}GW_{ising_type}.json"

            print(output_file_name)
            
            with open(output_file_name, 'w') as f:
                json.dump(result, f)
            #print(f'Exporting results to {output_file_name} ...')
            #export_results_to_json(result, output_file_name)
            print(' ... Done.')

