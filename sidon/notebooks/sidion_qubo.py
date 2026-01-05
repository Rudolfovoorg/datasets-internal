import json
import networkx as nx
import time
from copy import deepcopy
from dwave.system import DWaveSampler, FixedEmbeddingComposite



class Ising(nx.Graph):
    def __init__(self, *args, info={}, **kwargs):
        super().__init__(*args, **kwargs)
        self.info = info
    
    def serialize(self):
        return {'data': nx.node_link_data(self, edges='edges'), 'info': self.info}

    def deserialize(self, serialized_data):

        self.__init__(nx.node_link_graph(serialized_data['data'], edges='edges'), info=serialized_data['info'])

        return self
    
    def copy(self):
        return self.__copy__()
    
    def __copy__(self):
        cls = self.__class__
        result = cls.__new__(cls)
        result.__dict__.update(self.__dict__)
        return result
    

def get_topology(solver="Advantage2_system1.8"):
    try:
    
        sampler = DWaveSampler(solver=solver)
        active_qubits = sampler.nodelist
        active_couplers = sampler.edgelist
        properties = sampler.properties

        return {
            "active_qubits": active_qubits,
            "active_couplers": active_couplers,
            "properties": properties,
        }

    except Exception:
        raise ImportError("D-Wave system package not available or import failed.")
    
def load_data(case, isings):
    print(f'Loading data from: {case} ...')
    with open(f'{case}/graph_all.txt') as f:
        data = f.read() 
        items = (item.strip().replace('(', '[').replace(')', ']') for item in data.split('\n'))
        key_value = (item.split('=') for item in items if item)
        d = dict(((key.strip(), json.loads(value)) for key, value in key_value))

    for ising in isings:
        with open(f'{case}/{ising}.json') as f:
            J = ((eval(key), value) for key, value in json.load(f).items())
            d[ising] = dict(J)

    print('... Done.')

    return d

def get_G_ising(data, info={}):
 #   G = nx.Graph(weight='h', default=0.0)
 #   G.add_weighted_edges_from(((u, v, w) for (u, v), w in data.items() if w != 0) , weight='J')
 #   nx.set_node_attributes(G, 0.0, name="h")


    ising = Ising(info=info, weight='h', default=0.0)
    ising.add_weighted_edges_from(((u, v, w) for (u, v), w in data.items() if w != 0) , weight='J')
    nx.set_node_attributes(ising, 0.0, name="h")    

    return ising


def solve_ising(ising, sampler, **kwargs):
    time_b = time.time()
    sampleset = sampler.sample_ising(h=nx.get_node_attributes(ising, 'h'), 
                                     J=nx.get_edge_attributes(ising, 'J'), 
                                     **kwargs)
    running_time = time.time() - time_b
    
    return {'value': sampleset.first.energy, 'X': sampleset.first.sample, 'sampler_response': sampleset, 'ising': ising, 
            'time': running_time,
            'sampler': {'class': type(sampler).__name__, 'kwargs': kwargs}}

def export_results_to_json(results, filename, **additional_info):
    results = {
        'value': float(results['value']),
        'X': {node: int(value) for node, value in results['X'].items()},
        'sampler_response': results['sampler_response'].to_serializable(),
        'ising': results['ising'].serialize(),
        'time': results['time'],
        'sampler': results['sampler'],
        'additional_info': additional_info
    }

    with open(filename, 'w') as f:
        json.dump(results, f, indent=4)