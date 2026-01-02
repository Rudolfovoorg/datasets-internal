## k-cluster or densest subgraph Problem Instances

https://cedric.cnam.fr/~lamberta/Library/k-cluster.html

The instances are encoded as unconstrained binary optimization problem on a complement of a graph. Thus, the densest k-subgraph problem become an equivalent sparsest k-subgraph problem on a graph compment. The objective is modified as follows: 


$$x^TQx + offset = 1/2*x^TA(\overline G)x-\lambda(\sum_{x_i\in x} x_i-k) + \mu(\sum_{x_i\in x} x_i-k)^2$$


Where: $A$ is adjacency matrix, $offset = \lambda k + \mu k^2$.

Note, that optimal value for the densest k-subgraph problem imply optimal value for the sparsest k-subgraph problem as follows:


$$opt^{sparsest}_k = \frac{k(k-1)}{2}-opt^{densest}_k
$$.


The quadratic penalty term for problem instance $i$ is encoded ussing varying $\mu$ as follows:


$$
\mu^i_1 < \mu^i_2 < \mu^i_3
$$.

    .
    └── k-cluster
        ├── README.md           # this file            
        ├── 1 # problem instances encoded with $\mu^i_1$
        ├── 1 # problem instances encoded with $\mu^i_2$
        ├── 3 # problem instances encoded with $\mu^i_3$
        └── results 
            ├── 1 # problem instances encoded with $\mu^i_1$
            ├── 1 # problem instances encoded with $\mu^i_2$
            └── 3 # problem instances encoded with $\mu^i_3$

The matrix $Q$ is encoded in sparse COO format and can be loaded as follows:
```
import scipy as sp

Q_data = data["qubo"]["data"]
Q_row = data["qubo"]["row"]
Q_col = data["qubo"]["col"]
Q_shape = data["qubo"]["shape"]

Q = sp.sparse.coo_matrix((Q_data, (Q_row, Q_col), shape=Q_shape))
```
