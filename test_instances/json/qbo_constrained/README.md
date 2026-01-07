### A collection of internal benchmark instances for the sparsest k-subgraph problem

These problem instances were generated from the collection of benchmark instances for the k-cluster (densest) subgraph problem. See https://cedric.cnam.fr/~lamberta/Library/k-cluster.html

The k-cluster (densest) subgraph problem and the sparsest k-subgraph problem are trivially equivalent via a complement of a graph. I.e. an optimal k-cluster (densest k-subgraph) in some fraph $G$ is an optimal sparsest k-subgraph in the complement of $G: \overline G$. In this collection constraints are added in a way that the original optimal solution is preserved.

In the text bellow we use the following notation: 
- $A^D$ - an adjacency matrix for the k-cluster problem
- $A^S$ - an adjacency matrix for its equivallent sparsest k-subgraph problem

### Files


    .
    └── README.md           # this file      
    ├── linear
        ├── 1               # see linear/1 
        ├── 2               # see linear/2             
        ├── 3               # see linear/3             
        └── 4               # see linear/4 
    └── quadratic
        ├── 1               # see quadratic/1
        ├── 2               # see quadratic/2             
        ├── 3               # see quadratic/3             
        └── 4               # see quadratic/4 
### Linear
- linear/1: [($B_1$, $c_1$, $=$)] 
    - $Q:=\frac{1}{2}A^S$
    - $B_1:= \mathbf{1}^{1,n} $
    - $r_1:=(k,)$

- linear/2: [($B_1$, $c_1$, $=$)] 
    - $Q$ form `biqbin/k-cluster/1`
    - $B_1:= Rand(Z^0)^{m<n,n} $
    - $r_1\in R^{m, 1}$, s.t. $B_1 x^* =r_1$

- linear/3: [($B_1$, $c_1$, $=$)] 
    - $Q$ form `biqbin/k-cluster/3`
    - $B_1:= Rand(Z^0)^{m<n,n} $
    - $r_1\in R^{m, 1}$, s.t. $B_1 x^* =r_1$

- linear/4: [($B_1$, $c_1$, $\ge$)] 
    - $Q$ form `biqbin/k-cluster/3`
    - $B_1:= Rand(Z^0)^{m<n,n} $
    - $r_1\in R^{m, 1}$, s.t. $B_1 x^* \ge r_1$

### Quadratic
- quadratic/1: [($Q_1$, $r_1$, $\ge$)] 
    - $Q$ form `biqbin/k-cluster/1`
    - $Q_1:=\frac{1}{2}A^D$
    - $r_1:=\overline{opt}-\epsilon, \epsilon\in Rand(Z^+)$

- quadratic/2: [($Q_1$, $r_1$, $\le$)] 
    - $Q$ form `biqbin/k-cluster/1`
    - $Q_1:=\frac{1}{2}A^D$
    - $r_1:=\overline{opt}+\epsilon, \epsilon\in Rand(Z^+)$

- quadratic/1: [($Q_1$, $r_1$, $\ge$)] 
    - $Q$ form `biqbin/k-cluster/3`
    - $Q_1:=\frac{1}{2}A^S$
    - $r_1:=opt-\epsilon, \epsilon\in Rand(Z^+)$

- quadratic/1: [($Q_1$, $r_1$, $\le$)] 
    - $Q$ form `biqbin/k-cluster/3`
    - $Q_1:=\frac{1}{2}A^S$
    - $r_1:=opt+\epsilon, \epsilon\in Rand(Z^+)$

### Data
The matrices $Q, B_i, Q_i$ is encoded in sparse COO format and can be loaded as follows:
```
import json
import scipy as sp


((Q_data, (Q_row, Q_col)), Q_shape) = data['QBO']['Q']
Q = sp.sparse.coo_matrix((Q_data, (Q_row, Q_col)), Q_shape)

linear = [
    ( 
        sp.sparse.coo_matrix((Bi_data, (Bi_row, Bi_col)), shape=Bi_shape), 
        Ci, 
        sense
    ) for ((Bi_data, (Bi_row, Bi_col)), Bi_shape), Ci, sense in  data['QBO']['constraints']['linear']]

quadratic = [
    ( 
        sp.sparse.coo_matrix((Qi_data, (Qi_row, Qi_col)), shape=Qi_shape), 
        ri, 
        sense
    ) for ((Qi_data, (Qi_row, Qi_col)), Qi_shape), ri, sense in  data['QBO']['constraints']['quadratic']]

x = np.array(data['x'])
instance = data['instance']
optimum = data['optimum']
k = data['info']['k']
l = data['info']['lambda']
mu = data['info']['mu']
offset = data['QBO']['offset']


```
