# datasets-internal

This repository contains a curated collection of datasets and benchmark instances used in our computational experiments related to the QBIQ project (https://www.rudolfovo.eu/en/nacionalni-projekti/qbiq%3A-)
The data is organized by problem type, with each top-level folder corresponding to a specific class of instances.

Where applicable, we provide:
- A description of the problem setting,
- the structure and format of the data files,
- links to original data sources,
- citations and usage notes.

---

## Repository Structure

The repository is organized as follows:
datasets-internal/
- biqbin/k_cluster/      # BiqBin k-cluster benchmark instances
- k-cluster/    # Sparsest k-cluster datasets
- sidon/                 # Sidon set instances
- test_instances/         # Internal sparsest
- max_independent_set/    # MIS / Max Clique benchmark metadata

Each folder is documented in detail below.

---

## `biqbin/k_cluster/`

### Description

This folder contains benchmark instances for the **k-cluster (densest k-subgraph)** problem encoded as **quadratic unconstrained binary optimization (QUBO)** problem.  

The densest k-subgraph problem is formulated on the **complement graph**, making it equivalent to a **sparsest k-subgraph** problem. Each instance enforces the cardinality constraint via a quadratic penalty formulation, where the objective function is defined as follows:

$$\frac{1}{2} \mathbf{x}^\mathrm{T} A(\overline{G}) \mathbf{x} - \lambda \Big( \sum_{x_i \in \mathbf{x}} x_i - k \Big) + \mu \Big( \sum_{x_i \in \mathbf{x}} x_i - k \Big)^2$$ 

where 
$A(\overline{G})$ is the adjecency matrix of the complement of the graph and the offset is defined as $offset = \lambda k + \mu k^2$.

Original data source:  
https://cedric.cnam.fr/~lamberta/Library/k-cluster.html


### Instance difficulties

For each problem difficulty, **three difficulty levels** are provided, corresponding to different values of the quadratic penalty parameter $\mu$:

$$
\mu_i^{(1)} < \mu_i^{(2)} < \mu_i^{(3)}.
$$

These difficulty levels affect how strongly the cardinality constraint is enforced and are documented in more detail in the local `README.md` inside this folder.

### Instance format
All instances are provided as JSON files and follow the naming scheme:
`kcluster{n}\_{edge_density}\_{k}\_{instance_id}.json` where
- $n$ is the number of vertices,
- edge_density specifies the graph edge density,
- $k$ is the target cluster size,
- instance_id is the instance number.
  
Example: `kcluster40_025_10_3.json` corresonds to a graph on $n = 40$ vertices with edge density of $25 \%$, cluster size $k=10$ and is the instance number 3.

---

## `k-cluster/`

### Description
This folder contains benchmark instances for the **Sparsest (or Densest) k-Cluster** problem.

The **original data source** and the **file naming convention** are the same as in [`biqbin/k_cluster/`](./biqbin/k_cluster/).  
Unlike `biqbin/k_cluster/`, this folder contains the **original instances only** (i.e., without additional difficulty variants / penalty settings).

More details about the dataset organization and file contents are provided in the local `README.md` inside this folder.

---

## `sidon/`

### Description

This folder contains QUBO problem instances based on **Sidon sets**  
(see: https://en.wikipedia.org/wiki/Sidon_sequence).

The instances correspond to **random spin glass models** defined on an actual **D-Wave hardware topology**, with zero local fields $h$ and varying coupling strengths $J$.

- The data is stored in the `data/` directory.
- The folders  
  `data/Advantage2_system1.8/send_8` and  
  `data/original_miha_instances/send_8`  
  each contain **12 subfolders (`m1`–`m12`)**, corresponding to increasing instance difficulty.

Additional literature and background material are provided in the `literature/` folder.

Computational results and analysis are available in the `notebooks/` folder and are based on instances from  
`data/Advantage2_system1.8/send_8`.

---

## `test_instances/`

### Description
A collection of internal benchmark instances for the *Sparsest k-Subgraph* problem. The original sparsest $k$-subgraph instances are modified with the addition of linear and quadratic constraints in a way that the original optimal solution is preserved. Additionally, we include some QPLIB examples of linearly and quadratically constrained QUBOs with on graphs of size $n \leq 200$.

This collection includes:
- old and new internal .json files in the folders `json/qbo_constrained` and `json_old_format/qbo_constrained`,
- instances derived from the QPLIB dataset (examples 1976, 2047, 2055, 2060 and 2067), available at https://qplib.zib.de/instances.html.

A detailed description of how the constraints are generated is available in the local `README.md` inside this folder.

**Remarks**
- The internal json test instances are not publicly available elsewhere.
- They are intended for internal comparison and reproducibility.
- All internal json test instances were provided by dr. Roman Kužel of Rudolfovo (e-mail: roman.kuzel@rudolfovo.eu).

---

## `max_independent_set/`

**Description**  
This folder documents the *Maximum Independent Set (MIS)* / *Maximum Clique* benchmark instances used in our experiments.

**Important note**  
The original datasets are **not included** in this repository.  
They are publicly available from the Network Data Repository:

- https://networkrepository.com/networks.php

We therefore only document **which instances were used** and their **properties**.

**Instances used**

The table below lists the MIS instances used in our computations, together with their number of vertices $n$, number of edges $m$, and known optimal value $opt$.  
These instances coincide with those reported in the *Max Independent Set* section of the BiqCrunch results (https://biqcrunch.lipn.univ-paris13.fr/results): 

| Instance name   | n    | m     | opt |
|-----------------|------|-------|-----|
| hamming6-2      | 64   | 192   | 32  |
| hamming6-4      | 64   | 1312  | 4   |
| hamming8-2      | 256  | 1024  | 128 |
| johnson16-2-4   | 120  | 1680  | 8   |
| johnson8-4-4    | 70   | 560   | 14  |
| MANN_a9         | 45   | 72    | 16  |
| keller4         | 171  | 5100  | 11  |
| san200_0.7_1    | 200  | 5970  | 30  |
| san200_0.7_2    | 200  | 5970  | 18  |
| san200_0.9_1    | 200  | 1990  | 70  |
| san200_0.9_2    | 200  | 1990  | 60  |
| san200_0.9_3    | 200  | 1990  | 44  |
| brock200_1      | 1006 | 5066  | 21  |
| brock200_2      | 200  | 10024 | 12  |
| brock200_3      | 200  | 7852  | 15  |
| brock200_4      | 200  | 6811  | 17  |

**Required citation**  
If you use these instances, please cite:

> Ryan A. Rossi and Nesreen K. Ahmed.  
> *The Network Data Repository with Interactive Graph Analytics and Visualization.*  
> AAAI, 2015.

---

## Citation

If you use data from this repository in academic work, please:
1. Cite the original data sources listed for each folder.
2. Cite any accompanying publications describing the experimental setup.

---

## Notes and Future Extensions
- Each folder may include its own `README.md` with more detailed documentation.


