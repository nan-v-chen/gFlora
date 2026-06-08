# gFlora
graph convolution-based functional co-response group discovery (https://ieeexplore.ieee.org/abstract/document/10964872)

## 1. Installation
You can install the latest development version of **gFlora** from GitHub:
```r
install.packages("devtools")
devtools::install_github("nan-v-chen/gFlora")
```

## 2. A quick example
The input data includes
1) the topological abundance data (**M.csv**) and,
2) the functional variable (**fv.csv**).

For demonstration purposes, the example dataset provided in this package can be loaded directly using:
```r
data <- gFlora::load_example_data()

M <- data$M
y <- data$y
```
**It is worth noting that the example data were ONLY sampled for testing purposes, NOT a whole microbial community!**

To discover the functional group, simply run:
```r
out <- gFlora::discover(M, y, n_max = 5)
```

To repeat the discovery multiple times, collect the results, and draw an aggregated graph, run:
```r
res <- gFlora::discover_pipeline(M, y, n_max = 5, n_iter = 10)
```

## 3. Functions and Parameters
### `discover(M, y, pk=NULL, fit_f="nmax", y_type="continuous", fun_s="+", n_max=NULL, alpha=40, max_iter=500, max_fitness=1, run=100, pop_size=200, parallel=TRUE, monitor=FALSE)`
Discover a functional co-response group using a genetic algorithm that maximizes the association between the aggregated topological abundance and a functional variable.
#### Parameters
| Parameter | Type | Description |
|-----------|------|-------------|
| `M` | matrix | Topological abundance matrix of shape `(n_samples, n_taxa)`. |
| `y` | vector | Target functional variable. |
| `pk` | vector | Optional prior knowledge vector. Elements equal to 1 indicate taxa that must be included in the discovered group. |
| `fit_f` | character | Fitness function. `"nmax"` constrains the maximum group size, `"regul"` applies regularization. |
| `y_type` | character | Type of functional variable. One of `"continuous"`, `"binary"`, or `"categorical"`. |
| `fun_s` | character | Direction of association. `"+"` searches for positively associated groups and `"-"` searches for negatively associated groups. |
| `n_max` | integer | Maximum number of selected taxa. |
| `alpha` | numeric | Regularization parameter used when `fit_f="regul"`. |
| `max_iter` | integer | Maximum number of genetic algorithm iterations. |
| `max_fitness` | numeric | Early stopping threshold for fitness value. |
| `run` | integer | Number of generations without improvement before stopping. |
| `pop_size` | integer | Population size of the genetic algorithm. |
| `parallel` | logical | Whether to enable parallel computation. |
| `monitor` | logical | Whether to record fitness evolution during optimization. |
#### Returns
| Return | Type | Description |
|---------|------|-------------|
| `out` | list | A list containing the discovered functional co-response group and associated statistics. |

#### Elements
| Element | Type | Description |
|---------|------|-------------|
| `fitness` | numeric | Best fitness value obtained. |
| `x` | numeric vector | Binary selection vector indicating selected taxa. |
| `members` | character vector | Names of selected taxa. |
| `abundance` | numeric vector | Aggregated topological abundance of the selected group. |
| `performance` | numeric | Pearson correlation (`continuous`) or eta coefficient (`binary`/`categorical`). |
| `fitnesses` | data.frame | Fitness values recorded during optimization (if monitored). |

---

### `discover_pipeline(M, y, n_iter=10, pk=NULL, fit_f="nmax", y_type="continuous", fun_s="+", n_max=NULL, alpha=40, max_iter=500, max_fitness=1, run=100, pop_size=200, parallel=TRUE)`
Repeat functional co-response group discovery multiple times and construct an aggregated co-selection network of the discovered taxa.
#### Parameters
| Parameter | Type | Description |
|-----------|------|-------------|
| `M` | matrix | Topological abundance matrix of shape `(n_samples, n_taxa)`. |
| `y` | vector | Target functional variable. |
| `n_iter` | integer | Number of repeated discovery iterations. |
| `pk` | vector | Optional prior knowledge vector. |
| `fit_f` | character | Fitness function. `"nmax"` or `"regul"`. |
| `y_type` | character | Type of functional variable. One of `"continuous"`, `"binary"`, or `"categorical"`. |
| `fun_s` | character | Direction of association. `"+"` or `"-"`. |
| `n_max` | integer | Maximum number of selected taxa. |
| `alpha` | numeric | Regularization parameter used when `fit_f="regul"`. |
| `max_iter` | integer | Maximum number of genetic algorithm iterations. |
| `max_fitness` | numeric | Early stopping threshold for fitness value. |
| `run` | integer | Number of generations without improvement before stopping. |
| `pop_size` | integer | Population size of the genetic algorithm. |
| `parallel` | logical | Whether to enable parallel computation. |
#### Returns
| Return | Type | Description |
|---------|------|-------------|
| `res_list` | list | Results from all discovery iterations. Each element contains the iteration index, selected taxa, performance, and the full output of `discover()`. |
#### Attributes
| Attribute | Type | Description |
|-----------|------|-------------|
| `nodes` | data.frame | Aggregated node information used for visualization. |
| `edges` | data.frame | Aggregated edge information used for visualization. |
| `graph` | visNetwork | Interactive co-selection network generated from repeated discoveries. |

## 4. More instruction
To get the topological abundance data (**M.csv**) from the original OTU file (samples * taxa), run the following Python script:
```python
import pandas as pd
import numpy as np
from elemi import EleMi, row_clr, col_normalize

otu_row_normalized = row_clr(otu.astype(float).values, pseudo_switch=False, clr_switch=False)
otu_normalized = col_normalize(otu_row_normalized)
A = EleMi(otu_normalized, 0.1, 0.01)
adj = (A + A.T) / 2 + np.eye(A.shape[0])
adj = pd.DataFrame(adj）

D = np.diag(np.sum(adj, axis=1))
D_sqrt_inv = np.linalg.inv(np.sqrt(D))
L = D_sqrt_inv.dot(adj).dot(D_sqrt_inv)
M = otu_row_normalized.dot(L)
M = pd.DataFrame(M)
```
