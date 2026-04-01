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
out <- gFlora::discover(M, y, Nmax = 5)
```

## 3. More instruction
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
