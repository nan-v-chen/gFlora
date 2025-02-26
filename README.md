# gFlora
graph convolution-based functional co-response group discovery (https://arxiv.org/abs/2407.03897)

## 1. Installation
The code of gFlora mainly includes two parts: a) code in Python to infer the co-occurrence networks; b) code in R to run the genetic algorithm.

For both of them, you can download the source code and use them directly as functions following these steps:

### a) Code in Python
1. **Install required dependencies**

   To install the necessary dependencies, run the following command:
   ```sh
   pip install numpy pandas
   ```

   Or
   ```sh
   conda install numpy pandas
   ```
   
3. **Download the source code from GitHub**

   You can clone the whole repository with [Git](https://git-scm.com/):
   ```sh
   git clone https://github.com/nan-v-chen/gFlora.git
   ```
   
   Or manually download the zip and unzip it.
   
4. **Import the function in Python**

   ```python
   from EleMi import EleMi, row_clr, col_normalize
   ```

### b) Code in R
1. **Install required dependencies**

   Before using gFlora, make sure you have the necessary R packages installed. You can install them with:  
   ```r
   install.packages(c("GA", "doParallel", "data.table", "reshape2"))
   ```
   
2. **Download the source code from GitHub**

   Same as in a).2
   
3. **Load the function in R**

   ```r
   source("gFlora.R")
   ```

Now, you can simply use them as functions in Python and R.

## 2. A quick example
The input data can be found in **example data**, which includes 1) the abundance data (**otu.csv**) and 2) the functional variable (**fv.csv**). **It is worth noting that the example data were ONLY sampled for testing purposes, NOT a whole microbial community!**
1. **Co-occurrence network construction**

   First, load the abundance data from CSV file:
   ```python
   import pandas as pd
   import numpy as np
   from EleMi import EleMi, row_clr, col_normalize

   otu = pd.read_csv("example data/otu.csv", index_col=0)
   ```
   otu is shaped like 51*50: 51 samples as rows and 50 taxa as columns.

   Then, infer the co-occurrence network using EleMi to get the adjacency matrix:
   ```python
   otu_row_normalized = row_clr(otu.astype(float).values, pseudo_switch=False, clr_switch=False)
   otu_normalized = col_normalize(otu_row_normalized)
   A = EleMi(otu_normalized, 0.1, 0.01)
   adj = (A + A.T) / 2 + np.eye(A.shape[0])
   # reindex
   adj = pd.DataFrame(adj, index=otu.columns, columns=otu.columns)
   adj is a weighted adjacency matrix shaped like 50*50.
   ```
   
2. **Graph convolution**

    Do the graph convolution to get the topological abundance matrix:
   ```python
   D = np.diag(np.sum(adj, axis=1))
   D_sqrt_inv = np.linalg.inv(np.sqrt(D))
   L = D_sqrt_inv.dot(adj).dot(D_sqrt_inv)
   M = otu_row_normalized.dot(L)
   M = pd.DataFrame(M, index=otu.index, columns=otu.columns)

   # save M for later usage
   M.to_csv("example data/M.csv")
   ```
   M is the topological abundance matrix shaped the same as otu.
   
3. **Choose a proper group size with AIC**

   ```r
   source("gFlora.R")
   
   y <- read.csv("example data/fv.csv", row.names=1)[, 1]
   
   ks <- seq(from=2, to=30, by=1)
   for (x in 1:10) {
     print(sprintf("Iteration: %d", x))
     set.seed(x)
     M <- read.csv("example data/M.csv", row.names = 1)
     
     # compute AIC
     aic <- sapply(ks, function(k){
       print(k)
       out <- gFlora(M, y, k=k)
       print(out$performance)
       assemblage <- out$abundance
       aic_k <- AIC(lm(y~assemblage))+2*(k-1)
       print(aic_k)
       return(aic_k)
     })
     write.csv(aic, sprintf("example data/aic_%d.csv", x), row.names=ks)
   }
   
   aic_all <- read.csv("example data/aic_1.csv", row.names=1)
   for (x in 2:10){
     aic <- read.csv(sprintf("example data/aic_%d.csv", x), row.names=1)
     aic_all <- cbind(aic_all, aic)
   }
   write.csv(aic_all, "example data/aic_all.csv", row.names=TRUE)
   ```
   An example of the AIC curve can be seen shown below (not related to the example dataset):
   
   <img src="https://github.com/user-attachments/assets/efe913c2-c3d1-48d3-81d5-67ded6dadfed" alt="aic" height="280" width="400" align="center" />

   According to the elbow criterion, the optimal group size should be 7.

4. 







