# gFlora
graph convolution-based functional co-response group discovery (https://ieeexplore.ieee.org/abstract/document/10964872)

## 1. Installation
You can install the latest development version of **gFlora** from GitHub:
```r
install.packages("devtools")
devtools::install_github("n-v-chen/gFlora")
```

## 2. A quick example
The input data includes
1) the topological abundance data (**M.csv**) and,
2) the functional variable (**fv.csv**).

For demonstration purposes, the example dataset provided in this package can be loaded directly using:
```r
data <- load_example_data()

M <- data$M
y <- data$y
```
**It is worth noting that the example data were ONLY sampled for testing purposes, NOT a whole microbial community!**

To discover the functional group, simply run:
```r
out <- discover(M, y, Nmax = 5)
```
