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
2. **Download the source code from GitHub**

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

Now, you can simple use them as functions in Python and R.

## 2. A quick example








