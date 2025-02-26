#!/usr/bin/python3
# -*- coding: utf-8 -*-
# @Time    : 14-3-2024 15:38
# @Author  : Nan Chen
# @File    : graph_convolution.py

import pandas as pd
import numpy as np
from EleMi import EleMi, row_clr, col_normalize


otu = pd.read_csv("data/otu.csv", index_col=0)

# compute the ratio of 0 for every col (taxon)
zero_ratio = (otu == 0).astype(int).sum(axis=0) / otu.shape[0]
# filter out taxa whose ratio of 0 across samples are over 0.8
otu = otu[zero_ratio[zero_ratio < 0.8].index]

# compute the adjacency matrix using EleMi
otu_row_normalized = row_clr(otu.astype(float).values, pseudo_switch=False, clr_switch=False)
otu_normalized = col_normalize(otu_row_normalized)
A = EleMi(otu_normalized, 0.1, 0.01)
adj = (A + A.T) / 2 + np.eye(A.shape[0])
# reindex
adj = pd.DataFrame(adj, index=otu.columns, columns=otu.columns)

# compute M
D = np.diag(np.sum(adj, axis=1))
D_sqrt_inv = np.linalg.inv(np.sqrt(D))
L = D_sqrt_inv.dot(adj).dot(D_sqrt_inv)
M = otu_row_normalized.dot(L)
M = pd.DataFrame(M, index=otu.index, columns=otu.columns)
