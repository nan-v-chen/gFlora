#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@author:Nan Chen
@file: EleMi.py
@time: 2023/08/20
"""

import numpy as np
import logging


# use logging module for easy debug
logging.basicConfig(format='%(asctime)s %(levelname)8s: %(message)s', datefmt='%m-%d %H:%M:%S')
logger = logging.getLogger(__name__)
logger.setLevel('DEBUG')


def S(alpha, v):
    """
    shrinkage thresholding operator
    :param alpha:
    :param v:
    :return:
    """
    return np.sign(v) * max(0, abs(v) - alpha)


def EleMi(X, miu1, miu2):
    """
    Elastic Net Regularized Multi-Regression
    :param X:
    :param miu1:
    :param miu2:
    :return:
    """
    max_iter = 500
    p = X.shape[1]
    rho = 0.1
    rho_max = 1e10
    beta0 = 1.1
    Y = np.zeros((p, p))
    A = np.zeros((p, p))
    Z = np.zeros((p, p))       # C
    t = 0
    thre_conv = 1e-7
    thre_beta = 1e-4
    while t < max_iter:
        t = t + 1
        # logger.debug(t)
        A_old = A
        Z_old = Z
        temp1 = 2 * miu2 * np.eye(p) + X.T @ X + rho * np.eye(p)
        temp2 = X.T @ X + rho * Z + Y
        A = np.linalg.inv(temp1) @ temp2
        A = A - np.diag(np.diag(A))
        for m in range(p):
            for n in range(p):
                Z[m][n] = S(2*miu1/rho, A[m][n]-Y[m][n]/rho)
        if np.max(abs(Z-A)) < thre_conv:
            break
        else:
            Y = Y + rho * (Z - A)
            if max(np.max(np.abs(Z-Z_old)), np.max(np.abs(A-A_old))) > thre_beta:
                beta = beta0
            else:
                beta = 1
            rho = min(rho_max, beta * rho)
    return Z


def row_clr(X, pseudo=0, pseudo_switch=True, pre_norm="TSS", pre_norm_switch=True, clr_switch=True):
    """
    sample normalization
    :param X:
    :param pseudo:
    :param pseudo_switch:
    :param pre_norm:
    :param pre_norm_switch:
    :param clr_switch:
    :return:
    """
    X_normalized = X.copy()
    for row_id in range(X_normalized.shape[0]):
        if pseudo_switch:
            # try different pseudo
            if isinstance(pseudo, (int, float)) and pseudo > 0:
                zero_rep = pseudo
            else:
                # default
                zero_rep = 0.1 * np.min(X_normalized[row_id, :][X_normalized[row_id, :] != 0])
            X_normalized[row_id, :][X_normalized[row_id, :] == 0] = zero_rep
        if pre_norm_switch:
            # try different pre_norm to address different sequencing depth issue
            if pre_norm.upper() == "CSS":
                cumulative_sum = np.cumsum(X_normalized[row_id, :])
                scaling_factor = np.percentile(cumulative_sum, 0.9 * 100)
                X_normalized[row_id, :] = X_normalized[row_id, :] / scaling_factor
            else:
                # TSS (default)
                X_normalized[row_id, :] = X_normalized[row_id, :] / np.sum(X_normalized[row_id, :])
        if clr_switch:
            geometric_mean = np.exp(np.mean(np.log(X_normalized[row_id, :])))
            X_normalized[row_id, :] = np.log(X_normalized[row_id, :] / geometric_mean)
    return X_normalized


def col_normalize(X):
    """
    feature normalization
    :param X:
    :return:
    """
    X_normalized = X.copy()
    for col_id in range(X_normalized.shape[1]):
        X_normalized[:, col_id] = X_normalized[:, col_id] / np.linalg.norm(X_normalized[:, col_id])
    return X_normalized