""" This module collects some alternatives for modeling inter-regional variation

Each function creates a variance-covariance matrix suitable for using
in the negative binomial model, which is a positive (n,n)-matrix, with
n-2 regions plus a dimension for sex and another for time """

import numpy as np

def all_related_equally(n, sigma):
    C = np.eye(n)
    for ii in range(n-2):
        for jj in range(n-2):
            C[ii,jj] += 1.
    C *= sigma**2.
    return C

def regions_nested_in_superregions(n, sigma):
    C = np.eye(n)
    for ii in range(n-2):
        for jj in range(n-2):
            C[ii,jj] += .1
    for S in superregions:
        for ii in S:
            for jj in S:
                C[ii,jj] += 10.
    C *= sigma**2.

    C[n-2,n-2] = 5.
    C[n-1,n-1] = 5.
    
    return C

# indices correspond to order of regions in dismod3.settings.gbd_regions
superregions = [
    [15, 5, 9, 0, 12],
    [7, 8, 1],
    [17, 18, 19, 20],
    [14],
    [3],
    [4, 2, 16],
    [10, 11, 13, 6],
    ]
