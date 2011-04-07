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

    no_regions = False
    if no_regions:
        # no borrowing between regions
        C = np.eye(n) * 10.e6 * sigma**2.
        return C
    
    C = np.eye(n)
    for ii in range(n-2):
        for jj in range(n-2):
            C[ii,jj] += .1
    for S in superregions:
        for ii in S:
            for jj in S:
                C[ii,jj] += 10.

    C[n-2,n-2] = 11.
    C[n-1,n-1] = 11.

    #print "making prior for sex effect very uninformative"
    #C[n-1,n-1] = 1001.
    
    C *= sigma**2.

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

altsuperregions = [
    [15, 5, 9],
    [7, 8, 1],
    [17, 18, 19, 20],
    [14],
    [3],
    [0, 4, 2, 16],
    [10, 11, 12, 13, 6],
    ]
    
gbd_regions = {u'Asia Pacific, High Income': 0,
               u'Asia, Central': 1,
               u'Asia, East': 2,
               u'Asia, South': 3,
               u'Asia, Southeast': 4,
               u'Australasia': 5,
               u'Caribbean': 6,
               u'Europe, Central': 7,
               u'Europe, Eastern': 8,
               u'Europe, Western': 9,
               u'Latin America, Andean': 10,
               u'Latin America, Central': 11,
               u'Latin America, Southern': 12,
               u'Latin America, Tropical': 13,
               u'North Africa/Middle East': 14,
               u'North America, High Income': 15,
               u'Oceania': 16,
               u'Sub-Saharan Africa, Central': 17,
               u'Sub-Saharan Africa, East': 18,
               u'Sub-Saharan Africa, Southern': 19,
               u'Sub-Saharan Africa, West': 20}
