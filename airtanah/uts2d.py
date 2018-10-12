# -*- coding: utf-8 -*-
"""
Created on Wed Oct 10 07:42:16 2018

@author: tarum
"""
import numpy as np
import matplotlib.pyplot as plt
from math import log10


dik = dict(NX=101,
           NY=101,
           LX=1000,
           LY=1000,
           B=30,
           H0=30,
           K=13.33,
           T=400,
           Q0=2000
           )

NX, NY, LX, LY = [dik.get(x) for x in 'NX NY LX LY'.split()]

dx, dy = LX/(NX-1), LY/(NY-1)

Q0_, T_, K_, B_, H0_ = [dik.get(x) for x in
                  'Q0 T K B H0'.split()]

T_ = K_*B_

max_iter = 1e5
con_crit = 1e-5

# initialization
H = np.ones((NX, NY)) * H0_
Q = np.zeros_like(H)
X = np.linspace(0,LX,NX)
Y = np.linspace(0,LY,NY)
i1, j1 = 25,25
i2, j2 = 75,75
Q[i1, j1] = Q0_
Q[i1, j2] = Q0_
Q[i2, j1] = Q0_
Q[i2, j2] = Q0_
nH = np.ones((NX,NY)) * H0_


# Iteration
iteration = 0
err = 0
result = []

for _ in range(10000):
    iteration += 1
    err = 0
    nH = np.copy(H)
    for i in range(NX):
        for j in range(NY):
            if ((i==0) or (i==NX-1)) or ((j == 0) or (j == NY-1)):
                pass
            else:
                H[i,j] = (-Q[i, j]/T_ + nH[i+1,j] + nH[i-1,j] +
                                      + nH[i,j+1] + nH[i,j-1])/4
        pass

    err = abs(H-nH).sum()

#    

    result.append([iteration,log10(err),err])
    
    
#    if iteration % 10 == 0:
#        print(f'{iteration:10d}, {err:10.7f}')
        

#    
    
print(np.array(result))
