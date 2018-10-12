# -*- coding: utf-8 -*-
"""
Created on Tue Oct  9 19:40:53 2018

@author: tarum
"""

import numpy as np
import matplotlib.pyplot as plt

dik = dict(L=500,
           N=6,
           dx=100,
           K=0.1,
           B=10,
           T=1,
           q=0.0001,
           HR=5
           )

L, N, dx, K, B, T, q, HR = dik.values()

Ti = []
for i in range(N):
    if i < N // 2:
        Ti.append(2*T)
    else:
        Ti.append(1*T)
Ti = np.array(Ti)

def caly(T):
    return - q * dx**2 / T

y = []
for i in range(N-1):
    if i == 0:
        y.append(caly(Ti[i]))
    elif i == N-2:
        y.append(caly(Ti[i]) - HR)
    else:
        y.append(caly(Ti[i]))
y = np.array(y)        


def create_matrix(n,a=1,b=2,c=-1):
    mat = (np.eye(n, k=-1)*a +
           np.eye(n, k=0)*b +
           np.eye(n, k=1)*c)
    return mat

X = create_matrix(5,1,-2,1)
X[0,1] = 2
H = np.linalg.solve(X,y)
H = np.append(H, HR)

print('Solusi: {}'.format(H))


fig, ax = plt.subplots(1,1)
ax.plot(np.linspace(0,L,N),H)
#plt.show()


## Kecepatan Darcy
def calq(h2,h0,T,i=True):
    if i: return -T/B*(h2-h0)/(2*dx)
    else: return -T/B*(h2-h0)/(dx)
    
q = []
for i in range(N):
    if i == 0:
        print('A')
        q += [calq(H[i+1],H[i],Ti[i],False)]
    elif i == N-1:
        print('B')
        q += [calq(H[i],H[i-1],Ti[i],False)]
    else:
        print('C')
        q += [calq(H[i+1],H[i-1],Ti[i])]

print(q)