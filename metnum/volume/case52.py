# -*- coding: utf-8 -*-

import numpy as np
from numpy import exp
import matplotlib.pyplot as plt
from case42 import create_zmatrix, create_axis

def func_c1exact(x):
    x = np.array(x)
    return (2.7183-exp(x))/1.7183

def func_c2exact(x):
    return (1+(1-exp(25*x))/(7.2*10**10))

def set_var(L=1, rho=1, gamma=0.1, nodes=5):
    dx = L/nodes
    return {'L': L, 'rho': rho, 'gamma': gamma, 'nodes': nodes, 'dx': dx}

def set_case(u=0.1, D=0.5, gamma=0.1, rho=1, dx=0.2, phiA=1, phiB=0):
    F, D = rho*u, gamma/dx
    return {'u': u, 'F': F, 'D': D, 'gamma': gamma, 
            'rho': rho, 'dx': dx, 'phiA': phiA, 'phiB': phiB}

def positive_nodes(position, F, D, phiA, phiB):
    if position.lower() == 'mid':
        aW, aE = F + D, D
        SP = 0
        Su = 0
        aP = aW + aE + SP
        return aW, aP, aE, Su    
    elif position.lower() == 'left':
        aW, aE = 0, D
        SP = F + 2*D
        Su = (F+2*D)*phiA
        aP = aW + aE + SP
        return aW, aP, aE, Su
    elif position.lower() == 'right':
        aW, aE = F+D, 0
        SP = 2*D
        Su = 2*D*phiB
        aP = aW + aE + SP
        return aW, aP, aE, Su

def calc_matrix(mat_a, mat_d, node_set, val):
    """
    Calculate matrix
    """
    F, D, phiA, phiB = val[0], val[1], val[2], val[3]
    
    for i in range(0, nodes):
        for j in range(0, nodes):
            if i == j and i == 0:
                aW, aP, aE, Su = node_set('left', F, D, phiA, phiB)
                mat_a[i, j] = aP
                mat_a[i, j+1] = -aE
                mat_d[i] = Su
            elif i == j and (i > 0 and i < nodes-1):
                aW, aP, aE, Su = node_set('mid', F, D, phiA, phiB)
                mat_a[i, j-1] = -aW
                mat_a[i, j] = aP
                mat_a[i, j+1] = -aE
                mat_d[i] = Su
            elif i == j and (i == nodes-1):
                aW, aP, aE, Su = node_set('right', F, D, phiA, phiB)
                mat_a[i, j-1] = -aW
                mat_a[i, j] = aP
                mat_d[i] = Su

def prep_y(matrix, A, B):
    y_num = np.append(matrix, B)
    y_num = np.insert(y_num, 0, A)
    return y_num

# PROGRAM Here
if __name__ == '__main__':
    
    # Set Parameter (noda, L, rho, gamma)
    var = set_var(nodes=10)
    dx, nodes, L, rho = var['dx'], var['nodes'], var['L'], var['rho']
    gamma = var['gamma']
    
    # CASE 1 (u=0.1)
    dc1 = set_case(u=0.1, dx=dx, rho=rho, gamma=gamma)
    c1phiA, c1phiB = dc1['phiA'], dc1['phiB']
    c1u, c1F, c1D = dc1['u'], dc1['F'], dc1['D']

    
    # Solusi Case 1
    c1_axisx = create_axis(nodes, L)
    # Buat Matrix
    matc1a, matc1d = create_zmatrix(nodes)
    calc_matrix(matc1a, matc1d, positive_nodes, [c1F, c1D, c1phiA, c1phiB])
    matc1x = np.linalg.solve(matc1a, matc1d)
    c1_y = prep_y(matc1x, c1phiA, c1phiB)
    
    # Exact Solution
    c1_exactx = np.linspace(0, L, 50)
    c1_exacty = func_c1exact(c1_exactx)
    
    # plotting
    plt.figure(1)
    plt.title('Case 1: $u = 0.1$')
    plt.xlabel('Distance x (m)')
    plt.ylabel('$\phi$')
    plt.plot(c1_axisx, c1_y, 'bo', label='numerical')
    plt.plot(c1_exactx, c1_exacty, 'y--', label='analytical')
    plt.legend()
    plt.grid()
    
    # CASE 2 (u=2.5)
    dc2 = set_case(u=2.5, dx=dx, rho=rho, gamma=gamma)
    c2phiA, c2phiB = dc2['phiA'], dc2['phiB']
    c2u, c2F, c2D = dc2['u'], dc2['F'], dc2['D']

    # Solusi Case 2
    c2_axisx = create_axis(nodes, L)
    # Buat Matrix
    matc2a, matc2d = create_zmatrix(nodes)
    calc_matrix(matc2a, matc2d, positive_nodes, [c2F, c2D, c2phiA, c2phiB])
    matc2x = np.linalg.solve(matc2a, matc2d)
    c2_y = prep_y(matc2x, c2phiA, c2phiB)
    
    # Exact Solution
    c2_exactx = np.linspace(0, L, 50)
    c2_exacty = func_c2exact(c2_exactx)
    
    # plotting
    plt.figure(2)
    plt.title('Case 2: $u = 2.5$')
    plt.xlabel('Distance x (m)')
    plt.ylabel('$\phi$')
    plt.plot(c2_axisx, c2_y, 'bo', label='numerical')
    plt.plot(c2_exactx, c2_exacty, 'y--', label='analytical')
    plt.legend()
    plt.grid()