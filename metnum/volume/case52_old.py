# -*- coding: utf-8 -*-
# Metode Numerik dan Elemen Hingga
# Latihan Soal 5.2 berdasarkan Buku 
# An Introduction to Computational Fluid Dynamics The Finite Volume Method
# by H K Versteeg, W Malalasekera
# Example 5.2 pg. 147-149
# Code by UMA (github: github.com/taruma/belajar-tsa)

# =============================================================================
# IMPORT FUNCTION AND LIBRARY
# Note: Please download module uma_func.py and put in the same folder
# =============================================================================

import numpy as np
import matplotlib.pyplot as plt
import argparse
from numpy import exp
from uma_func import create_zmatrix, create_axis

def parse_command_line():
    """
    Parse the command line arguments and return the parse_args object.
    
    args:
        None
        
    returns:
        args: generated argparse object with all the passed command line arguments      
    """
    parser = argparse.ArgumentParser()

    # Optional Argument
    parser.add_argument('-n', '--nodes', metavar='nodes', type=int, 
                        default=5, help='nodes (positive integer)')
    parser.add_argument('-c1', '--case1', action='store_true', help='show only case 1')
    parser.add_argument('-c2', '--case2', action='store_true', help='show only case 2')
    parser.add_argument('-nf', '--nofigure', action='store_true', help='disable figure')
    
    return parser.parse_args()

def func_c1exact(x):
    """
    Return value of exact solution of case 1
    
    args:
        x: distance (m)
    
    returns:
        float: based on x
    """
    x = np.array(x)
    return (2.7183-exp(x))/1.7183

def func_c2exact(x):
    """
    Return value of exact solution of case 2
    
    args:
        x: distance (m)
    
    returns:
        float: based on x
    """
    return 1+(1-exp(25*x))/(7.2*10**10)

def set_var(L=1, rho=1, gamma=0.1, nodes=5):
    """
    Set variable/parameter
    
    args:
        L: length (m) (default=1)
        rho: density (kg/m3) (default=1)
        gamma: gamma (kg/(m.s)) (default=0.1)
        nodes: nodes (default=5)
    
    returns:
        dictionary: with key 'L, rho, gamma, nodes, dx'
    """
    dx = L/nodes
    return {'L': L, 'rho': rho, 'gamma': gamma, 'nodes': nodes, 'dx': dx}

def set_case(u=0.1, D=0.5, gamma=0.1, rho=1, dx=0.2, phiA=1, phiB=0):
    """
    Set variable/parameter for specific Case
    
    args:
        u: velocity (m/s) (default=0.1)
        L: length (m) (default=1)
        rho: density (kg/m3) (default=1)
        gamma: gamma (kg/(m.s)) (default=0.1)
        nodes: nodes (default=5)
        D: (default=0.5)
        phiA: (default=1)
        phiB: (default=0)
    
    returns:
        dictionary: with key 'u, F, D, gamma, rho, dx, phiA, phiB'
    """
    F, D = rho*u, gamma/dx
    return {'u': u, 'F': F, 'D': D, 'gamma': gamma, 
            'rho': rho, 'dx': dx, 'phiA': phiA, 'phiB': phiB}

def equation_nodes(direction, position, case):
    """
    Calculation based upwind scheme
    
    args:
        direction: direction of flow (Boolean) (default=True as positive)
        position: position of nodes ('left'/'mid'/'right')
        case: dictionary of set_case

    returns:
        list: value of aW, aP, aE, Su
    """
    keys = 'F,D,phiA,phiB'.split(',')
    F, D, phiA, phiB = map(case.get, keys)
    if direction:
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
    else:
        if position.lower() == 'mid':
            aW, aE = D, D-F
            SP = 0
            Su = 0
            aP = aW + aE + SP
            return aW, aP, aE, Su    
        elif position.lower() == 'left':
            aW, aE = 0, D-F
            SP = 2*D
            Su = phiA*(2*D)
            aP = aW + aE + SP
            return aW, aP, aE, Su
        elif position.lower() == 'right':
            aW, aE = D, 0
            SP = 2*D-F
            Su = phiB*(2*D-F)
            aP = aW + aE + SP
            return aW, aP, aE, Su
        
def calc_matrix(mat_a, mat_d, node_set, case):
    """
    Fill matrix A and D based on direction/equation of nodes
    
    args:
        mat_a: Matrix A
        mat_d: Matrix D
        node_set: function of matrix
        case: dictionary of case
        
    """
    F = case['F']
    direction = True if F >= 0 else False
    
    for i in range(0, nodes):
        for j in range(0, nodes):
            if i == j and i == 0:
                aW, aP, aE, Su = node_set(direction, 'left', case)
                mat_a[i, j] = aP
                mat_a[i, j+1] = -aE
                mat_d[i] = Su
            elif i == j and (i > 0 and i < nodes-1):
                aW, aP, aE, Su = node_set(direction, 'mid', case)
                mat_a[i, j-1] = -aW
                mat_a[i, j] = aP
                mat_a[i, j+1] = -aE
                mat_d[i] = Su
            elif i == j and (i == nodes-1):
                aW, aP, aE, Su = node_set(direction, 'right', case)
                mat_a[i, j-1] = -aW
                mat_a[i, j] = aP
                mat_d[i] = Su

def prep_y(matrix, case):
    """
    Return value of y including at point A and B
    
    args:
        matrix: Matrix Result from solver
        case: dictionary of case
    
    returns:
        list: value of y
    """
    keys = 'phiA,phiB'.split(',')
    A, B = map(case.get, keys)
    y_num = np.append(matrix, B)
    y_num = np.insert(y_num, 0, A)
    return y_num

def create_xy(func, L=1, n=50):
    """
    Return x, y of exact solution
    
    args:
        func: exact solution function
        L: length (m) (default=1)
        n: nodes (default=50)
        
    returns:
        list (x, y)
    """
    x = np.linspace(0, L, n)
    y = func(x)
    return x, y

def calc_case(dictvar, dictcase, case=1):
    """
    Calculating case and returns coordinate numerical and exact solution
    
    args:
        dictvar: dictionary of variable/parameter
        dictcase: dictionary of case
        case: number of case (default=1)
        
    returns:
        (x, y) numerical solution and (x, y) exact solution
    """
    
    
    varkeys = ('dx', 'nodes', 'L', 'rho', 'gamma')
    dx, nodes, L, rho, gamma = map(dictvar.get, varkeys)

    # Buat Matrix
    mat_a, mat_d = create_zmatrix(nodes)
    calc_matrix(mat_a, mat_d, equation_nodes, dictcase)
    mat_x = np.linalg.solve(mat_a, mat_d)    

    # Buat x, y numerical
    x_num = create_axis(nodes, L)
    y_num = prep_y(mat_x, dictcase)
    
    # Exact Solution
    if case == 1:
        func = func_c1exact
        x_exact, y_exact = create_xy(func, L=L)
    elif case == 2:
        func = func_c2exact
        x_exact, y_exact = create_xy(func, L=L)
    else:
        x_exact, y_exact = 0, 0
    
    return (x_num, y_num), (x_exact, y_exact)

def plot_this(num, x1, y1, x2, y2, title=''):
    """
    Plotting the result
    """   
    plt.figure(num)
    plt.title(title)
    plt.xlabel('Distance x (m)')
    plt.ylabel('$\phi$')
    plt.plot(x1, y1, 'bo', label='numerical')
    plt.plot(x2, y2, 'y--', label='analytical')
    plt.legend()
    plt.grid()
    #plt.show()
    
def calc_error(var, sol_num, func, title='Percentage Error',
               kolom=15, prec=6):
    """
    Calculate and print of Percentage Error
    """
    nodes, L = var['nodes'], var['L']
    
    node = np.arange(1, nodes+1)
    distance = create_axis(nodes, L)[1:-1]
    sol_exact = func(distance)
    sol_num = sol_num[1:-1]
    
    diff, percenterr = [], []
    for i, val in enumerate(sol_exact):
        dval = val - sol_num[i]
        diff.append(dval)
        percenterr.append(dval/val*100)
    
    lebar = kolom*6+7
    print('='*lebar)
    print('{:^{l:d}s}'.format(title, l=lebar))
    print('='*lebar)
    print('|{:^{k:d}s}|{:^{k:d}s}|{:^{k:d}s}|{:^{k:d}s}|{:^{k:d}s}|{:^{k:d}s}|'.format(
            'node', 'distance', 'numeric', 'exact', 'diff', 'error', k=kolom))
    print('-'*lebar)
    
    for i in range(0, nodes):
        print('| {:^{k:d}d} | {:> {k:d}.{p:d}f} | {:> {k:d}.{p:d}f} | {:> {k:d}.{p:d}f} | {:> {k:d}.{p:d}f} | {:> {k:d}.{p:d}f} |'.format(
                node[i], distance[i], sol_num[i], sol_exact[i], diff[i], percenterr[i], k=kolom-2, p=prec
                ))
    print('='*lebar)
    
def print_info(var, case, title='Case', kolom=20):
    def print_desc(desc, val, unit):
        print('| {:<{k:d}s} | {:> {k:d}f} | {:^{k:d}s} |'.format(
                desc, val, unit, k=kolom-2))
    
    lebar = kolom*3+4
    print('='*lebar)
    print('|{:^{l:d}s}|'.format(title, l=lebar-2))
    print('='*lebar)
    print('|{:^{k:d}s}|{:^{k:d}s}|{:^{k:d}s}|'.format('description', 'value', 'unit', k=kolom))
    print('-'*lebar)
    
    nodes, rho, gamma, dx, L = map(var.get, ('nodes', 'rho', 'gamma', 'dx', 'L'))
    u, F, D = map(case.get, ('u', 'F', 'D'))
    
    print_desc('\\rho', rho, 'kg/(m^3)')
    print_desc('\\gamma', gamma, 'kg/(m.s)')
    print_desc('L', L, 'm')
    print_desc('nodes', nodes, '-')
    print_desc('dx', dx, 'm')
    print('-'*lebar)
    print_desc('u', u, 'm/s')
    print_desc('F', F, '-')
    print_desc('D', D, '-')
    print('='*lebar)
    print()

# MAIN PROGRAM Here
if __name__ == '__main__':
    args = parse_command_line()
    var = set_var(nodes=args.nodes)
    keys = 'dx,rho,gamma,nodes'.split(',')
    dx, rho, gamma, nodes = map(var.get, keys)
    
    check = True if not ((args.case1) or (args.case2)) else False
    
    if args.case1 or check:
        case1 = set_case(u=0.1, dx=dx, rho=rho, gamma=gamma)
        print_info(var, case1, title='Case 1')
        (x1, y1), (x2, y2) = calc_case(var, case1, 1)
        calc_error(var, y1, func_c1exact, title='Percentage Error Case 1')
        plot_this(1, x1, y1, x2, y2, 'case 1 $u = 0.1$')
    
    print('-'*30)
    
    if args.case2 or check:
        case2 = set_case(u=2.5, dx=dx, rho=rho, gamma=gamma)
        print_info(var, case2, title='Case 2')
        (x1, y1), (x2, y2) = calc_case(var, case2, 2)
        calc_error(var, y1, func_c2exact, title='Percentage Error Case 2')
        plot_this(2, x1, y1, x2, y2, 'case 2 $u = 2.5$')
    
    if not args.nofigure:
        plt.show()