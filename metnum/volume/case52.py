# -*- coding: utf-8 -*-
# Metode Numerik dan Elemen Hingga
# Latihan Soal 5.2 berdasarkan Buku 
# An Introduction to Computational Fluid Dynamics The Finite Volume Method
# by H K Versteeg, W Malalasekera
# Example 5.2 pg. 147-149
# Code by UMA (github: github.com/taruma/belajar-tsa)
# Note: Run this with python command in CMD or use Spyder.

# =============================================================================
# IMPORT FUNCTION AND LIBRARY
# Note: Please download module uma_func.py and put in the same folder
# =============================================================================

import numpy as np
import matplotlib.pyplot as plt
import argparse
from numpy import exp
from met_uma import get_valdict, create_axis, create_dictionary, create_zmatrix

# =============================================================================
# DEFINE FUNCTION
# =============================================================================

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
    parser.add_argument('-n', '--nodes', metavar='nodes', type=int, default=5, help='nodes (positive integer)')
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
    Set general variable/parameter
    
    args:
        L: length (m) (default=1)
        rho: density (kg/m3) (default=1)
        gamma: gamma (kg/(m.s)) (default=0.1)
        nodes: nodes (default=5)
    
    returns:
        dictionary: with key 'L, rho, gamma, nodes, dx'
    """
    dx = L/nodes
    value = [L, rho, gamma, nodes, dx]
    dictionary = create_dictionary('L,rho,gamma,nodes,dx',value)
    return dictionary

def set_case(dictvar, u=0.1, phiA=1, phiB=0):
    """
    Set variable/parameter for specific Case
    
    args:
        dictvar: dictionary of variable/parameter
        u: velocity (m/s) (default=0.1)
        phiA: Value of Phi at A (default=1)
        phiB: Value of Phi at B (default=0)
    
    returns:
        dictionary: with key 'u, F, D, gamma, rho, dx, phiA, phiB'
    """
    L, rho, gamma, nodes, dx = get_valdict(dictvar, 'L,rho,gamma,nodes,dx')
    F, D = rho*u, gamma/dx
    val = [u, phiA, phiB, F, D, L, rho, gamma, nodes, dx]
    
    dictionary = create_dictionary('u,phiA,phiB,F,D,L,rho,gamma,nodes,dx', val)
    
    return dictionary

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
    F, D, phiA, phiB = get_valdict(case, 'F,D,phiA,phiB')
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
    F, nodes = get_valdict(case, 'F,nodes')
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

def solve_exact(func, case, n=50):
    """
    Creating axis of (x, y) for exact solution
    
    args:
        func: exact solution function
        L: length (m) (default=1)
        n: nodes (default=50)
        
    returns:
        Axis of (x, y)
    """
    L, = get_valdict(case, 'L')
    x = np.linspace(0, L, n)
    y = func(x)
    return x, y

def solve_num(case):
    """
    Creating axis of (x, y) for numerical solution
    
    args:
        case: dictionary of case
    
    return:
        Axis of (x, y)
    """
    mat_a, mat_d = create_zmatrix(case)
    calc_matrix(mat_a, mat_d, equation_nodes, case)
    mat_x = np.linalg.solve(mat_a, mat_d)
    
    axis = create_axis(case)
    phiA, phiB = get_valdict(case, 'phiA,phiB')
    ordinat = np.append(mat_x, phiB)
    ordinat = np.insert(ordinat, 0, phiA)
    
    sol_num = (axis, ordinat)
    return sol_num
    
def plot_this(ncase, numerical, exact, title=''):
    """
    Plotting the result
    
    args:
        ncase: index of case
        numerical: (x, y) of numerical solution
        exact: (x, y) of exact solution
        title: title (opt, def='')
    """
    exact_x, exact_y = exact
    num_x, num_y = numerical
    plt.figure(ncase)
    plt.title(title)
    plt.xlabel('Distance x (m)')
    plt.ylabel('$\phi$')
    plt.plot(num_x, num_y, 'bo', label='numerical')
    plt.plot(exact_x, exact_y, 'y--', label='analytical')
    plt.legend()
    plt.grid()
    
def calc_error(case, numerical, func, title='Percentage Error',
               kolom=15, prec=6):
    """
    Calculate and Print of Percentage Error
    
    args:
        case: dictionary of case
        numerical: list of (x, y) numerical solution
        func: function of exact solution
        title: string of title (optional, default="Percentage Error")
        kolom: integer of column space (optional, default=15)
        prec: integer of precision (optional, default=6)
    """
    nodes, L = get_valdict(case, 'nodes,L')
    y_num = numerical[1]
    
    node = np.arange(1, nodes+1)
    distance = create_axis(case)[1:-1]
    sol_exact = func(distance)
    y_num = y_num[1:-1]
    
    diff, percenterr = [], []
    for i, val in enumerate(sol_exact):
        dval = val - y_num[i]
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
                node[i], distance[i], y_num[i], sol_exact[i], diff[i], percenterr[i], k=kolom-2, p=prec
                ))
    print('='*lebar)
    
def print_info(case, title='Case', kolom=20):
    """
    Print Information of case
    
    args:
        case: dictionary of case
    """
    def print_desc(desc, val, unit):
        print('| {:<{k:d}s} | {:> {k:d}f} | {:^{k:d}s} |'.format(
                desc, val, unit, k=kolom-2))
    
    lebar = kolom*3+4
    print('='*lebar)
    print('|{:^{l:d}s}|'.format(title, l=lebar-2))
    print('='*lebar)
    print('|{:^{k:d}s}|{:^{k:d}s}|{:^{k:d}s}|'.format('description', 'value', 'unit', k=kolom))
    print('-'*lebar)
    
    nodes, rho, gamma, dx, L = get_valdict(case, 'nodes,rho,gamma,dx,L')
    u, F, D = get_valdict(case, 'u,F,D')
    
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

# =============================================================================
# MAIN PROGRAM HERE
# =============================================================================

if __name__ == '__main__':
    args = parse_command_line()
    var = set_var(nodes=args.nodes)
    
    check = True if not ((args.case1) or (args.case2)) else False
    
    if args.case1 or check:
        case1 = set_case(var, u=0.1)
        print_info(case1)
        exact = solve_exact(func_c1exact, case1)
        numerical = solve_num(case1)
        calc_error(case1, numerical, func_c1exact, title='Percentage Error Case 1')
        plot_this(1, numerical, exact, title='case 1 $u = 0.1$')

    print('-'*30)

    if args.case2 or check:
        case2 = set_case(var, u=2.5)
        print_info(case2)
        exact = solve_exact(func_c2exact, case2)
        numerical = solve_num(case2)
        calc_error(case2, numerical, func_c2exact, title='Percentage Error Case 1')
        plot_this(2, numerical, exact, title='case 2 $u = 2.5$')    
    
    if not args.nofigure:
        plt.show()