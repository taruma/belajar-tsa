# -*- coding: utf-8 -*-
# Metode Numerik dan Elemen Hingga
# Latihan Soal 4.2 berdasarkan Buku 
# An Introduction to Computational Fluid Dynamics The Finite Volume Method
# by H K Versteeg, W Malalasekera
# Example 4.2 pg. 121-125
# Code by UMA (github: github.com/taruma/belajar-tsa)
# Note: Run this with python command in CMD or use Spyder.

# =============================================================================
# IMPORT FUNCTION AND LIBRARY
# Note: Please download module uma_func.py and put in the same folder
# =============================================================================

import numpy as np
import matplotlib.pyplot as plt
import argparse
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
    parser.add_argument('-l', '--length', metavar='length', type=float, default=2, help='length (meter)')
    parser.add_argument('-k', '--conductivity', metavar='conductivity', type=float, default=0.5, help='constant thermal conductivity (W/m.K)')
    parser.add_argument('-q', '--heatgeneration', metavar='heatgeneration', type=float, default=1000, help='uniform heat generation (kW/m^3)')
    parser.add_argument('-TA', '--tempA', metavar='tempA', type=int, default=100, help='temperature at A (Celcius)')
    parser.add_argument('-TB', '--tempB', metavar='tempB', type=int, default=200, help='temperature at A (Celcius)')
    parser.add_argument('-n', '--nodes', metavar='nodes', type=int, default=5, help='nodes (positive integer)')
    parser.add_argument('-A', '--area', metavar='area', type=float, default=1, help='area (m^2)')
    parser.add_argument('-nf', '--nofigure', action='store_true', help='disable figure')
    parser.add_argument('-nd', '--nodetail', action='store_true', help='disable detail')
    return parser.parse_args()

def func_exact_sol(x, case):
    """
    Function of exact solution
    
    args:
        x: position (meter)

    return:
        result of f(x)
    """
    TB, TA, L, q, k = get_valdict(case, 'TB,TA,L,q,k')
    return ((TB-TA)/L + (q/(2*k)*(L-x)))*x + TA

def set_parameter(L=2, k=0.5, q=1000, TA=100, TB=200, nodes=5, A=1):
    """
    Set parameter/variable of case 4.2. Default value: 
    (L=2, k=0.5, q=1000, TA=100, TB=200, nodes=5, A=1)
    
    args:
        L: length (cm) 
        k: constant thermal conductivity (W/m.K)
        q: uniform heat generation (kW/m^3)
        TA: temperature at faces A (Celcius)
        TB: temperature at faces B (Celcius)
        nodes: nodes of control volume
        A: Area (m^2)
    
    return:
        dictionary of (L, k, q, TA, TB, nodes, A, dx)
    
    examples:
        set_parameter(nodes=10) 'Use default parameter except nodes'
    """
    # Set to meter unit
    L, q = L/100, q*1000
    dx = L/nodes
    val = [L, k, q, TA, TB, nodes, A, dx]
    
    dictparameter = create_dictionary('L,k,q,TA,TB,nodes,A,dx', val)
    return dictparameter

def equation_nodes(position, case):
    """
    Return value of aW, aE, aP, and Su based on position
    
    args:
        position: (string) 'left'/'mid'/'right'
        case: dictionary of case
        
    returns:
        list of [aW, aE, aP, Su]
    
    """
    TA, TB, k, A, dx, q = get_valdict(case, 'TA,TB,k,A,dx,q')
    if position.lower() == 'left':
        aW = 0
        aE = k*A/dx
        SP = -2*k*A/dx
        aP = aW+aE-SP
        Su = q*A*dx + 2*k*A*TA/dx
        return aW, aE, aP, Su
    elif position.lower() == 'mid':
        aW = k*A/dx
        aE = k*A/dx
        SP = 0
        aP = aW+aE-SP
        Su = q*A*dx
        return aW, aE, aP, Su        
    elif position.lower() == 'right':
        aW = k*A/dx
        aE = 0
        SP = -2*k*A/dx
        aP = aW+aE-SP
        Su = q*A*dx + 2*k*A*TB/dx
        return aW, aE, aP, Su        

def calc_matrix(mat_a, mat_d, case):
    """
    Calculate tridiagonal matrix. [A][X] = [D]
    
    args:
        mat_a: Matrix [A]
        mat_d: Matrix [D]
        case: dictionary of case
    """
    nodes, = get_valdict(case, 'nodes')
    for i in range(0, nodes):
        for j in range(0, nodes):
            if i == j and i == 0:
                aW, aE, aP, Su = equation_nodes('left', case)
                mat_a[i, j] = aP
                mat_a[i, j+1] = -aE
                mat_d[i] = Su
            elif i == j and (i > 0 and i < nodes-1):
                aW, aE, aP, Su = equation_nodes('mid', case)
                mat_a[i, j-1] = -aW
                mat_a[i, j] = aP
                mat_a[i, j+1] = -aE
                mat_d[i] = Su
            elif i == j and (i == nodes-1):
                aW, aE, aP, Su = equation_nodes('right', case)
                mat_a[i, j-1] = -aW
                mat_a[i, j] = aP
                mat_d[i] = Su

def solve_num(case):
    """
    Creating axis of (x, y) for numerical solution
    
    args:
        case: dictionary of case
        nodes: Nodes (default=100)
        
    return:
        Axis of (x, y)
    """
    mat_a, mat_d = create_zmatrix(case)
    calc_matrix(mat_a, mat_d, case)
    # create axis and ordinat numerical
    TA, TB = get_valdict(case, 'TA,TB')
    axis = create_axis(case)
    result = np.linalg.solve(mat_a, mat_d)
    ordinat = np.append(result, TB)
    ordinat = np.insert(ordinat, 0, TA)
    
    sol_num = (axis, ordinat)
    
    return sol_num

def solve_exact(case, nodes = 100):
    """
    Creating axis of (x, y) for exact solution
    
    args:
        case: dictionary of case
        nodes: Nodes (default=100)
        
    return:
        Axis of (x, y)
    """
    L, = get_valdict(case, 'L')
    axis = np.linspace(0, L, nodes)
    ordinat = func_exact_sol(axis, case)
    
    sol_exact = (axis, ordinat)
    
    return sol_exact
    
def plot_result(numerical, exact):
    """
    Plot comparison numerical vs analytical/exact solution. 
    plot.show() needed for show the figure. 
    
    args:
        numerical: list of (x, y) numerical solution
        exact: list of (x, y) exact solution
    """
    
    numx, numy = numerical
    exax, exay = exact
    plt.title('Comparison Numerical vs. Analytical Solution')
    plt.xlabel('Distance x (m)')
    plt.ylabel('Temperature $(^{o}C)$')
    plt.plot(numx, numy, 'ro', label='numerical')
    plt.plot(exax, exay, 'y--', label='analytical')
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
    sol_exact = func(distance, case)
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

def print_info(case, title='Case', kolom=25):
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
    
    L, k, q, TA, TB, nodes, dx = get_valdict(case, 'L,k,q,TA,TB,nodes,dx')
    
    print_desc('L (length)', L, 'm')
    print_desc('k (conductivity)', k, 'W/(m.K)')
    print_desc('q (heat generation)', q, 'kW/(m^3)')
    print_desc('TA (temperature at A)', TA, 'Celcius')
    print_desc('TB (temperature at B)', TB, 'Celcius')
    print_desc('nodes', nodes, '-')
    print_desc('dx (grid space)', dx, 'm')

    print('='*lebar)
    print()

# =============================================================================
# MAIN PROGRAM HERE
# =============================================================================
if __name__ == '__main__':
    args = parse_command_line()
    case42 = set_parameter(
            L=args.length,
            k=args.conductivity,
            q=args.heatgeneration,
            TA=args.tempA,
            TB=args.tempB,
            nodes=args.nodes,
            A=args.area)
    
    exact_solution = solve_exact(case42)
    numerical_solution = solve_num(case42)
    plot_result(numerical_solution, exact_solution)
    
    if not args.nodetail:
        print_info(case42, title='Example 4.2')
        calc_error(case42, numerical_solution, func_exact_sol)
    
    if not args.nofigure:
        plt.show()