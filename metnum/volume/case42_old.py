# -*- coding: utf-8 -*-
"""
Metode Numerik
Penyelesaian Example 4.2
"""

import numpy as np
import matplotlib.pyplot as plt
from decimal import Decimal as dec
import argparse

def parse_command_line():
    """
    Parse the command line arguments and return the parse_args object.
    
    There are 6 optional arguments.
    The help message generated by the parser should look like:
    
    usage: case42.py    [-L length] [-k conductivity] 
                        [-q heat_generation] 
                        [-TA temperature_A] [-TB temperature_B] 
                        [-n nodes] [-A Area]

    optional arguments:
      -h, --help            show this help message and exit
      -l, --length          length (meter)
      -k, --conductivity    constant thermal conductivity (W/m.K)
      -q, --heatgeneration  uniform heat generation (kW/m^3)
      -TA, --tempA          temperature at A (Celcius)
      -TB, --tempB          temperature at B (Celcius)
      -n, --nodes           nodes (positive integer)
      -A, --area            area (m^2)

    args:
        None
        
    returns:
        args: generated argparse object with all the passed command line arguments      
    """
    parser = argparse.ArgumentParser()

    # Optional Argument
    parser.add_argument('-l', '--length', metavar='length', type=float,
                        default=2, help='length (meter)')
    parser.add_argument('-k', '--conductivity', metavar='conductivity', type=float,
                        default=0.5, help='constant thermal conductivity (W/m.K)')
    parser.add_argument('-q', '--heatgeneration', metavar='heatgeneration', type=float,
                        default=1000, help='uniform heat generation (kW/m^3)')
    parser.add_argument('-TA', '--tempA', metavar='tempA', type=int,
                        default=100, help='temperature at A (Celcius)')
    parser.add_argument('-TB', '--tempB', metavar='tempB', type=int,
                        default=200, help='temperature at A (Celcius)')
    parser.add_argument('-n', '--nodes', metavar='nodes', type=int, 
                        default=5, help='nodes (positive integer)')
    parser.add_argument('-A', '--area', metavar='area', type=float,
                        default=1, help='area (m^2)')
    
    return parser.parse_args()

def add_dec(num1, num2):
    """
    Return addition of two decimal numbers
    
    args:
        num1: first number
        num2: second number
    
    returns:
        addition of two number
        
    examples:
        add_dec(0.1, 0.1/2) ==> 0.15 instead 0.150000...
    """
    return float(dec(str(num1)) + dec(str(num2)))

def create_axis(nodes, L):
    """
    Create axis of x from 0 based on nodes and dx
    
    args:
        nodes: nodes of control volume
        L: length (m)
    
    return:
        list of axis of x
    """
    dx = L/nodes
    axisx = [0]
    x = 0
    for i in range(1, nodes+1):
        if i == 1:
            x = add_dec(x, dx/2)
            axisx.append(x)
        else:
            x = add_dec(x, dx)
            axisx.append(x)
    axisx.append(add_dec(x, dx/2))
    return np.array(axisx)

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
        list of (L, k, q, TA, TB, nodes, A, dx)
    
    examples:
        set_parameter(nodes=10) 'Use default parameter except nodes'
    """
    # Set to meter unit
    L = L/100
    q = q*1000
    
    return L, k, q, TA, TB, nodes, A

def create_zmatrix(nodes):
    """
    Create two zeros-array based on nodes, one dimension and two dimensions
    
    args:
        nodes: nodes
    
    return:
        list of two numpy.array (1D, 2D)
    """
    return np.zeros([nodes, nodes]), np.zeros([nodes])

def middle_nodes():
    """
    Calculate value in middle nodes
    
    return:
        list of aW, aE, aP, Su
    """
    
    aW = k*A/dx
    aE = k*A/dx
    SP = 0
    aP = aW+aE-SP
    Su = q*A*dx
    return aW, aE, aP, Su
    
def left_nodes():
    """
    Calculate value in left nodes
    
    return:
        list of aW, aE, aP, Su
    """
    
    aW = 0
    aE = k*A/dx
    SP = -2*k*A/dx
    aP = aW+aE-SP
    Su = q*A*dx + 2*k*A*TA/dx
    return aW, aE, aP, Su    

def right_nodes():
    """
    Calculate value in right nodes
    
    return:
        list of aW, aE, aP, Su
    """
    
    aW = k*A/dx
    aE = 0
    SP = -2*k*A/dx
    aP = aW+aE-SP
    Su = q*A*dx + 2*k*A*TB/dx
    return aW, aE, aP, Su        

def calc_matrix(mat_a, mat_d):
    """
    Calculate matrix
    """
    for i in range(0, nodes):
        for j in range(0, nodes):
            if i == j and i == 0:
                aW, aE, aP, Su = left_nodes()
                mat_a[i, j] = aP
                mat_a[i, j+1] = -aE
                mat_d[i] = Su
            elif i == j and (i > 0 and i < nodes-1):
                aW, aE, aP, Su = middle_nodes()
                mat_a[i, j-1] = -aW
                mat_a[i, j] = aP
                mat_a[i, j+1] = -aE
                mat_d[i] = Su
            elif i == j and (i == nodes-1):
                aW, aE, aP, Su = right_nodes()
                mat_a[i, j-1] = -aW
                mat_a[i, j] = aP
                mat_d[i] = Su

def solve_matrix(mat_a, mat_d):
    """
    Solving matrix use numpy.linalg.solve
    [A][R] = [D]
    
    args:
        mat_a: matrix [A]
        mat_d: matrix [D]
    
    return:
        matrix [R]
    """
    result = np.linalg.solve(mat_a, mat_d)
    return result

def prep_y(result):
    """
    Creating axis of y for plotting

    args:
        result: result of solve_matrix
    
    return:
        array of y value with TA and TB
    """
    y_num = np.append(result, TB)
    y_num = np.insert(y_num, 0, TA)
    return y_num

def func_exact_sol(x):
    """
    Function of exact solution
    
    args:
        x: position (meter)

    return:
        result of f(x)
    """
    return ((TB-TA)/L + (q/(2*k)*(L-x)))*x + TA

def exact_sol(nodes = 100):
    """
    Creating axis of (x, y) for exact solution
    
    args:
        nodes: Nodes (default=100)
        
    return:
        Axis of (x, y) and nodes
    """
    exact = []
    axisx_exact = np.linspace(0, L, nodes)
    for x in axisx_exact:
        hasil = func_exact_sol(x)
        exact.append(hasil)
    return axisx_exact, exact, nodes
    
def plot_result():
    plt.title('Comparison Numerical vs. Analytical Solution')
    plt.xlabel('Distance x (m)')
    plt.ylabel('Temperature $(^{o}C)$')
    plt.plot(axisx, y_num, 'ro', label='numerical')
    plt.plot(axisx_exact, exact, 'y--', label='analytical')
    plt.legend()
    
    node_ann = 3
    plt.annotate('Numerical', xy=(axisx[node_ann], y_num[node_ann]), xytext=(
            axisx[node_ann], y_num[node_ann]-20),
            arrowprops=dict(arrowstyle='->', connectionstyle='arc3'),)
    
    mid_exact = int(nodes_exact/2)-10
    plt.annotate('Exact Solution', xy=(axisx_exact[mid_exact], exact[mid_exact]),
                 xytext=(axisx_exact[mid_exact], exact[mid_exact]-20),
                 arrowprops=dict(arrowstyle='->', connectionstyle='arc3'),)
    plt.grid()
    plt.show()    

if __name__ == '__main__':

    args = parse_command_line()
    L, k, q, TA, TB, nodes, A = set_parameter(
            L=args.length,
            k=args.conductivity,
            q=args.heatgeneration,
            TA=args.tempA,
            TB=args.tempB,
            nodes=args.nodes,
            A=args.area)

    
    dx = L/nodes
    
    print('Suhu di titik \tA = {:>5d} C\nSuhu di titik \tB = {:>5d} C\ndengan panjang \tL = {:>1.3f} m'.format(
            TA, TB, L))
    print('Jumlah Noda \t  = {:>5d}\tdengan dx = {:1.3f} m'.format(nodes, dx))
    
    axisx = create_axis(nodes, L)
    mat_a, mat_d = create_zmatrix(nodes)
    calc_matrix(mat_a, mat_d)
    result = solve_matrix(mat_a, mat_d)
    
    print('Matrix A = \n{}\n====='.format(mat_a))
    print('Matrix D = \n{}\n====='.format(mat_d))
    print(
    'Penyelesaian Matrix [AX = D], diperoleh matrix X = \n{}\n====='.format(result))

    y_num = prep_y(result)
    axisx_exact, exact, nodes_exact = exact_sol()
    
    result_exact = []
    for x in axisx:
        result_exact.append(func_exact_sol(x))
    
    for counter, val in enumerate(result_exact):
        hasil = (y_num[counter] - val)/val*100
        print('Percentage Error pada titik {:3d} = {:2.2f}'.format(counter, hasil))

    plot_result() 