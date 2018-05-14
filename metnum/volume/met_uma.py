# -*- coding: utf-8 -*-
# Personal Function developed by UMA for Numerical Method Class SI-5231
# Code by UMA (github: github.com/taruma/belajar-tsa)

from decimal import Decimal as dec
import numpy as np

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

def get_valdict(dictionary, string):
    """
    get value of dictionary based list of string
    
    args:
        dictionary: dictionary
        string: string with comma seperated no spaces
        
    return:
        list of value of string-key
    
    example:
        let DICT = {'a':3, 'ba':10, 'c':30}
        get_valdict(DICT, 'a,ba')
            return [3, 10]
        get_valdict(DICT, 'a')
            return [3]
        get_valdict(DICT, 'c,ba')
            return [30, 10]
    """
    return map(dictionary.get, string.split(','))

def create_axis(case):
    """
    Create axis of x from 0 based on nodes and dx
    
    args:
        case: dictionary of case
    
    return:
        list of axis of x
    """
    nodes, L = get_valdict(case, 'nodes,L')
    dx, axisx, x = L/nodes, [0], 0
    for i in range(1, nodes+1):
        if i == 1:
            x = add_dec(x, dx/2)
            axisx.append(x)
        else:
            x = add_dec(x, dx)
            axisx.append(x)
    axisx.append(add_dec(x, dx/2))
    return np.array(axisx)

def create_dictionary(string, value):
    """
    Creating dictionary based on string and value order
    
    args:
        string: string with comma seperated no spaces
        value: list of value to be assigned
        
    return:
        dictionary with {string-keys: value}
        
    example:
        create_dictionary('a,b,c,d', [1,10,20,300])
            return {'a':1, 'b':10, 'c':20, 'd':300}
        create_dictionary('a', [1])
            return {'a':1}
    """    

    dictionary = {}
    keys = string.split(',')
    for index, value in enumerate(value):
        dictionary[keys[index]] = value
    
    return dictionary

def create_zmatrix(case):
    """
    Create two zeros-array based on nodes, one dimension and two dimensions
    
    args:
        case: dictionary of case
    
    return:
        list of two numpy.array (2D, 1D)
    """
    nodes, = get_valdict(case, 'nodes')
    return np.zeros([nodes, nodes]), np.zeros([nodes])