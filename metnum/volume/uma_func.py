# -*- coding: utf-8 -*-

"""
My Personal Function by UMA
"""
from decimal import Decimal as dec

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
    return map(dictionary.get, string.split(','))
