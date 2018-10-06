# -*- coding: utf-8 -*-
"""
Created on Fri Oct  5 09:01:45 2018

@author: tarum
"""

# Ignore this block code.
# Define Function for printing result
def cetak(name, val, unit = "", pad=50):
    print('| {:>15s} = {:>10.5f} {:<{n}s}|'.format(name, val, unit, n = pad-15-10-7))
    
def new_line(pad = 50):
    print('='*pad)


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
    string = string.replace(' ', '')
    return map(dictionary.get, string.split(','))