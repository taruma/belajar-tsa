# -*- coding: utf-8 -*-
"""
Created on Fri Oct  5 18:32:45 2018

@author: tarum
"""

test = dict(a=3,
            b=4,
            c=5,
            d=7)

deret = [1,2,3,4,5]

def printdis(e=2, f=10, *args, **kwargs):
    print(args)
    print(kwargs)
    print(e, f)
    
printdis(2, **test)