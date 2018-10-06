# -*- coding: utf-8 -*-
"""
Created on Fri Oct  5 18:11:34 2018

@author: tarum
"""

import sys, os
sys.path.append(os.path.abspath('..'))

import umakit
from umakit.general import get_valdict

print(umakit.basicmath.luas_lingkaran(4))
tes = dict(a=3,
           b=5,
           c=8)

a = get_valdict(tes,'a,c')
print(list(a))