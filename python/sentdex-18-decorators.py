# -*- coding: utf-8 -*-
"""
Created on Fri Sep 28 09:51:08 2018

@author: tarum
"""
from functools import wraps


def add_wraping(item):
    @wraps(item)
    def wrapped_item():
        return 'a wrapped up box of {}'.format(str(item()))
    return wrapped_item

@add_wraping
def new_gpu():
    return 'a new Tesla P100 GPU'

print(new_gpu())
print(new_gpu.__name__)