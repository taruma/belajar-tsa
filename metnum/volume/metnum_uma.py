import numpy as np
from decimal import Decimal as dec

def add_dec(num1, num2):
    return float(dec(str(num1)) + dec(str(num2)))

def create_axis(nodes, dx):
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
    return axisx