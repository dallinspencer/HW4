# -*- coding: utf-8 -*-
"""
Created on Sat Feb 13 15:05:07 2021

@author: zanka
Physics 513R
Problem 1
"""
import numpy as np
from scipy.sparse import diags

def popMatrices(n, gamma):
    grid,h = np.linspace(0, 1, n,retstep=True)
    A = diags([1/h**2-gamma/2,2/h**2,1/h**2-gamma/2], [-1,0,1], shape=(n,n))
    b = np.zeros_like(grid) + h
    b[0] = h/2
    b[-1] = h/2
    
    return A, b