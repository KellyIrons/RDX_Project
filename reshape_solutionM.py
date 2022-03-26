#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 27 13:49:42 2022

@author: kirons
"""

def reshape_solutionM(solutionM):
    
    import numpy as np
    
    a = len(solutionM)
    b = int(np.shape(solutionM[0])[0])
    
    new_solutionM = np.zeros((a,b))
    
    for i in range(a):
        new_solutionM[i,:] = solutionM[i]
        
    return new_solutionM