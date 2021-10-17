#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 16 14:24:14 2021

@author: kirons
"""
def write_pathways(prodrates, Slist, Rlist):
    import numpy as np
    f = open('myfile.txt', 'w')
    
    for i in range(len(Slist)):
        f.write('Species: %s \n' % Slist[i]['name'])
        f.write('#:       Weight:        Equation: \n')
        
        rates = prodrates[i,:]
        
        #maxrate = max(abs(rates))
        
        # find the nonzero indices
        posind = rates > 0
        negind = rates < 0
        
        pos = np.sort(rates[posind])
        pos = pos[::-1]
        neg = np.sort(rates[negind])   
        
        
        if len(pos) > 0:
            totratepos = sum(pos)
        else:
            totratepos = 1
            
        if len(neg) > 0:
            totrateneg = sum(abs(neg))
        else:
            totrateneg = 1
        
        
        for a in range(len(pos)):
            reac = np.where(rates == pos[a])
            f.write('%2d       %1.4e     %s\n' %(reac[0]+1, pos[a]/totratepos, Rlist[int(reac[0])]['eq']))
        for a in range(len(neg)):
            reac = np.where(rates == neg[a])
            f.write('%2d      %1.4e     %s\n' %(reac[0]+1, neg[a]/totrateneg, Rlist[int(reac[0])]['eq']))
        f.write('\n')