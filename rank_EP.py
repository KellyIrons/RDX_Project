#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Rank_EP
"""

def rank_EP(EP, mechobj):   
    
    #Sort EP Vector from highest to lowest value
    import numpy as np
    J = len(EP)
    sort_ind = np.argsort(EP)
    sort_ind = np.flipud(sort_ind)
    ranked_EP = EP[sort_ind]
    
    print(sort_ind)
    
    species = mechobj['species_names']
    ranked_spec = []
    for i in range(J):
        
        index = sort_ind[i]
        #print(i)
        #print(index)
        print(species[index])
        ranked_spec.append(species[index])
        '''
        index = np.where(sort_ind == i)
        ranked_spec.append(species)
        '''
        
    
    J = len(EP)
    print('Rank       Species         EP')
    print('------------------------------')
    for i in range(J):
        string = '%d       %s       %e \n' % (i+1,ranked_spec[i], ranked_EP[i])
        print(string)
     

    
    
    return [ranked_EP, ranked_spec, sort_ind]