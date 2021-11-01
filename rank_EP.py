#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Rank_EP
"""

def rank_EP(EP, mechobj, reduction_type, Rlist):   
    
    #Sort EP Vector from highest to lowest value
    import numpy as np
    if reduction_type in ['species', 'S']:
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

        ranked = ranked_spec        

        print('Rank       Species         EP')
        print('------------------------------')
        for i in range(J):
            string = '%d       %s       %e \n' % (i+1,ranked_spec[i], ranked_EP[i])
            print(string)
     
    elif reduction_type in ['reactions', 'R']:
        R= len(EP)
        sort_ind = np.argsort(EP)
        sort_ind = np.flipud(sort_ind)
        ranked_EP = EP[sort_ind]

        ranked_reac = []
        for i in range(R):            
            index = sort_ind[i]
            ranked_reac.append(Rlist[index]['eq'])
            
        ranked = ranked_reac
            
        print('Rank    EP                Reaction')
        print('-------------------------------------------------')
        for i in range(R):
            string = '%d       %e      %s \n' % (i+1, ranked_EP[i], ranked_reac[i])
            print(string)

    
    
    return [ranked_EP, ranked, sort_ind]