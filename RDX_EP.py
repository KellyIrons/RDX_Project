#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 28 09:28:45 2021

@author: kirons
"""
def local_error_propagation(mechobj, mysample, targets, reduction_type):
#def local_error_propagation(mechobj, mysample, targets, reduction_type, alpha_norm_loc):
    """ Computes error propagation coefficients on one sample

    :param mechobj: Mechanism object
    :param mysample: sample on which error propagation is applied
    :param targets: list of reduction targets
    :param reduction_type: type of reduction performed (species or reactions)
    :param alpha_norm_loc: scaling coefficient for the given sample
    :return: list of error coefficients

    Created: 17/11/14 [PP]
    Last modified: 19/01/28 [QC]
    """

    #To call: EP = local_error_propagation(mymechobj, mysample, ['RDX'], 'species')
        #The brackets on RDX are important

    import numpy as np
    
    from RDX_compute_DIC import compute_DIC
    [DIC_spec, DIC_reac] = compute_DIC(mechobj, mysample, reduction_type)
    
    " TEST TO MAKE SURE ZEROS ARE SYMMETRIC"
    A = np.ceil(DIC_spec)
    B = np.transpose(A[0:53,:])
    C = A[0:53,:] - B[:,0:53]
    D = C == 0
    if not D.all():
        print('error!')
        
    E = DIC_spec <= 1
    if not E.all():
        print('error 2!')
    
    

    #myns = mechobj.ns
    #mynr = mechobj.nr
    #net = mechobj.network
    #ctmech = mechobj.ctmech
    #species_names = ctmech.species_names
    myns = mechobj['ns']
    mynr = mechobj['nr']
    species_names = mechobj['species_names']

    # Parameters
    #EPmin = 1e-7
    EPmin = 1E-15
    EPcomp = np.zeros(myns, 'd')
    EP_spec = np.zeros(myns, 'd')
    EP_reac = np.zeros(mynr, 'd')
    EP_spec_dict = {}
    
    #Added code - Find indices of target species
    target_indices = []
    for target in targets:
        if target not in ['HeatRelease', 'HeatRelease']:
            index = 0;
            while index <= myns and target != species_names[index]:
                index = index + 1
            target_indices.append(index)
            #print(target_indices)
    # Indices of targets
  #  target_indices = [species_names.index(target) for target in targets 
  #                          if target not in ['HeatRelease', 'HeatRelease']]
    #target_indices = [species_names.index(target) for target in targets 
    #                        if target not in ['HeatRelease', 'HeatRelease']]
    if 'HeatRelease' in targets or 'HeatRelease' in targets:
        target_indices.append(myns)

    DIC_spec, DIC_reac = compute_DIC(mechobj, mysample, reduction_type)


    # ------------
    # Error Propagation (EP)
    # ------------

    # Go through each target
    S_ind = mechobj['S_ind']
    for index_target_local, index_target_global in enumerate(target_indices):
        
        
        # Initialize working array - Works for HR as well
        #EPtmp = DIC_spec[index_target_global, :] * alpha_norm_loc[index_target_local]
        EPtmp = DIC_spec[index_target_global, :]

        # Initialize EP arrays
        array_up = np.zeros(myns, 'd')
        array_down = np.zeros(myns, 'd')

        # Initialize array_up
        array_up[:] = -1.0
        for i in range(myns):
            if i == index_target_global or EPtmp[i] < EPmin:
                continue
                

            array_up[i] = EPtmp[i]

        # Iterate until all relevant coefficients have been included

        flag = True
        while flag:
        

            # Init inner loop
            flag = False
            array_down[:] = -1.0

            # Loop over array_up
            for i in range(myns):
                #indj = net.indsj[net.indsi == i]
                indj = S_ind[i]['other_spec'] 
                # If coeff is positive, store coeff and continue

                if array_up[i] > 0.0:
                    coeff_up = array_up[i]

                    # Loop over all species
                    for j in indj:
                        coeff_down = DIC_spec[i, j] * coeff_up

                        # New coeff counts if i!=j and > EPmin
                        if i != j and coeff_down > EPmin:
                            flag = True
                            # Update EPtmp and array_down for next iteration
                            if coeff_down > EPtmp[j]:
                                EPtmp[j] = coeff_down
                                array_down[j] = coeff_down

            if list(EPcomp) == list(EPtmp):
                flag = False
            else:
                EPcomp = EPtmp

            #if targets[index_target_local] is not 'HeatRelease':
            if targets[index_target_local] != 'HeatRelease':
                EP_spec_dict[species_names[index_target_global]] = EPtmp
            else:
                EP_spec_dict['HeatRelease'] = EPtmp

            array_up[:] = array_down[:]
        
        # Adjust maximum coefficient for target/species pair
        EP_spec = np.fmax(EP_spec, EPtmp)

    if reduction_type in ['reactions', 'R']:
        EPtmp = np.zeros(mynr, 'd')

        for index_target_local, index_target_global in enumerate(target_indices):
            # Initialize working array - Works for HR as well
            R_target = np.zeros(mynr, 'd')
            for index_spec in range(myns):
                #if targets[index_target_local] is 'HeatRelease':
                if targets[index_target_local] == 'HeatRelease':
                    Rtmp = EP_spec_dict['HeatRelease'][index_spec] * DIC_reac[index_spec, :]
                else:
                    Rtmp = EP_spec_dict[species_names[index_target_global]][index_spec] \
                           * DIC_reac[index_spec, :]

                R_target = np.fmax(R_target, Rtmp)

            EPtmp[:] = R_target[:]

            # Adjust maximum coefficient for target/species pair
            EP_reac = np.fmax(EP_reac, EPtmp)

    if reduction_type in ['species', 'S']:
        EP = EP_spec
    else:
        EP = EP_reac

    #print(np.shape(EP_spec))
    return EP