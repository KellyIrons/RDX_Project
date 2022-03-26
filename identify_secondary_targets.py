"""
Further Reduction of RDX Mechanism - Identifying Secondary Targets
"""

def identify_secondary_targets(TIME, solutionM, targets, t, Tset, Rlist, Slist, lg, expmatrixF, coeffmatrixF, expmatrixB, coeffmatrixB, EP):
    
    
    '''
    Identify "locally important" species to include in targets for further reduction
    
    '''
    import numpy as np
    from my_mech_obj import my_mech_obj
    from RDX_compute_DIC import compute_DIC
    
    
    J = len(Slist)
    EP_max = 1E-2
    DIC_min = 1E-10
    N_min = 5
    g_mean = np.zeros((J))
    Ns = np.zeros((J))
    
    ind  = 0
    found = False
    while ind <= len(TIME) and not found:
        found = t <= TIME[ind]
        ind = ind + 1
    
    
    sample = solutionM[ind][:J]
    
    
    mechobj = my_mech_obj(sample, t, Tset, Rlist, Slist, lg, expmatrixF, coeffmatrixF, expmatrixB, coeffmatrixB)
    
    
    [DIC_spec, DIC_reac] = compute_DIC(mechobj, sample, 'species')
    
    
    # Check each species
    for i in range(J):
        # make sure it's not already in the targets
        if Slist[i]['name'] not in targets:
        
            # check that the EP isn't too high
            if EP[i] <= EP_max:
                
                # Calulcate geometric mean of the DIC values for this species
                specs = []
                prod = 1
                N = 0
                for j in range(J):
                    if DIC_spec[j,i] > DIC_min:   
                        if i != j:
                            N = N + 1
                            specs.append(DIC_spec[j,i])

                for species in specs:
                    prod = prod*(species**(1/N))
                #g_mean[i] = N*prod**(1/N)
                if N >= N_min:
                    g_mean[i] = N*prod
                    Ns[i] = N
            
                
    sort_ind = np.argsort(g_mean)[::-1]
    
    print('Rank     pecies         Weight       Connections')
    print('---------------------------------------------------')
                
    for i in range(len(sort_ind)):
        if g_mean[sort_ind[i]] > 0:
            string = '%d       %s       %e      %d \n' % (i+1, Slist[int(sort_ind[i])]['name'], g_mean[sort_ind[i]], Ns[sort_ind[i]] )
            print(string)
                
    
    return [g_mean, sort_ind]
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
    

    

