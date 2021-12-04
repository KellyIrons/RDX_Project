"""
RDX Overall EP
"""
def overall_EP(N, TIME, solutionM, Tset, Rlist, Slist, lg, expmatrixF, coeffmatrixF, expmatrixB, coeffmatrixB, targets, reduction_type):

    import numpy as np    

    from my_mech_obj import my_mech_obj
    from RDX_EP import local_error_propagation
    
    t_array = np.linspace(0,max(TIME), N)
    index = np.zeros((N))
    time = np.zeros((N))
    if reduction_type in ['species', 'S']:
        EP_all = np.zeros((len(Slist),N))
        DIC_ind = np.zeros((len(Slist)))
        DIC = np.zeros((len(Slist),len(Slist), N))
    
    else:
        EP_all = np.zeros((len(Rlist),N))
    
    ind_list = []
    ind_list_f = []
    coeffs = []
    #coeff_list = []
    #coeffs = np.zeros((len(Slist), 2, N))
    iterations = np.zeros((N))
    for i in range(N):
        a = 1
        while TIME[a] < t_array[i] and a <= len(TIME):
            a = a+1

        index[i] = a
        time[i] = TIME[int(a)]
 
    for i in range(len(Slist)):
        ind_list_f.append([])
 
    
    for i in range(N):
        
        t = time[i]
    
        sample = solutionM[int(index[i]),0:len(Slist)]
        samplesign = sample >= 0
        
        #if samplesign.all:
        if t > -100:
        
        #print(np.shape(sample))
            
            mechobj = my_mech_obj(sample, t, Tset, Rlist, Slist, lg, expmatrixF, coeffmatrixF, expmatrixB, coeffmatrixB)
        
            [EP_i, ind_listi, DIC_i, coeffs_i, iterations_i] = local_error_propagation(mechobj, sample, targets, reduction_type)
           # print(np.shape(EP_i))
    
            
            EP_all[:,i] = EP_i
            ind_list.append(ind_listi)
            DIC[:,:,i] = DIC_i[:-1,:]
            iterations[i] = iterations_i
        
            #coeffs[:,:,i] = coeffs_i[:,1:3]       
            coeffs.append(coeffs_i[:,1:])
        
    EP = np.max(EP_all, axis = 1)

    
    for i in range(len(Slist)):
        if sum(EP[i] == EP_all[i,:]) == 1:
            ind = np.where(EP[i] == EP_all[i,:])
            ind_list_f[i] = ind_list[int(ind[0])][i]
            DIC_ind[i] = int(ind[0])

            #coeff_list.append(coeffs[int(ind[0])][i][1:])
            
        elif sum(EP[i] == EP_all[i,:]) > 1:
            print('Multiple Matches: ', Slist[i]['name'])
            ind_list_f[i] = ind_list[N-1][i]
            DIC_ind[i] = int(N-1)
            #coeff_list.append(coeffs[N-1][i][1:])
        else:
            print('Error: ', Slist[i]['name'])
    

    return [EP, mechobj, ind_list_f, DIC, DIC_ind, coeffs, iterations]




