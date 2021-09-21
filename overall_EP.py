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
    EP_all = np.zeros((len(Slist),N))
    
    for i in range(N):
        a = 1
        while TIME[a] < t_array[i] and a <= len(TIME):
            a = a+1

        index[i] = a
        time[i] = TIME[int(a)]
 

    
    for i in range(N):
        
        t = time[i]
    
        sample = solutionM[int(index[i]),:]
        samplesign = sample >= 0
        
        #if samplesign.all:
        if t > -100:
        
        #print(np.shape(sample))
            
            mechobj = my_mech_obj(sample, t, Tset, Rlist, Slist, lg, expmatrixF, coeffmatrixF, expmatrixB, coeffmatrixB)
        
            EP_i = local_error_propagation(mechobj, sample, targets, reduction_type)
           # print(np.shape(EP_i))
    
            
            EP_all[:,i] = EP_i
        
    EP = np.max(EP_all, axis = 1)
    

    return [EP, mechobj]