'''
RDX QSS Species Finder
'''
def rdx_qss_finder(r,N,Tset, S):
    #r = ode result from RDX_Attempt script
    #N = desired number of Jacobian evaluations
    #Tset = set temperature
    #S = number of species
    
    #Set up
    import numpy as np
    from RDX_jac import rdx_jac
    
    #Pull apart r
    t = r.t; y = r.y
    tmin = min(t); tmax = max(t)
    t_array = np.linspace(tmin,tmax,N)
    index = np.zeros((1,N))
    index[0,0] = 1
    time = np.zeros((1,N))
    diag_vecs = np.zeros((S,N))
    
    
    for i in range(N-1):
        a = 1
        while t[a] < t_array[i+1] and a <= len(t):
            a = a+1
        index[0,i+1] = int(a)
        time[0,i+1] = t[int(a)]
        
     
    for i in range(N):
        
        #Set Up
        ti = t_array[i] #Current time
        ind_i = int(index[0,i])
        y0 = y[:,ind_i] #mass fractions at this time
        
        #Calculate the Jacobian at this time
        Jac = rdx_jac(ti, Tset, y0)
        
        #Take the diagonals
        [a,b] = np.shape(Jac) 
       # diag_vec = np.zeros((a,1))
        for j in range(a):
           # diag_vec[j,0] = np.abs((Jac[j,j])**(-1))
            diag_vecs[j,i] = np.abs((Jac[j,j])**(-1))
        
        
        
        print(i)

    return(index, time, diag_vecs)
    
    
    
    
'''

Jac = rdx_jac(0.5, 538)

det_Jac = np.linalg.det(Jac)

inv_Jac = np.linalg.inv(Jac)
'''