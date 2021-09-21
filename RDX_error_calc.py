"""
RDX Error Calculation
"""

def RDX_error_calc(TIME_o, solutionM_o, solutionM_r, TIME_r, N):
    
    
    #Measure the median percent difference in results at set time steps
        #Note: this draws heavily from reduction_analysis.py
    
    import numpy as np
    
    #Get the max time
    T = max(TIME_o)

    #Determine times at which to analyze results
    ts = np.zeros((1,N))
    frac = 1/N
    for i in range(N):
        ts[0,i] = frac*(i+1)
        
        #KELLY - see what this looks like, can I use linspace?
    
    old_i = np.zeros((1,N))
    new_i = np.zeros((1,N))
    
    #Find the indices of these points
    for t in range(N-1):      
        i = 0
        while ts[0,t] > TIME_o[i] and i <= len(TIME_o):
            i = i+1
        old_i[0,t] = i
        i = 0
        while ts[0,t] > TIME_r[i] and i <= len(TIME_r):
            i = i+1
        new_i[0,t] = i
    old_i[0, N-1] = len(TIME_o) - 1 
    new_i[0, N-1] = len(TIME_r) - 1
    
        
          
    #Calculate the percent difference
    
    [J,b] = np.shape(solutionM_o)
    results = np.zeros((1,N))

    #For each time step - can I vectorize this?
    for i in range(N):
        old_ind = int(old_i[0,i])
        new_ind = int(new_i[0,i])
        p_diff = np.zeros((1,J))
        #For each species
        for j in range(J):    
            #find the percent different between new and old
            p_diff = abs(100*(solutionM_r[new_ind,:] - solutionM_o[old_ind,:])*(1/solutionM_o[old_ind,:]))
            #save this in a matrix
            #p_diff[0,j] = percent_diff
        results[0,i] = np.median(p_diff)
        #take the average 
    
    error = max(results[0,:])    
        
    import matplotlib.pylab as plt
   
    '''
    #plt.plot(TIME_o, solutionM_o, TIME_r, solutionM_r)
    plt.plot(TIME_r, solutionM_r)
    plt.show()
    '''
    
    return(error)