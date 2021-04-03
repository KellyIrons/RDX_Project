"""
Experimenttal Reduction Analysis
"""

def composition_Analysis(old_t, old_M, new_t, new_M, num):
    
    #Measure the average percent difference in results at set time steps
    
    import numpy as np
    
    #Get the max time
    T = max(old_t)

    #Determine times at which to analyze results
    ts = np.zeros((1,num))
    frac = 1/num
    for i in range(num):
        ts[0,i] = frac*(i+1)*T
    
    old_i = np.zeros((1,num))
    new_i = np.zeros((1,num))
    
    #Find the indices of these points
    for t in range(num-1):      
        i = 0
        while ts[0,t] > old_t[i] and i <= len(old_t):
            i = i+1
        old_i[0,t] = i
        i = 0
        while ts[0,t] > new_t[i] and i <= len(new_t):
            i = i+1
        new_i[0,t] = i
    old_i[0, num-1] = len(old_t) - 1 
    new_i[0, num-1] = len(new_t) - 1
    
        
          
    #Calculate the percent difference
    
    [J,b] = np.shape(old_M)
    results = np.zeros((2,num))
    results[0,:] = ts
    #For each time step
    for i in range(num):
        old_ind = int(old_i[0,i])
        new_ind = int(new_i[0,i])
        p_diff = np.zeros((1,J))
        #For each species
        for j in range(J):    
            #find the percent different between new and old
            percent_diff = 100*(new_M[j,new_ind] - old_M[j,old_ind])*(1/old_M[j,old_ind])
            #save this in a matrix
            p_diff[0,j] = percent_diff
        results[1,i] = np.mean(p_diff)
        #take the average 
        
    return results
    