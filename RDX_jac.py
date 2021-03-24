'''
RDX Jacobian
'''

def rdx_jac(t, Tset):
    
    import numpy as np
    
    #Call related functions
    from RDXskeletalparser import RDXskeletalparser
    [Rlist, Slist, lg] = RDXskeletalparser()
    J = len(Slist)
    ng = len(lg)
    #Calculate the rate constant for each reaction, in forward and backward direction
    from rateConstantCalc2 import rateConstantCalc2
    [KF, KB] = rateConstantCalc2(Rlist, Slist, Tset)
    from stoichcalc import stoichcalc
    [expmatrixF, coeffmatrixF, expmatrixB, coeffmatrixB] = stoichcalc(Rlist,Slist)
    from rdx_zdot import RDX_zdot
     
    
    #Set Up
    #y = np.ones((J+ng,1)) #give each species an initial concentration of one
        #this is the correct size for y
    y = np.ones((J,1))
    dist = 1e-6; #define a disturbance in each concentration (might make an input later)
    #Jac = np.zeros((J+ng,J+ng))
    Jac = np.zeros((J,J))
    #print(Jac)

    
    #Calculate the original rate vector
    zdot_o = RDX_zdot(t,y, Rlist, Slist, lg, KF, KB, expmatrixF, coeffmatrixF, expmatrixB, coeffmatrixB, Tset)     
    
    
 
    #for j in range(J+ng):
    for j in range(J):
        #perturb the species of interest
        y[j,0] = y[j,0] + dist;
            
        #calculate the new rate vector for each species
        zdot_new =  RDX_zdot(t,y, Rlist, Slist, lg, KF, KB, expmatrixF, coeffmatrixF, expmatrixB, coeffmatrixB, Tset)
            
        #evaluate the Jacobian entry for each species with this perturbation
        #print(np.shape(Jac[:,j]))
        #Jac = (zdot_new - zdot_o)/dist
        Jac_vector = (zdot_new - zdot_o)/dist
        for i in range(J):#G
            Jac[i,j] = Jac_vector[i]
    
    
    
    return Jac


###QUESTIONS###
# - concentrations vs mass fractions
# - rates are time dependent, so I'm pretty sure the Jacobian is too? If so,
#   how does this affect the determination of what species can be approximated 
#   as QSS?
# - do I need Jacobian entries for the gas phas? I think so







