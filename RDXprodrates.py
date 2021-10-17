"""
RDX Production Rates
"""

def RDX_prodrates(Tset, Rlist, Slist, KF, KB, expmatrixF, coeffmatrixF, expmatrixB, coeffmatrixB, TIME_o, solutionM_o):
    
    
    import numpy as np
    
    Ti = 337.6 #initial temp - K
    rho = 1820 #liquid density of RDX - kg/m^3
    
    S = len(Slist)
    R = len(Rlist)
    
    ProdRates = np.zeros((S,R))
    Ts = len(TIME_o)  
    
    
    
    # for each time step
    for t in range(Ts-1):
        t_step = TIME_o[t+1]- TIME_o[t]
        
        Y = solutionM_o[t,0:len(Slist)]
        
        T = Tset - (Tset-Ti)*np.math.exp(-TIME_o[t]/0.045) #liquid CV temp - K
            
        #Calculate the rate constant for each reaction, in forward and backward direction
        [a,b] = np.shape(KF) #Take size of the K matrices
        if T < Tset: #If the temperature is less than the set temp ...
            count = 0 #initialize the counter
            Trow = -1000000000 #initialize the temperature of the row
            while T > Trow and count < a: #while the current temp is greater than the row temp and the row is valid
                Trow = KF[count, 0] #set Trow
                count = count + 1 #add one to the count
                #at stop, Trow(count) is the closest temp greater than T
            if count == 1:
                Kf1 = KF[count-1, 1:b]; Kf2 = KF[count-1, 1:b]
                Kb1 = KB[count-1, 1:b]; Kb2 = KB[count-1, 1:b]  
                T1 = 0; T2 = KF[count-1,0];
            else:
                Kf1 = KF[count-2, 1:b]; Kb1 = KB[count-1, 1:b]
                Kf2 = KF[count-1, 1:b]; Kb2 = KB[count-1, 1:b]
                T1 = KF[count-2,0]; T2 = KF[count-1,0]
            #Interpolate to get the most accurate K values
            Kf = ((Kf2-Kf1)/(T2-T1))*(T-T1) + Kf1
            Kb = ((Kb2-Kb1)/(T2-T1))*(T-T1) + Kb1
        else: #Once T = Tset, the rate constants do not change
            Kf = KF[a-1,1:b] #set the rate constants to be the last values (T = Tset)
            Kb = KB[a-1,1:b]
        
        for y in range(R):               
                
            # From RDXsim...
            pitermf = np.ones([1,R])
            pitermb = np.ones([1,R]) #initialize the forward and backward pi terms
            for a in range(len(Rlist[y]['reactants'])):#for each reactant in the forward direction
                i = expmatrixF[Rlist[y]['reacI'][a],y]
                j = float(Y[Rlist[y]['reacI'][a]])
                pitermf[0,y] = pitermf[0,y]*(j**i)
            for a in range(len(Rlist[y]['products'])): #for each reactant in the backward direction (aka the products)
                i = expmatrixB[Rlist[y]['prodI'][a],y]
                j = float(Y[Rlist[y]['prodI'][a]])
                pitermb[0,y] = pitermb[0,y]*(j**i)
                
        #Calculate the net production rate of each species
        Wf = Kf*coeffmatrixF*pitermf
        Wb = Kb*coeffmatrixB*pitermb
        W = Wf + Wb
       # print(Wf)

        # for each species
        for x in range(S):
            for y in range(R):
                MW = Slist[x]['MW']
               # ProdRates[x,y] = ProdRates[x,y] + W[x,y]*MW*(1/rho)
                ProdRates[x,y] = ProdRates[x,y] + W[x,y]*(1/MW)*rho*t_step
        
    return ProdRates
                
                
            
            
            
            
            
            