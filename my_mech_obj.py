# -*- coding: utf-8 -*-
"""
Mechanism Object Creator
"""

def my_mech_obj(sample, t, Tset):
    
    import numpy as np
    
    from RDXskeletalparser import RDXskeletalparser
    [Rlist, Slist, lg] = RDXskeletalparser()
    
    
    from stoichcalc import stoichcalc
    [expmatrixF, coeffmatrixF, expmatrixB, coeffmatrixB] = stoichcalc(Rlist,Slist)
    
    
    J = len(Slist)
    R = len(Rlist)
    
    nup = np.zeros((R,J))
    nur = np.zeros((R,J))
    
   # pos_I = coeffmatrixF > 0
   # neg_I = coeffmatrixF < 0
    coeffmatrix = np.transpose(coeffmatrixF)
    
    for i in range(R):
        for j in range(J):
            val = coeffmatrix[i,j]
            if val > 0:
                nup[i,j] = val;
            elif val < 0:
                nur[i,j] = -1*val;
        
            #nup[i,j] = coeffmatrixF[pos_I[i,j]]
            #nur[i,j] = coeffmatrixF[neg_I[i,j]]
            #nup[i,:] = coeffmatrixF[pos_I[i,:]]
   
            
    #Calculate omega

    massfc = sample
    mt = 0.0005;
    Vc = 2.7473e-07
    MW = np.zeros([1,J])
    for x in range(J):
        MW[0,x] = Slist[x]['MW']
    m = massfc*mt #m is a vector of the mass of each condensed species
    m = np.transpose(m) #transpose m
    X = (m/(MW))*(1/Vc)
    Ti = 337.6 #initial temp - K
    
    #Calculate the rate constant for each reaction, in forward and backward direction
    T = Tset - (Tset-Ti)*np.math.exp(-t/0.045) #liquid CV temp - K
    from rateConstantCalc2 import rateConstantCalc2
    [KF, KB] = rateConstantCalc2(Rlist, Slist, Tset)
    [a,b] = np.shape(KF)
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
    
    pitermf = np.ones([1,R])
    pitermb = np.ones([1,R]) #initialize the forward and backward pi terms
    for x in range(R): #for each reaction 
        for a in range(len(Rlist[x]['reactants'])):#for each reactant in the forward direction
            i = expmatrixF[Rlist[x]['reacI'][a],x]
            j = float(X[0,Rlist[x]['reacI'][a]])
            pitermf[0,x] = pitermf[0,x]*(j**i)
        for a in range(len(Rlist[x]['products'])): #for each reactant in the backward direction (aka the products)
            i = expmatrixB[Rlist[x]['prodI'][a],x]
            j = float(X[0,Rlist[x]['prodI'][a]])
            pitermb[0,x] = pitermb[0,x]*(j**i)
    #Calculate the net production rate of each species
    Wf = Kf*coeffmatrixF*pitermf
    Wb = Kb*coeffmatrixB*pitermb #DOUBLE CHECK THIS TOO
    omegaf = Kf*pitermf
    omegab = Kb*pitermb
        #Wf and Wb are matrices of the K*coeff*piterm term for each species and reaction
    #sum each row of Wf and Wb to get wf and wb for each species
    wf = Wf.sum(axis=1)
    wb = Wb.sum(axis=1)
    
    w = wf + wb #sum the forward and backward rates to get the net    
    omega = omegaf + omegab
    
    
    #Calculate net production and consumption rates
    '''
    rho = 1820 #liquid density of RDX - kg/m^3
    posw = np.zeros((1,J))
    negw = np.zeros((1,J))
    W = Wf + Wb #JxR
    for i in range(J):
        for j in range(R):
            if W[i,j] > 0:
                posw[0,i] = posw[0,i] + W[i,j]
                #negw[0,i] = abs(wb[i])
                
            else:
                negw[0,i] = negw[0,i] + abs(W[i,j])
                #posw[0,i] = wb[i]
    
    PA = posw*MW*(1/rho)
    CA = negw*MW*(1/rho)
    '''
    
    rho = 1820 #liquid density of RDX - kg/m^3
    pos = np.zeros((1,J))
    neg = np.zeros((1,J))
    for i in range(J): #for each species
        for j in range(R): #for each reaction
            if nup[j,i]*omega[0,j] > 0:
                pos[0,i] += omega[0,j]*coeffmatrixF[i,j]
            else:
                neg[0,i] += omega[0,j]*coeffmatrixF[i,j]
    PA = pos*MW*(1/rho)
    CA = neg*MW*(1/rho)
                
            

    #Find the reactions containing each species
    S_ind = []
    for s in range(J):
        sindbool = [];
        sind = [];
        other_spec = []
        for r in range(R):
            r_inds = Rlist[r]['reacI'] + Rlist[r]['prodI'] #combine product and species indices
            sindbool.append(s in r_inds)
            if s in r_inds:
                sind.append(r)
                other_spec += r_inds
                '''
                if s in Rlist[r]['reacI']:
                    other_spec += (Rlist[r]['prodI'])
                elif s in Rlist[r]['prodI']:
                    other_spec += (Rlist[r]['reacI'])
                else:
                    print('error')
                '''
        Other_spec = []
        for i in other_spec:
            if i not in Other_spec and i!=s:
                Other_spec.append(i)
            
        s_ind = {'indbool':sindbool, 'ind':sind, 'other_spec': Other_spec}
        S_ind.append(s_ind)
        
    mechobj = {'ns':J, 'nr':R, 'nup': nup, 'nur':nur, 'omega':omega, 'PA':PA, 'CA':CA, 'S_ind':S_ind}


    return mechobj
