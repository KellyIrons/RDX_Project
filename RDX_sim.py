"""
RDX Simulation
"""

def RDX_sim(Rlist, Slist, lg, Tset, KF, KB, expmatrixF, coeffmatrixF, expmatrixB, coeffmatrixB):
    import numpy as np
    import matplotlib.pylab as mat
    from scipy.integrate import solve_ivp 

    Tg = 473  #Temperature of gas region - K
    P = 101325 #Atmospheric pressure - Pa
    Vg = 5.3996E-4 #m^3 %note: this is a guess - may be more like 0.0028 if dimensions match reference [38] or 5.4E-4 [32]
    Ru = 8.314 #Gas constant - J/molK
    mN2dot = 1.4583e-06 #mass flow of purge gas (guess) - kg/s
    Ti = 337.6 #initial temp - K
    #Tset = 538 #set temp - K
    rho = 1820 #liquid density of RDX - kg/m^3
    samplemass = 0.0005 #this is probably correct
    Vc = samplemass/rho
    constantMoles = (P*Vg)/(Ru*Tg) #constant number of moles in gas phase
    
    ng = len(lg)
    A = np.zeros((1,ng))
    E =  np.zeros((1,ng))
    n = np.zeros((1,ng))
    gasI = []
    
    '''
    for x in range(ng):
        A[0,x] = lg[x]["A"]
        n[0,x] = lg[x]["n"]
        E[0,x] = lg[x]["E"]
    '''
    for x in range(ng):    
       gasI.append(lg[x]['index'])
        #NOTE: need to make some sort of correction at some point for gas 
        #   phase species that are removed
        
    R = len(Rlist)
    J = len(Slist)
    
    MW = np.zeros([1,J])
    #code to make MWs, A, E, n vectors
    for x in range(J):
        MW[0,x] = Slist[x]['MW']
        
        
    [a,b] = np.shape(KF)
    
    
    #Define Initial quantities
    YinitialRDX = 1 #initial mass fraction of RDX in liquid CV is one
    initialmass = samplemass
    initialVector = np.zeros((1,J+ng))
    initialVector = initialVector.tolist()
    initialVector = initialVector[0]
    initialVector[0] = YinitialRDX ####REDUCED
    initialVector[J] = 1
    y0, t0 = initialVector, 0


    def myrhs(t, y,):# arg1):

        
        #Fix some naming inconsistencies 
        minitial = initialmass;
        rstruct = Rlist; sstruct = Slist
    
        
        #Unpack y
        massfc = y[0:J] #mass fraction of each species in condensed phase
        molefg = y[(J):(J+ng)] #mole fraction of each species in gas phase
    
        
        #"Correct" mass and mole fractions
        sumMassFrac = np.sum(massfc)  #sum mass fractions
        massfc = massfc/sumMassFrac #correct mass fractions so they sum to one
        sumMoleFrac = np.sum(molefg) 
        if sumMoleFrac != 0: #sum mole fractions #G
            molefg = molefg/sumMoleFrac #correct mole fractions so they sum to one #G
    
        
        #Calculate the mass that has left the condensed phase  #G
        gasMoles = constantMoles*molefg #gasMoles is a vector of the number of moles of each gas phase species
        gasMoles = 	np.transpose(gasMoles) #transpose gasMoles
        gasMass = gasMoles*MW[0,gasI] #gasMass is a vector of the mass of each gas phase species
        massLoss = np.sum(gasMass, axis=1) - gasMass[0,ng-1] #total mass in gas phase (excluding purge gas)
            # this is ok for now, since the N2 doesn't evaporate atm
    
    
        #Calculate the concentrations of the condensed phase species
        massLoss = 0  ####### debugging
        mt = minitial - massLoss #total mass in condensed phase
        
        
        m = massfc*mt #m is a vector of the mass of each condensed species
        m = np.transpose(m) #transpose m
        X = (m/(MW))*(1/Vc) #X is a vector of the concentration of each condensed species
        #this section looks good
    
        ### Start by looking at the condensed phase ###
        
        #Determine the reaction temperature using eq (1)
        T = Tset - (Tset-Ti)*np.math.exp(-t/0.045) #liquid CV temp - K
        
       
    
        
        #Calculate the rate constant for each reaction, in forward and backward direction
        [a,b] = np.shape(KF) #Take size of the K matrices
        #print(a,b) #the matrices are one unit longer but I think that's ok
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
            
    
     
    
        #Calculate the liquid-gas evaporation rate constants
        #e = np.math.exp(1)
        klg = np.zeros((1,ng))
        for x in range(ng):
            klg[0,x] = lg[x]["A"]*(T**lg[x]["n"])*(np.e**(-lg[x]["E"]/(Ru*T)));

    
    
        kLG = np.zeros([1,J]) #initialize a vector of klg's for all species #GAS
        for x in range(ng):
            kLG[0,gasI[x]] = klg[0,x]
        #kLG[gasI] = klg #fill kLG for species that evaporate
            #reminder: gasI = [299, 293, 298, 20, 16, 18, 17, 19, 1, 320, 319, 291]
            #  or      gasI = [6,17,4,37,22,30,25,45,8,3,1,24] (REDUCED)
        
    
        # A quick equation explanation ...
        #
        #     The equation used to calculate the net production rate in both
        #     the forward and backward direction is:
        #         
        #         _N_            ___M___
        #         \               |  |
        #  wj =   /   vji' * ki * |  | [Xj]^(vji")
        #        / __             |  |
        #        i=1              j=1
        #     where N is the total number of reactions, M is the total number
        #     of species, ki is the rate coefficient for the ith reaction, [Xj]
        #     is the concentration of the jth speies, and vji' and vji" are the
        #     stoichiometric coefficients for reaction i and species j. 'pi
        #     term' references the series of products denoted in the above
        #     equation by the large pi symbol
        
    
        #Calculate the pi term for each species  
        pitermf = np.ones([1,R])
        pitermb = np.ones([1,R]) #initialize the forward and backward pi terms
        for x in range(R): #for each reaction %This iterates about 2000 times!!!
            for a in range(len(rstruct[x]['reactants'])):#for each reactant in the forward direction
                i = expmatrixF[rstruct[x]['reacI'][a],x]
                j = float(X[0,rstruct[x]['reacI'][a]])
                pitermf[0,x] = pitermf[0,x]*(j**i)
                #pitermf[x] = pitermf[x]*X[0,rstruct[x]['reacI'][a]]#^expmatrixF[rstruct[x]['reacI'][a],x]
               #pitermf[x] = pitermf[x]*X[rstruct[x]['reacI'][a]]#^expmatrixF[rstruct[x]['reacI'][a],x]; 
                #include the contribution of this reactant
               # pitermb[x] = pitermb[x]*X[rstruct[x]['prodI'][a]]^expmatrixB[rstruct[x]['prodI'[a],x]
            for a in range(len(rstruct[x]['products'])): #for each reactant in the backward direction (aka the products)
                i = expmatrixB[rstruct[x]['prodI'][a],x]
                j = float(X[0,rstruct[x]['prodI'][a]])
                pitermb[0,x] = pitermb[0,x]*(j**i)
                #j = float(pitermb[x]*X[0,rstruct[x]['prodI'][a]])
                #pitermb[x] = j**i           
                #pitermb[x] = pitermb[x]*X[rstruct[x]['prodI'][a]]^expmatrixB[rstruct[x]['prodI'[a],x]
                #pitermb[x] = pitermb[x]*X(rstruct(x).prodI(a))^expmatrixB(rstruct(x).prodI(a),x);
                #include the contribution of this reactant
                
                
        #Calculate the net production rate of each species
        Wf = Kf*coeffmatrixF*pitermf
        #print(Wf)
        Wb = Kb*coeffmatrixB*pitermb #DOUBLE CHECK THIS TOO
            #Wf and Wb are matrices of the K*coeff*piterm term for each species and reaction
        
        #sum each row of Wf and Wb to get wf and wb for each species
        wf = Wf.sum(axis=1)
        wb = Wb.sum(axis=1)
        
        w = wf + wb #sum the forward and backward rates to get the net %%%%%%%%%%%%% DOUBLE CHECK THIS KELLY !!!!!!
        w = np.transpose(w, axes=None) #transpose w  
        
        sumTerm = 0 #initialize sum term in eq 2
        
        for i in range(ng): #for each species 
            index = lg[i]['index'] #identify location of current species
            sumTerm = sumTerm + massfc[index]*kLG[0,i] #add this species's contribution
        
        #calculate eq (2) for each condensed species
        massfcdot= w*MW*(1/rho) - np.transpose(massfc, axes=None)*kLG + np.transpose(massfc, axes=None)*sumTerm
        
        ### Now look at the gas phase ###
        
        sumTerm = 0 #initialize sum term in eq 6 and 7 %This iterates 11 times!
        for i in range(ng-1): #for each gas phase species except N2
            #index = lg[i]['index']
            #molecW = 
            #sumTerm = sumTerm +(massfc[index]
            sumTerm = sumTerm + (mt*massfc[lg[i]['index']]*klg[0,i])/(lg[i]['MW']) #add the species's contribution to the sum term    
        sumTerm = sumTerm + (mN2dot/0.028014) #add N2 to the sum term
    
        #Calculate eq 6 for each gas phase speices
        molefgdot = np.zeros((1,ng))
        for n in range (ng-1):
            #if lg[n]['MW'] > 0: #added if 9/1
            molefgdot[0,n] = ((Ru*Tg)/(P*Vg*lg[n]['MW']))*(mt*massfc[lg[n]['index']]*klg[0,n] - sumTerm*molefg[n]*lg[n]['MW'])
            #else:
             #   molefgdot[0,n] = 0

            
        #Calculate eq 7 for the purge gas KELLY LOOK AT THIS
        # molefgdot[0,ng-1] = 2#((Ru*Tg)/(P*Vg*lg[ng]['MW']))*(mN2dot - sumTerm*molefg[ng]*lg[ng]['MW'])  
        molefgdot[0,ng-1] = ((Ru*Tg)/(P*Vg*lg[ng-1]['MW']))*(mN2dot - sumTerm*molefg[ng-1]*lg[ng-1]['MW'])
        ### assume nitrogen is the last gas phase species  
        
        ### Now combine the two mass flow rate vectors ###
        #massfs = np.concatenate((massfcdot, molefgdot))
        massfs = np.hstack((massfcdot, molefgdot))
        zdot = np.transpose(massfs, axes=None)  
        #zdot = (massfcdot)
    
        zdot = np.ndarray.tolist(zdot)
        #print(zdot)
        print(t)
        
    
        return zdot

    t1 = 2.0 #set the stop time

    r = solve_ivp(myrhs, (0,t1), y0, method="LSODA", vectorized = True, atol = 1e-10)
    
    
    
    TIME = r['t']
    solutionM = r['y']
    solutionM = solutionM.transpose()
    #solutionM = solutionM[:,0:J]
    #solutionM = solutionM[:,0:J-1]
    
    return [TIME, solutionM]




