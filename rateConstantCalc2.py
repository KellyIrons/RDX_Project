"""
RATE CONSTANT CALCULATION V2
"""

def rateConstantCalc2(RList, SList, Tset):
    
    
    
    import numpy as np
    
    #This function uses a reaction dict and a species dict to calculate a
    #matrix of the forward and backward rate constants of each reaction at
    #small temperature intervals from the initial temperature to the set
    #temperature
    
    #Define Constants
    rxn = len(RList)
    kB = 1.380649E-23 #Boltzmann constant - J/K 
    h = 6.62607015E-34 #Planck's constant - J*s  
    cref = 1 #standard state concentration - M^(n-1)
    Ru = 8.314 #gas constant - J/molK 
    Na = 6.022E23 #Avogadro's number 
    S = len(SList) #number of species 
    e = np.e
    

    
    sigF=1# np.ones((1,rxn))
    sigB=1#np.ones((1,rxn))
    
    #T = Tset - (Tset-Ti)*exp(-t/0.045);
    
    #Figure out at what time the set temperature is reached
    Ti = 337.6 #initial temp - K
    tarray = np.linspace(0,.75,150) #overestimate the time interval in which Tset is reached
    T = -100000000000000 #initialize the temperature
    count = 0 #initialize a counter
    #while the temperture is smaller than the set temp and the time is valid
    while (T < (Tset - 0.05)) and (count <= len(tarray)):
        t = tarray[count] #set the current time
        T = Tset - (Tset-Ti)*(e**(-t/0.045)) #calculate the current temp
        count = count + 1 #add one to the counter
    setTime = t #the setTime is the time where the set temp is reached
 
    stepSize = 0.001 #difference between temp values in Kf and Kb (can change)
    tarray = np.arange(0, setTime+stepSize, stepSize)
    Tarray = Tset - (Tset-Ti)*(e**(-tarray/0.045)) #calculate temp at each time

    Kf = np.zeros((len(Tarray), rxn))
    Kb = np.zeros((len(Tarray), rxn))   
    
    for i in range(len(Tarray)):
        
        T = Tarray[i]
        
        #Use equation at the very top of page 460 to get dynamic visosity
        visc = ((-0.002551883266)*(T**3) + 5.1897440080*(T**2) - 3531.527767*T
                + 807165.85)*(10**(-6)) #%dynamic viscosity - Pa*s 
        
        kefff = np.zeros((1,rxn))
        keffb = np.zeros((1,rxn))
        
        
        #Calculate diffusion constant for each species at this temp 
        for x in range(S): #for each species
            r = SList[x]['rad'] #unpack molecular radius of a given species 
            SList[x]['D'] = (kB*T)/(6*(np.pi)*visc*r) #calculate D and put into struct - m^2/s

        #Calculate the rate constant for each reaction
        for x in range(rxn): #for each reaction
            w = RList[x]['w'] #unpack w 
            dGf = RList[x]['dGf']; dGb = RList[x]['dGb'] #unpack dGf and dGb  
            reacI = RList[x]['reacI'] #unpack reactant index array 
            prodI = RList[x]['prodI'] #unpack product index array 
                
            #Calculate kkin (eq 4) 
            tautun = 1 + (1/24)*((h*w)/(kB*T))**2 #calculate tautun using eq 3 in  
                #"Quantum mechanics investigation on initial decomposition of
                #ammonia borane in glyme" by Chatterjee and Thynell      
            kkinf = sigF*tautun*((kB*T)/(h*cref))*(e**((-dGf)/(Ru*T))) #kinetic constant for forward direction
            kkinb = sigB*tautun*((kB*T)/(h*cref))*(e**((-dGb)/(Ru*T))) #kinetic constant for backward direction
    
    
                #Calculate diffusion rate constant and overall rate constant based on molecularity 
                    #this section is based off of eqs (5) and (3) in the paper, as well
                    #as pages 3-5 of the supplementary material of "Examination of the
                    #Mechanism of the Yield of N2O from Nitroxyl (HNO) in Solution
                    #Phase by Theoretical Calculations" by Zhang and Thynell
                
            molecL = np.sum(RList[x]['rcoeffs'])
            molecR = np.sum(RList[x]['pcoeffs']) #molecularity of left and right sides
                
            if (molecL == 1 and molecR == 1) or (molecL > 2) or (molecR) > 2: #one product and one reactant 
                kdiffc = 0; kdiffd = 0 #no diffusion 
                kefff[0,x] = kkinf #rate constant is based only on kinetics 
                keffb[0,x] = kkinb
            elif (molecL == 2) and (molecR == 1): #two reactants and one product
                if len(RList[x]['reacI']) == 2:
                    i1 = reacI[0]; i2 = reacI[1]; #unpack species data 
                    r1 = SList[i1]['rad']; r2 = SList[i2]['rad']; 
                    Dp = SList[i1]['D']; Dq = SList[i2]['D'] 
                elif  len(RList[x]['reacI'])== 1:
                    i1 = reacI[0] #unpack species data 
                    r1 = r1 = SList[i1]['rad']; r2 = r1; 
                    Dp = SList[i1]['D'];  Dq = Dp;
                else:
                    raise NameError('error 1')
                R = r1 + r2 #calculate encounter radius 
                V = (4/3)*(np.pi)*R**3 #calculate encounter volume 
                kdiffc = 4*(np.pi)*Na*R*(Dp+Dq); kdiffd = 4*(np.pi)*R*(Dp+Dq)*(1/V) #calculate eq (5) 
                kefff[0,x] = (kdiffc*kkinf)/(kdiffd + kkinf) #apply eq (3) for forward  
                keffb[0,x] = (kdiffd*kkinb)/(kdiffd + kkinf) #apply eq (3) for backward rxn 
            elif (molecL == 1) and (molecR == 2): #one reactant and two products 
                if len(RList[x]['prodI']) == 2:
                    i1 = prodI[0]; i2 = prodI[1] #unpack species data 
                    r1 = SList[i1]['rad']; r2 = SList[i2]['rad']; 
                    Dp = SList[i1]['D']; Dq = SList[i2]['D'] 
                elif len(RList[x]['prodI']) == 1:
                    i1 = prodI[0] #unpack species data 
                    r1 = r1 = SList[i1]['rad']; r2 = r1; 
                    Dp = SList[i1]['D'];  Dq = Dp; 
                else:
                    raise NameError('error 2')
                R = r1 + r2 #calculate encounter radius 
                V = (4/3)*(np.pi)*R**3 #calculate encounter volume                                 
                kdiffc = 4*(np.pi)*Na*R*(Dp+Dq); kdiffd = 4*(np.pi)*R*(Dp+Dq)*(1/V) #calculate eq (5) 
                kefff[0,x] = (kdiffd*kkinf)/(kdiffd + kkinb) #apply eq (3) for forward rxn     
                keffb[0,x] = (kdiffc*kkinb)/(kdiffd + kkinb) #apply eq (3) for backward rxn
            elif (molecL == 2) and (molecR == 2): #two products and reactants
                if len(RList[x]['reacI']) == 2:
                    i1 = reacI[0]; i2 = reacI[1]; #unpack species data 
                    r1 = SList[i1]['rad']; r2 = SList[i2]['rad']; 
                    Da = SList[i1]['D']; Db = SList[i2]['D'] 
                elif len(RList[x]['reacI']) == 1:
                    i1 = reacI[0] #unpack species data 
                    r1 = r1 = SList[i1]['rad']; r2 = r1
                    Da = SList[i1]['D'];  Db = Da
                else:
                    raise NameError('error 3')
                if len(RList[x]['prodI']) == 2:
                    i3 = prodI[0]; i4 = prodI[1] #unpack species data 
                    r3 = SList[i3]['rad']; r4 = SList[i4]['rad'];
                    Dc = SList[i3]['D']; Dd = SList[i4]['D'] 
                elif len(RList[x]['prodI']) == 1:
                    i3 = prodI[0]
                    r3 = SList[i3]['rad']; r4 = r3
                    Dc = SList[i3]['D'];  Dd = Dc                
                else:
                    raise NameError('error 4')
                R1 = r1 + r2; R2 = r3 + r4; #calculate encounter radii 
                V1 = (4/3)*(np.pi)*R1**3; V2 = (4/3)*(np.pi)*R2**3; #calculate encounter volumes 
                k1 = 4*(np.pi)*Na*R1*(Da+Db); kneg1 = 4*(np.pi)*R1*(Da+Db)*(1/V1) #calculate eq (5) for foward rxn 
                k2 = 4*(np.pi)*Na*R2*(Dc+Dd); kneg2 = 4*(np.pi)*R2*(Dc+Dd)*(1/V2) #alculate eq (5) for backward rxn 
                kefff[0,x] = (k1*kkinf*kneg2)/(kneg1*(kneg2+kkinb) + kkinf*kneg2) #apply eq (3) for forward rxn 
                keffb[0,x] = (kneg1*kkinb*k2)/(kneg1*(kneg2+kkinb) + kkinf*kneg2) #apply eq (3) for backward rxn 
                kdiffc = 0.5*(k1 + k2); kdiffd = 0.5*(kneg1 + kneg2);
            else:
                 raise NameError('ERROR!')
                    
        Kf[i,:] = kefff
        Kb[i,:] = keffb
        
        
    Tarraytrans = Tarray[np.newaxis]
            
    Kf = np.concatenate((Tarraytrans.T,Kf), axis=1)
    Kb = np.concatenate((Tarraytrans.T,Kb), axis=1)
    
    #print(Tarraytrans)
        
    return [Kf, Kb]
    
    
    
    
    
    