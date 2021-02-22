"""
STOICH CALC
"""

def stoichcalc(Rlist,Slist):
    
    #This function takes a species and reaction struct created by RDXparser.py
    #(or similar) and examines the chemical equations to create matrices of the
    #coefficients needed for writing the net molar production rate of each
    #species. The matrices generated are denote the exponent and coefficient,
    #for each species for each reaction in each direction. For example, the
    #equations "A + B <=> C" and "2 C <=> B" would generate the following
    #matrices:       R1  R2                           R1  R2
    # expmatrixF = A [ 1,  0;       coeffmatrixF = A [ -1,  0;
    #              B   1,  0;                      B   -1,  1;
    #              C   0,  2]                      C    1,  -2]
    #                R1  R2                           R1  R2
    # expmatrixB = A [ 0,  0;       coeffmatrixB = A [ 1,  0;
    #              B   0,  1;                      B   1,  -1;
    #              C   1,  0]                      C   -1,  2]
    #These matrices result in the molar production rates:
    #   wf(A) = -k1f*[A]*[B]                  wb(A) = k1b*[C]
    #   wf(B) = -k1*[A]*[B] + k2f*[C]^2       wb(B) = k1b*[C]] - k2b*[B]
    #   wf(C) = k1f*[A]*[B]  - 2*k2f*[C]^2     wb(C) = -k1b*[C] + 2*k2b*[B]
    
    
    import numpy as np
    
    #Define Constants (the shouldn't change) 
    S = len(Slist) #number of species 
    rxn = len(Rlist) #number of reactions 
    

    #Calcualte the coefficient matrices
    expmatrixF = np.zeros((S, rxn)) #initialize exp matrix in forward direction
    coeffmatrixF = np.zeros((S, rxn)) #initialize coeff matrix in forward direction
    for a in range(S): #for each species 
        for b in range(rxn): #for each reaction
            reactcount = 0 #inititialize counter of reactant coefficient
            prodcount = 0 #inititialize counter of product coefficient
            for i in range(len(Rlist[b]['reactants'])): #for each reactant
                #print(Slist[a]['name'])
                #print(Rlist[b]['reactants'])
                if Slist[a]['name'] == Rlist[b]['reactants'][i]:
                #if strcmp(q(a).name, p(b).reactants{i}) #if the current species is in the reactants
                    reactcount = reactcount + Rlist[b]['rcoeffs'][i] #add its coefficient to the counter
            for i in range(len(Rlist[b]['products'])): #for each product
                #if strcmp(q(a).name, p(b).products{i}) %if the current species is in the products
                if Slist[a]['name'] == Rlist[b]['products'][i]:
                    prodcount = prodcount + Rlist[b]['pcoeffs'][i]  #add its coefficient to the counter
            expmatrixF[a,b] = reactcount #set the exponent for this species and reaction
            coeffmatrixF[a,b] = prodcount - reactcount #set the coefficient for this species and reaction


    expmatrixB = np.zeros((S, rxn)) #initialize exp matrix in backward direction
    coeffmatrixB = np.zeros((S, rxn)) #initialize coeff matrix in backward direction
    for a in range(S): #for each species 
        for b in range(rxn): #for each reaction
            reactcount = 0 #inititialize counter of reactant coefficient
            prodcount = 0 #inititialize counter of product coefficient
            for i in range(len(Rlist[b]['products'])): #for each product
                #if strcmp(q(a).name, p(b).products{i}) %if the current species is in the products
                if Slist[a]['name'] == Rlist[b]['products'][i]:
                    reactcount = reactcount + Rlist[b]['pcoeffs'][i] #add its coefficient to the counter
            for i in range(len(Rlist[b]['reactants'])): #for each reactant
                #if strcmp(q(a).name, p(b).reactants{i}) %if the current species is in the reactants
                if Slist[a]['name'] == Rlist[b]['reactants'][i]:
                    prodcount = prodcount + Rlist[b]['rcoeffs'][i] #add its coefficient to the counter
            expmatrixB[a,b] = reactcount #set the exponent for this species and reaction
            coeffmatrixB[a,b] = prodcount - reactcount #set the coefficient for this species and reaction        




    return [expmatrixF, coeffmatrixF, expmatrixB, coeffmatrixB]