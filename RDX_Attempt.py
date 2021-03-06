"""
RDX Attempt 
"""

'''
Set Up
'''

#Import Packages
import numpy as np
#from scipy import integrate
import matplotlib.pylab as mat
from scipy.integrate import solve_ivp #ode CHANGED 


#Define Constants
######A NOTE: I can't pass these through the function header, but I don't think I need to
p = {}
p['Tg'] = 473  #Temperature of gas region - K
p['P'] = 101325 #Atmospheric pressure - Pa
p['Vg'] = 5.3996E-4 #m^3 %note: this is a guess - may be more like 0.0028 if dimensions match reference [38] or 5.4E-4 [32]
p['Ru'] = 8.314 #Gas constant - J/molK
p['mN2dot'] = 1.4583e-06 #mass flow of purge gas (guess) - kg/s
p['Ti'] = 337.6 #initial temp - K
p['Tset'] = 538 #set temp - K
p['roe'] = 1820 #liquid density of RDX - kg/m^3
p['samplemass'] = 0.0005 #this is probably correct
p['Vc'] = p['samplemass']/p['roe']
constantMoles = (p['P']*p['Vg'])/(p['Ru']*p['Tg']); #constant number of moles in gas phase

#Only some species are assumed to evaporate: NO2 [299], N2O [293], 
    #NO [298], HCN [20], CH2O [16], CO2 [18], CO [17], H2O [19], 
    #C_HONO [1], T_HONO [320], RDX [319], CH2NH [7], CH2NCHO [3], N2 [291],
    #TAZ [321], INT86A [266], INT221A [187], INT175A [151], INT131A [100], 
    #INT128A [82], CH2NCN [4], HCNCHO [21], HCOOH [24], NH2CHO [294]

#gasI = [299, 293, 298, 20, 16, 18, 17, 19, 1, 320, 319, 291] #indices of gas phase species
#ng = len(gasI)
gasI = [5,16,3,36,21,29,24,44,7,2,0,23]


#Parse Mechanism
#from RDXparser2 import RDXparser2
#[Rlist, Slist, lg] = RDXparser2()
from RDXskeletalparser import RDXskeletalparser
[Rlist, Slist, lg] = RDXskeletalparser()



#Rlist = Rlist[0]

#REDUCED
'''
Slist[0] = Slist[318]
Slist[1] = Slist[150]
Slist[2] = Slist[319]

lg[0] = lg[9]
lg[1] = lg[10]
lg[0]['index'] = 2
lg[1]['index'] = 0

del Rlist[1:500]
del Slist[3:321]
del lg[2:12]
#del Rlist[1:499]
'''

ng = len(lg)
A = np.zeros((1,ng))
E =  np.zeros((1,ng))
n = np.zeros((1,ng))

for x in range(ng):
    A[0,x] = lg[x]["A"]
    E[0,x] = lg[x]["E"]
    n[0,x] = lg[x]["n"]
    
p['A'] = A
p['E'] = E
p['n'] = n

#REDUCED
'''
Rlist[0]['reacI'] = [0]
Rlist[0]['prodI'] = [1,2]
'''


p['R'] = len(Rlist)
p['J'] = len(Slist)
J = p['J']


MW = np.zeros([1,J])
#code to make MWs, A, E, n vectors
for x in range(J):
    MW[0,x] = Slist[x]['MW']

#Calculate the rate constant for each reaction, in forward and backward direction
from rateConstantCalc2 import rateConstantCalc2
[KF, KB] = rateConstantCalc2(Rlist, Slist, p['Tset'])




#Define stoichiometric coefficent matrices
from stoichcalc import stoichcalc
[expmatrixF, coeffmatrixF, expmatrixB, coeffmatrixB] = stoichcalc(Rlist,Slist)

#Define Initial quantities
#XinitialRDX = 1 #initial mole fraction of RDX in liquid CV is one
#XinitialPurge = 1 #initial mole fraction of purge gas in gas CV is one
YinitialRDX = 1 #initial mass fraction of RDX in liquid CV is one
initialmass = p['samplemass']
#initialVector = np.zeros((1,319+ng))
#initialVector[0,318] = YinitialRDX
#initialVector = [1,0,0]
#initialVector = np.zeros((1,3))
#initialVector[0,0] = 1
initialVector = np.zeros((1,J+ng))
initialVector = initialVector.tolist()
initialVector = initialVector[0]
initialVector[0] = YinitialRDX ####REDUCED
y0, t0 = initialVector, 0



def myrhs(t, y):# arg1):
    
    #Unpack "storage" dict p
    roe = p['roe']; Vg = p['Vg']; mN2dot = p['mN2dot']; P = p['P']; Ru = p['Ru']; 
    J = p['J']; R = p['R'];  Tg = p['Tg']; Vc = p['Vc']; Tset = p['Tset']; Ti = p['Ti'];
    minitial = initialmass;
    A = p['A']; n = p['n']; E = p['E']; #constantMoles = p['constantMoles'];
    rstruct = Rlist; sstruct = Slist
    #MW = p['MW'];

    
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
    e = np.math.exp(1)
    klg = A*(T**n)*(e**(-E/(Ru*T)));
    


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
    massfcdot= w*MW*(1/roe) - np.transpose(massfc, axes=None)*kLG + np.transpose(massfc, axes=None)*sumTerm
    
  
    ### Now look at the gas phase ###
    
    sumTerm = 0 #initialize sum term in eq 6 and 7 %This iterates 11 times!
    for i in range(ng-1): #for each gas phase species except N2
        #index = lg[i]['index']
        #molecW = 
        #sumTerm = sumTerm +(massfc[index]
        sumTerm = sumTerm + (mt*massfc[lg[i]['index']]*klg[0,i])/(lg[i]['MW']) #add the species's contribution to the sum term    
    sumTerm = sumTerm + (mN2dot/0.028) #add N2 to the sum term

    #Calculate eq 6 for each gas phase speices
    molefgdot = np.zeros((1,ng))
    for n in range (ng-1):
        molefgdot[0,n] = ((Ru*Tg)/(P*Vg*lg[n]['MW']))*(mt*massfc[lg[n]['index']]*klg[0,n] - sumTerm*molefg[n]*lg[n]['MW'])
            
    #Calculate eq 7 for the purge gas KELLY LOOK AT THIS
    molefgdot[0,ng-1] = 2#((Ru*Tg)/(P*Vg*lg[ng]['MW']))*(mN2dot - sumTerm*molefg[ng]*lg[ng]['MW'])  
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


#r = ode(myrhs).set_integrator('vode', method='bdf', atol = 1e-7, rtol = 1e-6, nsteps = 10E5)#, max_step = 1e-6) #set the integrator
#r = ode(myrhs).set_integrator('vode', method='bdf', atol = 1e-7, rtol = 1e-6, nsteps = 10E5)Ω
#r.set_initial_value(y0, t0) #.set_f_params(2.0) #set initial values
t1 = 2.0 #set the stop time
dt = 0.010 #set the time interval output in results

r = solve_ivp(myrhs, (0,t1), y0, method="LSODA", vectorized = True)

time = r['t']
solutionM = r['y']
solutionM = solutionM.transpose()
solutionM = solutionM[:,0:J-1]

mat.plot(time, solutionM)
mat.grid('on')
mat.xlabel('Time [seconds]')
mat.ylabel('Mass Fraction')


#increase the number of steps I allow
#
''' CHANGED
num_steps = np.floor((t1 - t0)/dt) + 1
num_steps = int(num_steps)
tarray = np.zeros((num_steps, 1))
#answers = np.zeros((num_steps, 2)) 
answersL = len(y0)
answers = np.zeros((num_steps, answersL)) 
tarray[0] = t0

for x in range(answersL):
    #answers[0,x] = y0[0,x]
    answers[0,x] = y0[x]
'''

'''
#answers[0,0] = y0[0] ####
#answers[0,1] = y0[1]
#answers[0,2] = y0[2]
#answers[0,3] = y0[3]
'''

'''
CHANGED
k = 1
while r.successful() and r.t < t1 and k< num_steps: #Added the last and, not sure if correct
#while r.t < t1:
    #print(r.t+dt, r.integrate(r.t+dt))
    #print(r.t+dt) #debugging
    r.integrate(r.t+dt)
    # Store the results to plot later
    tarray[k] = r.t
    for x in range(answersL):
        answers[k,x] = r.y[x]
       #answers[k,x] = r.y[0,x] ####
    #print(r.t)

    #answers[k,0] = r.y[0] ####
    #answers[k,1] = r.y[1] ####
    #answers[k,2] = r.y[2]
    #answers[k,3] = r.y[3]

    k += 1

A = answers[:,1:3]

mat.plot(tarray, A)
mat.grid('on')
mat.xlabel('Time [seconds]')
mat.ylabel('Mass Fraction')
'''

