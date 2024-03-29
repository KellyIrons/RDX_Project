"""
RDX Reduction
"""
##### SET UP #################################################################


##### User Inputs #####

Tset = 538 # Kelvin - set temperature
N = 10 # 
R = 1 # number of species/reactions to remove at a time 
targets = ['RDX']
s_tols = [7, 10, 15] # percent - absolute error tolerance
r_tols = [7, 10, 15] 
s_max_tols = [10, 15, 25] # percent - maximum error tolerance that will be temporarily accepted
r_max_tols =  [10, 15, 25]
initial_species = 'RDX'
atol = 1e-12
graph = False
NS = 140

mech_type = "full"

##### Import Functions #####

from RDXfullparser import RDXfullparser
from rateConstantCalc2 import rateConstantCalc2
from stoichcalc import stoichcalc
from RDX_sim2 import RDX_sim2
from RDX_sim import RDX_sim
from overall_EP import overall_EP
from rank_EP import rank_EP
from remove_and_remake import remove_and_remake
from RDX_error_calc import RDX_error_calc
from integrity_check import integrity_check
from reshape_solutionM import reshape_solutionM
from identify_secondary_targets import identify_secondary_targets 
import matplotlib.pylab as plt
import copy


##### Initial Simulation Run #####

# Make the original mechanism
[Rlist, Slist, lg] = RDXfullparser()
J = len(Slist)
[KF, KB] = rateConstantCalc2(Rlist, Slist, Tset)
[expmatrixF, coeffmatrixF, expmatrixB, coeffmatrixB] = stoichcalc(Rlist,Slist)


# Run the simulation
print("Running the simulation...thank you for your patience...")
[TIME_o, solutionM_o] = RDX_sim2(Rlist, Slist, lg, Tset, KF, KB, expmatrixF, coeffmatrixF, expmatrixB, coeffmatrixB, mech_type, initial_species, atol)

# Plot results
plt.plot(TIME_o, solutionM_o)
plt.grid()
plt.ylim((0, 0.15))
plt.xlabel("Time [s]")
plt.ylabel("Mass Fraction")
plt.title("Original Mechanism")
plt.show()


# Initialize an error measure
error = 0
removed_species = [] 
removed_reactions = [] 
removed_lg = [] 
TIME_r = []
solutionM_r = []

TIME_f = TIME_r.copy()
solutionM_f = solutionM_r.copy()
#Slist_f = Slist.copy()
#Rlist_f = Rlist.copy()
#lg_f = lg.copy()
Slist_f = copy.deepcopy(Slist)
Rlist_f = copy.deepcopy(Rlist)
lg_f = copy.deepcopy(lg)
removed_species_f = removed_species.copy()
removed_reactions_f = removed_reactions.copy()
removed_lg_f = removed_lg.copy()
error_f = error
count_d = 0


# Intializing reduction framework parameters
i = 0
j = 0
species_done = False
reactions_done = False
if mech_type == "full":
    EPmin  = 1E-35
elif mech_type == "skeletal":
    EPmin = 1E-15


'''Manual Reduction'''
reduction_type = 'species'
    
 # Calculate EP
[EP, mechobj] = overall_EP(N, TIME_o, solutionM_o, Tset, Rlist, Slist, lg, expmatrixF, coeffmatrixF, expmatrixB, coeffmatrixB, targets, reduction_type, EPmin, graph)

[ranked_EP, ranked_spec, sort_ind] = rank_EP(EP, mechobj, reduction_type, Rlist)

sort_ind_o_s = sort_ind

[new_Rlist, new_Slist, new_lg, removed_species, removed_reactions, removed_lg, sort_ind] = remove_and_remake(sort_ind, 150, removed_species, removed_reactions, removed_lg, Rlist, Slist, lg, reduction_type)
[new_Slist, new_Rlist, new_lg, sort_ind, removed_species, removed_reactions, removed_lg] = integrity_check(new_Slist, new_Rlist, new_lg, sort_ind, removed_species, removed_reactions, removed_lg)
        
 # Recalculate contant matrices for new mechanism
[KF, KB] = rateConstantCalc2(new_Rlist, new_Slist, Tset)
[expmatrixF, coeffmatrixF, expmatrixB, coeffmatrixB] = stoichcalc(new_Rlist,new_Slist)
        
        
# Re-run the simulation
print("Running the simulation...thank you for your patience...")
[TIME_r, solutionM_r] = RDX_sim2(new_Rlist, new_Slist, new_lg, Tset, KF, KB, expmatrixF, coeffmatrixF, expmatrixB, coeffmatrixB, mech_type, initial_species, atol)
solutionM_o = reshape_solutionM(solutionM_o)
solutionM_r = reshape_solutionM(solutionM_r)

        

# Calculate the error
[error, results, ts] = RDX_error_calc(TIME_o, solutionM_o, solutionM_r, TIME_r, N, removed_species, removed_lg, J)

TIME_f = TIME_r.copy()
solutionM_f = solutionM_r.copy()

Slist_f = copy.deepcopy(new_Slist)
Rlist_f = copy.deepcopy(new_Rlist)
lg_f = copy.deepcopy(new_lg)
removed_species_f = removed_species.copy()
removed_reactions_f = removed_reactions.copy()
removed_lg_f = removed_lg.copy()
error_f = error
print("Saved Mechanism!")

Slist = copy.deepcopy(new_Slist)
Rlist = copy.deepcopy(new_Rlist)
lg = copy.deepcopy(new_lg)

plt.grid()
plt.xlabel("Time [s]")
plt.ylabel("Median % Error")
plt.title("Reduced Mechanism Error: %i Species & %i Reactions" % (len(new_Slist), len (new_Rlist)))
plt.show()

        
plt.plot(TIME_r, solutionM_r)
plt.grid()
plt.ylim((0, 0.15))
plt.xlabel("Time [s]")
plt.ylabel("Mass Fraction")
plt.title("Reduced Mechanism: %i Species & %i Reactions" % (len(new_Slist), len (new_Rlist)))
plt.show()
        
print("Error = %f %%" % error)



##### OVERALL LOOP ############################################################
count = 0
while not species_done or not reactions_done: # stop when BOTH conditions are true

    ns_o = len(Slist)
    nr_o = len(Rlist)

    if count <= (len(s_tols)-1):
        s_tol = s_tols[count]
        r_tol = r_tols[count]
        s_max_tol = s_max_tols[count]
        r_max_tol = r_max_tols[count]
    else:
        s_tol = s_tols[-1]
        r_tol = r_tols[-1]
        s_max_tol = s_max_tols[-1]
        r_max_tol = r_max_tols[-1]
    

    ### Species Loop ###
    
    ## Calculate Coefficients
    
    reduction_type = 'species'
    
    # Calculate EP
    [EP, mechobj] = overall_EP(N, TIME_o, solutionM_o, Tset, Rlist, Slist, lg, expmatrixF, coeffmatrixF, expmatrixB, coeffmatrixB, targets, reduction_type, EPmin, graph)
    
    
    # Rank EP
    [ranked_EP, ranked_spec, sort_ind] = rank_EP(EP, mechobj, reduction_type, Rlist)
    
    sort_ind_o_r = sort_ind
    
    
    ## Remove Species
    ns = [len(Slist)]
    nr = [len(Rlist)]
    count_s = 0
    
    while not species_done: 
        
        # "Save" the old reduced arrays
        if count_s > 0:
            '''
            Slist = new_Slist.copy()
            Rlist = new_Rlist.copy()
            lg = new_lg.copy()
            '''
            Slist = copy.deepcopy(new_Slist)
            Rlist = copy.deepcopy(new_Rlist)
            lg = copy.deepcopy(new_lg)
            
        if len(Slist) > 150:
            R = 5
        else:
            R = 1
    
        # Remove a set number of species and remake the mechanism
        [new_Rlist, new_Slist, new_lg, removed_species, removed_reactions, removed_lg, sort_ind] = remove_and_remake(sort_ind, R, removed_species, removed_reactions, removed_lg, Rlist, Slist, lg, reduction_type)
        
        # Integrity Check
        [new_Slist, new_Rlist, new_lg, sort_ind, removed_species, removed_reactions, removed_lg] = integrity_check(new_Slist, new_Rlist, new_lg, sort_ind, removed_species, removed_reactions, removed_lg)
        
        
        
        
        # Recalculate contant matrices for new mechanism
        [KF, KB] = rateConstantCalc2(new_Rlist, new_Slist, Tset)
        [expmatrixF, coeffmatrixF, expmatrixB, coeffmatrixB] = stoichcalc(new_Rlist,new_Slist)
        
        
        # Re-run the simulation
        print("Running the simulation...thank you for your patience...")
        if len(new_Slist) >= NS: 
            [TIME_r, solutionM_r] = RDX_sim2(new_Rlist, new_Slist, new_lg, Tset, KF, KB, expmatrixF, coeffmatrixF, expmatrixB, coeffmatrixB, mech_type, initial_species, atol)
            solutionM_r = reshape_solutionM(solutionM_r)
        else:
            [TIME_r, solutionM_r] = RDX_sim(new_Rlist, new_Slist, new_lg, Tset, KF, KB, expmatrixF, coeffmatrixF, expmatrixB, coeffmatrixB, mech_type, initial_species, atol)    
        
        # Calculae the error
        #[error, results, ts] = RDX_error_calc(TIME_o, solutionM_o, solutionM_r, TIME_r, N, removed_species, removed_lg, J)
        failed_run = TIME_r[-1] < 0.75*TIME_o[-1]
        if failed_run:
            atol = 1e-17
            print("Running the simulation...thank you for your patience...")
            [TIME_r, solutionM_r] = RDX_sim2(new_Rlist, new_Slist, new_lg, Tset, KF, KB, expmatrixF, coeffmatrixF, expmatrixB, coeffmatrixB, mech_type, initial_species, atol)
            solutionM_r = reshape_solutionM(solutionM_r)  
            atol = 1e-15
            failed_run = TIME_r[-1] < 0.75*TIME_o[-1]
            if failed_run:
                error = 1000
            else:
                # Calculae the error
                [error, results, ts] = RDX_error_calc(TIME_o, solutionM_o, solutionM_r, TIME_r, N, removed_species, removed_lg, J)
        else:
            # Calculae the error
            [error, results, ts] = RDX_error_calc(TIME_o, solutionM_o, solutionM_r, TIME_r, N, removed_species, removed_lg, J)       
    
        if error <= s_tol:
            TIME_f = TIME_r.copy()
            solutionM_f = solutionM_r.copy()
            '''
            Slist_f = new_Slist.copy()
            Rlist_f = new_Rlist.copy()
            lg_f = new_lg.copy()
            '''
            Slist_f = copy.deepcopy(new_Slist)
            Rlist_f = copy.deepcopy(new_Rlist)
            lg_f = copy.deepcopy(new_lg)
            removed_species_f = removed_species.copy()
            removed_reactions_f = removed_reactions.copy()
            removed_lg_f = removed_lg.copy()
            error_f = error
            print("Saved Mechanism!")
        
        plt.plot(ts[0,:], results[0,:])
        plt.grid()
        plt.xlabel("Time [s]")
        plt.ylabel("Median % Error")
        plt.title("Reduced Mechanism Error: %i Species & %i Reactions" % (len(new_Slist), len (new_Rlist)))
        plt.show()
        
        
        plt.plot(TIME_r, solutionM_r)
        plt.grid()
        plt.ylim((0, 0.15))
        plt.xlabel("Time [s]")
        plt.ylabel("Mass Fraction")
        plt.title("Reduced Mechanism: %i Species & %i Reactions" % (len(new_Slist), len (new_Rlist)))
        plt.show()
    
        i = i+1
        ns.append(len(new_Slist))
        nr.append(len(new_Rlist))
        
        print("Error = %f %%" % error)
        
        count_s = count_s + 1
    
        #if ns[i] < ns[i-1]: #leads to oo loop :(
        if error <= s_max_tol:
            species_done = False
            reactions_done = False
        else:
            '''
            Slist = Slist_f.copy()
            Rlist = Rlist_f.copy()
            lg = lg_f.copy()
            '''
            Slist = copy.deepcopy(Slist_f)
            Rlist =  copy.deepcopy(Rlist_f)
            lg = copy.deepcopy(lg_f)
            removed_species = removed_species_f.copy()
            removed_reactions = removed_reactions_f.copy()
            removed_lg = removed_lg_f.copy()
            error = error_f
            species_done = True
        
            
            print("Exiting Species Loop")
 
    
    ### Reaction Loop ###
    
    
    reduction_type = 'reactions'
    
    # Recalculate coefficent matrices using species reduced mechanism
    [KF, KB] = rateConstantCalc2(Rlist, Slist, Tset)
    [expmatrixF, coeffmatrixF, expmatrixB, coeffmatrixB] = stoichcalc(Rlist,Slist)
    
    # Calculate EP
    [EP, mechobj] = overall_EP(N, TIME_o, solutionM_o, Tset, Rlist, Slist, lg, expmatrixF, coeffmatrixF, expmatrixB, coeffmatrixB, targets, reduction_type, EPmin, graph)
    
    
    # Rank EP
    [ranked_EP, ranked_reac, sort_ind] = rank_EP(EP, mechobj, reduction_type, Rlist)
    
    
    ## Remove Reactions
    ns = [len(Slist)]
    nr = [len(Rlist)]
    count_r = 0
    
    while not reactions_done:# and error <= max_tol:
        

        # "Save" the old reduced arrays
        if count_r > 0:
            '''
            Slist = new_Slist
            Rlist = new_Rlist
            lg = new_lg
            '''
            Slist = copy.deepcopy(new_Slist)
            Rlist = copy.deepcopy(new_Rlist)
            lg = copy.deepcopy(new_lg)
        
        if len(Slist) > 150:
            R = 5
        else:
            R = 1
        
    
        # Remove a set number of species and remake the mechanism
        [new_Rlist, new_Slist, new_lg, removed_species, removed_reactions, removed_lg, sort_ind] = remove_and_remake(sort_ind, R, removed_species, removed_reactions, removed_lg, Rlist, Slist, lg, reduction_type)
        
        
        # Recalculate contant matrices for new mechanism
        [KF, KB] = rateConstantCalc2(new_Rlist, new_Slist, Tset)
        [expmatrixF, coeffmatrixF, expmatrixB, coeffmatrixB] = stoichcalc(new_Rlist,new_Slist)
        
        
        # Re-run the simulation
        print("Running the simulation...thank you for your patience...")
        #[TIME_r, solutionM_r] = RDX_sim2(new_Rlist, new_Slist, new_lg, Tset, KF, KB, expmatrixF, coeffmatrixF, expmatrixB, coeffmatrixB, mech_type, initial_species, atol)
        #solutionM_r = reshape_solutionM(solutionM_r)
        
        if len(new_Slist) >= NS: 
            [TIME_r, solutionM_r] = RDX_sim2(new_Rlist, new_Slist, new_lg, Tset, KF, KB, expmatrixF, coeffmatrixF, expmatrixB, coeffmatrixB, mech_type, initial_species, atol)
            solutionM_r = reshape_solutionM(solutionM_r)
        else:
            [TIME_r, solutionM_r] = RDX_sim(new_Rlist, new_Slist, new_lg, Tset, KF, KB, expmatrixF, coeffmatrixF, expmatrixB, coeffmatrixB, mech_type, initial_species, atol)    
        
        
        failed_run = TIME_r[-1] < 0.75*TIME_o[-1]
        if failed_run:
            atol = 1e-17
            print("Running the simulation...thank you for your patience...")
            [TIME_r, solutionM_r] = RDX_sim2(new_Rlist, new_Slist, new_lg, Tset, KF, KB, expmatrixF, coeffmatrixF, expmatrixB, coeffmatrixB, mech_type, initial_species, atol)
            solutionM_r = reshape_solutionM(solutionM_r)  
            atol = 1e-15
            failed_run = TIME_r[-1] < 0.75*TIME_o[-1]
            if failed_run:
                error = 1000
            else:
                # Calculae the error
                [error, results, ts] = RDX_error_calc(TIME_o, solutionM_o, solutionM_r, TIME_r, N, removed_species, removed_lg, J)
        else:
            # Calculae the error
            [error, results, ts] = RDX_error_calc(TIME_o, solutionM_o, solutionM_r, TIME_r, N, removed_species, removed_lg, J)
        
        if error <= r_tol:
            TIME_f = TIME_r
            solutionM_f = solutionM_r
            '''
            Slist_f = new_Slist
            Rlist_f = new_Rlist
            lg_f = new_lg
            '''
            Slist_f = copy.deepcopy(new_Slist)
            Rlist_f = copy.deepcopy(new_Rlist)
            lg_f = copy.deepcopy(new_lg)
            removed_species_f = removed_species.copy()
            removed_reactions_f = removed_reactions.copy()
            removed_lg_f = removed_lg.copy()
            error_f = error
            print("Saved Mechanism!")
        
        
        plt.plot(ts[0,:], results[0,:])
        plt.grid()
        plt.xlabel("Time [s]")
        plt.ylabel("Median %% Error")
        plt.title("Reduced Mechanism Error: %i Species & %i Reactions" % (len(new_Slist), len (new_Rlist)))
        plt.show()
        
        
        plt.plot(TIME_r, solutionM_r)
        plt.grid()
        plt.ylim((0, 0.15))
        plt.xlabel("Time [s]")
        plt.ylabel("Mass Fraction")
        plt.title("Reduced Mechanism: %i Species & %i Reactions" % (len(new_Slist), len (new_Rlist)))
        plt.show()
    
        j = j+1
        ns.append(len(new_Slist))
        nr.append(len(new_Rlist))
    
    
        print("Error = %f %%" % error)
        
        count_r = count_r + 1
        
        if error <= r_max_tol:
            species_done = False
            reactions_done = False
        else:
            reactions_done = True
            '''
            Slist = Slist_f
            Rlist = Rlist_f
            lg = lg_f
            '''
            Slist = copy.deepcopy(Slist_f)
            Rlist =  copy.deepcopy(Rlist_f)
            lg = copy.deepcopy(lg_f)
            removed_species = removed_species_f.copy()
            removed_reactions = removed_reactions_f.copy()
            removed_lg = removed_lg_f.copy()
            error = error_f
            [expmatrixF, coeffmatrixF, expmatrixB, coeffmatrixB] = stoichcalc(Rlist,Slist) #need to do this so species reduction works
            
            print("Exiting Reation Loop")

    count = count+1
    
    
    if species_done and reactions_done:
        print('KELLY LOOK!')
    #if len(Slist) == len(Slist_f) and len(Rlist) == len(Rlist_f):
        if count_d <= 3:
            # Further reduction ?
            t = 1.75
            reduction_type = 'species'
            [expmatrixF, coeffmatrixF, expmatrixB, coeffmatrixB] = stoichcalc(Rlist_f,Slist_f)
            [EP, mechobj] = overall_EP(N, TIME_f, solutionM_f, Tset, Rlist_f, Slist_f, lg_f, expmatrixF, coeffmatrixF, expmatrixB, coeffmatrixB, targets, reduction_type, EPmin, graph)
            [g_mean, sort_ind_2] = identify_secondary_targets(TIME_f, solutionM_f, targets, t, Tset, Rlist_f, Slist_f, lg_f, expmatrixF, coeffmatrixF, expmatrixB, coeffmatrixB, EP)  
            targets.append(Slist_f[sort_ind_2[0]]['name'])
            count_d = count_d + 1
            
            species_done = False
            reactions_done = False
            
    
    
    
    
    

# Finishing Steps


print('%i Species Remaining!' % len(Slist))

print('%i Reactions Remaining!' % len(Rlist))


plt.plot(TIME_f, solutionM_f[:,0:len(Slist)])
plt.grid()
plt.ylim((0, 0.15))
plt.xlabel("Time [s]")
plt.ylabel("Mass Fraction")
plt.title("Final Mechanism")
plt.show()

