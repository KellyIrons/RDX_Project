"""
RDX Reduction
"""
# General set up
Tset = 538
N = 10
R = 1
targets = ['RDX']
reduction_type = 'species'
tol = 5 # percent


from RDXskeletalparser import RDXskeletalparser
from rateConstantCalc2 import rateConstantCalc2
from stoichcalc import stoichcalc
from RDX_sim import RDX_sim
from overall_EP import overall_EP
from rank_EP import rank_EP
from remove_and_remake import remove_and_remake
from RDX_error_calc import RDX_error_calc

import matplotlib.pylab as plt


# Make the original mechanism
[Rlist, Slist, lg] = RDXskeletalparser()

[KF, KB] = rateConstantCalc2(Rlist, Slist, Tset)

[expmatrixF, coeffmatrixF, expmatrixB, coeffmatrixB] = stoichcalc(Rlist,Slist)


# Run the simulation
[TIME_o, solutionM_o] = RDX_sim(Rlist, Slist, lg, Tset, KF, KB, expmatrixF, coeffmatrixF, expmatrixB, coeffmatrixB)
plt.plot(TIME_o, solutionM_o)
plt.ylim((0, 0.15))
plt.show()


# Calculate EP
[EP, mechobj] = overall_EP(N, TIME_o, solutionM_o, Tset, Rlist, Slist, lg, expmatrixF, coeffmatrixF, expmatrixB, coeffmatrixB, targets, reduction_type)


# Rank EP
[ranked_EP, ranked_spec, sort_ind] = rank_EP(EP, mechobj, reduction_type, Rlist)


# Initialize an error measure
error = 0
removed = []
TIME_r = []
solutionM_r = []

# While the error is less than the error tolerance
while error <= tol:
    
    # "Save" the old reduced arrays
    if len(removed) > 0:
        TIME_f = TIME_r
        solutionM_f = solutionM_r
        Slist = new_Slist
        Rlist = new_Rlist
        lg = new_lg

    # Remove a set number of species and remake the mechanism
    [new_Rlist, new_Slist, new_lg, removed, sort_ind] = remove_and_remake(sort_ind, R, removed, Rlist, Slist, lg, reduction_type)
    print('%i Species Removed' % len(removed)) # Change this when I include both species and reaction reduction
    
    
    # Recalculate contant matrices for new mechanism
    [KF, KB] = rateConstantCalc2(new_Rlist, new_Slist, Tset)
    [expmatrixF, coeffmatrixF, expmatrixB, coeffmatrixB] = stoichcalc(new_Rlist,new_Slist)
    
    
    # Re-run the simulation
    [TIME_r, solutionM_r] = RDX_sim(new_Rlist, new_Slist, new_lg, Tset, KF, KB, expmatrixF, coeffmatrixF, expmatrixB, coeffmatrixB)
    
    # Calculae the error
    [error, results, ts] = RDX_error_calc(TIME_o, solutionM_o, solutionM_r, TIME_r, N, removed)
    
    
    import matplotlib.pylab as plt
    plt.plot(ts[0,:], results[0,:])
    plt.show()
    
    
    plt.plot(TIME_r, solutionM_r)
    plt.ylim((0, 0.15))
    plt.show()


removed = removed[:-1]
print('%i Species Removed' % len(removed)) # Change this when I include both species and reaction reduction 


plt.plot(TIME_f, solutionM_f[:,0:len(Slist)])
plt.ylim((0, 0.15))
plt.show()



