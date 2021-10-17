"""
BASIC PATHWAY ANALYSIS
"""

# General set up
Tset = 538


from RDXskeletalparser import RDXskeletalparser
from rateConstantCalc2 import rateConstantCalc2
from stoichcalc import stoichcalc
from RDX_sim import RDX_sim
from RDXprodrates import RDX_prodrates
from write_pathways import write_pathways


# Make the original mechanism
[Rlist, Slist, lg] = RDXskeletalparser()

[KF, KB] = rateConstantCalc2(Rlist, Slist, Tset)

[expmatrixF, coeffmatrixF, expmatrixB, coeffmatrixB] = stoichcalc(Rlist,Slist)


# Run the simulation
[TIME_o, solutionM_o] = RDX_sim(Rlist, Slist, lg, Tset, KF, KB, expmatrixF, coeffmatrixF, expmatrixB, coeffmatrixB)


# Analyze the pathways
prodrates = RDX_prodrates(Tset, Rlist, Slist, KF, KB, expmatrixF, coeffmatrixF, expmatrixB, coeffmatrixB, TIME_o, solutionM_o)

# Display the pathways
write_pathways(prodrates, Slist, Rlist)
                
        
        
    
    
    