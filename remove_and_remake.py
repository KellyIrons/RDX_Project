"""
RDX Remove Species and Remake Mechanism
"""
def remove_and_remake(sort_ind, N, removed):
    
    # sort_ind = ranked species indices from rank_EP
    # N = the number of species to be removed *in this iteration*
    # removed = list of the species indices already removed (can be empty)
    
    # Get the original Rlist and Slist --> later can make this a function input
    # to make it more efficient if this function is called multiple times
    from RDXskeletalparser import RDXskeletalparser
    [Rlist, Slist, lg] = RDXskeletalparser()
    
    # Initialize the outputs
    new_Rlist = Rlist
    new_Slist = Slist
    new_lg = lg
    
    
    # Identify the species to remove (take the N lowest)
    if len(removed) != 0:
        to_remove = sort_ind[-N+-len(removed):-len(removed)] 
    else:
        to_remove = sort_ind[-N:]
    #to_remove = sort_ind[-N - removed: -removed]
    
    
    # Go through each species and "remove" it from the mechanism
    ''' #OLD METHOD
    for i in reversed(range(len(to_remove))):
        
        ind = to_remove[i]

        # Remove from Slist
        new_Slist[ind]['MW'] = 0
        new_Slist[ind]['rad'] = 0
        new_Slist[ind]['D'] = 0
        
        # Remove from Rlist
        for r in range(len(new_Rlist)):
            if ind in new_Rlist[r]['reacI']:
                pos = new_Rlist[r]['reacI'].index(ind)
                new_Rlist[r]['rcoeffs'].pop(pos)
                new_Rlist[r]['reactants'].pop(pos)
                new_Rlist[r]['reacI'].remove(ind)
                
            elif ind in Rlist[r]['prodI']:
                pos = Rlist[r]['prodI'].index(ind)
                new_Rlist[r]['pcoeffs'].pop(pos)
                new_Rlist[r]['products'].pop(pos)
                new_Rlist[r]['prodI'].remove(ind)
        
        # Remove from lg
    
        for i in range(len(new_lg)):
            if ind == new_lg[i]['index']:
                new_lg[i]['A'] = 0
                new_lg[i]['n'] = 0
                new_lg[i]['E'] = 0
                new_lg[i]['MW'] = 0
     
        '''
        
        
        # NEW METHOD
    for i in reversed(range(len(to_remove))):
        
        # Identify the species to remove
        ind = to_remove[i]
        
        # Identify the reactions containing that species and set the 
        # appropriate reaction parameters to empty
        for r in range(len(new_Rlist)):
            if ind in new_Rlist[r]['reacI'] or ind in Rlist[r]['prodI']:
                new_Rlist[r]['reacI'] = []
                new_Rlist[r]['prodI'] = []
                new_Rlist[r]['reactants'] = []
                new_Rlist[r]['products'] = []
        
        
        removed.append(ind)
    
    return [new_Rlist, new_Slist, new_lg, removed]
