"""
RDX Remove Species and Remake Mechanism
"""
def remove_and_remake(sort_ind, N, removed, Rlist, Slist, lg, reduction_type):
    
    # sort_ind = ranked species indices from rank_EP
    # N = the number of species to be removed *in this iteration*
    # removed = list of the species indices already removed (can be empty)
    
    # Get the original Rlist and Slist --> later can make this a function input
    # to make it more efficient if this function is called multiple times
   # from RDXskeletalparser import RDXskeletalparser
   # [Rlist, Slist, lg] = RDXskeletalparser()
    
    # Initialize the outputs
    new_Rlist = Rlist.copy()
    new_Slist = Slist.copy()
    new_lg = lg.copy()
    
    
    # Identify the species to remove (take the N lowest)
    to_remove = sort_ind[-N:]
    
    if reduction_type in ['S', 'species']:

        for i in reversed(range(len(to_remove))):
            ind = to_remove[i]
            
            print(Slist[ind]['name'])
            
            #Remove the species
            del new_Slist[ind]
        
            count =  0
            for i in range(len(new_Rlist)):
                indices = new_Rlist[i-count]['reacI'] + new_Rlist[i-count]['prodI']
                if ind in indices:
                    del new_Rlist[i-count]
                    count = count + 1
                else: 
                    for j in range(len(new_Rlist[i-count]['reacI'])):
                        if new_Rlist[i-count]['reacI'][j] > ind:
                            new_Rlist[i-count]['reacI'][j] = new_Rlist[i-count]['reacI'][j] - 1
                    for j in range(len(new_Rlist[i-count]['prodI'])):
                        if new_Rlist[i-count]['prodI'][j] > ind:
                            new_Rlist[i-count]['prodI'][j] = new_Rlist[i-count]['prodI'][j] - 1
                            
            count = 0               
            for i in range(len(lg)):
                if lg[i-count]['index'] > ind:
                    lg[i-count]['index'] = lg[i-count]['index'] - 1
                            
            
            removed.append(ind)
            
            count = 0
            sort_ind = sort_ind[0:-1]
            for i in range(len(sort_ind)):
                if sort_ind[i] > ind:
                    sort_ind[i] = sort_ind[i] - 1
            
        
        new_lg = lg
        
    else:
        
        for i in reversed(range(len(to_remove))):
            ind = to_remove[i]
            
            print(Rlist[ind]['eq'])
            
            # Remove the reaction
            del new_Rlist[ind]
            
            # Perform integrity check
            # For each species in the reaction
                # Check to see if it is involved in other reactions
                    # if not, remove the species and any associated reactions
            
            
            

                            
            
            removed.append(ind)
            
            # Update indices in sort_ind
            count = 0
            sort_ind = sort_ind[0:-1]
            for i in range(len(sort_ind)):
                if sort_ind[i] > ind:
                    sort_ind[i] = sort_ind[i] - 1
        
        
    
    return [new_Rlist, new_Slist, new_lg, removed, sort_ind]




