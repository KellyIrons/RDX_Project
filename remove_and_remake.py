"""
RDX Remove Species and Remake Mechanism
"""
def remove_and_remake(sort_ind, R, removed_species, removed_reactions, removed_lg, Rlist, Slist, lg, reduction_type):
    
    # sort_ind = ranked species indices from rank_EP
    # N = the number of species to be removed *in this iteration*
    # removed = list of the species indices already removed (can be empty)
    
    # Get the original Rlist and Slist --> later can make this a function input
    # to make it more efficient if this function is called multiple times
   # from RDXskeletalparser import RDXskeletalparser
   # [Rlist, Slist, lg] = RDXskeletalparser()
    
   
    from species_removal import species_removal
    import copy
    
    # Initialize the outputs
    new_Rlist = copy.deepcopy(Rlist)
    new_Slist = copy.deepcopy(Slist)
    new_lg = copy.deepcopy(lg)
    
    
    # Identify the species to remove (take the N lowest)
    to_remove = sort_ind[-R:]
    
    if reduction_type in ['S', 'species']:

        for i in reversed(range(len(to_remove))):
            ind = to_remove[i]
            
            print('Removed Species: %s' % new_Slist[ind]['name'])
            
            #Remove the species
            del new_Slist[ind]
         
            [new_Rlist, new_lg, sort_ind, removed_reactions, removed_lg] = species_removal(new_Rlist, ind, new_lg, sort_ind, reduction_type, removed_reactions, removed_lg)
            
        
            sort_ind = sort_ind[0:-1]
            for j in range(len(sort_ind)):
                if sort_ind[j] > ind:
                    sort_ind[j] = sort_ind[j] - 1
            
            
            removed_species.append(ind)
            
        
    elif reduction_type in ['reactions', 'R']:
        
        for i in reversed(range(len(to_remove))):
            ind = to_remove[i]
            
            
            species = Rlist[ind]['reacI'] + Rlist[ind]['prodI'] # all the species in this reaction
            
            # Remove the reaction
            print('Removed Reaction: %s' % new_Rlist[ind]['eq'])
            del new_Rlist[ind]
            
            removed_reactions.append(ind)
            
            sort_ind = sort_ind[0:-1]
            for j in range(len(sort_ind)):
                if sort_ind[j] > ind:
                    sort_ind[j] = sort_ind[j] - 1
            
            # Perform integrity check
            # For each species in the reaction
                # Check to see if it is involved in other reactions
                    # if not, remove the species and any associated reactions
            
            match = False
            found_species1 = []
            found_species2 = []
            while not match:
                for j in species: #need to do this until nothing changes
                    found = False
                    r = 0
                    while not found and r <= (len(new_Rlist)-1):
                        r_spec = new_Rlist[r]['reacI'] + new_Rlist[r]['prodI']
                        found = j in r_spec
                        r += 1
                    if not found:
                        print('Removed Species: %s' % new_Slist[j]['name'])
                        removed_species.append(j)
                        del new_Slist[j]
                        species.remove(j) # Added 1/11
                        # This is where I am going wrong, not changing sort_ind unless not found
                        [new_Rlist, new_lg, sort_ind, removed_reactions, removed_lg] = species_removal(new_Rlist, j, new_lg, sort_ind, reduction_type, removed_reactions, removed_lg)
                        
                    else:
                        found_species2.append(j) #note: this is in the original indexing

                match = found_species1 == found_species2
                found_species1 = found_species2

        
    
    return [new_Rlist, new_Slist, new_lg, removed_species, removed_reactions, removed_lg, sort_ind]




