"""
Mechanism Integrity Check
"""

def integrity_check(new_Slist, new_Rlist, new_lg, sort_ind, removed_species, removed_reactions, removed_lg):
    
    from species_removal import species_removal
    
    # Initialize
    found = [False for x in range(len(new_Slist))]
    ind = 0  
    lost_species = []
    lost_ind = []
    reduction_type = "species"
    
    while not all(found):
            
        ind = 0  
        lost_species = []
        lost_ind = []
        
        # Check if each species is present in a reaction
        while not all(found) and ind < len(new_Rlist):
            
            species =  new_Rlist[ind]['reacI'] + new_Rlist[ind]['prodI']
            
            for spec in species:
                found[spec] = True
            
            ind = ind+1
            
        # Remove the species that aren't in any reactions
        if not all(found):
            num = found.count(False)
            
            for i in range(num):
                index = found.index(False)
                
                val = new_Slist.pop(index)
                
                print("Lost Species: %s" % val['name'])
                
                lost_species.append(val)
                lost_ind.append(index)   
                removed_species.append(index)
                [new_Rlist, new_lg, sort_ind, removed_reactions, removed_lg] = species_removal(new_Rlist, index, new_lg, sort_ind, reduction_type, removed_reactions, removed_lg)  
                    
            found = [False for x in range(len(new_Slist))]
                
        #removed_species.append(lost_ind)
                
        #for i in lost_ind:
        #    [new_Rlist, new_lg, sort_ind, removed_reactions, removed_lg] = species_removal(new_Rlist, i, new_lg, sort_ind, reduction_type, removed_reactions, removed_lg)       

    return [new_Slist, new_Rlist, new_lg, sort_ind, removed_species, removed_reactions, removed_lg]
    
    
    