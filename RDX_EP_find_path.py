"""
RDX_find_EP_path
"""

def RDX_find_EP_path(ind_list, mechobj, EP, DIC_ind, DIC_spec, targets, coeff_list):
    
    # This is going to be slowwwwww, try to speed it up later
    
    import numpy as np
    import itertools as itertools
    from itertools import permutations
    
    S = len(ind_list)
    
    #find indices of targets
    target_indices = []
    species_names = mechobj['species_names']
    for target in targets:
        index = 0;
        while index <= (S-1) and target != species_names[index]:
            index = index + 1
        target_indices.append(index)
    #print(target_indices)
    
    EP_paths = []
    
    for i in range(S):
        
        #For multiple targets, would have to figure out which one is the one 
        # that the EP coefficient comes from. For now, just take the one target
        
        # Start out by adding the initial point (target) and the final point (species i)
        this_path = [target_indices[0], i]    
        
        
        
        # If the list isn't empty, then the target and and the species have a path
        #  of more than one edge
        DIC = DIC_spec[:, :, int(DIC_ind[i])]
        coeff_list_i = coeff_list[:,:,int(DIC_ind[i])]
        
        if len(ind_list[i]) > 0:
            
            mid_nodes = len(ind_list[i])
            
            temp_ind_list = list(set(ind_list[i]))
            
            
            # Check for repeats in ind_list (this complicates things)
            counts = np.ones(len(temp_ind_list))
            if mid_nodes != len(temp_ind_list):
                for j in (range(len(temp_ind_list))):
                    
                    repeats = ind_list[i].count(temp_ind_list[j])
                    counts[j] = int(repeats)
           
            # Modify the ind list to remove duplicates
            mid_nodes = len(temp_ind_list) #need to remove redunancies if I get this part to work
            ind_list_i = np.copy(temp_ind_list)
           
            counts = np.hstack(([1], np.hstack((counts, [1]))))
           
            print(counts)
            
            # Total number of possible paths, given the node list
            perm = list(permutations(ind_list_i, mid_nodes))
            
            perms = []
            
            l = [False, True]
            dic_direction = list(itertools.product(l, repeat=len(perm)+1))
            d = len(dic_direction)
            
            EPs = np.zeros((len(perm), d))
            
            
            for p in range(len(perm)): # for each possible combination...
                # make the combination and store
                #print(p)
                perms.append([this_path[0]]+list(perm[p])+[this_path[1]])
                print(perms)
                
                for q in range(d):
                   # print(q)
                    coeff = 1
                    # calculate the EP for the combination and store
                    for j in range(len(perms[p])-1):
                        
                        '''
                        # If both the current species and the next species aren't repeated
                        if counts[j] == 1 and counts[j+1] == 1:
                            if dic_direction[q][j]:
                                coeff = coeff*DIC[perms[p][j], perms[p][j+1]]
                            else:
                                coeff = coeff*DIC[perms[p][j+1], perms[p][j]]
                        # If both species are repeated
                        elif counts[j] > 1 and counts[j+1] >1:
                            coeff = coeff_list[j][int(counts[j]-2)]*coeff_list[j+1][int(counts[j+1]-2)]
                            # not sure how necessary this is
                        # Only one of the species is repeated
                        else:
                            coeff = coeff
                        '''
                            
                        # If the next species isn't repeated    
                        if counts[j+1] == 1:
                            if dic_direction[q][j]:
                                coeff = coeff*DIC[perms[p][j], perms[p][j+1]]
                            else:
                                coeff = coeff*DIC[perms[p][j+1], perms[p][j]]
                        # If the next species is repeated
                        elif counts[j+1] > 1:
                            #coeff = coeff*coeff_list_i[perms[p][j+1]][int(counts[j+1]-2)]
                            coeff = coeff*coeff_list_i[j+1][int(counts[j+1]-2)]
                        else:
                            print('oh no!')
                                
                                
                        
                    
                    EPs[p,q] = coeff
                    #print(EPs[p,q])
            # Find the EP that matches the actual EP, that is it!
                # Question: would this also be the max EP? I think so
            ind = np.where(EPs == EP[i])
           # print(i)
            print(ind)
            print(EPs)
            print(EP[i])

            if len(ind[0]) == 1:
                EP_paths.append(perms[int(ind[0])]) # Kelly - check this later
                
            else:
                if np.all(ind[0] == ind[0][0]):
                    EP_paths.append(perms[int(ind[0][0])])
                else:
                    print('oh no: %i \n'% i)
                    
        else: 
            EP_paths.append(this_path)      
            
        
        print(EP_paths[i])
    
    return EP_paths