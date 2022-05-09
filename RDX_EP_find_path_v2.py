"""
RDX_find_EP_path v2
"""

def RDX_find_EP_path2(ind_list, mechobj, EP, DIC_ind, DIC_spec, targets, coeff_list):
    

    
    S = len(ind_list)
    
    #find indices of targets
    target_indices = []
    species_names = mechobj['species_names']
    for target in targets:
        index = 0;
        while index <= (S-1) and target != species_names[index]:
            index = index + 1
        target_indices.append(index)
    
    EP_paths = []
    
    for i in range(S):
        
        #For multiple targets, would have to figure out which one is the one 
        # that the EP coefficient comes from. For now, just take the one target
        
        # Start out by adding the initial point (target) and the final point (species i)
        this_path = [target_indices[0], i]    
        
        
        
        # If the list isn't empty, then the target and and the species have a path
        #  of more than one edge     
        if len(ind_list[i]) > 0:
            
            # Take the necessary data subsets
            DIC = DIC_spec[:, :, int(DIC_ind[i])]
            #coeff_list_i = coeff_list[:,:,int(DIC_ind[i])]
            coeff_list_i = coeff_list[int(DIC_ind[i])]
            
            # Initialize a list of possible paths
 #           complete_paths = []
            paths = [this_path]
            
            
            # Take the last entry of ind_list as the second to last node, at it to the path
            node = ind_list[i][-1]
            paths[0].insert(1,node)
            
            #print(paths)
            
            final_EP = []
            
            flag = True
            if coeff_list_i[node,0] == -1:
                flag = False
                this_EP = DIC[target_indices[0], node]*DIC[node,i]
                flag2 = EP[i] == this_EP
                #print(flag2)
                if flag2:
                    final_EP = paths[0]
                else:
                    print('Error!')

            # If the numer of iterations for this species in local_error_propagation was greater than 1..
            else:
                count = 0 # need to count the number of levels we go through
                flag3 = False
                while flag and not flag3:
                    # Add a level
                    count = count+1
                    if not flag3:
                        for j in range(len(paths)):    # for each possible path                            
                            # For each level, look at the adjacent species as possible next paths
                            next_node = mechobj['S_ind'][node]['other_spec'].copy()
                            for ind in paths[j]:
                                if ind in next_node:
                                    next_node.remove(ind)
                            for k in range(len(next_node)):
                                temp_path = paths[j].copy()
                                temp_path.insert(1, next_node[k])
                                if coeff_list_i[next_node[k],count-1] == -1 and not flag3 :
                                    this_EP = 1
                                    for ind in range(len(temp_path)-1):
                                        this_EP = this_EP*DIC[temp_path[ind], temp_path[ind+1]]
                                    if EP[i] == this_EP:
                                        final_EP = temp_path.copy()
                                        flag3 = True
                                        continue
                                else:
                                    paths.append(temp_path)


        else:
            final_EP = this_path

        EP_paths.append(final_EP)
            
        
       
        
       
    # Checking step
    for i in range(S):
        DIC_i = DIC_ind[i]
        DIC = DIC_spec[:, :, int(DIC_i)]
        path = EP_paths[i]
        
        coeff = 1
        for ind in range(len(path)-1):
            coeff = coeff*DIC[path[ind], path[ind+1]]
        
        flag = EP[i] == coeff
        
        if not flag:
            print('Error: species %i', i)
    
    
    
    return EP_paths