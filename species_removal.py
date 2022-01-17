"""
Species Sanity Check
"""

def species_removal(new_Rlist, ind, new_lg, sort_ind, reduction_type, removed_reactions, removed_lg):
    
    import copy
    
    removed_Rs = []
    count =  0
    for i in range(len(new_Rlist)):
        indices = new_Rlist[i-count]['reacI'] + new_Rlist[i-count]['prodI']
        if ind in indices:
            print('Removed Reaction: %s' % new_Rlist[i-count]['eq'])
            removed_reactions.append(i-count)
            del new_Rlist[i-count]
            if reduction_type in ['R', 'reactions']:
                removed_Rs.append(i-count)
                
                for j in range(len(sort_ind)):
                    if sort_ind[j] > (i-count):
                        sort_ind[j] = sort_ind[j] - 1
                    
            count = count + 1
     
        else: 
            for j in range(len(new_Rlist[i-count]['reacI'])):
                if new_Rlist[i-count]['reacI'][j] > ind:
                    new_Rlist[i-count]['reacI'][j] = new_Rlist[i-count]['reacI'][j] - 1
            for j in range(len(new_Rlist[i-count]['prodI'])):
                if new_Rlist[i-count]['prodI'][j] > ind:
                    new_Rlist[i-count]['prodI'][j] = new_Rlist[i-count]['prodI'][j] - 1
                            
    count = 0         
    lg_tmp = copy.deepcopy(new_lg)  
    for i in range(len(new_lg)):
        if new_lg[i]['index'] == ind:  # added 1/7
            del lg_tmp[i]
            count = count+1
            removed_lg.append(i)
            print("Gas Phase Species Detected!")
        elif new_lg[i]['index'] > ind:
            lg_tmp[i-count]['index'] = new_lg[i]['index'] - 1
    '''
    for i in range(len(new_lg)):
        if new_lg[i]['index'] == ind:  # added 1/7
            del new_lg[i]
            #count = count+1
            print("Gas Phase Species Detected!")
        elif new_lg[i]['index'] > ind:
            new_lg[i]['index'] = new_lg[i]['index'] - 1
    '''
                            
    new_lg = copy.deepcopy(lg_tmp)
    
    return [new_Rlist, new_lg, sort_ind, removed_reactions, removed_lg]