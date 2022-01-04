"""
Species Sanity Check
"""

def species_removal(new_Rlist, ind, new_lg, sort_ind, reduction_type):
    
    removed_Rs = []
    count =  0
    for i in range(len(new_Rlist)):
        indices = new_Rlist[i-count]['reacI'] + new_Rlist[i-count]['prodI']
        if ind in indices:
            print('Removed Reaction: %s' % new_Rlist[i-count]['eq'])
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
    for i in range(len(new_lg)):
        if new_lg[i-count]['index'] > ind:
            new_lg[i-count]['index'] = new_lg[i-count]['index'] - 1
                            
    
    
    return [new_Rlist, new_lg, sort_ind]