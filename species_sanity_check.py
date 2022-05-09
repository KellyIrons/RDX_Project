"""
Species Sanity Check
"""

def species_removal(new_Rlist, ind, new_lg, sort_ind):
    
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
    for i in range(len(new_lg)):
        if new_lg[i-count]['index'] > ind:
            new_lg[i-count]['index'] = new_lg[i-count]['index'] - 1
                            
            
    count = 0
    sort_ind = sort_ind[0:-1]
    for i in range(len(sort_ind)):
        if sort_ind[i] > ind:
            sort_ind[i] = sort_ind[i] - 1
    
    return [new_Rlist, new_lg, sort_ind]