'''
RDX SKELETAL MECHANISM PARSER v2
'''
def RDXskeletalparser():

    #Input Skeletal Reactions   #Currently the data type is a list
    skeletalReactions = ['RDX <=> INT175a + t_HONO','RDX + NO <=> ONDNTA + NO2','RDX + NO2 <=> INT221a + c_HONO',
    'RDX + CH2N <=> INT221a + CH2NH','INT221a <=> INT175a + NO2','INT175a <=> INT128a + t_HONO','INT128a <=> TAZ + t_HONO',
    'INT175a + NO2 <=> INT174a + c_HONO','INT175a + c_HONO <=> INT222a','INT222a <=> INT222a_RO','INT222a_RO <=> INT132a + N2O + NO2',
    'INT132a <=> INT86a + NO2','INT86a + ONNO2 <=> INT162b','INT162b <=> INT132d + CH2O','INT132d <=> CH2NHCHO + NO2 + N2',
    'CH2NHCHO + NO2 <=> CH2NH + c_HONO + CO','CH2NHCHO + NO2 <=> CH2NCHO + t_HONO','CH2NCHO + NO2 <=> CH2NCO + c_HONO',
    'CH2NCO + NO2 <=> INT102h','INT102h <=> CH2NNO + CO2','CH2NNO <=> CH2N + NO','INT175a + ONNO2 <=> INT251a',
    'INT251a <=> INT251a_RO','INT251a_RO <=> INT161a + N2O + NO2','INT161a <=> INT115a + NO2','INT175a <=> INT101a + CH2NNO2',
    'CH2NNO2 <=> CH2N + NO2','CH2NNO2 <=> HCN + t_HONO','CH2NNO2 + c_HONO <=> INT121a','INT121a <=> INT121b','INT121b <=> INT90a + HNO',
    'CH2NNO2 + ONNO2 <=> INT150a','INT150a <=> CH2O + INT120b','INT120b <=> N2O + ONONO','2 t_HONO <=> ONNO2 + H2O','2 NO <=> ONNO',
    '2 NO2 <=> ONONO2','2 NO2 <=> O2NNO2','NO + NO2 <=> ONNO2','NO + NO2 <=> ONONO','t_HONO <=> c_HONO','CH2NH + NO2 <=> CH2N + c_HONO',
    'CH2NH + NO2 <=> CH2N + t_HONO','CH2N + NO2 <=> HCN + c_HONO','HNO2 + CH2N <=> NO2 + CH2NH','t_HONO + c_HONO <=> ONNO2 + H2O',
    '2 c_HONO <=> ONONO + H2O','ONONO <=> ONNO2','CH2O + NO2 <=> c_HONO + HCO','HCOOH + NO2 <=> HOCO + c_HONO','HNO + NO2 <=> NO + c_HONO',
    'CH2N + NO <=> HCN + HNO','HOCO + NO <=> CO + t_HONO','HCO + NO2 <=> INT75c','INT75c <=> CO + t_HONO','HCO + NO <=> CO + HNO']


    #nSkeletal = len(skeletalReactions) #number of skeletal reactions 

    #Parse Mechanism
    from RDXfullparser import RDXfullparser
    [Rlist, Slist, lg] = RDXfullparser() 
    
    nFull = len(Rlist) #total number of reactions
    
    #Build the skeletal reaction struct from the full reaction struct
    S_Rlist = []
    for reaction in skeletalReactions: #each reaction is a string
        match = False
        count = 0
        #print(type(reaction))
        while count <= nFull and match == False:
            match = reaction == Rlist[count]['eq']
            #match = strcmp(string, p(count).eq);
            #print((Rlist[count]['eq']))
            count = count + 1;
        if match == True:
            #print(count)
            S_Rlist.append(Rlist[count-1])
        else:
            print("Reaction not found :(")
            
    #Build the skeletal species struct from the full species struct
    #for each reaction go through each reactant find the species in the 
    #species struct, if it hasn't already been found add it to the species 
    #struct
    
    indices = []
    for reaction in S_Rlist:
        #identtify all the involved species, add to a long list of species indices
        indices = indices + reaction['reacI'] + reaction['prodI']
    #remove the duplicates to get a list of unique species
    indices = list(dict.fromkeys(indices))
  
   
    S_Slist = []
    S_Slist = list(S_Slist)
 
    for i in indices:
        S_Slist.append(Slist[i])

    
    #Now we need to go through the reactions and update the species indices
    for reaction in S_Rlist:       
        newR = []
        newP = []
        for species in reaction['reacI']:
            #ind = S_Slist['name'].index(species)
            ind = indices.index(species)
            newR.append(ind)
        for species in reaction['prodI']:
            #ind = S_Slist.index(species)
            ind = indices.index(species)
            newP.append(ind)
        reaction['reacI'] = []
        reaction['prodI'] = []       
        reaction['reacI'] = newR
        reaction['prodI'] = newP
    
    #Finally, need to update the indices in lg
    for species in lg:
        oldS = species['index']
        ind = indices.index(oldS)
        species['index'] = ind
        species['MW'] = S_Slist[ind]['MW']
    S_lg = lg
    

    return [S_Rlist, S_Slist, S_lg]


