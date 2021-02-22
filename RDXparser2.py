'''
RDX MECHANISM PARSER v2
'''
def RDXparser2():

    ### Create a string of species names ###
    spec = " " #initialize a string of species names
    List = [];
    with open("RDXspecies.txt") as fp: #open the species file
        line = fp.readline() #read the first line
        while line: #while there is a line to read
            text = line.split(" ") #split the line at each space to get a list (text) of strings
            spec = spec + " " + text[0] #add the species name to the species string
            string = text[0]
            List.append(string)
            line = fp.readline() #read the next line
    
    ### Create a list of dictionaries of species data ###
    Slist = [] #initialize a list of dicts
    linenumber = 0 #initialize line counter
    with open("RDXthermo.txt") as fp: #open the thermo file
        line = fp.readline() #read the first line
        while line: #while there is a line to read
            linenumber = linenumber + 1 #add one to the line counter
            if (linenumber+3)%4 == 0: #if the current line is a multiple of four from the beginning
                text = line.split(" ") #split the line at each space to be able to access its components
                name = text[0] #the name of the species is the fist
                carbon = int(text[2].rstrip('H'))
                hydrogen = int(text[3].rstrip('N'))
                nitrogen = int(text[4].rstrip('O'))
                oxygen = int(text[5].rstrip('L'))
                MW = (carbon*12.011 + hydrogen*1.008 +nitrogen*14.007 + oxygen*15.999)/1000; #molec weight - kg/mol
                rad = (0.617*(1000*MW)**(1/3)) * 1E-10 #molecular radius - m
                ThermoD = {'name':name, 'MW':MW, 'rad':rad}
                Slist.append(ThermoD)
            line = fp.readline() #read the next line
    
    
    ### Create a list of dictionaries of reaction data ###
    Rlist = [] #initialize list of dicts
    linenumber = 0 #initialize line (reaction number) counter
    with open("RDXreactions.txt") as fp: #open the reaction file
        line = fp.readline() #read the first line
        while line: #while there is a line to read
            linenumber = linenumber + 1 #add one to the line counter
            name = linenumber-1 #id = string containing the reaction number
            text = line.split(" ") #split the line at each space to be able to access its components
            rxn = text[0].split("=") #the reaction is the first part of the text, split it at the = sign
            lh = rxn[0] #lh = string containing the characters left of the = (reactants)
            rh = rxn[1] #rh = string containing the characters right of the = (products)
            lhs = lh.split("+") #lhs = list of each species left of the =  (reactants)
            rhs = rh.split("+") #rhs = list of each species right of the = (products)
            rcoeffs = []
            pcoeffs = []
            reactants = []
            products = []
            reacI = []
            prodI = []
            for species in lhs: #for each reactant...
                LL = len(lhs) #LL = number of reactants
                count = 0 #initialize counter
                for x in range(LL): #for the number of reactants
                    if species == lhs[x]: #if a species is found in the reactants...
                        count = count + 1 #... increase the count of times it appears
                if count > 1: #if the species appears more than once
                        for x in range(count): #for the number of times the species appears...
                            lhs.remove(species) #...remove the species from the reactants
                        lhs.append("%d %s" %(count, species)) #add the species back in, with correct coefficient
            for species in rhs: #for each product ...
                LR = len(rhs) #LR = number of products
                count = 0 #initialize counter
                for x in range(LR): #for the number of products
                    if species == rhs[x]: #if a species is found in the reactants...
                        count = count + 1 #... increase the count of times it appears
                if count > 1: #if the species appears more than once
                    for x in range(count): #for the number of times the species appears...
                        rhs.remove(species) #...remove that species from the products
                    rhs.append("%d %s" %(count, species)) #add the species back in, with correct coefficient
            for species in lhs:
                species2 = species.split(' ')
                if len(species2) == 1:
                    rcoeffs.append(1)
                    reactants.append(species)
                elif len(species2) == 2:
                    rcoeffs.append(int(species2[0]))
                    reactants.append(species2[1])
            for species in rhs:
                species2 = species.split(' ')
                if len(species2) == 1:
                    pcoeffs.append(1)
                    products.append(species)
                elif len(species2) == 2:
                    pcoeffs.append(int(species2[0]))
                    products.append(species2[1])  
            for species in reactants:
                ind = List.index(species)
                reacI.append(ind)
            for species in products:
                ind = List.index(species)
                prodI.append(ind)
            LHS = " + ".join(lhs) #join reactants together with correct format
            RHS = " + ".join(rhs) #join products together with correct format
            ctirxn = " <=> ".join([LHS, RHS]) #join the reactants and products with correct formatting
            ImagFreq = (3E8)*(100*float(text[1])) #define ImagFreq
            #dHf = float(text[2]) #define dHf
            #dHb = float(text[3]) #define dHb
            dGf = 4.184*float(text[4]) #define dGf
            dGb = 4.184*float(text[5]) #define dGb
            ReactionD = {'id':name, 'w':ImagFreq, 'dGf':dGf, 'dGb':dGb, 'eq':ctirxn,
                         'reactants':reactants, 'rcoeffs':rcoeffs, 'products':products, 
                         'pcoeffs':pcoeffs, 'reacI':reacI, 'prodI':prodI} #ReactionD = dict with reaction info
            Rlist.append(ReactionD) #add this reaction's dict onto the list
            line = fp.readline() #read the next line
    
    index = [299, 293, 298, 20, 16, 18, 17, 19, 1, 320, 319, 291]
    A = [4.226E13, 2.7E2, 6.133E11, 2.7E2, 3.9E5, 2.7E3, 1E3, 3.9E5, 3.9E5, 3.9E5, 3.9E5, 0]
    n = [1,0,1,0,0,0,0,0,0,0,0,0]
    E = [3.755E4,6.2E3, 3.393E4,8.644E3,6.904E3,1.4E3,6.934E3,1.3856E4,1.35E4,1.4E4,2.467E4,0]
    
    lg = []
    for x in range(len(index)):
        index[x] = index[x] - 1
        E[x] = 4.184*E[x]
        D = { 'index':index[x], 'A':A[x], 'n':n[x], 'E':E[x]}
        lg.append(D)
        
    return [Rlist, Slist, lg]