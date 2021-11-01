"""
RDX Modified compute_DIC from ARCANE
"""
def compute_DIC(mechobj, mysample, reduction_type):
    """Compute Direct interaction coefficients (DIC) for a given sample

    :param mechobj: Mechanism object
    :param mysample: sample for which DIC is computed
    :param reduction_type: type of reduction performed (species or reactions)
    :return: number_of_species * number_of_species matrix of coefficients for type species
             number_of_species * number_of_reactions matrix of coefficients for type reactions

    Created: 17/11/14 [PP]
    Last modified: 18/04/03 [QC]
    """
    import numpy as np

    # Parameters
    myns = mechobj['ns']
    mynr = mechobj['nr']
    #net = mechobj.network NA
    #ctmech = mechobj.ctmech
    nup = mechobj['nup']
    nur = mechobj['nur']
    nup = np.transpose(nup)
    nur = np.transpose(nur)

    # Get reaction rates for each reaction
    #ctmech.TPY = float(mysample.T), float(mysample.P), [float(myY) for myY in mysample.Y]

    #omega = ctmech.net_rates_of_progress
    omega = mechobj['omega']
    

    # Get production and consumption rates for species
    #PA = ctmech.creation_rates
    #CA = ctmech.destruction_rates
    PA = mechobj['PA']
    CA = mechobj['CA']

    # Get enthalpy from production and consumption of species in reactions
    '''
    HR_species = species_heat_release(mechobj, mysample)
    HR_reactions = reactions_heat_release(mechobj, mysample)
    HR_prod_spec = sum([HR_spec for HR_spec in HR_species if HR_spec > 0])
    HR_cons_spec = sum([HR_spec for HR_spec in HR_species if HR_spec < 0])
    HR_prod_reac = sum([HR_reac for HR_reac in HR_reactions if HR_reac > 0])
    HR_cons_reac = sum([HR_reac for HR_reac in HR_reactions if HR_reac < 0])
    i_HR = myns
    '''


    # Workspace
    DIC_spec = np.zeros([myns + 1, myns], 'd')

    # Evaluate DIC(i,j)
    S_ind = mechobj['S_ind']
    for i in range(myns): #for each species
        # reactions containing species i
        #booli = net.indr[net.inds == i]
        #indj = net.indsj[net.indsi == i]
        #booli = S_ind[i]['indbool']
        #indj = S_ind[i]['ind'] #species indices
        booli = S_ind[i]['ind'] #reactions with species i
        indj = S_ind[i]['other_spec'] #species involved in those reactions

        for j in indj:
            #boolj = net.indr[net.inds == j]  # reactions containing species j
            #boolj = S_ind[j]['indbool']
            boolj = S_ind[j]['ind'] #reactions with species j
            
            indk = np.intersect1d(booli, boolj)  # reactions containing species i and j

            # Compute the DIC
            for k in indk:
                DIC_spec[i, j] += (nup[i, k] - nur[i, k]) * omega[0,k]
                #Go through each pair of species, go through each reaction containing
                #both of those species and compute the DIC
                
            #Notes
            # - indj are species indices --> take to be the species on the other
            #   side of reactions containing species i 
            # - indk are reaction indices
            # - booli and boolj must be intergers 
                

            
            # Normalize

            DIC_spec[i, j] = abs(DIC_spec[i, j]) / max(max(PA[0,i], CA[0,i]), 1e-60)
            #DIC_spec[i, j] = abs(DIC_spec[i, j]) / val

        # Heat release term
        #DIC_spec[i_HR, i] += abs(HR_species[i])

        # Normalize
        #DIC_spec[i_HR, i] = abs(DIC_spec[i_HR, i]) / max(max(abs(HR_prod_spec), abs(HR_cons_spec)), 1e-60)

    if reduction_type in ['reactions', 'R']:
        # Workspace
        DIC_reac = np.zeros([myns + 1, mynr], 'd')

        # Evaluate DIC(i,j)
        for i in range(myns):
            # reactions containing species i
            #indk = net.indr[net.inds == i]  # reactions containing species i
            indk = S_ind[i]['indbool']

            # Compute the DIC
            for k in indk:
                DIC_reac[i, k] = abs(nup[i, k] - nur[i, k]) * omega[k]

                # Normalize
                DIC_reac[i, k] = abs(DIC_reac[i, k]) / max(max(PA[0,i], CA[0,i]), 1e-60)

        #for i in range(mynr):

            # Heat release term
            #DIC_reac[i_HR, i] += abs(HR_reactions[i])

            # Normalize
            #DIC_reac[i_HR, i] = abs(DIC_reac[i_HR, i]) / max(max(abs(HR_prod_reac), abs(HR_cons_reac)), 1e-60)

    else:
        DIC_reac = None

    return DIC_spec, DIC_reac
