"""
keep_reducing_RDX.py
"""

def keep_reducing_RDX(TIME_o,solutionM_o, TIME_f, solutionM_f, Slist_f, Rlist_f, lg_f, removed_species_f, removed_reactions_f, removed_lg_f, error_f, targets):
        
    from RDXfullparser import RDXfullparser
    from rateConstantCalc2 import rateConstantCalc2
    from stoichcalc import stoichcalc
    from RDX_sim2 import RDX_sim2
    from RDX_sim import RDX_sim
    from overall_EP import overall_EP
    from rank_EP import rank_EP
    from remove_and_remake import remove_and_remake
    from RDX_error_calc import RDX_error_calc
    from integrity_check import integrity_check
    from reshape_solutionM import reshape_solutionM
    from identify_secondary_targets import identify_secondary_targets 
    import matplotlib.pylab as plt
    import copy
    
    
    
    Tset = 538 # Kelvin - set temperature
    N = 10 # 
    R = 1 # number of species/reactions to remove at a time 
    #targets = ['RDX']
    s_tols = [7, 10, 15] # percent - absolute error tolerance
    r_tols = [7, 10, 15] 
    s_max_tols = [10, 15, 25] # percent - maximum error tolerance that will be temporarily accepted
    r_max_tols =  [10, 15, 25]
    initial_species = 'RDX'
    atol = 1e-12
    graph = False
    
    mech_type = "full"
    
    
    Slist = copy.deepcopy(Slist_f)
    Rlist =  copy.deepcopy(Rlist_f)
    lg = copy.deepcopy(lg_f)
    removed_species = removed_species_f.copy()
    removed_reactions = removed_reactions_f.copy()
    removed_lg = removed_lg_f.copy()
    J = len(Slist)
    
    [KF, KB] = rateConstantCalc2(Rlist, Slist, Tset)
    [expmatrixF, coeffmatrixF, expmatrixB, coeffmatrixB] = stoichcalc(Rlist,Slist)
    targets_f = targets.copy()
    targets_all = targets.copy()
    
    i = 0
    j = 0
    species_done = False
    reactions_done = False
    if mech_type == "full":
        EPmin  = 1E-35
    elif mech_type == "skeletal":
        EPmin = 1E-15
    
    ##### OVERALL LOOP ############################################################
    count = 2
    count_d = 0
    while not species_done or not reactions_done: # stop when BOTH conditions are true
    
        ns_o = len(Slist)
        nr_o = len(Rlist)
    
        if count <= (len(s_tols)-1):
            s_tol = s_tols[count]
            r_tol = r_tols[count]
            s_max_tol = s_max_tols[count]
            r_max_tol = r_max_tols[count]
        else:
            s_tol = s_tols[-1]
            r_tol = r_tols[-1]
            s_max_tol = s_max_tols[-1]
            r_max_tol = r_max_tols[-1]
        
    
        ### Species Loop ###
        
        ## Calculate Coefficients
        
        reduction_type = 'species'
        
        # Calculate EP
        print(targets)
        [EP, mechobj] = overall_EP(N, TIME_o, solutionM_o, Tset, Rlist, Slist, lg, expmatrixF, coeffmatrixF, expmatrixB, coeffmatrixB, targets, reduction_type, EPmin, graph)
        
        
        # Rank EP
        [ranked_EP, ranked_spec, sort_ind] = rank_EP(EP, mechobj, reduction_type, Rlist)
        
        sort_ind_o_r = sort_ind
        
        
        ## Remove Species
        ns = [len(Slist)]
        nr = [len(Rlist)]
        count_s = 0
        
        while not species_done: 
            
            # "Save" the old reduced arrays
            if count_s > 0:
                '''
                Slist = new_Slist.copy()
                Rlist = new_Rlist.copy()
                lg = new_lg.copy()
                '''
                Slist = copy.deepcopy(new_Slist)
                Rlist = copy.deepcopy(new_Rlist)
                lg = copy.deepcopy(new_lg)
                
            if len(Slist) > 150:
                R = 5
            else:
                R = 1
        
            # Remove a set number of species and remake the mechanism
            [new_Rlist, new_Slist, new_lg, removed_species, removed_reactions, removed_lg, sort_ind] = remove_and_remake(sort_ind, R, removed_species, removed_reactions, removed_lg, Rlist, Slist, lg, reduction_type)
            
            # Integrity Check
            [new_Slist, new_Rlist, new_lg, sort_ind, removed_species, removed_reactions, removed_lg] = integrity_check(new_Slist, new_Rlist, new_lg, sort_ind, removed_species, removed_reactions, removed_lg)
            
            
            
            
            # Recalculate contant matrices for new mechanism
            [KF, KB] = rateConstantCalc2(new_Rlist, new_Slist, Tset)
            [expmatrixF, coeffmatrixF, expmatrixB, coeffmatrixB] = stoichcalc(new_Rlist,new_Slist)
            
            
            # Re-run the simulation
            print("Running the simulation...thank you for your patience...")
            if len(new_Slist) >= 140: 
                [TIME_r, solutionM_r] = RDX_sim2(new_Rlist, new_Slist, new_lg, Tset, KF, KB, expmatrixF, coeffmatrixF, expmatrixB, coeffmatrixB, mech_type, initial_species, atol)
                solutionM_r = reshape_solutionM(solutionM_r)
            else:
                [TIME_r, solutionM_r] = RDX_sim(new_Rlist, new_Slist, new_lg, Tset, KF, KB, expmatrixF, coeffmatrixF, expmatrixB, coeffmatrixB, mech_type, initial_species, atol)    
            
            # Calculae the error
            #[error, results, ts] = RDX_error_calc(TIME_o, solutionM_o, solutionM_r, TIME_r, N, removed_species, removed_lg, J)
            failed_run = TIME_r[-1] < 0.75*TIME_o[-1]
            if failed_run:
                atol = 1e-17
                print("Running the simulation...thank you for your patience...")
                [TIME_r, solutionM_r] = RDX_sim2(new_Rlist, new_Slist, new_lg, Tset, KF, KB, expmatrixF, coeffmatrixF, expmatrixB, coeffmatrixB, mech_type, initial_species, atol)
                solutionM_r = reshape_solutionM(solutionM_r)  
                atol = 1e-15
                failed_run = TIME_r[-1] < 0.75*TIME_o[-1]
                if failed_run:
                    error = 1000
                else:
                    # Calculae the error
                    [error, results, ts] = RDX_error_calc(TIME_o, solutionM_o, solutionM_r, TIME_r, N, removed_species, removed_lg, J)
            else:
                # Calculae the error
                [error, results, ts] = RDX_error_calc(TIME_o, solutionM_o, solutionM_r, TIME_r, N, removed_species, removed_lg, J)       
        
            if error <= s_tol:
                TIME_f = TIME_r.copy()
                solutionM_f = solutionM_r.copy()
                '''
                Slist_f = new_Slist.copy()
                Rlist_f = new_Rlist.copy()
                lg_f = new_lg.copy()
                '''
                Slist_f = copy.deepcopy(new_Slist)
                Rlist_f = copy.deepcopy(new_Rlist)
                lg_f = copy.deepcopy(new_lg)
                removed_species_f = removed_species.copy()
                removed_reactions_f = removed_reactions.copy()
                removed_lg_f = removed_lg.copy()
                error_f = error
                targets_f = targets.copy()
                print("Saved Mechanism!")
            
            plt.plot(ts[0,:], results[0,:])
            plt.grid()
            plt.xlabel("Time [s]")
            plt.ylabel("Median % Error")
            plt.title("Reduced Mechanism Error: %i Species & %i Reactions" % (len(new_Slist), len (new_Rlist)))
            plt.show()
            
            
            plt.plot(TIME_r, solutionM_r)
            plt.grid()
            plt.ylim((0, 0.15))
            plt.xlabel("Time [s]")
            plt.ylabel("Mass Fraction")
            plt.title("Reduced Mechanism: %i Species & %i Reactions" % (len(new_Slist), len (new_Rlist)))
            plt.show()
        
            i = i+1
            ns.append(len(new_Slist))
            nr.append(len(new_Rlist))
            
            print("Error = %f %%" % error)
            
            count_s = count_s + 1
        
            #if ns[i] < ns[i-1]: #leads to oo loop :(
            if error <= s_max_tol:
                species_done = False
                reactions_done = False
            else:
                '''
                Slist = Slist_f.copy()
                Rlist = Rlist_f.copy()
                lg = lg_f.copy()
                '''
                Slist = copy.deepcopy(Slist_f)
                Rlist =  copy.deepcopy(Rlist_f)
                lg = copy.deepcopy(lg_f)
                removed_species = removed_species_f.copy()
                removed_reactions = removed_reactions_f.copy()
                removed_lg = removed_lg_f.copy()
                error = error_f
                targets = targets_f.copy()
                species_done = True
            
                
                print("Exiting Species Loop")
     
        
        ### Reaction Loop ###
        
        
        reduction_type = 'reactions'
        
        # Recalculate coefficent matrices using species reduced mechanism
        [KF, KB] = rateConstantCalc2(Rlist, Slist, Tset)
        [expmatrixF, coeffmatrixF, expmatrixB, coeffmatrixB] = stoichcalc(Rlist,Slist)
        
        # Calculate EP
        [EP, mechobj] = overall_EP(N, TIME_o, solutionM_o, Tset, Rlist, Slist, lg, expmatrixF, coeffmatrixF, expmatrixB, coeffmatrixB, targets, reduction_type, EPmin, graph)
        
        
        # Rank EP
        [ranked_EP, ranked_reac, sort_ind] = rank_EP(EP, mechobj, reduction_type, Rlist)
        
        
        ## Remove Reactions
        ns = [len(Slist)]
        nr = [len(Rlist)]
        count_r = 0
    
        while not reactions_done:# and error <= max_tol:
            
    
            # "Save" the old reduced arrays
            if count_r > 0:
                '''
                Slist = new_Slist
                Rlist = new_Rlist
                lg = new_lg
                '''
                Slist = copy.deepcopy(new_Slist)
                Rlist = copy.deepcopy(new_Rlist)
                lg = copy.deepcopy(new_lg)
            
            if len(Slist) > 150:
                R = 5
            else:
                R = 1
            
        
            # Remove a set number of species and remake the mechanism
            [new_Rlist, new_Slist, new_lg, removed_species, removed_reactions, removed_lg, sort_ind] = remove_and_remake(sort_ind, R, removed_species, removed_reactions, removed_lg, Rlist, Slist, lg, reduction_type)
            
            
            # Recalculate contant matrices for new mechanism
            [KF, KB] = rateConstantCalc2(new_Rlist, new_Slist, Tset)
            [expmatrixF, coeffmatrixF, expmatrixB, coeffmatrixB] = stoichcalc(new_Rlist,new_Slist)
            
            
            # Re-run the simulation
            print("Running the simulation...thank you for your patience...")
            #[TIME_r, solutionM_r] = RDX_sim2(new_Rlist, new_Slist, new_lg, Tset, KF, KB, expmatrixF, coeffmatrixF, expmatrixB, coeffmatrixB, mech_type, initial_species, atol)
            #solutionM_r = reshape_solutionM(solutionM_r)
            
            if len(new_Slist) >= 140: 
                [TIME_r, solutionM_r] = RDX_sim2(new_Rlist, new_Slist, new_lg, Tset, KF, KB, expmatrixF, coeffmatrixF, expmatrixB, coeffmatrixB, mech_type, initial_species, atol)
                solutionM_r = reshape_solutionM(solutionM_r)
            else:
                [TIME_r, solutionM_r] = RDX_sim(new_Rlist, new_Slist, new_lg, Tset, KF, KB, expmatrixF, coeffmatrixF, expmatrixB, coeffmatrixB, mech_type, initial_species, atol)    
            
            
            failed_run = TIME_r[-1] < 0.75*TIME_o[-1]
            if failed_run:
                atol = 1e-17
                print("Running the simulation...thank you for your patience...")
                [TIME_r, solutionM_r] = RDX_sim2(new_Rlist, new_Slist, new_lg, Tset, KF, KB, expmatrixF, coeffmatrixF, expmatrixB, coeffmatrixB, mech_type, initial_species, atol)
                solutionM_r = reshape_solutionM(solutionM_r)  
                atol = 1e-15
                failed_run = TIME_r[-1] < 0.75*TIME_o[-1]
                if failed_run:
                    error = 1000
                else:
                    # Calculae the error
                    [error, results, ts] = RDX_error_calc(TIME_o, solutionM_o, solutionM_r, TIME_r, N, removed_species, removed_lg, J)
            else:
                # Calculae the error
                [error, results, ts] = RDX_error_calc(TIME_o, solutionM_o, solutionM_r, TIME_r, N, removed_species, removed_lg, J)
            
            if error <= r_tol:
                TIME_f = TIME_r
                solutionM_f = solutionM_r
                '''
                Slist_f = new_Slist
                Rlist_f = new_Rlist
                lg_f = new_lg
                '''
                Slist_f = copy.deepcopy(new_Slist)
                Rlist_f = copy.deepcopy(new_Rlist)
                lg_f = copy.deepcopy(new_lg)
                removed_species_f = removed_species.copy()
                removed_reactions_f = removed_reactions.copy()
                removed_lg_f = removed_lg.copy()
                error_f = error
                targets_f = targets.copy()
                print("Saved Mechanism!")
            
            
            plt.plot(ts[0,:], results[0,:])
            plt.grid()
            plt.xlabel("Time [s]")
            plt.ylabel("Median %% Error")
            plt.title("Reduced Mechanism Error: %i Species & %i Reactions" % (len(new_Slist), len (new_Rlist)))
            plt.show()
            
            
            plt.plot(TIME_r, solutionM_r)
            plt.grid()
            plt.ylim((0, 0.15))
            plt.xlabel("Time [s]")
            plt.ylabel("Mass Fraction")
            plt.title("Reduced Mechanism: %i Species & %i Reactions" % (len(new_Slist), len (new_Rlist)))
            plt.show()
        
            j = j+1
            ns.append(len(new_Slist))
            nr.append(len(new_Rlist))
        
        
            print("Error = %f %%" % error)
            
            count_r = count_r + 1
            
            if error <= r_max_tol:
                species_done = False
                reactions_done = False
            else:
                reactions_done = True
                '''
                Slist = Slist_f
                Rlist = Rlist_f
                lg = lg_f
                '''
                Slist = copy.deepcopy(Slist_f)
                Rlist =  copy.deepcopy(Rlist_f)
                lg = copy.deepcopy(lg_f)
                removed_species = removed_species_f.copy()
                removed_reactions = removed_reactions_f.copy()
                removed_lg = removed_lg_f.copy()
                error = error_f
                targets = targets_f.copy()
                [expmatrixF, coeffmatrixF, expmatrixB, coeffmatrixB] = stoichcalc(Rlist,Slist) #need to do this so species reduction works
                
                print("Exiting Reation Loop")
    
        count = count+1
        
        print(species_done)
        print(reactions_done)
        print(count_d)
        print(targets)
        
        
        if species_done and reactions_done:
            if count_d <= 3:
                # Further reduction ?
                t = 1.75
                reduction_type = 'species'
                [expmatrixF, coeffmatrixF, expmatrixB, coeffmatrixB] = stoichcalc(Rlist_f,Slist_f)
                [EP, mechobj] = overall_EP(N, TIME_f, solutionM_f, Tset, Rlist_f, Slist_f, lg_f, expmatrixF, coeffmatrixF, expmatrixB, coeffmatrixB, targets, reduction_type, EPmin, graph)
                [g_mean, sort_ind_2] = identify_secondary_targets(TIME_f, solutionM_f, targets, t, Tset, Rlist_f, Slist_f, lg_f, expmatrixF, coeffmatrixF, expmatrixB, coeffmatrixB, EP)  
                # Check that this soecies hasn't already been tried
                found = Slist_f[sort_ind_2[0]]['name'] not in targets_all
                ind = 0
                
                while not found:
                    print('Kelly you need to switch things up!')
                    ind = ind+1
                    found = Slist_f[sort_ind_2[ind]]['name'] not in targets_all
                
                targets.append(Slist_f[sort_ind_2[ind]]['name'])
                targets_all.append(Slist_f[sort_ind_2[ind]]['name'])
                print('ADDED SPECIES: ' + Slist_f[sort_ind_2[ind]]['name'])
                count_d = count_d + 1
                species_done = False
                reactions_done = False
        
        
        
        
        
    
    # Finishing Steps
    
    
    print('%i Species Remaining!' % len(Slist))
    
    print('%i Reactions Remaining!' % len(Rlist))
    
    
    plt.plot(TIME_f, solutionM_f[:,0:len(Slist)])
    plt.grid()
    plt.ylim((0, 0.15))
    plt.xlabel("Time [s]")
    plt.ylabel("Mass Fraction")
    plt.title("Final Mechanism")
    plt.show()


    return [TIME_f, solutionM_f, Slist_f, Rlist_f, lg_f, removed_species_f, removed_reactions_f, removed_lg_f, error_f, targets]
