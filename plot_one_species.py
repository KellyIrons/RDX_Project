"""
Plot One Species
"""

def plot_one_species(species, Slist, TIME_o, solutionM_o, TIME_f, solutionM_f, mechobj):
    import matplotlib.pylab as plt
        
    ''' 
    species_names = mechobj['species_names']
    for new_ind in species:
        old_ind = species_names.index(Slist[new_ind]['name'])
        plt.plot(TIME_o, solutionM_o[:,old_ind], 'b-', label = 'Original Mechanism')  
        plt.plot(TIME_f, solutionM_f[:,new_ind], 'r:', label = 'Reduced Mechanism')
        plt.xlabel('Time (s)')    
        plt.ylabel('Mass Fraction')
        plt.title('%s Results Comparison' % Slist[new_ind]['name'])
        plt.grid(True)
        plt.legend()
        plt.show()
    '''
    
    from RDXskeletalparser import RDXskeletalparser
    [old_Rlist, old_Slist, old_lg] = RDXskeletalparser()
    
    all_species_names = [x['name'] for x in old_Slist]
    
    

    for new_ind in species:
        old_ind = all_species_names.index(Slist[new_ind]['name'])
        plt.plot(TIME_o, solutionM_o[:,old_ind], 'b-', label = 'Original Mechanism')  
        plt.plot(TIME_f, solutionM_f[:,new_ind], 'r:', label = 'Reduced Mechanism')
        plt.xlabel('Time (s)')    
        plt.ylabel('Mass Fraction')
        plt.title('%s Results Comparison' % Slist[new_ind]['name'])
        plt.grid(True)
        plt.legend()
        plt.show()
    