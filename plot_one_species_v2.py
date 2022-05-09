"""
plot_one_species_v2
"""

def plot_one_species_v2(species, TIME_o, solutionM_o, TIME_r, solutionM_r, TIME_s, solutionM_s):
    import matplotlib.pylab as plt
    
    for ind in range(len(species)):
        #old_ind = all_species_names.index(Slist[new_ind]['name'])
        plt.plot(TIME_o, solutionM_o[:,ind], 'b-', label = 'Original Mechanism')  
        plt.plot(TIME_s, solutionM_s[:,ind], 'r:', label = 'Skeletal Mechanism')
        plt.plot(TIME_r, solutionM_r[:,ind], 'g--', label = 'Reduced Mechanism')
        plt.xlabel('Time (s)')    
        plt.ylabel('Mass Fraction')
        plt.title('%s Results Comparison' % species[ind])
        plt.grid(True)
        plt.legend()
        plt.show()