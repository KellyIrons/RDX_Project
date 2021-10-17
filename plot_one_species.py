"""
Plot One Species
"""

def plot_one_species(species, Slist, TIME_o, solutionM_o, TIME_f, solutionM_f):
    import matplotlib.pylab as plt
        
    for spec in species:
        plt.plot(TIME_o, solutionM_o[:,spec], 'b-', label = 'Original Mechanism')  
        plt.plot(TIME_f, solutionM_f[:,spec], 'r:', label = 'Reduced Mechanism')
        #plt.plot(TIME_r, solutionM_r[:,spec], 'b--')
        #plt.legend('Original Mechanism', 'Reduced Mechanism')
        plt.xlabel('Time (s)')    
        plt.ylabel('Mass Fraction')
        plt.title('%s Results Comparison' % Slist[spec]['name'])
        plt.grid(True)
        plt.legend()
        plt.show()