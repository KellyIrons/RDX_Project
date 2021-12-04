"""
RDX_Graph_Creation
"""

def RDX_Graph_Creation(TIME, solutionM,t, Tset, Rlist, Slist, lg, EP, EP_paths):
    import numpy as np
    import networkx as nx
    from stoichcalc import stoichcalc
    from my_mech_obj import my_mech_obj
    from RDX_compute_DIC import compute_DIC
    
    '''
    # Test
    G = nx.Graph()
    G.add_edge(1, 2)  # default edge data=1
    G.add_edge(2, 3, weight=0.9)  # specify edge data
    '''

    [expmatrixF, coeffmatrixF, expmatrixB, coeffmatrixB] = stoichcalc(Rlist,Slist)


    reduction_type = 'species'
     # Find the DIC Coefficients at one time

    a = 1
    while TIME[a] < t and a <= len(TIME):
        a = a+1

    index = a
    time = TIME[int(a)]
    
    sample = solutionM[int(index),0:len(Slist)]
    
    mechobj = my_mech_obj(sample, time, Tset, Rlist, Slist, lg, expmatrixF, coeffmatrixF, expmatrixB, coeffmatrixB)
    
    DIC_spec, DIC_reac = compute_DIC(mechobj, sample, reduction_type)
    
    
    G = nx.MultiDiGraph()
    
    added = [False for i in range(len(Slist))]
    
    if len(Slist) == len(EP):
        EP_opt=True
    
    
   # for i in range(len(Slist)):
   #     G.add_nodes_from([(i, {"name": Slist[i]['name']})])
    
    for i in range(len(Slist)):
        for j in range(len(Slist)):
                if i!=j and DIC_spec[i,j] > 1E-1:
                    if not added[i]:
                        #G.add_nodes_from([(i, {"name": Slist[i]['name'], 'EP':EP[i], 'path':[EP_paths[i]]})])
                        G.add_nodes_from([(i, {"name": Slist[i]['name'], 'EP':EP[i]})])
                        added[i] = True
                    if not added[j]:
                        #G.add_nodes_from([(j, {"name": Slist[j]['name'],'EP':EP[j], 'path':[EP_paths[i]]})])
                        G.add_nodes_from([(j, {"name": Slist[j]['name'],'EP':EP[j]})])
                        added[j] = True
                    G.add_edge(i, j, weight=DIC_spec[i,j]) # arrow from i to j


    for i in range(len(Slist)):
        nx.add_path(G, EP_paths[i], path_name = Slist[i]['name'], path_num = i)

    
     
    nx.write_graphml(G, 'myRDXgraphskeletal_wpaths.graphml')
    return [G, DIC_spec]



