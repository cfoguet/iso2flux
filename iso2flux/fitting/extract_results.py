import numpy as np
def extract_results(archi,optimal_variables_fname=None,archi_xlist_fname=None,champion_replace=True):
    optimal_solution=min(([isl.population.champion.f for isl in archi]))
    print optimal_solution
    n_islands=len(archi)
    n_pop=len(archi[0])
    n_parameters=len(archi[0].population.champion.x)
    champion_x_list=np.zeros((n_islands,n_parameters))
    archi_xlist=np.zeros(shape=(n_islands,n_pop,n_parameters))
    v_list=np.zeros(shape=(n_islands,n_pop,n_parameters))
    for n1,population in enumerate([isl.population for isl in archi]):
        for n2,ind in enumerate(population):
            for n3,value  in enumerate(ind.cur_x):
                archi_xlist[n1][n2][n3]=value
            for n3,velocity  in enumerate(ind.cur_v):
                v_list[n1][n2][n3]=velocity
            if champion_replace and n2 in population.get_best_idx(1):
               archi_xlist[n1][n2]==np.asarray(population.champion.x)
        if population.champion.f==optimal_solution:
           optimal_variables=np.asarray(population.champion.x)
        champion_x_list[n1]=np.asarray(population.champion.x)
    if optimal_variables_fname!=None:
       np.savetxt(optimal_variables_fname,optimal_variables)
    if archi_xlist_fname!=None:
       np.savetxt(archi_xlist_fname,archi_xlist)
    return  optimal_solution[0],optimal_variables, archi_xlist, champion_x_list, v_list
   


"""
def extract_results(archi):
    archi_xlist=[]
    optimal_solution=min(([isl.population.champion.f for isl in archi]))
    print optimal_solution
    for population in ([isl.population for isl in archi]):
        x_list=[]
        for ind in population:
            x_list.append(list(ind.best_x))
        archi_xlist.append(x_list)
        if population.champion.f==optimal_solution:
           optimal_variables=population.champion.x
    return  optimal_solution[0],list(optimal_variables), archi_xlist


a,b,c,d=extract_results(archi)
   
"""
