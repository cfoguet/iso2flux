from ..flux_functions.minimal_flux import create_minimal_fux_model,flux_minimization_fva
from variable_sampling import variable_sampling
from objfunc import objfunc
from get_variable_bounds import get_variable_bounds

def build_confidence_interval_dict(label_model,reaction_2_test,archi_x=None,lb_list=None,ub_list=None,dirc="max",best_flux_value=None,make_target_flux_independent=True,flux_penalty_dict={},max_flux=1e16,n_islands=6,pop_size=16):
    extra_constraints_dict={}
    variable_list_dict=label_model.variable_list_dict
    variable_list_dict[-1]={"type":"NULL"} 
    xtra_constraints_dict={}
    min_model,minimal_flux,non_milp_reactions =create_minimal_fux_model(label_model.constrained_model,fraction_of_optimum_objective=0.0,      fraction_of_flux_minimum=None,boundaries_precision=0.001,label_model=None,metabolite_list_file_name=None,flux_penalty_dict=flux_penalty_dict,maximum_flux=max_flux, mutually_exclusive_directionality_constraint=False)
    fva_reaction_2_test=flux_minimization_fva(min_model,solver=None,reaction_list=[reaction_2_test])
    #fva_reaction_2_test=flux_variability_analysis(label_model.constrained_model,reaction_list=[reaction_2_test]) #TODO change by a flux minimization FVA
    reaction_2_test_ub=fva_reaction_2_test[reaction_2_test]["maximum"]
    reaction_2_test_lb=fva_reaction_2_test[reaction_2_test]["minimum"]
    if best_flux_value!=None:
       if "max" in dirc:
           reaction_2_test_lb=max(reaction_2_test_lb,best_flux_value)
       elif "min" in dirc:
           reaction_2_test_ub=min(reaction_2_test_ub,best_flux_value)
    if lb_list==None or ub_list==None:
       extra_constraints_dict={reaction_2_test:{"lb":reaction_2_test_lb,"ub":reaction_2_test_ub}}
       lb_list,ub_list=get_variable_bounds(label_model,flux_penalty_dict,max_flux,extra_constraints_dict=extra_constraints_dict)
       
    if archi_x!=None:
       corrected_archi_x=copy.deepcopy(archi_x)
    else:
       corrected_archi_x=[]
       for n in range(n_islands):
           corrected_archi_x.append([])
           for n2 in range(pop_size):
               corrected_archi_x[n].append(variable_sampling(label_model,lb_list,ub_list,max_flux=max_flux,flux_penalty_dict=flux_penalty_dict,extra_constraints_dict=extra_constraints_dict))
    objfunc(label_model,corrected_archi_x[0][0]) #run it once to clear the result of any previous confidence interval work
    #print fva
    ###write_corrector 
    count_non_0=0
    for n,coef in enumerate(label_model.flux_solver_nullmnp[label_model.flux_solver_reaction_n_dict[reaction_2_test]]):
     if coef!=0:
       count_non_0+=1
       non_0_n=n
       non_0_coef=coef
       print [non_0_n,non_0_coef]
    if count_non_0>1 and make_target_flux_independent:
       lb_list[-1]=reaction_2_test_lb
       ub_list[-1]=reaction_2_test_ub
       n2modiffy=non_0_n
       coef2target=non_0_coef
       func="def correct_flux(nullmnp,free_fluxes_vector,desired_flux):\n new_flux=((desired_flux-np.dot(nullmnp[%s],free_fluxes_vector)+%s*free_fluxes_vector[%s])/(%s))\n free_fluxes_vector[%s]=new_flux" %(label_model.flux_solver_reaction_n_dict[reaction_2_test],coef2target,n2modiffy,coef2target,n2modiffy)
       exec(func)
       print func
       for n_pop,population in enumerate(corrected_archi_x):
           for n_ind, individual in enumerate(population):
               print len(individual)
               objfunc(label_model,individual,verbose=True)
               corrected_archi_x[n_pop][n_ind][-1]=label_model.reversible_flux_dict[reaction_2_test]
       variable_list_dict[-1]=({"type":"non_indepndent_flux","fun":correct_flux})
    target_flux_dict={"reaction":reaction_2_test,"dir":dirc,"lb":reaction_2_test_lb,"ub":reaction_2_test_ub,"factor":1}
    return target_flux_dict,corrected_archi_x, lb_list, ub_list


