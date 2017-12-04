from ..flux_functions.minimal_flux import create_minimal_fux_model,flux_minimization_fva
from ..misc.round_functions import round_up,round_down

def get_variable_bounds(label_model,flux_penalty_dict={},max_flux=1e6,extra_constraints_dict={},min_model=None,precision=4):
   #TODO make fva based on the max flux solution
   if min_model==None:
      min_model,minimal_overall_flux,non_milp_reactions =create_minimal_fux_model(label_model.constrained_model,fraction_of_optimum_objective=0.0, fraction_of_flux_minimum=None,boundaries_precision=0.001,label_model=None,metabolite_list_file_name=None,flux_penalty_dict=flux_penalty_dict,maximum_flux=max_flux,mutually_exclusive_directionality_constraint=False,extra_constraints_dict=extra_constraints_dict)
   lb_list=[]
   ub_list=[]
   variable_list_dict=label_model.variable_list_dict
   #fva=add_flux_limit_constraints(copy.deepcopy(label_model.constrained_model),fraction_of_optimum_objective=0, fraction_of_flux_minimum=10)
   for n, x1 in enumerate(variable_list_dict):
        vtype=variable_list_dict[n]["type"]
        if vtype=="flux":
           nflux=variable_list_dict[n]["n"]
           reference_reaction=label_model.flux_solver_independent_flux_dict[nflux].keys()[0]
           coef=label_model.flux_solver_independent_flux_dict[nflux][reference_reaction]
           fva=flux_minimization_fva(min_model,solver=None,reaction_list=[reference_reaction])
           #fva=flux_variability_analysis(label_model.constrained_model,reaction_list=[reference_reaction], fraction_of_optimum=0,tolerance_feasibility=label_model.lp_tolerance_feasibility)
           lb=round_down(min(fva[reference_reaction]["minimum"]/coef,fva[reference_reaction]["maximum"]/coef),precision)#round_up(min(fva[reference_reaction]["minimum"]/coef,fva[reference_reaction]["maximum"]/coef),6)   
           ub=round_up(max(fva[reference_reaction]["minimum"]/coef,fva[reference_reaction]["maximum"]/coef),precision)#round_down(max(fva[reference_reaction]["minimum"]/coef,fva[reference_reaction]["maximum"]/coef),6)   
           print reference_reaction, lb,ub        
           lb_list.append(lb)
           ub_list.append(ub)          
        elif vtype=="turnover":
           reaction_id=variable_list_dict[n]["reaction"]
           fva=flux_minimization_fva(min_model,reaction_list=[reaction_id])
           #The turnover will be added twice to the overall flux (for the reverse and forward reaction)
           minimal_flux_value=min(abs(fva[reaction_id]["minimum"]),abs(fva[reaction_id]["maximum"]))
           #minimal_flux_value+max_ub*2<minimal_overall_flux
           max_theoretical_ub=(max_flux-minimal_flux_value)/2
           lb=float(label_model.turnover_flux_dict[reaction_id]["lb"])
           ub=float(label_model.turnover_flux_dict[reaction_id]["ub"])
           ub=min(ub,max_theoretical_ub)
           lb_list.append(lb)
           ub_list.append(ub) 
           print "turnover",lb,ub
        else: #Space saved for confidence intervals
           lb_list.append(0.0)
           ub_list.append(1.0)
   return lb_list, ub_list  
