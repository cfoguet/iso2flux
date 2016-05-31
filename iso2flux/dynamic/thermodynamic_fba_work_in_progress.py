from cplex import Cplex,SparsePair
from cobra.solvers.cplex_solver import create_problem

def write_thermofba_problem(label_model,non_labelled_model,temperature=310,objective=[]):
    model=non_labelled_model
    old_objective=[]
    for x in model.reactions:
           if x.objective_coefficient!=0:
              old_objective.append([x,x.objective_coefficient])
              x.objective_coefficient=0
              
    lp=create_problem(model, quadratic_component=None)
    n_variables=len(lp.variables.get_names())       
    #We can restore the old objective to the model now. 
    for x in old_objective:
        x[0].objective_coefficient=x[1]
    
    alpha=100000000
    RT=8.3114*temperature
    
    variable_names=[]
    lower_bounds=[]
    upper_bounds=[]
    objective_coefficients=[]
    variable_kinds=[]
    
    
    constraint_sense = []
    constraint_names = []
    constraint_limits = []
    the_linear_expressions = []
    
    for metabolite in model.metabolites: 
        if metabolite.id not in ref_con_dict:
           continue  
        
        variable_names.append(metabolite.id+"_log_con")
        if metabolite.id in concentration_range_dict:
           lb=math.log(concentration_range_dict[metabolite.id][0]/1000)
           ub=math.log(concentration_range_dict[metabolite.id][1]/1000)
           lower_bounds.append(lb) #Aproximatelly LN of 0.001
           upper_bounds.append(ub)
        else:
           lower_bounds.append(-16.2) #Aproximatelly LN of 0.001
           upper_bounds.append(-3.68) #Aproximatelly LN of 25
        if metabolite.id in objective:
              objective_coefficients.append(1)
        else:
           objective_coefficients.append(0.0)
        variable_kinds.append(Cplex.variables.type.continuous)  
    
    for reaction in model.reactions:
        if reaction.id in deltaGibs_dict:
           variable_names.append("Gibs_"+reaction.id)
           lower_bounds.append(-100000)
           upper_bounds.append(100000)
           objective_coefficients.append(0.0)
           variable_kinds.append(Cplex.variables.type.continuous)
           
           
           variable_names.append("d_"+reaction.id)  
           lower_bounds.append(0)
           upper_bounds.append(1)
           objective_coefficients.append(0.0)
           variable_kinds.append(Cplex.variables.type.binary)
           
           
           constraint_sense.append("G")
           constraint_names.append(reaction.id+"_d1.1")
           constraint_limits.append(0.0)
           the_linear_expressions.append(SparsePair(ind=[reaction.id,"d_"+reaction.id],val=[1,alpha]))
           
           constraint_sense.append("L")
           constraint_names.append(reaction.id+"_d1.2")
           constraint_limits.append(alpha)
           the_linear_expressions.append(SparsePair(ind=[reaction.id,"d_"+reaction.id],val=[1,alpha]))
           
           constraint_sense.append("G")
           constraint_names.append(reaction.id+"_d2.1")
           constraint_limits.append(10) #Changed from 0 to prevent DeltaGibs of 0.
           the_linear_expressions.append(SparsePair(ind=["Gibs_"+reaction.id,"d_"+reaction.id],val=[-1,alpha]))
           
           constraint_sense.append("L")
           constraint_names.append(reaction.id+"_d2.2")
           constraint_limits.append(alpha)
           the_linear_expressions.append(SparsePair(ind=["Gibs_"+reaction.id,"d_"+reaction.id],val=[-1,alpha]))           
           
           
           
           #print reaction.id
           variable_list = []
           coefficient_list = []
           for metabolite in reaction._metabolites:
               if metabolite.id in ref_con_dict:
                  variable_list.append(metabolite.id+"_log_con")   
                  coefficient_list.append(reaction._metabolites[metabolite])
           variable_list.append("Gibs_"+reaction.id)
           coefficient_list.append(-1.0/RT)
           constraint_sense.append("E")
           constraint_names.append("Gibs Constraint")
           constraint_limits.append(-1*deltaGibs_dict[reaction.id]/(RT))
           the_linear_expressions.append(SparsePair(ind=variable_list,val=coefficient_list))   
           
           #Constraint that Gibs must be larger than 500 or smaller than -500 
           variable_names.append("B_Gibs_"+reaction.id)  
           lower_bounds.append(0)
           upper_bounds.append(1)
           objective_coefficients.append(0.0)
           variable_kinds.append(Cplex.variables.type.binary)
           
           if reaction.id not in Free_E_lower_bound:
              Free_E_lower_bound[reaction.id]=500
           
           constraint_sense.append("G")
           constraint_names.append(reaction.id+"_Gibs_absolute_value1")
           constraint_limits.append(Free_E_lower_bound[reaction.id])
           the_linear_expressions.append(SparsePair(ind=["Gibs_"+reaction.id,"B_Gibs_"+reaction.id],val=[1,alpha]))
           
           constraint_sense.append("G")
           constraint_names.append(reaction.id+"_Gibs_absolute_value2")
           constraint_limits.append(Free_E_lower_bound[reaction.id]-alpha)
           the_linear_expressions.append(SparsePair(ind=["Gibs_"+reaction.id,"B_Gibs_"+reaction.id],val=[-1,-alpha]))
           #Constraint that Gibs must be smaller than a given value in absolute value: 
           if reaction.id in Free_E_upper_bound: 
                                   
              constraint_sense.append("L")
              constraint_names.append(reaction.id+"_Gibs_absolute_value1b")
              constraint_limits.append(Free_E_upper_bound[reaction.id])
              the_linear_expressions.append(SparsePair(ind=["Gibs_"+reaction.id],val=[1]))
              
              constraint_sense.append("L")
              constraint_names.append(reaction.id+"_Gibs_absolute_value2b")
              constraint_limits.append(Free_E_upper_bound[reaction.id])
              the_linear_expressions.append(SparsePair(ind=["Gibs_"+reaction.id],val=[-1]))
               
    
    
    
    n_non_quad_variables=n_variables+ len(variable_names)
    n_total_variables =n_non_quad_variables+len(ref_con_dict)
    quadratic_objective=dok_matrix((n_total_variables,n_total_variables))
    for i,metabolite in enumerate(ref_con_dict):
           variable_names.append(metabolite+"_difference")  
           lower_bounds.append(-100)
           upper_bounds.append(100)
           objective_coefficients.append(0.0)
           variable_kinds.append(Cplex.variables.type.continuous)
 
           constraint_sense.append("E")
           constraint_names.append(metabolite+"_con_difference")
           constraint_limits.append(math.log(float(ref_con_dict[metabolite])/1000))
           the_linear_expressions.append(SparsePair(ind=[metabolite+"_log_con",metabolite+"_difference"],val=[1,-1]))       
    
           quadratic_objective[n_non_quad_variables+i,n_non_quad_variables+i]=1           
           
    lp.variables.add(obj=objective_coefficients,
                         lb=lower_bounds,
                         ub=upper_bounds,
                         names=variable_names,
                         types=variable_kinds)
    
    lp.linear_constraints.add(lin_expr=the_linear_expressions,
                                  rhs=constraint_limits,
                                  senses=constraint_sense,
                                  names=constraint_names)
    
    set_quadratic_objective(lp, quadratic_objective)
    problem_type = Cplex.problem_type.MIQP
    lp.set_problem_type(problem_type)
    lp.objective.set_sense(1) #Quadratic problem must be minimized
    
          
    return lp    

def get_concentrations(model):
    met_cons_dict={}
    for metabolite in model.metabolites:
        met_name=metabolite.id+"_log_con"
        if met_name  in model.solution.x_dict:
           conc=round(math.exp(model.solution.x_dict[met_name])*1000,4)
           #print met_name+"  "+str(conc)
           met_cons_dict[metabolite.id]=conc
    return met_cons_dict

  
            


def thermodynamic_fba(model):
    problem=thermo_problem(model)
    problem.solve()
    model.solution=format_solution(problem, model)
    con_dict=get_concentrations(model)
    flux_dict={}
    for x in model.reactions:
       if x.id in model.solution.x_dict:
          flux_dict[x.id]=round(model.solution.x_dict[x.id],9)
       else:
          flux_dict[x.id]=0
    return (con_dict,flux_dict)
