import re
import numpy as np
import copy
from warnings import warn
from itertools import izip
#from gurobipy import Model, LinExpr, GRB, QuadExpr
from six import string_types, iteritems
from time import time
from openpyxl import load_workbook
import math

from cobra.flux_analysis.variability import flux_variability_analysis
from cobra import Reaction, Metabolite
from cobra.core.Solution import Solution
from cobra.solvers import solver_dict

from convert_model_to_irreversible import convert_to_irreversible_with_indicators
from ..misc.round_functions import round_up, round_down 
from ..misc.read_spreadsheets import read_spreadsheets
from ..input_functions.read_metabolomics_data import read_metabolomics_data
from add_turnover_metabolites import add_turnover_metabolites

if __name__ == "__main__": 
    from iso2flux.flux_functions.convert_model_to_irreversible import convert_to_irreversible_with_indicators

import math


#From Gimme





def cplex_create_imat_problem(cobra_model,hex_reactions=[],lex_reactions=[],epsilon=1,lex_epsilon=0.0001,imat_lower_bound=0.0,imat_upper_bound=2000.0,objective="imat",**kwargs):
    """
    Formulates and imat problem for the cplex solver
    imat_upper_bound is deprectaed
    """
    #print("creating problem")
    # Process parameter defaults
    from cplex import Cplex, SparsePair
    from cplex.exceptions import CplexError
    from cobra.solvers.cplex_solver import parameter_defaults, set_parameter,  variable_kind_dict
    lp = Cplex()
    lp.set_log_stream(None)
    lp.set_error_stream(None)
    lp.set_warning_stream(None)
    lp.set_results_stream(None)
    
    the_parameters = parameter_defaults
    
    if kwargs:
        the_parameters = parameter_defaults.copy()
        the_parameters.update(kwargs)
    
    """for k, v in iteritems(the_parameters):
        set_parameter(lp, k, v)"""
    
    
    objective_coefficients=[]
    for x in cobra_model.reactions :
       if x.id==objective:
          objective_coefficients.append(1)
       else:
          objective_coefficients.append(0)
       if float(x.objective_coefficient)!=0.0 : 
          warn('previous objective removed')        
    lower_bounds = [float(x.lower_bound) for x in cobra_model.reactions]
    upper_bounds = [float(x.upper_bound) for x in cobra_model.reactions]
    variable_names = cobra_model.reactions.list_attr("id")
    variable_kinds=[]
    for x in cobra_model.reactions:
        if x.variable_kind=="continuous":
            variable_kinds.append(Cplex.variables.type.continuous)
        elif x.variable_kind=="integer":
            if x.lower_bound==0 and x.upper_bound==1:
               variable_kinds.append(Cplex.variables.type.binary)   
            else:
               variable_kinds.append(Cplex.variables.type.integer) 
    for x in hex_reactions:
       objective_coefficients.append(0)
       lower_bounds.append(int(0))
       upper_bounds.append(int(1))
       variable_names.append(x+"_hexs_var")
       variable_kinds.append(Cplex.variables.type.binary)  
    for x in lex_reactions:
       objective_coefficients.append(0)
       lower_bounds.append(int(0))
       upper_bounds.append(int(1))
       variable_names.append(x+"_lexs_var")
       variable_kinds.append(Cplex.variables.type.binary) 
       
    goal_coef=0.0
    if objective=="imat":
       goal_coef=1.0
    objective_coefficients.append(goal_coef)
    lower_bounds.append(int(imat_lower_bound))
    upper_bounds.append(len(cobra_model.reactions))
    variable_names.append("imat_objective")
    variable_kinds.append(Cplex.variables.type.integer) 
    
    #print lower_bounds
    
    lp.variables.add(obj=objective_coefficients,
    lb=lower_bounds,
    ub=upper_bounds,
    names=variable_names,
    types=variable_kinds) 
    constraint_sense = [] #"E"=equal, "G"=Greater, "L" lower
    constraint_names = []
    constraint_limits = []
    
    [(constraint_sense.append(x._constraint_sense),
      constraint_names.append(x.id),
      constraint_limits.append(float(x._bound)))
     for x in cobra_model.metabolites]
    #CFC added v+xH(Vmin-epsilon)>Vmin
    [(constraint_sense.append("G"),
     constraint_names.append(x+"_hexs"),
     #constraint_limits.append(float(cobra_model.reactions.get_by_id(x).lower_bound)))
      constraint_limits.append(0))
     for x in hex_reactions]
    #CFC added V + xL*Vmax<Vmax
    for x in lex_reactions:
         constraint_sense.append("L")
         constraint_names.append(x+"_lexs")
         reaction=cobra_model.reactions.get_by_id(x) 
         if  "reflection" in reaction.notes:
             reflection=cobra_model.reactions.get_by_id(reaction.notes["reflection"]) 
             constraint_limits.append(float(max(reaction.upper_bound,reflection.upper_bound)))
         else:
             constraint_limits.append(float(reaction.upper_bound))
      
    constraint_sense.append("E")
    constraint_names.append("imat_objective_cons")
    constraint_limits.append(0)
     
     
    the_linear_expressions = []
    #NOTE: This won't work with metabolites that aren't in any reaction
    for the_metabolite in cobra_model.metabolites:
        variable_list = []
        coefficient_list = []
        for the_reaction in the_metabolite._reaction:
            variable_list.append(the_reaction.id)
            coefficient_list.append(float(the_reaction._metabolites[the_metabolite]))
        the_linear_expressions.append(SparsePair(ind=variable_list,val=coefficient_list))
#Added by CFC
    imat_variables_list=[]
    for x in hex_reactions:
         imat_variables_list.append(x+"_hexs_var")
         if "reflection" in cobra_model.reactions.get_by_id(x).notes :
            variable_list = [x, cobra_model.reactions.get_by_id(x).notes["reflection"],x+"_hexs_var"]
            coefficient_list = [1,1,-epsilon]     
         else:
            variable_list = [x,x+"_hexs_var"]
            coefficient_list = [1,-epsilon]
         #coefficient_list = [1,(cobra_model.reactions.get_by_id(x).lower_bound)-epsilon]
         the_linear_expressions.append(SparsePair(ind=variable_list,val=coefficient_list))
    for x in lex_reactions:
         imat_variables_list.append(x+"_lexs_var")
         if "reflection" in cobra_model.reactions.get_by_id(x).notes:
            reflection_id=cobra_model.reactions.get_by_id(x).notes["reflection"]
            variable_list = [x,cobra_model.reactions.get_by_id(x).notes["reflection"],x+"_lexs_var"]
            coefficient_list = [1,1,(max(cobra_model.reactions.get_by_id(x).upper_bound,cobra_model.reactions.get_by_id(reflection_id).upper_bound)-lex_epsilon)]
         else: 
            variable_list = [x,x+"_lexs_var"]
            coefficient_list = [1,(cobra_model.reactions.get_by_id(x).upper_bound-lex_epsilon)]
         the_linear_expressions.append(SparsePair(ind=variable_list, val=coefficient_list))
    
    variable_list=["imat_objective"]
    coefficient_list=[1]
    for x in imat_variables_list:
        variable_list.append(x)
        coefficient_list.append(-1) 
    the_linear_expressions.append(SparsePair(ind=variable_list, val=coefficient_list))
    # Set objective to quadratic program
    lp.linear_constraints.add(lin_expr=the_linear_expressions,
                                  rhs=constraint_limits,
                                  senses=constraint_sense,
                                  names=constraint_names)
    #Set the problem type as cplex doesn't appear to do this correctly
    problem_type = Cplex.problem_type.MILP
    lp.set_problem_type(problem_type)
    return(lp)




def cplex_solve_imat(model,hex_reactions=[],lex_reactions=[],epsilon=1,lex_epsilon=0.001,imat_lower_bound=0.0,imat_upper_bound=2000.0,objective="imat",sense="maximize",lp_output=False,problem=None): 
   """
   Solves an imat problem using CPLEX
   model: model object
   hex_reactions: list of reaction IDs
          List of highly expressed reactions
   lex_reactions: list of reactions IDs
          List of lowly expressed reactions
   epsilon: float
            Minimum flux to consider a reaction active
   lex_espilon: float
            Maximum flux to consider a reaction inactive
   imal_lower_bound: int
          Minimal value of the iMAT objective function
   imat_upper_bound=int
          Deprecated
   objective: str
          Objective to be maximized or minimized. It can be either iMAT (default) to optimize the iMat objective function or a reaction ID
   sense: str
          Either "maximize" or "minimize" depending on wether the objective should be maximize or minimized
   lp_output: bool
          Used for debugging
   problem: cplex problem object
          Used for debugging
   
   """
   from cobra.solvers.cplex_solver import format_solution, solve_problem, get_status
   name="problem_"+objective+"_"+sense+".lp"
   if problem==None:
             problem=cplex_create_imat_problem(model, hex_reactions ,lex_reactions, epsilon,lex_epsilon,imat_lower_bound,objective=objective)
   #problem.write(name)
   #problem.read(name) canviat,mira si dona fallo
   problem.set_log_stream(None)
   problem.set_error_stream(None)
   problem.set_warning_stream(None)
   problem.set_results_stream(None)
   """problem.parameters.mip.limits.gomorypass.set(10000)
   problem.parameters.mip.limits.gomorycand.set(10)"""
   problem.parameters.timelimit.set(500)         
   for x in range(0,5):
          
          problem.status='failed'
          if sense=="maximize":
             sense_coef=-1
          else:
             sense_coef=1
          problem.objective.set_sense(sense_coef)
          if x==1:
             f=open("time_limit","a")
             f.write(name+"\n")
             f.close()
             problem.write(name)
             problem.read(name)
          if x==2:
             problem.parameters.mip.limits.gomorypass.set(10000)
             problem.parameters.mip.limits.gomorycand.set(10)
          if x>2:
             problem.parameters.mip.limits.gomorypass.reset()
             problem.parameters.mip.limits.gomorycand.reset()
             print("Tunning parameters")
             problem.parameters.tune_problem()
             print("Done")
          if x==4:
             problem.parameters.timelimit.set(3600)
          #problem.parameters.mip.tolerances.integrality.set(1e-07)
          problem.solve()
          status=get_status(problem)
          print(status)
          #status = solve_problem(problem,objective_sense=sense, lp_method=the_method)
          #print "Solution: "+status
          if status=="optimal":
             break
   
   model.solution=format_solution(problem, model)
   if lp_output:
      return (problem.solution.get_objective_value(), status,problem)
   else: 
      return (problem.solution.get_objective_value(), status)

def guorbi_create_problem_imat(cobra_model,hexs=[],lexs=[],epsilon=1,lex_epsilon=0.0001,objective="imat",quadratic_component=None,imat_lower_bound=0, **kwargs):
    print "gurobi"
    """
    Creates an IMAT problem for the sofware GUROBI      
    """
    from cobra.solvers.gurobi_solver import parameter_defaults, set_parameter, solve_problem, variable_kind_dict, sense_dict
    from gurobipy import Model, LinExpr, GRB, QuadExpr
    lp = Model("")
    
    the_parameters = parameter_defaults
    if kwargs:
        the_parameters = parameter_defaults.copy()
        the_parameters.update(kwargs)
    
    if "verbose" in the_parameters:
        set_parameter(lp, "verbose", the_parameters["verbose"])
    for k, v in iteritems(the_parameters):
        set_parameter(lp, k, v)
    
    # Create variables
    #TODO:  Speed this up
    """variable_list = [lp.addVar(float(x.lower_bound),
                               float(x.upper_bound),
                               0, #float(x.objective_coefficient)
                               variable_kind_dict[x.variable_kind],
                               str(i))"""
    """variable_list = [lp.addVar(float(x.lower_bound),
                               float(x.upper_bound),
                               0, #float(x.objective_coefficient)
                               variable_kind_dict[x.variable_kind],
                               x.id)
                     for i, x in enumerate(cobra_model.reactions)]"""
    variable_list=[]
    for i,x in  enumerate(cobra_model.reactions):
        if objective==x.id:
           objective_coefficient=1
        else:
           objective_coefficient=0 
        variable_list.append(lp.addVar(float(x.lower_bound),
                               float(x.upper_bound),
                               objective_coefficient,
                               variable_kind_dict[x.variable_kind],
                               x.id))
    reaction_to_variable = dict(zip(cobra_model.reactions,
                                    variable_list))
    nvariable=len(variable_list)
    for x in hexs:
        reaction_to_variable[x+"_hexs_var"]=lp.addVar(0,
                               1,
                               0,
                               variable_kind_dict['integer'],
                               str(x+"_hexs_var"))
        nvariable=+1 
    for x in lexs:
        reaction_to_variable[x+"_lexs_var"]=lp.addVar(0,
                               1,
                               0,
                               variable_kind_dict['integer'],
                               str(x+"_lexs_var"))
        nvariable=+1
    if objective=="imat":
       imat_objective_coefficient=1
    else:
       imat_objective_coefficient=0
    reaction_to_variable["imat_ojective"]=lp.addVar(imat_lower_bound,
                               len(hexs+lexs),
                               imat_objective_coefficient,
                               variable_kind_dict['integer'],
                               str(x+"_lexs_var"))
    nvariable=+1  
    # Integrate new variables
    lp.update()
    
    #Constraints are based on mass balance
    #Construct the lin expression lists and then add
    #TODO: Speed this up as it takes about .18 seconds
    #HERE
    for i, the_metabolite in enumerate(cobra_model.metabolites):
        constraint_coefficients = []
        constraint_variables = []
        for the_reaction in the_metabolite._reaction:
            constraint_coefficients.append(the_reaction._metabolites[the_metabolite])
            constraint_variables.append(reaction_to_variable[the_reaction])
        #Add the metabolite to the problem
        lp.addConstr(LinExpr(constraint_coefficients, constraint_variables),
                     sense_dict[the_metabolite._constraint_sense.upper()],
                     the_metabolite._bound,
                     str(i))
    nconstr=len(cobra_model.metabolites)
    imat_objective_constraint_coefficients=[1]
    imat_objective_constraint_variables=[reaction_to_variable["imat_ojective"]]
    for x in hexs:
        varname=x+"_hexs_var"
        the_reaction=cobra_model.reactions.get_by_id(x)
        if "reflection" not in the_reaction.notes:
           constraint_coefficients = [1,-epsilon]
           constraint_variables = [reaction_to_variable[the_reaction],reaction_to_variable[varname]]
           lp.addConstr(LinExpr(constraint_coefficients, constraint_variables),sense_dict["G"],0, str(nconstr))
        else:
           the_reverse_reaction=cobra_model.reactions.get_by_id(the_reaction.notes["reflection"])
           constraint_coefficients = [1,1,-epsilon]
           constraint_variables = [reaction_to_variable[the_reaction],reaction_to_variable[the_reverse_reaction],reaction_to_variable[varname]]
           lp.addConstr(LinExpr(constraint_coefficients, constraint_variables), sense_dict["G"], 0,  str(nconstr))
        nconstr=nconstr+1
        imat_objective_constraint_coefficients.append(-1)
        imat_objective_constraint_variables.append(reaction_to_variable[varname])
    for x in lexs:
        varname=x+"_lexs_var"
        the_reaction=cobra_model.reactions.get_by_id(x)
        if "reflection" not in the_reaction.notes:
           constraint_coefficients = [1,(the_reaction.upper_bound-(lex_epsilon))]
           #constraint_coefficients = [1,(the_reaction.upper_bound)]
           constraint_variables = [reaction_to_variable[the_reaction],reaction_to_variable[varname]]
           lp.addConstr(LinExpr(constraint_coefficients, constraint_variables), sense_dict["L"],the_reaction.upper_bound, str(nconstr))
        else:
           the_reverse_reaction=cobra_model.reactions.get_by_id(the_reaction.notes["reflection"])
           #constraint_coefficients = [1,1,(the_reaction.upper_bound)]
           constraint_coefficients = [1,1,(max(the_reaction.upper_bound,the_reverse_reaction.upper_bound)-(lex_epsilon))]
           constraint_variables = [reaction_to_variable[the_reaction],reaction_to_variable[the_reverse_reaction],reaction_to_variable[varname]]
           lp.addConstr(LinExpr(constraint_coefficients, constraint_variables), sense_dict["L"],max(the_reaction.upper_bound,the_reverse_reaction.upper_bound), str(nconstr))   
        nconstr=nconstr+1
        imat_objective_constraint_coefficients.append(-1)
        imat_objective_constraint_variables.append(reaction_to_variable[varname])
    
    lp.addConstr(LinExpr(imat_objective_constraint_coefficients,imat_objective_constraint_variables), sense_dict["E"],0, str(nconstr))
    nconstr=nconstr+1
    # Set objective to quadratic program
    if quadratic_component is not None:
        set_quadratic_objective(lp, quadratic_component)
    
    lp.update()
    return(lp)


def gurobi_format_solution_imat(lp, cobra_model, **kwargs):
    """Formats the solution of the Gurobi IMAT problem"""
    from cobra.solvers.gurobi_solver import get_status
    from cobra.core.Solution import Solution
    status = get_status(lp)
    x=[]
    x_dict={}
    if status not in ('optimal', 'time_limit'):
        the_solution = Solution(None, status=status)
        objective_value=0
        
    else:
        objective_value = lp.ObjVal
        x = [v.X for v in lp.getVars()]
        names = [v.VarName for v in lp.getVars()]
        x_dict=dict(izip(names,x))
    the_solution = Solution(objective_value, x=x, x_dict=x_dict, y=[],
                                y_dict={}, status=status)
    return(the_solution)

def gurobi_solve_imat(cobra_model,hexs=[],lexs=[],epsilon=1,lex_epsilon=0.0001,imat_lower_bound=0,objective="imat",sense="maximize", **kwargs):
    """
   Solves an imat problem using CPLEX
   model: model object
   hex_reactions: list of reaction IDs
          List of highly expressed reactions
   lex_reactions: list of reactions IDs
          List of lowly expressed reactions
   epsilon: float
            Minimum flux to consider a reaction active
   lex_espilon: float
            Maximum flux to consider a reaction inactive
   imal_lower_bound: int
          Minimal value of the iMAT objective function
   objective: str
          Objective to be maximized or minimized. It can be either iMAT (default) to optimize the iMat objective function or a reaction ID
   sense: str
          Either "maximize" or "minimize" depending on wether the objective should be maximize or minimized
   """
    from gurobipy import Model, LinExpr, GRB, QuadExpr
    from cobra.solvers.gurobi_solver import solve_problem
    for i in ["new_objective", "update_problem", "the_problem"]:
        if i in kwargs:
            raise Exception("Option %s removed" % i)
    if 'error_reporting' in kwargs:
        warn("error_reporting deprecated")
        kwargs.pop('error_reporting')  
    #Create a new problem
    lp= guorbi_create_problem_imat(cobra_model,hexs,lexs,epsilon,lex_epsilon,objective=objective,sense=sense,imat_lower_bound=imat_lower_bound, **kwargs)
    #Define objective sense
    if sense=="minimize":
       objective_sense=GRB.MINIMIZE
    else:
       objective_sense=GRB.MAXIMIZE
    lp.setAttr('ModelSense', objective_sense)
    ###Try to solve the problem using other methods if the first method doesn't work
    try:
        lp_method = kwargs['lp_method']
    except:
        lp_method = 0
    the_methods = [0, 2, 1]
    if lp_method in the_methods:
        the_methods.remove(lp_method)
    #Start with the user specified method
    the_methods.insert(0, lp_method)
    for the_method in the_methods:
        try:
            status = solve_problem(lp, lp_method=the_method)
        except:
            status = 'failed'
        if status == 'optimal':
            break
    cobra_model.solution=gurobi_format_solution_imat(lp, cobra_model) 
    print cobra_model.solution.f, status
    return (cobra_model.solution.f, status) 


def imat(cobra_model,hex_reactions=[],lex_reactions=[],epsilon=1,lex_epsilon=0.001,objective="imat",imat_lower_bound=0.0,sense="maximize",output=False,solver="cplex"):
    """
    Solves an imat problem
    cobra_model: model object
    hex_reactions: list of reaction IDs
          List of highly expressed reactions
    lex_reactions: list of reactions IDs
          List of lowly expressed reactions
    epsilon: float
            Minimum flux to consider a reaction active
    lex_espilon: float
            Maximum flux to consider a reaction inactive
    imal_lower_bound: int
          Minimal value of the iMAT objective function
    objective: str
          Objective to be maximized or minimized. It can be either iMAT (default) to optimize the iMat objective function or a reaction ID
    sense: str
          Either "maximize" or "minimize" depending on wether the objective should be maximize or minimized
    output: bool
          If True it will print an overvier of the result
    solver: string
          Either "cplex" or "gurobi"
    """
    if solver=="gurobi":
       objective, status=gurobi_solve_imat(cobra_model,hex_reactions,lex_reactions,epsilon,lex_epsilon,imat_lower_bound=imat_lower_bound,objective=objective,sense=sense) #TODO complete parameters
    elif solver=="cplex":
       objective, status=cplex_solve_imat(cobra_model,hex_reactions,lex_reactions,epsilon,lex_epsilon,imat_lower_bound,imat_upper_bound=2000.0,objective=objective,sense=sense)
    
    if output==True:
      f=open("exp_consistency_"+solver,"w")
      f.write("Highly expressed:\n")
      for x in hex_reactions:
          if "reflection" not in cobra_model.reactions.get_by_id(x).notes: 
             f.write(x+" "+str(cobra_model.solution.x_dict[x+"_hexs_var"])+" "+str(cobra_model.solution.x_dict[x])+"\n")
          else:
             f.write(x+" "+str(cobra_model.solution.x_dict[x+"_hexs_var"])+" "+str(cobra_model.solution.x_dict[x]-cobra_model.solution.x_dict[cobra_model.reactions.get_by_id(x).notes["reflection"]])+"\n") 
      f.write("Lowly expressed:\n")
      for x in lex_reactions: 
         if "reflection"  not in cobra_model.reactions.get_by_id(x).notes:
            f.write(x+" "+str(cobra_model.solution.x_dict[x+"_lexs_var"])+" "+str(cobra_model.solution.x_dict[x])+"\n")
         else: 
            f.write(x+" "+str(cobra_model.solution.x_dict[x+"_lexs_var"])+" "+str(cobra_model.solution.x_dict[x]-cobra_model.solution.x_dict[cobra_model.reactions.get_by_id(x).notes["reflection"]])+"\n") 
      f.close()
     
    return objective, status 


def imat_variability(model,hexs,lexs,epsilon,lex_epsilon,fraction_of_optimum=1.0,check_hex=True,all_reactions=False,solver="cplex"):
    """
    Performs a flux variabilaity analysis using iMAT
    model: model object
    hex_reactions: list of reaction IDs
          List of highly expressed reactions
    lex_reactions: list of reactions IDs
          List of lowly expressed reactions
    epsilon: float
            Minimum flux to consider a reaction active
    lex_espilon: float
            Maximum flux to consider a reaction inactive
    fraction of optimum: float
          Minimal fraction of the iMAT objective that must be achieved
    check_hex: bool
          If false only the variability of lowly expressed reactions will be evaluated. This can be used when the goal using iMAT is to remove lowly expressed reactions
    all_reactions: bool
          If False the only the variability of highly expressed and lowly expressed reactions will be evaluated. If True all reactions will be evluated
    solver: string
          Either "cplex" or "gurobi"
    """
    optimal, status=imat(model,hexs,lexs,epsilon,lex_epsilon,solver=solver)
    #print [optimal,fraction_of_optimum]
    fraction_optimal=int(optimal*fraction_of_optimum)
    fva={}
    count=0
    import datetime
    if all_reactions:
       fva_list=[]
       for x in model.reactions:
           if "_reverse" in x.id:
              continue
           if "TMS_" in x.id:
              continue
           if  "RGROUP_" in x.id:
              continue
           fva_list.append(x.id) 
    else:
       if check_hex:
          fva_list=hexs+lexs
       else:
          fva_list=lexs
       
       
    n=str(len(fva_list))
    for x  in fva_list:
           count=count+1
           fva_entry={} 
           max_val, max_status=imat(model,hex_reactions=hexs,lex_reactions=lexs,epsilon=epsilon,lex_epsilon=lex_epsilon,objective=x,imat_lower_bound=fraction_optimal,sense="maximize",output=False,solver=solver)
           if max_status!="optimal":
              max_val=model.reactions.get_by_id(x).upper_bound
           min_val, min_status=imat(model,hex_reactions=hexs,lex_reactions=lexs,epsilon=epsilon,lex_epsilon=lex_epsilon,objective=x,imat_lower_bound=fraction_optimal,sense="minimize",output=False,solver=solver)
           if min_status!="optimal":
              min_val=model.reactions.get_by_id(x).lower_bound
           #print str(count)+"/"+str(n)+" "+x+" "+str(max_val)+" "+max_status+" "+str(min_val)+" "+min_status
           if "reflection" in model.reactions.get_by_id(x).notes:
              reflection_id=model.reactions.get_by_id(x).notes["reflection"]
              reflection=model.reactions.get_by_id(reflection_id)
              max_val_r, max_status=imat(model,hex_reactions=hexs,lex_reactions=lexs,epsilon=epsilon,lex_epsilon=lex_epsilon,objective=reflection_id,imat_lower_bound=fraction_optimal,sense="maximize",output=False,solver=solver)
              if max_status!="optimal":
                 max_val_r=reflection.upper_bound
              min_val_r, min_status=imat(model,hex_reactions=hexs,lex_reactions=lexs,epsilon=epsilon,lex_epsilon=lex_epsilon,objective=reflection_id,imat_lower_bound=fraction_optimal,sense="minimize",output=False,solver=solver)
              if min_status!="optimal":
                 min_val_r=reflection.lower_bound
              if max_val_r>0:
                 min_val=-1*max_val_r
              if min_val_r>0:
                 max_val=-1*min_val_r
           fva_entry["maximum"]= max_val
           fva_entry["minimum"]= min_val          
           fva[x]=fva_entry
           f=open("imat_log","a")
           f.write(str(count)+"/"+n+" "+str(x)+" "+str(max_val)+" "+str(min_val)+"\n")
           f.close()
           print(str(count)+"/"+n+" "+str(x)+" "+str(max_val)+" "+str(min_val))
           
          
    
    #count=0   
    #for x in lexs:
       
           
    
    return fva


def create_gim3e_model(cobra_model,file_name="gene_expression_data.xlsx",metabolite_list=[],label_model=None,epsilon=0.0001,gene_method="average",gene_prefix="",gene_sufix="",low_expression_threshold=25,absent_gene_expression=100,percentile=True):
    """
    Creates a Gim3e model
    """
    reaction_expression_dict,genexpraw_dict,expression_dict=get_gene_exp(cobra_model,absent_gene_expression=absent_gene_expression,percentile=percentile,file_name=file_name,gene_method=gene_method,gene_prefix=gene_prefix,gene_sufix=gene_sufix)
    print reaction_expression_dict
    if percentile==True: #If percentile is set to True, low_expression_threshold is assumed to be a percintile
       value_list=[] 
       for x in genexpraw_dict:
           value_list+=genexpraw_dict[x]     
       low_expression_threshold=lower_percentile=np.percentile(value_list,low_expression_threshold)
    print low_expression_threshold
    penalty_dict={}
    for reaction_id in  reaction_expression_dict:
        gene_expression_value=reaction_expression_dict[reaction_id]
        if gene_expression_value< low_expression_threshold:
           gene_expression_penalty=round(-gene_expression_value+low_expression_threshold,6)
           penalty_dict[reaction_id]=gene_expression_penalty
    print penalty_dict    
    convert_to_irreversible_with_indicators( cobra_model,penalty_dict.keys(),metabolite_list=[], mutually_exclusive_directionality_constraint = True,label_model=label_model)
    add_turnover_metabolites(cobra_model, metabolite_id_list=metabolite_list, epsilon=epsilon,label_model=label_model)
    
    objective_reaction = Reaction('gim3e_objective')
    objective_reaction.objective_coefficient=-1
    gim3e_indicator = Metabolite('gim3e_indicator',formula='',name='',compartment='')
    objective_reaction.add_metabolites({gim3e_indicator:-1})
    cobra_model.add_reaction(objective_reaction)
    
    total_bound=0         
    for reaction_id in penalty_dict:
           gene_expression_penalty=penalty_dict[reaction_id]
           reaction=cobra_model.reactions.get_by_id(reaction_id)
           reaction.add_metabolites({gim3e_indicator:gene_expression_penalty})
           total_bound+=gene_expression_penalty*reaction.upper_bound
           if "reflection" in reaction.notes:
               reflection_id=reaction.notes["reflection"]
               reflection=cobra_model.reactions.get_by_id(reflection_id)
               reflection.add_metabolites({gim3e_indicator:gene_expression_penalty})
               total_bound+=gene_expression_penalty*reflection.upper_bound
    print total_bound
    objective_reaction.lower_bound=0.0
    objective_reaction.upper_bound=total_bound
    return penalty_dict


def get_expression(model,file_name="gene_expression_data.xlsx",gene_method="average",gene_prefix="",gene_sufix=""):
    """
    Reads the gene expression file
    model: cobra model object
    file_name: str
          Name of the file with the gene expression data. It should be either a CSV or a XLSX file with gene names in the first column and gene expression value in the second column
    gene_method: str
          Determines wich gene expression value should be given to a gene when there are multiples entries for the same gene. It can either be "average" (it will use the average of all values), "maximum" (it will use the maximum) or "minimum" (it will use the minimum)
    gene_prefix: str
          Prefix used by genes in the cobra model but not present in the gene expression file
    gene_sufix= str
          Sufix  used by genes in the cobra model but not present in the gene expression file. In Recon 1 alternative transtricpts are indicated by appending _AT1, _AT2 , AT3_ at the end of gene. If in the gene expression file alternative scripts are not indicated in that case _AT should be defined as Sufix
    """
    genexpraw_dict={}
    spreadsheet_dict=read_spreadsheets(file_names=file_name,csv_delimiter=',',more_than_1=False,tkinter_title="Chose a file")
    gene_expression_dict=gene_expression_dict={}
    #wb = load_workbook(file_name, read_only=True)
    #ws=wb.active
    for sheet in spreadsheet_dict:
      for row in spreadsheet_dict[sheet]:
        geneid=str(row[0])
        gene_expression_value=float(row[1])
        if geneid==None or geneid=="":
           continue
        if gene_expression_value==None or geneid=="":
           continue  
        gene_matches=model.genes.query(geneid)
        regular_expression=""
        if gene_prefix=="":
           regular_expression+="^"
        else:
           regular_expression+=gene_prefix 
        regular_expression+=geneid
        if gene_sufix=="":
           regular_expression+="$"
        elif gene_sufix==".":
           regular_expression+="\." 
        else:
           regular_expression+=gene_sufix  
        print regular_expression
        for gene in gene_matches:
            if re.search(regular_expression,gene.id)!=None: 
               if gene.id not in genexpraw_dict:
                  genexpraw_dict[gene.id]=[]
               genexpraw_dict[gene.id].append(gene_expression_value)
    for geneid in genexpraw_dict:
        list_of_values=genexpraw_dict[geneid]
        if gene_method=="average" or "mean":
           value=np.mean(list_of_values)
        elif gene_method in("max","maximum"):
           value=max(list_of_values)
        elif gene_method in ("min","minimum"):
           value=min(list_of_values)
        else:
            print "Error: Method not supoprted"
            return {}
        gene_expression_dict[geneid]=round(value,4)     
    print gene_expression_dict      
    return genexpraw_dict, gene_expression_dict





def get_gene_exp(model,absent_gene_expression=50,percentile=True,file_name="gene_expression_data.xlsx",gene_method="average",gene_prefix="",gene_sufix=""):
    """
    Assigns each reaction a expression value based on the gene_expression file and the GPR rules
    """
    genexpraw_dict,expression_dict=get_expression(model,file_name=file_name,gene_method=gene_method,gene_prefix=gene_prefix,gene_sufix=gene_sufix)
    if percentile==True: #If percentile is set to True, low_expression_threshold and high_expression_threshold
       value_list=[] 
       for x in genexpraw_dict:
           value_list+=genexpraw_dict[x]  
         
       absent_gene_expression=np.percentile(value_list,absent_gene_expression)
    f=open("ReactionExpression","w")     
    for gene in model.genes:
        if gene.id not in expression_dict:
           expression_dict[gene.id]=absent_gene_expression
    reaction_expression_dict={}
    for the_reaction in model.reactions:
        if "reflection" in the_reaction.notes:
            if the_reaction.notes["reflection"] in reaction_expression_dict:
               continue
        if the_reaction.gene_reaction_rule != "" : 
           the_gene_reaction_relation = copy.deepcopy(the_reaction.gene_reaction_rule)
           for the_gene in the_reaction.genes: 
               the_gene_re = re.compile('(^|(?<=( |\()))%s(?=( |\)|$))'%re.escape(the_gene.id))
               the_gene_reaction_relation = the_gene_re.sub(str(expression_dict[the_gene.id]), the_gene_reaction_relation) 
           #print the_gene_reaction_relation
           expression_value=evaluate_gene_expression_string( the_gene_reaction_relation)
           reaction_expression_dict[the_reaction.id]=expression_value
           f.write(the_reaction.id+" "+str(expression_value)+"\n")
           #print(the_reaction.id+" "+str(expression_value))    
    f.close()
    return (reaction_expression_dict,genexpraw_dict,expression_dict)


def gene_exp_classify_reactions(model,low_expression_threshold=25,high_expression_threshold=75,percentile=True,file_name="gene_expression_data.xlsx",gene_method="average",gene_prefix="",gene_sufix=""):
    
    reaction_expression_dict,genexpraw_dict,expression_dict=get_gene_exp(model,absent_gene_expression=50,percentile=True,file_name=file_name,gene_method=gene_method,gene_prefix=gene_prefix,gene_sufix=gene_sufix)
    if percentile==True: #If percentile is set to True, low_expression_threshold and high_expression_threshold
       value_list=[] 
       for x in genexpraw_dict:
           value_list+=genexpraw_dict[x]  
         
       high_expression_threshold=upper_percentile=np.percentile(value_list,high_expression_threshold)
       low_expression_threshold=lower_percentile=np.percentile(value_list,low_expression_threshold)
    print [high_expression_threshold,low_expression_threshold]
    hex_reactions=[]
    lex_reactions=[]
    iex_reactions=[] #Inconclusive reactions
    for reaction_id in reaction_expression_dict:
           expression_value=reaction_expression_dict[reaction_id]
           if expression_value > high_expression_threshold: 
              hex_reactions.append(reaction_id) 
           elif expression_value < low_expression_threshold:      
                lex_reactions.append(reaction_id)
           else :      
                iex_reactions.append(reaction_id)
    return (hex_reactions, lex_reactions, iex_reactions )  


number_finder = re.compile("[\d]+\.?[\d]*")
class gene_expression_score:
    """Based on the work by : #Schmidt BJ1, Ebrahim A, Metz TO, Adkins JN, Palsson B, Hyduke DR. GIM3E: condition-specific models of cellular metabolism developed from metabolomics and expression data Bioinformatics. 2013 Nov 15;29(22):2900-8. doi: 10.1093/bioinformatics/btt493. Epub 2013 Aug 23."""
    def __init__(self, value):
        self.str_value = str(value)
        self.value = float(value)
    
    def __add__(self, other):
        # addition is like AND
        return gene_expression_score(min(self.value, other.value))
    
    def __mul__(self, other):
        # multiplication is like OR
        return gene_expression_score(max(self.value, other.value))
    
    def __neg__(self): #CFC Added
        return gene_expression_score (-self.value) #CFC Added

def evaluate_gene_expression_string(gene_expression_string):
    """ penalty string will have:
        * 'or' statements which need to be converted to min
        * 'and' statements which need to be converted to max
    >>> evaluate_penalty("(1 and 2 and 3)")
    max(1, 2, 3)"""
    # if there are no ands or ors, we are done
    gene_expression_string = gene_expression_string.lower()  # don't want to deal with cases
    
    if "and" not in gene_expression_string and "or" not in gene_expression_string:
        return eval(gene_expression_string)
    # we will replace AND and OR with addition and multiplication
    # equivalent to min/max
    #gene_expression_string = gene_expression_string.replace("or", "+").replace("and", "*") #Changed from GIMME
    gene_expression_string = gene_expression_string.replace("or", "*").replace("and", "+")
    # replace the numbers with the custom class which have overloaded AND/OR
    values = [gene_expression_score(i) for i in number_finder.findall(gene_expression_string)]
    values_strings = tuple("values[%i]" % i for i in range(len(values)))
    gene_expression_string = number_finder.sub("%s", gene_expression_string)
    gene_expression_string = gene_expression_string % values_strings
    return eval(gene_expression_string).value




def integrate_omics_imat(metabolic_model,gene_expression_file,fraction_of_optimum=1,low_expression_threshold=25,high_expression_threshold=75,percentile=True,gene_method="average",gene_prefix="",gene_sufix="",metabolite_list_fname=None,epsilon=0.1,lex_epsilon=0.0001,imat_fraction_optimum=1,label_model=None,add_as_constraints=True,boundaries_precision=0.001,solver=None):
   """
   Functions that reads gene expression data and optionally qualtitaive metabolomics data and uses the iMAT algorythm to find the flux distribution that is most consistent with such data. Optionally,you can use the results to constraint the cobra model
   metabolic_model: COBRA model object
             Model where the gene expression data and metabolomics data will be integrated
   gene_expression_file: str
             Name of the XLSX or CSV file that describes gene expression. In the first row it should have the gene identifiers and the second column the gene expression value. It is important that the type of identifier used in the file is the same as the one used in the metabolic model. 
   fraction_of_optimum: float
             Fraction of the model objective reactions that must be fulfilled
   low_expression_threshold: float
             Threshold at which a gene is considered lowly expressed- If percentile is set to True this will be a percentile
   high_expression_threshold: float
             Threshold at which a gene is considered highly expressed- If percentile is set to True this will be a percentile
   percentile: bool
              If True the low_expression_threshold and absent_gene_expression will be percintiles
   gene_method: str
          Determines wich gene expression value should be given to a gene when there are multiples entries for the same gene. It can either be "average" (it will use the average of all values), "maximum" (it will use the maximum) or "minimum" (it will use the minimum)
   gene_prefix: str
          Prefix used by genes in the cobra model but not present in the gene expression file
   gene_sufix= str
          Sufix  used by genes in the cobra model but not present in the gene expression file. In Recon 1 alternative transtricpts are indicated by appending _AT1, _AT2 , AT3_ at the end of gene. If in the gene expression file alternative scripts are not indicated in that case "_AT" should be defined as Sufix
   metabolite_list_fname: string, optional
          The path of a XLSX or CSV file that indicates the metabolites that have been detected in the study conditions. Metabolites can be indicated either from metabolite name or metabolite ID. 
   label_model: label_model object, optional
          If a label_model object is provided the algorythm will ensure that all metabolites whose isotopologues are quantified can be produced 
   epsilon: float
          Minimum flux to consider a reaction active
   lex_epsilon: float
          Maximum flux to consider a reaction innactive. Alos, this is used for setting the lower bound on turnover reactions
   gim3e_fraction_optimum: float
          Fraction of gim3e optimum that must be fulfilled. Must have a value bewteen 0 and 1
   add_as_constraints: bool
          If True the gim3e results will be used to constraint the input model
   boundaries_precision: float 
          If add_as_constraints is set to True this determines the precision of the constraints that will be added
   """
   if percentile in (True,"true","True",1,"1","yes"):
      percentile=True
   else:
      percentile=False
   if solver==None:
      if "cplex" in solver_dict:
          solver="cplex"
      elif "gurobi" in solver_dict:
          solver="gurobi"
      else:
          raise Exception('Error: Integrate gene expression requires either Cplex or Gurobi installed')  
          
   precision=int(-1*(math.log10(boundaries_precision)))
   gene_expression_model=copy.deepcopy(metabolic_model)
   if metabolite_list_fname!=None and metabolite_list_fname!="": 
      metabolite_list=read_metabolomics_data(gene_expression_model,metabolite_list_fname)
   else:
      metabolite_list=[]
   for reaction in gene_expression_model.reactions: #TODO Add support for multiple objectives
       if reaction.objective_coefficient!=0:
          fva=flux_variability_analysis(gene_expression_model,fraction_of_optimum=fraction_of_optimum,reaction_list=[reaction])
          reaction.lower_bound=max(round_down(fva[reaction.id]["minimum"]-boundaries_precision/2.0,precision),reaction.lower_bound)
          reaction.upper_bound=min(round_up(fva[reaction.id]["maximum"]+boundaries_precision/2.0,precision),reaction.upper_bound)  
          reaction.objective_coefficient=0
          break
   
   hexs,lexs,iexs=gene_exp_classify_reactions(metabolic_model,low_expression_threshold=low_expression_threshold,high_expression_threshold=high_expression_threshold,percentile=percentile,file_name=gene_expression_file,gene_method=gene_method,gene_prefix=gene_prefix,gene_sufix=gene_sufix)
   print hexs,lexs
   convert_to_irreversible_with_indicators( gene_expression_model,hexs+lexs,metabolite_list=metabolite_list, mutually_exclusive_directionality_constraint = True,label_model=label_model)
   add_turnover_metabolites(gene_expression_model, metabolite_id_list=metabolite_list, epsilon=epsilon,label_model=label_model)
   print gene_expression_model.optimize(solver=solver)
   objective,status=imat( gene_expression_model,hex_reactions=hexs,lex_reactions=lexs,epsilon=epsilon,lex_epsilon=lex_epsilon,objective="imat",imat_lower_bound=0.0,sense="maximize",output=True,solver=solver)
   print objective,status
   fva=imat_variability(gene_expression_model,hexs,lexs,epsilon=epsilon,lex_epsilon=lex_epsilon,fraction_of_optimum=imat_fraction_optimum,check_hex=True,all_reactions=True,solver=solver)
   print fva
   if add_as_constraints==True:
     for reaction_id in fva: 
       reaction=metabolic_model.reactions.get_by_id(reaction_id)
       lower_bound=fva[reaction_id]["minimum"]  
       upper_bound=fva[reaction_id]["maximum"]
       #reaction.lower_bound=max(round_down(fva[reaction_id]["minimum"],precision),reaction.lower_bound)
       #reaction.upper_bound=min(round_up(fva[reaction_id]["maximum"],precision),reaction.upper_bound)
       reaction.upper_bound=min(round_up(upper_bound,precision),reaction.upper_bound)  
       reaction.lower_bound=max(round_down(fva[reaction_id]["minimum"],precision),reaction.lower_bound)
     """if label_model!=None:
        for group_reaction_id in label_groups_reactions_dict:
            group_reaction=metabolic_model.reactions.get_by_id(group_reaction_id)
            lower_bound=0
            upper_bound=0
            for group_memeber_id in label_groups_reactions_dict["group_reaction_id"]:
                group_member_reaction=metabolic_model.reactions.get_by_id(group_memeber_id)
                coef=label_groups_reactions_dict["group_reaction_id"]["group_memeber_id"]
                lower_bound+=coef*group_member_reaction.lower_bound
                upper_bound+=coef*group_member_reaction.lower_bound"""
     """if abs(lower_bound)<=lex_epsilon:
          reaction.lower_bound=max(0,reaction.lower_bound)
       else:
          reaction.lower_bound=max(round_down(fva[reaction_id]["minimum"],precision),reaction.lower_bound)
       if abs(upper_bound)<=lex_epsilon:
          reaction.upper_bound=min(0,reaction.upper_bound)
       else:
          reaction.upper_bound=min(round_up(upper_bound,precision),reaction.upper_bound)  
       reaction.objective_coefficient=0.0"""
   return hexs,lexs,[objective,int(objective*imat_fraction_optimum)],status, fva


def integrate_omics_gim3e(metabolic_model,gene_expression_file,fraction_of_optimum=1,low_expression_threshold=25,absent_gene_expression=50,percentile=True,gene_method="average",gene_prefix="",gene_sufix="",metabolite_list_fname=None,label_model=None,epsilon=0.0001,gim3e_fraction_optimum=0.75,add_as_constraints=True,boundaries_precision=0.001):
   """
   Functions that reads gene expression data and optionally qualtitaive metabolomics data and uses the GIM3E algorythm to find the flux distribution that is most consistent with such data. Optionally,you can use the results to constraint the cobra model.
   metabolic_model
             Model where the gene expression data and metabolomics data will be integrated
   fraction_of_optimum: float
             Fraction of the model objective reactions that must be fulfilled
   low_expression_threshold: float
             Threshold at which a gene is considered lowly expressed- If percentile is set to True this will be a percentile
   absent_gene_expression: Gene expression level given to genes that are not measured.  If percentile is set to True this will be a percentile
   percentile: bool
              If True the low_expression_threshold and absent_gene_expression will be percintiles
    gene_method: str
          Determines wich gene expression value should be given to a gene when there are multiples entries for the same gene. It can either be "average" (it will use the average of all values), "maximum" (it will use the maximum) or "minimum" (it will use the minimum)
    gene_prefix: str
          Prefix used by genes in the cobra model but not present in the gene expression file
    gene_sufix= str
          Sufix  used by genes in the cobra model but not present in the gene expression file. In Recon 1 alternative transtricpts are indicated by appending _AT1, _AT2 , AT3_ at the end of gene. If in the gene expression file alternative scripts are not indicated in that case "_AT" should be defined as Sufix
   metabolite_list_fname: string, optional
          The path of a XLSX or CSV file that indicates the metabolites that have been detected in the study conditions. Metabolites can be indicated either from metabolite name or metabolite ID. 
   label_model: label_model object, optional
          If a label_model object is provided the algorythm will ensure that all metabolites whose isotopologues are quantified can be produced 
   epsilon: float
          Maximum flux to consider a reaction inaactive. Used for setting the lower bound on turnover reactions
   gim3e_fraction_optimum: float
          Fraction of gim3e optimum that must be fulfilled. Must have a value bewteen 0 and 1
   add_as_constraints: bool
          If True the gim3e results will be used to constraint the input model
   boundaries_precision: float 
          If add_as_constraints is set to True this determines the precision of the constraints that will be added
   """
   if percentile in (True,"true","True",1,"1","yes"):
      percentile=True
   else:
      percentile=False
   precision=int(-1*(math.log10(boundaries_precision)))
   gene_expression_model=copy.deepcopy(metabolic_model)
   if metabolite_list_fname!=None and metabolite_list_fname!="": 
      metabolite_list=read_metabolomics_data(gene_expression_model,metabolite_list_fname)
   else:
      metabolite_list=[]
   if len(gene_expression_model.objective)>1:
      raise Exception("Multiple objectives not supported")
   for reaction in gene_expression_model.objective: #TODO Add support for multiple objectives
       if reaction.objective_coefficient!=0:
          fva=flux_variability_analysis(gene_expression_model,fraction_of_optimum=fraction_of_optimum,reaction_list=[reaction])
          reaction.lower_bound=max(round_down(fva[reaction.id]["minimum"]-boundaries_precision/2.0,precision),reaction.lower_bound)
          reaction.upper_bound=min(round_up(fva[reaction.id]["maximum"]+boundaries_precision/2.0,precision),reaction.upper_bound)  
          reaction.objective_coefficient=0
   penalty_dict=create_gim3e_model(gene_expression_model,file_name=gene_expression_file,metabolite_list=metabolite_list,label_model=label_model,epsilon=epsilon,gene_method=gene_method,gene_prefix=gene_prefix,gene_sufix=gene_sufix,low_expression_threshold=low_expression_threshold,absent_gene_expression=100,percentile=True)
   gene_expression_model.optimize()
   objective=gene_expression_model.solution.f    
   fva=flux_variability_analysis(gene_expression_model,fraction_of_optimum=1.0/gim3e_fraction_optimum)
   if add_as_constraints==True:
     for reaction_id in fva: 
       if reaction_id in metabolic_model.reactions:
          if "RGROUP_" in reaction_id:
              continue
          reaction=metabolic_model.reactions.get_by_id(reaction_id)
          lower_bound=fva[reaction_id]["minimum"]#round_down(fva[reaction_id]["minimum"],precision)  
          upper_bound=fva[reaction_id]["maximum"]#round_up(fva[reaction_id]["maximum"],precision)  
          if "reflection" in gene_expression_model.reactions.get_by_id(reaction_id).notes:
              reflection_id=gene_expression_model.reactions.get_by_id(reaction_id).notes["reflection"]
              reflection_lower_bound=fva[reflection_id]["minimum"]#round_down(fva[reflection_id]["minimum"],precision) 
              reflection_upper_bound=fva[reflection_id]["maximum"]#round_up(fva[reflection_id]["maximum"],precision)
              if reflection_lower_bound==0.0  and reflection_upper_bound==0.0:
                 lower_bound=lower_bound
                 upper_bound=upper_bound
              elif lower_bound==0.0  and upper_bound==0.0:
                 lower_bound=-reflection_upper_bound
                 upper_bound=-reflection_lower_bound
              else:
                 lower_bound=-reflection_upper_bound
                 upper_bound=upper_bound
           
          reaction.lower_bound=max(round_down(lower_bound,precision),reaction.lower_bound)
          reaction.upper_bound=min(round_up(upper_bound,precision),reaction.upper_bound)       
          reaction.objective_coefficient=0 #The new constraints enforce the optimalty of the solution, objective coefficient can be
   return penalty_dict,objective,fva
   
    





