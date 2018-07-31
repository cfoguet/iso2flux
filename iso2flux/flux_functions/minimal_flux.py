
import copy
import cobra
from cobra import Model, Reaction, Metabolite
from flux_variability_analysis import flux_variability_analysis
import math
from convert_model_to_irreversible import convert_to_irreversible_with_indicators
from cobra.manipulation.modify import convert_to_irreversible
from ..misc.round_functions import round_up, round_down 
from ..input_functions.read_metabolomics_data import read_metabolomics_data
from add_turnover_metabolites import add_turnover_metabolites


def create_minimal_fux_model(metabolic_model,fraction_of_optimum_objective=0.8, fraction_of_flux_minimum=2,mutually_exclusive_directionality_constraint=True,boundaries_precision=0.001,label_model=None,metabolite_list_file_name=None,flux_penalty_dict={},maximum_flux=None,extra_constraints_dict={}):
    print "creating irreversible model"
    model=metabolic_model.copy()
    for reaction_id in extra_constraints_dict:
        reaction=model.reactions.get_by_id(reaction_id)
        if "lb" in extra_constraints_dict[reaction_id]:
            reaction.lower_bound=round_down(extra_constraints_dict[reaction_id]["lb"],6)
        if "ub" in extra_constraints_dict[reaction_id]:
            reaction.upper_bound=round_up(extra_constraints_dict[reaction_id]["ub"],6)
    #original_objective=copy.deepcopy(model.objective)
    for reaction in model.reactions:
       if reaction.objective_coefficient!=0:
         fva=flux_variability_analysis(model,reaction_list=[reaction.id], fraction_of_optimum=fraction_of_optimum_objective)
         if reaction.objective_coefficient>0:
            reaction.lower_bound=min(fva[reaction.id]["maximum"],reaction.upper_bound)
            reaction.upper_bound=min(fva[reaction.id]["maximum"],reaction.upper_bound)
         elif reaction.objective_coefficient<0:
            reaction.lower_bound=max(fva[reaction.id]["minimum"],reaction.lower_bound)
            reaction.upper_bound=max(fva[reaction.id]["minimum"],reaction.lower_bound)
         reaction.lower_bound-=boundaries_precision/2.0
         reaction.upper_bound+=boundaries_precision/2.0
         reaction.objective_coefficient=0
    reversible_reactions=[]
    convert_to_irreversible(model)
    for reaction in model.reactions:
        if reaction.lower_bound<0:
           reversible_reactions.append(reaction.id)  
    #convert_to_irreversible_with_indicators(model,reversible_reactions,metabolite_list=[], mutually_exclusive_directionality_constraint = False)
    if metabolite_list_file_name!=None and metabolite_list_file_name!="":
       metabolite_id_list=read_metabolomics_data(model,metabolite_list_fname=metabolite_list_file_name)
    else:
       metabolite_id_list=[]
    add_turnover_metabolites(model, metabolite_id_list=metabolite_id_list, epsilon=1e-6,label_model=label_model)
    flux_reporter = Metabolite('flux_reporter_0',formula='',name='',compartment='0')
    for reaction in model.reactions:
        if "IRRMILP_" in reaction.id or "RATIO_" in reaction.id or "RGROUP_" in reaction.id or "TMS_" in reaction.id:
            continue
        if reaction.id in flux_penalty_dict:
           reporter_coef=flux_penalty_dict[reaction.id]
        elif "reflection" in reaction.notes: #If a reversible reaction has no penalty assigned test if the reverse has a penalty
           reflection=reaction.notes["reflection"]
           if reflection in flux_penalty_dict:
              reporter_coef=flux_penalty_dict[reflection]
           else:
              reporter_coef=1
              #print "Warning: reaction "+reaction.id+" assigned default weight of 1"
        else:
              reporter_coef=1
              #print "Warning: reaction "+reaction.id+" assigned default weight of 1"
        if reporter_coef!=0:
           reaction.add_metabolites({flux_reporter:reporter_coef})        
    total_flux_reaction = Reaction('total_flux')
    total_flux_reaction.name = 'total_flux'
    total_flux_reaction.subsystem = 'misc'
    total_flux_reaction.lower_bound = 0  
    total_flux_reaction.upper_bound = 1000*len(model.reactions)
    total_flux_reaction.add_metabolites({flux_reporter:-1}) 
    model.add_reaction(total_flux_reaction)
    total_flux_reaction.objective_coefficient=-1  
    minimal_flux=-1*model.optimize().f
    if maximum_flux!=None:
      total_flux_reaction.upper_bound = maximum_flux
    else:
      total_flux_reaction.upper_bound = minimal_flux*fraction_of_flux_minimum
    total_flux_reaction.objective_coefficient=0.0
    """for objective in original_objective:
        reaction=model.reactions.get_by_id(objective.id)
        reaction.lower_bound=(fva[reaction.id]["minimum"]-boundaries_precision/2.0)
        reaction.upper_bound=(fva[reaction.id]["maximum"]+boundaries_precision/2.0)"""
    milp_reactions=model.reactions.query("IRRMILP_")
    non_milp_reactions=[]
    for reaction in  model.reactions:
      if reaction not in milp_reactions:
         non_milp_reactions.append(reaction)
    return model,minimal_flux,non_milp_reactions 


def add_flux_limit_constraints(metabolic_model,fraction_of_optimum_objective=0.8, fraction_of_flux_minimum=2,boundaries_precision=0.001,label_model=None,metabolite_list_file_name=None):
    precision=int(-1*(math.log10(boundaries_precision)))
    irreversible_model,minimal_flux,non_milp_reactions=create_minimal_fux_model(metabolic_model,fraction_of_optimum_objective, fraction_of_flux_minimum,boundaries_precision=boundaries_precision,label_model=label_model,metabolite_list_file_name=metabolite_list_file_name)
    irreversible_fva=flux_variability_analysis(irreversible_model,reaction_list=non_milp_reactions,fraction_of_optimum=0,solver=solver)
    reversible_fva={}
    reverese_reactions=irreversible_model.reactions.query("_reverse") 
    for reaction in irreversible_model.reactions:
        if reaction in reverese_reactions or "IRRMILP_" in reaction.id:
           continue
        reversible_fva[reaction.id]=copy.deepcopy(irreversible_fva[reaction.id])
        if "reflection" in reaction.notes:
            reverse_reaction= reaction.notes["reflection"]
            reversible_fva[reaction.id]["minimum"]=-irreversible_fva[reverse_reaction]["maximum"]
            if irreversible_fva[reverse_reaction]["minimum"]!=0: 
               reversible_fva[reaction.id]["maximum"]=-irreversible_fva[reverse_reaction]["minimum"] 
        if reaction.id in metabolic_model.reactions:
           original_model_reaction=metabolic_model.reactions.get_by_id(reaction.id) 
           if original_model_reaction.lower_bound!=0:
              original_model_reaction.lower_bound=max(round_down(reversible_fva[reaction.id]["minimum"],precision),original_model_reaction.lower_bound)
           if original_model_reaction.upper_bound!=0:
              original_model_reaction.upper_bound=min(round_up(reversible_fva[reaction.id]["maximum"],precision),original_model_reaction.upper_bound)
           original_model_reaction.objective_coefficient=0.0
        
    return reversible_fva  






def flux_minimization_fva(irreversible_model,solver=None,reaction_list=[],lp_tolerance_feasibility=1e-6):
  reaction2test=[]
  if reaction_list!=[]:
       for reaction_id in reaction_list:
           reaction2test.append(reaction_id)  
           reaction=irreversible_model.reactions.get_by_id(reaction_id)
           if "reflection" in reaction.notes:
               reflection= reaction.notes["reflection"]
               if reflection not in reaction2test:
                  reaction2test.append(reflection)
       
  else:
     for reaction in irreversible_model.reactions:
         if  "IRRMILP_" not in reaction.id:
             reaction2test.append(reaction.id) 
  #print reaction2test
  #print len(reaction2test)
  #irreversible_fva=flux_variability_analysis(irreversible_model,reaction_list=reaction2test,fraction_of_optimum=0,solver=solver)
  irreversible_fva={}
  model_reaction_n_dict={}
  for n,reaction in enumerate(irreversible_model.reactions):
      model_reaction_n_dict[reaction.id]=n
  for reaction_id in reaction2test:
      reaction=irreversible_model.reactions.get_by_id(reaction_id)
      if "reflection" in reaction.notes :
               #print "------------------------------------"
               irreversible_model2=irreversible_model
               reaction=irreversible_model2.reactions.get_by_id(reaction_id)
               reflection_id= reaction.notes["reflection"]
               n=model_reaction_n_dict[reaction_id]
               n_reflection=model_reaction_n_dict[reflection_id]
               #print "REFLECTIONNNNNNNNNNNN",reaction,reflection_id
               reflection=irreversible_model2.reactions.get_by_id(reflection_id)
               #Test max
               irreversible_model2.reactions.get_by_id(reflection_id).objective_coefficient=-100.0 #make the reflection inacive 
               irreversible_model2.reactions.get_by_id(reaction_id).objective_coefficient=1.0
               #print irreversible_model2.reactions.get_by_id(reflection_id).objective_coefficient,irreversible_model.reactions.get_by_id(reaction_id).objective_coefficient
               max_value=irreversible_model2.optimize(solver=solver,tolerance_feasibility=lp_tolerance_feasibility).x_dict[reaction_id]
               solution_dict=irreversible_model2.optimize(solver=solver,tolerance_feasibility=lp_tolerance_feasibility).x_dict
               #print [[reaction_id,solution_dict[reaction_id]],[reflection_id,solution_dict[reflection_id]]]
               #reflection.objective_coefficient=1
               irreversible_model3=irreversible_model
               irreversible_model3.reactions.get_by_id(reflection_id).objective_coefficient=0 #make the reflection inacive 
               irreversible_model3.reactions.get_by_id(reaction_id).objective_coefficient=-1.0
               #print irreversible_model.optimize()
               #print irreversible_model3.reactions.get_by_id(reflection_id).objective_coefficient,irreversible_model.reactions.get_by_id(reaction_id).objective_coefficient
               min_value=irreversible_model3.optimize(solver=solver,tolerance_feasibility=lp_tolerance_feasibility).x_dict[reaction_id]
               solution_dict=irreversible_model3.optimize(solver=solver,tolerance_feasibility=lp_tolerance_feasibility).x_dict
               #print [[reaction_id,solution_dict[reaction_id]],[reflection_id,solution_dict[reflection_id]]]
               irreversible_fva[reaction_id]={"minimum":min_value,"maximum":max_value}
               #print "sdfghjk",reaction_id,max_value,min_value
               #print irreversible_fva
               #Restore
               irreversible_model.reactions.get_by_id(reflection_id).objective_coefficient=0
               irreversible_model.reactions.get_by_id(reaction_id).objective_coefficient=0
      else:
        reaction_fva=flux_variability_analysis(irreversible_model,reaction_list=[reaction_id],fraction_of_optimum=0,tolerance_feasibility=lp_tolerance_feasibility)
        irreversible_fva[reaction_id]=reaction_fva[reaction_id]
  reversible_fva={}
  reverese_reactions=irreversible_model.reactions.query("_reverse") 
  #print irreversible_fva
  for reaction_id in irreversible_fva:
        reaction=irreversible_model.reactions.get_by_id(reaction_id)
        if reaction in reverese_reactions or "IRRMILP_" in reaction.id:
           continue
        reversible_fva[reaction.id]=copy.deepcopy(irreversible_fva[reaction.id])
        if "reflection" in reaction.notes:
            reverse_reaction= reaction.notes["reflection"]
            if irreversible_fva[reverse_reaction]["maximum"]!=0:
               reversible_fva[reaction.id]["minimum"]=-irreversible_fva[reverse_reaction]["maximum"]
            if irreversible_fva[reverse_reaction]["minimum"]!=0: 
               reversible_fva[reaction.id]["maximum"]=-irreversible_fva[reverse_reaction]["minimum"] 
  #print reversible_fva  
  return reversible_fva  

