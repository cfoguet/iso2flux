
import copy
import cobra
from cobra import Model, Reaction, Metabolite
from cobra.flux_analysis.variability import flux_variability_analysis
import math
from convert_model_to_irreversible import convert_to_irreversible_with_indicators
from ..misc.round_functions import round_up, round_down 
from ..input_functions.read_metabolomics_data import read_metabolomics_data
from add_turnover_metabolites import add_turnover_metabolites

def create_minimal_fux_model(metabolic_model,fraction_of_optimum_objective=0.8, fraction_of_flux_minimum=2,boundaries_precision=0.001,label_model=None,metabolite_list_file_name=None):
    model=copy.deepcopy(metabolic_model)
    if len(model.objective)>1:
      raise ValueError('Error:Only one objective supported')
    original_objective=copy.deepcopy(model.objective)
    for reaction in model.objective:
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
    #convert_to_irreversible(model)
    for reaction in model.reactions:
        if reaction.lower_bound<0:
           reversible_reactions.append(reaction.id)  
    convert_to_irreversible_with_indicators(model,reversible_reactions,metabolite_list=[], mutually_exclusive_directionality_constraint = True)
    if metabolite_list_file_name!=None and metabolite_list_file_name!="":
       metabolite_id_list=read_metabolomics_data(model,metabolite_list_fname=metabolite_list_file_name)
    else:
       metabolite_id_list=[]
    add_turnover_metabolites(model, metabolite_id_list=metabolite_id_list, epsilon=1e-6,label_model=label_model)
    flux_reporter = Metabolite('flux_reporter_0',formula='',name='',compartment='0')
    reporter_coef=1
    for reaction in model.reactions:
        if "IRRMILP_" in reaction.id or "RATIO_" in reaction.id or "RGROUP_" in reaction.id or "TMS_" in reaction.id:
            continue
        reaction.add_metabolites({flux_reporter:reporter_coef})        
    total_flux_reaction = Reaction('total_flux')
    total_flux_reaction.name = 'total_flux'
    total_flux_reaction.subsystem = 'misc'
    total_flux_reaction.lower_bound = 0  
    total_flux_reaction.upper_bound = 1000*reporter_coef*len(model.reactions)
    total_flux_reaction.objective_coefficient=-1  
    total_flux_reaction.add_metabolites({flux_reporter:-1}) 
    model.add_reaction(total_flux_reaction)
    minimal_flux=-1*model.optimize().f
    total_flux_reaction.upper_bound = minimal_flux*fraction_of_flux_minimum
    total_flux_reaction.objective_coefficient=0.0
    for objective in original_objective:
        reaction=model.reactions.get_by_id(objective.id)
        reaction.lower_bound=(fva[reaction.id]["minimum"]-boundaries_precision/2.0)
        reaction.upper_bound=(fva[reaction.id]["maximum"]+boundaries_precision/2.0)
    milp_reactions=model.reactions.query("IRRMILP_")
    non_milp_reactions=[]
    for reaction in  model.reactions:
      if reaction not in milp_reactions:
         non_milp_reactions.append(reaction)
    return model,minimal_flux,non_milp_reactions 


def add_flux_limit_constraints(metabolic_model,fraction_of_optimum_objective=0.8, fraction_of_flux_minimum=2,solver=None,boundaries_precision=0.001,label_model=None,metabolite_list_file_name=None):
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


