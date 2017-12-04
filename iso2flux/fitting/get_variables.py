import numpy as np
from ..flux_functions.flux_variability_analysis import flux_variability_analysis

def get_variables(label_model):
    variable_vector=np.empty(0)
    variable_list_dict=[]
    for n in label_model.flux_solver_independent_flux_dict:
        #Check if its really variable or if its constant:
        reference_reaction=label_model.flux_solver_independent_flux_dict[n].keys()[0]
        fva=flux_variability_analysis(label_model.constrained_model,reaction_list=[reference_reaction], fraction_of_optimum=0,tolerance_feasibility=label_model.lp_tolerance_feasibility)
        flux_range=fva[reference_reaction]["maximum"]-fva[reference_reaction]["minimum"]
        if flux_range>1e-6:
           variable_list_dict.append({"type":"flux","n":n})
           variable_vector=np.append(variable_vector,label_model.flux_solver_free_fluxes[n])
    for reaction in label_model.turnover_flux_dict:
       if label_model.turnover_flux_dict[reaction]["lb"]!=label_model.turnover_flux_dict[reaction]["ub"]:
          variable_list_dict.append({"type":"turnover","reaction":reaction})
          variable_vector=np.append(variable_vector,label_model.turnover_flux_dict[reaction]["v"])
    
    variable_vector=np.append(variable_vector,0.0)
    variable_list_dict.append({"type":"NULL","reaction":None}) #The last position is saved for use when doing confidence intervals
    label_model.variable_vector=variable_vector
    label_model.variable_list_dict=variable_list_dict
    return variable_vector
