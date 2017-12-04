import copy
import numpy as np
from ..flux_functions.flux_variability_analysis import flux_variability_analysis

def build_flux_solver_variables(label_model):
    #Find out which are the independent fluxes from the coefficients in the nullnmp
    label_model.flux_solver_independent_flux_dict=independent_flux_dict={}
    for n_row, row in enumerate(label_model.flux_solver_nullmnp):
        reaction_id=label_model.flux_solver_n_reaction_dict[n_row]
        count_non_0=0
        for n_col, col in enumerate(row):
            if col!=0:
               count_non_0+=1
               coef=col
               non0_n_col=n_col
        if count_non_0==1:
           #print reaction_id, row
           if non0_n_col not in independent_flux_dict:
              independent_flux_dict[non0_n_col]={}
           independent_flux_dict[non0_n_col][reaction_id]=coef
    
    label_model.flux_solver_free_fluxes=free_fluxes_vector=np.zeros(len(label_model.flux_solver_nullmnp[0]))
    #Start the free fluxes vector with the solution of the constrained model. To prevent errors maje all reactions that are irreversible sliglghy positive
    model_copy=copy.deepcopy(label_model.constrained_model)
    for x in model_copy.reactions:
        if x.lower_bound==0:
           fva=flux_variability_analysis(model_copy,reaction_list=[x.id],fraction_of_optimum=0.0,tolerance_feasibility=1e-9)
           if fva[x.id]["maximum"]>1e-06:
              #print x.id
              x.lower_bound=1e-06
    solution=model_copy.optimize()
    solution_dict={}
    for n,reaction in enumerate(model_copy.reactions):
        solution_dict[reaction.id]=solution.x[n]
    label_model.flux_solver_free_fluxes=free_fluxes_vector=np.zeros(len(label_model.flux_solver_nullmnp[0]))
    for n in independent_flux_dict:
        reaction_id=independent_flux_dict[n].keys()[0]
        coef=independent_flux_dict[n][reaction_id]
        free_fluxes_vector[n]=solution_dict[reaction_id]/coef
