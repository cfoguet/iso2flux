import copy
import math
from ..misc.round_functions import round_up,round_down

def apply_parameters(label_model,parameter_dict=None,apply_flux_values=True,parameter_precision=None,parameter_list=[],model=None):
    """
    Applies the parameters to the COBRA model and the flux turnover dict
    label_model: label_model object
    parameter_dict: dict, optional
	Dict of parameters that should be applied. If none the parameters in the label_model.parameter dict will be taken
    apply_flux_values: bool, optional
        If false the parameters relating flux values wont be applied
    parameter_precision: float, optional
        Precision which will be used to set flux bounds for flux values. If none the default precision defined in the label model object will be used
    parameter_list: list, optional
	List of parameters in the parameter dict that should be applied. If empty all parameters of the dict will be applied
    model: COBRApy model object, optional
        Model where the parameters will be applied. If none is provided it will use the constrained_model in the label_model object
    """
    if parameter_dict=={}:
       return
    if parameter_precision==None:
       parameter_precision=label_model.parameter_precision
    precision=int(-1*(math.log10(parameter_precision)))
    if model==None: 
       model=label_model.constrained_model
    if parameter_dict==None:
       parameter_dict=label_model.parameter_dict
    else: #Check that that functions does not break any existing fuctions
       label_model.parameter_dict=parameter_dict
    if parameter_list==[]:
       parameter_to_apply=parameter_dict.keys()
    else:
       parameter_to_apply=parameter_list
       #print parameter_list
    for parameter in parameter_to_apply:
        if "ratio" in parameter_dict[parameter]:
            label_model.ratio_dict[parameter]=copy.copy(parameter_dict[parameter]["ratio"])
            for reaction in label_model.ratio_dict[parameter]:
                if label_model.ratio_dict[parameter][reaction]=="v":
                   label_model.ratio_dict[parameter][reaction]= parameter_dict[parameter]["v"]
        elif parameter_dict[parameter]["type"]=="flux value" and apply_flux_values==True:
             for reaction_id in parameter_dict[parameter]["reactions"]:
                 reaction=model.reactions.get_by_id(reaction_id)
                 reaction.lower_bound=round_down(parameter_dict[parameter]["v"],precision)
                 reaction.upper_bound=round_up(parameter_dict[parameter]["v"],precision)
                 if  "original_objective_coefficient" in parameter_dict[parameter]:
                      if parameter_dict[parameter]["original_objective_coefficient"]!=0:
                         reaction.original_objective_coefficient=0
                         reaction.lower_bound=max(reaction.lower_bound,parameter_dict[parameter]["lb"])
                         reaction.upper_bound=min(reaction.upper_bound,parameter_dict[parameter]["ub"])
                 if  reaction.lower_bound>=reaction.upper_bound:
                     reaction.upper_bound=reaction.lower_bound+parameter_precision
                 reaction.lower_bound=max(reaction.lower_bound,parameter_dict[parameter]["original_lb"])
                 reaction.upper_bound=min(reaction.upper_bound,parameter_dict[parameter]["original_ub"])
                 """if precision>0:
                    reaction.lower_bound=round(reaction.lower_bound,precision)
                    reaction.upper_bound=round(reaction.upper_bound,precision)""" 
        elif parameter_dict[parameter]["type"]=="turnover":
             for reaction_id in parameter_dict[parameter]["reactions"]:            
               label_model.turnover_flux_dict[reaction_id]["v"]=round(parameter_dict[parameter]["v"],precision)
