from ..flux_functions.apply_ratios import remove_ratio

def clear_parameters(label_model,parameter_dict=None,parameter_list=[], clear_ratios=True,clear_turnover=True,clear_fluxes=True,restore_objectives=True,delete_parameters=False):
    if parameter_dict==None:
       parameter_dict=label_model.parameter_dict
    if parameter_list==[]:
       parameter_list=parameter_dict.keys()
    for parameter in parameter_list:
        if  parameter_dict[parameter]["type"]=="flux value" and clear_fluxes:
            for reaction_id in label_model.parameter_dict[parameter]["reactions"]:
                     reaction_is_objective=False
                     reaction=label_model.constrained_model.reactions.get_by_id(reaction_id)
                     if "original_objective_coefficient" in label_model.parameter_dict[parameter] and not restore_objectives:
                         if label_model.parameter_dict[parameter]["original_objective_coefficient"]!=0:
                            reaction_is_objective=True
                     if not reaction_is_objective:
                        reaction.lower_bound=label_model.parameter_dict[parameter]["original_lb"]
                        reaction.upper_bound=label_model.parameter_dict[parameter]["original_ub"]
                     else:
                        reaction.lower_bound=label_model.parameter_dict[parameter]["lb"]
                        reaction.upper_bound=label_model.parameter_dict[parameter]["ub"]
                     if restore_objectives:
                        if "original_objective_coefficient" in label_model.parameter_dict[parameter]:
                            reaction.objective_coefficient=label_model.parameter_dict[parameter]["original_objective_coefficient"]
                     
        elif label_model.parameter_dict[parameter]["type"]=="turnover" and clear_turnover:
              for reaction_id in label_model.parameter_dict[parameter]["reactions"]:
                  label_model.turnover_flux_dict[reaction_id]["v"]=(label_model.turnover_flux_dict[reaction_id]["ub"]+label_model.turnover_flux_dict[reaction_id]["lb"])/2.0
        elif  label_model.parameter_dict[parameter]["type"]=="ratio" and clear_ratios:
               remove_ratio(label_model.constrained_model,parameter,label_model.ratio_dict)
    if delete_parameters==True:
       for parameter in parameter_list:
            del(label_model.parameter_dict[parameter])
               
