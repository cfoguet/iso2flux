from cobra import Model, Reaction, Metabolite
from ..misc.round_functions import round_up, round_down 
def define_reaction_group(model,reaction_dict,group_reaction_id=None,lower_bound=None,upper_bound=None,objective_coefficient=0):
    new_reaction_id="RGROUP"
    if group_reaction_id!=None:
       if "RGROUP" in group_reaction_id:
           new_reaction_id=group_reaction_id 
       else:
           new_reaction_id+="_"+group_reaction_id
    else:
       for reaction_id in reaction_dict:
           new_reaction_id+="_"+reaction_id
    if new_reaction_id in model.reactions:
       model.reactions.get_by_id(new_reaction_id).remove_from_model()
    new_reaction_name="Reaction Group:"
    for reaction_id in reaction_dict:
        if  reaction_dict[reaction_id]>0:
            new_reaction_name+="+"+reaction_id
        else:
            new_reaction_name+="-"+reaction_id
    metabolite = Metabolite("m"+new_reaction_id,formula='',name="mGROUP"+new_reaction_id,compartment='gr')
    group_reaction = Reaction(new_reaction_id)
    group_reaction.name = new_reaction_name
    group_reaction.subsystem = 'Reaction group'
    if upper_bound!=None:
       group_reaction.upper_bound=upper_bound 
    group_reaction.add_metabolites({metabolite:-1})
    if objective_coefficient==None:
        group_reaction.objective_coefficient=0
    model.add_reaction(group_reaction)
    group_reaction.objective_coefficient=objective_coefficient
    theoretical_lower_bound=0
    theoretical_upper_bound=0
    for reaction_id in reaction_dict:
        coef=reaction_dict[reaction_id]
        reaction=model.reactions.get_by_id(reaction_id)
        reaction.add_metabolites({metabolite:coef})
        if coef>=0:
           theoretical_upper_bound+=reaction.upper_bound
           theoretical_lower_bound+=reaction.lower_bound
        else:
           theoretical_upper_bound-=reaction.lower_bound
           theoretical_lower_bound-=reaction.upper_bound
    if lower_bound==None:
        group_reaction.lower_bound=min(round_down(theoretical_lower_bound,2),0)
    else:
        group_reaction.lower_bound=lower_bound
    if upper_bound==None:
        group_reaction.upper_bound=max(round_up(theoretical_upper_bound,2),1000)
    else:
        group_reaction.upper_bound=upper_bound


