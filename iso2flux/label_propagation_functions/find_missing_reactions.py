import copy
from get_reactants_dict import get_reactants_dict 
def find_missing_reactions(label_model):
   missing_reactions_dict={}
   model=copy.deepcopy(label_model.irreversible_metabolic_model)
   ex_reactions=model.reactions.query("EX_")
   for reaction in model.reactions:
       if reaction.id not in label_model.simplified_metabolic_model.reactions and (reaction not in ex_reactions) and (reaction.id not in label_model.reaction_merged_reactions_dict):
           if "reflection" in reaction.notes:
                 if reaction.notes["reflection"] in label_model.reaction_merged_reactions_dict:
                    continue
           substrates, products=get_reactants_dict(label_model,reaction.id)
           if reaction.lower_bound==0 and reaction.upper_bound==0:
             continue
           if substrates!={} and products!={}:
              missing_reactions_dict[reaction.id]=[substrates,products]
   for x in missing_reactions_dict:   print "missing "+x+" reaction"
   return missing_reactions_dict 

