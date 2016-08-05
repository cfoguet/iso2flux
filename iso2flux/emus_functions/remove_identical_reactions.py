from cobra import Model, Reaction, Metabolite
import copy
import re
from ..flux_functions.define_reaction_group import define_reaction_group

def remove_identical_reactions(label_model):
  """
  identifies identical reactions in the emu_models and replaces them with a single reaction 
  label_model: label_model object 
  """
  count=0
  for size in label_model.size_model_dict:
   emu_model2=copy.deepcopy(label_model.size_model_dict[size])
   analysed_reactions=[]
   matched_reactions=[]
   for original_reaction in emu_model2.reactions:
       matched_reactions_temp=[original_reaction]  
       for original_reaction2 in emu_model2.reactions:
           if original_reaction in analysed_reactions or original_reaction==original_reaction2:
              continue   
           if original_reaction.metabolites==original_reaction2.metabolites:
              matched_reactions_temp.append(original_reaction2)
              analysed_reactions.append(original_reaction2)
              print (original_reaction.id+" "+original_reaction2.id)
              count+=1
       analysed_reactions.append(original_reaction)    
       if len(matched_reactions_temp)>1:
          matched_reactions.append(matched_reactions_temp)
   
   for reactions_sets in matched_reactions:
       metabolites_dict=reactions_sets[0].metabolites
       new_reaction_name=""
       base_reactions=[]
       for reaction in reactions_sets:
           print reaction.id
           if isinstance(label_model.emu_reaction_dict[reaction.id],list):
               base_reactions+=label_model.emu_reaction_dict[reaction.id]
           else:
               base_reactions.append(label_model.emu_reaction_dict[reaction.id])
           print base_reactions
           """if reaction==reactions_sets[len(reactions_sets)-1]: #When the last reaction is reached the id is added to the new_id without the additional "_"
              new_reaction_name+=reaction.id
              break  
           base_reaction_id=label_model.emu_reaction_dict[reaction.id]
           new_reaction_name+=base_reaction_id+"_"""
       sub_id=""
       prod_id=""
       for metabolites in metabolites_dict:
           coef=metabolites_dict[metabolites]
           if coef<0:
              sub_id+=metabolites.id+"_"
           if coef>0:
              prod_id+=metabolites.id
       if sub_id=="":
          sub_id="input_"
       if prod_id=="":
          prod_id="output"
       new_reaction_id=sub_id+"to_"+prod_id 
       new_reaction = Reaction(new_reaction_id)
       new_reaction.name =  new_reaction_name
       new_reaction.subsystem = ''
       new_reaction.lower_bound = 0  # This is the default
       new_reaction.upper_bound = 1000.  # This is the default
       new_reaction.add_metabolites(metabolites_dict) 
       emu_model2.add_reaction(new_reaction)
       for reaction in base_reactions:
           label_model.reaction_emu_dict[reaction].append(new_reaction.id)
           if new_reaction.id not in label_model.emu_reaction_dict:  
              label_model.emu_reaction_dict[new_reaction.id]=[reaction]
           else:
              label_model.emu_reaction_dict[new_reaction.id].append(reaction)
       for reaction in reactions_sets:
           base_reaction=label_model.emu_reaction_dict[reaction.id]
           if isinstance(base_reaction,list):
              for x in base_reaction:
                 label_model.reaction_emu_dict[x].remove(reaction.id) 
           else:
              label_model.reaction_emu_dict[base_reaction].remove(reaction.id)
           del label_model.emu_reaction_dict[reaction.id]  
           reaction.remove_from_model()
   label_model.size_model_dict[size]=emu_model2
   
