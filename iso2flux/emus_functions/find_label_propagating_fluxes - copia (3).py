from cobra import Model, Reaction, Metabolite
import copy
import re
from ..flux_functions.define_reaction_group import define_reaction_group
def find_label_propagating_fluxes(label_model):
  reverse_re=re.compile("(.+)_reverse$")
  label_model.label_groups_reactions_dict={} #Add to save/load function
  label_model.reactions_propagating_label=[] #Add to save/load function
  group_reaction_candidates=[]
  for emu_reaction in label_model.emu_reaction_dict:
      associated_reactions=label_model.emu_reaction_dict[emu_reaction]
      if not isinstance (associated_reactions,list):
         associated_reactions=[associated_reactions]
      if len(associated_reactions)>1:
         group_reaction_candidates.append(emu_reaction) 
      else: 
         reaction_id=associated_reactions[0]
         if reaction_id in label_model.merged_reactions_reactions_dict:
            reaction_id=label_model.merged_reactions_reactions_dict[reaction_id][0]    
         if reverse_re.match(reaction_id)!=None:
            reaction_id=reverse_re.match(reaction_id).group(1)           
         if reaction_id not in label_model.reactions_propagating_label: 
               label_model.reactions_propagating_label.append(reaction_id)
         else:
              reaction_id +" already present"
  for emu_reaction in group_reaction_candidates:
         #print group_reaction_candidates
         associated_reactions=label_model.emu_reaction_dict[emu_reaction]
         reaction_id="LABEL_RGROUP_"
         local_label_reactions_dict={}
         for reaction in sorted(associated_reactions):
             if reaction in label_model.merged_reactions_reactions_dict:
                reaction=label_model.merged_reactions_reactions_dict[reaction][0]
             if reverse_re.match(reaction)==None:
                if reaction in label_model.reactions_propagating_label:
                   continue 
                   print "A"
                else:
                   reaction_id+=reaction+"_"
                   local_label_reactions_dict[reaction]=1
             else:
                forward_reaction=reverse_re.match(reaction).group(1)
                if forward_reaction in label_model.reactions_propagating_label:
                   continue
                   print "B"
                else:
                   reaction_id+=forward_reaction+"_"
                   local_label_reactions_dict[forward_reaction]=-1
         if len(local_label_reactions_dict)==1:
           reaction_id=local_label_reactions_dict.keys()[0]
           label_model.reactions_propagating_label.append(reaction_id)
         elif len(local_label_reactions_dict)>1:
           reaction_id=reaction_id[:-1]
           
           if reaction_id not in label_model.reactions_propagating_label:
             label_model.label_groups_reactions_dict[reaction_id]=local_label_reactions_dict 
             label_model.reactions_propagating_label.append(reaction_id) 
             define_reaction_group(label_model.metabolic_model,local_label_reactions_dict,group_reaction_id=reaction_id,lower_bound=None,upper_bound=None,objective_coefficient=0)
         else:
            print reaction_id+" already present"
  inverse_label_groups_reactions_dict={}
  for group in label_model.label_groups_reactions_dict:
      for reaction in label_model.label_groups_reactions_dict[group]:
          if reaction not in inverse_label_groups_reactions_dict:
             inverse_label_groups_reactions_dict[reaction]=[]
          inverse_label_groups_reactions_dict[reaction].append(group)
          if reaction in label_model.reactions_propagating_label:
             inverse_label_groups_reactions_dict[reaction].append(reactions)
  
              
             
   
