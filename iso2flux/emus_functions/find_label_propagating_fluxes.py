from cobra import Model, Reaction, Metabolite
import copy
import re
from ..flux_functions.define_reaction_group import define_reaction_group
def find_label_propagating_fluxes(label_model):
  """
  Automatically groups reactions of the metabolic models and defines the "reactions_propagating_label" variable. It requires that remove_identical_reactions has been prior to being executed or it will fail to group reactions
  label_model: label_model object 
  """
  reverse_re=re.compile("(.+)_reverse$")
  group_emu_dict={}
  emu_group_dict={}
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
  print group_reaction_candidates
  for emu_reaction in group_reaction_candidates:
         print emu_reaction
         #Find the size where the reaction belongs
         for size in label_model.size_model_dict:
          if emu_reaction in label_model.size_model_dict[size].reactions:
             emu_reaction_object=label_model.size_model_dict[size].reactions.get_by_id(emu_reaction)
             break 
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
                else:
                   reaction_id+=reaction+"_"
                   local_label_reactions_dict[reaction]=1
             else:
                forward_reaction=reverse_re.match(reaction).group(1)
                if forward_reaction in label_model.reactions_propagating_label:
                   continue
                else:
                   reaction_id+=forward_reaction+"_"
                   local_label_reactions_dict[forward_reaction]=-1
         if len(local_label_reactions_dict)==1:
           reaction_id=local_label_reactions_dict.keys()[0]
           label_model.reactions_propagating_label.append(reaction_id)
         elif len(local_label_reactions_dict)>1:
           reaction_id=reaction_id[:-1]
           print [reaction_id,local_label_reactions_dict]
           if reaction_id not in label_model.reactions_propagating_label:
             label_model.label_groups_reactions_dict[reaction_id]=local_label_reactions_dict 
             label_model.reactions_propagating_label.append(reaction_id) 
             define_reaction_group(label_model.metabolic_model,local_label_reactions_dict,group_reaction_id=reaction_id,lower_bound=None,upper_bound=None,objective_coefficient=0)
           else:
            #check if it is the inverse reaction
            if local_label_reactions_dict!=label_model.label_groups_reactions_dict[reaction_id]:
               reaction_id+="_reverse"
           if reaction_id not in group_emu_dict:
               group_emu_dict[reaction_id]=[]
           if emu_reaction not in emu_group_dict:
               emu_group_dict[emu_reaction]=[]
           group_emu_dict[reaction_id].append(emu_reaction)
           emu_group_dict[emu_reaction].append(reaction_id)
           #print reaction_id+" already present"
  #label_model.temp=group_emu_dict
  
  inverse_label_groups_reactions_dict={}
  #Find the reactions that only participate in one reaction group
  for group in label_model.label_groups_reactions_dict:
      for reaction in label_model.label_groups_reactions_dict[group]:
          if reaction not in inverse_label_groups_reactions_dict:
             inverse_label_groups_reactions_dict[reaction]=[]
          inverse_label_groups_reactions_dict[reaction].append(group)
          if reaction in label_model.reactions_propagating_label:
              inverse_label_groups_reactions_dict[reaction].append(reactions)
  reactions_in_only_one_group=[]
  for reaction in inverse_label_groups_reactions_dict:
      if len(inverse_label_groups_reactions_dict[reaction])==1:
         reactions_in_only_one_group.append(reaction)
  groups_with_unique_reactions=[]
  for group in label_model.label_groups_reactions_dict:
      unique_reactions=True
      for reaction in label_model.label_groups_reactions_dict[group]:
          if reaction not in reactions_in_only_one_group:
             unique_reactions=False
      if unique_reactions:
         groups_with_unique_reactions.append(group)
  #label_model.temp=group_emu_dict
  #label_model.temp1=groups_with_unique_reactions
  for group in groups_with_unique_reactions:
      for emu_reaction in group_emu_dict[group]:
          label_model.emu_reaction_dict[emu_reaction].append(group)
          for reaction in label_model.label_groups_reactions_dict[group]:
               if reaction in label_model.reaction_merged_reactions_dict:
                 reaction=label_model.reaction_merged_reactions_dict[reaction]
               print reaction
               print label_model.emu_reaction_dict[emu_reaction]
               print emu_reaction
               if reaction in  label_model.emu_reaction_dict[emu_reaction]:
                  label_model.emu_reaction_dict[emu_reaction].remove(reaction)
               elif reaction+"_reverse" in  label_model.emu_reaction_dict[emu_reaction]:
                    label_model.emu_reaction_dict[emu_reaction].remove(reaction+"_reverse")

      for reaction in label_model.label_groups_reactions_dict[group]:
          if reaction in label_model.reaction_merged_reactions_dict:
             reaction=label_model.reaction_merged_reactions_dict[reaction]
          del label_model.reaction_emu_dict[reaction]
      label_model.reaction_emu_dict[group]=group_emu_dict[group]
      if group+"_reverse" in group_emu_dict:
        for emu_reaction in group_emu_dict[group+"_reverse"]:
          label_model.emu_reaction_dict[emu_reaction].append(group+"_reverse")
          for reaction in label_model.label_groups_reactions_dict[group]:
               if reaction in label_model.reaction_merged_reactions_dict:
                 reaction=label_model.reaction_merged_reactions_dict[reaction]
               label_model.emu_reaction_dict[emu_reaction].remove(reaction)
        for reaction in label_model.label_groups_reactions_dict[group]:
          if reaction in label_model.reaction_merged_reactions_dict:
             reaction=label_model.reaction_merged_reactions_dict[reaction]
          del label_model.reaction_emu_dict[reaction+"_reverse"]
        label_model.reaction_emu_dict[group+"_reverse"]=group_emu_dict[group+"_reverse"]
  

         
  
              
             
   
