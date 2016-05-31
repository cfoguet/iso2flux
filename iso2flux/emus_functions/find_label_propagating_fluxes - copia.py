from cobra import Model, Reaction, Metabolite
import copy
import re
from ..flux_functions.define_reaction_group import define_reaction_group
from set_reflections import set_reflections

def find_label_propagating_fluxes(label_model):
  set_reflections(label_model.size_model_dict)     
  reverse_re=re.compile("(.+)_reverse$")
  label_model.label_groups_reactions_dict={} #Add to save/load function
  label_model.reactions_propagating_label=[] #Add to save/load function
  done_reflections=[]
  for emu_reaction in label_model.emu_reaction_dict:
      if emu_reaction in done_reflections:
         continue
      #Find the size where the reaction belongs
      for size in label_model.size_model_dict:
          if emu_reaction in label_model.size_model_dict[size].reactions:
             emu_reaction_object=label_model.size_model_dict[size].reactions.get_by_id(emu_reaction)
             break 
          
      associated_reactions=label_model.emu_reaction_dict[emu_reaction]
      if not isinstance (associated_reactions,list):
         associated_reactions=[associated_reactions]
      if len(associated_reactions)>1:
         reaction_id="LABEL_RGROUP_"
         local_label_reactions_dict={}
         for reaction in sorted(associated_reactions):
             label_model.emu_reaction_dict[emu_reaction].remove(reaction)
             label_model.reaction_emu_dict[reaction].remove(emu_reaction)
             if reaction in label_model.merged_reactions_reactions_dict:
                reaction=label_model.merged_reactions_reactions_dict[reaction][0]
             if reverse_re.match(reaction)==None:
                reaction_id+=reaction+"_"
                local_label_reactions_dict[reaction]=1
             else:
                forward_reaction=reverse_re.match(reaction).group(1)
                reaction_id+=forward_reaction+"_"
                local_label_reactions_dict[forward_reaction]=-1
         group_reaction_id=reaction_id=reaction_id[:-1]
         if group_reaction_id not in label_model.emu_reaction_dict[emu_reaction]:
            label_model.emu_reaction_dict[emu_reaction].append(group_reaction_id)
         if group_reaction_id not in label_model.reaction_emu_dict:
            label_model.reaction_emu_dict[group_reaction_id]=[]
         label_model.reaction_emu_dict[group_reaction_id].append(emu_reaction)
         if reaction_id not in label_model.reactions_propagating_label:
            label_model.label_groups_reactions_dict[reaction_id]=local_label_reactions_dict 
            label_model.reactions_propagating_label.append(reaction_id) 
            define_reaction_group(label_model.metabolic_model,local_label_reactions_dict,group_reaction_id=reaction_id,lower_bound=None,upper_bound=None,objective_coefficient=0)
         """else:
            print reaction_id+" already present"""
         if "reflection" in emu_reaction_object.notes:
           emu_reflection=emu_reaction_object.notes["reflection"]
           print ["A1",emu_reflection]
           reflection_associated_reactions=label_model.emu_reaction_dict[emu_reflection]
           print ["A2",reflection_associated_reactions]
           done_reflections.append(emu_reflection)
           if not isinstance (reflection_associated_reactions,list):
               reflection_associated_reactions=[reflection_associated_reactions]
           if len(reflection_associated_reactions)>1:
            print ["B"]
            for reaction in reflection_associated_reactions:
                add_group_reverse=False
                if reaction in label_model.merged_reactions_reactions_dict:
                   reaction=label_model.merged_reactions_reactions_dict[reaction][0]
                if reaction in local_label_reactions_dict:
                      add_group_reverse=True
                elif reverse_re.match(reaction)!=None:
                      if reverse_re.match(reaction).group(1) in local_label_reactions_dict:
                           add_group_reverse=True
                if add_group_reverse:
                      label_model.emu_reaction_dict[emu_reflection].remove(reaction)
                      label_model.reaction_emu_dict[reaction].remove(emu_reflection)
                      reverse_group_reaction_id=group_reaction_id+"_reverse"
                      if reverse_group_reaction_id not in label_model.emu_reaction_dict[emu_reflection]:
                         label_model.emu_reaction_dict[emu_reflection].append(reverse_group_reaction_id)
                      if reverse_group_reaction_id not in label_model.reaction_emu_dict:
                        label_model.reaction_emu_dict[reverse_group_reaction_id]=[]
                      label_model.reaction_emu_dict[reverse_group_reaction_id].append(emu_reflection)
                      
            
      else: 
         reaction_id=associated_reactions[0]
         if reaction_id in label_model.merged_reactions_reactions_dict:
            reaction_id=label_model.merged_reactions_reactions_dict[reaction_id][0]    
         if reverse_re.match(reaction_id)!=None:
            reaction_id=reverse_re.match(reaction_id).group(1)           
         if reaction_id not in label_model.reactions_propagating_label: 
               label_model.reactions_propagating_label.append(reaction_id)
         """else:
              reaction_id +" already present"""
  for reaction in label_model.reaction_emu_dict.keys(): #Delete reactions that no longuer have any emu associated
      if len(label_model.reaction_emu_dict[reaction])==0:
         del label_model.reaction_emu_dict[reaction] 
