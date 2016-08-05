from cobra import Model, Reaction, Metabolite
from get_emu_metabolite import get_emu_metabolite
import copy

def build_emu_model(label_model,emu0_dict):
  """
  A function that constructs the emu reaction newtork
  label_model: label_model object
  emu0_dict: emu dict, optional
      dictionary of the emus that must be produced by the emu network. If not provided it will take them from the variable emu0_dict from label_model 
  """
  if emu0_dict!={} or emu0_dict!=None:
     label_model.emu0_dict=emu0_dict
  else:
     label_model.emu0_dict.update(emu0_dict)
  simp_model=label_model.simplified_metabolic_model
  label_model.emu_dict=emu_dict=copy.copy(label_model.emu0_dict)
  label_model.emu_size_dict=emu_size_dict={}
  label_model.metabolite_emu_dict=metabolite_emu_dict={}
  label_model.emu_metabolite_dict=emu_metabolite_dict={}
  label_model.reaction_emu_dict=reaction_emu_dict={}
  label_model.emu_reaction_dict=emu_reaction_dict={}
  for emu in emu_dict:
      size=emu_dict[emu]["size"]
      if size not in emu_size_dict:
         emu_size_dict[size]=[emu]
      else:
         emu_size_dict[size].append(emu)
   
  emu_model=label_model.emu_model= Model('emu_model')
  still_working=True
  while still_working==True:
   still_working=False
   emu_size_dict_copy=copy.deepcopy(emu_size_dict)
   sorted_size=sorted(emu_size_dict_copy.keys(),reverse=True)
   for emu_size in sorted_size:
    emu_dict_copy=copy.deepcopy(emu_dict)   
    for emu_id in emu_dict_copy:
       emu=emu_dict_copy[emu_id]
       print emu
       if emu_dict_copy[emu_id]["done"]==True:
          continue
       else:
          still_working=True 
       met_id=emu["met_id"]
       met=simp_model.metabolites.get_by_id(met_id) 
       #print met_id
       carbons=[]
       for carbon in emu["carbons"]:
           carbons.append(carbon-1)
       carbons_list=[sorted(carbons)]
       isotopomer_object=label_model.id_isotopomer_object_dict[met_id]
       if isotopomer_object.symm==True:
          if emu["symm_carbons"]!=[]: #If no symmetry is supplied skip this (this occures when the symmetric emu is the same as the original one) 
             carbons=[]
             for carbon in emu["symm_carbons"]:
                 carbons.append(carbon-1)
             carbons_list.append(sorted(carbons))
       for carbons in carbons_list:
          #print carbons_list
          for reaction in met.reactions:
              if reaction.id not in reaction_emu_dict:
                 reaction_emu_dict[reaction.id]=[]
              #print reaction.id
              prod_coef=reaction.metabolites[met]
              if prod_coef<0: #we only want reactions that produce the metabolite 
                 continue
              possible_metabolites=[]
              if isotopomer_object.pool==True:
                 for ref_met in  isotopomer_object.ref_met:
                      possible_metabolites.append(ref_met.id)
              else:
                 possible_metabolites=[isotopomer_object.ref_met.id]
              #print possible_metabolites
              label_propagation=label_model.label_propagation_dict[reaction.id]
              label_propagation_products=label_propagation["prod_label_dict"]
              prod_emu,emu,flag=get_emu_metabolite(emu,label_model)
              for product_key in label_propagation["prod_id_dict"]:
                  metabolites_dict={prod_emu:1} #assume only one is produce for reaction. Check that that is correct
                  sub_string=""
                  prod_id=label_propagation["prod_id_dict"][product_key]
                  #print prod_id
                  if prod_id in possible_metabolites: #check if the metabolite is one of the ones associated to the isotopomer
                     propagation=label_propagation_products[product_key]
                     #print propagation
                     sub_n_dict={}
                     for n in carbons: 
                         if propagation[n]==["carb"] or propagation[n]==["carboxylation"]:
                            continue
                         sub_key=propagation[n][0]
                         sub_n=propagation[n][1]
                         sub_id=label_propagation["sub_id_dict"][sub_key]
                         iso_id=label_model.met_id_isotopomer_dict[sub_id].id
                         iso_object=label_model.id_isotopomer_object_dict[iso_id]
                         sub_symm=iso_object.symm    
                         if sub_key in sub_n_dict:
                            sub_n_dict[sub_key].append(sub_n)
                         else: 
                            sub_n_dict[sub_key]=[sub_n]
                     for sub_key in sub_n_dict:
                         sub_id=label_propagation["sub_id_dict"][sub_key]
                         iso_id=label_model.met_id_isotopomer_dict[sub_id].id
                         #sub_coef=reaction.metabolites[simp_model.metabolites.get_by_id(iso_id)]#unused 
                         n_list=sub_n_dict[sub_key]
                         temporary_dict={"met_id":iso_id,"carbons":sorted([n+1 for n in n_list]),"done":False,"emu_met":None,"size":len(n_list)}
                         #print iso_id
                         emu_msub,temporary_dict, already_present_flag=get_emu_metabolite(temporary_dict,label_model)
                         sub_string+="_"+emu_msub.id
                         if emu_msub in  metabolites_dict:
                            metabolites_dict[emu_msub]+=-1
                         else:
                            metabolites_dict[emu_msub]=-1
                         if emu_msub.id in label_model.emu0_dict: #Prevents dupplicates if during the firt iteration a metabolite define in emu0 dict is found
                            already_present_flag=True
                         if already_present_flag==False:
                            #metabolite_emu_dict={}
                            #emu_metabolite_dict[emu_msub]={}
                            size=len(n_list)
                            temporary_dict["size"]=size
                            emu_dict[emu_msub.id]=temporary_dict
                            if size in  emu_size_dict:
                               emu_size_dict[size].append(emu_msub.id)
                            else:
                               emu_size_dict[size]=[emu_msub.id] 
                         #emu_met=get_emu_metabolite(temporary_dict,simp_model,emu_model)
                     reaction_id=reaction.id+sub_string
                     if reaction_id in emu_model.reactions: 
                        reaction_id+="_bis"
                        #print reaction_id
                     new_reaction = Reaction(reaction_id)
                     new_reaction.name = ''
                     new_reaction.subsystem = ''
                     new_reaction.lower_bound = 0  # This is the default
                     new_reaction.upper_bound = 1000.  # This is the default
                     new_reaction.add_metabolites(metabolites_dict) 
                     #print [new_reaction.id,new_reaction.reaction]
                     emu_model.add_reaction(new_reaction)
                     reaction_emu_dict[reaction.id].append(new_reaction.id)
                     emu_reaction_dict[new_reaction.id]=reaction.id
       emu_dict[emu_id]["done"]=True
  label_model.emu_dict=emu_dict
  label_model.emu_size_dict=emu_size_dict
  label_model.metabolite_emu_dict=metabolite_emu_dict
  label_model.emu_metabolite_dict=emu_metabolite_dict
  label_model.reaction_emu_dict=reaction_emu_dict
  label_model.emu_reaction_dict=emu_reaction_dict
  return emu_model
