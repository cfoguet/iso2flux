from cobra import Model, Reaction, Metabolite
#from emu_add_input_reaction import emu_add_input_reaction
#from emu_add_output_reaction import emu_add_output_reaction

def emu_add_label_ouputs_inputs(label_model,excluded_reactions_id=[]):
   """
   Automatically adds input and output reactions to the emu models based on the stochiometry of the metabolic model.
   label_model: label_model object
   excluded_reactions_id: list of strings,optional
   		IDs of the reactions that should not be added as inputs or ouputs 
   """
   simp_model=label_model.simplified_metabolic_model
   size_model_dict=label_model.size_model_dict 
   #Add outputs that are reactions for different size models
   for size in label_model.emu_size_dict:
      reactions_to_add=[]
      print("adding emu output reactions "+str(size))
      for emu in label_model.emu_size_dict[size]:
       
       met_id=label_model.emu_dict[emu]["met_id"]
       if label_model.id_isotopomer_object_dict[met_id].input==True:
          continue
       print met_id
       met=simp_model.metabolites.get_by_id(met_id)
       met_emu=size_model_dict[size].metabolites.get_by_id(emu)
       for met_reaction in met.reactions:
           found_flag=0
           coef=met_reaction.metabolites[met]
           for reaction in met_emu.reactions:
               if label_model.emu_reaction_dict[reaction.id]==met_reaction.id:
                  found_flag+=1
               elif met_reaction.id in label_model.emu_reaction_dict[reaction.id] and isinstance(label_model.emu_reaction_dict[reaction.id],list):
                  found_flag+=1
           if found_flag<abs(coef): #if this is not used it will cause error with reactions with 2 identcal substrates
              output_coef=coef+found_flag
              print(emu+" "+met_reaction.id+" "+str(output_coef))
              reactions_to_add.append({"emu":emu,"met_reaction":met_reaction.id,"coef":output_coef})
      print reactions_to_add
      for reaction_to_add in reactions_to_add:
              emu=reaction_to_add["emu"]
              met_reaction= reaction_to_add["met_reaction"]
              coef=reaction_to_add["coef"]
              model=size_model_dict[size]           
              emu_reaction=Reaction(met_reaction+"_"+emu+"_output")
              emu_reaction.name=reaction.name+" "+emu+" (output)"
              emu_reaction.subsystem = ''
              emu_reaction.add_metabolites({size_model_dict[size].metabolites.get_by_id(emu): coef})
              model.add_reaction(emu_reaction)
              if emu_reaction.id in label_model.emu_reaction_dict:
                 label_model.emu_reaction_dict[emu_reaction.id].append(met_reaction)
              else:
                 label_model.emu_reaction_dict[emu_reaction.id]=[met_reaction]
              if met_reaction in label_model.reaction_emu_dict:
                 label_model.reaction_emu_dict[met_reaction].append(emu_reaction.id)
              else:
                 label_model.reaction_emu_dict[emt_reaction.id]=[emu_reaction.id]
   #Add Inputs not included in simp_model
   added_inputs={}
   added_outputs={}
   for metabolite in label_model.irreversible_metabolic_model.metabolites:
       if metabolite.id in label_model.met_id_isotopomer_dict: #Continue if the metabolite is not defined as isotopomer
          isotopomer_object=label_model.met_id_isotopomer_dict[metabolite.id]
          if isotopomer_object.id not in label_model.metabolite_emu_dict: #Continue if the metabolite is defined as isotopomer but is not found as emu
             continue 
          if isotopomer_object.input==True:
             continue 
       else:
          continue 
       for reaction in metabolite._reaction:
           if reaction.lower_bound==0 and reaction.upper_bound==0:
              if "reflection" not in reaction.notes:
                  continue
              else:
                 reflection_reaction_id=reaction.notes["reflection"]
                 reflection_reaction=label_model.irreversible_metabolic_model.reactions.get_by_id(reflection_reaction_id)
                 if reflection_reaction.lower_bound==0.0 and reflection_reaction.upper_bound==0.0:
                    continue
           if (reaction.id in excluded_reactions_id):
              continue
           if (reaction.id in label_model.reaction_merged_reactions_dict):
              continue
           if "reflection" in reaction.notes:
              if (reaction.notes["reflection"] in label_model.reaction_merged_reactions_dict):
                 continue
           if isotopomer_object.id in  added_outputs:
              if reaction.id in added_outputs[isotopomer_object.id]:
                 continue
           if isotopomer_object.id in  added_inputs:
              if reaction.id in added_inputs[isotopomer_object.id]:
                 continue
           if (reaction.id not in simp_model.reactions): # and (isotopomer_object.input==False):
              coef=reaction.metabolites[metabolite]
              net_coef=0 #added coefficients of all the members of the same pool particpating in the reaction
              if isotopomer_object.pool==True:
                 metabolite_list=isotopomer_object.ref_met
              else:
                 metabolite_list=[metabolite]
              for met in  metabolite_list: 
                     if met in reaction.metabolites:
                        net_coef+=reaction.metabolites[met] 
              if net_coef!=0:
                 if net_coef<0:
                    reaction_type="output"
                    if isotopomer_object.id not in added_outputs:    
                         added_outputs[isotopomer_object.id]=[reaction.id]
                    else:
                         added_outputs[isotopomer_object.id].append(reaction.id) 
                 else:
                    reaction_type="input"
                    if isotopomer_object.id not in added_inputs:    
                         added_inputs[isotopomer_object.id]=[reaction.id]
                    else:
                         added_inputs[isotopomer_object.id].append(reaction.id) 
                 print("adding "+reaction.id+" as "+reaction_type+" ("+isotopomer_object.id+")") 
                 #new_reaction=Reaction(reaction.id+"_"+isotopomer_object.id+"_"+reaction_type)
                 new_reaction=Reaction(reaction.id+"_"+isotopomer_object.id+"_"+reaction_type)
                 new_reaction.name=reaction.name
                 new_reaction.subsystem = reaction.subsystem
                 emu_metabolite=label_model.simplified_metabolic_model.metabolites.get_by_id(isotopomer_object.id)
                 new_reaction.add_metabolites({emu_metabolite:net_coef})  
                 label_model.simplified_metabolic_model.add_reaction(new_reaction)
                 if emu_metabolite.id in label_model.metabolite_emu_dict:
                    emus=label_model.metabolite_emu_dict[emu_metabolite.id]
                 else:
                    continue
                 for emu in emus:
                      symmetric=False
                      if "symm_carbons" in  label_model.emu_dict[emu]:
                          if label_model.emu_dict[emu]["symm_carbons"]!=[]:
                             symmetric=True
                      size=label_model.emu_dict[emu]["size"]
                      emu_reaction=Reaction(reaction.id+"_"+emu+"_"+reaction_type)
                      emu_reaction.name=reaction.name+" ("+reaction_type+")"
                      emu_reaction.subsystem = ''
                      if not symmetric:
                             emu_reaction.add_metabolites({label_model.size_model_dict[size].metabolites.get_by_id(emu): net_coef})
                      else:
                             emu_reaction.add_metabolites({label_model.size_model_dict[size].metabolites.get_by_id(emu): 2*net_coef})
                      emu_reaction.lower_bound=0
                      label_model.size_model_dict[size].add_reaction(emu_reaction)
                      if reaction.id in label_model.reaction_emu_dict: 
                            label_model.reaction_emu_dict[reaction.id].append(emu_reaction.id)
                      else:
                            label_model.reaction_emu_dict[reaction.id]=[emu_reaction.id]      
                      label_model.emu_reaction_dict[emu_reaction.id]=reaction.id 
