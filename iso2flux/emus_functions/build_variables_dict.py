def build_variables_dict(label_model,force_balance=True):
   label_model.size_variable_dict={}
   label_model.size_inverse_variable_dict={}
   for size in label_model.emu_size_dict:
       label_model.size_variable_dict[size]={}
       label_model.size_inverse_variable_dict[size]={}
       sorted_emus=sorted(label_model.emu_size_dict[size])
       n=0  
       for emu_id in sorted_emus:
           isotopomer_id=label_model.emu_dict[emu_id]["met_id"]
           isotopomer_object=label_model.id_isotopomer_object_dict[isotopomer_id]
           if isotopomer_object.input==True: #If it is an input metabolite do not add it to the variable list
              continue  
           local_emu_dict=label_model.emu_dict[emu_id]
           sorted_mi=sorted(local_emu_dict["mid"].keys())
           for mi in sorted_mi:
               mi_metabolite=local_emu_dict["mid"][mi]
               if mi_metabolite not in label_model.size_expanded_model_dict[size].metabolites:
                  print (mi_metabolite+" not present")
                  continue
               if force_balance==True and mi==0: #If we force that all mass isotopomer fractions must add up to 1 then whe can skip m0
                  continue
               label_model.size_variable_dict[size][mi_metabolite]=n
               label_model.size_inverse_variable_dict[size][n]=mi_metabolite
               n+=1
  
   return label_model.size_variable_dict, label_model.size_inverse_variable_dict
