def set_reflections(emu_model_dict):
   print("setting reflections...")
   count=0
   if not isinstance(emu_model_dict,dict):
      emu_model_dict={1:emu_model_dict}
   for x in emu_model_dict:
      emu_model=emu_model_dict[x]
      for original_reaction in emu_model.reactions:
          for original_reaction2 in emu_model.reactions:
              if "reflection" in original_reaction.notes or "reflection" in original_reaction2.notes:
                 continue  
              reverse_metabolite_dict={}
              for metabolite in original_reaction.metabolites:
                  reverse_metabolite_dict[metabolite]=-1*original_reaction.metabolites[metabolite]
              
              if reverse_metabolite_dict==original_reaction2.metabolites:
                 if "input" in original_reaction.id or "input" in  original_reaction2.id:
                    print(original_reaction.id+" "+original_reaction2.id)
                    continue
                 original_reaction.notes["reflection"]=original_reaction2.id
                 original_reaction2.notes["reflection"]= original_reaction.id  
                 count+=1
   return count
