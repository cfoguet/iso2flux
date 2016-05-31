from cobra import Model, Reaction, Metabolite
import copy
def split_model(label_model):
   size_model_dict=label_model.size_model_dict
   reaction_size_dict={}
   for reaction in label_model.emu_model.reactions:
      for emu_metabolite in reaction.metabolites:
           size=label_model.emu_dict[emu_metabolite.id]["size"]
           coef=reaction.metabolites[emu_metabolite]
           if coef>0:#Classification is based on product size
              if size not in reaction_size_dict:
                 reaction_size_dict[size]=[]
              reaction_size_dict[size].append(reaction.id)
   for size in reaction_size_dict:
       new_model=copy.deepcopy(label_model.emu_model)
       new_model.id='emu model (size'+str(size)+")"
       for dict_size in reaction_size_dict:
           if dict_size==size: #Remove all reactions except the ones with the same size
              continue
           for reaction_id in reaction_size_dict[dict_size]:
                if reaction_id in  new_model.reactions:
                   new_model.reactions.get_by_id(reaction_id).remove_from_model(remove_orphans=True) 
       size_model_dict[size]=new_model
   return size_model_dict

