from cobra import Model, Reaction, Metabolite

def emu_add_output_reaction(label_model,reaction,isotopomer_object_id,coef):
    if coef>0:
       print("Error: emu_add_output_reaction: positive coefficient for output reaction")
    new_reaction=Reaction(reaction.id+"_"+isotopomer_object_id+"_output")
    new_reaction.name=reaction.name
    new_reaction.subsystem = reaction.subsystem
    metabolite=label_model.simplified_metabolic_model.metabolites.get_by_id(isotopomer_object_id)
    new_reaction.add_metabolites({metabolite:coef})  
    label_model.simplified_metabolic_model.add_reaction(new_reaction)
    
    emus=label_model.metabolite_emu_dict[ metabolite.id]
    for emu in emus:
            size=label_model.emu_dict[emu]["size"]
            emu_reaction=Reaction(reaction.id+"_"+emu+"_output")
            emu_reaction.name=reaction.name+" (output)"
            emu_reaction.subsystem = ''
            emu_reaction.add_metabolites({label_model.size_model_dict[size].metabolites.get_by_id(emu): coef})
            emu_reaction.lower_bound=0
            label_model.size_model_dict[size].add_reaction(emu_reaction)
            if reaction.id in label_model.reaction_emu_dict: 
                  label_model.reaction_emu_dict[reaction.id].append(emu_reaction.id)
            else:
                  label_model.reaction_emu_dict[reaction.id]=[emu_reaction.id]      
            label_model.emu_reaction_dict[emu_reaction.id]=reaction.id 

