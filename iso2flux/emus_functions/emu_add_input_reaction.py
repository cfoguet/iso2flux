from cobra import Model, Reaction, Metabolite
def emu_add_input_reaction(label_model,reaction):
    product_coef_dict={}
    for metabolite in reaction.metabolites:
        same_pool=False
        if metabolite in label_model.metabolite_isotopomers_dict:
           coef=reaction.metabolites[metabolite]
           isotopomer_object=label_model.metabolite_isotopomers_dict[metabolite]
           pool=isotopomer_object.pool  
           if coef>0:
              if pool==True: 
                 for pool_metabolite in isotopomer_object.ref_met:
                     if pool_metabolite in reaction.metabolites and pool_metabolite!=metabolite:
                        if reaction.metabolites[pool_metabolite]<0:
                           same_pool=True
              if same_pool==False:
                 product_id=label_model.metabolite_isotopomers_dict[metabolite].id
                 product=label_model.simplified_metabolic_model.metabolites.get_by_id(product_id)
                 if product in product_coef_dict:
                    product_coef_dict[product]+=coef
                 else:
                    product_coef_dict[product]=coef      
    if len(product_coef_dict)>1:
       split_flag=True
       label_model.split_reactions_dict[reaction.id]=[]
    else:
       split_flag=False              
    for metabolite in product_coef_dict:
        coef=product_coef_dict[metabolite]
        if  metabolite.id not in label_model.metabolite_emu_dict:
            print(metabolite.id+ " not found")
            continue 
        emus=label_model.metabolite_emu_dict[ metabolite.id]
        if split_flag==True:
           new_reaction=Reaction(reaction.id+"_"+metabolite.id)
           new_reaction.name=reaction.name+" ("+metabolite.id+")"
           new_reaction.subsystem = reaction.subsystem
           new_reaction.add_metabolites({metabolite:coef})  
           label_model.simplified_metabolic_model.add_reaction(new_reaction)
           label_model.split_reactions_dict[reaction.id].append(new_reaction.id)
        else:
           new_reaction=Reaction(reaction.id)
           new_reaction.name=reaction.name
           new_reaction.subsystem = reaction.subsystem
           new_reaction.add_metabolites({metabolite:coef})  
           label_model.simplified_metabolic_model.add_reaction(new_reaction)
        for emu in emus:
            size=label_model.emu_dict[emu]["size"]
            model=label_model.size_model_dict[size]           
            emu_reaction=Reaction(reaction.id+"_"+emu+"_input")
            emu_reaction.name=reaction.name+" (input)"
            emu_reaction.subsystem = ''
            emu_reaction.add_metabolites({model.metabolites.get_by_id(emu): coef})
            emu_reaction.lower_bound=0
            model.add_reaction(emu_reaction)
            
            #print expanded_emu_reaction.reaction
            if reaction.id in label_model.reaction_emu_dict: 
               label_model.reaction_emu_dict[reaction.id].append(emu_reaction.id)
            else:
               label_model.reaction_emu_dict[reaction.id]=[emu_reaction.id]      
            label_model.emu_reaction_dict[emu_reaction.id]=reaction.id
