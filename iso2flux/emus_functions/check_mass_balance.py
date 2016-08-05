def check_mass_balance(label_model):
 """
 checks if the mass is consdereved in the reactions of the emu model network. This function is currently not used in iso2flux
 """ 
 for reaction in label_model.emu_model.reactions:
    mass_balance=0
    
    for emu_metabolite in reaction.metabolites:
        size=label_model.emu_dict[emu_metabolite.id]["size"]
        coef=reaction.metabolites[emu_metabolite]
        mass_balance+=size*coef
    if mass_balance!=0:
       print (reaction.id+" "+str(mass_balance))

