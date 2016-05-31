def check_mass_balance(label_model):
 for reaction in label_model.emu_model.reactions:
    mass_balance=0
    
    for emu_metabolite in reaction.metabolites:
        size=label_model.emu_dict[emu_metabolite.id]["size"]
        coef=reaction.metabolites[emu_metabolite]
        mass_balance+=size*coef
    if mass_balance!=0:
       print (reaction.id+" "+str(mass_balance))

