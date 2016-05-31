def remove_reflection_co2_reactions(model,output=False): #remove the link between forward and reverse reactions in reactions where CO2 or HCO3 participates
    co2_list=[]
    for metabolite in model.metabolites:
        if metabolite.formula=="CO2" or metabolite.formula=="HCO3":
           co2_list.append(metabolite)
    reversible_co2_reactions=[] 
    for co2 in co2_list:
       for reaction in co2.reactions:
           if "reflection" in reaction.notes:
              if output:
                 print("removed "+ reaction.id +" reflection")
              reversible_co2_reactions.append(reaction.id)
              del(reaction.notes["reflection"])
    return reversible_co2_reactions
       
