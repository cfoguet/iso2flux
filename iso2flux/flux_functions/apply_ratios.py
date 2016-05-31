from cobra import Model, Reaction, Metabolite
def apply_ratios(constrained_model,ratio_dict):
     for ratio in ratio_dict:
         if "RATIO_"+ratio in constrained_model.reactions:
            #print "Reaction"+" RATIO_"+ratio+ "already present"
            ratio_reaction=constrained_model.reactions.get_by_id("RATIO_"+ratio)
            for flux in ratio_dict[ratio]:
                flux_metabolite=constrained_model.metabolites.get_by_id(ratio+"_"+flux)
                coef=ratio_reaction.metabolites[flux_metabolite]
                if coef!=-ratio_dict[ratio][flux]:
                   deltacoef=-ratio_dict[ratio][flux]-coef
                   ratio_reaction.add_metabolites({flux_metabolite:deltacoef}) 
         else:
            metabolite_dict={}
            for flux in ratio_dict[ratio]:
                reaction=constrained_model.reactions.get_by_id(flux)
                metabolite=Metabolite(ratio+"_"+flux,
                   name="",
                   compartment="r") #ratio
                reaction.add_metabolites({metabolite:1}) 
                coef=ratio_dict[ratio][flux]
                metabolite_dict[metabolite]=-coef
                #print metabolite_dict
                
            ratio_reaction = Reaction("RATIO_"+ratio)
            ratio_reaction.name = "RATIO_"+ratio
            ratio_reaction.subsystem = 'Flux ratio'
            ratio_reaction.lower_bound = 0
            ratio_reaction.upper_bound = 10000.  
            ratio_reaction.add_metabolites(metabolite_dict)
            constrained_model.add_reaction(ratio_reaction)


def remove_ratio(constrained_model,ratio_id,ratio_dict):
    for flux in ratio_dict[ratio_id]:
             reaction=constrained_model.reactions.get_by_id(flux) 
             ratio_metabolite_id=ratio_id+"_"+flux
             ratio_metabolite=constrained_model.metabolites.get_by_id(ratio_metabolite_id)
             reaction.add_metabolites({ratio_metabolite:-1}) 
    ratio_reaction=constrained_model.reactions.get_by_id("RATIO_"+ratio_id)
    ratio_reaction.remove_from_model()
    del(ratio_dict[ratio_id])
 
          
