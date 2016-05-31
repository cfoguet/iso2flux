from cobra import Model, Reaction, Metabolite
def add_exchange_reaction(met_id,model,lb=0,ub=1000):
    met=model.metabolites.get_by_id(met_id)
    reaction = Reaction(met.id)
    reaction.lower_bound=lb
    reaction.upper_bound=ub
    reaction.add_metabolites({met: -1.0})
    model.add_reaction(reaction)

