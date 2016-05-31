def showr(model):
    for x in model.reactions: 
        print (x.id+" "+x.reaction+" lb="+str(x.lower_bound)+" ub="+str(x.upper_bound))
