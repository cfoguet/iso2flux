import cobra
from ..misc.read_spreadsheets import read_spreadsheets
from cobra import Model, Reaction, Metabolite
def create_cobra_model_from_file(fn):
    data_rows_dict=read_spreadsheets(file_names=fn,csv_delimiter=',',more_than_1=False)
    model = Model('Constraint based model')
    candidate_metabolites=[]
    for data in data_rows_dict:
        for row in data_rows_dict[data]:
          #print row
          if len(row)>1:
            if row[1]==None:
               continue
            if "->" in row[1] or "=>" in row[1] or "<-" in row[1]: #Is a reaction:
               rid=row[0]
               new_reaction=Reaction(str(rid))
               print row[1]
               try:
                  model.add_reaction(new_reaction)
                  new_reaction.build_reaction_from_string(str(row[1]))
               except:
                  print rid+": not added"
               try:
                  new_reaction.name=str(row[2])
               except:
                   print rid+": name not added"
               try:
                  new_reaction.lower_bound=float(row[3])
               except:
                   print rid+": lb not added"
               try:
                  new_reaction.upper_bound=float(row[4])  
               except:
                  print rid+": ub not added"
               try:
                  new_reaction.objective_coefficient=float(row[5])   
               except:
                  print rid+": Objective coefficient not added"
               try:
                  new_reaction.gene_reaction_rule=str(row[6])   
               except:
                  print rid+" gene rules not added"
            else: #Probably is a metabolite add them to the list of candidate metabolites
                 candidate_metabolites.append(row)
    #Check if all candidate metabolites are efectivelly metabolites
    for row in  candidate_metabolites:
        try:
          str(row[0])
        except:
          continue
        if str(row[0]) in model.metabolites:
           metabolite=model.metabolites.get_by_id(str(row[0]))
           try:
             metabolite.name=str(row[1])
           except:
             print row[0]+": name not added"
           try:
             metabolite.compartment=str(row[3])
           except:
             print row[0]+": compartment not added"
           try:
             metabolite.formula=str(row[2])
           except:
             print row[0]+": formula not added"
     
    return model
