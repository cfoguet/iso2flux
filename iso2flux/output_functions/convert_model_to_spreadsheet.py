from ..misc.write_spreadsheet import write_spreadsheet



def convert_model_to_spreadsheet(model,fn):
    sheet_row_data_dict={"model":[]} 
    sheet_row_data_dict["model"].append(["Reaction id","Reaction","Reaction Name",	"Lower Bound",	"Upper Bound","Objective coefficient",	"Gene rules"])
    for reaction in model.reactions:
        sheet_row_data_dict["model"].append([reaction.id,reaction.reaction,reaction.name,	reaction.lower_bound,	reaction.upper_bound,reaction.objective_coefficient,	reaction._gene_reaction_rule]) 
    sheet_row_data_dict["model"].append([])
    sheet_row_data_dict["model"].append(["Metabolite id","Metabolite Name","Formula",	"Compartment"])
    for metabolite in model.metabolites:
        sheet_row_data_dict["model"].append([metabolite.id,metabolite.name,metabolite.formula,	metabolite.compartment]) 
    if not any(x in fn.lower() for x in ["xlsx","xlsm","xltx","xltm","csv","txt" ]): 
       fn+=".xlsx" 
    write_spreadsheet(file_name=fn,sheet_row_data_dict=sheet_row_data_dict,sheet_order=None)    




