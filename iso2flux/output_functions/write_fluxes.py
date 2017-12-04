import openpyxl 

from ..misc.write_spreadsheet import write_spreadsheet
from ..fitting.objfunc import objfunc

def export_flux_results(label_model,variables,fn="output.xlsx"):
  a,b=objfunc(label_model,variables)
  print a
  """a,chi_dict_not_rsm,simulation_not_rsm=get_objective_function(label_model,output=False,rsm="never")
  a,chi_dict_with_rsm,simulation_with_rsm=get_objective_function(label_model,output=False,rsm="always")"""
  sheet_row_data_dict={}
  row=["ID","Name","Stoichiometry","Value"]
  sheet_row_data_dict["Reaction fluxes"]=[row]
  for x in sorted(label_model.reversible_flux_dict,key=lambda v: v.upper()):  
         if "_RATIO" in x:
            continue
         if x in label_model.constrained_model.reactions:
             reaction=label_model.constrained_model.reactions.get_by_id(x)
             sheet_row_data_dict["Reaction fluxes"].append([reaction.id,reaction.name,reaction.reaction,label_model.reversible_flux_dict[x]])
  write_spreadsheet(file_name=fn,sheet_row_data_dict=sheet_row_data_dict)    

