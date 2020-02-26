import openpyxl 

from ..misc.write_spreadsheet import write_spreadsheet
from ..fitting.objfunc import objfunc

def export_flux_results(label_model,variables,fn="output.xlsx",flux_sd_dict={},reversible=True):
  a,b=objfunc(label_model,variables)
  print a
  """a,chi_dict_not_rsm,simulation_not_rsm=get_objective_function(label_model,output=False,rsm="never")
  a,chi_dict_with_rsm,simulation_with_rsm=get_objective_function(label_model,output=False,rsm="always")"""
  sheet_row_data_dict={}
  row=["ID","Name","Stoichiometry","Value"]
  if flux_sd_dict not in (None,{}):
     row.append("local SD")
  sheet_row_data_dict["Reaction fluxes"]=[row]
  if reversible:
     for x in sorted(label_model.reversible_flux_dict,key=lambda v: v.upper()):  
         if "_RATIO" in x:
            continue
         if x in label_model.constrained_model.reactions:
             reaction=label_model.constrained_model.reactions.get_by_id(x)
             row=[reaction.id,reaction.name,reaction.reaction,label_model.reversible_flux_dict[x]]
             if x in flux_sd_dict:
                row.append(flux_sd_dict[x])   
             sheet_row_data_dict["Reaction fluxes"].append(row)
  else:
     for x in sorted(label_model.flux_dict,key=lambda v: v.upper()):  
             row=[x,label_model.flux_dict[x]]
             if x in flux_sd_dict:
                row.append(flux_sd_dict[x])   
             sheet_row_data_dict["Reaction fluxes"].append(row) 
  write_spreadsheet(file_name=fn,sheet_row_data_dict=sheet_row_data_dict)    

