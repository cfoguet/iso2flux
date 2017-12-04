from ..flux_functions.flux_variability_analysis import flux_variability_analysis
import openpyxl 
import math

from ..misc.write_spreadsheet import write_spreadsheet

def write_fva(model,fn="reaction_fluxes.xlsx",fva=None,fraction=1,remove0=True,change_threshold=0.001,mode="full",lp_tolerance_feasibility=1e-6,flux_precision=1e-3,reaction_list=None):
    if reaction_list==None:
       reaction_list=[x.id for x in model.reactions]
    precision=max(int(-1*(math.log10(flux_precision))),6)
    if fva==None:
       fva=flux_variability_analysis(model,fraction_of_optimum=fraction,tolerance_feasibility=lp_tolerance_feasibility)
    if mode=="full":
       row=["ID","Name","Stoichiometry","Minimum","Maximum"]
    else:
       row=["ID","Name","Flux"]
    sheet_row_data_dict={"fva":[row]}   
    for x in sorted(fva,key=lambda v: v.upper()):
         if x not in reaction_list:
            continue  
         if "_RATIO" in x:
            continue
         if abs(fva[x]["maximum"])>0.000001 or abs(fva[x]["minimum"])>0.000001 or remove0==False:
            row=[]
            reaction=model.reactions.get_by_id(x)
            row.append(reaction.id)
            row.append(reaction.name)
            minimum=fva[x]["minimum"]
            maximum=fva[x]["maximum"]
            if mode=="full":
               row.append(reaction.reaction)
               row.append(round(minimum,precision))
               row.append(round(maximum,precision))
            else:
               if (maximum-minimum)>change_threshold:
                  string=("%s < "+reaction.id+" <%s ")%( round(minimum,precision),round(maximum,precision))
               else:
                  string=str(round(maximum,precision)) 
               row.append(string)
            sheet_row_data_dict["fva"].append(row)
            #f.write(reaction.name+";"+reaction.subsystem+" ;"+reaction.reaction+";"+str(fva[x]["minimum"])+";"+str(fva[x]["maximum"])+"\n")
    print 1
    write_spreadsheet(file_name=fn,sheet_row_data_dict=sheet_row_data_dict,force_extension=True)



