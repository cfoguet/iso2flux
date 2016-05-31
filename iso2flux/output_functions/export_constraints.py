import openpyxl 
import Tkinter
from ..misc.write_spreadsheet import write_spreadsheet

def export_constraints(model,fn=None,ratio_dict={}):
  if fn==None:
     try:
       tk=Tkinter.Tk()
       tk.withdraw()
       fn = tkFileDialog.asksaveasfilename(parent=tk,title="Save constraints as...",filetypes=[("xlsx","*.xlsx"),("csv","*.csv")])
       if ".xlsx" not in fn and ".csv" not in fn:
          fn+=".csv"
       print fn
       tk.destroy()
     except:
       fn="constraints.xlsx"
  sheet_row_data_dict={"constraints":[["Constraints","Reaction id","Lower bound","Upper bound","Objective coefficient","Reaction name"]]}
  
  
  sorted_reaction_id=sorted([x.id for x in model.reactions],key=lambda v: v.upper())
  for reaction_id in sorted_reaction_id:
      row=[None]*5
      reaction=model.reactions.get_by_id(reaction_id)
      row[0]=reaction.id
      row[1]=reaction.lower_bound
      row[2]=reaction.upper_bound
      row[3]=reaction.objective_coefficient
      row[4]=reaction.name
      sheet_row_data_dict["constraints"].append(row)
  if ratio_dict!={}:
      sheet_row_data_dict["constraints"].append([])
      sheet_row_data_dict["constraints"].append(["Ratios"])
      for ratio in ratio_dict:
         #sheet['A'+str(n)]=ratio
         if len(ratio_dict[ratio])>2:
            raise Exception('Error: Export constraints does not suport ratios with more than 2 elements')
         ratio_value=None
         for reaction in ratio_dict[ratio]:
             if ratio_value==None:
                ratio_value=ratio_dict[ratio][reaction]
             elif ratio_dict[ratio]!=1:
                ratio_value=ratio_dict[ratio][reaction]
         sheet_row_data_dict["constraints"].append(ratio,ratio_value)
  write_spreadsheet(file_name=fn,sheet_row_data_dict=sheet_row_data_dict,sheet_order=None)    
