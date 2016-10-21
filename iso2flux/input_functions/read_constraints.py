import re
import copy
import Tkinter
import tkFileDialog
import csv
from openpyxl import load_workbook
from ..flux_functions.apply_ratios import apply_ratios, remove_ratio
from ..flux_functions.define_reaction_group import define_reaction_group
from ..misc.round_functions import round_up,round_down
import math

def read_flux_constraints(metabolic_model,ratio_dict={},file_name=None,create_copies=True,delimiter=",",boundary_precision=1e-6):
    precision=int(-1*(math.log10(boundary_precision)))
    if file_name==None:
      tk=Tkinter.Tk()
      tk.withdraw()
      loaded_file = tkFileDialog.askopenfile(title='Choose a constraint file',filetypes=[("xlsx","*.xlsx"),("csv","*.csv"),('All files','*.*')]) 
      file_name=loaded_file.name
      tk.destroy()
      if len(file_name)==0:
         return model,ratio_dict
    if create_copies:
       model=copy.deepcopy(metabolic_model)
       ratio_dict={} 
    else:
       model=metabolic_model
    ratio_re=re.compile("^(.+)/(.+)")
    row_list=[]
    if ".xlsx"in file_name:
        wb = load_workbook(file_name, read_only=True,data_only=True)
        ws=wb.active
        for n,row_xlsx in enumerate(ws.rows):
            row=[]
            for cell in row_xlsx:   
                row.append(cell.value)
            row_list.append(row)
    else:
        csv_file=open(file_name)
        csv_reader=csv.reader(csv_file,delimiter=delimiter)
        row_list=list(csv_reader)
        csv_file.close()
        print row_list  
    for n,row in enumerate(row_list):
        row_len=len(row)
        if row_len==0:
           continue
        if row[0]!=None and row[0]!="":
           if "#" in row[0]:
               continue
           contraint_id=row[0].replace(" ","") #Remove empty spaces
           if contraint_id in model.reactions:
              reaction_id=str(contraint_id)
              reaction=model.reactions.get_by_id(reaction_id)
              if row_len>1:
                 if row[1]!=None and row[1]!="":
                    try:
                      reaction.lower_bound=round(float(row[1]),precision)
                      print reaction_id+": lower bound set to "+ str(reaction.lower_bound)
                    except:
                       print reaction_id+": lower bound not valid"
              if row_len>2:
                 if row[2]!=None and row[2]!="":
                    try:
                      reaction.upper_bound=round(float(row[2]),precision)
                      print reaction_id+": upper bound set to "+ str(reaction.upper_bound)
                    except:
                       print reaction_id+": upper bound not valid"
              if row_len>3:
                 if row[3]!=None and row[3]!="":
                    try:
                      reaction.objective_coefficient=round(float(row[3]),precision)  
                      print reaction_id+": objective coefficient set to "+ str(reaction.objective_coefficient)
                    except:
                      print reaction_id+": objective coefficient not valid"
                 #print model.optimize()
           elif ratio_re.match(contraint_id)!=None:#Check if it is a ratio
                if row_len>0:
                   if row[1]!=None and row[1]!="":
                      reaction_id1=str(ratio_re.match(contraint_id).group(1))
                      reaction_id2=str(ratio_re.match(contraint_id).group(2)) 
                      ratio_id=contraint_id
                      inverse_ratio_id=str(reaction_id2+"/"+reaction_id1)
                      if str(row[1]).lower()=="remove" or str(row[1]).lower()=="clear":
                         if ratio_id in ratio_dict:
                             remove_ratio(model,ratio_id,ratio_dict)
                             print "Ratio "+ratio_id+" removed"
                         elif inverse_ratio_id in ratio_dict:
                             remove_ratio(model,inverse_ratio_id,ratio_dict)
                             print "Ratio "+inverse_ratio_id+" removed"
                      else:   
                       try:
                        ratio_value=float(row[1])
                        if inverse_ratio_id in ratio_dict:
                           ratio_value=1.0/ratio_value
                           ratio_dict[inverse_ratio_id]={reaction_id2:ratio_value,reaction_id1:1.0}
                        else:
                           ratio_dict[ratio_id]={reaction_id1:ratio_value,reaction_id2:1.0}
                        apply_ratios(model,ratio_dict)
                        print "Ratio "+ratio_id+" added/updated" 
                       except:
                        print "Ratio "+contraint_id+" not valid"
           elif "+" in contraint_id or "-" in contraint_id:
               add_reaction_group=True
               reaction_group_dict={}
               contraint_id=contraint_id.replace("-","+-1*")
               reaction_group_elements=contraint_id.split("+")
               if len(reaction_group_elements)>1:
                  for reaction_string in reaction_group_elements:
                      if "*" in reaction_string:
                        try:
                         str_coef,reaction_id=reaction_string.split("*")
                         coef=float(str_coef)
                        except:
                         add_reaction_group=False
                      else:
                         coef=1      
                         reaction_id=reaction_string
                      if reaction_id in model.reactions:
                         reaction_group_dict[str(reaction_id)]=coef
                      else:
                         print "Reaction "+reaction_id+" not valid"
                         add_reaction_group=False
                  if add_reaction_group:
                    if row_len>1:
                       if row[1]!=None and row[1]!="":
                         try:
                           lower_bound=round(float(row[1]),precision)
                           print contraint_id+": lower bound set to "+ str(reaction.lower_bound)
                         except:
                            print contraint_id+": lower bound not valid"
                            lower_bound=None
                       else:
                          lower_bound=None  
                    else:
                          lower_bound=None
                    if row_len>2:
                       if row[2]!=None and row[2]!="":
                          try:
                            upper_bound=round((row[2]),precision)
                            print contraint_id+": upper bound set to "+ str(reaction.upper_bound)
                          except:
                             print contraint_id+": upper bound not valid"
                             upper_bound=None
                       else:
                             upper_bound=None 
                    else:
                         upper_bound=None
                    if row_len>3: 
                       if row[3]!=None and row[3]!="":
                          try:
                            objective_coefficient=round(float(row[3]),precision)  
                            print contraint_id+": objective coefficient set to "+ str(float(row[3]))
                          except:
                            print contraint_id+": objective coefficient not valid"
                            objective_coefficient=0
                       else:
                          objective_coefficient=0
                    else:
                        objective_coefficient=0
                    
                    print reaction_group_dict  
                    define_reaction_group(model,reaction_dict=reaction_group_dict,group_reaction_id=None,lower_bound=lower_bound,upper_bound=upper_bound,objective_coefficient=objective_coefficient)    
    return model,ratio_dict   

