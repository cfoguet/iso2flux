import csv
from openpyxl import load_workbook
import re
import Tkinter
import tkFileDialog
from ..misc.read_spreadsheets import read_spreadsheets
def read_experimental_mid(label_model,file_names=None,emu0_dict={},experimental_dict={},minimum_sd=0.01,delimiter=','):
    """
    label_model: label_model object
    file_name: list of string
      name of the file/s containing the initial label of substrates and the measured isotopologues at the end of the experiments. If its an .xlx files diferent condition (i.e labelled Glucose, labelled glutamine) should be in different tabs, if it is a CSV file should be in separed files and those files should be loaded at the same time. If it is not provided a pop will open to select the files 
    emu0_dict: dict,optional
	dict of emus that have to be simulated but are not present in the quantfied isotopologues
    experimental_dict: dict,optional
        dict of experimental measuments not present in the file
    minimum_sd: float,optional
        minimum value of standard deviations, lower experimentak standard deviations will be set to this value
   delimiter: string,optional
        separator used if a CSV file is loaded
    """  
    label_model.minimum_sd=minimum_sd 
    label_model.data_name_emu_dict={}
    mi_re=re.compile("^m[0-9]")
    #label_pattern_re=re.compile("^([0,1][,]+)") 
    label_pattern_re=re.compile("^([0,1])")
    condition_rows_dict={}
    """if file_names==None:
      tk=Tkinter.Tk()
      tk.withdraw()
      loaded_files = tkFileDialog.askopenfiles(title='Choose a experimental label file',filetypes=[("xlsx","*.xlsx"),("csv","*.csv"),('All files','*.*')]) 
      file_names=[x.name for x in loaded_files]
      tk.destroy()
      print file_names  
    if not isinstance(file_names,list):
       file_names=[file_names]
    for file_name in file_names:
        if "xlsx" in file_name:
           wb = load_workbook(file_name, read_only=True)
           for ws in wb.worksheets:
                  condition=ws.title
                  #
                  if condition=="Sheet1" or condition=="Hoja1":
                     condition="control"
                  condition_rows_dict[condition]=[]
                  for xlsx_row in ws.rows:
                      row=[]
                      for cell in xlsx_row:
                          row.append(cell.value)
                      condition_rows_dict[condition].append(row)
        else:
           csv_file=open(file_name)
           csv_reader=csv.reader(csv_file,delimiter=delimiter)
           row_list=list(csv_reader)
           condition=(file_name.split(".")[0]).split("/")[-1] #Get the name of the file minus the extension and the path
           condition_rows_dict[condition]=row_list
           csv_file.close() 
           print row_list"""
    condition_rows_dict=read_spreadsheets(file_names=file_names,csv_delimiter=',',more_than_1=True,tkinter_title="Choose a experimental label file") 
    for condition in condition_rows_dict:
         inital_label_dict={}
         for n,row in enumerate(condition_rows_dict[condition]):
          print n
          label_pattern_flag=False
          mi_flag=False
          #Identify wether the row contains information about initial label or mid distribution
          if len(row)>2:
                if (row[1]!=None and row[1]!="") and (row[2]!=None and row[2]!=""): 
                   if label_pattern_re.match(row[1])!=None: #Is an initial label pattern
                      try: 
                        float(row[2])
                        label_pattern_flag=True
                      except:
                        label_pattern_flag=False  
                if len(row)>4:
                    if (row[3]!=None and row[3]!="") and (row[4]!=None and row[4]!=""):
                       if mi_re.match(row[3])!=None:
                          try: 
                           float(row[4])
                           mi_flag=True
                          except:
                           mi_flag=False
          if mi_flag: 
            if row[0]!=None and row[0]!="": #define emu
               metabolite_id=row[1]
               carbon_range=row[2].split("-")
               if len(carbon_range)>1: 
                  carbons=[x for x in range(int(carbon_range[0]),int(carbon_range[1])+1)]
               else:
                  carbons=[int(carbon_range[0])]
               local_emu0_dict={}
               if metabolite_id in label_model.met_id_isotopomer_dict:
                  iso_object=label_model.met_id_isotopomer_dict[metabolite_id]
               else:
                  print (metabolite_id+" not defined as isotopomer")
               emuid="emu_"+iso_object.id+"_"
               local_emu0_dict["done"]=False
               local_emu0_dict["size"]=len(carbons)
               local_emu0_dict["met_id"]=iso_object.id
               local_emu0_dict["carbons"]=carbons
               carbon_range_string=""
               for carbon in sorted(carbons):  
                   carbon_range_string+=str(carbon)
               if iso_object.symm==True:
                  #build a symetryc dic
                  symm_dict={}
                  forward_range=range(1,iso_object.n+1)
                  reverse_range=range(iso_object.n,0,-1)
                  for x in forward_range:
                      symm_dict[x]=reverse_range[x-1]
                  #print symm_dict
                  symm_carbons=[]
                  for carbon in carbons:
                      symm_carbons.append(symm_dict[carbon])
                  #print symm_carbons
                  symm_carbons=sorted(symm_carbons) 
                  local_emu0_dict["symm_carbons"]=symm_carbons
                  if symm_carbons!=carbons: #Check if they are not equal: 
                     #Identfy the lower range, which will be written first in the metabolite id
                     symm_carbon_range_string="" 
                     for carbon in symm_carbons:
                         symm_carbon_range_string+=str(carbon)
                      
                     if symm_carbons[0]<carbons[0]:
                        emuid+=symm_carbon_range_string+"_and_"+carbon_range_string
                     else:
                        emuid+=carbon_range_string+"_and_"+symm_carbon_range_string
                  else:
                     local_emu0_dict["symm_carbons"]=[] 
                     emuid+=carbon_range_string
               else:
                     emuid+=carbon_range_string
               if emuid not in emu0_dict:  
                  emu0_dict[emuid]= local_emu0_dict
               if emuid not in label_model.data_name_emu_dict:
                  label_model.data_name_emu_dict[emuid]=row[0]
               if len(row)>6:
                  if "/sm" in str(row[6]).lower() or "true" in str(row[6]).lower() or "yes" in str(row[6]).lower():
                     if emuid not in label_model.rsm_list:
                        label_model.rsm_list.append(emuid)
            if condition not in experimental_dict:
               experimental_dict[condition]={}
            if emuid not in experimental_dict[condition]:
               experimental_dict[condition][emuid]={}
            print [row[3],row[4]]
            if (row[3]!=None and row[3]!="") and (row[4]!=None,row[4]!=""):
               mi=row[3].lower()
               """if "/sm" in row[3].value.lower():
                  mi=mi.replace("/sm","")
                  if emuid not in label_model.rsm_list:
                     label_model.rsm_list.append(emuid)"""
               n_mi=int(mi.replace("m",""))
               print n_mi
               """if emuid in label_model.rsm_list and n_mi==0:
                  continue"""
               mean=float(row[4])
               if len(row)>5:
                  if row[5]==None:
                     sd=minimum_sd
                  elif float(row[5])<minimum_sd: 
                     sd=minimum_sd
                  else:
                     sd=float(row[5])
               else:
                     sd=minimum_sd
               experimental_dict[condition][emuid][n_mi]={"m":mean,"sd":sd}
          elif label_pattern_flag:
               metabolite_id=row[0]
               string_pattern=row[1].split(",") 
               pattern=[int(x) for x in string_pattern]
               abundance=float(row[2])
               if metabolite_id  not in inital_label_dict:
                  inital_label_dict[metabolite_id]=[]
               inital_label_dict[metabolite_id].append([pattern,abundance])
         for metabolite_id in inital_label_dict:
            label_pattern=inital_label_dict[metabolite_id] 
            label_model.add_initial_label(metabolite_id,label_pattern,condition=condition,total_concentration=1)
    label_model.experimental_dict=experimental_dict
    label_model.emu0_dict=label_model.emu_dict=emu0_dict
    return emu0_dict,label_model.experimental_dict 

