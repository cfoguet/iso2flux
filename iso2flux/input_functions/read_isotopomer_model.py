from openpyxl import load_workbook
import re

from ..classes.isotopomer import isotopomer #Import the label_model class
from ..label_propagation_functions.add_label_reactions import add_label_reactions
from ..misc.read_spreadsheets import read_spreadsheets
#def read_experimental_mid(label_model,file_name,emu0_dict={},experimental_dict={},minimum_sd=0.01): 


def read_isotopomer_model(label_model,file_name,header=True):
    """
    reads the label propagation model from an xlsx file
    label_model: label_model object
    file_name: string
	name of the xlsx or CSVs files that contains the model. If its a xlsx file it should have to sheets, one containing the word "metabolite"s and and one containing the word "propagation" in their tittle. Alternatevly it can be 2 CSV files, one containing "metabolites" in the name and one containing "propagation" in the name 
		The metabolites sheet should have the folowing infromation in that order: 
			Metabolite/s (ID of the metabolites that are assumed to share label distribution separated by ",")	N carbons (Optional,number of carbons of the Metabolites/s)	Symmetry (Optional,TRUE or FALSE)	Constant(Optional,TRUE if the label distribution is not variable)
                The propagation sheet should have the following information:
			Reaction (Reaction/s id)	Substrates(Optional, ID Substrate1(c1,c2,c3,ci) +ID Substrate2(ci+1,ci+2,...) +ID SubstrateN ...)	Products (Optional, ID Product1(c1,c2,c3,ci) +ID Product2(ci+1,ci+2,...) +ID ProductN ...)
   header: bool,optional:
     indicates if the first row in the file is a header


    """
    sheet_rows_dict=read_spreadsheets(file_names=file_name,csv_delimiter=',',more_than_1=True,tkinter_title="Choose a Label Metabolites/Label propagation file(s)") 
    #wb = load_workbook(file_name, read_only=True)
    for data in sheet_rows_dict:
        print data
        if "metabolite" in data:
           for n,row in enumerate(sheet_rows_dict[data]):
               if n==0 and header:
                  continue
               if row[0]==None:
                  continue
               reference_metabolites=str(row[0].replace(" ","")).split(",")
               if row[1]!=None:
                   ncarbons=int(row[1])
               else:
                   ncarbons=-1
               if row[2]!=None:
                  if row[2]==True:
                     symmetric=True
                  elif str(row[2]).lower()=="true" or str(row[2]).lower()=="1" or str(row[2]).lower()=="yes":
                     symmetric=True
                  else:
                     symmetric=False
               else:
                    symmetric=False
               if row[3]!=None:
                  if row[3]==True:
                     label_input=True
                  elif str(row[3]).lower()=="true" or str(row[3]).lower()=="1" or str(row[3]).lower()=="yes":
                     label_input=True
                  else:
                     label_input=False
               else: 
                    label_input=False 
               #print[references_metabolites,ncarbons,symmetric]
               iso=isotopomer(reference_metabolite_id=reference_metabolites,label_model=label_model,ncarbons=ncarbons,symmetric=symmetric,label_input=label_input,iso_id=None)
               print iso.id
        elif "propagation" in data:
           label_propagation_ws=sheet_rows_dict[data]
    remove_produced_inputs(label_model) #Remove inputs that are products of irreversible reactions
    label_re=re.compile("(.+)[(](.+)[)]")
    for n_row,row in enumerate(label_propagation_ws):
               if n_row==0 and header:
                  continue
               if row[0]==None:
                  continue
               reaction_id=str(row[0].replace(" ",""))
               if "," in reaction_id:
                  reaction_id=reaction_id.split(",")
               label_id_dict={}
               substrates_dict={}
               products_dict={}
               label_propagation={}
               print reaction_id
               if row[1]!=None and row[2]!=None:
                  #Get substrates
                  substrates_list=row[1].replace(" ","").split("+")
                  #print substrates_list
                  for n_sub,substrate in enumerate(substrates_list):
                      match=label_re.match(substrate)
                      metabolite_id=match.group(1) 
                      substrates_dict["sub"+str(n_sub)]=metabolite_id
                      label_id_list=match.group(2).split(",")
                      for n_label,label_id in enumerate(label_id_list):
                          label_id_dict[label_id]=["sub"+str(n_sub),n_label]
                  products_list=row[2].replace(" ","").split("+")
                  #print products_list
                  #Get products
                  for n_prod, product in enumerate(products_list):     
                      match=label_re.match(product)
                      metabolite_id=match.group(1) 
                      products_dict["prod"+str(n_prod)]=metabolite_id
                      label_id_list=match.group(2).split(",")
                      label_propagation["prod"+str(n_prod)]=[]
                      for n_label,label_id in enumerate(label_id_list):
                          if label_id.lower()=="carb": # or label_id.lower()=="co2":
                             label_propagation["prod"+str(n_prod)].append(["carb"])
                          else:
                             substrate_id=label_id_dict[label_id][0]
                             substrate_n=label_id_dict[label_id][1]
                             label_propagation["prod"+str(n_prod)].append([substrate_id,substrate_n])
               #print [label_propagation,products_dict,substrates_dict]
               add_label_reactions(label_model,reaction_id,label_propagation=label_propagation,products_dict=products_dict,substrates_dict=substrates_dict)




def remove_produced_inputs(label_model):
    for met_id in label_model.met_id_isotopomer_dict:
        if label_model.met_id_isotopomer_dict[met_id].input==True:
           metabolite=label_model.irreversible_metabolic_model.metabolites.get_by_id(met_id)
           for reaction in metabolite.reactions:
               if "reflection" in reaction.notes:
                   continue
               coef=reaction.metabolites[metabolite]
               if coef>0.0:
                  print reaction
                  reaction.add_metabolites({metabolite:-coef})
               

#read_isotopomer_model(label_model,"isotopomer_model.xlsx",header=True)






