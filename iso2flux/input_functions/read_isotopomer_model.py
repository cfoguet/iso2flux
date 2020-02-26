from openpyxl import load_workbook
import re
import cobra
from ..classes.isotopomer import isotopomer #Import the label_model class
from ..label_propagation_functions.add_label_reactions import add_label_reactions
from ..misc.read_spreadsheets import read_spreadsheets
from ..label_propagation_functions.find_missing_reactions import find_missing_reactions
#def read_experimental_mid(label_model,file_name,emu0_dict={},experimental_dict={},minimum_sd=0.01): 


def read_isotopomer_model(label_model,file_name,header=True):
    """
    reads the label propagation model from an xlsx file
    label_model: label_model object
    file_name: string
	name of the xlsx or CSVs files that contains the model. It should define the metabolites that can labelled and the label propagation rules of reactions. 
		The metabolites  should have the folowing infromation in that order: 
			Metabolite/s (ID of the metabolites that are assumed to share label distribution separated by ",")	N carbons (Optional,number of carbons of the Metabolites/s)	Symmetry (Optional,TRUE or FALSE)	Constant(Optional,TRUE if the label distribution is not variable)
                The propagation sheet should have the following information:
			Reaction (Reaction/s id)	Substrates(Optional, ID Substrate1(c1,c2,c3,ci) +ID Substrate2(ci+1,ci+2,...) +ID SubstrateN ...)	Products (Optional, ID Product1(c1,c2,c3,ci) +ID Product2(ci+1,ci+2,...) +ID ProductN ...)
   header: bool,deprectaed


    """
    label_model.rsm_metabolite_id_list=[]
    label_model.label_rules_file=file_name
    sheet_rows_dict=read_spreadsheets(file_names=file_name,csv_delimiter=',',more_than_1=False,tkinter_title="Choose label propagation rules") 
    #wb = load_workbook(file_name, read_only=True)
    metabolites_rows=[]
    reactions_rows=[]
    for data in sheet_rows_dict:
        for row in sheet_rows_dict[data]:
           try:
            if row[0].split(",")[0] in label_model.metabolic_model.metabolites: #Check if you are looking at the list of metabolites
               metabolites_rows.append(row)
               print row
            elif row[0].split(",")[0] in label_model.metabolic_model.reactions:
               reactions_rows.append(row)
           except:
               pass 
    for n,row in enumerate(metabolites_rows):
               print row
               if row[0]==None:
                  continue
               reference_metabolites=str(row[0].replace(" ","")).split(",")
               if len(row)>1:
                  if row[1]!=None:
                      try:
                         ncarbons=int(row[1])
                      except:
                         ncarbons=-1
                  else:
                      ncarbons=-1
               if len(row)>2:
                 if row[2]!=None:
                   if row[2]==True:
                      symmetric=True
                   elif str(row[2]).lower()=="true" or str(row[2]).lower()=="1" or str(row[2]).lower()=="yes":
                      symmetric=True
                   else:
                      symmetric=False
                 else:
                    symmetric=False
               if len(row)>3:
                if row[3]!=None:
                   if row[3]==True:
                      label_input=True
                   elif str(row[3]).lower()=="true" or str(row[3]).lower()=="1" or str(row[3]).lower()=="yes":
                      label_input=True
                   else:
                      label_input=False
                else:
                   label_input=False
               else: 
                     label_input=False
               if len(row)>4:
                  if row[4]!=None:
                     if "/sm" in str(row[4]).lower():
                        for metabolite_id in reference_metabolites:
                            label_model.rsm_metabolite_id_list.append(metabolite_id) 
                     elif row[4]==True or str(row[4]).lower()=="true" or str(row[4]).lower()=="1" or str(row[4]).lower()=="yes": #Does it have a large unlabelled pool?
                        for metabolite_id in reference_metabolites:
                            EX_present=False
                            EX_reaction_id="EX_"+metabolite_id
                            if EX_reaction_id in label_model.metabolic_model.reactions:
                                label_model.reactions_with_forced_turnover.append(EX_reaction_id)
                                EX_present=True
                                EX_reaction=label_model.irreversible_metabolic_model.reactions.get_by_id(EX_reaction_id)
                                if "reflection" not in EX_reaction.notes:
                                    EX_reaction.lower_bound=-1000
                                    EX_reaction.upper_bound=1000
                                break
                        if  not EX_present:
                            #Add it to the irreversible_model
                            metabolite=label_model.irreversible_metabolic_model.metabolites.get_by_id(reference_metabolites[0])
                            reaction = cobra.core.Reaction('EX_'+metabolite.id)
                            reaction.name = 'Exchange of '+metabolite.id
                            reaction.subsystem = 'Exchange reaction'
                            reaction.lower_bound=-1000
                            reaction.upper_bound=1000
                            reaction.add_metabolites({metabolite:-1})
                            label_model.irreversible_metabolic_model.add_reaction(reaction)
                            #Add it to the metabolic model
                            metabolite=label_model.metabolic_model.metabolites.get_by_id(reference_metabolites[0])
                            reaction = cobra.core.Reaction('EX_'+metabolite.id)
                            reaction.name = 'Exchange of '+metabolite.id
                            reaction.subsystem = 'Exchange reaction'
                            reaction.lower_bound=0
                            reaction.upper_bound=0
                            reaction.add_metabolites({metabolite:-1})
                            label_model.metabolic_model.add_reaction(reaction)
                            
                            label_model.reactions_with_forced_turnover.append(reaction.id) 
               
               iso=isotopomer(reference_metabolite_id=reference_metabolites,label_model=label_model,ncarbons=ncarbons,symmetric=symmetric,label_input=label_input,iso_id=None)
               print iso.id      
    cobra.manipulation.convert_to_irreversible(label_model.irreversible_metabolic_model)#Convert any Exchange we migh have added to irreversible
    remove_produced_inputs(label_model) #Remove inputs that are products of irreversible reactions
    label_re=re.compile("(.+)[(](.+)[)]")
    for n_row,row in enumerate(reactions_rows):
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
               if len(row)>2:
                 if (row[1] not in (None,""," ")) and (row[2] not in (None,""," ")):
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
               print [reaction_id,substrates_dict,substrates_dict]
               if len(row)>3:
                  if "turn" in str(row[3]).lower() and not isinstance(reaction_id,list):
                      label_model.reactions_with_forced_turnover.append(reaction_id)
                      reaction=label_model.irreversible_metabolic_model.reactions.get_by_id(reaction_id)
                      if reaction.upper_bound>0:
                         reaction.lower_bound=-reaction.upper_bound 
                      else:
                         reaction.lower_bound=-1000
               cobra.manipulation.convert_to_irreversible(label_model.irreversible_metabolic_model)       
               add_label_reactions(label_model,reaction_id,label_propagation=label_propagation,products_dict=products_dict,substrates_dict=substrates_dict)
    
    add_missing_uni_uni_reactions(label_model)



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

def add_missing_uni_uni_reactions(label_model):
    missing_reactions_dict=find_missing_reactions(label_model,verbose=False)
    print missing_reactions_dict
    for missing_reaction in missing_reactions_dict:
        substrates=missing_reactions_dict[missing_reaction][0]
        products=missing_reactions_dict[missing_reaction][1]
        if len(products)==1 and len(substrates)==1: #Identfy uni-uni reactions
           substrate_object= label_model.met_id_isotopomer_dict[substrates[0]]
           product_object= label_model.met_id_isotopomer_dict[products[0]]
           if product_object.n==substrate_object.n: #Check if they have the same number of carbons
              add_label_reactions(label_model,missing_reaction)
    #Check if they are any reactions still missing              
#read_isotopomer_model(label_model,"isotopomer_model.xlsx",header=True)






