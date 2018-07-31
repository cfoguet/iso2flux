import csv
from openpyxl import load_workbook
import re
import Tkinter
import tkFileDialog
from ..misc.read_spreadsheets import read_spreadsheets
import numpy

""" Names of the columns in the fluxmic files, modify them if the file format is changed"""
labelled_substrate_col_name="tracer molecule"
label_pattern_col_name="labelled positions"
lab_sub_abundance_col_name="abundance"
lab_sub_abundance_percentage=True #Set it to true if the abundance of labelled substrate are in percentatges
replicate_col_name="Parameter Value[Replicate]"
injection_col_name="Parameter Value[injection]"
metabolite_name_col_name="Metabolite name"
atomic_positions_to_the_parent_molecule_col_name="atomic positions to the parent molecule/metabolite name"
isotopologue_col_name="isotopologue"
isotopologue_fraction_col_name="isotologue abundance relative concentration"
isotopologue_fraction_abundance_percentage=True #Set it to true if the isotopologue are in percentatges
CHEBI_identifier_col_name="CHEBI identifier"
incubation_time_colname="Factor Value[Incubation time]"


tracer_regular_expression=re.compile("(.*)\[(.+)-(13C|C13)(.*)\]-(.+)")
isotopologue_regular_expression=re.compile(".*13c(.+)")
if __name__=="__main__":
   from iso2flux.misc.read_spreadsheets import read_spreadsheets

def read_metabolights(label_model,file_name,name_id_dict={},selected_condition="Ctr",selected_time=None,minimum_sd=0.01,rsm=True):
   label_model.minimum_sd=minimum_sd
   #TODO account for multiple time units
   emu0_dict={}
   label_model.experimental_dict={}
   label_model.data_name_emu_dict={}
   
   #Get name equivalency
   for x in label_model.metabolic_model.metabolites:
       if x.name==None or x.name=="":
          x.name=x.id
       if x.name.lower() not in  name_id_dict:
          name_id_dict[x.name.lower()]=[]
       name_id_dict[x.name.lower()].append(x.id)
   label_rules_dict=read_spreadsheets(file_names=label_model.label_rules_file,csv_delimiter=',',more_than_1=True,tkinter_title="Choose a experimental label file") 
   for data in label_rules_dict:
        for row in label_rules_dict[data]:
            #print row
            if row[0]==None:
               continue
            try: 
             if row[0].split(",")[0] in label_model.metabolic_model.metabolites:
               if len(row)>5:
                  if row[5] not in ["",None]:
                     metabolite_id=label_model.metabolic_model.metabolites.get_by_id(row[0].split(",")[0]).id
                     name_id_dict[str(row[5].lower())]=[metabolite_id]
            except:
                pass
   print name_id_dict
   sheet_rows_dict=read_spreadsheets(file_names=file_name,csv_delimiter=',',more_than_1=False,tkinter_title="") 
   for sheet in sheet_rows_dict:
       #Identfy the contents of each col
       col_n_dict={}
       rows=sheet_rows_dict[sheet]
       for n,header in enumerate(rows[0]):
           if header!=None:
              col_n_dict[header.lower()]=n #Remove caps
       #print col_n_dict
       #Update this if the headers are changed
       labelled_substrate_emuid_isotopologue_replicate_injection_dict={}
       #n_condition=col_n_dict["conditions"]
       n_substrate=col_n_dict[labelled_substrate_col_name.lower()]
       n_lab_pattern_substrate=col_n_dict[label_pattern_col_name.lower()]  
       #n_lab_sub_abundance=col_n_dict["abundance [value]"]
       n_lab_sub_abundance=col_n_dict[lab_sub_abundance_col_name.lower()]
       #n_lab_sub_abundance_units=col_n_dict["abundance [units]"]
       #n_time=col_n_dict["incubation time [value]"]
       n_replicate=col_n_dict[replicate_col_name.lower()]
       n_injection=col_n_dict[injection_col_name.lower()]
       n_metabolite_name=col_n_dict[metabolite_name_col_name.lower()]
       n_carbon_range=col_n_dict[atomic_positions_to_the_parent_molecule_col_name.lower()]
       n_isotopologue=col_n_dict[isotopologue_col_name.lower()]
       n_isotopologue_abundance=col_n_dict[isotopologue_fraction_col_name.lower()]
       try:
         n_CHEBI_identifier=col_n_dict[CHEBI_identifier_col_name.lower]  
         CHEBI_identifier=str(row[n_CHEBI_identifier])
       except:
         pass
       n_time=col_n_dict[incubation_time_colname.lower()]
       #n_isotopologue_units=col_n_dict["isotopologue [units]"]
       #n_isotopologue_units=col_n_dict["isotopologue [units]"]
       time_row_dict={}
       for n,row in enumerate(rows): #Select the columns matching the selected time. If no time is given the maximum time will be takem
           if n==0:
              continue
           if str(row[n_time]).rstrip().lstrip()=="":
              continue
           try:
              time=float(row[n_time])
           except:
              time=0.0
           #print time
           if time not in time_row_dict:
              time_row_dict[time]=[]
           time_row_dict[time].append(row)
           
       if selected_time==None:
          selected_time=max(time_row_dict)
          print "selected time ", selected_time
       if selected_time not in time_row_dict:
          raise Exception ("Incubation time "+str(selected_time)+" "+"is not a valid option. Valid options are:"+str(time_row_dict.keys()))              
       for n,row in enumerate(time_row_dict[selected_time]):
           #print row
           #condition=row[n_condition]
           metabolite_name=str(row[n_metabolite_name])
           unrpocessed_carbon_range=row[n_carbon_range]
           str_isotopologue=row[n_isotopologue]
           replicate=str(row[n_replicate])
           injection=str(row[n_injection])
           #print [replicate,injection]
           substrate=row[n_substrate]
           
           new_substrate=""
           for potential_substrate in substrate.split("/"):
             tracer_match=tracer_regular_expression.match(potential_substrate)
             try:
               potential_substrate=tracer_match.group(1)+tracer_match.group(5)
               tracer_pattern_temp=tracer_match.group(2) #To be used in future version
             except:
               pass 
                
               #potential substrate
               
             if potential_substrate.lower() in name_id_dict:
                new_substrate+=name_id_dict[potential_substrate.lower()][0]+"/" #Assume that so far all metabolites are the same label pool regardles of the compartment
             else: #If the name was not found try the CHEBY_ID
                try:
                  if CHEBI_identifier in name_id_dict:
                     new_substrate+=name_id_dict[potential_substrate.lower()][0]+"/"
                except:
                     if potential_substrate.rstrip()!="":
                        print "Warning: Substrate "+str(potential_substrate)+ " was not found in the constraint based model and will be ignored"
                     continue
           #print new_substrate
           substrate=new_substrate[:-1]
           abundance=row[n_lab_sub_abundance]
           pattern=row[n_lab_pattern_substrate]
           #print [substrate,n_substrate,abundance,pattern]
           
           isotopologue_abundance_str=row[n_isotopologue_abundance]
           if  any(x==None or x=="" or x==" " or x=="NA" for x in [abundance,pattern,substrate,metabolite_name,unrpocessed_carbon_range,str_isotopologue,replicate,isotopologue_abundance_str]):
               #print "continue",[abundance,pattern,substrate,metabolite_name,unrpocessed_carbon_range,str_isotopologue,replicate,isotopologue_abundance_str]
               continue
               #print "aaaaaaaaaaAAaa"
           #print "bbbbbbbb"
           #Obsolete
           #isotopologue=str(row[n_isotopologue].lower().replace("m",""))
           #print (str(row[n_isotopologue].lower())),"aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"            
           try:
              
              isotopologue=isotopologue_regular_expression.match(str(row[n_isotopologue].lower())).group(1)
           except:
             continue
           #print [n_isotopologue,row[n_isotopologue],isotopologue]
           isotopologue_abundance=max(float(row[n_isotopologue_abundance]),0)
           if isotopologue_fraction_abundance_percentage:
              isotopologue_abundance=float(isotopologue_abundance)/100.0
           labelled_substrate=str(substrate)+"$/$"+str(pattern)+"$/$"+str(abundance) 
           # labelled_substrate,"BbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbB"
           #Define Emu
           if metabolite_name.lower() in name_id_dict:
              metabolite_id=name_id_dict[metabolite_name.lower()][0] #Assume that so far all metabolites are the same label pool regardles of the compartment
           else: #If the name was not found try the CHEBY_ID
              try:
                if CHEBI_identifier in name_id_dict:
                   metabolite_id=name_id_dict[CHEBI_identifier.lower()][0]
              except:
                   if metabolite_name.rstrip()!="":
                      print "Warning: Metabolite "+str(metabolite_name)+ "was not found in the constraint based model and will be ignored"
                   continue
           carbon_range=unrpocessed_carbon_range.lower().replace("c","").split("-")
           if len(carbon_range)>1: 
                  carbons=[x for x in range(int(carbon_range[0]),int(carbon_range[1])+1)]
           else:
                  carbons=[int(carbon_range[0])]
           local_emu0_dict={}
           if metabolite_id in label_model.met_id_isotopomer_dict:
                  iso_object=label_model.met_id_isotopomer_dict[metabolite_id]
           else:
                  print (metabolite_id+" not defined as isotopomer")
           #print row
           emuid="emu_"+iso_object.id+"_"
           #print emuid
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
           if emuid not in label_model.rsm_list and rsm==True:
                        label_model.rsm_list.append(emuid) 
           if emuid not in emu0_dict:  
                  emu0_dict[emuid]= local_emu0_dict
                  if emuid not in label_model.data_name_emu_dict:
                     label_model.data_name_emu_dict[emuid]=metabolite_name+"_"+unrpocessed_carbon_range
                  if metabolite_id in label_model.rsm_metabolite_id_list:
                     label_model.rsm_list.append(emuid)
           if labelled_substrate not in labelled_substrate_emuid_isotopologue_replicate_injection_dict:
              labelled_substrate_emuid_isotopologue_replicate_injection_dict[labelled_substrate]={}
           if emuid not in labelled_substrate_emuid_isotopologue_replicate_injection_dict[labelled_substrate]:
              labelled_substrate_emuid_isotopologue_replicate_injection_dict[labelled_substrate][emuid]={}
           if isotopologue not in labelled_substrate_emuid_isotopologue_replicate_injection_dict[labelled_substrate][emuid]:
              labelled_substrate_emuid_isotopologue_replicate_injection_dict[labelled_substrate][emuid][isotopologue]={}
           if replicate not in labelled_substrate_emuid_isotopologue_replicate_injection_dict[labelled_substrate][emuid][isotopologue]:  
              labelled_substrate_emuid_isotopologue_replicate_injection_dict[labelled_substrate][emuid][isotopologue][replicate]=[]
           labelled_substrate_emuid_isotopologue_replicate_injection_dict[labelled_substrate][emuid][isotopologue][replicate].append(isotopologue_abundance)
   emus_to_remove=[]
   for labelled_substrate in labelled_substrate_emuid_isotopologue_replicate_injection_dict:
       initial_label=labelled_substrate.split("$/$")
       substrate_name=initial_label[0]
       string_pattern=initial_label[1]
       #abundance=float(initial_label[2]) 
       #print initial_label
       condition_name=(substrate_name+"_"+str(string_pattern)+"_"+str(abundance)).replace(" ","")
       substrate_id=substrate_name
       #print substrate_name
       substrate_list=substrate_id.split("/")
       #print substrate_list
       pattern_list=string_pattern.split("/")
       abundance_list=[float(x) for x in str(abundance).split("/")]
       if lab_sub_abundance_percentage:
              abundance_list=[x/100.0 for x in abundance_list]
       #print substrate_list
       for n_tracer,individual_substrate in enumerate(substrate_list):
           pattern=[int(x) for x in pattern_list[n_tracer].split(",") ]
           label_model.add_initial_label(individual_substrate,[[pattern,abundance_list[n_tracer]]],condition=condition_name,total_concentration=1)
       label_model.experimental_dict[condition_name]={}
       #Remove the measruments from the subtrate if thexy exist as they can lead to error 
       for emuid in emu0_dict.keys():
           if emu0_dict[emuid]["met_id"]in substrate_id.split("/"):
              emus_to_remove.append(emuid)
              del(emu0_dict[emuid])
          
       for emuid in labelled_substrate_emuid_isotopologue_replicate_injection_dict[labelled_substrate]:
           if emuid in emus_to_remove:
              continue
           label_model.experimental_dict[condition_name][emuid]={}
           #print labelled_substrate_emuid_isotopologue_replicate_injection_dict
           for mi in labelled_substrate_emuid_isotopologue_replicate_injection_dict[labelled_substrate][emuid]:
               if float(mi)>emu0_dict[emuid]["size"] or float(mi)<0:
                  #print mi,emu0_dict[emuid]["size"],"ccccccccccccccccccccccccccccccccccccccc"
                  continue
               data_dict=labelled_substrate_emuid_isotopologue_replicate_injection_dict[labelled_substrate][emuid][mi] 
               replicates_list=[numpy.mean(data_dict[replicates]) for replicates in data_dict]
               #print replicates_list
               #print data_dict
               mean=numpy.mean(replicates_list)
               sd=numpy.std(replicates_list)
               #sd=max(numpy.std(replicates_list),minimum_sd)
               #print [emuid,mi]
               label_model.experimental_dict[condition_name][emuid][int(mi)]={"m":mean,"sd":sd}
   #Remove the measruments from the subtrate if thexy exist as they can lead to error
   for  emuid in emus_to_remove:
        for condition in label_model.experimental_dict:
            if emuid in label_model.experimental_dict[condition]:
               del(label_model.experimental_dict[condition][emuid]) 
   label_model.emu0_dict=label_model.emu_dict=emu0_dict
   return emu0_dict,label_model.experimental_dict 
