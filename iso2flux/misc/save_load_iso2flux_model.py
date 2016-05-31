import json
import copy
import tkFileDialog
import Tkinter 
from ..emus_functions.set_equations_variables import set_equations_variables
from ..classes.label_model import Label_model #Import the label_model class
import sys
import os
import cobra

def save_iso2flux_model(label_model,name="project",write_sbml=True,ask_sbml_name=False,gui=False):
   if gui:
      tk=Tkinter.Tk()
      tk.withdraw()
      project_name = tkFileDialog.asksaveasfilename(parent=tk,title="Save project as...",filetypes=[("iso2flux",".iso2flux")])
      if ".iso2flux" not in project_name:
        project_name+=".iso2flux"
      if write_sbml and ask_sbml_name:
         sbml_name=tkFileDialog.asksaveasfilename(parent=tk,title="Save reference SBML model as...",filetypes=[("sbml",".sbml"),("xml",".xml")])
         if ".sbml" not in sbml_name:
           sbml_name+=".sbml"
      tk.destroy()
   else:
      project_name=name+".iso2flux"
   project_dict={}
   #project_dict["condition_size_yy_dict"]=label_model.condition_size_yy_dict
   #project_dict["condition_size_yy0_dict"]=label_model.condition_size_yy0_dict
   project_dict["eqn_dir"]=label_model.eqn_dir
   project_dict["reaction_n_dict"]=label_model.reaction_n_dict
   project_dict["merged_reactions_reactions_dict"]=label_model.merged_reactions_reactions_dict
   project_dict["force_balance"]=label_model.force_balance
   project_dict["ratio_dict"]=label_model.ratio_dict
   project_dict["parameter_dict"]=label_model.parameter_dict
   project_dict["turnover_flux_dict"]=label_model.turnover_flux_dict
   project_dict["size_variable_dict"]=label_model.size_variable_dict
   project_dict["emu_dict"]=label_model.emu_dict
   project_dict["rsm_list"]=label_model.rsm_list
   project_dict["emu_size_dict"]=label_model.emu_size_dict
   project_dict["initial_label"]=label_model.initial_label
   project_dict["experimental_dict"]=label_model.experimental_dict
   project_dict["input_m0_list"]=label_model.input_m0_list
   project_dict["input_n_dict"]=label_model.input_n_dict
   project_dict["data_name_emu_dict"]=label_model.data_name_emu_dict
   project_dict["lp_tolerance_feasibility"]=label_model.lp_tolerance_feasibility
   project_dict["parameter_precision"]=label_model.parameter_precision
   project_dict["metabolite_id_isotopomer_id_dict"]=label_model.metabolite_id_isotopomer_id_dict
   project_dict["isotopomer_id_metabolite_id_dict"]=label_model.isotopomer_id_metabolite_id_dict
   project_dict["reactions_propagating_label"]=label_model.reactions_propagating_label
   project_dict["label_groups_reactions_dict"]=label_model.label_groups_reactions_dict
   if write_sbml:
      if ask_sbml_name==False or gui==False:
         sbml_name=project_name[:-9]+".sbml"
      cobra.io.write_sbml_model(label_model.metabolic_model, sbml_name)
   with open(project_name, 'w') as fp:
         json.dump(project_dict, fp)



def load_iso2flux_model(project_file="project.iso2flux",sbml_name="project_metabolic_model.sbml",ask_sbml_name=False,gui=False):
   if gui:
      tk=Tkinter.Tk()
      tk.withdraw()
      loaded_file = tkFileDialog.askopenfile(title='Select iso2flux file',filetypes=[("iso2flux",".iso2flux")]) 
      project_file=loaded_file.name
      if ask_sbml_name==True:
         loaded_file = tkFileDialog.askopenfile(title='Select reference SBML model',filetypes=[("sbml",".sbml"),("xml",".xml")]) 
         sbml_name=loaded_file.name
      tk.destroy()
   if ask_sbml_name==False:
      sbml_name=project_file[:-9]+".sbml"    
   metabolic_model=cobra.io.read_sbml_model(sbml_name)
   with open(project_file, 'r') as fp:
         project_dict=json.load(fp)
   loaded_label_model=Label_model(metabolic_model) 
   loaded_label_model.constrained_model=copy.deepcopy(self.metabolic_model)
   loaded_label_model.eqn_dir=os.getcwd()+"/equations"#project_dict["eqn_dir"]
   loaded_label_model.reaction_n_dict=project_dict["reaction_n_dict"]
   loaded_label_model.merged_reactions_reactions_dict=project_dict["merged_reactions_reactions_dict"]
   loaded_label_model.force_balance=project_dict["force_balance"]
   loaded_label_model.ratio_dict=project_dict["ratio_dict"]
   loaded_label_model.parameter_dict=project_dict["parameter_dict"]
   loaded_label_model.turnover_flux_dict=project_dict["turnover_flux_dict"]
   size_variable_dict =project_dict["size_variable_dict"]
   #convert indices back to int and generate size_inverse_variable_dict
   loaded_label_model.size_inverse_variable_dict={}
   loaded_label_model.size_variable_dict={}
   for size in size_variable_dict:
       int_size=int(size)
       loaded_label_model.size_variable_dict[int_size]=size_variable_dict[size]
       loaded_label_model.size_inverse_variable_dict[int_size]={}
       for mi in size_variable_dict[size]:
           n=size_variable_dict[size][mi]
           loaded_label_model.size_inverse_variable_dict[int_size][n]=mi
      
   emu_dict=project_dict["emu_dict"]
   loaded_label_model.emu_dict=copy.deepcopy(emu_dict)
   #convert indices back to int
   for emu in emu_dict:
       for n in emu_dict[emu]["mid"]:
           loaded_label_model.emu_dict[emu]["mid"][int(n)]=emu_dict[emu]["mid"][n]
   #
   loaded_label_model.rsm_list=project_dict["rsm_list"]
   emu_size_dict=project_dict["emu_size_dict"]
   #convert indices back to int
   loaded_label_model.emu_size_dict={}
   for size in emu_size_dict:
       int_size=int(size)
       loaded_label_model.emu_size_dict[int_size]=emu_size_dict[size]
   #
   loaded_label_model.initial_label=project_dict["initial_label"]
   experimental_dict=project_dict["experimental_dict"]
   #convert string back to int
   loaded_label_model.experimental_dict={}
   for condition in experimental_dict:
       loaded_label_model.experimental_dict[condition]={}
       for emu in experimental_dict[condition]:
           loaded_label_model.experimental_dict[condition][emu]={}
           for n in experimental_dict[condition][emu]:
               loaded_label_model.experimental_dict[condition][emu][int(n)]=experimental_dict[condition][emu][n]
   loaded_label_model.input_m0_list=project_dict["input_m0_list"]
   loaded_label_model.input_n_dict=project_dict["input_n_dict"]
   loaded_label_model.data_name_emu_dict=project_dict["data_name_emu_dict"]
   loaded_label_model.lp_tolerance_feasibility=project_dict["lp_tolerance_feasibility"]
   loaded_label_model.parameter_precision=project_dict["parameter_precision"]
   loaded_label_model.metabolite_id_isotopomer_id_dict=project_dict["metabolite_id_isotopomer_id_dict"]
   loaded_label_model.isotopomer_id_metabolite_id_dict=project_dict["isotopomer_id_metabolite_id_dict"]
   loaded_label_model.reactions_propagating_label=project_dict["reactions_propagating_label"]
   loaded_label_model.project_dict=project_dict["reactions_propagating_label"]
   loaded_label_model.project_dict=project_dict["label_groups_reactions_dict"]
   
   set_equations_variables(loaded_label_model,force_balance=loaded_label_model.force_balance)
   sys.path.insert(0, loaded_label_model.eqn_dir)
   from get_equations import get_equations
   loaded_label_model.size_emu_c_eqn_dict={}
   get_equations(loaded_label_model.size_emu_c_eqn_dict)
   del loaded_label_model.irreversible_metabolic_model #Delete irreversible model as it not used if the model is already built
   return loaded_label_model
   


