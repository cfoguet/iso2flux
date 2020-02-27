# -*- coding: utf-8 -*-
"""
This script works similar to the p13cmfa.py script but instead of minimizing fluxes it minimizes and maximizes each flux individually to estimate the variation interval for each flux. 
--iso2flux_model_file=,-I (mandatory):  path to the iso2flux model. 
--output_prefix=","-o= (optional): Used to define a prefix that will be added to all output files. It can be used both to name the outputs and to select the directory where the outputs will be saved (directory must already exist)
--number_of_processes=,-n(optional): Number of islands (processes) that will be used in the optimization. Should not be larger than the number of processor threads. Default is 4.  
--population_size=(optional),-p(optional):  Size of the population in each island. Default is 20. 
--generations_per_cycle=,-g(optional):  Number of generations that each island will evolve before exchanging individuals with other islands. Default is 200.
--max_cycles_without_improvement=,-m (optional): Maximum number of cycles without a significant improvement of the objective function. When this number is reached the algorithm will stop. Default is 9 
--absolute, -a (optional) If this flag is used the tolerance will be used as an absolute tolerance
--tolerance_of_label_objective=,-t(optional): Tolerance of the primary  13C MFA objective in the interval estimation. If the absolute flag is not used the maximum primary objective allowed value will be the optimal value of the primary objective plus the tolerance_of_label_objective. The optimal value of the primary objective will be taken from the iso2flux_model_file if solve_iso2flux_label.py  has been run with this model. Alternatively, 13C MFA will be run to find the optimal value of the primary objective. If the absolute flag is used, the tolerance_of_label_objective will be the absolute maximum value allowed for the primary objective. Default is 3.84
"""
import tkFileDialog
import Tkinter
import sys, getopt
import numpy as np
import copy
import os



# -*- coding: utf-8 -*-
from iso2flux import fitting
from iso2flux.misc.save_load_iso2flux_model import load_iso2flux_model, save_iso2flux_model
from iso2flux.flux_functions.minimal_flux import create_minimal_fux_model,flux_minimization_fva
from iso2flux.fitting.variable_sampling import variable_sampling
from iso2flux.fitting.get_variable_bounds import get_variable_bounds
from iso2flux.fitting.objfunc import objfunc
from iso2flux.fitting.extract_results import extract_results
from iso2flux.output_functions.write_fluxes import export_flux_results
from iso2flux.output_functions.export_label_results import export_label_results
from iso2flux.output_functions.write_fva import write_fva
from iso2flux.fitting.optimize import optimize,minimize_fluxes, flux_variation, define_isoflux_problem
from iso2flux.flux_functions.build_flux_penalty_dict import build_flux_penalty_dict, add_gene_expression_to_flux_penalty,write_flux_penalty_dict,read_flux_penalty_dict_from_file
from iso2flux.misc.read_spreadsheets import read_spreadsheets


try:
 argv=sys.argv[1:]
 
 opts, args = getopt.getopt(argv,"i:f:n:p:g:m:o:t:ar:",["iso2flux_model_file=","flux_penalty_file=","number_of_processes=","population_size=","generations_per_cycle=","max_cycles_without_improvement=","output_prefix=","tolerance_of_objective=","absolute","reaction_list_file="])
except getopt.GetoptError as err:
   # print help information and exit:
   print str(err)  # will print something like "option -a not recognized":
   #print "wrong parameters"#'test.py -i <inputfile> -o <outputfile>'
   sys.exit(2)

#file_name="hipoxia.iso2flux"
pop_size=20
n_gen=50
number_of_processes=6
output_prefix=None
flux_penalty_file=None
objective_tolerance=3.84
relative_tolerance=True
initial_flux_estimation=None
flux_penalty_dict={}
reaction_list_file=None
max_cycles_without_improvement=9
file_name=None

for opt, arg in opts:
         print [opt,arg]
         if opt in ("--iso2flux_model_file","-i"):
             file_name=arg
         elif opt in ("--number_of_processes","-n"):
             number_of_processes=int(arg)            
         elif opt in ("--output_prefix","-o"):
             output_prefix=arg
         elif opt in ("--population_size","-p"):
             pop_size=int(arg)
         elif opt in ("--generations_per_cycle","-g"):
             n_gen=int(arg)
         elif opt in ("--flux_penalty_file","-f"):
             flux_penalty_file=arg
         elif opt in ("--tolerance_of_objective","-t"):
             objective_tolerance=float(arg)
         elif opt in ("--absolute","-a"):
             relative_tolerance=False
         elif opt in ("--reaction_list_file","-r"):
             reaction_list_file=arg
         elif opt in ("max_cycles_without_improvement","-m"):
             max_cycles_without_improvement=int(arg)



print file_name

if file_name==None:
    tk=Tkinter.Tk()
    tk.withdraw()
    loaded_file = tkFileDialog.askopenfile(title='Select iso2flux file',filetypes=[("iso2flux",".iso2flux")]) 
    file_name=loaded_file.name
    tk.destroy()

print "file_name",file_name


if ".iso2flux" not in file_name:
    file_name+=".iso2flux"


if output_prefix==None:
   output_prefix=file_name.split("/")[-1].replace(".iso2flux","") 
     




label_model=load_iso2flux_model(file_name)
#lb_list,ub_list=get_variable_bounds(label_model,max_flux=1e6)

print flux_penalty_file
"""
default_name_flux_dict=file_name.replace(".iso2flux","_flux_penalty.csv")
if flux_penalty_file==None:
   flux_penalty={} 
   #Try to see if the flux_penalty dict saved with the default name exists
else:
   print "loading "+flux_penalty_file
   flux_penalty_dict=read_flux_penalty_dict_from_file(flux_penalty_file)
"""

"""

path_to_problem = os.path.dirname(fitting.__file__)+"/pygmo_problem.py"
f=open(path_to_problem,"r")
problem_script=f.read()
f.close()
exec(problem_script)
label_model.iso2flux_problem=iso2flux_problem

"""

#Get the best fit from iso2flux model, if does not exist run a simulation to get it
reference_parameters=[]

if flux_penalty_dict not in (None,{}): 
   if label_model.best_p13cmfa_variables!=None: 
      reference_parameters.append(label_model.best_p13cmfa_variables)   
      best_variables=label_model.best_p13cmfa_variables
      a,objective_dict=objfunc(label_model,label_model.best_p13cmfa_variables,flux_penalty_dict=flux_penalty_dict,flux_weight=1)
   else:
     print "Relative tolerance with flux minimization requires having done p13cmfa"
     sys.exit(2)
   if label_model.best_label_variables!=None:
      best_variables=label_model.best_label_variables
      reference_parameters.append(best_variables)
   if relative_tolerance:
      a,objective_dict=objfunc(label_model,label_model.best_p13cmfa_variables,flux_penalty_dict=flux_penalty_dict,flux_weight=1)
      #max_chi=objective_dict["chi2_score"]/objective_tolerance
      max_chi=label_model.label_tolerance
      max_flux=objective_dict["flux_score"]+objective_tolerance
   else:
      max_chi=label_model.label_tolerance
      max_flux=objective_tolerance
else:
 max_flux=1e6
 label_model.best_label_variables
 if label_model.best_label_variables!=None:
      best_variables=label_model.best_label_variables
      reference_parameters.append(best_variables)
      a,objective_dict=objfunc(label_model,best_variables,flux_penalty_dict=flux_penalty_dict,flux_weight=1)
      print objective_dict
 if label_model.best_p13cmfa_variables!=None: 
      reference_parameters.append(label_model.best_p13cmfa_variables)   
      """best_variables=label_model.best_p13cmfa_variables
      a,objective_dict=objfunc(label_model,label_model.best_p13cmfa_variables)
      print objective_dict"""
 if relative_tolerance and label_model.best_label_variables!=None:
    a,objective_dict=objfunc(label_model,best_variables,flux_penalty_dict=flux_penalty_dict,flux_weight=1)
    max_chi=objective_dict["chi2_score"]+objective_tolerance
    print max_chi
 elif relative_tolerance:
     if label_model.label_tolerance!=None:
        max_chi=label_model.label_tolerance
     else:
        print "Relative tolerance  requires having solved the model first"
        sys.exit(2)
   #label_problem_parameters={"label_weight":0.0001,"target_flux_dict":None,"max_chi":max_chi,"verbose":True}
 else:
   max_chi=objective_tolerance
   print max_chi

print max_chi,  len(reference_parameters) 

print max_chi


#label_problem_parameters={"label_weight":0.000,"target_flux_dict":None,"max_chi":max_chi,"max_flux":1e6,"flux_penalty_dict":{},"verbose":True,"flux_weight":0.00000,"flux_unfeasible_penalty":100,"label_unfeasible_penalty":100}
  
print flux_penalty_dict

def find_flux_grups(label_model,reaction_list=None,irreversible_flag=False):
    reference_flux_group_dict={}
    reaction_reference_flux_group_dict={}
    grouped_reactions=[]
    nullm=label_model.flux_solver_nullmnp
    corr_matrix=np.corrcoef(nullm)
    for n_flux,row in enumerate(corr_matrix):
        reaction1=label_model.flux_solver_n_reaction_dict[n_flux]
        if "LABEL_RGROUP_" in reaction1:
            continue
        if "RATIO_" in reaction1:
            continue
        if reaction1 in grouped_reactions:
            continue
        grouped_reactions.append(reaction1)
        reference_flux_group_dict[reaction1]=[reaction1]
        reaction_reference_flux_group_dict[reaction1]=reaction1
        for n_flux2,row2 in enumerate(corr_matrix):
            reaction2=label_model.flux_solver_n_reaction_dict[n_flux2]
            if reaction1==reaction2:
               continue
            if reaction2 in grouped_reactions:
                continue
            coef=corr_matrix[n_flux,n_flux2]
            if abs(coef)>0.99999:
               print coef,reaction1,reaction2  
               grouped_reactions.append(reaction2)                          
               reference_flux_group_dict[reaction1].append(reaction2)
               reaction_reference_flux_group_dict[reaction2]=reaction1
    if reaction_list!=None: #TODO make it work with custom list
       if irreversible_flag:
          reference_flux_group_dict=reaction_list
          reaction_reference_flux_group_dict=reaction_list
       else:
         reduced_reference_flux_group_dict={}
         reduced_reaction_reference_flux_dict={}
         select_reaction_groups=[]
         print reaction_list
         for reaction_id in reaction_list:
           select_reaction_groups.append(reaction_reference_flux_group_dict[reaction_id])
         print select_reaction_groups
         for reaction_id in set(select_reaction_groups):
             reduced_reference_flux_group_dict[reaction_id]=reference_flux_group_dict[reaction_id]
             print reduced_reference_flux_group_dict[reaction_id]
             for reaction_id2 in reduced_reference_flux_group_dict[reaction_id]:
               reduced_reaction_reference_flux_dict[reaction_id2]=reaction_id
         
         print reduced_reference_flux_group_dict
         reference_flux_group_dict=reduced_reference_flux_group_dict
         reaction_reference_flux_group_dict=reduced_reaction_reference_flux_dict
          
       
    return   reference_flux_group_dict,reaction_reference_flux_group_dict


irreversible_flag=False
if reaction_list_file!=None:
   reaction_list=[]
   sheets=read_spreadsheets(file_names=reaction_list_file)
   for sheet in sheets:
    for row in sheets[sheet]:
        for element in row:
            print element
            if "_forward" in element or "_reverse" in element:
               irreversible_flag=True
            if str(element).replace("_forward","").replace("_reverse","") in label_model.constrained_model.reactions:
               reaction_list.append(str(element)) 
   print reaction_list
else:
   reaction_list=None

reference_flux_group_dict,reaction_reference_flux_group_dict=find_flux_grups(label_model,reaction_list=reaction_list,irreversible_flag=irreversible_flag)
print reference_flux_group_dict,reaction_reference_flux_group_dict


label_problem_parameters={"label_weight":0.00000,"target_flux_dict":None,"max_chi":max_chi,"max_flux":max_flux,"flux_penalty_dict":flux_penalty_dict,"verbose":True,"flux_weight":0.00000,"flux_unfeasible_penalty":10,"label_unfeasible_penalty":2}


iso2flux_problem=define_isoflux_problem(label_model)
flux_interval_dict,irreversible_flux_interval_dict=flux_variation(label_model,iso2flux_problem,reference_flux_group_dict,reference_parameters,label_problem_parameters,max_chi=max_chi,max_flux=max_flux,flux_penalty_dict=flux_penalty_dict ,pop_size=pop_size ,n_gen=n_gen ,n_islands=number_of_processes ,max_cycles_without_improvement=max_cycles_without_improvement ,stop_criteria_relative=0.005 ,max_iterations=5,log_name="confidence.txt")

#variation_dict2=flux_variation(label_model,["rib5p_dem","fbp"],best_variables,label_problem_parameters,max_chi=max_chi,max_flux=max_flux,flux_penalty_dict=flux_penalty_dict ,pop_size=20 ,n_gen=20 ,n_islands=2 ,evolves_per_cycle=8 ,max_evolve_cycles=20 ,stop_criteria_relative=0.005 ,max_iterations=10,log_name="confidence.txt")

write_fva(label_model.constrained_model,fn=output_prefix+"_flux_interval.csv",fva=flux_interval_dict,fraction=1,remove0=False,change_threshold=1e-6,mode="full",lp_tolerance_feasibility=1e-6,flux_precision=1e-3,reaction_list=reaction_reference_flux_group_dict)


