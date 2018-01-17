"""
This script creates a iso2flux 13C propagation model object and solves is it. It is equivalent to running create_iso2flux_model and solve_iso2flux_label. 
"""

#p_dict={'reactions_with_forced_turnover': [], 'annealing_cycle_time_limit': 8600, 'confidence_max_absolute_perturbation': 10, 'turnover_exclude_EX': True, 'annealing_n_processes': 3, 'annealing_p0': 0.4, 'identify_free_parameters_add_turnover': True, 'minimum_sd': 0.01, 'annealing_max_perturbation': 1, 'turnover_upper_bound': 20, 'confidence_perturbation': 0.1, 'annealing_m': 2000, 'annealing_n': 10, 'annealing_relative_max_sample': 0.4, 'confidence_min_absolute_perturbation': 0.05, 'annealing_pf': 0.0001, 'confidence_significance': 0.95, 'identify_free_parameters_change_threshold': 0.005, 'parameter_precision': 0.0001, 'fraction_of_optimum': 0, 'lp_tolerance_feasibility': 1e-09, 'identify_free_parameters_n_samples': 200, 'annealing_relative_min_sample': 0.25, 'annealing_iterations': 2,"gene_expression_mode":"gim3e", "gene_expression_low_expression_threshold":25,"gene_expression_high_expression_threshold":75,"gene_expression_percentile":True,"gene_expression_gene_method":"avearge", "gene_expression_gene_prefix":"","gene_expression_gene_sufix":"_AT","gene_expression_epsilon":1, "gene_expression_lex_epsilon":1e-6,"gene_expression_fraction_optimum":1, "gene_expression_absent_gene_expression_value":50}

try:
  import iso2flux
except:
 pass



import cobra
from  warnings import warn
import numpy as np
#from scipy.integrate import odeint
#import subprocess
import copy
import math
import json
import tkFileDialog
import Tkinter
import sys, getopt
import os
import time
#from cobra.core.arraybasedmodel import ArrayBasedModel
#iso2flux imports
from iso2flux.classes.label_model import Label_model #Import the label_model class
from iso2flux.label_propagation_functions.find_missing_reactions import find_missing_reactions
from iso2flux.emus_functions.solver import solver #To do, move it to emu equations

from iso2flux.output_functions.model_to_excel import model_to_excel
from iso2flux.output_functions.showr import showr
from iso2flux.output_functions.export_label_results import export_label_results
from iso2flux.output_functions.write_fva import write_fva
#from iso2flux.output_functions.print_emu_results import  print_emu_results



from iso2flux.misc.check_steady_state import check_steady_state
from iso2flux.misc.check_simulated_fractions import check_simulated_fractions
#from iso2flux.flux_functions.expression_analysis import integrate_omics_imat,integrate_omics_gim3e
from iso2flux.misc.save_load_iso2flux_model import save_iso2flux_model,load_iso2flux_model
from iso2flux.misc.round_functions import round_up, round_down 
#from iso2flux.misc.convert_model_to_irreversible import convert_to_irreversible_with_indicators

#from iso2flux.gui.gui import launch_gui

from iso2flux.input_functions.read_experimental_mid import read_experimental_mid
from iso2flux.input_functions.read_isotopomer_model import read_isotopomer_model
from iso2flux.input_functions.read_constraints import read_flux_constraints 
from iso2flux.input_functions.read_isotoflux_settings import read_isotoflux_settings
from iso2flux.input_functions.read_metabolights import read_metabolights
from iso2flux.misc.read_spreadsheets import read_spreadsheets
from iso2flux.input_functions.create_cobra_model_from_file import create_cobra_model_from_file
from iso2flux.doc.open_manual import open_manual
from iso2flux.flux_functions.check_feasibility import check_feasibility
from iso2flux.flux_functions.get_fluxes import get_fluxes


from iso2flux.flux_functions.compute_nullmatrix import compute_nullmatrix
from iso2flux.flux_functions.minimal_flux import create_minimal_fux_model,flux_minimization_fva
from iso2flux.flux_functions.build_flux_penalty_dict import build_flux_penalty_dict
from iso2flux.fitting.build_flux_solver_variables import build_flux_solver_variables
from iso2flux.emus_functions.c13solver import c13solver, flux_solver, check_bounds 
from iso2flux.fitting.get_variables import get_variables
from iso2flux.fitting.variable_sampling import variable_sampling
from iso2flux.fitting.get_variable_bounds import get_variable_bounds
from iso2flux.fitting.objfunc import objfunc
from iso2flux.fitting.build_confidence_interval_dict import build_confidence_interval_dict
from iso2flux.fitting.extract_results import extract_results
from iso2flux.fitting.get_objective_function import get_objective_function as get_xi2
from iso2flux.output_functions.write_fluxes import export_flux_results


from iso2flux.flux_functions.minimal_flux import create_minimal_fux_model,flux_minimization_fva
from iso2flux.fitting.variable_sampling import variable_sampling
from iso2flux.fitting.get_variable_bounds import get_variable_bounds
from iso2flux.fitting.objfunc import objfunc
from iso2flux.fitting.extract_results import extract_results
from iso2flux.output_functions.write_fluxes import export_flux_results
from iso2flux.output_functions.export_label_results import export_label_results
from iso2flux.fitting.optimize import optimize,flux_variation,define_isoflux_problem

#################Function used to select the fluxes to be computed for confidence intervals

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

############################

#Default Value for the parameters
model=None
constraints_file=None
mid_data_name=None
iso_model_file=None
null_matrix=None
eqn_dir=None
max_t=20
validate=False
pop_size=20
n_gen=400
number_of_processes=1
output_prefix="Iso2Flux"
max_cycles_without_improvement=8
compute_intervals=False

try:
 argv=sys.argv[1:]
 opts, args = getopt.getopt(argv,"e:l:c:s:f:o:t:q:t:vn:p:g:m:i",["experimental_data_file=","label_propagation_rules=","constraint_based_model=","settings_file=","flux_constraints=","output_name=","eqn_dir=","max_reversible_turnover=","validate","number_of_processes=","population_size=","generations_per_cycle=","max_cycles_without_improvement=","compute_confidence_intervals"])
 #opts, args = getopt.getopt(argv,"n:p:g:m:",["number_of_processes=","population_size=","generations_per_cycle=","max_cycles_without_improvement="])

except getopt.GetoptError as err:
   # print help information and exit:
   print str(err)  # will print something like "option -a not recognized":
   #print "wrong parameters"#'test.py -i <inputfile> -o <outputfile>'
   sys.exit(2)

output_name=None
for opt, arg in opts:
         print [opt,arg]
         if opt in ("--experimental_data_file","-e"):
             mid_data_name=arg
         elif opt in ("--label_propagation_rules=","-l"):
             iso_model_file=arg
             if "]" in iso_model_file:
                iso_model_file.replace("[","").replace("]","").split(",")       
         elif opt in ("--constraint_based_model=","-c"):
              if ".sbml" in arg.lower() or ".xml" in arg.lower():
                  model=cobra.io.read_sbml_model(arg)
              else:
                  model=create_cobra_model_from_file(arg) 
         elif opt in ("flux_constraints=","-f"):
              constraints_file=arg
         elif opt in ("--output_prefix=","-o"):
              output_name=arg
         elif opt in ("--eqn_dir=","-q"):
              eqn_dir=arg
         elif opt in ("--max_reversible_turnover=","-t"):
              max_t=arg
         elif opt in ("-v","--validate"):
              validate=True
         elif opt in ("--number_of_processes=","-n"):
             number_of_processes=int(arg)            
         elif opt in ("--output_prefix=","-o"):
             output_prefix=arg
         elif opt in ("--population_size=","-p"):
             pop_size=int(arg)
         elif opt in ("--generations_per_cycle=","-g"):
             n_gen=int(arg)
         elif opt in ("--max_cycles_without_improvement=","-m"):
             max_cycles_without_improvement=int(arg)
         elif opt in ("--compute_confidence_intervals","-i"):
             compute_intervals=True



if mid_data_name==None:
        raise Exception ("'--experimental_data_file' required") 
        #sys.exit(2)

if iso_model_file==None:
        iso_model_file="simple_label_model.xlsx"     

if model==None:
         model=cobra.io.read_sbml_model("simple_model.sbml")

if constraints_file!=None:
        print constraints_file
        model,ratio_dict=read_flux_constraints(model,ratio_dict={},file_name=constraints_file)


print [iso_model_file]     
     

"""
fraction_of_optimum=p_dict["fraction_of_optimum"]
reactions_with_forced_turnover=p_dict["reactions_with_forced_turnover"]
"""

#label_model=Label_model(model,lp_tolerance_feasibility=p_dict["lp_tolerance_feasibility"],parameter_precision=p_dict["parameter_precision"],reactions_with_forced_turnover=reactions_with_forced_turnover,make_all_reversible=False,moddify_false_reversible_reactions=False) #initialize the label model class   
label_model=Label_model(model,make_all_reversible=False,moddify_false_reversible_reactions=False) #initialize the label model class   
read_isotopomer_model(label_model,iso_model_file,header=True)
missing_reactions_list= find_missing_reactions(label_model).keys()
if len(missing_reactions_list)>0:
   raise Exception ("Some reactions lack label propagation rules: "+str(missing_reactions_list)) 

#emu_dict0,label_model.experimental_dict =read_experimental_mid(label_model,mid_data_name,emu0_dict={},experimental_dict={},minimum_sd=p_dict["minimum_sd"])
try:
   emu_dict0,label_model.experimental_dict =read_metabolights(label_model,mid_data_name,selected_condition=0,selected_time=0,minimum_sd=0.01,rsm=False)
except:
   if  "[" in mid_data_name and "]" in mid_data_name:
      try:
        split_name=mid_data_name.replace("[","").replace("]","")
        
        emu_dict0,label_model.experimental_dict =read_experimental_mid(label_model,split_name.split(","),emu0_dict={},experimental_dict={},minimum_sd=0.01)
      except:
        emu_dict0,label_model.experimental_dict =read_experimental_mid(label_model,split_name,emu0_dict={},experimental_dict={},minimum_sd=0.01) 
   else:
      emu_dict0,label_model.experimental_dict =read_experimental_mid(label_model,mid_data_name,emu0_dict={},experimental_dict={},minimum_sd=0.01) 


print emu_dict0


if output_name==None: 
   try:
      tk=Tkinter.Tk()
      tk.withdraw()
      output_name = tkFileDialog.asksaveasfilename(parent=tk,title="Save project as...",filetypes=[("iso2flux",".iso2flux")])
      
   except: 
     output_name=mid_data_name.replace(".csv","").replace(".xlsx","")


if ".iso2flux" in output_name:
        output_name=output_name.replace(".iso2flux","")

if eqn_dir==None:
   label_model.eqn_dir=output_name.split("/")[-1]+"_equations"
else:
   label_model.eqn_dir=eqn_dir

label_model.build(emu_dict0,force_balance=not validate,recompile_c_code=True,remove_impossible_emus=True,isotopic_steady_state=True,excluded_outputs_inputs=[],turnover_upper_bound=max_t,clear_intermediate_data=True,turnover_exclude_EX=True)


#check_steady_state(label_model,only_initial_m0=True,threshold=1e-9)#Check initial dy for steady state deviations
#check_simulated_fractions(label_model)

if null_matrix==None:
   fname=output_name+"_nullmatrix.txt" #TODO allow to name it
   compute_nullmatrix(label_model,remove_inactive_reactions=True)
   #np.savetxt(fname, label_model.flux_solver_nullmnp)
else:
   compute_nullmatrix(label_model,remove_inactive_reactions=True,fname=null_matrix) 
   #label_model.flux_solver_nullmnp=np.loadtxt("nullmatrix.txt")

build_flux_solver_variables(label_model)

get_variables(label_model)

if not validate:
   flux_penalty_dict=build_flux_penalty_dict(label_model,base_penalty=1,fn=output_name+"_flux_penalty.csv")
   label_model.best_label_variables=None
   save_iso2flux_model(label_model,name=output_name,write_sbml=True,ask_sbml_name=False,gui=False)
   #sys.exit(2)
else:
   lb_list,ub_list=get_variable_bounds(label_model,flux_penalty_dict={},max_flux=1e6)
   sampled_variables=variable_sampling(label_model,lb_list,ub_list,maximum_flux=1e6,flux_penalty_dict={},n_pop=1000,n_processes=1)
   for x in sampled_variables:
       objfunc(label_model,x)
       #check_simulated_fractions(label_model)
       passed_flag=check_steady_state(label_model,only_initial_m0=True,threshold=1e-6,precision=10,fn=output_name+"_validation.txt",fn_mode="w")
       if passed_flag==False:
          "Summary saved in "+output_name+"_validation.txt"
          break
   if passed_flag==True:
      print "Validation passed"
   sys.exit(2)


label_problem_parameters={"label_weight":1,"target_flux_dict":None,"max_chi":1e6,"max_flux":1e6,"flux_penalty_dict":{},"verbose":True,"flux_weight":0.0,"label_unfeasible_penalty":1,"flux_unfeasible_penalty":10}


iso2flux_problem=define_isoflux_problem(label_model)
optimal_solution,optimal_variables=optimize(label_model,iso2flux_problem,pop_size = pop_size,n_gen = n_gen,n_islands=number_of_processes ,max_cycles_without_improvement=max_cycles_without_improvement,stop_criteria_relative=0.005,initial_archi_x=[],lb_list=[],ub_list=[],max_flux=1e6,label_problem_parameters=label_problem_parameters)

label_model.best_chi2=optimal_solution
export_flux_results(label_model,optimal_variables,fn=output_prefix+"_fluxes.csv")
objfunc(label_model,optimal_variables)
export_label_results(label_model,fn=output_prefix+"_label.csv",show_chi=True,show_emu=False,show_fluxes=False)
np.savetxt(output_prefix+"_variables.txt",optimal_variables)
reference_parameters=[optimal_variables]
label_model.best_label_variables=optimal_variables

save_iso2flux_model(label_model,name=output_prefix+".iso2flux",write_sbml=True,ask_sbml_name=False,gui=False)
if not compute_intervals:
   output_model=label_model.constrained_model.copy()
   for reaction in output_model.reactions:
        if reaction.id in label_model.reversible_flux_dict:
           reaction.lower_bound=max(round_down(label_model.reversible_flux_dict[reaction.id],6),reaction.lower_bound)
           reaction.upper_bound=min(round_up(label_model.reversible_flux_dict[reaction.id],6),reaction.upper_bound)


if compute_intervals:
    irreversible_flag=False
    objective_tolerance=3.84
    max_chi=label_model.best_chi2+objective_tolerance
    reference_flux_group_dict,reaction_reference_flux_group_dict=find_flux_grups(label_model,reaction_list=None,irreversible_flag=irreversible_flag)
    label_problem_parameters={"label_weight":0.00000,"target_flux_dict":None,"max_chi":max_chi,"max_flux":1e6,"flux_penalty_dict":{},"verbose":True,"flux_weight":0.00000,"flux_unfeasible_penalty":10,"label_unfeasible_penalty":2}
    flux_interval_dict,irreversible_flux_interval_dict=flux_variation(label_model,iso2flux_problem,reference_flux_group_dict,reference_parameters,label_problem_parameters,max_chi=max_chi,max_flux=1e6,flux_penalty_dict=flux_penalty_dict ,pop_size=pop_size ,n_gen=n_gen ,n_islands=number_of_processes ,max_cycles_without_improvement=max_cycles_without_improvement ,stop_criteria_relative=0.005 ,max_iterations=3,log_name="confidence.txt")
    #variation_dict2=flux_variation(label_model,["rib5p_dem","fbp"],best_variables,label_problem_parameters,max_chi=max_chi,max_flux=max_flux,flux_penalty_dict=flux_penalty_dict ,pop_size=20 ,n_gen=20 ,n_islands=2 ,evolves_per_cycle=8 ,max_evolve_cycles=20 ,stop_criteria_relative=0.005 ,max_iterations=10,log_name="confidence.txt")
    write_fva(label_model.constrained_model,fn=output_prefix+"_flux_interval.csv",fva=flux_interval_dict,fraction=1,remove0=False,change_threshold=1e-6,mode="full",lp_tolerance_feasibility=1e-6,flux_precision=1e-3,reaction_list=reaction_reference_flux_group_dict)
    #Add output to a SBML model
    output_model=label_model.constrained_model.copy()
    for reaction in output_model.reactions:
        if reaction.id in flux_interval_dict:
           reaction.lower_bound=max(round_down(flux_interval_dict[reaction.id]["minimum"],6),reaction.lower_bound)
           reaction.upper_bound=min(round_up(flux_interval_dict[reaction.id]["maximum"],6),reaction.upper_bound)


cobra.io.write_sbml_model(output_model,output_prefix+"_constrained_model.xml",use_fbc_package=False)
