"""
This script creates a iso2flux 13C propagation model object. This object is used as input for the other scripts. It takes the following command line arguments: 
--experimental_data_file=, -e (Mandatory): Path to the file describing the labelled substrates used and the quantified isotopologue fractions in metabolic products. See here. 
--label_propagation_rules=,-l (Mandatory) : Path to the file describing the label propagation rules. See here  
--constraint_based_model=,-c (Mandatory) Path to the constraint based model that will be used in the analysis. See here 
--flux_constraints=,-f (Optional):Path the file defining additional constraints  for the  constraint_based_model. See here
--max_reversible_turnover=,-t (Optional) :  Maximum reversible reaction turnover allowed. Reversible reaction turnover is defined as the flux that is common to the forward (Jf) and reverse (Jr) reactions in a reversible reaction. Default is 20.  
--output_prefix=,-o (Optional): Used to define a prefix that will be added to all output files. It can be used both to name the outputs and to select the directory where the outputs will be saved (directory must already exist). If left empty the name of the experimental_data_file will be used. 
--eqn_dir=,-q (Optional): Used to define the directory where iso2flux equations will be saved. 
 -v,--validate (Optional): If this flag is used iso2flux will validate if the label propagation rules are properly coupled to the constraint based model. If validation fails it is usually a sign that there is an error in the label propagation rules file.  This mode will not create normal iso2flux inputs. To create a iso2flux model run the script without this flag. 
The main outputs (when not run with the -v flag) will be an iso2flux model (.iso2flux) file and a set of files associated to it.  As a general rule, the files generated should not be manually edited or altered. The exception is the flux_penalty.csv file which can be used to set the flux minimization weights used by the p13cmfa script. 

"""

#p_dict={'reactions_with_forced_turnover': [], 'annealing_cycle_time_limit': 8600, 'confidence_max_absolute_perturbation': 10, 'turnover_exclude_EX': True, 'annealing_n_processes': 3, 'annealing_p0': 0.4, 'identify_free_parameters_add_turnover': True, 'minimum_sd': 0.01, 'annealing_max_perturbation': 1, 'turnover_upper_bound': 20, 'confidence_perturbation': 0.1, 'annealing_m': 2000, 'annealing_n': 10, 'annealing_relative_max_sample': 0.4, 'confidence_min_absolute_perturbation': 0.05, 'annealing_pf': 0.0001, 'confidence_significance': 0.95, 'identify_free_parameters_change_threshold': 0.005, 'parameter_precision': 0.0001, 'fraction_of_optimum': 0, 'lp_tolerance_feasibility': 1e-09, 'identify_free_parameters_n_samples': 200, 'annealing_relative_min_sample': 0.25, 'annealing_iterations': 2,"gene_expression_mode":"gim3e", "gene_expression_low_expression_threshold":25,"gene_expression_high_expression_threshold":75,"gene_expression_percentile":True,"gene_expression_gene_method":"avearge", "gene_expression_gene_prefix":"","gene_expression_gene_sufix":"_AT","gene_expression_epsilon":1, "gene_expression_lex_epsilon":1e-6,"gene_expression_fraction_optimum":1, "gene_expression_absent_gene_expression_value":50}

try:
  import iso2flux
except:
 pass


import sympy
import numpy
import cobra
import random
from cobra import Model, Reaction, Metabolite
from cobra.manipulation.modify import convert_to_irreversible
from iso2flux.flux_functions.flux_variability_analysis import flux_variability_analysis
from iso2flux.flux_functions.minimal_flux import add_flux_limit_constraints
from  warnings import warn
import numpy as np
#from scipy.integrate import odeint
#import subprocess
from timeit import timeit
import copy
import math
import json
import tkFileDialog
import Tkinter
import sys, getopt
import os

from scipy.optimize import broyden1
from scipy.optimize import broyden2
from scipy.optimize import fsolve
from scipy.integrate import odeint
from scipy.optimize import broyden1
from scipy.optimize import broyden2
from scipy.optimize import fsolve
from scipy.integrate import odeint
import json 
import cobra
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

from iso2flux.fitting.anneal import annealing
from iso2flux.fitting.get_objective_function import get_objective_function
from iso2flux.fitting.apply_parameters import apply_parameters
from iso2flux.fitting.clear_parameters import clear_parameters 
from iso2flux.fitting.save_parameters import save_parameters
from iso2flux.fitting.load_parameters import load_parameters
from iso2flux.fitting.identify_free_parameters import identify_free_parameters
from iso2flux.fitting.confidence_intervals import estimate_confidence_intervals
from iso2flux.fitting.sampling import sampling
from iso2flux.fitting.clear_parameters import clear_parameters

from iso2flux.misc.check_steady_state import check_steady_state
from iso2flux.misc.check_simulated_fractions import check_simulated_fractions
#from iso2flux.flux_functions.expression_analysis import integrate_omics_imat,integrate_omics_gim3e
from iso2flux.flux_functions.minimal_flux import add_flux_limit_constraints
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
#To do, move it to emu equations
#from iso2flux.dynamic.interpolate_timecourse import interpolate_timecourse
#from iso2flux.dynamic.interpolate_timecourse import interpolate_timecourse
"""mid_data_name="CA_34_label_reduit.xlsx" 
iso_model_file="isotopomer_model_reduit2.xlsx"""
#Default Value of the parameters
model=None
constraints_file=None
mid_data_name=None
iso_model_file=None
output_prefix=""
null_matrix=None
eqn_dir=None
max_t=20
validate=False
try:
 argv=sys.argv[1:]
 opts, args = getopt.getopt(argv,"e:l:c:s:f:o:t:n:q:t:v",["experimental_data_file=","label_propagation_rules=","constraint_based_model=","settings_file=","flux_constraints=","output_name=","nullmatrix_file=","eqn_dir=","max_reversible_turnover=","validate"])
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
         elif opt in ("--nullmatix","-n"):
              null_matrix=arg #Must be a numpytext
         elif opt in ("--eqn_dir=","-q"):
              eqn_dir=arg
         elif opt in ("--max_reversible_turnover=","-t"):
              max_t=arg
         elif opt in ("-v","--validate"):
              validate=True



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

label_model.build(emu_dict0,force_balance=not validate,recompile_c_code=True,remove_impossible_emus=True,isotopic_steady_state=True,excluded_outputs_inputs=[],turnover_upper_bound=max_t,clear_intermediate_data=False,turnover_exclude_EX=True)


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
   sys.exit(2)
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
      "Validation passed"

