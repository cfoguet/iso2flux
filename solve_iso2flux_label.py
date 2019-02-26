# -*- coding: utf-8 -*-
"""
This script will run standard 13C MFA on a iso2flux model and provide the optimal flux distribution. Iso2flux uses the PyGMO generalized island-model. Under this model there are several “islands” (each running on a processor thread) and on each island, there is a population of solutions vectors evolves. After a n number of generations, the best solutions migrate across islands. Some of the parameters of the script are used to define the characteristics this island-model. The parameters that can be taken by the script are: 
--iso2flux_model_file=,-I (mandatory):  path to the iso2flux model that will be solved. Iso2flux model can be created with the create_iso2flux_model.py and has  “.iso2flux” as extension. 
--output_prefix=","-o= (optional): Used to define a prefix that will be added to all output files. It can be used both to name the outputs and to select the directory where the outputs will be saved (directory must already exist)
--number_of_processes=,-n(optional): Number of islands (processes) that will be used in the optimization. Should not be larger than the number of processor threads. Default is 4.  
--population_size=(optional),-p(optional):  Size of the population in each island. Default is 20. 
--generations_per_cycle=,-g(optional):  Number of generations that each island will evolve before exchanging individuals with other islands. Default is 200.
--max_cycles_without_improvement=,-m (optional): Maximum number of cycles without a significant improvement of the objective function. When this number is reached the algorithm will stop. Default is 9 

"""

import tkFileDialog
import Tkinter
import sys, getopt
import numpy as np
import copy
import os
import time


from iso2flux import fitting
from iso2flux.misc.save_load_iso2flux_model import load_iso2flux_model, save_iso2flux_model
from iso2flux.flux_functions.minimal_flux import create_minimal_fux_model,flux_minimization_fva
from iso2flux.fitting.variable_sampling import variable_sampling
from iso2flux.fitting.get_variable_bounds import get_variable_bounds
from iso2flux.fitting.objfunc import objfunc
from iso2flux.fitting.extract_results import extract_results
from iso2flux.output_functions.write_fluxes import export_flux_results
from iso2flux.output_functions.export_label_results import export_label_results
from iso2flux.fitting.optimize import optimize,define_isoflux_problem
from iso2flux.fitting.covariance import get_std_deviation

try:
 argv=sys.argv[1:]
 opts, args = getopt.getopt(argv,"i:n:p:g:m:o:s:",["iso2flux_model_file=","number_of_processes=","population_size=","generations_per_cycle=","max_cycles_without_improvement=","output_prefix=","--max_flux_for_sampling="])
except getopt.GetoptError as err:
   # print help information and exit:
   print str(err)  # will print something like "option -a not recognized":
   sys.exit(2)


pop_size=20
n_gen=200
number_of_processes=6
output_prefix=None
max_cycles_without_improvement=9
file_name=None
max_flux_for_sampling=1e6

for opt, arg in opts:
         print [opt,arg]
         if opt in ("--iso2flux_model_file=","-i"):
             file_name=arg
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
         elif opt in ("--max_flux_for_sampling=","-s"):
             max_cycles_without_improvement=float(arg)




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
     


f=open(output_prefix+"_log.txt","w")
f.write("started at :"+time.strftime('%l:%M%p %Z on %b %d, %Y')+"\n")
f.close()

label_model=load_iso2flux_model(file_name)
#lb_list,ub_list=get_variable_bounds(label_model,max_flux=1e6)





label_problem_parameters={"label_weight":1,"target_flux_dict":None,"max_chi":1e6,"max_flux":1e6,"flux_penalty_dict":{},"verbose":True,"flux_weight":0.0,"label_unfeasible_penalty":1,"flux_unfeasible_penalty":10}


iso2flux_problem=define_isoflux_problem(label_model)
optimal_solution,optimal_variables=optimize(label_model,iso2flux_problem,pop_size = pop_size,n_gen = n_gen,n_islands=number_of_processes ,max_cycles_without_improvement=max_cycles_without_improvement,stop_criteria_relative=0.005,initial_archi_x=[],lb_list=[],ub_list=[],max_flux=1e6,label_problem_parameters=label_problem_parameters,migrate="one_direction",max_flux_sampling=max_flux_for_sampling)

label_model.best_chi2=optimal_solution

try:
  flux_sd_dict, hessian,inverse_hessian,covariance=get_std_deviation(label_model,optimal_variables,initial_step=1e-3)
except:
   flux_sd_dict={}

export_flux_results(label_model,optimal_variables,fn=output_prefix+"_fluxes.csv",flux_sd_dict=flux_sd_dict)
objfunc(label_model,optimal_variables)
export_label_results(label_model,fn=output_prefix+"_label.csv",show_chi=True,show_emu=True,show_fluxes=False)
np.savetxt(output_prefix+"_variables.txt",optimal_variables)
label_model.best_label_variables=optimal_variables

save_iso2flux_model(label_model,name=output_prefix+".iso2flux",write_sbml=True,ask_sbml_name=False,gui=False)
f=open(output_prefix+"log.txt","w")
f.write("finished at :"+time.strftime('%l:%M%p %Z on %b %d, %Y')+"\n")
