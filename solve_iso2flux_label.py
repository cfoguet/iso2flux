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


try:
 argv=sys.argv[1:]
 opts, args = getopt.getopt(argv,"i:n:p:g:m:o:",["iso2flux_model_file=","number_of_processes=","population_size=","generations_per_cycle=","max_cycles_without_improvement=","output_prefix="])
except getopt.GetoptError as err:
   # print help information and exit:
   print str(err)  # will print something like "option -a not recognized":
   sys.exit(2)


pop_size=20
n_gen=200
number_of_processes=6
output_prefix=None
max_cycles_without_improvement=9



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
         elif opt in ("max_cycles_without_improvement=","-m"):
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
     


f=open(output_prefix+"_log.txt","w")
f.write("started at :"+time.strftime('%l:%M%p %Z on %b %d, %Y')+"\n")
f.close()

label_model=load_iso2flux_model(file_name)
#lb_list,ub_list=get_variable_bounds(label_model,max_flux=1e6)





label_problem_parameters={"label_weight":1,"target_flux_dict":None,"max_chi":1e6,"max_flux":1e6,"flux_penalty_dict":{},"verbose":True,"flux_weight":0.0,"label_unfeasible_penalty":1,"flux_unfeasible_penalty":10}


iso2flux_problem=define_isoflux_problem(label_model)
optimal_solution,optimal_variables=optimize(label_model,iso2flux_problem,pop_size = pop_size,n_gen = n_gen,n_islands=number_of_processes ,max_cycles_without_improvement=max_cycles_without_improvement,stop_criteria_relative=0.005,initial_archi_x=[],lb_list=[],ub_list=[],max_flux=1e6,label_problem_parameters=label_problem_parameters)

label_model.best_chi2=optimal_solution
export_flux_results(label_model,optimal_variables,fn=output_prefix+"_fluxes.csv")
objfunc(label_model,optimal_variables)
export_label_results(label_model,fn=output_prefix+"_label.csv",show_chi=True,show_emu=True,show_fluxes=False)
np.savetxt(output_prefix+"_variables.txt",optimal_variables)
label_model.best_label_variables=optimal_variables

save_iso2flux_model(label_model,name=output_prefix+".iso2flux",write_sbml=True,ask_sbml_name=False,gui=False)
f=open(output_prefix+"log.txt","w")
f.write("finished at :"+time.strftime('%l:%M%p %Z on %b %d, %Y')+"\n")
