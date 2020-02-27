# -*- coding: utf-8 -*-
import os
import getopt
import sys
import subprocess
import shlex
###Flux functions
import cobra
import iso2flux
from iso2flux.input_functions.create_cobra_model_from_file import create_cobra_model_from_file
##Script Dir
scripts_path = os.path.dirname(os.path.realpath(iso2flux.__file__))+"/scripts/"



try:
 argv=sys.argv[1:]
 print argv
 opts, args = getopt.getopt(argv,"e:l:c:s:f:o:t:q:t:vn:p:g:m:iu:rz:a:",["experimental_data_file=","label_propagation_rules=","constraint_based_model=","settings_file=","flux_constraints=","output_name=","eqn_dir=","max_reversible_turnover=","validate","number_of_processes=","population_size=","generations_per_cycle=","max_cycles_without_improvement=","compute_confidence_intervals","--incubation_time=","--run_p13cfma=","--xi_tolerance=","--gene_expression_file=","--max_flux_for_sampling="])
 #opts, args = getopt.getopt(argv,"n:p:g:m:",["number_of_processes=","population_size=","generations_per_cycle=","max_cycles_without_improvement="])

except getopt.GetoptError as err:
   # print help information and exit:
   print str(err)  # will print something like "option -a not recognized":
   #print "wrong parameters"#'test.py -i <inputfile> -o <outputfile>'
   sys.exit(2)

compute_intervals=False
run_p13cfma=False
constraints_file=None
n_process_str=""
n_pop_str=""
n_gen_str=""
xi_tolerance_str=" -t 3.84" #ChiSquare 1 degree of freedom 95%
integrate_gene_expression_flag=False
max_flux_for_sampling=""
compute_intervals_flag=False


for opt, arg in opts:
         print [opt,arg]
         if opt in ("--experimental_data_file","-e"):
             mid_data_name=arg
         elif opt in ("--label_propagation_rules","-l"):
             iso_model_file=arg
             if "]" in iso_model_file:
                iso_model_file.replace("[","").replace("]","").split(",")       
         elif opt in ("--constraint_based_model","-c"):
             model_name=arg
             try:
                  model=cobra.io.read_sbml_model(arg)
             except:
                  model=create_cobra_model_from_file(arg) 
         elif opt in ("flux_constraints","-f"):
              constraints_file=arg
         elif opt in ("--output_prefix","-o"):
              pass
              #output_name=output_prefix=arg
         elif opt in ("--eqn_dir","-q"):
              eqn_dir=arg
         elif opt in ("--max_reversible_turnover","-t"):
              max_t=arg
         elif opt in ("-v","--validate"):
              validate=True
         elif opt in ("--number_of_processes","-n"):
             #number_of_processes=int(arg)
             n_process_str=" -n "+arg            
         elif opt in ("--population_size","-p"):
             #pop_size=int(arg)
             n_pop_str=" -p "+arg
         elif opt in ("--generations_per_cycle","-g"):
             n_gen=int(arg)
             n_gen_str=" -g "+arg
         elif opt in ("--max_cycles_without_improvement","-m"): #Unused thus far
             max_cycles_without_improvement=int(arg)
         elif opt in ("--compute_confidence_intervals","-i"):
             compute_intervals_flag=True
         elif opt in ("--incubation_time","-u"):
              try:
                 incubation_time=float(arg)
              except: 
                 raise Exception ("Incubation time "+str(arg)+" could not be converted into float")
         elif opt in ("--run_p13cfma","-r"):
              run_p13cfma=True
         elif opt in ("--xi_tolerance","-x"):
             xi_tolerance_str=" -t "+arg
         elif opt in ("--gene_expression_file","-z"):
             gene_expression_file=arg
             integrate_gene_expression_flag=True
         elif opt in ("--max_flux_for_sampling","-a"):
             max_flux_for_sampling=" --max_flux_for_sampling="+str(min_sampling_flux)
             

#######
if max_flux_for_sampling=="":
   min_possible_flux=cobra.flux_analysis.pfba(model).f
   min_sampling_flux=max(min_possible_flux*5,1) #As a default set the max flux for sampling to 5 times the maximum value
   max_flux_for_sampling=" --max_flux_for_sampling="+str(min_sampling_flux)



create_model="python "+ scripts_path+"create_iso2flux_model.py -e "+mid_data_name+"  -c "+model_name+" -l "+iso_model_file+" -o iso2flux "
#create_model="create_iso2flux_model.py -e "+mid_data_name+"  -c "+model_name+" -l "+iso_model_file+" -o iso2flux "

if constraints_file!=None:
   create_model+=" -f"+ constraints_file


solve_model="python "+ scripts_path+"solve_iso2flux_label.py -i iso2flux -o iso2flux"+n_process_str+n_pop_str+n_gen_str+max_flux_for_sampling
#solve_model="solve_iso2flux_label.py -i iso2flux -o iso2flux"+n_process_str+n_pop_str+n_gen_str+max_flux_for_sampling

print solve_model

if run_p13cfma:
    p13cmfa="python"+ scripts_path+"p13cmfa.py -i iso2flux -o iso2flux"+n_process_str+n_pop_str+n_gen_str+ xi_tolerance_str+max_flux_for_sampling
    #p13cmfa="p13cmfa.py -i iso2flux -o iso2flux"+n_process_str+n_pop_str+n_gen_str+ xi_tolerance_str+max_flux_for_sampling
    print p13cmfa
else:
    p13cmfa=""
    
if integrate_gene_expression_flag:
   integrate_gene_expression="python"+ scripts_path+"integrate_gene_expression.py -i iso2flux -f iso2flux_flux_penalty.csv --low_expression_threshold 100 --gene_expression_file  "+ gene_expression_file
   #integrate_gene_expression="integrate_gene_expression.py -i iso2flux -f iso2flux_flux_penalty.csv --low_expression_threshold 100 --gene_expression_file  "+ gene_expression_file
else:
   integrate_gene_expression="" 
   
if compute_intervals_flag:
   get_intervals="python"+ scripts_path+"get_intervals.py -i iso2flux -o iso2flux"+n_process_str+n_pop_str+n_gen_str+ xi_tolerance_str
   #get_intervals="get_intervals.py -i iso2flux -o iso2flux"+n_process_str+n_pop_str+n_gen_str+ xi_tolerance_str
else:
    get_intervals=""
    



f=open("command_list.txt","w")
f.write(create_model+"\n"+solve_model+"\n"+integrate_gene_expression+"\n"+p13cmfa+"\n"+get_intervals)
f.close()



subprocess.call(shlex.split(create_model)) 

   
if integrate_gene_expression!="":
   subprocess.call(shlex.split(integrate_gene_expression)) 



subprocess.call(shlex.split(solve_model))


if p13cmfa!="":
   subprocess.call(shlex.split(p13cmfa))


if get_intervals!="":
   subprocess.call(shlex.split(get_intervals))


 
