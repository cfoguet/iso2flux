"""A script that performs 13C Metabolic Flux Analysis with the software Iso2flux. It will find the flux distribution most consistent with experimental tracer data and (if not disabled) also estimate confidence interval for fluxes. 

Options:

-e , --experimental_data_file=    Name of file(s) containing the C13 patterns. If more than one file must be entered this can be defined by putting the name of the files in a string separated by a coma (e.g. "file1,file2")

-l , --label_propagation_rules= Name of the file describing the label propagation rules. 

-c,--constrained_based_model=  Name of the file (sbml, xslx or csv) describing the constraint based model that will be used
-s,--settings_file=     Name of the file (xlsx or csv) defining additional settings for Iso2flux (Optional)

f-,--flux_constraints_file=   Name of the file containing additional constraints for the constraint based model (Optional)

-o,--output_prefix=         Prefix appended to all the output files (Optional)

-q,--quick_analysis When this flag is used it disables the confidence interval analysis (Optional)
-w,--working_directory= Name of the working directory (Optional). If none is defined it will use the one where the script is run
-g --gene_expression_file=  Name of the file (xlsx or csv) indicating the gene expression in the conditions of study (Optional)
-m --metabolomics_file=     Name of the file (xlsx or csv) indicating the metabolites that have been detected in the conditions of study (Optional). It will only be used if a gene expression file is provided.
-t --targetted_fluxes       Name of the files (xlsx or csv) indicating the list of fluxes whose confidence intervals will be computed (Optional). If none is provided confidence intervals will be computed for all fluxes.
"""

p_dict={'reactions_with_forced_turnover': [], 'annealing_cycle_time_limit': 1800, 'confidence_max_absolute_perturbation': 10, 'turnover_exclude_EX': True, 'annealing_n_processes': 3, 'annealing_p0': 0.4, 'identify_free_parameters_add_turnover': True, 'minimum_sd': 0.01, 'annealing_max_perturbation': 1, 'turnover_upper_bound': 10, 'confidence_perturbation': 0.1, 'annealing_m': 20, 'annealing_n': 10, 'annealing_relative_max_sample': 0.4, 'confidence_min_absolute_perturbation': 0.05, 'annealing_pf': 0.0001, 'confidence_significance': 0.95, 'identify_free_parameters_change_threshold': 0.005, 'parameter_precision': 0.0001, 'fraction_of_optimum': 0, 'lp_tolerance_feasibility': 1e-09, 'identify_free_parameters_n_samples': 200, 'annealing_relative_min_sample': 0.25, 'annealing_iterations': 2,"gene_expression_mode":"gim3e", "gene_expression_low_expression_threshold":25,"gene_expression_high_expression_threshold":75,"gene_expression_percentile":True,"gene_expression_gene_method":"avearge", "gene_expression_gene_prefix":"","gene_expression_gene_sufix":"_AT","gene_expression_epsilon":1, "gene_expression_lex_epsilon":1e-6,"gene_expression_fraction_optimum":1, "gene_expression_absent_gene_expression_value":50}
if __name__ == "__main__":
     try:
        import iso2flux
     except:
         pass
     import cobra
     from cobra import Model, Reaction, Metabolite
     from cobra.manipulation.modify import convert_to_irreversible
     from cobra.flux_analysis.variability import flux_variability_analysis
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
     from iso2flux.input_functions.create_cobra_model_from_file import create_cobra_model_from_file
     from iso2flux.misc.check_steady_state import check_steady_state
     from iso2flux.misc.check_simulated_fractions import check_simulated_fractions
     from iso2flux.flux_functions.expression_analysis import integrate_omics_imat,integrate_omics_gim3e
     from iso2flux.flux_functions.minimal_flux import add_flux_limit_constraints
     from iso2flux.misc.save_load_iso2flux_model import save_iso2flux_model,load_iso2flux_model
     from iso2flux.misc.round_functions import round_up, round_down 
     #from iso2flux.misc.convert_model_to_irreversible import convert_to_irreversible_with_indicators
     
     from iso2flux.gui.gui import launch_gui
     
     from iso2flux.input_functions.read_experimental_mid import read_experimental_mid
     from iso2flux.input_functions.read_isotopomer_model import read_isotopomer_model
     from iso2flux.input_functions.read_constraints import read_flux_constraints 
     from iso2flux.input_functions.read_isotoflux_settings import read_isotoflux_settings
     from iso2flux.input_functions.read_metabolights import read_metabolights
     from iso2flux.misc.read_spreadsheets import read_spreadsheets
     from iso2flux.doc.open_manual import open_manual
     from iso2flux.flux_functions.check_feasibility import check_feasibility
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
     quick_analysis=False
     gene_expression_file=None
     metabolomics_file=None
     fluxes_to_evaluate_file=None
     try:
      argv=sys.argv[1:]
      #Aargv=("-e hypoxia_label.xlsx -l isotopomer_model_anusha.xlsx -c reference_model.sbml -o esborram_ -f constraints_hipoxia.xlsx -t esborraml.csv").split()
      #argv=("-e output_Midcor_input_iso2flux-2.csv -l isotopomer_model_reduit2.xlsx -p parameters.csv -s simple_model.sbml -t 18 -f Ctr -o example_").split()
      opts, args = getopt.getopt(argv,"e:l:c:s:f:o:w:g:m:qt:",["experimental_data_file=","label_propagation_rules=","constraint_based_model=","settings_file=","flux_constraints=","output_prefix=","working_directory=","gene_expression_file=","metabolomics_file=","quick_analysis","targetted_fluxes="])
      #opts, args = getopt.getopt(argv,"e:l:s:p:c:o:t:f:w:g:m:q",["experimental_data_file=","label_propagation_files=","sbml_model=","parameters_file=","constraints_file=","output_prefix=","time=","factor=","working_directory=","gene_expression_file=","metabolomics_file=","quick_analysis"])
      #opts, args = getopt.getopt(sys.argv[1:],"hi:o:",["ifile=","ofile="])
     except getopt.GetoptError as err:
        # print help information and exit:
        print str(err)  # will print something like "option -a not recognized":
        #print "wrong parameters"#'test.py -i <inputfile> -o <outputfile>'
        sys.exit(2)
     for opt, arg in opts:
         print [opt,arg]
         if opt in ("--experimental_data_file","-e"):
             mid_data_name=arg
         elif opt in ("--label_propagation_rules=","-l"):
             iso_model_file=arg
             if "]" in iso_model_file:
                iso_model_file.replace("[","").replace("]","").split(",")       
         elif opt in ("--settings_file==","-s"):
              try:
                 p_dict.update(read_isotoflux_settings(arg))
              except: 
                 print "settings could not be loaded"
         elif opt in ("--constraint_based_model=","-c"):
              if ".sbml" in arg.lower() or ".xml" in arg.lower():
                  model=cobra.io.read_sbml_model(arg)
              else:
                  model=create_cobra_model_from_file(arg) 
         elif opt in ("flux_constraints=","-f"):
              constraints_file=arg
         elif opt in ("--output_prefix=","-o"):
              output_prefix=arg
         elif opt in ("--quick_analysis","-q"):
              quick_analysis=True
         elif opt in ("--working_directory=","-w"):
              os.chdir(arg)
         elif opt in ("--gene_expression_file=","-g"):
              gene_expression_file=arg
         elif opt in ("--metabolomics_file","-m"):
              metabolomics_file=arg
         elif opt in ("--targetted_fluxes","-t"):
              fluxes_to_evaluate_file=arg
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
        """try:
           model,ratio_dict=read_flux_constraints(model,ratio_dict={},file_name=constraints_file)
        except:
           print "Constraints could not be loaded"
        """
     print [iso_model_file]     
     
     """
     try:
       mid_data_name=sys.argv[1]
     except:
       quit()
     print mid_data_name
     try:
       iso_model_file=sys.argv[2]
     except:
       quit()
     print iso_model_file
     if len(sys.argv)>3:
        print sys.argv[3]
        p_dict.update(read_isotoflux_settings(sys.argv[3]))
     if len(sys.argv)>4:
        print sys.argv[4]
        model=cobra.io.read_sbml_model((sys.argv[4]))     
     if len(sys.argv)>5: 
        try:
          model,ratio_dict=read_flux_constraints(model,ratio_dict={},file_name=sys.argv[5])
        except: 
          print "Wrong constraints file" 
     """
     
     
     fraction_of_optimum=p_dict["fraction_of_optimum"]
     
     
        
     #model,ratio_dict=read_flux_constraints(model,file_name="normoxia_constraints.xlsx")
     fva=flux_variability_analysis(model,fraction_of_optimum=1)
     write_fva(model,fn=output_prefix+"initial_solution_space.csv",fraction=p_dict["fraction_of_optimum"],remove0=False,change_threshold=0.001,mode="full",lp_tolerance_feasibility=p_dict["lp_tolerance_feasibility"])
     #Remember which reactions can be negative to add them as turnover later as negative lower bound may be removed when minimizing flux
     reactions_with_forced_turnover=p_dict["reactions_with_forced_turnover"]
     """for reaction in model.reactions:
         if "EX_" not in reaction.id and reaction.lower_bound<0:
            reactions_with_forced_turnover.append(reaction.id) 
     reversible_fva=add_flux_limit_constraints(model,fraction_of_optimum_objective=0, fraction_of_flux_minimum=10,boundaries_precision=0.001,solver="cplex")
     write_fva(model,fn="normoxia_minimum_fluxes.xlsx",fraction=1,remove0=False,change_threshold=0.001,mode="full",lp_tolerance_feasibility=1e-6)"""
     label_model=Label_model(model,lp_tolerance_feasibility=p_dict["lp_tolerance_feasibility"],parameter_precision=p_dict["parameter_precision"],reactions_with_forced_turnover=reactions_with_forced_turnover) #initialize the label model class   
     read_isotopomer_model(label_model,iso_model_file,header=True)
     missing_reactions_list= find_missing_reactions(label_model).keys()
     if len(missing_reactions_list)>0:
        raise Exception ("Some reactions lack label propagation rules: "+str(missing_reactions_list)) 
     
     #emu_dict0,label_model.experimental_dict =read_experimental_mid(label_model,mid_data_name,emu0_dict={},experimental_dict={},minimum_sd=p_dict["minimum_sd"])
     try:
        emu_dict0,label_model.experimental_dict =read_metabolights(label_model,mid_data_name,selected_condition=factor,selected_time=time,minimum_sd=p_dict["minimum_sd"],rsm=False)
     except:
        emu_dict0,label_model.experimental_dict =read_experimental_mid(label_model,mid_data_name,emu0_dict={},experimental_dict={},minimum_sd=p_dict["minimum_sd"])
     print emu_dict0
     label_model.build(emu_dict0,force_balance=True,recompile_c_code=True,remove_impossible_emus=True,isotopic_steady_state=True,excluded_outputs_inputs=[],turnover_upper_bound=p_dict["turnover_upper_bound"],clear_intermediate_data=False,turnover_exclude_EX=p_dict['turnover_exclude_EX'])
     """label_model.turnover_flux_dict["glc6p_pdif"]={"ub":1000,"lb":0,"v":500}
     label_model.turnover_flux_dict["pyr_pdif"]={"ub":1000,"lb":0,"v":500}
     label_model.turnover_flux_dict["gludxm"]={"ub":1000,"lb":0,"v":500}"""
     check_steady_state(label_model,only_initial_m0=True,threshold=1e-9)#Check initial dy for steady state deviations
     check_simulated_fractions(label_model)
     if gene_expression_file!=None:
        if str(p_dict["gene_expression_mode"]).lower()=="imat":
           hexs,lexs,objective,status, genefva=integrate_omics_imat(label_model.constrained_model,gene_expression_file,fraction_of_optimum=fraction_of_optimum,low_expression_threshold=p_dict["gene_expression_low_expression_threshold"],high_expression_threshold=p_dict["gene_expression_high_expression_threshold"],percentile=p_dict["gene_expression_percentile"],gene_method=p_dict["gene_expression_gene_method"],gene_prefix=p_dict["gene_expression_gene_prefix"],gene_sufix=p_dict["gene_expression_gene_sufix"],metabolite_list_fname=metabolomics_file,epsilon=p_dict["gene_expression_epsilon"],lex_epsilon=p_dict["gene_expression_lex_epsilon"],imat_fraction_optimum=p_dict["gene_expression_fraction_optimum"],label_model=label_model,add_as_constraints=True,boundaries_precision=p_dict["parameter_precision"],solver=None) 
           
           #hexs,lexs,objective,status, genefva=integrate_omics_imat(model,file_name=gene_expression_file,fraction_of_optimum=fraction_of_optimum,low_expression_threshold=25,high_expression_threshold=75,percentile=True,gene_method="max",gene_prefix="",gene_sufix="AT",metabolite_list_fname=None,epsilon=0.1,lex_epsilon=0.0001,imat_fraction_optimum=1,label_model=None,add_as_constraints=True,boundaries_precision=0.01,solver=None) 
     #penalty_dict,objective,genefva=integrate_omics_gim3e(model,file_name="condition1.csv",fraction_of_optimum=p_dict["fraction_of_optimum"],low_expression_threshold=25,absent_gene_expression=50,percentile=True,gene_method="max",gene_prefix="",gene_sufix="AT",metabolite_list_fname=None,label_model=None,epsilon=0.0001,gim3e_fraction_optimum=0.75,add_as_constraints=True,boundaries_precision=0.001) 
        elif str(p_dict["gene_expression_mode"]).lower()=="gim3e" or str(p_dict["gene_expression_mode"]).lower()=="gimme": 
           penalty_dict,objective,genefva=integrate_omics_gim3e(label_model.constrained_model,gene_expression_file,fraction_of_optimum=p_dict["fraction_of_optimum"],low_expression_threshold=p_dict["gene_expression_low_expression_threshold"],absent_gene_expression=p_dict["gene_expression_absent_gene_expression_value"],percentile=p_dict["gene_expression_percentile"],gene_method=p_dict["gene_expression_gene_method"],gene_prefix=p_dict["gene_expression_gene_prefix"],gene_sufix=p_dict["gene_expression_gene_sufix"],metabolite_list_fname=metabolomics_file,label_model=label_model,epsilon=p_dict["gene_expression_lex_epsilon"],gim3e_fraction_optimum=p_dict["gene_expression_fraction_optimum"],add_as_constraints=True,boundaries_precision=p_dict["parameter_precision"])  
        write_fva(label_model.constrained_model,fn=output_prefix+p_dict["gene_expression_mode"]+"_solution.csv",fraction=p_dict["fraction_of_optimum"],remove0=False,change_threshold=0.001,mode="full",lp_tolerance_feasibility=p_dict["lp_tolerance_feasibility"])
     
     
     
     #launch_gui(label_model)
     
     a,b=solver(label_model,mode="fsolve")
     get_objective_function(label_model,output=True)
     def test1(): solver(label_model,mode="fsolve")
     
     timeit(test1,number=1)
     #model=label_model.constrained_model=copy.deepcopy(label_model.metabolic_model)
     #best_parameters=parameters=identify_free_parameters(label_model,parameter_dict={},fraction_of_optimum=fraction_of_optimum,change_threshold=p_dict["identify_free_parameters_change_threshold"],add_turnover=p_dict["identify_free_parameters_add_turnover"],parameter_precision=label_model.parameter_precision,n_samples=p_dict["identify_free_parameters_n_samples"])
     relative_max_sample=p_dict["annealing_relative_max_sample"]
     relative_min_sample=p_dict["annealing_relative_min_sample"]
     print p_dict
     best_parameters={}
     f_best=99999999
     for x in range(0,p_dict["annealing_iterations"]):
         parameters=parameters=identify_free_parameters(label_model,parameter_dict={},fraction_of_optimum=fraction_of_optimum,change_threshold=p_dict["identify_free_parameters_change_threshold"],add_turnover=p_dict["identify_free_parameters_add_turnover"],parameter_precision=label_model.parameter_precision,n_samples=p_dict["identify_free_parameters_n_samples"]) 
         print parameters
         max_random_sample=int(relative_max_sample*len(parameters))
         min_random_sample=int(relative_min_sample*len(parameters))
         parameters,best_flux_dict,f=annealing(label_model,max_random_sample=max_random_sample,min_random_sample=min_random_sample,n=p_dict["annealing_n"],m=p_dict["annealing_m"],p0=p_dict["annealing_p0"],pf=p_dict["annealing_pf"],parameter_precision=label_model.parameter_precision,max_perturbation=p_dict["annealing_max_perturbation"],      fraction_of_optimum=fraction_of_optimum,parameter_dict=parameters,n_processes=p_dict["annealing_n_processes"],output=False)
         if f<f_best:
            f_best=f
            best_parameters=copy.deepcopy(parameters)
         clear_parameters(self.label_model,parameter_dict=None,parameter_list=[], clear_ratios=True,clear_turnover=True,clear_fluxes=True,restore_objectives=True) #Clear previous parameters
         self.label_model.parameter_dict={} #Clear previous parameters
     #Select parameters to evaluate
     parameters_to_evaluate=[]
     export_label_results(label_model,fn=(output_prefix+"best_label.csv"),show_chi=True)
     write_fva(label_model.constrained_model,fn=output_prefix+"best_fluxes.csv",fraction=1,remove0=False,change_threshold=0.001,mode="full",lp_tolerance_feasibility=1e-6)
     if quick_analysis==False: 
        if fluxes_to_evaluate_file!=None:   
            data_rows_dict=read_spreadsheets(file_names=fluxes_to_evaluate_file,csv_delimiter=",")
            for data in data_rows_dict:
              for row in data_rows_dict[data]:
                  for element in row:
                      if element in parameters_to_evaluate:
                         continue 
                      if element in label_model.constrained_model.reactions:
                         parameters_to_evaluate.append(element)
                      elif "/" in element:
                          element_list=element.split("/")
                          if element_list[0] in label_model.constrained_model.reactions and element_list[1] in label_model.constrained_model.reactions:
                             parameters_to_evaluate.append(element)
        print parameters_to_evaluate
        parameter_confidence_interval_dict,flux_confidence_interval_dict,chi_parameters_sets_dict,constrained_model=estimate_confidence_intervals(label_model,significance=p_dict["confidence_significance"],perturbation=p_dict["confidence_perturbation"],min_absolute_perturbation=p_dict["confidence_min_absolute_perturbation"],max_absolute_perturbation=p_dict["confidence_max_absolute_perturbation"],parameter_precision=label_model.parameter_precision,best_parameter_dict=best_parameters,parameter_list=parameters_to_evaluate,fraction_of_optimum=fraction_of_optimum,relative_max_random_sample=relative_max_sample, relative_min_random_sample= relative_min_sample,annealing_n=p_dict["annealing_n"],annealing_m=p_dict["annealing_m"],annealing_p0=p_dict["annealing_p0"],annealing_pf=p_dict["annealing_pf"],output=True,annealing_n_processes=p_dict["annealing_n_processes"],annealing_cycle_time_limit=p_dict["annealing_cycle_time_limit"], annealing_cycle_max_attempts=5,annealing_iterations=p_dict["annealing_iterations"],fname=output_prefix+"_confidence.csv",sbml_name=output_prefix+"_13C_constrained_model.sbml")
        export_label_results(label_model,fn=(output_prefix+"best_label.csv"),show_chi=True)
        write_fva(label_model.constrained_model,fn=output_prefix+"best_fluxes.csv",fraction=1,remove0=False,change_threshold=0.001,mode="full",lp_tolerance_feasibility=1e-6)
        label_model.constrained_model=constrained_model
