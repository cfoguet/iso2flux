from scipy.stats import chi2
import math
import random
import copy
import json
import cobra
from cobra.flux_analysis.variability import flux_variability_analysis
from get_objective_function import get_objective_function
from anneal import annealing, coordinate_descent
from apply_parameters import apply_parameters
from clear_parameters import clear_parameters
from ..emus_functions.solver import solver
from ..misc.round_functions import round_up,round_down
from ..flux_functions.check_feasibility import check_feasibility
from multiprocessing import Pool
import openpyxl
import json

def error_dump(label_model,e_type=None):
    import time
    import json
    import cobra
    dump_dict={"parameter_dict":label_model.parameter_dict,"flux_dict":label_model.flux_dict,"turnover_dict":label_model.turnover_flux_dict,"e_type":e_type}
    fname="zzzzzz"+str(time.time())
    with open(fname+".json", 'w') as fp:
         json.dump(dump_dict, fp)
    cobra.io.write_sbml_model(label_model.constrained_model, fname+"sbml")

def estimate_confidence_intervals(label_model,significance=0.95,n_perturbations=100,minimum_perturbation=0.1,max_perturbation=1,parameter_precision=0.01,best_parameter_dict=None,parameter_list=None,fraction_of_optimum=0.9,relative_max_random_sample=0.5, relative_min_random_sample= 0.25,annealing_n=50,annealing_m=100,annealing_p0=0.4,annealing_pf=0.001,output=True,annealing_n_processes=1,annealing_cycle_time_limit=1800, annealing_cycle_max_attempts=5,annealing_iterations=4,fname="confidence.json"):
   if best_parameter_dict==None:
      best_parameter_dict=label_model.parameter_dict
   if parameter_list==None or parameter_list==[]:
      parameter_list=best_parameter_dict.keys()   
   precision=int(-1*(math.log10(parameter_precision)))
   max_random_sample=int(relative_max_random_sample*len(best_parameter_dict))
   min_random_sample=int(relative_min_random_sample*len(best_parameter_dict))
   #chi_parameters_sets_dict={}
   parameter_confidence_interval_dict={}
   flux_confidence_interval_dict={}
   parameter_value_parameters_sets_dict={}
   build_confidence_dicts(parameter_confidence_interval_dict,parameter_value_parameters_sets_dict,best_parameter_dict)
   """for flux in label_model.flux_dict:
       flux_confidence_interval_dict[flux]={"lb":label_model.flux_dict[flux],"ub":label_model.flux_dict[flux]}
       if (flux+"_reverse") in label_model.flux_dict:
          net_flux=label_model.flux_dict[flux]-label_model.flux_dict[flux+"_reverse"]
          flux_confidence_interval_dict["net_"+flux]={"lb":net_flux,"ub":net_flux}"""
   print flux_confidence_interval_dict
   apply_parameters(label_model,best_parameter_dict,parameter_precision=parameter_precision)
   a,b=solver(label_model)
   best_objective,b,c=get_objective_function(label_model,output=False)
   delta_chi = chi2.isf(q=1-significance, df=1)
   signficance_threshold=delta_chi+best_objective
   print signficance_threshold
   original_objectives_bounds={}
   for reaction in label_model.constrained_model.objective: 
              original_objectives_bounds[reaction.id]={}
              original_objectives_bounds[reaction.id]["lb"]=reaction.lower_bound
              original_objectives_bounds[reaction.id]["ub"]=reaction.upper_bound
              original_objectives_bounds[reaction.id]["obj_coef"]=reaction.objective_coefficient
              fva=flux_variability_analysis(label_model.constrained_model,reaction_list=[reaction], fraction_of_optimum=fraction_of_optimum,tolerance_feasibility=label_model.lp_tolerance_feasibility)
              reaction.lower_bound=max(round_down(fva[reaction.id]["minimum"],precision),reaction.lower_bound)
              reaction.upper_bound=min(round_up(fva[reaction.id]["maximum"],precision),reaction.upper_bound)
              reaction.objective_coefficient=0
   
   flux_parameter_list=[]
   
   for parameter in best_parameter_dict:
       if best_parameter_dict[parameter]["type"]=="flux value":
          flux_parameter_list.append(parameter)
   feasability_process = Pool(processes=1)
   for parameter in parameter_list:
       if output:
          with open("estimate_confidence_interval_output.txt", "a") as myfile:
               myfile.write("///////"+parameter+"(lb="+str(best_parameter_dict[parameter]["lb"])+" ub="+str(best_parameter_dict[parameter]["ub"])+ ")\n")
       #chi_parameters_sets_dict[parameter]={}
       variation_range= best_parameter_dict[parameter]["ub"]-best_parameter_dict[parameter]["lb"] 
       delta_parameter=min(max(variation_range/n_perturbations,minimum_perturbation),max_perturbation)
       print delta_parameter
       #Find the highest/lowest value found on previous simulations
       n=1
       sign=1
       is_flux_value=parameter in flux_parameter_list
       min_parameter_dict=parameter_value_parameters_sets_dict[parameter]["lb_parameter_dict"]
       max_parameter_dict=parameter_value_parameters_sets_dict[parameter]["ub_parameter_dict"]
       lb_list=[]
       ub_list=[]
       #lb_list=[best_parameter_dict[parameter]["lb"]]
       #ub_list=[best_parameter_dict[parameter]["ub"]]
       if is_flux_value==True:
          clear_parameters(label_model,parameter_dict=best_parameter_dict,parameter_list=flux_parameter_list, clear_ratios=False,clear_turnover=False,clear_fluxes=True,restore_objectives=False) #Clear all parameters
          #Get real upper and lower bound for the parameters
          for reaction_id in best_parameter_dict[parameter]["reactions"]: 
                fva=flux_variability_analysis(label_model.constrained_model,fraction_of_optimum=0,reaction_list=[reaction_id],tolerance_feasibility=label_model.lp_tolerance_feasibility)
                lb_list.append(fva[reaction_id]["minimum"])
                ub_list.append(fva[reaction_id]["maximum"])
       else:
          lb_list=[best_parameter_dict[parameter]["lb"]]
          ub_list=[best_parameter_dict[parameter]["ub"]]
       parameter_lb=max(lb_list)
       parameter_ub=min(ub_list)
       #Max value
       test_value=parameter_ub
       apply_parameters(label_model,max_parameter_dict)
       parameter_dict,f_best=evaluate_parameter(label_model,copy.deepcopy(max_parameter_dict),flux_parameter_list,parameter,parameter_lb,parameter_ub,test_value,signficance_threshold,feasability_process,parameter_precision,max_perturbation,annealing_n,annealing_m,annealing_p0,annealing_pf,max_random_sample,min_random_sample,annealing_n_processes,annealing_cycle_time_limit, annealing_cycle_max_attempts,annealing_iterations)
       if output:
             with open("estimate_confidence_interval_output.txt", "a") as myfile:
               myfile.write(parameter+" "+"v="+str(test_value)+" chi="+str(f_best)+"\n")
       if f_best>=signficance_threshold:  
         lower_rejected_value=test_value
         higher_accepted_value= max_parameter_dict[parameter]["v"]
         higher_accepted_value_parameter_dict=copy.deepcopy(max_parameter_dict)
         iteration=0
         while(n<=10000):
          iteration+=1
          test_value=(lower_rejected_value+higher_accepted_value)/2.0
          parameter_dict,f_best=evaluate_parameter(label_model,copy.deepcopy(higher_accepted_value_parameter_dict),flux_parameter_list,parameter,parameter_lb,parameter_ub,test_value,signficance_threshold,feasability_process,parameter_precision,max_perturbation,annealing_n,annealing_m,annealing_p0,annealing_pf,max_random_sample,min_random_sample,annealing_n_processes,annealing_cycle_time_limit, annealing_cycle_max_attempts,annealing_iterations)
          if f_best>=signficance_threshold:
             lower_rejected_value=test_value
          else:
              higher_accepted_value_parameter_dict=copy.deepcopy(parameter_dict)
              higher_accepted_value=test_value
              build_flux_confidence_interval_dict(label_model,flux_confidence_interval_dict)
              build_confidence_dicts(parameter_confidence_interval_dict,parameter_value_parameters_sets_dict,parameter_dict)   
          if output:
             with open("estimate_confidence_interval_output.txt", "a") as myfile:
               myfile.write(parameter+" "+"v="+str(test_value)+" chi="+str(f_best)+" higher_accepted_value="+str(higher_accepted_value)+" lower_rejected_value="+str(lower_rejected_value)+"\n")
          if minimum_perturbation>=abs((higher_accepted_value-lower_rejected_value)):
             break
       else:
            build_flux_confidence_interval_dict(label_model,flux_confidence_interval_dict)
            build_confidence_dicts(parameter_confidence_interval_dict,parameter_value_parameters_sets_dict,parameter_dict)
       #Min value
       apply_parameters(label_model,min_parameter_dict)
       test_value=parameter_lb
       parameter_dict,f_best=evaluate_parameter(label_model,copy.deepcopy(max_parameter_dict),flux_parameter_list,parameter,parameter_lb,parameter_ub,test_value,signficance_threshold,feasability_process,parameter_precision,max_perturbation,annealing_n,annealing_m,annealing_p0,annealing_pf,max_random_sample,min_random_sample,annealing_n_processes,annealing_cycle_time_limit, annealing_cycle_max_attempts,annealing_iterations)
       if output:
             with open("estimate_confidence_interval_output.txt", "a") as myfile:
               myfile.write(parameter+" "+"v="+str(test_value)+" chi="+str(f_best)+"\n")
       if f_best>=signficance_threshold:  
         higher_rejected_value=test_value
         lower_accepted_value= max_parameter_dict[parameter]["v"]
         lower_accepted_value_parameter_dict=copy.deepcopy(min_parameter_dict)
         iteration=0
         while(n<=10000):
          iteration+=1
          test_value=(higher_rejected_value+lower_accepted_value)/2.0
          parameter_dict,f_best=evaluate_parameter(label_model,copy.deepcopy(lower_accepted_value_parameter_dict),flux_parameter_list,parameter,parameter_lb,parameter_ub,test_value,signficance_threshold,feasability_process,parameter_precision,max_perturbation,annealing_n,annealing_m,annealing_p0,annealing_pf,max_random_sample,min_random_sample,annealing_n_processes,annealing_cycle_time_limit, annealing_cycle_max_attempts,annealing_iterations)
          if f_best>=signficance_threshold:
             higher_rejected_value=test_value
          else:
              lower_accepted_value_parameter_dict=copy.deepcopy(parameter_dict)
              lower_accepted_value=test_value
              build_flux_confidence_interval_dict(label_model,flux_confidence_interval_dict)
              build_confidence_dicts(parameter_confidence_interval_dict,parameter_value_parameters_sets_dict,parameter_dict)   
          if output:
             with open("estimate_confidence_interval_output.txt", "a") as myfile:
               myfile.write(parameter+" "+"v="+str(test_value)+" chi="+str(f_best)+" higher_rejected_value="+str(lower_accepted_value)+" higher_rejected_value="+str(lower_rejected_value)+"\n")
          if minimum_perturbation>=abs(higher_rejected_value-lower_accepted_value):
             break
       else:
            build_flux_confidence_interval_dict(label_model,flux_confidence_interval_dict)
            build_confidence_dicts(parameter_confidence_interval_dict,parameter_value_parameters_sets_dict,parameter_dict)   
   for reaction_id in original_objectives_bounds:
            reaction=label_model.constrained_model.reactions.get_by_id(reaction_id)
            reaction.lower_bound=original_objectives_bounds[reaction_id]["lb"]
            reaction.upper_bound=original_objectives_bounds[reaction_id]["ub"]
            reaction.objective_coefficient=original_objectives_bounds[reaction_id]["obj_coef"]
   apply_parameters(label_model,best_parameter_dict,parameter_precision=parameter_precision)
   feasability_process.close()
   if "xlsx" in fname:
      save_flux_confidence_interval(flux_confidence_interval_dict,significance=significance,fn=fname) 
   elif "json" in  fname:
      save_confidence_interval_json(flux_confidence_interval_dict,parameter_confidence_interval_dict,fn=fname)     
   return parameter_confidence_interval_dict,flux_confidence_interval_dict,parameter_value_parameters_sets_dict 



def evaluate_parameter(label_model,parameter_dict,flux_parameter_list,parameter,parameter_lb,parameter_ub,parameter_new_value,signficance_threshold,feasability_process,parameter_precision,max_perturbation,annealing_n,annealing_m,annealing_p0,annealing_pf,max_random_sample,min_random_sample,annealing_n_processes,annealing_cycle_time_limit, annealing_cycle_max_attempts,annealing_iterations):
   parameter_dict=copy.deepcopy(parameter_dict)
   lb_list=[parameter_dict[parameter]["lb"]]
   ub_list=[parameter_dict[parameter]["ub"]]
   is_flux_value=parameter in flux_parameter_list
   if is_flux_value==True:
             clear_parameters(label_model,parameter_dict=parameter_dict,parameter_list=flux_parameter_list, clear_ratios=False,clear_turnover=False,clear_fluxes=True,restore_objectives=False) #Clear all parameters
             """#Get real upper and lower bound for the parameters
             for reaction_id in parameter_dict[parameter]["reactions"]: 
                status,feasability_process=check_feasibility(label_model.constrained_model,tolerance_feasibility=label_model.lp_tolerance_feasibility,time_limit=60,pool=feasability_process)
                fva=flux_variability_analysis(label_model.constrained_model,fraction_of_optimum=0,reaction_list=[reaction_id],tolerance_feasibility=label_model.lp_tolerance_feasibility)
                lb_list.append(fva[reaction_id]["minimum"])
                ub_list.append(fva[reaction_id]["maximum"])
   parameter_lb=max(lb_list)
   parameter_ub=min(ub_list)"""
   parameter_dict[parameter]["v"]=parameter_new_value
   print [parameter,parameter_dict[parameter]]
   apply_parameters(label_model,parameter_dict,apply_flux_values=True,parameter_precision=parameter_precision,parameter_list=[parameter])
   if is_flux_value==True:
      parameter_backup=copy.deepcopy(parameter_dict)
      for attempt in range(0,10):
          retry_flag=False
          random_parameter_sample=random.sample(flux_parameter_list, len(flux_parameter_list))
          for flux_value in random_parameter_sample:
                  if flux_value==parameter:
                     continue
                  lb=-999999
                  ub=999999
                  for reaction_id in parameter_dict[flux_value]["reactions"]:
                      status,feasability_process=check_feasibility(label_model.constrained_model,tolerance_feasibility=label_model.lp_tolerance_feasibility,time_limit=60,pool=feasability_process)
                      if status =="infeasible":
                         retry_flag=True 
                         clear_parameters(label_model,parameter_dict=parameter_dict,parameter_list=flux_parameter_list, clear_ratios=False ,clear_turnover=False ,clear_fluxes=True, restore_objectives=False)
                         parameter_dict=copy.deepcopy(parameter_backup)
                         apply_parameters(label_model,parameter_dict,apply_flux_values=True,parameter_precision=parameter_precision,parameter_list=[parameter])
                         break
                      fva=flux_variability_analysis(label_model.constrained_model,fraction_of_optimum=0,reaction_list=[reaction_id],tolerance_feasibility=label_model.lp_tolerance_feasibility)
                      ub=min(fva[reaction_id]["maximum"],ub)
                      lb=max(fva[reaction_id]["minimum"],lb)
                  if retry_flag==True:
                     break 
                  value=parameter_dict[flux_value]["v"] #previous value
                  parameter_dict[flux_value]["v"]=min(max(lb,value),ub)
                  apply_parameters(label_model,parameter_dict,apply_flux_values=True,parameter_precision=parameter_precision,parameter_list=[flux_value])
          #print model.optimize()
          if  retry_flag==False:
              break #If no errors where encountered no need for further attempts
   apply_parameters(label_model,parameter_dict,apply_flux_values=True,parameter_precision=parameter_precision)
   a,b=solver(label_model)
   f_best,b,c=get_objective_function(label_model,output=False)
   best_flux_dict=copy.deepcopy(label_model.flux_dict)
   print ["FBEST",f_best]
   if f_best>=signficance_threshold:
             print "coordinated descent"
             parameters_to_fit=copy.copy(parameter_dict)
             del parameters_to_fit[parameter]
             f_best, new_parameters,best_flux_dict=coordinate_descent(label_model,mode="fsolve",hot_start=False, check_feasability=True,parameter_precision=parameter_precision,parameter_to_be_fitted=parameters_to_fit,max_perturbation=max_perturbation,perturbation=1.2,fba_mode="fba",parameter_dict=parameter_dict)  
             best_flux_dict=label_model.flux_dict
             print [parameter_new_value,f_best]
             if f_best<=signficance_threshold:
                parameter_dict=new_parameters
             else:
                for x in range(0,annealing_iterations): #Try several times to make sure the result is above the signficance threeshol
                    parameter_dict,best_flux_dict,f_best=annealing(label_model,n=annealing_n,m=annealing_m,p0=annealing_p0,pf=annealing_pf,max_random_sample=max_random_sample,min_random_sample=min_random_sample,save_parameters=False,mode="fsolve",hot_start=False, check_feasability=True,fraction_of_optimum=0,parameter_precision=parameter_precision,parameter_to_be_fitted=parameters_to_fit,max_perturbation=max_perturbation,gui=None,fba_mode="fba", break_threshold=signficance_threshold,parameter_dict=parameter_dict,n_processes=annealing_n_processes,cycle_time_limit=annealing_cycle_time_limit, cycle_max_attempts=annealing_cycle_max_attempts,output=False)
                    if f_best<signficance_threshold:
                        break
   return parameter_dict, f_best




def build_flux_confidence_interval_dict(label_model,flux_confidence_interval_dict):
     fva=flux_variability_analysis(label_model.constrained_model, fraction_of_optimum=0,tolerance_feasibility=label_model.lp_tolerance_feasibility)
     for flux in label_model.flux_dict:
              forward_reaction=False
              reverse_reaction=False
              flux_value=label_model.flux_dict[flux]
              if flux in fva:
                 fva_max=fva[flux]["maximum"]
                 fva_min=fva[flux]["minimum"]
              if (flux+"_reverse") in label_model.flux_dict:
                 forward_reaction=True
                 #net_flux=flux_value-flux_dict[flux+"_reverse"]
              elif ("_reverse") in flux:
                 reverse_reaction=True
                 #net_flux=flux_value-flux_dict[flux+"_reverse"]
              if flux not in  flux_confidence_interval_dict:
                 if forward_reaction:
                     flux_confidence_interval_dict[flux]={"lb":fva_min,"ub":fva_max}
                     flux_confidence_interval_dict[flux+"_forward"]={"lb":flux_value,"ub":flux_value} 
                 elif reverse_reaction:
                     flux_confidence_interval_dict[flux]={"lb":flux_value,"ub":flux_value}
                 else: 
                     flux_confidence_interval_dict[flux]={"lb":fva_min,"ub":fva_max}
              else: 
                 if forward_reaction:   
                   flux_confidence_interval_dict[flux]["lb"]=min(flux_confidence_interval_dict[flux]["lb"],fva_min)       
                   flux_confidence_interval_dict[flux]["ub"]=max(flux_confidence_interval_dict[flux]["ub"],fva_max)       
                   
                   flux_confidence_interval_dict[flux+"_forward"]["lb"]=min(flux_confidence_interval_dict[flux+"_forward"]["lb"],flux_value)       
                   flux_confidence_interval_dict[flux+"_forward"]["ub"]=max(flux_confidence_interval_dict[flux+"_forward"]["ub"],flux_value) 
                 elif reverse_reaction:
                   flux_confidence_interval_dict[flux]["lb"]=min(flux_confidence_interval_dict[flux]["lb"],flux_value)       
                   flux_confidence_interval_dict[flux]["ub"]=max(flux_confidence_interval_dict[flux]["ub"],flux_value) 
                 else:     
                   flux_confidence_interval_dict[flux]["lb"]=min(flux_confidence_interval_dict[flux]["lb"],fva_min)       
                   flux_confidence_interval_dict[flux]["ub"]=max(flux_confidence_interval_dict[flux]["ub"],fva_max) 


#def build_parameter_value_parameters_sets_dict(parameter_value_parameters_sets_dict,parameter_dict)
def build_confidence_dicts(parameter_confidence_interval_dict,parameter_value_parameters_sets_dict,parameter_dict):
    for parameter in parameter_dict:
        parameter_value=parameter_dict[parameter]["v"]
        if parameter not in parameter_value_parameters_sets_dict or  parameter not in parameter_confidence_interval_dict:
           parameter_dict_copy=copy.deepcopy(parameter_dict)
           parameter_confidence_interval_dict[parameter]={"upper_limit":parameter_value,"lower_limit":parameter_value}
           parameter_value_parameters_sets_dict[parameter]={"value_lb":parameter_value,"value_ub":parameter_value,"lb_parameter_dict":parameter_dict_copy,"ub_parameter_dict":parameter_dict_copy}
        else:
           if parameter_value<parameter_value_parameters_sets_dict[parameter]["value_lb"]:
              parameter_dict_copy=copy.deepcopy(parameter_dict)
              parameter_confidence_interval_dict[parameter]["lower_limit"]=parameter_value
              parameter_value_parameters_sets_dict[parameter]["value_lb"]=parameter_value 
              parameter_value_parameters_sets_dict[parameter]["lb_parameter_dict"]=parameter_dict_copy
           elif parameter_value>parameter_value_parameters_sets_dict[parameter]["value_ub"]:
              parameter_dict_copy=copy.deepcopy(parameter_dict)
              parameter_value_parameters_sets_dict[parameter]["value_ub"]=parameter_value
              parameter_confidence_interval_dict[parameter]["upper_limit"]=parameter_value
              parameter_value_parameters_sets_dict[parameter]["ub_parameter_dict"]=parameter_dict_copy 

       
    


def save_flux_confidence_interval(flux_confidence_interval_dict={},significance=0.95,fn="confidence.xlsx"):
  wb = openpyxl.Workbook()
  sheet = wb.active
  sheet.title = "Fluxes"
  sheet['A1']="Confidence level:"
  sheet['B1']=significance
  sheet['A2']="Reaction id"
  sheet['B2']="Lower bound"
  sheet['C2']="Upper bound"
  #sheet['D2']="Objective coefficient"
  sorted_reactions=sorted([x for x in flux_confidence_interval_dict],key=lambda v: v.upper())
  n=2
  for reaction_id in sorted_reactions:
      n+=1
      sheet['A'+str(n)]=reaction_id
      sheet['B'+str(n)]=flux_confidence_interval_dict[reaction_id]["lb"]
      sheet['C'+str(n)]=flux_confidence_interval_dict[reaction_id]["ub"]
  
  wb.save(fn)


def save_confidence_interval_json(flux_confidence_interval_dict,parameter_confidence_interval_dict,fn="confidence.json"):
     with open(fn, 'w') as fp:
          json.dump({flux_confidence_interval_dict,parameter_confidence_interval_dict}, fp)
 
