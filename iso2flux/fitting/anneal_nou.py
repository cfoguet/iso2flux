import numpy as np
import matplotlib.pyplot as plt
import random
import math
import copy
import time
from multiprocessing import Pool
from cobra.flux_analysis.variability import flux_variability_analysis
from get_objective_function import get_objective_function
from apply_parameters import apply_parameters
from clear_parameters import clear_parameters
from ..emus_functions.solver import solver
from ..misc.round_functions import round_up,round_down
from ..flux_functions.check_feasibility import check_feasibility

if __name__=="__main__":
   from iso2flux.emus_functions.solver import solver
   from iso2flux.misc.round_functions import round_up,round_down
   from iso2flux.flux_functions.check_feasibility import check_feasibility
   from iso2flux.fitting.clear_parameters import clear_parameters

def error_dump(label_model,e_type=None,extra_id="zzzzzz"):
    import time
    import json
    import cobra
    dump_dict={"parameter_dict":label_model.parameter_dict,"flux_dict":label_model.flux_dict,"turnover_dict":label_model.turnover_flux_dict,"e_type":e_type}
    fname=extra_id+str(time.time())
    with open(fname+".json", 'w') as fp:
         json.dump(dump_dict, fp)
    cobra.io.write_sbml_model(label_model.constrained_model, fname+"sbml")

def coordinate_descent(label_model,parameter_dict=None,parameter_to_be_fitted=None,parameter_precision=0.001,perturbation=0.01,max_perturbation=0.1,force_flux_value_bounds=False,mode="fsolve",fba_mode="fba",debug=True,is_subprocess=False):
  """
  function that performs a coordinated descent optimization of the model parameters
  label_model: label_model object
  parameter_dict: dict,optional:
       dict of the current model parameters
  mode: string,optional
       prefered solver used by the emu solver function
  parameter_precision: float, optional
       maximum variation tolerated for parameters
  parameter_to_be_fitted: list, optional 
       list with the parameters that should be fitted 
  perturbation: float,optional
       fraction of the parameter range (lower bound - upper bound) that will be addded or substracted from the parameters value at each iteration.  	
  max_perturbation: float, optional
       maximum absolute variation allowed to parameter value at each step. The formula for computing the value added or substracted at each step is change=min((ub-lb)*perturbation,max_perturbation)
  fba_mode: string, optional
       mode for the linear programing solver (either fba or pfba). In normal use of the program fba should be used
  force_flux_value_bounds: bool,optional
       if set to True it will enforce the bounds of the flux value parameters defined in the parameter dict. If set to false it will only enforce the bounds for flux value parameters thar were objective reactions. Seting this parameter to True can potentially result on reaching unfeasible solutions.
  debug: bool, optional 
  """
  print "coordinate_descent"
  #error_dump(label_model,e_type="coordinated")
  precision=int(-1*(math.log10(parameter_precision)))
  model=label_model.constrained_model
  if parameter_dict==None:
     parameter_dict=label_model.parameter_dict
  if parameter_to_be_fitted==None:
     parameter_to_be_fitted=parameter_dict.keys()    
  #print "precision:"+str(precision)
  apply_parameters(label_model,parameter_dict=parameter_dict,apply_flux_values= True,parameter_precision=parameter_precision,parameter_list=parameter_to_be_fitted) #remove me
  solver(label_model,mode=mode,fba_mode=fba_mode)
  fc,b,c=get_objective_function(label_model)
  f_best=fc
  best_flux_dict=copy.deepcopy(label_model.flux_dict)
  print f_best
  working_parameters_sample=random.sample(parameter_to_be_fitted, len(parameter_to_be_fitted))
  flux_value_parameters=[]
  for parameter in working_parameters_sample:
      if parameter_dict[parameter]["type"]=="flux value":
               flux_value_parameters.append(parameter)
               """for reaction_id in parameter_dict[parameter]["reactions"]:
                   reaction=model.reactions.get_by_id(reaction_id)
                   reaction.lower_bound=round_down(parameter_dict[parameter]["lb"],precision)   
                   reaction.upper_bound=round_up(parameter_dict[parameter]["ub"],precision)"""
  
  if is_subprocess==False:
     feasability_process = Pool(processes=1)
  #print label_model.turnover_flux_dict 
  #print label_model.flux_dict
  for parameter in working_parameters_sample:
      lb_list=[]
      ub_list=[]  
      if parameter not in flux_value_parameters or force_flux_value_bounds:
         lb_list.append(parameter_dict[parameter]["lb"])
         ub_list.append(parameter_dict[parameter]["ub"])
      if parameter in flux_value_parameters:
          #reaction_list=[]
          for reaction_id in parameter_dict[parameter]["reactions"]:
              #reaction_list.append(reaction_id)
              """lower_bound=parameter_dict[parameter]["original_lb"]
              upper_bound=parameter_dict[parameter]["original_ub"]"""
              if  parameter_dict[parameter]["original_objective_coefficient"]!=0:
                   lower_bound=parameter_dict[parameter]["lb"]
                   upper_bound=parameter_dict[parameter]["ub"] 
              else:
                   lower_bound=parameter_dict[parameter]["original_lb"]
                   upper_bound=parameter_dict[parameter]["original_ub"] 
              reaction=model.reactions.get_by_id(reaction_id)
              reaction.lower_bound=lower_bound
              reaction.upper_bound=upper_bound
              fva=flux_variability_analysis(label_model.constrained_model,fraction_of_optimum=0,reaction_list=[reaction_id],tolerance_feasibility=label_model.lp_tolerance_feasibility) #TODO ADD Exception handler
              lb_list.append(fva[reaction_id]["minimum"])
              ub_list.append(fva[reaction_id]["maximum"])
      ub=min(ub_list)
      lb=max(lb_list)
      """if ub-lb<parameter_precision and parameter in flux_value_parameters:
         parameter_dict[parameter]["v"]=model.optimize().solution.dict[reaction.id]
         apply_parameters(label_model,parameter_dict,apply_flux_values= True,parameter_precision=parameter_precision,parameter_list=[parameter])
         if parameter in flux_value_parameters and debug: 
            reaction=model.reactions.get_by_id(parameter)
            print ["A",model.optimize(tolerance_feasibility=label_model.lp_tolerance_feasibility),fva,[reaction.lower_bound,reaction.upper_bound]] 
      continue """   
      #value=parameter_dict[parameter]["v"]
      sign=1
      n_perturbation=0
      stop_flag=False
      # "perturbation "+str(n_perturbation) +" "+parameter+" "+str(parameter_to_be_fitted[parameter]["v"])
      while True:
            accept=True
            #print model.optimize(tolerance_feasibility=label_model.lp_tolerance_feasibility)
            n_perturbation+=1
            #ADD FVA to check for bound
            previous_value=parameter_dict[parameter]["v"]
            if debug:
               print ["A",ub,lb,ub-lb]
            change=min((ub-lb)*perturbation,max_perturbation)
            new_value=previous_value+sign*change
            #print parameter_dict[parameter]["v"]
            parameter_dict[parameter]["v"]=min(max(new_value,lb),ub) 
            #parameter_dict[parameter]["v"]=min(max(new_value,round_up(lb,precision)),round_down(ub,precision))
            #parameter_dict[parameter]["v"]=min(max(new_value,round_up(lb,precision)),round_down(ub,precision)) #Canviat##########################################################################333333333
            if debug:
               print [new_value,change,previous_value]
            #print parameter+" "+str(parameter_to_be_fitted[parameter]["v"])+" "+str(n_perturbation)+" "+str(sign)
            apply_parameters(label_model,parameter_dict,apply_flux_values= True,parameter_precision=parameter_precision,parameter_list=[parameter])
            if parameter in flux_value_parameters:
               if is_subprocess==False: 
                  status,feasability_process=check_feasibility(label_model.constrained_model,tolerance_feasibility=label_model.lp_tolerance_feasibility,time_limit=60,pool=feasability_process)
               else:
                  status=model.optimize(tolerance_feasibility=label_model.lp_tolerance_feasibility).status 
               if status=="infeasible":
                  accept=False
            if accept==True:
               # print "perturbation "+str(n_perturbation) +" "+parameter+" "+str(parameter_to_be_fitted[parameter]["v"])
               solver(label_model)
               fc,b,c=get_objective_function(label_model)
               #print fc
               """if model.optimize().status!="optimal": 
               return"""
               #print [fc,f_best]
               #if parameter_to_be_fitted[parameter]["v"]==lb or  parameter_to_be_fitted[parameter]["v"]==ub:
               #   print "MIN/MAX " +parameter 
               #   stop_flag=True  
               if fc+0.00001>f_best:
                  accept=False
               else:
                  accept=True
            if not accept: 
               sign*=-1
               if debug:
                  print ["refused",parameter,fc,parameter_dict[parameter]["v"]]
               parameter_dict[parameter]["v"] =previous_value
               apply_parameters(label_model,parameter_dict,apply_flux_values= True,parameter_precision=parameter_precision,parameter_list=[parameter])
               if parameter in flux_value_parameters and debug: 
                  reaction=model.reactions.get_by_id(parameter) 
                  print ["restored",model.optimize(tolerance_feasibility=label_model.lp_tolerance_feasibility)  ,fva,[reaction.lower_bound,reaction.upper_bound]] 
               #solver(label_model)
               #f2,b,c=get_objective_function(label_model)
               """if f2!=f_best:
                  print sign
                  print "previous"+str(previous_value)
                  print parameter
                  print parameter_to_be_fitted[parameter]
                  print [f2,f_best]
                  print label_model.flux_dict
                  print label_model.turnover_flux_dict
                  return"""
               #print " restoring "+str(n_perturbation) +" "+parameter+" "+str(parameter_to_be_fitted[parameter]["v"])
               stop_flag=True
            else:
               if debug:
                  print ["accepted",parameter,fc,parameter_dict[parameter]["v"]]
               f_best=fc
               best_flux_dict=copy.deepcopy(label_model.flux_dict)
               """print "BEST["
               label_model.best_flux_dict=copy.deepcopy(label_model.flux_dict)
               label_model.best_parameter_dict=copy.deepcopy(parameter_to_be_fitted)
               print parameter
               print parameter_to_be_fitted[parameter] 
               print [fc,f_best]
               print label_model.flux_dict
               print label_model.turnover_flux_dict
               print "]BEST"""
            if stop_flag==True:
               """print "stop"""
               if n_perturbation>1:
                  break
  
  """print "Final["
  print parameter_to_be_fitted["aldo_turnover"]
  print parameter_to_be_fitted["tk1_turnover"]
  ftemp4,b,c=get_objective_function(label_model)
  apply_parameters(label_model,parameter_to_be_fitted,apply_flux_values= True,parameter_precision=parameter_precision,parameter_list=[])
  solver(label_model)
  #print label_model.flux_dict
  ftemp5,b,c=get_objective_function(label_model)
  print [f_best,ftemp4]
  print label_model.turnover_flux_dict
  print [label_model.flux_dict["aldo"],label_model.flux_dict["tk1"]]
  print "]FINAL"""
  #REMOVE ME
  print f_best
  """for parameter in flux_value_parameters:
      reaction=model.reactions.get_by_id(parameter) 
      print [parameter,parameter_dict[parameter]["v"],[reaction.lower_bound,reaction.upper_bound]]""" 
  apply_parameters(label_model,parameter_dict,parameter_precision=parameter_precision)
  solver(label_model,mode=mode)
  fi,b,c=get_objective_function(label_model)
  """print fi
  for parameter in flux_value_parameters:
      reaction=model.reactions.get_by_id(parameter) 
      print [parameter,parameter_dict[parameter]["v"],[reaction.lower_bound,reaction.upper_bound]]"""
  if is_subprocess==False: 
     feasability_process.close()
  return f_best, parameter_dict,best_flux_dict



##################################################
# Simulated Annealing
##################################################
def annealing(label_model,parameter_dict=None,parameter_to_be_fitted=None,parameter_precision=None,n=50,m=50,p0=0.7,pf=0.001,max_random_sample=None,min_random_sample=None,max_perturbation=1,min_perturbation=0.1,break_threshold=0,n_processes=1,cycle_time_limit=3600, cycle_max_attempts=10,fraction_of_optimum=0,mode="fsolve",fba_mode="fba",force_flux_value_bounds=False,output=True,gui=None):
  """function that performs a coordinated descent optimization of the model parameters
  label_model: label_model object
  parameter_dict: dict,optional:
       dict of the current model parameters
  mode: string,optional
       prefered solver used by the emu solver function
  parameter_precision: float, optional
       maximum variation tolerated for parameters
  parameter_to_be_fitted: list, optional 
       list with the parameters that should be fitted 
  n: int,optional
     number of simulated annealing cycles
  m: int, optional
     number of iterations for cycle
  p0: float, optional
     probability of accepting a non-optimal solution in the first cycle
  pf: float, optional
     probability of accepting a non-optimal solution in the last cycle
  max_random_sample: int, optional
     maximum number of parameters perturbed at each iteration
  min_random_sample: int, optional
     minimal number of parameters perturbed at each iteration
  max_perturbation: float, optional
       maximum absolute variation allowed to parameter value at each step. 
  min_perturbation: float, optional 
      minimimum possible absolute variation allowed to parameter value at each step. Prevents parameters close to 0 from not varying 
  break_threshold: float, optional
      Chi Square value at wich the algorythm will be stopped (used for computing confidence intervals)
  n_processes: float,optional
      number of paralel process that will be launched at each cycle. A value larger than 1 allows to take advantage of multiple CPU threads. 
  cycle_time_limit: float,optional
      Cycles that last more than this value will be interrupted and restasted. Used to prevent the program for hanging if one of the process encounters and error and hangs 
  cycle_max_attempts: int, optional:
      Max Number of times a cycle will be restarted after the cycle time limit has been reached before stoping the program
  fraction_of_optimum: float, optional
      The minimum/maximum fraction of the reaction objective that all valid parameters should satisfy. If the parameters have been added automatcally this is not necessary as the automatic process removes the objective and converts it into a parameter
  fba_mode: string, optional
       mode for the linear programing solver (either fba or pfba). In normal use of the program fba should be used
  force_flux_value_bounds: bool,optional
       if set to True it will enforce the bounds of the flux value parameters defined in the parameter dict. If set to false it will only enforce the bounds for flux value parameters thar were objective reactions. Seting this parameter to True can potentially result on reaching unfeasible solutions.
  output: bool,optional
       if set to True a fig named annealing.png will be created. This figure represent the Chi square obtained at each cycle
  gui: tkinter widged, optional
       if GUI is active this will be used update the label figure."""
  if parameter_precision==None:
     parameter_precision=label_model.parameter_precision
  precision=int(-1*(math.log10(parameter_precision)))
  model=label_model.constrained_model
  #print "precision:"+str(precision)
  if parameter_dict==None:
     parameter_dict=label_model.parameter_dict 
  if parameter_to_be_fitted==None:
     parameter_to_be_fitted=parameter_dict.keys()
  if len(model.objective)>1:
      raise ValueError('Error:Only one objective supported')
  original_objectives_bounds={}
  for reaction in model.objective: 
              original_objectives_bounds[reaction.id]={}
              original_objectives_bounds[reaction.id]["lb"]=reaction.lower_bound
              original_objectives_bounds[reaction.id]["ub"]=reaction.upper_bound
              original_objectives_bounds[reaction.id]["obj_coef"]=reaction.objective_coefficient
              fva=flux_variability_analysis(model,reaction_list=[reaction], fraction_of_optimum=fraction_of_optimum,tolerance_feasibility=label_model.lp_tolerance_feasibility)
              reaction.lower_bound=max(round_down(fva[reaction.id]["minimum"],precision),reaction.lower_bound)
              reaction.upper_bound=min(round_up(fva[reaction.id]["maximum"],precision),reaction.upper_bound)
              reaction.objective_coefficient=0
              """if reaction.id in parameter_dict:
                 parameter_dict[reaction.id]["original_lb"]=fva[reaction.id]["minimum"]
                 parameter_dict[reaction.id]["original_ub"]=fva[reaction.id]["maximum"]"""
  if max_random_sample==None:
     max_random_sample=len(parameter_to_be_fitted)
  if min_random_sample==None:
     min_random_sample=1
  """if fraction_of_optimum!=0:
     check_feasability==True"""
  # n=Number of cycles
  # m=Number of trials per cycle
  # p0=Probability of accepting worse solution at the start
  # pf Probability of accepting worse solution at the end  
  # Number of accepted solutions
  na = 0.0
  # Initial temperature
  t1 = -1.0/math.log(p0)
  # Final temperature
  tf = -1.0/math.log(pf)
  # Fractional reduction every cycle
  frac = (tf/t1)**(1.0/(n-1.0))
  # Initialize x
  working_parameters=copy.deepcopy(parameter_dict)
  #parameter_list=[]
  #Identfy flux value type parameters and store they default lower and upper_bounds
  flux_value_parameter_list=[]
  #original_reaction_bounds_dict={}
  for parameter in working_parameters:
      if "type" in working_parameters[parameter]: # and check_feasability==True:
         if working_parameters[parameter]["type"]=="flux value":
            flux_value_parameter_list.append(parameter)
            #parameter_lb= working_parameters[parameter]["lb"]
            #parameter_ub= working_parameters[parameter]["ub"]
            #for reaction in working_parameters[parameter]["reactions"]:
                #original_reaction_bounds_dict[reaction]={"lb":max(parameter_dict[parameter]["original_lb"],parameter_lb),"ub":min(parameter_dict[parameter]["original_ub"],parameter_ub)}"""
  #http://apmonitor.com/me575/index.php/Main/SimulatedAnnealing
  # Current best results so far
  current_parameters=copy.deepcopy(parameter_dict)
  best_parameters=copy.deepcopy(current_parameters)
  best_flux_dict=copy.deepcopy(label_model.flux_dict) 
  apply_parameters(label_model,apply_flux_values=True,parameter_precision=parameter_precision,parameter_dict=current_parameters)
  solver(label_model,mode=mode,fba_mode=fba_mode)
  fc,b,c=get_objective_function(label_model)
  f_best=fc
  #f_best,best_parameters=coordinate_descent2(label_model,mode="fsolve",hot_start=hot_start,parameter_precision=parameter_precision,parameter_dict=best_parameters,max_perturbation=max_perturbation, original_reaction_bounds_dict=original_reaction_bounds_dict,perturbation=1.105,fba_mode=fba_mode,parameter_to_be_fitted=parameter_to_be_fitted)
  print f_best
  """fc,current_parameters=coordinate_descent(label_model,mode="fsolve",hot_start=hot_start,parameter_precision=parameter_precision,parameter_to_be_fitted=current_parameters,max_perturbation=max_perturbation,original_reaction_bounds_dict=original_reaction_bounds_dict,perturbation=1.2)"""
  f_best=fc 
  print f_best
  fs = np.zeros(n+1)
  fs[0] = fc
  # Current temperature
  t = t1
  # DeltaE Average
  DeltaE_avg = 0.0
  #cycle_arg_dict={"label_model":label_model,"m":m,"i":0,"max_perturbation":max_perturbation,"min_perturbation":min_perturbation,"current_parameters":current_parameters,"parameter_to_be_fitted":parameter_to_be_fitted, "min_random_sample":min_random_sample,"max_random_sample":max_random_sample,"parameter_precision":parameter_precision,"DeltaE_avg":DeltaE_avg,"fc":fc, "f_best":f_best,"t":t,"flux_value_parameter_list":flux_value_parameter_list,"mode":mode,"fba_mode":fba_mode,"na":na,"best_parameters":best_parameters,"best_flux_dict":best_flux_dict,"force_flux_value_bounds":force_flux_value_bounds}
  if n_processes>1:
     pool = Pool(processes=n_processes)
  for i in range(n):
    cycle_arg_dict={"label_model":label_model,"m":m,"i":i,"max_perturbation":max_perturbation,"min_perturbation":min_perturbation,"current_parameters":current_parameters,"parameter_to_be_fitted":parameter_to_be_fitted, "min_random_sample":min_random_sample,"max_random_sample":max_random_sample,"parameter_precision":parameter_precision,"DeltaE_avg":DeltaE_avg,"fc":fc, "f_best":f_best,"t":t,"flux_value_parameter_list":flux_value_parameter_list,"mode":mode,"fba_mode":fba_mode,"na":na,"best_parameters":best_parameters,"best_flux_dict":best_flux_dict,"force_flux_value_bounds":force_flux_value_bounds}
    print 'Cycle: ' + str(i) + ' with Temperature: ' + str(t)+" best fit:"+str(fc)
    #print label_model.constrained_model.optimize()
    #cycle_arg_dict["i"]=i
    #pool = Pool(processes=4)
    if n_processes>1:
       task_completed=False
       n_attempt=0
       time_limit=cycle_time_limit
       while not(task_completed):
           n_attempt+=1
           task_start = time.time()   # start time
           task = pool.map_async(cycle, [cycle_arg_dict]* n_processes)
           timeout=False
           while not(task.ready()):
                 if task._number_left==1 and cycle_time_limit==time_limit: #When all task are done except one the remaining task should finish soon if it has not crashed
                    time_limit=2*(time.time() - task_start)
                 if (time.time() - task_start) > time_limit: # check maximum time (user def.)
                    print "timeout, restaring process"
                    pool.terminate()                    # kill old pool
                    pool = Pool(processes=n_processes)      # create new pool
                    timeout = True                         # redo computation
                    break                               # break loop, (not finished)
           if timeout == False:
              results = task.get() 
              task_completed=True
           if n_attempt>cycle_max_attempts:
              raise Exception("Annealing can't complete cycles within timelimit")
       print [time.time() - task_start]
       #results=pool.map(cycle, [cycle_arg_dict]* n_processes)
       new_best=9999999
       for n_process,result in enumerate(results):
        print [n,"f_best",result["f_best"]]
        if result["f_best"]<new_best:
           new_best=result["f_best"]
           n_best=n_process
        elif result["f_best"]==new_best and result["fc"]<results[n_best]["fc"]:
           n_best=n_process
       #cycle_arg_dict=copy.deepcopy(results[n_best])
       best_result=copy.deepcopy(results[n_best])
       
    else:
       best_result=cycle_arg_dict=cycle(cycle_arg_dict)
    
    """best_parameters=cycle_arg_dict["best_parameters"]
    current_parameters=cycle_arg_dict["current_parameters"]
    f_best=cycle_arg_dict["f_best"]
    fc=cycle_arg_dict["fc"]
    best_flux_dict=cycle_arg_dict["best_flux_dict"]"""
    best_parameters=best_result["best_parameters"]
    current_parameters=best_result["current_parameters"]
    fc=f_best=best_result["f_best"]
    #fc=best_result["fc"]
    best_flux_dict=best_result["best_flux_dict"]
    DeltaE_avg=best_result["DeltaE_avg"]
    apply_parameters(label_model,current_parameters,parameter_precision=parameter_precision)
    solver(label_model,mode=mode)
    fi,b,c=get_objective_function(label_model)
    #apply_parameters(label_model,current_parameters,parameter_precision=parameter_precision)
    #solver(label_model,mode=mode)
    #fi,b,c=get_objective_function(label_model)    
    print [fi,b,c]
    #error_dump(label_model,e_type=None,extra_id="coord")
    #return
    if gui!=None:
          apply_parameters(label_model,best_parameters,parameter_precision=parameter_precision)
          solver(label_model,mode=mode)
          gui.update_label() 
          gui.root.update_idletasks()
          fi,b,c=get_objective_function(label_model)
          
    if best_result["f_best"]<=break_threshold:
                  if n_processes>1:
                     pool.close()
                  return best_parameters,best_flux_dict, f_best
    print ["f_best",best_result["f_best"]]      
    #cycle_arg_dict=results[0] #TODO GET THE BEST FIT
    #cycle_arg_dict=cycle(cycle_arg_dict)
    print ["f_best2",f_best]
    #print ["f",fc,f_best]
    # Record the best x values at the end of every cycle
    """x[i+1][0] = copy.deepcopy(current_parameters)
    x[i+1][0]["obj"] = fc
    fs = fc"""
    fs[i+1] = fc
    # Lower the temperature for next cycle
    """if save_parameters:
       parameter_list.append(current_parameters)"""
    t = frac * t
    #cycle_arg_dict["t"]=t
  if output:
     plt.plot(fs,'r.-')
     plt.savefig('annealing.png')
     parameter_dict=best_parameters
     plt.clf() #Clears buffer
  if n_processes>1:
     pool.close()
  best_parameters=best_result["best_parameters"]
  best_flux_dict=best_result["best_flux_dict"]
  label_model.constrained_model=best_result["label_model"].constrained_model
  print 'Cycle: ' + str(n) + ' with Temperature: ' + str(t)+" best fit:"+str(fc)
  #f_best,best_parameters,best_flux_dict=coordinate_descent(label_model,mode="fsolve",parameter_precision=parameter_precision,parameter_dict=best_parameters,max_perturbation=max_perturbation,perturbation=0.01,fba_mode=fba_mode,parameter_to_be_fitted=parameter_to_be_fitted,force_flux_value_bounds=force_flux_value_bounds)
  print("best Chi:"+str(f_best))
  apply_parameters(label_model,best_parameters,parameter_precision=parameter_precision)
  solver(label_model,mode=mode)
  fi,b,c=get_objective_function(label_model)
  print fi
  for reaction_id in original_objectives_bounds:
            reaction=model.reactions.get_by_id(reaction_id)
            reaction.lower_bound=original_objectives_bounds[reaction_id]["lb"]
            reaction.upper_bound=original_objectives_bounds[reaction_id]["ub"]
            reaction.objective_coefficient=original_objectives_bounds[reaction_id]["obj_coef"]
  """if save_parameters:
     return parameter_list"""
  return best_parameters,best_flux_dict, f_best


























def cycle(cycle_arg_dict):
    """
    function that does a cycle of simulated annealing
    cycle_arg_dict: dict
        dict with all the parameters for running a cycle of simulated annealing
    """
    debug=False
    local_cycle_arg_dict=copy.deepcopy(cycle_arg_dict)
    #print local_cycle_arg_dict
    #return local_cycle_arg_dict
    #current_parameters=local_cycle_arg_dict["current_parameters"]  
    #parameter_to_be_fitted=local_cycle_arg_dict["parameter_to_be_fitted"]
    #min_random_sample=local_cycle_arg_dict["min_random_sample"]
    #max_random_sample=local_cycle_arg_dict["max_random_sample"]
    #parameter_precision=local_cycle_arg_dict["parameter_precision"]
    #print parameter_precision
    #DeltaE_avg=local_cycle_arg_dict["DeltaE_avg"]
    #fc=local_cycle_arg_dict["fc"]
    #f_best=local_cycle_arg_dict["f_best"]
    #t=local_cycle_arg_dict["t"]
    #i=local_cycle_arg_dict["i"]
    #m=local_cycle_arg_dict["m"]
    #mode=local_cycle_arg_dict["mode"]
    #original_reaction_bounds_dict=local_cycle_arg_dict["original_reaction_bounds_dict"]
    #flux_value_parameter_list=local_cycle_arg_dict["flux_value_parameter_list"]
    #fba_mode=local_cycle_arg_dict["fba_mode"]
    #gui=local_cycle_arg_dict["gui"]
    #break_threshold=local_cycle_arg_dict["break_threshold"]
    #na=local_cycle_arg_dict["na"]
    #min_random_sample=local_cycle_arg_dict["min_random_sample"]
    #max_random_sample=local_cycle_arg_dict["max_random_sample"]
    
    original_f_best=local_cycle_arg_dict["f_best"]
    model=local_cycle_arg_dict["label_model"].constrained_model
    precision=int(-1*(math.log10(local_cycle_arg_dict["parameter_precision"])))
    print ["initialfc", local_cycle_arg_dict["fc"]]
    for j in range(local_cycle_arg_dict["m"]):
        #print "j"+str(j)
        working_parameters=copy.deepcopy(local_cycle_arg_dict["current_parameters"]  )   
        sample_size = random.randint(local_cycle_arg_dict["min_random_sample"],local_cycle_arg_dict["max_random_sample"]) 
        working_parameters_sample=random.sample(local_cycle_arg_dict["parameter_to_be_fitted"], sample_size)
        #Remove constraints for the flux values parameters we will fit
        sampled_flux_value_parameters=[]
        failure_flag=False
        counter=0
        for parameter in working_parameters_sample: 
            if parameter in local_cycle_arg_dict["flux_value_parameter_list"]:
               counter+=1
               reaction=model.reactions.get_by_id(parameter)
               ###print [counter,reaction.id,reaction.lower_bound,reaction.upper_bound,model.optimize()]
               #backup_model=copy.deepcopy(model)
               clear_parameters(local_cycle_arg_dict["label_model"],parameter_dict=working_parameters,parameter_list=[parameter], default_turnover=0, clear_ratios=False,clear_turnover=False,clear_fluxes=True,restore_objectives=False)
               sampled_flux_value_parameters.append(parameter)
               """if model.optimize(tolerance_feasibility=local_cycle_arg_dict["label_model"].lp_tolerance_feasibility).status!="optimal":
                  ###print [reaction.id,reaction.lower_bound,reaction.upper_bound,model.optimize(tolerance_feasibility=local_cycle_arg_dict["label_model"].lp_tolerance_feasibility)] 
                  error_dump(local_cycle_arg_dict["label_model"],"not optimal 1","clear2")
                  #local_cycle_arg_dict["label_model"].constrained_model=backup_model
                  #error_dump(local_cycle_arg_dict["label_model"],"not optimal 1","clear1") 
                  #return """
               """ 
               lb= working_parameters[parameter]["lb"]
               ub= working_parameters[parameter]["ub"]
               original_lb= working_parameters[parameter]["original_lb"]
               original_ub= working_parameters[parameter]["original_ub"]
               if  working_parameters[parameter]["original_objective_coefficient"]!=0:
                   lower_bound=lb
                   upper_bound=ub
               else:
                   lower_bound=original_lb#round_down(lb,precision)
                   upper_bound=original_ub#round_up(ub,precision)
               #lower_bound=max(round_down(lb,precision),original_lb)
               #upper_bound=min(round_up(ub,precision),original_ub)
               for reaction_id in working_parameters[parameter]["reactions"]:
                   reaction=model.reactions.get_by_id(reaction_id)
                   reaction.lower_bound=lower_bound
                   reaction.upper_bound=upper_bound"""
            else:
                   rand=random.uniform(-local_cycle_arg_dict["current_parameters"]  [parameter]["max_d"], local_cycle_arg_dict["current_parameters"]  [parameter]["max_d"])
                   lb= working_parameters[parameter]["lb"]
                   ub= working_parameters[parameter]["ub"]
                   """if local_cycle_arg_dict["current_parameters"][parameter]["v"]-lb>local_cycle_arg_dict["max_perturbation"]:
                      perturbation=local_cycle_arg_dict["max_perturbation"]
                   else:
                      perturbation=local_cycle_arg_dict["current_parameters"][parameter]["v"]-lb"""
                   """delta_parameter=max(abs(working_parameters[parameter]["v"]),local_cycle_arg_dict["min_perturbation"])*rand #If the parameter is to close to 0 we want to use the value of min_perturbation instead of the parameter value
                   delta_parameter=max(min(delta_parameter*rand,local_cycle_arg_dict["max_perturbation"]),-local_cycle_arg_dict["max_perturbation"]) #Ensure that the perturbation is not larger than max_perturbation
                   #print delta_parameter        
                   value=working_parameters[parameter]["v"]+delta_parameter"""
                   rand=random.uniform(1-local_cycle_arg_dict["current_parameters"]  [parameter]["max_d"], 1+local_cycle_arg_dict["current_parameters"]  [parameter]["max_d"])
                   lb= working_parameters[parameter]["lb"]
                   ub= working_parameters[parameter]["ub"]
                   if local_cycle_arg_dict["current_parameters"][parameter]["v"]-lb>local_cycle_arg_dict["max_perturbation"]:
                      perturbation=local_cycle_arg_dict["max_perturbation"]
                   elif local_cycle_arg_dict["current_parameters"][parameter]["v"]-lb<local_cycle_arg_dict["min_perturbation"]:
                      perturbation=local_cycle_arg_dict["min_perturbation"]
                   else:
                      perturbation=local_cycle_arg_dict["current_parameters"][parameter]["v"]-lb
                   relative_value=(perturbation)*rand
                   value=(lb+relative_value)
                   value=max(min(value,ub),lb)
                   working_parameters[parameter]["v"] =value 
                   apply_parameters(local_cycle_arg_dict["label_model"],working_parameters,apply_flux_values= True,parameter_precision=local_cycle_arg_dict["parameter_precision"],parameter_list=[parameter])
        if len(sampled_flux_value_parameters)<sample_size:
           #print sampled_flux_value_parameters 
           #If ratios or other parameters are being fitted check if they make the solution unfeasiable
           if model.optimize(tolerance_feasibility=local_cycle_arg_dict["label_model"].lp_tolerance_feasibility).status!="optimal": 
              print ["not optimal 1",j,model.optimize(tolerance_feasibility=local_cycle_arg_dict["label_model"].lp_tolerance_feasibility)]
              failure_flag=True
              error_dump(local_cycle_arg_dict["label_model"],"not optimal 1","not_optimal_1")
              local_cycle_arg_dict["label_model"].error_parameters=local_cycle_arg_dict["current_parameters"]
              #return
        if failure_flag==False:
           #print "start"
           #print "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" 
           #clear_parameters(local_cycle_arg_dict["label_model"],parameter_dict=working_parameters,parameter_list=sampled_flux_value_parameters, default_turnover=0, clear_ratios=False,clear_turnover=False,clear_fluxes=True,restore_objectives=False)    
           for parameter in sampled_flux_value_parameters:
                  lb_list=[]
                  ub_list=[]
                  if local_cycle_arg_dict["force_flux_value_bounds"]:
                     lb_list.append(working_parameters[parameter]["lb"])
                     ub_list.append(working_parameters[parameter]["ub"])
                  #original_lb_list=[]
                  #original_ub_list=[]
                  
                  for reaction in local_cycle_arg_dict["current_parameters"][parameter]["reactions"]:
                      try:
                         fva=flux_variability_analysis(model,reaction_list=[reaction], fraction_of_optimum=0,tolerance_feasibility=local_cycle_arg_dict["label_model"].lp_tolerance_feasibility)
                         lb_list.append(fva[reaction]["minimum"])
                         ub_list.append(fva[reaction]["maximum"])
                         #original_lb_list.append(local_cycle_arg_dict["original_reaction_bounds_dict"][reaction]["lb"])
                         #original_ub_list.append(local_cycle_arg_dict["original_reaction_bounds_dict"][reaction]["ub"])
                         #print [model.optimize(),parameter,fva]
                      except:
                        print "errror 1"
                        error_dump(local_cycle_arg_dict["label_model"],"error1")
                        failure_flag=True
                        #return 
                  if failure_flag==True:
                     print "not optimal 1"
                     error_dump(local_cycle_arg_dict["label_model"],"error1")
                     apply_parameters(local_cycle_arg_dict["label_model"], local_cycle_arg_dict["current_parameters"],apply_flux_values=True,parameter_precision=local_cycle_arg_dict["parameter_precision"])
                     #return 
                  else:
                     lb=max(lb_list)
                     ub=min(ub_list)
                     #original_lb=max(original_lb_list)
                     #original_ub=min(original_ub_list)
                     if abs(lb-ub)<local_cycle_arg_dict["parameter_precision"]:
                        value=((lb+ub)/2.0)
                        #value=min(max(value,lb+0.5*local_cycle_arg_dict["parameter_precision"]),ub-0.5*local_cycle_arg_dict["parameter_precision"]) #Make sure value will not be oustide original lower and upper               
                        value=min(max(value,lb),ub)
                        working_parameters[parameter]["v"]=value
                        apply_parameters(local_cycle_arg_dict["label_model"],working_parameters,apply_flux_values= True,parameter_precision=local_cycle_arg_dict["parameter_precision"],parameter_list=[parameter])
                        #Esborram
                        if debug:
                           reaction=local_cycle_arg_dict["label_model"].constrained_model.reactions.get_by_id(parameter)
                           print [fva,[reaction.lower_bound,reaction.upper_bound],[lb,ub]]
                           print model.optimize(tolerance_feasibility=local_cycle_arg_dict["label_model"].lp_tolerance_feasibility)
                        
                        # [model.reactions.get_by_id(reaction_id),model.reactions.get_by_id(reaction_id).lower_bound,model.reactions.get_by_id(reaction_id).upper_bound]  
                     else:
                        """rand=random.uniform(-local_cycle_arg_dict["current_parameters"][parameter]["max_d"], local_cycle_arg_dict["current_parameters"][parameter]["max_d"])
                        relative_value=(local_cycle_arg_dict["current_parameters"][parameter]["v"]-lb)#*rand
                        delta_parameter=max(relative_value,local_cycle_arg_dict["min_perturbation"])*rand #If the parameter is to close to 0 we want to use the value of min_perturbation instead of the parameter value
                        delta_parameter=max(min(delta_parameter*rand,local_cycle_arg_dict["max_perturbation"]),-local_cycle_arg_dict["max_perturbation"]) #Ensure that the perturbation is not larger than max_perturbation
                        value=(lb+relative_value+delta_parameter)"""
                        rand=random.uniform(1-local_cycle_arg_dict["current_parameters"][parameter]["max_d"], 1+local_cycle_arg_dict["current_parameters"][parameter]["max_d"])
                        relative_value=(local_cycle_arg_dict["current_parameters"][parameter]["v"]-lb)
                        if relative_value<local_cycle_arg_dict["min_perturbation"]:
                           relative_value=local_cycle_arg_dict["min_perturbation"]
                        relative_value*=rand 
                        value=(lb+relative_value)
                        #value=min(max(value,lb+0.5*local_cycle_arg_dict["parameter_precision"]),ub-0.5*local_cycle_arg_dict["parameter_precision"])
                        value=min(max(value,lb),ub)
                        #value=min(max(value,round_up(lb,precision)),round_down(ub,precision)) #############################Canviat######################################################################################33
                        #print [parameter,delta_parameter,relative_value]
                        #value=min(max(value,lb+0.5*local_cycle_arg_dict["parameter_precision"]),ub-0.5*local_cycle_arg_dict["parameter_precision"])
                        #value=min(max(value,lb),ub)
                        #value=round(max(min(value,ub-0.5*local_cycle_arg_dict["parameter_precision"]),lb+0.5*local_cycle_arg_dict["parameter_precision"]),precision)
                        #value=round(max(min(value,ub),lb),precision)
                        #3value=max(value,(original_lb+0.5*local_cycle_arg_dict["parameter_precision"]))
                        #value=min(value,(original_ub-0.5*local_cycle_arg_dict["parameter_precision"]))
                        working_parameters[parameter]["v"] =value 
                        apply_parameters(local_cycle_arg_dict["label_model"],working_parameters,apply_flux_values= True,parameter_precision=local_cycle_arg_dict["parameter_precision"],parameter_list=[parameter])
                        if debug:
                           reaction=model.reactions.get_by_id(parameter)
                           print [fva,[reaction.lower_bound,reaction.upper_bound],[lb,ub]]
                           print model.optimize(tolerance_feasibility=local_cycle_arg_dict["label_model"].lp_tolerance_feasibility)
                        """for reaction_id in current_parameters[parameter]["reactions"]:  
                         model.reactions.get_by_id(reaction_id).lower_bound= round(working_parameters[parameter]["v"]-0.5*parameter_precision,precision+1) 
                         model.reactions.get_by_id(reaction_id).upper_bound= round(working_parameters[parameter]["v"]+0.5*parameter_precision,precision+1)"""
                     #print [model.reactions.get_by_id(reaction_id),model.reactions.get_by_id(reaction_id).lower_bound,model.reactions.get_by_id(reaction_id).upper_bound]  
        if failure_flag==False: #Do not run the simulation if the mdoel was unfeasible when ratios where applied
            #fi,working_parameters=coordinate_descent(label_model,mode="fsolve",hot_start=hot_start,parameter_precision=parameter_precision,parameter_to_be_fitted=working_parameters,max_perturbation=max_perturbation,original_reaction_bounds_dict=original_reaction_bounds_dict,perturbation=1.2)
            #print [fi, len(working_parameters)]
            status,b=solver(local_cycle_arg_dict["label_model"],mode=local_cycle_arg_dict["mode"],fba_mode=local_cycle_arg_dict["fba_mode"])
            fi,b,c=get_objective_function(local_cycle_arg_dict["label_model"])
            #print status
            #print label_model.constrained_model.optimize()
        else:
           fi=99999999
           #print "Code 1" 
        DeltaE = abs(fi-local_cycle_arg_dict["fc"])
        #print fi-local_cycle_arg_dict["fc"]
        if (fi>local_cycle_arg_dict["fc"]):
            # Initialize DeltaE_avg if a worse solution was found
            #   on the first iteration
            if (local_cycle_arg_dict["i"]==0 and j==0): 
                local_cycle_arg_dict["DeltaE_avg"] = min(DeltaE,50)
                #print "set avg"
                """if DeltaE>100:
                   DeltaE_avg=100"""
            
            # objective function is worse
            # generate probability of acceptance
            try :
              p = math.exp(-DeltaE/(local_cycle_arg_dict["DeltaE_avg"] * local_cycle_arg_dict["t"]))
              
            except:
              p=1
            # determine whether to accept worse point
            if (random.random()<p):
                # accept the worse solution
                accept = True
                if debug:
                   print [j,fi,fi-local_cycle_arg_dict["fc"],"accept",p]
            else:
                # don't accept the worse solution
                accept = False
        else:
            # objective function is lower, automatically accept
            accept = True
            if debug:
               print [j,fi,fi-local_cycle_arg_dict["fc"],"accept"]
            #check if it is the new best result
            
        if (accept==True):
            local_cycle_arg_dict["current_parameters"]=copy.deepcopy(working_parameters)
            ###error_dump(local_cycle_arg_dict["label_model"],"not optimal 1","accept_no_apply") 
            ###apply_parameters(local_cycle_arg_dict["label_model"], local_cycle_arg_dict["current_parameters"],apply_flux_values=True,parameter_precision=local_cycle_arg_dict["parameter_precision"])
            ###print ["Accept",model.optimize(tolerance_feasibility=local_cycle_arg_dict["label_model"].lp_tolerance_feasibility)]
            if model.optimize(tolerance_feasibility=local_cycle_arg_dict["label_model"].lp_tolerance_feasibility).status!="optimal":
               error_dump(local_cycle_arg_dict["label_model"],"not optimal 1","accept_apply")
            #print label_model.flux_dict
            #apply_parameters(label_model,working_parameters,apply_flux_values=True,parameter_precision=parameter_precision)
            #solver(label_model)
            #print label_model.flux_dict
            #ftemp,b,c=get_objective_function(label_model)
            #print "accept"
            #print model.optimize()
            """if ftemp!=fi:
               print [ftemp,fi]
               print label_model.flux_dict
               return"""
            #fi,working_parameters=coordinate_descent(label_model,mode="fsolve",hot_start=hot_start,parameter_precision=parameter_precision,parameter_to_be_fitted=working_parameters,max_perturbation=max_perturbation,original_reaction_bounds_dict=original_reaction_bounds_dict,perturbation=1.2)
            #print [fi, len(working_parameters)]
            if fi<local_cycle_arg_dict["f_best"]:
               best_parameters=copy.deepcopy(working_parameters)
               best_flux_dict=copy.deepcopy(local_cycle_arg_dict["label_model"].flux_dict) 
               local_cycle_arg_dict["f_best"]=fi
               local_cycle_arg_dict["best_parameters"]=best_parameters
               local_cycle_arg_dict["best_flux_dict"]=best_flux_dict
               print("New best fit is %s"%(fi))
            #Remove me
            #print label_model.flux_dict
            """for x in label_model.constrained_model.reactions:
                if x.upper_bound!=1000:
                   print [x.id,x.upper_bound,x.lower_bound]"""
            """solver(label_model)
            ftemp0,b,c=get_objective_function(label_model)
            apply_parameters(label_model,working_parameters,apply_flux_values=True,parameter_precision=parameter_precision)
            solver(label_model)
            #print label_model.flux_dict
            ftemp,b,c=get_objective_function(label_model)"""
            #print [fi,ftemp0,ftemp]
            #print working_parameters
            """for x in label_model.constrained_model.reactions:
                if x.upper_bound!=1000:
                   print [x.id,x.upper_bound,x.lower_bound]"""
            #Remove me    
            #fi2,b=get_objective_function(label_model) #Remova
            #print[fi,fi2]
            # update currently accepted solution
            #local_cycle_arg_dict["current_parameters"]=copy.deepcopy(working_parameters)
            #print current_parameters
            #apply_parameters(label_model,current_parameters,apply_flux_values=True)
            local_cycle_arg_dict["fc"]= fi
            
            # increment number of accepted solutions
            local_cycle_arg_dict["na"]+=1
            na=local_cycle_arg_dict["na"]
            #print na
            # update DeltaE_avg
            local_cycle_arg_dict["DeltaE_avg"] = (local_cycle_arg_dict["DeltaE_avg"] * (na-1.0) +  DeltaE) / na
            #print [na,i,=local_cycle_arg_dict["t"],DeltaE,DeltaE_avg]
        else:
            if debug:
               print [j,fi,fi-local_cycle_arg_dict["fc"],"reject"]
            #If we rejected the simulation we need to restore the flux bound to the value of current_parameters
            #print "Solution Rejected"
            #print "rejected" 
            #apply_parameters(local_cycle_arg_dict["label_model"], local_cycle_arg_dict["current_parameters"],apply_flux_values=True,parameter_precision=local_cycle_arg_dict["parameter_precision"],parameter_list=working_parameters_sample)
            apply_parameters(local_cycle_arg_dict["label_model"], local_cycle_arg_dict["current_parameters"],apply_flux_values=True,parameter_precision=local_cycle_arg_dict["parameter_precision"],parameter_list=working_parameters_sample)
            ###print ["reject",model.optimize(tolerance_feasibility=local_cycle_arg_dict["label_model"].lp_tolerance_feasibility)]
            """for parameter in sampled_flux_value_parameters: 
                   for reaction_id in current_parameters[parameter]["reactions"]:
                       reaction=model.reactions.get_by_id(reaction_id)
                       reaction.lower_bound=round(current_parameters[parameter]["v"]-0.5*parameter_precision,precision+1)
                       reaction.upper_bound=round(current_parameters[parameter]["v"]+0.5*parameter_precision,precision+1)"""
    #apply_parameters(local_cycle_arg_dict["label_model"], local_cycle_arg_dict["best_parameters"],apply_flux_values=True,parameter_precision=local_cycle_arg_dict["parameter_precision"]) 
    f_best1=f_best2=9999999
    if original_f_best!=local_cycle_arg_dict["f_best"]:
       f_best1, parameter_dict1,best_flux_dict1=coordinate_descent(local_cycle_arg_dict["label_model"],parameter_dict=local_cycle_arg_dict["best_parameters"],parameter_to_be_fitted=local_cycle_arg_dict["parameter_to_be_fitted"],parameter_precision=local_cycle_arg_dict["parameter_precision"],perturbation=0.01,max_perturbation=local_cycle_arg_dict["max_perturbation"],force_flux_value_bounds=False,mode="fsolve",fba_mode="fba",debug=False,is_subprocess=True)
    if local_cycle_arg_dict["fc"]!=local_cycle_arg_dict["f_best"]:
       f_best2, parameter_dict2,best_flux_dict2=coordinate_descent(local_cycle_arg_dict["label_model"],parameter_dict=local_cycle_arg_dict["current_parameters"],parameter_to_be_fitted=local_cycle_arg_dict["parameter_to_be_fitted"],parameter_precision=local_cycle_arg_dict["parameter_precision"],perturbation=0.01,max_perturbation=local_cycle_arg_dict["max_perturbation"],force_flux_value_bounds=False,mode="fsolve",fba_mode="fba",debug=False,is_subprocess=True)
    if  f_best1<f_best2:
       print ["coordinated1"]
       local_cycle_arg_dict["f_best"]=f_best1
       local_cycle_arg_dict["fc"]=f_best1
       local_cycle_arg_dict["best_parameters"]=copy.deepcopy(parameter_dict1)
       local_cycle_arg_dict["best_flux_dict"]=copy.deepcopy(best_flux_dict1)
       local_cycle_arg_dict["current_parameters"]=copy.deepcopy(parameter_dict1)
    elif f_best2<f_best1 and f_best2<original_f_best:
       print ["coordinated2"]
       local_cycle_arg_dict["f_best"]=f_best2
       local_cycle_arg_dict["fc"]=f_best2
       local_cycle_arg_dict["best_parameters"]=copy.deepcopy(parameter_dict2)
       local_cycle_arg_dict["best_flux_dict"]=copy.deepcopy(best_flux_dict2)
       local_cycle_arg_dict["current_parameters"]=copy.deepcopy(parameter_dict2)
    """print ["current_f_best",local_cycle_arg_dict["f_best"]]
    if original_f_best!=local_cycle_arg_dict["f_best"]:
       f_best1, parameter_dict1,best_flux_dict1=coordinate_descent(local_cycle_arg_dict["label_model"],parameter_dict=local_cycle_arg_dict["best_parameters"],parameter_to_be_fitted=local_cycle_arg_dict["parameter_to_be_fitted"],parameter_precision=local_cycle_arg_dict["parameter_precision"],perturbation=0.01,max_perturbation=local_cycle_arg_dict["max_perturbation"],force_flux_value_bounds=False,mode="fsolve",fba_mode="fba",debug=False,is_subprocess=True)
       local_cycle_arg_dict["f_best"]=f_best1
       local_cycle_arg_dict["fc"]=f_best1
       local_cycle_arg_dict["best_parameters"]=copy.deepcopy(parameter_dict1)
       local_cycle_arg_dict["best_flux_dict"]=copy.deepcopy(best_flux_dict1)
       print ["best",f_best1]
       #local_cycle_arg_dict["current_parameters"]=copy.deepcopy(parameter_dict)
    else:
      f_best1=original_f_best   
    f_best2, parameter_dict2,best_flux_dict2=coordinate_descent(local_cycle_arg_dict["label_model"],parameter_dict=local_cycle_arg_dict["current_parameters"],parameter_to_be_fitted=local_cycle_arg_dict["parameter_to_be_fitted"],parameter_precision=local_cycle_arg_dict["parameter_precision"],perturbation=0.01,max_perturbation=local_cycle_arg_dict["max_perturbation"],force_flux_value_bounds=False,mode="fsolve",fba_mode="fba",debug=False,is_subprocess=True)
    print ["current",f_best2]
    if f_best2<f_best1:
       print ["accept current",f_best2] 
       local_cycle_arg_dict["f_best"]=f_best2
       local_cycle_arg_dict["fc"]=f_best2
       local_cycle_arg_dict["best_parameters"]=copy.deepcopy(parameter_dict2)
       local_cycle_arg_dict["best_flux_dict"]=copy.deepcopy(best_flux_dict2)
    local_cycle_arg_dict["current_parameters"]=copy.deepcopy(local_cycle_arg_dict["best_parameters"])"""
    return local_cycle_arg_dict

