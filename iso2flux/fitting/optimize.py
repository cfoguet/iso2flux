import numpy as np
import copy
import pygmo
#from PyGMO import algorithm, island,archipelago, topology
import gc
from objfunc import objfunc
from get_variable_bounds import get_variable_bounds
from extract_results import extract_results
from ..flux_functions.minimal_flux import create_minimal_fux_model,flux_minimization_fva
from variable_sampling import variable_sampling
from multiprocessing import Pool
"""
except: 
  from iso2flux.flux_functions.minimal_flux import create_minimal_fux_model,flux_minimization_fva
  from iso2flux.fitting.variable_sampling import variable_sampling
"""
 

import time

from random import shuffle



def migrate_one_direction(fs,xs,n_bests,verbose=False):
     #Migrate only in one direction
     if verbose:
        print "-------------------\nMigration:"
     champions_f= [min(f) for f in fs]
     champions_x=[]
     for n,n_best in enumerate(n_bests):
         champions_x.append(xs[n][n_best]) 
     variables_sets=xs
     n_archis=range(0,len(variables_sets))
     #For each island check if it can migrate to the next island
     for n_archi in n_archis:
            best_individual_f=champions_f[n_archi]
            best_individual_x=list(champions_x[n_archi])
            #Get the worst element of the adjacents islands
            worst_f_left=None
            n_left_island=n_archis[n_archi-1]
            for n_left,f_left in enumerate(fs[n_left_island]):
                if worst_f_left==None or worst_f_left>f_left:
                   worst_n_left=n_left
                   worst_f_left=f_left
            
            #If the best element of the island is better than the worst of the adjacent islands replace them
            if best_individual_f<worst_f_left:
               if verbose:
                  print [n_archi,"left",[n_archi,"->",n_left_island],[best_individual_f,"->",worst_f_left]] 
               fs[n_left_island][worst_n_left]=best_individual_f
               variables_sets[n_left_island][worst_n_left]=copy.deepcopy(best_individual_x)
               
            """#Do the same for the right island
            if (n_archi+1)==len(variables_sets):
               n_right_island=0
            else:
               n_right_island=n_archis[n_archi+1]
            worst_f_right=None
            for n_right,f_right in enumerate(fs[n_right_island]):
                if worst_f_right==None or worst_f_right>f_right:
                   worst_n_right=n_right
                   worst_f_right=f_right
            
            #If the best element of the island is better than the worst of the adjacent islands replace them
            if best_individual_f<worst_f_right:
               if verbose: 
                  print [n_archi,"right",[n_archi,"->",n_right_island],[best_individual_f,"->",worst_f_right]] 
               fs[n_right_island][worst_n_right]=best_individual_f
               variables_sets[n_right_island][worst_n_right]=copy.deepcopy(best_individual_x)"""
     #print variables_sets
     print "-------------------"    
     return variables_sets


def migrate_ring(fs,xs,n_bests,verbose=False):
     if verbose:
        print "-------------------\nMigration:"
     champions_f= [min(f) for f in fs]
     champions_x=[]
     for n,n_best in enumerate(n_bests):
         champions_x.append(xs[n][n_best]) 
     variables_sets=xs
     n_archis=range(0,len(variables_sets))
     #For each island check if it can migrate to the next island
     for n_archi in n_archis:
            best_individual_f=champions_f[n_archi]
            best_individual_x=list(champions_x[n_archi])
            #Get the worst element of the adjacents islands
            worst_f_left=None
            n_left_island=n_archis[n_archi-1]
            for n_left,f_left in enumerate(fs[n_left_island]):
                if worst_f_left==None or worst_f_left>f_left:
                   worst_n_left=n_left
                   worst_f_left=f_left
            
            #If the best element of the island is better than the worst of the adjacent islands replace them
            if best_individual_f<worst_f_left:
               if verbose:
                  print [n_archi,"left",[n_archi,"->",n_left_island],[best_individual_f,"->",worst_f_left]] 
               fs[n_left_island][worst_n_left]=best_individual_f
               variables_sets[n_left_island][worst_n_left]=copy.deepcopy(best_individual_x)
               
            
            #Do the same for the right island
            if (n_archi+1)==len(variables_sets):
               n_right_island=0
            else:
               n_right_island=n_archis[n_archi+1]
            worst_f_right=None
            for n_right,f_right in enumerate(fs[n_right_island]):
                if worst_f_right==None or worst_f_right>f_right:
                   worst_n_right=n_right
                   worst_f_right=f_right
            
            #If the best element of the island is better than the worst of the adjacent islands replace them
            if best_individual_f<worst_f_right:
               if verbose: 
                  print [n_archi,"right",[n_archi,"->",n_right_island],[best_individual_f,"->",worst_f_right]] 
               fs[n_right_island][worst_n_right]=best_individual_f
               variables_sets[n_right_island][worst_n_right]=copy.deepcopy(best_individual_x)
     #print variables_sets
     print "-------------------"    
     return variables_sets

"""def migrate_ring(archi):
     champions_f=archi.get_champions_f()
     champions_x=archi.get_champions_x()
     n_archis=range(0,len(archi))
     f_sets=[]
     variables_sets=[]
     if n_archis==1:
        variables_sets=[archi[0].get_population().get_x()]
        return variables_sets
     for island in archi:
            island_variable_set=[]
            islands_x=island.get_population().get_x()
            island_f= island.get_population().get_f() 
            f_sets.append(island_f)
            #variables_sets.append(islands_x)     
            for x in islands_x:
                island_variable_set.append(x)
            variables_sets.append(island_variable_set)
     #For each island check if it can migrate to the next island
     for n_archi in n_archis:
            best_individual_f=champions_f[n_archi]
            best_individual_x=list(champions_x[n_archi])
            #Get the worst element of the adjacents islands
            worst_f_left=None
            n_left_island=n_archis[n_archi-1]
            for n_left,f_left in enumerate(f_sets[n_left_island]):
                if worst_f_left==None or worst_f_left>f_left:
                   worst_n_left=n_left
                   worst_f_left=f_left
            
            #If the best element of the island is better than the worst of the adjacent islands replace them
            if best_individual_f<worst_f_left:
               f_sets[n_left_island][worst_n_left]=best_individual_f
               variables_sets[n_left_island][worst_n_left]=copy.deepcopy(best_individual_x)
            
            #Do the same for the right island
            if (n_archi+1)==len(archi):
               n_right_island=0
            else:
               n_right_island=n_archis[n_archi+1]
            worst_f_right=None
            for n_right,f_right in enumerate(f_sets[n_right_island]):
                if worst_f_right==None or worst_f_right>f_right:
                   worst_n_right=n_right
                   worst_f_right=f_right
            
            #If the best element of the island is better than the worst of the adjacent islands replace them
            if best_individual_f<worst_f_right:
               f_sets[n_right_island][worst_n_right]=best_individual_f
               variables_sets[n_right_island][worst_n_right]=copy.deepcopy(best_individual_x)
     #print variables_sets          
     return variables_sets
"""

def evolve_process(pop): 
      algo = pygmo.algorithm(pygmo.sade(gen = 300)) 
      pop=algo.evolve(pop)
      return [pop.get_f(),pop.get_x(),pop.best_idx()] 


def optimize(label_model,iso2flux_problem,pop_size = 25,n_gen = 500,n_islands=6,max_evolve_cycles=999,max_cycles_without_improvement=10,stop_criteria_relative=0.01,stop_criteria_absolute=-1e6,initial_archi_x=[],lb_list=[],ub_list=[],flux_penalty_dict=None,max_flux=None,label_problem_parameters={},min_model=None,extra_constraint_dict={},log_file_name="optimize_log.txt",max_flux_sampling=None,migrate="ring"):
  log_file=open(log_file_name,"w")
  log_file.write(time.strftime("%c")+": starting optimization\n")
  log_file.close()  
  if flux_penalty_dict==None or flux_penalty_dict=={}:
          flux_penalty_dict={}
          if "flux_penalty_dict" in  label_problem_parameters:
               flux_penalty_dict=label_problem_parameters["flux_penalty_dict"]
  if  max_flux==None:
           max_flux=1e6
           if "max_flux" in  label_problem_parameters:
               max_flux=label_problem_parameters["max_flux"]
  input_flag=False
  if min_model==None:
       min_model,minimal_flux,non_milp_reactions =create_minimal_fux_model(copy.deepcopy(label_model.constrained_model),fraction_of_optimum_objective=0.0,      fraction_of_flux_minimum=None,boundaries_precision=0.001,label_model=None,metabolite_list_file_name=None,flux_penalty_dict=flux_penalty_dict,maximum_flux=max_flux, mutually_exclusive_directionality_constraint=False,extra_constraints_dict=extra_constraint_dict) 
  else:
       min_model.reactions.get_by_id("total_flux").upper_bound=max_flux
       min_model.reactions.get_by_id("total_flux").lower_bound=0
  
  #label_problem=label_model.iso2flux_problem(dim=len(label_model.variable_vector))
  with min_model:
    if lb_list==[] or ub_list==[]:     
       lb_list,ub_list=get_variable_bounds(label_model,flux_penalty_dict=flux_penalty_dict,max_flux=max_flux,min_model=min_model)
       for n,x in enumerate(lb_list): #If upper and lower bound are identical program will get stuck creating problems
           if lb_list[n]-ub_list[n]<1e-9:
              ub_list[n]+=1e-6
  prob=pygmo.problem(iso2flux_problem(lb_list,ub_list,parameter_dict=label_problem_parameters,verbose=False))
  if "target_flux_unfeasible_penalty" not in label_problem_parameters:
      label_problem_parameters["target_flux_unfeasible_penalty"]=10
  uncont=pygmo.unconstrain(prob = prob, method = "weighted",weights=[label_problem_parameters["label_unfeasible_penalty"],label_problem_parameters["flux_unfeasible_penalty"],label_problem_parameters["target_flux_unfeasible_penalty"],label_problem_parameters["target_flux_unfeasible_penalty"]])
  #algo = pygmo.algorithm(pygmo.sade(gen = n_gen))
  algo = pygmo.algorithm(pygmo.sade(gen = n_gen)) 
  
  #variables_sets=variable_sampling(label_model,lb_list,ub_list,maximum_flux=max_flux,flux_penalty_dict=flux_penalty_dict,n_pop=pop_size*n_islands,n_processes=1)
  if max_flux_sampling==None:
     max_flux_sampling= max_flux  
  sampled_variables=variable_sampling(label_model,lb_list,ub_list,maximum_flux=max_flux_sampling,flux_penalty_dict=flux_penalty_dict,n_pop=pop_size*n_islands,n_processes=1,extra_constraints_dict=extra_constraint_dict)
  variables_sets=[]
  counter=0
  for n_island in range(0,n_islands):
      
      island_variables=[]
      for n_pop in range(0,pop_size):
          island_variables.append(sampled_variables[counter])
          counter+=1
      variables_sets.append(island_variables)
  count_no_progress=0
  pool = Pool(processes=n_islands)
  for ncycle in range(0,max_evolve_cycles):
        
        pop_list=[]
        for n_island in range(0,n_islands):
            pop = pygmo.population(uncont, size = pop_size)
            for n_pop in range(0,pop_size):
                   pop.set_x(n_pop,variables_sets[n_island][n_pop])
            pop_list.append(pop)
        champions=[x.champion_f[0] for x in pop_list]
        if ncycle==0:
           previous_best_solution=None
           for n_island, obj in enumerate(champions):
            f=obj
            if f<previous_best_solution or n_island==0:
                n_optimal=n_island
                previous_best_solution=f
                previous_best_variables=pop_list[n_optimal].champion_x
        if True:
               print len(pop_list)
               results=pool.map(evolve_process,pop_list)
               #pool.close()
               #variables_sets=[x[1] for x in results]
               fs=[[x2[0] for x2 in x1[0]] for x1 in results]
               xs=[x[1] for x in results]
               best_ns=[x[2] for x in results]
               champions=[min(x) for x in fs]
               optimal_solution=min(champions)
               for n_champ,champion in enumerate(champions):
                   if champion==optimal_solution:
                      optimal_variables=xs[n_champ][best_ns[n_champ]]
                      
               #variables_sets=migrate_ring(fs,xs,best_ns)
        log_file=open(log_file_name,"a")
        log_file.write("\n"+time.strftime("%c")+":\n")       
        for n_pop,n_best in enumerate(best_ns):
            x=xs[n_pop][n_best]
            obj, obj_dict=objfunc(label_model,x,verbose=False,label_weight=label_problem_parameters["label_weight"],flux_weight=label_problem_parameters["flux_weight"],max_chi=label_problem_parameters["max_chi"],target_flux_dict=label_problem_parameters["target_flux_dict"],max_flux=label_problem_parameters["max_flux"],flux_penalty_dict=label_problem_parameters["flux_penalty_dict"])
            
            flux_obj=obj_dict["flux_score"]
            label_obj=obj_dict["chi2_score"]
            fltarget_obj=obj_dict["fltarget_obj"]
            output="Island "+str(n_pop)+ ": obj="+str(round(obj,3))
            output+=" Chi2="+str(round(label_obj,3))
            if label_problem_parameters["max_chi"]<99999:
                output+="/"+str(round(label_problem_parameters["max_chi"],3))
            if flux_obj>0:
                output+=" flux="+str(round(flux_obj,3))+"/"+str(label_problem_parameters["max_flux"])
            target_flux_dict=label_problem_parameters["target_flux_dict"]
            if fltarget_obj>0 and target_flux_dict!={}:
               irreversible_model_flag=target_flux_dict.get("irreversible_model")
               if irreversible_model_flag in (None,False):
                  flux=label_model.reversible_flux_dict[target_flux_dict["reaction"]]
               else:
                  flux=label_model.flux_dict[target_flux_dict["reaction"].replace("_forward","")]
               output+=" "+target_flux_dict["reaction"]+"="+str(round(flux,3))+" "+target_flux_dict["dir"]
            print output
            log_file.write(output+"\n")
        """for n_island, obj in enumerate(champions):
            f=obj[0]
            if f<optimal_solution:
                n_optimal=n_island
                optimal_solution=f"""
        #print archi.get_champions_f()
        print "previous best objective: "+str(round(previous_best_solution,3))+"; current best objective: "+str(round(optimal_solution,3))
        log_file.write("previous best objective: "+str(round(previous_best_solution,3))+"; current best objective: "+str(round(optimal_solution,3))+"\n")
        #optimal_variables=champions_x[n_optimal]
        #variables_sets=migrate_ring(archi)
        #variables_sets=migrate_ring(fs,xs,best_ns)
        #shuffle(variables_sets)
        #print previous_best_solution,optimal_solution
        #print optimal_variables
        log_file.write(str(optimal_variables)+"\n")
        log_file.close()
        if previous_best_solution-optimal_solution<max(previous_best_solution*stop_criteria_relative,1e-4) or optimal_solution<stop_criteria_absolute:
           if previous_best_solution<optimal_solution:
              optimal_variables=previous_best_variables
              optimal_solution=previous_best_solution
           count_no_progress+=1
           
           print "no signficant improvement"+ str(count_no_progress)+"/"+str(max_cycles_without_improvement)
           if count_no_progress>max_cycles_without_improvement:
               print "convergence reached at cycle "+str(ncycle)
               break
        else:
           #print [previous_best_solution,optimal_solution]
           print "improvement"
           count_no_progress=0
           previous_best_solution=optimal_solution
           previous_best_variables=optimal_variables
        if optimal_solution<1e-6: #Mostly used when computing intervals
           break
        if migrate=="one_direction":
            variables_sets=migrate_one_direction(fs,xs,best_ns,False)
        else:       
           variables_sets=migrate_ring(fs,xs,best_ns)
  
  obj, obj_dict=objfunc(label_model,optimal_variables,verbose=False,max_chi=label_problem_parameters["max_chi"],target_flux_dict=label_problem_parameters["target_flux_dict"],max_flux=label_problem_parameters["max_flux"],flux_penalty_dict=label_problem_parameters["flux_penalty_dict"])     
  if (obj_dict["chi2_score"]-label_problem_parameters["max_chi"])>0.0001 or (obj_dict["flux_score"]-label_problem_parameters["max_flux"])>0.0001:
     optimal_solution=1e16
     print obj_dict
  if label_problem_parameters["target_flux_dict"] not in [None,{}]:
     if irreversible_model_flag in (None,False):
              flux=label_model.reversible_flux_dict[target_flux_dict["reaction"]]
     else:
              flux=label_model.flux_dict[target_flux_dict["reaction"].replace("_forward","")]
     if flux-label_problem_parameters["target_flux_dict"]["ub"]>1e-6:
        optimal_solution=1e16
     if flux-label_problem_parameters["target_flux_dict"]["lb"]<1e-6:
        optimal_solution=1e16
  pool.close() 
  return optimal_solution,optimal_variables


 
"""               
def fake():               
        print "Cycle "+str(ncycle)
        input_flag=False
        counter=0
        if ncycle==0:
           archi=pygmo.archipelago() #empty archi
           for n_island in range(0,n_islands): 
               #print "island"
               pop = pygmo.population(uncont, size = 0)
               for n_pop in range(0,pop_size):
                   #print n_pop
                   pop.push_back(variables_sets[n_island][n_pop])
                   #pop.push_back(variables_sets[counter])
                   counter+=1
               archi.push_back(pop=pop,algo=algo)
           previous_best_solution=None
           for n_island, obj in enumerate(archi.get_champions_f()):
            f=obj[0]
            if f<previous_best_solution or n_island==0:
                n_optimal=n_island
                previous_best_solution=f
           previous_best_variables=archi.get_champions_x()[n_optimal]
        print  gc.garbage
        if input_flag: yourvar = input('Choose a number:1 ') 
        #return archi
        archi.evolve()
        if input_flag: yourvar = input('Choose a number:2 ') 
        #archi.wait()
        archi.wait_check()
        if input_flag: yourvar = input('Choose a number:3 ') 
        #print aaaaaaaaaaaaaaaa
        if input_flag: yourvar = input('Choose a number:3.1 ') 
        optimal_solution=1e12
        print "champions:"
        champions=archi.get_champions_f()
        champions_x=archi.get_champions_x()
        if input_flag: yourvar = input('Choose a number:3.2 ') 
        for n,x in enumerate(champions_x):
            obj, obj_dict=objfunc(label_model,x,verbose=False,max_chi=label_problem_parameters["max_chi"],target_flux_dict=label_problem_parameters["target_flux_dict"],max_flux=label_problem_parameters["max_flux"],flux_penalty_dict=label_problem_parameters["flux_penalty_dict"])
            
            flux_obj=obj_dict["flux_score"]
            label_obj=obj_dict["chi2_score"]
            fltarget_obj=obj_dict["fltarget_obj"]
            output="Island "+str(n)+ ": obj="+str(round(champions[n][0],3))
            output+=" Chi2="+str(round(label_obj,3))
            if flux_obj>0:
                output+=" flux="+str(round(flux_obj,3))
            target_flux_dict=label_problem_parameters["target_flux_dict"]
            if fltarget_obj>0 and target_flux_dict!={}:
               output+=" "+target_flux_dict["reaction"]+"="+str(round(label_model.reversible_flux_dict[target_flux_dict["reaction"]],3))+" "+target_flux_dict["dir"]
            print output
        if input_flag: yourvar = input('Choose a number:3.3 ') 
        for n_island, obj in enumerate(champions):
            f=obj[0]
            if f<optimal_solution:
                n_optimal=n_island
                optimal_solution=f
        #print archi.get_champions_f()
        
        print "previous best objective: "+str(round(previous_best_solution,3))+"; current best objective: "+str(round(optimal_solution,3))
        optimal_variables=champions_x[n_optimal]
        variables_sets=migrate_ring(archi)
        if input_flag: yourvar = input('Choose a number: 4') 
        #shuffle(variables_sets)
        #print previous_best_solution,optimal_solution
        if previous_best_solution-optimal_solution<max(previous_best_solution*stop_criteria_relative,1e-4) or optimal_solution<stop_criteria_absolute:
           if previous_best_solution<optimal_solution:
              optimal_variables=previous_best_variables
              optimal_solution=previous_best_solution
           count_no_progress+=1
           
           print "no signficant improvement"+ str(count_no_progress)+"/"+str(max_cycles_without_improvement)
           if count_no_progress>max_cycles_without_improvement:
               print "convergence reached at cycle "+str(ncycle)
               break
        else:
           #print [previous_best_solution,optimal_solution]
           print "improvement"
           count_no_progress=0
           previous_best_solution=optimal_solution
           previous_best_variables=optimal_variables
        if optimal_solution<1e-6: #Mostly used when computing intervals
           break
        if input_flag: yourvar = input('Choose a number: 5') 
        print "------------------------------------"
  
  obj, obj_dict=objfunc(label_model,optimal_variables,verbose=False,max_chi=label_problem_parameters["max_chi"],target_flux_dict=label_problem_parameters["target_flux_dict"],max_flux=label_problem_parameters["max_flux"],flux_penalty_dict=label_problem_parameters["flux_penalty_dict"])     
  if (obj_dict["chi2_score"]-label_problem_parameters["max_chi"])>0.0001 or (obj_dict["flux_score"]-label_problem_parameters["max_flux"])>0.0001:
     optimal_solution=1e99
     print obj_dict
  if label_problem_parameters["target_flux_dict"] not in [None,{}]:
     if label_model.reversible_flux_dict[self.target_flux_dict["reaction"]]-label_problem_parameters["target_flux_dict"]["ub"]>1e-6:
        optimal_solution=1e99
     if label_model.reversible_flux_dict[self.target_flux_dict["reaction"]]-label_problem_parameters["target_flux_dict"]["lb"]<1e-6:
        optimal_solution=1e99
  return optimal_solution,optimal_variables 

"""








"""

def optimize(label_model,pop_size = 25,n_gen = 500,n_islands=6,max_evolve_cycles=999,max_cycles_without_improvement=10,stop_criteria_relative=0.01,stop_criteria_absolute=-1e6,initial_archi_x=[],lb_list=[],ub_list=[],flux_penalty_dict=None,max_flux=None,label_problem_parameters={},min_model=None,extra_constraint_dict={}):
  if flux_penalty_dict==None or flux_penalty_dict=={}:
          flux_penalty_dict={}
          if "flux_penalty_dict" in  label_problem_parameters:
               flux_penalty_dict=label_problem_parameters["flux_penalty_dict"]
  if  max_flux==None:
           max_flux=1e6
           if "max_flux" in  label_problem_parameters:
               max_flux=label_problem_parameters["max_flux"]
  if min_model==None:
       min_model,minimal_flux,non_milp_reactions =create_minimal_fux_model(copy.deepcopy(label_model.constrained_model),fraction_of_optimum_objective=0.0,      fraction_of_flux_minimum=None,boundaries_precision=0.001,label_model=None,metabolite_list_file_name=None,flux_penalty_dict=flux_penalty_dict,maximum_flux=max_flux, mutually_exclusive_directionality_constraint=False,extra_constraints_dict=extra_constraint_dict) 
  else:
       min_model.reactions.get_by_id("total_flux").upper_bound=max_flux
       min_model.reactions.get_by_id("total_flux").lower_bound=0
  
  #label_problem=label_model.iso2flux_problem(dim=len(label_model.variable_vector))
  with min_model:
    if lb_list==[] or ub_list==[]:     
       lb_list,ub_list=get_variable_bounds(label_model,flux_penalty_dict=flux_penalty_dict,max_flux=max_flux,min_model=min_model)
       for n,x in enumerate(lb_list): #If upper and lower bound are identical program will get stuck creating problems
           if lb_list[n]-ub_list[n]<1e-9:
              ub_list[n]+=1e-6
  prob=pygmo.problem(iso2flux_problem(label_model,lb_list,ub_list,parameter_dict=label_problem_parameters,verbose=False))
  if "target_flux_unfeasible_penalty" not in label_problem_parameters:
      label_problem_parameters["target_flux_unfeasible_penalty"]=10
  uncont=pygmo.unconstrain(prob = prob, method = "weighted",weights=[label_problem_parameters["label_unfeasible_penalty"],label_problem_parameters["flux_unfeasible_penalty"],label_problem_parameters["target_flux_unfeasible_penalty"],label_problem_parameters["target_flux_unfeasible_penalty"]])
  #algo = pygmo.algorithm(pygmo.sade(gen = n_gen))
  algo = pygmo.algorithm(pygmo.de(gen = n_gen)) 
  
  #variables_sets=variable_sampling(label_model,lb_list,ub_list,maximum_flux=max_flux,flux_penalty_dict=flux_penalty_dict,n_pop=pop_size*n_islands,n_processes=1)
  sampled_variables=variable_sampling(label_model,lb_list,ub_list,maximum_flux=max_flux,flux_penalty_dict=flux_penalty_dict,n_pop=pop_size*n_islands,n_processes=1,extra_constraints_dict=extra_constraint_dict)
  a=iso2flux_problem(label_model,lb_list,ub_list,parameter_dict=label_problem_parameters)
  a.fitness(sampled_variables[0])
  variables_sets=[]
  counter=0
  for n_island in range(0,n_islands):
      
      island_variables=[]
      for n_pop in range(0,pop_size):
          island_variables.append(sampled_variables[counter])
          counter+=1
      variables_sets.append(island_variables)
  #print    variables_sets               
  #jsjsjjs
  count_no_progress=0
  for ncycle in range(0,max_evolve_cycles):
        
        print "Cycle "+str(ncycle)
        input_flag=False
        counter=0
        if ncycle==0:
           archi=pygmo.archipelago() #empty archi
           for n_island in range(0,n_islands): 
               #print "island"
               pop = pygmo.population(uncont, size = 0)
               for n_pop in range(0,pop_size):
                   #print n_pop
                   pop.push_back(variables_sets[n_island][n_pop])
                   #pop.push_back(variables_sets[counter])
                   counter+=1
               archi.push_back(pop=pop,algo=algo)
           previous_best_solution=None
           for n_island, obj in enumerate(archi.get_champions_f()):
            f=obj[0]
            if f<previous_best_solution or n_island==0:
                n_optimal=n_island
                previous_best_solution=f
           previous_best_variables=archi.get_champions_x()[n_optimal]
        print  gc.garbage
        if input_flag: yourvar = input('Choose a number:1 ') 
        #return archi
        archi.evolve()
        if input_flag: yourvar = input('Choose a number:2 ') 
        #archi.wait()
        archi.wait_check()
        if input_flag: yourvar = input('Choose a number:3 ') 
        #print aaaaaaaaaaaaaaaa
        if input_flag: yourvar = input('Choose a number:3.1 ') 
        optimal_solution=1e12
        print "champions:"
        champions=archi.get_champions_f()
        champions_x=archi.get_champions_x()
        if input_flag: yourvar = input('Choose a number:3.2 ') 
        for n,x in enumerate(champions_x):
            obj, obj_dict=objfunc(label_model,x,verbose=False,max_chi=label_problem_parameters["max_chi"],target_flux_dict=label_problem_parameters["target_flux_dict"],max_flux=label_problem_parameters["max_flux"],flux_penalty_dict=label_problem_parameters["flux_penalty_dict"])
            
            flux_obj=obj_dict["flux_score"]
            label_obj=obj_dict["chi2_score"]
            fltarget_obj=obj_dict["fltarget_obj"]
            output="Island "+str(n)+ ": obj="+str(round(champions[n][0],3))
            output+=" Chi2="+str(round(label_obj,3))
            if flux_obj>0:
                output+=" flux="+str(round(flux_obj,3))
            target_flux_dict=label_problem_parameters["target_flux_dict"]
            if fltarget_obj>0 and target_flux_dict!={}:
               output+=" "+target_flux_dict["reaction"]+"="+str(round(label_model.reversible_flux_dict[target_flux_dict["reaction"]],3))+" "+target_flux_dict["dir"]
            print output
        if input_flag: yourvar = input('Choose a number:3.3 ') 
        for n_island, obj in enumerate(champions):
            f=obj[0]
            if f<optimal_solution:
                n_optimal=n_island
                optimal_solution=f
        #print archi.get_champions_f()
        
        print "previous best objective: "+str(round(previous_best_solution,3))+"; current best objective: "+str(round(optimal_solution,3))
        optimal_variables=champions_x[n_optimal]
        variables_sets=migrate_ring(archi)
        if input_flag: yourvar = input('Choose a number: 4') 
        #shuffle(variables_sets)
        #print previous_best_solution,optimal_solution
        if previous_best_solution-optimal_solution<max(previous_best_solution*stop_criteria_relative,1e-4) or optimal_solution<stop_criteria_absolute:
           if previous_best_solution<optimal_solution:
              optimal_variables=previous_best_variables
              optimal_solution=previous_best_solution
           count_no_progress+=1
           
           print "no signficant improvement"+ str(count_no_progress)+"/"+str(max_cycles_without_improvement)
           if count_no_progress>max_cycles_without_improvement:
               print "convergence reached at cycle "+str(ncycle)
               break
        else:
           #print [previous_best_solution,optimal_solution]
           print "improvement"
           count_no_progress=0
           previous_best_solution=optimal_solution
           previous_best_variables=optimal_variables
        if optimal_solution<1e-6: #Mostly used when computing intervals
           break
        if input_flag: yourvar = input('Choose a number: 5') 
        print "------------------------------------"
  
  obj, obj_dict=objfunc(label_model,optimal_variables,verbose=False,max_chi=label_problem_parameters["max_chi"],target_flux_dict=label_problem_parameters["target_flux_dict"],max_flux=label_problem_parameters["max_flux"],flux_penalty_dict=label_problem_parameters["flux_penalty_dict"])     
  if (obj_dict["chi2_score"]-label_problem_parameters["max_chi"])>0.0001 or (obj_dict["flux_score"]-label_problem_parameters["max_flux"])>0.0001:
     optimal_solution=1e99
     print obj_dict
  if label_problem_parameters["target_flux_dict"] not in [None,{}]:
     if label_model.reversible_flux_dict[self.target_flux_dict["reaction"]]-label_problem_parameters["target_flux_dict"]["ub"]>1e-6:
        optimal_solution=1e99
     if label_model.reversible_flux_dict[self.target_flux_dict["reaction"]]-label_problem_parameters["target_flux_dict"]["lb"]<1e-6:
        optimal_solution=1e99
  return optimal_solution,optimal_variables 

"""




def find_flux_grups(label_model,reaction_list=None):
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
    if reaction_list!=None:
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







def define_isoflux_problem(label_model):
  
  class iso2flux_problem:
    """
    iso2flux_problem
    """
    
    def __init__(self,lb_list,ub_list,allow_unfeasible_solutions=False,verbose=True,parameter_dict={}):#, dim=len(label_model.variable_vector),i_dim=0,n_obj=1,c_dim=0,c_ineq_dim=0,c_tol=1e-6):
        # First we call the constructor of the base class telling PyGMO
        # what kind of problem to expect ('dim' dimensions, 1 objective, 0 contraints etc.)
        self.allow_unfeasible_solutions=allow_unfeasible_solutions
        self.label_weight=1
        self.flux_weight=1
        self.target_flux_dict=None
        self.max_chi=1e16
        self.max_flux=1e16
        self.flux_penalty_dict={}
        self.verbose=verbose
        self.fmin=1e36
        self.label_unfeasible_penalty=1e9
        self.flux_unfeasible_penalty=1e6
        self.lb=lb_list
        self.ub=ub_list
        self.n_constraints=4
        self.constraints_bounds={}
        #a,b=find_flux_grups(label_model)
        if allow_unfeasible_solutions:
           for reaction_id in a:
              n=label_model.flux_solver_reaction_n_dict[reaction_id]
              reaction=label_model.constrained_model.reactions.get_by_id(reaction_id)
              ub=reaction.upper_bound
              lb=reaction.lower_bound
              if ub>=1000 and lb<=-1000:
                 print "omit",reaction_id,[lb,ub]
                 continue
              self.constraints_bounds[n]={"lb":lb,"ub":ub}
           self.n_constraints+=2*len(self.constraints_bounds)
        if "label_weight" in parameter_dict:
             self.label_weight=float(parameter_dict["label_weight" ])
        if "flux_weight" in parameter_dict:
             self.flux_weight=float(parameter_dict["flux_weight" ])
        if "target_flux_dict" in parameter_dict:
             self.target_flux_dict=parameter_dict["target_flux_dict" ]
        if "max_chi" in parameter_dict: 
             self.max_chi=parameter_dict["max_chi"]
        if "max_flux" in parameter_dict: 
             self.max_flux=parameter_dict["max_flux"]
        """if "verbose" in parameter_dict:
              self.verbose=parameter_dict["verbose"]"""
        if "flux_penalty_dict" in parameter_dict:
              self.flux_penalty_dict=parameter_dict["flux_penalty_dict"]
        if "label_unfeasible_penalty" in parameter_dict:
            self.label_unfeasible_penalty=parameter_dict["label_unfeasible_penalty"]
        if "flux_unfeasible_penalty" in parameter_dict:
            self.flux_unfeasible_penalty=parameter_dict["flux_unfeasible_penalty"]
    def get_nic(self):
        
        return self.n_constraints
        
    def get_nec(self): #TODO find if any reactions are equalities
        
        return 0 
         
    # Reimplement the virtual method that defines the objective function.
    def set_bounds(self,lb_list,ub_list):
        self.lb=lb_list
        self.ub=ub_list
        
    def gradient(self, x):
        return pygmo.estimate_gradient_h(lambda x: self.fitness(x), x)
        
    def get_bounds(self):
        return (self.lb,self.ub)
            
    def fitness(self, x):
        #print "running ob"
        # Compute the sphere function
        #try:
        #f,obj_dict = objfunc(label_model,x,simulate_infeasible=allow_unfeasible_solutions,label_weight=self.label_weight,mode="fsolve",target_flux_dict=self.target_flux_dict,max_chi=self.max_chi,max_flux=self.max_flux,flux_penalty_dict=self.flux_penalty_dict,verbose=False,flux_weight=self.flux_weight,label_unfeasible_penalty=self.label_unfeasible_penalty,flux_unfeasible_penalty=self.flux_unfeasible_penalty)
        f,obj_dict = objfunc(label_model,x,simulate_infeasible=self.allow_unfeasible_solutions,label_weight=self.label_weight,mode="fsolve",target_flux_dict=self.target_flux_dict,max_chi=1e12,max_flux=1e12,flux_penalty_dict=self.flux_penalty_dict,verbose=False,flux_weight=self.flux_weight,label_unfeasible_penalty=self.label_unfeasible_penalty,flux_unfeasible_penalty=self.flux_unfeasible_penalty)
        if f<self.fmin and self.verbose:
             self.fmin=f
             #print obj_dict
             flux_obj=obj_dict["flux_score"]
             label_obj=obj_dict["chi2_score"]
             fltarget_obj=obj_dict["fltarget_obj"]
             output="obj="+str(round(f,3))
             output+=" Chi2="+str(round(label_obj,3))
             if flux_obj>0:
                output+=" flux="+str(round(flux_obj,3))
             if fltarget_obj>0:
                irreversible_model_flag=self.target_flux_dict.get("irreversible_model")
                if irreversible_model_flag in (None,False):
                   output+=" "+self.target_flux_dict["reaction"]+"="+str(round(label_model.reversible_flux_dict[self.target_flux_dict["reaction"]],3))+" "+self.target_flux_dict["dir"]+" (lb:"+str(round(self.target_flux_dict["lb"],3))+" ub:"+str(round(self.target_flux_dict["ub"],3))+")"
                else:
                   output+=" "+self.target_flux_dict["reaction"]+"="+str(round(label_model.flux_dict[self.target_flux_dict["reaction"].replace("_forward","")],3))+" "+self.target_flux_dict["dir"]+" (lb:"+str(round(self.target_flux_dict["lb"],3))+" ub:"+str(round(self.target_flux_dict["ub"],3))+")"
             
             print output
        else:
           pass
           #print "not best: "+str(obj_dict) 
        # Note that we return a tuple with one element only. In PyGMO the objective functions
        # return tuples so that multi-objective optimization is also possible.
        #print f
        return_list=[f]
        constrains_list=[0]*self.n_constraints
        constrains_list[0]=obj_dict["chi2_score"]-self.max_chi
        constrains_list[1]=obj_dict["flux_score"]-self.max_flux
        counter=2
        if self.target_flux_dict not in [{},None]:
           irreversible_model_flag=self.target_flux_dict.get("irreversible_model")
           if irreversible_model_flag in (None,False):
              flux=label_model.reversible_flux_dict[self.target_flux_dict["reaction"]]
           else:
              flux=label_model.flux_dict[self.target_flux_dict["reaction"].replace("_forward","")]
           constrains_list[2]=flux-self.target_flux_dict["ub"]
           constrains_list[3]=self.target_flux_dict["lb"]-flux
        if self.allow_unfeasible_solutions:
           for n in self.constraints_bounds:
               lb=self.constraints_bounds[n]["lb"]
               ub=self.constraints_bounds[n]["ub"]
               constrains_list[counter]=label_model.flux_solution[n]-ub
               counter+=1
               constrains_list[counter]=lb-label_model.flux_solution[n]
               counter+=1
        
           
        #print return_list, max(constrains_list)       
        #return_list+constraints_bounds
        return return_list+constrains_list  
        #except:
        #print "Unkown error in objfun"
        #return (99999999999999999, )
  return   iso2flux_problem
        



"""

def optimize(label_model,pop_size = 25,n_gen = 500,n_islands=6,evolves_per_cycle=10,max_evolve_cycles=8,stop_criteria_relative=0.01,stop_criteria_absolute=-1e6,initial_archi_x=[],lb_list=[],ub_list=[],flux_penalty_dict=None,max_flux=None,label_problem_parameters={},min_model=None):

  if flux_penalty_dict==None or flux_penalty_dict=={}:
          flux_penalty_dict={}
          if "flux_penalty_dict" in  label_problem_parameters:
               flux_penalty_dict=label_problem_parameters["flux_penalty_dict"]
  if  max_flux==None:
           max_flux=1e6
           if "max_flux" in  label_problem_parameters:
               max_flux=label_problem_parameters["max_flux"]
  if min_model==None:
       min_model,minimal_flux,non_milp_reactions =create_minimal_fux_model(copy.deepcopy(label_model.constrained_model),fraction_of_optimum_objective=0.0,      fraction_of_flux_minimum=None,boundaries_precision=0.001,label_model=None,metabolite_list_file_name=None,flux_penalty_dict=flux_penalty_dict,maximum_flux=max_flux, mutually_exclusive_directionality_constraint=False,extra_constraints_dict={}) 
  else:
       min_model.reactions.get_by_id("total_flux").upper_bound=max_flux
  label_problem=label_model.iso2flux_problem(dim=len(label_model.variable_vector))
  with min_model:
    if lb_list==[] or ub_list==[]:
       
       lb_list,ub_list=get_variable_bounds(label_model,flux_penalty_dict=flux_penalty_dict,max_flux=max_flux,min_model=min_model)
    print   lb_list, ub_list,max_flux
    for n,x in enumerate(lb_list): #If upper and lower bound are identical program will get stuck creating problems
        if lb_list[n]-ub_list[n]<1e-9:
           ub_list[n]+=1e-6
    label_problem.set_bounds(lb_list,ub_list)
    label_problem.set_parameters(label_problem_parameters)
    archi_x=copy.deepcopy(initial_archi_x)
    algo_1 = algorithm.jde(gen = n_gen, xtol=1e-9, ftol=1e-9)
    
    total_flux_reaction=min_model.reactions.get_by_id("total_flux")
    total_flux_reaction.objective_coefficient=-1
    
    #return min_model,lb_list,ub_list
    #min_model=None
    sampled_variables=variable_sampling(label_model,lb_list,ub_list,maximum_flux=max_flux,flux_penalty_dict=flux_penalty_dict,n_pop=pop_size*n_islands,n_processes=1)
    for ncycle in range(0,max_evolve_cycles):
        print "starting new cycle"
        archi = archipelago(topology = topology.ring()) #Empty archipelago 
        print "created empty archipelago"
        counter=0
        for n in range(n_islands):
          isl = island(algo_1,label_problem,pop_size)
          print "created island"
          for n2 in range(pop_size):
               if archi_x==[]:
                  print "sampling pop "+str(n2)+" in island "+str(n)
                  print max_flux
                  #sample_variables=variable_sampling(label_model,lb_list,ub_list,flux_penalty_dict=flux_penalty_dict,max_flux=max_flux,flux_min_model=min_model)
                  
                  isl.set_x(n2,sampled_variables[counter])#variable_sampling(label_model,lb_list,ub_list,max_flux=max_flux,flux_penalty_dict=flux_penalty_dict,flux_min_model=min_model))
                  counter+=1
               else:
                  print "end of cycle"
                  isl.set_x(n2,archi_x[n][n2])
          archi.push_back(isl)
        print [isl.population.champion.f for isl in archi]
        if ncycle==0:
           previous_best_objective=min([isl.population.champion.f for isl in archi])[0]
        print [isl.population.champion.f for isl in archi]
        print evolves_per_cycle
        archi.evolve(evolves_per_cycle)
        archi.join()
        best_objective,optimal_variables,archi_x,champion_x_list,v_list=extract_results(archi)
        print [previous_best_objective,best_objective]
        if previous_best_objective-best_objective<max(previous_best_objective*stop_criteria_relative,1e-4) or best_objective<stop_criteria_absolute:
           print "convergence reached at cycle "+str(ncycle)
           break
        else:
           previous_best_objective=best_objective
    return best_objective,optimal_variables,archi_x

"""

def minimize_fluxes(label_model,iso2flux_problem,label_problem_parameters,max_chi=999999,flux_penalty_dict={} ,pop_size=20 ,n_gen=500 ,n_islands=6 ,max_cycles_without_improvement=8 ,max_evolve_cycles=20 ,stop_criteria_relative=0.005 ,max_iterations=10,  initial_flux_estimation=1000,log_name="log.txt",migrate="ring",max_flux_sampling=None):
    if max_flux_sampling==None:
       max_flux_sampling= initial_flux_estimation  
    min_model,minimal_flux,non_milp_reactions =create_minimal_fux_model(copy.deepcopy(label_model.constrained_model),fraction_of_optimum_objective=0.0,      fraction_of_flux_minimum=None,boundaries_precision=0.001,label_model=None,metabolite_list_file_name=None,flux_penalty_dict=flux_penalty_dict,maximum_flux=max_flux_sampling, mutually_exclusive_directionality_constraint=False,extra_constraints_dict={}) 
    best_flux=max_flux=initial_flux_estimation
    best_variables=None
    f=open(log_name, "a")
    f.write(time.strftime("%c")+"//"+"Starting flux minimization\nmax_chi="+str(max_chi)+"\n")
    f.close()
    label_problem_parameters["max_chi"]=min(label_problem_parameters["max_chi"],max_chi)
    for iteration in range(0,max_iterations):
        label_problem_parameters["max_flux"]=max_flux
        #label_problem_parameters={"label_weight":0.001,"target_flux_dict":None,"max_chi":max_chi,"max_flux":max_flux,"flux_penalty_dict":flux_penalty_dict,"verbose":True,"flux_weight":1,"flux_unfeasible_penalty":1e6,"label_unfeasible_penalty":1e3}
        flux_objective,flux_optimal_variables=optimize(label_model,iso2flux_problem,pop_size = pop_size,n_gen = n_gen,n_islands=n_islands,max_evolve_cycles=max_evolve_cycles,max_cycles_without_improvement=max_cycles_without_improvement, stop_criteria_relative=stop_criteria_relative,initial_archi_x=[],lb_list=[],ub_list=[],flux_penalty_dict=flux_penalty_dict, max_flux=max_flux,label_problem_parameters=label_problem_parameters,min_model=min_model,migrate=migrate)
        #a,objective_dict=objfunc(label_model,flux_optimal_variables,flux_penalty_dict=flux_penalty_dict,flux_weight=1)
        #flux_objective=objective_dict["flux_score"]
        f=open(log_name, "a")
        f.write(time.strftime("%c")+"//Flux Objective="+str(flux_objective)+" parameters are:\n"+str(flux_optimal_variables)+"\n")
        f.close() 
        
        if best_flux-flux_objective>stop_criteria_relative*best_flux:
           print best_flux,flux_objective,stop_criteria_relative
           best_flux=max_flux=flux_objective  
           best_variables=flux_optimal_variables
        else: 
           if best_flux>flux_objective:
              best_flux=max_flux=flux_objective  
              best_variables=flux_optimal_variables              
           print best_flux,flux_objective,stop_criteria_relative
           break
        if abs(minimal_flux-flux_objective)<1e-6:
           print "minimal flux reached, breaking"
           break
    f=open(log_name, "a")
    f.write(time.strftime("%c")+"//"+"Flux minimization done\n\n")
    f.close()
    return best_flux,best_variables




def update_flux_interval_dict(flux_interval_dict,parameters=[],max_chi=99999,max_flux=99999,label_model=None,irreversible_flux_dict=False):
    label_unfeasible_penalty=1e9
    flux_unfeasible_penalty=1e9
    if len(parameters)>0:
       obj, obj_dict=objfunc(label_model,parameters,max_chi=max_chi,max_flux=max_flux,label_unfeasible_penalty=1e9,flux_unfeasible_penalty=1e9)
          
       print obj
       if obj>=label_unfeasible_penalty:
          print "not updating fluxes"
          return
    if irreversible_flux_dict:
       flux_dict=label_model.flux_dict
       print flux_dict
       
    else:
       flux_dict=label_model.reversible_flux_dict
    if flux_interval_dict=={}:
       for flux in flux_dict:
           value=flux_dict[flux]
           flux_interval_dict[flux]={"maximum":value,"minimum":value}
           
    else:
        for flux in flux_interval_dict:
            flux_interval_dict[flux]["maximum"]=max(flux_interval_dict[flux]["maximum"],flux_dict[flux])
            flux_interval_dict[flux]["minimum"]=min(flux_interval_dict[flux]["minimum"],flux_dict[flux])



def flux_variation(label_model,iso2flux_problem,fluxes_to_evaluate,reference_variables,label_problem_parameters,max_chi=999999,max_flux=999999,flux_penalty_dict={} ,pop_size=20 ,n_gen=500 ,n_islands=6 ,max_cycles_without_improvement=10 ,stop_criteria_relative=0.005 ,max_iterations=10,log_name="log.txt",flux_interval_dict={},irreversible_flux_interval_dict={},max_flux_sampling=None,migrate="ring"):
    
    min_model,minimal_flux,non_milp_reactions =create_minimal_fux_model(copy.deepcopy(label_model.constrained_model),fraction_of_optimum_objective=0.0,      fraction_of_flux_minimum=None,boundaries_precision=0.001,label_model=None,metabolite_list_file_name=None,flux_penalty_dict=flux_penalty_dict,maximum_flux=max_flux, mutually_exclusive_directionality_constraint=False,extra_constraints_dict={}) 
    """best_flux=max_flux=initial_flux_estimation
    best_variables=None"""
    f=open(log_name, "w")
    f.write(time.strftime("%c")+"//"+"Starting Variation estimation\nmax_chi="+str(max_chi)+"\n")
    f.close()
    for parameters in reference_variables:
        update_flux_interval_dict(flux_interval_dict,parameters=parameters,max_chi=max_chi,max_flux=max_flux,label_model=label_model,irreversible_flux_dict=False)
        update_flux_interval_dict(irreversible_flux_interval_dict,parameters=parameters,max_chi=max_chi,max_flux=max_flux,label_model=label_model,irreversible_flux_dict=True)
    label_problem_parameters["max_chi"]=max_chi
    label_problem_parameters["max_flux"]=max_flux
    for reaction_id in fluxes_to_evaluate:
        print reaction_id
        if "_forward" in reaction_id or "_reverse" in reaction_id:
           irreversible_model_flag=True 
           base_reaction=reaction_id.replace("_reverse","").replace("_forward","")
           fva_reaction_2_test=flux_minimization_fva(min_model,solver=None,reaction_list=[base_reaction])
           original_minimum=fva_reaction_2_test[base_reaction]["minimum"]
           original_maximum=fva_reaction_2_test[base_reaction]["maximum"]
           fva_reaction_2_test[reaction_id]={}
           if  "_reverse" in reaction_id:
              fva_reaction_2_test[reaction_id]["minimum"]=max(original_minimum,0)
              fva_reaction_2_test[reaction_id]["maximum"]=max(original_maximum,0)+label_model.turnover_flux_dict[base_reaction]["ub"]  
           elif "_forward" in reaction_id:
              fva_reaction_2_test[reaction_id]["minimum"]=max(original_minimum,0)
              fva_reaction_2_test[reaction_id]["maximum"]=max(original_maximum,0)+label_model.turnover_flux_dict[base_reaction]["ub"]  
        else:
           irreversible_model_flag=False
           fva_reaction_2_test=flux_minimization_fva(min_model,solver=None,reaction_list=[reaction_id])
           print fva_reaction_2_test
        #fva_reaction_2_test=flux_variability_analysis(label_model.constrained_model,reaction_list=[reaction_2_test]) #TODO change by a flux minimization FVA
        reaction_2_test_ub=fva_reaction_2_test[reaction_id]["maximum"]
        reaction_2_test_lb=fva_reaction_2_test[reaction_id]["minimum"]
        task_list={"max":"maximum","min":"minimum"}
        print irreversible_flux_interval_dict
        for task in task_list:
          for iteration in range(0,max_iterations):
                target_flux_dict={"reaction":reaction_id,"dir":task,"lb":reaction_2_test_lb,"ub":reaction_2_test_ub,"factor":1,"irreversible_model":irreversible_model_flag}
                print flux_interval_dict
                print reaction_2_test_ub
                extra_constraint_dict={}
                if task=="max":
                 if not irreversible_model_flag:
                   if abs(flux_interval_dict[reaction_id]["maximum"]-reaction_2_test_ub)<1e-5:
                        break
                   maximum=flux_interval_dict[reaction_id]["maximum"]
                   extra_constraint_dict={reaction_id:{"lb":maximum}}
                   target_flux_dict["lb"]=maximum
                 else:
                   if abs(irreversible_flux_interval_dict[reaction_id.replace("_forward","")]["maximum"]-reaction_2_test_ub)<1e-5:
                        break
                   maximum=irreversible_flux_interval_dict[reaction_id.replace("_forward","")]["maximum"]
                   extra_constraint_dict={}
                   target_flux_dict["lb"]=maximum  
                else:
                  if not irreversible_model_flag:
                   if abs(flux_interval_dict[reaction_id]["minimum"]-reaction_2_test_lb)<1e-5:
                     break
                   minimum=flux_interval_dict[reaction_id]["minimum"]
                   extra_constraint_dict={reaction_id:{"ub":minimum}}
                   target_flux_dict["ub"]=minimum
                  else:
                   if abs(irreversible_flux_interval_dict[reaction_id.replace("_forward","")]["minimum"]-reaction_2_test_lb)<1e-5:
                     break
                   minimum=irreversible_flux_interval_dict[reaction_id.replace("_forward","")]["minimum"]
                   extra_constraint_dict={}
                   target_flux_dict["ub"]=minimum
                label_problem_parameters["target_flux_dict"]=target_flux_dict
                if not irreversible_model_flag:
                   if "max"==task:
                      ub=target_flux_dict["ub"]
                      best_objective=(ub-label_model.reversible_flux_dict[reaction_id]) #Minimize the distance between the upper bound and the flux value
                   elif "min"==task:
                      lb=target_flux_dict["lb"]
                      best_objective=(label_model.reversible_flux_dict[reaction_id]-lb)
                if irreversible_model_flag:
                   if "max"==task:
                      ub=target_flux_dict["ub"]
                      best_objective=(ub-label_model.flux_dict[reaction_id.replace("_forward","")]) #Minimize the distance between the upper bound and the flux value
                   elif "min"==task:
                      lb=target_flux_dict["lb"]
                      best_objective=(label_model.flux_dict[reaction_id.replace("_forward","")]-lb)                   
                   extra_constraint_dict={} #TODO use constraints with the irreversible fluxes                 
                """if task=="max":
                  if abs(flux_interval_dict[reaction_id]["maximum"]-reaction_2_test_ub)<1e-5:
                     break
                  if label_model.constrained_model.reactions.get_by_id(reaction_id).lower_bound>=0: #If the reaction is not reversible
                    maximum=flux_interval_dict[reaction_id]["maximum"]
                    min_model.reactions.get_by_id(reaction_id).lower_bound=max(reaction_2_test_lb,flux_interval_dict[reaction_id]["maximum"])
                  else: 
                    min_model.reactions.get_by_id(reaction_id).lower_bound=max(max(reaction_2_test_lb,flux_interval_dict[reaction_id]["maximum"]),0)
                    min_model.reactions.get_by_id(reaction_id+"_reverse").upper_bound=max(min(-1*reaction_2_test_lb,-1*flux_interval_dict[reaction_id]["maximum"]),0)
                elif task=="min":
                  if abs(flux_interval_dict[reaction_id]["minimum"]-reaction_2_test_lb)<1e-5:
                     break
                  if label_model.constrained_model.reactions.get_by_id(reaction_id).lower_bound>=0: #If the reaction is not reversible
                    min_model.reactions.get_by_id(reaction_id).upper_bound=min(reaction_2_test_ub,flux_interval_dict[reaction_id]["minimum"])
                  else: 
                    min_model.reactions.get_by_id(reaction_id).upper_bound=max(min(reaction_2_test_ub,flux_interval_dict[reaction_id]["minimum"]),0)
                    min_model.reactions.get_by_id(reaction_id+"_reverse").lower_bound=max(max(-1*reaction_2_test_ub,-1*flux_interval_dict[reaction_id]["minimum"]),0)"""                     
                f=open(log_name, "a")
                f.write("starting "+str(target_flux_dict)+"\n")
                objective,optimal_variables=optimize(label_model,iso2flux_problem,pop_size = pop_size,n_gen = n_gen,n_islands=n_islands,max_cycles_without_improvement=max_cycles_without_improvement, stop_criteria_relative=stop_criteria_relative,initial_archi_x=[],lb_list=[],ub_list=[],flux_penalty_dict=flux_penalty_dict, max_flux=max_flux,label_problem_parameters=label_problem_parameters,extra_constraint_dict=extra_constraint_dict,max_flux_sampling=max_flux_sampling,migrate=migrate)
                objfunc(label_model,optimal_variables)
                update_flux_interval_dict(flux_interval_dict,optimal_variables,label_problem_parameters["max_chi"],label_problem_parameters["max_flux"],label_model=label_model)
                update_flux_interval_dict(irreversible_flux_interval_dict,optimal_variables,label_problem_parameters["max_chi"],label_problem_parameters["max_flux"],label_model=label_model,irreversible_flux_dict=True)
                #return flux_interval_dict
                f=open(log_name, "a")
                try:
                  f.write(time.strftime("%c")+"//"+reaction_id+task+" "+str(flux_interval_dict[reaction_id][task_list[task]])+" parameters are:\n"+str(optimal_variables)+"\n")
                except:
                  print irreversible_flux_interval_dict
                  f.write(time.strftime("%c")+"//"+reaction_id+task+" "+str(irreversible_flux_interval_dict[reaction_id.replace("_forward","")][task_list[task]])+" parameters are:\n"+str(optimal_variables)+"\n")
                f.write(str(flux_interval_dict)+"\n")
                f.close()
                if best_objective-objective>stop_criteria_relative*best_objective:
                   best_objective=objective
                   best_variables=optimal_variables
                else: 
                   if best_objective>objective:
                      best_objective=objective
                      best_variables=optimal_variables
                   break
                if abs(objective)<1e-6:
                   print "Bound reached, breaking"
                   break
    
    return flux_interval_dict,irreversible_flux_interval_dict


"""def update_flux_interval_dict(flux_interval_dict,parameters=[],max_chi=99999,max_flux=99999,label_model=None):
    label_unfeasible_penalty=1e9
    flux_unfeasible_penalty=1e9
    if len(parameters)>0:
       obj, obj_dict=objfunc(label_model,parameters,max_chi=max_chi,max_flux=max_flux,label_unfeasible_penalty=1e9,flux_unfeasible_penalty=1e9)
       print obj
       if obj>=label_unfeasible_penalty:
          print "not updating fluxes"
          return
    if flux_interval_dict=={}:
       for flux in label_model.reversible_flux_dict:
           value=label_model.reversible_flux_dict[flux]
           flux_interval_dict[flux]={"maximum":value,"minimum":value}
    else:
        for flux in flux_interval_dict:
            flux_interval_dict[flux]["maximum"]=max(flux_interval_dict[flux]["maximum"],label_model.reversible_flux_dict[flux])
            flux_interval_dict[flux]["minimum"]=min(flux_interval_dict[flux]["minimum"],label_model.reversible_flux_dict[flux])

"""

"""def flux_variation(label_model,iso2flux_problem,fluxes_to_evaluate,reference_variables,label_problem_parameters,max_chi=999999,max_flux=999999,flux_penalty_dict={} ,pop_size=20 ,n_gen=500 ,n_islands=6 ,max_cycles_without_improvement=10 ,stop_criteria_relative=0.005 ,max_iterations=10,log_name="log.txt"):
    
    min_model,minimal_flux,non_milp_reactions =create_minimal_fux_model(copy.deepcopy(label_model.constrained_model),fraction_of_optimum_objective=0.0,      fraction_of_flux_minimum=None,boundaries_precision=0.001,label_model=None,metabolite_list_file_name=None,flux_penalty_dict=flux_penalty_dict,maximum_flux=max_flux, mutually_exclusive_directionality_constraint=False,extra_constraints_dict={}) 
    f=open(log_name, "w")
    f.write(time.strftime("%c")+"//"+"Starting Variation estimation\nmax_chi="+str(max_chi)+"\n")
    f.close()
    flux_interval_dict={}
    for parameters in reference_variables:
        update_flux_interval_dict(flux_interval_dict,parameters=parameters,max_chi=max_chi,max_flux=max_flux,label_model=label_model)
    label_problem_parameters["max_chi"]=max_chi
    label_problem_parameters["max_flux"]=max_flux
    for reaction_id in fluxes_to_evaluate:
        fva_reaction_2_test=flux_minimization_fva(min_model,solver=None,reaction_list=[reaction_id])
        #fva_reaction_2_test=flux_variability_analysis(label_model.constrained_model,reaction_list=[reaction_2_test]) #TODO change by a flux minimization FVA
        reaction_2_test_ub=fva_reaction_2_test[reaction_id]["maximum"]
        reaction_2_test_lb=fva_reaction_2_test[reaction_id]["minimum"]
        task_list={"max":"maximum","min":"minimum"}
        for task in task_list:
          target_flux_dict={"reaction":reaction_id,"dir":task,"lb":reaction_2_test_lb,"ub":reaction_2_test_ub,"factor":1}
          if "max"==task:
             ub=target_flux_dict["ub"]
             best_objective=(ub-label_model.reversible_flux_dict[reaction_id]) #Minimize the distance between the upper bound and the flux value
          elif "min"==task:
             lb=target_flux_dict["lb"]
             best_objective=(label_model.reversible_flux_dict[reaction_id]-lb)  
          label_problem_parameters["target_flux_dict"]=target_flux_dict        
          for iteration in range(0,max_iterations):
              
                print flux_interval_dict
                print reaction_2_test_ub
                extra_constraint_dict={}
                if task=="max":
                   if abs(flux_interval_dict[reaction_id]["maximum"]-reaction_2_test_ub)<1e-5:
                     break
                   maximum=flux_interval_dict[reaction_id]["maximum"]
                   extra_constraint_dict={reaction_id:{"lb":maximum}}
                else:
                   if abs(flux_interval_dict[reaction_id]["minimum"]-reaction_2_test_lb)<1e-5:
                     break
                   minimum=flux_interval_dict[reaction_id]["minimum"]
                   extra_constraint_dict={reaction_id:{"ub":minimum}}
                
                
                #print "AAAAAAAAAAAAAAAA",reaction_id,min_model.reactions.get_by_id(reaction_id).upper_bound,min_model.reactions.get_by_id(reaction_id).lower_bound
                objective,optimal_variables=optimize(label_model,iso2flux_problem,pop_size = pop_size,n_gen = n_gen,n_islands=n_islands,max_cycles_without_improvement=max_cycles_without_improvement, stop_criteria_relative=stop_criteria_relative,initial_archi_x=[],lb_list=[],ub_list=[],flux_penalty_dict=flux_penalty_dict, max_flux=max_flux,label_problem_parameters=label_problem_parameters,extra_constraint_dict=extra_constraint_dict)
                objfunc(label_model,optimal_variables)
                update_flux_interval_dict(flux_interval_dict,optimal_variables,label_problem_parameters["max_chi"],label_problem_parameters["max_flux"],label_model=label_model)
                #return flux_interval_dict
                f=open(log_name, "a")
                f.write(time.strftime("%c")+"//"+reaction_id+task+" "+str(flux_interval_dict[reaction_id][task_list[task]])+" parameters are:\n"+str(optimal_variables)+"\n"+str(flux_interval_dict))
                f.close()
                if best_objective-objective>stop_criteria_relative*best_objective:
                   best_objective=objective
                   best_variables=optimal_variables
                else: 
                   if best_objective>objective:
                      best_objective=objective
                      best_variables=optimal_variables
                   break
                if abs(objective)<1e-6:
                   print "Bound reached, breaking"
                   break
    
    return flux_interval_dict

"""


"""
path_to_problem = os.path.dirname(fitting.__file__)+"/pygmo_problem.py"
f=open(path_to_problem,"r")
problem_script=f.read()
f.close()
exec(problem_script)
label_model.iso2flux_problem=iso2flux_problem

best_objective,optimal_variables,archi_x=optimize(label_model,pop_size = 30,n_gen = 200,n_islands=2,evolves_per_cycle=4,max_evolve_cycles=20,stop_criteria_relative=0.01,initial_archi_x=[],lb_list=[],ub_list=[])

label_problem_parameters={"label_weight":0.0001,"flux_weight":0.0001}
flux_penalty_dict=build_flux_penalty_dict(label_model)
flux_variation(label_model,reference_flux_group_dict,optimal_variables,label_problem_parameters,max_chi=0.1,max_flux=10,flux_penalty_dict=flux_penalty_dict ,pop_size=20 ,n_gen=500 ,n_islands=2 ,evolves_per_cycle=8 ,max_evolve_cycles=20 ,stop_criteria_relative=0.005 ,max_iterations=10,log_name="confidence.txt")
"""

