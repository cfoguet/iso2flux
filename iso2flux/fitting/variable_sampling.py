from ..flux_functions.minimal_flux import create_minimal_fux_model,flux_minimization_fva
import random
import numpy as np
import copy
from ..misc.round_functions import round_up, round_down 
from cobra.flux_analysis import sample
from ..flux_functions.flux_variability_analysis import flux_variability_analysis
from cobra import Reaction
import json

def variable_sampling(label_model,lb_list,ub_list,maximum_flux=1e6,flux_penalty_dict={},n_pop=40,n_processes=6,extra_constraints_dict={},debug=False):
  with label_model.constrained_model as model:
     for reaction_id in extra_constraints_dict:
        reaction=model.reactions.get_by_id(reaction_id)
        if "lb" in extra_constraints_dict[reaction_id]:
            reaction.lower_bound=round_down(extra_constraints_dict[reaction_id]["lb"],6)
        if "ub" in extra_constraints_dict[reaction_id]:
            reaction.upper_bound=round_up(extra_constraints_dict[reaction_id]["ub"],6)
     fva=flux_variability_analysis(model,fraction_of_optimum=0)
     reactions2add=[]
     turnover_reaction_list=[]
     for reaction_id in fva:
         if maximum_flux==None:
            continue
         #print reaction_id
         reaction=model.reactions.get_by_id(reaction_id)
         if reaction.lower_bound<0:
            if reaction_id not in label_model.turnover_flux_dict:
               reaction.upper_bound=min(round_up(fva[reaction_id]["maximum"],8),reaction.upper_bound)
               reaction.lower_bound=max(round_down(fva[reaction_id]["minimum"],8),reaction.lower_bound)
            elif reaction_id in label_model.turnover_flux_dict and reaction.upper_bound>0:
               reaction.upper_bound=min(round_up(fva[reaction_id]["maximum"],8)+label_model.turnover_flux_dict[reaction_id]["ub"],reaction.upper_bound)
               reaction.lower_bound=max(round_down(fva[reaction_id]["minimum"],8)-label_model.turnover_flux_dict[reaction_id]["ub"],reaction.lower_bound)
               if reaction.lower_bound>=0:
                  #create a dummy reverse reaction
                  reverse_reaction = Reaction(reaction_id + "_reverse")
                  reverse_reaction.name="turnover"
                  reactions2add.append(reverse_reaction)
                  reverse_reaction.upper_bound =label_model.turnover_flux_dict[reaction_id]["ub"]*2
                  turnover_reaction_list.append(reaction_id + "_reverse")
            elif reaction_id in label_model.turnover_flux_dict and reaction.upper_bound<=0:
                 #create a reverse reaction
                 reverse_reaction = Reaction(reaction_id + "_reverse")
                 #reverse_reaction.name="turnover"
                 reaction_dict ={}
                 for k in reaction._metabolites:
                     reaction_dict[k]=reaction._metabolites[k] *-1
                 reverse_reaction.add_metabolites(reaction_dict)
                 reverse_reaction.lower_bound = max(0, -reaction.upper_bound)
                 reverse_reaction.upper_bound = -reaction.lower_bound
                 reactions2add.append(reverse_reaction)
                 #Create a dummy forward reaction
                 reaction.lower_bound =0
                 reaction.name="turnover"
                 reaction.upper_bound = label_model.turnover_flux_dict[reaction_id]["ub"]*2 #to represent that turnover has double effect ont fluxes
                 reaction.add_metabolites(reaction_dict) #Adding this removes all metabolites creating an empty reaction
                 turnover_reaction_list.append(reaction_id)
          
         elif reaction.upper_bound>0 and reaction_id in label_model.turnover_flux_dict:
              #create a dummy reverse reaction
              reverse_reaction = Reaction(reaction_id + "_reverse")
              reverse_reaction.name="turnover"
              reactions2add.append(reverse_reaction)
              reverse_reaction.upper_bound =label_model.turnover_flux_dict[reaction_id]["ub"]*2
              turnover_reaction_list.append(reaction_id + "_reverse")
         elif reaction.upper_bound==0 and reaction.lower_bound==0 and reaction_id in label_model.turnover_flux_dict :
              reverse_reaction = Reaction(reaction_id + "_reverse")
              reverse_reaction.name="turnover"
              reactions2add.append(reverse_reaction)
              reverse_reaction.upper_bound =label_model.turnover_flux_dict[reaction_id]["ub"]*2
              turnover_reaction_list.append(reaction_id + "_reverse")
     if maximum_flux==None:
            min_model=model.copy()
            fluxes_dict=sampling(min_model,n=n_pop*n_processes,processes=n_processes,objective=None)
     else:        
        model.add_reactions(reactions2add)
        min_model,minimal_flux,non_milp_reactions =create_minimal_fux_model(model,fraction_of_optimum_objective=0, boundaries_precision=0.001,label_model=None,metabolite_list_file_name=None,flux_penalty_dict=flux_penalty_dict,maximum_flux=maximum_flux,extra_constraints_dict={})
        fluxes_dict=sampling(min_model,n=n_pop*n_processes,processes=n_processes,objective="total_flux")
  if debug:
    f=open("fluxes_dict.json","w")
    json.dump(fluxes_dict,f)
    f.close()
     
  variable_sets=[]
  #print label_model.turnover_flux_dict  
  for flux_dict in fluxes_dict:
     variables=[]
     for n,variable in enumerate(label_model.variable_list_dict):
        lb=lb_list[n]
        ub=ub_list[n]
        if variable["type"]=="NULL":
           value=0
        elif variable["type"]=="flux":
             reaction=label_model.flux_solver_independent_flux_dict[variable["n"]].keys()[0]
             forward_flux=0
             reverse_flux=0
             reverse_id=reaction+"_reverse"
             if reverse_id in flux_dict and reverse_id not in turnover_reaction_list:
                reverse_flux=flux_dict[reverse_id]
             if reaction in flux_dict and reaction not in turnover_reaction_list:
                forward_flux=flux_dict[reaction] 
                
             coef=label_model.flux_solver_independent_flux_dict[variable["n"]][reaction]
             value=(forward_flux-reverse_flux)/coef
             #variables.append((forward_flux-reverse_flux)/coef)
        elif variable["type"]=="turnover":
             if maximum_flux==None:
                value=random.uniform(lb_list[n],ub_list[n])
             else:
               reaction=variable["reaction"]
               reverse_id=reaction+"_reverse"
               if reverse_id in turnover_reaction_list:
                  value=(flux_dict[reverse_id]/2.0)
               elif reaction in  turnover_reaction_list:
                  value=(flux_dict[reaction]/2.0)
               else:
                if reverse_id not in flux_dict:
                   flux_dict[reverse_id]=0.0
                value=min(abs(flux_dict[reaction]),abs(flux_dict[reverse_id]))
        #print variable, value,lb,ub
        value=max(min(value,ub),lb)
        variables.append(value)
     variable_sets.append(variables)
  return variable_sets       


def sampling(model,n=100,processes=6,objective=None):
    reaction_ids=[x.id for x in model.reactions]
    #print model.reactions.get_by_id(objective).lower_bound
    result_matrix = sample(model, n,processes=processes).as_matrix()
    flux_dict_list=[]
    for row in result_matrix:
        flux_dict={}
        for n_flux,flux in enumerate(row):
            flux_dict[reaction_ids[n_flux]]=flux
        flux_dict_list.append(flux_dict)
        """if objective!=None:
           print flux_dict[objective]"""
    return flux_dict_list

"""
def variable_sampling(label_model,lb_list,ub_list,max_flux=1e6, flux_penalty_dict={},extra_constraints_dict={},flux_min_model=None):
    random.seed(a=None) #not sure if necessary
    variables=np.zeros(len(label_model.variable_list_dict))
    #order=np.random.choice(range(0,len(variables)), size=len(variables), replace=False, p=None)
    #order=[62, 75, 82, 15, 51, 40, 74, 43, 25,  4, 79,  8, 73, 80, 42, 24, 53,57, 67, 17, 14, 41, 52, 12, 33, 13, 35, 72, 32, 34, 31, 70, 29, 78,65, 11, 64, 71,  6,  3, 69, 76, 60, 37, 66, 10, 68, 26,  1, 48, 28,       30,  7, 22, 27, 50, 63, 59,  0, 21, 39, 36,  9, 19, 38, 47, 23, 45,       55, 20, 49,  5, 54, 46, 77, 44, 16, 18, 56, 58,  2, 61, 81]
    
    flux_variables=[]
    for n in range(0,len(variables)):
        vtype=label_model.variable_list_dict[n]["type"]
        if vtype=="flux":
           flux_variables.append(n)
        elif vtype=="turnover":
           variables[n]=random.uniform(lb_list[n],ub_list[n])
    
    copy_model=copy.deepcopy(label_model.constrained_model)
    if flux_min_model==None:
       print "creating model"
       min_model,minimal_flux,non_milp_reactions =create_minimal_fux_model(copy_model,fraction_of_optimum_objective=0.0,      fraction_of_flux_minimum=None,boundaries_precision=0.001,label_model=None,metabolite_list_file_name=None,flux_penalty_dict=flux_penalty_dict,maximum_flux=max_flux, mutually_exclusive_directionality_constraint=True,extra_constraints_dict=extra_constraints_dict)
       total_flux_reaction=min_model.reactions.get_by_id("total_flux")
       #total_flux_reaction.objective_coefficient=-1  
    else:
        min_model=flux_min_model
    #print max_flux
    #print min_model.reactions[-1].upper_bound
    random_order=np.random.choice(range(0,len(variables)), size=len(variables), replace=False, p=None)
    total_flux_reaction=min_model.reactions.get_by_id("total_flux")
    #total_flux_reaction.objective_coefficient=-1  
    for n in random_order: 
        vtype=label_model.variable_list_dict[n]["type"]
        if vtype=="flux": 
           #print "flux"
           nflux=label_model.variable_list_dict[n]["n"]
           reference_reaction_id=label_model.flux_solver_independent_flux_dict[nflux].keys()[0]
           coef=label_model.flux_solver_independent_flux_dict[nflux][reference_reaction_id]
           fva=flux_minimization_fva(min_model,solver=None,reaction_list=[reference_reaction_id],lp_tolerance_feasibility=1e-8)       
           #fva=flux_variability_analysis(copy_model,reaction_list=[reference_reaction_id], fraction_of_optimum=0,tolerance_feasibility=label_model.lp_tolerance_feasibility)
           lb=min(fva[reference_reaction_id]["minimum"]/coef,fva[reference_reaction_id]["maximum"]/coef)  
           ub=max(fva[reference_reaction_id]["minimum"]/coef,fva[reference_reaction_id]["maximum"]/coef)
           value=random.uniform(max(lb_list[n],lb),min(ub_list[n],ub))
           #print  lb,ub,value, coef
           if abs(lb-ub)>5e-3:
              reference_reaction=min_model.reactions.get_by_id(reference_reaction_id) 
              if "reflection" in reference_reaction.notes:
                  if value*coef<0:
                     reference_reaction.lower_bound=0
                     reference_reaction.upper_bound=0
                     reflection_reaction=min_model.reactions.get_by_id(reference_reaction.notes["reflection"])
                     reflection_reaction.lower_bound=round_down(-value*coef,3)
                     reflection_reaction.upper_bound=round_up(-value*coef,3)
                  else:
                     reference_reaction.lower_bound=round_down(value*coef,3)
                     reference_reaction.upper_bound=round_up(value*coef,3)
                     reflection_reaction=min_model.reactions.get_by_id(reference_reaction.notes["reflection"])
                     reflection_reaction.lower_bound=0
                     reflection_reaction.upper_bound=0
              else:
                     reference_reaction.lower_bound=round_down(value*coef,3)
                     reference_reaction.upper_bound=round_up(value*coef,3)
              #reference_reaction=copy_model.reactions.get_by_id(reference_reaction_id)
              #reference_reaction.lower_bound=round_down(value*coef,6)
              #reference_reaction.upper_bound=round_up(value*coef,6)
        elif  vtype=="turnover":
              #variables[n]=0; continue
              #total_flux_reaction=min_model.reactions.get_by_id("total_flux")
              total_flux_reaction.objective_coefficient=-1   
              #current_minimal_flux=min_model.optimize(tolerance_feasibility=label_model.lp_tolerance_feasibility).x_dict["total_flux"]
              current_minimal_flux=min_model.optimize(tolerance_feasibility=label_model.lp_tolerance_feasibility).x_dict["total_flux"]
              #print current_minimal_flux
              total_flux_reaction.objective_coefficient=0
              #current_flux_space=max_flux-current_minimal_flux #How much room do we have to increase flux
              reaction=label_model.variable_list_dict[n]["reaction"]
              base_penalty=0
              if reaction in flux_penalty_dict:
                 base_penalty+=flux_penalty_dict[reaction]
              if reaction+"_reverse" in flux_penalty_dict: #Assume that the reverse reactions is the name of forward reaction+_reverse
                 base_penalty+=flux_penalty_dict[reaction+"_reverse"]
              elif reaction in flux_penalty_dict: 
                 base_penalty+=flux_penalty_dict[reaction]
              original_ub=ub_list[n]
              original_lb=lb_list[n]
              if base_penalty==0: 
                 #print reaction
                 value=random.uniform(original_lb,original_ub)
                 variables[n]=max(min(value,ub_list[n]),lb_list[n])
              else:
                 #print(total_flux_reaction.upper_bound,current_minimal_flux,base_penalty)
                 current_upper_bound=(total_flux_reaction.upper_bound-current_minimal_flux)/float(base_penalty)
                 value=random.uniform(original_lb,min(original_ub,current_upper_bound)) 
                 variables[n]=max(min(value,ub_list[n]),lb_list[n])
                 total_flux_reaction.upper_bound=max(float(round_down(total_flux_reaction.upper_bound-base_penalty*variables[n],3)),current_minimal_flux) #Decrease the maximum flux to account for these new reaction
                 #print ("abxchdf",reaction,base_penalty,[original_ub,current_upper_bound],variables[n],current_upper_bound,current_minimal_flux,total_flux_reaction.upper_bound)
              #print current_minimal_flux,min_model.optimize(tolerance_feasibility=label_model.lp_tolerance_feasibility).x_dict["total_flux"] 
    #Once the model has been constrained solver and assugb the solution to the flux variables 
    min_model.reactions.get_by_id("total_flux").objective_coefficient=-1   
    min_model_solution=min_model.optimize(tolerance_feasibility=label_model.lp_tolerance_feasibility) 
   
    #print min_model.optimize(tolerance_feasibility=label_model.lp_tolerance_feasibility).x_dict["total_flux"] 
    #f=open("esborramt.txt","w")  
    #print len(flux_variables)
    model_reaction_n_dict={}
    for n,reaction in enumerate(min_model.reactions):
        model_reaction_n_dict[reaction.id]=n
    for n in flux_variables:
           nflux=label_model.variable_list_dict[n]["n"]
           reference_reaction_id=label_model.flux_solver_independent_flux_dict[nflux].keys()[0]
           n_reaction=model_reaction_n_dict[reference_reaction_id]
           reference_reaction=min_model.reactions.get_by_id(reference_reaction_id)
           coef=label_model.flux_solver_independent_flux_dict[nflux][reference_reaction_id]
           if "reflection" in reference_reaction.notes:
               reflection_id=reference_reaction.notes["reflection"]
               n_reflection=model_reaction_n_dict[reflection_id]
               value_reflection=min_model_solution.x[n_reflection]/coef
               value_forward=min_model_solution.x[n_reaction]/coef
               if abs(value_reflection)>abs(value_forward):
                  value=-value_reflection
               else:
                  value=value_forward
               #f.write(reference_reaction_id+str([n,coef,[value_forward,value_reflection],value,max(min(value,ub_list[n]),lb_list[n])])+"\n")
               #print  value_forward,value_reflection,value,max(min(value,ub_list[n]),lb_list[n])
           else:
               value=min_model_solution.x[n_reaction]/coef
               #f.write(reference_reaction_id+str([n,coef,value,max(min(value,ub_list[n]),lb_list[n])])+"\n")
           variables[n]=max(min(value,ub_list[n]),lb_list[n])
           
    #f.close()
    return variables

"""

"""
def variable_sampling(label_model,lb_list,ub_list,max_flux=1e6, flux_penalty_dict={},extra_constraints_dict={}):
    random.seed(a=None) #not sure if necessary
    variables=np.zeros(len(label_model.variable_list_dict))
    #order=np.random.choice(range(0,len(variables)), size=len(variables), replace=False, p=None)
    #order=[62, 75, 82, 15, 51, 40, 74, 43, 25,  4, 79,  8, 73, 80, 42, 24, 53,57, 67, 17, 14, 41, 52, 12, 33, 13, 35, 72, 32, 34, 31, 70, 29, 78,65, 11, 64, 71,  6,  3, 69, 76, 60, 37, 66, 10, 68, 26,  1, 48, 28,       30,  7, 22, 27, 50, 63, 59,  0, 21, 39, 36,  9, 19, 38, 47, 23, 45,       55, 20, 49,  5, 54, 46, 77, 44, 16, 18, 56, 58,  2, 61, 81]
    
    flux_variables=[]
    for n in range(0,len(variables)):
        vtype=label_model.variable_list_dict[n]["type"]
        if vtype=="flux":
           flux_variables.append(n)
        elif vtype=="turnover":
           variables[n]=random.uniform(lb_list[n],ub_list[n])
    
    copy_model=copy.deepcopy(label_model.constrained_model)
    min_model,minimal_flux,non_milp_reactions =create_minimal_fux_model(copy_model,fraction_of_optimum_objective=0.0,      fraction_of_flux_minimum=None,boundaries_precision=0.001,label_model=None,metabolite_list_file_name=None,flux_penalty_dict=flux_penalty_dict,maximum_flux=max_flux, mutually_exclusive_directionality_constraint=False,extra_constraints_dict=extra_constraints_dict)
    flux_variables_order=np.random.choice(flux_variables, size=len(flux_variables), replace=False, p=None)
    for n in flux_variables_order: #First randomly constraint the variables of the constraint based model
           nflux=label_model.variable_list_dict[n]["n"]
           reference_reaction_id=label_model.flux_solver_independent_flux_dict[nflux].keys()[0]
           coef=label_model.flux_solver_independent_flux_dict[nflux][reference_reaction_id]
           fva=flux_minimization_fva(min_model,solver=None,reaction_list=[reference_reaction_id],lp_tolerance_feasibility=1e-8)       
           #fva=flux_variability_analysis(copy_model,reaction_list=[reference_reaction_id], fraction_of_optimum=0,tolerance_feasibility=label_model.lp_tolerance_feasibility)
           lb=min(fva[reference_reaction_id]["minimum"]/coef,fva[reference_reaction_id]["maximum"]/coef)  
           ub=max(fva[reference_reaction_id]["minimum"]/coef,fva[reference_reaction_id]["maximum"]/coef)
           value=random.uniform(max(lb_list[n],lb),min(ub_list[n],ub))
           print fva
           print  lb,ub,value, coef
           if abs(lb-ub)>1e-3:
              print "a"
              reference_reaction=min_model.reactions.get_by_id(reference_reaction_id) 
              if "reflection" in reference_reaction.notes:
                  if value<0:
                     reference_reaction.lower_bound=0
                     reference_reaction.upper_bound=0
                     reflection_reaction=min_model.reactions.get_by_id(reference_reaction.notes["reflection"])
                     reflection_reaction.lower_bound=round_down(-value*coef,4)
                     reflection_reaction.upper_bound=round_up(-value*coef,4)
                  else:
                     reference_reaction.lower_bound=round_down(value*coef,4)
                     reference_reaction.upper_bound=round_up(value*coef,4)
                     reflection_reaction=min_model.reactions.get_by_id(reference_reaction.notes["reflection"])
                     reflection_reaction.lower_bound=0
                     reflection_reaction.upper_bound=0
              else:
                     reference_reaction.lower_bound=round_down(value*coef,4)
                     reference_reaction.upper_bound=round_up(value*coef,4)
              #reference_reaction=copy_model.reactions.get_by_id(reference_reaction_id)
              #reference_reaction.lower_bound=round_down(value*coef,6)
              #reference_reaction.upper_bound=round_up(value*coef,6)
    #Once the model has been constrained solver and assugb the solution to the flux variables 
    min_model.reactions.get_by_id("total_flux").objective_coefficient=-1   
    min_model.optimize(tolerance_feasibility=label_model.lp_tolerance_feasibility)   
    #f=open("esborramt.txt","w")   
    for n in flux_variables:
           nflux=label_model.variable_list_dict[n]["n"]
           reference_reaction_id=label_model.flux_solver_independent_flux_dict[nflux].keys()[0]
           reference_reaction=min_model.reactions.get_by_id(reference_reaction_id)
           coef=label_model.flux_solver_independent_flux_dict[nflux][reference_reaction_id]
           if "reflection" in reference_reaction.notes:
               reflection_id=reference_reaction.notes["reflection"]
               value_reflection=min_model.solution.x_dict[reflection_id]/coef
               value_forward=min_model.solution.x_dict[reference_reaction_id]/coef
               if abs(value_reflection)>abs(value_forward):
                  value=-value_reflection
               else:
                  value=value_forward
               #f.write(reference_reaction_id+str([n,coef,[value_forward,value_reflection],value,max(min(value,ub_list[n]),lb_list[n])])+"\n")
               #print  value_forward,value_reflection,value,max(min(value,ub_list[n]),lb_list[n])
           else:
               value=min_model.solution.x_dict[reference_reaction_id]/coef
               #f.write(reference_reaction_id+str([n,coef,value,max(min(value,ub_list[n]),lb_list[n])])+"\n")
           variables[n]=max(min(value,ub_list[n]),lb_list[n])
           
    #f.close()
    return variables  

flux_penalty_dict_normoxia


a,min_model=variable_sampling(label_model,lb_list,ub_list,max_flux=800, flux_penalty_dict=flux_penalty_dict_normoxia)
objfunc(label_model,a,flux_penalty_dict=flux_penalty_dict_normoxia,label_weight=0)
label_model.turnover_flux_dict

total=0
for x in label_model.turnover_flux_dict:
    total+=label_model.turnover_flux_dict[x]["v"]*2

total


for flux in flux_penalty_dict:
     if flux in label_model.flux_dict:
                   flux_score+=abs(label_model.flux_dict[flux])*flux_penalty_dict[flux]
     else: 
         print flux


flux_score=0
for flux in flux_penalty_dict:
     if flux in label_model.flux_dict:
         flux_score+=abs(label_model.flux_dict[flux])*flux_penalty_dict[flux]
     else: 
         print flux

flux_score

for x in label_model.flux_dict:
    if abs(label_model.flux_dict[x]-min_model.solution.x_dict[x])>0.01:
       if not abs(label_model.flux_dict[x]-min_model.solution.x_dict[x+"_reverse"])>0.01:
          print x,label_model.flux_dict[x],min_model.solution.x_dict[x]



flux_score=0
min_model.optimize()
for flux in flux_penalty_dict:
     if "reverse" in flux:
         continue
     if flux in min_model.solution.x_dict:
         flux_score+=abs(min_model.solution.x_dict[flux])*flux_penalty_dict[flux]
         if flux+"_reverse" in min_model.solution.x_dict:
            flux_score+=abs(min_model.solution.x_dict[flux+"_reverse"])*flux_penalty_dict[flux]
     else: 
         print flux

flux_score

flux_score=0
for flux in label_model.flux_dict:
                if flux in flux_penalty_dict:
                   flux_score+=abs(label_model.flux_dict[flux])*flux_penalty_dict[flux]

for reaction in model.reactions:
        if "IRRMILP_" in reaction.id or "RATIO_" in reaction.id or "RGROUP_" in reaction.id or "TMS_" in reaction.id:
            continue
        if reaction.id in flux_penalty_dict:
           reporter_coef=flux_penalty_dict[reaction.id]
        elif "reflection" in reaction.notes: #If a reversible reaction has no penalty assigned test if the reverse has a penalty
           reflection=reaction.notes["reflection"]
           if reflection in flux_penalty_dict:
              reporter_coef=flux_penalty_dict[reflection]
           else:
              reporter_coef=1
              print "Warning: reaction "+reaction.id+" assigned default weight of 1"
        else:
              reporter_coef=1
              print "Warning: reaction "+reaction.id+" assigned default weight of 1"

"""
