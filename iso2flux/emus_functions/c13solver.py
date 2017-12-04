import numpy as np
from scipy.optimize import broyden1
from scipy.optimize import broyden2
from scipy.optimize import fsolve
from scipy.integrate import odeint
import json 
import cobra
import time
import copy
from ..flux_functions.apply_ratios import apply_ratios
from ..flux_functions.get_fluxes import get_fluxes

#from ..iso2flux.flux_functions.check_bounds import check_bounds

def check_bounds(solution_dict,model,label_model_label_reactions={},penalty_factor=100000, verbose=False,label_reaction_tolerance=1e-8,non_label_reaction_tolerance=1e-6):
    penalty=0
    string=""
    for x in model.reactions:
        if x.id in label_model_label_reactions:
           tolerance=label_reaction_tolerance
        else:
           tolerance=non_label_reaction_tolerance
        reaction_id=x.id
        if reaction_id not in solution_dict:
           solution_dict[reaction_id]=0
        reaction_value=solution_dict[reaction_id]
        if reaction_value>(x.upper_bound+tolerance): 
           if verbose:
              string+="upper boun violated for "+x.id+str([reaction_value,x.upper_bound,x.lower_bound])+"\n"
              print("upper boun violated for "+x.id+str([reaction_value,x.upper_bound,x.lower_bound]))
           penalty+=penalty_factor*(reaction_value-x.upper_bound) 
        if reaction_value<(x.lower_bound-tolerance):
           if verbose:
              string+="lowerd bound violated for "+x.id+str([reaction_value,x.upper_bound,x.lower_bound])+"\n"
              print("lowerd bound violated for "+x.id+str([reaction_value,x.upper_bound,x.lower_bound]))
           penalty+=penalty_factor*(x.lower_bound-reaction_value) 
    #f=open("check_bounds.txt","w")
    #f.write(string)
    #f.close()
    return penalty

def flux_solver(nullmnp,free_fluxes_vector,n_reaction_dict=None,label_model=None):
   if label_model!=None:
      flux_solution=label_model.compute_fluxes(nullmnp,free_fluxes_vector)#np.dot(nullmnp,free_fluxes_vector)
   else:
      flux_solution=np.dot(nullmnp,free_fluxes_vector)
   solution_dict={}
   if n_reaction_dict!=None:
      for x in n_reaction_dict:
          solution_dict[n_reaction_dict[x]]=flux_solution[x]
   return flux_solution, solution_dict


#flux_solver(nullmnp,free_fluxes_vector,n_reaction_dict=n_reaction_dict)



def check_dy(label_model,condition_size_yy_dict,flux_list,condition):
    sumt=0
    #label_model.active_condition=condition
    for size in label_model.size_emu_c_eqn_dict:
        dy=label_model.size_emu_c_eqn_dict[size](condition_size_yy_dict[condition][size],0,condition_size_yy_dict,flux_list,label_model.condition_initial_label_yy_dict, condition)
        #print dy
        for x in dy:
            sumt+=abs(x)  
    return sumt



def c13solver(label_model,mode="fsolve",model=None,check_unity=True,hot_start=False,simulate_infeasible=False,simulate_label=True):
      """
      Function that computes the steady state isotopologue fractions by solving the emu balances model. 
     
      label_model: label_model object
      mode: str, optional
	Either "fsolve" (default) or "ODE". Indicates how should the problem be solved
      fba_mode: str,optional
        Either "fba" (flux balance analysis,default) or "pfba" (parsimonious). Indicates wich method should be used to generate the set of setady state fluxes from the COBRA model
      model: COBRApy model object, optional
           COBRApy model that will be used to compute the steady state flux distribution. If none is provided it will use the constrained_model in the label_model object
      check_unity: bool, optional
	If set to True it will check that isotoplogues fractions for a given EMU are never larger than 1
      hot_start: bool, deprectated
        currently unnused
           
      """
      penalty=0
      label_model.solver_flag="success"
      
      
      label_model.flux_solution,label_model.flux_dict= flux_solver(label_model.flux_solver_nullmnp,label_model.flux_solver_free_fluxes,label_model.flux_solver_n_reaction_dict,label_model)
      label_model.reversible_flux_dict=copy.deepcopy(label_model.flux_dict)
      if simulate_infeasible==False:
         penalty=check_bounds(label_model.flux_dict,label_model.constrained_model,label_model.reaction_emu_dict)
         if penalty>0:
            flux_list=get_fluxes(label_model,model,precision=12)
            return "fail" , label_model.condition_size_yy_dict,penalty
      if simulate_label==False:
         
         flux_list=get_fluxes(label_model,model,precision=12) #This is needed to compute reversible fluxes, which must be considered for total flux, when label weight set to 0
         return "success" , label_model.condition_size_yy_dict,penalty
      if hot_start:
         yy0=label_model.condition_size_yy_dict
      else:
         yy0=label_model.condition_size_yy0_dict
      
      for n_precision in xrange(12,20):
       #print n_precision
       flux_list=get_fluxes(label_model,model,precision=n_precision)
       condition_size_yy_dict=copy.deepcopy(label_model.condition_size_yy0_dict)
       for condition in label_model.initial_label:
         label_model.active_condition=condition
         if mode=="fsolve":
            for size in label_model.size_emu_c_eqn_dict:
                condition_size_yy_dict[condition][size] = fsolve(label_model.size_emu_c_eqn_dict[size] , yy0[condition][size],(0,condition_size_yy_dict,flux_list,label_model.condition_initial_label_yy_dict, condition),maxfev=len(yy0[condition][size])*1000000,xtol=1e-8)
                if check_unity:
                   max_value=condition_size_yy_dict[condition][size].max()
                   min_value=condition_size_yy_dict[condition][size].min() 
                   if max_value>1.001 or min_value<-0.001:
                      label_model.solver_flag="fail"
                      print mode+" failed unity check"+str([max_value,min_value])
                      print list(label_model.variable_vector)
                      print n_precision
                      #error_dump(label_model,e_type="fsolve error") 
                      break
                      
         elif mode=="ode":
            t=np.linspace(0,5000,2) 
            nt=len(t)-1
            for size in label_model.size_emu_c_eqn_dict:
                yy=odeint(label_model.size_emu_c_eqn_dict[size], yy0[condition][size], t,(condition_size_yy_dict,flux_list,label_model.condition_initial_label_yy_dict, condition),mxstep=5000)
                condition_size_yy_dict[condition][size]=yy[nt]
                dy=check_dy(label_model,condition_size_yy_dict,flux_list,condition)
                #total=sum(dy)
                if check_unity:
                   max_value=condition_size_yy_dict[condition][size].max() 
                   min_value=condition_size_yy_dict[condition][size].min() 
                   if max_value>1.001 or min_value<-0.001:
                      label_model.solver_flag="fail"
                      print mode+" failed unity check"+str([max_value,min_value])
                      #error_dump(label_model,e_type="solver error") 
                      break
         dy=check_dy(label_model,condition_size_yy_dict,flux_list,condition)
         if dy>1e-4:
            label_model.solver_flag="fail"
            print mode+" failed with an error of "+str(dy)
         if label_model.solver_flag=="fail":
            print "trying ode"
            if mode!="ode": #If fsolve failes try ode
               label_model.solver_flag, label_model.condition_size_yy_dict,penalty=c13solver(label_model,mode="ode")
            else:
               break
                
       if label_model.solver_flag!="fail":
           break #Break the precision loop if it worked
           #error_dump(label_model,e_type="solver error")  
      label_model.condition_size_yy_dict=condition_size_yy_dict        
      return  label_model.solver_flag, label_model.condition_size_yy_dict,penalty
