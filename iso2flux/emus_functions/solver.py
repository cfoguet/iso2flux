from scipy.optimize import broyden1
from scipy.optimize import broyden2
from scipy.optimize import fsolve
from scipy.integrate import odeint
import json 
import cobra
import time

from cobra.flux_analysis.parsimonious import optimize_minimal_flux
from ..flux_functions.apply_ratios import apply_ratios
from ..flux_functions.get_fluxes import get_fluxes

if __name__=="__main__":
    from iso2flux.flux_functions.apply_ratios import apply_ratios
    from iso2flux.flux_functions.get_fluxes import get_fluxes

import numpy as np
t=np.linspace(0,1000,2)
nt=len(t)-1 #position where the final solution will be found

def check_dy(label_model,condition):
    sumt=0
    label_model.active_condition=condition
    for size in label_model.size_emu_c_eqn_dict:
        dy=label_model.size_emu_c_eqn_dict[size](label_model.condition_size_yy_dict[condition][size],0,label_model)
        for x in dy:
            sumt+=abs(x)  
    return sumt 

def error_dump(label_model,e_type=None):
    #Fucntion used for debbuging
    dump_dict={"parameter_dict":label_model.parameter_dict,"flux_dict":label_model.flux_dict,"turnover_dict":label_model.turnover_flux_dict,"e_type":e_type}
    fname="zzzzzz"+str(time.time())
    with open(fname+".json", 'w') as fp:
         json.dump(dump_dict, fp)
    cobra.io.write_sbml_model(label_model.constrained_model, fname+"sbml")

def solver(label_model,mode="fsolve",fba_mode="fba",model=None,check_unity=True,hot_start=False):
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
    label_model.solver_flag="success"
    if model==None:
       model=label_model.constrained_model
    apply_ratios(model,label_model.ratio_dict)
    try:
      if fba_mode.lower()=="fba":
             model.optimize(tolerance_feasibility=label_model.lp_tolerance_feasibility)
      elif fba_mode.lower()=="pfba":
             #convert_model_to_irreversible(label_model.constrained_model)
             optimize_minimal_flux(model, already_irreversible=False)
             #revert_model_to_reversible(label_model.constrained_model, update_solution=True)
    except:
             label_model.solver_flag="fail"    
             error_dump(label_model,e_type="error in optimize") 
             return  label_model.solver_flag, label_model.condition_size_yy_dict
    if model.solution.status =="optimal":
      if hot_start:
         yy0=label_model.condition_size_yy_dict
      else:
         yy0=label_model.condition_size_yy0_dict
      for n_precision in xrange(10,16):
       get_fluxes(label_model,model,precision=n_precision)
       for condition in label_model.initial_label:
         label_model.active_condition=condition
         if mode=="fsolve":
            for size in label_model.size_emu_c_eqn_dict:
                label_model.condition_size_yy_dict[condition][size] = fsolve(label_model.size_emu_c_eqn_dict[size] , yy0[condition][size],(0,label_model))
                if check_unity:
                   max_value=label_model.condition_size_yy_dict[condition][size].max()
                   min_value=label_model.condition_size_yy_dict[condition][size].min() 
                   if max_value>1.001 or min_value<-0.001:
                      label_model.solver_flag="fail"
                      print mode+" failed unity check"+str([max_value,min_value])
                      error_dump(label_model,e_type="fsolve error") 
                      break
                      
         elif mode=="ode":
            for size in label_model.size_emu_c_eqn_dict:
                yy=odeint(label_model.size_emu_c_eqn_dict[size], yy0[condition][size], t,(label_model,),mxstep=5000)
                label_model.condition_size_yy_dict[condition][size]=yy[nt]
                dy=label_model.size_emu_c_eqn_dict[size](yy[nt],0,label_model)
                total=sum(dy)
                if total>1e-6: #If it is not a strict steady state use fsolve to refine the results
                   label_model.condition_size_yy_dict[condition][size] =  fsolve(label_model.size_emu_c_eqn_dict[size] , yy0[condition][size],(0,label_model))
                if check_unity:
                   max_value=label_model.condition_size_yy_dict[condition][size].max() 
                   min_value=label_model.condition_size_yy_dict[condition][size].min() 
                   if max_value>1.001 or min_value<-0.001:
                      label_model.solver_flag="fail"
                      print mode+" failed unity check"+str([max_value,min_value])
                      #error_dump(label_model,e_type="solver error") 
                      break
         dy=check_dy(label_model,condition)
         if dy>1e-4:
            label_model.solver_flag="fail"
            print mode+" failed with an error of "+str(dy)
         if label_model.solver_flag=="fail":
            if mode!="ode": #If fsolve failes try ode
               label_model.solver_flag, label_model.condition_size_yy_dict=solver(label_model,mode="ode")
            else:
               break
                
       if label_model.solver_flag!="fail":
           break #Break the precision loop if it worked
           error_dump(label_model,e_type="solver error")               
    else:
      #remove 
      label_model.solver_flag="fail"    
      error_dump(label_model,e_type="unfeasible") 
    return  label_model.solver_flag, label_model.condition_size_yy_dict
