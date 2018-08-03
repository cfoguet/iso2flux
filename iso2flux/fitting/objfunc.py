from ..emus_functions.c13solver import c13solver 
from get_objective_function import get_objective_function as get_xi2
def objfunc(label_model,variables,label_weight=1,simulate_label=True,mode="fsolve",target_flux_dict=None,max_chi=1e16,max_flux=1e16,flux_penalty_dict={},flux_weight=1,verbose=False,abort_when_above_max_flux=False,label_unfeasible_penalty=1e9,flux_unfeasible_penalty=1e6,simulate_infeasible=False,update_variables=True):
    #label_model.variable_vector=variables
    abort=False #When True Label, or fluxes wont be simulated
    reaction=""
    variable_list_dict=label_model.variable_list_dict
    flux_obj=label_obj=chi2_score=flux_score=0
    fltarget_obj=0
    obj_dict={}
    if update_variables:
     for n, x1 in enumerate(variables):
        vtype=variable_list_dict[n]["type"]
        if vtype=="flux":
           label_model.flux_solver_free_fluxes[variable_list_dict[n]["n"]]=round(x1,12)
        elif vtype=="turnover":
           label_model.turnover_flux_dict[variable_list_dict[n]["reaction"]]["v"]=round(max(x1,0),12) 
        elif  vtype=="non_indepndent_flux":
           correct_flux_fun=variable_list_dict[n]["fun"]
           correct_flux_fun(label_model.flux_solver_nullmnp,label_model.flux_solver_free_fluxes,x1)
           #print "r0"+str(variables[-1])
    succes_flag,b,penalty=c13solver(label_model,simulate_label=simulate_label,mode=mode,simulate_infeasible=simulate_infeasible)
    if succes_flag=="fail": 
       obj=round((penalty+1),3)*1e6
       flux_obj=1e3
       flux_score=1e3
       label_obj=1e3
       chi2_score=1e3
       #print obj
    else:
       for flux in label_model.flux_dict:
                if flux in flux_penalty_dict:
                   flux_score+=abs(label_model.flux_dict[flux])*flux_penalty_dict[flux]
       #obj+=sum([abs(label_model.flux_dict[x]) for x in label_model.flux_dict])*flux_penalty
       #print flux_score
       if flux_score>max_flux:
            flux_obj=(flux_score+1)*flux_unfeasible_penalty
            if abort_when_above_max_flux:
               abort=True
       else:
            flux_obj=flux_score*flux_weight
       if simulate_label and not abort:
	    chi2_score,b,c=get_xi2(label_model,output=False) #TODO rename function
            label_obj=chi2_score*label_weight
            """if chi2_score>max_chi:
               label_obj=(chi2_score+1)*label_unfeasible_penalty
            else:
               label_obj=chi2_score*label_weight"""
       else:
           chi2_score=1e3
           flux_score=1e3
           flux_obj=1e3
           label_obj=1e3
       if target_flux_dict!=None:
          #print "r "+str(label_model.reversible_flux_dict[reaction])
          #print variables[-4:]
          reaction=target_flux_dict["reaction"]
          direction=target_flux_dict["dir"] #maximum if tis must be maximized, minimum if it must be minimized
          factor=target_flux_dict["factor"] #weitght given to the flux value
          irreversible_model_flag=target_flux_dict["irreversible_model"]
          if irreversible_model_flag in (None,False):
           if "max" in direction.lower():
             ub=target_flux_dict["ub"]
             fltarget_obj+=(ub-label_model.reversible_flux_dict[reaction])*factor #Minimize the distance between the upper bound and the flux value
           elif "min" in direction.lower():
             lb=target_flux_dict["lb"]
             fltarget_obj+=(label_model.reversible_flux_dict[reaction]-lb)*factor
           else:
               raise ValueError('Field dir in target_flux_dict must be either max or min')
          else:
            reaction=reaction.replace("_forward","")
            if "max" in direction.lower():
             ub=target_flux_dict["ub"]
             fltarget_obj+=(ub-label_model.flux_dict[reaction])*factor #Minimize the distance between the upper bound and the flux value
            elif "min" in direction.lower():
             lb=target_flux_dict["lb"]
             fltarget_obj+=(label_model.flux_dict[reaction]-lb)*factor
       obj=flux_obj+label_obj+fltarget_obj
    obj_dict={"flux_obj":flux_obj,"label_obj":label_obj,"fltarget_obj":fltarget_obj,"chi2_score":chi2_score,"flux_score":flux_score}
    if verbose and label_obj>0:
       output="Chi2="+str(round(label_obj,3))
       if flux_obj>0:
          output+=" flux="+str(round(flux_obj,3))
       if fltarget_obj>0:
          output+=" "+reaction+"="+str(round(label_model.reversible_flux_dict[reaction],3))+" "+direction.lower()
       print output
    #print obj_dict        
    return obj, obj_dict 
