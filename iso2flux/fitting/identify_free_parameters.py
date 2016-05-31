import random
import copy
from apply_parameters import apply_parameters 
from clear_parameters import clear_parameters 
from ..flux_functions.apply_ratios import apply_ratios
from cobra.flux_analysis.variability import flux_variability_analysis
import math
from ..misc.round_functions import round_up,round_down

def identify_free_fluxes(label_model,parameter_dict={},fraction_of_optimum=0,change_threshold=0.1,parameter_precision=0.01,max_d=0.1,key_reactions=[],original_fva=None,restore_model=False,debug=False):
    precision=int(-1*(math.log10(parameter_precision)))
    #print label_model.constrained_model.reactions.get_by_id("biomass").lower_bound
    reaction_list=[]
    original_model=copy.deepcopy(label_model.constrained_model)
    #apply_parameters(label_model,parameter_dict=parameter_dict) #Paremeters will be applied with the default precision im label_model
    #original_objective_list=[]
    #model=label_model.constrained_model
    for reaction in label_model.reactions_propagating_label:#label_model.reaction_n_dict: #Find reactions that propagate label
        if reaction in label_model.merged_reactions_reactions_dict:
           for merged_reaction in label_model.merged_reactions_reactions_dict[reaction]:
               if "_reverse" not in merged_reaction:
                  reaction_list.append(merged_reaction) 
        elif "_reverse" not in reaction and "RATIO_" not in reaction:
           reaction_list.append(reaction)
    for reaction in key_reactions:
        if reaction in label_model.constrained_model.reactions and reaction not in reaction_list:
           reaction_list.append(reaction) 
    #Classify reactions 
    if original_fva==None:
       original_fva=flux_variability_analysis(label_model.constrained_model,fraction_of_optimum=0,tolerance_feasibility=label_model.lp_tolerance_feasibility, reaction_list=reaction_list)
    #original_mfva=flux_variability_analysis(label_model.constrained_model,reaction_list=reaction_list, fraction_of_optimum=0,tolerance_feasibility=label_model.lp_tolerance_feasibility)
    additional_parameters_list=[]
    first_time=True
    done_reactions=[]
    for i in xrange(0,10000):
      #print model.optimize()
      free_ex_reactions_list=[]
      free_reactions_list=[]
      priortizided_reactions_list=[]
      mfva=flux_variability_analysis(label_model.constrained_model,reaction_list=reaction_list, fraction_of_optimum=0,tolerance_feasibility=label_model.lp_tolerance_feasibility)
      for reaction in mfva:
          """if reaction not in label_model.reaction_emu_dict and reaction not in label_model.reaction_merged_reactions_dict:
             print "Error"
             continue"""
          maximum=mfva[reaction]["maximum"]
          minimum=mfva[reaction]["minimum"]
          if abs(maximum-minimum)>change_threshold:
             if reaction in key_reactions:
                priortizided_reactions_list.append(reaction)
                """elif "EX_" in reaction:
                free_ex_reactions_list.append(reaction)
                #print free_ex_reactions_list"""
             else:
                free_reactions_list.append(reaction)
      #print free_reactions_list 
      #print free_ex_reactions_list
      if  free_reactions_list==[] and priortizided_reactions_list==[]:
          break
      elif priortizided_reactions_list!=[]:
           #Find the reaction with more variation allowed
           max_variation=0
           for reaction in priortizided_reactions_list:
               variation=mfva[reaction]["maximum"]-mfva[reaction]["minimum"]
               if variation>=max_variation:
                   max_variation=variation
                   reaction_id=reaction
                   maximum=mfva[reaction_id]["maximum"]
                   minimum=mfva[reaction_id]["minimum"]  
                   original_maximum=round_down(original_fva[reaction_id]["maximum"],precision)
                   original_minimum=round_up(original_fva[reaction_id]["minimum"],precision)
           """elif free_ex_reactions_list!=[]:
           #Find the reaction with more variation allowed
           max_variation=0
           for reaction in free_ex_reactions_list:
               variation=mfva[reaction]["maximum"]-mfva[reaction]["minimum"]
               if variation>=max_variation:
                   max_variation=variation
                   reaction_id=reaction
                   maximum=mfva[reaction_id]["maximum"]
                   minimum=mfva[reaction_id]["minimum"]  
                   original_maximum=round(original_mfva[reaction_id]["maximum"],precision)
                   original_minimum=round(original_mfva[reaction_id]["minimum"],precision)"""
      elif free_reactions_list!=[]:
               max_variation=0
               for reaction in free_reactions_list:
                 variation=mfva[reaction]["maximum"]-mfva[reaction]["minimum"]
                 if variation>=max_variation:
                   max_variation=variation
                   reaction_id=reaction 
                   maximum=mfva[reaction_id]["maximum"]
                   minimum=mfva[reaction_id]["minimum"]  
                   original_maximum=original_fva[reaction_id]["maximum"]#round_down(original_fva[reaction_id]["maximum"],precision)
                   original_minimum=original_fva[reaction_id]["minimum"]#round_up(original_fva[reaction_id]["minimum"],precision)
      print ["max variation",max_variation]
      reaction=label_model.constrained_model.reactions.get_by_id(reaction_id)
      value=random.uniform(minimum,maximum)
      #value=round((4*minimum+1*maximum)/5,precision) #Weighted average
      parameter_dict[reaction_id]={"v":value,"lb":original_minimum,"ub":original_maximum ,"type":"flux value","reactions":[reaction_id],"max_d":max_d,"original_lb":original_model.reactions.get_by_id(reaction_id).lower_bound,"original_ub":original_model.reactions.get_by_id(reaction_id).upper_bound,"original_objective_coefficient":0.0}
      apply_parameters(label_model,parameter_dict=parameter_dict,parameter_precision=parameter_precision,parameter_list=[reaction_id])
      done_reactions.append(reaction_id)
      if debug:
         print [reaction_id,mfva[reaction_id],[reaction.lower_bound,reaction.upper_bound]]
      """reaction.lower_bound=value
      reaction.upper_bound=value"""
      #print([reaction.id,mfva[reaction_id],value,model.optimize()])
      #print(reaction_id+" "+str(parameter_dict[reaction_id]))
    print ("%s free fluxes found"%(len(done_reactions)))
    """if add_turnover==True:
       for turnover in label_model.turnover_flux_dict:
        if turnover not in excluded_turnovers or (turnover+"_turnover") in label_model.parameter_dict: 
           parameter_dict[turnover+"_turnover"]={"v":label_model.turnover_flux_dict[turnover],"lb":0,"ub":turnover_upper_bound,"max_d":max_d,"type":"turnover","reactions":[turnover]}"""
    if restore_model==True:   
        label_model.constrained_model=original_model
    return parameter_dict

def identify_free_parameters(label_model,parameter_dict={},n_samples=50,add_to_model=True,parameter_precision=1e-5,change_threshold=1e-3,fraction_of_optimum=0,max_d=0.1,key_reactions=[],add_turnover=True,excluded_turnovers=[],turnover_upper_bound=None,debug=True):
 free_parameters=copy.deepcopy(parameter_dict)
 precision=int(-1*(math.log10(parameter_precision)))
 apply_parameters(label_model,parameter_dict=free_parameters)
 for reaction in label_model.constrained_model.objective: 
     fva=flux_variability_analysis(label_model.constrained_model,reaction_list=[reaction], fraction_of_optimum=fraction_of_optimum,tolerance_feasibility=label_model.lp_tolerance_feasibility)
     print fva
     minimum=round_down(fva[reaction.id]["minimum"],precision)#max(round_up(fva[reaction.id]["minimum"],precision),reaction.lower_bound)
     maximum=round_up(fva[reaction.id]["maximum"],precision)#min(round_down(fva[reaction.id]["maximum"],precision),reaction.upper_bound)
     value=((minimum+maximum)/2) 
     #original_objective_list.append(reaction)
     free_parameters[reaction.id]={"v":value,"lb":minimum,"ub":maximum ,"type":"flux value","reactions":[reaction.id],"max_d":max_d,"original_lb":reaction.lower_bound,"original_ub":reaction.upper_bound,"original_objective_coefficient":reaction.objective_coefficient}
     reaction.lower_bound=minimum
     reaction.upper_bound=maximum
     reaction.objective_coefficient=0
     #print [model.optimize(),minimum,maximum]
 original_fva=flux_variability_analysis(label_model.constrained_model,fraction_of_optimum=0,tolerance_feasibility=label_model.lp_tolerance_feasibility)
 #print free_parameters
 apply_parameters(label_model,parameter_dict=free_parameters,parameter_precision=parameter_precision)
 #print model.optimize()
 for i in range(0,n_samples):
  print "sample "+str(i)
  try:
   #original_model=copy.deepcopy(label_model.constrained_model)
   free_parameters=identify_free_fluxes(label_model,parameter_dict=free_parameters,fraction_of_optimum=fraction_of_optimum ,change_threshold=change_threshold, parameter_precision=parameter_precision,max_d=max_d,key_reactions=key_reactions,original_fva=original_fva,debug=debug)
   flux_value_parameter_list=[]
   for parameter in free_parameters:
        local_parameter_dict=free_parameters[parameter]
        if local_parameter_dict["type"]=="flux value":
           flux_value_parameter_list.append(parameter)
           #Save the original lower and upper_bounds and turnover 
           #Generate the samples
   for parameter in flux_value_parameter_list:
       for reaction_id in free_parameters[parameter]["reactions"]:
           reaction=label_model.constrained_model.reactions.get_by_id(reaction_id)
           reaction.lower_bound=max(free_parameters[parameter]["original_lb"],free_parameters[parameter]["lb"])
           reaction.upper_bound=min(free_parameters[parameter]["original_ub"],free_parameters[parameter]["ub"])
   working_flux_value_parameter_list=random.sample(flux_value_parameter_list, len(flux_value_parameter_list))
   for reaction_id in working_flux_value_parameter_list: 
    
            fva=flux_variability_analysis(label_model.constrained_model, reaction_list=[reaction_id],fraction_of_optimum=0,tolerance_feasibility=label_model.lp_tolerance_feasibility)
            #if fva[reaction_id]["maximum"]-fva[reaction_id]["minimum"]>parameter_precision:
            new_value=random.uniform(fva[reaction_id]["minimum"],fva[reaction_id]["maximum"])
            free_parameters[reaction_id]["v"]=new_value
            apply_parameters(label_model,parameter_dict=free_parameters,parameter_precision=parameter_precision,parameter_list=[reaction_id])
            reaction=label_model.constrained_model.reactions.get_by_id(reaction_id)
            #print[fva,new_value,[reaction.lower_bound,reaction.upper_bound],label_model.constrained_model.optimize()]
            
  except:
            print "Error caught"
            continue
 if add_turnover==True:
       """if turnover_upper_bound==None: 
          turnover_upper_bound=0
          for reaction_id in label_model.reaction_n_dict:
              if reaction_id in original_fva:
                  turnover_upper_bound=max(original_fva[reaction_id]["maximum"],turnover_upper_bound)
          turnover_upper_bound=round_up(turnover_upper_bound,0)""" 
       for turnover in label_model.turnover_flux_dict:
        if turnover not in excluded_turnovers and (turnover+"_turnover") not in free_parameters and label_model.turnover_flux_dict[turnover]["ub"]!=label_model.turnover_flux_dict[turnover]["lb"]:
           ub=label_model.turnover_flux_dict[turnover]["ub"] 
           lb=label_model.turnover_flux_dict[turnover]["lb"]
           v=label_model.turnover_flux_dict[turnover]["v"] 
           free_parameters[turnover+"_turnover"]={"v":v,"lb":lb,"ub":ub,"max_d":max_d,"type":"turnover","reactions":[turnover]}
 if add_to_model==True:
     label_model.parameter_dict=free_parameters
     apply_parameters(label_model,parameter_dict=free_parameters)
 else:
     clear_parameters(label_model,parameter_dict=copy.deepcopy(free_parameters))
 return free_parameters   





