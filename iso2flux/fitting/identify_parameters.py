import random
import copy
from apply_parameters import apply_parameters 
from ..flux_functions.apply_ratios import apply_ratios
from cobra.flux_analysis.variability import flux_variability_analysis
import math
from ..misc.round_functions import round_up,round_down

def identify_parameters(label_model,add_to_parameter_dict=False,fraction_of_optimum=0,change_threshold=0.1,parameter_precision=0.01,max_d=0.1,add_turnover=True,turnover_upper_bound=100,excluded_turnovers=[],key_reactions=[],debug=False):
    parameter_dict={}
    """if change_threshold<parameter_precision:
       parameter_precision=change_threshold/10"""
    precision=int(-1*(math.log10(parameter_precision)))
    #print label_model.constrained_model.reactions.get_by_id("biomass").lower_bound
    reaction_list=[]
    original_model=copy.deepcopy(label_model.constrained_model)
    apply_parameters(label_model,parameter_list=[]) #Paremeters will be applied with the default precision im label_model
    original_objective_list=[]
    for reaction in label_model.constrained_model.objective: 
            fva=flux_variability_analysis(label_model.constrained_model,reaction_list=[reaction], fraction_of_optimum=fraction_of_optimum,tolerance_feasibility=label_model.lp_tolerance_feasibility)
            minimum=max(round_down(fva[reaction.id]["minimum"],precision),reaction.lower_bound)
            maximum=min(round_up(fva[reaction.id]["maximum"],precision),reaction.upper_bound)
            value=round((minimum+maximum)/2,precision) 
            reaction.objective_coefficient=0
            original_objective_list.append(reaction)
            parameter_dict[reaction.id]={"v":value,"lb":minimum,"ub":maximum ,"type":"flux value","reactions":[reaction.id],"max_d":max_d,"original_lb":reaction.lower_bound,"original_ub":reaction.upper_bound,"original_objective_coefficient":reaction.objective_coefficient}
            apply_parameters(label_model,parameter_dict=parameter_dict,parameter_precision=parameter_precision,parameter_list=[reaction.id])
    #model=label_model.constrained_model
    for reaction in label_model.reaction_n_dict: #Find reactions that propagate label
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
    print reaction_list
    original_mfva=flux_variability_analysis(label_model.constrained_model,reaction_list=reaction_list, fraction_of_optimum=0,tolerance_feasibility=label_model.lp_tolerance_feasibility)
    additional_parameters_list=[]
    first_time=True
    for i in xrange(0,10000):
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
                   original_maximum=round_up(original_mfva[reaction_id]["maximum"],precision)
                   original_minimum=round_down(original_mfva[reaction_id]["minimum"],precision)
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
                   original_maximum=round_up(original_mfva[reaction_id]["maximum"],precision)
                   original_minimum=round_down(original_mfva[reaction_id]["minimum"],precision)
      print ["max variation",max_variation]
      reaction=label_model.constrained_model.reactions.get_by_id(reaction_id)
      value=round((4*minimum+1*maximum)/5,precision) #Weighted average
      parameter_dict[reaction_id]={"v":value,"lb":original_minimum,"ub":original_maximum ,"type":"flux value","reactions":[reaction_id],"max_d":max_d,"original_lb":original_model.reactions.get_by_id(reaction_id).lower_bound,"original_ub":original_model.reactions.get_by_id(reaction_id).upper_bound}
      apply_parameters(label_model,parameter_dict=parameter_dict,parameter_precision=parameter_precision,parameter_list=[reaction_id])
      if debug:
         print [reaction_id,mfva[reaction_id],[reaction.lower_bound,reaction.upper_bound]]
      """reaction.lower_bound=value
      reaction.upper_bound=value"""
      #print([reaction.id,mfva[reaction_id],value,model.optimize()])
      #print(reaction_id+" "+str(parameter_dict[reaction_id]))
    print ("%s free fluxes found"%(len(parameter_dict))) 
    print  sorted(parameter_dict.keys() )  
    if add_turnover==True:
       for turnover in label_model.turnover_flux_dict:
        if turnover not in excluded_turnovers or (turnover+"_turnover") in label_model.parameter_dict: 
           parameter_dict[turnover+"_turnover"]={"v":label_model.turnover_flux_dict[turnover],"lb":0,"ub":turnover_upper_bound,"max_d":max_d,"type":"turnover","reactions":[turnover]}
    
    if add_to_parameter_dict==True:
       apply_parameters(label_model,parameter_dict=parameter_dict) #Paremeters will be applied with the default precision im label_model
       """for reaction in original_objective_list:
           original_reaction=original_model.reactions.get_by_id(reaction.id)
           reaction.objective_coefficient=0.0#=original_reaction.objective_coefficient
           if fraction_of_optimum==0: 
              reaction.lower_bound=original_reaction.lower_bound
              reaction.upper_bound=original_reaction.upper_bound"""
       label_model.parameter_dict.update(parameter_dict)   
    else:
       label_model.constrained_model=original_model
    return parameter_dict    

