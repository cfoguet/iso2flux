import random
import copy
from apply_parameters import apply_parameters 
from cobra.flux_analysis.variability import flux_variability_analysis
from ..flux_functions.apply_ratios import apply_ratios
from identify_free_parameters import identify_free_fluxes
from get_objective_function import get_objective_function
from apply_parameters import apply_parameters
from ..emus_functions.solver import solver
import math
from ..misc.round_functions import round_up,round_down

def sampling(label_model,n=100,fraction_of_optimum=0,output_emu_list=None,max_turnover=100,parameter_dict={},fba_mode="fba",parameter_precision=0.001,gui=None,change_threshold=0.01):
    print parameter_precision
    precision=int(-1*(math.log10(parameter_precision)))
    #print ["hello",label_model.constrained_model.reactions.get_by_id("biomass").lower_bound]
    original_model=copy.deepcopy(label_model.constrained_model)
    original_turnover=copy.copy(label_model.turnover_flux_dict)
    apply_parameters(label_model,parameter_dict)
    #model=label_model.constrained_model
    #try:
    if fraction_of_optimum!=0:
       for reaction in label_model.constrained_model.objective: 
            fva=flux_variability_analysis(label_model.constrained_model,reaction_list=[reaction], fraction_of_optimum=fraction_of_optimum,tolerance_feasibility=label_model.lp_tolerance_feasibility)
            reaction.lower_bound=round_down(fva[reaction.id]["minimum"],precision)
            reaction.upper_bound=round_up(fva[reaction.id]["maximum"],precision)
            reaction.objective_coefficient=0
    if output_emu_list==None:
       output_emu_list=label_model.data_name_emu_dict.keys()
    #Prepare output
    output_dict={}
    for condition in label_model.condition_initial_label_yy_dict:
        output_dict[condition]={}
        for emu in output_emu_list:
            size=label_model.emu_dict[emu]["size"]
            if size not in output_dict[condition]:
               output_dict[condition][size]={}
            for mi in label_model.emu_dict[emu]["mid"]:
                output_dict[condition][size][label_model.emu_dict[emu]["mid"][mi]]=[]    
    #identfy the turnover fluxes that are already part of parameters as they will be excluded 
    apply_ratios(label_model.constrained_model,label_model.ratio_dict)
    #free_parameters=identify_free_fluxes(label_model,parameter_dict={},fraction_of_optimum=0,change_threshold=change_threshold,parameter_precision=parameter_precision,max_d=0.1,key_reactions=[],original_fva=None,debug=False)
    """if len (free_parameters)==0:
       return output_dict"""
    #print "%s free parameters found"%(len(free_parameters))
    #turnover_parameter_list=[]
    """flux_value_parameter_list=[]
    for parameter in free_parameters:
        parameter_dict=free_parameters[parameter]
        if parameter_dict["type"]=="turnover":
           turnover_parameter_list.append(parameter)
        elif parameter_dict["type"]=="flux value":
           flux_value_parameter_list.append(parameter)"""
    #Save the original lower and upper_bounds and turnover 
    #Generate the samples
    free_parameters={}
    for i in range(0,n):
         free_parameters=identify_free_fluxes(label_model,parameter_dict=free_parameters,fraction_of_optimum=0,change_threshold=change_threshold,parameter_precision=parameter_precision,max_d=0.1,key_reactions=[],original_fva=None,debug=False)
         for reaction_id in free_parameters:
            reaction=label_model.constrained_model.reactions.get_by_id(reaction_id)
            original_reaction=original_model.reactions.get_by_id(reaction_id)
            reaction.lower_bound=original_reaction.lower_bound
            reaction.upper_bound=original_reaction.upper_bound
        
         working_flux_value_parameter_list=random.sample(free_parameters, len(free_parameters))
         #try:
         print free_parameters
         for reaction_id in working_flux_value_parameter_list: 
            fva=flux_variability_analysis(label_model.constrained_model, reaction_list=[reaction_id],fraction_of_optimum=0,tolerance_feasibility=label_model.lp_tolerance_feasibility)
            if fva[reaction_id]["maximum"]-fva[reaction_id]["minimum"]>parameter_precision:
               new_value=random.uniform(fva[reaction_id]["minimum"],fva[reaction_id]["maximum"])
               free_parameters[reaction_id]["v"]=new_value
               print [reaction_id,new_value]
               apply_parameters(label_model,parameter_dict=free_parameters,parameter_precision=parameter_precision,parameter_list=[reaction_id])
         for turnover in label_model.turnover_flux_dict:
             if (turnover+"_turnover") not in parameter_dict: 
                lb=label_model.turnover_flux_dict[turnover]["lb"]
                ub=label_model.turnover_flux_dict[turnover]["lb"]
                new_value=round(random.uniform(lb,ub),precision) 
                label_model.turnover_flux_dict[turnover]["v"]=new_value
         #apply_parameters(label_model,parameter_dict=free_parameters,parameter_list=turnover_parameter_list,parameter_precision=parameter_precision)
         print label_model.constrained_model.optimize()
         a,b=solver(label_model,mode="fsolve",fba_mode=fba_mode)
         print a
         #Store output
         for condition in output_dict:
            for size in output_dict[condition]:
                for mi in output_dict[condition][size]:
                    if mi in label_model.size_variable_dict[size]:
                       n_emu=label_model.size_variable_dict[size][mi]
                       mi_value=round(label_model.condition_size_yy_dict[condition][size][n_emu],4)
                       output_dict[condition][size][mi].append(mi_value)
         if gui!=None:
           #get_objective_function(label_model)
           gui.update_label_sampling()
           gui.root.update_idletasks()
           a,b,c=get_objective_function(label_model,force_balance=label_model.force_balance,output=False)               
            
         print "sample %s of %s..."%((i+1),n)
         #Restore original values for flux constraints
         #except:
          # continue
    label_model.turnover_flux_dict=original_turnover
    label_model.constrained_model=original_model
    #print ["hello",label_model.constrained_model.reactions.get_by_id("biomass").lower_bound]
    #apply_parameters(label_model,parameter_dict=None,parameter_precision=parameter_precision)
    a,b=solver(label_model,mode="fsolve",fba_mode=fba_mode)
    get_objective_function(label_model)
    #except:
    """print "Error"
    label_model.turnover_flux_dict=original_turnover
    label_model.constrained_model=original_model"""
    #print ["hello",label_model.constrained_model.reactions.get_by_id("biomass").lower_bound]
    return output_dict     





#label_model.constrained_model=copy.deepcopy(label_model.metabolic_model)
