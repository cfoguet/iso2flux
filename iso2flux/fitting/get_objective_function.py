def get_objective_function(label_model,output=False,condition_size_yy_dict=None,rsm="dynamic",force_balance=True):
    """
    Computes computes square deviation normalized by standard deviation between the experimentally determined and simulated isotoplogues fractions
    label_model: label_model object
    output:boolean, optional
		If true it will print in the terminal the square deviation associated to each measurment 
    condition_size_yy_dict: dict,optional
		Provides the simulated isotopologue distributions. If it is None it will take the distribution stored in the label_model object
    rsm: string,optional:
         Indicates how the fucntion should deal with experimental measurments that have the m/Sm flag. It can either be "dynamic" (only compute m/Sm when experimental m0>simulated m0), always (always compute m/Sm) or never (never compute m/Sm). Default is dynamic	  
    """
    force_balance=label_model.force_balance
    if condition_size_yy_dict==None:
       condition_size_yy_dict=label_model.condition_size_yy_dict
    label_model.objective_function_dict={}
    label_model.condition_simulation_results_dict={}
    label_model.obj=0
    for condition in label_model.condition_size_yy_dict:
        label_model.condition_simulation_results_dict[condition]={}
        mz_dict={}
        if label_model.solver_flag!="success":
           return 99999999, {},{}
        label_model.objective_function_dict[condition]={}
        for emu in  label_model.experimental_dict[condition]:
            do_rsm=False
            total_label=0
            size=label_model.emu_dict[emu]["size"]
            for n in label_model.experimental_dict[condition][emu]:
                mi=label_model.emu_dict[emu]["mid"][n]
                """if n==0 and emu in label_model.rsm_list:
                   continue"""  
                if mi in label_model.size_variable_dict[size]: #Do not attempt to get simulated m0 if balance is forced
                   n_mi=label_model.size_variable_dict[size][mi]
                   sim=label_model.condition_size_yy_dict[condition][size][n_mi]
                   if emu not in label_model.condition_simulation_results_dict[condition]:
                      label_model.condition_simulation_results_dict[condition][emu]={}
                   label_model.condition_simulation_results_dict[condition][emu][n]=sim
                   if n!=0:
                      total_label+=sim
                else:
                   if emu not in label_model.condition_simulation_results_dict[condition]:
                      label_model.condition_simulation_results_dict[condition][emu]={}
                   label_model.condition_simulation_results_dict[condition][emu][n]=0.0 
            if force_balance==True:
               label_model.condition_simulation_results_dict[condition][emu][0]=1-total_label
            if emu in label_model.rsm_list and rsm!="never": #TODO change it to be a property of experimental data dict
               if label_model.experimental_dict[condition][emu][0]["m"]>label_model.condition_simulation_results_dict[condition][emu][0] or rsm=="always": #Only do m/SM ratio when simulated m0 is higuer than experimental m0
                  do_rsm=True
                  for n in label_model.experimental_dict[condition][emu]:
                      
                      if do_rsm==True and n==0:
                         continue
                      if n>0:
                        if total_label==0 or n not in label_model.condition_simulation_results_dict[condition][emu]:
                           label_model.condition_simulation_results_dict[condition][emu][n]=0
                        else:      
                           label_model.condition_simulation_results_dict[condition][emu][n]/=total_label  
            for n in label_model.experimental_dict[condition][emu]:#label_model.condition_simulation_results_dict[condition][emu]:
                if do_rsm==True and n==0:
                   continue
                sim=label_model.condition_simulation_results_dict[condition][emu][n]
                if do_rsm==True:
                   exp_m=label_model.experimental_dict[condition][emu][n]["m"]/(1-label_model.experimental_dict[condition][emu][0]["m"])
                   exp_sd=max(label_model.experimental_dict[condition][emu][n]["sd"]/(1-label_model.experimental_dict[condition][emu][0]["m"]),label_model.minimum_sd)   
                else:
                   exp_m=label_model.experimental_dict[condition][emu][n]["m"]
                   exp_sd=max(label_model.experimental_dict[condition][emu][n]["sd"],label_model.minimum_sd) 
                chi=round(pow((sim-exp_m)/exp_sd,2),3)
                if emu not in label_model.objective_function_dict[condition]:
                   label_model.objective_function_dict[condition][emu]={}
                label_model.objective_function_dict[condition][emu][n]=chi
                label_model.obj+=chi
                if output==True:
                   if emu in label_model.data_name_emu_dict:
                      name=label_model.data_name_emu_dict[emu]
                   else:
                      name=emu
                   string=name+" m"+str(n)
                   if do_rsm:
                      string+="/ms"
                   string+=": "+str(sim)+" "+str(exp_m)+"+-"+str(exp_sd)+" chi:"+str(chi)
                   print string
    return label_model.obj,label_model.objective_function_dict,label_model.condition_simulation_results_dict
      
