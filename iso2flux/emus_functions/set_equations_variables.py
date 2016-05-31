import copy
import numpy as np
def set_equations_variables(label_model,force_balance=True,set_initial_label=True):
    label_model.condition_size_timecourse_dict={}
    label_model.condition_size_transposed_timecourse_dict={}
    label_model.condition_size_yy0_dict={}
    
    n=0
    if label_model.input_n_dict=={}:
      label_model.input_m0_list=[]
      for isotopomer_id in sorted(label_model.id_isotopomer_object_dict.keys()):
        if label_model.id_isotopomer_object_dict[isotopomer_id].input==True and isotopomer_id in label_model.metabolite_emu_dict:
           for emu in sorted(label_model.metabolite_emu_dict[isotopomer_id]):
               for n_mi in label_model.emu_dict[emu]["mid"]:
                   mi=label_model.emu_dict[emu]["mid"][n_mi]
                   label_model.input_n_dict[mi]=n
                   if n_mi==0:
                      label_model.input_m0_list.append(mi)
                   n+=1
    for condition in label_model.initial_label:
        label_model.condition_size_timecourse_dict[condition]={}
        label_model.condition_size_transposed_timecourse_dict[condition]={}
        label_model.condition_initial_label_yy_dict[condition]=np.zeros(len(label_model.input_n_dict),dtype=np.float64)
        for mi0 in label_model.input_m0_list:
            n=label_model.input_n_dict[mi0]
            label_model.condition_initial_label_yy_dict[condition][n]=1
        label_model.condition_size_yy0_dict[condition]={}
        for size in label_model.emu_size_dict:
            label_model.condition_size_yy0_dict[condition][size]=np.zeros(shape=(len(label_model.size_variable_dict[size])),dtype=np.float64)
            if force_balance==False:
              for emu in label_model.emu_size_dict[size]:
                  mi0=label_model.emu_dict[emu]["mid"][0]
                  if mi0 in label_model.size_variable_dict[size]: #Used to check if Force Balance=True
                     n=label_model.size_variable_dict[size][mi0]
                     label_model.condition_size_yy0_dict[condition][size][n]=1
        for size in range(1,max(label_model.emu_size_dict.keys())): #If some size is not present create empty variables for it
            if size not in label_model.condition_size_yy0_dict[condition]:
               label_model.condition_size_yy0_dict[condition][size]=np.zeros(shape=(1),dtype=np.float64)
        label_model.condition_size_yy_dict=copy.deepcopy(label_model.condition_size_yy0_dict)
        if set_initial_label==True:
          for mi in label_model.initial_label[condition]:
             value=label_model.initial_label[condition][mi]
             n=label_model.input_n_dict[mi]
             label_model.condition_initial_label_yy_dict[condition][n]=value
    #label_model.size_dy_dict=copy.deepcopy(label_model.size_yy0_dict)
