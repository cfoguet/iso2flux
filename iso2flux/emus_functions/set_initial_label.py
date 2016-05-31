import numpy as np
def set_initial_label(metabolite_id,label_model,label_patterns=[[[1,1,0,0,0,0],0.5]],condition="condition",total_concentration=1):
    isotopomer_object=label_model.met_id_isotopomer_dict[metabolite_id]
    isotopomer_object.input=True
    isotopomer_id=label_model.met_id_isotopomer_dict[metabolite_id].id
    labelled_positions=[]
    if condition not in label_model.initial_label:
       label_model.initial_label[condition]={}
    for label_pattern in label_patterns:
       labelled_ratio=label_pattern[1]
       for n, x in enumerate(label_pattern[0]):
           if x==1:
              labelled_positions.append(n+1)
       for emu in label_model.metabolite_emu_dict[isotopomer_id]:
           emu_carbons=label_model.emu_dict[emu]["carbons"]
           emu_size=label_model.emu_dict[emu]["size"]
           weight=0
           for position in labelled_positions:
               if position in emu_carbons:
                  weight+=1
           if weight>0:
              emu_mi=label_model.emu_dict[emu]["mid"][weight]
              #emu_mi0=label_model.emu_dict[emu]["mid"][0]
              label_model.initial_label[condition][emu_mi]=total_concentration*labelled_ratio
              #emu_fractions_dict[emu_mi0]=total_concentration-emu_fractions_dict[emu_mi]
    #Set unlabelled ratio
    for emu in label_model.metabolite_emu_dict[isotopomer_id]: 
        total_labelled=0
        for weight in label_model.emu_dict[emu]["mid"]:
            emu_mi=label_model.emu_dict[emu]["mid"][weight]
            if weight==0:
               continue
            elif emu_mi in label_model.initial_label[condition]:
               total_labelled+=label_model.initial_label[condition][emu_mi]
        emu_mi0=label_model.emu_dict[emu]["mid"][0]
        label_model.initial_label[condition][emu_mi0]=total_concentration-total_labelled   
    return label_model.initial_label
