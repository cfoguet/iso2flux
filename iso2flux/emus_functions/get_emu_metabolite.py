from cobra import  Metabolite

def get_emu_metabolite(emu_dict,label_model):
    """
    Creates a cobra metabolite object for the new emu. Additionally it will add additional information if the emu is symmetryc
    label_model: label_model object
    emu_dict: dict
        dict that contains the ID of the metabolite (or group of metabolites) and the carbon range
    """
    #print emu_dict
    met_id=emu_dict["met_id"]
    carbons=emu_dict["carbons"]
    iso_object=label_model.id_isotopomer_object_dict[met_id]
    emuid="emu_"+met_id+"_"
    carbon_range_string=""
    for carbon in carbons:  
        carbon_range_string+=str(carbon)
    if iso_object.symm==True:
       #build a symetryc dic
       symm_dict={}
       forward_range=range(1,iso_object.n+1)
       reverse_range=range(iso_object.n,0,-1)
       for x in forward_range:
           symm_dict[x]=reverse_range[x-1]
       #print symm_dict
       symm_carbons=[]
       for carbon in carbons:
           symm_carbons.append(symm_dict[carbon])
       #print symm_carbons
       symm_carbons=sorted(symm_carbons) 
       emu_dict["symm_carbons"]=symm_carbons
       if symm_carbons!=carbons: #Check if they are not equal: 
          #Identfy the lower range, which will be written first in the metabolite id
          symm_carbon_range_string="" 
          for carbon in symm_carbons:
              symm_carbon_range_string+=str(carbon)
           
          if symm_carbons[0]<carbons[0]:
             emuid+=symm_carbon_range_string+"_and_"+carbon_range_string
            
          else:
             emuid+=carbon_range_string+"_and_"+symm_carbon_range_string
       else:
          emu_dict["symm_carbons"]=[]
          emuid+=carbon_range_string
    else:
          emuid+=carbon_range_string 
    if emuid in label_model.emu_model.metabolites:
       emu_met=label_model.emu_model.metabolites.get_by_id(emuid)
       present_flag=True
    else:
       comp=label_model.simplified_metabolic_model.metabolites.get_by_id(emu_dict["met_id"]).compartment
       emu_met=Metabolite(emuid,formula='',name=emuid, compartment=comp)
       present_flag=False
       met_id=emu_dict["met_id"]
       label_model.emu_metabolite_dict[emu_met]=met_id
       if met_id not in label_model.metabolite_emu_dict:
          label_model.metabolite_emu_dict[met_id]=[]
       label_model.metabolite_emu_dict[met_id].append(emu_met.id) 
    return(emu_met,emu_dict,present_flag)
