from cobra import Model, Reaction, Metabolite
import copy
def get_reactants_dict(label_model,reaction_id,check_count=True):
    model=label_model.irreversible_metabolic_model 
    substrate_dict={}
    product_dict={}
    isotopomer_substrate_coef_dict={}
    isotopomer_product_coef_dict={}
    coef_dict={} 
    isotopomer_coef_dict={}
    if not isinstance(reaction_id,list):
       reaction_id_list=[reaction_id]
    else:
       reaction_id_list=reaction_id
    for reaction_id in reaction_id_list:
        reaction=model.reactions.get_by_id(reaction_id)
        for metabolite in reaction.metabolites:
            if metabolite.id not in label_model.metabolite_id_isotopomer_id_dict:
               continue 
            coef=reaction.metabolites[metabolite]
            if metabolite.id in coef_dict:
               coef_dict[metabolite.id]+=coef
            else:
               coef_dict[metabolite.id]=coef
    #print coef_dict
    for metabolite in coef_dict:
        coef=coef_dict[metabolite]
        isotopomer=label_model.metabolite_id_isotopomer_id_dict[metabolite]
        
        if coef<0:
           substrate_dict[metabolite]=metabolite
           if isotopomer in isotopomer_substrate_coef_dict:
              isotopomer_substrate_coef_dict[isotopomer]-=coef
           else:
              isotopomer_substrate_coef_dict[isotopomer]=-coef
        elif coef>0:
           product_dict[metabolite]=metabolite
           if isotopomer in isotopomer_product_coef_dict:
              isotopomer_product_coef_dict[isotopomer]+=coef
           else:
              isotopomer_product_coef_dict[isotopomer]=coef
    count=0
    for isotopomer in isotopomer_substrate_coef_dict:
        if isotopomer in isotopomer_product_coef_dict:
           coef_sub=isotopomer_substrate_coef_dict[isotopomer]
           coef_prod=isotopomer_product_coef_dict[isotopomer]
           if coef_sub==coef_prod: 
              if len(label_model.isotopomer_id_metabolite_id_dict[isotopomer])>1:
               for met in label_model.isotopomer_id_metabolite_id_dict[isotopomer]:
                   if met in substrate_dict:
                      del(substrate_dict[met])
                   if met in product_dict:
                      del(product_dict[met])
           if coef_sub>coef_prod: #Input reactions
             count+=1
             product_dict={}
           if coef_sub<coef_prod: #Output reaction
             count+=1
             substrate_dict={}
      
    if count>1 and check_count: #If several "inputs or outputs" are detected they migh not really be inputs and outputs and the produc_dict and substrate_dict will need to be build manually
       raise("Exception: Substrates and products could not be identfied uniquelly, please enter substrates_dict and product_dict manually")  
    return substrate_dict, product_dict

