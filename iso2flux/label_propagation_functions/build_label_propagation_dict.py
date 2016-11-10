def build_label_propagation_dict(label_model,reaction,label_propagation,substrates_dict,products_dict):
    product_propagation_dict=label_propagation
    substrate_propagation_dict={}
    print product_propagation_dict 
    for substrate_key in  substrates_dict:
        substrate=substrates_dict[substrate_key] 
        substrate_propagation_dict[substrate_key]=[]
        isotopomer_object=label_model.met_id_isotopomer_dict[substrate]
        for x in xrange(isotopomer_object.n):
            substrate_propagation_dict[substrate_key].append(["carb"]) #This ensure that if the postion does not appear in the products is treated as decarboxylation
    for product_key in products_dict:
        for n,source in enumerate(product_propagation_dict[product_key]):
            if source==["carb"] or source==["carboxylation"]:
               continue
            substrate_key=source[0]
            substrate_carbon=source[1]
            substrate_propagation_dict[substrate_key][substrate_carbon]=[product_key,n]
    label_model.label_propagation_dict[reaction.id]={"sub_id_dict":substrates_dict,"prod_id_dict":products_dict,"prod_label_dict":product_propagation_dict,"sub_label_dict":substrate_propagation_dict} 
    if "reflection" in reaction.notes:
       label_model.label_propagation_dict[reaction.notes["reflection"]]={"sub_id_dict":products_dict, "prod_id_dict":substrates_dict, "prod_label_dict":substrate_propagation_dict,  "sub_label_dict":product_propagation_dict}   
