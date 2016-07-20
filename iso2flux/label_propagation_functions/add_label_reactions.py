from cobra import Model, Reaction, Metabolite
from get_reactants_dict import get_reactants_dict
from build_label_propagation_dict import build_label_propagation_dict
import copy
from ..flux_functions.define_reaction_group import define_reaction_group

def add_label_reactions(label_model,reaction_id,label_propagation={},substrates_dict={},products_dict={},inv=False,split_reaction=False,split_reactions_dict={}):
    model=label_model.irreversible_metabolic_model
    simp_model=label_model.simplified_metabolic_model
    if inv==False:
      if isinstance(reaction_id,list) and len(reaction_id)==1:
	    reaction_id=reaction_id[0]
      if not isinstance(reaction_id,list): 
          if reaction_id not in model.reactions:
             print "Reaction "+reaction_id+" not found"
             return
          reference_reaction=model.reactions.get_by_id(reaction_id)
          if reference_reaction.lower_bound==0.0 and reference_reaction.upper_bound==0.0 and reference_reaction not in label_model.reactions_with_forced_turnover: #If the reaction is inactive do not add it
             if "reflection" not in reference_reaction.notes:
                 print "Reaction "+reaction_id+" ignored as it is inactive"
                 return
             else:
               reflection_reaction_id=reference_reaction.notes["reflection"]
               reflection_reaction=model.reactions.get_by_id(reflection_reaction_id)
               if reflection_reaction.lower_bound==0.0 and reflection_reaction.upper_bound==0.0:
                  print "Reaction "+reaction_id+" ignored as it is inactive"
                  return
      else: #If a list of reactions is added the reactions will be merged
          reference_reaction=process_reaction_list(label_model,reaction_id)
          if reference_reaction==None:
             return
          
          
          
    else: 
         raise ValueError("Error: invisible reactions not supported yet")
    if substrates_dict=={} or products_dict=={}:
       generated_substrates_dict,generated_products_dict=get_reactants_dict(label_model,reaction_id)
    if substrates_dict=={}:
       substrates_dict=generated_substrates_dict
    if products_dict=={}:
          if not inv:
             products_dict=generated_products_dict
          else:
             products_dict=substrates_dict
    if products_dict=={}:
          #raise ValueError("Error: No product found. Input/Output reactions are added autimatically") 
          #label_output_reaction(model,reaction_id)
          return
    if substrates_dict=={}:
          #print("Error: No substrate found.Input/Output reactions are added autimatically")
          #label_input_reaction(model,reaction_id)
          return
          #input
              
    substrates_object_dict={}
    products_object_dict={}
    products=[] 
    reaction_list=[]
    reverse_reaction_list=[]
    simp_metabolites_dict={}
    for substrate_key in substrates_dict:
           if inv==True:
              continue
           isotopomer_object=label_model.met_id_isotopomer_dict[substrates_dict[substrate_key]]
           #substrates_object_dict[substrate_key]=isotopomer_object
           isotopomer_object_id=isotopomer_object.id
           if isotopomer_object_id in simp_model.metabolites:
              metabolite=simp_model.metabolites.get_by_id(isotopomer_object_id)
           else:
              metabolite = Metabolite(isotopomer_object_id,formula='',name=isotopomer_object_id,compartment=isotopomer_object.comp)
           if metabolite in simp_metabolites_dict:
                  simp_metabolites_dict[metabolite]+=-1
           else:
                  simp_metabolites_dict[metabolite]=-1         
    for product_key in  products_dict:
           isotopomer_object=label_model.met_id_isotopomer_dict[products_dict[product_key]]
           products_object_dict[product_key]=isotopomer_object
           isotopomer_object_id=isotopomer_object.id
           if isotopomer_object_id in simp_model.metabolites:
              metabolite=simp_model.metabolites.get_by_id(isotopomer_object_id)
           else:
              metabolite = Metabolite(isotopomer_object_id,formula='',name=isotopomer_object_id,compartment=isotopomer_object.comp)
           if metabolite in simp_metabolites_dict:
                  simp_metabolites_dict[metabolite]+=1
           else:
                  simp_metabolites_dict[metabolite]=1
    if reference_reaction.id in simp_model.reactions: #If a previous entry of the reaction has been added,  remove it
       old_reaction=simp_model.reactions.get_by_id(reference_reaction.id)
       old_reaction.remove_from_model()
    reaction_simp = Reaction(reference_reaction.id)
    reaction_simp.name = reference_reaction.name
    reaction_simp.subsystem = ''
    reaction_simp.lower_bound = 0  # This is the default
    reaction_simp.upper_bound = 1000.  # This is the default
    reaction_simp.add_metabolites(simp_metabolites_dict) 
    if reaction_simp not in simp_model.reactions:
       simp_model.add_reaction(reaction_simp)
    if "reflection" in reference_reaction.notes:
       simp_reverse_metabolites_dict={}
       for met in simp_metabolites_dict:
           simp_reverse_metabolites_dict[met]=-1*simp_metabolites_dict[met]
       
       reverse_reaction_simp = Reaction(reference_reaction.notes["reflection"])
       reverse_reaction_simp.name = reference_reaction.notes["reflection"]
       reverse_reaction_simp.subsystem = ''
       reverse_reaction_simp.lower_bound = 0  # This is the default
       reverse_reaction_simp.upper_bound = 1000.  # This is the default
       reverse_reaction_simp.add_metabolites(simp_reverse_metabolites_dict)
       if reverse_reaction_simp not in simp_model.reactions: 
          simp_model.add_reaction(reverse_reaction_simp)
          reverse_reaction_simp.notes["reflection"]=reaction_simp.id
          reaction_simp.notes["reflection"]=reverse_reaction_simp.id
            
    #if the reaction is monosubstrate and no label propagation is indicated, assume no change ocurs in label (ie 10101->10101) 
    if label_propagation=={} and len(substrates_dict)==1 and len(products_dict)==1:
       substrate_key=substrates_dict.keys()[0]
       product_key=products_dict.keys()[0]
       product=products_object_dict[product_key]
       product_label=[]
       for n in xrange(0,product.n):
           product_label.append([substrate_key,n])
       label_propagation={product_key: product_label}
     
    build_label_propagation_dict(label_model,reference_reaction,label_propagation,substrates_dict,products_dict)

def process_reaction_list(label_model,reaction_list):
	#Find out if they are sequential reactions or paralel reactions
	reference_reactants_list=[]
	is_reaction_group=True
	for n,reaction_id in enumerate(reaction_list):
		reactant_list=[]
		substrates,products=get_reactants_dict(label_model,reaction_id)
		for reactant in substrates.keys()+products.keys():
                        isotopomer=label_model.met_id_isotopomer_dict[reactant]
                        if isotopomer.input==False:
			   reactant_list.append(isotopomer.id)
			
		print reactant_list
		
		if reference_reactants_list==[]:
		   reference_reactants_list=sorted(reactant_list)
		if reference_reactants_list!=sorted(reactant_list):
		   is_reaction_group=False
		   break
	if is_reaction_group:
		reaction_group_dict={}
		reversible=False
		group_id="LABEL_RGROUP_"
		stochiometry={}
		reverse_stochiometry={}
		for n,reaction_id in enumerate(reaction_list):
			group_id+=reaction_id+"_"
			reaction_group_dict[reaction_id]=1
			reference_reaction=label_model.irreversible_metabolic_model.reactions.get_by_id(reaction_id)
			if stochiometry=={}:
				for metabolite in reference_reaction.metabolites:
					stochiometry[metabolite]=reference_reaction.metabolites[metabolite]
					reverse_stochiometry[metabolite]=-reference_reaction.metabolites[metabolite]
			if "reflection" in reference_reaction.notes:
				reverse_reaction_id=reference_reaction.notes["reflection"]
				if reverse_reaction_id in label_model.irreversible_metabolic_model.reactions:
				   reverse_reaction=label_model.irreversible_metabolic_model.reactions.get_by_id(reverse_reaction_id)
				   reversible=True
				   reverse_reaction.remove_from_model()
			reference_reaction.remove_from_model()
                group_id=group_id[:-1]
		define_reaction_group(label_model.metabolic_model,reaction_group_dict,group_reaction_id=group_id,lower_bound=None,upper_bound=None,objective_coefficient=0)
                reference_reaction=Reaction(group_id)
                reference_reaction.add_metabolites(stochiometry)
                if reversible:
			reversible_id=group_id+"_reverse"
			reverse_reaction=Reaction(reversible_id)
			reverse_reaction.add_metabolites(reverse_stochiometry)
			reverse_reaction.notes["reflection"]=group_id 
			label_model.irreversible_metabolic_model.add_reaction(reverse_reaction)
			reference_reaction.notes["reflection"]=reversible_id
                label_model.irreversible_metabolic_model.add_reaction(reference_reaction)
	else:
	  merged_reaction_id=""
          for reaction_to_add in reaction_list:
            merged_reaction_id+=reaction_to_add+"_"
            merged_reaction_id=merged_reaction_id[:-1] #Remove "_" at the end
            label_model.merged_reactions_reactions_dict[merged_reaction_id]=reaction_list
            reference_reaction=Reaction(merged_reaction_id)
            reversible_flag=True
            reference_reaction.name="" 
            reversible_merged_reaction_list=[]
            for reaction_to_add in reaction_list:
              label_model.reaction_merged_reactions_dict[reaction_to_add]=merged_reaction_id
              reaction=label_model.irreversible_metabolic_model.reactions.get_by_id(reaction_to_add)
              reference_reaction.name+=reaction.name+" "  
              if "reflection" not in reaction.notes and reaction_to_add not in label_model.reversible_co2_reactions: #If of one of the merged reactions is irrvesible the new reaction is irreversible
                 reversible_flag=False
              else:
                 reversible_merged_reaction_list.append(reaction)
              if reaction.lower_bound==0.0 and reaction.upper_bound==0.0:  #If one of the reactions merged is inactive do not add it
                 return None
            if reversible_flag==True:
             #reflection=Reaction(merged_reaction_id+"_reverse")
             reference_reaction.notes[reflection]=merged_reaction_id+"_reverse"
            else:
             print reversible_merged_reaction_list
             print "Merged reaction removed reflections"
             #Make al reactions of the merged reactions irreversible if one is irreversible
             for forward_reaction in reversible_merged_reaction_list:
                 reflection=label_model.irreversible_metabolic_model.reactions.get_by_id(forward_reaction.id+"_reverse")
                 reflection.remove_from_model()
                 if "reflection" in forward_reaction.notes:
                 	del forward_reaction.notes["reflection"]
	return reference_reaction			   
		       
