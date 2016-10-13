from cobra import Model, Reaction, Metabolite
from cobra.flux_analysis.variability import flux_variability_analysis
import copy
def remove_impossible_emus(label_model):
   """
   Removes emus that cannot with labelled substrates used in the experiment
    
   """
   #size_expanded_model2_dict={} 
   possible_isotopomers=[]
   impossible_isotopomers=[]
   f=open("impossible_isotopologues","w")
   #If some metabolites are marked as input but have no label described assume they all are m0 and add them to initial label
   for isotopomer in label_model.isotopomer_object_list:
       if isotopomer.input==True:
          if isotopomer.id not in label_model.metabolite_emu_dict:
             continue
          emus=label_model.metabolite_emu_dict[isotopomer.id]
          for condition in label_model.initial_label:
              present=False
              for emu in emus:
                for mi in label_model.emu_dict[emu]["mid"]:
                  if label_model.emu_dict[emu]["mid"][mi] in label_model.initial_label[condition]:
                     present=True
                     break
                if not present:
                 for mi in label_model.emu_dict[emu]["mid"]:
                     if mi==0:
                        label_model.initial_label[condition][label_model.emu_dict[emu]["mid"][mi]]=1
                     else:
                        label_model.initial_label[condition][label_model.emu_dict[emu]["mid"][mi]]=0
   present_initial_label=[]
   for condition in label_model.initial_label:  #For simplicity the label model is shared across all conditions
        for label in label_model.initial_label[condition]:
          if label_model.initial_label[condition][label]>0.0:
             if label not in present_initial_label:
                present_initial_label.append(label)                            
   for model_size in sorted(label_model.size_expanded_model_dict.keys()): #needs to be done in order of size
      expanded_emu_model=label_model.size_expanded_model_dict[model_size]
      label_exchange_reactions=[]
      for metabolite in expanded_emu_model.metabolites: 
          reaction = Reaction("EX_"+metabolite.id)
          reaction.name = 'EX'+metabolite.name
          reaction.subsystem = 'Test'
          reaction.lower_bound = 0
          reaction.upper_bound = 1000.  # This is the default
          reaction.add_metabolites({metabolite:-1})        
          expanded_emu_model.add_reaction(reaction) 
          label_exchange_reactions.append(reaction)
      for condition in label_model.initial_label:  #For simplicity the label model is shared across all conditions
        for label in label_model.initial_label[condition]:
             ex_id="EX_"+label
             if ex_id in expanded_emu_model.reactions:
                if label in present_initial_label:
                   expanded_emu_model.reactions.get_by_id(ex_id).lower_bound=-1000
                else: 
                    expanded_emu_model.reactions.get_by_id(ex_id).lower_bound=0
                    expanded_emu_model.reactions.get_by_id(ex_id).upper_bound=0
        """for isotopomer in label_model.isotopomer_object_list:
            if isotopomer.input==True:
             if isotopomer.id in label_model.metabolite_emu_dict:
               for emu in label_model.metabolite_emu_dict[isotopomer.id]:
                   mid0=label_model.emu_dict[emu]["mid"][0]
                   ex_id="EX_"+mid0
                   if ex_id in expanded_emu_model.reactions:
                      expanded_emu_model.reactions.get_by_id(ex_id).lower_bound=-1000"""
      if model_size==1.0:
         import cobra
         cobra.io.write_sbml_model(expanded_emu_model,"size1.sbml")
      
      for isotopomer in possible_isotopomers: #If a isotopomer is produced in a smaller size it should appear as present in this size
          ex_id="EX_"+isotopomer
          if ex_id in expanded_emu_model.reactions:
             #print ex_id 
             expanded_emu_model.reactions.get_by_id(ex_id).lower_bound=-1000
      fva=flux_variability_analysis(expanded_emu_model,reaction_list=label_exchange_reactions,fraction_of_optimum=0)
      for x in label_exchange_reactions: 
          if fva[x.id]["maximum"]>0.00001 or fva[x.id]["minimum"]<-0.00001 : 
             if x.id[3:] not in possible_isotopomers: #prevent repeated entries
                possible_isotopomers.append(x.id[3:])
             #print x.id
          else: 
             if x.id[3:] not in impossible_isotopomers: #prevent repeated entries
                impossible_isotopomers.append(x.id[3:])
             #print x.id
      #Remove reactions where those isotopomers participate: 
      reactions_to_remove=[]
      for x in  label_exchange_reactions:
          x.remove_from_model(remove_orphans=True)
      impossible_isotopomers_object_list=[]    
      for isotopomer_id in impossible_isotopomers:
          if isotopomer_id not in expanded_emu_model.metabolites:
             continue  
          isotopomer=expanded_emu_model.metabolites.get_by_id(isotopomer_id)
          for x in isotopomer._reaction:
              if x not in reactions_to_remove and x.metabolites[isotopomer]<1:
                 reactions_to_remove.append(x)
          isotopomer.remove_from_model(method='subtractive')
      for reaction in expanded_emu_model.reactions:
           if len(reaction.metabolites)==0 and reaction not in reactions_to_remove:
              reactions_to_remove.append(reaction) 
      for reaction in reactions_to_remove:
         print "removing "+reaction.id
         f.write("removing "+reaction.id+"\n")
         emu_reaction=label_model.expanded_reaction_emu_dict[reaction.id]
         label_model.emu_reaction_expanded_dict[emu_reaction].remove(reaction.id)
         del label_model.expanded_reaction_emu_dict[reaction.id]               
         reaction.remove_from_model(remove_orphans=True)
      #size_expanded_model2_dict[model_size]=expanded_emu_model
      #labek_model=size_expanded_model_dict[model_size]=expanded_emu_model 
              
   for emu in impossible_isotopomers:
       f.write(emu+"\n")
   f.close()
   #Remove the reactions where only m0 is transmited and replace them with input reactions
   remove_m0_reactions(label_model,impossible_isotopomers,possible_isotopomers)
   return label_model.size_expanded_model_dict



def remove_m0_reactions(label_model,impossible_isotopomers,possible_isotopomers):
    #This will remove all reactions that have only m0 substrates and products
    only_m0_list=[]
    for emu in label_model.emu_dict:
        only_m0=True
        for mi in label_model.emu_dict[emu]["mid"]:
            """if label_model.id_isotopomer_object_dict[label_model.emu_dict[emu]["met_id"]].input==True:
               only_m0=False
               break"""
            if mi>0 and label_model.emu_dict[emu]["mid"][mi] in possible_isotopomers:
               only_m0=False
               break
        if only_m0:
           only_m0_list.append(emu)
    print only_m0_list
    reactions_to_remove=[]
    for emu in sorted(only_m0_list):
        mi0_id=label_model.emu_dict[emu]["mid"][0]
        for size in label_model.size_expanded_model_dict:
          expanded_model=label_model.size_expanded_model_dict[size]
          if mi0_id in expanded_model.metabolites:
             mi0=expanded_model.metabolites.get_by_id(mi0_id)
             mi0.remove_from_model(method='subtractive')
        """for reaction in mi0.reactions:
            coef=reaction.metabolites[mi0]
            if len(reaction.metabolites)==1: #If its an output or input reacton we can remove it
               if reaction not in reactions_to_remove:
                  reactions_to_remove.append(reaction)
            else:
                  reaction.add_metabolites({mi0:-1*coef}) #If the reaction consumes the metabolite we will remove the metabolite from the reaction"""
    for size in label_model.size_expanded_model_dict:
        for reaction in label_model.size_expanded_model_dict[size].reactions:
            if len(reaction.metabolites)==0:
                reactions_to_remove.append(reaction)         
    f=open("remove_m0_reactions","w")
    for emu in only_m0_list:
       f.write(emu+"\n")  
    for reaction in reactions_to_remove:
       f.write(reaction.id+" "+reaction.reaction+"\n")
    f.close()
    for reaction in reactions_to_remove:
         print "removing"+reaction.id
         emu_reaction=label_model.expanded_reaction_emu_dict[reaction.id]
         del label_model.emu_reaction_expanded_dict[emu_reaction]
         del label_model.expanded_reaction_emu_dict[reaction.id]
         metabolic_reactions=copy.copy(label_model.emu_reaction_dict[emu_reaction])
         if not isinstance(metabolic_reactions,list):
            metabolic_reactions=[metabolic_reactions]
         for metabolic_reaction in metabolic_reactions:
             label_model.reaction_emu_dict[metabolic_reaction].remove(emu_reaction)
             if len(label_model.reaction_emu_dict[metabolic_reaction])==0:
                del label_model.reaction_emu_dict[metabolic_reaction]
         del label_model.emu_reaction_dict[emu_reaction]
         label_model.reaction_emu_dict             
         reaction.remove_from_model(remove_orphans=True)

             
               
        
           
            
