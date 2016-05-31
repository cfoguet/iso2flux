from cobra import Model, Reaction, Metabolite
import copy
from set_reflections import set_reflections
def expand_emu_models(label_model):
   for model_size in label_model.size_model_dict:
      expanded_emu_model=model = Model('expanded_emu_model (size'+str(model_size)+")")
      mid_emu_dict={}
      mid_weigh_dict={}
      for metabolite in label_model.size_model_dict[model_size].metabolites:
           size=label_model.emu_dict[metabolite.id]["size"]
           local_met_dict={}
           for n in range(0,size+1):
               met=metabolite.id+"_m"+str(n)
               local_met_dict[n]=met
               mid_emu_dict[met]=metabolite
               mid_weigh_dict[n]=met
               label_model.expanded_emu_dict[met]=metabolite.id  
           label_model.emu_dict[metabolite.id]["mid"]=local_met_dict
      for original_reaction in label_model.size_model_dict[model_size].reactions:
          new_reactions=[]
          subs=[]
          prod=None
          for metabolite in original_reaction.metabolites: 
              coef=original_reaction.metabolites[metabolite]
              if coef==1 or (coef>0 and "input" in original_reaction.id) :
                 prod=metabolite.id
              elif coef==-1:
                 subs.append(metabolite.id)
              elif coef==-2:
                 subs=[metabolite.id,metabolite.id]
              elif coef<0 and ("output" in original_reaction.id):
                 subs.append(metabolite.id)  
              else: #The coeficient can only be different from 1 or 2 in input or output reactions   
                 raise ValueError("Error: wrong coeffient in "+original_reaction.id)
          if subs==[]: #input reactions, assumed to generate only m0
             print("adding "+original_reaction.id+" as input")
             if label_model.emu_dict[prod]["mid"][0] in expanded_emu_model.metabolites:
                prod_m=expanded_emu_model.metabolites.get_by_id(label_model.emu_dict[prod]["mid"][0])
             else:
                prod_m=Metabolite(label_model.emu_dict[prod]["mid"][0],formula='', name=label_model.emu_dict[prod]["mid"][0], compartment= label_model.size_model_dict[model_size].metabolites.get_by_id( prod ).compartment )           
             metabolites_dict={prod_m:coef}
             new_reaction_id=original_reaction.id+"_m0"
             new_reaction = Reaction(new_reaction_id)
             new_reaction.name =  original_reaction.name+"m0 " 
             new_reaction.subsystem = ''
             new_reaction.lower_bound = 0  # This is the default
             new_reaction.upper_bound = 1000.  # This is the default
             new_reaction.add_metabolites(metabolites_dict) 
             expanded_emu_model.add_reaction(new_reaction)
             new_reactions.append(new_reaction.id)
             print new_reaction.reaction
             
             
          else: #if the reaction has substrate
              print subs   
              sub1_dict=label_model.emu_dict[subs[0]]["mid"]
              if prod!=None:
                 prod_dict=label_model.emu_dict[prod]["mid"]
                 print prod_dict
              for m in sub1_dict:
                  if sub1_dict[m] in expanded_emu_model.metabolites:
                     sub1_m=expanded_emu_model.metabolites.get_by_id(sub1_dict[m])
                  else:
                     sub1_m=Metabolite(sub1_dict[m],formula='',name=sub1_dict[m], compartment=label_model.size_model_dict[model_size].metabolites.get_by_id(subs[0]).compartment) 
                  if len(subs)==2:
                     sub2_dict=label_model.emu_dict[subs[1]]["mid"] 
                     for m2 in sub2_dict:
                         if sub2_dict[m2] in expanded_emu_model.metabolites:
                           sub2_m=expanded_emu_model.metabolites.get_by_id(sub2_dict[m2])
                         else:
                           sub2_m=Metabolite(sub2_dict[m2],formula='',name=sub2_dict[m2], compartment=label_model.size_model_dict[model_size].metabolites.get_by_id(subs[1]).compartment) 
                         
                         weigh=m+m2
                         print ["weigh",weigh]
                         if prod_dict[weigh] in expanded_emu_model.metabolites:
                            prod_m=expanded_emu_model.metabolites.get_by_id(prod_dict[weigh])
                         else:
                            prod_m=Metabolite(prod_dict[weigh],formula='', name=prod_dict[weigh], compartment= label_model.size_model_dict[model_size].metabolites.get_by_id( prod ).compartment )                       
                            
                         local_met_dict[n]=met
                         mid_emu_dict[met]=metabolite
                         mid_weigh_dict[n]=met
                         metabolites_dict={sub1_m:-1,sub2_m:-1,prod_m:1}
                         new_reaction_id=original_reaction.id+"_m"+str(m)+"_m"+str(m2) 
                         new_reaction = Reaction(new_reaction_id)
                         new_reaction.name =  original_reaction.name+"_m"+str(m)+"_m"+str(m2) 
                         new_reaction.subsystem = ''
                         new_reaction.lower_bound = 0  # This is the default
                         new_reaction.upper_bound = 1000.  # This is the default
                         new_reaction.add_metabolites(metabolites_dict) 
                         expanded_emu_model.add_reaction(new_reaction)
                         new_reactions.append(new_reaction.id)
                  elif len(subs)==1 and prod!=None:
                     weigh=m
                     if prod_dict[weigh] in expanded_emu_model.metabolites:
                            prod_m=expanded_emu_model.metabolites.get_by_id(prod_dict[weigh])
                     else:
                            prod_m=Metabolite(prod_dict[weigh],formula='', name=prod_dict[weigh], compartment= label_model.size_model_dict[model_size].metabolites.get_by_id( prod ).compartment )                                 
                     metabolites_dict={sub1_m:-1,prod_m:1}
                     new_reaction_id=original_reaction.id+"_m"+str(m)
                     new_reaction = Reaction(new_reaction_id)
                     new_reaction.name =  original_reaction.name+"_m"+str(m) 
                     new_reaction.subsystem = ''
                     new_reaction.lower_bound = 0  # This is the default
                     new_reaction.upper_bound = 1000.  # This is the default
                     new_reaction.add_metabolites(metabolites_dict) 
                     expanded_emu_model.add_reaction(new_reaction)
                     new_reactions.append(new_reaction.id)
                     
                  elif prod==None: #output
                     weigh=m                                 
                     metabolites_dict={sub1_m:coef}
                     new_reaction_id=original_reaction.id+"_m"+str(m)
                     new_reaction = Reaction(new_reaction_id)
                     new_reaction.name =  original_reaction.name+"_m"+str(m)
                     new_reaction.subsystem = ''
                     new_reaction.lower_bound = 0  # This is the default
                     new_reaction.upper_bound = 1000.  # This is the default
                     new_reaction.add_metabolites(metabolites_dict) 
                     expanded_emu_model.add_reaction(new_reaction)
                     new_reactions.append(new_reaction.id)
          for new_reaction in new_reactions: 
            if original_reaction.id in label_model.emu_reaction_expanded_dict:
                label_model.emu_reaction_expanded_dict[original_reaction.id].append(new_reaction)
            else:
                label_model.emu_reaction_expanded_dict[original_reaction.id]=[new_reaction]
            label_model.expanded_reaction_emu_dict[new_reaction]=original_reaction.id
                            
                     
      label_model.size_expanded_model_dict[model_size]=expanded_emu_model
   set_reflections(label_model.size_expanded_model_dict)     
   return(label_model.size_expanded_model_dict)
