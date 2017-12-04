from get_f_fluxes import get_f_fluxes
import copy
def get_fluxes(label_model,met_con_dict={},split_reactions_dict={},steady_state=True,model=None,precision=10):
    met_con_dict={}
    if model!=None:
       label_model.flux_dict=copy.deepcopy(model.solution.x_dict)
       original_fluxes=model.solution.x_dict
    else:
       original_fluxes=copy.deepcopy(label_model.reversible_flux_dict)
    reaction_list=[]
    for reaction_id in label_model.reaction_n_dict:#label_model.reaction_emu_dict:
         if reaction_id in label_model.metabolic_model.reactions:
            reaction_list.append(reaction_id)
         elif "_reverse" in reaction_id and reaction_id[:-8] not in label_model.reaction_n_dict:
             reaction_list.append(reaction_id[:-8]) #Even If only the reverse reaction is present in reaction_n_dict add the original reaction to the flux dict
         """elif reaction_id in label_model.merged_reactions_reactions_dict:
              for merged_reaction in label_model.merged_reactions_reactions_dict[reaction_id]:
                  reaction_list.append(label_model.metabolic_model.reactions.get_by_id(merged_reaction))"""
    for reaction_id in reaction_list:
         #reaction=model.reactions.get_by_id(reaction_id)
         if reaction_id not in  label_model.flux_dict: #Zero flux reactions are removed when computing the nullmatrix, but they can still carry flux through turnover 
            label_model.flux_dict[reaction_id]=0
         label_model.flux_dict[reaction_id]=round(label_model.flux_dict[reaction_id],precision)
         
         if label_model.flux_dict[reaction_id]<0:
            label_model.flux_dict[reaction_id+"_reverse"]=-label_model.flux_dict[reaction_id]
            label_model.flux_dict[reaction_id]=0
         elif reaction_id+"_reverse" in label_model.reaction_n_dict:
            label_model.flux_dict[reaction_id+"_reverse"]=0
         
         
    for reaction_id in label_model.turnover_flux_dict: 
        """if reaction in label_model.reaction_merged_reactions_dict:
           #If reaction is part of a merged reaction that is irreversible do not apply turnover
           merged_reaction_id=label_model.reaction_merged_reactions_dict[reaction]
           if "reflection" not in label_model.simplified_metabolic_model.reactions.get_by_id(merged_reaction_id).notes:
              continue"""
        #net_flux=model.solution.x_dict[reaction]
        net_flux=round(original_fluxes[reaction_id],precision)
        if net_flux>=0:
           label_model.flux_dict[reaction_id+"_reverse"]=label_model.turnover_flux_dict[reaction_id]["v"]
           label_model.flux_dict[reaction_id]=net_flux+label_model.turnover_flux_dict[reaction_id]["v"]
        else:
            label_model.flux_dict[reaction_id]=label_model.turnover_flux_dict[reaction_id]["v"]
            label_model.flux_dict[reaction_id+"_reverse"]=-net_flux+label_model.flux_dict[reaction_id]
    flux_list=get_f_fluxes(label_model,met_con_dict=met_con_dict,steady_state=steady_state)   
    return flux_list#label_model.flux_dict  
