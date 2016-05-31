import numpy as np
def get_f_fluxes(label_model,met_con_dict={},steady_state=True):
    label_model.flux_list=np.zeros(shape=len(label_model.reaction_n_dict),dtype=np.float64)
    """for reaction in label_model.split_reactions_dict:
        flux=label_model.flux_dict[reaction]
        for split_reaction in label_model.split_reactions_dict[reaction]:
            label_model.flux_dict[split_reaction]=flux"""
    for reaction_id in label_model.reaction_n_dict:
        reaction_n=label_model.reaction_n_dict[reaction_id]
        f_flux=label_model.flux_dict[reaction_id]
        if not steady_state:
           raise ValueError("Error: non steady state not supported yet")           
           """reaction=model.reactions.get_by_id(reaction_id)
           for metabolites in reaction._metabolites:
               if metabolites in metabolite_isotopomers_dict:
                  isotopomer_object=metabolite_isotopomers_dict[metabolites]
                  if not isotopomer_object.pool:
                     con=met_con_dict[metabolites.id]
                  else:
                     con=0
                     for met in isotopomer_object.ref_met:
                         con+=met_con_dict[met.id]*vol_dict[met.compartment]
                     con/=vol_dict[isotopomer_object.comp]
                  coef=reaction._metabolites[metabolites]
                  if isotopomer_object.comp!="e" and isotopomer_object.comp!="ec":  
                     if coef<=0 and coef.is_integer():  
                        f_flux*=pow(con,coef)
                     elif coef<=0:
                        f_flux/=con"""
        label_model.flux_list[reaction_n]=f_flux
    return label_model.flux_list
