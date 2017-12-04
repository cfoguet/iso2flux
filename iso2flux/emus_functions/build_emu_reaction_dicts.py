def build_emu_reaction_dicts(label_model):
    "Builds the dict of the emu reactions"
    "label_model: label_model object"
    label_model.reaction_n_dict={}
    label_model.n_reaction_dict={}
    sorted_reactions=sorted(label_model.reaction_emu_dict.keys())
    for n,reaction in enumerate( sorted_reactions):
        if "_reverse" in reaction:
            if reaction[:-8] in label_model.merged_reactions_reactions_dict:
               reaction_merged=sorted(label_model.merged_reactions_reactions_dict[reaction[:-8]])[0]+"_reverse"
               label_model.reaction_n_dict[reaction_merged]=n
               label_model.n_reaction_dict[n]=reaction_merged
               continue
        if reaction not in label_model.merged_reactions_reactions_dict:
           label_model.reaction_n_dict[reaction]=n
           label_model.n_reaction_dict[n]=reaction
        else: #If it is a merged reaction take the flux of the first one (a merged reaction should only be declared when the member reactions will always have the same flux)
           reaction_merged=sorted(label_model.merged_reactions_reactions_dict[reaction])[0]
           label_model.reaction_n_dict[reaction_merged]=n
           label_model.n_reaction_dict[n]=reaction_merged
    return label_model.reaction_n_dict, label_model.n_reaction_dict
