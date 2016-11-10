from cobra.core.Reaction import Reaction
from cobra.core.Metabolite import Metabolite
from numpy import abs
import copy

def add_turnover_metabolites(cobra_model, metabolite_id_list=[], epsilon=1e-6,label_model=None):
    #Function based on the work by : """Schmidt BJ1, Ebrahim A, Metz TO, Adkins JN, Palsson B, Hyduke DR. GIM3E: condition-specific models of cellular metabolism developed from metabolomics and expression data Bioinformatics. 2013 Nov 15;29(22):2900-8. doi: 10.1093/bioinformatics/btt493. Epub 2013 Aug 23."""
    
    
    """ NOTE: Model must first be converted to irreversible!
    This entry creates a corresponding turnover metabolite
    that ensures flux through the metabolite of interest.
    
    Arguments:
     cobra_model: the model to be updated.
     metabolite_id_list: list of model metabolites for
                         which to add a turnover metabolites. If one element of the list contains several metabolites it will add a common turnover metabolites to all the metabolites in the list forcing at least one of the metabolites to be present. Example: [["pyr_c","pyr_m","pyr_e"]] will force pyruvate to be present regardless of the compartment 
     epsilon: minimal flux to force through turnover metabolites.
      recommend 1.01 X solver_tolerance    
    
    
    """
    
    maxupper_bound=max([x.upper_bound for x in cobra_model.reactions])
    # Set the minimum flux for metabolites equal to some factor larger than the solver's tolerance
    the_min_flux = epsilon
    turnover_metabolites = []
    sink_reactions = []
    full_metabolite_id_list=copy.copy(metabolite_id_list)
    single_metabolites=[]
    group_metabolites=[]
    #If a label model object  is provided make sure all experimentally measured metabolites (or at least one of the metabolites in the pool) is produced
    if label_model!=None:
       emus=[]
       for condition in label_model.experimental_dict:
           for emu in label_model.experimental_dict[condition]:
               if emu not in emus:
                  emus.append(emu)
       measured_metabolite_dict={}
       print emus
       for emu in emus:
           iso_id=label_model.emu_dict[emu]["met_id"]
           #isotopomer_object=label_model.id_isotopomer_object_dict[iso_id]
           if len(label_model.isotopomer_id_metabolite_id_dict[iso_id])>1:
              group_metabolites.append(label_model.isotopomer_id_metabolite_id_dict[iso_id])
           else:
              single_metabolites.append(str(label_model.isotopomer_id_metabolite_id_dict[iso_id][0]))
    for metabolite in metabolite_id_list:
        if isinstance(metabolite,list) and len(metabolite)>1:
           group_metabolites.append(metabolite)
        elif isinstance(metabolite,list):
           single_metabolites.apppend(metabolite[0])
        elif metabolite not in single_metabolites:  
            single_metabolites.apppend(metabolite)
    print group_metabolites
    for metabolites in  group_metabolites:
              met_id=""
              for metabolite in metabolites:
                  met_id+=metabolite+"_"
              met_id=str(met_id[:-1])
              #print isotopomer_object.ref_met
              sum_abs_source_reaction_bounds = 0
              v_metabolite = Metabolite("TM_" + met_id)
              if v_metabolite in cobra_model.metabolites:
                 continue
              #ref_met_id=[x for x in label_model.id_isotopomer_object_dict[iso_id]]
              for metabolite_id in metabolites:
                  r_metabolite = cobra_model.metabolites.get_by_id(metabolite_id)
                  the_reaction_id_list = [x.id for x in r_metabolite.reactions]
                  for the_reaction_id in the_reaction_id_list: 
                      the_reaction = cobra_model.reactions.get_by_id(the_reaction_id)
                      coefficient = (the_reaction.metabolites[r_metabolite])
                      the_reaction.add_metabolites({v_metabolite: coefficient})
                      sum_abs_source_reaction_bounds += maxupper_bound
              for the_reaction in v_metabolite.reactions:
                  if the_reaction.metabolites[v_metabolite]<0:
                     the_reaction.add_metabolites({v_metabolite:-2*the_reaction.metabolites[v_metabolite]})
              sink_reaction = Reaction("TMS_" +met_id)
              sink_reaction.add_metabolites({v_metabolite:-2})
              sink_reaction.lower_bound = the_min_flux
              sink_reaction.upper_bound = sum_abs_source_reaction_bounds / 2.
              turnover_metabolites.append(v_metabolite)
              sink_reactions.append(sink_reaction)
    print single_metabolites
    for metabolite_id in single_metabolites:
        print metabolite_id
        v_metabolite = Metabolite("TM_" + str(metabolite_id))
        if v_metabolite in cobra_model.metabolites:
                 continue
        # Now for reactions.  We include all reactions 
        # that create or consume the real metabolite.
        # These reactions therefore also drive creation of the
        # turnover metabolite, and we need to add a reaction
        # that is the sink.  By constraining this single
        # sink reaction, we ensure flux through the real reactions
        # to create the real metabolite.
        r_metabolite = cobra_model.metabolites.get_by_id(metabolite_id)
        sum_abs_source_reaction_bounds = 0
        
        the_reaction_id_list = [x.id for x in r_metabolite.reactions]
        for the_reaction_id in the_reaction_id_list:
            the_reaction = cobra_model.reactions.get_by_id(the_reaction_id)
            coefficient = abs(the_reaction.metabolites[r_metabolite])
            the_reaction.add_metabolites({v_metabolite:coefficient})
            #v_metabolite._reaction.add(the_reaction)
            # The model should be irreversible, and the upper bound
            # should be greater than the lower bound and positive.
            # Use 1000 for each reaction to be conservative.
            # E.g. might flip reactions back on later, don't
            # use the current upper bound for the reactions.
            #CFC Changed sum_abs_source_reaction_bounds += 1000.     
            sum_abs_source_reaction_bounds += maxupper_bound   
        # Add the sink reaction for the turnover metabolite
        # Since both creation and consumption of
        # the real metabolite drive the turnover,
        # we require 2 units of metabolite per
        # 1 unit of flux sink so this matches
        # the turnover through the real metabolite.
        sink_reaction = Reaction("TMS_" + metabolite_id)
        sink_reaction.add_metabolites({v_metabolite:-2})
        
        # Ensure a positive flux through the sink.
        # and maximum equal to maximal needed to
        # sustain reactions at their maximum.
        sink_reaction.lower_bound = the_min_flux
        sink_reaction.upper_bound = sum_abs_source_reaction_bounds / 2.
        v_metabolite._reaction.add(sink_reaction)        
        turnover_metabolites.append(v_metabolite)
        sink_reactions.append(sink_reaction)
    
    
    print sink_reactions
    #cobra_model.add_metabolites(turnover_metabolites)
    cobra_model.add_reactions(sink_reactions)
