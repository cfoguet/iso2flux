import copy
from warnings import warn
from itertools import izip
#from gurobipy import Model, LinExpr, GRB, QuadExpr
from cobra.core.Solution import Solution
from six import string_types, iteritems
import re
from time import time
#from ilabel.label_propagation_functions.get_reactants_dict import get_reactants_dict 
from cobra.flux_analysis.variability import flux_variability_analysis
import re
import numpy as np
#from ..label_propagation_functions.get_reactants_dict import get_reactants_dict 
#From Gimme

def convert_to_irreversible_with_indicators(cobra_model,reaction_id_list,metabolite_list, mutually_exclusive_directionality_constraint = False,label_model=None):
    #Function modified from the work by : """Schmidt BJ1, Ebrahim A, Metz TO, Adkins JN, Palsson B, Hyduke DR. GIM3E: condition-specific models of cellular metabolism developed from metabolomics and expression data Bioinformatics. 2013 Nov 15;29(22):2900-8. doi: 10.1093/bioinformatics/btt493. Epub 2013 Aug 23."""
    """Will break all of the reversible reactions into two separate irreversible
     reactions with different directions.  This function call modified from
     a version in the core cobra to facilitate the MILP formulation and
     include gene_reaction_rules with the reverse reaction
   
     Arguments:
      cobra_model: A model object which will be modified in place.
      mutually_exclusive_directionality_constraint: Boolean.  If True, turnover 
       reactions are constructed to serve as MILP constraints to prevent loops.
      
     Returns:
      None, cobra_model is modified in place
    
    
    
    
     """
    reactions_to_add = []
    from cobra.core.Reaction import Reaction
    from cobra.core import Metabolite
    reactions_to_make_irreversible=[]
    for x in reaction_id_list:
        reactions_to_make_irreversible.append(cobra_model.reactions.get_by_id(x))
    """for x in lexs:
        reactions_to_make_irreversible.append(cobra_model.reactions.get_by_id(x))"""
    
    #If a label model object  is provided make sure all experimentally measured metabolites (or at least one of the metabolites in the pool) is produced
    full_metabolite_list=copy.copy(metabolite_list)
    print full_metabolite_list
    if label_model!=None:
       emus=[]
       for condition in label_model.experimental_dict:
           for emu in label_model.experimental_dict[condition]:
               if emu not in emus:
                  emus.append(emu)
       measured_metabolite_dict={}
       
       for emu in emus:
           iso_id=str(label_model.emu_dict[emu]["met_id"])
           #print label_model.id_isotopomer_object_dict
           #isotopomer_object=label_model.id_isotopomer_object_dict[iso_id]
           metabolites=label_model.isotopomer_id_metabolite_id_dict[iso_id]
           print [iso_id,label_model.isotopomer_id_metabolite_id_dict[iso_id]]
           if isinstance(metabolites,list):
              for metabolite in metabolites:
                  full_metabolite_list.append(metabolite)
           else:
              full_metabolite_list.append(metabolites)
    
    for metabolites in full_metabolite_list:
       print metabolites
       if not isinstance(metabolites,list):
          metabolites=[metabolites]
       for metabolite in metabolites:
          print metabolite
          the_metabolite=cobra_model.metabolites.get_by_id(metabolite)
          for x in the_metabolite.reactions:
             if x not in reactions_to_make_irreversible:
              reactions_to_make_irreversible.append(x)    
                  
    for reaction in reactions_to_make_irreversible:
        # Potential artifact because a reaction might run backwards naturally
        # and this would result in adding an empty reaction to the
        # model in addition to the reverse reaction.
        if reaction.lower_bound < 0:
            #reverse_reaction = Reaction(reaction.id + "_reverse")
            reverse_reaction = reaction.copy()
            reverse_reaction.id = reaction.id + "_reverse"
            reverse_reaction.lower_bound = max(0,-1*reaction.upper_bound)
            reverse_reaction.upper_bound = reaction.lower_bound * -1.
            reaction.lower_bound = 0
            if reaction.upper_bound<0:
               reaction.upper_bound=0
            # Make the directions aware of each other
            reaction.notes["reflection"] = reverse_reaction.id
            reverse_reaction.notes["reflection"] = reaction.id
            reaction_dict = {}
            current_metabolites = [x for x in reaction.metabolites]
            for the_metabolite in current_metabolites:
                reaction_dict[the_metabolite] = -2 * reaction.get_coefficient(the_metabolite.id)
            reverse_reaction.add_metabolites(reaction_dict)
            reactions_to_add.append(reverse_reaction)
            # Also: GPRs should already copy
            # reverse_reaction.gene_reaction_rule = reaction.gene_reaction_rule
            # reverse_reaction._genes = reaction._genes
            
            if mutually_exclusive_directionality_constraint:
                # A continuous reaction bounded by 0., 1.
                # Serves as a source for the indicator metabolites
                tmp_source = Reaction('IRRMILP_direction_constraint_source_for_%s_and_%s'
                                                           %(reaction.id,
                                                             reverse_reaction.id))
                tmp_source.upper_bound = 1.
                tmp_source.lower_bound = 0.
                # The reverse indicator reaction is
                # an integer-valued reaction bounded by 0,1
                # that activates flux to the reverse reaction
                # and deactivates the forward reaction only when it is
                # turned on to 1
                tmp_indicator = Reaction('IRRMILP_reverse_indicator_for_%s_and_%s'
                                                           %(reaction.id,
                                                             reverse_reaction.id))
                tmp_indicator.upper_bound = 1
                tmp_indicator.lower_bound = 0
                tmp_indicator.variable_kind = 'integer'                    
                flux_constraint_forward = Metabolite(id = 
                     'IRRMILP_direction_constraint_for_%s'%reaction.id)
                flux_constraint_reverse = Metabolite(id = 
                     'IRRMILP_direction_constraint_for_%s'%reverse_reaction.id)
                flux_constraint_reverse._constraint_sense = 'G'
                flux_constraint_reverse._bound = 0.
                
                tmp_source.add_metabolites({flux_constraint_forward: 1})
                
                tmp_indicator.add_metabolites({flux_constraint_forward: -1,
                                      flux_constraint_reverse: 1})
                if reaction.upper_bound != 0:
                        reaction.add_metabolites({flux_constraint_forward: -1./reaction.upper_bound})
                else:
                    # could put 1.01 X the tolerance here,
                    # This is arbitrary.  Use 0.001
                    # since 1000 is a typical upper bound
                    reaction.add_metabolites({flux_constraint_forward: -0.001})
                if reverse_reaction.upper_bound != 0:
                    reverse_reaction.add_metabolites({flux_constraint_reverse: -1./reverse_reaction.upper_bound})
                else:
                    reverse_reaction.add_metabolites({flux_constraint_reverse: -0.001})
                reactions_to_add.append(tmp_indicator)
                reactions_to_add.append(tmp_source)
    cobra_model.add_reactions(reactions_to_add)

