from cobra.flux_analysis.variability import flux_variability_analysis
from cobra.core.ArrayBasedModel import ArrayBasedModel
import sympy
import numpy as np

def compute_nullmatrix(label_model,remove_inactive_reactions=True):
       reaction_n_dict={}
       n_reaction_dict={}
       if remove_inactive_reactions:
         model=label_model.constrained_model
         fva=flux_variability_analysis(model,fraction_of_optimum=0)
         reactions_to_remove=[]
         print "The following reactions will be omitted as they cannot carry flux:" 
         for x in fva: #Remove inactive reactions
           if fva[x]["minimum"]==0.0 and 0.0==fva[x]["maximum"]:
              print x,fva[x]
              reactions_to_remove.append(model.reactions.get_by_id(x))
         for r in reactions_to_remove:
           r.remove_from_model()
       rgroup_dict={}
       list_group_metabolites=[]
       n=0
       for reaction in model.reactions:
           if "LABEL_RGROUP_" in reaction.id:
               group_dict={}
               group_metabolite=reaction.metabolites.keys()[0] #Group reactions only have one metabolite
               #Get the coefficients of reactions in the list
               for group_member_reaction in group_metabolite.reactions:
                   if reaction!=group_member_reaction:
                      group_dict[group_member_reaction]=group_member_reaction.metabolites[group_metabolite]
               list_group_metabolites.append(group_metabolite)
               rgroup_dict[reaction]=group_dict
           else:    
               reaction_n_dict[reaction.id]=n
               n_reaction_dict[n]=reaction.id
               n+=1
       print rgroup_dict 
       array_based=ArrayBasedModel(model)
       sm=array_based.S.toarray()
       sm_list=[]
       for n_row,metabolite in enumerate(model.metabolites):
           if metabolite in list_group_metabolites:
              continue
           for n_col,reaction in enumerate(model.reactions):
               if n_col not in  n_reaction_dict:
                  continue
               coef=sm[n_row,n_col]
               if coef.is_integer():
                  sm_list.append(int(coef))
               else:
                  sm_list.append(sympy.Rational(coef))
       
       smm=sympy.Matrix(len(model.metabolites)-len(list_group_metabolites),len(reaction_n_dict),sm_list)
       
       print "computing null matrix..."
       nullms=smm.nullspace(True) #Simplify True
       nullmnp=np.transpose(np.array(nullms,dtype=np.float64))
       #original_nullmnp=copy.deepcopy(nullmnp)
       #Adding group reactions:
       for group_reaction in rgroup_dict:
           n=len(n_reaction_dict)
           reaction_n_dict[group_reaction.id]=n
           n_reaction_dict[n]=group_reaction.id
           new_row=np.zeros((1,len(nullmnp[0])))
           for grouped_reaction in rgroup_dict[group_reaction]:
               grouped_reaction_n=reaction_n_dict[grouped_reaction.id]
               grouped_reaction_coef=rgroup_dict[group_reaction][grouped_reaction]
               for n_coef,null_coef in enumerate(nullmnp[grouped_reaction_n]):
                   new_row[0][n_coef]+=null_coef*grouped_reaction_coef
           nullmnp = np.append(nullmnp,new_row,axis=0)
           #Create empty numpy array to append
       label_model.flux_solver_nullmnp=nullmnp    
       label_model.flux_solver_n_reaction_dict=n_reaction_dict
       label_model.flux_solver_reaction_n_dict=reaction_n_dict
       return label_model.flux_solver_nullmnp
