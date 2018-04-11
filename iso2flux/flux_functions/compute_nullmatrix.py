from flux_variability_analysis import flux_variability_analysis
try:
  from cobra.core.arraybasedmodel import ArrayBasedModel
except:
  from cobra.core.ArrayBasedModel import ArrayBasedModel

import sympy
import numpy as np
import copy

from types import FunctionType
from sympy.simplify import simplify as _simplify, signsimp, nsimplify


def nullspace(matrix, simplify=False):
        """Returns list of vectors (Matrix objects) that span nullspace of self
        """
        from sympy.matrices import zeros
        
        simpfunc = simplify if isinstance(
            simplify, FunctionType) else _simplify
        reduced, pivots = matrix.rref(simplify=simpfunc)
        
        basis = []
        # create a set of vectors for the basis
        for i in range(matrix.cols - len(pivots)):
            basis.append(zeros(matrix.cols, 1))
        # contains the variable index to which the vector corresponds
        basiskey, cur = [-1]*len(basis), 0
        for i in range(matrix.cols):
            if i not in pivots:
                basiskey[cur] = i
                cur += 1
        for i in range(matrix.cols):
            if i not in pivots:  # free var, just set vector's ith place to 1
                basis[basiskey.index(i)][i, 0] = 1
            else:               # add negative of nonpivot entry to corr vector
                for j in range(i + 1, matrix.cols):
                    line = pivots.index(i)
                    v = reduced[line, j]
                    if simplify:
                        v = simpfunc(v)
                    if v:
                        if j in pivots:
                            # XXX: Is this the correct error?
                            raise NotImplementedError(
                                "Could not compute the nullspace of matrix self.")
                        basis[basiskey.index(j)][i, 0] = -v
        return [matrix._new(b) for b in basis]


def compute_nullmatrix(label_model,remove_inactive_reactions=True,fname=None):
       reaction_n_dict={}
       n_reaction_dict={}
       reaction0_list=[]
       model=copy.deepcopy(label_model.constrained_model)
       if remove_inactive_reactions:
         
         fva=flux_variability_analysis(model,fraction_of_optimum=0)
         reactions_to_remove=[]
         print "The following reactions will be omitted as they cannot carry flux:" 
         for x in fva: #Remove inactive reactions
           if fva[x]["minimum"]==0.0 and 0.0==fva[x]["maximum"]:
              print x,fva[x]
              reactions_to_remove.append(model.reactions.get_by_id(x))
              reaction0_list.append(x)
         for r in reactions_to_remove:
           r.remove_from_model()
         metabolites2remove=[]
         for x in model.metabolites: 
             if len(x.reactions)==0:
                 metabolites2remove.append(x)
         for x in metabolites2remove:
                 print x.id+" removed"
                 x.remove_from_model()
       rgroup_dict={}
       list_group_metabolites=[]
       n=0
       for reaction in  model.reactions:# sorted([x.id for x in model.reactions ]):
           #print str(n)+reaction_id
           #reaction=model.reactions.get_by_id(reaction_id)
           if "LABEL_RGROUP_" in reaction.id:
               group_dict={}
               group_metabolite=reaction.metabolites.keys()[0] #Group reactions only have one metabolite
               #Get the coefficients of reactions in the list
               for group_member_reaction in group_metabolite.reactions:
                   if reaction!=group_member_reaction:
                      group_dict[group_member_reaction.id]=group_member_reaction.metabolites[group_metabolite]
               list_group_metabolites.append(group_metabolite)
               rgroup_dict[reaction.id]=group_dict
           else:    
               reaction_n_dict[reaction.id]=n
               n_reaction_dict[n]=reaction.id
               n+=1
       print rgroup_dict 
       if fname==None:
         array_based=ArrayBasedModel(model)
         sm=array_based.S.toarray()
         print len(sm),len(sm[0])
         sm_list=[]
         for n_row,metabolite in enumerate(model.metabolites):#in enumerate(sorted([x.id for x in model.metabolites])):
           #metabolite=model.metabolites.get_by_id(metabolite_id)
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
         #print len(smm),len(smm[0])
         print "computing null matrix..."
         nullms=nullspace(smm,True)
         #nullms=smm.nullspace(True) #Simplify True
         nullmnp=np.transpose(np.array(nullms,dtype=np.float64))
         #original_nullmnp=copy.deepcopy(nullmnp)
         #Adding group reactions:
       else: 
         nullmnp=np.loadtxt(fname)
       for group_reaction in sorted(rgroup_dict.keys()):
           n=len(n_reaction_dict)
           reaction_n_dict[group_reaction]=n
           n_reaction_dict[n]=group_reaction
           if fname!=None: #If it has been loaded this is not necessary
              continue
           new_row=np.zeros((1,len(nullmnp[0])))
           for grouped_reaction in rgroup_dict[group_reaction]:
               print rgroup_dict[group_reaction]
               grouped_reaction_n=reaction_n_dict[grouped_reaction]
               grouped_reaction_coef=rgroup_dict[group_reaction][grouped_reaction]
               for n_coef,null_coef in enumerate(nullmnp[grouped_reaction_n]):
                   new_row[0][n_coef]+=null_coef*grouped_reaction_coef
           nullmnp = np.append(nullmnp,new_row,axis=0)
       for reaction0 in reaction0_list:
           n=len(n_reaction_dict)
           reaction_n_dict[reaction0]=n
           n_reaction_dict[n]=reaction0
           if fname!=None: #If it has been loaded this is not necessary
              continue
           new_row=np.zeros((1,len(nullmnp[0])))
           nullmnp = np.append(nullmnp,new_row,axis=0)
           #Create empty numpy array to append
       label_model.flux_solver_nullmnp=nullmnp    
       label_model.flux_solver_n_reaction_dict=n_reaction_dict
       label_model.flux_solver_reaction_n_dict=reaction_n_dict
       label_model.compute_fluxes=get_compute_fluxes_function(label_model)
       return label_model.flux_solver_nullmnp



def get_compute_fluxes_function(label_model):
 nullmnp=label_model.flux_solver_nullmnp
 n_reaction_dict=label_model.flux_solver_n_reaction_dict
 string="def compute_fluxes(nullmnp,ff): \n res=np.zeros(%s)\n" %(len(n_reaction_dict))
 for n,x in enumerate(n_reaction_dict):
    flux_expression=" res[%s]="%(n)
    abs_sum=sum([abs(x) for x in nullmnp[n]])
    if abs_sum==0:
       continue
    for c, coef in enumerate(nullmnp[n]):
        if coef==0:
           continue
        if coef.is_integer():
           if coef==1:
              flux_expression+="ff[%s]+" % (c) 
           elif coef==-1:
              flux_expression+="-ff[%s]+" % (c) 
           else:
              flux_expression+=str(coef)+"*ff[%s]+" % (c)
        else:
           value="{:.20f}".format(coef)
           flux_expression+=value+"*ff[%s]+" % (c)
           #flux_expression+="nullmnp[%s,%s]*ff[%s]+" % (n,c,c)
    flux_expression=flux_expression[:-1]+"\n" #Remove last plus sign
    string+=flux_expression
    #print flux_expression       
 
 string+=" return(res)"
 exec(string)
 #f=open("compute_fluxes.txt","w")
 #f.write(string)
 #f.close()
 return  compute_fluxes
