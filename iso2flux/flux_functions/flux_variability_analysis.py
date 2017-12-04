import cobra
from cobra.flux_analysis.variability import flux_variability_analysis as cobra_flux_variability_analysis

#If cobra version is 0.6 or newer it converts the FVA results from pandas.DataFrame to dict

if cobra.__version__ >="0.6": 
   def flux_variability_analysis(model,fraction_of_optimum=0,tolerance_feasibility=1e-6,reaction_list=None):
       fva={}
       if reaction_list!=None:
          if isinstance(reaction_list[0], basestring):
             reaction_list=[model.reactions.get_by_id(x) for x in reaction_list] 
       pandas_fva=cobra_flux_variability_analysis(model,fraction_of_optimum=fraction_of_optimum,tolerance_feasibility=tolerance_feasibility,reaction_list=reaction_list) 
       for reaction in pandas_fva.index:
           fva[reaction]={"maximum":pandas_fva.loc[reaction]["maximum"],"minimum":pandas_fva.loc[reaction]["minimum"]}
           
       return fva
else:
     flux_variability_analysis=cobra_flux_variability_analysis   
