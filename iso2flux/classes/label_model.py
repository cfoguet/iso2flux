import copy
import cobra
from cobra import Model, Reaction, Metabolite
from ..flux_functions.remove_reflection_co2_reactions import remove_reflection_co2_reactions
import os
import sys

from ..emus_functions.build_emu_model import build_emu_model
from ..emus_functions.check_mass_balance import check_mass_balance
from ..emus_functions.remove_identical_reactions import remove_identical_reactions
from ..emus_functions.split_model import split_model
from ..emus_functions.emu_add_label_ouputs_inputs import emu_add_label_ouputs_inputs
from ..emus_functions.expand_emu_models import expand_emu_models
from ..emus_functions.set_initial_label import set_initial_label
from ..emus_functions.remove_impossible_emus import remove_impossible_emus as rm_impossible_emus 
from ..emus_functions.build_variables_dict import build_variables_dict
from ..emus_functions.build_emu_reaction_dicts import build_emu_reaction_dicts
from ..emus_functions.find_label_propagating_fluxes import find_label_propagating_fluxes
from ..emus_functions.write_emus_equations import write_emus_equations
from ..emus_functions.compile_c_code import compile_c_code
from ..emus_functions.set_equations_variables import set_equations_variables
from ..emus_functions.solver import solver #To do, move it to emu equations
from ..emus_functions.set_reflections import set_reflections

class Label_model:
      """Class used to store and process all information rellevant for building the emu model and simulating label propagation"""
      def __init__(self,model,split_co2_reactions=False,reactions_with_forced_turnover=[],lp_tolerance_feasibility=1e-8,parameter_precision=0.0001):
          """
          Initializes the label model class
          
          model: COBRApy model object
		reference metabolic model
          split_co2_reactions: bool, optional. 
		If set to TRUE the forward and reverse reactions where CO2 is consumed or produced will be considered as independent reactions when it comes to label propagations allowing to define a separate label propagation for the forward and reverse reactions. This allows to have reversible carboxylation/decarboxylation reactions without explicitlty simulating label in CO2 
          reactions_with_forced_turnover: list of strings,optional.
		 List of reactions where label propagation will be simulated forward and backward regardless of wether flux bounds allow negative reactions fluxes
          lp_tolerance_feasibility: float, optional. 
		Tolerance feasibility for the linear programing solvers 
          parameter_precision: float optional. 
		Default Parameter precision used on parameter fitting. 
          """
          self.minimum_sd=0.01
          self.parameter_precision=parameter_precision
          self.lp_tolerance_feasibility=lp_tolerance_feasibility
          self.met_id_isotopomer_dict={} # metabolite -> isotopomer object dict
          self.metabolite_isotopomers_dict={} 
          self.label_propagation_dict={}
          self.isotopomer_object_list=[]
          self.id_isotopomer_object_dict={}
          self.metabolite_id_isotopomer_id_dict={}
          self.isotopomer_id_metabolite_id_dict={}
          #EMU dicts
          self.reaction_merged_reactions_dict={}
          self.merged_reactions_reactions_dict={} 
          #models
          self.metabolic_model=None
          self.constrained_model=None #same as metabolic model but with ratios added
          self.irreversible_metabolic_model=None
          self.simplified_metabolic_model=Model('simplified_metabolic_model')
          self.emu_model=None #Initial emu model 
          #emu_model2=None #Emu model withouth redudant reactions
          self.size_model_dict={}
          self.size_expanded_model_dict={}
          #size_expanded_model2_dict={} #Expanded model_dict withouth the impossible reactions
          #Experimental dicts
          self.labelled_substrates=[]
          self.initial_label={}
          self.condition_initial_label_yy_dict={}
          self.rsm_list=[]
          #Flux dict
          self.ratio_dict={} #dictionary of the ratios 
          self.turnover_flux_dict={}#dictionary of the turnover fluxes
          self.flux_dict={}
          self.flux_list=[]         
          #Solver variables
          self.condition_size_yy0_dict={}
          self.condition_size_yy_dict={} 
          #Active condition
          self.active_condition=None
          #Variables used in dynamic simulations
          self.vol_dict={}
          self.con_n_dict={}
          self.np_con=None
          self.size_emu_c_ode_dict={}
          self.split_reactions_dict={} #unused but keept for compatibility
          self.parameter_dict={}
          self.metabolic_model=model
          self.eqn_dir="equations"
          self.label_groups_reactions_dict={} #Add to save/load function
          self.reactions_propagating_label=[] #Add to save/load function
          self.p_dict={}
          self.emu_dict0={}
          self.emu_dict={} #Dictionary of all emus in the model
          self.emu_size_dict={} #Ditctionady of all emus sorted by size
          self.reaction_emu_dict={} #Dictionary metabolic reactions -> emu reactions
          self.emu_reaction_dict={} #Dictionady emu reactions -> reactions 
          self.metabolite_emu_dict={} #Dictionary metabolite -> emu
          self.emu_metabolite_dict={} #Dictionadyr emu -> metabolites
          self.expanded_reaction_emu_dict={}
          self.emu_reaction_expanded_dict={}
          self.expanded_emu_dict={}
          self.size_variable_dict={}
          self.size_inverse_variable_dict={}
          self.input_n_dict={}
          self.reaction_n_dict={}
          self.n_reaction_dict={}
          self.size_emu_c_eqn_dict={} #Dict of the size-> equations systems
          self.force_balance=True 
          self.experimental_dict={}
          self.input_m0_list={}
          self.data_name_emu_dict={}
          #Create irreversible model
          self.reactions_with_forced_turnover=reactions_with_forced_turnover
          self.irreversible_metabolic_model=copy.deepcopy(model)
          for reaction_id in reactions_with_forced_turnover:
              reaction=self.irreversible_metabolic_model.reactions.get_by_id(reaction_id)
              reaction.lower_bound=-1000
              reaction.upper_bound=1000
          cobra.manipulation.convert_to_irreversible(self.irreversible_metabolic_model) #Reactions with a negative lower_bound will be split into 2
          if split_co2_reactions==True:
             self.reversible_co2_reactions=remove_reflection_co2_reactions(self.irreversible_metabolic_model)
          else:
             self.reversible_co2_reactions=[]
          #self.reactions_with_forced_turnover=reactions_with_forced_turnover
          """for reaction_id in reactions_with_forced_turnover:
              reaction=irreversible_metabolic_model.reactions.get_by_id(reaction_id)
              reaction.lower_bound=-1000
              reaction.upper_bound=1000
          cobra.manipulation.convert_to_irreversible(self.irreversible_metabolic_model) #Reactions with a negative lower_bound will be split into 2
          if split_co2_reactions==True:
             remove_reflection_co2_reactions(self.irreversible_metabolic_model)"""
          
      def add_initial_label(self,metabolite_id,label_patterns,condition="condition0",total_concentration=1.0):
          """
          adds the initial labeling patters and abundances for one metabolite. Metabolites not defined using this function will be assumed to be unlabeled. 
          
          metabolite_id: string. 
		Id of the metabolite in the cobra model
          label_patterns: list of lists. 
		Labeling pattern in the initial metabolites and its abundance. It is not necessary to define unlabelled metabolites. As an example, for Glcuose 50% [1,2-13C]-Glucose and and 1% [1-13C]-Glucose this field should be [[[1,1,0,0,0,0],0.5],[[1,0,0,0,0,0],0.01]]
          condition: string,optional. 
		Name of the condition (i.e "labelled glucose,condition 1").
          total_concentration: 
		float,optional. Currently unused, should be left at default value 
          """
          label_dict={}
          label_dict["met_id"]=metabolite_id
          label_dict["label_patterns"]=label_patterns
          label_dict["condition"]=condition
          label_dict["total_concentration"]=total_concentration
                    
          self.labelled_substrates.append(label_dict)
          
      def build(self,emu0_dict={},force_balance=True,recompile_c_code=True,remove_impossible_emus=True,isotopic_steady_state=True,excluded_outputs_inputs=[],default_turnover=None,turnover_upper_bound=None, turnover_exclude_EX=True,clear_intermediate_data=True):
          """
          Builds the emu label propagation model
          
          emu_dict0: dict,optional. 
		Definition of the initial emus, typycally those that are measured experimentally. Should follow the structure and nomenclature of emus. For instance for Lactate it could be {'emu_lac': {'met_id': 'lac_c','done': False, 'carbons': [1, 2, 3], 'size': 3}. It can be generated automatically using the "read_experimental_mid" function. If experiemntal data has been loaded using read_experimental_mid function this variable is not necessary 
          force_balance:bool,optional. 
		If set to true will force isotopologues fractions to add up to 1 by writing the m0 fraction as a function of the others fractions. Performance is greatly improved forcing balance as the number of variables is reduced. However, seting force balance to false helpful to check for errors or missing reaction in the label propagation model. 
          recompile_c_code: bool, optional. 
		If set to False the compilation setp will be skipped. 
          remove_impossible_emus: bool,optional. 
		If True isotopologues that cannot be generated with defined intitial label will be removed
          isotopic_steady_state: bool,optional. 
		Currently unused.
          excluded_outputs_inputs: list of strings, optional. 
		Reaction IDs found in this list won't be automatically added as inputs or ouputs in the model
          default_turnover: float,deprecated.
          turnover_upper_bound: float, optional
                Maximum turnover allowed (the max absolute difference between the forward and and reverse flux of a reaction). If left at None it will automatically assign the turnover as the maximum upper bound or absolute lower bound 
          turnover_exclude_EX: bool,optional. 
		If set to true turnovers will not be added to Exchange reactions (EX_).  
          clear_intermediate_data: bool,optional. 
		If true intemediary models and dicts generated when building the emu models will be deleted 
          """ 
          self.emu_dict0={}
          self.emu_dict={} #Dictionary of all emus in the model
          self.emu_size_dict={} #Ditctionady of all emus sorted by size
          self.reaction_emu_dict={} #Dictionary metabolic reactions -> emu reactions
          self.emu_reaction_dict={} #Dictionady emu reactions -> reactions 
          self.metabolite_emu_dict={} #Dictionary metabolite -> emu
          self.emu_metabolite_dict={} #Dictionadyr emu -> metabolites
          self.expanded_reaction_emu_dict={}
          self.emu_reaction_expanded_dict={}
          self.expanded_emu_dict={}
          self.size_variable_dict={}
          self.size_inverse_variable_dict={}
          self.input_n_dict={}
          self.reaction_n_dict={}
          self.n_reaction_dict={}
          self.size_emu_c_eqn_dict={} #Dict of the size-> equations systems
          self.force_balance=force_balance 
          build_emu_model(self,emu0_dict)  
          split_model(self)
          emu_add_label_ouputs_inputs(self,excluded_outputs_inputs)
          remove_identical_reactions(self) 
          expand_emu_models(self)
          if turnover_exclude_EX!=True and turnover_exclude_EX in ("true",1,"yes","True","Yes"):
             turnover_exclude_EX=True
          else: 
             turnover_exclude_EX=False
          for initial_label in self.labelled_substrates:
              l0=initial_label
              set_initial_label(l0["met_id"],self,l0["label_patterns"],l0["condition"],l0["total_concentration"])
          if remove_impossible_emus:
             rm_impossible_emus(self)
          #set_reflections(self.size_expanded_model_dict)     
          find_label_propagating_fluxes(self)
          build_variables_dict(self,force_balance=force_balance)
          build_emu_reaction_dicts(self)
          set_equations_variables(self,force_balance=force_balance)
          write_emus_equations(self,c_code=True,force_balance=force_balance)
          #return
          compile_c_code(self,recompile=recompile_c_code)
          #original_directory=os.getcwd()
          #print original_directory
          sys.path.insert(0, self.eqn_dir)
          #os.chdir(self.eqn_dir)
          #print os.getcwd()
          import get_equations
          """try:
            get_equations=reload(get_equations)
          except:
            pass"""
          get_equations.get_equations(self.size_emu_c_eqn_dict)
          #os.chdir(original_directory)  
          self.constrained_model=copy.deepcopy(self.metabolic_model)
          """if default_turnover==None:
             default_turnover=0"""   
          if turnover_upper_bound==None: 
             turnover_upper_bound=0
             for reaction_id in self.reaction_n_dict:
              if reaction_id in self.metabolic_model.reactions:
                  reaction=self.metabolic_model.reactions.get_by_id(reaction_id)
                  turnover_upper_bound=max(abs(reaction.upper_bound),abs(reaction.lower_bound),turnover_upper_bound)
          for reaction in self.reaction_emu_dict:
               if reaction+"_reverse" in self.reaction_emu_dict:#self.simplified_metabolic_model.reactions :
                  if reaction in self.merged_reactions_reactions_dict:
                     reaction=self.merged_reactions_reactions_dict[reaction][0]
                  print [self.reactions_with_forced_turnover]
                  if ("EX_" in reaction and turnover_exclude_EX) and (reaction not in self.reactions_with_forced_turnover):
                      self.turnover_flux_dict[reaction]={"v":0,"lb":0,"ub":0}
                  else:
                      self.turnover_flux_dict[reaction]={"v":turnover_upper_bound/2.0,"lb":0,"ub":turnover_upper_bound}
                      #Add turnover for inputs and outputs         
          if clear_intermediate_data:
             del self.irreversible_metabolic_model
             del self.size_model_dict
             del self.size_expanded_model_dict
             del self.label_propagation_dict
             del self.isotopomer_object_list
             del self.simplified_metabolic_model
             del self.emu_model


