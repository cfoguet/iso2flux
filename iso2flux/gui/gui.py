from Tkinter import *
import copy
import re
import json
import tkFileDialog
import math
import numpy
import os
import sys
import time
import random
#Ilabel imports
import cobra
from cobra.flux_analysis.variability import flux_variability_analysis
from ..flux_functions.apply_ratios import apply_ratios, remove_ratio #Functions to apply ratios of fluxes
from ..flux_functions.define_reaction_group import define_reaction_group
from cobra.solvers import solver_dict

from ..emus_functions.solver import solver #To do, move it to emu equations

from ..output_functions.model_to_excel import model_to_excel
from ..output_functions.showr import showr
from ..output_functions.write_fva import write_fva
from ..output_functions.print_emu_results import  print_emu_results
from ..output_functions.convert_model_to_spreadsheet import  convert_model_to_spreadsheet
from ..output_functions.export_constraints import export_constraints
from ..output_functions.export_label_results import export_label_results
from ..input_functions.read_constraints import read_flux_constraints 

#from ilabel.fitting.anneal import annealing
from ..fitting.get_objective_function import get_objective_function
from ..fitting.apply_parameters import apply_parameters
from ..fitting.save_parameters import save_parameters
from ..fitting.load_parameters import load_parameters
from ..fitting.load_parameters import load_parameters
from ..fitting.identify_free_parameters import identify_free_parameters
from ..fitting.sampling import sampling
from ..fitting.anneal import annealing
from ..fitting.clear_parameters import clear_parameters
from ..flux_functions.expression_analysis import integrate_omics_imat,integrate_omics_gim3e
from ..misc.save_load_iso2flux_model import save_iso2flux_model,load_iso2flux_model
from ..fitting.confidence_intervals import estimate_confidence_intervals

from ..flux_functions.minimal_flux import add_flux_limit_constraints

from ..doc.open_manual import open_manual
from ..classes.label_model import Label_model #Import the label_model class

from cobra.flux_analysis.parsimonious import optimize_minimal_flux

from ..label_propagation_functions.find_missing_reactions import find_missing_reactions
from ..input_functions.read_experimental_mid import read_experimental_mid
from ..input_functions.read_isotopomer_model import read_isotopomer_model
from ..input_functions.read_constraints import read_flux_constraints 
from ..input_functions.read_isotoflux_settings import read_isotoflux_settings
from ..input_functions.create_cobra_model_from_file import create_cobra_model_from_file
from ..misc.write_spreadsheet import write_spreadsheet
from ..misc.read_spreadsheets import read_spreadsheets
from ..misc.check_steady_state import check_steady_state
from ..misc.check_simulated_fractions import check_simulated_fractions
known_extensions=["xlsx","xlsm","xltx","xltm","csv","txt" ]

if __name__=="__main__":
   #from ilabel.fitting.anneal import annealing
   #from iso2flux.fitting.get_objective_function import get_objective_function
   from iso2flux.fitting.apply_parameters import apply_parameters
   from iso2flux.fitting.save_parameters import save_parameters
   from iso2flux.fitting.load_parameters import load_parameters
   from iso2flux.fitting.load_parameters import load_parameters
   from iso2flux.fitting.identify_free_parameters import identify_free_parameters
   from iso2flux.fitting.sampling import sampling
   from iso2flux.fitting.anneal import annealing
   from iso2flux.fitting.clear_parameters import clear_parameters
   from iso2flux.flux_functions.expression_analysis import integrate_omics_imat,integrate_omics_gim3e
   from iso2flux.misc.save_load_iso2flux_model import save_iso2flux_model,load_iso2flux_model
   from iso2flux.fitting.confidence_intervals import estimate_confidence_intervals
   
   from iso2flux.flux_functions.minimal_flux import add_flux_limit_constraints
   
   from iso2flux.doc.open_manual import open_manual

class GUI:
      label_window=False
      ratio_inverse_name_dict={}
      def add_menubar(self):
          menubar = Menu(self.root)
          filemenu = Menu(menubar, tearoff=0)
          filemenu.add_command(label="Reset", command=self.reset_model)
          filemenu.add_command(label="Load", command=self.load_model) #TODO Swicth to load model
          filemenu.add_command(label="Save", command=self.save_model) #TODO Swicth to save model
          #TODO Save constraints/Load constraints
          menubar.add_cascade(label="File", menu=filemenu)
          if self.label_model.emu_dict!={}: #Create the fitting menu only if a label model exists
             fittingmenu = Menu(menubar, tearoff=0)
             fittingmenu.add_command(label="Automatically add parameters", command=self.create_add_parameters_window)
             fittingmenu.add_command(label="View parameters", command=self.create_view_parameters_window)
             fittingmenu.add_command(label="Fit parameters", command=self.create_fitting_window)
             fittingmenu.add_command(label="Sampling", command=self.create_sampling_window)
             fittingmenu.add_command(label="Calculate confidence intervals", command=self.create_confidence_interval_window)
             menubar.add_cascade(label="Fitting", menu=fittingmenu)
             options = Menu(menubar, tearoff=0)
           
          add_constrain_menu = Menu(menubar, tearoff=0)
          add_constrain_menu.add_command(label="Load constraints", command=self.load_constraints)
          add_constrain_menu.add_command(label="Restrict max fluxes", command=self.create_limit_max_fluxes_window)
          #add_constrain_menu.add_command(label="Integrate gene expression")
          menubar.add_cascade(label="Add constraints", menu=add_constrain_menu)
          
          gene_expression_menu = Menu(menubar, tearoff=0)
          if "cplex" in solver_dict or "gurobi" in solver_dict: #Only add the optin if there are the righ solvers to run the program
              gene_expression_menu.add_command(label="iMAT", command=self.create_imat_window)
          gene_expression_menu.add_command(label="GIM3E", command=self.create_gim3e_window)
          add_constrain_menu.add_cascade(label="Integrate gene expression",menu=gene_expression_menu)
          
          export_menu = Menu(menubar, tearoff=0)           
          export_menu.add_command(label="Export Fluxes", command=self.export_fluxes)
          export_menu.add_command(label="Export Constraints", command=self.export_constraints)
          export_menu.add_command(label="Export model", command=self.export_model) 
          if self.label_model.emu_dict!={}: #Create the fitting menu only if a label model exists
             export_menu.add_command(label="Export label simulation results", command=self.export_label)
          menubar.add_cascade(label="Export", menu=export_menu)
          
          add_help_menu = Menu(menubar, tearoff=0)
          add_help_menu.add_command(label="Open Manual", command=open_manual)
          #add_help_menu.add_command(label="Integrate gene expression", command=self.create_integrate_gene_expression_window)
          menubar.add_cascade(label="Help", menu=add_help_menu)
          
          #menubar.add_cascade(label="Options", menu=options)
          self.root.config(menu=menubar)
      
      def export_constraints(self):
          export_constraints(self.label_model.constrained_model,fn=None,ratio_dict=self.label_model.ratio_dict)
      
      def load_constraints(self):
          #self.reset_all_bounds()
          clear_parameters(self.label_model,parameter_dict=None,parameter_list=[], clear_ratios=False,clear_turnover=False,clear_fluxes=True) 
          self.label_model.constrained_model,self.label_model.ratio_dict=read_flux_constraints(self.label_model.constrained_model,ratio_dict=self.label_model.ratio_dict,file_name=None,create_copies=True) 
          self.update_ratio_list()
      
      def reset_model(self):
          self.reset_all_bounds()
          self.reset_all_turnover()
          clear_parameters(self.label_model,parameter_dict=None,parameter_list=[], clear_ratios=False,clear_turnover=False,clear_fluxes=True) 
          for ratio_id in self.label_model.ratio_dict.keys():
              remove_ratio(self.label_model.constrained_model,ratio_id,self.label_model.ratio_dict)
          
      def add_search_bar(self):
          search_labelframe=LabelFrame(self.root,text="Search")
          search_labelframe.pack(side=TOP,expand=False)
          self.search_entry=Entry(search_labelframe)
          self.search_entry.pack(side=LEFT,expand=False)
          self.add_filters_bar(search_labelframe)
          
      def add_filters_bar(self,root):
          self.show_variable = IntVar()
          self.show_variable.set(0)
          filters_labelframe=Frame(root)
          filters_labelframe.pack(side=LEFT,expand=False)
          show_frame=Frame(filters_labelframe)
          show_frame.pack(side=TOP)
          if self.label_model.emu_dict!={}: #Create the only if a label model exists
             radiobutton1 = Radiobutton(show_frame, text="Show all reactions", variable=self.show_variable, value=0)
             radiobutton1.pack( side = LEFT)
             radiobutton2 = Radiobutton(show_frame, text="Show only label reactions", variable=self.show_variable, value=1)
             radiobutton2.pack( side = LEFT)
             radiobutton3 = Radiobutton(show_frame, text="Show only metabolic reactions", variable=self.show_variable, value=2)
             radiobutton3.pack( side = LEFT)
          
      
      
      
      def create_limit_max_fluxes_window(self):
          top=self.limit_max_fluxes_window= Toplevel() 
          fraction_of_flux_minimum_frame=Frame(top)
          fraction_of_flux_minimum_frame.pack(side=TOP)
          Label(fraction_of_flux_minimum_frame,text="Maximum flux/min flux").pack(side=LEFT)
          self.fraction_of_flux_minimum_entry=Entry(fraction_of_flux_minimum_frame)
          self.fraction_of_flux_minimum_entry.delete(0, END)
          self.fraction_of_flux_minimum_entry.insert(0, 1) 
          self.fraction_of_flux_minimum_entry.pack(side=LEFT)
          metabolomics_file_frame=Frame(top)
          metabolomics_file_frame.pack(side=TOP)
          Label(metabolomics_file_frame,text="metabolomics file").pack(side=LEFT)
          self.metabolomics_file_entry=Entry(metabolomics_file_frame)
          self.metabolomics_file_entry.delete(0, END)
          self.metabolomics_file_entry.insert(0, "") 
          self.metabolomics_file_entry.pack(side=LEFT)
          button2 = Button(metabolomics_file_frame, text="Open",command=self.get_metabolomics_file)
          button2.pack( side = LEFT)  
          button = Button(top, text="Add contsraints",command=self.limit_max_fluxes)
          button.pack( side = TOP)
          
          
      def get_gene_expression_file(self):
          loaded_file = tkFileDialog.askopenfile(title='Choose a file',filetypes=[("xlsx","*.xlsx"),("csv","*.csv"),('All files','*.*')]) 
          file_name=loaded_file.name    
          self.gene_expression_file_entry.delete(0, END)
          self.gene_expression_file_entry.insert(0, file_name)
          
      def get_metabolomics_file(self):
          loaded_file = tkFileDialog.askopenfile(title='Choose a file',filetypes=[("xlsx",".xlsx"),("csv",".csv")]) 
          file_name=loaded_file.name    
          self.metabolomics_file_entry.delete(0, END)
          self.metabolomics_file_entry.insert(0, file_name)  
          
                
      
      def create_imat_window(self):
          top=self.integrate_gene_expression_window= Toplevel() 
          top.title="iMAT"
          low_expression_threshold_frame=Frame(top)
          low_expression_threshold_frame.pack(side=TOP)
          Label(low_expression_threshold_frame,text="Low expression percentile").pack(side=LEFT)
          self.low_expression_threshold_entry=Entry(low_expression_threshold_frame)
          self.low_expression_threshold_entry.delete(0, END)
          self.low_expression_threshold_entry.insert(0, 25) 
          self.low_expression_threshold_entry.pack(side=LEFT)
          
          high_expression_threshold_frame=Frame(top)
          high_expression_threshold_frame.pack(side=TOP)
          Label(high_expression_threshold_frame,text="High expression percentile").pack(side=LEFT)
          self.high_expression_threshold_entry=Entry(high_expression_threshold_frame)
          self.high_expression_threshold_entry.delete(0, END)
          self.high_expression_threshold_entry.insert(0, 75) 
          self.high_expression_threshold_entry.pack(side=LEFT)
          
          hex_epsilon_frame=Frame(top)
          hex_epsilon_frame.pack(side=TOP)
          Label(hex_epsilon_frame,text="Epsilon").pack(side=LEFT)
          self.hex_epsilon_entry=Entry(hex_epsilon_frame)
          self.hex_epsilon_entry.delete(0, END)
          self.hex_epsilon_entry.insert(0, 1) 
          self.hex_epsilon_entry.pack(side=LEFT)
          """lex_epsilon_frame=Frame(top)
          lex_epsilon_frame.pack(side=TOP)
          Label(lex_epsilon_frame,text="lex epsilon").pack(side=LEFT)
          self.lex_epsilon_entry=Entry(lex_epsilon_frame)
          self.lex_epsilon_entry.delete(0, END)
          self.lex_epsilon_entry.insert(0, 1e-6) 
          self.lex_epsilon_entry.pack(side=LEFT)"""
          
          imat_fraction_optimum_frame=Frame(top)
          imat_fraction_optimum_frame.pack(side=TOP)
          Label(imat_fraction_optimum_frame,text="Fraction of optimal imat objective").pack(side=LEFT)
          self.imat_fraction_optimum_entry=Entry(imat_fraction_optimum_frame)
          self.imat_fraction_optimum_entry.delete(0, END)
          self.imat_fraction_optimum_entry.insert(0, 1) 
          self.imat_fraction_optimum_entry.pack(side=LEFT)
                    
          gene_expression_file_frame=Frame(top)
          gene_expression_file_frame.pack(side=TOP)
          Label(gene_expression_file_frame,text="Gene expression file").pack(side=LEFT)
          self.gene_expression_file_entry=Entry(gene_expression_file_frame)
          self.gene_expression_file_entry.delete(0, END)
          self.gene_expression_file_entry.insert(0, "") 
          self.gene_expression_file_entry.pack(side=LEFT)
          button2 = Button(gene_expression_file_frame, text="Open",command=self.get_gene_expression_file)
          button2.pack( side = LEFT)  
          
          metabolomics_file_frame=Frame(top)
          metabolomics_file_frame.pack(side=TOP)
          Label(metabolomics_file_frame,text="Metabolomics file").pack(side=LEFT)
          self.metabolomics_file_entry=Entry(metabolomics_file_frame)
          self.metabolomics_file_entry.delete(0, END)
          self.metabolomics_file_entry.insert(0, "") 
          self.metabolomics_file_entry.pack(side=LEFT)
          button2 = Button(metabolomics_file_frame, text="Open",command=self.get_metabolomics_file)
          button2.pack( side = LEFT)  
          
          button2 = Button(top, text="Constrain using gene expression",command=self.imat)
          button2.pack( side = TOP)
      
      def imat(self):
           low_expression_threshold=float(self.low_expression_threshold_entry.get())
           high_expression_threshold=float(self.high_expression_threshold_entry.get())
           hex_epsilon=float(self.hex_epsilon_entry.get())
           lex_epsilon=self.label_model.parameter_precision#float(self.lex_epsilon_entry.get())
           imat_fraction_optimum=float(self.imat_fraction_optimum_entry.get())
           excel_name=self.gene_expression_file_entry.get()
           metabolite_list_fname=self.metabolomics_file_entry.get()
           """if metabolite_list_fname=="":
              metabolite_list_fname=None
           if "G_" in self.label_model.metabolic_model.genes[0].id:
              prefix="G_"
           else:
              prefix=""
           if "_AT" in self.label_model.metabolic_model.genes[0].id:
              sufix="_AT"
           else:
              sufix=""
           """
           prefix=self.label_model.p_dict["gene_expression_gene_prefix"]
           sufix=self.label_model.p_dict["gene_expression_gene_sufix"]   
           hexs,lexs,objective,status, genefva= integrate_omics_imat(self.label_model.constrained_model,excel_name,fraction_of_optimum=self.fraction_of_optimum,low_expression_threshold=low_expression_threshold, high_expression_threshold=high_expression_threshold,percentile=True,gene_method="average",metabolite_list_fname=metabolite_list_fname, epsilon=hex_epsilon,lex_epsilon=lex_epsilon,imat_fraction_optimum=imat_fraction_optimum,label_model=self.label_model,add_as_constraints=True,solver=None,gene_prefix=prefix,gene_sufix=sufix)
           print hexs
           print lexs 
           print objective
           reaction_id=self.active_reaction_get_bounds
           reaction=self.label_model.constrained_model.reactions.get_by_id(reaction_id)
           self.boundaries_lb_spinbox.delete(0,"end")
           self.boundaries_lb_spinbox.insert(INSERT,str(reaction.lower_bound))
           self.boundaries_ub_spinbox.delete(0,"end")
           self.boundaries_ub_spinbox.insert(INSERT,str(reaction.upper_bound))
           #print self.label_model.constrained_model.optimize()
           self.fraction_of_optimum=0
           self.integrate_gene_expression_window.destroy()
           self.run_fva()
           
      def create_gim3e_window(self):
          top=self.integrate_gene_expression_window= Toplevel() 
          top.title="gim3e"
          low_expression_threshold_frame=Frame(top)
          low_expression_threshold_frame.pack(side=TOP)
          Label(low_expression_threshold_frame,text="Low expression percentile").pack(side=LEFT)
          self.low_expression_threshold_entry=Entry(low_expression_threshold_frame)
          self.low_expression_threshold_entry.delete(0, END)
          self.low_expression_threshold_entry.insert(0, 25) 
          self.low_expression_threshold_entry.pack(side=LEFT)
          
          gim3e_fraction_optimum_frame=Frame(top)
          gim3e_fraction_optimum_frame.pack(side=TOP)
          Label(gim3e_fraction_optimum_frame,text="Fraction of optimal gim3e objective").pack(side=LEFT)
          self.gim3e_fraction_optimum_entry=Entry(gim3e_fraction_optimum_frame)
          self.gim3e_fraction_optimum_entry.delete(0, END)
          self.gim3e_fraction_optimum_entry.insert(0, 1) 
          self.gim3e_fraction_optimum_entry.pack(side=LEFT)
          
          
          gene_expression_file_frame=Frame(top)
          gene_expression_file_frame.pack(side=TOP)
          Label(gene_expression_file_frame,text="Gene expression file").pack(side=LEFT)
          self.gene_expression_file_entry=Entry(gene_expression_file_frame)
          self.gene_expression_file_entry.delete(0, END)
          self.gene_expression_file_entry.insert(0, "") 
          self.gene_expression_file_entry.pack(side=LEFT)
          button2 = Button(gene_expression_file_frame, text="Open",command=self.get_gene_expression_file)
          button2.pack( side = LEFT)
          
          metabolomics_file_frame=Frame(top)
          metabolomics_file_frame.pack(side=TOP)
          Label(metabolomics_file_frame,text="Metabolomics file").pack(side=LEFT)
          self.metabolomics_file_entry=Entry(metabolomics_file_frame)
          self.metabolomics_file_entry.delete(0, END)
          self.metabolomics_file_entry.insert(0, "") 
          self.metabolomics_file_entry.pack(side=LEFT)
          button3 = Button(metabolomics_file_frame, text="Open",command=self.get_metabolomics_file)
          button3.pack( side = LEFT)   
          
          button4 = Button(top, text="Constrain using gene expression",command=self. gim3e)
          button4.pack( side = TOP)  
           
          
      def gim3e(self):
           low_expression_threshold=float(self.low_expression_threshold_entry.get())
           gim3e_fraction_optimum=float(self.gim3e_fraction_optimum_entry.get())
           excel_name=self.gene_expression_file_entry.get()
           metabolite_list_fname=self.metabolomics_file_entry.get()
           """if metabolite_list_fname=="":
              metabolite_list_fname=None
           if "G_" in self.label_model.metabolic_model.genes[0].id:
              prefix="G_"
           else:
              prefix=""
           if "_AT" in self.label_model.metabolic_model.genes[0].id:
              sufix="_AT"
           else:
              sufix=""  
           """
           prefix=self.label_model.p_dict["gene_expression_gene_prefix"]
           sufix=self.label_model.p_dict["gene_expression_gene_sufix"]  
           penalty_dict,objective,fva=integrate_omics_gim3e(self.label_model.constrained_model,excel_name,fraction_of_optimum=self.fraction_of_optimum,low_expression_threshold=low_expression_threshold,absent_gene_expression=50,percentile=True,gene_method="average",metabolite_list_fname=metabolite_list_fname,label_model=self.label_model,epsilon=0.0001,gim3e_fraction_optimum=gim3e_fraction_optimum,add_as_constraints=True,boundaries_precision=0.01,gene_prefix=prefix,gene_sufix=sufix)
           print penalty_dict
           print objective
           print fva
           reaction_id=self.active_reaction_get_bounds
           reaction=self.label_model.constrained_model.reactions.get_by_id(reaction_id)
           self.boundaries_lb_spinbox.delete(0,"end")
           self.boundaries_lb_spinbox.insert(INSERT,str(reaction.lower_bound))
           self.boundaries_ub_spinbox.delete(0,"end")
           self.boundaries_ub_spinbox.insert(INSERT,str(reaction.upper_bound))
           #print self.label_model.constrained_model.optimize()
           self.integrate_gene_expression_window.destroy()
           self.fraction_of_optimum=0
           self.run_fva()
      def limit_max_fluxes(self):
           print "constraining fluxes"
           fraction_of_flux_minimum=max(float(self.fraction_of_flux_minimum_entry.get()),1)
           metabolomics_file=self.metabolomics_file_entry.get()
           
           if "cplex" in solver_dict:
               solver="cplex"
           elif "gurobi" in solver_dict:
               solver="gurobi"
           else:
               solver=None #Cobra will select the default solver
           add_flux_limit_constraints(self.label_model.constrained_model,fraction_of_optimum_objective=self.fraction_of_optimum, fraction_of_flux_minimum=fraction_of_flux_minimum,solver=solver,label_model=self.label_model,metabolite_list_file_name=metabolomics_file)
           self.limit_max_fluxes_window.destroy()
           reaction_id=self.active_reaction_get_bounds
           reaction=self.label_model.constrained_model.reactions.get_by_id(reaction_id)
           self.boundaries_lb_spinbox.delete(0,"end")
           self.boundaries_lb_spinbox.insert(INSERT,str(reaction.lower_bound))
           self.boundaries_ub_spinbox.delete(0,"end")
           self.boundaries_ub_spinbox.insert(INSERT,str(reaction.upper_bound))
           #print self.label_model.constrained_model.optimize()
           self.run_fva()
                    
      def create_sampling_window(self):
          top=self.sampling_window= Toplevel()
          nframe=Frame(top)
          nframe.pack(side=TOP)
          Label(nframe,text="Enter number of samples:").pack(side=LEFT)
          self.sampling_n_entry=Entry(nframe)
          self.sampling_n_entry.delete(0, END)
          self.sampling_n_entry.insert(0, 100) 
          self.sampling_n_entry.pack(side=LEFT)
          sampling_max_turnoverframe=Frame(top)
          sampling_max_turnoverframe.pack(side=TOP)
          """Label(sampling_max_turnoverframe,text="max turnover").pack(side=LEFT)
          self.sampling_max_turnover_entry=Entry(sampling_max_turnoverframe)
          self.sampling_max_turnover_entry.delete(0, END)
          self.sampling_max_turnover_entry.insert(0, 1000) 
          self.sampling_max_turnover_entry.pack(side=LEFT)"""
          
          """sampling_change_threshold=Frame(top)
          sampling_change_threshold.pack(side=TOP)
          Label(sampling_change_threshold,text="change_threshold").pack(side=LEFT)
          self.sampling_change_threshold_entry=Entry(sampling_change_threshold)
          self.sampling_change_threshold_entry.delete(0, END)
          self.sampling_change_threshold_entry.insert(0, 0.1) 
          self.sampling_change_threshold_entry.pack(side=LEFT)"""  
     
          button01 = Button(top, text="Start",command=self.sampling)
          button01.pack( side = TOP)          
          
          
      
      def sampling(self):
          if self.label_window==False:
             self.create_figure_window()
          self.root.update_idletasks()
          n=int(self.sampling_n_entry.get())
          #max_turnover=float(self.sampling_max_turnover_entry.get())
          self.sampling_window.destroy()
          #sampling(self.label_model,n=n,fraction_of_optimum=self.fraction_of_optimum,output_emu_list=None,max_turnover=50,fba_mode="fba",parameter_precision= sel,gui=self,change_threshold=0.1)
          original_parameters=copy.deepcopy(self.label_model.parameter_dict)
          sampling(self.label_model,n=n,fraction_of_optimum=self.fraction_of_optimum,output_emu_list=None,max_turnover=1000,fba_mode=self.fba_mode.get(),parameter_precision= self.change_threshold/100.0,gui=self,change_threshold=self.change_threshold,parameter_dict=self.label_model.parameter_dict)
          #print ["bye",self.label_model.constrained_model.reactions.get_by_id("biomass").lower_bound]
          #self.label_model.parameter_dict=original_parameters
          apply_parameters(self.label_model,original_parameters)
      """def fba_options(self):
          top=Toplevel()
          top.title="Select FBA mode"
          Label(top,text="Select FBA mode").pack(side=TOP)  
          fba_mode_radio_button0 = Radiobutton(top, text="FBA", variable=self.fba_mode, value="FBA")
          fba_mode_radio_button0.pack( side = TOP)
          fba_mode_radio_button1 = Radiobutton(top, text="pFBA", variable=self.fba_mode, value="pFBA")
          fba_mode_radio_button1.pack( side = TOP) 
          Button(top, text="Apply",command=self.run_fva).pack( side = TOP)
          Button(top, text="Close",command=top.destroy).pack( side = TOP)"""  
      def create_add_parameters_window(self):
          self.add_parameters_window=top= Toplevel()
          top.title="Automatically add parameters"
          self.add_fluxes_flag = IntVar()
          self.add_turnover_flag=IntVar()
          self.add_fluxes_flag.set(1)
          self.add_turnover_flag.set(1)
          text=Label(top,text="Automatically add:")
          text.pack(side = TOP)
          C1 = Checkbutton(top, text = "Free fluxes as parameters", variable = self.add_fluxes_flag,onvalue = 1, offvalue = 0)
          C2 = Checkbutton(top, text = "Reversible turnover fluxes as parameters", variable = self.add_turnover_flag,onvalue = 1, offvalue = 0)
          C1.pack()
          C2.pack()
          abutton=Button(top, text="Accept",command=self.add_parameters)
          abutton.pack( side = TOP)
          cbutton=Button(top, text="Cancel",command=top.destroy)
          cbutton.pack( side = TOP)
      
      def add_parameters(self):
          self.add_parameters_window.destroy() 
          if self.add_fluxes_flag.get()==1:
             self.add_all_parameters()             
          if self.add_turnover_flag.get()==1:
             self.add_all_turnover_as_parameter()
          if self.waiting_for_parameters==True:
             self.create_fitting_window()
          
      def create_fitting_window(self):
          if len(self.label_model.parameter_dict)==0:
                 self.waiting_for_parameters=True
                 self.create_add_parameters_window()
                 return
          top=self.fitting_window= Toplevel()
          top.title="Parameter fitting"
          nframe=Frame(top)
          nframe.pack(side=TOP)
          Label(nframe,text="n cycles").pack(side=LEFT)
          self.nframeentry=Entry(nframe)
          self.nframeentry.delete(0, END)
          self.nframeentry.insert(0, self.annealing_n) 
          self.nframeentry.pack(side=LEFT)
          mframe=Frame(top)
          mframe.pack(side=TOP)
          Label(mframe,text="simulations per cycle").pack(side=LEFT)
          self.mframeentry=Entry(mframe)
          self.mframeentry.delete(0, END)
          self.mframeentry.insert(0, self.annealing_m) 
          self.mframeentry.pack(side=LEFT)
          p0frame=Frame(top) 
          p0frame.pack(side=TOP)
          Label(p0frame,text="p0").pack(side=LEFT)
          self.p0frameentry=Entry(p0frame)
          self.p0frameentry.delete(0, END)
          self.p0frameentry.insert(0, self.annealing_p0) 
          self.p0frameentry.pack(side=LEFT)
          pfframe=Frame(top) 
          pfframe.pack(side=TOP)
          Label(pfframe,text="pf").pack(side=LEFT)
          self.pfframeentry=Entry(pfframe)
          self.pfframeentry.delete(0, END)
          self.pfframeentry.insert(0, self.annealing_pf) 
          self.pfframeentry.pack(side=LEFT)
          
          """max_dframe=Frame(top) 
          max_dframe.pack(side=TOP)
          Label(max_dframe,text="max relative perturbation").pack(side=LEFT)
          self.max_dframeentry=Entry(max_dframe)
          self.max_dframeentry.delete(0, END)
          self.max_dframeentry.insert(0, 0.1) 
          self.max_dframeentry.pack(side=LEFT)"""
          
          max_perturbationframe=Frame(top) 
          max_perturbationframe.pack(side=TOP)
          """Label(max_perturbationframe,text="max parameter perturbation").pack(side=LEFT)
          self.max_perturbationframeentry=Entry(max_perturbationframe)
          self.max_perturbationframeentry.delete(0, END)
          self.max_perturbationframeentry.insert(0, 1) 
          self.max_perturbationframeentry.pack(side=LEFT)"""
          
          n_process_frame=Frame(top) 
          n_process_frame.pack(side=TOP)
          Label(n_process_frame,text="number of paralel processes").pack(side=LEFT)
          self.n_process_entry=Entry(n_process_frame)
          self.n_process_entry.delete(0, END)
          self.n_process_entry.insert(0, self.annealing_n_processes) 
          self.n_process_entry.pack(side=LEFT)
          
          annealing_n_iterations=Frame(top) 
          annealing_n_iterations.pack(side=TOP)
          Label(annealing_n_iterations,text="Number of iterations of the annealing algorythm").pack(side=LEFT)
          self.annealing_n_iterations_entry=Entry(annealing_n_iterations)
          self.annealing_n_iterations_entry.delete(0, END)
          self.annealing_n_iterations_entry.insert(0, self.annealing_n_iterations) 
          self.annealing_n_iterations_entry.pack(side=LEFT)
          
          Label(top,text="n parameters is "+str(len(self.label_model.parameter_dict))).pack(side=TOP)
          
          max_sampleframe=Frame(top) 
          max_sampleframe.pack(side=TOP)
          Label(max_sampleframe,text="max size of parameter sample").pack(side=LEFT)
          self.max_sampleframeentry=Entry(max_sampleframe)
          self.max_sampleframeentry.delete(0, END)
          self.max_sampleframeentry.insert(0, int(len(self.label_model.parameter_dict)*self.relative_max_random_sample)) 
          self.max_sampleframeentry.pack(side=LEFT)
          
          min_sampleframe=Frame(top) 
          min_sampleframe.pack(side=TOP)
          Label(min_sampleframe,text="min size of parameter sample").pack(side=LEFT)
          self.min_sampleframeentry=Entry(min_sampleframe)
          self.min_sampleframeentry.delete(0, END)
          self.min_sampleframeentry.insert(0, int(len(self.label_model.parameter_dict)*self.relative_min_random_sample)) 
          self.min_sampleframeentry.pack(side=LEFT)
          
          
          
          
          button01 = Button(top, text="Start",command=self.fitting)
          button01.pack( side = TOP)
          
      def create_confidence_interval_window(self):
          top=self.confidence_window= Toplevel()
          top.title="Estimate confidence interval"
          
          confidence_parameters_frame=LabelFrame(top,text="Parameters:")
          confidence_parameters_frame.pack(side=TOP)
          
          signficance_level_frame=Frame(confidence_parameters_frame)
          signficance_level_frame.pack(side=TOP)
          Label(signficance_level_frame,text="Significance Level").pack(side=LEFT)
          self.signficance_level_frame_entry=Entry(signficance_level_frame)
          self.signficance_level_frame_entry.delete(0, END)
          self.signficance_level_frame_entry.insert(0, self.label_model.p_dict['confidence_significance']*100) 
          self.signficance_level_frame_entry.pack(side=LEFT)
          
          setp_size_frame=Frame(confidence_parameters_frame)
          setp_size_frame.pack(side=TOP)
          Label(setp_size_frame,text="Step size (%)").pack(side=LEFT)
          self.setp_size_entry=Entry(setp_size_frame)
          self.setp_size_entry.delete(0, END)
          self.setp_size_entry.insert(0, 10) 
          self.setp_size_entry.pack(side=LEFT)
          
          
          annealing_parameters_frame=LabelFrame(top,text="Annealing specific parameters:")
          annealing_parameters_frame.pack(side=TOP)
          nframe=Frame(annealing_parameters_frame)
          nframe.pack(side=TOP)
          Label(nframe,text="n cycles").pack(side=LEFT)
          self.nframe_confidence_entry=Entry(nframe)
          self.nframe_confidence_entry.delete(0, END)
          self.nframe_confidence_entry.insert(0, self.annealing_n) 
          self.nframe_confidence_entry.pack(side=LEFT)
          mframe=Frame(annealing_parameters_frame)
          mframe.pack(side=TOP)
          Label(mframe,text="simulations per cycle").pack(side=LEFT)
          self.mframe_confidence_entry=Entry(mframe)
          self.mframe_confidence_entry.delete(0, END)
          self.mframe_confidence_entry.insert(0, self.annealing_m) 
          self.mframe_confidence_entry.pack(side=LEFT)
          p0frame=Frame(annealing_parameters_frame) 
          p0frame.pack(side=TOP)
          Label(p0frame,text="p0").pack(side=LEFT)
          self.p0frame_confidence_entry=Entry(p0frame)
          self.p0frame_confidence_entry.delete(0, END)
          self.p0frame_confidence_entry.insert(0, self.annealing_p0) 
          self.p0frame_confidence_entry.pack(side=LEFT)
          pfframe=Frame(annealing_parameters_frame) 
          pfframe.pack(side=TOP)
          Label(pfframe,text="pf").pack(side=LEFT)
          self.pfframe_confidence_entry=Entry(pfframe)
          self.pfframe_confidence_entry.delete(0, END)
          self.pfframe_confidence_entry.insert(0, self.annealing_pf) 
          self.pfframe_confidence_entry.pack(side=LEFT)
          
          """max_dframe=Frame(top) 
          max_dframe.pack(side=TOP)
          Label(max_dframe,text="max relative perturbation").pack(side=LEFT)
          self.max_dframeentry=Entry(max_dframe)
          self.max_dframeentry.delete(0, END)
          self.max_dframeentry.insert(0, 0.1) 
          self.max_dframeentry.pack(side=LEFT)"""
          
          #max_perturbationframe=Frame(top) 
          #max_perturbationframe.pack(side=TOP)
          """Label(max_perturbationframe,text="max parameter perturbation").pack(side=LEFT)
          self.max_perturbationframeentry=Entry(max_perturbationframe)
          self.max_perturbationframeentry.delete(0, END)
          self.max_perturbationframeentry.insert(0, 1) 
          self.max_perturbationframeentry.pack(side=LEFT)"""
          
          n_process_frame=Frame(annealing_parameters_frame) 
          n_process_frame.pack(side=TOP)
          Label(n_process_frame,text="number of paralel processes").pack(side=LEFT)
          self.n_process_confidence_entry=Entry(n_process_frame)
          self.n_process_confidence_entry.delete(0, END)
          self.n_process_confidence_entry.insert(0, self.annealing_n_processes) 
          self.n_process_confidence_entry.pack(side=LEFT)
          
          annealing_n_iterations=Frame(annealing_parameters_frame) 
          annealing_n_iterations.pack(side=TOP)
          Label(annealing_n_iterations,text="Number of iterations of the annealing algorythm").pack(side=LEFT)
          self.annealing_n_iterations_entry=Entry(annealing_n_iterations)
          self.annealing_n_iterations_entry.delete(0, END)
          self.annealing_n_iterations_entry.insert(0, self.annealing_n_iterations) 
          self.annealing_n_iterations_entry.pack(side=LEFT)
          
          
          Label(annealing_parameters_frame,text="n parameters is "+str(len(self.label_model.parameter_dict))).pack(side=TOP)
          
          max_sampleframe=Frame(annealing_parameters_frame) 
          max_sampleframe.pack(side=TOP)
          Label(max_sampleframe,text="max size of parameter sample").pack(side=LEFT)
          self.max_sampleframe_confidence_entry=Entry(max_sampleframe)
          self.max_sampleframe_confidence_entry.delete(0, END)
          self.max_sampleframe_confidence_entry.insert(0, int(len(self.label_model.parameter_dict)*self.relative_max_random_sample)) 
          self.max_sampleframe_confidence_entry.pack(side=LEFT)
          
          min_sampleframe=Frame(annealing_parameters_frame) 
          min_sampleframe.pack(side=TOP)
          Label(min_sampleframe,text="min size of parameter sample").pack(side=LEFT)
          self.min_sampleframe_confidence_entry=Entry(min_sampleframe)
          self.min_sampleframe_confidence_entry.delete(0, END)
          self.min_sampleframe_confidence_entry.insert(0, int(len(self.label_model.parameter_dict)*self.relative_min_random_sample)) 
          self.min_sampleframe_confidence_entry.pack(side=LEFT)
          
          self.confidence_sbml = IntVar()
          self.confidence_constraint = IntVar()
          check_confidence_sbml=Checkbutton(top, text = "Save as Constrained Model", variable = self.confidence_sbml,onvalue = 1, offvalue = 0)
          check_confidence_constraint=Checkbutton(top, text = "Add results as constraints to the active model", variable = self.confidence_constraint,onvalue = 1, offvalue = 0)
          evaluate_all_button = Button(top, text="Evaluate all fluxes",command=self.confidence_start_warning)
          evaluate_specific_button = Button(top, text="Evaluate specific fluxes & ratios",command=self.confidence_select_fluxes)
          check_confidence_sbml.pack(side=TOP)
          check_confidence_constraint.pack(side=TOP)
          evaluate_all_button.pack(side=TOP)
          evaluate_specific_button.pack(side=TOP) 
          """select_parameters_frame=LabelFrame(top,text="Select Parameters")    
          select_parameters_frame.pack(side=TOP)
          
          all_parameter_frame=LabelFrame(select_parameters_frame,text="Parameters")   
          all_parameter_frame.pack(side=LEFT,fill=Y,expand="yes") 
          
          confidence_parameter_selector_listbox_frame=Frame(all_parameter_frame)   
          confidence_parameter_selector_listbox_frame.pack(side=TOP) 
          parameter_selector_scrollbar = Scrollbar(confidence_parameter_selector_listbox_frame)
          parameter_selector_scrollbar.pack( side = RIGHT, fill=Y,expand="yes")
          self.confidence_parameter_selector_listbox=confidence_parameter_selector_listbox = Listbox(confidence_parameter_selector_listbox_frame, yscrollcommand =parameter_selector_scrollbar.set , height = 3, width=14)
          parameter_selector_scrollbar['command'] = confidence_parameter_selector_listbox.yview
          confidence_parameter_selector_listbox['yscrollcommand'] = parameter_selector_scrollbar.set
          confidence_parameter_selector_listbox.pack( side = RIGHT, fill = BOTH ,expand="yes")
          parameter_selector_scrollbar.config( command = confidence_parameter_selector_listbox.yview )
          for parameter in sorted(self.label_model.parameter_dict,key=lambda v: v.upper()):
                 confidence_parameter_selector_listbox.insert(END,parameter)
          remove_parameter_button = Button(all_parameter_frame, text="Add selected",command=self.confidence_add_parameter)
          remove_parameter_button.pack( side = TOP)
          remove_parameter_button = Button(all_parameter_frame, text="Add all",command=self.confidence_add_all_parameters)
          remove_parameter_button.pack( side = TOP)
          
                
          added_parameter_frame=LabelFrame(select_parameters_frame,text="Selected parameters")    
          added_parameter_frame.pack(side=LEFT,fill=Y,expand="yes")
          
          added_parameter_listbox_frame=Frame(added_parameter_frame)   
          added_parameter_listbox_frame.pack(side=TOP) 
           
          added_parameter_selector_scrollbar = Scrollbar(added_parameter_listbox_frame)
          added_parameter_selector_scrollbar.pack( side = RIGHT, fill=Y,expand="yes")
          self.confidence_added_parameter_selector_listbox=confidence_added_parameter_selector_listbox = Listbox(added_parameter_listbox_frame, yscrollcommand = parameter_selector_scrollbar.set , height = 3, width=14)
          added_parameter_selector_scrollbar['command'] = confidence_added_parameter_selector_listbox.yview
          confidence_added_parameter_selector_listbox['yscrollcommand'] = added_parameter_selector_scrollbar.set
          confidence_added_parameter_selector_listbox.pack( side = RIGHT, fill = BOTH ,expand="yes")
          added_parameter_selector_scrollbar.config( command = confidence_added_parameter_selector_listbox.yview )
          remove_parameter_button = Button(added_parameter_frame, text="Remove selected",command=self.confidence_remove_selected)
          remove_parameter_button.pack( side = TOP)"""
          
          self.confidence_selected_parameters=[]
          
          #button01 = Button(top, text="Start",command=self.confidence_start_warning)
          #button01.pack( side = TOP)
          
      """def confidence_add_parameter(self):
          selection=self.confidence_parameter_selector_listbox.curselection()
          if selection!=():
             parameter=self.confidence_parameter_selector_listbox.get(selection[0])
             if parameter not in self.confidence_selected_parameters:
                self.confidence_selected_parameters.append(parameter)
                self.confidence_added_parameter_selector_listbox.insert(END,parameter)
      
      def confidence_add_all_parameters(self):
          for parameter in self.label_model.parameter_dict:
              if parameter not in self.confidence_selected_parameters:
                self.confidence_selected_parameters.append(parameter)
                self.confidence_added_parameter_selector_listbox.insert(END,parameter)
      
      def confidence_remove_selected(self):
          selection=self.confidence_added_parameter_selector_listbox.curselection()
          if selection!=():
             parameter=self.confidence_added_parameter_selector_listbox.get(selection[0])
             self.confidence_selected_parameters.remove(parameter)
             self.confidence_added_parameter_selector_listbox.delete(selection)"""
      
      def confidence_select_fluxes(self):
          data_rows_dict=read_spreadsheets(file_names=None,csv_delimiter=',',more_than_1=True,tkinter_title="Chose the file with reaction ID or ratios that should be evaluated")
          for data in data_rows_dict:
              for row in data_rows_dict[data]:
                  for element in row:
                      if element==None:
                         continue
                      if element in self.confidence_selected_parameters:
                         continue 
                      if element in self.label_model.constrained_model.reactions:
                         self.confidence_selected_parameters.append(element)
                      elif "/" in element:
                          element_list=element.split("/")
                          if element_list[0] in self.label_model.constrained_model.reactions and element_list[1] in self.label_model.constrained_model.reactions:
                             self.confidence_selected_parameters.append(element)
          self.confidence_start_warning()  
      
          
      def confidence_start_warning(self):
          self.confidence_warning_window=top=Toplevel()      
          self.confidence_warning_window.title("Warning")
          if self.fitted_parameters==True:
             Label(top,text="Calcuating confidence intervals\nmay take a long time.\nAre your sure you want to proceed?").pack(side=TOP)  
          else:
             Label(top,text='Calcuating confidence intervals requires having excetued\nthe "Fit Parameters" function beforehand to be signficartive\n"Fit Parameters has not been exectued in this session"\nAre your sure you want to proceed?').pack(side=TOP)
          yes_button= Button(top, text="Yes",command=self.calculate_confidence_intervals)
          no_button= Button(top, text="No",command=top.destroy)
          yes_button.pack(side=TOP)
          no_button.pack(side=TOP)
      
      def calculate_confidence_intervals(self):
          parameter_list=self.confidence_selected_parameters
          signficance=float(self.signficance_level_frame_entry.get())/100.0
          n_perturbations=200
          setp_size=float(self.setp_size_entry.get())/100.0
          minimum_perturbation=self.change_threshold
          n_parameters=len(self.label_model.parameter_dict)
          self.relative_max_random_sample=relative_max_random_sample=float(self.max_sampleframe_confidence_entry.get())/n_parameters
          self.relative_min_random_sample=relative_min_random_sample=float(self.min_sampleframe_confidence_entry.get())/n_parameters
          self.annealing_n=annealing_n=int(self.nframe_confidence_entry.get())
          self.annealing_m=annealing_m=int(self.mframe_confidence_entry.get())
          self.annealing_p0=annealing_p0=float(self.p0frame_confidence_entry.get())
          self.annealing_pf=annealing_pf=float(self.pfframe_confidence_entry.get())
          self.annealing_n_iterations=int(self.annealing_n_iterations_entry.get())
          output=True
          self.annealing_n_processes=annealing_n_processes=int(self.n_process_confidence_entry.get())
          annealing_cycle_time_limit=self.label_model.p_dict['annealing_cycle_time_limit']#(annealing_m*50)
          annealing_cycle_max_attempts=5
          self.create_waiting_window() 
          self.confidence_warning_window.destroy()
          self.root.update_idletasks()
          fname = tkFileDialog.asksaveasfilename(parent=self.root,title="Save results as...",filetypes=[("xlsx","*.xlsx"),("csv","*.csv")])
          if not any(x in fname.lower() for x in ["xlsx","xlsm","xltx","xltm","csv","txt" ]): 
             fname+=".xlsx" 
          if self.confidence_sbml.get()==1:
                sbml_name = tkFileDialog.asksaveasfilename(parent=self.root,title="Save SBML as...",filetypes=[("SBML","*.sbml"),("XML","*.xml")])
          else:
                sbml_name=None 
          parameter_confidence_interval_dict , flux_confidence_interval_dict, chi_parameters_sets_dict,constrained_model =estimate_confidence_intervals(self.label_model,significance=signficance,perturbation=setp_size,min_absolute_perturbation=0.01,max_absolute_perturbation=self.label_model.p_dict["confidence_max_absolute_perturbation"],parameter_precision=self.parameter_precision,best_parameter_dict=self.label_model.parameter_dict,parameter_list=parameter_list,fraction_of_optimum=self.fraction_of_optimum,relative_max_random_sample=relative_max_random_sample, relative_min_random_sample= relative_min_random_sample,annealing_n=annealing_n,annealing_m=annealing_m,annealing_p0=annealing_p0,annealing_pf=annealing_pf,output=output,annealing_n_processes=annealing_n_processes,annealing_cycle_time_limit=annealing_cycle_time_limit, annealing_cycle_max_attempts=annealing_cycle_max_attempts,fname=fname,sbml_name=sbml_name,annealing_iterations=self.annealing_n_iterations)
          if self.confidence_constraint.get()==1:
             self.label_model.constrained_model=constrained_model
          self.waiting.destroy()
          self.confidence_window.destroy()
              
             
      def create_waiting_window(self):
          self.waiting= Toplevel()
          self.waiting.title("Working...")
          frame=Frame(self.waiting)
          frame.pack() 
          Label(frame,text="\n\nWorking, please wait...\n\n").pack(side=LEFT)
          #text_frame=Frame(self.waiting).pack()           
          #Label(text_frame,text="Please wait..").pack()
          
          
      def fitting(self):
          self.create_waiting_window()
          self.create_figure_window()
          #waiting= Toplevel()
          #waiting.title("About this application...")
          #text_frame=Frame(waiting).pack()
          #Label(text_frame,text="Fitting in progress..").pack()
          
          self.root.update_idletasks()
          #waiting.update()
          #waiting.update_idletasks()
          #label_model.parameter_dict={}
          #max_d=float(self.max_dframeentry.get())
          self.annealing_p0=p0=float(self.p0frameentry.get())
          self.annealing_pf=pf=float(self.pfframeentry.get())
          self.annealing_n=n=int(self.nframeentry.get())
          self.annealing_m=m=int(self.mframeentry.get())
          fraction_of_optimum=float(self.fraction_of_optimum)
          max_perturbation=self.label_model.p_dict['annealing_max_perturbation']#2#float(self.max_perturbationframeentry.get())
          max_sample=min(int(self.max_sampleframeentry.get()),len(self.label_model.parameter_dict))
          min_sample=min(int(self.min_sampleframeentry.get()),len(self.label_model.parameter_dict))
          self.relative_max_random_sample=float(max_sample)/len(self.label_model.parameter_dict)
          self.relative_min_random_sample=float(min_sample)/len(self.label_model.parameter_dict)
          self.annealing_n_processes=n_processes=int(self.n_process_entry.get())
          self.annealing_n_iterations=n_iterations=int(self.annealing_n_iterations_entry.get())
          cycle_time_limit=(m*50)
          print [max_sample,min_sample]
          """max_d=0.1
          p0=0.4
          pf=0.01
          n=10
          m=10
          fraction_of_optimum=1.0
          max_perturbation=1.0"""
          """if self.label_model.parameter_dict=={}:
             identify_parameters(self.label_model,add_to_parameter_dict=True,fraction_of_optimum=fraction_of_optimum,change_threshold=0.01,max_d=max_d,add_turnover=True"""
          #print [max_d,p0,pf,n,m,max_perturbation]
          self.fitting_window.destroy()
          f_best=999999999999
          backup_parameters=copy.deepcopy(self.label_model.parameter_dict)
          f_list=[]
          for x in range(0,self.annealing_n_iterations):
             if x>0:
              clear_parameters(self.label_model,parameter_dict=None,parameter_list=[], clear_ratios=True,clear_turnover=False,clear_fluxes=True,restore_objectives=False) #Clear previous parameters
              for parameter in self.label_model.parameter_dict: 
                lb_list=[]
                ub_list=[]
                if self.label_model.parameter_dict[parameter]["type"]=="flux value":
                   for reaction_id in  self.label_model.parameter_dict[parameter]["reactions"]:
                       fva=flux_variability_analysis(self.label_model.constrained_model, reaction_list=[reaction_id],fraction_of_optimum=0,tolerance_feasibility=self.label_model.lp_tolerance_feasibility)
                       lb_list.append(fva[reaction_id]["minimum"])
                       ub_list.append(fva[reaction_id]["maximum"])
                   new_value=random.uniform(max(lb_list),min(ub_list))
                   self.label_model.parameter_dict[parameter]["v"]=new_value
                elif self.label_model.parameter_dict[parameter]["type"]=="turnover": 
                     for reaction_id in  self.label_model.parameter_dict[parameter]["reactions"]:       
                         lb_list.append(self.label_model.turnover_flux_dict[reaction_id]["lb"])
                         ub_list.append(self.label_model.turnover_flux_dict[reaction_id]["ub"])
                     new_value=round(random.uniform(max(lb_list),min(ub_list)),self.precision) 
                     self.label_model.parameter_dict[parameter]["v"]=new_value
                apply_parameters(self.label_model,parameter_list=[parameter],parameter_precision=self.parameter_precision)
             #self.label_model.parameter_dict=copy.deepcopy(backup_parameters)
             
             parameters,best_flux_dict, f=annealing(self.label_model,max_random_sample=max_sample,min_random_sample=min_sample,n=n,m=m,p0=p0,pf=pf,parameter_precision=self.parameter_precision,max_perturbation=max_perturbation,fraction_of_optimum=fraction_of_optimum,gui=self,n_processes=n_processes,cycle_time_limit=cycle_time_limit)
             if f<f_best:
                f_best=f
                best_parameters=copy.deepcopy(parameters)
             f_list.append(f)
             logfile=open(self.log_name, "a")
             logfile.write(time.strftime("%c")+("//Annealing iteration %s of %s completed. Best Chi achieved is %s\n")%(x+1,self.annealing_n_iterations,f))
             logfile.close()
          self.label_model.parameter_dict=best_parameters
          apply_parameters(self.label_model,parameter_precision=self.parameter_precision)
          #best_parameters,best_flux_dict, f_best=annealing(self.label_model,max_random_sample=max_sample,min_random_sample=min_sample,n=n,m=m,p0=p0,pf=pf,parameter_precision=self.parameter_precision,max_perturbation=max_perturbation,fraction_of_optimum=fraction_of_optimum,gui=self,n_processes=n_processes,cycle_time_limit=cycle_time_limit)
          parameters=self.label_model.parameter_dict=best_parameters
          for parameter in parameters:
             if parameters[parameter]["type"]=="turnover":
                 for reaction_id in parameters[parameter]["reactions"]:
                     if reaction_id not in self.modified_turnovers:
                        self.modified_turnovers.append(reaction_id)
          print f_list
          print "Annealing completed, best Chi is %s"%(f_best)
          #self.highlight_modified(self.turnover_selector_listbox,self.modified_turnovers) 
          self.get_bounds_objective() 
          self.run_fva()
          self.waiting.destroy()
          self.update_label()
          self.fitted_parameters=True
      
      def load_model(self):
          new_label_model=load_iso2flux_model(project_file="",sbml_name="",ask_sbml_name=False,gui=True)#load_iso2flux_model(project_file="",metabolic_model_name="",gui=True)
          self.root.after_cancel(self.get_ratio_selection_after)
          self.root.after_cancel(self.get_bounds_objective_after)
          self.root.after_cancel(self.update_search_after)
          #self.root.after_cancel(self.get_selected_turnover_after)
          self.root.after_cancel(self.get_selected_reaction_name_after)
          self.root.destroy()
          launch_gui(new_label_model)
          
          """loaded_file = tkFileDialog.askopenfile(title='Choose a file',filetypes=[("parameters",".par")]) 
          fileName=loaded_file.name
          if len(fileName)>0:
             with open(fileName, 'r') as fp:
                  loaded_data=json.load(fp)
             self.label_model.ratio_dict=loaded_data[1]
             self.label_model.turnover_flux_dict=loaded_data[2]
             self.update_ratio_list()
             apply_ratios(self.label_model.constrained_model,self.label_model.ratio_dict)
             for reaction_id in loaded_data[0]:
                 reaction=self.label_model.constrained_model.reactions.get_by_id(reaction_id)
                 reaction.upper_bound=loaded_data[0][reaction_id]["ub"]
                 reaction.lower_bound=loaded_data[0][reaction_id]["lb"]
                 reaction.objective_coefficient=loaded_data[0][reaction_id]["obj"]
             self.label_model.parameter_dict=loaded_data[3]
             self.run_fva()
             #Force to the update the display
             self.force_bounds_update=True
             self.get_bounds_objective()
             self.force_bounds_update=False"""
             
      
      def save_model(self):
          project_name=save_iso2flux_model(self.label_model,name="project",write_sbml=True,gui=True)
          self.label_model.project_name=project_name
          self.log_name=project_name[:-9]+"_log.txt"
          self.root.title(self.label_model.project_name)
          f=open(self.log_name, "w")
          f.write(time.strftime("%c")+"//"+"Instance created\n")
          """fileName = tkFileDialog.asksaveasfilename(parent=self.root,title="Save the parameters as...",filetypes=[("parameters",".par")])
          if len(fileName)>0:
              reaction_dict={}
              for reaction in self.label_model.constrained_model.reactions:
                  if "RATIO_" in reaction.id:
                     continue
                  reaction_dict[reaction.id]={"lb":reaction.lower_bound,"ub":reaction.upper_bound,"obj":reaction.objective_coefficient}
              store_paremeters=[reaction_dict,self.label_model.ratio_dict,self.label_model.turnover_flux_dict,self.label_model.parameter_dict]  
              with open(fileName, 'w') as fp:
                   json.dump(store_paremeters, fp)
          print fileName"""
      
      
      
      def get_ratio_selection(self):
          previous_selected_ratio=self.selected_ratio
          self.get_ratio_selection_after=self.root.after(200, self.get_ratio_selection)
          selection_tupple=self.ratio_selector_listbox.curselection()
          if selection_tupple!=(): 
             print ["st",selection_tupple]
             self.selected_ratio,=selection_tupple
          if previous_selected_ratio!=self.selected_ratio:
             ratio_string=self.ratio_selector_listbox.get(self.ratio_selector_listbox.curselection())
             ratio_id=self.ratio_regular_expression.match(ratio_string).group(1)
             if ratio_id not in self.label_model.ratio_dict:
                ratio_id=self.ratio_inverse_name_dict[ratio_id]
             print ratio_id
             if len(self.label_model.ratio_dict[ratio_id])>2:
                #raise("Ratio with more than 2 elements not supported in GUI mode")
                print "Ratio with more than 2 elements not supported in GUI mode"
             reaction_list=self.label_model.ratio_dict[ratio_id].keys()
             if self.label_model.ratio_dict[ratio_id][reaction_list[0]]==1: 
                self.selected_r_ratio1=reaction_list[1]
                self.selected_r_ratio2=reaction_list[0]
             elif self.label_model.ratio_dict[ratio_id][reaction_list[1]]==1:
                self.selected_r_ratio1=reaction_list[0]
                self.selected_r_ratio2=reaction_list[1]
             else:
                 print "At least one coefficient must be one in GUI mode"
                #raise("At least one coefficient must be one in GUI mode")
             #print[self.inverse_llistafluxos_dict[self.selected_r_ratio1],self.inverse_llistafluxos_dict[self.selected_r_ratio2]]
             self.ratios_spinbox.delete(0,"end")
             self.ratios_spinbox.insert(INSERT,str(self.label_model.ratio_dict[ratio_id][self.selected_r_ratio1])) 
      
      
      def update_bounds(self):
          self.get_bounds_objective()
          changed=False
          #self.root.after(50, self.update_bounds)
          reaction_id=self.active_reaction_get_bounds
          reaction=self.label_model.constrained_model.reactions.get_by_id(reaction_id)
          try:
            new_reaction_lb=float(self.boundaries_lb_spinbox.get())
          except:
            new_reaction_lb=self.boundaries_lb_spinbox.get()
            new_reaction_lb=float(new_reaction_lb.replace(",","."))
          try:
            new_reaction_ub=float(self.boundaries_ub_spinbox.get())
          except:
            new_reaction_ub=self.boundaries_ub_spinbox.get()
            new_reaction_ub=float(new_reaction_ub.replace(",","."))
                
            
          """if self.fraction_of_optimum!=new_fraction_of_optimum:
                self.fraction_of_optimum=new_fraction_of_optimum
                changed=True
          if reaction.lower_bound!=new_reaction_lb:
                reaction.lower_bound=new_reaction_lb
                changed=True
          if reaction.upper_bound!=new_reaction_ub:
                reaction.upper_bound=new_reaction_ub 
                changed=True
          if reaction.objective_coefficient!=self.varoptimizacion.get():
                reaction.objective_coefficient=self.varoptimizacion.get()
                changed=True"""
          if reaction_id not in self.modified_reactions:
                 self.modified_reactions.append(reaction_id)
                 #self.boundaries_reaction_selector_listbox.itemconfig(self.variable13b1, {'fg': 'green'}) 
          reaction.lower_bound=new_reaction_lb
          reaction.lower_bound=new_reaction_lb
          reaction.upper_bound=new_reaction_ub
          reaction.upper_bound=new_reaction_ub 
          reaction.objective_coefficient=self.varoptimizacion.get()
          if self.active_reaction_get_bounds in self.label_model.turnover_flux_dict:
             try:
               turnover_lb=float(self.boundaries_lb_turnover_spinbox.get())
             except:
                turnover_lb=float(self.boundaries_lb_turnover_spinbox.get().replace(",","."))
             try:
               turnover_ub=float(self.boundaries_ub_turnover_spinbox.get())
             except:
                turnover_ub=float(self.boundaries_ub_turnover_spinbox.get().replace(",","."))
             self.label_model.turnover_flux_dict[self.active_reaction_get_bounds]["lb"]=turnover_lb
             self.label_model.turnover_flux_dict[self.active_reaction_get_bounds]["ub"]=turnover_ub
             self.label_model.turnover_flux_dict[self.active_reaction_get_bounds]["v"]=(turnover_lb+turnover_ub)/2
             print self.label_model.turnover_flux_dict[self.active_reaction_get_bounds]
          self.highlight_modified(self.boundaries_reaction_selector_listbox,self.modified_reactions)
          self.waiting_for_fva=True
          self.run_fva()
          
      def reset_bound(self):
          reaction_id=self.active_reaction_get_bounds
          orginal_reaction=self.backup_constrained_model.reactions.get_by_id(reaction_id)
          lb=orginal_reaction.lower_bound
          ub=orginal_reaction.upper_bound
          obj=orginal_reaction.objective_coefficient
          reaction_to_reset=self.label_model.constrained_model.reactions.get_by_id(reaction_id)
          reaction_to_reset.lower_bound=lb
          reaction_to_reset.uppee_bound=ub
          reaction_to_reset.objective_coefficient=obj
          self.boundaries_lb_spinbox.delete(0,"end")
          self.boundaries_lb_spinbox.insert(INSERT,str(lb))
          self.boundaries_ub_spinbox.delete(0,"end")
          self.boundaries_ub_spinbox.insert(INSERT,str(ub))
          self.varoptimizacion.set(int(obj))
          if reaction_id in self.label_model.parameter_dict:
             del self.label_model.parameter_dict[reaction_id]
          if reaction_id in self.label_model.turnover_flux_dict:
             if reaction_id+"_turnover" in self.label_model.parameter_dict:
                del self.label_model.parameter_dict[reaction_id+"_turnover"]
             self.label_model.turnover_flux_dict[reaction_id]=copy.deepcopy(self.backup_turnover_flux_dict[reaction_id])
             self.boundaries_lb_turnover_spinbox.delete(0,"end")
             self.boundaries_lb_turnover_spinbox.insert(INSERT,str(self.label_model.turnover_flux_dict[reaction_id]["lb"]))
             self.boundaries_ub_turnover_spinbox.delete(0,"end")
             self.boundaries_ub_turnover_spinbox.insert(INSERT,str(self.label_model.turnover_flux_dict[reaction_id]["ub"]))
          self.run_fva()
          if reaction_to_reset.id in self.label_model.parameter_dict:
             del self.label_model.parameter_dict[reaction_to_reset.id]
          if reaction_to_reset.id in self.modified_reactions:
             self.modified_reactions.remove(reaction_to_reset.id)
          #search_term=self.search_entry.get()
          #self.boundaries_reaction_selector_listbox.delete(0, END)
          #reaction_list=self.search_function(search_term) 
  
          self.highlight_modified(self.boundaries_reaction_selector_listbox,self.modified_reactions)                    
          """for n,reaction in enumerate(sorted(reaction_list,key=lambda v: v.upper())):
              self.boundaries_reaction_selector_listbox.insert(END,reaction)
              if reaction in self.modified_reactions:
                  self.boundaries_reaction_selector_listbox.itemconfig(n, {'fg': 'green'})"""
          
          
      def reset_all_bounds(self):
          self.label_model.constrained_model=copy.deepcopy(self.backup_constrained_model)
          """for reaction_id in self.modified_reactions:
              #self.boundaries_reaction_selector_listbox.itemconfig(self.inverse_llistafluxos_dict[reaction_id], {'fg': 'black'}) 
              reaction=self.backup_constrained_model.reactions.get_by_id(reaction_id)
              lb=reaction.lower_bound
              ub=reaction.upper_bound
              obj=reaction.objective_coefficient
              reaction_to_reset=self.label_model.constrained_model.reactions.get_by_id(reaction_id)
              reaction_to_reset.lower_bound=lb
              reaction_to_reset.upper_bound=ub
              reaction_to_reset.objective_coefficient=obj
              if reaction_to_reset.id==self.active_reaction_get_bounds:
                 self.boundaries_lb_spinbox.delete(0,"end")
                 self.boundaries_lb_spinbox.insert(INSERT,str(lb))
                 self.boundaries_ub_spinbox.delete(0,"end")
                 self.boundaries_ub_spinbox.insert(INSERT,str(ub))
                 self.varoptimizacion.set(int(obj))"""
          self.run_fva()
          self.modified_reactions=[]
          self.label_model.parameter_dict={}
          search_term=self.search_entry.get()
          #reaction_list=[x.id for x in self.label_model.metabolic_model.reactions.query(search_term)]
          #self.boundaries_reaction_selector_listbox.delete(0, END)
          self.highlight_modified(self.boundaries_reaction_selector_listbox,self.modified_reactions)
          print "Boundaries reset"
          
      
      def add_ratio(self):
          reaction_id_a=self.selected_r_ratio1
          reaction_id_b=self.selected_r_ratio2
          ratio_id1=reaction_id_a+"/"+reaction_id_b
          ratio_id2=reaction_id_b+"/"+reaction_id_a
          if ratio_id2 in self.label_model.ratio_dict:
             remove_ratios(self.label_model.constrained_model,ratio_id2,self.label_model.ratio_dict)
          self.ratio_inverse_name_dict[ratio_id2]=ratio_id1
          try: 
             ratio_value=float(self.ratios_spinbox.get())
          
          except:
             ratio_value=(self.ratios_spinbox.get())
             ratio_value=float(ratio_value.replace(",","."))
             print ratio_value 
          self.label_model.ratio_dict[ratio_id1]={reaction_id_a:ratio_value,reaction_id_b:1.0}
          self.update_ratio_list()
          apply_ratios(self.label_model.constrained_model,self.label_model.ratio_dict)
          print "Ratio_added"
          
          self.run_fva()
      #123456789
      """def get_r1_r2(self): # Select Ratio
          self.root.after(200, self.get_r1_r2)
          selection_tupple=self.flux_ratio_selector1_listbox.curselection()
          previous_variable31b1_a=self.variable31b1_a
          previous_variable31b1_b=self.variable31b1_b  
          if selection_tupple!=(): 
             self.variable31b1_a,=selection_tupple
             new_r1_ratio=self.flux_ratio_selector1_listbox.get(self.variable31b1_a) 
             if previous_variable31b1_a!=self.variable31b1_a:
                self.selected_r_ratio1 =new_r1_ratio
          selection_tupple=self.flux_ratio_selector2_listbox.curselection()
          if selection_tupple!=(): 
             self.variable31b1_b,=selection_tupple
             new_r2_ratio=self.flux_ratio_selector2_listbox.get(self.variable31b1_b) 
             if previous_variable31b1_b!=self.variable31b1_b:
                self.flux_ratio_selector2_listbox.itemconfig(self.variable31b1_b, {'fg': 'blue'}) 
                self.flux_ratio_selector2_listbox.itemconfig(previous_variable31b1_b, {'fg': 'black'})
                if self.selected_r_ratio2!=None:
                   self.flux_ratio_selector2_listbox.itemconfig(self.inverse_llistafluxos_dict[self.selected_r_ratio2], {'fg': 'black'})  
                self.selected_r_ratio2 =new_r2_ratio
          if self.selected_r_ratio1!=None and self.selected_r_ratio2!=None:
             #sel2 = self.llistafluxos_dict[self.variable31b1_a] + " / " +self.llistafluxos_dict[self.variable31b1_b] + " = " + str(self.ratios_spinbox.get()) + "; "
             #lab31b1.config(text = sel2)
             self.sublabelframe_ratios_spinbox.config(text = self.selected_r_ratio1 + " / " +self.selected_r_ratio2)"""
      
      
      def delete_ratio(self):
          ratio_string=self.ratio_selector_listbox.get(self.ratio_selector_listbox.curselection())
          ratio_id=self.ratio_regular_expression.match(ratio_string).group(1)
          print ratio_id
          if ratio_id not in self.label_model.ratio_dict:
             ratio_id=self.ratio_inverse_name_dict[ratio_id]
          remove_ratio(self.label_model.constrained_model,ratio_id,self.label_model.ratio_dict)
          self.update_ratio_list()
          self.run_fva()
     
      
      def update_ratio_list(self):
         self.ratio_selector_listbox.delete(0, "end")
         for ratio_id in self.label_model.ratio_dict:
            ratio=1.0
            for flux in self.label_model.ratio_dict[ratio_id]:
                if self.label_model.ratio_dict[ratio_id][flux]!=1.0:
                   ratio=self.label_model.ratio_dict[ratio_id][flux]
            string=ratio_id+"="+str(ratio)  
            self.ratio_selector_listbox.insert(END, string)
      
      def run_fva(self): #TODO Rename function
          self.waiting_for_fva=False
          self.run_label_simulation() 
          print self.fba_mode.get() 
          try:
            if self.fba_mode.get().lower()=="pfba":
               #convert_model_to_irreversible(self.label_model.constrained_model)
               optimize_minimal_flux(self.label_model.constrained_model, already_irreversible=False)
               #revert_model_to_reversible(self.label_model.constrained_model, update_solution=True)
            else:
               self.label_model.constrained_model.optimize(tolerance_feasibility=self.label_model.lp_tolerance_feasibility)
            print self.label_model.constrained_model.solution
            optimal_solution=self.label_model.constrained_model.solution
            self.solution_dict=optimal_solution.x_dict
            self.fva=flux_variability_analysis(self.label_model.constrained_model,fraction_of_optimum=self.fraction_of_optimum,tolerance_feasibility=self.label_model.lp_tolerance_feasibility)
            print "FVA updated"
          except:
            self.fva={}
            print "solution not optimal"
          self.update_fva_listbox()
       
      def update_fva_listbox(self):
          search_term=self.search_entry.get()
          reaction_list=self.search_function(search_term)                         
          self.fva_listbox.delete(0,"end")
          self.optimal_solution_listbox.delete(0, "end")
          self.free_fluxes_n_dict={}
          self.n_free_fluxes_dict={}
          n=0
          for flux in sorted(reaction_list,key=lambda v: v.upper()):
              if "RATIO_" in flux:
                 continue
              """optimal_flux=round(solution_dict[flux],3)
              string=(flux+"=%s ")%( optimal_flux) 
              self.optimal_solution_listbox.insert(END, string)"""
              try:
                max_flux=self.fva[flux]["maximum"]
                min_flux=self.fva[flux]["minimum"]
                """if max_flux==min_flux:
                   string=("    "+flux+"=%s ")%( min_flux)
                   else:
                   string=("%s < "+flux+" <%s ")%( min_flux, max_flux)"""
                if (round(min_flux,4)!=round(max_flux,4) and (max_flux-min_flux)>1.1*self.change_threshold) and flux not in self.label_model.parameter_dict:#if (max_flux-min_flux)>1.5*self.change_threshold and flux not in self.label_model.parameter_dict:
                 string=("%s < "+flux+" <%s ")%( round(min_flux,3), round(max_flux,3))
                 self.fva_listbox.insert(END, string)
                 self.free_fluxes_n_dict[flux]=n
                 self.n_free_fluxes_dict[n]=flux
                 n+=1
                else:
                  optimal_flux=round((max_flux+min_flux)/2,4)
                  string=(flux+"=%s ")%( optimal_flux) 
                  self.optimal_solution_listbox.insert(END, string)
              except:
                  string=flux+"=NAN"
                  self.optimal_solution_listbox.insert(END, string)


      def run_label_simulation(self):
          self.waiting_for_label_simulation=False 
          if self.label_window==True:
             a1,b=solver(self.label_model,mode="fsolve",fba_mode=self.fba_mode.get().lower())
             #a2,b,c=get_objective_function(self.label_model,force_balance=self.label_model.force_balance,output=False)
             self.update_label()
           
      
      
      
      def get_bounds_objective(self): #Assign bounds and objective based on reaction
          #print variable13b1
          #active_reaction=self.variable13b1
          self.get_bounds_objective_after=self.root.after(100, self.get_bounds_objective)
          selection_tupple=self.boundaries_reaction_selector_listbox.curselection()
          if selection_tupple!=(): 
             position,=selection_tupple
             new_active_reaction=self.boundaries_reaction_selector_listbox.get(position)
          else:
             new_active_reaction=self.active_reaction_get_bounds
          #print selection_tupple
          #sel1=str(variable13b1)
          if new_active_reaction!=self.active_reaction_get_bounds or self.force_bounds_update==True:
             """if self.llistafluxos_dict[previous_variable13b1] not in self.modified_reactions:
                self.boundaries_reaction_selector_listbox.itemconfig(previous_variable13b1, {'fg':'black'})
             if self.llistafluxos_dict[self.variable13b1] not in self.modified_reactions:
                self.boundaries_reaction_selector_listbox.itemconfig(self.variable13b1, {'fg':'blue'})"""
             self.active_reaction_get_bounds=new_active_reaction
             reaction=self.label_model.constrained_model.reactions.get_by_id(self.active_reaction_get_bounds)
             #Set current value for lower bound
             self.boundaries_lb_spinbox.delete(0,"end")
             lower_bound_var=reaction.lower_bound
             self.boundaries_lb_spinbox.insert(INSERT,str(lower_bound_var))
             #Set current value for upper bound
             self.boundaries_ub_spinbox.delete(0,"end")
             upper_bound_var=reaction.upper_bound
             self.boundaries_ub_spinbox.insert(INSERT,str(upper_bound_var))
             #Set current objective
             obj=reaction.objective_coefficient
             print reaction.objective_coefficient
             self.varoptimizacion.set(int(obj))
             if self.label_model.emu_dict!={}: #Create the fitting menu only if a label model exists
                if self.active_reaction_get_bounds in self.label_model.parameter_dict:
                   self.parameter_button.config(text="Remove from parameters")
                else:
                   self.parameter_button.config(text="Make parameter")
             text=(self.active_reaction_get_bounds)
             if len(text)>=2*(self.average_width+1):
                text=text[:(2*(self.average_width))]+".."
             self.sublabelframe_boundaries_rid.config(text= text)
             if self.active_reaction_get_bounds in self.label_model.turnover_flux_dict:
                self.turnover_filler.pack_forget()
                if not self.frame_turnover_boundaries.winfo_ismapped():
                   self.frame_turnover_boundaries.pack(side=TOP)
                   self.add_parameter_turnover_button.pack(side=TOP)
                if self.active_reaction_get_bounds+"_turnover" in self.label_model.parameter_dict:
                   lower_bound=upper_bound=self.label_model.parameter_dict[self.active_reaction_get_bounds+"_turnover"]["v"]
                   self.add_parameter_turnover_button.config(text="Remove from parameters")
                else:
                   lower_bound=self.label_model.turnover_flux_dict[self.active_reaction_get_bounds]["lb"]
                   upper_bound=self.label_model.turnover_flux_dict[self.active_reaction_get_bounds]["ub"]
                   self.add_parameter_turnover_button.config(text="Make parameter")
                self.boundaries_lb_turnover_spinbox.delete(0,"end")
                self.boundaries_lb_turnover_spinbox.insert(INSERT,lower_bound)
                self.boundaries_ub_turnover_spinbox.delete(0,"end")
                self.boundaries_ub_turnover_spinbox.insert(INSERT,upper_bound)
             elif self.frame_turnover_boundaries.winfo_ismapped():
                self.frame_turnover_boundaries.pack_forget()
                self.add_parameter_turnover_button.pack_forget()
                self.turnover_filler.pack(side=TOP)
          
          #lab13b1.config(text = sel1)
          #print (boundaries_lb_spinbox.get())
          
          
      """def poll13b2(self):
          
          selection_tupple=self.boundaries_reaction_selector_listbox.curselection()
          if selection_tupple!=(): 
             self.variable13b2,=selection_tupple
          #print "poll13b2"
          self.root.after(200, self.poll13b2)
          sel12 = " del " + (self.boundaries_reaction_selector_listbox.get(self.variable13b1)) + "; "
          #lab13b2.config(text = sel12)"""
      
      
      
      """def input_frame(self):
              labelframe_input = LabelFrame(self.root, text="INPUT")
              labelframe_input.pack(fill=None, expand=False)
              button01 = Button(labelframe_input, text="Import SBML", fg="red")
              button01.pack( side = LEFT)
              button02 = Button(labelframe_input, text="Import Excell", fg="brown")
              button02.pack( side = LEFT )"""
      
      def reaction_selector(self,frame):
          reaction_selector = Frame(frame)
          reaction_sub_selector = Frame(reaction_selector)
          #entry_search=Entry(reaction_selector,width=14)
          #entry_search.pack(side=TOP)
          #average_width=int(numpy.mean([len(x.id) for x in label_model.constrained_model.reactions]))
          reaction_sub_selector.pack(side=TOP,fill=Y,expand="yes")
          reaction_selector_scrollbar = Scrollbar(reaction_sub_selector)
          reaction_selector_xscrollbar = Scrollbar(reaction_sub_selector,orient=HORIZONTAL)
          reaction_selector_scrollbar.pack( side = RIGHT, fill=Y,expand="yes")
          reaction_selector_listbox = Listbox(reaction_sub_selector, yscrollcommand =reaction_selector_scrollbar.set , height = 3, width=self.average_width+7)
          reaction_selector_scrollbar['command'] = reaction_selector_listbox.yview
          reaction_selector_xscrollbar['command'] = reaction_selector_listbox.xview
          reaction_selector_listbox['yscrollcommand'] = reaction_selector_scrollbar.set
          reaction_selector_listbox['xscrollcommand'] = reaction_selector_xscrollbar.set
          reaction_selector_listbox.pack( side = TOP, fill = BOTH ,expand="yes")
          reaction_selector_xscrollbar.pack( side = TOP, fill=X,expand=False) 
          reaction_selector_scrollbar.config( command = reaction_selector_listbox.yview )
          reaction_selector_xscrollbar.config( command = reaction_selector_listbox.xview )
          reaction_list=[x.id for x in self.label_model.metabolic_model.reactions]
          for reaction in sorted(reaction_list,key=lambda v: v.upper()):
                 reaction_selector_listbox.insert(END,reaction)
          return(reaction_selector,reaction_selector_listbox) 
      
      def add_boundaries(self):
         "Flux Boundaries and objective"
         labelframe_flux_boundaries_master = LabelFrame(self.root, text="Set flux bounds and optimization goals")
         labelframe_flux_boundaries_master.pack(fill="y", expand=False)
         labelframe_boundaries_objective_turnover = LabelFrame(labelframe_flux_boundaries_master, text="")
         labelframe_boundaries_objective_turnover.pack(fill="y", expand=True)
         
         labelframe_boundaries_buttons = Frame(labelframe_boundaries_objective_turnover)
         apply_button=Button(labelframe_boundaries_buttons, text="Apply", command = self.update_bounds)
         reset_bounds_button=Button(labelframe_boundaries_buttons, text="Reset", command = self.reset_bound)
         #reset_all_bounds_button = Button(labelframe_boundaries_buttons, text="Reset All", command = self.reset_all_bounds)
         
         
         apply_button.pack(side = TOP)
         reset_bounds_button.pack(side=TOP)
         #reset_all_bounds_button.pack( side = TOP)
         self.frame_boundaries = Frame(labelframe_boundaries_objective_turnover)
         frame_flux_boundaries=Frame(self.frame_boundaries)
         frame_flux_boundaries.pack(side=TOP,fill="y", expand=True)
         if self.label_model.emu_dict!={}: #Create the fitting menu only if a label model exists
            self.parameter_button=add_parameter_button = Button(self.frame_boundaries, text="Make parameter",command=self.add_parameter)
            add_parameter_button.pack(side=TOP)
         sublabelframe_boundaries_lb = Frame(frame_flux_boundaries)
         self.boundaries_lb_spinbox = Spinbox(sublabelframe_boundaries_lb, from_=-99999, to=99999, increment = 0.01, format="%.2f", width = 7)
         reaction_selector,self.boundaries_reaction_selector_listbox=self.reaction_selector(labelframe_boundaries_objective_turnover)
         self.active_reaction_get_bounds=self.boundaries_reaction_selector_listbox.get(0)
         #self.boundaries_reaction_active_search=""
         #average_width=int(numpy.mean([len(x.id) for x in label_model.constrained_model.reactions] ))
         self.sublabelframe_boundaries_rid = Label(frame_flux_boundaries,width=2*(self.average_width+1))
         sublabelframe_boundaries_ub = Frame( frame_flux_boundaries)
         #Label(sublabelframe_boundaries_ub, text="<").pack( side = LEFT)
         self.boundaries_ub_spinbox = Spinbox(sublabelframe_boundaries_ub, from_=-99999, to=99999, increment = 0.01, format="%.2f", width = 7)
         upper_bound_initial_reaction=self.label_model.constrained_model.reactions.get_by_id(self.boundaries_reaction_selector_listbox.get(0)).upper_bound
         self.boundaries_ub_spinbox.delete(0,"end")
         self.boundaries_ub_spinbox.insert(INSERT,upper_bound_initial_reaction)
         self.boundaries_ub_spinbox.pack( side = LEFT)
         
         lower_bound_initial_reaction=self.label_model.constrained_model.reactions.get_by_id(self.boundaries_reaction_selector_listbox.get(0)).lower_bound
         self.boundaries_lb_spinbox.delete(0,"end")
         self.boundaries_lb_spinbox.insert(INSERT,lower_bound_initial_reaction)
         self.boundaries_lb_spinbox.pack( side = LEFT)
         
         sublabelframe_optimization = LabelFrame(labelframe_boundaries_objective_turnover,text="Objective")
         
         self.varoptimizacion = IntVar()
         self.varoptimizacion.set(0)
         radiobutton411 = Radiobutton(sublabelframe_optimization, text="none", variable=self.varoptimizacion, value=0)
         radiobutton411.pack( side = TOP)
         radiobutton412 = Radiobutton(sublabelframe_optimization, text="minimize", variable=self.varoptimizacion, value=-1)
         radiobutton412.pack( side = TOP)
         radiobutton413 = Radiobutton(sublabelframe_optimization, text="maximize", variable=self.varoptimizacion, value=1)
         radiobutton413.pack( side = TOP)
         
         greater_label1=Label(frame_flux_boundaries, text=" > ")
         greater_label2=Label(frame_flux_boundaries, text=" > ")
         
         reaction_selector.pack(side = LEFT,fill=Y,expand=True)
         self.frame_boundaries.pack(side=LEFT)
         Label(sublabelframe_boundaries_lb, text="").pack( side = LEFT)
         sublabelframe_boundaries_lb.pack(side = LEFT)
         greater_label1.pack(side = LEFT)
         self.sublabelframe_boundaries_rid.pack(side=LEFT)
         greater_label2.pack(side = LEFT)
         sublabelframe_boundaries_ub.pack(side = LEFT)
         sublabelframe_optimization.pack(side=LEFT)
         labelframe_boundaries_buttons.pack(side = LEFT)
         
         self.frame_turnover_boundaries=frame_turnover_boundaries=Frame(self.frame_boundaries)
         frame_turnover_boundaries.pack(side=TOP)
         self.add_parameter_turnover_button= Button(self.frame_boundaries, text="Make parameter",command=self.add_turnover_as_parameter)
         self.turnover_filler=Label(self.frame_boundaries,text="\n\n")
         self.turnover_filler.pack(side=TOP)
         self.add_parameter_turnover_button.pack(side=TOP)
         
         sublabelframe_turnover_boundaries_lb = Frame(frame_turnover_boundaries)
         self.boundaries_lb_turnover_spinbox = Spinbox(sublabelframe_turnover_boundaries_lb, from_=0, to=99999, increment = 0.01, format="%.2f", width = 7)
         self.boundaries_lb_turnover_spinbox.pack( side = LEFT)
         
         self.sublabelframe_turnover_rid = Label(frame_turnover_boundaries)
         self.sublabelframe_turnover_rid.config(text ="Reversible turnover")
         
         sublabelframe_turnover_ub = Frame( frame_turnover_boundaries)
         self.boundaries_ub_turnover_spinbox = Spinbox(sublabelframe_turnover_ub, from_=-99999, to=99999, increment = 0.01, format="%.2f", width = 7)
         self.boundaries_ub_turnover_spinbox.pack( side = LEFT)
         
         greater_turnover_label1=Label(frame_turnover_boundaries, text=" > ")
         greater_turnover_label2=Label(frame_turnover_boundaries, text=" > ")
         
         sublabelframe_turnover_boundaries_lb.pack(side=LEFT)
         greater_turnover_label1.pack(side=LEFT)
         self.sublabelframe_turnover_rid.pack(side=LEFT)
         greater_turnover_label2.pack(side=LEFT)
         sublabelframe_turnover_ub.pack(side=LEFT)
         self.sublabelframe_boundaries_rid.config(text =self.active_reaction_get_bounds)
         if self.boundaries_reaction_selector_listbox.get(0) in self.label_model.turnover_flux_dict:
            self.turnover_filler.pack_forget()
            lower_bound=self.label_model.turnover_flux_dict[self.boundaries_reaction_selector_listbox.get(0)]["lb"]
            upper_bound=self.label_model.turnover_flux_dict[self.boundaries_reaction_selector_listbox.get(0)]["ub"]
            self.boundaries_lb_turnover_spinbox.delete(0,"end")
            self.boundaries_lb_turnover_spinbox.insert(INSERT,lower_bound)
            self.boundaries_ub_turnover_spinbox.delete(0,"end")
            self.boundaries_ub_turnover_spinbox.insert(INSERT,upper_bound)
            if self.boundaries_reaction_selector_listbox.get(0)+"_turnover" in self.label_model.parameter_dict:
               self.add_parameter_turnover_button.config(text="Remove from parameters")
         else:
            self.frame_turnover_boundaries.pack_forget()
            self.add_parameter_turnover_button.pack_forget()
         if  self.boundaries_reaction_selector_listbox.get(0) in self.label_model.parameter_dict:
          self.parameter_button.config(text="Remove from parameters")
          
          
                
      def update_search(self):
          self.update_search_after=self.root.after(200,self.update_search)
          ###Reaction bounds
          search_term=self.search_entry.get()
          if self.active_search!=search_term or self.active_filter!=self.show_variable.get():
             self.active_search=search_term
             self.active_filter=self.show_variable.get()
             
             reaction_list=self.search_function(search_term)                       
             self.boundaries_reaction_selector_listbox.delete(0, END)
             self.flux_ratio_selector1_listbox.delete(0, END)
             #self.turnover_selector_listbox.delete(0,END)
             n_turnover=0
             for n,reaction in enumerate(sorted(reaction_list,key=lambda v: v.upper())):
                 self.boundaries_reaction_selector_listbox.insert(END,reaction)
                 self.flux_ratio_selector1_listbox.insert(END,reaction)
                 if reaction in self.modified_reactions:
                    self.boundaries_reaction_selector_listbox.itemconfig(n, {'fg': 'green'})
                 """if reaction in self.label_model.turnover_flux_dict:
                    self.turnover_selector_listbox.insert(END,reaction)
                    if reaction in self.modified_turnovers:
                       self.turnover_selector_listbox.itemconfig(n_turnover, {'fg': 'green'})
                    n_turnover+=1"""
             self.update_fva_listbox()
      
      def search_function(self,search_term):  
             reaction_list=[]
             lower_case_search_term=search_term.lower()
             #reaction_list=[x.id for x in self.label_model.metabolic_model.reactions.query(search_term)]
             #Add reactions whose name matches
             for reaction_id in self.id_name_dict:
                 reaction_name=self.id_name_dict[reaction_id]
                 lower_case_reaction_id=reaction_id.lower()
                 lower_case_reaction_name=reaction_name.lower()
                 if (lower_case_search_term in lower_case_reaction_id) or (lower_case_search_term in lower_case_reaction_name)  and reaction_id not in reaction_list:
                    if self.show_variable.get()==1 and reaction_id not in self.label_model.reactions_propagating_label:
                       continue
                    elif self.show_variable.get()==2 and ("LABEL_RGROUP_" in reaction_id):
                       continue
                    reaction_list.append(reaction_id)
             return reaction_list
      
      def add_ratios_interface(self):
        root=self.root
        labelframe_fluxratios = LabelFrame(root, text="Relative fluxes values (R1/R2)")
        labelframe_fluxratios.pack(fill=Y, expand=False)
        sublabelframe_flux_ratio_selector1,self.flux_ratio_selector1_listbox=self.reaction_selector(labelframe_fluxratios)  
        self.search_term_active_ratio=""  
        #sublabelframe_flux_ratio_selector1 = LabelFrame(labelframe_fluxratios,text="R1 ")
        #sublabelframe_flux_ratio_selector1.pack(side = LEFT)
        self.sublabelframe_ratios_spinbox = LabelFrame(labelframe_fluxratios,text="ratio")
        self.ratios_spinbox=Spinbox(self.sublabelframe_ratios_spinbox, from_=-1000, to=1000, increment = 1, format="%.2f", width = 7)
        self.ratios_spinbox.delete(0,"end")
        self.ratios_spinbox.insert(INSERT,1)
        #self.flux_ratios_scale = Scale( self.sublabelframe_ratios_scale, resolution=0.01, from_= 1, to = 100, orient = HORIZONTAL)
        #self.flux_ratios_scale.set(1)
        self.ratios_spinbox.pack( side = TOP)
        #sublabelframe_flux_ratio_selector2,self.flux_ratio_selector2_listbox,self.ratio_entry2=self.reaction_selector(labelframe_fluxratios) 
        add_ratio_button = Button(self.sublabelframe_ratios_spinbox, text="Add/Update\nRatio",command=self.add_ratio)
        sublabelframe_flux_ratio_selector1.pack(side = LEFT,fill=Y, expand=True)
        #sublabelframe_flux_ratio_selector2.pack(side = LEFT,fill=Y, expand=True)"""
        self.sublabelframe_ratios_spinbox.pack(side=LEFT)
        add_ratio_button.pack( side = TOP)
        subframe_ratio_selector = Frame(labelframe_fluxratios)
        subframe_ratio_selector.pack(side = RIGHT, fill=Y)
        ratio_selector_scrollbar = Scrollbar(subframe_ratio_selector)
        self.ratio_selector_listbox = Listbox(subframe_ratio_selector, yscrollcommand = ratio_selector_scrollbar.set , height = 3, width=14)
        ratio_selector_scrollbar['command'] = self.ratio_selector_listbox.yview
        self.ratio_selector_listbox['yscrollcommand'] = ratio_selector_scrollbar.set
        self.ratio_selector_listbox.pack( side = LEFT, fill=Y, expand=True )
        ratio_selector_scrollbar.pack( side = LEFT, fill=Y)
        ratio_selector_scrollbar.config( command = self.ratio_selector_listbox.yview )
        del_ratio_button = Button(subframe_ratio_selector, text="Delete",command=self.delete_ratio)
        del_ratio_button.pack( side = LEFT)
        set_r1 = Button(sublabelframe_flux_ratio_selector1, text="Select as R1",command=self.select_r1)
        set_r2 = Button(sublabelframe_flux_ratio_selector1, text="Select as R2",command=self.select_r2)
        set_r1.pack(side=TOP)
        set_r2.pack(side=TOP)
      
      def select_r1(self):
          selection_tupple=self.flux_ratio_selector1_listbox.curselection()
          if selection_tupple!=(): 
             selection,=selection_tupple
	     self.selected_r_ratio1=self.flux_ratio_selector1_listbox.get(selection)
             self.sublabelframe_ratios_spinbox.config(text = self.selected_r_ratio1 + " / " +self.selected_r_ratio2)    
      
      def select_r2(self):
          selection_tupple=self.flux_ratio_selector1_listbox.curselection()
          if selection_tupple!=(): 
             selection,=selection_tupple
	     self.selected_r_ratio2=self.flux_ratio_selector1_listbox.get(selection)
             self.sublabelframe_ratios_spinbox.config(text = self.selected_r_ratio1 + " / " +self.selected_r_ratio2)   
        
      def add_fraction_optimum_interface(self):
            labelframe_fraction_optimum = Frame(self.root)
            labelframe_fraction_optimum.pack(side=TOP)
            Label(labelframe_fraction_optimum,text="\nfraction of optimum:").pack(side=LEFT)
            self.fraction_optimum_spinbox = Spinbox( labelframe_fraction_optimum, increment=0.05, from_= 0, to = 1,width=5)
            self.fraction_optimum_spinbox.delete(0,"end")
            self.fraction_optimum_spinbox.insert(INSERT,float(self.fraction_of_optimum))
            self.fraction_optimum_spinbox.pack( side = LEFT)
            button=Button(labelframe_fraction_optimum, text="Apply",command=self.update_fraction_of_optimum)
            button.pack(side = LEFT)
      
      def update_fraction_of_optimum(self):
          try:
            self.fraction_of_optimum=float(self.fraction_optimum_spinbox.get())
          except:
            self.fraction_of_optimum=(self.fraction_optimum_spinbox.get())
            self.fraction_of_optimum=float(new_fraction_of_optimum.replace(",","."))
          self.run_fva()
      
      def add_flux_display(self):
         labelframe6 = LabelFrame(self.root, text="Reaction fluxes")
         #export_to_excel_button = Button(labelframe6, text="Export  Excel",command=self.export_fluxes )
         #export_to_excel_button.pack(side = TOP)
         #self.flux_display_search=Entry(labelframe6,width=14)
         #self.flux_display_search.pack(side=TOP)
         #self.search_term_flux_active=""
         #average_width=int(numpy.mean([len(x.id) for x in label_model.constrained_model.reactions] ))
         labelframe6.pack(fill=Y, expand=False)
         label_frame_determined_fluxes = LabelFrame(labelframe6,text="Fluxes with a single solution")
         label_frame_determined_fluxes.pack(side = LEFT,fill=Y,expand=True)
         scrollbar612 = Scrollbar(label_frame_determined_fluxes)
         scrollbar612.pack( side = RIGHT, fill=Y,expand=True)
         scrollbar613 = Scrollbar(label_frame_determined_fluxes,orient=HORIZONTAL)
         self.optimal_solution_listbox = Listbox(label_frame_determined_fluxes, yscrollcommand = scrollbar612.set,xscrollcommand = scrollbar613.set, width= self.average_width+20 )
         scrollbar612.config( command = self.optimal_solution_listbox.yview ) 
         scrollbar613.config( command = self.optimal_solution_listbox.xview ) 
         self.optimal_solution_listbox.pack( side = TOP, fill=Y,expand=True)
         scrollbar613.pack( side = TOP, fill=X,expand=False) 
         
         label_frame_fva = LabelFrame(labelframe6,text="Fluxes with a solution range")
         frame_fva_list=Frame(label_frame_fva)
         frame_fva_list.pack(side = TOP,fill=Y,expand=True)
         label_frame_fva.pack(side = LEFT,fill=Y,expand=True)
         scrollbar611 = Scrollbar(frame_fva_list)
         scrollbar611.pack( side = RIGHT, fill=Y,expand=True)
         scrollbar614 = Scrollbar(frame_fva_list,orient=HORIZONTAL)
         self.fva_listbox = Listbox(frame_fva_list, yscrollcommand = scrollbar611.set,xscrollcommand = scrollbar614.set ,width=self.average_width+20)
         self.fva_listbox.pack( side = TOP, fill=Y,expand=True)
         scrollbar614.pack( side = TOP, fill=X,expand=False) 
         scrollbar611.config( command = self.fva_listbox.yview ) 
         scrollbar614.config( command = self.fva_listbox.xview ) 
         #add_parameter_button = Button(label_frame_fva, text="Add selected flux as parameter",command=self.add_parameter)
         #add_parameter_button.pack( side = TOP)
         #add_parameter_button = Button(label_frame_fva, text="automatically\nadd parameters",command=self.add_all_parameters)
         #add_parameter_button.pack( side = TOP)
         
         """for flux in sorted(self.fva):
             if "RATIO_" in flux:
                 continue
             max_flux=round(fva[flux]["maximum"],2)
             min_flux=round(fva[flux]["minimum"],2)
             string=("%s < "+flux+" <%s ")%( min_flux, max_flux)
             fva_listbox.insert(END, string)"""
      
      
      def add_parameter(self):
             reaction_id=self.active_reaction_get_bounds
             reaction=self.label_model.constrained_model.reactions.get_by_id(reaction_id)
             if reaction_id not in self.label_model.parameter_dict:
                fva=flux_variability_analysis(self.label_model.constrained_model,fraction_of_optimum=self.fraction_of_optimum,reaction_list=[reaction_id],tolerance_feasibility=self.label_model.lp_tolerance_feasibility)
                minimum=fva[reaction_id]["minimum"] 
                maximum=fva[reaction_id]["maximum"]
                value=round((19*minimum+1*maximum)/20,5) #Weighted average
                model=self.label_model.constrained_model
                self.label_model.parameter_dict[reaction_id]={"v":value,"lb":model.reactions.get_by_id(reaction_id).lower_bound,"ub":model.reactions.get_by_id(reaction_id).upper_bound, "type":"flux value","reactions":[reaction_id] ,"max_d":0.1,"original_lb":model.reactions.get_by_id(reaction_id).lower_bound,"original_ub":model.reactions.get_by_id(reaction_id).upper_bound,"original_objective_coefficient":model.reactions.get_by_id(reaction_id).objective_coefficient}
                apply_parameters(self.label_model,parameter_precision=self.parameter_precision,parameter_list=[reaction_id])
                self.parameter_button.config(text="Remove from parameters")
             else:
                 clear_parameters(self.label_model,parameter_dict=None,parameter_list=[reaction_id], clear_ratios=False,clear_turnover=False,clear_fluxes=True) 
                 del self.label_model.parameter_dict[reaction_id]
                 self.parameter_button.config(text="Make parameter")
             reaction.objective_coefficient=0.0
             self.boundaries_lb_spinbox.delete(0,"end")
             self.boundaries_lb_spinbox.insert(INSERT,str(reaction.lower_bound))
             self.boundaries_ub_spinbox.delete(0,"end")
             self.boundaries_ub_spinbox.insert(INSERT,str(reaction.upper_bound))
             self.varoptimizacion.set(0)
             if reaction.id not in self.modified_reactions:
                         self.modified_reactions.append(reaction.id)
             self.highlight_modified(self.boundaries_reaction_selector_listbox,self.modified_reactions)
             self.run_fva()
      
      
      """def create_add_all_parameters_windows(self): 
          top=self.add_all_parameters_windows=Toplevel()  
          threshold_change_frame=Frame(top)
          threshold_change_frame.pack(side=TOP)
          Label(threshold_change_frame,text="Change threshold").pack(side=LEFT)
          self.threshold_change_entry=Entry(threshold_change_frame)
          self.threshold_change_entry.delete(0, END)
          self.threshold_change_entry.insert(0, 0.1) 
          self.threshold_change_entry.pack(side=LEFT)
          button01 = Button(top, text="Accept",command=self.add_all_parameters)
          button01.pack( side = TOP) """         
                   
      
      def add_all_parameters(self):
          #self.change_threshold=max(float( self.threshold_change_entry.get()),self.parameter_precision)
          parameters=identify_free_parameters(self.label_model,parameter_dict={},n_samples=self.label_model.p_dict["identify_free_parameters_n_samples"],add_to_model=True,parameter_precision=self.parameter_precision,change_threshold=self.change_threshold,fraction_of_optimum=self.fraction_of_optimum,max_d=0.1,key_reactions=[],add_turnover=False,excluded_turnovers=[],turnover_upper_bound=1000,debug=True)
          reaction=self.label_model.constrained_model.reactions.get_by_id(self.active_reaction_get_bounds)
          self.boundaries_lb_spinbox.delete(0,"end")
          self.boundaries_lb_spinbox.insert(INSERT,str(reaction.lower_bound))
          self.boundaries_ub_spinbox.delete(0,"end")
          self.boundaries_ub_spinbox.insert(INSERT,str(reaction.upper_bound))
          for parameter in parameters:
              if "flux value" in parameters[parameter]["type"]:
                  for reaction in parameters[parameter]["reactions"]:
                      if reaction not in self.modified_reactions:
                         self.modified_reactions.append(reaction)
          self.highlight_modified(self.boundaries_reaction_selector_listbox,self.modified_reactions)
          self.run_fva()
      
      
      def export_fluxes(self):
           fileName = tkFileDialog.asksaveasfilename(parent=self.root,title="Save fluxes as...",filetypes=[("xlsx","*.xlsx"),("csv","*.csv")])
           #if "xlsx" not in fileName.lower() and "csv" not in fileName.lower()  and "txt" not in fileName.lower() and "xlsm" not in fileName.lower() and "xltx" not in fileName.lower() and "xltm" :
           if not any(x in fileName.lower() for x in known_extensions): 
                fileName+=".xlsx"  
           if len(fileName)>0:
                write_fva(self.label_model.constrained_model, fn=fileName,fraction=self.fraction_of_optimum,remove0=False,change_threshold=10*self.parameter_precision, mode="full", lp_tolerance_feasibility=self.label_model.lp_tolerance_feasibility)
      
          
      def export_model(self):
           fileName = tkFileDialog.asksaveasfilename(parent=self.root,title="Save model with constrains as...",filetypes=[("sbml","*.sbml"),("xml","*.xml"),("xlsx","*.xlsx"),("csv","*.csv")])
           if ".sbml" in fileName.lower() or ".xml" in fileName.lower():
              cobra.io.write_sbml_model(self.label_model.constrained_model, fileName)
           else:
              convert_model_to_spreadsheet(self.label_model.constrained_model,fileName) 
      
      def export_label(self):
           fileName = tkFileDialog.asksaveasfilename(parent=self.root,title="Save label simulation results as...",filetypes=[("xlsx",".xlsx"),("csv",".csv")])
           if not any(x in fileName.lower() for x in known_extensions): 
                fileName+=".xlsx"  
           if len(fileName)>0:
              a1,b=solver(self.label_model,mode="fsolve",fba_mode=self.fba_mode.get().lower())
              a2,b,c=get_objective_function(self.label_model,force_balance=self.label_model.force_balance,output=True)
              export_label_results(self.label_model,fn=fileName,show_chi=True)
      """def add_file_output_interface(self):
        labelframe7 = LabelFrame(self.root, text="")
        labelframe7.pack(fill=None, expand=False)
        #button71 = Button(labelframe7, text="Export SBML", fg="red")
        #button71.pack( side = LEFT)
        button72 = Button(labelframe7, text="Export Excel",command=self.export_fluxes )
        button72.pack( side = LEFT )"""
           
      def label_button(self):
          if self.label_model.emu_dict!={}: #Create the button only if a label model exists
             button01 = Button(self.root, text="Simulate Label",command=self.create_figure_window)
             button01.pack( side = TOP)
          
      """def add_turnover_display(self):
          add_turnover_frame=LabelFrame(self.root,text="Set flux turnover")
          add_turnover_frame.pack(fill=None, expand=False) 
          turnover_selector = Frame(add_turnover_frame)
          turnover_selector.pack(side = "left", fill=Y,expand=True)
          turnover_selector_scrollbar = Scrollbar(turnover_selector)
          turnover_selector_scrollbar.pack( side = LEFT, fill=Y,expand=True)
          turnover_selector_listbox = Listbox(turnover_selector, yscrollcommand =turnover_selector_scrollbar.set , height = 3, width=14)
          turnover_selector_scrollbar['command'] = turnover_selector_listbox.yview
          turnover_selector_listbox['yscrollcommand'] = turnover_selector_scrollbar.set
          turnover_selector_listbox.pack( side = LEFT, fill = Y, expand=True )
          turnover_selector_scrollbar.config( command = turnover_selector_listbox.yview )
          for flux in sorted(self.label_model.turnover_flux_dict,key=lambda v: v.upper()):
             turnover_selector_listbox.insert(END, flux)
          self.selected_turnover=turnover_selector_listbox.get(0)
          self.turnover_selector_listbox=turnover_selector_listbox
          self.turnover_spinbox=Spinbox(add_turnover_frame, from_=0, to=1000, increment = 0.1, format="%.2f", width = 7)
          self.turnover_spinbox.delete(0,"end")
          value0=self.label_model.turnover_flux_dict[sorted(self.label_model.turnover_flux_dict,key=lambda v: v.upper())[0]]
          self.turnover_spinbox.insert(INSERT,value0)
          self.get_selected_turnover() 
          #self.flux_ratios_scale = Scale( self.sublabelframe_ratios_scale, resolution=0.01, from_= 1, to = 100, orient = HORIZONTAL)
          #self.flux_ratios_scale.set(1)
          self.selected_turnover_label=Label(turnover_selector,text=self.selected_turnover+" turnover",width=18)
          equal_label=Label(turnover_selector,text="=")
          self.selected_turnover_label.pack(side=LEFT)
          equal_label.pack(side=LEFT)
          self.turnover_spinbox.pack( side = "left")
          apply_button= Button(add_turnover_frame, text="Apply",command=self.update_turnover)
          apply_button.pack(side=TOP)
          resetbutton= Button(add_turnover_frame, text="Reset",command=self.reset_turnover)
          resetbutton.pack( side = TOP )
          add_turnover_as_parameter_button= Button(add_turnover_frame, text="Add as parameter",command=self.add_turnover_as_parameter)
          add_turnover_as_parameter_button.pack( side = TOP )
          #reset_all_button= Button(add_turnover_frame, text="Reset all",command=self.reset_all_turnover)
          #reset_all_button.pack( side = TOP )
          #add_all_turnover_as_parameter= Button(add_turnover_frame, text="Add all as parameters",command=self.add_all_turnover_as_parameter)
          #add_all_turnover_as_parameter.pack( side = TOP )"""
      
      def reset_turnover(self):
          flux=self.selected_turnover
          self.label_model.turnover_flux_dict[flux]=copy.deepcopy(self.backup_turnover_flux_dict[flux])
          if (flux+"_turnover") in self.label_model.parameter_dict:
              self.label_model.parameter_dict[flux+"_turnover"]["v"]=self.label_model.turnover_flux_dict[flux]["v"]
          self.turnover_spinbox.delete(0,"end")
          self.turnover_spinbox.insert(INSERT,self.backup_turnover_flux_dict[flux])
          #self.turnover_selector_listbox.itemconfig(self.selected_turnover, {'fg': 'blue'})
          """if flux in self.modified_turnovers:
             self.modified_turnovers.remove(flux)
          self.highlight_modified(self.turnover_selector_listbox,self.modified_turnovers)  """
                    
      def reset_all_turnover(self):
          #self.turnover_selector_listbox
          print "resiting turnovers"
          self.label_model.turnover_flux_dict=copy.copy(self.backup_turnover_flux_dict)
          #self.turnover_selector_listbox.delete(0, END)  
          for flux in sorted(self.label_model.turnover_flux_dict,key=lambda v: v.upper()):
              if flux+"_turnover" in self.label_model.parameter_dict:
                 self.label_model.parameter_dict[flux+"_turnover"]["v"]=self.label_model.turnover_flux_dict[flux]["v"]
          """#self.selected_turnover=0
          self.turnover_spinbox.delete(0,"end")
          value0=self.label_model.turnover_flux_dict[sorted(self.label_model.turnover_flux_dict,key=lambda v: v.upper())[0]]
          self.turnover_spinbox.insert(INSERT,value0)
          self.modified_turnovers=[]  
          self.highlight_modified(self.turnover_selector_listbox,self.modified_turnovers)"""  
      
      
      def add_turnover_as_parameter(self):
          turnover=self.active_reaction_get_bounds
          if turnover+"_turnover" not in self.label_model.parameter_dict:
             self.label_model.parameter_dict[turnover+"_turnover"]={"v":self.label_model.turnover_flux_dict[turnover]["v"],"lb":self.label_model.turnover_flux_dict[turnover]["lb"],"ub":self.label_model.turnover_flux_dict[turnover]["ub"],"max_d":0.1,"type":"turnover","reactions":[turnover]} 
             apply_parameters(self.label_model,parameter_precision=self.parameter_precision,parameter_list=[turnover+"_turnover"])
             self.add_parameter_turnover_button.config(text="Remove from parameter")
          else:
             del self.label_model.parameter_dict[turnover+"_turnover"]
             self.add_parameter_turnover_button.config(text="Make parameter")
      
      def add_all_turnover_as_parameter(self):
          #Automatically find a good the upper bound for turnovers
          turnover_upper_bound=0
          for reaction_id in self.label_model.reaction_n_dict:
              if reaction_id in self.fva:
                  turnover_upper_bound=max(self.fva[reaction_id]["maximum"],turnover_upper_bound)
          for turnover in self.label_model.turnover_flux_dict:
             self.label_model.parameter_dict[turnover+"_turnover"]={"v":self.label_model.turnover_flux_dict[turnover]["v"],"lb":self.label_model.turnover_flux_dict[turnover]["lb"],"ub":self.label_model.turnover_flux_dict[turnover]["ub"],"max_d":0.1,"type":"turnover","reactions":[turnover]} 
          apply_parameters(self.label_model,parameter_precision=self.parameter_precision,parameter_list=[turnover+"_turnover"])
          
           
      
      """def get_selected_turnover(self):
          self.get_selected_turnover_after=self.root.after(50, self.get_selected_turnover)
          previous_slected_turnover =self.selected_turnover
          selection_tupple=self.turnover_selector_listbox.curselection()
          if selection_tupple!=(): 
             position,=selection_tupple
             self.selected_turnover=self.turnover_selector_listbox.get(position)
             if previous_slected_turnover!=self.selected_turnover:
                #self.turnover_selector_listbox.itemconfig(self.selected_turnover, {'fg': 'blue'})
                #if  self.turnover_selector_listbox.get(previous_slected_turnover) in self.modified_turnovers:
                #    self.turnover_selector_listbox.itemconfig(previous_slected_turnover, {'fg': 'green'})
                #else:
                #    self.turnover_selector_listbox.itemconfig(previous_slected_turnover, {'fg': 'black'})
                turnover_id=self.selected_turnover
                turnover_value=self.label_model.turnover_flux_dict[turnover_id]["v"]  
                self.turnover_spinbox.delete(0,"end")
                self.turnover_spinbox.insert(INSERT,str(turnover_value))
                self.selected_turnover_label.config(text=self.selected_turnover+" turnover")
                print self.selected_turnover"""
      
      
      """def update_turnover(self):
          #self.root.after(200, self.update_turnover)
          selected_turnover=self.selected_turnover
          print selected_turnover
          try:
            new_value=float(self.turnover_spinbox.get())
          except:
            new_value=(self.turnover_spinbox.get())
            new_value=float(new_value.replace(",","."))
          if selected_turnover not in self.modified_turnovers:
            self.modified_turnovers.append(selected_turnover)
          self.label_model.turnover_flux_dict[selected_turnover]["v"]=new_value
          print self.label_model.turnover_flux_dict
          self.highlight_modified(self.turnover_selector_listbox,self.modified_turnovers)  
          #self.turnover_selector_listbox.itemconfig(self.selected_turnover, {'fg': 'green'})"""
            
      
      
      def create_fig(self,root,condition="glc"):
          a,b,simulation_not_rsm=get_objective_function(self.label_model,force_balance=self.label_model.force_balance,output=False,rsm="never")
          a,b,simulation_with_rsm=get_objective_function(self.label_model,force_balance=self.label_model.force_balance,output=False,rsm="always")
          #ventanafigura = Tk()
          figure_height=sum([30+25*len(self.label_model.experimental_dict[condition][emu]) for emu in self.label_model.experimental_dict[condition]])
          for emu in  self.label_model.rsm_list:
              if emu in self.label_model.experimental_dict[condition]:
                 figure_height+=25*len(self.label_model.experimental_dict[condition][emu])
          fig = Canvas(root, width=280, height=figure_height,scrollregion=(0,0,0, figure_height))
          pos=0
          for emu_id in self.label_model.experimental_dict[condition]:
             if emu_id in self.label_model.data_name_emu_dict:
                title=self.label_model.data_name_emu_dict[emu_id]
             else:
                title=emu_id
             base=60
             pos+=20
             fig.create_text(base, pos+10, text='0%')
             fig.create_text(base+50, pos, text=title)
             fig.create_text(base+(100*2), pos+10, text='100%')
             fig.create_line(base, pos+20, base+(100*2), pos+20)
             size=self.label_model.emu_dict[emu_id]["size"]
             pos+=10
             #self.emu_canvas_dict[condition][emu_id]=fig
             self.simulated_points_object_dict[condition][emu_id]={}
             
             for n in sorted(self.label_model.emu_dict[emu_id]["mid"]):
                #mid=label_model.emu_dict[emu_id]["mid"][n]
                if n not in  self.label_model.experimental_dict[condition][emu_id]:
                    continue
                """if n==0 and emu_id in self.label_model.rsm_list:
                    continue"""
                #print "a"
                pos+=20
                string="m"+str(n)
                mean=self.label_model.experimental_dict[condition][emu_id][n]["m"]*100 #Percentatge
               #print "b"
                sd=self.label_model.experimental_dict[condition][emu_id][n]["sd"]*100 #Percentatge
                #print "c"
                fig.create_text(base-30, pos+5, text=string)
                #print "d"
                fig.create_rectangle(base, pos, base+(mean*2), pos+10, fill="gray")
                fig.create_line(max(base+(mean*2)-(sd*2),base), pos+5, min(base+(mean*2)+(sd*2),base+210), pos+5)
                if n in simulation_not_rsm[condition][emu_id]:
                   #n_emu=label_model.size_variable_dict[size][mid]
                   #print "e"
                   value=round(simulation_not_rsm[condition][emu_id][n],2)
                   print value
                   oval_id=fig.create_oval((base+value*100*2)-3, pos+5-2, (base+value*100*2)+3, pos+5+2, fill="red")
                   self.simulated_points_object_dict[condition][emu_id][n]=[oval_id,pos] 
                   """for xx in range(0, len(argument2[x])-3): 
                   fig.create_oval((base+argument2[x][xx+3]*2)-3, pos+5-2, (base+argument2[x][xx+3]*2)+3, pos+5+2, fill="red")"""
             if emu_id in self.label_model.rsm_list:
                    base=60
                    pos+=20
                    fig.create_text(base, pos+10, text='0%')
                    fig.create_text(base+50, pos, text=title+"(m/Sm)")
                    fig.create_text(base+(100*2), pos+10, text='100%')
                    fig.create_line(base, pos+20, base+(100*2), pos+20)
                    size=self.label_model.emu_dict[emu_id]["size"]
                    pos+=10
                    #self.emu_canvas_dict[condition][emu_id]=fig
                    #self.simulated_points_object_dict[condition][emu_id]={}
                    
                    for n in sorted(self.label_model.emu_dict[emu_id]["mid"]):
                       #mid=label_model.emu_dict[emu_id]["mid"][n]
                       if n not in  self.label_model.experimental_dict[condition][emu_id]:
                           continue
                       if n==0 :
                           continue
                       #print "a"
                       pos+=20
                       string="m"+str(n)+"/Sm"
                       mean=self.label_model.experimental_dict[condition][emu_id][n]["m"]/(1-self.label_model.experimental_dict[condition][emu_id][0]["m"])*100
                       sd=self.label_model.experimental_dict[condition][emu_id][n]["sd"]/(1-self.label_model.experimental_dict[condition][emu_id][0]["m"])*100 #Percentatge
                       fig.create_text(base-30, pos+5, text=string)
                       #print "d"
                       fig.create_rectangle(base, pos, base+(mean*2), pos+10, fill="gray")
                       fig.create_line(max(base+(mean*2)-(sd*2),base), pos+5, min(base+(mean*2)+(sd*2),base+210), pos+5)
                       #fig.create_line(base+(mean*2)-(sd*2), pos+5, base+(mean*2)+(sd*2), pos+5)
                       if n in simulation_with_rsm[condition][emu_id]:
                          #n_emu=label_model.size_variable_dict[size][mid]
                          #print "e"
                          value=round(simulation_with_rsm[condition][emu_id][n],2)#round(self.label_model.condition_simulation_results_dict[condition][emu_id][n],2)
                          print value
                          value=max(0,value)
                          oval_id=fig.create_oval((base+value*100*2)-3, pos+5-2, (base+value*100*2)+3, pos+5+2, fill="red")
                          self.simulated_points_object_dict[condition][emu_id]["ms"+str(n)]=[oval_id,pos] 
                          """for xx in range(0, len(argument2[x])-3): 
                          fig.create_oval((base+argument2[x][xx+3]*2)-3, pos+5-2, (base+argument2[x][xx+3]*2)+3, pos+5+2, fill="red")"""
             
             pos+=10
          
          print self.simulated_points_object_dict
          return fig
      
      def update_label(self):
        #print self.label_model.flux_dict
        a,b,simulation_not_rsm=get_objective_function(self.label_model,force_balance=self.label_model.force_balance,output=False,rsm="never")
        a,b,simulation_with_rsm=get_objective_function(self.label_model,force_balance=self.label_model.force_balance,output=False,rsm="always")
        for condition in self.simulated_points_object_dict:
            base=60
            fig=self.condition_canvas_dict[condition]
            if self.label_window==True:
              for emu_id in self.label_model.experimental_dict[condition]:
                for n in sorted(self.simulated_points_object_dict[condition][emu_id]):
                    if "ms" not in str(n):
                        new_value=round(simulation_not_rsm[condition][emu_id][n],3)
                    else:
                        actual_n=int(n.replace("ms",""))
                        new_value=round(simulation_with_rsm[condition][emu_id][actual_n],4)
                    new_value=max(0,new_value)
                    oval_id=self.simulated_points_object_dict[condition][emu_id][n][0]
                    pos=self.simulated_points_object_dict[condition][emu_id][n][1]
                    fig.delete(oval_id)
                    oval_id=fig.create_oval((base+new_value*100*2)-3, pos+5-2, (base+new_value*100*2)+3, pos+5+2, fill="red")
                    self.simulated_points_object_dict[condition][emu_id][n]=[oval_id,pos]
                    #print emu_id+" updatded"+str(new_value)  
      
      
      
      def update_label_sampling(self):
        a,b,simulation_not_rsm=get_objective_function(self.label_model,force_balance=self.label_model.force_balance,output=False,rsm="never")
        a,b,simulation_with_rsm=get_objective_function(self.label_model,force_balance=self.label_model.force_balance,output=False,rsm="always")
        for condition in self.simulated_points_object_dict:
            base=60
            fig=self.condition_canvas_dict[condition]
            if self.label_window==True:
              for emu_id in self.label_model.experimental_dict[condition]:
                for n in sorted(self.simulated_points_object_dict[condition][emu_id]):
                    if "ms" not in str(n):
                        new_value=round(simulation_not_rsm[condition][emu_id][n],3)
                    else:
                        actual_n=int(n.replace("ms",""))
                        new_value=max(round(simulation_with_rsm[condition][emu_id][actual_n],3),0)
                    pos=self.simulated_points_object_dict[condition][emu_id][n][1]
                    oval_id=fig.create_oval((base+new_value*100*2)-3, pos+5-2, (base+new_value*100*2)+3, pos+5+2, fill="blue")
                    self.sampling_points_condition_dict[condition].append(oval_id)           
      
      
      def create_figure_window(self):
          if self.label_window==True:
             return 
          #TODO create a single canvas and make it scrolable
          self.condition_canvas_dict={}
          self.simulated_points_object_dict={}
          self.sampling_points_condition_dict={}
          a,b=solver(self.label_model,mode="fsolve",fba_mode=self.fba_mode.get().lower())
          #a,b,simulation_not_rsm=get_objective_function(self.label_model,force_balance=self.label_model.force_balance,output=False,rsm="never")
          #a,b,simulation_with_rsm=get_objective_function(self.label_model,force_balance=self.label_model.force_balance,output=False,rsm="always")
          print a
          self.windows_condtion_dict={}
          for condition in self.label_model.experimental_dict:
              self.simulated_points_object_dict[condition]={}
              self.sampling_points_condition_dict[condition]=[]
              top = Toplevel()
              top.title(condition)
              self.windows_condtion_dict[condition]=top
              nrow=0
              ncol=0
              fig=self.create_fig(top,condition)
              scrollbar=Scrollbar(top,orient=VERTICAL)
              scrollbar.pack(side=RIGHT,fill=Y)
              scrollbar.config(command=fig.yview)
              fig.pack(fill=Y) 
              fig.config(yscrollcommand=scrollbar.set)
              self.condition_canvas_dict[condition]=fig 
              top.protocol("WM_DELETE_WINDOW",self.close_all_label_windows) #Closes all windows if one is closed
              """for emu in self.label_model.condition_simulation_results_dict[condition]:
                  if nrow>4:
                     ncol+=1
                     nrow=0
                  fig=self.create_fig(self.label_model,top,emu_id=emu,condition=condition)
                  #fig.pack()
                  fig.grid(row=nrow,column=ncol) #Replace with grid
                  nrow+=1"""
          #self.update_label()
          self.label_window=True
          
      
      def close_all_label_windows(self):
          for condition in self.windows_condtion_dict:
              window=self.windows_condtion_dict[condition]  
              window.destroy()  
              self.label_window=False            
          
      
      def restore_original_parameters(self):
          self.label_model.constrained_model=copy.deepcopy(self.backup_constrained_model)
          #self.label_model.parameter_dict={}
          self.label_model.turnover_flux_dict=copy.deepcopy(self.backup_turnover_flux_dict)
          self.label_model.ratio_dict=copy.deepcopy(self.backup_ratio_dict)
          self.root.destroy()
          
      def create_view_parameters_window(self):
          self.view_parameters_window=top = Toplevel()
          top.title("View parameters")
          list_frame=Frame(top)
          list_frame.pack(side=TOP)
          self.parameters_scrollbar = Scrollbar(list_frame)
          self.parameters_xscrollbar = Scrollbar(list_frame,orient=HORIZONTAL)
          self.parameters_scrollbar.pack( side = RIGHT, fill=Y,expand=True)
          self.parameters_listbox = Listbox(list_frame, yscrollcommand = self.parameters_scrollbar.set,xscrollcommand = self.parameters_xscrollbar.set, width=self.average_width+20 )
          self.parameters_scrollbar.config( command = self.parameters_listbox.yview )
          self.parameters_xscrollbar.config( command = self.parameters_listbox.xview ) 
          self.parameters_listbox.pack( side = TOP, fill=Y,expand=True) 
          self.parameters_xscrollbar.pack( side = TOP, fill=X,expand=False) 
          self.precision=int(-1*(math.log10(self.parameter_precision)))
          for parameter in sorted(self.label_model.parameter_dict,key=lambda v: v.upper()):
              string=parameter+"="+str(round(self.label_model.parameter_dict[parameter]["v"],self.precision))
              self.parameters_listbox.insert(END,string)
          button = Button(top, text="Clear Parameters", command=self.clear_all_parameters)
          button.pack(side=TOP)
      
      def clear_all_parameters(self):
          self.fitted_parameters=False
          parameters=self.label_model.parameter_dict
          for parameter in parameters:
              if parameters[parameter]["type"]=="flux value":
                 for reaction_id in parameters[parameter]["reactions"]:
                     if reaction_id in self.modified_reactions:
                        self.modified_reactions.remove(reaction_id)
                        """self.boundaries_reaction_selector_listbox.itemconfig(self.inverse_llistafluxos_dict[reaction_id], {'fg': 'green'})"""
              elif parameters[parameter]["type"]=="turnover":
                 for reaction_id in parameters[parameter]["reactions"]:
                     if reaction_id in self.modified_turnovers:
                        self.modified_turnovers.remove(reaction_id)
          clear_parameters(self.label_model,parameter_dict=None,parameter_list=[], clear_ratios=False,clear_turnover=False,clear_fluxes=True) 
          self.label_model.parameter_dict={}
          self.run_fva() 
          self.view_parameters_window.destroy()  
          #self.highlight_modified(self.turnover_selector_listbox,self.modified_turnovers)  
          self.highlight_modified(self.boundaries_reaction_selector_listbox,self.modified_reactions)
          self.update_fva_listbox()
      
      def highlight_modified(self,listbox,modified_list):
          for n, reaction in enumerate(listbox.get(0, END)):
              print reaction
              if reaction in modified_list:
                 listbox.itemconfig(n,{'fg': 'green'})
              else:
                 listbox.itemconfig(n,{'fg': 'black'})
      
      
      def add_reaction_name_label(self):
          self.reaction_name_label=Label(self.root,text="Selected reaction:")
          self.reaction_name_label.pack(side=TOP,expand=False)
      
      def get_selected_reaction_name(self):
          self.get_selected_reaction_name_after=self.root.after(200, self.get_selected_reaction_name)
          reaction_id=None
          if self.boundaries_reaction_selector_listbox.curselection()!=():
             position,=self.boundaries_reaction_selector_listbox.curselection()
             reaction_id=self.boundaries_reaction_selector_listbox.get(position)
             print [reaction_id,position]
          elif  self.flux_ratio_selector1_listbox.curselection()!=():
                print self.flux_ratio_selector1_listbox.curselection
                position,=self.flux_ratio_selector1_listbox.curselection()
                reaction_id=self.flux_ratio_selector1_listbox.get(position)
                print [reaction_id,position]
          elif self.optimal_solution_listbox.curselection()!=():
               position,=self.optimal_solution_listbox.curselection()
               flux_string=self.optimal_solution_listbox.get(position)
               reaction_id=self.determined_flux_regular_expression.match(flux_string).group(1)
               print [reaction_id,position]
          elif self.fva_listbox.curselection()!=(): 
               position,=self.fva_listbox.curselection()
               flux_string=self.fva_listbox.get(position)
               reaction_id=self.free_flux_regular_expression.match(flux_string).group(1)
               print [reaction_id,position]
          if reaction_id!=None:
             self.selected_reaction_name=self.label_model.metabolic_model.reactions.get_by_id(reaction_id).name
             if len(self.selected_reaction_name)>65:
                self.selected_reaction_name=self.selected_reaction_name[:62]+"..."
          self.reaction_name_label.config(text="Selected reaction: "+self.selected_reaction_name)
      
      """def view_label_propagating_reactions(self):
          label_propagating_reactions=LabelFrame(self.root,text="Label propagating reactions")
          label_propagating_reactions.pack(fill=True, expand=True) 
          label_propagating_selector = Frame(add_turnover_frame)
          label_propagating_selector.pack(side = "left", fill=Y,expand=True)
          label_propagating_selector_scrollbar = Scrollbar(label_propagating_selector)
          label_propagating_selector_scrollbar.pack( side = LEFT, fill=Y,expand=True)
          label_propagating_selector_listbox = Listbox(label_propagating_selector, yscrollcommand =label_propagating_selector_scrollbar.set , height = 3, width=14)
          label_propagating_selector_scrollbar['command'] = label_propagating_selector_listbox.yview
          label_propagating_selector_listbox['yscrollcommand'] = label_propagating_selector_scrollbar.set
          label_propagating_selector_listbox.pack( side = LEFT, fill = Y, expand=True )
          label_propagating_selector_scrollbar.config( command = label_propagating_selector_listbox.yview )
          reactions_and_turnover_list=[] 

          turnover_selector_listbox.insert(END, flux)
          self.selected_turnover=turnover_selector_listbox.get(0)
          self.turnover_selector_listbox=turnover_selector_listbox
          self.turnover_spinbox=Spinbox(add_turnover_frame, from_=0, to=1000, increment = 0.1, format="%.2f", width = 7)
          self.turnover_spinbox.delete(0,"end")
          value0=self.label_model.turnover_flux_dict[sorted(self.label_model.turnover_flux_dict,key=lambda v: v.upper())[0]]
          self.turnover_spinbox.insert(INSERT,value0)
          self.get_selected_turnover() 
          #self.flux_ratios_scale = Scale( self.sublabelframe_ratios_scale, resolution=0.01, from_= 1, to = 100, orient = HORIZONTAL)
          #self.flux_ratios_scale.set(1)
          self.selected_turnover_label=Label(turnover_selector,text=self.selected_turnover+" turnover",width=18)
          equal_label=Label(turnover_selector,text="=")
          self.selected_turnover_label.pack(side=LEFT)
          equal_label.pack(side=LEFT)
          self.turnover_spinbox.pack( side = "left")
          apply_button= Button(add_turnover_frame, text="Apply",command=self.update_turnover)
          apply_button.pack(side=TOP)
          resetbutton= Button(add_turnover_frame, text="Reset",command=self.reset_turnover)
          resetbutton.pack( side = TOP )
          add_turnover_as_parameter_button= Button(add_turnover_frame, text="Add as parameter",command=self.add_turnover_as_parameter)
          add_turnover_as_parameter_button.pack( side = TOP )
          #reset_all_button= Button(add_turnover_frame, text="Reset all",command=self.reset_all_turnover)
          #reset_all_button.pack( side = TOP )
          #add_all_turnover_as_parameter= Button(add_turnover_frame, text="Add all as parameters",command=self.add_all_turnover_as_parameter)
          #add_all_turnover_as_parameter.pack( side = TOP )
      
      def create_label_propagating_reactions.listbox(self,listbox):
          reactions_and_turnover_dict=[]
          for label_propagating_reaction in self.label_model.label_propagating_reactions:
              minimum=round(self.fva[label_propagating_reaction]["minimum"],3)
              maximum=round(self.fva[label_propagating_reaction]["maximum"],3)
              reactions_and_turnover_dict[label_propagating_reaction]=str(minimum)+"< "+label_propagating_reaction+" <"+str(maximum)
          for turnover in label_model.turnover_flux_dict:
                  
                            
          for flux in sorted(self.label_model.turnover_flux_dict,key=lambda v: v.upper()):"""
          
          
      def __init__(self,root,label_model): 
          self.ratio_regular_expression = re.compile('(.+)=(.*)') #Regular expression used to obtain ratio from
          self.free_flux_regular_expression=re.compile(".+< (.+) <.+")
          self.determined_flux_regular_expression=re.compile("(.+)=.+")
          self.root = root
          self.root.protocol("WM_DELETE_WINDOW",self.restore_original_parameters)
          self.label_model=label_model
          self.parameter_precision=label_model.parameter_precision
          self.precision=int(-1*(math.log10(self.parameter_precision)))
          print self.label_model.p_dict
          self.change_threshold=self.label_model.p_dict["identify_free_parameters_change_threshold"]
          root.title(self.label_model.project_name)
          self.log_name=self.label_model.project_name[:-9]+"_log.txt"
          self.label_model.constrained_model=label_model.constrained_model
          self.backup_constrained_model=copy.deepcopy(label_model.constrained_model) #Used to restore bounds and objectives
          self.backup_turnover_flux_dict=copy.deepcopy(label_model.turnover_flux_dict) #Used to restore bounds and objectives
          self.fba_mode = StringVar()
          self.fba_mode.set("FBA")
          self.backup_ratio_dict=copy.deepcopy(label_model.ratio_dict)
          self.modified_reactions=[]
          self.modified_turnovers=[]
          self.selected_r_ratio1=""
          self.selected_r_ratio2=""
          self.selected_turnover=0
          self.rsm="never" #"Never" dont display m/Sm, always: always display m/SM
          #self.variable13b1=0
          #self.variable13b2=0
          #self.variable31b1_a=0
          #self.variable31b1_b=1
          #self.variable42b1=0
          self.selected_ratio=-1
          self.waiting_for_fva=False
          self.waiting_for_label_simulation=True
          self.reseting_bounds=False
          self.fraction_of_optimum=self.label_model.p_dict["fraction_of_optimum"]
          self.force_bounds_update=False 
          self.waiting_for_parameters=False
          self.fitted_parameters=True
          self.llistafluxos=[]
          #Default annealing parameters
          self.annealing_n=self.label_model.p_dict["annealing_n"]
          self.annealing_m=self.label_model.p_dict["annealing_m"]
          self.annealing_p0=self.label_model.p_dict["annealing_p0"]
          self.annealing_pf=self.label_model.p_dict["annealing_pf"]
          self.annealing_n_processes=self.label_model.p_dict["annealing_n_processes"]
          self.annealing_n_iterations=self.label_model.p_dict["annealing_iterations"]
          self.relative_max_random_sample=self.label_model.p_dict["annealing_relative_max_sample"]
          self.relative_min_random_sample=self.label_model.p_dict["annealing_relative_min_sample"]
          #Build reaction name_dict:
          self.id_name_dict={}
          for reaction in self.label_model.metabolic_model.reactions:
              if reaction.name==None:
                 reaction.name=reaction.id 
              self.id_name_dict[reaction.id]=reaction.name
          self.average_width=int(numpy.mean([len(x) for x in self.id_name_dict ])) 
          """for reaction in label_model.metabolic_model.reactions:
                 self.llistafluxos.append(reaction.id)
          self.llistafluxos2=self.llistafluxos=sorted(self.llistafluxos,key=lambda v: v.upper())
          self.llistafluxos_dict={}
          self.inverse_llistafluxos_dict={}
          for n,x in enumerate(self.llistafluxos):
             self.llistafluxos_dict[n]=x
             self.inverse_llistafluxos_dict[x]=n"""      
          self.add_menubar()
          self.add_search_bar()
          #self.add_filters_bar()
          #self.input_frame()
          self.add_boundaries()
          self.add_ratios_interface()
          self.add_fraction_optimum_interface()
          #self.add_turnover_display()
          #add_turnover_frame=Frame(root,width=768, height=576)
          #add_turnover_frame.pack(fill="y", expand="yes") 
          self.add_flux_display()
          #self.add_file_output_interface() 
          #poll13b2()
          self.get_bounds_objective() 
          #self.update_bounds()
          #poll42b1()
          #self.get_r1_r2()
          self.update_ratio_list()
          self.run_fva()
          #self.update_turnover()  
          self.label_button()
          self.add_reaction_name_label()
          self.selected_reaction_name=""
          self.get_selected_reaction_name()
          self.get_ratio_selection()
          self.active_search=""
          self.active_filter=0
          self.update_search()
 
      

def launch_gui(label_model):
    """print 
    if __name__ == '__main__':"""
    root = Tk()
    gui = GUI(root,label_model)
    root.mainloop()

class build_model_gui:
    def get_file(self,root,text="",default_file=None,open_file_widget=None,row=0):
        file_frame=root  
        #file_frame.grid(side=TOP)
        Label(file_frame,text=text).grid(row=row,column=0)
        file_entry=Entry(file_frame,width=60)
        file_entry.grid(row=row,column=1) 
        if default_file!=None:
           file_entry.insert(0,default_file)     
        button = Button(file_frame,text="Browse",command=open_file_widget)
        button.grid(row=row,column=2)
        return file_entry
    
    def get_working_directory(self):
        wd=tkFileDialog.askdirectory(title="Select working directory:",parent=self.root)  
        os.chdir(wd)
        self.wd_entry.delete(0, END)
        self.wd_entry.insert(0,wd) 
       
    def get_sbml_model(self):
        loaded_file = tkFileDialog.askopenfile(title='Choose a constraint based model',filetypes=[("sbml","*.sbml"),("xml","*.xml"),("xlsx","*.xlsx"),("CSV","*.csv"),("all files","*")],parent=self.root)
        self.sbml_entry.delete(0, END) 
        self.sbml_entry.insert(0,str(loaded_file.name)) 
    
    def get_label_propagation_rules(self):
        loaded_file = tkFileDialog.askopenfile(title='Choose a set of label propagation rules',filetypes=[("xlsx","*.xlsx"),("csv","*.csv"),('All files','*.*')],parent=self.root)
        self.label_rules_entry.delete(0, END) 
        self.label_rules_entry.insert(0,str(loaded_file.name))
        
    def get_exp_data(self):
        loaded_files = tkFileDialog.askopenfiles(title='Choose experimental isotopologues data',filetypes=[("xlsx","*.xlsx"),("csv","*.csv"),('All files','*.*')],parent=self.root)
        file_names=[x.name for x in loaded_files]
        self.e_data_entry.delete(0, END) 
        for file_name in file_names:
           self.e_data_entry.insert(0,file_name+",")
        self.e_data_entry.delete(len(self.e_data_entry.get())-1) #delete last coma
        
    def get_flux_constraints(self):
        loaded_file = tkFileDialog.askopenfile(title='Chose constraints file',filetypes=[("xlsx","*.xlsx"),("csv","*.csv"),('All files','*.*')],parent=self.root)
        self.constraints_entry.delete(0, END) 
        self.constraints_entry.insert(0,str(loaded_file.name))
        
    def get_settings(self):
        loaded_file = tkFileDialog.askopenfile(title='Chose Settings file',filetypes=[("xlsx","*.xlsx"),("csv","*.csv"),('All files','*.*')],parent=self.root)
        self.settings_entry.delete(0, END) 
        self.settings_entry.insert(0,str(loaded_file.name))
    def new_model(self):
        self.store_inputs()
        if ""==self.sbml_entry.get():# or ""==self.e_data_entry.get or ""==self.label_rules_entry.get():
             print "A constrained model must be defined"
             return
        p_dict={'reactions_with_forced_turnover': [], 'annealing_cycle_time_limit': 1800, 'confidence_max_absolute_perturbation': 10, 'turnover_exclude_EX': True, 'annealing_n_processes': 3, 'annealing_p0': 0.4, 'identify_free_parameters_add_turnover': True, 'minimum_sd': 0.01, 'annealing_max_perturbation': 1, 'turnover_upper_bound': 100, 'confidence_perturbation': 0.1, 'annealing_m': 1000, 'annealing_n': 20, 'annealing_relative_max_sample': 0.3, 'confidence_min_absolute_perturbation': 0.05, 'annealing_pf': 0.0001, 'confidence_significance': 0.95, 'identify_free_parameters_change_threshold': 0.001, 'parameter_precision': 0.0001, 'fraction_of_optimum': 1, 'lp_tolerance_feasibility': 1e-09, 'identify_free_parameters_n_samples': 400, 'annealing_relative_min_sample': 0.1, 'annealing_iterations': 3,"gene_expression_mode":"imat", "gene_expression_low_expression_threshold":25,"gene_expression_high_expression_threshold":75,"gene_expression_percentile":True,"gene_expression_gene_method":"avearge", "gene_expression_gene_sufix":"_AT","gene_expression_gene_prefix":"","gene_expression_epsilon":1, "gene_expression_lex_epsilon":1e-6,"gene_expression_fraction_optimum":1, "gene_expression_absent_gene_expression_value":50}
        #p_dict={'reactions_with_forced_turnover': [], 'annealing_cycle_time_limit': 1800, 'confidence_max_absolute_perturbation': 10, 'turnover_exclude_EX': True, 'annealing_n_processes': 4, 'annealing_p0': 0.4, 'identify_free_parameters_add_turnover': True, 'minimum_sd': 0.01, 'annealing_max_perturbation': 1, 'turnover_upper_bound': 100, 'confidence_perturbation': 0.1, 'annealing_m': 1000, 'annealing_n': 10, 'annealing_relative_max_sample': 0.35, 'confidence_min_absolute_perturbation': 0.05, 'annealing_pf': 0.0001, 'confidence_significance': 0.95, 'identify_free_parameters_change_threshold': 0.005, 'parameter_precision': 0.0001, 'fraction_of_optimum': 0, 'lp_tolerance_feasibility': 1e-09, 'identify_free_parameters_n_samples': 200, 'annealing_relative_min_sample': 0.2, 'annealing_iterations': 2,"gene_prefix":"gene","gene_sufix":"_AT"}
        model_file=self.sbml_entry.get()
        if ".sbml" in model_file.lower() or ".xml" in model_file.lower():
            model=cobra.io.read_sbml_model(self.sbml_entry.get())
        else:
            model=create_cobra_model_from_file(model_file) 
        try:
          read_flux_constraints(model,ratio_dict={},file_name=self.constraints_entry.get(),create_copies=False)
        except:
          print "Constraints loaded"
          pass
        try: 
          p_dict.update(read_isotoflux_settings(settings_entry.get()))
          print "Settings loaded"
        except:
          pass 
        self.project_name = tkFileDialog.asksaveasfilename(title="Save project as...",filetypes=[("iso2flux",".iso2flux")])
        if ".iso2flux" not in self.project_name:
            self.project_name+=".iso2flux"
        log_name=self.project_name[:-9]+"_log.txt"
        f=open(log_name, "w")
        f.write(time.strftime("%c")+"//"+"Instance created")
        self.label_model=label_model=Label_model(model,lp_tolerance_feasibility=p_dict["lp_tolerance_feasibility"],parameter_precision=p_dict["parameter_precision"],reactions_with_forced_turnover=p_dict["reactions_with_forced_turnover"])
        self.label_model.p_dict=p_dict
        self.label_model.eqn_dir=self.project_name[:-9]+"_equations"
        try:  
          read_isotopomer_model(label_model,self.label_rules_entry.get())
          missing_reactions_list= find_missing_reactions(label_model).keys()
          if len(missing_reactions_list)>0:
             raise Exception ("Some reactions lack label propagation rules: "+str(missing_reactions_list)) 
          #loaded_file = tkFileDialog.askopenfile(title='Choose experimental measuments file',filetypes=[("xlsx",".xlsx")]) 
          #fileName=loaded_file.name
          e_data_names=self.e_data_entry.get().replace("[","").replace("]","").split(",")
          emu_dict0,label_model.experimental_dict =read_experimental_mid(label_model,e_data_names,emu0_dict={},experimental_dict={},minimum_sd=p_dict["minimum_sd"])
          label_model.build(emu_dict0,force_balance=True,recompile_c_code=True,remove_impossible_emus=True,isotopic_steady_state=True,excluded_outputs_inputs=[],turnover_upper_bound=p_dict["turnover_upper_bound"],clear_intermediate_data=True,turnover_exclude_EX=p_dict['turnover_exclude_EX'])
          self.finish_create_new_model() 
        except Exception, e:
           print e 
           top=Toplevel()
           input_missing_label=Label(top,text="Warning:\nIso2flux was unable to build the label propagation\nmodel due to missing or wrong inputs.\nThis will prevent simulating 13C propagation.\n Do you wish to continue?")
           input_missing_label.pack(side=TOP)
           yes_button=Button(top,text="Yes",command=self.finish_create_new_model) 
           yes_button.pack(side=TOP)
           no_button=Button(top,text="No",command=top.destroy) 
           no_button.pack(side=TOP)
           label_model.constrained_model=copy.deepcopy(label_model.metabolic_model)
        #save_iso2flux_model(label_model,name=self.project_name,write_sbml=True,gui=False)
        #self.root.destroy()
    
    def finish_create_new_model(self):
        save_iso2flux_model(self.label_model,name=self.project_name,write_sbml=True,gui=False)
        self.root.destroy()  
        
    def validate_model(self):
        e=None
        steady_state_flag=False
        missing_dict={}
        self.root.withdraw() 
        self.store_inputs()
        if ""==self.sbml_entry.get() or ""==self.e_data_entry.get or ""==self.label_rules_entry.get():
             print "Mandatory Input missing"
             return
        with open("validation_results.txt","w") as f:
                 f.write("Validation summary:\n") 
        p_dict={'reactions_with_forced_turnover': [], 'annealing_cycle_time_limit': 1800, 'confidence_max_absolute_perturbation': 10, 'turnover_exclude_EX': True, 'annealing_n_processes': 4, 'annealing_p0': 0.4, 'identify_free_parameters_add_turnover': True, 'minimum_sd': 0.01, 'annealing_max_perturbation': 1, 'turnover_upper_bound': 100, 'confidence_perturbation': 0.1, 'annealing_m': 1000, 'annealing_n': 10, 'annealing_relative_max_sample': 0.35, 'confidence_min_absolute_perturbation': 0.05, 'annealing_pf': 0.0001, 'confidence_significance': 0.95, 'identify_free_parameters_change_threshold': 0.005, 'parameter_precision': 0.0001, 'fraction_of_optimum': 1.0, 'lp_tolerance_feasibility': 1e-09, 'identify_free_parameters_n_samples': 200, 'annealing_relative_min_sample': 0.2, 'annealing_iterations': 2}
        model_file=self.sbml_entry.get()
        if ".sbml" in model_file.lower() or ".xml" in model_file.lower():
            model=cobra.io.read_sbml_model(self.sbml_entry.get())
        else:
            model=create_cobra_model_from_file(model_file) 
        try:
          read_flux_constraints(model,ratio_dict={},file_name=self.constraints_entry.get(),create_copies=False)
        except:
          print "Constraints loaded"
          pass
        try: 
          p_dict.update(read_isotoflux_settings(settings_entry.get()))
          print "Settings loaded"
        except:
          pass
        try:
          if model.optimize().status=="infeasible":
             feasible=False
             with open("validation_results.txt","a") as f:
                 f.write("Error:Constraint based model is infeasible with the current constraints"+"\n") 
          else:
             feasible=True    
             self.label_model=label_model=Label_model(model,lp_tolerance_feasibility=p_dict["lp_tolerance_feasibility"],parameter_precision=p_dict["parameter_precision"],reactions_with_forced_turnover=p_dict["reactions_with_forced_turnover"])
             self.label_model.project_name="validation"
             read_isotopomer_model(label_model,self.label_rules_entry.get())
             missing_dict=find_missing_reactions(label_model,fn="validation_results.txt",fn_mode="a")
             #print self.label_model.flux_dict
             if len(missing_dict)==0:
              #loaded_file = tkFileDialog.askopenfile(title='Choose experimental measuments file',filetypes=[("xlsx",".xlsx")]) 
              #fileName=loaded_file.name
              e_data_names=self.e_data_entry.get().replace("[","").replace("]","").split(",")
              emu_dict0,label_model.experimental_dict =read_experimental_mid(label_model,e_data_names,emu0_dict={},experimental_dict={},minimum_sd=p_dict["minimum_sd"])
              label_model.build(emu_dict0,force_balance=False,recompile_c_code=True,remove_impossible_emus=True,isotopic_steady_state=True,excluded_outputs_inputs=[],turnover_upper_bound=p_dict["turnover_upper_bound"],clear_intermediate_data=False,turnover_exclude_EX=p_dict['turnover_exclude_EX'])
              steady_state_flag=check_steady_state(label_model,only_initial_m0=True,threshold=1e-9,fn="validation_results.txt",fn_mode="a")#Check initial dy for steady state deviations
               
              #check_simulated_fractions(label_model)
              #self.label_model.p_dict=p_dict
              #save_iso2flux_model(label_model,name="project",write_sbml=True,gui=True)
              #self.root.destroy()
             else:
               steady_state_flag=False
        except Exception, e:
            print e
            with open("validation_results.txt","a") as f:
                 f.write(str(e)+"\n")       
        
        top=Toplevel()
        top.protocol("WM_DELETE_WINDOW",self.restart_script)
        if steady_state_flag and len(missing_dict)==0 and e==None and feasible==True:
           validation_summary=Label(top,text="Succesfully validated\n")
        else:  
           validation_summary=Label(top,text="Validation failed,\nsee validation_results.txt for details\n")
           #print self.label_model.flux_dict
        validation_summary.pack(side=TOP)
        ok_button=Button(top,text="Ok",command=self.restart_script) 
        ok_button.pack(side=TOP)
        
    def restart_script(self):
       os.execv(sys.executable, ['python'] + sys.argv)    
    def store_inputs(self):
        with open("inputs.temp", 'w') as fp:
             json.dump({"sbml":self.sbml_entry.get(),"label_rules":self.label_rules_entry.get(),"e_data":self.e_data_entry.get(),"constraints":self.constraints_entry.get(),"settings":self.settings_entry.get()}, fp)
    def load_inputs(self):
        try:
          with open("inputs.temp", 'r') as fp:
                  loaded_inputs=json.load(fp)
          self.sbml_entry.delete(0, END) 
          self.sbml_entry.insert(0,str(loaded_inputs["sbml"]))
          self.label_rules_entry.delete(0, END) 
          self.label_rules_entry.insert(0,str(loaded_inputs["label_rules"]))
          self.e_data_entry.delete(0, END) 
          self.e_data_entry.insert(0,str(loaded_inputs["e_data"]))
          self.constraints_entry.delete(0, END) 
          self.constraints_entry.insert(0,str(loaded_inputs["constraints"]))
          self.settings_entry.delete(0, END) 
          self.settings_entry.insert(0,str(loaded_inputs["settings"]))
        except:
          pass
        
    def __init__(self,root):
        self.root=root
        root.title("Iso2Flux")
        main_frame=Frame(root)
        #self.wd_entry=self.get_file(root,text="Working directory",default_file=os.getcwd(),open_file_widget=self.get_working_directory,row=)
        self.sbml_entry=self.get_file(root,text="Constraint based model",default_file=None,open_file_widget=self.get_sbml_model,row=0)
        self.label_rules_entry=self.get_file(root,text="Label propagation rules",default_file=None,open_file_widget=self.get_label_propagation_rules,row=1)
        self.e_data_entry=self.get_file(root,text="13C labelling patterns",default_file=None,open_file_widget=self.get_exp_data,row=2)
        self.constraints_entry=self.get_file(root,text="Flux constraints (Optional)",default_file=None,open_file_widget=self.get_flux_constraints,row=3)
        self.settings_entry=self.get_file(root,text="Advanced settings(Optional)",default_file=None,open_file_widget=self.get_settings,row=4)
        button1=Button(root,text="Create Iso2flux Instance",command=self.new_model)
        button2=Button(root,text="Validate Inputs",command=self.validate_model)
        button1.grid(row=6,column=1)
        button2.grid(row=7,column=1)
        self.load_inputs()


"""def launch_build_model_gui():
    tk.destroy()
    global label_model
    root = Tk()
    gui = build_model_gui(root)
    root.mainloop()
    label_model=gui.label_model


def load_model():
    tk.withdraw()
    global label_model
    label_model=load_iso2flux_model(gui=True)
    tk.destroy()"""


