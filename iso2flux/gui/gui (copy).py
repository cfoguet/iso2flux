from Tkinter import *
import copy
import re
import json
import tkFileDialog
import math
import numpy

#Ilabel imports
from cobra.flux_analysis.variability import flux_variability_analysis
from ..flux_functions.apply_ratios import apply_ratios, remove_ratio #Functions to apply ratios of fluxes
from ..flux_functions.define_reaction_group import define_reaction_group


from ..emus_functions.solver import solver #To do, move it to emu equations

from ..output_functions.model_to_excel import model_to_excel
from ..output_functions.showr import showr
from ..output_functions.write_fva import write_fva
from ..output_functions.print_emu_results import  print_emu_results

#from ilabel.fitting.anneal import annealing
from ..fitting.get_objective_function import get_objective_function
from ..fitting.apply_parameters import apply_parameters
from ..fitting.save_parameters import save_parameters
from ..fitting.load_parameters import load_parameters
from ..fitting.load_parameters import load_parameters
from ..fitting.identify_parameters import identify_parameters
from ..fitting.sampling import sampling
from ..fitting.anneal import annealing
from ..fitting.clear_parameters import clear_parameters
from ..misc.expression_analysis import integrate_omics

from ..misc.minimal_flux import add_flux_limit_constraints

from cobra.flux_analysis.parsimonious import optimize_minimal_flux


class GUI:
      label_window=False
      ratio_inverse_name_dict={}
      def add_menubar(self):
          menubar = Menu(self.root)
          filemenu = Menu(menubar, tearoff=0)
          filemenu.add_command(label="Load parameters", command=self.load_parameters)
          filemenu.add_command(label="Save parameters", command=self.save_parameters)
          menubar.add_cascade(label="File", menu=filemenu)
          fittingmenu = Menu(menubar, tearoff=0)
          #fittingmenu.add_command(label="Fit fluxes", command=self.create_fitting_window)
          fittingmenu.add_command(label="Fit fluxes", command=self.create_fitting_window)
          fittingmenu.add_command(label="Sampling", command=self.create_sampling_window)
          fittingmenu.add_command(label="View parameters", command=self.create_view_parameters_window)
          menubar.add_cascade(label="Fitting", menu=fittingmenu)
          options = Menu(menubar, tearoff=0)
           
          add_constrain_menu = Menu(menubar, tearoff=0)
          add_constrain_menu.add_command(label="Restrict max fluxes", command=self.create_limit_max_fluxes_window)
          add_constrain_menu.add_command(label="Integrate gene expression", command=self.create_integrate_gene_expression_window)
          menubar.add_cascade(label="Add constraints", menu=add_constrain_menu)
          
          #menubar.add_cascade(label="Options", menu=options)
          self.root.config(menu=menubar)
      
      def create_limit_max_fluxes_window(self):
          top=self.limit_max_fluxes_window= Toplevel() 
          fraction_of_flux_minimum_frame=Frame(top)
          fraction_of_flux_minimum_frame.pack(side=TOP)
          Label(fraction_of_flux_minimum_frame,text="relative max flux").pack(side=LEFT)
          self.fraction_of_flux_minimum_entry=Entry(fraction_of_flux_minimum_frame)
          self.fraction_of_flux_minimum_entry.delete(0, END)
          self.fraction_of_flux_minimum_entry.insert(0, 1) 
          self.fraction_of_flux_minimum_entry.pack(side=LEFT)
          button = Button(top, text="Add contsraints",command=self.limit_max_fluxes)
          button.pack( side = TOP)
          
          
      def get_gene_expression_file(self):
          loaded_file = tkFileDialog.askopenfile(title='Choose a file',filetypes=[("xlsx",".xlsx")]) 
          file_name=loaded_file.name    
          self.gene_expression_file_entry.delete(0, END)
          self.gene_expression_file_entry.insert(0, file_name) 
                
      
      def create_integrate_gene_expression_window(self):
          top=self.integrate_gene_expression_window= Toplevel() 
          low_expression_threshold_frame=Frame(top)
          low_expression_threshold_frame.pack(side=TOP)
          Label(low_expression_threshold_frame,text="low expression percentile").pack(side=LEFT)
          self.low_expression_threshold_entry=Entry(low_expression_threshold_frame)
          self.low_expression_threshold_entry.delete(0, END)
          self.low_expression_threshold_entry.insert(0, 25) 
          self.low_expression_threshold_entry.pack(side=LEFT)
          
          high_expression_threshold_frame=Frame(top)
          high_expression_threshold_frame.pack(side=TOP)
          Label(high_expression_threshold_frame,text="high expression percentile").pack(side=LEFT)
          self.high_expression_threshold_entry=Entry(high_expression_threshold_frame)
          self.high_expression_threshold_entry.delete(0, END)
          self.high_expression_threshold_entry.insert(0, 75) 
          self.high_expression_threshold_entry.pack(side=LEFT)
          
          hex_epsilon_frame=Frame(top)
          hex_epsilon_frame.pack(side=TOP)
          Label(hex_epsilon_frame,text="hex epsilon").pack(side=LEFT)
          self.hex_epsilon_entry=Entry(hex_epsilon_frame)
          self.hex_epsilon_entry.delete(0, END)
          self.hex_epsilon_entry.insert(0, 1) 
          self.hex_epsilon_entry.pack(side=LEFT)
          
          lex_epsilon_frame=Frame(top)
          lex_epsilon_frame.pack(side=TOP)
          Label(lex_epsilon_frame,text="lex epsilon").pack(side=LEFT)
          self.lex_epsilon_entry=Entry(lex_epsilon_frame)
          self.lex_epsilon_entry.delete(0, END)
          self.lex_epsilon_entry.insert(0, 0.001) 
          self.lex_epsilon_entry.pack(side=LEFT)
                    
          gene_expression_file_frame=Frame(top)
          gene_expression_file_frame.pack(side=TOP)
          Label(gene_expression_file_frame,text="gene_expression_file").pack(side=LEFT)
          self.gene_expression_file_entry=Entry(gene_expression_file_frame)
          self.gene_expression_file_entry.delete(0, END)
          self.gene_expression_file_entry.insert(0, "") 
          self.gene_expression_file_entry.pack(side=LEFT)
          button2 = Button(gene_expression_file_frame, text="Open",command=self.get_gene_expression_file)
          button2.pack( side = LEFT)  
          
          button2 = Button(top, text="Constrain using gene expression",command=self.integrate_gene_expression)
          button2.pack( side = TOP)
          
      def integrate_gene_expression(self):
           low_expression_threshold=float(self.low_expression_threshold_entry.get())
           high_expression_threshold=float(self.high_expression_threshold_entry.get())
           hex_epsilon=float(self.hex_epsilon_entry.get())
           lex_epsilon=float(self.lex_epsilon_entry.get())
           excel_name=self.gene_expression_file_entry.get()
           integrate_omics(self.constrained_model,excel_name,fraction_of_optimum=self.fraction_of_optimum,low_expression_threshold=low_expression_threshold, high_expression_threshold=high_expression_threshold,percentile=True,gene_method="average",metabolite_id_list=[], epsilon=hex_epsilon,lex_epsilon=lex_epsilon,imat_fraction_optimum=1,maxupper_bound=1000,label_model=None,add_as_constraints=True)
           
           
           reaction_id=self.boundaries_reaction_selector_listbox.get(self.variable13b1)
           reaction=self.constrained_model.reactions.get_by_id(reaction_id)
           self.boundaries_lb_spinbox.delete(0,"end")
           self.boundaries_lb_spinbox.insert(INSERT,str(reaction.lower_bound))
           self.boundaries_ub_spinbox.delete(0,"end")
           self.boundaries_ub_spinbox.insert(INSERT,str(reaction.upper_bound))
           print self.constrained_model.optimize()
           self.integrate_gene_expression_window.destroy()
           self.run_fva()  
          
      def limit_max_fluxes(self):
           print "constraining fluxes"
           fraction_of_flux_minimum=float(self.fraction_of_flux_minimum_entry.get())
           add_flux_limit_constraints(self.constrained_model,fraction_of_optimum_objective=self.fraction_of_optimum, fraction_of_flux_minimum=fraction_of_flux_minimum,solver="cplex")
           self.limit_max_fluxes_window.destroy()
           reaction_id=self.boundaries_reaction_selector_listbox.get(self.variable13b1)
           reaction=self.constrained_model.reactions.get_by_id(reaction_id)
           self.boundaries_lb_spinbox.delete(0,"end")
           self.boundaries_lb_spinbox.insert(INSERT,str(reaction.lower_bound))
           self.boundaries_ub_spinbox.delete(0,"end")
           self.boundaries_ub_spinbox.insert(INSERT,str(reaction.upper_bound))
           print self.constrained_model.optimize()
           self.run_fva()
                    
      def create_sampling_window(self):
          top=self.sampling_window= Toplevel()
          nframe=Frame(top)
          nframe.pack(side=TOP)
          Label(nframe,text="n").pack(side=LEFT)
          self.sampling_n_entry=Entry(nframe)
          self.sampling_n_entry.delete(0, END)
          self.sampling_n_entry.insert(0, 100) 
          self.sampling_n_entry.pack(side=LEFT)
          sampling_max_turnoverframe=Frame(top)
          sampling_max_turnoverframe.pack(side=TOP)
          Label(sampling_max_turnoverframe,text="max turnover").pack(side=LEFT)
          self.sampling_max_turnover_entry=Entry(sampling_max_turnoverframe)
          self.sampling_max_turnover_entry.delete(0, END)
          self.sampling_max_turnover_entry.insert(0, 1000) 
          self.sampling_max_turnover_entry.pack(side=LEFT)  
          
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
          max_turnover=float(self.sampling_max_turnover_entry.get())
          self.sampling_window.destroy()
          #sampling(self.label_model,n=n,fraction_of_optimum=self.fraction_of_optimum,output_emu_list=None,max_turnover=50,fba_mode="fba",parameter_precision= sel,gui=self,change_threshold=0.1)
          sampling(self.label_model,n=n,fraction_of_optimum=self.fraction_of_optimum,output_emu_list=None,max_turnover=max_turnover,fba_mode=self.fba_mode.get(),parameter_precision= self.parameter_precision,gui=self,change_threshold=self.change_threshold)
          
      
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
          
          
      def create_fitting_window(self):
          top=self.fitting_window= Toplevel()
          nframe=Frame(top)
          nframe.pack(side=TOP)
          Label(nframe,text="n cycles").pack(side=LEFT)
          self.nframeentry=Entry(nframe)
          self.nframeentry.delete(0, END)
          self.nframeentry.insert(0, 50) 
          self.nframeentry.pack(side=LEFT)
          mframe=Frame(top)
          mframe.pack(side=TOP)
          Label(mframe,text="simulations per cycle").pack(side=LEFT)
          self.mframeentry=Entry(mframe)
          self.mframeentry.delete(0, END)
          self.mframeentry.insert(0, 50) 
          self.mframeentry.pack(side=LEFT)
          p0frame=Frame(top) 
          p0frame.pack(side=TOP)
          Label(p0frame,text="p0").pack(side=LEFT)
          self.p0frameentry=Entry(p0frame)
          self.p0frameentry.delete(0, END)
          self.p0frameentry.insert(0, 0.4) 
          self.p0frameentry.pack(side=LEFT)
          pfframe=Frame(top) 
          pfframe.pack(side=TOP)
          Label(pfframe,text="pf").pack(side=LEFT)
          self.pfframeentry=Entry(pfframe)
          self.pfframeentry.delete(0, END)
          self.pfframeentry.insert(0, 0.01) 
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
          Label(max_perturbationframe,text="max parameter perturbation").pack(side=LEFT)
          self.max_perturbationframeentry=Entry(max_perturbationframe)
          self.max_perturbationframeentry.delete(0, END)
          self.max_perturbationframeentry.insert(0, 1) 
          self.max_perturbationframeentry.pack(side=LEFT)
          
          n_process_frame=Frame(top) 
          n_process_frame.pack(side=TOP)
          Label(n_process_frame,text="number of paralel processes").pack(side=LEFT)
          self.n_process_entry=Entry(n_process_frame)
          self.n_process_entry.delete(0, END)
          self.n_process_entry.insert(0, 1) 
          self.n_process_entry.pack(side=LEFT)
          
          Label(top,text="n parameters is "+str(len(self.label_model.parameter_dict))).pack(side=TOP)
          
          max_sampleframe=Frame(top) 
          max_sampleframe.pack(side=TOP)
          Label(max_sampleframe,text="max size of parameter sample").pack(side=LEFT)
          self.max_sampleframeentry=Entry(max_sampleframe)
          self.max_sampleframeentry.delete(0, END)
          self.max_sampleframeentry.insert(0, len(self.label_model.parameter_dict)) 
          self.max_sampleframeentry.pack(side=LEFT)
          
          min_sampleframe=Frame(top) 
          min_sampleframe.pack(side=TOP)
          Label(min_sampleframe,text="min size of parameter sample").pack(side=LEFT)
          self.min_sampleframeentry=Entry(min_sampleframe)
          self.min_sampleframeentry.delete(0, END)
          self.min_sampleframeentry.insert(0, len(self.label_model.parameter_dict)) 
          self.min_sampleframeentry.pack(side=LEFT)
          
          
          
          
          button01 = Button(top, text="Start",command=self.fitting)
          button01.pack( side = TOP)
          
      def create_waiting_window(self):
          self.waiting= Toplevel()
          self.waiting.title("Fitting in progress")
          frame=Frame(self.waiting)
          frame.pack() 
          Label(frame,text="\n\nPlease wait parameter fitting in progress...\n\n").pack(side=LEFT)
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
          max_d=float(self.max_dframeentry.get())
          p0=float(self.p0frameentry.get())
          pf=float(self.pfframeentry.get())
          n=int(self.nframeentry.get())
          m=int(self.mframeentry.get())
          fraction_of_optimum=float(self.fraction_of_optimum)
          max_perturbation=float(self.max_perturbationframeentry.get())
          max_sample=min(int(self.max_sampleframeentry.get()),len(self.label_model.parameter_dict))
          min_sample=min(int(self.min_sampleframeentry.get()),len(self.label_model.parameter_dict))
          n_processes=int(self.n_process_entry.get())
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
          print [max_d,p0,pf,n,m,max_perturbation]
          self.fitting_window.destroy()
          best_parameters,best_flux_dict, f_best=annealing(self.label_model,max_random_sample=max_sample,min_random_sample=min_sample,n=n,m=m,p0=p0,pf=pf,check_feasability=True,parameter_precision=self.parameter_precision,max_perturbation=max_perturbation,fraction_of_optimum=fraction_of_optimum,gui=self,n_processes=n_processes)
          parameters=self.label_model.parameter_dict=best_parameters
          for parameter in parameters:
              if parameters[parameter]["type"]=="flux value":
                 for reaction_id in parameters[parameter]["reactions"]:
                     if reaction_id not in self.modified_reactions:
                        self.modified_reactions.append(reaction_id)
                        """self.boundaries_reaction_selector_listbox.itemconfig(self.inverse_llistafluxos_dict[reaction_id], {'fg': 'green'})"""
              elif parameters[parameter]["type"]=="turnover":
                 for reaction_id in parameters[parameter]["reactions"]:
                     if reaction_id not in self.modified_turnovers:
                        self.modified_turnovers.append(reaction_id)
          for n,flux in enumerate(sorted(self.label_model.turnover_flux_dict,key=lambda v: v.upper())):
              if flux in self.modified_turnovers:
                 self.turnover_selector_listbox.itemconfig(n, {'fg': 'green'})
          self.get_bounds_objective() 
          self.run_fva()
          self.waiting.destroy()
      
      def load_parameters(self):
          loaded_file = tkFileDialog.askopenfile(title='Choose a file') 
          fileName=loaded_file.name
          if len(fileName)>0:
             with open(fileName, 'r') as fp:
                  loaded_data=json.load(fp)
             self.label_model.ratio_dict=loaded_data[1]
             self.label_model.turnover_flux_dict=loaded_data[2]
             self.update_ratio_list()
             apply_ratios(self.constrained_model,self.label_model.ratio_dict)
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
             self.force_bounds_update=False
             
      
      def save_parameters(self):
          fileName = tkFileDialog.asksaveasfilename(parent=self.root,title="Save the parameters as...")
          if len(fileName)>0:
              reaction_dict={}
              for reaction in self.constrained_model.reactions:
                  if "RATIO_" in reaction.id:
                     continue
                  reaction_dict[reaction.id]={"lb":reaction.lower_bound,"ub":reaction.upper_bound,"obj":reaction.objective_coefficient}
              store_paremeters=[reaction_dict,self.label_model.ratio_dict,self.label_model.turnover_flux_dict,self.label_model.parameter_dict]  
              with open(fileName, 'w') as fp:
                   json.dump(store_paremeters, fp)
          print fileName
      
      
      
      def get_ratio_selection(self):
          previous_selected_ratio=self.selected_ratio
          self.root.after(200, self.get_ratio_selection)
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
          reaction_id=self.boundaries_reaction_selector_listbox.get(self.variable13b1)
          reaction=self.constrained_model.reactions.get_by_id(reaction_id)
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
                 self.boundaries_reaction_selector_listbox.itemconfig(self.variable13b1, {'fg': 'green'}) 
          reaction.lower_bound=new_reaction_lb
          reaction.lower_bound=new_reaction_lb
          reaction.objective_coefficient=self.varoptimizacion.get()
          print [reaction.id, reaction.objective_coefficient]
          self.waiting_for_fva=True
          self.run_fva()
          
      def reset_bound(self):
          reaction_id=self.boundaries_reaction_selector_listbox.get(self.variable13b1)
          orginal_reaction=self.backup_constrained_model.reactions.get_by_id(reaction_id)
          lb=orginal_reaction.lower_bound
          ub=orginal_reaction.upper_bound
          obj=orginal_reaction.objective_coefficient
          reaction_to_reset=self.constrained_model.reactions.get_by_id(reaction_id)
          reaction_to_reset.lower_bound=lb
          reaction_to_reset.uppee_bound=ub
          reaction_to_reset.objective_coefficient=obj
          if reaction_id in self.label_model.parameter_dict:
             del self.label_model.parameter_dict[reaction_to_reset_id]
          
          self.boundaries_lb_spinbox.delete(0,"end")
          self.boundaries_lb_spinbox.insert(INSERT,str(lb))
          self.boundaries_ub_spinbox.delete(0,"end")
          self.boundaries_ub_spinbox.insert(INSERT,str(ub))
          self.varoptimizacion.set(int(obj))
          self.run_fva()
          if reaction_to_reset.id in self.label_model.parameter_dict:
             del self.label_model.parameter_dict[reaction_to_reset.id]
          if reaction_to_reset.id in self.modified_reactions:
             self.modified_reactions.remove(reaction_to_reset.id)
          search_term=self.boundaries_reaction_search_entry.get()
          self.boundaries_reaction_selector_listbox.delete(0, END)
          reaction_list=[x.id for x in self.label_model.metabolic_model.reactions.query(search_term)]
          for n,reaction in enumerate(sorted(reaction_list,key=lambda v: v.upper())):
              self.boundaries_reaction_selector_listbox.insert(END,reaction)
              if reaction in self.modified_reactions:
                  self.boundaries_reaction_selector_listbox.itemconfig(n, {'fg': 'green'})
          
          
      def reset_all_bounds(self):
          for reaction_id in self.modified_reactions:
              if reaction_id not in self.constrained_model.reactions:
                 continue
              #self.boundaries_reaction_selector_listbox.itemconfig(self.inverse_llistafluxos_dict[reaction_id], {'fg': 'black'}) 
              reaction=self.backup_constrained_model.reactions.get_by_id(reaction_id)
              lb=reaction.lower_bound
              ub=reaction.upper_bound
              obj=reaction.objective_coefficient
              reaction_to_reset=self.constrained_model.reactions.get_by_id(reaction_id)
              reaction_to_reset.lower_bound=lb
              reaction_to_reset.upper_bound=ub
              reaction_to_reset.objective_coefficient=obj
              if reaction_to_reset.id==self.boundaries_reaction_selector_listbox.get(self.variable13b1):
                 self.boundaries_lb_spinbox.delete(0,"end")
                 self.boundaries_lb_spinbox.insert(INSERT,str(lb))
                 self.boundaries_ub_spinbox.delete(0,"end")
                 self.boundaries_ub_spinbox.insert(INSERT,str(ub))
                 self.varoptimizacion.set(int(obj))
          self.run_fva()
          self.modified_reactions=[]
          self.label_model.parameter_dict={}
          search_term=self.boundaries_reaction_search_entry.get()
          reaction_list=[x.id for x in self.label_model.metabolic_model.reactions.query(search_term)]
          self.boundaries_reaction_selector_listbox.delete(0, END)
          for n,reaction in enumerate(sorted(reaction_list,key=lambda v: v.upper())):
              self.boundaries_reaction_selector_listbox.insert(END,reaction)
              if reaction in self.modified_reactions:
                  self.boundaries_reaction_selector_listbox.itemconfig(n, {'fg': 'green'})
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
          apply_ratios(self.constrained_model,self.label_model.ratio_dict)
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
          remove_ratio(self.constrained_model,ratio_id,self.label_model.ratio_dict)
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
               #convert_model_to_irreversible(self.constrained_model)
               optimize_minimal_flux(self.constrained_model, already_irreversible=False)
               #revert_model_to_reversible(self.constrained_model, update_solution=True)
            else:
               self.constrained_model.optimize()
            print self.constrained_model.solution
            optimal_solution=self.constrained_model.solution
            self.solution_dict=optimal_solution.x_dict
            self.fva=flux_variability_analysis(self.constrained_model,fraction_of_optimum=self.fraction_of_optimum)
            print "FVA updated"
          except:
            self.fva={}
            print "solution not optimal"
          self.update_fva_listbox()
       
      def update_fva_listbox(self):
          search_term=self.flux_display_search.get()
          reaction_list=[x.id for x in self.label_model.metabolic_model.reactions.query(search_term)]   
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
              max_flux=round(self.fva[flux]["maximum"],4)
              min_flux=round(self.fva[flux]["minimum"],4)
              """if max_flux==min_flux:
                 string=("    "+flux+"=%s ")%( min_flux)
              else:
                 string=("%s < "+flux+" <%s ")%( min_flux, max_flux)"""
              if (max_flux-min_flux)>=10*self.parameter_precision and flux not in self.label_model.parameter_dict:
                 string=("%s < "+flux+" <%s ")%( min_flux, max_flux)
                 self.fva_listbox.insert(END, string)
                 self.free_fluxes_n_dict[flux]=n
                 self.n_free_fluxes_dict[n]=flux
                 n+=1
              else:
                  optimal_flux=round(self.solution_dict[flux],5)
                  string=(flux+"=%s ")%( optimal_flux) 
                  self.optimal_solution_listbox.insert(END, string)
         
      
      def run_label_simulation(self):
          self.waiting_for_label_simulation=False 
          if self.label_window==True:
             a1,b=solver(self.label_model,mode="fsolve",fba_mode=self.fba_mode.get().lower())
             a2,b,c=get_objective_function(self.label_model,force_balance=self.label_model.force_balance,output=False)
             print a1,a2
             self.update_label()
           
      
      
      
      def get_bounds_objective(self): #Assign bounds and objective based on reaction
          #print variable13b1
          #active_reaction=self.variable13b1
          self.root.after(100, self.get_bounds_objective)
          selection_tupple=self.boundaries_reaction_selector_listbox.curselection()
          if selection_tupple!=(): 
             self.variable13b1,=selection_tupple
             new_active_reaction=self.boundaries_reaction_selector_listbox.get(self.variable13b1)
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
             reaction=self.constrained_model.reactions.get_by_id(self.active_reaction_get_bounds)
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
          text=(self.active_reaction_get_bounds)
          self.sublabelframe_boundaries_rid .config(text= text)
          
          #lab13b1.config(text = sel1)
          #print (boundaries_lb_spinbox.get())
          
          
      def poll13b2(self):
          
          selection_tupple=self.boundaries_reaction_selector_listbox.curselection()
          if selection_tupple!=(): 
             self.variable13b2,=selection_tupple
          #print "poll13b2"
          self.root.after(200, self.poll13b2)
          sel12 = " del " + (self.boundaries_reaction_selector_listbox.get(self.variable13b1)) + "; "
          #lab13b2.config(text = sel12)
      
      
      
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
          entry_search=Entry(reaction_selector,width=14)
          entry_search.pack(side=TOP)
          reaction_sub_selector.pack(side=TOP,expand="yes")
          reaction_selector_scrollbar = Scrollbar(reaction_sub_selector)
          reaction_selector_scrollbar.pack( side = RIGHT, fill=Y,expand="yes")
          reaction_selector_listbox = Listbox(reaction_sub_selector, yscrollcommand =reaction_selector_scrollbar.set , height = 3, width=14)
          reaction_selector_scrollbar['command'] = reaction_selector_listbox.yview
          reaction_selector_listbox['yscrollcommand'] = reaction_selector_scrollbar.set
          reaction_selector_listbox.pack( side = RIGHT, fill = BOTH ,expand="yes")
          reaction_selector_scrollbar.config( command = reaction_selector_listbox.yview )
          reaction_list=[x.id for x in self.label_model.metabolic_model.reactions]
          for reaction in sorted(reaction_list,key=lambda v: v.upper()):
                 reaction_selector_listbox.insert(END,reaction)
          return(reaction_selector,reaction_selector_listbox,entry_search) 
      
      def add_boundaries(self):
         "Flux Boundaries and objective"
         labelframe_flux_boundaries_master = LabelFrame(self.root, text="Set reaction bounds and optimization goal")
         labelframe_flux_boundaries_master.pack(fill="y", expand=False)
         labelframe_boundaries = LabelFrame(labelframe_flux_boundaries_master, text="")
         labelframe_boundaries.pack(fill="y", expand="yes")
         
         labelframe_boundaries_buttons = Frame(labelframe_boundaries)
         apply_button=Button(labelframe_boundaries_buttons, text="Apply", command = self.update_bounds)
         reset_bounds_button=Button(labelframe_boundaries_buttons, text="Reset", command = self.reset_bound)
         reset_all_bounds_button = Button(labelframe_boundaries_buttons, text="Reset All", command = self.reset_all_bounds)
         
         
         apply_button.pack(side = TOP)
         reset_bounds_button.pack(side=TOP)
         reset_all_bounds_button.pack( side = TOP)
         
         sublabelframe_boundaries_lb = Frame(labelframe_boundaries)
         self.boundaries_lb_spinbox = Spinbox(sublabelframe_boundaries_lb, from_=-1000, to=1000, increment = 0.01, format="%.2f", width = 7)
         reaction_selector,self.boundaries_reaction_selector_listbox,self.boundaries_reaction_search_entry=self.reaction_selector(labelframe_boundaries)
         self.active_reaction_get_bounds=self.boundaries_reaction_selector_listbox.get(0)
         self.boundaries_reaction_active_search=""
         self.sublabelframe_boundaries_rid = Label(labelframe_boundaries,width=10)
         self.sublabelframe_boundaries_rid.config(text ="")
         sublabelframe_boundaries_ub = Frame(labelframe_boundaries)
         #Label(sublabelframe_boundaries_ub, text="<").pack( side = LEFT)
         self.boundaries_ub_spinbox = DoubleVar()
         self.boundaries_ub_spinbox = Spinbox(sublabelframe_boundaries_ub, from_=-1000, to=1000, increment = 0.01, format="%.2f", width = 7)
         upper_bound_initial_reaction=self.constrained_model.reactions.get_by_id(self.boundaries_reaction_selector_listbox.get(0)).upper_bound
         self.boundaries_ub_spinbox.delete(0,"end")
         self.boundaries_ub_spinbox.insert(INSERT,upper_bound_initial_reaction)
         self.boundaries_ub_spinbox.pack( side = LEFT)
         
         lower_bound_initial_reaction=self.constrained_model.reactions.get_by_id(self.boundaries_reaction_selector_listbox.get(0)).lower_bound
         self.boundaries_lb_spinbox.delete(0,"end")
         self.boundaries_lb_spinbox.insert(INSERT,lower_bound_initial_reaction)
         self.boundaries_lb_spinbox.pack( side = LEFT)
         
         sublabelframe_optimization = LabelFrame(labelframe_boundaries,text="Objective")
         
         self.varoptimizacion = IntVar()
         self.varoptimizacion.set(0)
         radiobutton411 = Radiobutton(sublabelframe_optimization, text="none", variable=self.varoptimizacion, value=0)
         radiobutton411.pack( side = TOP)
         radiobutton412 = Radiobutton(sublabelframe_optimization, text="minimize", variable=self.varoptimizacion, value=-1)
         radiobutton412.pack( side = TOP)
         radiobutton413 = Radiobutton(sublabelframe_optimization, text="maximize", variable=self.varoptimizacion, value=1)
         radiobutton413.pack( side = TOP)
         
         greater_label1=Label(labelframe_boundaries, text=" > ")
         greater_label2=Label(labelframe_boundaries, text=" > ")
         
         reaction_selector.pack(side = LEFT,fill=Y,expand=True)
         Label(sublabelframe_boundaries_lb, text="").pack( side = LEFT)
         sublabelframe_boundaries_lb.pack(side = LEFT)
         greater_label1.pack(side = LEFT)
         self.sublabelframe_boundaries_rid.pack(side=LEFT)
         greater_label2.pack(side = LEFT)
         sublabelframe_boundaries_ub.pack(side = LEFT)
         sublabelframe_optimization.pack(side=LEFT)
         labelframe_boundaries_buttons.pack(side = LEFT)
         
      
      def search_function(self):
          self.root.after(200,self.search_function)
          ###Reaction bounds
          search_term=self.boundaries_reaction_search_entry.get()
          if self.boundaries_reaction_active_search!=search_term:
             self.boundaries_reaction_active_search=search_term
             self.boundaries_reaction_selector_listbox.delete(0, END)
             reaction_list=[x.id for x in self.label_model.metabolic_model.reactions.query(search_term)]
             for n,reaction in enumerate(sorted(reaction_list,key=lambda v: v.upper())):
                 self.boundaries_reaction_selector_listbox.insert(END,reaction)
                 if reaction in self.modified_reactions:
                    self.boundaries_reaction_selector_listbox.itemconfig(n, {'fg': 'green'})
          search_term_ratio=self.ratio_entry1.get()
          if self.search_term_active_ratio!= search_term_ratio:
             self.search_term_active_ratio=search_term_ratio
             self.flux_ratio_selector1_listbox.delete(0, END)
             reaction_list=[x.id for x in self.label_model.metabolic_model.reactions.query(search_term_ratio)]
             for reaction in sorted(reaction_list,key=lambda v: v.upper()):
                 self.flux_ratio_selector1_listbox.insert(END,reaction)
          search_term_flux=self.flux_display_search.get()
          if search_term_flux!= self.search_term_flux_active:
             self.search_term_flux_active=search_term_flux
             self.update_fva_listbox()
      
      def add_ratios_interface(self):
        root=self.root
        labelframe_fluxratios = LabelFrame(root, text="Relative values (R1/R2)")
        labelframe_fluxratios.pack(fill=Y, expand=False)
        sublabelframe_flux_ratio_selector1,self.flux_ratio_selector1_listbox,self.ratio_entry1=self.reaction_selector(labelframe_fluxratios)  
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
        ratio_selector_scrollbar.pack( side = RIGHT, fill=Y)
        self.ratio_selector_listbox = Listbox(subframe_ratio_selector, yscrollcommand = ratio_selector_scrollbar.set , height = 3, width=14)
        ratio_selector_scrollbar['command'] = self.ratio_selector_listbox.yview
        self.ratio_selector_listbox['yscrollcommand'] = ratio_selector_scrollbar.set
        self.ratio_selector_listbox.pack( side = LEFT, fill=Y, expand=True )
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
            self.fraction_optimum_spinbox.insert(INSERT,"1.0")
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
         self.flux_display_search=Entry(labelframe6,width=14)
         self.flux_display_search.pack(side=TOP)
         self.search_term_flux_active=""
         labelframe6.pack(fill=Y, expand=False)
         label_frame_determined_fluxes = LabelFrame(labelframe6,text="Uniquelly determined fluxes")
         label_frame_determined_fluxes.pack(side = LEFT,fill=Y,expand=True)
         scrollbar612 = Scrollbar(label_frame_determined_fluxes)
         scrollbar612.pack( side = RIGHT, fill=Y,expand=True)
         self.optimal_solution_listbox = Listbox(label_frame_determined_fluxes, yscrollcommand = scrollbar612.set )
         scrollbar612.config( command = self.optimal_solution_listbox.yview ) 
         self.optimal_solution_listbox.pack( side = RIGHT, fill=Y,expand=True) 
         
         label_frame_fva = LabelFrame(labelframe6,text="Free Fluxes")
         frame_fva_list=Frame(label_frame_fva)
         frame_fva_list.pack(side = TOP,fill=Y,expand=True)
         label_frame_fva.pack(side = LEFT,fill=Y,expand=True)
         scrollbar611 = Scrollbar(frame_fva_list)
         scrollbar611.pack( side = RIGHT, fill=Y,expand=True)
         self.fva_listbox = Listbox(frame_fva_list, yscrollcommand = scrollbar611.set )
         self.fva_listbox.pack( side = TOP, fill=Y,expand=True)
         scrollbar611.config( command = self.fva_listbox.yview ) 
         add_parameter_button = Button(label_frame_fva, text="add as\nparameter",command=self.add_parameter)
         add_parameter_button.pack( side = TOP)
         add_parameter_button = Button(label_frame_fva, text="automatically\nadd parameters",command=self.add_all_parameters)
         add_parameter_button.pack( side = TOP)
         
         
         """for flux in sorted(self.fva):
             if "RATIO_" in flux:
                 continue
             max_flux=round(fva[flux]["maximum"],2)
             min_flux=round(fva[flux]["minimum"],2)
             string=("%s < "+flux+" <%s ")%( min_flux, max_flux)
             fva_listbox.insert(END, string)"""
      
      
      def add_parameter(self):
          model=self.constrained_model
          selection_tupple=self.fva_listbox.curselection()
          if selection_tupple!=(): 
             n=selection_tupple[0]
             reaction_id=self.n_free_fluxes_dict[n]
             fva=flux_variability_analysis(model,fraction_of_optimum=self.fraction_of_optimum,reaction_list=[reaction_id])
             reaction=model.reactions.get_by_id(reaction_id)
             minimum=fva[reaction_id]["minimum"] 
             maximum=fva[reaction_id]["maximum"]
             value=round((19*minimum+1*maximum)/20,5) #Weighted average
             self.label_model.parameter_dict[reaction_id]={"v":value,"lb":model.reactions.get_by_id(reaction_id).lower_bound,"ub":model.reactions.get_by_id(reaction_id).upper_bound, "type":"flux value","reactions":[reaction_id] ,"max_d":0.1,"original_lb":model.reactions.get_by_id(reaction_id).lower_bound,"original_ub":model.reactions.get_by_id(reaction_id).upper_bound,"original_objective_coefficient":model.reactions.get_by_id(reaction_id).objective_coefficient}
             apply_parameters(self.label_model,parameter_precision=self.parameter_precision,parameter_list=[reaction_id])
             model.reactions.get_by_id(reaction_id).objective_coefficient=0.0
             if reaction_id==self.boundaries_reaction_selector_listbox.get(self.variable13b1):
                 self.boundaries_lb_spinbox.delete(0,"end")
                 self.boundaries_lb_spinbox.insert(INSERT,str(reaction.lower_bound))
                 self.boundaries_ub_spinbox.delete(0,"end")
                 self.boundaries_ub_spinbox.insert(INSERT,str(reaction.upper_bound))
                 self.varoptimizacion.set(0)
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
          parameters=identify_parameters( self.label_model,add_to_parameter_dict=True ,fraction_of_optimum=self.fraction_of_optimum ,parameter_precision=self.parameter_precision ,change_threshold=self.change_threshold,max_d=0.1,add_turnover=False)
          reaction=self.constrained_model.reactions.get_by_id(self.boundaries_reaction_selector_listbox.get(self.variable13b1))
          self.boundaries_lb_spinbox.delete(0,"end")
          self.boundaries_lb_spinbox.insert(INSERT,str(reaction.lower_bound))
          self.boundaries_ub_spinbox.delete(0,"end")
          self.boundaries_ub_spinbox.insert(INSERT,str(reaction.upper_bound))
          self.run_fva()
      
      
      def export_fluxes(self):
          write_fva(self.constrained_model,fn="flux_distribution",fraction=self.fraction_of_optimum,remove0=False,change_threshold=10*self.parameter_precision,mode="reduced")
          
            
      def add_file_output_interface(self):
        labelframe7 = LabelFrame(self.root, text="")
        labelframe7.pack(fill=None, expand=False)
        #button71 = Button(labelframe7, text="Export SBML", fg="red")
        #button71.pack( side = LEFT)
        button72 = Button(labelframe7, text="Export Excel",command=self.export_fluxes )
        button72.pack( side = LEFT )
           
      def label_button(self):
          button01 = Button(self.root, text="Simulate Label",command=self.create_figure_window)
          button01.pack( side = TOP)
          
      def add_turnover_display(self):
          add_turnover_frame=LabelFrame(self.root,text="Set reaction turnover")
          add_turnover_frame.pack(fill=None, expand=False) 
          turnover_selector = Frame(add_turnover_frame)
          turnover_selector.pack(side = "left", fill=Y,expand=True)
          turnover_selector_scrollbar = Scrollbar(turnover_selector)
          turnover_selector_scrollbar.pack( side = RIGHT, fill=Y,expand=True)
          turnover_selector_listbox = Listbox(turnover_selector, yscrollcommand =turnover_selector_scrollbar.set , height = 3, width=14)
          turnover_selector_scrollbar['command'] = turnover_selector_listbox.yview
          turnover_selector_listbox['yscrollcommand'] = turnover_selector_scrollbar.set
          turnover_selector_listbox.pack( side = RIGHT, fill = Y, expand=True )
          turnover_selector_scrollbar.config( command = turnover_selector_listbox.yview )
          for flux in sorted(self.label_model.turnover_flux_dict,key=lambda v: v.upper()):
             turnover_selector_listbox.insert(END, flux)
          self.turnover_selector_listbox=turnover_selector_listbox
          self.turnover_spinbox=Spinbox(add_turnover_frame, from_=0, to=1000, increment = 0.1, format="%.2f", width = 7)
          self.turnover_spinbox.delete(0,"end")
          value0=self.label_model.turnover_flux_dict[sorted(self.label_model.turnover_flux_dict,key=lambda v: v.upper())[0]]
          self.turnover_spinbox.insert(INSERT,value0)
          self.get_selected_turnover() 
          #self.flux_ratios_scale = Scale( self.sublabelframe_ratios_scale, resolution=0.01, from_= 1, to = 100, orient = HORIZONTAL)
          #self.flux_ratios_scale.set(1)
          self.turnover_spinbox.pack( side = "left")
          resetbutton= Button(add_turnover_frame, text="Reset",command=self.reset_turnover)
          resetbutton.pack( side = TOP )
          add_turnover_as_parameter_button= Button(add_turnover_frame, text="Add as parameter",command=self.add_turnover_as_parameter)
          add_turnover_as_parameter_button.pack( side = TOP )
          reset_all_button= Button(add_turnover_frame, text="Reset all",command=self.reset_all_turnover)
          reset_all_button.pack( side = TOP )
          add_all_turnover_as_parameter= Button(add_turnover_frame, text="Add all as parameters",command=self.add_all_turnover_as_parameter)
          add_all_turnover_as_parameter.pack( side = TOP )
      
      def reset_turnover(self):
          flux=self.turnover_selector_listbox.get(self.selected_turnover)
          self.label_model.turnover_flux_dict[flux]=copy.copy(self.backup_turnover_flux_dict[flux])
          if (flux+"_turnover") in self.label_model.parameter_dict:
              self.label_model.parameter_dict[flux+"_turnover"]["v"]=self.label_model.turnover_flux_dict[flux]
          self.turnover_spinbox.delete(0,"end")
          self.turnover_spinbox.insert(INSERT,self.backup_turnover_flux_dict[flux])
          self.turnover_selector_listbox.itemconfig(self.selected_turnover, {'fg': 'blue'})
          if flux in self.modified_turnovers:
             self.modified_turnovers.remove(flux)  
                    
      def reset_all_turnover(self):
          self.turnover_selector_listbox
          print "resiting turnovers"
          self.label_model.turnover_flux_dict=copy.copy(self.backup_turnover_flux_dict)
          self.turnover_selector_listbox.delete(0, END)  
          for flux in sorted(self.label_model.turnover_flux_dict,key=lambda v: v.upper()):
              self.turnover_selector_listbox.insert(END, flux)
              if flux+"_turnover" in self.label_model.parameter_dict:
                 self.label_model.parameter_dict[flux+"_turnover"]["v"]=self.label_model.turnover_flux_dict[flux]
          self.selected_turnover=0
          self.turnover_spinbox.delete(0,"end")
          value0=self.label_model.turnover_flux_dict[sorted(self.label_model.turnover_flux_dict,key=lambda v: v.upper())[0]]
          self.turnover_spinbox.insert(INSERT,value0)
          self.turnover_selector_listbox.itemconfig(0, {'fg': 'blue'})
          self.modified_turnovers=[]  
      
      
      def add_turnover_as_parameter(self):
          turnover=self.turnover_selector_listbox.get(self.selected_turnover)
          self.label_model.parameter_dict[turnover+"_turnover"]={"v":self.label_model.turnover_flux_dict[turnover],"lb":0,"ub":500,"max_d":0.1,"type":"turnover","reactions":[turnover]} 
          apply_parameters(self.label_model,parameter_precision=self.parameter_precision,parameter_list=[turnover+"_turnover"])
      
      def add_all_turnover_as_parameter(self):
          for turnover in self.label_model.turnover_flux_dict:
             self.label_model.parameter_dict[turnover+"_turnover"]={"v":self.label_model.turnover_flux_dict[turnover],"lb":0,"ub":1000,"max_d":0.1,"type":"turnover","reactions":[turnover]} 
          apply_parameters(self.label_model,parameter_precision=self.parameter_precision,parameter_list=[turnover+"_turnover"])
          
           
      
      def get_selected_turnover(self):
          self.root.after(50, self.get_selected_turnover)
          previous_slected_turnover =self.selected_turnover
          selection_tupple=self.turnover_selector_listbox.curselection()
          if selection_tupple!=(): 
             self.selected_turnover,=selection_tupple
             if previous_slected_turnover!=self.selected_turnover:
                self.turnover_selector_listbox.itemconfig(self.selected_turnover, {'fg': 'blue'})
                if  self.turnover_selector_listbox.get(previous_slected_turnover) in self.modified_turnovers:
                    self.turnover_selector_listbox.itemconfig(previous_slected_turnover, {'fg': 'green'})
                else:
                    self.turnover_selector_listbox.itemconfig(previous_slected_turnover, {'fg': 'black'})
                turnover_id=self.turnover_selector_listbox.get(self.selected_turnover)
                turnover_value=self.label_model.turnover_flux_dict[turnover_id]  
                self.turnover_spinbox.delete(0,"end")
                self.turnover_spinbox.insert(INSERT,str(turnover_value))
      
      
      def update_turnover(self):
          self.root.after(200, self.update_turnover)
          selected_turnover=self.turnover_selector_listbox.get(self.selected_turnover)
          current_value=self.label_model.turnover_flux_dict[selected_turnover]
          try:
            new_value=float(self.turnover_spinbox.get())
          except:
            new_value=(self.turnover_spinbox.get())
            new_value=float(new_value.replace(",","."))
          if current_value!=new_value:
             if selected_turnover not in self.modified_turnovers:
                self.modified_turnovers.append(selected_turnover)
             self.label_model.turnover_flux_dict[selected_turnover]=new_value
             self.turnover_selector_listbox.itemconfig(self.selected_turnover, {'fg': 'green'})
             if self.waiting_for_label_simulation==False:
                self.root.after(400,self.run_label_simulation)
                self.waiting_for_label_simulation=True
      
      
      def create_fig(self,root,condition="glc"):
          #ventanafigura = Tk()
          figure_height=sum([30+25*len(self.label_model.experimental_dict[condition][emu]) for emu in self.label_model.experimental_dict[condition]])
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
                #print "a"
                pos+=20
                mean=self.label_model.experimental_dict[condition][emu_id][n]["m"]*100 #Percentatge
                #print "b"
                sd=self.label_model.experimental_dict[condition][emu_id][n]["sd"]*100 #Percentatge
                #print "c"
                string="m"+str(n)
                if emu_id in self.label_model.rsm_list:
                   string+="/Sm" 
                fig.create_text(base-30, pos+5, text=string)
                #print "d"
                fig.create_rectangle(base, pos, base+(mean*2), pos+10, fill="gray")
                fig.create_line(base+(mean*2)-(sd*2), pos+5, base+(mean*2)+(sd*2), pos+5)
                if n in self.label_model.condition_simulation_results_dict[condition][emu_id]:
                   #n_emu=label_model.size_variable_dict[size][mid]
                   #print "e"
                   value=round(self.label_model.condition_simulation_results_dict[condition][emu_id][n],2)
                   print value
                   oval_id=fig.create_oval((base+value*100*2)-3, pos+5-2, (base+value*100*2)+3, pos+5+2, fill="red")
                   self.simulated_points_object_dict[condition][emu_id][n]=[oval_id,pos] 
                   """for xx in range(0, len(argument2[x])-3): 
                   fig.create_oval((base+argument2[x][xx+3]*2)-3, pos+5-2, (base+argument2[x][xx+3]*2)+3, pos+5+2, fill="red")"""
             pos+=10
          
          
          return fig
      
      def update_label(self):
        for condition in self.simulated_points_object_dict:
            print "label_updatded"
            base=60
            fig=self.condition_canvas_dict[condition]
            if self.label_window==True:
              for emu_id in self.label_model.experimental_dict[condition]:
                for n in sorted(self.simulated_points_object_dict[condition][emu_id]):
                    new_value=round(self.label_model.condition_simulation_results_dict[condition][emu_id][n],3)
                    oval_id=self.simulated_points_object_dict[condition][emu_id][n][0]
                    pos=self.simulated_points_object_dict[condition][emu_id][n][1]
                    fig.delete(oval_id)
                    oval_id=fig.create_oval((base+new_value*100*2)-3, pos+5-2, (base+new_value*100*2)+3, pos+5+2, fill="red")
                    self.simulated_points_object_dict[condition][emu_id][n]=[oval_id,pos] 
      
      
      
      def update_label_sampling(self):
        for condition in self.simulated_points_object_dict:
            base=60
            fig=self.condition_canvas_dict[condition]
            if self.label_window==True:
              for emu_id in self.label_model.experimental_dict[condition]:
                for n in sorted(self.simulated_points_object_dict[condition][emu_id]):
                    new_value=round(self.label_model.condition_simulation_results_dict[condition][emu_id][n],3)
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
          a,b,c=get_objective_function(self.label_model,force_balance=self.label_model.force_balance,output=False)
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
          self.constrained_model=copy.deepcopy(self.backup_constrained_model)
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
          self.parameters_scrollbar.pack( side = RIGHT, fill=Y,expand=True)
          self.parameters_listbox = Listbox(list_frame, yscrollcommand = self.parameters_scrollbar.set )
          self.parameters_scrollbar.config( command = self.parameters_listbox.yview ) 
          self.parameters_listbox.pack( side = RIGHT, fill=Y,expand=True) 
          for parameter in sorted(self.label_model.parameter_dict,key=lambda v: v.upper()):
              string=parameter+"="+str(self.label_model.parameter_dict[parameter]["v"])
              self.parameters_listbox.insert(END,string)
          button = Button(top, text="Clear Parameters", command=self.clear_all_parameters)
          button.pack(side=TOP)
      
      def clear_all_parameters(self):
          clear_parameters(self.label_model,parameter_dict=None,parameter_list=[], clear_ratios=False,clear_turnover=False,clear_fluxes=True) 
          self.label_model.parameter_dict={}
          self.run_fva() 
          self.view_parameters_window.destroy()  
      
      
      
      def __init__(self,root,label_model): 
          self.ratio_regular_expression = re.compile('(.+)=(.*)') #Regular expression used to obtain ratio from
          self.root = root
          self.root.protocol("WM_DELETE_WINDOW",self.restore_original_parameters)
          self.label_model=label_model
          self.parameter_precision=0.0001
          self.change_threshold=self.parameter_precision*2
          root.title("iso2flux GUI")
          self.constrained_model=label_model.constrained_model
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
          self.variable13b1=0
          self.variable13b2=0
          self.variable31b1_a=0
          self.variable31b1_b=1
          self.variable42b1=0
          self.selected_ratio=-1
          self.waiting_for_fva=False
          self.waiting_for_label_simulation=True
          self.reseting_bounds=False
          self.fraction_of_optimum=1
          self.force_bounds_update=False 
          self.llistafluxos=[]
          """for reaction in label_model.metabolic_model.reactions:
                 self.llistafluxos.append(reaction.id)
          self.llistafluxos2=self.llistafluxos=sorted(self.llistafluxos,key=lambda v: v.upper())
          self.llistafluxos_dict={}
          self.inverse_llistafluxos_dict={}
          for n,x in enumerate(self.llistafluxos):
             self.llistafluxos_dict[n]=x
             self.inverse_llistafluxos_dict[x]=n"""      
          self.add_menubar()
          #self.input_frame()
          self.add_boundaries()
          self.add_ratios_interface()
          self.add_fraction_optimum_interface()
          self.add_turnover_display()
          #add_turnover_frame=Frame(root,width=768, height=576)
          #add_turnover_frame.pack(fill="y", expand="yes") 
          self.add_flux_display()
          self.add_file_output_interface() 
          #poll13b2()
          self.get_bounds_objective() 
          #self.update_bounds()
          #poll42b1()
          #self.get_r1_r2()
          self.update_ratio_list()
          self.run_fva()
          self.update_turnover()  
          self.label_button()
          self.get_ratio_selection()
          self.search_function()
 
      

def launch_gui(label_model):
    """print 
    if __name__ == '__main__':"""
    root = Tk()
    gui = GUI(root,label_model)
    root.mainloop()


