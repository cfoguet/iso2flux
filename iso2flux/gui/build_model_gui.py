import os
import cobra
import json
import tkFileDialog
import Tkinter
import iso2flux
from iso2flux.classes.label_model import Label_model #Import the label_model class
from iso2flux.label_propagation_functions.find_missing_reactions import find_missing_reactions
from iso2flux.misc.save_load_iso2flux_model import save_iso2flux_model,load_iso2flux_model
from iso2flux.input_functions.read_experimental_mid import read_experimental_mid
from iso2flux.input_functions.read_isotopomer_model import read_isotopomer_model
from iso2flux.input_functions.read_constraints import read_flux_constraints 
from Tkinter import *
from iso2flux.input_functions.read_isotoflux_settings import read_isotoflux_settings


class build_model_gui:
    def get_file(self,root,text="",default_file=None,open_file_widget=None,row=0):
        file_frame=root  
        #file_frame.grid(side=TOP)
        Label(file_frame,text=text).grid(row=row,column=0)
        file_entry=Entry(file_frame)
        file_entry.grid(row=row,column=1) 
        if default_file!=None:
           file_entry.insert(0,default_file)     
        button = Button(file_frame,text="Select",command=open_file_widget)
        button.grid(row=row,column=2)
        return file_entry
    
    def get_working_directory(self):
        wd=tkFileDialog.askdirectory(title="Select working directory:",parent=self.root)  
        os.chdir(wd)
        self.wd_entry.delete(0, END)
        self.wd_entry.insert(0,wd) 
       
    def get_sbml_model(self):
        loaded_file = tkFileDialog.askopenfile(title='Choose a SBML model',filetypes=[("sbml",".sbml"),("xml",".xml")],parent=self.root)
        self.sbml_entry.delete(0, END) 
        self.sbml_entry.insert(0,str(loaded_file.name)) 
    
    def get_label_propagation_rules(self):
        loaded_file = tkFileDialog.askopenfile(title='Choose a set of label propagation rules',filetypes=[("xlsx","*.xlsx"),("CSV","*.csv")],parent=self.root)
        self.label_rules_entry.delete(0, END) 
        self.label_rules_entry.insert(0,str(loaded_file.name))
        
    def get_exp_data(self):
        loaded_files = tkFileDialog.askopenfiles(title='Choose epxerimental isotopologues data',filetypes=[("xlsx","*.xlsx"),("CSV","*.csv")],parent=self.root)
        file_names=[x.name for x in loaded_files]
        self.e_data_entry.delete(0, END) 
        for file_name in file_names:
           self.e_data_entry.insert(0,file_name+",")
        self.e_data_entry.delete(len(self.e_data_entry.get())-1) #delete last coma
        
    def get_flux_constraints(self):
        loaded_file = tkFileDialog.askopenfile(title='Chose constraints file',filetypes=[("xlsx","*.xlsx"),("CSV","*.csv")],parent=self.root)
        self.constraints_entry.delete(0, END) 
        self.constraints_entry.insert(0,str(loaded_file.name))
        
    def get_settings(self):
        loaded_file = tkFileDialog.askopenfile(title='Chose Settings file',filetypes=[("xlsx","*.xlsx"),("CSV","*.csv")],parent=self.root)
        self.settings_entry.delete(0, END) 
        self.settings_entry.insert(0,str(loaded_file.name))
    def new_model(self):
        if ""==self.sbml_entry.get() or ""==self.e_data_entry.get or ""==self.label_rules_entry.get():
             print "Mandatory Input missing"
             return
        p_dict={'reactions_with_forced_turnover': [], 'annealing_cycle_time_limit': 1800, 'confidence_max_absolute_perturbation': 10, 'turnover_exclude_EX': True, 'annealing_n_processes': 4, 'annealing_p0': 0.4, 'identify_free_parameters_add_turnover': True, 'minimum_sd': 0.01, 'annealing_max_perturbation': 1, 'turnover_upper_bound': 100, 'confidence_perturbation': 0.1, 'annealing_m': 1000, 'annealing_n': 10, 'annealing_relative_max_sample': 0.35, 'confidence_min_absolute_perturbation': 0.05, 'annealing_pf': 0.0001, 'confidence_significance': 0.95, 'identify_free_parameters_change_threshold': 0.005, 'parameter_precision': 0.0001, 'fraction_of_optimum': 0, 'lp_tolerance_feasibility': 1e-09, 'identify_free_parameters_n_samples': 200, 'annealing_relative_min_sample': 0.2, 'annealing_iterations': 2}
        model=cobra.io.read_sbml_model(self.sbml_entry.get())
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
        self.label_model=label_model=Label_model(model,lp_tolerance_feasibility=p_dict["lp_tolerance_feasibility"],parameter_precision=p_dict["parameter_precision"],reactions_with_forced_turnover=p_dict["reactions_with_forced_turnover"])
        read_isotopomer_model(label_model,self.label_rules_entry.get())
        find_missing_reactions(label_model)
        #loaded_file = tkFileDialog.askopenfile(title='Choose experimental measuments file',filetypes=[("xlsx",".xlsx")]) 
        #fileName=loaded_file.name
        e_data_names=self.e_data_entry.get().replace("[","").replace("]","").split(",")
        emu_dict0,label_model.experimental_dict =read_experimental_mid(label_model,e_data_names,emu0_dict={},experimental_dict={},minimum_sd=p_dict["minimum_sd"])
        label_model.build(emu_dict0,force_balance=True,recompile_c_code=True,remove_impossible_emus=True,isotopic_steady_state=True,excluded_outputs_inputs=[],turnover_upper_bound=p_dict["turnover_upper_bound"],clear_intermediate_data=False,turnover_exclude_EX=p_dict['turnover_exclude_EX'])
        
        save_iso2flux_model(label_model,name="project",write_sbml=True,gui=True)
        self.root.destroy()
        
    def __init__(self,root):
        self.root=root
        root.title="Iso2flux"
        main_frame=Frame(root)
        #self.wd_entry=self.get_file(root,text="Working directory",default_file=os.getcwd(),open_file_widget=self.get_working_directory,row=)
        self.sbml_entry=self.get_file(root,text="SBML model",default_file=None,open_file_widget=self.get_sbml_model,row=0)
        self.label_rules_entry=self.get_file(root,text="Label propagation rules",default_file=None,open_file_widget=self.get_label_propagation_rules,row=1)
        self.e_data_entry=self.get_file(root,text="Experimental measurements",default_file=None,open_file_widget=self.get_exp_data,row=2)
        self.constraints_entry=self.get_file(root,text="Flux constraints (Optional)",default_file=None,open_file_widget=self.get_flux_constraints,row=3)
        self.settings_entry=self.get_file(root,text="Advanced settings(Optional)",default_file=None,open_file_widget=self.get_settings,row=4)
        button=Button(root,text="Create Iso2flux model",command=self.new_model)
        button.grid(row=6,column=1)



import Tkinter
import os
import tkFileDialog
import iso2flux
from iso2flux.misc.save_load_iso2flux_model import save_iso2flux_model,load_iso2flux_model
from iso2flux.gui.gui import build_model_gui, launch_gui

def launch_build_model_gui():
    tk.destroy()
    global label_model
    root = Tkinter.Tk()
    gui = build_model_gui(root)
    root.mainloop()
    label_model=gui.label_model


def load_model():
    tk.withdraw()
    global label_model
    label_model=load_iso2flux_model(gui=True)
    tk.destroy()
    

if __name__ == "__main__": 
   tk=Tkinter.Tk()
   tk.withdraw()
   wd=tkFileDialog.askdirectory(title="Select working directory")
   os.chdir(wd)
   tk.deiconify()
   tk.title("Iso2flux")
   text=Tkinter.Label(tk,text="Welcome to Iso2Flux, select one option:   ")
   text.pack(side="top")
   load_button=Tkinter.Button(tk,text="Load model",command=load_model)
   create_button=Tkinter.Button(tk,text="Create new model",command=launch_build_model_gui)
   load_button.pack()
   create_button.pack()
   tk.mainloop()
   launch_gui(label_model)

