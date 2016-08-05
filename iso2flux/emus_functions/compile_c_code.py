import sys
import subprocess
import os
#from ..dynamic.interpolate_timecourse import interpolate_timecourse

#label_model=None

def compile_c_code(label_model,recompile=True,steady_state=True):
    """
    Compiles the emu balances into C using cython. 

    label_model: label_model object
    recompile: bool, optional
            if set to false the emu balances won't be compiled and existing ones will be used
    steady_state: bool,optional
            currently unnused
    """
    if recompile==False:
       return
    os.chdir(label_model.eqn_dir)
    for size in label_model.emu_size_dict:
      if steady_state==True:
           theproc = subprocess.Popen([sys.executable, 'emu_equations_size%s_setup.py'%(size), 'build_ext','--inplace'])
           theproc.communicate()
    os.chdir("..")            

           
    
