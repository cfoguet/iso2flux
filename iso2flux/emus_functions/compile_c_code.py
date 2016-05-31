import sys
import subprocess
import os
#from ..dynamic.interpolate_timecourse import interpolate_timecourse

#label_model=None

def compile_c_code(label_model,recompile=True,steady_state=True):
    os.chdir("equations")
    for size in label_model.emu_size_dict:
      if steady_state==True:
        if recompile==True:
           theproc = subprocess.Popen([sys.executable, 'emu_equations_size%s_setup.py'%(size), 'build_ext','--inplace'])
           theproc.communicate()
    os.chdir("..")            

           
    
