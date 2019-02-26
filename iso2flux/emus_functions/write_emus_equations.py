import os

def write_emus_equations(label_model,c_code=True,force_balance=True,steady_state=True):
    """
    Writes the emu balances into files .pyx that so that they can be compiled using cython.
    label_model: label_model object
    c_code: bool, deprectaed:
        Currently not used
    force_balance: bool,optional
	If True it will force the isotologues of a given emu to add up to 1 resulting on increased performance
    steady_state: bool, deprecated 
        Currently not used
         
    """
    external_metabolites_list=[]
    eqn_dir=label_model.eqn_dir
    if not os.path.exists(eqn_dir):
       os.makedirs(eqn_dir)
    os.chdir( eqn_dir )
    #label_model.eqn_dir= os.getcwd() 
    f=open("__init__.py","w")
    f.write("from . import *")
    f.close()
    function_str="" 
    for size in sorted(label_model.size_expanded_model_dict.keys()):
        #function_str+="import emu_equations_size%s\ntry:\n  emu_equations_size%s=reload(emu_equations_size%s)\nexcept:\n  pass\n"%(size,size,size)
        #function_str+="def f%s(yy,t=0,label_model=None):\n  condition=label_model.active_condition\n  dy=emu_equations_size%s.emu_equations_size%s(yy,"%(size,size,size)
        function_str+="from emu_equations_size%s import emu_equations_size%s\n"%(size,size)
        function_str+="def f%s(yy,t,condition_size_yy_dict,flux_list,condition_initial_label_yy_dict,condition):\n  dy=emu_equations_size%s(yy,"%(size,size)  
        #function_str+="def f%s(yy,t=0,label_model=None):\n  condition=label_model.active_condition\n  dy=emu_equations_size%s(yy,"%(size,size)  
        string="import cython\nfrom libc.stdlib cimport malloc, free\nimport numpy as np\ncimport numpy as np\n" 
        #string_d="cdef double dy_size%s[%s];\n"%(size,len(label_model.size_variable_dict[size]))
        #for sizes in  xrange(size,0,-1):#string_d+="cdef double yy_size%s[%s];\n"%(sizes,len(label_model.size_variable_dict[sizes]))
        if steady_state==True:
           string+="def emu_equations_size%s("%(size)
        else:
           string+="def emu_odes_size%s("%(size)
        
        for sizes in  range(size,0,-1):
               if sizes!=size:
                  function_str+="condition_size_yy_dict[condition][%s],"%(sizes)
               if steady_state==True or sizes==size:
                  string+="np.ndarray[np.float64_t, ndim=1] np_size%s_yy ,"%(sizes)
               else:
                  string+="np.ndarray[np.float64_t, ndim=2] np_transposed_time_course_size%s ,"%(sizes)
        string+="np.ndarray[np.float64_t, ndim=1] flux_values,"  #"def emu_equations(np_yy,flux_values,yy0)
        function_str+="flux_list,"
        if steady_state==True:
              string+="np.ndarray[np.float64_t, ndim=1] np_input_yy0):\n"
              function_str+="condition_initial_label_yy_dict[condition])\n  return dy\n\n"
        else:
              string+="np.ndarray[np.float64_t, ndim=1] np_input_yy0,np.ndarray[np.float64_t, ndim=1] np_con, np.ndarray[np.float64_t, ndim=1] np_t, double t):\n"
              function_str+="label_model.condition_initial_label_yy_dict[condition],label_model.np_con,np_t,t)\n"
        string+="  cdef np.ndarray[np.float64_t, ndim=1] dy_size%s=np.zeros(%s);\n"%(size,len(label_model.size_variable_dict[size]))
        #Write the flux associations
        if steady_state==False:
           string+="  #volumes\n" 
           string+="  cdef double vol_e=%s;\n"%(label_model.vol_dict["e"])
           string+="  cdef double vol_c=%s;\n"%(label_model.vol_dict["c"])
           string+="  cdef double vol_m=%s;\n"%(label_model.vol_dict["m"])
           string+="  cdef double vol_ec=%s;\n"%(label_model.vol_dict["e"]+label_model.vol_dict["c"])
           string+="  cdef double vol_cm=%s;\n"%(label_model.vol_dict["c"]+label_model.vol_dict["m"])
           string+="  cdef double vol_ecm=%s;\n"%(label_model.vol_dict["c"]+label_model.vol_dict["m"]+label_model.vol_dict["e"])
           string+="  #Changes in extracelular metabolites concentrations\n"
           for met_id in label_model.met_id_isotopomer_dict:
               if label_model.met_id_isotopomer_dict[met_id].input==True:
                  continue
               met=label_model.metabolic_model.metabolites.get_by_id(met_id)
               if met.compartment=="e":
                  isotopomer_object=label_model.met_id_isotopomer_dict[met_id]
                  external_metabolites_list.append(isotopomer_object.id)
                  string+="  cdef double con_"+isotopomer_object.id+"=np_con[%s]+(flux_values[%s]"%(label_model.con_n_dict[isotopomer_object.id],label_model.reaction_n_dict["EX_"+met_id])
                  if "EX_"+met_id+"_reverse" in label_model.reaction_n_dict:
                     n=label_model.reaction_n_dict["EX_"+met_id+"_reverse"]
                     string+="-flux_values[%s]"%(n)
                  string+=")*t/vol_"+isotopomer_object.comp+";\n"
           string+="  #Interpolate value of odes of lower size\n"
           for local_size in label_model.lower_size_dict[size]:
               for mi in label_model.lower_size_dict[size][local_size]:
                   n_variable=label_model.size_variable_dict[local_size][mi]
                   string+=("  cdef double "+mi+"=np.interp(t, np_t, np_transposed_time_course_size%s[%s]);\n")%(local_size,n_variable) 
        string+="  #flux associations\n"
        input_reactions=label_model.size_model_dict[size].reactions.query("_input") 
        for emu_reaction in label_model.size_model_dict[size].reactions:
              if emu_reaction.id not in label_model.emu_reaction_dict:
                 continue
              if force_balance==True and emu_reaction in input_reactions:
                 continue
              """if steady_state==False and "EX_" in emu_reaction.id:
                 continue #If we are working outside of stedy state EX_reactions are not necessary"""
              comment="#"
              if c_code:
                 string+="  cdef double "
              else: 
                 string+="  "
              string+="J_"+emu_reaction.id.replace(")","_").replace("(","_").replace("[","_").replace("]","_").replace(".","_")+"="
              list_of_reactions=[]
              if isinstance(label_model.emu_reaction_dict[emu_reaction.id],list):
                 for reaction in label_model.emu_reaction_dict[emu_reaction.id]:
                     list_of_reactions.append(reaction) 
              else:
                  list_of_reactions=[label_model.emu_reaction_dict[emu_reaction.id]]
              for reaction in list_of_reactions:
                  if reaction not in label_model.merged_reactions_reactions_dict:
                     if "_reverse" in reaction:
                         if reaction[:-8] in label_model.merged_reactions_reactions_dict:
                            reaction=sorted(label_model.merged_reactions_reactions_dict[reaction[:-8]])[0]+"_reverse"
                     reaction_n=label_model.reaction_n_dict[reaction]
                  else:
                     reaction_n=label_model.reaction_n_dict[sorted(label_model.merged_reactions_reactions_dict[reaction])[0]] #We take the first reaction as defined in emu_build_reaction_dict 
                  string+="+flux_values[%s] "%(reaction_n)
                  comment+=" +J_"+reaction.replace(")","_").replace("(","_").replace("[","_").replace("]","_").replace(".","_")
              string+=";"+comment+"\n"
        if force_balance==True:
           mi0_list=[]
           string+="  #force isotopomer balance\n"
           for emu in label_model.emu_dict:
               local_emu_dict=label_model.emu_dict[emu]
               isotopomer_object=label_model.id_isotopomer_object_dict[local_emu_dict["met_id"]]
               if isotopomer_object.input==True:
                  continue
               mi0=local_emu_dict["mid"][0]
               mi0_list.append(mi0)
               local_size=local_emu_dict["size"]
               if mi0 in label_model.size_expanded_model_dict[size].metabolites: #only add the balance for the emus that participate
                  string+="  cdef double  "+mi0.replace(")","_").replace("(","_").replace("[","_").replace("]","_").replace(".","_")+"=1" #consider moving the declaration at the begining
                  for n in local_emu_dict["mid"]:
                      if n==0:
                         continue 
                      else:
                         mi=local_emu_dict["mid"][n]
                         if mi in label_model.size_expanded_model_dict[size].metabolites:
                            if steady_state==False and local_size<size: #If we are doing dynamic simulations lower size emus need to be indicated by variable name
                               string+=" -"+mi
                            else:  
                               string+=" -np_size%s_yy[%s]"%(local_size,label_model.size_variable_dict[local_size][mi])
                  string+=";\n"          
        ###############
        string+="  #compute emu fluxes\n"
        #string_d+="cdef double "
        reverse_reactions=[]
        for reaction in label_model.size_expanded_model_dict[size].reactions:
            if reaction.id in reverse_reactions:
               #print(reaction.id) 
               continue
            product=None
            """if steady_state==False and "EX_" in reaction.id:
                     continue #If we are working outside of stedy state EX_reactions are not necessary"""
            #string_d_flux="J_"+reaction.id+", "
            string_flux ="  cdef double J_"+reaction.id.replace(")","_").replace("(","_").replace("[","_").replace("]","_").replace(".","_")+"= J_"+label_model.expanded_reaction_emu_dict[reaction.id].replace(")","_").replace("(","_").replace("[","_").replace("]","_").replace(".","_")
            all_m0_flag=True #A flag to check if all substrates are m0
            for metabolite in reaction.metabolites:
                print reaction.metabolites
                print reaction.id
                print reaction.reaction
                emu=label_model.expanded_emu_dict[metabolite.id]
                local_size=label_model.emu_dict[emu]["size"] 
                coef=float(reaction.metabolites[metabolite])
                if coef>0:
                   product=metabolite
                   product_emu=emu 
                   product_size=local_size
                elif coef<0:
                     if (metabolite.id!=label_model.emu_dict[emu]["mid"][0] or force_balance==False) and metabolite.id not in label_model.input_n_dict:
                        all_m0_flag=False
                        print [emu,local_size,size]
                        if steady_state==False and local_size<size:
                           variable_name="*"+metabolite.id
                        else: 
                           variable_name="*np_size%s_yy[%s]"%(local_size,label_model.size_variable_dict[local_size][metabolite.id])
                     elif metabolite.id in label_model.input_n_dict:
                          all_m0_flag=False
                          variable_name="*np_input_yy0[%s]"%(label_model.input_n_dict[metabolite.id])
                     else:
                        variable_name="*"+metabolite.id.replace(")","_").replace("(","_").replace("[","_").replace("]","_").replace(".","_")  
                     if coef.is_integer(): #coef.is_integer checks if coefficient is a whole number.Coefs will be integer in all cases except output_reactions
                        string_flux+=variable_name
                        if coef==2:
                           string_flux+=variable_name
                     else:
                        string_flux+=variable_name
                        
            if "reflection" in reaction.notes:
                      if reaction.notes["reflection"] in  label_model.expanded_reaction_emu_dict:
                        reverse_reactions.append(reaction.notes["reflection"])
                        if product==None:
                           variable_name=""
                        elif (product.id!=label_model.emu_dict[product_emu]["mid"][0] or force_balance==False) and product.id not in label_model.input_n_dict:
                           variable_name="*np_size%s_yy[%s]"%(product_size,label_model.size_variable_dict[product_size][product.id])
                        elif product.id in label_model.input_n_dict:
                            variable_name="*np_input_yy0[%s]"%(label_model.input_n_dict[product.id])
                        else:
                           variable_name="*"+product.id
                        print [reaction.id,reaction.notes["reflection"]]
                        string_flux+=(" -J_"+label_model.expanded_reaction_emu_dict[reaction.notes["reflection"]].replace(")","_").replace("(","_").replace("[","_").replace("]","_").replace(".","_")+variable_name)  
            if all_m0_flag==True and force_balance==True: #If we do not simulate m0 we do not need to write a raection that only consumes and produces m0
               continue                 
            string_flux+=";\n"
            string+=string_flux
            #string_d+=string_d_flux
        #string_d+="\n"
        #string_d=string_d[:-3]+";\n"
        #print len(reverse_reactions)
        string+="  \n  \n  #ODES\n"
        for metabolite in label_model.size_expanded_model_dict[size].metabolites:
            emu=label_model.expanded_emu_dict[metabolite.id]
            met_id=label_model.emu_dict[emu]["met_id"]
            local_size=label_model.emu_dict[emu]["size"] 
            if force_balance==True and label_model.emu_dict[emu]["mid"][0]==metabolite.id:
               continue
            if metabolite.id not in label_model.size_variable_dict[size]: #if the metabolite has a lower size or is constant ignore it
               continue
            string+="  dy_size%s[%s]=("%(size,label_model.size_variable_dict[size][metabolite.id])
            for reaction in  metabolite.reactions:
               """if steady_state==False and "EX_" in reaction.id:
                     continue #If we are working outside of stedy state EX_reactions are not necessary"""
               if reaction.id in reverse_reactions: #if it is a reverse reaction ignore it, it will be addded with the forward reaction
                   continue
               coef=reaction.metabolites[metabolite]
               if coef==1:
                  string+=" +J_"+reaction.id.replace(")","_").replace("(","_").replace("[","_").replace("]","_").replace(".","_")
               elif coef==-1:
                  string+=" -J_"+reaction.id.replace(")","_").replace("(","_").replace("[","_").replace("]","_").replace(".","_")
               elif coef>0:
                    string+=" +"+str(coef)+"*J_"+reaction.id.replace(")","_").replace("(","_").replace("[","_").replace("]","_").replace(".","_")
               else: #negative coefficients different than -1
                    string+=" "+str(coef)+"*J_"+reaction.id.replace(")","_").replace("(","_").replace("[","_").replace("]","_").replace(".","_")
            string+=")"
            if steady_state==False:
               string+="/(vol_"+metabolite.compartment+"*"
               if met_id in external_metabolites_list:
                  string+="con_"+met_id+")"
               else:
                  string+="np_con[%s])"%(label_model.con_n_dict[met_id])
            string+=";\n"
        string+="  \n  return(dy_size%s)"%(size)  
        #string+="  \n  return([x for x in dy_size%s])"%(size) 
        
        if c_code:
           if steady_state==True:
              f2=open("emu_equations_size%s_setup.py"%(size),"w")
              setup="try:\n    from setuptools import setup\n    from setuptools import Extension\nexcept ImportError:\n    from distutils.core import setup\n    from distutils.extension import Extension\nfrom Cython.Distutils import build_ext\nimport numpy as np\n\next_modules = [Extension('emu_equations_size%s', ['emu_equations_size%s.pyx'])]\nsetup(name = 'emu_equations_size%s',  cmdclass = {'build_ext': build_ext},  include_dirs = [np.get_include()],  ext_modules = ext_modules)"%(size,size,size)
              f2.write(setup)
              f2.close()
              f=open("emu_equations_size%s.pyx"%(size),"w")
           else:
              f=open("emu_odes_size%s.pyx"%(size),"w")
        else:
           f=open("emu_equations_size%s_py.py"%(size),"w")
        f.write(string)
        f.close()
    function_str+="def get_equations(size_emu_c_eqn_dict):\n"
    for size in label_model.size_expanded_model_dict:
        function_str+="  size_emu_c_eqn_dict[%s]=f%s\n"%(size,size)
    f=open("get_equations.py","w")
    f.write(function_str)
    f.close() 
    os.chdir("..")




