from ..flux_functions.get_fluxes import get_fluxes
from ..emus_functions.set_equations_variables import set_equations_variables
def check_steady_state(label_model,only_initial_m0=True,threshold=1e-6,precision=10,fn=None,fn_mode="a",verbose=False):
    passed_flg=True
    if fn!=None:
       output_file=open(fn,fn_mode)
       output_file.write("Steady state test:\n")
    label_model.constrained_model.optimize()
    get_fluxes(label_model,model=None,precision=precision)
    if only_initial_m0:
       set_equations_variables(label_model,force_balance=label_model.force_balance,set_initial_label=False)
    for condition in label_model.initial_label:
        #print "condition:"+condition
        label_model.active_condition=condition
        for size in label_model.size_emu_c_eqn_dict:
            label_model.condition_size_yy_dict[condition][size]=label_model.condition_size_yy0_dict[condition][size]
            soln=label_model.size_emu_c_eqn_dict[size](label_model.condition_size_yy_dict[condition][size],0,label_model.condition_size_yy_dict,label_model.flux_list,label_model.condition_initial_label_yy_dict, condition)
            for n,x in enumerate(soln):
                if abs(x)>threshold:
                   passed_flg=False  
                   print(str(n)+" "+label_model.size_inverse_variable_dict[size][n]+" "+str(x))
                   if fn!=None:
                       output_file.write(str(n)+" "+label_model.size_inverse_variable_dict[size][n]+" "+str(x)+"\n")
    if only_initial_m0:
       set_equations_variables(label_model,force_balance=label_model.force_balance,set_initial_label=True)
       if  passed_flg==False:
           if fn!=None:
              output_file.write("failed\n\n")
              output_file.close() 
           print "test failed"
       else:
           if verbose:
              print "test passed"
           if fn!=None:
              output_file.write("passed\n\n")
              output_file.close()  
    return passed_flg
