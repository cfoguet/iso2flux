from ..flux_functions.get_fluxes import get_fluxes
from ..emus_functions.set_equations_variables import set_equations_variables
def check_steady_state(label_model,only_initial_m0=True,threshold=1e-6):
    passed_flg=True
    label_model.constrained_model.optimize()
    get_fluxes(label_model)
    if only_initial_m0:
       set_equations_variables(label_model,force_balance=label_model.force_balance,set_initial_label=False)
    for condition in label_model.initial_label:
        print "condition:"+condition
        label_model.active_condition=condition
        for size in label_model.size_emu_c_eqn_dict:
            label_model.condition_size_yy_dict[condition][size]=label_model.condition_size_yy0_dict[condition][size]
            soln=label_model.size_emu_c_eqn_dict[size](label_model.condition_size_yy0_dict[condition][size],0,label_model)
            for n,x in enumerate(soln):
                if abs(x)>threshold:
                   passed_flg=False  
                   print(str(n)+" "+label_model.size_inverse_variable_dict[size][n]+" "+str(x))
    if only_initial_m0:
       set_equations_variables(label_model,force_balance=label_model.force_balance,set_initial_label=True)
       if  passed_flg==False:
           print "test failed"
       else:
           print "test passed"
