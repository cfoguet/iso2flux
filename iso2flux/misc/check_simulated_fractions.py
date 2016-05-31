from ..emus_functions.solver import solver #To do, move it to emu equations
def check_simulated_fractions(label_model): #Checks that the simulation can be succesfully run and that there are not negative or larger than 1 fractions
  a,b=solver(label_model,mode="fsolve")
  """if a!="success":
     print "solver failed"
     return"""
  for condition in label_model.initial_label:
    error=False
    label_model.active_condition=condition
    for size in label_model.condition_size_yy0_dict[condition]:
        for n,x in enumerate(label_model.condition_size_yy_dict[condition][size]):
          if abs(x)>1.01:
           print(str(n)+" "+label_model.size_inverse_variable_dict[size][n]+" "+str(x))
           print("Warning: Value larger than 1 found on solution")
           error=True 
          if abs(round(x,6))<0:
           print(str(n)+" "+label_model.size_inverse_variable_dict[size][n]+" "+str(x))
           print("Warning: Negative value found on solution")
           error=True
    if error==True:
       print("Test failed")
    else:
       print("Test succesful")  
