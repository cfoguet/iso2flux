
import numpy as np

from objfunc import objfunc

def get_derivative(label_model,reference_vector,n_row,original_objective,optimal_variables,step=1e-2,min_step=1e-6):
 #print "-----------------------------"   
 #print original_objective
 previous_derivative=1e99
 for n in range(0,100):
    step=step/2
    #Increase vector
    new_vector_up=reference_vector.copy()
    new_vector_up[n_row]+=+step
    label_model.flux_solver_free_fluxes=new_vector_up
    new_objective_up,c=objfunc(label_model,optimal_variables,update_variables=False,simulate_infeasible=True)
    #print step,new_objective_up,label_model.flux_solver_free_fluxes 
    first_derivative_up=(new_objective_up-original_objective)/step
    #print first_derivative_up,previous_derivative
    if abs(first_derivative_up-previous_derivative)<abs(first_derivative_up*0.05):
       #print "break",first_derivative_up-previous_derivative,first_derivative_up*0.05
       break
    previous_derivative=first_derivative_up
    if step<1e-6:
       break
 return (first_derivative_up+previous_derivative)/2




def get_second_derivative(label_model,reference_vector,n_row,n_col,optimal_variables,original_derivative,step=1e-2,min_step=1e-6):
 print "-----------------------------"   
 print original_derivative
 previous_derivative=1e99
 for n in range(0,50):
    step=step/2
    #Increase vector
    new_vector_up=reference_vector.copy()
    new_vector_up[n_col]+=+step
    label_model.flux_solver_free_fluxes=new_vector_up
    new_objective,c=objfunc(label_model,optimal_variables,update_variables=False,simulate_infeasible=True)
    print optimal_variables
    derivative_up=get_derivative(label_model,new_vector_up,n_row,new_objective,optimal_variables,step=1e-3)
    print step,derivative_up,derivative_up-original_derivative,label_model.flux_solver_free_fluxes 
    second_derivative_up=(derivative_up-original_derivative)/step
    print second_derivative_up
    if abs(second_derivative_up-previous_derivative)<abs(second_derivative_up*0.05):
       #print "break",first_derivative_up-previous_derivative,first_derivative_up*0.05
       break
    previous_derivative=second_derivative_up
    if step<min_step:
       break
 return (previous_derivative+second_derivative_up)/2  


def get_std_deviation(label_model,optimal_variables,initial_step=1e-3):
    original_objective,b=objfunc(label_model,optimal_variables)
    original_solution=label_model.flux_solver_free_fluxes.copy()
    n_free_fluxes=len(label_model.flux_solver_free_fluxes)
    hessian=np.empty((n_free_fluxes,n_free_fluxes), dtype=np.float64)
    for n_row in range(0,n_free_fluxes):
        first_derivative=get_derivative(label_model,original_solution,n_row,original_objective,optimal_variables,step=initial_step)
        for n_col in range(0,n_free_fluxes):
            second_derivative=get_second_derivative(label_model,original_solution,n_row,n_col,optimal_variables,original_derivative=first_derivative,step=initial_step)
            """if abs(second_derivative)<1e-10:
               second_derivative=0 """
            hessian[n_row][n_col]=hessian[n_row][n_col]=max(min(round(second_derivative,10),100000000),-100000000)#round(second_derivative,10)#max(min(round(second_derivative,10),100000),-100000)
    label_model.flux_solver_free_fluxes=original_solution
    inverse_hessian=np.linalg.inv(hessian)
    try:
       inverse_hessian=np.linalg.inv(hessian)
    except:
       inverse_hessian=inverse_hessian=np.linalg.pinv(hessian)
    nan_free_fluxes=[]
    #Remove negative variances from the main diagonal
    for n_row in range(0,n_free_fluxes):
        if inverse_hessian[n_row][n_row]<0:
           inverse_hessian[n_row][n_row]=max(inverse_hessian[n_row][n_row],0)
           for reaction_id in label_model.flux_solver_independent_flux_dict[n_row]:
               nan_free_fluxes.append(reaction_id)
    covariance=np.dot(np.dot(label_model.flux_solver_nullmnp,inverse_hessian),np.transpose(label_model.flux_solver_nullmnp))
    flux_std_error_dict={}
    for n_row in range(0,len(covariance)):
        value=pow(covariance[n_row][n_row],0.5)
        if label_model.flux_solver_n_reaction_dict[n_row] in nan_free_fluxes:
           value="nan"
        print n_row,label_model.flux_solver_n_reaction_dict[n_row], value
        flux_std_error_dict[label_model.flux_solver_n_reaction_dict[n_row]]=value
    
    return flux_std_error_dict, hessian,inverse_hessian,covariance

