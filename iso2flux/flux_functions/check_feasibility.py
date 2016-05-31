from multiprocessing import Pool
import time
import copy
from cobra.io import write_sbml_model

def get_status(input_dict={"model":None,"tolerance_feasibility":1e-6}):
    try:
       return input_dict["model"].optimize(tolerance_feasibility=input_dict["tolerance_feasibility"]).status
    except:
       return infeasible

def check_feasibility(metabolic_model,time_limit=60,tolerance_feasibility=1e-6,pool=None,debug=True):
    #Some solvers may crash python with unfeasible parameters, by creating a new process we prevent that error from crashing the main process
    #model_copy=copy.deepcopy(metabolic_model)
    if pool==None:
       pool = Pool(processes=1)
    task_start = time.time()   # start time
    task = pool.map_async(get_status,[{"model":metabolic_model,"tolerance_feasibility":tolerance_feasibility}])
    while not(task.ready()): 
                 #print task._number_left
                 if (time.time() - task_start) > time_limit: # check maximum time (user def.)               print "timeout"
                    pool.terminate()                    # kill old pool
                    timeout = True                         # redo computation
                    pool = Pool(processes=1)
                    if debug:
                       print "infeasible"
                       write_sbml_model(metabolic_model, "feasibility.sbml")   
                    return "infeasible",pool
    status=task.get()[0]
    #pool.close()
    if debug:
       print status
    return status,pool
