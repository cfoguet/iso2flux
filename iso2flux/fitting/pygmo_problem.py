from PyGMO.problem import base
class iso2flux_problem(base):
    """
    iso2flux_problem
    """
    
    def __init__(self, dim=len(label_model.variable_vector),i_dim=0,n_obj=1,c_dim=0,c_ineq_dim=0,c_tol=1e-6):
        # First we call the constructor of the base class telling PyGMO
        # what kind of problem to expect ('dim' dimensions, 1 objective, 0 contraints etc.)
        super(iso2flux_problem,self).__init__(dim)#, i_dim, n_obj, c_dim, c_ineq_dim, c_tol)
        self.label_weight=1
        self.flux_weight=1
        self.target_flux_dict=None
        self.max_chi=1e16
        self.max_flux=1e16
        self.flux_penalty_dict={}
        self.verbose=None
        self.fmin=1e36
        self.label_unfeasible_penalty=1e9
        self.flux_unfeasible_penalty=1e6
    # Reimplement the virtual method that defines the objective function.
    def  set_parameters(self,parameter_dict):
         if "label_weight" in parameter_dict:
             self.label_weight=float(parameter_dict["label_weight" ])
         if "flux_weight" in parameter_dict:
             self.flux_weight=float(parameter_dict["flux_weight" ])
         if "target_flux_dict" in parameter_dict:
             self.target_flux_dict=parameter_dict["target_flux_dict" ]
         if "max_chi" in parameter_dict: 
             self.max_chi=parameter_dict["max_chi"]
         if "max_flux" in parameter_dict: 
             self.max_flux=parameter_dict["max_flux"]
         if "verbose" in parameter_dict:
              self.verbose=parameter_dict["verbose"]
         if "flux_penalty_dict" in parameter_dict:
              self.flux_penalty_dict=parameter_dict["flux_penalty_dict"]
         if "label_unfeasible_penalty" in parameter_dict:
            self.label_unfeasible_penalty=parameter_dict["label_unfeasible_penalty"]
         if "flux_unfeasible_penalty" in parameter_dict:
            self.flux_unfeasible_penalty=parameter_dict["flux_unfeasible_penalty"]    
    def _objfun_impl(self, x):
       #print "running ob"
       # Compute the sphere function
       try:
        f,obj_dict = objfunc(label_model,x,label_weight=self.label_weight,mode="fsolve",target_flux_dict=self.target_flux_dict,max_chi=self.max_chi,max_flux=self.max_flux,flux_penalty_dict=self.flux_penalty_dict,verbose=False,flux_weight=self.flux_weight,label_unfeasible_penalty=self.label_unfeasible_penalty,flux_unfeasible_penalty=self.flux_unfeasible_penalty)
        if f<self.fmin and self.verbose:
             self.fmin=f
             #print obj_dict
             flux_obj=obj_dict["flux_score"]
             label_obj=obj_dict["chi2_score"]
             fltarget_obj=obj_dict["fltarget_obj"]
             output="obj="+str(round(f,3))
             output+=" Chi2="+str(round(label_obj,3))
             if flux_obj>0:
                output+=" flux="+str(round(flux_obj,3))
             if fltarget_obj>0:
                output+=" "+self.target_flux_dict["reaction"]+"="+str(round(label_model.reversible_flux_dict[self.target_flux_dict["reaction"]],3))+" "+self.target_flux_dict["dir"]
             
             print output
        else:
           pass
           #print "not best: "+str(obj_dict) 
        # Note that we return a tuple with one element only. In PyGMO the objective functions
        # return tuples so that multi-objective optimization is also possible.
        #print f
        return (f, )
       except:
        print "Unkown error in objfun"
        return (99999999999999999, )
    """def _compute_constraints_impl(self,x):
        solution=[]
        for function in constraint_functions:
            solution.append(-function(x))
        print sum(solution)
        return tuple(solution)"""
