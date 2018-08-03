import math
class isotopomer:
    """isotopomer is a class for holding information regarding
    a isotopomers
    """
    #Class used to define metabolites where label propgation is to be simulated
    def __init__(self,label_model,reference_metabolite_id,ncarbons=None,symmetric=False,label_input=False,iso_id=None):
        """
        label_model: label model object
        reference_metabolite_id: string or a list of strings. 
		COBRA model metabolite/s id where label propgation will be simulated. If a metabolite list is entered all the metabolites will be assumed to share the same label distibution, to be part of the same pool 
        ncarbons: int, optional. 
		Number of carbons for the metabolite/s. If this field is absent the number of carbons will be taken from the SBML model
        symmetric: bool,optional. 
		True indicates that the metabolite is symmetrical
        label_input: bool,optional. 
		True indicates that the label in this metabolite/s is constant
        iso_id: string,optional. 
		Name given to the isotopomer object. 
        """  
        model=label_model.irreversible_metabolic_model
        pool=False 
        if ncarbons==None:
           ncarbons=-1
        if isinstance(reference_metabolite_id, list): #Check if there is a pool of metabolites
           if len(reference_metabolite_id)!=1:
               pool=True
           else:
               reference_metabolite_id=reference_metabolite_id[0]
                     
        if not pool:
               reference_metabolite=model.metabolites.get_by_id(reference_metabolite_id)
               if ncarbons<=0: #If no number of carbons is provided in the initialization read from the cobra metabolite formula 
                  if "C" in reference_metabolite.elements:
                      ncarbons=reference_metabolite.elements["C"] 
                  else: #If the formula does not contain Carbons raise and error
                      raise ValueError("Error: cannot determine the number of Carbons for "+reference_metabolite_id) 
        else: 
               reference_metabolite=[]
               reference_metabolite_c=[] #The number of carbibs will be stored here from all pool membres to check for inconsistencies 
               for x in reference_metabolite_id:
                   individual_met=model.metabolites.get_by_id(x)
                   reference_metabolite.append(individual_met)
                   if ncarbons<=0:
                      if "C" in individual_met.elements:
                          ncarbons=individual_met.elements["C"]        
                          reference_metabolite_c.append(ncarbons)
                   else: 
                      reference_metabolite_c=[ncarbons] 
               if max(reference_metabolite_c)!=  min(reference_metabolite_c):
                  raise ValueError("Error: Inconsitent the number of Carbons in pool")
                     
                          
        if self not in label_model.isotopomer_object_list: label_model.isotopomer_object_list.append(self)
        self.n=ncarbons
        self.len=int(math.pow(2.0,ncarbons)) #Used for full isotopomer simulations. Obsolete with emus
        
        self.input=label_input #Defines a metabolite as a label input, ie labelled glucose
        
        #Define the sarting position of the isotopomers of this metabolite in the variable array (similar to isodyn) 
        """if global_n_isotopomer_dict=={}:
           self.starting=0
        else:
           self.starting=len(global_n_isotopomer_dict)"""
        self.ref_met=reference_metabolite
        #Define compartment of isotopomers (used in dynamic simulations)      
        if not pool:
           self.comp=reference_metabolite.compartment
           self.name=reference_metabolite.name
           self.id=reference_metabolite.id
           self.pool=False
           label_model.metabolite_isotopomers_dict[reference_metabolite]=self
           label_model.met_id_isotopomer_dict[reference_metabolite.id]=self 
           label_model.metabolite_id_isotopomer_id_dict[reference_metabolite.id]=self.id
           label_model.isotopomer_id_metabolite_id_dict[self.id]=reference_metabolite.id
        else: 
           compartment_list=[]
           self.pool=True 
           self.name="pool of"
           self.id="pool"
           for x in reference_metabolite: 
               if x.compartment not in compartment_list:
                  compartment_list.append(x.compartment) 
               self.name+=" "+str(x.name)
               if len(self.id)<50:
                  self.id+="_"+x.id
               elif "_etc" not in self.id:
                  self.id+="_etc"
               label_model.metabolite_isotopomers_dict[x]=self
               label_model.met_id_isotopomer_dict[x.id]=self
               
           if "e" in compartment_list:
               if "c" in compartment_list:
                   if "m" in compartment_list:
                       self.comp="ecm"
                   else:
                       self.comp="ec"
               else: 
                   self.comp="e"
           elif "c" in compartment_list:
               if "m" in compartment_list:
                   self.comp="cm"
               else:
                   self.comp="c" 
           else:
               self.comp="m"    
        if iso_id!=None: #If a id is given, overwrite default id
           self.id=iso_id
        #Build the metabolite_id, isotopomer_id dicts:
        if pool==True:
           for x in reference_metabolite:
               label_model.metabolite_id_isotopomer_id_dict[x.id]=self.id
           label_model.isotopomer_id_metabolite_id_dict[self.id]=[x.id for x in reference_metabolite]
        else:
           label_model.metabolite_id_isotopomer_id_dict[reference_metabolite.id]=self.id
           label_model.isotopomer_id_metabolite_id_dict[self.id]=[reference_metabolite.id]
        
        label_model.id_isotopomer_object_dict[self.id]=self
        self.dict={}
        self.n_dict={}
        self.met_n_dict={}
        self.symm=symmetric
