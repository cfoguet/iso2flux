"""
This script uses gene expression data to set add additional minimization weights to reactions with gene expression bellow a given threshold. 
--iso2flux_model_file=,-i (mandatory):  path to the iso2flux model that will be solved. Iso2flux model can be created with the create_iso2flux_model.py and has  “.iso2flux” as extension. 
--gene_expression_file=,-g(optional): path to the file with gene expression data. See here
--output_flux_penalty_file=,-o(optional): Used to define a prefix that will be added to the output flux penalty file. It can be used both to name the output and to select the directory where the output will be saved (directory must already exist)
--reference_flux_penalty_file=,-f(optional): path to the reference flux_penalty_file, gene expression penalties will use the values defined in this file as a starting point. Default is the default reference_flux_penalty_file associated to the iso2flux_model loaded.   
--gene_prefix==,-p (optional): prefix of that the constraint based model uses for genes (e.g gene_) that is not present in the gene_expression_file. Default is empty (no prefix )  
--gene_sufix==,-s (optional): suffix of that the constraint based model uses for genes (_at) that is not present in the gene_expression_file. Default is empty (no suffix )  
--remove_penalty_for_spontaneous,-r (Optional): If this flag is set the penalty for reactions not associated to a gene will be 0 

"""
import tkFileDialog
import Tkinter
import sys, getopt
from iso2flux.misc.save_load_iso2flux_model import load_iso2flux_model
from iso2flux.flux_functions.minimal_flux import create_minimal_fux_model,flux_minimization_fva
from iso2flux.flux_functions.build_flux_penalty_dict import build_flux_penalty_dict, add_gene_expression_to_flux_penalty,write_flux_penalty_dict,read_flux_penalty_dict_from_file

try:
 argv=sys.argv[1:]
 opts, args = getopt.getopt(argv,"i:g:l:f:p:s:o:",["iso2flux_model_file=","gene_expression_file=","low_expression_threshold=","reference_flux_penalty_file=","gene_prefix=","gene_sufix=","output_flux_penalty_file="])
 #opts, args = getopt.getopt(argv,"e:l:s:p:c:o:t:f:w:g:m:q",["experimental_data_file=","label_propagation_files=","sbml_model=","parameters_file=","constraints_file=","output_prefix=","time=","factor=","working_directory=","gene_expression_file=","metabolomics_file=","quick_analysis"])
 #opts, args = getopt.getopt(sys.argv[1:],"hi:o:",["ifile=","ofile="])
except getopt.GetoptError as err:
   # print help information and exit:
   print str(err)  # will print something like "option -a not recognized":
   #print "wrong parameters"#'test.py -i <inputfile> -o <outputfile>'
   sys.exit(2)

iso2flux_file_name=None
gene_expression_file=None
output_flux_penalty_file=None
reference_flux_penalty_file=None
gene_prefix=""
gene_sufix=""

for opt, arg in opts:
         print [opt,arg]
         if opt in ("--iso2flux_model_file=","-i"):
             iso2flux_file_name=arg
         elif opt in ("--gene_expression_file=","-g"):
             gene_expression_file=arg       
         elif opt in ("--output_flux_penalty_file=","-o"):
             output_flux_penalty_file=arg
         elif opt in ("--low_expression_threshold=","-l"):
             low_expression_threshold=float(arg)
         elif opt in ("--reference_flux_penalty_file=","-f"):
            reference_flux_penalty_file=arg
         elif opt in ("--gene_prefix==","-p"):
            gene_prefix=arg
         elif opt in ("--gene_sufix==","-s"):
            gene_sufix=arg




if iso2flux_file_name==None:
    tk=Tkinter.Tk()
    tk.withdraw()
    loaded_file = tkFileDialog.askopenfile(title='Select iso2flux file',filetypes=[("iso2flux",".iso2flux")]) 
    iso2flux_file_name=loaded_file.name
    tk.destroy()

if ".iso2flux" not in iso2flux_file_name:
    iso2flux_file_name+=".iso2flux"


label_model=load_iso2flux_model(iso2flux_file_name)

if gene_expression_file==None:
    tk=Tkinter.Tk()
    tk.withdraw()
    loaded_file = tkFileDialog.askopenfile(title='Select gene expression file',filetypes=[("csv","*.csv"),("xlsx","*.xlsx"),('All files','*.*')]) 
    gene_expression_file=loaded_file.name
    tk.destroy()

default_name=iso2flux_file_name.replace(".iso2flux","_flux_penalty.csv")
if output_flux_penalty_file==None:
   output_flux_penalty_file=default_name

if reference_flux_penalty_file==None:
   #Try to see if the flux_penalty dict saved with the default name exists
   try:
     #default_name=iso2flux_file_name.replace(".iso2flux","_flux_penalty.csv")   
     flux_penalty_dict=read_flux_penalty_dict_from_file(default_name)
   except:
     print default_name+"_not found. Creating new flux penalties"
     flux_penalty_dict=build_flux_penalty_dict(label_model,base_penalty=1,fn=None)

else:
   print "loading "+reference_flux_penalty_file
   flux_penalty_dict=read_flux_penalty_dict_from_file(reference_flux_penalty_file)

gene_expression_options_dict={"file_name":gene_expression_file,"absent_gene_expression":50,"percentile":True,"low_expression_threshold":low_expression_threshold,"gene_method":"average","gene_prefix":gene_prefix,"gene_sufix":gene_sufix}
flux_penalty_dict=add_gene_expression_to_flux_penalty(label_model.constrained_model,flux_penalty_dict,gene_expression_options_dict)
print "new flux penalty file saved as:"
write_flux_penalty_dict(output_flux_penalty_file,flux_penalty_dict,irreversible_metabolic_model=None)


