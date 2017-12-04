import re
import copy
import numpy as np
from ..misc.read_spreadsheets import read_spreadsheets
from ..misc.write_spreadsheet import write_spreadsheet
from ..misc.round_functions import round_up, round_down
from ..fitting.objfunc import objfunc

number_finder = re.compile("[\d]+\.?[\d]*")
class gene_expression_score:
    """Based on the work by : #Schmidt BJ1, Ebrahim A, Metz TO, Adkins JN, Palsson B, Hyduke DR. GIM3E: condition-specific models of cellular metabolism developed from metabolomics and expression data Bioinformatics. 2013 Nov 15;29(22):2900-8. doi: 10.1093/bioinformatics/btt493. Epub 2013 Aug 23."""
    def __init__(self, value):
        self.str_value = str(value)
        self.value = float(value)
    
    def __add__(self, other):
        # addition is like AND
        return gene_expression_score(min(self.value, other.value))
    
    def __mul__(self, other):
        # multiplication is like OR
        return gene_expression_score(max(self.value, other.value))
    
    def __neg__(self): #CFC Added
        return gene_expression_score (-self.value) #CFC Added

def evaluate_gene_expression_string(gene_expression_string):
    """Based on the work by : #Schmidt BJ1, Ebrahim A, Metz TO, Adkins JN, Palsson B, Hyduke DR. GIM3E: condition-specific models of cellular metabolism developed from metabolomics and expression data Bioinformatics. 2013 Nov 15;29(22):2900-8. doi: 10.1093/bioinformatics/btt493. Epub 2013 Aug 23."""
    """ penalty string will have:
        * 'or' statements which need to be converted to min
        * 'and' statements which need to be converted to max
    >>> evaluate_penalty("(1 and 2 and 3)")
    max(1, 2, 3)"""
    # if there are no ands or ors, we are done
    gene_expression_string = gene_expression_string.lower()  # don't want to deal with cases
    
    if "and" not in gene_expression_string and "or" not in gene_expression_string:
        return eval(gene_expression_string)
    # we will replace AND and OR with addition and multiplication
    # equivalent to min/max
    #gene_expression_string = gene_expression_string.replace("or", "+").replace("and", "*") #Changed from GIMME
    gene_expression_string = gene_expression_string.replace("or", "*").replace("and", "+")
    # replace the numbers with the custom class which have overloaded AND/OR
    values = [gene_expression_score(i) for i in number_finder.findall(gene_expression_string)]
    values_strings = tuple("values[%i]" % i for i in range(len(values)))
    gene_expression_string = number_finder.sub("%s", gene_expression_string)
    gene_expression_string = gene_expression_string % values_strings
    return eval(gene_expression_string).value


def get_expression(model,file_name="gene_expression_data.xlsx",gene_method="average",gene_prefix="",gene_sufix=""):
    """
    Reads the gene expression file
    model: cobra model object
    file_name: str
          Name of the file with the gene expression data. It should be either a CSV or a XLSX file with gene names in the first column and gene expression value in the second column
    gene_method: str
          Determines wich gene expression value should be given to a gene when there are multiples entries for the same gene. It can either be "average" (it will use the average of all values), "maximum" (it will use the maximum) or "minimum" (it will use the minimum)
    gene_prefix: str
          Prefix used by genes in the cobra model but not present in the gene expression file
    gene_sufix= str
          Sufix  used by genes in the cobra model but not present in the gene expression file. In Recon 1 alternative transtricpts are indicated by appending _AT1, _AT2 , AT3_ at the end of gene. If in the gene expression file alternative scripts are not indicated in that case _AT should be defined as Sufix
    """
    genexpraw_dict={}
    spreadsheet_dict=read_spreadsheets(file_names=file_name,csv_delimiter=',',more_than_1=False,tkinter_title="Chose a file")
    gene_expression_dict=gene_expression_dict={}
    #wb = load_workbook(file_name, read_only=True)
    #ws=wb.active
    for sheet in spreadsheet_dict:
      for row in spreadsheet_dict[sheet]:
        geneid=str(row[0])
        gene_expression_value=float(row[1])
        if geneid in (None,"","NA","---"):
           continue
        if gene_expression_value==None or geneid=="":
           continue  
        gene_matches=model.genes.query(geneid)
        regular_expression=""
        if gene_prefix=="":
           regular_expression+="^"
        else:
           regular_expression+=gene_prefix 
        regular_expression+=geneid
        if gene_sufix=="":
           regular_expression+="$"
        elif gene_sufix==".":
           regular_expression+="\." 
        else:
           regular_expression+=gene_sufix  
        print regular_expression
        for gene in gene_matches:
            if re.search(regular_expression,gene.id)!=None: 
               if gene.id not in genexpraw_dict:
                  genexpraw_dict[gene.id]=[]
               genexpraw_dict[gene.id].append(gene_expression_value)
    for geneid in genexpraw_dict:
        list_of_values=genexpraw_dict[geneid]
        if gene_method=="average" or "mean":
           value=np.mean(list_of_values)
        elif gene_method in("max","maximum"):
           value=max(list_of_values)
        elif gene_method in ("min","minimum"):
           value=min(list_of_values)
        else:
            print "Error: Method not supoprted"
            return {}
        gene_expression_dict[geneid]=round(value,4)     
    print gene_expression_dict      
    return genexpraw_dict, gene_expression_dict

def get_gene_exp(model,absent_gene_expression=50,percentile=True,file_name="gene_expression_data.xlsx",gene_method="average",gene_prefix="",gene_sufix=""):
    """
    Assigns each reaction a expression value based on the gene_expression file and the GPR rules
    """
    genexpraw_dict,expression_dict=get_expression(model,file_name=file_name,gene_method=gene_method,gene_prefix=gene_prefix,gene_sufix=gene_sufix)
    if percentile==True: #If percentile is set to True, low_expression_threshold and high_expression_threshold
       value_list=[] 
       for x in genexpraw_dict:
           value_list+=genexpraw_dict[x]  
         
       absent_gene_expression=np.percentile(value_list,absent_gene_expression)
    f=open("ReactionExpression","w")     
    for gene in model.genes:
        if gene.id not in expression_dict:
           expression_dict[gene.id]=absent_gene_expression
    reaction_expression_dict={}
    for the_reaction in model.reactions:
        """if "reflection" in the_reaction.notes:
            if the_reaction.notes["reflection"] in reaction_expression_dict:
               continue"""
        if the_reaction.gene_reaction_rule != "" : 
           the_gene_reaction_relation = copy.deepcopy(the_reaction.gene_reaction_rule)
           for the_gene in the_reaction.genes: 
               the_gene_re = re.compile('(^|(?<=( |\()))%s(?=( |\)|$))'%re.escape(the_gene.id))
               the_gene_reaction_relation = the_gene_re.sub(str(expression_dict[the_gene.id]), the_gene_reaction_relation) 
           #print the_gene_reaction_relation
           expression_value=evaluate_gene_expression_string( the_gene_reaction_relation)
           reaction_expression_dict[the_reaction.id]=expression_value
           f.write(the_reaction.id+" "+str(expression_value)+"\n")
           #print(the_reaction.id+" "+str(expression_value))    
    f.close()
    return (reaction_expression_dict,genexpraw_dict)


def build_flux_penalty_dict(label_model=None,base_penalty=1,gene_expression_options_dict={},fn=None):
    flux_penalty_dict={}
    if label_model.flux_dict=={}:
       objfunc(label_model,label_model.variable_vector)
    for reaction in label_model.constrained_model.reactions:
        flux=reaction.id
        if "LABEL_RGROUP" in flux:
           flux_penalty_dict[flux]=0
        else:
           flux_penalty_dict[flux]=base_penalty
        if reaction.lower_bound<0:
           flux_penalty_dict[flux+"_reverse"]=flux_penalty_dict[flux]
    for flux in label_model.flux_dict: #just in case there is a reaction that exist as a flux but not as real reaction in constrained model
        if flux not in flux_penalty_dict:
          if "LABEL_RGROUP" in flux:
           flux_penalty_dict[flux]=0
          else:
           flux_penalty_dict[flux]=base_penalty
                 
    if gene_expression_options_dict!={}:
       flux_penalty_dict=add_gene_expression_to_flux_penalty(label_model.constrained_model,flux_penalty_dict,gene_expression_options_dict)
       """
       absent_gene_expression=gene_expression_options_dict["absent_gene_expression"]
       percentile=gene_expression_options_dict["percentile"]
       file_name=gene_expression_options_dict["file_name"]
       gene_method=gene_expression_options_dict["gene_method"]
       gene_prefix=gene_expression_options_dict["gene_prefix"]
       gene_sufix=gene_expression_options_dict["gene_sufix"]
       low_expression_threshold=gene_expression_options_dict["low_expression_threshold"]
       reaction_expression_dict,genexpraw_dict=get_gene_exp(label_model.constrained_model,absent_gene_expression=absent_gene_expression, percentile=percentile,file_name=file_name,gene_method=gene_method,gene_prefix=gene_prefix,gene_sufix=gene_sufix)
       print reaction_expression_dict
       if percentile==True: #If percentile is set to True, low_expression_threshold is assumed to be a percintile
          value_list=[] 
          for x in genexpraw_dict:
              value_list+=genexpraw_dict[x]     
          low_expression_threshold=lower_percentile=np.percentile(value_list,low_expression_threshold)
       print low_expression_threshold
       for reaction_id in  reaction_expression_dict:
           gene_expression_value=reaction_expression_dict[reaction_id]
           if gene_expression_value< low_expression_threshold:
              gene_expression_penalty=round_down(-gene_expression_value+low_expression_threshold,3)
              flux_penalty_dict[reaction_id]+=gene_expression_penalty
              if reaction_id+"_reverse" in flux_penalty_dict: 
                 flux_penalty_dict[reaction_id+"_reverse"]+=gene_expression_penalty"""
    if fn!=None:
       write_flux_penalty_dict(fn,flux_penalty_dict,irreversible_metabolic_model=label_model.irreversible_metabolic_model)
       """sheet_row_data_dict={"flux_penalty":[["Reaction id","Penalty","Reaction name","Reaction stoichiometry"]]}
       for reaction_id in sorted(flux_penalty_dict):
           if "LABEL_RGROUP_" in reaction_id:
              continue
           row=[reaction_id,flux_penalty_dict[reaction_id]]
           if label_model.irreversible_metabolic_model!=None:
             if reaction_id in label_model.irreversible_metabolic_model.reactions:
                 reaction=label_model.irreversible_metabolic_model.reactions.get_by_id(reaction_id)
                 row.append(reaction.name)
                 row.append(reaction.reaction)
           
           sheet_row_data_dict["flux_penalty"].append(row)
       write_spreadsheet(file_name=fn,sheet_row_data_dict=sheet_row_data_dict,sheet_order=None)"""
    return flux_penalty_dict

def add_gene_expression_to_flux_penalty(model,flux_penalty_dict,gene_expression_options_dict):
       absent_gene_expression=gene_expression_options_dict["absent_gene_expression"]
       percentile=gene_expression_options_dict["percentile"]
       file_name=gene_expression_options_dict["file_name"]
       gene_method=gene_expression_options_dict["gene_method"]
       gene_prefix=gene_expression_options_dict["gene_prefix"]
       gene_sufix=gene_expression_options_dict["gene_sufix"]
       low_expression_threshold=gene_expression_options_dict["low_expression_threshold"]
       reaction_expression_dict,genexpraw_dict=get_gene_exp(model,absent_gene_expression=absent_gene_expression, percentile=percentile,file_name=file_name,gene_method=gene_method,gene_prefix=gene_prefix,gene_sufix=gene_sufix)
       print reaction_expression_dict
       if percentile==True: #If percentile is set to True, low_expression_threshold is assumed to be a percintile
          value_list=[] 
          for x in genexpraw_dict:
              value_list+=genexpraw_dict[x]     
          low_expression_threshold=lower_percentile=np.percentile(value_list,low_expression_threshold)
       print low_expression_threshold
       for reaction_id in  reaction_expression_dict:
           gene_expression_value=reaction_expression_dict[reaction_id]
           if gene_expression_value< low_expression_threshold:
              gene_expression_penalty=round_down(-gene_expression_value+low_expression_threshold,3)
              flux_penalty_dict[reaction_id]+=gene_expression_penalty
              if reaction_id+"_reverse" in flux_penalty_dict: 
                 flux_penalty_dict[reaction_id+"_reverse"]+=gene_expression_penalty
       
       return flux_penalty_dict


def read_flux_penalty_dict_from_file(fn):
    flux_penalty_dict={}
    spreadsheet_dict=read_spreadsheets(file_names=fn,csv_delimiter=',',more_than_1=False,tkinter_title="Chose a file")
    flux_penalty_sheet=spreadsheet_dict[spreadsheet_dict.keys()[0]]
    for row in flux_penalty_sheet:
        try:
          flux_penalty_dict[str(row[0])]=float(row[1])
        except:
          continue
    return flux_penalty_dict


def write_flux_penalty_dict(fn,flux_penalty_dict,irreversible_metabolic_model=None):
    sheet_row_data_dict={"flux_penalty":[["Reaction id","Penalty","Reaction name","Reaction stoichiometry"]]}
    for reaction_id in sorted(flux_penalty_dict):
           if "LABEL_RGROUP_" in reaction_id:
              continue
           row=[reaction_id,flux_penalty_dict[reaction_id]]
           if irreversible_metabolic_model!=None:
             if reaction_id in irreversible_metabolic_model.reactions:
                 reaction=irreversible_metabolic_model.reactions.get_by_id(reaction_id)
                 row.append(reaction.name)
                 row.append(reaction.reaction)
           
           sheet_row_data_dict["flux_penalty"].append(row)
    write_spreadsheet(file_name=fn,sheet_row_data_dict=sheet_row_data_dict,sheet_order=None)
    
        
    
