from ..misc.read_spreadsheets import read_spreadsheets

def read_metabolomics_data(metabolic_model,metabolite_list_fname=None):
    data_dict=read_spreadsheets(file_names=metabolite_list_fname,csv_delimiter=',',more_than_1=False,tkinter_title="Chose a file")
    #Build a name id dict
    name_id_dict={}
    for metabolite in metabolic_model.metabolites:
        if metabolite.name.lower() not in name_id_dict:
           name_id_dict[metabolite.name.lower()]=[metabolite.id]
        else:
           name_id_dict[metabolite.name.lower()].append(metabolite.id)
    metabolite_list=[]
    for sheet in data_dict:
        for row in data_dict[sheet]:
          metabolite_found=False 
          local_metabolite_list=[]
          for cell in row:
            value_list=str(cell).split(",") 
            for value in value_list:
               if value=="":
                  continue
               print value
               if value in metabolic_model.metabolites:
                  local_metabolite_list.append(value)
               elif value.lower() in name_id_dict:
                  for met_id in name_id_dict[value.lower()]:
                      local_metabolite_list.append(met_id)
          if local_metabolite_list!=[]:
             metabolite_list.append(local_metabolite_list)
    return  metabolite_list   
