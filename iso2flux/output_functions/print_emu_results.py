import openpyxl 
def print_emu_results(label_model,to_file=False,fn="emu_simulation",only_emu0=True):
  for condition in label_model.initial_label:
    print condition
    for size in label_model.size_emu_c_eqn_dict:
        for n,x in enumerate(label_model.condition_size_yy_dict[condition][size]):
            emu_mid=label_model.size_inverse_variable_dict[size][n]
            emu=label_model.expanded_emu_dict[emu_mid]
            if emu in label_model.emu_dict0 or only_emu0==False:
               print(str(n)+" "+emu_mid+" "+str(round(x,8)))
    if to_file==True:
       wb = openpyxl.Workbook()
       sheet = wb.active
       sheet.title = 'Fluxes'
       n=0
       for flux in label_model.flux_dict:
           n+=1
           sheet['A'+str(n)]=flux
           sheet['B'+str(n)]=str(round(label_model.flux_dict[flux],8))
       for size in label_model.size_emu_c_eqn_dict:
           sheet=wb.create_sheet(title="emu size %s"%(size))
           for n,x in enumerate(label_model.condition_size_yy_dict[condition][size]):
                sheet['A'+str(n+1)]=label_model.size_inverse_variable_dict[size][n]
                sheet['B'+str(n+1)]=str(round(x,8))
       wb.save(fn+"_"+condition+".xlsx")
