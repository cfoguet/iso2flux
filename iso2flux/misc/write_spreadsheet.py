import string
import openpyxl 
import csv

alphabet=list(string.ascii_uppercase)

def write_spreadsheet(file_name,sheet_row_data_dict,sheet_order=None,force_extension=False):
    if sheet_order==None or sheet_order==[]:
       sheet_order=sorted(sheet_row_data_dict)
    print file_name
    if ".xlsx" in file_name:
       wb = openpyxl.Workbook()
       for n_sheet,sheet_id in enumerate(sheet_order):
           sheet_title = (str(n_sheet)+"_"+sheet_id[:27] + '..') if len(sheet_id) > 31 else sheet_id
           if n_sheet==0:
              sheet = wb.active
              sheet.title=sheet_title
           else:
                          
              sheet=wb.create_sheet(title=sheet_title)
           for n_row, row in enumerate(sheet_row_data_dict[sheet_id]):
               for n_col,data in enumerate(sheet_row_data_dict[sheet_id][n_row]):
                    sheet[alphabet[n_col]+str(n_row+1)]=data
       wb.save(file_name)
       return
    else:
        if ".txt" not in file_name.lower() and ".csv" not in file_name.lower() and force_extension: 
           file_name+=".csv" 
        csvFile = open(file_name, 'w')
        outputWriter = csv.writer(csvFile)
        for n_sheet,sheet_id in enumerate(sheet_order):
            if len(sheet_order)>1:
               outputWriter.writerow(['###'+sheet_id+"###"])
            for n_row, row in enumerate(sheet_row_data_dict[sheet_id]):
                outputWriter.writerow(row)
            if len(sheet_order)>1: 
                outputWriter.writerow(["######"])
        csvFile.close()
        return       
 
#write_spreadsheet("essborram.csv",sheet_row_data_dict={"hola":[[1,2],[3,4]],"adeu":[["a","b"],["c","d"]]},sheet_order=["hola","adeu"])
