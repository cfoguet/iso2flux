from ..misc.read_spreadsheets import read_spreadsheets  
def read_isotoflux_settings(fn):
    sheet_rows_dict=read_spreadsheets(file_names=fn,csv_delimiter=',',more_than_1=True,tkinter_title="Chose a file")
    parameter_dict={}
    for sheet in sheet_rows_dict:
        for row in sheet_rows_dict[sheet]:
             try:
                key=row[0]
                value=str(row[1])
                print [key,value]
                if value.lower()=="true" or value.lower()=="yes":
                   value=True
                elif value.lower()=="false" or value.lower()=="no":
                   value=False
                elif "," in value.lower() or ("[" in value.lower() or "]" in value.lower()):
                   value=value.replace("]","").replace("[","").split(",")
                   if value==[""]:
                      value=[]
                   
                elif "." in value:
                   value=float(value)
                else:
                  value=int(value)
                parameter_dict[key]=value
             except:
               print "error in parameter "+key
         
    return parameter_dict
