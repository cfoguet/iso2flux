def model_to_excel(model,name="model.xlsx"):
    import openpyxl 
    wb = openpyxl.Workbook()
    sheet = wb.get_active_sheet()
    sheet.title = 'Metabolites'
    sheet['A1']="Metabolite id"
    sheet['B1']="Metabolite name"
    sheet['C1']="Metabolite compartment"
    sheet['D1']="Metabolite formula"
    for n, metabolite in enumerate(model.metabolites):
        nrow=n+2
        sheet['A'+str(nrow)]=metabolite.id
        sheet['B'+str(nrow)]=metabolite.name
        sheet['C'+str(nrow)]=metabolite.compartment
        #sheet['D'+str(nrow)]=metabolite.formula.id    
    
    sheet=wb.create_sheet(title="Reactions")
    sheet['A1']="Reaction id"
    sheet['B1']="Reaction name"
    sheet['C1']="Reaction subsystem"
    sheet['D1']="Reaction"
    for n, reaction in enumerate(model.reactions):
        nrow=n+2
        sheet['A'+str(nrow)]=reaction.id
        sheet['B'+str(nrow)]=reaction.name
        sheet['C'+str(nrow)]=reaction.subsystem
        sheet['D'+str(nrow)]=reaction.reaction
    
    wb.save(name)
