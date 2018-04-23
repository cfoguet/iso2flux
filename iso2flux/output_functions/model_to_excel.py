import openpyxl

def model_to_excel(model,name="model.xlsx"):
         
    wb = openpyxl.Workbook()
    sheet = wb.active
    sheet.title="Reactions"
    sheet['A1']="Reaction id"
    sheet['B1']="Reaction"
    sheet['C1']="Reaction name"
    sheet['D1']="Lower bound"
    sheet['E1']="Upper bound"
    sheet['F1']="Objective coefficient"
    sheet['G1']="Gene reaction rule"
    sheet['H1']="Subsystem"
    for n, reaction in enumerate(model.reactions):
        nrow=n+2
        sheet['A'+str(nrow)]=reaction.id
        sheet['B'+str(nrow)]=reaction.reaction
        sheet['C'+str(nrow)]=reaction.name
        sheet['D'+str(nrow)]=reaction.lower_bound
        sheet['E'+str(nrow)]=reaction.upper_bound
        sheet['F'+str(nrow)]=reaction.objective_coefficient
        sheet['G'+str(nrow)]=reaction.gene_reaction_rule
        sheet['H'+str(nrow)]=reaction.subsystem
        
    sheet=wb.create_sheet(title="Metabolites")
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
    
    
    
    wb.save(name)
