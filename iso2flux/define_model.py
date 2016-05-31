import cobra
from cobra import Model, Reaction, Metabolite
from cobra.manipulation.modify import convert_to_irreversible
from cobra.flux_analysis.variability import flux_variability_analysis


#warning an exit for adn, dttp must be added
def showr(model):
    for x in model.reactions: 
        print (x.id+" "+x.reaction+" lb="+str(x.lower_bound)+" ub="+str(x.upper_bound))


from cobra.flux_analysis.variability import flux_variability_analysis
def write_fva(model,fraction=1,remove0=True):
    #ToDO:Export to excel
    fva=flux_variability_analysis(model,fraction_of_optimum=fraction)
    f=open("fva.csv","w")
    f.write("Name;subsystem;stoichiometry;minflux;maxflux\n")
    for x in fva:
         if abs(fva[x]["maximum"])>0.000001 or abs(fva[x]["minimum"])>0.000001 or remove0==False:
            reaction=label_model.reactions.get_by_id(x) 
            f.write(reaction.name+";"+reaction.subsystem+" ;"+reaction.reaction+";"+str(fva[x]["minimum"])+";"+str(fva[x]["maximum"])+"\n")
    f.close()

   
model = Model('metabolic_model')

#AA and extracelular metabolites
ala_c = Metabolite('ala_c',formula='C3H7NO2',name='Alanine',compartment='c')        
ala_e = Metabolite('ala_e',formula='C3H7NO2',name='Alanine',compartment='e') 
arg_c = Metabolite('arg_c',formula='',name='Arginine',compartment='c')        
arg_e = Metabolite('arg_e',formula='',name='Arginine',compartment='e') 
asn_c = Metabolite('asn_c',formula='',name='Asparagine',compartment='c')  
asn_m = Metabolite('asn_m',formula='',name='Asparagine',compartment='m')              
asn_e = Metabolite('asn_e',formula='',name='Asparagine',compartment='e')  
asp_c = Metabolite('asp_c',formula='',name='Aspartate',compartment='c')        
asp_e = Metabolite('asp_e',formula='',name='Aspartate',compartment='e') 
asp_m = Metabolite('asp_m',formula='',name='Aspartate',compartment='m') 
cys_c = Metabolite('cys_c',formula='',name='Cysteine',compartment='c')        
cys_e = Metabolite('cys_e',formula='',name='Cysteine',compartment='e') 
#Warning
glc_c = Metabolite('glc_c',formula='C6H12O6',name='Glucose',compartment='c')        
glc_e = Metabolite('glc_e',formula='C6H12O6',name='Glucose',compartment='e')  
gln_c = Metabolite('gln_c',formula='',name='Glutamine',compartment='c')        
gln_m = Metabolite('gln_m',formula='',name='Glutamine',compartment='m')        
gln_e = Metabolite('gln_e',formula='',name='Glutamine',compartment='e') 
glu_c = Metabolite('glu_c',formula='',name='Glutamate',compartment='c')        
glu_e = Metabolite('glu_e',formula='',name='Glutamate',compartment='e')  
glu_m = Metabolite('glu_m',formula='',name='Glutamate',compartment='m')  
gly_c = Metabolite('gly_c',formula='',name='Glycine',compartment='c')        
gly_e = Metabolite('gly_e',formula='',name='Glycine',compartment='e')  
gly_m = Metabolite('gly_m',formula='C2H5NO2',name='Glycine',compartment='m')
his_c = Metabolite('his_c',formula='',name='Histidne',compartment='c')        
his_e = Metabolite('his_e',formula='',name='Histidne',compartment='e')  
ile_c = Metabolite('ile_c',formula='',name='Isoleucine',compartment='c')   
ile_m = Metabolite('ile_m',formula='',name='Isoleucine',compartment='m')             
ile_e = Metabolite('ile_e',formula='',name='Isoleucine',compartment='e')
lac_c = Metabolite('lac_c',formula='C3H5O3',name='Lactate',compartment='c')        
lac_e = Metabolite('lac_e',formula='C3H5O3',name='Lactate',compartment='e')
leu_c = Metabolite('leu_c',formula='',name='Leucine',compartment='c')
leu_m = Metabolite('leu_m',formula='',name='Leucine',compartment='m')                
leu_e = Metabolite('leu_e',formula='',name='Leucine',compartment='e')
lys_c = Metabolite('lys_c',formula='',name='Lysine',compartment='c')        
lys_e = Metabolite('lys_e',formula='',name='Lysine',compartment='e')  
lys_m = Metabolite('lys_m',formula='',name='Lysine',compartment='m')  
met_c = Metabolite('met_c',formula='',name='Methionine',compartment='c')        
met_e = Metabolite('met_e',formula='',name='Methionine',compartment='e') 
phe_c = Metabolite('phe_c',formula='',name='Phenylalanine',compartment='c')        
phe_e = Metabolite('phe_e',formula='',name='Phenylalanine',compartment='e')  
pro_c = Metabolite('pro_c',formula='',name='Proline',compartment='c') 
pro_m = Metabolite('pro_m',formula='',name='Proline',compartment='m')               
pro_e = Metabolite('pro_e',formula='',name='Proline',compartment='e') 
ser_c = Metabolite('ser_c',formula='',name='Serine',compartment='c')        
ser_e = Metabolite('ser_e',formula='',name='Serine',compartment='e') 
thr_c = Metabolite('thr_c',formula='',name='Threonine',compartment='c')        
thr_e = Metabolite('thr_e',formula='',name='Threonine',compartment='e') 
trp_c = Metabolite('trp_c',formula='',name='Tryptophan',compartment='c')        
trp_e = Metabolite('trp_e',formula='',name='Tryptophan',compartment='e')
tyr_c = Metabolite('tyr_c',formula='',name='Tyrosine',compartment='c')        
tyr_e = Metabolite('tyr_e',formula='',name='Tyrosine',compartment='e') 
tyr_m = Metabolite('tyr_m',formula='',name='Tyrosine',compartment='m') 
val_c = Metabolite('val_c',formula='',name='Valine',compartment='c') 
val_m = Metabolite('val_m',formula='',name='Valine',compartment='m')               
val_e = Metabolite('val_e',formula='',name='Valine',compartment='e')   
prot_c = Metabolite('prot_e',formula='',name='Protein',compartment='c')
dna_c = Metabolite('prot_e',formula='',name='DNA',compartment='c')
rna_c = Metabolite('prot_e',formula='',name='RNA',compartment='c')
#palm_e = Metabolite('palm_e',formula='',name='Palmitate',compartment='e')
#palm_c = Metabolite('palm_c',formula='',name='Palmitate',compartment='c')
malcoa_c= Metabolite('malcoa_c',formula='',name='Malonyl coenzyme A ',compartment='c')   
#Glycolysis and Glycogen
glc6p_a_c = Metabolite('glc6p_a_c',formula='C6',name='Glucose 6-phosphate (pool A)',compartment='c')
glc1p_c = Metabolite('glc1p_c',formula='C6',name='Glucose 1-phosphate',compartment='c')
glygn_e = Metabolite('glygn_e',formula='C6',name='Glycogen',compartment='e')
udpglc_c = Metabolite('udpglc_c',formula='C6',name='Uridine diphosphate glucose',compartment='c')
glc6p_b_c = Metabolite('glc6p_b_c',formula='C6',name='Glucose 6-phosphate (pool B)' ,compartment='c')
fru6p_c = Metabolite('fru6p_c',formula='C6',name='Fructose 6-phosphate' ,compartment='c')         
fru16bp_c=Metabolite('fru16bp_c',formula='C6',name='Fructose 1-6 biphosphate', compartment='c')
gap_c=Metabolite('gap_c',formula='C3',name='Glyceraldehyde 3-phosphate',  compartment='c')
dhap_c=Metabolite('dhap_c',formula='C3',name='Dihydroxyacetone phsophate',compartment='c')
bpg13_c=Metabolite("bpg13_c",formula='C3',name='13-bisphosphoglycerate',compartment='c')
pg3_c=Metabolite("3pg_c",formula='C3',  name='3-phosphoglycerate', compartment='c')
pg2_c=Metabolite("2pg_c",formula='C3',  name='2-phosphoglycerate', compartment='c')
pep_c=Metabolite("pep_c",formula='C3',  name='Posphoenolpyruvate', compartment='c')                                                              
pyr_a_c = Metabolite('pyr_a_c',formula='C3',name='Pyruvate (pool A)',compartment='c')
pyr_b_c = Metabolite('pyr_b_c',formula='C3',name='Pyruvate (pool b)',compartment='c')
pyr_m = Metabolite('pyr_m',formula='C3',name='Pyruvate',compartment='m')

#Pentose phosphate pathway
pglclac_c=Metabolite('pglclac_c',formula='C6',name='6-phosphogluconolactone',     compartment='c')
glucon6p_c=Metabolite('glucon6p_c',formula='C6',name='D-gluconate 6-phosphate',     compartment='c')
rul5p_c=Metabolite('rul5p_c',formula='C5',name='Ribulose 5-phosphate',     compartment='c')
rib5p_c=Metabolite('rib5p_c',formula='C5',name='Ribose 5-phosphate', compartment='c')
xyl5p_c=Metabolite('xyl5p_c', formula='C5', name='Xylulose 5-phosphate',    compartment='c')
e4p_c=Metabolite('e4p_c',formula='C4',name='Erythrose 4-phosphate', compartment='c')                                                      
sed7p_c=Metabolite('sed7p_c',formula='C7', name='Sedulose 7-phosphate',  compartment='c')
#pools [glc6p_a_c,glc1p_c,udpglc_c]

#TCA and shuttle

akg_c=Metabolite('akg_c',formula='C5',name='alpha-Ketoglutarate',compartment='c')
akg_m=Metabolite('akg_m',formula='C5',name='alpha-Ketoglutarate',compartment='m')

coa_m=Metabolite('coa_m',formula='',name='Coenzyme A',compartment='m')
coa_c=Metabolite('coa_c',formula='',name='Coenzyme A',compartment='c')
nad_m=Metabolite('nad_m',formula='',name='NAD',compartment='m')
nadh_m=Metabolite('nadh_m',formula='',name='NAD',compartment='m')
accoa_c=Metabolite('accoa_c',formula='',name='acetyl-CoA',compartment='c')
accoa_m=Metabolite('accoa_m',formula='',name='acetyl-CoA',compartment='m')
cit_c=Metabolite('cit_c',formula='',name='Citrate',compartment='c')
cit_m=Metabolite('cit_m',formula='',name='Citrate',compartment='m')
icit_m=Metabolite('icit_m',formula='',name='isocitrate',compartment='m')
icit_c=Metabolite('icit_c',formula='',name='isocitrate',compartment='c')
succ_m=Metabolite('succ_m',formula='',name='alpha-Ketoglutarate',compartment='m')
succoa_m=Metabolite('succoa_m',formula='',name='alpha-Ketoglutarate',compartment='m')
mal_m=Metabolite('mal_m',formula='C4',name='L-Malate',compartment='m')
oaa_c=Metabolite('oaa_c',formula='',name='Oxaloacetate',compartment='c')
oaa_m=Metabolite('oaa_m',formula='',name='Oxaloacetate',compartment='m')
fum_m=Metabolite('fum_m',formula='C4',name='fumarate',compartment='m')
fum_c=Metabolite('fum_c',formula='C4',name='fumarate',compartment='c')

mal_m=Metabolite('mal_m',formula='C4',name='L-Malate',compartment='m')
mal_c=Metabolite('mal_c',formula='C4',name='L-Malate',compartment='c')

#Beta oxidation

hdca_c = Metabolite('hdca_c',formula='C16H31O2',name='Hexadecanoate (n-C16:0)',compartment='c')
hdca_e = Metabolite('hdca_e',formula='C16H31O2',name='Hexadecanoate (n-C16:0)',compartment='e')
pmtcoa_c= Metabolite('pmtcoa_c',formula='C37H62N7O17P3S',name='Palmitoyl-CoA (n-C16:0CoA)',compartment='c')
crn_c = Metabolite('crn_c',formula='C7H15NO3',name='L-Carnitine',compartment='c')
crn_e = Metabolite('crn_e',formula='C7H15NO3',name='L-Carnitine',compartment='e')
crn_m = Metabolite('crn_m',formula='C7H15NO3',name='L-Carnitine',compartment='m')
pmtcrn_c = Metabolite('pmtcrn_c',formula='C23H45NO4',name='L-Palmitoylcarnitine',compartment='c')
pmtcrn_m = Metabolite('pmtcrn_m',formula='C23H45NO4',name='L-Palmitoylcarnitine',compartment='m')
pmtcoa_m = Metabolite('pmtcoa_m',formula='C37H62N7O17P3S',name='Palmitoyl-CoA (n-C16:0CoA)',compartment='m')
occoa_m = Metabolite('occoa_m',formula='C29H46N7O17P3S',name='Octanoyl-CoA (n-C8:0CoA)',compartment='m')


#Energy metabolism
pi_e = Metabolite('pi_e',formula='',name='Inorganic phosphate',compartment='e')
pi_c = Metabolite('pi_c',formula='',name='Inorganic phosphate',compartment='c')
pi_m = Metabolite('pi_m',formula='',name='Inorganic phosphate',compartment='m')
ppi_c = Metabolite('ppi_c',formula='',name='Pyrophosphate',compartment='c')
                          

atp_c = Metabolite('atp_c',formula='',name='ATP',compartment='c')
#atp_novo_c = Metabolite('atp_novo_c',formula='C5',name='ATP de novo synthesis',compartment='c')                                                    
adp_c = Metabolite('adp_c',formula='',name='ADP',compartment='c')
amp_c = Metabolite('amp_c',formula='',name='AMP',compartment='c')

atp_m = Metabolite('atp_m',formula='',name='ATP',compartment='m')                          
adp_m = Metabolite('adp_m',formula='',name='ADP',compartment='m') 


utp_c = Metabolite('utp_c',formula='',name='UTP',compartment='c')                          
udp_c = Metabolite('udp_c',formula='',name='UDP',compartment='c')      

gtp_c = Metabolite('gtp_c',formula='',name='GTP',compartment='c')                          
gdp_c = Metabolite('gdp_c',formula='',name='GDP',compartment='c')
gtp_m = Metabolite('gtp_m',formula='',name='GTP',compartment='m')                          
gdp_m = Metabolite('gdp_m',formula='',name='GDP',compartment='m')


nadph_c = Metabolite('nadph_c',formula='',name='NADPH', compartment='c')
nadp_c = Metabolite('nadp_c',formula='',name='NADP', compartment='c')
nadph_m = Metabolite('nadph_m',formula='',name='NADPH', compartment='m')
nadp_m = Metabolite('nadp_m',formula='',name='NADP', compartment='m')

nadh_c = Metabolite('nadh_c',formula='',name='NADH', compartment='c')
nad_c = Metabolite('nad_c',formula='',name='NAD', compartment='c') 
nadh_e = Metabolite('nadh_c',formula='',name='NADH', compartment='e')
nad_e = Metabolite('nad_c',formula='',name='NAD', compartment='e')   

fadh2_m=Metabolite('fadh2_m',formula='',name='Flavin adenine dinucleotide reduced',compartment='m')
fad_m=Metabolite('fad_m',formula='',name='Flavin adenine dinucleotide oxidized',compartment='m')
q10_m=Metabolite('q10_m',formula='',name='Ubiquinone-10',compartment='m')
q10h2_m=Metabolite('q10h2_m',formula='Ubiquinol-10',name='Coenzyme A',compartment='m')
ficytc_m=Metabolite('ficytc_m',formula='',name='Ferricytochrome c',compartment='m')
focytc_m=Metabolite('focytc_m',formula='',name='Ferrocytochrome c',compartment='m')
etfrd_m=Metabolite('etfrd_m',formula='',name='Electron transfer flavoprotein reduced',compartment='m')
etfox_m=Metabolite('etfox_m',formula='',name='Electron transfer flavoprotein oxidized',compartment='m')


#dNTP
datp_c = Metabolite('datp_c',formula='',name='dATP',compartment='c')                          
dadp_c = Metabolite('dadp_c',formula='',name='dADP',compartment='c')
damp_c = Metabolite('damp_c',formula='',name='dAMP',compartment='c')

dttp_e = Metabolite('dttp_e',formula='',name='dTTP',compartment='e')                          
dttp_c = Metabolite('dttp_c',formula='',name='dTTP',compartment='c')                          
dtdp_c = Metabolite('dtdp_c',formula='',name='dTDP',compartment='c')   
dtmp_c = Metabolite('dtmp_c',formula='',name='dTMP',compartment='c')   

dctp_c = Metabolite('dctp_c',formula='',name='dCTP',compartment='c')                          
dcdp_c = Metabolite('dcdp_c',formula='',name='dCDP',compartment='c') 
dcmp_c = Metabolite('dcmp_c',formula='',name='dCTMP',compartment='c')

dgtp_c = Metabolite('dgtp_c',formula='',name='dGTP',compartment='c')                          
dgdp_c = Metabolite('dgdp_c',formula='',name='dGDP',compartment='c') 
dgmp_c = Metabolite('dgmp_c',formula='',name='dGMP',compartment='c') 

dutp_c = Metabolite('dutp_c',formula='',name='dUTP',compartment='c')                          
dudp_c = Metabolite('dudp_c',formula='',name='dUDP',compartment='c') 
dump_c = Metabolite('dump_c',formula='',name='dUMP',compartment='c') 

#Fatty Acid Biosynthesis 
chol_c = Metabolite('chol_c',formula='C5H14NO',name='Choline',compartment='c')
chol_e = Metabolite('chol_e',formula='C5H14NO',name='Choline',compartment='e')
cholp_c = Metabolite('cholp_c',formula='C5H13NO4P',name='Choline phosphate',compartment='c')
dag_hs_c = Metabolite('dag_hs_c',formula='C3H6OFULLRCO2FULLR2CO2',name='diacylglycerol (homo sapiens)',compartment='c')
pchol_hs_c=Metabolite('pchol_hs_c',formula='C8H18NO4PFULLRCO2FULLR2CO2',name='Phosphatidylcholine (homo sapiens)',compartment='c')
pchol_hs_e=Metabolite('pchol_hs_e',formula='C8H18NO4PFULLRCO2FULLR2CO2',name='Phosphatidylcholine (homo sapiens)',compartment='e')
cdpchol_c = Metabolite('cdpchol_c',formula='C14H25N4O11P2',name='CDPcholine',compartment='c')
rtotalcoa_c = Metabolite('rtotalcoa_c',formula='CO2FULLRC21H31N7O15P3S',name='R total Coenzyme A',compartment='c')
glyc3p_c = Metabolite('glyc3p_c',formula='C3H7O6P',name='Glycerol 3-phosphate',compartment='c')
alpa_hs_c = Metabolite('alpa_hs_c',formula='C3H6O5PFULLRCO2',name='lysophosphatidic acid (homo sapiens)',compartment='c')
odecoa_c = Metabolite('odecoa_c',formula='C39H64N7O17P3S',name='Octadecenoyl-CoA (n-C18:1CoA)',compartment='c')
stcoa_c = Metabolite('stcoa_c',formula='C39H66N7O17P3S',name='Stearoyl-CoA (n-C18:0CoA)',compartment='c')
hdcoa_c = Metabolite('hdcoa_c',formula='C37H60N7O17P3S',name='Hexadecenoyl-CoA (n-C16:1CoA)',compartment='c')
pa_hs_c = Metabolite('pa_hs_c',formula='C3H5O4PFULLRCO2FULLR2CO2',name='phosphatidic acid (homo sapiens)',compartment='c')
dag_hs_c = Metabolite('dag_hs_c',formula='C3H6OFULLRCO2FULLR2CO2',name='diacylglycerol (homo sapiens)',compartment='c')
ps_hs_c = Metabolite('ps_hs_c',formula='C6H11NO6PFULLRCO2FULLR2CO2',name='phosphatidylserine (homo sapiens)',compartment='c')
ps_hs_m = Metabolite('ps_hs_m',formula='C6H11NO6PFULLRCO2FULLR2CO2',name='phosphatidylserine (homo sapiens)',compartment='m')
pe_hs_c = Metabolite('pe_hs_c',formula='C5H12NO4PFULLRCO2FULLR2CO2',name='phosphatidylethanolamine (homo sapiens)',compartment='c')
pe_hs_m = Metabolite('pe_hs_m',formula='C5H12NO4PFULLRCO2FULLR2CO2',name='phosphatidylethanolamine (homo sapiens)',compartment='c')
pe_hs_e = Metabolite('pe_hs_e',formula='C5H12NO4PFULLRCO2FULLR2CO2',name='phosphatidylethanolamine (homo sapiens)',compartment='e')
ocdcea_e = Metabolite('ocdcea_e',formula='C18H33O2',name='octadecenoate (n-C18:1)',compartment='e')
ocdcea_c = Metabolite('ocdcea_c',formula='C18H33O2',name='octadecenoate (n-C18:1)',compartment='c')
ocdca_c = Metabolite('ocdca_c',formula='C18H35O2',name='octadecanoate (n-C18:0)',compartment='c')
ocdca_e = Metabolite('ocdca_e',formula='C18H35O2',name='octadecanoate (n-C18:0)',compartment='e')
#AA metabolism
p3hp_c = Metabolite('p3hp_c',formula='C3',name='3-phospho-hydroxypyruvate', compartment='c')
pser_c=Metabolite('pser_c',formula='C3',name='3-phospho-L-serine', compartment='c')
hpyr_c=Metabolite('hpyr_c',formula='C3',name='Hydroxy-pyruvate', compartment='c')
glyc_c=Metabolite('glyc_c',formula='C3',name='Glycerate', compartment='c')

amet_c=Metabolite('amet_c',formula='C15H23N6O5S',name='S-Adenosyl-L-methionine', compartment='c')
ahcys_c=Metabolite('ahcys_c',formula='C14H20N6O5S',name='"S-Adenosyl-L-homocysteine"', compartment='c')
methyl_a_c=Metabolite('methyl_a_c',formula='',name='demethylated acceptor', compartment='c')
ch3_methyl_a_c=Metabolite('ch3_methyl_a_c',formula='',name='methylated acceptor', compartment='c')
ch3_methyl_a_c=Metabolite('ch3_methyl_a_c',formula='',name='methylated acceptor', compartment='c')
hcys_c=Metabolite('hcys_c',formula='C4H9NO2S',name='Homocysteine', compartment='c')
adn_c=Metabolite('adn_c',formula='C10H13N5O4',name='Adenosine', compartment='c')
adn_e=Metabolite('adn_e',formula='C10H13N5O4',name='Adenosine', compartment='e')
ade_c = Metabolite('ade_c',formula='C5H5N5',name='Adenine',compartment='c')
ade_e = Metabolite('ade_e',formula='C5H5N5',name='Adenine',compartment='e')

glu5p_m=Metabolite('glu5p_m',formula='',name='L-Glutamate 5-phosphate', compartment='m')
glu5sa_m=Metabolite('glu5sa_m',formula='',name='L-Glutamate 5-semialdehyde', compartment='m')
n1pyr5c_m=Metabolite('1pyr5c_m',formula='',name='1-Pyrroline-5-carboxylate', compartment='m')

orn_c=Metabolite('orn_c',formula='',name='Ornithine', compartment='c')
orn_m=Metabolite('orn_m',formula='',name='Ornithine', compartment='m')
urea_c=Metabolite('urea_c',formula='',name='Urea', compartment='c')
urea_e=Metabolite('urea_e',formula='',name='Urea', compartment='e')

urcan_c=Metabolite('urcan_c',formula='C6H5N2O2',name='Urocanate', compartment='c')
n4izp_c=Metabolite('4izp_c',formula='4-Imidazolone-5-propanoate',name='Urea', compartment='c')
forglu_c=Metabolite('forglu_c',formula="C6H8N2O4",name='N-Formimidoyl-L-glutamate', compartment='c')
ptrc_c=Metabolite('ptrc_c',formula='',name='Putrescine', compartment='c')
ametam_c = Metabolite('ametam_c',formula='C14H24N6O3S',name='S-Adenosylmethioninamine',compartment='c')
n5mdru1p_c = Metabolite('5mdru1p_c',formula='C6H11O7PS',name='5-Methylthio-5-deoxy-D-ribulose 1-phosphate',compartment='c')
n5mta_c = Metabolite('5mta_c',formula='C11H15N5O3S',name='5-Methylthioadenosine',compartment='c')
n5mdr1p_c = Metabolite('5mdr1p_c',formula='C6H11O7PS',name='5-Methylthio-5-deoxy-D-ribose 1-phosphate',compartment='c')
dkmpp_c = Metabolite('dkmpp_c',formula='C6H9O6PS',name='2,3-diketo-5-methylthio-1-phosphopentane',compartment='c')
n2kmb_c = Metabolite('2kmb_c',formula='C5H7O3S',name='2-keto-4-methylthiobutyrate',compartment='c')
for_c = Metabolite('for_c',formula='CH1O2',name='Formate',compartment='c')

spmd_c = Metabolite("spmd_c",formula='C7H22N3',name='Spermidine',compartment='c')
sprm_c = Metabolite('sprm_c',formula='C10H30N4',name='Spermine',compartment='c')

n3mop_m = Metabolite('3mop_m',formula='C6H9O3',name='(S)-3-Methyl-2 oxopentanoate',compartment='m')
n2mbcoa_m = Metabolite('2mbcoa_m',formula='C26H40N7O17P3S',name='2-Methylbutanoyl-CoA',compartment='m')
n2mb2coa_m = Metabolite('2mb2coa_m',formula='C26H38N7O17P3S',name='trans-2-Methylbut-2-enoyl-CoA',compartment='m')
n3hmbcoa_m = Metabolite('3hmbcoa_m',formula='C26H40N7O18P3S',name='(S)-3-Hydroxy-2-methylbutyryl-CoA',compartment='m')
n2maacoa_m = Metabolite('2maacoa_m',formula='C26H38N7O18P3S',name='2-Methyl-3-acetoacetyl-CoA',compartment='m')
ppcoa_c = Metabolite('ppcoa_c',formula='C24H36N7O17P3S',name='Propanoyl-CoA',compartment='c')
ppcoa_m = Metabolite('ppcoa_m',formula='C24H36N7O17P3S',name='Propanoyl-CoA',compartment='m')
mmcoa_s_m = Metabolite('mmcoa_s_m',formula='C25H35N7O19P3S',name='(S)-Methylmalonyl-CoA',compartment='m')
mmcoa_r_m = Metabolite('mmcoa_r_m',formula='C25H35N7O19P3S',name='(R)-Methylmalonyl-CoA',compartment='m')

n4mop_m = Metabolite('4mop_m',formula='C6H9O3',name='4-Methyl-2-oxopentanoate',compartment='m')
ivcoa_m = Metabolite('ivcoa_m',formula='C26H40N7O17P3S',name='Isovaleryl-CoA',compartment='m')
n3mgcoa_m = Metabolite('3mgcoa_m',formula='C27H37N7O19P3S',name='3-Methylglutaconyl-CoA',compartment='m')
n3mb2coa_m = Metabolite('3mb2coa_m',formula='C26H38N7O17P3S',name='3-Methylbut-2-enoyl-CoA',compartment='m')
hmgcoa_m = Metabolite('hmgcoa_m',formula='C27H39N7O20P3S',name='Hydroxymethylglutaryl-CoA',compartment='m')
acac_m = Metabolite('acac_m',formula='C4H5O3',name='Acetoacetate',compartment='m')
aacoa_m = Metabolite('aacoa_m',formula='C25H36N7O18P3S',name='Acetoacetyl-CoA',compartment='m')

saccrp_m = Metabolite('saccrp_m',formula='C11H19N2O6',name='L-Saccharopine',compartment='m')
nl2aadp6sa_m = Metabolite('l2aadp6sa_m',formula='C6H11NO3',name='L-2-Aminoadipate 6-semialdehyde',compartment='m')
nl2aadp_m = Metabolite('l2aadp_m',formula='C6H10NO4',name='L-2-Aminoadipate',compartment='m')
nl2aadp_c = Metabolite('l2aadp_c',formula='C6H10NO4',name='L-2-Aminoadipate',compartment='c')
n2oxoadp_c = Metabolite('2oxoadp_c',formula='C6H6O5',name='2-Oxoadipate',compartment='c')
n2oxoadp_m = Metabolite('2oxoadp_m',formula='C6H6O5',name='2-Oxoadipate',compartment='m')
glutcoa_m = Metabolite('glutcoa_m',formula='C26H37N7O19P3S',name='Glutaryl-CoA',compartment='m')
b2coa_m = Metabolite('b2coa_m',formula='C25H36N7O17P3S',name='Crotonoyl-CoA',compartment='m')
n3hbcoa_m = Metabolite('3hbcoa_m',formula='C25H38N7O18P3S',name='(S)-3-Hydroxybutanoyl-CoA',compartment='m')

n34hpp_c=Metabolite('34hpp_c',formula='C9H7O4',name='3-(4-Hydroxyphenyl)pyruvate',compartment='c')
hgentis_c= Metabolite('hgentis_c',formula='C8H7O4',name='Homogentisate',compartment='c')
n4mlacac_c= Metabolite('4mlacac_c',formula='C8H6O6',name='4-Maleylacetoacetate',compartment='c')
n4fumacac_c= Metabolite('4fumacac_c',formula='C8H6O6',name='4-Fumarylacetoacetate',compartment='c')
acac_c = Metabolite('acac_c',formula='C4H5O3',name='Acetoacetate',compartment='c')

n2obut_c = Metabolite('2obut_c',formula='C4H5O3',name='2-Oxobutanoate',compartment='c')

lfmkynr_c = Metabolite('lfmkynr_c',formula='C11H12N2O4',name='L-Formylkynurenine',compartment='c')
lkynr_c = Metabolite('lkynr_c',formula='C10H12N2O3',name='L-Kynurenine',compartment='c')
hlkynr_c = Metabolite('hlkynr_c',formula='C10H12N2O4',name='3-Hydroxy-L-kynurenine',compartment='c')
n3hanthrn_c = Metabolite('3hanthrn_c',formula='C7H7NO3',name='3-Hydroxyanthranilate',compartment='c')
cmusa_c = Metabolite('cmusa_c',formula='C7H6NO5',name='2-Amino-3-carboxymuconate semialdehyde',compartment='c')
am6sa_c = Metabolite('am6sa_c',formula='C6H7NO3',name='2-Aminomuconate 6-semialdehyde',compartment='c')
amuco_c = Metabolite('amuco_c',formula='C6H6NO4',name='2-Aminomuconate',compartment='c')

n3sala_c = Metabolite('3sala_c',formula='C3H5NO4S',name='3-Sulfino-L-alanine',compartment='c')
n3snpyr_c = Metabolite('3snpyr_c',formula='C3H2O5S',name='3-Sulfinopyruvate',compartment='c')
so3_c = Metabolite('so3_c',formula='O3S',name='Sulfite',compartment='c')
so4_e = Metabolite('so4_e',formula='O4S',name='Sulfate',compartment='e')
so4_c = Metabolite('so4_c',formula='O4S',name='Sulfate',compartment='c')

cyst_c= Metabolite('cyst_dash_l_c',formula='C7H14N2O4S',name='L-Cystathionine',compartment='c')

n2amac_c = Metabolite('2amac_c',formula='C3H5NO2',name='2-Aminoacrylate',compartment='c')
#folate metabolism: 
thf_c = Metabolite('thf_c',formula='',name='5,6,7,8-Tetrahydrofolate', compartment='c')
mlthf_c=Metabolite('mlthf_c',formula='C3',name='5,10-methylenetetrahydrofolate', compartment='c')
dhf_c = Metabolite('dhf_c',formula='',name="7,8-Dihydrofolate", compartment='c')
m5thf_c = Metabolite('5mthf_c',formula='',name="5-Methyltetrahydrofolate", compartment='c')
hcys_c = Metabolite('hcys_c',formula='',name="L-Homocysteine", compartment='c')

n3mob_c = Metabolite('3mob_c',formula='C5H7O3',name='3-Methyl-2-oxobutanoate',compartment='c')
n3mob_m = Metabolite('3mob_m',formula='C5H7O3',name='3-Methyl-2-oxobutanoate',compartment='m')
ibcoa_m = Metabolite('ibcoa_m',formula='C25H38N7O17P3S',name='Isobutyryl-CoA',compartment='m')
n2mp2coa_m = Metabolite('2mp2coa_m',formula='C25H36N7O17P3S',name='2-Methylprop-2-enoyl-CoA',compartment='m')
n3hibutcoa_m = Metabolite('3hibutcoa_m',formula='C25H38N7O18P3S',name='(S)-3-Hydroxyisobutyryl-CoA',compartment='m')
n3hmp_m = Metabolite('3hmp_m',formula='C4H7O3',name='3-Hydroxy-2-methylpropanoate',compartment='m')
n2mop_m = Metabolite('2mop_m',formula='C4H5O3',name='2-Methyl-3-oxopropanoate',compartment='m')

thbpt_c = Metabolite('thbpt_c',formula='C9H15N5O3',name='Tetrahydrobiopterin',compartment='c')
thbpt4acam_c = Metabolite('thbpt4acam_c',formula='C9H15N5O4',name='Tetrahydrobiopterin-4a-carbinolamine',compartment='c')
dhbpt_c = Metabolite('dhbpt_c',formula='C9H13N5O3',name='6,7-Dihydrobiopterin',compartment='c')
n5forthf_c = Metabolite('5forthf_c',formula='C20H22N8O6',name='5-Formiminotetrahydrofolate',compartment='c')
methf_c = Metabolite('methf_c',formula='C20H20N7O6',name='5,10-Methenyltetrahydrofolate',compartment='c')
n10fthf_c = Metabolite('10fthf_c',formula='C20H21N7O7',name='10-Formyltetrahydrofolate',compartment='c')


#Misc
nh4_e=Metabolite('nh4_e',formula='NH4',name='NH4', compartment='e')
nh4_c=Metabolite('nh4_c',formula='NH4',name='NH4', compartment='c')
nh4_m=Metabolite('nh4_m',formula='NH4',name='NH4', compartment='m')

o2_e=Metabolite('o2_e',formula='O2',name='O2', compartment='e')
o2_c=Metabolite('o2_c',formula='O2',name='O2', compartment='c')
o2_m=Metabolite('o2_m',formula='O2',name='O2', compartment='m')

co2_e=Metabolite('co2_e',formula='CO2',name='CO2', compartment='e')
co2_c=Metabolite('co2_c',formula='CO2',name='CO2', compartment='c')
co2_m=Metabolite('co2_m',formula='CO2',name='CO2', compartment='m')

h2o_e=Metabolite('h2o_e',formula='H2O',name='Water', compartment='e')
h2o_c=Metabolite('h2o_c',formula='H2O',name='Water', compartment='c')
h2o_m=Metabolite('h2o_m',formula='H2O',name='Water', compartment='m')

h_e=Metabolite('h_e',formula='H',name='Proton', compartment='e')
h_c=Metabolite('h_c',formula='H',name='Proton', compartment='c')
h_m=Metabolite('h_m',formula='H',name='Proton', compartment='m')

hco3_e=Metabolite('hco3_e',formula='HCO3',name='HCO3', compartment='e')
hco3_c=Metabolite('hco3_c',formula='HCO3',name='HCO3', compartment='c')
hco3_m=Metabolite('hco3_m',formula='HCO3',name='HCO3', compartment='m')

#Nucleotide Biosynthesis 
prpp_c = Metabolite('prpp_c',formula='C5H8O14P3',name='5-Phospho-alpha-D-ribose 1-diphosphate',compartment='c')
pram_c = Metabolite('pram_c',formula='C5H11NO7P',name='5-Phospho-beta-D-ribosylamine',compartment='c')
gar_c = Metabolite('gar_c',formula='C7H14N2O8P',name='N1-(5-Phospho-D-ribosyl)glycinamide',compartment='c')
fgam_c = Metabolite('fgam_c',formula='C8H13N2O9P',name='N2-Formyl-N1-(5-phospho-D-ribosyl)glycinamide',compartment='c')
fpram_c = Metabolite('fpram_c',formula='C8H15N3O8P',name='2-(Formamido)-N1-(5-phospho-D-ribosyl)acetamidine',compartment='c')
air_c = Metabolite('air_c',formula='C8H12N3O7P',name='5-amino-1-(5-phospho-D-ribosyl)imidazole',compartment='c')
n5aizc_c = Metabolite('5aizc_c',formula='C9H11N3O9P',name='5-amino-1-(5-phospho-D-ribosyl)imidazole-4-carboxylate',compartment='c')
n25aics_c = Metabolite('25aics_c',formula='C13H15N4O12P',name='(S)-2-[5-Amino-1-(5-phospho-D-ribosyl)imidazole-4-carboxamido]succinate',compartment='c')
aicar_c = Metabolite('aicar_c',formula='C9H13N4O8P',name='5-Amino-1-(5-Phospho-D-ribosyl)imidazole-4-carboxamide',compartment='c')
fprica_c = Metabolite('fprica_c',formula='C10H13N4O9P',name='5-Formamido-1-(5-phospho-D-ribosyl)imidazole-4-carboxamide',compartment='c')
imp_c = Metabolite('imp_c',formula='C10H11N4O8P',name='IMP',compartment='c')
dcamp_c = Metabolite('dcamp_c',formula='C14H14N5O11P',name='N6-(1,2-Dicarboxyethyl)-AMP',compartment='c')
xmp_c = Metabolite('xmp_c',formula='C10H11N4O9P',name='Xanthosine 5-phosphate',compartment='c')
gmp_c = Metabolite('gmp_c',formula='C10H12N5O8P',name='GMP',compartment='c')
ctp_c = Metabolite('ctp_c',formula='C9H12N3O14P3',name='CTP',compartment='c')
cdp_c = Metabolite('cdp_c',formula='C9H12N3O11P2',name='CDP',compartment='c')
cmp_c = Metabolite('cmp_c',formula='C9H12N3O8P',name='CMP',compartment='c')
ump_c = Metabolite('ump_c',formula='C9H11N2O9P',name='UMP',compartment='c')
dttp_c = Metabolite('dttp_c',formula='C10H13N2O14P3',name='dTTP',compartment='c')
dttp_m = Metabolite('dttp_m',formula='C10H13N2O14P3',name='dTTP',compartment='m')
dttp_n = Metabolite('dttp_n',formula='C10H13N2O14P3',name='dTTP',compartment='n')
dudp_c = Metabolite('dudp_c',formula='C9H11N2O11P2',name='dUDP',compartment='c')
dudp_m = Metabolite('dudp_m',formula='C9H11N2O11P2',name='dUDP',compartment='m')
dudp_n = Metabolite('dudp_n',formula='C9H11N2O11P2',name='dUDP',compartment='n')
dump_c = Metabolite('dump_c',formula='C9H11N2O8P',name='dUMP',compartment='c')
dump_m = Metabolite('dump_m',formula='C9H11N2O8P',name='dUMP',compartment='m')
dump_n = Metabolite('dump_n',formula='C9H11N2O8P',name='dUMP',compartment='n')
duri_c = Metabolite('duri_c',formula='C9H12N2O5',name='Deoxyuridine',compartment='c')
duri_e = Metabolite('duri_e',formula='C9H12N2O5',name='Deoxyuridine',compartment='e')
duri_m = Metabolite('duri_m',formula='C9H12N2O5',name='Deoxyuridine',compartment='m')
duri_n = Metabolite('duri_n',formula='C9H12N2O5',name='Deoxyuridine',compartment='n')
dutp_c = Metabolite('dutp_c',formula='C9H11N2O14P3',name='dUTP',compartment='c')
dutp_m = Metabolite('dutp_m',formula='C9H11N2O14P3',name='dUTP',compartment='m')
dutp_n = Metabolite('dutp_n',formula='C9H11N2O14P3',name='dUTP',compartment='n')
trdox_c = Metabolite('trdox_c',formula='X',name='Oxidized thioredoxin',compartment='c')
trdrd_c = Metabolite('trdrd_c',formula='XH2',name='Reduced thioredoxin',compartment='c')
orot5p_c = Metabolite('orot5p_c',formula='C10H10N2O11P',name='Orotidine 5-phosphate',compartment='c')
orot_c = Metabolite('orot_c',formula='C5H3N2O4',name='Orotate',compartment='c')
cbasp_c = Metabolite('cbasp_c',formula='C5H6N2O5',name='N-Carbamoyl-L-aspartate',compartment='c')
dhor_c = Metabolite('dhor_c',formula='C5H5N2O4',name='(S)-Dihydroorotate',compartment='c')
uri_c = Metabolite('uri_c',formula='C9H12N2O6',name='Uridine',compartment='c')
cbp_c = Metabolite('cbp_c',formula='CH2NO5P',name='Carbamoyl phosphate',compartment='c')

#Cholesterol Biosynthesis 
hmgcoa_x = Metabolite('hmgcoa_x',formula='C27H39N7O20P3S',name='Hydroxymethylglutaryl-CoA',compartment='x')
h_x = Metabolite('h_x',formula='H',name='H+',compartment='x')
h_r = Metabolite('h_r',formula='H',name='H+',compartment='r')
coa_x = Metabolite('coa_x',formula='C21H32N7O16P3S',name='Coenzyme A',compartment='x')
nadp_x = Metabolite('nadp_x',formula='C21H25N7O17P3',name='Nicotinamide adenine dinucleotide phosphate',compartment='x')
nadph_x = Metabolite('nadph_x',formula='C21H26N7O17P3',name='Nicotinamide adenine dinucleotide phosphate - reduced',compartment='x')
mev_r_x = Metabolite('mev_dash_r_x',formula='C6H11O4',name='(R)-Mevalonate',compartment='x')
n5pmev_x = Metabolite('5pmev_x',formula='C6H10O7P',name='(R)-5-Phosphomevalonate',compartment='x')
adp_x = Metabolite('adp_x',formula='C10H12N5O10P2',name='ADP',compartment='x')
atp_x = Metabolite('atp_x',formula='C10H12N5O13P3',name='ATP',compartment='x')
n5dpmev_x = Metabolite('5dpmev_x',formula='C6H10O10P2',name='(R)-5-Diphosphomevalonate',compartment='x')
co2_x = Metabolite('co2_x',formula='CO2',name='CO2',compartment='x')
co2_r = Metabolite('co2_r',formula='CO2',name='CO2',compartment='r')
pi_x = Metabolite('pi_x',formula='HO4P',name='Phosphate',compartment='x')
ppi_x = Metabolite('ppi_x',formula='HO7P2',name='Diphosphate',compartment='x')
ppi_r = Metabolite('ppi_r',formula='HO7P2',name='Diphosphate',compartment='r')
ipdp_x = Metabolite('ipdp_x',formula='C5H9O7P2',name='Isopentenyl diphosphate',compartment='x')
dmpp_x = Metabolite('dmpp_x',formula='C5H9O7P2',name='Dimethylallyl diphosphate',compartment='x')
grdp_x = Metabolite('grdp_x',formula='C10H17O7P2',name='Geranyl diphosphate',compartment='x')
frdp_x = Metabolite('frdp_x',formula='C15H25O7P2',name='Farnesyl diphosphate',compartment='x')
frdp_r = Metabolite('frdp_r',formula='C15H25O7P2',name='Farnesyl diphosphate',compartment='r')
nadp_r = Metabolite('nadp_r',formula='C21H25N7O17P3',name='Nicotinamide adenine dinucleotide phosphate',compartment='r')
nadph_r = Metabolite('nadph_r',formula='C21H26N7O17P3',name='Nicotinamide adenine dinucleotide phosphate - reduced',compartment='r')
ssq23epx_r = Metabolite('ssq23epx_r',formula='C30H50O',name='(S)-Squalene-2,3-epoxide',compartment='r')
o2_r = Metabolite('o2_r',formula='O2',name='O2',compartment='r')
h2o_r = Metabolite('h2o_r',formula='H2O',name='H2O',compartment='r')
sql_r = Metabolite('sql_r',formula='C30H50',name='Squalene',compartment='r')
lanost_r = Metabolite('lanost_r',formula='C30H50O',name='Lanosterol',compartment='r')
n44mctr_r = Metabolite('44mctr_r',formula='C29H46O',name='4,4-dimethylcholesta-8,14,24-trienol',compartment='r')
for_r = Metabolite('for_r',formula='CH1O2',name='Formate',compartment='r')
n44mzym_r = Metabolite('44mzym_r',formula='C29H48O',name='4,4-dimethylzymosterol',compartment='r')
n4mzym_int1_r = Metabolite('4mzym_int1_r',formula='C29H46O3',name='4-Methylzymosterol intermediate 1',compartment='r')
n4mzym_int2_r = Metabolite('4mzym_int2_r',formula='C28H44O',name='4-Methylzymosterol intermediate 2',compartment='r')
zym_int2_r = Metabolite('zym_int2_r',formula='C27H42O',name='zymosterol intermediate 2',compartment='r')
zymst_r = Metabolite('zymst_r',formula='C27H44O',name='zymosterol',compartment='r')
zymstnl_r = Metabolite('zymstnl_r',formula='C27H46O',name='Zymostenol',compartment='r')
lthstrl_r = Metabolite('lthstrl_r',formula='C27H46O',name='5alpha-Cholest-7-en-3beta-ol',compartment='r')
n7dhchsterol_r = Metabolite('7dhchsterol_r',formula='C27H44O',name='7-Dehydrocholesterol',compartment='r')
chsterol_r = Metabolite('chsterol_r',formula='C27H46O',name='Cholesterol',compartment='r')
chsterol_e = Metabolite('chsterol_e',formula='C27H46O',name='Cholesterol',compartment='e')
chsterol_m = Metabolite('chsterol_m',formula='C27H46O',name='Cholesterol',compartment='m')
chsterol_c = Metabolite('chsterol_c',formula='C27H46O',name='Cholesterol',compartment='c')
hmgcoa_c = Metabolite('hmgcoa_c',formula='C27H39N7O20P3S',name='Hydroxymethylglutaryl-CoA',compartment='c')
aacoa_c = Metabolite('aacoa_c',formula='C25H36N7O18P3S',name='Acetoacetyl-CoA',compartment='c')
nad_r = Metabolite('nad_r',formula='C21H26N7O14P2',name='Nicotinamide adenine dinucleotide',compartment='r')
nadh_r = Metabolite('nadh_r',formula='C21H27N7O14P2',name='Nicotinamide adenine dinucleotide - reduced',compartment='r')
 #Inositol Phosphate synthesis 
cdpdag_hs_c = Metabolite('cdpdag_hs_c',formula='C12H17N3O11P2FULLRCO2FULLR2CO2',name='CDP diacylglycerol (homo sapiens)',compartment='c')
mi1p_d_c = Metabolite('mi1p_d_c',formula='C6H11O9P',name='1D-myo-Inositol 1-phosphate',compartment='c')
inost_c = Metabolite('inost_c',formula='C6H12O6',name='myo-Inositol',compartment='c')
pail_hs_c = Metabolite('pail_hs_c',formula='C9H16O9PFULLRCO2FULLR2CO2',name='phosphatidylinositol (homo sapiens)',compartment='c')
####### clpn_hs_c and pglyc_hs_c sysnthesis
pgp_hs_c = Metabolite('pgp_hs_c',formula='C6H11O9P2FULLRCO2FULLR2CO2',name='phosphatidyl glycerol phosphate (homo sapiens)',compartment='c')
pglyc_hs_c = Metabolite('pglyc_hs_c',formula='C6H12O6PFULLRCO2FULLR2CO2',name='phosphatidylglycerol (homo sapiens)',compartment='c')
clpn_hs_c = Metabolite('clpn_hs_c',formula='C9H16O9P2FULLRCO2FULLR2CO2FULLRCO2FULLR2CO2',name='cardiolipin (homo sapiens)',compartment='c')


########### Sphingomyelin synthesis
n3dsphgn_c = Metabolite('3dsphgn_c',formula='C18H38NO2',name='3-Dehydrosphinganine',compartment='c')
sphgn_c = Metabolite('sphgn_c',formula='C18H40NO2',name='Sphinganine',compartment='c')
dhcrm_hs_c = Metabolite('dhcrm_hs_c',formula='C18H38NO2FULLRCO',name='dihydroceramide (homo sapiens)',compartment='c')
crm_hs_c = Metabolite('crm_hs_c',formula='C18H36NO2FULLRCO',name='ceramide (homo sapiens)',compartment='c')
sphmyln_hs_c = Metabolite('sphmyln_hs_c',formula='C23H48N2O5PFULLRCO',name='sphingomyelin (homo sapiens)',compartment='c')

#Urea Cycle 
cbp_m = Metabolite('cbp_m',formula='CH2NO5P',name='Carbamoyl phosphate',compartment='m')
citr_c = Metabolite('citr_c',formula='C6H13N3O3',name='L-Citrulline',compartment='c')
citr_m = Metabolite('citr_m',formula='C6H13N3O3',name='L-Citrulline',compartment='m')
argsuc_c = Metabolite('argsuc_c',formula='C10H17N4O6',name='N(omega)-(L-Arginino)succinate',compartment='c')


#pools [glc6p_a_c,glc1p_c,udpglc_c,glucon6p_c, pglclac_c]
#pools[pyr_b_c,pyr_m]
      
#Transport reactions 
reaction = Reaction('o2_etr')
reaction.name = 'O2 extracellular transport'
reaction.subsystem = 'misc'
reaction.lower_bound = -1000  
reaction.upper_bound = 1000.  
reaction.add_metabolites({o2_e:-1, o2_c:1}) 
model.add_reaction(reaction)

reaction = Reaction('o2_mtr')
reaction.name = 'O2 mitochondrial transport'
reaction.subsystem = 'misc'
reaction.lower_bound = -1000  
reaction.upper_bound = 1000.  
reaction.add_metabolites({o2_c:-1, o2_m:1}) 
model.add_reaction(reaction)

reaction = Reaction('co2_etr')
reaction.name = 'CO2 extracellular transport'
reaction.subsystem = 'misc'
reaction.lower_bound = -1000  
reaction.upper_bound = 1000.  
reaction.add_metabolites({co2_e:-1, co2_c:1}) 
model.add_reaction(reaction)

reaction = Reaction('so4_etr')
reaction.name = 'SO4 extracellular transport'
reaction.subsystem = 'misc'
reaction.lower_bound = -1000  
reaction.upper_bound = 1000.  
reaction.add_metabolites({so4_e:-1, so4_c:1}) 
model.add_reaction(reaction)

reaction = Reaction('co2_mtr')
reaction.name = 'CO2 mitochondrial transport'
reaction.subsystem = 'misc'
reaction.lower_bound = -1000  
reaction.upper_bound = 1000.  
reaction.add_metabolites({co2_c:-1, co2_m:1}) 
model.add_reaction(reaction)

reaction = Reaction('h2o_etr')
reaction.name = 'H2O extracellular transport'
reaction.subsystem = 'misc'
reaction.lower_bound = -1000  
reaction.upper_bound = 1000.  
reaction.add_metabolites({h2o_e:-1, h2o_c:1}) 
model.add_reaction(reaction)

reaction = Reaction('h2o_mtr')
reaction.name = 'H2O mitochondrial transport'
reaction.subsystem = 'misc'
reaction.lower_bound = -1000  
reaction.upper_bound = 1000.  
reaction.add_metabolites({h2o_c:-1, h2o_m:1}) 
model.add_reaction(reaction)

reaction = Reaction('h_etr')
reaction.name = 'Proton extracellular transport'
reaction.subsystem = 'misc'
reaction.lower_bound = -1000  
reaction.upper_bound = 1000.  
reaction.add_metabolites({h_e:-1, h_c:1}) 
model.add_reaction(reaction)

reaction = Reaction('nh4_etr')
reaction.name = 'NH4 extracellular transport'
reaction.subsystem = 'misc'
reaction.lower_bound = -1000  
reaction.upper_bound = 1000.  
reaction.add_metabolites({nh4_e:-1, nh4_c:1}) 
model.add_reaction(reaction)

reaction = Reaction('nh4_mtr')
reaction.name = 'NH4 mitochondrial transport'
reaction.subsystem = 'misc'
reaction.lower_bound = -1000  
reaction.upper_bound = 1000.  
reaction.add_metabolites({nh4_c:-1, nh4_m:1}) #incomplete stochiometry
model.add_reaction(reaction)

reaction = Reaction('lysmtr')
reaction.name = 'Lysine mitochondrial transport via ornithine carrier'
reaction.subsystem = 'Transport, Mitochondrial'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({lys_c: -1.0, h_m: -1.0, h_c: 1.0, lys_m: 1.0})
model.add_reaction(reaction)

reaction = Reaction('ala_etr')
reaction.name = 'Alanine extracellular transport'
reaction.subsystem = 'AA'
reaction.lower_bound = -1000  
reaction.upper_bound = 1000.  
reaction.add_metabolites({ala_e:-1, ala_c:1}) 
model.add_reaction(reaction)

reaction = Reaction('arg_etr')
reaction.name = 'Arginine extracellular transport'
reaction.subsystem = 'AA'
reaction.lower_bound = -1000  
reaction.upper_bound = 1000.  
reaction.add_metabolites({arg_e:-1, arg_c:1}) 
model.add_reaction(reaction)

reaction = Reaction('asn_etr')
reaction.name = 'Asparagine extracellular transport'
reaction.subsystem = 'AA'
reaction.lower_bound = -1000  
reaction.upper_bound = 1000.  
reaction.add_metabolites({asn_e:-1, asn_c:1}) 
model.add_reaction(reaction)

reaction = Reaction('asp_etr')
reaction.name = 'Aspartate extracellular transport'
reaction.subsystem = 'AA'
reaction.lower_bound = -1000  
reaction.upper_bound = 1000.  
reaction.add_metabolites({asp_e:-1, asp_c:1}) 
model.add_reaction(reaction)

reaction = Reaction('cys_etr')
reaction.name = 'Aspartate extracellular transport'
reaction.subsystem = 'AA'
reaction.lower_bound = -1000  
reaction.upper_bound = 1000.  
reaction.add_metabolites({cys_e:-1, cys_c:1}) 
model.add_reaction(reaction)

#Warning!!!!
reaction = Reaction('glc_etr')
reaction.name = 'Glucose extracellular transport'
reaction.subsystem = 'Glycolysis'
reaction.lower_bound = -1000  
reaction.upper_bound = 1000.  
reaction.add_metabolites({glc_e:-1, glc_c:1}) 
model.add_reaction(reaction)

reaction = Reaction('gln_etr')
reaction.name = 'Glutamine extracellular transport'
reaction.subsystem = 'AA'
reaction.lower_bound = -1000  
reaction.upper_bound = 1000.  
reaction.add_metabolites({gln_e:-1, gln_c:1}) 
model.add_reaction(reaction)

reaction = Reaction('ile_mtr')
reaction.name = 'Isoleucine extracellular transport'
reaction.subsystem = 'AA'
reaction.lower_bound = -1000  
reaction.upper_bound = 1000.  
reaction.add_metabolites({ile_c:-1, ile_m:1}) 
model.add_reaction(reaction)

reaction = Reaction('leu_mtr')
reaction.name = 'Leucine extracellular transport'
reaction.subsystem = 'AA'
reaction.lower_bound = -1000  
reaction.upper_bound = 1000.  
reaction.add_metabolites({leu_c:-1, leu_m:1}) 
model.add_reaction(reaction)

reaction = Reaction('glu_etr')
reaction.name = 'Glutamate extracellular transport'
reaction.subsystem = 'AA'
reaction.lower_bound = -1000  
reaction.upper_bound = 1000.  
reaction.add_metabolites({glu_e:-1, glu_c:1}) 
model.add_reaction(reaction)

reaction = Reaction('gly_etr')
reaction.name = 'Glycine extracellular transport'
reaction.subsystem = 'AA'
reaction.lower_bound = -1000  
reaction.upper_bound = 1000.  
reaction.add_metabolites({gly_e:-1, gly_c:1}) 
model.add_reaction(reaction)

reaction = Reaction('his_etr')
reaction.name = 'Histidine extracellular transport'
reaction.subsystem = 'AA'
reaction.lower_bound = -1000  
reaction.upper_bound = 1000.  
reaction.add_metabolites({his_e:-1, his_c:1}) 
model.add_reaction(reaction)

reaction = Reaction('ile_etr')
reaction.name = 'Isoleucine extracellular transport'
reaction.subsystem = 'AA'
reaction.lower_bound = -1000  
reaction.upper_bound = 1000.  
reaction.add_metabolites({ile_e:-1, ile_c:1}) 
model.add_reaction(reaction)

reaction = Reaction('lac_etr')
reaction.name = 'Lactate extracellular secretion'
reaction.subsystem = 'Glycolysis'
reaction.lower_bound = -1000  
reaction.upper_bound = 1000.  
reaction.add_metabolites({lac_e:1, lac_c:-1}) 
model.add_reaction(reaction)

reaction = Reaction('leu_etr')
reaction.name = 'Leucine extracellular transport'
reaction.subsystem = 'AA'
reaction.lower_bound = -1000  
reaction.upper_bound = 1000.  
reaction.add_metabolites({leu_e:-1, leu_c:1}) 
model.add_reaction(reaction)

reaction = Reaction('lys_etr')
reaction.name = 'Lysine extracellular transport'
reaction.subsystem = 'AA'
reaction.lower_bound = -1000  
reaction.upper_bound = 1000.  
reaction.add_metabolites({lys_e:-1, lys_c:1}) 
model.add_reaction(reaction)

reaction = Reaction('met_etr')
reaction.name = 'Methionine extracellular transport'
reaction.subsystem = 'AA'
reaction.lower_bound = -1000  
reaction.upper_bound = 1000.  
reaction.add_metabolites({met_e:-1, met_c:1}) 
model.add_reaction(reaction)

reaction = Reaction('phe_etr')
reaction.name = 'Phenylalanine extracellular transport'
reaction.subsystem = 'AA'
reaction.lower_bound = -1000  
reaction.upper_bound = 1000.  
reaction.add_metabolites({phe_e:-1, phe_c:1}) 
model.add_reaction(reaction)

reaction = Reaction('pro_etr')
reaction.name = 'Proline extracellular transport'
reaction.subsystem = 'AA'
reaction.lower_bound = -1000  
reaction.upper_bound = 1000.  
reaction.add_metabolites({pro_e:-1, pro_c:1}) 
model.add_reaction(reaction)

reaction = Reaction('ser_etr')
reaction.name = 'Serine extracellular transport'
reaction.subsystem = 'AA'
reaction.lower_bound = -1000  
reaction.upper_bound = 1000.  
reaction.add_metabolites({ser_e:-1, ser_c:1}) 
model.add_reaction(reaction)

reaction = Reaction('thr_etr')
reaction.name = 'Threonine extracellular transport'
reaction.subsystem = 'AA'
reaction.lower_bound = -1000  
reaction.upper_bound = 1000.  
reaction.add_metabolites({thr_e:-1, thr_c:1}) 
model.add_reaction(reaction)

reaction = Reaction('trp_etr')
reaction.name = 'Tryptophan extracellular transport'
reaction.subsystem = 'AA'
reaction.lower_bound = -1000  
reaction.upper_bound = 1000.  
reaction.add_metabolites({trp_e:-1, trp_c:1}) 
model.add_reaction(reaction)

reaction = Reaction('tyr_etr')
reaction.name = 'Tyrosine extracellular transport'
reaction.subsystem = 'AA'
reaction.lower_bound = -1000  
reaction.upper_bound = 1000.  
reaction.add_metabolites({tyr_e:-1, tyr_c:1}) 
model.add_reaction(reaction)

reaction = Reaction('val_etr')
reaction.name = 'Valine extracellular transport'
reaction.subsystem = 'AA'
reaction.lower_bound = -1000  
reaction.upper_bound = 1000.  
reaction.add_metabolites({val_e:-1, val_c:1}) 
model.add_reaction(reaction)



"""reaction = Reaction('adn_etr')
reaction.name = 'Adenosine extracellular transport'
reaction.subsystem = 'misc'
reaction.lower_bound = -1000  
reaction.upper_bound = 1000.  
reaction.add_metabolites({adn_e:-1, adn_c:1}) 
model.add_reaction(reaction)"""

"""reaction = Reaction('ade_etr')
reaction.name = 'Adenine extracellular transport'
reaction.subsystem = 'misc'
reaction.lower_bound = -1000  
reaction.upper_bound = 1000.  
reaction.add_metabolites({ade_e:-1, ade_c:1}) 
model.add_reaction(reaction)"""

reaction = Reaction('dttp_etr')
reaction.name = 'dTTP extracellular transport'
reaction.subsystem = 'misc'
reaction.lower_bound = -1000  
reaction.upper_bound = 1000.  
reaction.add_metabolites({dttp_e:-1, dttp_c:1}) 
model.add_reaction(reaction)

reaction = Reaction('hco3ec')
reaction.name = 'HCO3 equilibration reaction'
reaction.subsystem = 'misc'
reaction.lower_bound = -1000  
reaction.upper_bound = 1000.  
reaction.add_metabolites({co2_c:-1, h2o_c:-1,hco3_c:1,h_c:1}) 
model.add_reaction(reaction)

#Warning: according to Recon1 HCO3 equilibration reaction should be mitochondrial, according to Recon2 extracelular

reaction = Reaction('hco3_mtr')
reaction.name = 'dTTP extracellular transport'
reaction.subsystem = 'misc'
reaction.lower_bound = -1000  
reaction.upper_bound = 1000.  
reaction.add_metabolites({hco3_c:-1, hco3_m:1}) 
model.add_reaction(reaction)


reaction = Reaction('pyr_mtr')
reaction.name = 'Pyruvate mitochondrial transport'
reaction.subsystem = 'misc'
reaction.lower_bound = -1000  
reaction.upper_bound = 1000.  
reaction.add_metabolites({pyr_b_c:-1,h_c:-1, pyr_m:1,h_m:1}) 
model.add_reaction(reaction)

reaction = Reaction('pro_mtr')
reaction.name = 'Proline mitochondrial transport'
reaction.subsystem = 'misc'
reaction.lower_bound = -1000  
reaction.upper_bound = 1000.  
reaction.add_metabolites({pro_c:-1,pro_m:1}) 
model.add_reaction(reaction)


reaction = Reaction('gln_mtr')
reaction.name = 'L-glutamine transport via electroneutral transporter'
reaction.subsystem = 'misc'
reaction.lower_bound = -1000  
reaction.upper_bound = 1000.  
reaction.add_metabolites({gln_c:-1,gln_m:1}) 
model.add_reaction(reaction)

reaction = Reaction('glu_mtr')
reaction.name = 'L-glutamate reversible transport via proton symport, mitochondrial'
reaction.subsystem = 'misc'
reaction.lower_bound = -1000  
reaction.upper_bound = 1000.  
reaction.add_metabolites({glu_c:-1,h_c:-1,glu_m:1,h_m:1}) 
model.add_reaction(reaction)  

reaction = Reaction('val_mtr')
reaction.name = 'Valine reversible mitochondrial transport'
reaction.subsystem = 'Transport, Mitochondrial'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({val_m: 1.0, val_c: -1.0})
model.add_reaction(reaction) 

reaction = Reaction('mal_mtr')
reaction.name = 'malate transport, mitochondrial'
reaction.subsystem = 'misc'
reaction.lower_bound = -1000  
reaction.upper_bound = 1000.  
reaction.add_metabolites({mal_c:-1,pi_m:-1,mal_m:1,pi_c:1}) 
model.add_reaction(reaction)    

reaction = Reaction('cit_mtr')
reaction.name = 'citrate transport, mitochondrial'
reaction.subsystem = 'misc'
reaction.lower_bound = -1000  
reaction.upper_bound = 1000.  
reaction.add_metabolites({cit_c:-1,mal_m:-1,cit_m:1,mal_c:1}) 
model.add_reaction(reaction)    

 

reaction = Reaction('urea_etr')
reaction.name = 'urea transport, extracellular'
reaction.subsystem = 'misc'
reaction.lower_bound = -1000  
reaction.upper_bound = 1000.  
reaction.add_metabolites({urea_e:-1,urea_c:1}) 
model.add_reaction(reaction)   


#Glycolysis

reaction = Reaction('hex1')
reaction.name = 'Hexokinase'
reaction.subsystem = 'Glycolysis/Gluconeogenesis'
reaction.lower_bound = 0.  
reaction.upper_bound = 1000.  
reaction.add_metabolites({glc_c: -1.0,
                          atp_c: -1.0,h_c: 1.0,
                          glc6p_a_c: 1.0,
                          adp_c: 1.0})  
model.add_reaction(reaction)

reaction = Reaction('glc6p_pdif')
reaction.name = 'Difussion between Glucose 6-phosphate pools'
reaction.subsystem = 'Glycolysis/Gluconeogenesis'
reaction.lower_bound = -1000.  
reaction.upper_bound = 1000.  
reaction.add_metabolites({glc6p_a_c: -1.0,
                          glc6p_b_c: 1.0})  
model.add_reaction(reaction)

reaction = Reaction('pgi')
reaction.name = 'Glucose-6-phosphate isomerase'
reaction.subsystem = 'Glycolysis/Gluconeogenesis'
reaction.lower_bound = -1000  
reaction.upper_bound = 1000.  
reaction.add_metabolites({glc6p_b_c: -1.0,
                           fru6p_c:1.0})  
model.add_reaction(reaction)



reaction = Reaction('pfk')
reaction.name = 'Phosphofructokinase'
reaction.subsystem = 'Glycolysis/Gluconeogenesis'
reaction.lower_bound = 0.  
reaction.upper_bound = 1000.  
reaction.add_metabolites({fru6p_c: -1.0,
                          atp_c: -1.0,
                          fru16bp_c: 1.0,
                          adp_c: 1.0,h_c: 1.0})        
model.add_reaction(reaction)

reaction = Reaction('fbp')
reaction.name = 'fructose-bisphosphatase'
reaction.subsystem = 'Glycolysis/Gluconeogenesis'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({h2o_c: -1.0, pi_c: 1.0, fru16bp_c: -1.0, fru6p_c: 1.0})
model.add_reaction(reaction)

reaction = Reaction('g6p')
reaction.name = 'Glucose-phosphatase'
reaction.subsystem = 'Glycolysis/Gluconeogenesis'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({h2o_c: -1.0, pi_c: 1.0, glc6p_a_c: -1.0, glc_c: 1.0})
model.add_reaction(reaction)

reaction = Reaction('fba')
reaction.name = 'Aldolase'
reaction.subsystem = 'Glycolysis/Gluconeogenesis'
reaction.lower_bound = -1000
reaction.upper_bound = 1000.  
reaction.add_metabolites({gap_c: 1.0,dhap_c: 1.0,fru16bp_c: -1.0})        
model.add_reaction(reaction)

reaction = Reaction('tpi')
reaction.name = 'Triose-phosphate isomerase'
reaction.subsystem = 'Glycolysis/Gluconeogenesis'
reaction.lower_bound = -1000
reaction.upper_bound = 1000.  
reaction.add_metabolites({gap_c: 1.0,dhap_c: -1.0})        
model.add_reaction(reaction)

reaction = Reaction('gapd')
reaction.name = 'Glyceraldehyde phosphate dehydrogenase'
reaction.subsystem = 'Glycolysis/Gluconeogenesis'
reaction.lower_bound = -1000
reaction.upper_bound = 1000.  
reaction.add_metabolites({gap_c: -1.0,nad_c: -1.0,pi_c: -1.0,bpg13_c:1.0,nadh_c:1.0,h_c: 1.0})        
model.add_reaction(reaction)


reaction = Reaction('pgk')
reaction.name = 'Phosphoglycerate Kinase'
reaction.subsystem = 'Glycolysis/Gluconeogenesis'
reaction.lower_bound = -1000
reaction.upper_bound = 1000.  
reaction.add_metabolites({bpg13_c:-1.0,adp_c:-1.0,atp_c:1.0,pg3_c:1.0})        
model.add_reaction(reaction)

reaction = Reaction('pgm')
reaction.name = 'Phosphoglycerate mutase'
reaction.subsystem = 'Glycolysis/Gluconeogenesis'
reaction.lower_bound = -1000
reaction.upper_bound = 1000.  
reaction.add_metabolites({pg3_c:-1.0,pg2_c:1.0})        
model.add_reaction(reaction)

reaction = Reaction('eno')
reaction.name = 'Enolase'
reaction.subsystem = 'Glycolysis/Gluconeogenesis'
reaction.lower_bound = -1000
reaction.upper_bound = 1000.  
reaction.add_metabolites({pg2_c:-1.0,pep_c:1.0, h2o_c:1})        
model.add_reaction(reaction)

reaction = Reaction('pyk')
reaction.name = 'Pyruvate Kinase'
reaction.subsystem = 'Glycolysis/Gluconeogenesis'
reaction.lower_bound = 0
reaction.upper_bound = 1000.  
reaction.add_metabolites({pep_c:-1.0,adp_c:-1.0,atp_c:1.0,pyr_a_c:1.0,h_c: -1.0})        
model.add_reaction(reaction)

reaction = Reaction('ldh_l')
reaction.name = 'L-lactate dehydrogenase'
reaction.subsystem = 'Glycolysis/Gluconeogenesis'
reaction.lower_bound = -1000
reaction.upper_bound = 1000.  
reaction.add_metabolites({pyr_a_c:-1.0,nadh_c:-1.0,lac_c:1.0,nad_c:1.0,h_c: -1.0})        
model.add_reaction(reaction)


#Glycogen 
reaction = Reaction('pgmt')
reaction.name = 'Phosphoglucomutase'
reaction.subsystem = 'Glycogen'
reaction.lower_bound =-1000.  
reaction.upper_bound = 1000.  
reaction.add_metabolites({glc6p_a_c:-1, glc1p_c:1}) 
model.add_reaction(reaction)

reaction = Reaction('galui')
reaction.name = 'UTP-glucose-1-phosphate uridylyltransferase (irreversible)'
reaction.subsystem = 'Starch and Sucrose Metabolism'
reaction.lower_bound =-1000.  
reaction.upper_bound = 1000.  
reaction.add_metabolites({glc1p_c:-1, utp_c:-1,h_c:-1,udpglc_c:1}) 
model.add_reaction(reaction)

reaction = Reaction('gs')
reaction.name = 'Glycogen synthase'
reaction.subsystem = 'Starch and Sucrose Metabolism'
reaction.lower_bound =0  
reaction.upper_bound = 1000  
reaction.add_metabolites({udpglc_c:-1,glygn_e:1,ppi_c:1,udp_c:1}) 
model.add_reaction(reaction)

reaction = Reaction('gp')
reaction.name = 'Glycogen phosphorylase'
reaction.subsystem = 'Starch and Sucrose Metabolism'
reaction.lower_bound =0  
reaction.upper_bound = 1000  
reaction.add_metabolites({glc1p_c:1,glygn_e:-1,ppi_c:-1}) 
model.add_reaction(reaction)

reaction = Reaction('pyr_pdif')
reaction.name = 'Difussion between pyruvate pools'
reaction.subsystem = 'Glycolysis'
reaction.lower_bound = -1000  
reaction.upper_bound = 1000.  
reaction.add_metabolites({pyr_a_c: -1.0,
                          pyr_b_c: 1.0})  
model.add_reaction(reaction)

#Pentose phosphate pathway

reaction = Reaction('g6pdh2r')
reaction.name = 'Glucose 6-Phosphate dehydrogenase'
reaction.subsystem = 'Pentose Phosphate Pathway'
reaction.lower_bound = 0  
reaction.upper_bound = 1000.  
reaction.add_metabolites({glc6p_b_c: -1.0,nadp_c:-1.0,
                           pglclac_c:1.0,nadph_c: 1.0,h_c: 1.0})  
model.add_reaction(reaction)

reaction = Reaction('pgl')
reaction.name = '6-phosphogluconolactonase'
reaction.subsystem = 'Pentose Phosphate Pathway'
reaction.lower_bound = 0  
reaction.upper_bound = 1000.  
reaction.add_metabolites({pglclac_c: -1.0,h2o_c:-1.0,glucon6p_c:1,
                            h_c: 1.0})  
model.add_reaction(reaction)

reaction = Reaction('gnd')
reaction.name = "6-phosphogluconate dehydrogenase"
reaction.subsystem = 'Pentose Phosphate Pathway'
reaction.lower_bound = 0  
reaction.upper_bound = 1000.  
reaction.add_metabolites({glucon6p_c: -1.0,nadp_c:-1.0,nadph_c:1.0,
                            co2_c:1.0,rul5p_c:1})  
model.add_reaction(reaction)

reaction = Reaction('rpi')
reaction.name = 'Ribose phosphate isomerase'
reaction.subsystem = 'Pentose Phosphate Pathway'
reaction.lower_bound = -1000
reaction.upper_bound = 1000.  
reaction.add_metabolites({rul5p_c:-1.0,rib5p_c:1.0})        
model.add_reaction(reaction)

reaction = Reaction('rpe')
reaction.name = 'Ribulose 5 phosphate epimerase'
reaction.subsystem = 'Pentose Phosphate Pathway'
reaction.lower_bound = -1000
reaction.upper_bound = 1000.  
reaction.add_metabolites({rul5p_c:-1.0,xyl5p_c:1.0})        
model.add_reaction(reaction)

reaction = Reaction('tkt1')
reaction.name = 'Transketolase 1'
reaction.subsystem = 'Pentose Phosphate Pathway'
reaction.lower_bound = -1000
reaction.upper_bound = 1000.  
reaction.add_metabolites({xyl5p_c:-1.0,rib5p_c:-1.0,gap_c:1,sed7p_c:1})        
model.add_reaction(reaction)



reaction = Reaction('tkt2')
reaction.name = 'Transketolase 2'
reaction.subsystem = 'Pentose Phosphate Pathway'
reaction.lower_bound = -1000
reaction.upper_bound = 1000.  
reaction.add_metabolites({gap_c:1.0,fru6p_c:1.0,e4p_c:-1,xyl5p_c:-1})        
model.add_reaction(reaction)

"""reaction = Reaction('tk3')
reaction.name = 'Transketolase 3'
reaction.subsystem = 'Test'
reaction.lower_bound = 0
reaction.upper_bound = 0 
reaction.add_metabolites({rib5p_c:-1.0,fru6p_c:-1.0,e4p_c:1,sed7p_c:1})        
model.add_reaction(reaction)"""

reaction = Reaction('tala')
reaction.name = 'Transaldolase'
reaction.subsystem = 'Pentose Phosphate Pathway'
reaction.lower_bound = -1000
reaction.upper_bound = 1000.  
reaction.add_metabolites({fru6p_c:1.0,e4p_c:1,sed7p_c:-1,gap_c:-1.0})        
model.add_reaction(reaction)

"""reaction = Reaction('nadphase')
reaction.name = 'NADPH hydrolase'
reaction.subsystem = 'Test'
reaction.lower_bound = 0
reaction.upper_bound = 1000.  
reaction.add_metabolites({nadph_c:-1.0,nadp_c:1.0})        
model.add_reaction(reaction)"""



#Energy metabolism
reaction = Reaction('ndpk1')
reaction.name = 'nucleoside-diphosphate kinase (ATP:GDP)'
reaction.subsystem = 'Energy'
reaction.lower_bound =-1000.  
reaction.upper_bound = 1000.  
reaction.add_metabolites({atp_c:-1, adp_c:1,gdp_c:-1, gtp_c:1}) 
model.add_reaction(reaction)

reaction = Reaction('ndpk1m')
reaction.name = 'nucleoside-diphosphate kinase (ATP:GDP)'
reaction.subsystem = 'Energy'
reaction.lower_bound =-1000.  
reaction.upper_bound = 1000.  
reaction.add_metabolites({atp_m:-1, adp_m:1,gdp_m:-1, gtp_m:1}) 
model.add_reaction(reaction)

reaction = Reaction('ndpk2')
reaction.name = 'nucleoside-diphosphate kinase (ATP:UDP)'
reaction.subsystem = 'Energy'
reaction.lower_bound =-1000.  
reaction.upper_bound = 1000.  
reaction.add_metabolites({atp_c:-1, adp_c:1,udp_c:-1, utp_c:1}) 
model.add_reaction(reaction)



reaction = Reaction('atptm')
reaction.name = 'ADP/ATP transporter, mitochondrial'
reaction.subsystem = 'Transport, Mitochondrial'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({adp_m: 1.0, atp_m: -1.0, atp_c: 1.0, adp_c: -1.0})
model.add_reaction(reaction)

reaction = Reaction('pi_etr')
reaction.name = 'Phosphate extracelular carrier'
reaction.subsystem = 'Energy'
reaction.lower_bound =-1000.  
reaction.upper_bound = 1000.  
reaction.add_metabolites({pi_c:1, pi_e:-1}) 
model.add_reaction(reaction)

reaction = Reaction('pi_mtr')
reaction.name = 'phosphate transporter, mitochondrial'
reaction.subsystem = 'Energy'
reaction.lower_bound =-1000.  
reaction.upper_bound = 1000.  
reaction.add_metabolites({pi_c:-1,h_c:-1, pi_m:1,h_m:1}) 
model.add_reaction(reaction)

reaction = Reaction('ppa')
reaction.name = 'Inorganic diphosphatase"'
reaction.subsystem = 'Energy'
reaction.lower_bound =0  
reaction.upper_bound = 1000.  
reaction.add_metabolites({ppi_c:-1, pi_c:2}) 
model.add_reaction(reaction)

reaction = Reaction('atpase')
reaction.name = 'ATP hydrolase'
reaction.subsystem = 'Energy'
reaction.lower_bound = 0
reaction.upper_bound = 1000.  
reaction.add_metabolites({atp_c:-1.0,adp_c:1.0,pi_c:1.0,h_c:1})        
model.add_reaction(reaction)

reaction = Reaction('adk1')
reaction.name = 'Adenylate kinase'
reaction.subsystem = 'Energy'
reaction.lower_bound = -1000
reaction.upper_bound = 1000.  
reaction.add_metabolites({amp_c:-1.0,atp_c:-1.0,adp_c:2.0})        
model.add_reaction(reaction)


##Protein metabolism: 


#Serine metabolism
#c3PG + Glu + cNAD ↔ Ser + αKG + cNADH
reaction = Reaction('pgcd')
reaction.name = 'Phosphoglycerate dehydrogenase'
reaction.subsystem = 'Glycine, Serine, and Threonine Metabolism'
reaction.lower_bound = 0
reaction.upper_bound = 1000.  
reaction.add_metabolites({pg3_c:-1.0,nad_c:-1.0,p3hp_c:1,nadh_c:1.0,h_c:1.0})        
model.add_reaction(reaction)

reaction = Reaction('psert')
reaction.name = 'Phosphoserine transaminase'
reaction.subsystem = 'Glycine, Serine, and Threonine Metabolism'
reaction.lower_bound = 0
reaction.upper_bound = 1000.  
reaction.add_metabolites({glu_c:-1.0,p3hp_c:-1.0,pser_c:1.0,akg_c:1.0})        
model.add_reaction(reaction)

reaction = Reaction('psp_l')
reaction.name = 'Phosphoserine phosphatase'
reaction.subsystem = 'Glycine, Serine, and Threonine Metabolism'
reaction.lower_bound = 0
reaction.upper_bound = 1000.  
reaction.add_metabolites({ser_c:1.0,pi_c:1.0,pser_c:-1.0,h2o_c:-1.0})        
model.add_reaction(reaction)

#Ser + Pyr2 + cNADPH + ATP ↔ c3PG + Ala + cNADP + ADP
#Serine degradation 1 (active in tumor cells) 
#TO DO, change toe perixosomal
reaction = Reaction('sptix')
reaction.name = 'serine-pyruvate-transaminase'
reaction.subsystem = 'Serine metabolsim'
reaction.lower_bound = 0
reaction.upper_bound = 1000.  
reaction.add_metabolites({ser_c:-1,pyr_b_c:-1,ala_c:1,hpyr_c:1})        
model.add_reaction(reaction)

reaction = Reaction('hpyrr')
reaction.name = 'Hydroxypyruvate reductase (NADPH)'
reaction.subsystem = 'Serine metabolsim'
reaction.lower_bound = 0
reaction.upper_bound = 1000.  
reaction.add_metabolites({hpyr_c:-1,nadph_c:-1,h_c:-1,glyc_c:1,nadp_c:1})        
model.add_reaction(reaction)

reaction = Reaction('glyck2')
reaction.name = 'Glycerate kinase'
reaction.subsystem = 'Serine metabolsim'
reaction.lower_bound = 0
reaction.upper_bound = 1000.  
reaction.add_metabolites({glyc_c:-1.0,atp_c:-1.0,pg2_c:1.0,h_c:1.0,adp_c:1.0})        
model.add_reaction(reaction)

#Serine degradation 2 (active in somatic cells) 
reaction = Reaction('serhl')
reaction.name = 'L-Serine hydro-lyase'
reaction.subsystem = 'Glycine, Serine, and Threonine Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({n2amac_c: 1.0, ser_c: -1.0, h2o_c: 1.0})
model.add_reaction(reaction)

reaction = Reaction('2amachyd')
reaction.name = ' 2-Aminoacrylate hydrolysis'
reaction.subsystem = 'Glycine, Serine, and Threonine Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({n2amac_c: -1.0, h2o_c: -1.0, pyr_b_c: 1.0, nh4_c: 1.0})
model.add_reaction(reaction)



reaction = Reaction('alata_l')
reaction.name = 'Alanine transaminase'
reaction.subsystem = 'alanine metabolsim'
reaction.lower_bound = -1000
reaction.upper_bound = 1000  
reaction.add_metabolites({ala_c:-1.0,akg_c:-1.0,pyr_b_c:1.0,glu_c:1.0})        
model.add_reaction(reaction)
"""Nota: En uniprot hay referencia a que la última reacción (glycerate kinase) tiene como reactivo 3PG y no 2PG, lo mismo que metacyc y otras bases de datos. Para el libro “Metabolism at a glance” es 2PG y toda la via es mitocondrial. Para el RECON1 es toda la via citoplasmática y puede ser 2PG o 3PG"""

#Glycine, Serine, and Threonine Metabolism
reaction = Reaction('ghmt2r')
reaction.name = 'glycine hydroxymethyltransferase, reversible'
reaction.subsystem = 'Glycine, Serine, and Threonine Metabolism'
reaction.lower_bound = -1000
reaction.upper_bound = 1000  
reaction.add_metabolites({thf_c:-1.0,ser_c:-1.0,mlthf_c:1.0,gly_c:1.0})        
model.add_reaction(reaction)

#c510mTHF-> DHF
reaction = Reaction('tmds')
reaction.name = 'Thymidylate synthase '
reaction.subsystem = 'Glycine, Serine, and Threonine Metabolism'
reaction.lower_bound = 0 #Assumed irreversibles, recon 1 considers it reveresible
reaction.upper_bound = 1000  
reaction.add_metabolites({mlthf_c:-1.0,dump_c:-1.0,dtmp_c:1.0,dhf_c:1.0})        
model.add_reaction(reaction)



reaction = Reaction('dtmpk')
reaction.name = 'dTMP kinase'
reaction.subsystem = 'Glycine, Serine, and Threonine Metabolism'
reaction.lower_bound = -1000
reaction.upper_bound = 1000  
reaction.add_metabolites({atp_c:-1.0,dtmp_c:-1.0,dtdp_c:1.0,adp_c:1.0})        
model.add_reaction(reaction)

"""reaction = Reaction('dtdpk')
reaction.name = 'dTDP kinase'
reaction.subsystem = 'Glycine, Serine, and Threonine Metabolism'
reaction.lower_bound = -1000
reaction.upper_bound = 1000  
reaction.add_metabolites({atp_c:-1.0,dtdp_c:-1.0,dttp_c:1.0,adp_c:1.0})        
model.add_reaction(reaction)"""
#Warning, add exit of dttp to biomass

#To do: cytoplasmatic fumarasa, add exit of dttp to biomass

reaction = Reaction('dhfr')
reaction.name = 'dihydrofolate reductase (NADPH)'
reaction.subsystem = 'Glycine, Serine, and Threonine Metabolism'
reaction.lower_bound = 0
reaction.upper_bound = 1000.  
reaction.add_metabolites({dhf_c:-1,nadph_c:-1,h_c:-1,thf_c:1,nadp_c:1})        
model.add_reaction(reaction)


reaction = Reaction('mthfr3')
reaction.name = '5,10-methylenetetrahydrofolatereductase (NADPH)'
reaction.subsystem = 'Glycine, Serine, and Threonine Metabolism'
reaction.lower_bound = 0
reaction.upper_bound = 1000.  
reaction.add_metabolites({mlthf_c:-1,nadph_c:-1,h_c:-1,m5thf_c:1,nadp_c:1})        
model.add_reaction(reaction)
#cytoplasmatic fumarasa
#Alerta el pedro no considera el NADPH Dihydrofolate reductase, i Methylenetetrahydrofolate reductase

reaction = Reaction('mets')
reaction.name = 'methionine synthase'
reaction.subsystem = 'Glycine, Serine, and Threonine Metabolism'
reaction.lower_bound = 0
reaction.upper_bound = 1000.  
reaction.add_metabolites({m5thf_c:-1,hcys_c:-1,thf_c:1,met_c:1})        
model.add_reaction(reaction)



"""
#TODO IMPROVE
Nota simplificacion
R29
Adomet ↔ HCys
A variety of methyl transferases
Adomet + demethylated acceptor => AdoHCys + methylated acceptor + H+
Adenosylhomocysteinase, encoded by AHCYL1, AHCY and AHCYL2, NAD+ as cofactor
S-adenosyl-L-homocysteine (AdoHcys) + H2O ↔ L-homocysteine (HCys) + adenosine"""

reaction = Reaction('ametmtr')
reaction.name = 'variety of methyl transferases'
reaction.subsystem = ''
reaction.lower_bound = 0
reaction.upper_bound = 1000.  
reaction.add_metabolites({amet_c:-1,methyl_a_c:-1,ahcys_c:1,ch3_methyl_a_c:1,h_c:1})        
model.add_reaction(reaction)

reaction = Reaction('methyl_a_demet')
reaction.name = 'methylated acceptor demethylator'
reaction.subsystem = 'artificial'
reaction.lower_bound = -1000
reaction.upper_bound = 1000.  
reaction.add_metabolites({methyl_a_c:1,ch3_methyl_a_c:-1})        
model.add_reaction(reaction)

reaction = Reaction('ahc')
reaction.name = 'Adenosylhomocysteinase'
reaction.subsystem = ''
reaction.lower_bound = 0
reaction.upper_bound = 1000.  
reaction.add_metabolites({ahcys_c:-1,h2o_c:-1,hcys_c:1,adn_c:1})        
model.add_reaction(reaction)

#TCA 
reaction = Reaction('pdhm')
reaction.name = 'pyruvate dehydrogenase'
reaction.subsystem = 'Glycolysis/Gluconeogenesis'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({pyr_m: -1.0, accoa_m: 1.0, nadh_m: 1.0, nad_m: -1.0, coa_m: -1.0, co2_m: 1.0})
model.add_reaction(reaction)

reaction = Reaction('pcm')
reaction.name = 'pyruvate carboxylase'
reaction.subsystem = 'Pyruvate Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({pyr_m: -1.0, hco3_m: -1.0, atp_m: -1.0, oaa_m: 1.0, pi_m: 1.0, adp_m: 1.0, h_m: 1.0})
model.add_reaction(reaction)

reaction = Reaction('pepck')
reaction.name = 'Phosphoenolpyruvate carboxykinase (GTP)'
reaction.subsystem = 'Glycolysis/Gluconeogenesis'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({pep_c: 1.0, gdp_c: 1.0, co2_c: 1.0, gtp_c: -1.0, oaa_c: -1.0})
model.add_reaction(reaction)


reaction = Reaction('csm')
reaction.name = 'citrate synthase'
reaction.subsystem = 'Citric Acid Cycle'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({cit_m: 1.0, h2o_m: -1.0, accoa_m: -1.0, oaa_m: -1.0, coa_m: 1.0, h_m: 1.0})
model.add_reaction(reaction)

reaction = Reaction('acontm')
reaction.name = 'Aconitate hydratase, mitochondrial'
reaction.subsystem = 'Citric Acid Cycle'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({cit_m: -1.0, icit_m: 1.0})
model.add_reaction(reaction)

reaction = Reaction('acont')
reaction.name = 'Aconitate hydratase, cytosolic'
reaction.subsystem = 'Citric Acid Cycle'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({cit_c: -1.0, icit_c: 1.0})
model.add_reaction(reaction)

reaction = Reaction('icdhxm')
reaction.name = 'Isocitrate dehydrogenase (NAD+),mitochondrial'
reaction.subsystem = 'Citric Acid Cycle'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({nad_m: -1.0, nadh_m: 1.0, icit_m: -1.0, akg_m: 1.0, co2_m: 1.0})
model.add_reaction(reaction)

reaction = Reaction('icdhyrm')
reaction.name = 'Isocitrate dehydrogenase (NADP+), mitochondrial'
reaction.subsystem = 'Citric Acid Cycle'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({co2_m: 1.0, nadph_m: 1.0, nadp_m: -1.0, akg_m: 1.0, icit_m: -1.0})
model.add_reaction(reaction)


reaction = Reaction('icdhy')
reaction.name = 'isocitrate dehydrogenase (NADP), cytosolic'
reaction.subsystem = 'Citric Acid Cycle'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({akg_c: 1.0, icit_c: -1.0, nadp_c: -1.0, co2_c: 1.0, nadph_c: 1.0})
model.add_reaction(reaction)


reaction = Reaction('akgdm')
reaction.name = ' 2-oxoglutarate dehydrogenase'
reaction.subsystem = 'Citric Acid Cycle'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({succoa_m: 1.0, akg_m: -1.0, nadh_m: 1.0, nad_m: -1.0, coa_m: -1.0, co2_m: 1.0})
model.add_reaction(reaction)

reaction = Reaction('sucoas1m')
reaction.name = 'Succinate--CoA ligase (GDP-forming)'
reaction.subsystem = 'Citric Acid Cycle'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({succoa_m: 1.0, gdp_m: 1.0, pi_m: 1.0, coa_m: -1.0, succ_m: -1.0, gtp_m: -1.0})
model.add_reaction(reaction)

reaction = Reaction('sucd1m')
reaction.name = 'succinate dehydrogenase'
reaction.subsystem = 'Citric Acid Cycle'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({fad_m: -1.0, fadh2_m: 1.0, succ_m: -1.0, fum_m: 1.0})
model.add_reaction(reaction)

reaction = Reaction('nadh2_u10m')
reaction.name = 'NADH dehydrogenase, mitochondrial'
reaction.subsystem = 'Oxidative Phosphorylation'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({q10h2_m: 1.0, nadh_m: -1.0, h_c: 4.0, nad_m: 1.0, h_m: -5.0, q10_m: -1.0})
model.add_reaction(reaction)

reaction = Reaction('cyor_u10m')
reaction.name = 'ubiquinol-6 cytochrome c reductase, Complex III'
reaction.subsystem = 'Oxidative Phosphorylation'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({h_m: -2.0, ficytc_m: -2.0, h_c: 4.0, q10h2_m: -1.0, focytc_m: 2.0, q10_m: 1.0})
model.add_reaction(reaction)

"""
reaction including small ross production
reaction = Reaction('cy00m3')
reaction.name = 'cytochrome c oxidase, mitochondrial Complex IV'
reaction.subsystem = 'Oxidative Phosphorylation'
reaction.lower_bound = 0
reaction.upper_bound = 1000.  
reaction.add_metabolites({ficytc_m:4,focytc_m:-4,h_c:4,h_m:-7.92,o2_m:-1,h2o_m:1.96,o2s_m:1})        
model.add_reaction(reaction)"""

reaction = Reaction('cyoom3')
reaction.name = 'cytochrome c oxidase, mitochondrial Complex IV'
reaction.subsystem = 'Oxidative Phosphorylation'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({h2o_m: 2, o2_m: -1.0, ficytc_m: 4.0, h_c: 4.0,  h_m: -8, focytc_m: -4.0})
model.add_reaction(reaction)

reaction = Reaction('etf')
reaction.name = 'electron transfer flavoprotein'
reaction.subsystem = 'Fatty acid oxidation'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({fad_m: 1.0, fadh2_m: -1.0, etfox_m: -1.0, etfrd_m: 1.0})
model.add_reaction(reaction)

reaction = Reaction('etfqo')
reaction.name = 'Electron transfer flavoprotein-ubiquinone oxidoreductase'
reaction.subsystem = 'Fatty acid oxidation'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({q10h2_m: 1.0, etfox_m: 1.0, etfrd_m: -1.0, q10_m: -1.0})
model.add_reaction(reaction)

reaction = Reaction('atps4m')
reaction.name = 'ATP synthase (four protons for one ATP)'
reaction.subsystem = 'Oxidative Phosphorylation'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({atp_m: 1.0, h2o_m: 1.0, h_c: -4.0, adp_m: -1.0, h_m: 3.0, pi_m: -1.0})
model.add_reaction(reaction)

reaction = Reaction('fumm')
reaction.name = 'fumarase, mitochondrial'
reaction.subsystem = 'Citric Acid Cycle'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({mal_m: 1.0, h2o_m: -1.0, fum_m: -1.0})
model.add_reaction(reaction)

reaction = Reaction('mdhm')
reaction.name = 'malate dehydrogenase, mitochondrial'
reaction.subsystem = 'Citric Acid Cycle'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({nad_m: -1.0, oaa_m: 1.0, h_m: 1.0, mal_m: -1.0, nadh_m: 1.0})
model.add_reaction(reaction)

reaction = Reaction('mdhc')
reaction.name = "malate dehydrogenase, cytosolic"
reaction.subsystem = 'Oxidative Phosphorylation'
reaction.lower_bound = -1000
reaction.upper_bound = 1000.  
reaction.add_metabolites({mal_c:-1,nad_c:-1,oaa_c:1,nadh_c:1,h_c:1})        
model.add_reaction(reaction)

reaction = Reaction('glnm')
reaction.name = "glutaminase (mitochondrial)"
reaction.subsystem = 'Glutamate metabolism'
reaction.lower_bound = 0
reaction.upper_bound = 1000.  
reaction.add_metabolites({gln_m:-1,h2o_m:-1,nh4_m:1,glu_m:1})        
model.add_reaction(reaction)

reaction = Reaction('gludxm')
reaction.name = "glutamate dehydrogenase (NAD),mitochondrial"
reaction.subsystem = 'Glutamate metabolism'
reaction.lower_bound = -1000
reaction.upper_bound = 1000.  
reaction.add_metabolites({glu_m:-1,nad_m:-1,h2o_m:-1,akg_m:1,nadh_m:1,nh4_m:1,h_m:1})        
model.add_reaction(reaction)

reaction = Reaction('gludym')
reaction.name = "glutamate dehydrogenase (NADP), mitochondrial"
reaction.subsystem = 'Glutamate metabolism'
reaction.lower_bound = -1000
reaction.upper_bound = 1000.  
reaction.add_metabolites({glu_m:-1,nadp_m:-1,h2o_m:-1,akg_m:1,nadph_m:1,nh4_m:1,h_m:1})        
model.add_reaction(reaction)

reaction = Reaction('glns')
reaction.name = 'glutamine synthetase'
reaction.subsystem = 'Glutamate metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({atp_c: -1.0, gln_c: 1.0, h_c: 1.0, adp_c: 1.0, glu_c: -1.0, pi_c: 1.0, nh4_c: -1.0})
model.add_reaction(reaction)

reaction = Reaction('asptam')
reaction.name = "aspartate transaminase, mitochondrial"
reaction.subsystem = 'Malate aspratate shuttle'
reaction.lower_bound = -1000.
reaction.upper_bound = 1000.  
reaction.add_metabolites({akg_m:-1,asp_m:-1,glu_m:1,oaa_m:1})        
model.add_reaction(reaction)

reaction = Reaction('asptac')
reaction.name = "aspartate transaminase, cytosolic"
reaction.subsystem = 'Malate aspratate shuttle'
reaction.lower_bound = -1000.
reaction.upper_bound = 1000.  
reaction.add_metabolites({akg_c:-1,asp_c:-1,glu_c:1,oaa_c:1})        
model.add_reaction(reaction)

reaction = Reaction('akgmaltm')
reaction.name = "alpha-ketoglutarate/malate transporter"
reaction.subsystem = 'Malate aspratate shuttle'
reaction.lower_bound = -1000.
reaction.upper_bound = 1000.  
reaction.add_metabolites({akg_m:-1,mal_c:-1,akg_c:1,mal_m:1})        
model.add_reaction(reaction)

reaction = Reaction('aspglum')
reaction.name = "aspartate-glutamate mitochondrial shuttle"
reaction.subsystem = 'Malate aspratate shuttle'
reaction.lower_bound = -1000.
reaction.upper_bound = 1000.  
reaction.add_metabolites({asp_m:-1,glu_c:-1,asp_c:1,glu_m:1})        
model.add_reaction(reaction)

reaction = Reaction('asns1')
reaction.name = "asparagine synthase (glutamine-hydrolysing)"
reaction.subsystem = "Alanine and Aspartate Metabolism"
reaction.lower_bound = 0
reaction.upper_bound = 1000.  
reaction.add_metabolites({asp_c:-1,atp_c:-1,gln_c:-1,h2o_c:-1,amp_c:1,asn_c:1,glu_c:1,h_c:1,ppi_c:1})        
model.add_reaction(reaction)

reaction = Reaction('acitl')
reaction.name = "ATP-Citrate lyase"
reaction.subsystem = "Faty acid synthesis"
reaction.lower_bound = 0
reaction.upper_bound = 1000.  
reaction.add_metabolites({atp_c:-1,cit_c:-1,coa_c:-1,accoa_c:1,adp_c:1,oaa_c:1,pi_c:1})        
model.add_reaction(reaction)

reaction = Reaction('accoac')
reaction.name = "acetyl-CoA carboxylase"
reaction.subsystem = "Faty acid synthesis"
reaction.lower_bound = 0
reaction.upper_bound = 1000.  
reaction.add_metabolites({accoa_c:-1,atp_c:-1,hco3_c:-1,adp_c:1,h_c:1,malcoa_c:1,pi_c:1})        
model.add_reaction(reaction)

reaction = Reaction('kas8')
reaction.name = 'b-ketoacyl synthetase (palmitate, n-C16:0)'
reaction.subsystem = 'Fatty acid elongation'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({h2o_c: 6.0, malcoa_c: -7.0, nadph_c: -14.0, nadp_c: 14.0, coa_c: 8.0, h_c: -20.0, accoa_c: -1.0, co2_c: 7.0, hdca_c: 1.0})
model.add_reaction(reaction)

reaction = Reaction('me1m')
reaction.name = "malic enzyme (NAD), mitochondrial"
reaction.subsystem = "Pyruvate Metabolism"
reaction.lower_bound = 0
reaction.upper_bound = 1000.  
reaction.add_metabolites({mal_m:-1,nad_m:-1,co2_m:1,nadh_m:1,pyr_m:1})  
model.add_reaction(reaction)

reaction = Reaction('me2')
reaction.name = "malic enzyme (NADP), cytosolic"
reaction.subsystem = "Pyruvate Metabolism"
reaction.lower_bound = 0
reaction.upper_bound = 1000.  
reaction.add_metabolites({mal_c:-1,nadp_c:-1,co2_c:1,nadph_c:1,pyr_b_c:1})  
model.add_reaction(reaction)

reaction = Reaction('glu5km')
reaction.name = "glutamate 5-kinase (m)"
reaction.subsystem = "Urea cycle/amino group metabolism"
reaction.lower_bound = 0
reaction.upper_bound = 1000.  
reaction.add_metabolites({atp_m:-1,glu_m:-1,adp_m:1,glu5p_m:1})  
model.add_reaction(reaction)

reaction = Reaction('g5sdym')
reaction.name = "glutamate-5-semialdehyde dehydrogenase (m)"
reaction.subsystem = "Urea cycle/amino group metabolism"
reaction.lower_bound = 0
reaction.upper_bound = 1000.  
reaction.add_metabolites({glu5p_m:-1,h_m:-1,nadph_m:-1,glu5sa_m:1,nadp_m:1,pi_m:1})  
model.add_reaction(reaction)


reaction = Reaction('g5sadrm')
reaction.name = "L-glutamate 5-semialdehyde dehydratase, reversible, mitochondrial"
reaction.subsystem = "Urea cycle/amino group metabolism"
reaction.lower_bound = -1000
reaction.upper_bound = 1000.  
reaction.add_metabolites({glu5sa_m:-1,n1pyr5c_m:1,h2o_m:1,h_m:1})  
model.add_reaction(reaction)


reaction = Reaction('p5crxm')
reaction.name = "pyrroline-5-carboxylate reductase (m)"
reaction.subsystem = "Urea cycle/amino group metabolism"
reaction.lower_bound = 0
reaction.upper_bound = 1000.  
reaction.add_metabolites({n1pyr5c_m:-1,h_m:-1,nadh_m:-1,nad_m:1,pro_m:1})  
model.add_reaction(reaction)




reaction = Reaction('argn')
reaction.name = "arginase"
reaction.subsystem = "Urea cycle/amino group metabolism"
reaction.lower_bound = 0
reaction.upper_bound = 1000.  
reaction.add_metabolites({arg_c:-1,h2o_c:-1,orn_c:1,urea_c:1})  
model.add_reaction(reaction)

reaction = Reaction('ornt3m')
reaction.name = "ornithine mitochondrial transport via proton antiport)"
reaction.subsystem = "Urea cycle/amino group metabolism"
reaction.lower_bound = -1000
reaction.upper_bound = 1000.  
reaction.add_metabolites({h_c:-1,orn_m:-1,orn_c:1,h_m:1})  
model.add_reaction(reaction) #Check if it can work in anaeerobic conditions

reaction = Reaction('orntarm')
reaction.name = "ornithine transaminase reversible (m)"
reaction.subsystem = "Urea cycle/amino group metabolism"
reaction.lower_bound = -1000
reaction.upper_bound = 1000.  
reaction.add_metabolites({akg_m:-1,orn_m:-1,glu5sa_m:1,glu_m:1})  
model.add_reaction(reaction)

reaction = Reaction('p5cdm')
reaction.name = "1-pyrroline-5-carboxylate dehydrogenase, mitochondrial"
reaction.subsystem = "Urea cycle/amino group metabolism"
reaction.lower_bound = 0
reaction.upper_bound = 1000.  
reaction.add_metabolites({n1pyr5c_m:-1,h2o_m:-1,nad_m:-1,glu_m:1,h_m:1,nadh_m:1})  
model.add_reaction(reaction)

#ProlineCatabolism
reaction = Reaction('prod2m')
reaction.name = 'Proline dehydrogenase (m)'
reaction.subsystem = 'Arginine and Proline Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({fad_m: -1.0, h_m: 1.0, fadh2_m: 1.0, pro_m: -1.0, n1pyr5c_m: 1.0})
model.add_reaction(reaction)


reaction = Reaction('hisd')
reaction.name = "histidase"
reaction.subsystem = "Histidine Metabolism"
reaction.lower_bound = 0
reaction.upper_bound = 1000.  
reaction.add_metabolites({his_c:-1,nh4_c:1,urcan_c:1})  
model.add_reaction(reaction)

reaction = Reaction('urcn')
reaction.name = "urocanase"
reaction.subsystem = "Histidine Metabolism"
reaction.lower_bound = 0
reaction.upper_bound = 1000.  
reaction.add_metabolites({urcan_c:-1,h2o_c:-1,n4izp_c:1})  
model.add_reaction(reaction)

reaction = Reaction('izpn')
reaction.name = "imidazolonepropionase"
reaction.subsystem = "Histidine Metabolism"
reaction.lower_bound = 0
reaction.upper_bound = 1000.  
reaction.add_metabolites({n4izp_c:-1,h2o_c:-1,forglu_c:1,h_c:1})  
model.add_reaction(reaction)

reaction = Reaction('glufortx')
reaction.name = 'Glutamate formimidoyltransferase'
reaction.subsystem = 'Histidine Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({thf_c: -1.0, glu_c: 1.0, forglu_c: -1.0, h_c: -1.0, n5forthf_c: 1.0})
model.add_reaction(reaction)

reaction = Reaction('ftcd')
reaction.name = 'formimidoyltransferase cyclodeaminase'
reaction.subsystem = 'Folate Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({nh4_c: 1.0, methf_c: 1.0, h_c: -2.0, n5forthf_c: -1.0})
model.add_reaction(reaction)

reaction = Reaction('mthfc')
reaction.name = 'methenyltetrahydrofolate cyclohydrolase'
reaction.subsystem = 'Folate Metabolism'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({h2o_c: -1.0, methf_c: -1.0, h_c: 1.0, n10fthf_c: 1.0})
model.add_reaction(reaction)

reaction = Reaction('fthfl')
reaction.name = 'formate-tetrahydrofolate ligase'
reaction.subsystem = 'Folate Metabolism'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({adp_c: 1.0, atp_c: -1.0, thf_c: -1.0, for_c: -1.0, pi_c: 1.0, n10fthf_c: 1.0})
model.add_reaction(reaction)


reaction = Reaction('orndc')
reaction.name = "Ornithine Decarboxylase"
reaction.subsystem = "Urea cycle/amino group metabolism"
reaction.lower_bound = 0
reaction.upper_bound = 1000.  
reaction.add_metabolites({h_c:-1,orn_c:-1,co2_c:1,ptrc_c:1})  
model.add_reaction(reaction)

reaction = Reaction('metat')
reaction.name = 'methionine adenosyltransferase'
reaction.subsystem = 'Methionine Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({h2o_c: -1.0, atp_c: -1.0, met_c: -1.0, ppi_c: 1.0, amet_c: 1.0, pi_c: 1.0})
model.add_reaction(reaction)

reaction = Reaction('admdc')
reaction.name = 'adenosylmethionine decarboxylase'
reaction.subsystem = 'Arginine and Proline Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({amet_c: -1.0, co2_c: 1.0, h_c: -1.0, ametam_c: 1.0})
model.add_reaction(reaction) 




reaction = Reaction('mtap')
reaction.name = ' 5-methylthioadenosine phosphorylase'
reaction.subsystem = 'Arginine and Proline Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({n5mta_c: -1.0, ade_c: 1.0, pi_c: -1.0, n5mdr1p_c: 1.0})
model.add_reaction(reaction)

reaction = Reaction('adpt')
reaction.name = 'adenine phosphoribosyltransferase'
reaction.subsystem = 'Salvage Pathway'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({ppi_c: 1.0, ade_c: -1.0, amp_c: 1.0, prpp_c: -1.0})
model.add_reaction(reaction)

reaction = Reaction('mtri')
reaction.name = ' 5-methylthioribose-1-phosphate isomerase'
reaction.subsystem = 'Arginine and Proline Metabolism'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({n5mdr1p_c: -1.0, n5mdru1p_c: 1.0})
model.add_reaction(reaction)

reaction = Reaction('mdrpd')
reaction.name = ' 5-Methylthio-5-deoxy-D-ribulose 1-phosphate dehydratase'
reaction.subsystem = 'Arginine and Proline Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({dkmpp_c: 1.0, h2o_c: 1.0, n5mdru1p_c: -1.0})
model.add_reaction(reaction)

#Warning some of the reactions of the methionine salvage pathway differ from recon1 # http://humancyc.org/HUMAN/NEW-IMAGE?type=PATHWAY&object=PWY-7527&detail-level=2# 

reaction = Reaction('dkmppd')
reaction.name = ' 2,3-diketo-5-methylthio-1-phosphopentane degradation reaction'
reaction.subsystem = 'Arginine and Proline Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({dkmpp_c: -1.0, h2o_c: -1.0, o2_c: -1.0, n2kmb_c: 1.0, h_c: 2.0, for_c: 1.0, pi_c: 1.0})
model.add_reaction(reaction)

reaction = Reaction('unk3')
reaction.name = ' 2-keto-4-methylthiobutyrate transamination'
reaction.subsystem = 'Arginine and Proline Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({akg_c: 1.0, glu_c: -1.0, n2kmb_c: -1.0, met_c: 1.0})
model.add_reaction(reaction)

reaction = Reaction('fdh')
reaction.name = 'formate dehydrogenase'
reaction.subsystem = 'Fatty Acid Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({nadh_c: 1.0, co2_c: 1.0, for_c: -1.0, nad_c: -1.0})
model.add_reaction(reaction)

reaction = Reaction('spms')
reaction.name = 'spermidine synthase'
reaction.subsystem = 'Urea cycle/amino group metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({ametam_c: -1.0, spmd_c: 1.0, n5mta_c: 1.0, h_c: 1.0, ptrc_c: -1.0})
model.add_reaction(reaction)

reaction = Reaction('sprms')
reaction.name = 'spermine synthase'
reaction.subsystem = 'Urea cycle/amino group metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({ametam_c: -1.0, spmd_c: -1.0, n5mta_c: 1.0, sprm_c: 1.0, h_c: 1.0})
model.add_reaction(reaction)

"""Omited from the model for now: 
cAcCoA + Spd ↔ AcSpd + cCoA
AcCoA + Spn ↔ AcSpn + CoA

Spermidine/spermine N1-acetyltransferase (SSAT), encoded by SAT1 (peroxisomal) (no SAT2)
acetyl-CoA + Spd ↔ N1-acetylspermidine (AcSpd) + coenzyme A + H+ 
acetyl-CoA + Spm ↔ N1-acetylspermine (AcSpn) + coenzyme A + H+

Spn ↔ Sdp



Spermine oxidase (SMO), encoded by SMOX, FAD as cofactor
Spn + O2 + H2O ↔ Spd + 3-aminopropanal + H2O2 

AcSpd →

AcSpn →

In addition, export and import systems exist for polyamines [2], using both endocytic and solute carrier transport mechanisms as shown in HCT116 cells [3]. Thus, in HCT116 cells transport into cells has been shown that occurs via caveolar endocytosis – negatively regulated by caveolin-1 – and also export has been identified that occurs via a polyamine/arginine exchange transporter, including the export of acetylated polyamines [3]. SLC3A2 gene encoded one of the subunits of the heterodimeric solute carrier transporter, which also catalyzed the uptake of putrescine by a reverse reaction in certain conditions.

AcSpd ↔ Put

AcSpn ↔ Spd

N1-acetylpolyamine oxidase (APAO), encoded by PAOX
AcSpd + O2 + H2O ↔ Put + 3-acetamidopropanal + H2O2 
AcSpm + O2 + H2O ↔ Spd + 3-acetamidopropanal + H2O2 
"""

#add isoleucine mitochondrial transporter

reaction = Reaction('iletam')
reaction.name = 'isoleucine transaminase, mitochondrial'
reaction.subsystem = 'Valine, Leucine, and Isoleucine Metabolism'
reaction.lower_bound=0
reaction.upper_bound=1000.0
reaction.add_metabolites({ile_m: -1.0, glu_m: 1.0, akg_m: -1.0, n3mop_m: 1.0})
model.add_reaction(reaction)

reaction = Reaction('oivd3m')
reaction.name = ' 2-oxoisovalerate dehydrogenase (acylating; 3-methyl-2-oxopentanoate), mitochondrial'
reaction.subsystem = 'Valine, Leucine, and Isoleucine Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({n3mop_m: -1.0, n2mbcoa_m: 1.0, nadh_m: 1.0, nad_m: -1.0, coa_m: -1.0, co2_m: 1.0})
model.add_reaction(reaction)

reaction = Reaction('acoad10m')
reaction.name = 'acyl-CoA dehydrogenase (2-methylbutanoyl-CoA), mitochondrial'
reaction.subsystem = 'Valine, Leucine, and Isoleucine Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({fad_m: -1.0, fadh2_m: 1.0, n2mbcoa_m: -1.0, n2mb2coa_m: 1.0})
model.add_reaction(reaction)

reaction = Reaction('ecoah9m')
reaction.name = ' 2-Methylprop-2-enoyl-CoA (2-Methylbut-2-enoyl-CoA), mitochondrial'
reaction.subsystem = 'Valine, Leucine, and Isoleucine Metabolism'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({h2o_m: -1.0, n3hmbcoa_m: 1.0, n2mb2coa_m: -1.0})
model.add_reaction(reaction)

reaction = Reaction('hacd9m')
reaction.name = ' 3-hydroxyacyl-CoA dehydrogenase (2-Methylacetoacetyl-CoA), mitochondrial'
reaction.subsystem = 'Valine, Leucine, and Isoleucine Metabolism'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({nad_m: -1.0, h_m: 1.0, n2maacoa_m: 1.0, n3hmbcoa_m: -1.0, nadh_m: 1.0})
model.add_reaction(reaction)

reaction = Reaction('acact10m')
reaction.name = 'acetyl-CoA C-acetyltransferase, mitochondrial'
reaction.subsystem = 'Valine, Leucine, and Isoleucine Metabolism'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({accoa_m: 1.0, ppcoa_m: 1.0, coa_m: -1.0, n2maacoa_m: -1.0})
model.add_reaction(reaction)

reaction = Reaction('ppcoacm')
reaction.name = 'Propionyl-CoA carboxylase, mitochondrial'
reaction.subsystem = 'Valine, Leucine, and Isoleucine Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({hco3_m: -1.0, atp_m: -1.0, adp_m: 1.0, pi_m: 1.0, ppcoa_m: -1.0, h_m: 1.0, mmcoa_s_m: 1.0})
model.add_reaction(reaction)

reaction = Reaction('mmem')
reaction.name = 'methylmalonyl-CoA epimerase/racemase'
reaction.subsystem = 'Valine, Leucine, and Isoleucine Metabolism'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({mmcoa_s_m: 1.0, mmcoa_r_m: -1.0})
model.add_reaction(reaction)

reaction = Reaction('mmmm')
reaction.name = 'methylmalonyl-CoA mutase'
reaction.subsystem = 'Valine, Leucine, and Isoleucine Metabolism'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({succoa_m: 1.0, mmcoa_r_m: -1.0})
model.add_reaction(reaction)

reaction = Reaction('leutam')
reaction.name = 'leucine transaminase, mitochondrial'
reaction.subsystem = 'Valine, Leucine, and Isoleucine Metabolism'
reaction.lower_bound=0
reaction.upper_bound=1000.0
reaction.add_metabolites({glu_m: 1.0, n4mop_m: 1.0, leu_m: -1.0, akg_m: -1.0})
model.add_reaction(reaction)

reaction = Reaction('oivd1m')
reaction.name = ' 2-oxoisovalerate dehydrogenase (acylating; 4-methyl-2-oxopentaoate), mitochondrial'
reaction.subsystem = 'Valine, Leucine, and Isoleucine Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({n4mop_m: -1.0, co2_m: 1.0, nad_m: -1.0, coa_m: -1.0, ivcoa_m: 1.0, nadh_m: 1.0})
model.add_reaction(reaction)

reaction = Reaction('acoad8m')
reaction.name = 'acyl-CoA dehydrogenase (isovaleryl-CoA), mitochondrial'
reaction.subsystem = 'Valine, Leucine, and Isoleucine Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({fad_m: -1.0, fadh2_m: 1.0, n3mb2coa_m: 1.0, ivcoa_m: -1.0})
model.add_reaction(reaction)

reaction = Reaction('mcccrm')
reaction.name = 'methylcrotonoyl-CoA carboxylase, mitochondrial'
reaction.subsystem = 'Valine, Leucine, and Isoleucine Metabolism'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({n3mgcoa_m: 1.0, hco3_m: -1.0, atp_m: -1.0, n3mb2coa_m: -1.0, pi_m: 1.0, adp_m: 1.0, h_m: 1.0})
model.add_reaction(reaction)

reaction = Reaction('mgchrm')
reaction.name = 'methylglutaconyl-CoA hydratase (reversible), mitochondrial'
reaction.subsystem = 'Valine, Leucine, and Isoleucine Metabolism'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({h2o_m: -1.0, n3mgcoa_m: -1.0, hmgcoa_m: 1.0})
model.add_reaction(reaction)

reaction = Reaction('hmglm')
reaction.name = 'hydroxymethylglutaryl-CoA lyase'
reaction.subsystem = 'Cholesterol Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({accoa_m: 1.0, acac_m: 1.0, hmgcoa_m: -1.0})
model.add_reaction(reaction)

reaction = Reaction('ocoat1m')
reaction.name = ' 3-oxoacid CoA-transferase'
reaction.subsystem = 'Valine, Leucine, and Isoleucine Metabolism'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({succoa_m: -1.0, succ_m: 1.0, acac_m: -1.0, aacoa_m: 1.0})
model.add_reaction(reaction)

reaction = Reaction('acact1rm')
reaction.name = 'acetyl-CoA C-acetyltransferase, mitochondrial'
reaction.subsystem = 'Tryptophan metabolism'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({accoa_m: -2.0, coa_m: 1.0, aacoa_m: 1.0})
model.add_reaction(reaction)

reaction = Reaction('saccd3m')
reaction.name = 'saccharopine dehydrogenase (NADP, L-lysine forming), mitochondrial'
reaction.subsystem = 'Lysine Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({lys_m: -1.0, nadph_m: -1.0, h2o_m: 1.0, nadp_m: 1.0, akg_m: -1.0, saccrp_m: 1.0, h_m: -1.0})
model.add_reaction(reaction)

reaction = Reaction('saccd4m')
reaction.name = 'saccharopine dehydrogenase (NADP, L-glutamate forming), mitochondrial'
reaction.subsystem = 'Lysine Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({nadph_m: 1.0, h2o_m: -1.0, nadp_m: -1.0, saccrp_m: -1.0, nl2aadp6sa_m: 1.0, glu_m: 1.0, h_m: 1.0})
model.add_reaction(reaction)

reaction = Reaction('aasad3m')
reaction.name = 'L-aminoadipate-semialdehyde dehydrogenase (NADH), mitochondrial'
reaction.subsystem = 'Lysine Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({h2o_m: -1.0, nl2aadp6sa_m: -1.0, nadh_m: 1.0, nad_m: -1.0, h_m: 2.0, nl2aadp_m: 1.0})
model.add_reaction(reaction)

reaction = Reaction('2amadptm')
reaction.name = 'L-2-aminoadipate shuttle (cytosol/mitochondria)'
reaction.subsystem = 'Lysine Metabolism'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({akg_c: 1.0, nl2aadp_c: -1.0, nl2aadp_m: 1.0, akg_m: -1.0})
model.add_reaction(reaction)

reaction = Reaction('aatai')
reaction.name = ' 2-aminoadipate transaminase, irreversible'
reaction.subsystem = 'Lysine Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({akg_c: -1.0, nl2aadp_c: -1.0, n2oxoadp_c: 1.0, glu_c: 1.0})
model.add_reaction(reaction)

reaction = Reaction('2oxoadptm')
reaction.name = ' 2-oxoadipate shuttle (cytosol/mitochondria)'
reaction.subsystem = 'Tryptophan metabolism'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({akg_c: 1.0, n2oxoadp_m: 1.0, n2oxoadp_c: -1.0, akg_m: -1.0})
model.add_reaction(reaction)

reaction = Reaction('2oxoadoxm')
reaction.name = ' 2-Oxoadipate:lipoamde 2-oxidoreductase(decarboxylating and acceptor-succinylating) (mitochondria)'
reaction.subsystem = 'Lysine Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({n2oxoadp_m: -1.0, nadh_m: 1.0, nad_m: -1.0, coa_m: -1.0, glutcoa_m: 1.0, co2_m: 1.0})
model.add_reaction(reaction)

reaction = Reaction('glutcoadhm')
reaction.name = 'glutaryl-CoA dehydrogenase (mitochondria)'
reaction.subsystem = 'Tryptophan metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({b2coa_m: 1.0, fadh2_m: 1.0, fad_m: -1.0, co2_m: 1.0, h_m: -1.0, glutcoa_m: -1.0})
model.add_reaction(reaction)

reaction = Reaction('ecoah1m')
reaction.name = ' 3-hydroxyacyl-CoA dehydratase (3-hydroxybutanoyl-CoA) (mitochondria)'
reaction.subsystem = 'Tryptophan metabolism'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({h2o_m: 1.0, b2coa_m: 1.0, n3hbcoa_m: -1.0})
model.add_reaction(reaction)

reaction = Reaction('hacd1m')
reaction.name = ' 3-hydroxyacyl-CoA dehydrogenase (acetoacetyl-CoA) (mitochondria)'
reaction.subsystem = 'Tryptophan metabolism'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({nad_m: 1.0, h_m: -1.0, n3hbcoa_m: 1.0, nadh_m: -1.0, aacoa_m: -1.0})
model.add_reaction(reaction)

#degradacio threonina diferent pedro

reaction = Reaction('thrd_l')
reaction.name = 'L-threonine deaminase'
reaction.subsystem = 'Glycine, Serine, and Threonine Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({thr_c: -1.0, n2obut_c: 1.0, nh4_c: 1.0})
model.add_reaction(reaction)

reaction = Reaction('obdhc')
reaction.name = ' 2-Oxobutanoate dehydrogenase, cytosolic'
reaction.subsystem = 'Glycine, Serine, and Threonine Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({nadh_c: 1.0, nad_c: -1.0, n2obut_c: -1.0, coa_c: -1.0, ppcoa_c: 1.0, co2_c: 1.0})
model.add_reaction(reaction)

#Modified with respect to recon2
reaction = Reaction('ppcoamtr')
reaction.name = 'Postulated transport reaction for propanoyl-CoA'
reaction.subsystem = 'Glycine, Serine, and Threonine Metabolism'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({ppcoa_c: -1.0, ppcoa_m: 1.0,coa_c:1,coa_m:-1})
model.add_reaction(reaction)

#Tyrosine degradation
reaction = Reaction('tyrta')
reaction.name = 'tyrosine transaminase'
reaction.subsystem = 'Tyrosine metabolism'
reaction.lower_bound=0 #Assumed irreversible but in recon1 is reverseible 
reaction.upper_bound=1000.0
reaction.add_metabolites({akg_c: -1.0, glu_c: 1.0, n34hpp_c: 1.0, tyr_c: -1.0})
model.add_reaction(reaction)

reaction = Reaction('34hppor')
reaction.name = ' 4-Hydroxyphenylpyruvate:oxygen oxidoreductase'
reaction.subsystem = 'Tyrosine metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({hgentis_c: 1.0, n34hpp_c: -1.0, co2_c: 1.0, o2_c: -1.0})
model.add_reaction(reaction)

reaction = Reaction('hgntor')
reaction.name = 'Homogentisate:oxygen 1,2-oxidoreductase (decyclizing)'
reaction.subsystem = 'Tyrosine metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({hgentis_c: -1.0, n4mlacac_c: 1.0, h_c: 1.0, o2_c: -1.0})
model.add_reaction(reaction)

reaction = Reaction('macaci')
reaction.name = 'maleylacetoacetate isomerase'
reaction.subsystem = 'Tyrosine metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({n4fumacac_c: 1.0, n4mlacac_c: -1.0})
model.add_reaction(reaction)

reaction = Reaction('fumac')
reaction.name = 'fumarylacetoacetase'
reaction.subsystem = 'Tyrosine metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({acac_c: 1.0, n4fumacac_c: -1.0, h2o_c: -1.0, fum_c: 1.0, h_c: 1.0})
model.add_reaction(reaction)

reaction = Reaction('fum')
reaction.name = 'fumarase'
reaction.subsystem = 'Citric Acid Cycle'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({h2o_c: -1.0, fum_c: -1.0, mal_c
: 1.0})
model.add_reaction(reaction)

#Fate 1 of acetoactetate: mitochondrial acoA
reaction = Reaction('acact2m')
reaction.name = 'Acetoacetate mitochondrial transport via H+ symport'
reaction.subsystem = 'Transport, Mitochondrial'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({acac_c: -1.0, h_m: 1.0, h_c: -1.0, acac_m: 1.0})
model.add_reaction(reaction)
#Fate 2 of aceatoactetate cytosolic acoa, this is omited to prevent redundancy (acetyl coa can still be exported to the cytosol to cytrate liase with the same enrgetic cost)  

#tryptophan metabolism: 
reaction = Reaction('trpo2')
reaction.name = 'L-Tryptophan:oxygen 2,3-oxidoreductase (decyclizing)'
reaction.subsystem = 'Tryptophan metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({lfmkynr_c: 1.0, o2_c: -1.0, trp_c: -1.0})
model.add_reaction(reaction)

reaction = Reaction('fkynh')
reaction.name = 'N-Formyl-L-kynurenine amidohydrolase'
reaction.subsystem = 'Tryptophan metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({lfmkynr_c: -1.0, h2o_c: -1.0, lkynr_c: 1.0, for_c: 1.0, h_c: 1.0})
model.add_reaction(reaction)

reaction = Reaction('kyn3ox')
reaction.name = 'kynurenine 3-monooxygenase'
reaction.subsystem = 'Tryptophan metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({o2_c: -1.0, h2o_c: 1.0, nadp_c: 1.0, nadph_c: -1.0, hlkynr_c: 1.0, h_c: -1.0, lkynr_c: -1.0})
model.add_reaction(reaction)

reaction = Reaction('hkynh')
reaction.name = ' 3-Hydroxy-L-kynurenine hydrolase'
reaction.subsystem = 'Tryptophan metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({n3hanthrn_c: 1.0, ala_c: 1.0, hlkynr_c: -1.0, h2o_c: -1.0})
model.add_reaction(reaction)

reaction = Reaction('3hao')
reaction.name = ' 3-hydroxyanthranilate 3,4-dioxygenase'
reaction.subsystem = 'Tryptophan metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({n3hanthrn_c: -1.0, o2_c: -1.0, cmusa_c: 1.0, h_c: 1.0})
model.add_reaction(reaction)

reaction = Reaction('pclad')
reaction.name = 'picolinic acid decarboxylase'
reaction.subsystem = 'Tryptophan metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({cmusa_c: -1.0, co2_c: 1.0, am6sa_c: 1.0, h_c: -1.0})
model.add_reaction(reaction)

reaction = Reaction('am6sad')
reaction.name = 'aminomuconate-semialdehyde dehydrogenase'
reaction.subsystem = 'Tryptophan metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({h2o_c: -1.0, amuco_c: 1.0, am6sa_c: -1.0, h_c: 2.0, nad_c: -1.0, nadh_c: 1.0})
model.add_reaction(reaction)

reaction = Reaction('amcoxo')
reaction.name = ' 2-aminomuconate reductase'
reaction.subsystem = 'Tryptophan metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({h2o_c: -1.0, nadp_c: 1.0, nadph_c: -1.0, n2oxoadp_c: 1.0, h_c: -1.0, amuco_c: -1.0, nh4_c: 1.0})
model.add_reaction(reaction)

#Cysteine metabolism
reaction = Reaction('cyso')
reaction.name = 'cysteine oxidase'
reaction.subsystem = 'Taurine and hypotaurine metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({o2_c: -1.0, cys_c: -1.0, h_c: 2.0, n3sala_c: 1.0})
model.add_reaction(reaction)

reaction = Reaction('3salatai')
reaction.name = ' 3-sulfino-alanine transaminase (irreversible)'
reaction.subsystem = 'Cysteine Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({akg_c: -1.0, glu_c: 1.0, n3snpyr_c: 1.0, h_c: -1.0, n3sala_c: -1.0})
model.add_reaction(reaction)

reaction = Reaction('3spyrsp')
reaction.name = ' 3-sulfinopyruvate hydrolase (spotaneous reaction)'
reaction.subsystem = 'Cysteine Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({h2o_c: -1.0, n3snpyr_c: -1.0, pyr_b_c: 1.0, so3_c: 1.0, h_c: 1.0})
model.add_reaction(reaction)

reaction = Reaction('sulfox')
reaction.name = 'sulfite oxidase'
reaction.subsystem = 'Cysteine Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({h2o_c: -1.0, so3_c: -1.0, ficytc_m: -2.0, so4_c: 1.0, h_c: 2.0, focytc_m: 2.0})
model.add_reaction(reaction)

#Methionine 
reaction = Reaction('cysts')
reaction.name = 'cystathionine beta-synthase'
reaction.subsystem = 'Methionine Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({ser_c: -1.0, hcys_c: -1.0, h2o_c: 1.0, cyst_c: 1.0})
model.add_reaction(reaction)

reaction = Reaction('cystgl')
reaction.name = 'cystathionine g-lyase'
reaction.subsystem = 'Cysteine Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({h2o_c: -1.0, n2obut_c: 1.0, nh4_c: 1.0, cyst_c: -1.0, cys_c: 1.0})
model.add_reaction(reaction)

reaction = Reaction('valtam')
reaction.name = 'valine transaminase, mitochondiral'
reaction.subsystem = 'Valine, Leucine, and Isoleucine Metabolism'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({n3mob_m: 1.0, glu_m: 1.0, akg_m: -1.0, val_m: -1.0})
model.add_reaction(reaction)

"""reaction = Reaction('valta')
reaction.name = 'valine transaminase'
reaction.subsystem = 'Valine, Leucine, and Isoleucine Metabolism'
reaction.lower_bound=0 #assumed irreversible but in recon1 is reversibles 
reaction.upper_bound=1000.0
reaction.add_metabolites({akg_c: -1.0, glu_c: 1.0, val_c: -1.0, n3mob_c: 1.0})
model.add_reaction(reaction)"""

"""
reaction = Reaction('valtam')
reaction.name = 'valine transaminase, mitochondiral'
reaction.subsystem = 'Valine, Leucine, and Isoleucine Metabolism'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({'3mob_m': 1.0, 'glu_dash_l_m': 1.0, 'akg_m': -1.0, 'val_dash_l_m': -1.0})
model.add_reaction(reaction)

reaction = Reaction('oivd2m')
reaction.name = ' 2-oxoisovalerate dehydrogenase (acylating; 3-methyl-2-oxobutanoate), mitochondrial'
reaction.subsystem = 'Valine, Leucine, and Isoleucine Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({'3mob_m': -1.0, 'ibcoa_m': 1.0, 'nadh_m': 1.0, 'nad_m': -1.0, 'coa_m': -1.0, 'co2_m': 1.0})
model.add_reaction(reaction)

reaction = Reaction('acoad9m')
reaction.name = 'acyl-CoA dehydrogenase (isobutyryl-CoA), mitochondrial'
reaction.subsystem = 'Valine, Leucine, and Isoleucine Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({'fad_m': -1.0, 'ibcoa_m': -1.0, 'fadh2_m': 1.0, '2mp2coa_m': 1.0})
model.add_reaction(reaction)

reaction = Reaction('ecoah12m')
reaction.name = ' 3-hydroxyacyl-CoA dehydratase (3-hydroxyisobutyryl-CoA) (mitochondria)'
reaction.subsystem = 'Valine, Leucine, and Isoleucine Metabolism'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({'3hibutcoa_m': 1.0, 'h2o_m': -1.0, '2mp2coa_m': -1.0})
model.add_reaction(reaction)

reaction = Reaction('3hbcoahlm')
reaction.name = ' 3-hydroxyisobutyryl-CoA hydrolase, mitochondrial'
reaction.subsystem = 'Valine, Leucine, and Isoleucine Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({'3hmp_m': 1.0, '3hibutcoa_m': -1.0, 'coa_m': 1.0, 'h_m': 1.0, 'h2o_m': -1.0})
model.add_reaction(reaction)

reaction = Reaction('hibdm')
reaction.name = ' 3-hydroxyisobutyrate dehydrogenase, mitochondrial'
reaction.subsystem = 'Valine, Leucine, and Isoleucine Metabolism'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({'3hmp_m': -1.0, 'nad_m': -1.0, 'h_m': 1.0, 'nadh_m': 1.0, '2mop_m': 1.0})
model.add_reaction(reaction)

reaction = Reaction('mmsad1m')
reaction.name = 'methylmalonate-semialdehyde dehydrogenase'
reaction.subsystem = 'Valine, Leucine, and Isoleucine Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({'2mop_m': -1.0, 'ppcoa_m': 1.0, 'nadh_m': 1.0, 'nad_m': -1.0, 'coa_m': -1.0, 'co2_m': 1.0})
model.add_reaction(reaction)

"""

"""reaction = Reaction('3mobt2im')
reaction.name = ' 3-methyl-2-oxobutanoate mitochondrial transport via proton symport'
reaction.subsystem = 'Transport, Mitochondrial'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({n3mob_m: 1.0, h_m: 1.0, h_c: -1.0, n3mob_c: -1.0})
model.add_reaction(reaction)"""

reaction = Reaction('oivd2m')
reaction.name = ' 2-oxoisovalerate dehydrogenase (acylating; 3-methyl-2-oxobutanoate), mitochondrial'
reaction.subsystem = 'Valine, Leucine, and Isoleucine Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({n3mob_m: -1.0, ibcoa_m: 1.0, nadh_m: 1.0, nad_m: -1.0, coa_m: -1.0, co2_m: 1.0})
model.add_reaction(reaction)

reaction = Reaction('acoad9m')
reaction.name = 'acyl-CoA dehydrogenase (isobutyryl-CoA), mitochondrial'
reaction.subsystem = 'Valine, Leucine, and Isoleucine Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({fad_m: -1.0, ibcoa_m: -1.0, fadh2_m: 1.0, n2mp2coa_m: 1.0})
model.add_reaction(reaction)

reaction = Reaction('ecoah12m')
reaction.name = ' 3-hydroxyacyl-CoA dehydratase (3-hydroxyisobutyryl-CoA) (mitochondria)'
reaction.subsystem = 'Valine, Leucine, and Isoleucine Metabolism'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({n3hibutcoa_m: 1.0, h2o_m: -1.0, n2mp2coa_m: -1.0})
model.add_reaction(reaction)

reaction = Reaction('3hbcoahlm')
reaction.name = ' 3-hydroxyisobutyryl-CoA hydrolase, mitochondrial'
reaction.subsystem = 'Valine, Leucine, and Isoleucine Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({n3hmp_m: 1.0, n3hibutcoa_m: -1.0, coa_m: 1.0, h_m: 1.0, h2o_m: -1.0})
model.add_reaction(reaction)

reaction = Reaction('hibdm')
reaction.name = ' 3-hydroxyisobutyrate dehydrogenase, mitochondrial'
reaction.subsystem = 'Valine, Leucine, and Isoleucine Metabolism'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({n3hmp_m: -1.0, nad_m: -1.0, h_m: 1.0, nadh_m: 1.0, n2mop_m: 1.0})
model.add_reaction(reaction)

"""reaction = Reaction('mmtsadm')
reaction.name = 'malonate-semialdehyde dehydrogenase (acetylating), mitochondrial'
reaction.subsystem = 'Valine, Leucine, and Isoleucine Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({n2mop_m: -1.0, nadh_m: 1.0, nad_m: -1.0, coa_m: -1.0, h_m: 1.0, mmcoa_r_m: 1.0})
model.add_reaction(reaction)"""

reaction = Reaction('mmsad1m')
reaction.name = 'methylmalonate-semialdehyde dehydrogenase'
reaction.subsystem = 'Valine, Leucine, and Isoleucine Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({n2mop_m: -1.0, ppcoa_m: 1.0, nadh_m: 1.0, nad_m: -1.0, coa_m: -1.0, co2_m: 1.0})
model.add_reaction(reaction)

reaction = Reaction('phethptox2')
reaction.name = 'L-Phenylalanine,tetrahydrobiopterin:oxygen oxidoreductase (4-hydroxylating)'
reaction.subsystem = 'Tyr, Phe, Trp Biosynthesis'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({thbpt_c: -1.0, o2_c: -1.0, phe_c: -1.0, thbpt4acam_c: 1.0, tyr_c: 1.0})
model.add_reaction(reaction)

reaction = Reaction('thbpt4acamdase')
reaction.name = 'Tetrahydrobiopterin-4a-carbinolamine dehydratase'
reaction.subsystem = 'Tetrahydrobiopterin'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({h2o_c: 1.0, thbpt4acam_c: -1.0, dhbpt_c: 1.0})
model.add_reaction(reaction)

reaction = Reaction('dhpr')
reaction.name = ' 6,7-dihydropteridine reductase'
reaction.subsystem = 'Tetrahydrobiopterin'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({thbpt_c: 1.0, nadh_c: -1.0, h_c: -1.0, nad_c: 1.0, dhbpt_c: -1.0})
model.add_reaction(reaction)

#Beta oxidation: 

reaction = Reaction('hdcatr')
reaction.name = 'fatty acid transport via diffusion'
reaction.subsystem = 'Transport, Extracellular'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({hdca_e: -1.0, hdca_c: 1.0})
model.add_reaction(reaction)

reaction = Reaction('facoal160i')
reaction.name = 'C160 fatty acid activation'
reaction.subsystem = 'Fatty acid activation'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({atp_c: -1.0, ppi_c: 1.0, pmtcoa_c: 1.0, coa_c: -1.0, amp_c: 1.0, hdca_c: -1.0})
model.add_reaction(reaction)

reaction = Reaction('c160cpt1')
reaction.name = 'carnitine O-palmitoyltransferase'
reaction.subsystem = 'Carnitine shuttle'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({crn_c: -1.0, pmtcrn_c: 1.0, coa_c: 1.0, pmtcoa_c: -1.0})
model.add_reaction(reaction)

reaction = Reaction('c160crnt')
reaction.name = 'C160 transport into the mitochondria'
reaction.subsystem = 'Carnitine shuttle'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({pmtcrn_c: -1.0, pmtcrn_m: 1.0})
model.add_reaction(reaction)

reaction = Reaction('c160cpt2')
reaction.name = 'C160 transport into the mitochondria'
reaction.subsystem = 'Carnitine shuttle'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({coa_m: -1.0, crn_m: 1.0, pmtcoa_m: 1.0, pmtcrn_m: -1.0})
model.add_reaction(reaction)

reaction = Reaction('crntim')
reaction.name = 'L-carnitine transport out of mitochondria via diffusion'
reaction.subsystem = 'Transport, Mitochondrial'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({crn_c: 1.0, crn_m: -1.0})
model.add_reaction(reaction)

reaction = Reaction('faoxc16080m')
reaction.name = 'Beta oxidation of long chain fatty acid'
reaction.subsystem = 'Fatty acid oxidation'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({h2o_m: -4.0, fad_m: -4.0, accoa_m: 4.0, occoa_m: 1.0, nadh_m: 4.0, nad_m: -4.0, coa_m: -4.0, h_m: 4.0, fadh2_m: 4.0, pmtcoa_m: -1.0})
model.add_reaction(reaction)

reaction = Reaction('faoxc80')
reaction.name = 'Beta oxidation of med/long chain fatty acid'
reaction.subsystem = 'Fatty acid oxidation'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({h2o_m: -3.0, fad_m: -3.0, accoa_m: 4.0, occoa_m: -1.0, nadh_m: 3.0, nad_m: -3.0, coa_m: -3.0, h_m: 3.0, fadh2_m: 3.0})
model.add_reaction(reaction)

#####Nucleotodies Biosynthesis
#http://biocyc.org/HUMAN/NEW-IMAGE?type=PATHWAY&object=PWY-841&detail-level=3
#http://biocyc.org/HUMAN/NEW-IMAGE?type=PATHWAY&object=PWY0-162 

reaction = Reaction('prpps')
reaction.name = 'phosphoribosylpyrophosphate synthetase'
reaction.subsystem = 'Pentose Phosphate Pathway'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({prpp_c: 1.0, rib5p_c: -1.0, amp_c: 1.0, atp_c: -1.0, h_c: 1.0})
model.add_reaction(reaction)

#Consumes prpp_c,10fthf_c
#IMP Biosinthesys
reaction = Reaction('gluprt')
reaction.name = 'glutamine phosphoribosyldiphosphate amidotransferase'
reaction.subsystem = 'IMP Biosynthesis'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({pram_c: 1.0, h2o_c: -1.0, ppi_c: 1.0, gln_c: -1.0, glu_c: 1.0, prpp_c: -1.0})
model.add_reaction(reaction)

reaction = Reaction('pragsr')
reaction.name = 'phosphoribosylglycinamide synthase'
reaction.subsystem = 'IMP Biosynthesis'
reaction.lower_bound=0 #Assumed irreversible
reaction.upper_bound=1000.0
reaction.add_metabolites({pram_c: -1.0, atp_c: -1.0, h_c: 1.0, adp_c: 1.0, gar_c: 1.0, pi_c: 1.0, gly_c: -1.0})
model.add_reaction(reaction)

reaction = Reaction('garft')
reaction.name = 'phosphoribosylglycinamide formyltransferase'
reaction.subsystem = 'IMP Biosynthesis'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({gar_c: -1.0, n10fthf_c: -1.0, fgam_c: 1.0, h_c: 1.0, thf_c: 1.0})
model.add_reaction(reaction)

reaction = Reaction('prfgs')
reaction.name = 'phosphoribosylformylglycinamidine synthase'
reaction.subsystem = 'IMP Biosynthesis'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({h2o_c: -1.0, atp_c: -1.0, fgam_c: -1.0, gln_c: -1.0, h_c: 1.0, adp_c: 1.0, fpram_c: 1.0, glu_c: 1.0, pi_c: 1.0})
model.add_reaction(reaction)

reaction = Reaction('prais')
reaction.name = 'phosphoribosylaminoimidazole synthase'
reaction.subsystem = 'IMP Biosynthesis'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({atp_c: -1.0, fpram_c: -1.0, h_c: 2.0, adp_c: 1.0, air_c: 1.0, pi_c: 1.0})
model.add_reaction(reaction)

reaction = Reaction('aircr')
reaction.name = 'phosphoribosylaminoimidazole carboxylase'
reaction.subsystem = 'IMP Biosynthesis'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({air_c: -1.0, co2_c: -1.0, h_c: 1.0, n5aizc_c: 1.0})
model.add_reaction(reaction)

reaction = Reaction('prascs')
reaction.name = 'phosphoribosylaminoimidazolesuccinocarboxamide synthase'
reaction.subsystem = 'IMP Biosynthesis'
reaction.lower_bound=0 #Converted to irreversible, originally was reversible
reaction.upper_bound=1000.0
reaction.add_metabolites({atp_c: -1.0, n5aizc_c: -1.0, n25aics_c: 1.0, h_c: 1.0, adp_c: 1.0, pi_c: 1.0, asp_c: -1.0})
model.add_reaction(reaction)

reaction = Reaction('adsl2')
reaction.name = 'adenylosuccinate lyase'
reaction.subsystem = 'IMP Biosynthesis'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({aicar_c: 1.0, fum_c: 1.0, n25aics_c: -1.0})
model.add_reaction(reaction)


reaction = Reaction('aicart')
reaction.name = 'phosphoribosylaminoimidazolecarboxamide formyltransferase'
reaction.subsystem = 'IMP Biosynthesis'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({aicar_c: -1.0, fprica_c: 1.0, n10fthf_c: -1.0, thf_c: 1.0})
model.add_reaction(reaction)


reaction = Reaction('impc')
reaction.name = 'IMP cyclohydrolase'
reaction.subsystem = 'IMP Biosynthesis'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({h2o_c: -1.0, fprica_c: 1.0, imp_c: -1.0})
model.add_reaction(reaction)


#AMP Biosynthesis 
reaction = Reaction('adss')
reaction.name = 'adenylosuccinate synthase'
reaction.subsystem = 'Nucleotides'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({gdp_c: 1.0, gtp_c: -1.0, h_c: 2.0, pi_c: 1.0, asp_c: -1.0, imp_c: -1.0, dcamp_c: 1.0})
model.add_reaction(reaction)

reaction = Reaction('adsl1')
reaction.name = 'adenylosuccinate lyase'
reaction.subsystem = 'Nucleotides'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({fum_c: 1.0, amp_c: 1.0, dcamp_c: -1.0})
model.add_reaction(reaction)

#GMP Biosynyhesis 

reaction = Reaction('impd')
reaction.name = 'IMP dehydrogenase'
reaction.subsystem = 'Nucleotides'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({xmp_c: 1.0, h2o_c: -1.0, h_c: 1.0, nad_c: -1.0, nadh_c: 1.0, imp_c: -1.0})
model.add_reaction(reaction)

reaction = Reaction('gmps2')
reaction.name = 'GMP synthase'
reaction.subsystem = 'Nucleotides'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({xmp_c: -1.0, h2o_c: -1.0, gmp_c: 1.0, atp_c: -1.0, ppi_c: 1.0, gln_c: -1.0, h_c: 2.0, glu_c: 1.0, amp_c: 1.0})
model.add_reaction(reaction)


"""reaction = Reaction('tmds')
reaction.name = 'thymidylate synthase'
reaction.subsystem = 'Nucleotides'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({dtmp_c: 1.0, dump_c: -1.0, dhf_c: 1.0, mlthf_c: -1.0})
model.add_reaction(reaction)"""

#Pyrimidine Biosythesis 

reaction = Reaction('cbps')
reaction.name = 'carbamoyl-phosphate synthase (glutamine-hydrolysing)'
reaction.subsystem = 'Pyrimidine Biosynthesis'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({h2o_c: -1.0, atp_c: -2.0, hco3_c: -1.0, cbp_c: 1.0, gln_c: -1.0, h_c: 2.0, adp_c: 2.0, glu_c: 1.0, pi_c: 1.0})
model.add_reaction(reaction)

"""reaction = Reaction('aspctr')
reaction.name = 'aspartate carbamoyltransferase (reversible)'
reaction.subsystem = 'Pyrimidine Biosynthesis'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({cbasp_c: 1.0, cbp_c: -1.0, pi_c: 1.0, asp_c: -1.0, h_c: 1.0})
model.add_reaction(reaction)"""

#Converted to irreversible
reaction = Reaction('aspct')
reaction.name = 'aspartate carbamoyltransferase'
reaction.subsystem = 'Pyrimidine Biosynthesis'
reaction.lower_bound=0
reaction.upper_bound=1000.0
reaction.add_metabolites({cbasp_c: 1.0, cbp_c: -1.0, pi_c: 1.0, asp_c: -1.0, h_c: 1.0})
model.add_reaction(reaction)

reaction = Reaction('dhorts')
reaction.name = 'dihydroorotase'
reaction.subsystem = 'Pyrimidine Biosynthesis'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({cbasp_c: 1.0, h2o_c: -1.0, h_c: 1.0, dhor_c: -1.0})
model.add_reaction(reaction)

reaction = Reaction('dhord9')
reaction.name = 'dihydoorotic acid dehydrogenase (quinone10)'
reaction.subsystem = 'Pyrimidine Biosynthesis'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({dhor_c: -1.0, q10h2_m: 1.0, orot_c: 1.0, q10_m: -1.0})
model.add_reaction(reaction)


reaction = Reaction('orpt')
reaction.name = 'orotate phosphoribosyltransferase'
reaction.subsystem = 'Pyrimidine Biosynthesis'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({ppi_c: -1.0, prpp_c: 1.0, orot_c: 1.0, orot5p_c: -1.0})
model.add_reaction(reaction)

#ADD URIDINE INPUT?
reaction = Reaction('orpt')
reaction.name = 'uridine kinase (ATP:Uridine)'
reaction.subsystem = 'Pyrimidine Biosynthesis'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({ump_c: 1.0, uri_c: -1.0, atp_c: -1.0, adp_c: 1.0, h_c: 1.0})
model.add_reaction(reaction)

reaction = Reaction('ompdc')
reaction.name = 'orotidine-5-phosphate decarboxylase'
reaction.subsystem = 'Pyrimidine Biosynthesis'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({ump_c: 1.0, co2_c: 1.0, h_c: -1.0, orot5p_c: -1.0})
model.add_reaction(reaction)

reaction = Reaction('umpk')
reaction.name = 'UMP kinase'
reaction.subsystem = 'Nucleotides'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({ump_c: -1.0, udp_c: 1.0, atp_c: -1.0, adp_c: 1.0})
model.add_reaction(reaction)

reaction = Reaction('ctps2')
reaction.name = 'CTP synthase (glutamine)'
reaction.subsystem = 'Pyrimidine Biosynthesis'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({utp_c: -1.0, atp_c: -1.0, h2o_c: -1.0, gln_c: -1.0, h_c: 2.0, adp_c: 1.0, ctp_c: 1.0, glu_c: 1.0, pi_c: 1.0})
model.add_reaction(reaction)

#ribonucleoside-diphosphate reductase

reaction = Reaction('rndr1')
reaction.name = 'ribonucleoside-diphosphate reductase (ADP)'
reaction.subsystem = 'Nucleotides'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({dadp_c: 1.0, trdox_c: 1.0, h2o_c: 1.0, adp_c: -1.0, trdrd_c: -1.0})
model.add_reaction(reaction)

reaction = Reaction('rndr3')
reaction.name = 'ribonucleoside-diphosphate reductase (CDP)'
reaction.subsystem = 'Nucleotides'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({cdp_c: -1.0, trdox_c: 1.0, trdrd_c: -1.0, h2o_c: 1.0, dcdp_c: 1.0})
model.add_reaction(reaction)

reaction = Reaction('rndr2')
reaction.name = 'ribonucleoside-diphosphate reductase (GDP)'
reaction.subsystem = 'Nucleotides'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({h2o_c: 1.0, dgdp_c: 1.0, gdp_c: -1.0, trdox_c: 1.0, trdrd_c: -1.0})
model.add_reaction(reaction)

reaction = Reaction('rndr4')
reaction.name = 'ribonucleoside-diphosphate reductase (UDP)'
reaction.subsystem = 'Nucleotides'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({trdrd_c: -1.0, h2o_c: 1.0, trdox_c: 1.0, udp_c: -1.0, dudp_c: 1.0})
model.add_reaction(reaction)

reaction = Reaction('ndp8')
reaction.name = 'nucleoside-diphosphatase (dUDP)'
reaction.subsystem = 'Nucleotides'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({h2o_c: -1.0, pi_c: 1.0, dudp_c: -1.0, dump_c: 1.0, h_c: 1.0})
model.add_reaction(reaction)

reaction = Reaction('trdr')
reaction.name = 'thioredoxin reductase (NADPH)'
reaction.subsystem = 'Nucleotides'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({trdox_c: -1.0, nadp_c: 1.0, nadph_c: -1.0, h_c: -1.0, trdrd_c: 1.0})
model.add_reaction(reaction)

#Nucleotides Kinasses 
reaction = Reaction('ndpk5')
reaction.name = 'nucleoside-diphosphate kinase (ATP:dGDP)'
reaction.subsystem = 'Nucleotides'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({dgdp_c: -1.0, dgtp_c: 1.0, atp_c: -1.0, adp_c: 1.0})
model.add_reaction(reaction)

reaction = Reaction('ndpk7')
reaction.name = 'nucleoside-diphosphate kinase (ATP:dCDP)'
reaction.subsystem = 'Nucleotides'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({dctp_c: 1.0, atp_c: -1.0, adp_c: 1.0, dcdp_c: -1.0})
model.add_reaction(reaction)

reaction = Reaction('ndpk4')
reaction.name = 'nucleoside-diphosphate kinase (ATP:dTDP)'
reaction.subsystem = 'Nucleotides'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({dtdp_c: -1.0, dttp_c: 1.0, atp_c: -1.0, adp_c: 1.0})
model.add_reaction(reaction)

reaction = Reaction('ndpk3')
reaction.name = 'nucleoside-diphosphate kinase (ATP:CDP)'
reaction.subsystem = 'Nucleotides'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({ctp_c: 1.0, cdp_c: -1.0, atp_c: -1.0, adp_c: 1.0})
model.add_reaction(reaction)

reaction = Reaction('ndpk8')
reaction.name = 'nucleoside-diphosphate kinase (ATP:dADP)'
reaction.subsystem = 'Nucleotides'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({dadp_c: -1.0, atp_c: -1.0, adp_c: 1.0, datp_c: 1.0})
model.add_reaction(reaction)

reaction = Reaction('gk1')
reaction.name = 'guanylate kinase (GMP:ATP)'
reaction.subsystem = 'Nucleotides'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({gmp_c: -1.0, gdp_c: 1.0, atp_c: -1.0, adp_c: 1.0})
model.add_reaction(reaction)

"""reaction = Reaction('uridk2m')
reaction.name = 'uridylate kinase (dUMP), mitochondrial'
reaction.subsystem = 'Nucleotides'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({atp_m: -1.0, dump_m: -1.0, adp_m: 1.0, dudp_m: 1.0})
model.add_reaction(reaction)"""

reaction = Reaction('ocdcatr')
reaction.name = 'fatty acid transport via diffusion'
reaction.subsystem = 'Transport, Extracellular'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({ocdca_e: -1.0, ocdca_c: 1.0})
model.add_reaction(reaction)





#Phosphatydyl choline 
reaction = Reaction('choltu')
reaction.name = 'Choline uniport'
reaction.subsystem = 'Transport, Extracellular'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({chol_c: 1.0,chol_e: -1.0})
model.add_reaction(reaction)

reaction = Reaction('cholk')
reaction.name = 'Choline kinase'
reaction.subsystem = 'Glycerophospholipid Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({chol_c: -1.0,cholp_c: 1.0,h_c: 1.0,adp_c: 1.0,atp_c: -1.0})
model.add_reaction(reaction)

reaction = Reaction('chlpctd')
reaction.name = 'choline phosphate cytididyltransferase'
reaction.subsystem = 'Glycerophospholipid Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({ctp_c: -1.0,cholp_c: -1.0,ppi_c: 1.0,cdpchol_c: 1.0,h_c: -1.0})
model.add_reaction(reaction)

reaction = Reaction('ceptc')
reaction.name = 'choline phosphotransferase'
reaction.subsystem = 'Glycerophospholipid Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({pchol_hs_c: 1.0,dag_hs_c: -1.0,cdpchol_c: -1.0,h_c: 1.0,cmp_c: 1.0})
model.add_reaction(reaction)

#Modified from recon1, removed energy, requirement
reaction = Reaction('pct')
reaction.name = 'phosphatidylcholine transporter'
reaction.subsystem = 'Transport, Extracellular'
reaction.lower_bound=0
reaction.upper_bound=1000.0
reaction.add_metabolites({pchol_hs_c: -1.0, pchol_hs_e: 1.0, atp_c: -1.0, h_c: 1.0, adp_c: 1.0, pi_c: 1.0})
model.add_reaction(reaction)



reaction = Reaction('cytk1')
reaction.name = 'cytidylate kinase (CMP)'
reaction.subsystem = 'Nucleotides'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({cdp_c: 1.0,adp_c: 1.0,atp_c: -1.0,cmp_c: -1.0})
model.add_reaction(reaction)

reaction = Reaction('rtot1')
reaction.name = 'R total flux'
reaction.subsystem = 'R Group Synthesis'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({rtotalcoa_c: 1.0, pmtcoa_c: -0.186,hdcoa_c:-0.074,odecoa_c:-0.629,stcoa_c:-0.111})
model.add_reaction(reaction)

reaction = Reaction('facoal181i')
reaction.name = 'C18:1 fatty acid activation'
reaction.subsystem = 'Fatty acid activation'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({ocdcea_c: -1.0, odecoa_c: 1.0, atp_c: -1.0, ppi_c: 1.0, coa_c: -1.0, amp_c: 1.0})
model.add_reaction(reaction)

reaction = Reaction('facoal180i')
reaction.name = 'C180 fatty acid activation'
reaction.subsystem = 'Fatty acid activation'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({atp_c: -1.0, ppi_c: 1.0, stcoa_c: 1.0, coa_c: -1.0, amp_c: 1.0, ocdca_c: -1.0})
model.add_reaction(reaction)

#Stoichiometric coefficents based on relative abundance of aminoacids from Kashif Sheikh et al 2005, C16:1 0.066,C16:0 0.165,C18:1 0.559,C18:0 0.099

"""
Modified. from recon 1, original version 
reaction = Reaction('rtot1')
reaction.name = 'R total flux'
reaction.subsystem = 'R Group Synthesis'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({'rtotalcoa_c': 1.0, 'r1coa_hs_c': -1.0})
model.add_reaction(reaction)


"""



"""reaction = Reaction('artplm1')
reaction.name = 'R group to palmitate conversion'
reaction.subsystem = 'R Group Synthesis'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({pmtcoa_c: 1.0,rtotalcoa_c: -1.0})
model.add_reaction(reaction)
Removed from Recon1
"""


reaction = Reaction('gpam_hs')
reaction.name = 'glycerol-3-phosphate acyltransferase'
reaction.subsystem = 'Triacylglycerol Synthesis'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({rtotalcoa_c: -1.0,coa_c: 1.0,glyc3p_c: -1.0,alpa_hs_c: 1.0})
model.add_reaction(reaction)

reaction = Reaction('agpat1')
reaction.name = ' 1-acylglycerol-3-phosphate O-acyltransferase 1'
reaction.subsystem = 'Triacylglycerol Synthesis'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({pa_hs_c: 1.0,rtotalcoa_c: -1.0,coa_c: 1.0,alpa_hs_c: -1.0})
model.add_reaction(reaction)

"""
Modified with respect to recon1 
reaction = Reaction('agpat1')
reaction.name = ' 1-acylglycerol-3-phosphate O-acyltransferase 1'
reaction.subsystem = 'Triacylglycerol Synthesis'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({pa_hs_c: 1.0,rtotal2coa_c: -1.0,coa_c: 1.0,alpa_hs_c: -1.0})
model.add_reaction(reaction)
"""

reaction = Reaction('ppap')
reaction.name = 'phosphatidic acid phosphatase'
reaction.subsystem = 'Triacylglycerol Synthesis'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({pa_hs_c: -1.0,h2o_c: -1.0,pi_c: 1.0,dag_hs_c: 1.0})
model.add_reaction(reaction)

reaction = Reaction('g3pd1')
reaction.name = 'glycerol-3-phosphate dehydrogenase (NAD)'
reaction.subsystem = 'Glycerophospholipid Metabolism'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({glyc3p_c: -1.0,nadh_c: 1.0,dhap_c: 1.0,h_c: 1.0,nad_c: -1.0})
model.add_reaction(reaction)

"""reaction = Reaction('rtot_2')
reaction.name = 'R total flux 2 position'
reaction.subsystem = 'R Group Synthesis'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({'rtotal2coa_c': 2.0, 'r2coa_hs_c': -1.0, 'r4coa_hs_c': -1.0})
model.add_reaction(reaction)
Removed from Recon1
"""
reaction = Reaction('desat18_3')
reaction.name = 'stearoyl-CoA desaturase (n-C18:0CoA -> n-C18:1CoA)'
reaction.subsystem = 'Fatty acid elongation'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({o2_c: -1.0, h2o_c: 2.0, odecoa_c: 1.0, stcoa_c: -1.0, h_c: -1.0, nad_c: 1.0, nadh_c: -1.0})
model.add_reaction(reaction)

reaction = Reaction('fas180coa')
reaction.name = 'fatty-acyl-CoA synthase (n-C18:0CoA)'
reaction.subsystem = 'Fatty acid elongation'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({h2o_c: 1.0, malcoa_c: -1.0, nadph_c: -2.0, pmtcoa_c: -1.0, stcoa_c: 1.0, nadp_c: 2.0, coa_c: 1.0, h_c: -3.0, co2_c: 1.0})
model.add_reaction(reaction)

reaction = Reaction('desat16_2')
reaction.name = 'palmitoyl-CoA desaturase (n-C16:0CoA -> n-C16:1CoA)'
reaction.subsystem = 'Fatty acid elongation'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({o2_c: -1.0, h2o_c: 2.0, hdcoa_c: 1.0, pmtcoa_c: -1.0, h_c: -1.0, nad_c: 1.0, nadh_c: -1.0})
model.add_reaction(reaction)

reaction = Reaction('ocdceatr')
reaction.name = 'fatty acid transport via diffusion'
reaction.subsystem = 'Transport, Extracellular'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({ocdcea_c: 1.0, ocdcea_e: -1.0})
model.add_reaction(reaction)


#pe_hs_c

reaction = Reaction('pssa1_hs')
reaction.name = 'Phosphatidylserine synthase homo sapiens'
reaction.subsystem = 'Glycerophospholipid Metabolism'
reaction.lower_bound=0 #Made irrevesible for simplicity
reaction.upper_bound=1000.0
reaction.add_metabolites({chol_c: 1.0,pchol_hs_c: -1.0,ser_c: -1.0,ps_hs_c: 1.0})
model.add_reaction(reaction)

reaction = Reaction('psflipm')
reaction.name = 'phosphatidylserine flippase'
reaction.subsystem = 'Transport, Mitochondrial'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({h2o_c: -1.0,atp_c: -1.0,ps_hs_c: -1.0,ps_hs_m: 1.0,h_c: 1.0,adp_c: 1.0,pi_c: 1.0})
model.add_reaction(reaction)

reaction = Reaction('psdm_hs')
reaction.name = 'Phosphatidylserine decarboxylase'
reaction.subsystem = 'Glycerophospholipid Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({pe_hs_m: 1.0,ps_hs_m: -1.0,h_m: -1.0,co2_m: 1.0})
model.add_reaction(reaction)

reaction = Reaction('pe_hstm')
reaction.name = 'phosphatidylethanolamine scramblase'
reaction.subsystem = 'Transport, Mitochondrial'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({pe_hs_m: 1.0,pe_hs_c: -1.0})
model.add_reaction(reaction)

##Used to test the model, remove in the full version
reaction = Reaction('peflip')
reaction.name = 'phosphatidylethanolamine flippase'
reaction.subsystem = 'Transport, Extracellular'
reaction.lower_bound=-1000 #Switched
reaction.upper_bound=1000.0
reaction.add_metabolites({h2o_c: -1.0, atp_c: -1.0, h_c: 1.0, adp_c: 1.0, pe_hs_e: -1.0, pi_c: 1.0, pe_hs_c: 1.0})
model.add_reaction(reaction)


##################Cholsterol 

# Lipid metabolism Hypoxia http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3820259/
#Buscar fonts de NADPH I NADH al reticle, atp ppi al peroxisome
"""for x in model.reactions:
    if x.subsystem=="Cholesterol Metabolism":
       print (x.id+" "+str(fva[x.id]))"""

reaction = Reaction('acact1r')
reaction.name = 'acetyl-CoA C-acetyltransferase'
reaction.subsystem = 'Tryptophan metabolism'
reaction.lower_bound=0 #Made irreversible for simplicity
reaction.upper_bound=1000.0
reaction.add_metabolites({aacoa_c: 1.0, coa_c: 1.0, accoa_c: -2.0})
model.add_reaction(reaction)

reaction = Reaction('hmgcoasi')
reaction.name = 'Hydroxymethylglutaryl CoA synthase (ir)'
reaction.subsystem = 'Cholesterol Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({h2o_c: -1.0, aacoa_c: -1.0, hmgcoa_c: 1.0, coa_c: 1.0, h_c: 1.0, accoa_c: -1.0})
model.add_reaction(reaction)

reaction = Reaction('hmgcoatx')
reaction.name = 'Hydroxymethylglutaryl-CoA reversible peroxisomal transport'
reaction.subsystem = 'Transport, Peroxisomal'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({hmgcoa_c: -1.0,hmgcoa_x: 1.0})
model.add_reaction(reaction)

reaction = Reaction('hmgcoarx')
reaction.name = 'Hydroxymethylglutaryl CoA reductase (ir)'
reaction.subsystem = 'Cholesterol Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({hmgcoa_x: -1.0,h_x: -2.0,coa_x: 1.0,nadp_x: 2.0,mev_r_x: 1.0,nadph_x: -2.0})
model.add_reaction(reaction)

reaction = Reaction('mevk1x')
reaction.name = 'mevalonate kinase (atp)'
reaction.subsystem = 'Cholesterol Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({mev_r_x: -1.0,n5pmev_x: 1.0,adp_x: 1.0,h_x: 1.0,atp_x: -1.0})
model.add_reaction(reaction)

reaction = Reaction('pmevkx')
reaction.name = 'phosphomevalonate kinase'
reaction.subsystem = 'Cholesterol Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({n5dpmev_x: 1.0,n5pmev_x: -1.0,atp_x: -1.0,adp_x: 1.0})
model.add_reaction(reaction)

reaction = Reaction('dpmvdx')
reaction.name = 'diphosphomevalonate decarboxylase'
reaction.subsystem = 'Cholesterol Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({pi_x: 1.0,co2_x: 1.0,adp_x: 1.0,ipdp_x: 1.0,n5dpmev_x: -1.0,atp_x: -1.0})
model.add_reaction(reaction)

reaction = Reaction('ipddix')
reaction.name = 'isopentenyl-diphosphate D-isomerase'
reaction.subsystem = 'Cholesterol Metabolism'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({ipdp_x: -1.0,dmpp_x: 1.0})
model.add_reaction(reaction)

reaction = Reaction('dmattx')
reaction.name = 'dimethylallyltranstransferase'
reaction.subsystem = 'Cholesterol Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({ipdp_x: -1.0,ppi_x: 1.0,dmpp_x: -1.0,grdp_x: 1.0})
model.add_reaction(reaction)

reaction = Reaction('grttx')
reaction.name = 'geranyltranstransferase'
reaction.subsystem = 'Cholesterol Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({ipdp_x: -1.0,ppi_x: 1.0,grdp_x: -1.0,frdp_x: 1.0})
model.add_reaction(reaction)



reaction = Reaction('frdptr')
reaction.name = 'lipid, flip-flop intracellular transport'
reaction.subsystem = 'Transport, Peroxisomal'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({frdp_r: 1.0,frdp_x: -1.0})
model.add_reaction(reaction)

reaction = Reaction('sqlsr')
reaction.name = 'Squalene synthase'
reaction.subsystem = 'Cholesterol Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({h_r: -1.0,nadp_r: 1.0,frdp_r: -2.0,ppi_r: 2.0,sql_r: 1.0,nadph_r: -1.0})
model.add_reaction(reaction)

reaction = Reaction('sqler')
reaction.name = 'Squalene epoxidase, endoplasmic reticular (NADP)'
reaction.subsystem = 'Cholesterol Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({h_r: -1.0,ssq23epx_r: 1.0,o2_r: -1.0,nadp_r: 1.0,h2o_r: 1.0,nadph_r: -1.0,sql_r: -1.0})
model.add_reaction(reaction)



reaction = Reaction('lnstlsr')
reaction.name = 'lanosterol synthase'
reaction.subsystem = 'Cholesterol Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({lanost_r: 1.0,ssq23epx_r: -1.0})
model.add_reaction(reaction)

reaction = Reaction('lns14dmr')
reaction.name = 'cytochrome P450 lanosterol 14-alpha-demethylase'
reaction.subsystem = 'Cholesterol Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({n44mctr_r: 1.0,h_r: -2.0,for_r: 1.0,o2_r: -3.0,nadp_r: 3.0,h2o_r: 4.0,nadph_r: -3.0,lanost_r: -1.0})
model.add_reaction(reaction)

reaction = Reaction('c14strr')
reaction.name = 'C-14 sterol reductase'
reaction.subsystem = 'Cholesterol Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({h_r: -1.0,n44mctr_r: -1.0,nadp_r: 1.0,n44mzym_r: 1.0,nadph_r: -1.0})
model.add_reaction(reaction)

reaction = Reaction('c4stmo1r')
reaction.name = 'C-4 sterol methyl oxidase (4,4-dimethylzymosterol)'
reaction.subsystem = 'Cholesterol Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({h_r: -3.0,n44mzym_r: -1.0,o2_r: -3.0,nadp_r: 3.0,h2o_r: 4.0,n4mzym_int1_r: 1.0,nadph_r: -3.0})
model.add_reaction(reaction)

reaction = Reaction('c3stdh1r')
reaction.name = 'C-3 sterol dehydrogenase (4-methylzymosterol)'
reaction.subsystem = 'Cholesterol Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({h_r: 1.0,nad_r: -1.0,nadh_r: 1.0,n4mzym_int2_r: 1.0,n4mzym_int1_r: -1.0,co2_r: 1.0})
model.add_reaction(reaction)

reaction = Reaction('c4stmo2r')
reaction.name = 'C-4 methyl sterol oxidase'
reaction.subsystem = 'Cholesterol Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({h_r: 1.0,nad_r: -1.0,co2_r: 1.0,n4mzym_int2_r: -1.0,o2_r: -1.0,zym_int2_r: 1.0,nadh_r: 1.0})
model.add_reaction(reaction)

reaction = Reaction('c3stkr2r')
reaction.name = 'C-3 sterol keto reductase (zymosterol)'
reaction.subsystem = 'Cholesterol Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({zym_int2_r: -1.0,nadp_r: 1.0,h_r: -1.0,nadph_r: -1.0,zymst_r: 1.0})
model.add_reaction(reaction)

reaction = Reaction('dhcr241r')
reaction.name = ' 24-dehydrocholesterol reductase [Precursor]'
reaction.subsystem = 'Cholesterol Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({nadh_r: -1.0,nad_r: 1.0,h_r: -1.0,zymst_r: -1.0,zymstnl_r: 1.0})
model.add_reaction(reaction)

"""
Modified from recon1, original reaction is as follows: 
reaction = Reaction('dhcr241r')
reaction.name = ' 24-dehydrocholesterol reductase [Precursor]'
reaction.subsystem = 'Cholesterol Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({fadh2_r: -1.0,fad_r: 1.0,zymst_r: -1.0,zymstnl_r: 1.0})
model.add_reaction(reaction)

"""


reaction = Reaction('ebp2r')
reaction.name = ' 3-beta-hydroxysteroid-delta(8),delta(7)-isomerase'
reaction.subsystem = 'Cholesterol Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({lthstrl_r: 1.0,zymstnl_r: -1.0})
model.add_reaction(reaction)

reaction = Reaction('lsto2r')
reaction.name = 'Lathosterol oxidase'
reaction.subsystem = 'Cholesterol Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({h_r: -1.0,n7dhchsterol_r: 1.0,o2_r: -1.0,lthstrl_r: -1.0,nadp_r: 1.0,h2o_r: 2.0,nadph_r: -1.0})
model.add_reaction(reaction)

reaction = Reaction('dhcr72r')
reaction.name = ' 7-dehydrocholesterol reductase'
reaction.subsystem = 'Cholesterol Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({n7dhchsterol_r: -1.0,nadp_r: 1.0,h_r: -1.0,chsterol_r: 1.0,nadph_r: -1.0})
model.add_reaction(reaction)

reaction = Reaction('chsterolt2')
reaction.name = 'cholesterol intracellular transport'
reaction.subsystem = 'Transport, Mitochondrial'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({chsterol_m: 1.0,chsterol_r: -1.0})
model.add_reaction(reaction)

reaction = Reaction('chsterolt3')
reaction.name = 'cholesterol intracellular transport'
reaction.subsystem = 'Transport, Mitochondrial'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({chsterol_m: 1.0,chsterol_c: -1.0})
model.add_reaction(reaction)

#Cholesterol input from medium
#Taken from recon2 as it is not present in recon1
reaction = Reaction('chsterol_etr')
reaction.name = 'cholesterol vesicular transport'
reaction.subsystem = 'Transport'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({chsterol_e: 1.0,chsterol_c: -1.0})
model.add_reaction(reaction)




reaction = Reaction('ppitr')
reaction.name = 'Diphosphate transporter, endoplasmic reticulum'
reaction.subsystem = 'Transport, Endoplasmic Reticular'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({ppi_c: -1.0,ppi_r: 1.0})
model.add_reaction(reaction)

reaction = Reaction('nadphtru')
reaction.name = 'NADPH transporter, endoplasmic reticulum'
reaction.subsystem = 'Transport, Endoplasmic Reticular'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({nadph_c: -1.0,nadph_r: 1.0})
model.add_reaction(reaction)

reaction = Reaction('nadptru')
reaction.name = 'NADP transporter, endoplasmic reticulum'
reaction.subsystem = 'Transport, Endoplasmic Reticular'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({nadp_r: -1.0,nadp_c: 1.0})
model.add_reaction(reaction)

reaction = Reaction('nadhtru')
reaction.name = 'NADH transporter, endoplasmic reticulum'
reaction.subsystem = 'Transport, Endoplasmic Reticular'
reaction.lower_bound=-1000 #Made reversible compared to recon1
reaction.upper_bound=1000.0
reaction.add_metabolites({nadh_r: 1.0,nadh_c: -1.0})
model.add_reaction(reaction)

reaction = Reaction('nadtru')
reaction.name = 'NAD transporter, endoplasmic reticulum'
reaction.subsystem = 'Transport, Endoplasmic Reticular'
reaction.lower_bound=-1000 #Made reversible compared to recon1
reaction.upper_bound=1000.0
reaction.add_metabolites({nad_c: 1.0,nad_r: -1.0})
model.add_reaction(reaction)


reaction = Reaction('coatp')
reaction.name = 'coenzyme A transport, peroxisomal'
reaction.subsystem = 'Transport, Peroxisomal'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({coa_c: -1.0, coa_x: 1.0})
model.add_reaction(reaction)

reaction = Reaction('o2ter')
reaction.name = 'O2  transport, endoplasmic reticulum'
reaction.subsystem = 'Transport, Endoplasmic Reticular'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({o2_r: 1.0, o2_c: -1.0})
model.add_reaction(reaction)

reaction = Reaction('htx')
reaction.name = 'H transporter, peroxisome'
reaction.subsystem = 'Transport, Peroxisomal'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({h_c: -1.0, h_x: 1.0})
model.add_reaction(reaction)

reaction = Reaction('nadphtxu')
reaction.name = 'NADPH transporter, peroxisome'
reaction.subsystem = 'Transport, Peroxisomal'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({nadph_c: -1.0, nadph_x: 1.0})
model.add_reaction(reaction)

reaction = Reaction('nadptxu')
reaction.name = 'NADP transporter, peroxisome'
reaction.subsystem = 'Transport, Peroxisomal'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({nadp_x: -1.0, nadp_c: 1.0})
model.add_reaction(reaction)

reaction = Reaction('atptx')
reaction.name = 'ATP transporter, peroxisomal'
reaction.subsystem = 'Transport, Peroxisomal'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({atp_x: 1.0, atp_c: -1.0})
model.add_reaction(reaction)

reaction = Reaction('adptx')
reaction.name = 'ADP transporter, peroxisomal'
reaction.subsystem = 'Transport, Peroxisomal'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({adp_x: 1.0, adp_c: -1.0})
model.add_reaction(reaction)

reaction = Reaction('co2tp')
reaction.name = 'CO2 peroxisomal transport'
reaction.subsystem = 'Transport, Peroxisomal'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({co2_x: 1.0, co2_c: -1.0})
model.add_reaction(reaction)

reaction = Reaction('pitx')
reaction.name = 'Phosphate transporter, peroxisome'
reaction.subsystem = 'Transport, Peroxisomal'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({pi_x: 1.0, pi_c: -1.0})
model.add_reaction(reaction)

reaction = Reaction('ppitx')
reaction.name = 'Diphosphate transporter, peroxisome'
reaction.subsystem = 'Transport, Peroxisomal'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({ppi_x: 1.0, ppi_c: -1.0})
model.add_reaction(reaction)

reaction = Reaction('htr')
reaction.name = 'H transporter, endoplasmic reticulum'
reaction.subsystem = 'Transport, Endoplasmic Reticular'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({h_r: 1.0, h_c: -1.0})
model.add_reaction(reaction)

reaction = Reaction('h2oter')
reaction.name = 'H2O endoplasmic reticulum transport'
reaction.subsystem = 'Transport, Endoplasmic Reticular'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({h2o_c: -1.0, h2o_r: 1.0})
model.add_reaction(reaction)

reaction = Reaction('fortr')
reaction.name = 'FOR transporter, endoplasmic reticulum'
reaction.subsystem = 'Transport, Endoplasmic Reticular'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({for_c: -1.0, for_r: 1.0})
model.add_reaction(reaction)

reaction = Reaction('co2ter')
reaction.name = 'CO2 endoplasmic reticular transport via diffusion'
reaction.subsystem = 'Transport, Endoplasmic Reticular'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({co2_r: 1.0, co2_c: -1.0})
model.add_reaction(reaction)

################Inositol Biosynthesis

reaction = Reaction('cds')
reaction.name = 'phosphatidate cytidylyltransferase'
reaction.subsystem = 'Glycerophospholipid Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({ctp_c: -1.0, pa_hs_c: -1.0, ppi_c: 1.0, h_c: -1.0, cdpdag_hs_c: 1.0})
model.add_reaction(reaction)

reaction = Reaction('mi1ps')
reaction.name = 'myo-Inositol-1-phosphate synthase'
reaction.subsystem = 'Inositol Phosphate Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({glc6p_b_c: -1.0, mi1p_d_c: 1.0})
model.add_reaction(reaction)

reaction = Reaction('mi1pp')
reaction.name = 'myo-inositol 1-phosphatase'
reaction.subsystem = 'Inositol Phosphate Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({h2o_c: -1.0, pi_c: 1.0, inost_c: 1.0, mi1p_d_c: -1.0})
model.add_reaction(reaction)

reaction = Reaction('cdiptr')
reaction.name = 'phosphatidylinositol synthase (Homo sapiens)'
reaction.subsystem = 'Glycerophospholipid Metabolism'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({pail_hs_c: 1.0, inost_c: -1.0, h_c: 1.0, cdpdag_hs_c: -1.0, cmp_c: 1.0})
model.add_reaction(reaction)

####### clpn_hs_c and pglyc_hs_c sysnthesis
reaction = Reaction('pgppt')
reaction.name = 'phosphatidyl-CMP: glycerophosphate phosphatidyltransferase'
reaction.subsystem = 'Glycerophospholipid Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({pgp_hs_c: 1.0, glyc3p_c: -1.0, h_c: 1.0, cdpdag_hs_c: -1.0, cmp_c: 1.0})
model.add_reaction(reaction)

reaction = Reaction('pgpp_hs')
reaction.name = 'Phosphatidylglycerol phosphate phosphatase (homo sapiens)'
reaction.subsystem = 'Glycerophospholipid Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({pgp_hs_c: -1.0, h2o_c: -1.0, pglyc_hs_c: 1.0, pi_c: 1.0})
model.add_reaction(reaction)

reaction = Reaction('cls_hs')
reaction.name = 'cardiolipin synthase (homo sapiens)'
reaction.subsystem = 'Glycerophospholipid Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({clpn_hs_c: 1.0, pglyc_hs_c: -1.0, h_c: 1.0, cdpdag_hs_c: -1.0, cmp_c: 1.0})
model.add_reaction(reaction)

########### Sphingomyelin synthesis
reaction = Reaction('serpt')
reaction.name = 'serine C-palmitoyltransferase'
reaction.subsystem = 'Sphingolipid Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({n3dsphgn_c: 1.0, pmtcoa_c: -1.0, ser_c: -1.0, coa_c: 1.0, h_c: -1.0, co2_c: 1.0})
model.add_reaction(reaction)

reaction = Reaction('3dsphr')
reaction.name = ' 3-Dehydrosphinganine reductase'
reaction.subsystem = 'Sphingolipid Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({n3dsphgn_c: -1.0, nadp_c: 1.0, nadph_c: -1.0, h_c: -1.0, sphgn_c: 1.0})
model.add_reaction(reaction)

reaction = Reaction('dsat')
reaction.name = 'dihydrosphingosine N-acyltransferase'
reaction.subsystem = 'Sphingolipid Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({dhcrm_hs_c: 1.0, rtotalcoa_c: -1.0, coa_c: 1.0, h_c: 1.0, sphgn_c: -1.0})
model.add_reaction(reaction)


reaction = Reaction('dhcrd2')
reaction.name = 'dihydroceramide desaturase'
reaction.subsystem = 'Sphingolipid Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({dhcrm_hs_c: -1.0,fadh2_m: 1.0, crm_hs_c: 1.0, fad_m: -1.0,h_c: 1.0}) #Original reaction used cytosolic
model.add_reaction(reaction)
"""
modified with respect of recon1, original reaction is
reaction = Reaction('dhcrd2')
reaction.name = 'dihydroceramide desaturase'
reaction.subsystem = 'Sphingolipid Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({dhcrm_hs_c: -1.0, fadh2_c: 1.0, crm_hs_c: 1.0, fad_c: -1.0})
model.add_reaction(reaction)"""

reaction = Reaction('sms')
reaction.name = 'Sphingomyelin synthase (homo sapiens)'
reaction.subsystem = 'Sphingolipid Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({pchol_hs_c: -1.0, dag_hs_c: 1.0, sphmyln_hs_c: 1.0, crm_hs_c: -1.0})
model.add_reaction(reaction)


reaction = Reaction('biomass')
reaction.name = 'Biomass'
reaction.subsystem = ''
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({utp_c: -0.053446, lys_c: -0.59211, thr_c: -0.31269, pro_c: -0.41248, h2o_c: -20.6508, ser_c: -0.39253, h_c: 20.6508, ctp_c: -0.039036, pail_hs_c: -0.023315, val_c: -0.35261, trp_c: -0.013306, ala_c: -0.50563, pglyc_hs_c: -0.002914, chsterol_c: -0.020401, atp_c: -20.7045, asp_c: -0.35261, leu_c: -0.54554, gln_c: -0.326, dgtp_c: -0.009898, pe_hs_c: -0.055374, adp_c: 20.6508, tyr_c: -0.15967, arg_c: -0.35926, glc6p_b_c: -0.27519, ile_c: -0.28608, gly_c: -0.53889, cys_c: -0.046571, met_c: -0.15302, phe_c: -0.25947, dttp_c: -0.013091, asn_c: -0.27942, his_c: -0.12641, glu_c: -0.38587, pchol_hs_c: -0.15446, sphmyln_hs_c: -0.017486, ps_hs_c: -0.005829, clpn_hs_c: -0.011658, dctp_c: -0.009442, gtp_c: -0.036117, datp_c: -0.013183, pi_c: 20.6508}) #gDW of Biomass produced
model.add_reaction(reaction)

#warning artifical reaction, remove
reaction = Reaction('spms_dm')
reaction.name = 'Spermine demmand'
reaction.subsystem = ''
reaction.lower_bound=0
reaction.upper_bound=1000.0
reaction.add_metabolites({sprm_c:-1})
model.add_reaction(reaction)

reaction = Reaction('spmd_dm')
reaction.name = 'Spermidine demmand'
reaction.subsystem = ''
reaction.lower_bound=0
reaction.upper_bound=1000.0
reaction.add_metabolites({spmd_c:-1})
model.add_reaction(reaction)

"""dump_c
reaction = Reaction('dump_c_entry')
reaction.name = 'dump_c_entry'
reaction.subsystem = 'REMOVE ME'
reaction.lower_bound=0
reaction.upper_bound=1000.0
reaction.add_metabolites({dump_c:1})
model.add_reaction(reaction)

reaction = Reaction('atp_synthesis')
reaction.name = 'atp_synthesis'
reaction.subsystem = 'REMOVE ME'
reaction.lower_bound=0
reaction.upper_bound=1000.0
reaction.add_metabolites({rib5p_c:-1,atp_novo_c:1,atp_c:-2,adp_c:2,h_c:2})
model.add_reaction(reaction)"""

reaction = Reaction('adnk1')
reaction.name = 'adenosine kinase'
reaction.subsystem = 'Nucleotides'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({atp_c: -1.0, adn_c: -1.0, amp_c: 1.0, adp_c: 1.0, h_c: 1.0})
model.add_reaction(reaction)

"""reaction = Reaction('atp_exchange')
reaction.name = 'atp_exchange'
reaction.subsystem = 'REMOVE ME'
reaction.lower_bound=-1000
reaction.upper_bound=1000.0
reaction.add_metabolites({atp_novo_c:-1,atp_c:1})
model.add_reaction(reaction)"""

"""reaction = Reaction('atp_removal')
reaction.name = 'atp_removal'
reaction.subsystem = 'REMOVE ME'
reaction.lower_bound=0
reaction.upper_bound=1000.0
reaction.add_metabolites({atp_novo_c:-1})
model.add_reaction(reaction)"""

#warning artifical reaction, remove

#Biomass synthesis 
"""
reaction = Reaction('prot_syn')
reaction.name = 'Protein synthesis'
reaction.subsystem = 'AA'
reaction.lower_bound = 0
reaction.upper_bound = 1000.  
reaction.add_metabolites({ala_c:-0.088,  arg_c:-0.055, asn_c:-0.042, asp_c:-0.053,gln_c:-0.047,glu_c:-0.057,gly_c:-0.079,  his_c:-0.021,pro_c:-0.046, ser_c: -0.063,  thr_c:-0.057,  trp_c:-0.006,  cys_c:-0.021,  lys_c:-0.084,   ile_c:-0.048,  leu_c:-0.083,  tyr_c:-0.027,  phe_c:-0.032,  met_c:-0.020,  val_c:-0.061, atp_c:-4.3,adp_c:4.3,pi_c:4.3,h_c: 4.3,prot_c:1}) #Reaction produces 1mmol of protein
        
model.add_reaction(reaction)
"""


"""reaction = Reaction('prot_deg')
reaction.name = 'Protein degradation'
reaction.subsystem = 'AA'
reaction.lower_bound = 0
reaction.upper_bound = 1000.  
reaction.add_metabolites({ala_c:0.088,  arg_c:0.055, asn_c:0.042, asp_c:0.053,gln_c:0.047,glu_c:0.057,gly_c:0.079,  his_c:0.021,pro_c:0.046, ser_c:0.063,  thr_c:0.057,  trp_c:0.006,  cys_c:0.021,  lys_c:0.084,   ile_c:0.048,  leu_c:0.083,  tyr_c:0.027,  phe_c:0.032,  met_c:0.020,  val_c:0.061, prot_c:-1})        
model.add_reaction(reaction) #Reaction degrades 1mmol of protein

reaction = Reaction('dna_syn')
reaction.name = 'DNA_syn'
reaction.subsystem = 'REMOVE ME'
reaction.lower_bound=0
reaction.upper_bound=1000.0
reaction.add_metabolites({dctp_c:-0.2,dgtp_c:-0.2,datp_c:-0.3,dttp_c:-0.3,atp_c:-1.372,adp_c:1.372,pi_c:3.372,h_c: 1.372,dna_c:1}) #1 mM of DNA
model.add_reaction(reaction)


reaction = Reaction('rna_syn')
reaction.name = 'RNA_syn'
reaction.subsystem = 'REMOVE ME'
reaction.lower_bound=0
reaction.upper_bound=1000.0
reaction.add_metabolites({ctp_c:-0.3,gtp_c:-0.34,atp_c:-0.58,utp_c:-0.18,adp_c:0.4,pi_c:2.4,h_c:0.4,rna_c:1})
model.add_reaction(reaction)"""

#Urea cycle and Arginine synthesis 
reaction = Reaction('cbpsam')
reaction.name = 'carbamoyl-phosphate synthase (ammonia) (mitochondria'
reaction.subsystem = 'Glutamate metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({hco3_m: -1.0, atp_m: -2.0, cbp_m: 1.0, nh4_m: -1.0, pi_m: 1.0, adp_m: 2.0, h_m: 2.0})
model.add_reaction(reaction)

reaction = Reaction('ocbtm')
reaction.name = 'ornithine carbamoyltransferase, irreversible'
reaction.subsystem = 'Urea cycle/amino group metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({orn_m: -1.0, h_m: 1.0, cbp_m: -1.0, pi_m: 1.0, citr_m: 1.0})
model.add_reaction(reaction)

reaction = Reaction('ornt4m')
reaction.name = 'ornithine mitochondrial transport exchange with citruline'
reaction.subsystem = 'Transport, Mitochondrial'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({citr_m: 1.0, citr_c: -1.0, orn_c: 1.0, h_c: -1.0, orn_m: -1.0, h_m: 1.0})
model.add_reaction(reaction)

reaction = Reaction('argss')
reaction.name = 'argininosuccinate synthase'
reaction.subsystem = 'Alanine and Aspartate Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({atp_c: -1.0, ppi_c: 1.0, citr_c: -1.0, argsuc_c: 1.0, amp_c: 1.0, asp_c: -1.0, h_c: 1.0})
model.add_reaction(reaction)

reaction = Reaction('argsl')
reaction.name = 'argininosuccinate lyase'
reaction.subsystem = 'Alanine and Aspartate Metabolism'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({arg_c: 1.0, argsuc_c: -1.0, fum_c: 1.0})
model.add_reaction(reaction)

#Asparagine Degradation 

reaction = Reaction('asntm')
reaction.name = 'L-asparagine transport, mitochondrial'
reaction.subsystem = 'Transport, Mitochondrial'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({asn_c: -1.0, asn_m: 1.0})
model.add_reaction(reaction)

reaction = Reaction('asnnm')
reaction.name = 'L-asparaginase (mitochondrial)'
reaction.subsystem = 'Alanine and Aspartate Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({asn_m: -1.0, asp_m: 1.0, nh4_m: 1.0, h2o_m: -1.0})
model.add_reaction(reaction)

#To allow the degradation of glycine through glycine
methf_m = Metabolite('methf_m',formula='C20H20N7O6',name='5,10-Methenyltetrahydrofolate',compartment='m')
n10fthf_m = Metabolite('10fthf_m',formula='C20H21N7O7',name='10-Formyltetrahydrofolate',compartment='m')
thf_m = Metabolite('thf_m',formula='C19H21N7O6',name='5,6,7,8-Tetrahydrofolate',compartment='m')
lpam_m = Metabolite('lpam_m',formula='C8H15NOS2',name='Lipoamide',compartment='m')
dhlam_m = Metabolite('dhlam_m',formula='C8H17NOS2',name='Dihydrolipoamide',compartment='m')
alpam_m = Metabolite('alpam_m',formula='C9H21N2OS2',name='S-aminomethyldihydrolipoamide',compartment='m')
mlthf_m = Metabolite('mlthf_m',formula='C20H21N7O6',name='5,10-Methylenetetrahydrofolate',compartment='m')

reaction = Reaction('gcc2am')
reaction.name = 'glycine-cleavage complex (lipoamide), mitochondrial'
reaction.subsystem = 'Glycine, Serine, and Threonine Metabolism'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({lpam_m: -1.0, gly_m: -1.0, h_m: -1.0, alpam_m: 1.0, co2_m: 1.0})
model.add_reaction(reaction)

reaction = Reaction('gcc2bim')
reaction.name = 'glycine-cleavage system (lipoamide) irreversible, mitochondrial'
reaction.subsystem = 'Glycine, Serine, and Threonine Metabolism'
reaction.lower_bound=0.0
reaction.upper_bound=1000.0
reaction.add_metabolites({mlthf_m: 1.0, thf_m: -1.0, nh4_m: 1.0, alpam_m: -1.0, dhlam_m: 1.0})
model.add_reaction(reaction)

reaction = Reaction('gcc2cm')
reaction.name = 'glycine-cleavage complex (lipoamide), mitochondrial'
reaction.subsystem = 'Glycine, Serine, and Threonine Metabolism'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({nad_m: -1.0, lpam_m: 1.0, h_m: 1.0, nadh_m: 1.0, dhlam_m: -1.0})
model.add_reaction(reaction)

reaction = Reaction('thftm')
reaction.name = ' 5,6,7,8-Tetrahydrofolate transport, diffusion, mitochondrial'
reaction.subsystem = 'Transport, Mitochondrial'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({thf_m: 1.0, thf_c: -1.0})
model.add_reaction(reaction)


reaction = Reaction('mthfd2m')
reaction.name = 'methylenetetrahydrofolate dehydrogenase (NAD), mitochondrial'
reaction.subsystem = 'Folate Metabolism'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({nad_m: -1.0, mlthf_m: -1.0, methf_m: 1.0, nadh_m: 1.0})
model.add_reaction(reaction)

reaction = Reaction('mthfcm')
reaction.name = 'methenyltetrahydrifikate cyclohydrolase, mitochondrial'
reaction.subsystem = 'Folate Metabolism'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({methf_m: -1.0, h_m: 1.0, h2o_m: -1.0, n10fthf_m: 1.0})
model.add_reaction(reaction)

reaction = Reaction('10fthftm')
reaction.name = ' 10-Formyltetrahydrofolate mitochondrial transport via diffusion'
reaction.subsystem = 'Transport, Mitochondrial'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({n10fthf_m: 1.0, n10fthf_c: -1.0})
model.add_reaction(reaction)

reaction = Reaction('glytm')
reaction.name = 'glycine passive transport to mitochondria'
reaction.subsystem = 'Transport, Mitochondrial'
reaction.lower_bound=-1000.0
reaction.upper_bound=1000.0
reaction.add_metabolites({gly_m: 1.0, gly_c: -1.0})
model.add_reaction(reaction)
########





#TOFIX
#TDODO
"""
Argine cannot be produced
Asparragine cannot be consumed (other than protein synthesis) 
Glcyine cannot be consumed other than protein synthesis
"""

media_metabolites=["h_e","pi_e","o2_e","co2_e","chol_e","h2o_e"]
for x in model.metabolites:
    if x.compartment=="e":
       reaction = Reaction('EX_'+x.id)
       reaction.name = 'Exchange of '+x.name
       reaction.subsystem = 'Exchange reaction'
       if x.id in media_metabolites:
          reaction.lower_bound=-1000
       else:
          reaction.lower_bound=0
       reaction.upper_bound=1000.0
       reaction.add_metabolites({x:-1})
       model.add_reaction(reaction)


#Remove in full version
#Remove in full version
#model.reactions.get_by_id("dna_syn").objective_coefficient=0

#Remove in full version
#Remove in full version
cobra.io.write_sbml_model(model, "test_model.sbml")


model.reactions.get_by_id("EX_glc_e").lower_bound=-1
model.reactions.get_by_id("EX_lac_e").lower_bound=0.5
model.reactions.get_by_id("biomass").objective_coefficient=1
model.reactions.get_by_id("admdc").lower_bound=0.0001





#cobra.io.read_sbml_model("Recon1Cobra.v01.xml")
model.reactions.get_by_id("EX_glc_e").lower_bound=-1
model.reactions.get_by_id("EX_ala_e").lower_bound=-0.1
model.reactions.get_by_id("EX_arg_e").lower_bound=-0.1
model.reactions.get_by_id("EX_asn_e").lower_bound=-0.1
model.reactions.get_by_id("EX_asp_e").lower_bound=-0.1
model.reactions.get_by_id("EX_cys_e").lower_bound=-0.1
model.reactions.get_by_id("EX_gln_e").lower_bound=-0.1
model.reactions.get_by_id("EX_glu_e").lower_bound=-0.1
model.reactions.get_by_id("EX_gly_e").lower_bound=-0.1
model.reactions.get_by_id("EX_his_e").lower_bound=-0.1
model.reactions.get_by_id("EX_ile_e").lower_bound=-0.1
model.reactions.get_by_id("EX_leu_e").lower_bound=-0.1
model.reactions.get_by_id("EX_lys_e").lower_bound=-0.1
model.reactions.get_by_id("EX_met_e").lower_bound=-0.1
model.reactions.get_by_id("EX_phe_e").lower_bound=-0.1
model.reactions.get_by_id("EX_pro_e").lower_bound=-0.1
model.reactions.get_by_id("EX_ser_e").lower_bound=-0.1
model.reactions.get_by_id("EX_thr_e").lower_bound=-0.1
model.reactions.get_by_id("EX_trp_e").lower_bound=-0.1
model.reactions.get_by_id("EX_tyr_e").lower_bound=-0.1
model.reactions.get_by_id("EX_val_e").lower_bound=-0.1
model.reactions.get_by_id("EX_hdca_e").lower_bound=-0.1
model.reactions.get_by_id("EX_chol_e").lower_bound=-0.1
model.reactions.get_by_id("EX_hdca_e").lower_bound=-0.1
model.reactions.get_by_id("EX_ocdcea_e").lower_bound=-0.1


fva=flux_variability_analysis(model,fraction_of_optimum=0.0)

model.reactions.get_by_id("EX_ala_e").lower_bound=0
model.optimize()
model.reactions.get_by_id("EX_arg_e").lower_bound=0
model.optimize()
model.reactions.get_by_id("EX_asn_e").lower_bound=0
model.optimize()
model.reactions.get_by_id("EX_asp_e").lower_bound=0
model.optimize()
model.reactions.get_by_id("EX_cys_e").lower_bound=0
model.optimize()
model.reactions.get_by_id("EX_gly_e").lower_bound=-0.1
model.optimize()
model.reactions.get_by_id("EX_his_e").lower_bound=0
model.optimize()
model.reactions.get_by_id("EX_ile_e").lower_bound=0
model.optimize()
model.reactions.get_by_id("EX_leu_e").lower_bound=0
model.optimize()
model.reactions.get_by_id("EX_lys_e").lower_bound=0
model.optimize()
model.reactions.get_by_id("EX_met_e").lower_bound=0
model.optimize()
model.reactions.get_by_id("EX_phe_e").lower_bound=-0.1
model.optimize()
model.reactions.get_by_id("EX_pro_e").lower_bound=-0.1
model.optimize()
model.reactions.get_by_id("EX_ser_e").lower_bound=0
model.optimize()
model.reactions.get_by_id("EX_thr_e").lower_bound=0
model.optimize()
model.reactions.get_by_id("EX_trp_e").lower_bound=0
model.optimize()
model.reactions.get_by_id("EX_tyr_e").lower_bound=0
model.optimize()
model.reactions.get_by_id("EX_val_e").lower_bound=0
model.optimize()
model.reactions.get_by_id("EX_hdca_e").lower_bound=-0.1
model.optimize()

fva=flux_variability_analysis(model,fraction_of_optimum=0.0)


model.reactions.get_by_id("gs").lower_bound=0.001
model.reactions.get_by_id("admdc").lower_bound=0.0001
model.reactions.get_by_id("EX_glc_e").lower_bound=-0.282
model.reactions.get_by_id("EX_glc_e").upper_bound=-0.282
model.reactions.get_by_id("EX_gln_e").upper_bound=-0.235
model.reactions.get_by_id("EX_gln_e").lower_bound=-0.235
model.reactions.get_by_id("EX_glu_e").upper_bound=0.061
model.reactions.get_by_id("EX_glu_e").lower_bound=0.061
model.reactions.get_by_id("EX_lac_e").lower_bound=1.558
model.reactions.get_by_id("EX_lac_e").upper_bound=1.558

model.reactions.get_by_id("dna_syn").objective_coefficient=1
#model.reactions.get_by_id("pfk").lower_bound=0.1 #TEMP
#model.reactions.get_by_id("pdhm").lower_bound=0.1 #TEMP
#model.reactions.get_by_id("EX_ile_e").objective_coefficient=1
model.reactions.get_by_id("EX_ala_e").lower_bound=-4
#model.reactions.get_by_id("EX_arg_e").lower_bound=-0.5*factor
#model.reactions.get_by_id("EX_asn_e").lower_bound=-0.5*factor
model.reactions.get_by_id("EX_asp_e").lower_bound=-4
#model.reactions.get_by_id("EX_cys_e").lower_bound=-1*factor
#model.reactions.get_by_id("EX_gln_e").lower_bound=1
#model.reactions.get_by_id("EX_glu_e").lower_bound=1
model.reactions.get_by_id("EX_gly_e").lower_bound=-1
#model.reactions.get_by_id("EX_his_e").lower_bound=-0.5*factor
model.reactions.get_by_id("EX_ile_e").lower_bound=-2
#model.reactions.get_by_id("EX_leu_e").lower_bound=-0.5*factor
#model.reactions.get_by_id("EX_lys_e").lower_bound=-0.5*factor
#model.reactions.get_by_id("EX_met_e").lower_bound=-0.5*factor
#model.reactions.get_by_id("EX_phe_e").lower_bound=-0.5*factor
#model.reactions.get_by_id("EX_pro_e").lower_bound=-0.5*factor
model.reactions.get_by_id("EX_ser_e").lower_bound=-2
#model.reactions.get_by_id("EX_thr_e").lower_bound=-0.5*factor
#model.reactions.get_by_id("EX_trp_e").lower_bound=-0.5*factor
#model.reactions.get_by_id("EX_tyr_e").lower_bound=-0.5*factor
#model.reactions.get_by_id("EX_val_e").lower_bound=-0.5*factor
model.reactions.get_by_id("EX_hdca_e").lower_bound=-1


#TEMP
model.reactions.get_by_id("gluxm").lower_bound=0
model.reactions.get_by_id("icdhyrm").lower_bound=0
model.reactions.get_by_id("me1m").lower_bound=0
model.reactions.get_by_id("me1m").upper_bound=0
model.reactions.get_by_id("me2").lower_bound=0
model.reactions.get_by_id("me2").upper_bound=0
model.reactions.get_by_id("gp").upper_bound=0

model.optimize()

for x in fva:
    if abs(fva[x]["maximum"])<0.001 and abs(fva[x]["minimum"])<0.001:
       print(x+" "+str(fva[x]))

for x in fva:
    if "EX_" in x:
       print(x+" "+str(fva[x]))


for x in model.reactions:
    if x.metabolites=={} or x.metabolites==None:
       print x 



model=cobra.io.read_sbml_model("Recon1Cobra.v01.xml")
model.reactions.get_by_id("biomass").objective_coefficient=1
model.reactions.get_by_id("EX_glc_LPAREN_e_RPAREN_").lower_bound=-1
model.reactions.get_by_id("EX_ala_DASH_L_LPAREN_e_RPAREN_").lower_bound=-0.1
model.reactions.get_by_id("EX_arg_DASH_L_LPAREN_e_RPAREN_").lower_bound=-0.1
model.reactions.get_by_id("EX_asn_DASH_L_LPAREN_e_RPAREN_").lower_bound=-0.1
model.reactions.get_by_id("EX_asp_DASH_L_LPAREN_e_RPAREN_").lower_bound=-0.1
model.reactions.get_by_id("EX_cys_DASH_L_LPAREN_e_RPAREN_").lower_bound=-0.1
model.reactions.get_by_id("EX_gln_DASH_L_LPAREN_e_RPAREN_").lower_bound=-0.1
model.reactions.get_by_id("EX_glu_DASH_L_LPAREN_e_RPAREN_").lower_bound=-0.1
model.reactions.get_by_id("EX_glyc_LPAREN_e_RPAREN_").lower_bound=-0.1
model.reactions.get_by_id("EX_his_DASH_L_LPAREN_e_RPAREN_").lower_bound=-0.1
model.reactions.get_by_id("EX_ile_DASH_L_LPAREN_e_RPAREN_").lower_bound=-0.1
model.reactions.get_by_id("EX_leu_DASH_L_LPAREN_e_RPAREN_").lower_bound=-0.1
model.reactions.get_by_id("EX_lys_DASH_L_LPAREN_e_RPAREN_").lower_bound=-0.1
model.reactions.get_by_id("EX_met_DASH_L_LPAREN_e_RPAREN_").lower_bound=-0.1
model.reactions.get_by_id("EX_phe_DASH_L_LPAREN_e_RPAREN_").lower_bound=-0.1
model.reactions.get_by_id("EX_pro_DASH_L_LPAREN_e_RPAREN_").lower_bound=-0.1
model.reactions.get_by_id("EX_ser_DASH_L_LPAREN_e_RPAREN_").lower_bound=-0.1
model.reactions.get_by_id("EX_thr_DASH_L_LPAREN_e_RPAREN_").lower_bound=-0.1
model.reactions.get_by_id("EX_trp_DASH_L_LPAREN_e_RPAREN_").lower_bound=-0.1
model.reactions.get_by_id("EX_tyr_DASH_L_LPAREN_e_RPAREN_").lower_bound=-0.1
model.reactions.get_by_id("EX_val_DASH_L_LPAREN_e_RPAREN_").lower_bound=-0.1
model.reactions.get_by_id("EX_hdca_LPAREN_e_RPAREN_").lower_bound=-0.1
model.reactions.get_by_id("EX_chol_LPAREN_e_RPAREN_").lower_bound=-0.1



for x in model.reactions.query("EX_"):
    if model.solution.x_dict[x.id]<0:
       print x

fva=flux_variability_analysis(model,fraction_of_optimum=1)
for x in fva:
    if

