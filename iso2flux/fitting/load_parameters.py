import json

def load_parameters(name="parameter_dict.json",label_model=None):
    with open(name, 'r') as fp:
         parameter_dict=json.load(fp)
    if label_model!=None:
       label_model.parameter_dict= parameter_dict
    return parameter_dict

