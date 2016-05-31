import json
def save_parameters(label_model,name="parameter_dict.json"):
    with open(name, 'w') as fp:
         json.dump(label_model.parameter_dict, fp)


