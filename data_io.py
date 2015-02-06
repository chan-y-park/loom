import logging
import json
import pdb


def save_spectral_network_data(data_list, file_object, **kwargs):
    json_data = {}
    for data in data_list:
        json_data['phase'] = data['phase']
        json_data['hit_table'] = data['hit_table'].get_json_data()
        json_data['s_walls'] = [s_wall.get_json_data()
                                for s_wall in data['s_walls']]
        json_data['joints'] = [joint.get_json_data()
                               for joint in data['joints']]

    json.dump(json_data, file_object, **kwargs)
    

def load_spectral_network_data(file_object, **kwargs):
    data = json.load(file_object, **kwargs)
