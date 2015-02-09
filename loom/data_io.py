import logging
import json
import pdb

from spectral_network import SpectralNetwork
from s_wall import Joint, SWall


#def save_spectral_network_data(data, file_object, **kwargs):
#    json_data = {}
#    json_data['phase'] = data['phase']
#    json_data['hit_table'] = data['hit_table'].get_json_data()
#    json_data['s_walls'] = [s_wall.get_json_data()
#                            for s_wall in data['s_walls']]
#    json_data['joints'] = [joint.get_json_data()
#                           for joint in data['joints']]
#    json.dump(json_data, file_object, **kwargs)
#    
#
#def load_spectral_network_data(file_object, sw_curve, sw_diff, config_data,
#                               **kwargs):
#    json_data = json.load(file_object, **kwargs)
#    spectral_network = SpectralNetwork(sw_curve, sw_diff, json_data['phase'],
#                                       config_data,)
#    spectral_network.hit_table.load_json_data(json_data['hit_table'])
#    for joint_data in json_data['joints']:
#        a_joint = Joint()
#        spectral_network.joints.append(a_joint.set_json_data(joint_data))
#
#    for s_wall_data in json_data['s_walls']:
#        an_s_wall = SWall()
#        spectral_network.s_walls.append(an_s_wall.set_json_data(s_wall_data))
#
#    # NOTE: Loaded data indicates parents as labels,
#    # but they should be changed to the corresponding instances.
