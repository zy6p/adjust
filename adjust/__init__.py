import numpy as np
import pandas as pd


class Daoxian:
    def __init__(self, info1='', info2='', info3=''):
        self.info1 = info1
        self.info2 = info2
        self.info3 = info3
        self.edges = []
        self.corners = []
        self.known_poi = []
        self.error = {"edge": [], "corner": []}
        self.params = []

    def check_info(self):
        # if self.info1.endswith('\n'):
        #     self.info1 = self.info1[0:-1]
        edges = self.info1.split('\n')
        edges_observation = {"one": [], "two": [], "value": []}
        for edge in edges:
            if len(edge) == 0:
                pass
            elif len(edge) == 11:
                edges_observation['one'].append(int(edge[1]))
                edges_observation['two'].append(int(edge[6]))
                edges_observation['value'].append(0)
            elif len(edge) > 11:
                edges_observation['one'].append(int(edge[1]))
                edges_observation['two'].append(int(edge[6]))
                edges_observation['value'].append(float(edge[8:]))
        self.edges = pd.DataFrame(edges_observation)

        corners = self.info2.split('\n')
        corners_observation = {"one": [], "two": [], "three": [], "value": []}
        for corner in corners:
            if len(corner) == 0:
                pass
            elif len(corner) > 15:
                corners_observation['one'].append(int(corner[1]))
                corners_observation['two'].append(int(corner[3]))
                corners_observation['three'].append(int(corner[5]))
                value_str = corner[7:].split('°')[1]
                corner_value = np.pi / 180 * (float(corner[7:].split('°')[0]) +
                                              float(value_str[0:2]) / 60 + float(value_str[3:-1]) / 3600)
                corners_observation['value'].append(corner_value)
        self.corners = pd.DataFrame(corners_observation)

        known = self.info3.split('\n')
        known = [i for i in known if len(i) > 2]
        self.error['edge'] = eval(known[2][5:])
        self.error['corner'] = np.pi / 180 * float(known[3][5:7])
        known_poi = {"poi_name": [], "poi_x": [], "poi_y": []}
        for poi in known[0:2]:
            poi_info = poi.split('=')
            poi_name = int(poi_info[0][-3])
            poi_x = float(poi_info[1].split('(')[0])
            poi_y = float(poi_info[1].split('(')[0])
            known_poi['poi_name'].append(poi_name)
            known_poi['poi_x'].append(poi_x)
            known_poi['poi_y'].append(poi_y)
        self.known_poi = pd.DataFrame(known_poi)
        return True

    def init_params(self, params=None):
        if params is None:
            ox = [269.011, 163.475, 152.66, 104.35, 172.45, 278.44, 164.65, 193.77, 284.43]
            oy = [292.728, 292.866, 264.52, 183.83, 123.89, 124.98, 253.76, 209.33, 242.22]
            params = pd.DataFrame(ox, oy, columns=['ox', 'oy'])
        self.params = params

    def cal(self):




        return True
