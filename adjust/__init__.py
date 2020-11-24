import numpy as np
import pandas as pd


class Daoxian:
    def __init__(self, info1='', info2='', info3=''):
        self.rho = 206265
        self.info1 = info1
        self.info2 = info2
        self.info3 = info3
        self.edges = pd.DataFrame({"node1": [], "node2": [], "value": []})
        self.corners = pd.DataFrame({"node1": [], "node2": [], "node3": [], "value": []})
        self.known_poi = pd.DataFrame({"poi_name": [], "poi_x": [], "poi_y": []})
        self.error = {"edge": [], "corner": []}
        self.params = pd.DataFrame({'ox': [], 'oy': []})
        self.p_weight_matrix = []

    def check_info(self):
        # if self.info1.endswith('\n'):
        #     self.info1 = self.info1[0:-1]
        edges = self.info1.split('\n')
        edges_observation = {"node1": [], "node2": [], "value": []}
        for edge in edges:
            if len(edge) == 0:
                pass
            elif len(edge) == 11:
                edges_observation['node1'].append(int(edge[1]))
                edges_observation['node2'].append(int(edge[6]))
                edges_observation['value'].append(0)
            elif len(edge) > 11:
                edges_observation['node1'].append(int(edge[1]))
                edges_observation['node2'].append(int(edge[6]))
                edges_observation['value'].append(float(edge[8:]))
        self.edges = pd.DataFrame(edges_observation)

        corners = self.info2.split('\n')
        corners_observation = {"node1": [], "node2": [], "node3": [], "value": []}
        for corner in corners:
            if len(corner) == 0:
                pass
            elif len(corner) > 15:
                corners_observation['node1'].append(int(corner[1]))
                corners_observation['node2'].append(int(corner[3]))
                corners_observation['node3'].append(int(corner[5]))
                value_str = corner[7:].split('°')[1]
                corner_value = 3600 * (float(corner[7:].split('°')[0]) +
                                       float(value_str[0:2]) / 60 + float(value_str[3:-1]) / 3600)
                corners_observation['value'].append(corner_value)
        self.corners = pd.DataFrame(corners_observation)

        known = self.info3.split('\n')
        known = [i for i in known if len(i) > 2]
        self.error['edge'] = eval(known[2][5:])
        self.error['corner'] = eval(known[3][5:7])
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
            # index = [i for i in range(len(ox))]
            for i in self.known_poi['poi_name']:
                ox[i] = 0
                oy[i] = 0
            params = pd.DataFrame({'ox': ox, 'oy': oy})
        self.params = params
        return True

    def init_error(self):
        p_weight_edge = [np.power(self.error['corner'] / self.error['edge'] / i / 100, 2) for i in self.edges['value']
                         if i > 0]
        p_weight_matrix = [1 for i in range(len(self.corners))]
        p_weight_matrix.extend(p_weight_edge)
        self.p_weight_matrix = np.diag(p_weight_matrix)

    def cal(self):
        b_coefficient = np.zeros(np.shape(self.p_weight_matrix))
        params = self.params

        edge_function = self.edges
        edge_function = edge_function[~edge_function['value'].isin([0])]

        edge_function_index = [edge_function['node1']]

        l_observation_corner = [i for i in self.corners['value']]
        l_observation_edge = [i for i in self.edges['value'] if i > 0]
        l_observation_corner.extend(l_observation_edge)
        l_observation = np.array(l_observation_corner).reshape((-1, 1))

        x_params = [i for i in self.params['ox'] if i != 0]
        x_params.extend([i for i in self.params['oy'] if i != 0])

        for i in range(100):
            _all_xy_coordinate = self.params
            for poi in range(len(self.known_poi)):
                _all_xy_coordinate['ox'][self.known_poi['poi_name'][poi]] = self.known_poi['poi_x'][poi]
                _all_xy_coordinate['oy'][self.known_poi['poi_name'][poi]] = self.known_poi['poi_y'][poi]

            _delta_x_jk = _all_xy_coordinate['ox'][edge_function['node1']].values - _all_xy_coordinate['ox'][
                edge_function['node2']].values
            _delta_y_jk = _all_xy_coordinate['oy'][edge_function['node1']].values - _all_xy_coordinate['ox'][
                edge_function['node2']].values
            _s_jk = np.power(np.power(_delta_x_jk, 2) + np.power(_delta_y_jk, 2), 0.5)
            l_observation_residual = l_observation
        # for corner

        # for edge

        return True
