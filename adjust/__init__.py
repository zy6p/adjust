import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import copy


class Daoxian:
    def __init__(self, info1='', info2='', info3=''):
        self.rho = 1
        self.info1 = info1
        self.info2 = info2
        self.info3 = info3
        self.edges = pd.DataFrame({"node1": [], "node2": [], "value": []})
        self.corners = pd.DataFrame({"node1": [], "node2": [], "node3": [], "value": []})
        self.known_poi = pd.DataFrame({"poi_name": [], "poi_x": [], "poi_y": []})
        self.error = {"edge": [], "corner": []}
        self.params = pd.DataFrame({'ox': [], 'oy': [], 'o_index': [], 'ox_index': [], 'oy_index': []})
        self.params_xy = {"oxy": [], 'oxy_index': []}
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
                pass
                # edges_observation['node1'].append(int(edge[1]))
                # edges_observation['node2'].append(int(edge[6]))
                # edges_observation['value'].append(0)
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
                corner_value = np.pi / 180.0 * (float(corner[7:].split('°')[0]) +
                                                float(value_str[0:2]) / 60 + float(value_str[3:-1]) / 3600)
                # corner_value = 3600 * (float(corner[7:].split('°')[0]) +
                #                        float(value_str[0:2]) / 60 + float(value_str[3:-1]) / 3600)
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
            poi_y = float(poi_info[2].split('(')[0])
            known_poi['poi_name'].append(poi_name)
            known_poi['poi_x'].append(poi_x)
            known_poi['poi_y'].append(poi_y)
        self.known_poi = pd.DataFrame(known_poi)
        return True

    def init_params(self, params=None):
        if params is None:
            ox = [284.43, 269.011, 163.475, 152.66, 97.35, 184.45, 278.44, 192.65, 178.77]
            oy = [242.22, 292.728, 292.866, 264.52, 155.83, 90.89, 124.98, 225.76, 225.33]
            ox_index = [i for i in range(len(ox))]
            oy_index = [i for i in range(len(ox) - 2, len(ox) * 2 - 2)]
            for ii in range(len(self.known_poi)):
                poi = self.known_poi['poi_name'][ii]
                ox[poi] = 0  # self.known_poi['poi_x'][ii]
                oy[poi] = 0  # self.known_poi['poi_y'][ii]
                ox_index = [i if i < poi else i - 1 for i in ox_index]
                oy_index = [i if i < poi + 7 else i - 1 for i in oy_index]
                ox_index[poi] = -1
                oy_index[poi] = -1
            oxy, oxy_index = copy.deepcopy(ox), copy.deepcopy(ox_index)
            oxy.extend(oy)
            oxy = [i for i in oxy if i != 0]
            oxy_index.extend(oy_index)
            oxy_index = [i for i in oxy_index if i != -1]
            params = pd.DataFrame({'ox': ox, 'oy': oy, 'ox_index': ox_index, 'oy_index': oy_index})
            self.params_xy = {"oxy": oxy, 'oxy_index': oxy_index}
            self.params = params
        return True

    def init_weight(self):
        p_weight_edge = [np.power(self.error['corner'] / self.error['edge'] / i / 1000, 2)
                         for i in self.edges['value']]
        p_weight_matrix = [1 for i in range(len(self.corners))]
        p_weight_matrix.extend(p_weight_edge)
        self.p_weight_matrix = np.diag(p_weight_matrix)

    def cal(self):
        n_num_edge = len(self.edges)
        n_num_corner = len(self.corners)
        n_num = n_num_edge + n_num_corner
        x_params = np.array(self.params_xy['oxy']).reshape((-1, 1))
        t_num = len(x_params)
        t_num_half = int(t_num / 2)
        known_num = len(self.known_poi)

        b_coefficient = np.zeros((n_num, t_num))

        edge_function = self.edges
        # edge_function = edge_function[~edge_function['value'].isin([0])]
        # edge_function_index = [edge_function['node1']]

        corner_function = self.corners

        l_observation_corner = [i for i in self.corners['value']]
        l_observation_edge = [i for i in self.edges['value']]
        l_observation_corner.extend(l_observation_edge)
        l_observation = np.array(l_observation_corner).reshape((-1, 1))

        l_observation_residual = l_observation

        for iteration in range(100):
            _all_xy_coordinate = {'ox': x_params[0:t_num_half], 'oy': x_params[t_num_half:]}
            for poi in range(known_num):
                _all_xy_coordinate['ox'] = np.insert(_all_xy_coordinate['ox'], self.known_poi['poi_name'][poi],
                                                     self.known_poi['poi_x'][poi])
                _all_xy_coordinate['oy'] = np.insert(_all_xy_coordinate['oy'], self.known_poi['poi_name'][poi],
                                                     self.known_poi['poi_y'][poi])
            # for edge
            _edge_delta_x_jk = _all_xy_coordinate['ox'][edge_function['node1']] - _all_xy_coordinate[
                'ox'][edge_function['node2']]
            _edge_delta_y_jk = _all_xy_coordinate['oy'][edge_function['node1']] - _all_xy_coordinate[
                'oy'][edge_function['node2']]
            _edge_s_jk = np.power(np.power(_edge_delta_x_jk, 2) + np.power(_edge_delta_y_jk, 2), 0.5)
            for ii in range(n_num_edge):
                b_coefficient[ii + n_num_corner][self.params['ox_index'][self.edges['node1'][ii]]] =\
                    - _edge_delta_x_jk[ii] / _edge_s_jk[ii]
                b_coefficient[ii + n_num_corner][self.params['oy_index'][self.edges['node1'][ii]]] =\
                    - _edge_delta_y_jk[ii] / _edge_s_jk[ii]
                b_coefficient[ii + n_num_corner][self.params['ox_index'][self.edges['node2'][ii]]] =\
                    _edge_delta_x_jk[ii] / _edge_s_jk[ii]
                b_coefficient[ii + n_num_corner][self.params['oy_index'][self.edges['node2'][ii]]] =\
                    _edge_delta_y_jk[ii] / _edge_s_jk[ii]
                l_observation[ii + n_num_corner] = l_observation[ii + n_num_corner] - _edge_s_jk[ii]
            # for corner
            _corner_delta_x_jk = _all_xy_coordinate['ox'][corner_function['node2']] - _all_xy_coordinate['ox'][
                corner_function['node3']]
            _corner_delta_y_jk = _all_xy_coordinate['oy'][corner_function['node2']] - _all_xy_coordinate['oy'][
                corner_function['node3']]
            _corner_delta_x_jh = _all_xy_coordinate['ox'][corner_function['node2']] - _all_xy_coordinate['ox'][
                corner_function['node1']]
            _corner_delta_y_jh = _all_xy_coordinate['oy'][corner_function['node2']] - _all_xy_coordinate['oy'][
                corner_function['node1']]
            _corner_s_jk_power2 = np.power(_corner_delta_x_jk, 2) + np.power(_corner_delta_y_jk, 2)
            _corner_s_jh_power2 = np.power(_corner_delta_x_jh, 2) + np.power(_corner_delta_y_jh, 2)
            _corner_a_jk = self.rho * _corner_delta_y_jk / _corner_s_jk_power2
            _corner_b_jk = - self.rho * _corner_delta_x_jk / _corner_s_jk_power2
            _corner_a_jh = self.rho * _corner_delta_y_jh / _corner_s_jh_power2
            _corner_b_jh = - self.rho * _corner_delta_x_jh / _corner_s_jh_power2
            for ii in range(n_num_corner):
                b_coefficient[ii][self.params['ox_index'][self.corners['node1'][ii]]] = _corner_a_jh[ii]
                b_coefficient[ii][self.params['oy_index'][self.corners['node1'][ii]]] = _corner_b_jh[ii]
                b_coefficient[ii][self.params['ox_index'][self.corners['node2'][ii]]] = _corner_a_jk[ii] - _corner_a_jh[ii]
                b_coefficient[ii][self.params['oy_index'][self.corners['node2'][ii]]] = _corner_b_jk[ii] - _corner_b_jh[ii]
                b_coefficient[ii][self.params['ox_index'][self.corners['node3'][ii]]] = - _corner_a_jk[ii]
                b_coefficient[ii][self.params['oy_index'][self.corners['node3'][ii]]] = - _corner_b_jk[ii]
                l_observation_residual[ii] = l_observation[ii] - np.arctan(_corner_delta_y_jk[ii] / _corner_delta_x_jk[
                    ii]) + np.arctan(_corner_delta_y_jh[ii] / _corner_delta_x_jh[ii])
            # l_observation_residual = l_observation + np.dot(b_coefficient, x_params)
            x_params_correct = np.dot(np.dot(np.dot(np.linalg.inv(np.dot(np.dot(
                b_coefficient.T, self.p_weight_matrix), b_coefficient)), b_coefficient.T), self.p_weight_matrix),
                l_observation_residual)
            x_params += x_params_correct
            plt.plot(_all_xy_coordinate['oy'], _all_xy_coordinate['ox'])
            plt.show()
            print(x_params)

        print(x_params)


        return True
