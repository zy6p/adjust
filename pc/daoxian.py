import copy

import numpy as np
import pandas as pd

from .inverseadjust import DxwInverseAdjust


class Daoxian:
    def __init__(self, info1='', info2='', info3='', weight_ratio='1'):
        self.weight_ratio = eval(weight_ratio)  # deg 1 ss = length 1 mm
        self.info1 = info1.replace('\r\n', '\n')

        self.info2 = info2.replace('\r\n', '\n')
        self.info3 = info3.replace('\r\n', '\n')
        self.edges = pd.DataFrame({"node1": [], "node2": [], "value": []})
        self.corners = pd.DataFrame({"node1": [], "node2": [], "node3": [], "value": []})
        self.known_poi = pd.DataFrame({"poi_name": [], "poi_x": [], "poi_y": []})
        self.error = {"edge": [], "corner": []}
        self.params = pd.DataFrame(
            {'origin_x': [], 'origin_y': [], 'o_index': [], 'origin_x_index': [], 'origin_y_index': []})
        self.params_xy = {"origin_xy": [], 'origin_xy_index': []}
        self.p_weight_matrix = []
        self.intermediates = {'temp': []}

    def check_info(self):
        self.convert_info1_to_edge()

        self.convert_info2_to_corner()

        self.convert_info3_to_error()
        return True

    def convert_info3_to_error(self):
        known = self.info3.split('\n')
        known = [i for i in known if len(i) > 2]
        self.error['edge'] = eval(known[2].split('：')[-1])
        self.error['corner'] = eval(known[3].split('：')[-1][0:-1])
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

    def convert_info2_to_corner(self):
        corners = self.info2.split('\n')
        corners_observation = {"node1": [], "node2": [], "node3": [], "value": [], "deg_ddd": [], "deg_mm": [],
                               "deg_ss": []}
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
                corners_observation['value'].append(corner_value)
                corners_observation['deg_ddd'].append(float(corner[7:].split('°')[0]))
                corners_observation['deg_mm'].append(float(value_str[0:2]))
                corners_observation['deg_ss'].append(float(value_str[3:-1]))
        self.corners = pd.DataFrame(corners_observation)

    def convert_info1_to_edge(self):
        edges = self.info1.split('\n')
        edges_observation = {"node1": [], "node2": [], "value": []}
        for edge in edges:
            if len(edge) == 0:
                pass
            elif edge.endswith('已知值'):
                pass
            elif len(edge) > 10:
                edges_observation['node1'].append(int(edge[1]))
                edges_observation['node2'].append(int(edge[6]))
                edges_observation['value'].append(float(edge[8:]))
        self.edges = pd.DataFrame(edges_observation)

    # TODO: This calculation method is sensitive to the initial value, so need to calculate the initial value first.
    def init_params(self, params=None):
        if params is None:
            origin_x = [2.870990000000000e+02, 1.919990000000000e+02, 1.406959751853701e+02, 1.492484041733704e+02,
                        2.045984049476477e+02, 2.819792529204340e+02, 3.275893976516628e+02, 2.154350873574551e+02,
                        2.556998807055271e+02]
            origin_y = [2.551870000000000e+02, 2.924290000000000e+02, 2.341834543352806e+02, 1.456793271234235e+02,
                        81.853698181164750, 1.103393419062568e+02, 1.943975064036699e+02, 1.362573869669805e+02,
                        2.137389646671174e+02]
            origin_x_index = [i for i in range(len(origin_x))]
            origin_y_index = [i for i in range(len(origin_x) - 2, len(origin_x) * 2 - 2)]
            for ii in range(len(self.known_poi)):
                poi = self.known_poi['poi_name'][ii]
                origin_x[poi] = 0  # self.known_poi['poi_x'][ii]
                origin_y[poi] = 0  # self.known_poi['poi_y'][ii]
                origin_x_index = [i if i < poi else i - 1 for i in origin_x_index]
                origin_y_index = [i if i < poi + 7 else i - 1 for i in origin_y_index]
                origin_x_index[poi] = -1
                origin_y_index[poi] = -1
            origin_xy, origin_xy_index = copy.deepcopy(origin_x), copy.deepcopy(origin_x_index)
            origin_xy.extend(origin_y)
            origin_xy = [i for i in origin_xy if i != 0]
            origin_xy_index.extend(origin_y_index)
            origin_xy_index = [i for i in origin_xy_index if i != -1]
            params = pd.DataFrame({'origin_x': origin_x, 'origin_y': origin_y, 'origin_x_index': origin_x_index,
                                   'origin_y_index': origin_y_index})
            self.params_xy = {"origin_xy": origin_xy, 'origin_xy_index': origin_xy_index}
            self.params = params
        return True

    def init_weight(self):
        p_weight_edge = [np.power(self.error['corner'] / (i * self.error['edge'] / self.weight_ratio), 2)
                         for i in self.edges['value']]
        p_weight_matrix = [1 for i in range(len(self.corners))]
        p_weight_matrix.extend(p_weight_edge)
        self.p_weight_matrix = np.diag(p_weight_matrix)
        return True

    def cal(self):
        dxw_inverse_adjust = DxwInverseAdjust(self)
        x_params = np.array(self.params_xy['origin_xy']).reshape((-1, 1))
        b_coefficient = np.zeros((dxw_inverse_adjust.n_num, dxw_inverse_adjust.t_num))

        l_observation_corner = [i for i in self.corners['value']]
        l_observation_edge = [i for i in self.edges['value']]
        l_observation_corner.extend(l_observation_edge)
        l_observation = np.array(l_observation_corner).reshape((-1, 1))

        l_observation_residual = l_observation * 2.0

        for iteration in range(20):
            _all_xy_coordinate = {'origin_x': x_params[0:dxw_inverse_adjust.t_num_half],
                                  'origin_y': x_params[dxw_inverse_adjust.t_num_half:]}
            for poi in range(dxw_inverse_adjust.known_num):
                _all_xy_coordinate['origin_x'] = np.insert(_all_xy_coordinate['origin_x'],
                                                           self.known_poi['poi_name'][poi],
                                                           self.known_poi['poi_x'][poi])
                _all_xy_coordinate['origin_y'] = np.insert(_all_xy_coordinate['origin_y'],
                                                           self.known_poi['poi_name'][poi],
                                                           self.known_poi['poi_y'][poi])
            # for edge ?
            _edge_delta_x_jk = _all_xy_coordinate['origin_x'][self.edges['node2']] - _all_xy_coordinate[
                'origin_x'][self.edges['node1']]
            _edge_delta_y_jk = _all_xy_coordinate['origin_y'][self.edges['node2']] - _all_xy_coordinate[
                'origin_y'][self.edges['node1']]
            _edge_s_jk = np.power(np.power(_edge_delta_x_jk, 2) + np.power(_edge_delta_y_jk, 2), 0.5)
            _edge_c_jk = _edge_delta_x_jk / _edge_s_jk
            _edge_d_jk = _edge_delta_y_jk / _edge_s_jk
            for ii in range(dxw_inverse_adjust.n_num_edge):
                if self.params['origin_x_index'][self.edges['node1'][ii]] >= 0:
                    b_coefficient[ii + dxw_inverse_adjust.n_num_corner][
                        self.params['origin_x_index'][self.edges['node1'][ii]]] = - 1000 * _edge_c_jk[ii]
                    b_coefficient[ii + dxw_inverse_adjust.n_num_corner][
                        self.params['origin_y_index'][self.edges['node1'][ii]]] = - 1000 * _edge_d_jk[ii]
                if self.params['origin_x_index'][self.edges['node2'][ii]] >= 0:
                    b_coefficient[ii + dxw_inverse_adjust.n_num_corner][
                        self.params['origin_x_index'][self.edges['node2'][ii]]] = 1000 * _edge_c_jk[ii]
                    b_coefficient[ii + dxw_inverse_adjust.n_num_corner][
                        self.params['origin_y_index'][self.edges['node2'][ii]]] = 1000 * _edge_d_jk[ii]
                l_observation_residual[ii + dxw_inverse_adjust.n_num_corner] = 1000 * (
                        l_observation[ii + dxw_inverse_adjust.n_num_corner] - _edge_s_jk[ii])
            # for corner √
            _corner_delta_x_jk = _all_xy_coordinate['origin_x'][self.corners['node2']] - \
                                 _all_xy_coordinate['origin_x'][self.corners['node3']]
            _corner_delta_y_jk = _all_xy_coordinate['origin_y'][self.corners['node2']] - \
                                 _all_xy_coordinate['origin_y'][self.corners['node3']]
            _corner_delta_x_jh = _all_xy_coordinate['origin_x'][self.corners['node2']] - \
                                 _all_xy_coordinate['origin_x'][self.corners['node1']]
            _corner_delta_y_jh = _all_xy_coordinate['origin_y'][self.corners['node2']] - \
                                 _all_xy_coordinate['origin_y'][self.corners['node1']]
            _corner_s_jk_power2 = np.power(_corner_delta_x_jk, 2) + np.power(_corner_delta_y_jk, 2)
            _corner_s_jh_power2 = np.power(_corner_delta_x_jh, 2) + np.power(_corner_delta_y_jh, 2)
            _corner_a_jk = _corner_delta_y_jk / _corner_s_jk_power2
            _corner_b_jk = - _corner_delta_x_jk / _corner_s_jk_power2
            _corner_a_jh = _corner_delta_y_jh / _corner_s_jh_power2
            _corner_b_jh = - _corner_delta_x_jh / _corner_s_jh_power2
            _corner_jh = np.arctan2(- _corner_delta_y_jh, - _corner_delta_x_jh)
            _corner_jk = np.arctan2(- _corner_delta_y_jk, - _corner_delta_x_jk)
            for ii in range(dxw_inverse_adjust.n_num_corner):
                if self.params['origin_x_index'][self.corners['node1'][ii]] >= 0:
                    b_coefficient[ii][self.params['origin_x_index'][self.corners['node1'][ii]]] = - 206065 * \
                                                                                                  _corner_a_jh[ii]
                    b_coefficient[ii][self.params['origin_y_index'][self.corners['node1'][ii]]] = - 206065 * \
                                                                                                  _corner_b_jh[ii]
                if self.params['origin_x_index'][self.corners['node2'][ii]] >= 0:
                    b_coefficient[ii][self.params['origin_x_index'][self.corners['node2'][ii]]] = 206065 * (
                            - _corner_a_jk[ii] + _corner_a_jh[ii])
                    b_coefficient[ii][self.params['origin_y_index'][self.corners['node2'][ii]]] = 206065 * (
                            - _corner_b_jk[ii] + _corner_b_jh[ii])
                if self.params['origin_x_index'][self.corners['node3'][ii]] >= 0:
                    b_coefficient[ii][self.params['origin_x_index'][self.corners['node3'][ii]]] = 206065 * _corner_a_jk[
                        ii]
                    b_coefficient[ii][self.params['origin_y_index'][self.corners['node3'][ii]]] = 206065 * _corner_b_jk[
                        ii]
                l_observation_residual[ii] = 206265 * (l_observation[ii] - np.mod(_corner_jk[ii] - _corner_jh[ii],
                                                                                  2 * np.pi))
            x_params_correct = np.dot(np.dot(np.dot(np.linalg.inv(np.dot(np.dot(
                b_coefficient.T, self.p_weight_matrix), b_coefficient)), b_coefficient.T), self.p_weight_matrix),
                l_observation_residual)
            x_params += x_params_correct
            # plt.plot(_all_xy_coordinate['origin_y'], _all_xy_coordinate['origin_x'])
            # plt.show()

        # cal weight intermediates result
        v_residual = - l_observation_residual
        sigma_final = np.sqrt(np.dot(np.dot(v_residual.T, self.p_weight_matrix), v_residual) / (
                dxw_inverse_adjust.n_num - dxw_inverse_adjust.t_num))[0, 0]
        q_xx_covariance = np.linalg.inv(np.dot(np.dot(b_coefficient.T, self.p_weight_matrix), b_coefficient))
        q_ll_covariance = np.dot(np.dot(b_coefficient, q_xx_covariance), b_coefficient.T)
        edge_sigma_max = np.sqrt(
            np.max(np.diag(q_ll_covariance)[- dxw_inverse_adjust.n_num_edge:])) * sigma_final * self.weight_ratio
        corner_sigma_max = np.sqrt(np.max(np.diag(q_ll_covariance)[0:dxw_inverse_adjust.n_num_corner])) * sigma_final
        _q_point_covariance = np.diag(q_xx_covariance).reshape((-1, 2))
        q_point_covariance = np.power(_q_point_covariance[:, 0] + _q_point_covariance[:, 1], 0.5) * sigma_final
        for poi in range(dxw_inverse_adjust.known_num):
            q_point_covariance = np.insert(q_point_covariance, self.known_poi['poi_name'][poi], 0)
        self.params['point_sigma'] = q_point_covariance
        self.error['point_sigma_max'] = max(q_point_covariance)
        self.error['edge_sigma_max'] = edge_sigma_max
        self.error['corner_sigma_max'] = corner_sigma_max
        self.error['sigma_final'] = sigma_final

        self.edges['edge_adjust'] = _edge_s_jk
        self.edges['edge_adjust_residual_mm'] = - l_observation_residual[- dxw_inverse_adjust.n_num_edge:]
        self.edges['edge_post_test_error'] = np.sqrt(
            np.diag(q_ll_covariance)[- dxw_inverse_adjust.n_num_edge:]) * self.weight_ratio * sigma_final
        self.corners['corner_adjust_rad'] = self.corners['value'].values.reshape((-1, 1)) - l_observation_residual[
                                                                                            0:dxw_inverse_adjust.n_num_corner] / 206265
        _corner_adjust_deg = self.corners['corner_adjust_rad'] * 180 / np.pi
        self.corners['corner_adjust_deg_ddd'] = np.floor(_corner_adjust_deg.values)
        self.corners['corner_adjust_deg_mm'] = np.floor(
            (_corner_adjust_deg.values - self.corners['corner_adjust_deg_ddd'].values) * 60)
        self.corners['corner_adjust_deg_ss'] = (_corner_adjust_deg.values - self.corners[
            'corner_adjust_deg_ddd'].values) * 3600 - self.corners['corner_adjust_deg_mm'].values * 60
        self.corners['corner_adjust_residual_deg_ss'] = - l_observation_residual[0:dxw_inverse_adjust.n_num_corner]
        self.corners['corner_post_test_error'] = np.sqrt(
            np.diag(q_ll_covariance)[0:dxw_inverse_adjust.n_num_corner]) * sigma_final

        self.intermediates['b_coefficient'] = pd.DataFrame(b_coefficient).to_html()
        self.intermediates['l_observation_residual'] = pd.DataFrame(l_observation_residual).to_html()
        self.intermediates['v_residual'] = pd.DataFrame(v_residual).to_html()
        self.intermediates['P_weight_matrix'] = pd.DataFrame(self.p_weight_matrix).to_html()
        # self.intermediates['x_params'] = x_params
        self.params['adjust_x'] = _all_xy_coordinate['origin_x']
        self.params['adjust_y'] = _all_xy_coordinate['origin_y']
        return True

    def html_of_params(self):
        return self.params.to_html(columns=['adjust_x', 'adjust_y', 'point_sigma'])

    def html_of_edges(self):
        return self.edges.to_html()

    def html_of_corners(self):
        return self.corners.to_html()
