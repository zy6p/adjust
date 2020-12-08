import numpy as np
import pandas as pd
from .adjustfather import AdjustFather


class InverseAdjust(AdjustFather):
    def __init__(self):
        self.t_num = 1


class DxwInverseAdjust(InverseAdjust):
    def __init__(self, dxw):
        # super().__init__()
        self.n_num_corner = len(dxw.corners)
        self.n_num_edge = len(dxw.edges)
        self.n_num = self.n_num_edge + self.n_num_corner
        self.t_num = len(dxw.params_xy['origin_xy'])
        self.known_num = len(dxw.known_poi)
        self.t_num_half = int(self.t_num / 2)

