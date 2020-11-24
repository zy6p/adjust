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
        self.error = []
        self.params = []

    def check_info(self):
        if self.info1.endswith('\n'):
            self.info1 = self.info1[0:-1]
        edges = self.info1.split('\n')
        edge_observation = {"begin": [], "end": [], "value": []}
        for edge in edges:
            if len(edge) == 0:
                pass
            elif len(edge) == 11:
                pass
        return True

    def init_params(self, params=[]):
        self.params = params

    def cal(self):
        return True


