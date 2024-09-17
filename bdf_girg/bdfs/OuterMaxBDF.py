from bdfs.BDF import BDF
from itertools import product
import numpy as np


def create_outer_max(k1, k2):
    if k1 is not None and k2 is not None:
        return OuterMax(k1, k2)
    if k1 is not None:
        return k1
    return k2


class OuterMax(BDF):
    k1 = None
    k2 = None
    opt_length = []
    opt_bdf = []

    def __init__(self, k1, k2):
        super().__init__(k1.dimension + k2.dimension)
        self.k1 = k1
        self.k2 = k2
        self.opt_length = [None] * self.get_depth_vol()
        self.opt_bdf = [None] * self.get_depth_vol()

    def get_bdf_value(self, vector):
        return max(self.k1.get_bdf_value(vector),
                   self.k2.get_bdf_value(vector))

    def get_depth_vol(self):
        return self.k1.get_depth_vol() + self.k2.get_depth_vol()

    def get_depth_com(self):
        return self.k1.get_depth_com() + self.k2.get_depth_com()

    def get_length_vol(self):
        return self.k1.get_length_vol() * self.k2.get_length_vol()

    def get_length_com(self):
        return self.k1.get_length_com() * self.k2.get_length_com()

    def get_dimensions(self):
        return self.k1.get_dimensions() + self.k2.get_dimensions()

    def get_optimal_bdf(self, n):
        if n == self.get_depth_vol():
            return 0, None
        if self.opt_length[n] is not None:
            return self.opt_length[n], self.opt_bdf[n]

        d_v1 = self.k1.get_depth_vol()
        d_v2 = self.k2.get_depth_vol()
        minimal_length = float('inf')
        minimal_bdf = None

        for i in range(max(0, n - d_v2), min(n, d_v1) + 1):
            k1_len, k1_bdf = self.k1.get_optimal_bdf(i)
            k2_len, k2_bdf = self.k2.get_optimal_bdf(n - i)
            length = max(k1_len, 1) * max(k2_len, 1)
            if minimal_length > length:
                minimal_length = length
                minimal_bdf = create_outer_max(k1_bdf, k2_bdf)
        self.opt_length[n] = minimal_length
        self.opt_bdf[n] = minimal_bdf
        return minimal_length, minimal_bdf

    def get_min_max_form(self):
        return {a + b for a, b in product(self.k1.get_min_max_form(), self.k2.get_min_max_form())}

    def get_volume_poly(self):
        return np.polymul(self.k1.get_volume_poly(), self.k2.get_volume_poly())

    def get_simplified_poly(self, inter_depth):
        return np.polymul(self.k1.get_simplified_poly(inter_depth), self.k2.get_simplified_poly(inter_depth))

    def __str__(self):
        return "max({},{})".format(str(self.k1), str(self.k2))
