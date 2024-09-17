from bdfs.BDF import BDF
import numpy as np

class OuterMin(BDF):
    k1 = None
    k2 = None
    opt_length = []
    opt_bdf = []

    def __init__(self, k1, k2):
        super().__init__(k1.dimension + k2.dimension)
        self.k1 = k1
        self.k2 = k2
        self.opt_length = [None] * (self.get_depth_vol())
        self.opt_bdf = [None] * (self.get_depth_vol())

    def get_bdf_value(self, vector):
        return min(self.k1.get_bdf_value(vector),
                   self.k2.get_bdf_value(vector))

    def get_depth_vol(self):
        return min(self.k1.get_depth_vol(), self.k2.get_depth_vol())

    def get_depth_com(self):
        return max(self.k1.get_depth_com(), self.k2.get_depth_com())

    def get_length_vol(self):
        if self.k1.get_depth_vol() > self.k2.get_depth_vol():
            return self.k2.get_length_vol()
        elif self.k1.get_depth_vol() < self.k2.get_depth_vol():
            return self.k1.get_length_vol()
        return self.k1.get_length_vol() + self.k2.get_length_vol()

    def get_length_com(self):
        return self.k1.get_length_com() + self.k2.get_length_com()

    def get_dimensions(self):
        return self.k1.get_dimensions() + self.k2.get_dimensions()

    def get_optimal_bdf(self, n):
        if n == self.get_depth_vol():
            return 0, None
        if self.opt_length[n] is not None:
            return self.opt_length[n], self.opt_bdf[n]

        d_v1 = self.k1.get_depth_vol()
        d_v2 = self.k2.get_depth_vol()
        if d_v1 >= d_v2:
            k1_len, k1_bdf = self.k1.get_optimal_bdf(n + d_v1 - d_v2)
            k2_len, k2_bdf = self.k2.get_optimal_bdf(n)
            self.opt_length[n] = k1_len + k2_len
            self.opt_bdf[n] = OuterMin(k1_bdf, k2_bdf)
        else:
            k1_len, k1_bdf = self.k1.get_optimal_bdf(n)
            k2_len, k2_bdf = self.k2.get_optimal_bdf(n + d_v2 - d_v1)
            self.opt_length[n] = k1_len + k2_len
            self.opt_bdf[n] = OuterMin(k1_bdf, k2_bdf)

        return self.opt_length[n], self.opt_bdf[n]

    def get_min_max_form(self):
        return self.k1.get_min_max_form() | self.k2.get_min_max_form()

    def get_volume_poly(self):
        p1 = self.k1.get_volume_poly()
        p2 = self.k2.get_volume_poly()

        # First compute sum, then the intersection
        sum_poly = np.polyadd(p1, p2)
        prod_poly = np.polymul(p1, p2)
        return np.polysub(sum_poly, prod_poly)

    def get_simplified_poly(self, inter_depth):
        p1 = self.k1.get_simplified_poly(inter_depth-1)
        p2 = self.k2.get_simplified_poly(inter_depth-1)

        sum_poly = np.polyadd(p1, p2)
        if inter_depth > 0:
            prod_poly = np.polymul(p1, p2)
        else:
            prod_poly = [0]
        return np.polysub(sum_poly, prod_poly)

    def __str__(self):
        return "min({},{})".format(str(self.k1), str(self.k2))
