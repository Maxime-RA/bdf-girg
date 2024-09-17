from bdfs.BDF import BDF


class D(BDF):
    i = 0

    def __init__(self, i):
        super().__init__(1)
        self.i = i

    def get_bdf_value(self, vector):
        return vector[self.i]

    def get_depth_vol(self):
        return 1

    def get_depth_com(self):
        return 1

    def get_length_vol(self):
        return 1

    def get_length_com(self):
        return 1

    def get_dimensions(self):
        return [self.i]

    def get_optimal_bdf(self, n):
        if n == 0:
            return 1, self
        if n == 1:
            return 0, None
        raise AttributeError("1d-bdf cannot be shortened more")

    def get_min_max_form(self):
        return {(self.i,)}

    def get_volume_poly(self):
        return [2, 0]

    def get_simplified_poly(self, precision):
        return self.get_volume_poly()

    def __str__(self):
        return str(self.i)
