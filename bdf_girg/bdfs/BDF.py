from abc import ABC, abstractmethod



class BDF(ABC):
    dimension = 0

    def __init__(self, dimension):
        self.dimension = dimension

    @abstractmethod
    def get_bdf_value(self, vector):
        pass

    @abstractmethod
    def get_depth_vol(self):
        pass

    @abstractmethod
    def get_depth_com(self):
        pass

    @abstractmethod
    def get_length_vol(self):
        pass

    @abstractmethod
    def get_length_com(self):
        pass

    @abstractmethod
    def get_dimensions(self):
        pass

    @abstractmethod
    def get_optimal_bdf(self, n):
        pass

    @abstractmethod
    def get_min_max_form(self):
        pass

    @abstractmethod
    def get_volume_poly(self):
        pass

    @abstractmethod
    def get_simplified_poly(self, precision):
        pass

    def __str__(self):
        return "Abstract BDF"
