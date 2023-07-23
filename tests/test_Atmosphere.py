"""
@file
Test atmosphere module.
"""

import os
import numpy as np
import unittest
from tiempo2.Atmosphere import Atmosphere

class TestAtmosphere(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        prefix_filename = "test" 
        path_data = os.path.join(os.path.dirname(__file__), "resources")
        pwv0 = 1. 
        pwvgrid = 0.2 
        max_windspeed = 10 
        obs_duration = 420
        dish_radius = 5 # m
        max_num_strips = 1
        x_length_strip = 69

        cls.atm = Atmosphere(prefix_filename, path_data, pwv0, pwvgrid, max_windspeed, obs_duration, dish_radius, max_num_strips, x_length_strip)

if __name__ == "__main__":
    import nose2
    nose2.main()

