"""
@file
Test atmosphere module.
"""

import os
import numpy as np
import unittest
import tiempo2.Atmosphere as TAtm

#class TestAtmosphere(unittest.TestCase):
#    pass
class TestAtmosphere(unittest.TestCase):
    def test_prepAtmospherePWV(self):
        atmD = {
                "path" : "/home/arend/Projects/Simulations/tiempo2/thies/aris/",
                "filename" : "sample00.dat",
                "PWV0"      : 1,
                }

        telD = {
                "Dtel"  : 10
                }

        test = TAtm.prepAtmospherePWV(atmD, telD)

        print(test)

    @classmethod
    def tearDownClass(self):
        pass

if __name__ == "__main__":
    import nose2
    nose2.main(verbosity=10)

