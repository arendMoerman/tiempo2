"""!
@file
Test TiEMPO2 interface.
"""

import os
import numpy as np
import unittest
from tiempo2.Interface import Interface, FieldError

import InputDicts as ind

class TestInterface(unittest.TestCase):
    def setUp(self):
        self.interface = Interface()

    def test_setSourceDict(self):
        source = {}

        self.assertIsNone(self.interface.sourceDict) 
        self.assertRaises(FieldError, self.interface.setSourceDict, source)
        
        self.interface.setSourceDict(ind.SZsource)

        self.assertIsNotNone(self.interface.sourceDict) 

    def test_setTelescopeDict(self):
        telescope = {}

        self.assertIsNone(self.interface.telescopeDict) 
        self.assertRaises(FieldError, self.interface.setTelescopeDict, telescope)
        
        self.interface.setTelescopeDict(ind.telescope)

        self.assertIsNotNone(self.interface.telescopeDict) 
if __name__ == "__main__":
    import nose2
    nose2.main()


