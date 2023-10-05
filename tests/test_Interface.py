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
        
        self.interface.setSourceDict(ind.loadSZsource)
        self.assertIsNotNone(self.interface.sourceDict) 

    def test_setTelescopeDict(self):
        telescope = {}

        self.assertIsNone(self.interface.telescopeDict) 
        self.assertRaises(FieldError, self.interface.setTelescopeDict, telescope)
        
        self.interface.setTelescopeDict(ind.telescope)

        self.assertIsNotNone(self.interface.telescopeDict) 
    
    def test_setInstrumentDict(self):
        instrument = {}

        self.assertIsNone(self.interface.instrumentDict) 
        self.assertRaises(FieldError, self.interface.setInstrumentDict, instrument)
        
        self.interface.setInstrumentDict(ind.instrument)

        self.assertIsNotNone(self.interface.instrumentDict) 
    
    def test_setAtmosphereDict(self):
        atmosphere = {}

        self.assertIsNone(self.interface.atmosphereDict) 
        self.assertRaises(FieldError, self.interface.setAtmosphereDict, atmosphere)
        
        self.interface.setAtmosphereDict(ind.atmosphere)

        self.assertIsNotNone(self.interface.atmosphereDict) 
    
    def test_setObservationDict(self):
        observation = {}

        self.assertIsNone(self.interface.observationDict) 
        self.assertRaises(FieldError, self.interface.setObservationDict, observation)
        
        self.interface.setObservationDict(ind.observation)

        self.assertIsNotNone(self.interface.observationDict) 

    def test_runSimulation(self):
        self.interface.setObservationDict(ind.observation)
        self.interface.setAtmosphereDict(ind.atmosphere)
        self.interface.setInstrumentDict(ind.instrument)
        self.interface.setSourceDict(ind.SZsource)
        self.interface.setTelescopeDict(ind.telescope)
        
        out = self.interface.runSimulation()


if __name__ == "__main__":
    import nose2
    nose2.main()


