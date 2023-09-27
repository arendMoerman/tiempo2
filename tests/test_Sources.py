"""!
@file
File for testing the sources used in TiEMPO2.
"""

import os
import shutil
import numpy as np
import unittest
from nose2.tools import params
import tiempo2.Sources as TSources
import InputDicts as ind

class TestSources(unittest.TestCase):
    @params(("MJy", True), ("SI", False), ("mK", True))
    def test_generateSZMaps(self, unit, save):
        SZ, CMB, Az, El = TSources.generateSZMaps(ind.SZsource, 
                                                  convolve_beam=True,
                                                  ret_unit=unit,
                                                  save=save, 
                                                  sourceName="test", 
                                                  path=os.getcwd(), 
                                                  telescopeDict=ind.telescope)

        self.assertEqual(SZ.shape, (ind.SZsource.get("nAz"), ind.SZsource.get("nEl"), 
                                    ind.SZsource.get("freqs").size))
        self.assertEqual(CMB.shape, (ind.SZsource.get("nAz"), ind.SZsource.get("nEl"), 
                                    ind.SZsource.get("freqs").size))
        
        self.assertEqual(Az.shape, (ind.SZsource.get("nAz"),))
        self.assertEqual(El.shape, (ind.SZsource.get("nEl"),))
    
        if save:
            SZ_l, CMB_l, Az_l, El_l, freqs_l = TSources.loadSZMaps(ind.loadSZsource)

            for SZi, SZ_li in zip(SZ.ravel(), SZ_l.ravel()):
                self.assertEqual(SZi, SZ_li)
            
            for CMBi, CMB_li in zip(CMB.ravel(), CMB_l.ravel()):
                self.assertEqual(CMBi, CMB_li)
            
            for Azi, Az_li in zip(Az.ravel(), Az_l.ravel()):
                self.assertEqual(Azi, Az_li)
            
            for Eli, El_li in zip(El.ravel(), El_l.ravel()):
                self.assertEqual(Eli, El_li)
            
            for freqsi, freqs_li in zip(ind.SZsource.get("freqs"), freqs_l):
                self.assertEqual(freqsi*1e9, freqs_li)

    def test_parallelConvolve(self):
        SZ = np.ones((3,3,3))
        k = np.ones(3)
        threadId = 0

        Az, El = np.mgrid[-1:1:3j, -1:1:3j]

        args = (k, SZ, threadId) 
        SZ_test = TSources._parallelConvolve(args, ind.telescope.get("Dtel")/2, Az, El)
        self.assertEqual(SZ.shape, SZ_test.shape)

    @classmethod
    def tearDownClass(self):
        shutil.rmtree(os.path.join(os.getcwd(), "test"))

if __name__ == "__main__":
    import nose2
    nose2.main()

