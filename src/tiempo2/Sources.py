"""!
@file
This file contains the possible source signals for a TiEMPO2 simulation.

Currently, can only use MockSZ.
"""

from dataclasses import dataclass, field

import MockSZ.Models as MModel
import MockSZ.Conversions as MConv

import numpy as np

@dataclass
class MockSZ:
    """!
    Class to represent the Sunyaev-Zeldovich effect on-sky.

    In essence, this class wraps main functionality of MockSZ. 
    It generates an isothermal-beta model and simulates the tSZ and kSZ effects, as function of angle on-sky.
    As output, it generates two datacubes (one for tSZ and kSZ) with the following layout: (dAzimuth, dElevation, frequency).
    Here, dAzimuth is the azimuthal angle from cluster center, dElevation elevation angle from cluster center and frequency the frequency of radiation.
    
    @param Te Electron temperature of cluster plasma, in KeV.
    @param ne0 Central electron number density in 1 / cm**3.
    @param rc Cluster core radius in kpc.
    @param beta Structure constant of cluster in isothermal-beta model.
    @param Da Cluster angular diameter distance in Mpc.
    @param v_pec Peculiar cluster velocity in km /s.
    """

    Te      : float
    ne0     : float
    rc      : float
    beta    : float
    Da      : float
    v_pec   : float = field(default = 0.)

    def generateMaps(self, freq_Hz, lims_Az, lims_El, nAz, nEl, ret_unit="MJy"):
        """!
        Generate SZ maps on-sky over a rage of frequencies.

        Note that the lower and upper azimuth and elevation angles are with respect to cluster center, i.e., centered around 0.

        @param freq_Hz Range of frequencies over which to generate maps, in Hertz.
        @param lims_Az Lower and upper limits on azimuth range in degrees.
        @param lims_El Lower and upper limits on elevation range in degrees.
        @param nAz Number of cells along azimuth axis.
        @param nEl Number of cells along elevation axis.
        @param ret_unit Return maps in this unit. Choose between standard W m**-2 Hz**-1 sr**-1 ('SI'), MJy sr**-1 ('MJy') or mK ('mK').

        @returns tSZ Datacube containing tSZ intensity sky maps across the frequency range.
        @returns kSZ Datacube containing kSZ intensity sky maps across the frequency range.
        """

        isob = MModel.IsoBetaModel(self.Te, self.ne0, self.rc, self.beta, self.Da, v_pec=self.v_pec)

        Az, El = np.mgrid[lims_Az[0]:lims_Az[1]:nAz * 1j, lims_El[0]:lims_El[1]:nEl * 1j]
        theta = np.sqrt(Az**2 + El**2)
        
        if ret_unit == "MJy":
            tSZ = MConv.SI_JySr(isob.tSZMap(theta, freq_Hz)) * 1e-6 # MJy / sr
            kSZ = MConv.SI_JySr(isob.kSZMap(theta, freq_Hz)) * 1e-6 # MJy / sr

        elif ret_unit == "SI":
            tSZ = isob.tSZMap(theta, freq_Hz)
            kSZ = isob.kSZMap(theta, freq_Hz)

        elif ret_unit == "mK":
            tSZ = MConv.SI_Temp(isob.tSZMap(theta, freq_Hz), freq_Hz) * 1e3 # mK
            kSZ = MConv.SI_Temp(isob.kSZMap(theta, freq_Hz), freq_Hz) * 1e3 # mK

        return tSZ, kSZ
