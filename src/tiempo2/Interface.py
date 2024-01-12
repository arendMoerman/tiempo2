import math
import os
import time

import numpy as np
import matplotlib.pyplot as pt
import scipy.fft as fft
import scipy.signal as signal

import tiempo2.Filterbank as TFilter
import tiempo2.Sources as TSource
import tiempo2.InputChecker as TCheck
import tiempo2.Atmosphere as TAtm
import tiempo2.BindCPU as TBCPU

import psutil
import logging
from tiempo2.CustomLogger import CustomLogger

logging.getLogger(__name__)

MEMFRAC = 0.5
MEMBUFF = MEMFRAC * psutil.virtual_memory().total

class FieldError(Exception):
    """!
    Field error. Raised when a required field is not specified in an input dictionary. 
    """
    pass

class InitialError(Exception):
    """!
    Initial error. Raised when attempting to run a simulation without submitting all required dictionaries. 
    """
    pass

class Interface(object):
    sourceDict = None
    atmosphereDict = None
    telescopeDict = None
    instrumentDict = None
    observationDict = None

    count = 0
    
    clog_mgr = CustomLogger(os.path.basename(__file__))
    clog = clog_mgr.getCustomLogger()
    
    def setSourceDict(self, sourceDict):
        errlist = TCheck.checkSourceDict(sourceDict)

        if not errlist:
            self.sourceDict = sourceDict
            self.count += 1

        else:
            errstr = f"Errors encountered in source dictionary in fields :{errlist}."
            raise FieldError(errstr)

    def setTelescopeDict(self, telescopeDict):
        errlist = TCheck.checkTelescopeDict(telescopeDict)

        if not errlist:
            self.telescopeDict = telescopeDict
            self.count += 1

        else:
            errstr = f"Errors encountered in telescope dictionary in fields :{errlist}."
            raise FieldError(errstr)
    
    def setInstrumentDict(self, instrumentDict):
        errlist = TCheck.checkInstrumentDict(instrumentDict)

        if not errlist:
            self.instrumentDict = instrumentDict
            self.count += 1

        else:
            errstr = f"Errors encountered in instrument dictionary in fields :{errlist}."
            raise FieldError(errstr)
        
    
    def setAtmosphereDict(self, atmosphereDict):
        errlist = TCheck.checkAtmosphereDict(atmosphereDict)

        if not errlist:
            self.atmosphereDict = atmosphereDict
            self.count += 1

        else:
            errstr = f"Errors encountered in atmosphere dictionary in fields :{errlist}."
            raise FieldError(errstr)
    
    def setObservationDict(self, observationDict):
        errlist = TCheck.checkObservationDict(observationDict)

        if not errlist:
            self.observationDict = observationDict
            self.count += 1

        else:
            errstr = f"Errors encountered in observation dictionary in fields :{errlist}."
            raise FieldError(errstr)

    def runSimulation(self, device="CPU"):
        """!
        Start and run a simulation.

        This is the main routine of TiEMPO2 and should be called after filling all dictionaries.
        Generate or load a source datacube, atmosphere PWV screen, transmission curves.

        Will return exception if not all dictionaries have been set.
        """
        self.clog.info("\033[1;32m*** STARTING TiEMPO2 SIMULATION ***")
        
        start = time.time()

        if self.count < 5:
            raise InitialError
            exit(1)
       
        #### INITIALISING OBSERVATION PARAMETERS ####
        # Calculate number of time evaluations
        self.observationDict["nTimes"] = math.ceil(self.observationDict["t_obs"] * self.instrumentDict["freq_sample"])
        self.observationDict["time_range"] = np.arange(self.observationDict["nTimes"]) / self.instrumentDict["freq_sample"]

        #### INITIALISING INSTRUMENT PARAMETERS ####
        # Generate filterbank
        if isinstance(self.instrumentDict.get("eta_filt"), float):
            self.instrumentDict["eta_filt"] *= np.ones(self.instrumentDict.get("n_freqs"))

        if self.instrumentDict.get("R"):
            self.instrumentDict["freqs_filt"] = self.instrumentDict.get("freq_0") * (1 + 1 / self.instrumentDict.get("R"))**np.arange(self.instrumentDict.get("n_freqs"))
            self.instrumentDict["filterbank"] = TFilter.generateFilterbankFromR(self.instrumentDict, self.sourceDict)
            pt.plot(self.instrumentDict["filterbank"].T)
            pt.show()
        else:
            pass
            # Here we need to call function that reads a filterbank matrix from Louis files.

        #### INITIALISING ASTRONOMICAL SOURCE ####
        # Load or generate source
        if self.sourceDict.get("type") == "SZ":
            SZ, CMB, Az, El = TSource.generateSZMaps(self.sourceDict, telescopeDict=self.telescopeDict)
            I_nu = SZ.T + CMB.T
            freqs = self.sourceDict.get("freqs_src")
        
        elif self.sourceDict.get("type") == "atmosphere":
            I_nu = np.array([0])
            Az = np.array([0])
            El = np.array([0])
            freqs = self.sourceDict.get("freqs_src")

        else:
            SZ, CMB, Az, El, freqs = TSource.loadSZMaps(self.sourceDict)
            I_nu = SZ.T + CMB.T
        
        _sourceDict = {
                "type"      : self.sourceDict.get("type"),
                "Az"        : Az,
                "El"        : El,
                "I_nu"      : I_nu,
                "freqs_src" : freqs*1e9
                }
        
        
        #### INITIALISING ATMOSPHERE PARAMETERS ####
        PWV_atm, nx, ny = TAtm.generateAtmospherePWV(self.atmosphereDict, self.telescopeDict, self.clog)  
        eta_atm, freqs_atm, pwv_curve = TAtm.readAtmTransmissionText()        
      
        #PWV_atm = np.ones((nx, ny)) * 0.1

        fig, ax = pt.subplots(1,1)
        ax.imshow(PWV_atm, aspect='auto')
        pt.show()

        # At t=0, x=y=0 is in middle
        x_atm = (np.arange(0, nx) - ny/2)*self.atmosphereDict.get("dx")
        y_atm = (np.arange(0, ny) - ny/2)*self.atmosphereDict.get("dy")
        
        # Check if atmosphere screen is long enough for given time and windspeed
        length_req = self.atmosphereDict.get("v_wind") * self.observationDict.get("t_obs")

        if length_req > np.max(x_atm):
            t_obs_new = np.floor(np.max(x_atm) / self.atmosphereDict.get("v_wind"))
            self.clog.warning(f"Atmospheric screen too small for windspeed of {self.atmosphereDict.get('v_wind')} m/s and observation time of {self.observationDict.get('t_obs')} s. Reducing observation time to {t_obs_new} s.")

            self.observationDict["t_obs"] = t_obs_new

        _atmDict = {
                "Tatm"      : self.atmosphereDict.get("Tatm"),
                "v_wind"    : self.atmosphereDict.get("v_wind"),
                "h_column"  : self.atmosphereDict.get("h_column"),
                "x_atm"     : x_atm,
                "y_atm"     : y_atm,
                "nx"        : nx,
                "ny"        : ny,
                "PWV"       : PWV_atm,
                "freqs_atm" : freqs_atm * 1e9,
                "nfreqs_atm": freqs_atm.size,
                "PWV_atm"   : pwv_curve,
                "eta_atm"   : eta_atm
                }

        #### INITIALISING TELESCOPE PARAMETERS ####
        if isinstance(self.telescopeDict.get("eta_ap"), float):
            self.telescopeDict["eta_ap"] *= np.ones(freqs.size)

        self.telescopeDict["dAz_chop"] /= 3600

        self.clog.info("Starting observation.")

        if device == "CPU":
            res = TBCPU.runTiEMPO2(self.instrumentDict, self.telescopeDict, 
                        _atmDict, _sourceDict, self.observationDict)
        
        elif device == "GPU":
            res = TBCPU.runTiEMPO2_CUDA(self.instrumentDict, self.telescopeDict, 
                        _atmDict, _sourceDict, self.observationDict)
        
        end = time.time()        
        
        self.clog.info("\033[1;32m*** FINISHED TiEMPO2 SIMULATION ***")
        self.clog.info(f"\033[1;32m*** Elapsed time: {end-start:.2f} seconds ***")
        
        return res, self.observationDict["time_range"], self.instrumentDict["freqs_filt"]
    
    def calcSignalPSD(self, output, axis=0):
        """!
        Calculate signal PSD of a simulation output.

        @param output Output structure obtained from a callc to runSimulation.
        @param axis Axis over which to calculate signal PSD. 0 is time axis (constant channel), 1 is frequency axis (constant timeslice).

        @returns signal_psd Signal PSD.
        @returns freq_psd Frequencies at which the signal PSD is defined, in Hertz.
        """

        if axis == 0:
            dstep = 1 / self.instrumentDict["freq_sample"]
            N = self.observationDict["nTimes"]
        else:
            dstep = self.instrumentDict["freqs_filt"] / self.instrumentDict["R"]
            N = self.instrumentDict["n_freqs"]

        signal_psd = 2 * dstep / N * np.absolute(fft.fftshift(fft.fft(output["signal"], axis=axis)))**2
        freq_signal = fft.fftshift(fft.fftfreq(N, dstep))
        
        return signal_psd, freq_signal

    def avgDirectSubtract(self, output):
        """!
        Calculate a spectrum by averaging over the time-domain signal.
        
        Atmosphere removal is done by direct subtraction.

        @param output Output object obtained from a simulation.

        @returns spectrum Averaged and direct-subtracted spectrum.
        """


        N = output.get("signal")[output.get("flag") == 0]
        A = output.get("signal")[output.get("flag") == 1]
        B = output.get("signal")[output.get("flag") == -1]
        subtract = (np.mean(A, axis=0) + np.mean(B, axis=0)) / 2

        spectrum = np.mean(N, axis=0) - subtract

        return spectrum
