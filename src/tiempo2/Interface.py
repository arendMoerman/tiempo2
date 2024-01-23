import math
import os
import time
import copy 
import numpy as np
import matplotlib.pyplot as pt
import scipy.fft as fft
import scipy.signal as signal

import tiempo2.Filterbank as TFilter
import tiempo2.Sources as TSource
import tiempo2.InputChecker as TCheck
import tiempo2.Atmosphere as TAtm
import tiempo2.BindCPU as TBCPU
import tiempo2.RemoveATM as TRemove

import psutil
import logging
from tiempo2.CustomLogger import CustomLogger

logging.getLogger(__name__)

MEMFRAC = 0.5
MEMBUFF = MEMFRAC * psutil.virtual_memory().total
T_OBS_BUFF = 0.8

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
    # These are for storing blueprint dicts
    __sourceDict = None
    __atmosphereDict = None
    __telescopeDict = None
    __instrumentDict = None
    __observationDict = None

    # These ones are actually used
    sourceDict = None
    atmosphereDict = None
    telescopeDict = None
    instrumentDict = None
    observationDict = None
    
    initialised = False

    count = 0
    
    clog_mgr = CustomLogger(os.path.basename(__file__))
    clog = clog_mgr.getCustomLogger()
    
    def setSourceDict(self, sourceDict):
        errlist = TCheck.checkSourceDict(sourceDict)

        if not errlist:
            self.__sourceDict = sourceDict
            self.count += 1

        else:
            errstr = f"Errors encountered in source dictionary in fields :{errlist}."
            raise FieldError(errstr)

    def setTelescopeDict(self, telescopeDict):
        errlist = TCheck.checkTelescopeDict(telescopeDict)

        if not errlist:
            self.__telescopeDict = telescopeDict
            self.count += 1

        else:
            errstr = f"Errors encountered in telescope dictionary in fields :{errlist}."
            raise FieldError(errstr)
    
    def setInstrumentDict(self, instrumentDict):
        errlist = TCheck.checkInstrumentDict(instrumentDict)

        if not errlist:
            self.__instrumentDict = instrumentDict
            self.count += 1

        else:
            errstr = f"Errors encountered in instrument dictionary in fields :{errlist}."
            raise FieldError(errstr)
        
    
    def setAtmosphereDict(self, atmosphereDict):
        errlist = TCheck.checkAtmosphereDict(atmosphereDict)

        if not errlist:
            self.__atmosphereDict = atmosphereDict
            self.count += 1

        else:
            errstr = f"Errors encountered in atmosphere dictionary in fields :{errlist}."
            raise FieldError(errstr)
    
    def setObservationDict(self, observationDict):
        errlist = TCheck.checkObservationDict(observationDict)

        if not errlist:
            self.__observationDict = observationDict
            self.count += 1

        else:
            errstr = f"Errors encountered in observation dictionary in fields :{errlist}."
            raise FieldError(errstr)

    def initialise(self):
        if self.count < 5:
            raise InitialError
            exit(1)
        
        self.observationDict = copy.deepcopy(self.__observationDict)
        self.sourceDict = copy.deepcopy(self.__sourceDict)
        self.instrumentDict = copy.deepcopy(self.__instrumentDict)
        self.atmosphereDict = copy.deepcopy(self.__atmosphereDict)
        self.telescopeDict = copy.deepcopy(self.__telescopeDict)
        
        #### INITIALISING OBSERVATION PARAMETERS ####
        # Calculate number of time evaluations
        self.observationDict["nTimes"] = math.ceil(self.observationDict["t_obs"] * self.instrumentDict["freq_sample"])
        self.observationDict["time_range"] = np.arange(self.observationDict["nTimes"]) / self.instrumentDict["freq_sample"]
        
        #### INITIALISING ASTRONOMICAL SOURCE ####
        # Load or generate source
        if self.sourceDict.get("type") == "SZ":
            SZ, Az, El = TSource.generateSZMaps(self.sourceDict, telescopeDict=self.telescopeDict)
            I_nu = SZ.T
            freqs = self.sourceDict.get("freqs_src")
        
        elif self.sourceDict.get("type") == "GalSpec":
            I_nu, Az, El, freqs = TSource.generateGalSpecMaps(self.sourceDict, telescopeDict=self.telescopeDict)
            self.sourceDict["freqs_src"] = freqs 
        elif self.sourceDict.get("type") == "atmosphere":
            I_nu = np.array([0])
            Az = np.array([0])
            El = np.array([0])
            freqs = self.sourceDict.get("freqs_src")

        else:
            SZ, Az, El, freqs = TSource.loadSZMaps(self.sourceDict)
            I_nu = SZ
        
        self.sourceDict["Az"] = Az / 3600
        self.sourceDict["El"] = El / 3600
        self.sourceDict["I_nu"] = I_nu
        self.sourceDict["freqs_src"] = freqs * 1e9

        #### INITIALISING INSTRUMENT PARAMETERS ####
        # Generate filterbank
        
        if isinstance(self.instrumentDict.get("eta_filt"), float):
            self.instrumentDict["eta_filt"] *= np.ones(self.instrumentDict.get("n_freqs"))

        if self.instrumentDict.get("R"):
            self.instrumentDict["freqs_filt"] = self.instrumentDict.get("freq_0") * (1 + 1 / self.instrumentDict.get("R"))**np.arange(self.instrumentDict.get("n_freqs")) * 1e9
            self.instrumentDict["filterbank"] = TFilter.generateFilterbankFromR(self.instrumentDict, self.sourceDict)
        else:
            pass
            # Here we need to call function that reads a filterbank matrix from Louis files.
        
        #### INITIALISING ATMOSPHERE PARAMETERS ####
        
        PWV_atm, nx, ny = TAtm.generateAtmospherePWV(self.atmosphereDict, self.telescopeDict, self.clog)  
        eta_atm, freqs_atm, pwv_curve = TAtm.readAtmTransmissionText()        
      
        # At t=0, x=y=0 is in middle
        x_atm = (np.arange(0, nx) - ny/2)*self.atmosphereDict.get("dx")
        y_atm = (np.arange(0, ny) - ny/2)*self.atmosphereDict.get("dy")
        
        # Check if atmosphere screen is long enough for given time and windspeed
        length_req = self.atmosphereDict.get("v_wind") * self.observationDict.get("t_obs")

        if length_req > np.max(x_atm):
            t_obs_new = np.floor(T_OBS_BUFF * np.max(x_atm) / self.atmosphereDict.get("v_wind"))
            self.clog.warning(f"Atmospheric screen too small for windspeed of {self.atmosphereDict.get('v_wind')} m/s and observation time of {self.observationDict.get('t_obs')} s. Reducing observation time to {t_obs_new} s.")

            self.observationDict["t_obs"] = t_obs_new

        self.atmosphereDict["x_atm"] = x_atm
        self.atmosphereDict["y_atm"] = y_atm
        self.atmosphereDict["nx"] = nx
        self.atmosphereDict["ny"] = ny
        self.atmosphereDict["PWV"] = PWV_atm
        self.atmosphereDict["freqs_atm"] = freqs_atm * 1e9
        self.atmosphereDict["nfreqs_atm"] = freqs_atm.size
        self.atmosphereDict["PWV_atm"] = pwv_curve
        self.atmosphereDict["eta_atm"] = eta_atm

        #### INITIALISING TELESCOPE PARAMETERS ####
        if isinstance(self.telescopeDict.get("eta_ap"), float):
            self.telescopeDict["eta_ap"] *= np.ones(freqs.size)

        self.telescopeDict["dAz_chop"] /= 3600
        self.initialised = True

    def getSourceSignal(self, Az_point, El_point):
        res = TBCPU.getSourceSignal(self.instrumentDict, self.telescopeDict, self.sourceDict, Az_point/3600, El_point/3600)
        return res, self.instrumentDict["freqs_filt"]

    def getNEP(self, PWV_value):
        res = TBCPU.getNEP(self.instrumentDict, self.telescopeDict, self.atmosphereDict, self.sourceDict, PWV_value)
        return res, self.instrumentDict["freqs_filt"]

    def runSimulation(self, device="CPU"):
        """!
        Start and run a simulation.

        This is the main routine of TiEMPO2 and should be called after filling all dictionaries.
        Generate or load a source datacube, atmosphere PWV screen, transmission curves.

        Will return exception if not all dictionaries have been set.
        """

        self.clog.info("\033[1;32m*** STARTING TiEMPO2 SIMULATION ***")
        
        start = time.time()
       
        self.clog.info("Starting observation.")

        if device == "CPU":
            res = TBCPU.runTiEMPO2(self.instrumentDict, self.telescopeDict, 
                        self.atmosphereDict, self.sourceDict, self.observationDict)
        
        elif device == "GPU":
            res = TBCPU.runTiEMPO2_CUDA(self.instrumentDict, self.telescopeDict, 
                        self.atmosphereDict, self.sourceDict, self.observationDict)
        
        end = time.time()        
        
        self.clog.info("\033[1;32m*** FINISHED TiEMPO2 SIMULATION ***")
        
        if self.observationDict.get("get_t_diag") and (device == "GPU"):
            self.clog.info("\033[1;32m*** TIME DIAGNOSTICS ***")
            self.clog.info(f"\033[1;32m*** TOTAL    : {end-start:.2f} seconds ***")
            self.clog.info(f"\033[1;32m*** H2D      : {res.get('t_diag')[0]:.2f} seconds ***")
            self.clog.info(f"\033[1;32m*** KERNEL   : {res.get('t_diag')[1]:.2f} seconds ***")
            self.clog.info(f"\033[1;32m*** D2H      : {res.get('t_diag')[2]:.2f} seconds ***")
        
        elif self.observationDict.get("get_t_diag"):
            self.clog.info("\033[1;32m*** TIME DIAGNOSTICS ***")
            self.clog.info(f"\033[1;32m*** TOTAL    : {end-start:.2f} seconds ***")
            self.clog.info(f"\033[1;32m*** THREAD   : {res.get('t_diag'):.2f} seconds ***")

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

    def rebinSignal(self, output, freqs_old, nbins_add, final_bin=True):
        """!
        Rebin a simulation result into a coarser bin size.
        
        @param final_bin If number of old bins divided by nbins_add is not an integer, wether to rebin final new bin with less bins, or add extra bins to second-to-last bin.
        """

        shape_old  = output.get("signal").shape
        nbins_old = shape_old[1]
        
        if (not isinstance(nbins_add, int)) or (nbins_add == 1):
            self.clog.error(f"Rebin number must be an integer and larger than 1.")
            exit(1)

        if final_bin:
            nbins_new = math.ceil(nbins_old / nbins_add)
        else:
            nbins_new = math.floor(nbins_old / nbins_add)

        signal_new = np.zeros((shape_old[0], nbins_new))
        freqs_new = np.zeros(nbins_new)

        for nbin in range(nbins_new):
            start_bin = nbin * nbins_add
            if nbin == nbins_new - 1:
                signal_new[:,nbin] = np.mean(output.get("signal")[:,start_bin:], axis=1)
                freqs_new[nbin] = np.mean(freqs_old[start_bin:])

            else:
                signal_new[:,nbin] = np.mean(output.get("signal")[:,start_bin:start_bin+nbins_add], axis=1)
                freqs_new[nbin] = np.mean(freqs_old[start_bin:start_bin+nbins_add])

        output_binned = copy.deepcopy(output)
        output_binned["signal"] = signal_new

        return output_binned, freqs_new
    
    def avgDirectSubtract(self, output):
        self.clog.info("Applying averaging and ABBA subtraction")
        return TRemove.avgDirectSubtract(output)
