import math
import os
import shutil
import time
import copy 
import sys
import numpy as np
import matplotlib.pyplot as pt
import scipy.fft as fft
import scipy.signal as signal
from scipy.interpolate import griddata
from scipy.stats import binned_statistic_2d

import tiempo2.Filterbank as TFilter
import tiempo2.Sources as TSource
import tiempo2.InputChecker as TCheck
import tiempo2.Atmosphere as TAtm
import tiempo2.BindCPU as TBCPU
import tiempo2.RemoveATM as TRemove
import tiempo2.Parallel as TParallel
import tiempo2.Analyse as TAnalyse

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
    """!
    Interface for TiEMPO2 simulations

    Attributes:
        __<type>Dict    :   Storage for input dictionaries. These dictionaries are raw copies of user input dicts.
                            When the input dictionary passes the input test, it will be copied to this dictionary.
        <type>Dict      :   Storage for input dictionaries. 
                            These are evaluated versions of the user input, containing references to actual data.
        <type>_set      :   Flags for signifying that a dictionary has been validly set. 
                            Will only be set to 'True' if the dictionary passes the input test.
        clog_mgr        :   Logger manager. Top-level wrapper for actual logger.
                            Manager handles meta information, such as pwd.
        clog            :   Custom logger for communicating information, warnings, and errors to user.
        outPath         :   Path to directory where to store simulation output. Defaults to current working directory.
    """

    __sourceDict        = None
    __atmosphereDict    = None
    __telescopeDict     = None
    __instrumentDict    = None
    __observationDict   = None

    sourceDict          = None
    atmosphereDict      = None
    telescopeDict       = None
    instrumentDict      = None
    observationDict     = None

    source_set          = False
    atmosphere_set      = False
    telescope_set       = False
    instrument_set      = False
    observation_set     = False
    
    initialisedSetup    = False
    initialisedObserve  = False

    clog_mgr = CustomLogger(os.path.basename(__file__))
    clog = clog_mgr.getCustomLogger()

    outPath = os.getcwd() # Output path defaults to cwd    
    
    c = 2.99792458e8

    def __init__(self, verbose=True, outPath=None):
        """!
        Create an interface object for TiEMPO2.

        @param verbose Set verbosity of logger. If False, will not log any output to screen.
            If True, will log simulation information, warnings, and errors to screen. 
            Default is True.
        @param outPath Path to where the interface object will write output. 
            Default is the current working directory.
        """

        self.outPath = outPath if outPath is not None else self.outPath
        if not verbose:
            self.clog.setLevel(logging.CRITICAL)
    
    def setLoggingVerbosity(self, verbose=True):
        """!
        Set or change the verbosity of the TiEMPO2 logger.
        
        @param verbose Give verbose output. Default is True.

        @ingroup settersgetters
        """

        if not verbose:
            self.clog.setLevel(logging.CRITICAL)
        
        else:
            self.clog.setLevel(logging.INFO)

    def setOutPath(self, outPath):
        """!
        Set or change path to output directory.

        @param outPath Path to output directory.

        @ingroup settersgetters
        """

        self.outPath = outPath
    

    def setSourceDict(self, sourceDict):
        """!
        Set a source dictionary.

        @param sourceDict A dictionary specifying the source to simulate.

        @ingroup inputdicts
        """

        errlist = TCheck.checkSourceDict(sourceDict)

        if not errlist:
            self.__sourceDict = sourceDict
            self.source_set = True

        else:
            errstr = f"Errors encountered in source dictionary in fields :{errlist}."
            raise FieldError(errstr)

    def setTelescopeDict(self, telescopeDict):
        """!
        Set a telescope dictionary.

        @param telescopeDict A dictionary specifying the telescope for the simulation.

        @ingroup inputdicts
        """

        errlist = TCheck.checkTelescopeDict(telescopeDict)

        if not errlist:
            self.__telescopeDict = telescopeDict
            self.telescope_set = True

        else:
            errstr = f"Errors encountered in telescope dictionary in fields :{errlist}."
            raise FieldError(errstr)
    
    def setInstrumentDict(self, instrumentDict):
        """!
        Set an instrument dictionary.

        @param instrumentDict A dictionary specifying the instrument to simulate.

        @ingroup inputdicts
        """

        errlist = TCheck.checkInstrumentDict(instrumentDict)

        if not errlist:
            self.__instrumentDict = instrumentDict
            self.instrument_set = True

        else:
            errstr = f"Errors encountered in instrument dictionary in fields :{errlist}."
            raise FieldError(errstr)
        
    
    def setAtmosphereDict(self, atmosphereDict):
        """!
        Set an atmosphere dictionary.
        Note that this dictionary only specifies some parameters, such as windspeed, scale height, etc. Actual atmosphere screens need to be generated using ARIS. 
        Also, the directory containing the ARIS screens, together with the actual name of the screen files, is specified here.

        @param atmosphereDict A dictionary specifying the atmosphere to simulate.

        @ingroup inputdicts
        """

        errlist = TCheck.checkAtmosphereDict(atmosphereDict)

        if not errlist:
            self.__atmosphereDict = atmosphereDict
            self.atmosphere_set = True

        else:
            errstr = f"Errors encountered in atmosphere dictionary in fields :{errlist}."
            raise FieldError(errstr)

    def prepareAtmosphere(self):
        """!
        Prepare the ARIS screen chunks for usage in TiEMPO2 simulations.
        A telescope and atmosphere dictionary need to be set, otherwise an error is raised.
        
        This function loads the ARIS screens, convolves them with an apodized Gaussian, the size of the telescope aperture and an edge taper of -10 dB (this will be adjustable in future updates). 
        Then, the prepared ARIS screens are saved in the same directory, under the 'prepd' subdirectory. 
        The original names are stripped, and only the screen index is used as name. The extension type is changed from .dat into .datp.

        For each simulation setup (disregarding source) and ARIS screen set, this function only needs to be run once. 
        For subsequent simulations the prepared screens can be reused.
        
        @ingroup initialise
        """

        if not (self.telescope_set and
                self.atmosphere_set):
            
            self.clog.error("a telescope and atmosphere dictionary MUST be set before calling prepareAtmosphere()")
            raise InitialError
            sys.exit()

        TAtm.prepAtmospherePWV(self.__atmosphereDict, self.__telescopeDict, self.clog)

    # NUMBER PARAMETER IS TEMPORARY
    def initSetup(self, use_ARIS=True, number=1, start=0):
        """!
        Initialise a TiEMPO2 setup. 

        The idea is that a TiEMPO2 setup can be done sequentially, the first step being a setup of the terrestrial part.
        This means that the first setup should consist of everything on Earth: instrument, telescope, and atmosphere.
        Therefore, this is the first initialisation that should be performed.

        This function will check wether or not the instrument, telescope, and atmosphere are set.
        If not, an InitialError will be raised.

        @param use_ARIS Whether to load an ARIS screen or not. 
            Some functions of TiEMPO2 do not require an ARIS screen to be loaded. 
            Setting this parameter to False could reduce total memory footprint in these cases.
            Default is True (load the ARIS screen).
        @param number Number of ARIS chunks to concatenate and load into memory.
        @param start ARIS chunk to start with. 
        
        @ingroup initialise
        """

        if not (self.instrument_set and 
                self.telescope_set and
                self.atmosphere_set):
            
            self.clog.error("an instrument, telescope, and atmosphere dictionary MUST be set before calling initSetup()")
            raise InitialError
            sys.exit()
       
        # Deepcopy raw dicts into proper input dicts, so to not have to interact with raw templates.
        self.instrumentDict = copy.deepcopy(self.__instrumentDict)
        self.telescopeDict = copy.deepcopy(self.__telescopeDict)
        self.atmosphereDict = copy.deepcopy(self.__atmosphereDict)

        #### SI CONVERSION ####
        f0_ch = self.instrumentDict["f0_ch"] * 1e9
        f0_src = self.instrumentDict["f0_src"] * 1e9
        f1_src = self.instrumentDict["f1_src"] * 1e9
        nf_src = self.instrumentDict["nf_src"]

        f_src = np.linspace(f0_src, f1_src, nf_src)

        self.instrumentDict["f_src"] = f_src 

        # Generate range of channel frequencies. Only for Python-side usage
        
        #### INITIALISING INSTRUMENT PARAMETERS ####
        # Generate filterbank
        if self.instrumentDict.get("R"):
            if isinstance(f0_ch, float) or isinstance(f0_ch, int):
                idx_ch_arr = np.arange(self.instrumentDict["nf_ch"])
                self.instrumentDict["f_ch_arr"] = f0_ch * (1 + 1 / self.instrumentDict["R"])**idx_ch_arr
            else:
                self.instrumentDict["f_ch_arr"] = f0_ch
            self.instrumentDict["filterbank"] = TFilter.generateFilterbankFromR(self.instrumentDict)
        else:
            self.instrumentDict["f_ch_arr"] = np.linspace(100,200,nf_src) 
            # Here we need to call function that reads a filterbank matrix from Louis files.

        #pt.plot(self.instrumentDict["filterbank"].T)
        #pt.show()
        
        #### INITIALISING ATMOSPHERE PARAMETERS ####

        #### INITIALISING TELESCOPE PARAMETERS ####
        if isinstance(self.telescopeDict.get("eta_ap_ON"), float):
            self.telescopeDict["eta_ap_ON"] *= np.ones(nf_src)
        
        if isinstance(self.telescopeDict.get("eta_ap_OFF"), float):
            self.telescopeDict["eta_ap_OFF"] *= np.ones(nf_src)

        if self.telescopeDict["s_rms"] is not None:
            self.telescopeDict["s_rms"] *= 1e-6 # Convert um to m

            eta_surf = np.exp(-(4 * np.pi * self.telescopeDict["s_rms"] * f_src / self.c)**2)

            self.telescopeDict["eta_ap_ON"] *= eta_surf 
            self.telescopeDict["eta_ap_OFF"] *= eta_surf 
        
        self.telescopeDict["dAz_chop"] /= 3600
        self.telescopeDict["Ax"] /= 3600
        self.telescopeDict["Axmin"] /= 3600
        self.telescopeDict["Ay"] /= 3600
        self.telescopeDict["Aymin"] /= 3600
        
        #### END INITIALISATION ####
        self.initialisedSetup = True

    def initSource(self, pointing=None):
        """!
        Initialise a TiEMPO2 source. 

        The idea is that a TiEMPO2 setup can be done sequentially, the first step being a setup of the terrestrial part.
        This means that the second setup should consist of everything in space: the astronomical source.
        Therefore, this is the second initialisation that should be performed.

        @param pointing Pointing center of telescope, w.r.t. center of astronomical source, in arcseconds.
        
        @ingroup initialise
        """

        pointing = np.zeros(2) if pointing is None else pointing
        self.sourceDict = copy.deepcopy(self.__sourceDict)
        
        trace_src = None 

        #### INITIALISING ASTRONOMICAL SOURCE ####
        # Load or generate source
        if self.sourceDict.get("type") == "SZ":
            # UPDATE FOLLOWING FOR ARBITRARY Az-El CENTERS!!!!!!!!!
            if self.telescopeDict["scantype"] == 0 and self.telescopeDict["chop_mode"] == 2:
                trace_src = np.zeros((2,3))
                trace_src[0,:] = np.array([-1, 0, 1]) * self.telescopeDict["dAz_chop"]*3600 + pointing[0]
                trace_src[1,:] = np.array([0, 0, 0]) + pointing[1]
            
            elif self.telescopeDict["scantype"] == 0 and self.telescopeDict["chop_mode"] == 0:
                trace_src = np.zeros((2,3))
                trace_src[0,:] = np.array([pointing[0]])
                trace_src[1,:] = np.array([pointing[1]])
            
            I_nu, Az, El = TSource.generateSZMaps(self.sourceDict, self.instrumentDict, self.clog, telescopeDict=self.telescopeDict, trace_src=trace_src)

        elif self.sourceDict.get("type") == "GalSpec":
            I_nu, Az, El, freqs = TSource.generateGalSpecMaps(self.sourceDict, self.instrumentDict, self.telescopeDict)
        
        elif self.sourceDict.get("type") == "background":
            if self.telescopeDict["scantype"] == 0 and self.telescopeDict["chop_mode"] == 2:
                I_nu = np.zeros(3 * self.instrumentDict["nf_src"])
            Az = np.zeros(self.instrumentDict["nf_src"])
            El = np.zeros(self.instrumentDict["nf_src"])
        
        elif self.sourceDict.get("type") == "blackbody":
            I_nu = np.zeros(self.instrumentDict["nf_src"])
            Az = np.zeros(self.instrumentDict["nf_src"])
            El = np.zeros(self.instrumentDict["nf_src"])

        else:
            I_nu, Az, El, freqs = TSource.loadSZMaps(self.sourceDict)
        
        self.sourceDict["Az_src"] = Az / 3600
        self.sourceDict["El_src"] = El / 3600
        self.sourceDict["I_nu"] = I_nu

    def getSourceSignal(self, Az_point, El_point, PWV_value=-1, ON=True):
        """!
        Get astronomical signal without atmospheric noise, but with all efficiencies, atmospheric transmission and filterbank.

        @param Az_point Azimuth point on-sky in arcsec where source should be evaluated.
        @param El_point Elevation point on-sky in arcsec where source should be evaluated.
        @param PWV_value PWV value for atmospheric transmission. Defaults to -1 (no atmosphere).
        @param ON Whether to evaluate source in ON path (default), or OFF. Makes a difference when different eta_ap have been defined for each path.

        @returns Transmitted signal (SI) and its frequency range (Hz).

        @ingroup auxilliarymethods
        """
        
        trace_src = np.array([[Az_point], [El_point]]) 
        
        SZ, Az, El = TSource.generateSZMaps(self.sourceDict, self.instrumentDict, self.clog, telescopeDict=self.telescopeDict, trace_src=trace_src)
        SZ = np.squeeze(SZ) 
        res = TBCPU.getSourceSignal(self.instrumentDict, self.telescopeDict, self.atmosphereDict, SZ, PWV_value, ON)

        return res, self.instrumentDict["f_ch_arr"]
    
    def getChopperCalibration(self, Tcal):
        """!
        Get power emitted by a blackbody source in front of cryostat.
        For calibration purposes.
 
        @ingroup auxilliarymethods       
        """

        res = TBCPU.getChopperCalibration(self.instrumentDict, Tcal)
        return res
    
    def getEtaAtm(self, PWV_value):
        """!
        Get atmosphere transmission at source frequencies given a PWV.

        @param PWV_value PWV value for atmospheric transmission. 

        @returns Atmospheric transmission and its frequency range (Hz).
 
        @ingroup auxilliarymethods       
        """

        res = TBCPU.getEtaAtm(self.instrumentDict, PWV_value)
        return res, self.instrumentDict["f_src"]

    def getNEP(self, PWV_value):
        """!
        Calculate Noise Equivalent Power (NEP) from the atmosphere.

        @param PWV_value PWV value at which to calculate atmospheric NEP.

        @returns NEP (SI) and its frequency range (Hz).
 
        @ingroup auxilliarymethods       
        """

        res = TBCPU.getNEP(self.instrumentDict, self.telescopeDict, self.atmosphereDict, PWV_value)
        return res, self.instrumentDict["f_ch_arr"]

    def runSimulation(self, t_obs, device="CPU", nThreads=None, verbosity=1, outpath="./out/", overwrite=False):
        """!
        Run a TiEMPO2 simulation.

        This is the main routine of TiEMPO2 and should be called after filling all dictionaries and running the self.initialise() method.
        
        @param t_obs Total observation time in seconds.
        @param device Whether to run on CPU (device='CPU') or GPU (device='GPU'). Default is 'CPU'.
        @param nThreads Number of threads to use. Default is 1. Only relevant when device='CPU'.
        @param verbosity Level of verbosity for simulation.
            0           : no extra output w.r.t. logger.
            1 (default) : show execution times of important routines.
            2           : show execution times of important routines and memory transactions.
        @param outpath Directory to store output in.
        @param overwrite Whether to overwrite existing output directories. 
            If False (default), a prompt will appear to either overwrite or specify new directory.
        """

        nThreads = 1 if nThreads is None else nThreads

        if not self.initialisedSetup:
            self.clog.error("initSetup() MUST be called before running a simulation!")
            raise InitialError
            sys.exit()
        
        #### INITIALISING OBSERVATION PARAMETERS ####
        # Calculate number of time evaluations
        # Note that, in order to simplify TLS noise calculations, we make nTimes even
        nTimes = math.ceil(t_obs * self.instrumentDict["f_sample"])

        if nTimes % 2 == 1:
            nTimes -= 1
        
        t_range = 1 / self.instrumentDict["f_sample"] * np.arange(nTimes)

        outpath_succes = False
        while not outpath_succes:
            if not os.path.exists(outpath):
                os.makedirs(outpath)
                outpath_succes = True

            elif overwrite:
                shutil.rmtree(outpath)
                os.makedirs(outpath)
                outpath_succes = True

            elif not overwrite:
                self.clog.warning(f"Output path {outpath} exists! Overwrite or change path?")
                choice = input("\033[93mOverwrite (y/n)? > ").lower()
                if choice == "y":
                    shutil.rmtree(outpath)
                else:
                    outpath = input("\033[93mSpecify new output path: > ")


        self.clog.info("\033[1;32m*** STARTING TiEMPO2 SIMULATION ***")
        
        start = time.time()

        if device == "CPU":
            res = TBCPU.runTiEMPO2(self.instrumentDict, self.telescopeDict, 
                        self.atmosphereDict, self.sourceDict, nTimes, nThreads)
        
        elif device == "GPU":
            TBCPU.runTiEMPO2_CUDA(self.instrumentDict, self.telescopeDict, 
                        self.atmosphereDict, self.sourceDict, nTimes, outpath)
        
        end = time.time()        
        
        self.clog.info("\033[1;32m*** FINISHED TiEMPO2 SIMULATION ***")

    def calcW2K(self, nPWV, nThreads=None, verbosity=1):
        """!
        Calculate a Watt-to-Kelvin calibration table..

        @param nPWV Number of PWV points to consider.
        @param nThreads Number of threads to use. Default is 1.
        @param verbosity Level of verbosity for simulation.
            0           : no extra output w.r.t. logger.
            1 (default) : show execution times of important routines.
            2           : show execution times of important routines and memory transactions.

        @ingroup auxilliarymethods
        """

        nThreads = 1 if nThreads is None else nThreads
        
        if not self.initialisedSetup:
            self.clog.error("initSetup() MUST be called before running a W2K calibration!")
            raise InitialError
            sys.exit()

        self.clog.info("\033[1;32m*** STARTING TiEMPO2 W2K CALIBRATION ***")
        
        start = time.time()

        res = TBCPU.calcW2K(self.instrumentDict, self.telescopeDict, 
                            self.atmosphereDict, nPWV, nThreads)
        
        a = np.zeros(self.instrumentDict["nf_ch"])
        b = np.zeros(self.instrumentDict["nf_ch"])
        for k in range(self.instrumentDict["nf_ch"]):
            _a, _b = np.polyfit(res["power"][:,k], res["temperature"][:,k], 1)
            a[k] = _a
            b[k] = _b

        res["a"] = a
        res["b"] = b

        end = time.time()        
        
        self.clog.info("\033[1;32m*** FINISHED TiEMPO2 W2K CALIBRATION ***")

        return res, self.instrumentDict["f_ch_arr"]
    
    def Watt2Kelvin(self, output, w2k):
        """!
        Convert the signal in output from a Watt to a Kelvin temperature scale.

        This function updates the content of the 'signal' field of the output dictionary, 
        so if you want to keep the signal in Watts as well, make sure to (deep)copy it to a different array first.

        @param output An output dictionary containing the signal.
        @param w2k A Watt-to-Kelvin (w2k) dictionary.

        @ingroup auxilliarymethods
        """

        n_ch = output["signal"].shape[0]
        
        output["signal"] = w2k["a"] * output["signal"] + w2k["b"]

    def calcSignalPSD(self, outpath, num_threads=1, nperseg=256):
        """!
        Calculate signal PSD of a simulation output.

        @param output Output structure obtained from a callc to runSimulation.
        @param x Array over which output signal is defined.
        @param axis Axis over which to calculate signal PSD. 0 is time axis (constant channel), 1 is frequency axis (constant timeslice).

        @returns signal_psd Signal PSD.
        @returns freq_psd Frequencies at which the signal PSD is defined, in Hertz.
 
        @ingroup auxilliarymethods       
        """

        if not self.instrument_set:
            
            self.clog.error("the instrument with which the data was simulated must be set in the interface!")
            raise InitialError
            sys.exit()
        
        self.clog.info("Calculating signal PSD.")
        
        Pxx, f_Pxx = TParallel.parallel_job(outpath, 
                                      num_threads = num_threads, 
                                      job = TAnalyse.calcSignalPSD, 
                                      args_list = [self.instrumentDict["f_sample"], nperseg])
        
        arr_sizes = np.array([x.size for x in f_Pxx])

        tot_Pxx = np.nanmean(Pxx[arr_sizes == np.max(arr_sizes)], axis=0)
        tot_f_Pxx = np.nanmean(f_Pxx, axis=0)

        return tot_Pxx, tot_f_Pxx

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
    
    def avgDirectSubtract(self, outpath, resolution=0, num_threads=1, var_method="PSD"):
        """!
        Apply full time-averaging and direct atmospheric subtraction.

        @param outpath Path to where output from a simulation is saved.
        @param resolution How far to average and subtract. The following options are available:
            0: Reduce signal to a single spectrum, fully averaged and subtracted according to the chopping/nodding scheme.
            1: Reduce signal by averaging over ON-OFF chop positions.
        @param num_threads Number of threads to use. If this number is more than necessary, it will automatically scale down.

        @returns red_signal Reduced signal.
        @returns red_Az Reduced Azimuth array. Only returned if resolution == 1.
        @returns red_El Reduced Elevation array. Only returned if resolution == 1.
        """
        
        if not self.instrument_set:
            self.clog.error("the instrument with which the data was simulated must be set in the interface!")
            raise InitialError
            sys.exit()
        
        if resolution == 0:
            self.clog.info("Applying time-averaging and direct subtraction.")
            
            avg_l, var_l, N_l = TParallel.parallel_job(outpath, 
                                          num_threads = num_threads, 
                                          job = TAnalyse.avgDirectSubtract, 
                                          args_list = [self.instrumentDict["f_sample"], var_method])
            
            avg = np.nansum((N_l - 1) * avg_l.T, axis=1) / (np.nansum(N_l) - len(N_l))
            var = np.nansum((N_l - 1) * var_l.T, axis=1) / (np.nansum(N_l) - len(N_l))**2

            return avg, var

        elif resolution == 1:
            self.clog.info("Averaging and subtracting over ON-OFF pairs.")
            red_signal, red_Az, red_El = TRemove.avgDirectSubtract_chop(output)
            
            return red_signal, red_Az, red_El

    def getExposureTime(self, red_Az, red_El, nAz_grid, nEl_grid):
        points = np.array([red_Az, red_El])

        min_Az = np.min(red_Az)
        max_Az = np.max(red_Az)
        min_El = np.min(red_El)
        max_El = np.max(red_El)

        grid_Az, grid_El = np.mgrid[min_Az:max_Az:nAz_grid*1j, min_El:max_El:nEl_grid*1j]
        statObj = binned_statistic_2d(red_Az, red_El, np.ones(red_Az.size), bins=[nAz_grid, nEl_grid], statistic='count') 
        exposure_time = statObj.statistic / self.instrumentDict.get("f_sample")
    
        return exposure_time, grid_Az, grid_El

    def regridRedSignal(self, red_signal, red_Az, red_El, nAz_grid, nEl_grid, stack=False, idx=None):
        points = np.array([red_Az, red_El])

        min_Az = np.min(red_Az)
        max_Az = np.max(red_Az)
        min_El = np.min(red_El)
        max_El = np.max(red_El)

        grid_Az, grid_El = np.mgrid[min_Az:max_Az:nAz_grid*1j, min_El:max_El:nEl_grid*1j]

        if idx is None:
            if stack:
                grid_signal = binned_statistic_2d(red_Az, red_El, np.mean(red_signal, axis=1), bins=[nAz_grid, nEl_grid]).statistic 

            else:
                grid_signal = np.zeros((nAz_grid, nEl_grid, red_signal.shape[1]))
                for i in range(red_signal.shape[1]):
                    grid_signal[:,:,i] = binned_statistic_2d(red_Az, red_El, red_signal[:,i], bins=[nAz_grid, nEl_grid]).statistic 

        else:
            grid_signal = binned_statistic_2d(red_Az, red_El, red_signal[:,idx], bins=[nAz_grid, nEl_grid]).statistic 
        
        return grid_signal, grid_Az, grid_El














