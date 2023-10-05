import numpy as np
import matplotlib.pyplot as pt

import tiempo2.Filterbank as TFilter
import tiempo2.Sources as TSource
import tiempo2.InputChecker as TCheck
import tiempo2.Atmosphere as TAtm
import tiempo2.BindCPU as TBCPU

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

    def runSimulation(self):
        """!
        Start and run a simulation.

        This is the main routine of TiEMPO2 and should be called after filling all dictionaries.
        Generate or load a source datacube, atmosphere PWV screen, transmission curves.

        Will return exception if not all dictionaries have been set.
        """

        if self.count < 5:
            raise InitialError
            exit(1)

        # Generate filterbank
        if self.instrumentDict.get("R"):
            self.instrumentDict["filterbank"] = TFilter.generateFilterbankFromR(self.instrumentDict, self.sourceDict)

        else:
            pass
            # Here we need to call function that reads a filterbank matrix

        # Load or generate source
        if self.sourceDict.get("type") == "SZ":
            SZ, CMB, Az, El = TSource.generateSZMaps(self.sourceDict, telescopeDict=self.telescopeDict)
            freqs = self.sourceDict.get("freqs_src")
        
        else:
            SZ, CMB, Az, El, freqs = TSource.loadSZMaps(self.sourceDict)
        
        _sourceDict = {
                "Az"        : Az,
                "El"        : El,
                "I_nu"      : SZ.T + CMB.T,
                "freqs_src" : freqs*1e9
                }
        
        PWV_atm, nx, ny = TAtm.generateAtmospherePWV(self.atmosphereDict, self.telescopeDict)  
        eta_atm, freqs_atm, pwv_curve = TAtm.readAtmTransmission()        
        
        # At t=0, x=y=0 is in middle
        x_atm = (np.arange(0, nx) - ny/2)*self.atmosphereDict.get("dx")
        y_atm = (np.arange(0, ny) - ny/2)*self.atmosphereDict.get("dx")
        
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

        if isinstance(self.telescopeDict.get("eta_ap"), float):
            self.telescopeDict["eta_ap"] *= np.ones(freqs.size)

        self.telescopeDict["dAz_chop"] /= 3600

        res = TBCPU.runTiEMPO2(self.instrumentDict, self.telescopeDict, 
                        _atmDict, _sourceDict, self.observationDict)
        
        return res

