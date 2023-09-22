import numpy as np

import tiempo2.Sources as TSource
import tiempo2.InputChecker as TCheck

class FieldError(Exception):
    """!
    Field error. Raised when a required field is not specified in an input dictionary. 
    """
    pass

class Interface(object):
    sourceDict = None
    atmosphereDict = None
    telescopeDict = None
    instrumentDict = None
    observationDict = None

    def setSourceDict(self, sourceDict):
        errlist = TCheck.checkSourceDict(sourceDict)

        if not errlist:
            self.sourceDict = sourceDict

        else:
            errstr = f"Errors encountered in source dictionary in fields :{errlist}."
            raise FieldError(errstr)

    def setTelescopeDict(self, telescopeDict):
        errlist = TCheck.checkTelescopeDict(telescopeDict)

        if not errlist:
            self.telescopeDict = telescopeDict

        else:
            errstr = f"Errors encountered in source dictionary in fields :{errlist}."
            raise FieldError(errstr)
    
    def generateSources(self):
     
        #Here, need to convert input dict to a source dict as in the Structs.py file.
        if sourceDict.get("mode") == "make":
            SZ, CMB, Az, El = TSource.generateSZMaps(sourceDict)            
