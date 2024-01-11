"""!
@file
Structs that are passed to the C++ backend.
""" 

import ctypes
import numpy as np

class Instrument(ctypes.Structure):
    """!
    Struct representing the simulated instrument.
    """

    _fields_ = [("freqs_filt", ctypes.POINTER(ctypes.c_double)),
                ("nfreqs_filt", ctypes.c_int),
                ("R", ctypes.c_int),
                ("eta_inst", ctypes.c_double),
                ("freq_sample", ctypes.c_double),
                ("filterbank", ctypes.POINTER(ctypes.c_double)),
                ("delta", ctypes.c_double),
                ("eta_pb", ctypes.c_double)]

class Telescope(ctypes.Structure):
    """!
    Struct representing the simulated telescope.
    """

    _fields_ = [("Ttel", ctypes.c_double),
                ("Tgnd", ctypes.c_double),
                ("Dtel", ctypes.c_double),
                ("chop_mode", ctypes.c_int),
                ("dAz_chop", ctypes.c_double),
                ("freq_chop", ctypes.c_double),
                ("freq_nod", ctypes.c_double),
                ("eta_ap", ctypes.POINTER(ctypes.c_double)),
                ("eta_mir", ctypes.c_double),
                ("eta_fwd", ctypes.c_double)]

class Atmosphere(ctypes.Structure):
    """!
    Struct representing the simulated atmosphere.
    """

    _fields_ = [("Tatm", ctypes.c_double),
                ("v_wind", ctypes.c_double),
                ("h_column", ctypes.c_double),
                ("x_atm", ctypes.POINTER(ctypes.c_double)),
                ("y_atm", ctypes.POINTER(ctypes.c_double)),
                ("nx", ctypes.c_int),
                ("ny", ctypes.c_int),
                ("PWV", ctypes.POINTER(ctypes.c_double)),
                ("freqs_atm", ctypes.POINTER(ctypes.c_double)),
                ("nfreqs_atm", ctypes.c_int),
                ("PWV_atm", ctypes.POINTER(ctypes.c_double)),
                ("nPWV_atm", ctypes.c_int),
                ("eta_atm", ctypes.POINTER(ctypes.c_double))]

class Source(ctypes.Structure):
    """!
    Struct representing simulated astronomical source.
    """

    _fields_ = [("present", ctypes.c_int),
                ("Az", ctypes.POINTER(ctypes.c_double)),
                ("nAz", ctypes.c_int),
                ("El", ctypes.POINTER(ctypes.c_double)),
                ("nEl", ctypes.c_int),
                ("I_nu", ctypes.POINTER(ctypes.c_double)),
                ("freqs_src", ctypes.POINTER(ctypes.c_double)),
                ("nfreqs_src", ctypes.c_int)]

class SimParams(ctypes.Structure):
    """!
    Struct representing simulation parameters.
    """

    _fields_ = [("t_obs", ctypes.c_double),
                ("nTimes", ctypes.c_int),
                ("nThreads", ctypes.c_int),
                ("t0", ctypes.c_double),
                ("use_noise", ctypes.c_int)]

class Output(ctypes.Structure):
    """!
    Struct used as output container.
    """

    _fields_ = [("signal", ctypes.POINTER(ctypes.c_double)),
                ("Az", ctypes.POINTER(ctypes.c_double)),
                ("El", ctypes.POINTER(ctypes.c_double)),
                ("flag", ctypes.POINTER(ctypes.c_int))]

# FLOAT
class CuInstrument(ctypes.Structure):
    """!
    Struct representing the simulated instrument.
    """

    _fields_ = [("freqs_filt", ctypes.POINTER(ctypes.c_float)),
                ("nfreqs_filt", ctypes.c_int),
                ("R", ctypes.c_int),
                ("eta_inst", ctypes.c_float),
                ("freq_sample", ctypes.c_float),
                ("filterbank", ctypes.POINTER(ctypes.c_float)),
                ("delta", ctypes.c_float),
                ("eta_pb", ctypes.c_float)]

class CuTelescope(ctypes.Structure):
    """!
    Struct representing the simulated telescope.
    """

    _fields_ = [("Ttel", ctypes.c_float),
                ("Tgnd", ctypes.c_float),
                ("Dtel", ctypes.c_float),
                ("chop_mode", ctypes.c_int),
                ("dAz_chop", ctypes.c_float),
                ("freq_chop", ctypes.c_float),
                ("freq_nod", ctypes.c_float),
                ("eta_ap", ctypes.POINTER(ctypes.c_float)),
                ("eta_mir", ctypes.c_float),
                ("eta_fwd", ctypes.c_float)]

class CuAtmosphere(ctypes.Structure):
    """!
    Struct representing the simulated atmosphere.
    """

    _fields_ = [("Tatm", ctypes.c_float),
                ("v_wind", ctypes.c_float),
                ("h_column", ctypes.c_float),
                ("x_atm", ctypes.POINTER(ctypes.c_float)),
                ("y_atm", ctypes.POINTER(ctypes.c_float)),
                ("nx", ctypes.c_int),
                ("ny", ctypes.c_int),
                ("PWV", ctypes.POINTER(ctypes.c_float)),
                ("freqs_atm", ctypes.POINTER(ctypes.c_float)),
                ("nfreqs_atm", ctypes.c_int),
                ("PWV_atm", ctypes.POINTER(ctypes.c_float)),
                ("nPWV_atm", ctypes.c_int),
                ("eta_atm", ctypes.POINTER(ctypes.c_float))]

class CuSource(ctypes.Structure):
    """!
    Struct representing simulated astronomical source.
    """

    _fields_ = [("present", ctypes.c_int),
                ("Az", ctypes.POINTER(ctypes.c_float)),
                ("nAz", ctypes.c_int),
                ("El", ctypes.POINTER(ctypes.c_float)),
                ("nEl", ctypes.c_int),
                ("I_nu", ctypes.POINTER(ctypes.c_float)),
                ("freqs_src", ctypes.POINTER(ctypes.c_float)),
                ("nfreqs_src", ctypes.c_int)]

class CuSimParams(ctypes.Structure):
    """!
    Struct representing simulation parameters.
    """

    _fields_ = [("t_obs", ctypes.c_float),
                ("nTimes", ctypes.c_int),
                ("nThreads", ctypes.c_int),
                ("t0", ctypes.c_float),
                ("use_noise", ctypes.c_int)]

class CuOutput(ctypes.Structure):
    """!
    Struct used as output container.
    """

    _fields_ = [("signal", ctypes.POINTER(ctypes.c_float)),
                ("Az", ctypes.POINTER(ctypes.c_float)),
                ("El", ctypes.POINTER(ctypes.c_float)),
                ("flag", ctypes.POINTER(ctypes.c_int))]
