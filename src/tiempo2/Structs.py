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
                ("eta_misc", ctypes.c_double),
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
                ("eta_ap_ON", ctypes.POINTER(ctypes.c_double)),
                ("eta_ap_OFF", ctypes.POINTER(ctypes.c_double)),
                ("eta_mir", ctypes.c_double),
                ("eta_fwd", ctypes.c_double),
                ("scantype", ctypes.c_int),
                ("Ax", ctypes.c_double),
                ("Axmin", ctypes.c_double),
                ("Ay", ctypes.c_double),
                ("Aymin", ctypes.c_double),
                ("wx", ctypes.c_double),
                ("wxmin", ctypes.c_double),
                ("wy", ctypes.c_double),
                ("wymin", ctypes.c_double),
                ("phix", ctypes.c_double),
                ("phiy", ctypes.c_double)]

class Atmosphere(ctypes.Structure):
    """!
    Struct representing the simulated atmosphere.
    """

    _fields_ = [("Tatm", ctypes.c_double),
                ("v_wind", ctypes.c_double),
                ("h_column", ctypes.c_double),
                ("x0", ctypes.c_double),
                ("dx", ctypes.c_double),
                ("nx", ctypes.c_int),
                ("y0", ctypes.c_double),
                ("dy", ctypes.c_double),
                ("ny", ctypes.c_int),
                ("f0", ctypes.c_double),
                ("df", ctypes.c_double),
                ("nf", ctypes.c_int),
                ("PWV0", ctypes.c_double),
                ("dPWV", ctypes.c_double),
                ("nPWV", ctypes.c_int),
                ("PWV", ctypes.POINTER(ctypes.c_double)),
                ("eta_atm", ctypes.POINTER(ctypes.c_double))]

class Source(ctypes.Structure):
    """!
    Struct representing simulated astronomical source.
    """

    _fields_ = [("Az0", ctypes.c_double),
                ("dAz", ctypes.c_double),
                ("nAz", ctypes.c_int),
                ("El0", ctypes.c_double),
                ("dEl", ctypes.c_double),
                ("nEl", ctypes.c_int),
                ("f0", ctypes.c_double),
                ("df", ctypes.c_double),
                ("nf", ctypes.c_int),
                ("I_nu", ctypes.POINTER(ctypes.c_double)),
                ("nI_nu", ctypes.c_int)]

class SimParams(ctypes.Structure):
    """!
    Struct representing simulation parameters.
    """

    _fields_ = [("t_obs", ctypes.c_double),
                ("nTimes", ctypes.c_int),
                ("nThreads", ctypes.c_int)]

class Output(ctypes.Structure):
    """!
    Struct used as output container.
    """

    _fields_ = [("signal", ctypes.POINTER(ctypes.c_double)),
                ("Az", ctypes.POINTER(ctypes.c_double)),
                ("El", ctypes.POINTER(ctypes.c_double)),
                ("flag", ctypes.POINTER(ctypes.c_int)),
                ("t_thread", ctypes.c_double)]

# FLOAT
class CuInstrument(ctypes.Structure):
    """!
    Struct representing the simulated instrument.
    """

    _fields_ = [("freqs_filt", ctypes.POINTER(ctypes.c_float)),
                ("nfreqs_filt", ctypes.c_int),
                ("R", ctypes.c_int),
                ("eta_inst", ctypes.c_float),
                ("eta_misc", ctypes.c_float),
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
                ("eta_ap_ON", ctypes.POINTER(ctypes.c_float)),
                ("eta_ap_OFF", ctypes.POINTER(ctypes.c_float)),
                ("eta_mir", ctypes.c_float),
                ("eta_fwd", ctypes.c_float),
                ("scantype", ctypes.c_int),
                ("Ax", ctypes.c_float),
                ("Axmin", ctypes.c_float),
                ("Ay", ctypes.c_float),
                ("Aymin", ctypes.c_float),
                ("wx", ctypes.c_float),
                ("wxmin", ctypes.c_float),
                ("wy", ctypes.c_float),
                ("wymin", ctypes.c_float),
                ("phix", ctypes.c_float),
                ("phiy", ctypes.c_float)]


class CuAtmosphere(ctypes.Structure):
    """!
    Struct representing the simulated atmosphere.
    """

    _fields_ = [("Tatm", ctypes.c_float),
                ("v_wind", ctypes.c_float),
                ("h_column", ctypes.c_float),
                ("x0", ctypes.c_float),
                ("dx", ctypes.c_float),
                ("nx", ctypes.c_int),
                ("y0", ctypes.c_float),
                ("dy", ctypes.c_float),
                ("ny", ctypes.c_int),
                ("f0", ctypes.c_float),
                ("df", ctypes.c_float),
                ("nf", ctypes.c_int),
                ("PWV0", ctypes.c_float),
                ("dPWV", ctypes.c_float),
                ("nPWV", ctypes.c_int),
                ("PWV", ctypes.POINTER(ctypes.c_float)),
                ("eta_atm", ctypes.POINTER(ctypes.c_float))]

class CuSource(ctypes.Structure):
    """!
    Struct representing simulated astronomical source.
    """

    _fields_ = [("Az0", ctypes.c_float),
                ("dAz", ctypes.c_float),
                ("nAz", ctypes.c_int),
                ("El0", ctypes.c_float),
                ("dEl", ctypes.c_float),
                ("nEl", ctypes.c_int),
                ("f0", ctypes.c_float),
                ("df", ctypes.c_float),
                ("nf", ctypes.c_int),
                ("I_nu", ctypes.POINTER(ctypes.c_float)),
                ("nI_nu", ctypes.c_int)]

class CuSimParams(ctypes.Structure):
    """!
    Struct representing simulation parameters.
    """

    _fields_ = [("t_obs", ctypes.c_float),
                ("nTimes", ctypes.c_int),
                ("nThreads", ctypes.c_int)]

class CuOutput(ctypes.Structure):
    """!
    Struct used as output container.
    """

    _fields_ = [("signal", ctypes.POINTER(ctypes.c_float)),
                ("Az", ctypes.POINTER(ctypes.c_float)),
                ("El", ctypes.POINTER(ctypes.c_float)),
                ("flag", ctypes.POINTER(ctypes.c_int)),
                ("t_diag", ctypes.POINTER(ctypes.c_float))]
