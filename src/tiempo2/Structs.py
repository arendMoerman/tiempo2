"""!
@file
Structs that are passed to the C++ backend.
""" 

import ctypes
import numpy as np

class ArrSpec(ctypes.Structure):
    """!
    Struct used for passing array specifications a0, da, na.
    Used when it is not necessary to pass full dictionaries.
    """
    _fields_ = [("start", ctypes.c_double),
                ("step", ctypes.c_double),
                ("num", ctypes.c_int)]

class CuArrSpec(ctypes.Structure):
    """!
    Struct used for passing array specifications a0, da, na.
    Used when it is not necessary to pass full dictionaries.
    """

    _fields_ = [("start", ctypes.c_float),
                ("step", ctypes.c_float),
                ("num", ctypes.c_int)]

class Instrument(ctypes.Structure):
    """!
    Struct representing the simulated instrument.
    """

    _fields_ = [("nf_ch", ctypes.c_int),
                ("f_spec", ArrSpec),
                ("eta_inst", ctypes.c_double),
                ("eta_misc", ctypes.c_double),
                ("f_sample", ctypes.c_double),
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
                ("x_spec", ArrSpec),
                ("y_spec", ArrSpec),
                ("PWV", ctypes.POINTER(ctypes.c_double))]

class Source(ctypes.Structure):
    """!
    Struct representing simulated astronomical source.
    """

    _fields_ = [("Az_spec", ArrSpec),
                ("El_spec", ArrSpec),
                ("I_nu", ctypes.POINTER(ctypes.c_double)),
                ("nI_nu", ctypes.c_int)]

class CalOutput(ctypes.Structure):
    """!
    Struct used as output for a power to temperature conversion database.
    """

    _fields_ = [("power", ctypes.POINTER(ctypes.c_double)),
                ("temperature", ctypes.POINTER(ctypes.c_double))]

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

    _fields_ = [("nf_ch", ctypes.c_int),
                ("f_spec", CuArrSpec),
                ("eta_inst", ctypes.c_float),
                ("eta_misc", ctypes.c_float),
                ("f_sample", ctypes.c_float),
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
                ("x_spec", CuArrSpec),
                ("y_spec", CuArrSpec),
                ("PWV", ctypes.POINTER(ctypes.c_float))]

class CuSource(ctypes.Structure):
    """!
    Struct representing simulated astronomical source.
    """

    _fields_ = [("Az_spec", CuArrSpec),
                ("El_spec", CuArrSpec),
                ("I_nu", ctypes.POINTER(ctypes.c_float)),
                ("nI_nu", ctypes.c_int)]

class CuSimParams(ctypes.Structure):
    """!
    Struct representing simulation parameters.
    """

    _fields_ = [("t_obs", ctypes.c_float),
                ("nTimes", ctypes.c_int),
                ("nThreads", ctypes.c_int)]

class CuCalOutput(ctypes.Structure):
    """!
    Struct used as output for a power to temperature conversion database.
    """

    _fields_ = [("power", ctypes.POINTER(ctypes.c_float)),
                ("temperature", ctypes.POINTER(ctypes.c_float))]

class CuOutput(ctypes.Structure):
    """!
    Struct used as output container.
    """

    _fields_ = [("signal", ctypes.POINTER(ctypes.c_float)),
                ("Az", ctypes.POINTER(ctypes.c_float)),
                ("El", ctypes.POINTER(ctypes.c_float)),
                ("flag", ctypes.POINTER(ctypes.c_int)),
                ("t_diag", ctypes.POINTER(ctypes.c_float))]

