import numpy as np
import matplotlib.pyplot as pt
import os
import struct

_type2fmt = {
        "float"     : [4, 'f'],
        "double"    : [8, 'd'],
        }

NFIELDS = 4

def unpack_output(path, out_t, clog):
    """!
    Process binary files written with results from simulation into output jsons for further processing.
    """
    ntimes = 0
    nfreq = 0

    nchunks = len(os.listdir(path)) // NFIELDS
    if not isinstance(nchunks, int):
        clog.error(f"Cannot load data from {path}: directory corrupt or not valid TiEMPO2 output directory!")
        exit(1)

    for chunk in range(nchunks):
        clog.info(f"Loading chunk {chunk}")
        with open(os.path.join(path, f"{chunk}flag.out"), 'rb') as fh:
            data = fh.read()
            ntimes = len(data)//4
            data_l = struct.unpack(f"@{ntimes}i", data)
        
        with open(os.path.join(path, f"{chunk}signal.out"), 'rb') as fh:
            data = fh.read()
            nsig = len(data)//_type2fmt[out_t][0]
            data_l = struct.unpack(f"@{nsig}{_type2fmt[out_t][1]}", data)
            nfreq = nsig // ntimes
            np.save("test", np.array(data_l).reshape((nfreq, ntimes)).T)
    return np.array(data_l).reshape((nfreq, ntimes)).T
