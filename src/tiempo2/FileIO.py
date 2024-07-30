import numpy as np
import matplotlib.pyplot as pt
import os
import struct

def unpack_output(path, out_t):
    ntimes = 0
    with open(os.path.join(path, "0flag.out"), 'rb') as fh:
        data = fh.read()
        ntimes = len(data)//4
        data_l = struct.unpack(f"<{ntimes}i", data)
    
    with open(os.path.join(path, "0signal.out"), 'rb') as fh:
        data = fh.read()
        n = len(data)//4
        data_l = struct.unpack(f"<{len(data)//4}f", data)
        ny = n // ntimes
        np.save("test", np.array(data_l).reshape((ny, ntimes)).T)
    return np.array(data_l).reshape((ny, ntimes)).T
