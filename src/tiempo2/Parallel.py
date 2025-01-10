"""!
@file This file contains private functions for parallel execution of jobs on simulated data.
It is not meant to be used by users, only internally.
"""

from multiprocessing import Pool
import numpy as np
import os
from functools import partial

NFIELDS = 4

def get_num_chunks(path):
    return len(os.listdir(path)) // NFIELDS

def parallel_job(result_path, num_threads, job, args_list):
    """!
    Perform a job on a simulation result in parallel.
    This is used for reduction operations, so the functions passed should return exactly two values: average and standard deviation.
    Note that this method does not calculate number of cores. This is purely user-specified.

    @param result_path Path to simulation results.
    @param num_threads Number of CPU threads to use.
    @param job Function handle of function to apply to data.
    @param args_list List of extra arguments to pass to function.
    @param clog Custom logger for catching log.

    @returns Output of job.
    """

    num_chunks = get_num_chunks(result_path)

    # Check: is num_chunks > num_threads?
    if num_chunks < num_threads:
        num_threads = num_chunks

    print(num_threads)
    chunks = np.array_split(np.arange(num_chunks), num_threads)
    chunks = [chu for chu in chunks]
    
    args = zip(chunks, np.arange(0, num_threads))

    _func = partial(job, result_path=result_path, xargs=args_list)

    pool    = Pool(num_threads)
    
    sum_l, var_l, N_l = zip(*pool.map(_func, args))

    avg     = np.nansum(np.array(sum_l), axis=0) / np.nansum(np.array(N_l))
    var     = np.nansum(np.array(var_l), axis=0) / num_threads**2

    return avg, var

