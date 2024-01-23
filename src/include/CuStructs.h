/*!
 * \file
 * \brief Data structures for receiving data from Python interface.
 **/

#ifndef __CuStructs_h
#define __CuStructs_h

struct CuInstrument;
struct CuTelescope;
struct CuAtmosphere;
struct CuSource;
struct CuSimParams;
struct CuOutput;

// Ctypes structs - for communicating with python
struct CuInstrument {
    float *freqs_filt; /**< Array with frequencies in Hertz.*/
    int nfreqs_filt;    /**< Number of elements in freqs.*/
    int R;              /**< Resolving power of instrument: R = f / df.*/
    float eta_inst;    /**< Instrument efficiency.*/
    float eta_ant;     /**< Antenna efficiency.*/
    float freq_sample; /**< Readout frequency of instrument in Hertz.*/
    float *filterbank; /**< Array with filterbank matrix, flattened.*/
    float delta;       /**< Superconducting bandgap energy in Joules.*/
    float eta_pb;      /**< Pair breaking efficiency of superconductor.*/
};

struct CuTelescope {
    float Ttel;        /**< Telescope temperature in Kelvin.*/
    float Tgnd;        /**< Ground temperature in Kelvin.*/
    float Dtel;        /**< Primary aperture diameter in meters.*/
    int chop_mode;      /**< Chopping mode. 0 is 'none', 1 is 'direct', 2 is 'abba'.*/
    float dAz_chop;    /**< Azimuthal separation between chopping paths.*/
    float freq_chop;   /**< Chopping frequency in Hertz. If < 0, no chopping.*/
    float freq_nod;    /**< Nodding frequency in Hertz.*/
    float *eta_ap;     /**< Array of aperture efficiencies, as function of frequency (set by Instrument). Size is nfreqs of instrument.*/
    float eta_mir;     /**< Mirror reflection efficiency.*/
    float eta_fwd;     /**< Telescope forward efficiency.*/
};

struct CuAtmosphere {
    float Tatm;        /**< Temperature of atmosphere in Kelvin.*/
    float v_wind;      /**< Max windspeed in meters per second.*/
    float h_column;    /**< Reference column height of atmosphere, in meters.*/
    float *x_atm;      /**< Array of floats, representing x-coordinates of atmospheric screen.*/
    float *y_atm;      /**< Array of floats, representing y-coordinates of atmospheric screen.*/
    int nx;             /**< Number of elements in x.*/
    int ny;             /**< Number of elements in y.*/
    float *PWV;        /**< Flat array containing smoothed PWV values at x and y.*/
    float *freqs_atm;  /**< Frequency range over which eta_atm is defined in Hertz.*/
    int nfreqs_atm;     /**< Number of elements in freqs_atm.*/
    float *PWV_atm;    /** < Array containing PWV values in millimeter, over which eta_atm is defined.*/
    int nPWV_atm;       /** < Number of elements in PWV_atm.*/
    float *eta_atm;    /**< Flat array containing all atmospheric transmissions curves over frequency and PWV. Size is nfreqs_atm * nPWV_atm.*/
};

struct CuSource {
    int present;        /**< Whether or not a source is actually present.*/
    float *Az;         /**< Array containing azimuth angles of source on-sky, relative to source center.*/
    int nAz;            /**< Number of elements in Az.*/
    float *El;         /**< Array containing elevation angles on-sky, relative to source center.*/
    int nEl;            /**< Number of elements in El.*/
    float *I_nu;       /**< Flat array of specific intensities, indexed as [i * ni + k * ni * nj + j].
                          Here, i is axis 0 (Az), j axis 1 (El) and k axis 2 (frequency). Note that size of I_nu = nAz * nEl * nfreqs.*/
    float *freqs_src;  /**< Frequencies over which the source is defined, in Hertz.*/
    int nfreqs_src;     /**< Number of frequencies in freq_src.*/
};

struct CuSimParams {
    float t_obs;       /**< Observation time in seconds.*/
    int nTimes;         /**< Total number of time calculations.*/
    int nThreads;       /**< Number of threads to use for computation.*/
    float t0;          /**< Starting time of simulation. If not given will default to 0.*/
    int OFF_empty;   /**< Wether to use no source or interpolate source during OFF chopping.*/
    int use_noise;      /**< Whether to add photon noise. For real life situations, this should be 1. Only set to 0 for debug purposes.*/
};

struct CuOutput {
    float *signal;     /**< Timestream of output signal, time = slow axis, frequency = fast axis.*/
    float *Az;         /**< Timestream of Azimuth angle, in degrees.*/
    float *El;         /**< Timestream of Elevation angle, in degrees.*/
    int *flag;          /**< Timestream of flags specifying chop/nod position. 0 for ON, 1 for OFF-RIGHT, 2 for OFF-LEFT.*/           
    float *t_diag;      /**< Calculation/allocation times, for diagnostics.*/
};

#endif
