/*!
 * \file
 * \brief Data structures for receiving data from Python interface.
 **/

#ifndef __Structs_h
#define __Structs_h

struct Instrument;
struct Telescope;
struct Atmosphere;
struct Source;
struct SimParams;
struct Output;

// Ctypes structs - for communicating with python
struct Instrument {
    double *freqs_filt; /**< Array with frequencies in Hertz.*/
    int nfreqs_filt;    /**< Number of elements in freqs.*/
    int R;              /**< Resolving power of instrument: R = f / df.*/
    double eta_inst;    /**< Instrument efficiency.*/
    double eta_ant;    /**< Antenna radiation efficiency.*/
    double freq_sample; /**< Readout frequency of instrument in Hertz.*/
    double *filterbank; /**< Array with filterbank matrix, flattened.*/
    double delta;       /**< Superconducting bandgap energy in Joules.*/
    double eta_pb;      /**< Pair breaking efficiency of superconductor.*/
};

struct Telescope {
    double Ttel;        /**< Telescope temperature in Kelvin.*/
    double Tgnd;        /**< Ground temperature in Kelvin.*/
    double Dtel;        /**< Primary aperture diameter in meters.*/
    int chop_mode;      /**< Chopping mode. 0 is 'none', 1 is 'direct', 2 is 'abba'.*/
    double dAz_chop;    /**< Azimuthal separation between chopping paths.*/
    double freq_chop;   /**< Chopping frequency in Hertz. If < 0, no chopping.*/
    double freq_nod;    /**< Nodding frequency in Hertz.*/
    double *eta_ap;     /**< Array of aperture efficiencies, as function of frequency (set by Instrument). Size is nfreqs of instrument.*/
    double eta_mir;     /**< Mirror reflection efficiency.*/
    double eta_fwd;     /**< Telescope forward efficiency.*/
};

struct Atmosphere {
    double Tatm;        /**< Temperature of atmosphere in Kelvin.*/
    double v_wind;      /**< Max windspeed in meters per second.*/
    double h_column;    /**< Reference column height of atmosphere, in meters.*/
    double *x_atm;      /**< Array of doubles, representing x-coordinates of atmospheric screen.*/
    double *y_atm;      /**< Array of doubles, representing y-coordinates of atmospheric screen.*/
    int nx;             /**< Number of elements in x.*/
    int ny;             /**< Number of elements in y.*/
    double *PWV;        /**< Flat array containing smoothed PWV values at x and y.*/
    double *freqs_atm;  /**< Frequency range over which eta_atm is defined in Hertz.*/
    int nfreqs_atm;     /**< Number of elements in freqs_atm.*/
    double *PWV_atm;    /** < Array containing PWV values in millimeter, over which eta_atm is defined.*/
    int nPWV_atm;       /** < Number of elements in PWV_atm.*/
    double *eta_atm;    /**< Flat array containing all atmospheric transmissions curves over frequency and PWV. Size is nfreqs_atm * nPWV_atm.*/
};

struct Source {
    int present;        /**< Whether or not a source is actually present.*/
    double *Az;         /**< Array containing azimuth angles of source on-sky, relative to source center.*/
    int nAz;            /**< Number of elements in Az.*/
    double *El;         /**< Array containing elevation angles on-sky, relative to source center.*/
    int nEl;            /**< Number of elements in El.*/
    double *I_nu;       /**< Flat array of specific intensities, indexed as [i * ni + k * ni * nj + j].
                          Here, i is axis 0 (Az), j axis 1 (El) and k axis 2 (frequency). Note that size of I_nu = nAz * nEl * nfreqs.*/
    double *freqs_src;  /**< Frequencies over which the source is defined, in Hertz.*/
    int nfreqs_src;     /**< Number of frequencies in freq_src.*/
};

struct SimParams {
    double t_obs;       /**< Observation time in seconds.*/
    int nTimes;         /**< Total number of time calculations.*/
    int nThreads;       /**< Number of threads to use for computation.*/
    double t0;          /**< Starting time of simulation. If not given will default to 0.*/
    int OFF_empty;   /**< Wether to use no source or interpolate source during OFF chopping.*/
    int use_noise;      /**< Whether to add photon noise. For real life situations, this should be 1. Only set to 0 for debug purposes.*/
};

struct Output {
    double *signal;     /**< Timestream of output signal, time = slow axis, frequency = fast axis.*/
    double *Az;         /**< Timestream of Azimuth angle, in degrees.*/
    double *El;         /**< Timestream of Elevation angle, in degrees.*/
    int *flag;          /**< Timestream of flags specifying chop/nod position. 0 for ON, 1 for OFF-RIGHT, 2 for OFF-LEFT.*/     
    double t_thread;   /**< Time spent calculating in thread.*/
};

// Local structs - for use internally
struct Effs {
    double eta_tot_chain; /**< Total constant efficiency after atmosphere.*/
    double eta_tot_gnd;   /**< Total constant efficiency after groundi, including ground emissivity.*/
    double eta_tot_mir;   /**< Total constant efficiency after mirrors, including mirror emissivity.*/
};

#endif
