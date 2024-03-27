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
    double eta_misc;    /**< Miscellaneous constant efficiencies.*/
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
    double *eta_ap_ON;  /**< Array of aperture efficiencies in ON position, as function of frequency (set by Instrument). Size is nfreqs of instrument.*/
    double *eta_ap_OFF; /**< Array of aperture efficiencies in OFF position, as function of frequency (set by Instrument). Size is nfreqs of instrument.*/
    double eta_mir;     /**< Mirror reflection efficiency.*/
    double eta_fwd;     /**< Telescope forward efficiency.*/
    int scantype;
    double Ax;
    double Axmin;
    double Ay;
    double Aymin;
    double wx;
    double wxmin;
    double wy;
    double wymin;
    double phix;
    double phiy;
};

struct Atmosphere {
    double Tatm;        /**< Temperature of atmosphere in Kelvin.*/
    double v_wind;      /**< Max windspeed in meters per second.*/
    double h_column;    /**< Reference column height of atmosphere, in meters.*/
    
    double x0;          /**< Starting x-coordinate of atmospheric screen.*/
    double dx;          /**< Stepsize of x-coordinate in atmospheric screen.*/
    int nx;            /**< Number of x-coordinates in atmospheric screen.*/
    
    double y0;          /**< Starting y-coordinate of atmospheric screen.*/
    double dy;          /**< Stepsize of y-coordinate in atmospheric screen.*/
    int ny;            /**< Number of y-coordinates in atmospheric screen.*/
    
    double f0;          /**< Starting frequency of ATM-model.*/
    double df;          /**< Stepsize of frequencies of ATM-model.*/
    int nf;            /**< Number of frequencies in ATM-model.*/
    
    double PWV0;        /**< Starting PWV value of ATM-model.*/
    double dPWV;        /**< Stepsize of PWV values in ATM-model.*/
    int nPWV;          /**< Number of PWV values in ATM-model.*/
    
    double *PWV;        /**< Flat array containing smoothed PWV values at x and y.*/
    double *eta_atm;    /**< Flat array containing all atmospheric transmissions curves over frequency and PWV. Size is nfreqs_atm * nPWV_atm.*/
};

struct Source {
    double Az0;         /**< Starting value of source azimuth range.*/
    double dAz;         /**< Stepsize of values in source azimuth range.*/
    int nAz;           /**< Number of values in source azimuth range.*/
    
    double El0;         /**< Starting value of source elevation range.*/
    double dEl;         /**< Stepsize of values in source elevation range.*/
    int nEl;           /**< Number of values in source elevation range.*/
    
    double f0;         /**< Starting value of source frequency range.*/
    double df;         /**< Stepsize of values in source frequency range.*/
    int nf;           /**< Number of values in source frequency range.*/
    
    double *I_nu;       /**< Flat array of specific intensities.*/
    int nI_nu;         /**< Number of source intensities.*/
};

struct SimParams {
    double t_obs;       /**< Observation time in seconds.*/
    int nTimes;         /**< Total number of time calculations.*/
    int nThreads;       /**< Number of threads to use for computation.*/
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
