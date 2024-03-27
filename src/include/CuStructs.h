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
    float eta_misc;    /**< Miscellaneous constant efficiencies.*/
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
    float *eta_ap_ON;  /**< Array of aperture efficiencies in ON position, as function of frequency (set by Instrument). Size is nfreqs of instrument.*/
    float *eta_ap_OFF; /**< Array of aperture efficiencies in OFF position, as function of frequency (set by Instrument). Size is nfreqs of instrument.*/
    float eta_mir;     /**< Mirror reflection efficiency.*/
    float eta_fwd;     /**< Telescope forward efficiency.*/
    int scantype;
    float Ax;
    float Axmin;
    float Ay;
    float Aymin;
    float wx;
    float wxmin;
    float wy;
    float wymin;
    float phix;
    float phiy;
};

struct CuAtmosphere {
    float Tatm;        /**< Temperature of atmosphere in Kelvin.*/
    float v_wind;      /**< Max windspeed in meters per second.*/
    float h_column;    /**< Reference column height of atmosphere, in meters.*/
    
    float x0;          /**< Starting x-coordinate of atmospheric screen.*/
    float dx;          /**< Stepsize of x-coordinate in atmospheric screen.*/
    int nx;            /**< Number of x-coordinates in atmospheric screen.*/
    
    float y0;          /**< Starting y-coordinate of atmospheric screen.*/
    float dy;          /**< Stepsize of y-coordinate in atmospheric screen.*/
    int ny;            /**< Number of y-coordinates in atmospheric screen.*/
    
    float f0;          /**< Starting frequency of ATM-model.*/
    float df;          /**< Stepsize of frequencies of ATM-model.*/
    int nf;            /**< Number of frequencies in ATM-model.*/
    
    float PWV0;        /**< Starting PWV value of ATM-model.*/
    float dPWV;        /**< Stepsize of PWV values in ATM-model.*/
    int nPWV;          /**< Number of PWV values in ATM-model.*/
    
    float *PWV;        /**< Flat array containing smoothed PWV values at x and y.*/
    float *eta_atm;    /**< Flat array containing all atmospheric transmissions curves over frequency and PWV. Size is nfreqs_atm * nPWV_atm.*/
};

struct CuSource {
    float Az0;         /**< Starting value of source azimuth range.*/
    float dAz;         /**< Stepsize of values in source azimuth range.*/
    int nAz;           /**< Number of values in source azimuth range.*/
    
    float El0;         /**< Starting value of source elevation range.*/
    float dEl;         /**< Stepsize of values in source elevation range.*/
    int nEl;           /**< Number of values in source elevation range.*/
    
    float f0;         /**< Starting value of source frequency range.*/
    float df;         /**< Stepsize of values in source frequency range.*/
    int nf;           /**< Number of values in source frequency range.*/
    
    float *I_nu;       /**< Flat array of specific intensities.*/
    int nI_nu;         /**< Number of source intensities.*/
};

struct CuSimParams {
    float t_obs;       /**< Observation time in seconds.*/
    int nTimes;         /**< Total number of time calculations.*/
    int nThreads;       /**< Number of threads to use for computation.*/
};

struct CuOutput {
    float *signal;     /**< Timestream of output signal, time = slow axis, frequency = fast axis.*/
    float *Az;         /**< Timestream of Azimuth angle, in degrees.*/
    float *El;         /**< Timestream of Elevation angle, in degrees.*/
    int *flag;          /**< Timestream of flags specifying chop/nod position. 0 for ON, 1 for OFF-RIGHT, 2 for OFF-LEFT.*/           
    float *t_diag;      /**< Calculation/allocation times, for diagnostics.*/
};

#endif
