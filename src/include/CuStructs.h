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
struct CuOutput;
struct CuArrSpec;

struct CuArrSpec {
    float start;        /**< Starting value of array.*/
    float step;         /**< Stepsize of array.*/
    int num;            /**< Number of array elements.*/
};

// Ctypes structs - for communicating with python
struct CuInstrument {
    int nf_ch;          /**< Number of elements in freqs.*/
    
    struct CuArrSpec f_spec;

    float eta_inst;    /**< Instrument efficiency.*/
    float eta_misc;    /**< Miscellaneous constant efficiencies.*/
    float f_sample; /**< Readout frequency of instrument in Hertz.*/
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

    struct CuArrSpec x_spec;
    struct CuArrSpec y_spec;
    struct CuArrSpec f_spec;
    struct CuArrSpec PWV_spec;
    
    float *PWV;        /**< Flat array containing smoothed PWV values at x and y.*/
    float *eta_atm;    /**< Flat array containing all atmospheric transmissions curves over frequency and PWV. Size is nfreqs_atm * nPWV_atm.*/
};

struct CuSource {
    struct CuArrSpec Az_spec;
    struct CuArrSpec El_spec;
    
    float *I_nu;       /**< Flat array of specific intensities.*/
    int nI_nu;         /**< Number of source intensities.*/
};

struct CuOutput {
    float *signal;     /**< Timestream of output signal, time = slow axis, frequency = fast axis.*/
    float *Az;         /**< Timestream of Azimuth angle, in degrees.*/
    float *El;         /**< Timestream of Elevation angle, in degrees.*/
    int *flag;          /**< Timestream of flags specifying chop/nod position. 0 for ON, 1 for OFF-RIGHT, 2 for OFF-LEFT.*/           
    float *t_diag;      /**< Calculation/allocation times, for diagnostics.*/
};


#endif
