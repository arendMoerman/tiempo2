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
struct Output;
struct CalOutput;
struct ArrSpec;

struct ArrSpec {
    double start;
    double step;
    int num;
};

// Ctypes structs - for communicating with python
struct Instrument {
    int nf_ch;          /**< Number of elements in freqs.*/
    
    struct ArrSpec f_spec;
    
    double eta_inst;    /**< Instrument efficiency.*/
    double eta_misc;    /**< Miscellaneous constant efficiencies.*/
    double f_sample; /**< Readout frequency of instrument in Hertz.*/
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
    
    struct ArrSpec x_spec;
    struct ArrSpec y_spec;
    struct ArrSpec f_spec;
    struct ArrSpec PWV_spec;
    
    double *PWV;        /**< Flat array containing smoothed PWV values at x and y.*/
    double *eta_atm;    /**< Flat array containing all atmospheric transmissions curves over frequency and PWV. Size is nfreqs_atm * nPWV_atm.*/
};

struct Source {
    struct ArrSpec Az_spec;
    struct ArrSpec El_spec;
    
    double *I_nu;       /**< Flat array of specific intensities.*/
    int nI_nu;         /**< Number of source intensities.*/
};

struct Output {
    double *signal;     /**< Timestream of output signal, time = slow axis, frequency = fast axis.*/
    double *Az;         /**< Timestream of Azimuth angle, in degrees.*/
    double *El;         /**< Timestream of Elevation angle, in degrees.*/
    int *flag;          /**< Timestream of flags specifying chop/nod position. 0 for ON, 1 for OFF-RIGHT, 2 for OFF-LEFT.*/     
    double t_thread;   /**< Time spent calculating in thread.*/
};

struct CalOutput {
    double *power;       /**< Power in Watt as function of filter index (axis 0) and PWV (axis 1).*/
    double *temperature; /**< LOS brightness temperature in Kelvin as function of filter index (axis 0) and PWV (axis 1).*/
};

// Local structs - for use internally
struct Effs {
    double eta_tot_chain; /**< Total constant efficiency after atmosphere.*/
    double eta_tot_gnd;   /**< Total constant efficiency after groundi, including ground emissivity.*/
    double eta_tot_mir;   /**< Total constant efficiency after mirrors, including mirror emissivity.*/
};

#endif
