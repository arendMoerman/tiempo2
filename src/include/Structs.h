/*!
 * \file
 * \brief Data structures for receiving data from Python interface.
 **/

#ifndef __Structs_h
#define __Structs_h

template<typename T>
struct Instrument;

template<typename T>
struct Telescope;

template<typename T>
struct Atmosphere;

template<typename T>
struct Source;

template<typename T>
struct Output;

template<typename T>
struct CalOutput;

template<typename T>
struct ArrSpec;

template<typename T>
struct ArrSpec {
    T start;
    T step;
    int num;
};

// Ctypes structs - for communicating with python
template<typename T>
struct Instrument {
    int nf_ch;          /**< Number of elements in freqs.*/
    
    struct ArrSpec<T> f_spec;

    //T R;            /**< Resolving power.*/
    //T *f_ch;        /**< Array with channel frequencies.*/
    //int order;      /**< Order of Lorentzian filters.*/
    
    //T *eta_filt;     /**<Peak height of filter, for each channel.*/
    T eta_inst;    /**< Instrument efficiency.*/
    T eta_misc;    /**< Miscellaneous constant efficiencies.*/
    T f_sample; /**< Readout frequency of instrument in Hertz.*/
    T *filterbank; /**< Array with filterbank matrix, flattened.*/
    T delta;       /**< Superconducting bandgap energy in Joules.*/
    T eta_pb;      /**< Pair breaking efficiency of superconductor.*/
};

template<typename T>
struct Telescope {
    T Ttel;        /**< Telescope temperature in Kelvin.*/
    T Tgnd;        /**< Ground temperature in Kelvin.*/
    T Dtel;        /**< Primary aperture diameter in meters.*/
    int chop_mode;      /**< Chopping mode. 0 is 'none', 1 is 'direct', 2 is 'abba'.*/
    T dAz_chop;    /**< Azimuthal separation between chopping paths.*/
    T freq_chop;   /**< Chopping frequency in Hertz. If < 0, no chopping.*/
    T freq_nod;    /**< Nodding frequency in Hertz.*/
    T *eta_ap_ON;  /**< Array of aperture efficiencies in ON position, as function of frequency (set by Instrument). Size is nfreqs of instrument.*/
    T *eta_ap_OFF; /**< Array of aperture efficiencies in OFF position, as function of frequency (set by Instrument). Size is nfreqs of instrument.*/
    T eta_mir;     /**< Mirror reflection efficiency.*/
    T eta_fwd;     /**< Telescope forward efficiency.*/
    int scantype;
    T El0;
    T Ax;
    T Axmin;
    T Ay;
    T Aymin;
    T wx;
    T wxmin;
    T wy;
    T wymin;
    T phix;
    T phiy;
};

template<typename T>
struct Atmosphere {
    T Tatm;        /**< Temperature of atmosphere in Kelvin.*/
    T v_wind;      /**< Max windspeed in meters per second.*/
    T h_column;    /**< Reference column height of atmosphere, in meters.*/
    T dx;          /**< Gridsize along x axis in meters.*/
    T dy;          /**< Gridsize along y axis in meters.*/
    char* path;    /**< Path to prepd folder.*/
};

template<typename T>
struct Source {
    struct ArrSpec<T> Az_spec;
    struct ArrSpec<T> El_spec;
    
    T *I_nu;       /**< Flat array of specific intensities.*/
    int nI_nu;         /**< Number of source intensities.*/
};

template<typename T>
struct Output {
    T *signal;     /**< Timestream of output signal, time = slow axis, frequency = fast axis.*/
    T *Az;         /**< Timestream of Azimuth angle, in degrees.*/
    T *El;         /**< Timestream of Elevation angle, in degrees.*/
    int *flag;          /**< Timestream of flags specifying chop/nod position. 0 for ON, 1 for OFF-RIGHT, 2 for OFF-LEFT.*/     
    T *t_diag;      /**< Calculation/allocation times, for diagnostics.*/
    T t_thread;   /**< Time spent calculating in thread.*/
};

template<typename T>
struct CalOutput {
    T *power;       /**< Power in Watt as function of filter index (axis 0) and PWV (axis 1).*/
    T *temperature; /**< LOS brightness temperature in Kelvin as function of filter index (axis 0) and PWV (axis 1).*/
};

// Local structs - for use internally
template<typename T>
struct Effs {
    T eta_tot_chain; /**< Total constant efficiency after atmosphere.*/
    T eta_tot_gnd;   /**< Total constant efficiency after groundi, including ground emissivity.*/
    T eta_tot_mir;   /**< Total constant efficiency after mirrors, including mirror emissivity.*/
};

#endif
