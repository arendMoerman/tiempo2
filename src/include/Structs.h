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

struct Instrument {
    double *freqs;      /**< Array with frequencies in Hertz.*/
    int nfreqs;         /**< Number of elements in freqs.*/
    int R;              /**< Resolving power of instrument: R = f / df.*/
    double eta_inst;    /**< Instrument efficiency. Defined to contain all efficiencies of the instrument and coupling to telescope.*/
    double freq_sample; /**< Readout frequency of instrument in Hertz.*/
};

struct Telescope {
    double Ttel;        /**< Telescope temperature in Kelvin.*/
    double Tgnd;        /**< Ground temperature in Kelvin.*/
    double Dtel;       /**< Primary aperture diameter in meters.*/
    double dAz_chop;    /**< Azimuthal separation between chopping paths.*/
    double freq_chop;   /**< Chopping frequency in Hertz. If < 0, no chopping.*/
    double *eta_ap;     /**< Array of aperture efficiencies, as function of frequency (set by Instrument). Size is nfreqs of instrument.*/
    double eta_mir;     /**< Mirror reflection efficiency.*/
    double eta_fwd;     /**< Telescope forward efficiency.*/
};

struct Atmosphere {
    double Tatm;        /**< Temperature of atmosphere in Kelvin.*/
    double vel_w;       /**< Max windspeed in meters per second.*/
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
    double *Az;         /**< Array containing azimuth angles of source on-sky, relative to source center.*/
    int nAz;            /**< Number of elements in Az.*/
    double *El;         /**< Array containing elevation angles on-sky, relative to source center.*/
    int nEl;            /**< Number of elements in El.*/
    double *I_nu;       /**< Flat array of specific intensities, indexed as [i * ni + k * ni * nj + j].
                          Here, i is axis 0 (Az), j axis 1 (El) and k axis 2 (frequency). Note that size of I_nu = nAz * nEl * nfreqs.*/
};

struct SimParams {
    double t_obs;       /**< Observation time in seconds.*/
    int nThreads;       /**< Number of threads to use for computation.*/
};

#endif
