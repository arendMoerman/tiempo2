#include "InterfaceCUDA.h"

/*! \file Kernels.cu
    \brief Definitions of CUDA kernels for TiEMPO2.

    author: Arend Moerman
*/

// PHYSICAL CONSTANTS
__constant__ float cPI;                     // Pi
__constant__ float cCL;                     // speed of light
__constant__ float cHP;                     // Planck constant
__constant__ float cKB;                     // Boltzmann constant

// OBSERVATION-INSTRUMENT PARAMETERS
__constant__ float const_effs[CEFFSSIZE];   // Contains constant efficiencies:chain, gnd, mir, pb 
__constant__ float cdt;                     // Timestep
__constant__ float cfreq_chop;              // Chopping frequency
__constant__ float cfreq_nod;               // Nodding frequency
__constant__ float cfreq_sample;            // Sampling frequency of readout
__constant__ float cdAz_chop;               // Chopping throw
__constant__ float cdelta;                  // Bandgap energy of MKID
__constant__ int cnt;                       // Number of time evals
__constant__ int cnf_filt;                  // Number of filter freqs
__constant__ int cchop_mode;                // What chopping scheme to use

// ATMOSPHERE PARAMETERS
__constant__ float ch_column;               // Column height
__constant__ float cv_wind;                 // Windspeed

__constant__ float cx0_atm;                 /**< Start x-value of atmosphere screen.*/
__constant__ float cdx_atm;                 /**< Stepsize of x-array in atmosphere screen.*/
__constant__ int cnx_atm;                   /**< Number of x points in screen.*/

__constant__ float cy0_atm;                 /**< Start y-value of atmosphere screen.*/
__constant__ float cdy_atm;                 /**< Stepsize of y-array in atmosphere screen.*/
__constant__ int cny_atm;                   /**< Number of y points in screen.*/

__constant__ float cf0_atm;                 /**< Start frequency for ATM-model.*/
__constant__ float cdf_atm;                 /**< Stepsize of frequencies in ATM-model.*/
__constant__ int cnf_atm;                   /**< Number of frequencies in ATM-model.*/

__constant__ float cPWV0_atm;               /**< Start PWV value for ATM-model.*/
__constant__ float cdPWV_atm;               /**< Stepsize for PWV in ATM-model.*/
__constant__ int cnPWV_atm;                 /**< Number of PWV values in ATM-model.*/

// SOURCE PARAMETERS
__constant__ float cAz0_src;                 /**< Start azimuth of source.*/
__constant__ float cdAz_src;                 /**< Stepsize of azimuth.*/
__constant__ int cnAz_src;                   /**< Number of elements in source azimuth.*/

__constant__ float cEl0_src;                 /**< Start elevation of source.*/
__constant__ float cdEl_src;                 /**< Stepsize of elevation.*/
__constant__ int cnEl_src;                   /**< Number of elements in source elevation.*/

__constant__ float cf0_src;                  /**< Start frequency of source.*/
__constant__ float cdf_src;                 /**< Stepsize of array.*/
__constant__ int cnf_src;                   /**< Number of elements in source frequencies.*/

// SCAN PARAMETERS
__constant__ int cscantype;
__constant__ float cAx;
__constant__ float cAxmin;
__constant__ float cAy;
__constant__ float cAymin;
__constant__ float cwx;
__constant__ float cwxmin;
__constant__ float cwy;
__constant__ float cwymin;
__constant__ float cphix;
__constant__ float cphiy;

// TEXTURE MEMORY
texture<float, cudaTextureType1D, cudaReadModeElementType> tex_filterbank;
texture<float, cudaTextureType1D, cudaReadModeElementType> tex_eta_ap_ON;
texture<float, cudaTextureType1D, cudaReadModeElementType> tex_eta_ap_OFF;

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }

/**
  Check CUDA API error status of call.
 
  Wrapper for finding errors in CUDA API calls.
 
  @param code The errorcode returned from failed API call.
  @param file The file in which failure occured.
  @param line The line in file in which error occured.
  @param abort Exit code upon error.
 */
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort = true)
{
   if (code != cudaSuccess)
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

/**
  Struct for passing an array to kernel.

  Eliminates the need to pass long arrays that can be described by start, step and number.
 */
typedef struct
{
    float start, step;      /**< Start of array and stepsize.*/
    int number;             /**< Number of elements in array.*/
} arr_data;

/**
  Write a CUDA array to a file, for debugging.

  @param array Pointer to array of type T.
  @param s_array Size of array.
  @param name_txt Name of file to write array to. Name is appended with '.txt' by the function itself.
 */
template <typename T>
__host__ void writeArray(T *array, int s_array, std::string name_txt) {
    
    T *h_array = new T[s_array];
    gpuErrchk( cudaMemcpy(h_array, array, s_array * sizeof(T), cudaMemcpyDeviceToHost) );
    
    std::ofstream myfile (name_txt + ".txt");
    if (myfile.is_open())
    {
        for(int count = 0; count < s_array; count ++){
            myfile << h_array[count] << "\n" ;
        }

        myfile.close();
    }
    else std::cout << "Unable to open file";
    delete[] h_array;
}

/**
  Calculate Planckian distribution.

  Used for calculating blackbody intensities of atmosphere, ground and telescope.

  @param T Temperature of blackbody, in Kelvin.
  @param nu Frequency at which to evaluate blackbody, in Hertz.

  @returns Blackbody intensity.
 */
__host__ float getPlanck(float T, float nu)
{
    float CL = 2.9979246e8; // m s^-1
    float HP = 6.62607015e-34;
    float KB = 1.380649e-23;
    
    float prefac = 2 * HP * nu*nu*nu / (CL*CL);
    float dist = 1 / (exp(HP*nu / (KB*T)) - 1);
    
    return prefac * dist;
}

/**
  Calculate sign of number.

  Used in determining chop-nod state.

  @param val Value of number.
 */
__device__ __inline__ void sgn(float val, int &out) {
    out = (float(0) < val) - (val < float(0));
}

/**
  Convert angle in arcseconds to radian.

  @param ang Angle in arcseconds.
 */
__device__ __inline__ float as2rad(float ang) {
    return ang / 3600 / 180 * cPI;
}

/**
  Calculate new Azimuth-Elevation co-ordinate, accoding to chop position.

  This function just "scans" a single point, so seems sort of pointless. 
  Still implemented for completeness.

  @param center Az-El co-ordinate of point to observe, w.r.t. source Az-El.
  @param out Container for storing output Az-El co-ordinate.
  @param chop Whether chopper is in A (false) or B (true).
  @param sep Angular throw between chop A and B, in degrees.
 */
__device__ __inline__ void scanPoint(AzEl* center, AzEl* out, bool chop, float sep = 0.) {
    float offset = 0.;
    
    if (chop) {
        offset = sep;
    }    
    out->Az = center->Az + offset;
    out->El = center->El;
}

__device__ __inline__ void scanDaisy(AzEl* center, AzEl* out, float t, bool chop, float sep = 0.) {
    float offset = 0.;
    
    if (chop) {
        offset = sep;
    }    
    
    out->Az = center->Az + offset + cAx*sinf(cwx*t)*cosf(cwx*t + as2rad(cphix)) + cAxmin*sinf(cwxmin*t)*cosf(cwxmin*t + as2rad(cphix));
    out->El = center->El + cAy*sinf(cwy*t)*sinf(cwy*t + as2rad(cphiy)) + cAymin*sinf(cwymin*t)*sinf(cwymin*t + as2rad(cphiy)) - cAy;
}

/**
  Convert an Az-El co-ordinate to a projected x-y co-ordinate on the atmosphere.

  @param angles Az-El co-ordinate to convert.
  @param out Container for storing the calculated x-y point.
 */
__device__ __inline__ void convertAnglesToSpatialAtm(AzEl* angles, xy_atm* out) {
    
    float coord = tanf(cPI * angles->Az / 180.) * ch_column;
    
    out->xAz = coord;
    coord = tanf(cPI * angles->El / 180.) * ch_column;
    out->yEl = coord;
}

__device__ __inline__ void getABBA_posflag(float &t_start, AzEl *center, AzEl *pointing, int &flagout) {
    int n_chop;
    int n_nod;
    int position;

    bool chop_flag;

    float is_in_lower_half;
    int nod_flag;

    n_chop = floorf(t_start * cfreq_chop);
    n_nod = floorf(t_start * cfreq_nod);
    
    chop_flag = (n_chop % 2 != 0); // If even (false), ON. Odd (true), OFF.
    nod_flag = -1 + 2 * (n_nod % 2 != 0); // If even (false), AB. Odd (true), BA.
    
    is_in_lower_half = (t_start - n_nod / cfreq_nod) - (1 / cfreq_nod / 2);
    sgn(is_in_lower_half, position);
    position *= nod_flag;
    
    scanPoint(center, pointing, chop_flag, position * cdAz_chop);
    flagout = chop_flag * position + (1 - chop_flag) * (1 - position);
}

__device__ __inline__ void getONOFF_posflag(float &t_start, AzEl *center, AzEl *pointing, int &flagout) {
    int n_chop;
    bool chop_flag;

    n_chop = floorf(t_start * cfreq_chop);
    
    chop_flag = (n_chop % 2 != 0); // If even (false), ON. Odd (true), OFF.
    if(cscantype == 0) {scanPoint(center, pointing, chop_flag, cdAz_chop);}
    else if(cscantype == 1) {scanDaisy(center, pointing, t_start, chop_flag, cdAz_chop);}
    flagout = chop_flag;
}

__device__ __inline__ void getnochop_posflag(float &t_start, AzEl *center, AzEl *pointing, int &flagout) {
    if(cscantype == 0) {scanPoint(center, pointing, 0, cdAz_chop);}
    else if(cscantype == 1) {scanDaisy(center, pointing, t_start, 0, cdAz_chop);}
    flagout = 0;
}

/**
  Initialize CUDA.
 
  Instantiate program and populate constant memory.
 
  @param instrument CuInstrument object containing instrument to be simulated.
  @param telescope CuTelescope object containing telescope to be simulated.
  @param simparams CuSimParams object containing simulation parameters.
  @param source CuSource object containing source definitions.
  @param atmosphere CuAtmosphere object containing atmosphere parameters.
  @param nThreads Number of CUDA threads per block.
 
  @return BT Array of two dim3 objects, containing number of blocks per grid and number of threads per block.
 */
__host__ void initCUDA(CuInstrument *instrument, CuTelescope *telescope, CuSimParams *simparams, CuSource *source, CuAtmosphere *atmosphere) {
    //int nBlocks = ceilf((float)simparams->nTimes / nThreads);

    // Calculate nr of blocks per grid and nr of threads per block
    //dim3 nrb(nBlocks); dim3 nrt(nThreads);


    float PI = 3.1415926; /* pi */
    float CL = 2.9979246e8; // m s^-1
    float HP = 6.62607015e-34;
    float KB = 1.380649e-23;

    // Pack constant array
    float _con[CEFFSSIZE] = {instrument->eta_inst * instrument->eta_misc * telescope->eta_fwd * telescope->eta_mir * 0.5,
        instrument->eta_inst * instrument->eta_misc * (1 - telescope->eta_fwd) * telescope->eta_mir * 0.5,
        instrument->eta_inst * instrument->eta_misc * (1 - telescope->eta_mir) * 0.5, 
        instrument->eta_pb};

    float dt = 1. / instrument->freq_sample;
    
    // PHYSICAL CONSTANTS
    gpuErrchk( cudaMemcpyToSymbol(cPI, &PI, sizeof(float)) );
    gpuErrchk( cudaMemcpyToSymbol(cCL, &CL, sizeof(float)) );
    gpuErrchk( cudaMemcpyToSymbol(cHP, &HP, sizeof(float)) );
    gpuErrchk( cudaMemcpyToSymbol(cKB, &KB, sizeof(float)) );
     
    // OBSERVATION-INSTRUMENT PARAMETERS
    gpuErrchk( cudaMemcpyToSymbol(const_effs, &_con, CEFFSSIZE * sizeof(float)) );
    gpuErrchk( cudaMemcpyToSymbol(cdt, &dt, sizeof(float)) );
    gpuErrchk( cudaMemcpyToSymbol(cfreq_chop, &(telescope->freq_chop), sizeof(float)) );
    gpuErrchk( cudaMemcpyToSymbol(cfreq_nod, &(telescope->freq_nod), sizeof(float)) );
    gpuErrchk( cudaMemcpyToSymbol(cfreq_sample, &(instrument->freq_sample), sizeof(float)) );
    gpuErrchk( cudaMemcpyToSymbol(cdAz_chop, &(telescope->dAz_chop), sizeof(float)) );
    gpuErrchk( cudaMemcpyToSymbol(cdelta, &(instrument->delta), sizeof(float)) );
    gpuErrchk( cudaMemcpyToSymbol(cnt, &(simparams->nTimes), sizeof(int)) );
    gpuErrchk( cudaMemcpyToSymbol(cnf_filt, &(instrument->nfreqs_filt), sizeof(int)) );
    gpuErrchk( cudaMemcpyToSymbol(cchop_mode, &(telescope->chop_mode), sizeof(int)) );
    
    // ATMOSPHERE PARAMETERS
    gpuErrchk( cudaMemcpyToSymbol(ch_column, &(atmosphere->h_column), sizeof(float)) );
    gpuErrchk( cudaMemcpyToSymbol(cv_wind, &(atmosphere->v_wind), sizeof(float)) );
    
    gpuErrchk( cudaMemcpyToSymbol(cx0_atm, &(atmosphere->x0), sizeof(float)) );
    gpuErrchk( cudaMemcpyToSymbol(cdx_atm, &(atmosphere->dx), sizeof(float)) );
    gpuErrchk( cudaMemcpyToSymbol(cnx_atm, &(atmosphere->nx), sizeof(int)) );

    gpuErrchk( cudaMemcpyToSymbol(cy0_atm, &(atmosphere->y0), sizeof(float)) );
    gpuErrchk( cudaMemcpyToSymbol(cdy_atm, &(atmosphere->dy), sizeof(float)) );
    gpuErrchk( cudaMemcpyToSymbol(cny_atm, &(atmosphere->ny), sizeof(int)) );
    
    gpuErrchk( cudaMemcpyToSymbol(cf0_atm, &(atmosphere->f0), sizeof(float)) );
    gpuErrchk( cudaMemcpyToSymbol(cdf_atm, &(atmosphere->df), sizeof(float)) );
    gpuErrchk( cudaMemcpyToSymbol(cnf_atm, &(atmosphere->nf), sizeof(int)) );
    
    gpuErrchk( cudaMemcpyToSymbol(cPWV0_atm, &(atmosphere->PWV0), sizeof(float)) );
    gpuErrchk( cudaMemcpyToSymbol(cdPWV_atm, &(atmosphere->dPWV), sizeof(float)) );
    gpuErrchk( cudaMemcpyToSymbol(cnPWV_atm, &(atmosphere->nPWV), sizeof(int)) );

    // SOURCE PARAMETERS
    gpuErrchk( cudaMemcpyToSymbol(cAz0_src, &(source->Az0), sizeof(float)) );
    gpuErrchk( cudaMemcpyToSymbol(cdAz_src, &(source->dAz), sizeof(float)) );
    gpuErrchk( cudaMemcpyToSymbol(cnAz_src, &(source->nAz), sizeof(int)) );
    
    gpuErrchk( cudaMemcpyToSymbol(cEl0_src, &(source->El0), sizeof(float)) );
    gpuErrchk( cudaMemcpyToSymbol(cdEl_src, &(source->dEl), sizeof(float)) );
    gpuErrchk( cudaMemcpyToSymbol(cnEl_src, &(source->nEl), sizeof(int)) );
    
    gpuErrchk( cudaMemcpyToSymbol(cf0_src, &(source->f0), sizeof(float)) );
    gpuErrchk( cudaMemcpyToSymbol(cdf_src, &(source->df), sizeof(float)) );
    gpuErrchk( cudaMemcpyToSymbol(cnf_src, &(source->nf), sizeof(int)) );

    // SCAN PARAMETERS
    gpuErrchk( cudaMemcpyToSymbol(cscantype, &(telescope->scantype), sizeof(int)) );
    gpuErrchk( cudaMemcpyToSymbol(cAx, &(telescope->Ax), sizeof(float)) );
    gpuErrchk( cudaMemcpyToSymbol(cAxmin, &(telescope->Axmin), sizeof(float)) );
    gpuErrchk( cudaMemcpyToSymbol(cAy, &(telescope->Ay), sizeof(float)) );
    gpuErrchk( cudaMemcpyToSymbol(cAymin, &(telescope->Aymin), sizeof(float)) );
    gpuErrchk( cudaMemcpyToSymbol(cwx, &(telescope->wx), sizeof(float)) );
    gpuErrchk( cudaMemcpyToSymbol(cwxmin, &(telescope->wxmin), sizeof(float)) );
    gpuErrchk( cudaMemcpyToSymbol(cwy, &(telescope->wy), sizeof(float)) );
    gpuErrchk( cudaMemcpyToSymbol(cwymin, &(telescope->wymin), sizeof(float)) );
    gpuErrchk( cudaMemcpyToSymbol(cphix, &(telescope->phix), sizeof(float)) );
    gpuErrchk( cudaMemcpyToSymbol(cphiy, &(telescope->phiy), sizeof(float)) );
}

/**
  Obtain timestream for chopping states, PWV from atmosphere, and rng for Gaussian draw.
  Each timestep gets a new seed for the Gaussian, just to randomise it even harder.
  
  @param state Array of curand states. Should be initialised and sized to total number of threads in grid.
  @param seed Integer describing the seed of the generator.
 */
__global__ void get_chop_pwv_rng(float *center, float *PWV_screen, float *PWV_out, 
                                 int *flagout, float *azout, float *elout,
                                 curandState *state, unsigned long long int seed = 0) {
    if (!seed) {
        seed = clock64();
    }

    int idx = blockIdx.x * blockDim.x + threadIdx.x; 

    if (idx < cnt) {
        // FLOATS
        float time;         // Timepoint for thread in simulation.
        float _PWV_out;     // Container for storing interpolated PWV values.
        float t, u;         // Parameters needed for interpolation.

        // INTEGERS
        int flag;           // Flag for storing chop position.
        int ix, iy;         // Indices for interpolation along an x- and y axis.

        // CUSTOM STRUCTS
        AzEl pointing;      // Struct for storing the current pointing (w.r.t. center_p).
        AzEl center_p;      // Struct for storing the central pointing.
        xy_atm point_atm;   // Struct for storing projected pointing x-y coordinates on atmosphere screen.
        
        center_p.Az = center[0];
        center_p.El = center[1];

        time = idx * cdt;
        
        if(cchop_mode == 0) {getnochop_posflag(time, &center_p, &pointing, flag);}
        else if(cchop_mode == 1) {getONOFF_posflag(time, &center_p, &pointing, flag);}
        else if(cchop_mode == 2) {getABBA_posflag(time, &center_p, &pointing, flag);}

        convertAnglesToSpatialAtm(&pointing, &point_atm);

        // Add wind to this - currently only along x-axis and pretty manual
        point_atm.xAz = point_atm.xAz + cv_wind * time;
        
        ix = floorf((point_atm.xAz - cx0_atm) / cdx_atm);
        iy = floorf((point_atm.yEl - cy0_atm) / cdy_atm);
        
        t = (point_atm.xAz - (cx0_atm + cdx_atm * ix)) / cdx_atm;
        u = (point_atm.yEl - (cy0_atm + cdy_atm * iy)) / cdy_atm;
        
        _PWV_out = (1-t)*(1-u)*PWV_screen[ix * cny_atm + iy];
        _PWV_out += t*(1-u)*PWV_screen[(ix + 1) * cny_atm + iy];
        _PWV_out += (1-t)*u*PWV_screen[ix * cny_atm + iy + 1];
        _PWV_out += t*u*PWV_screen[(ix + 1) * cny_atm + iy + 1];
    
        __syncthreads();
        curand_init(seed, idx, 0, &state[idx]);
        azout[idx] = pointing.Az;
        elout[idx] = pointing.El;
        flagout[idx] = flag;
        
        PWV_out[idx] = _PWV_out;
    }
}

/**
  Calculate power and NEP in each channel.
  This kernel is optimised for position switching observations.

  @param I_atm Array containing blackbody intensity of atmosphere, in SI units.
  @param I_gnd Array containing blackbody intensity of ground, in SI units.
  @param I_tel Array containing blackbody intensity of telescope, in SI units.
  @param I_CMB Array containing blackbody intensity of CMB, in SI units.
  @param sigout Array for storing output power, for each channel, for each time, in SI units.
  @param nepout Array for storing output NEP, for each channel, for each time, in SI units.
  @param flagout Array for storing wether beam is in chop A or B, in nod AB or BA.
  @param PWV_trace Array containing PWV value of atmosphere as seen by telescope over observation, in millimeters.
  @param eta_atm Array with transmission parameters as function of PWV and frequency.
  @param source Array containing source intensity at three pointings, as function of frequency, in SI units.
 */
__global__ void calcPowerNEP_PS(float *I_atm, float *I_gnd, float *I_tel, float *I_CMB,
        float *sigout, float *nepout, int *flagout,
        float *PWV_trace, float *eta_atm, float *source) {
    
    int idx = blockIdx.x * blockDim.x + threadIdx.x; 
    int idy = blockIdx.y * blockDim.y + threadIdx.y; 

    extern __shared__ float sh_I_nu[]; // Shared array for storing source intensities.

    if (idx < cnt && idy < cnf_src) {
        // FLOATS
        float eta_atm_interp;   // Interpolated eta_atm, over frequency and PWV
        float freq;             // Bin frequency
        float I_nu;             // Specific intensity of source.
        float PSD_nu;           // Local variable for storing PSD.
        float eta_kj;           // Filter efficiency for bin j, at channel k.
        float PWV_tr;           // Local variable for storing PWV value at time.
        float eta_ap;           // Local variable for storing aperture efficiency
        float sigfactor;        // Factor for calculating power. Perform outside of channel loop for speed.
        float nepfactor1;       // Factor 1 for calculating NEP. Perform outside of channel loop for speed.
        float nepfactor2;       // Factor 2 for calculating NEP. Perform outside of channel loop for speed.

        // INTEGERS
        int flag;               // Flag for storing chop position.
        int idx_point;          // Index for pointing of source (0 is OFF-B, 1 is ON-AB, 2 is OFF-A)

        // Reusable symbols for interpolation stuff - listed separately for readability
        int ix, iy;
        float t, u;
        
        // Reuse three time indices for the three pointings.
        if(threadIdx.x < 3) {
            sh_I_nu[threadIdx.x*blockDim.y + threadIdx.y] = source[threadIdx.x*cnf_src + idy];
        }

        __syncthreads();

        // Transfer global memory into local variables.
        PWV_tr = PWV_trace[idx];
        flag = flagout[idx];

        // Determine idx in source from chopping flag
        if(flag==0 or flag==2) {idx_point = 1;}
        else if(flag==1) {idx_point = 2;}
        else {idx_point = 0;}
        
        freq = cf0_src + cdf_src * idy;
            
        ix = floorf((PWV_tr - cPWV0_atm) / cdPWV_atm);
        iy = floorf((freq - cf0_atm) / cdf_atm);
        
        t = (PWV_tr - (cPWV0_atm + cdPWV_atm * ix)) / cdPWV_atm;
        u = (freq - (cf0_atm + cdf_atm * iy)) / cdf_atm;
        
        eta_atm_interp = (1-t)*(1-u)*eta_atm[ix * cnf_atm + iy];
        eta_atm_interp += t*(1-u)*eta_atm[(ix + 1) * cnf_atm + iy];
        eta_atm_interp += (1-t)*u*eta_atm[ix * cnf_atm + iy + 1];
        eta_atm_interp += t*u*eta_atm[(ix + 1) * cnf_atm + iy + 1];

        I_nu = sh_I_nu[idx_point * blockDim.y + threadIdx.y];

        if(flag == 0 or flag == -1) {
            eta_ap = tex1Dfetch(tex_eta_ap_ON, idy); 
        }

        else {
            eta_ap = tex1Dfetch(tex_eta_ap_OFF, idy);
        }

        PSD_nu = eta_ap * eta_atm_interp * const_effs[0] * I_nu
            + ( const_effs[0] * (1 - eta_atm_interp) * I_atm[idy] 
            + const_effs[1] * I_gnd[idy]
            + const_effs[2] * I_tel[idy]) 
            * cCL*cCL / (freq*freq);

        sigfactor = PSD_nu * cdf_src;
        nepfactor1 = PSD_nu * (cHP * freq + 2 * cdelta / const_effs[3]) * cdf_src;
        nepfactor2 = PSD_nu * PSD_nu * cdf_src;

        #pragma unroll 
        for(int k=0; k<cnf_filt; k++) {
            eta_kj = tex1Dfetch( tex_filterbank, k*cnf_src + idy);
            atomicAdd(&sigout[k*cnt + idx], eta_kj * sigfactor); 
            atomicAdd(&nepout[k*cnt + idx], eta_kj * (nepfactor1 + nepfactor2 * eta_kj)); 
        }
    }
}

/**
  Main simulation kernel. This is where the magic happens.

  @param I_atm Array containing blackbody intensity of atmosphere, in SI units.
  @param I_gnd Array containing blackbody intensity of ground, in SI units.
  @param I_tel Array containing blackbody intensity of telescope, in SI units.
  @param I_CMB Array containing blackbody intensity of CMB, in SI units.
  @param sigout Array for storing output power, for each channel, for each time, in SI units.
  @param flagout Array for storing wether beam is in chop A or B, in nod AB or BA.
  @param eta_ap Array containing aperture efficiencies, for each bin frequency.
  @param PWV_screen Array containing PWV value of atmosphere, over the range described by x_atm and y_atm, in millimeters.
  @param eta_atm Array with transmission parameters as fuiunction of freqs_atm and PWV_atm.
  @param filterbank Array containing filterbank of instrument.
  @param source Array containing source intensity, as function of azsrc, elsrc and freqs_src, in SI units.
  @param state Array with states for drawing random Gaussian values for noise calculations.
 */
__global__ void runSimulation(float *I_atm, float *I_gnd, float *I_tel, float *I_CMB,
        float *sigout, float *nepout, float *azout, float *elout, int *flagout,
        float *PWV_trace, float *eta_atm, float *source) {
    
    int idx = blockIdx.x * blockDim.x + threadIdx.x; 
    int idy = blockIdx.y * blockDim.y + threadIdx.y; 

    if (idx < cnt && idy < cnf_src) {
        float eta_atm_interp;   // Interpolated eta_atm, over frequency and PWV
        float freq;             // Bin frequency
        float I_nu;             // Specific intensity of source.
        float PSD_nu;           // Local variable for storing PSD.
        float eta_kj;           // Filter efficiency for bin j, at channel k.
        float PWV_tr;           // Local variable for storing PWV value at time.
        float eta_ap;           // Local variable for storing aperture efficiency
        float sigfactor;        // Factor for calculating power. Perform outside of channel loop for speed.
        float nepfactor1;       // Factor 1 for calculating NEP. Perform outside of channel loop for speed.
        float nepfactor2;       // Factor 2 for calculating NEP. Perform outside of channel loop for speed.

        // Reusable symbols for interpolation stuff
        int x0y0, x1y0, x0y1, x1y1, ix, iy;
        float t, u;
            
        AzEl pointing;

        pointing.Az = azout[idx];
        pointing.El = elout[idx];

        int iAz = floorf((pointing.Az - cAz0_src) / cdAz_src);
        int iEl = floorf((pointing.El - cEl0_src) / cdEl_src);

        PWV_tr = PWV_trace[idx];
        int flag = flagout[idx];

        float Az_src_max = cAz0_src + cdAz_src * (cnAz_src - 1);
        float El_src_max = cEl0_src + cdEl_src * (cnEl_src - 1);

        bool offsource = ((pointing.Az < cAz0_src) or (pointing.Az > Az_src_max)) or 
                         ((pointing.El < cEl0_src) or (pointing.El > El_src_max));

        freq = cf0_src + cdf_src * idy;
            
        ix = floorf((PWV_tr - cPWV0_atm) / cdPWV_atm);
        iy = floorf((freq - cf0_atm) / cdf_atm);
        
        t = (PWV_tr - (cPWV0_atm + cdPWV_atm * ix)) / cdPWV_atm;
        u = (freq - (cf0_atm + cdf_atm * iy)) / cdf_atm;
        
        eta_atm_interp = (1-t)*(1-u)*eta_atm[ix * cnf_atm + iy];
        eta_atm_interp += t*(1-u)*eta_atm[(ix + 1) * cnf_atm + iy];
        eta_atm_interp += (1-t)*u*eta_atm[ix * cnf_atm + iy + 1];
        eta_atm_interp += t*u*eta_atm[(ix + 1) * cnf_atm + iy + 1];

        if(offsource) {I_nu = I_CMB[idy];}
        
        else {
            x0y0 = iAz * cnf_src + iEl * cnf_src * cnAz_src;
            x1y0 = (iAz + 1) * cnf_src + iEl * cnf_src * cnAz_src;
            x0y1 = iAz * cnf_src + (iEl+1) * cnf_src * cnAz_src;
            x1y1 = (iAz+1) * cnf_src + (iEl+1) * cnf_src * cnAz_src;
            
            t = (pointing.Az - (cAz0_src + cdAz_src*iAz)) / cdAz_src;
            u = (pointing.El - (cEl0_src + cdEl_src*iEl)) / cdEl_src;
            
            I_nu = (1-t)*(1-u) * source[x0y0 + idy];
            I_nu += t*(1-u) * source[x1y0 + idy];
            I_nu += (1-t)*u * source[x0y1 + idy];
            I_nu += t*u * source[x1y1 + idy];
        }

        if(flag == 0 or flag == -1) {
            eta_ap = tex1Dfetch(tex_eta_ap_ON, idy); 
        }

        else {
            eta_ap = tex1Dfetch(tex_eta_ap_OFF, idy);
        }

        PSD_nu = eta_ap * eta_atm_interp * const_effs[0] * I_nu
            + ( const_effs[0] * (1 - eta_atm_interp) * I_atm[idy] 
            + const_effs[1] * I_gnd[idy]
            + const_effs[2] * I_tel[idy]) 
            * cCL*cCL / (freq*freq);
        //printf("%d %d %d %d %d %d\n", start_x0y0, start_x1y0, start_x0y1, start_x1y1, _i_src, _j_src);

        sigfactor = PSD_nu * cdf_src;
        nepfactor1 = PSD_nu * (cHP * freq + 2 * cdelta / const_effs[3]) * cdf_src;
        nepfactor2 = PSD_nu * PSD_nu * cdf_src;

        #pragma unroll 
        for(int k=0; k<cnf_filt; k++) {
            eta_kj = tex1Dfetch( tex_filterbank, k*cnf_src + idy);
            atomicAdd(&sigout[k*cnt + idx], eta_kj * sigfactor); 
            atomicAdd(&nepout[k*cnt + idx], eta_kj * (nepfactor1 + nepfactor2 * eta_kj)); 
        }
    }
}

__global__ void getSignal(float *sigout, float *nepout, curandState *state) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x; 
    if (idx < cnt) {
        curandState localState = state[idx];
        //state[idx] = localState;
        float sqrt_samp = sqrtf(0.5 * cfreq_sample); // Constant term needed for noise calculation
        float sigma_k, P_k;


        for(int k=0; k<cnf_filt; k++) {
            sigma_k = sqrtf(2 * nepout[k*cnt + idx]) * sqrt_samp;
            P_k = sigma_k * curand_normal(&localState);

            state[idx] = localState;

            sigout[k*cnt + idx] += P_k;
        }
    }
}

/**
  Run a TiEMPO2 simulation using CUDA.
 
  This function is exposed to the ctypes interface and can be called from Python..
 
  @param instrument CuInstrument object containing instrument to be simulated.
  @param telescope CuTelescope object containing telescope to be simulated.
  @param atmosphere CuAtmosphere object containing atmosphere parameters.
  @param source CuSource object containing source definitions.
  @param simparams CuSimParams object containing simulation parameters.
  @param output CuOutput object for storing simulation output.
 */
void runTiEMPO2_CUDA(CuInstrument *instrument, CuTelescope *telescope, CuAtmosphere *atmosphere, CuSource *source, 
        CuSimParams *simparams, CuOutput *output) {
    // FLOATS
    float *d_sigout;        // Device pointer for output power array
    float *d_nepout;        // Device pointer for output NEP array
    float *d_azout;         // Device pointer for output Azimuth array 
    float *d_elout;         // Device pointer for output Elevation array
    float *d_I_nu;          // Device pointer for source intensities
    
    // INTEGERS
    int *d_flagout;         // Device pointer for output chopping flags
    int nffnt;              // Number of filter frequencies times number of time evaluations
    int numSMs;             // Number of streaming multiprocessors on GPU
    int nThreads1D;         // Number of threads in 1D block
    int nThreads2D;         // Number of threads per dimension of 2D block
    int nBlocks1D;          // Number of 1D blocks, in terms of number of SMs
    int nBlocks2Dx;         // Number of 2D blocks, in terms of number of SMs, along the x-dimension
    int nBlocks2Dy;         // Number of 2D blocks, in terms of number of SMs, along the y-dimension
    int nShared;            // Number of elements to allocate in shared memory

    // OTHER DECLARATIONS
    dim3 blockSize1D;       // Size of 1D block (same as nThreads1D, but dim3 type)
    dim3 gridSize1D;        // Number of 1D blocks per grid
    dim3 blockSize2D;       // Size of 2D block, along x and y dimension
    dim3 gridSize2D;        // Number of 2D blocks per grid
    Timer timer;            // Timer class for timing kernel invocations

    nffnt = instrument->nfreqs_filt * simparams->nTimes;
    
    // Setting preferred memory priorities and getting amount of SMs
    cudaDeviceGetAttribute(&numSMs, cudaDevAttrMultiProcessorCount, 0);
    //gpuErrchk( cudaDeviceSetLimit(cudaLimitMallocHeapSize, 128*1024*1024) );

    gpuErrchk( cudaFuncSetCacheConfig(get_chop_pwv_rng, cudaFuncCachePreferL1) );
    gpuErrchk( cudaFuncSetCacheConfig(runSimulation, cudaFuncCachePreferL1) );
    gpuErrchk( cudaFuncSetCacheConfig(calcPowerNEP_PS, cudaFuncCachePreferL1) );
    gpuErrchk( cudaFuncSetCacheConfig(getSignal, cudaFuncCachePreferL1) );
    
    nThreads1D = 256;
    nBlocks1D = ceilf((float)simparams->nTimes / nThreads1D / numSMs);

    blockSize1D = nThreads1D;
    gridSize1D = nBlocks1D*numSMs;
    
    nThreads2D = 16;
    nBlocks2Dx = ceilf((float)simparams->nTimes / nThreads2D / numSMs);
    nBlocks2Dy = ceilf((float)source->nf / nThreads2D / numSMs);

    blockSize2D = dim3(nThreads2D, nThreads2D);
    gridSize2D = dim3(nBlocks2Dx * numSMs, nBlocks2Dy * numSMs);

    timer.start();
    
    initCUDA(instrument, telescope, simparams, source, atmosphere);
    
    // Allocate output arrays
    gpuErrchk( cudaMalloc((void**)&d_sigout, nffnt * sizeof(float)) );
    gpuErrchk( cudaMalloc((void**)&d_nepout, nffnt * sizeof(float)) );
    gpuErrchk( cudaMalloc((void**)&d_azout, simparams->nTimes * sizeof(float)) );
    gpuErrchk( cudaMalloc((void**)&d_elout, simparams->nTimes * sizeof(float)) );
    gpuErrchk( cudaMalloc((void**)&d_flagout, simparams->nTimes * sizeof(int)) );

    // Allocate PWV screen now, delete CUDA allocation after first kernel call
    float *dPWV_screen;
    
    int nPWV_screen = atmosphere->nx * atmosphere->ny;
    
    gpuErrchk( cudaMalloc((void**)&dPWV_screen, nPWV_screen * sizeof(float)) );
    gpuErrchk( cudaMemcpy(dPWV_screen, atmosphere->PWV, nPWV_screen * sizeof(float), cudaMemcpyHostToDevice) );
    
    float *PWV_out;
    gpuErrchk( cudaMalloc((void**)&PWV_out, simparams->nTimes * sizeof(float)) );
    
    float pointing_center[2] = {0., 0.};
    float *dpointing_center;
    gpuErrchk( cudaMalloc((void**)&dpointing_center, 2 * sizeof(float)) );
    gpuErrchk( cudaMemcpy(dpointing_center, pointing_center, 2 * sizeof(float), cudaMemcpyHostToDevice) );
    
    curandState *devStates;
    gpuErrchk( cudaMalloc((void **)&devStates, simparams->nTimes * sizeof(curandState)) );

    get_chop_pwv_rng<<<gridSize1D, blockSize1D>>>(dpointing_center, dPWV_screen, PWV_out, d_flagout, d_azout, d_elout, devStates);

    float freq;    // Frequency, used for initialising background sources.

    // Allocate and copy blackbodies
    float *I_atm = new float[source->nf];
    float *I_gnd = new float[source->nf];
    float *I_tel = new float[source->nf];
    float *I_CMB = new float[source->nf];

    for(int j=0; j<source->nf; j++)
    {
        freq = source->f0 + source->df * j;
        
        I_atm[j] = getPlanck(atmosphere->Tatm, freq); 
        I_gnd[j] = getPlanck(telescope->Tgnd, freq); 
        I_tel[j] = getPlanck(telescope->Ttel, freq);
        I_CMB[j] = getPlanck(2.725, freq);
    }
   
    gpuErrchk( cudaFree(dPWV_screen) );


    float *dI_atm, *dI_gnd, *dI_tel, *dI_CMB;
    gpuErrchk( cudaMalloc((void**)&dI_atm, source->nf * sizeof(float)) );
    gpuErrchk( cudaMalloc((void**)&dI_gnd, source->nf * sizeof(float)) );
    gpuErrchk( cudaMalloc((void**)&dI_tel, source->nf * sizeof(float)) );
    gpuErrchk( cudaMalloc((void**)&dI_CMB, source->nf * sizeof(float)) );

    gpuErrchk( cudaMemcpy(dI_atm, I_atm, source->nf * sizeof(float), cudaMemcpyHostToDevice) );
    gpuErrchk( cudaMemcpy(dI_gnd, I_gnd, source->nf * sizeof(float), cudaMemcpyHostToDevice) );
    gpuErrchk( cudaMemcpy(dI_tel, I_tel, source->nf * sizeof(float), cudaMemcpyHostToDevice) );
    gpuErrchk( cudaMemcpy(dI_CMB, I_CMB, source->nf * sizeof(float), cudaMemcpyHostToDevice) );
    
    
    // Allocate and copy telescope arrays
    float *deta_ap_ON, *deta_ap_OFF;
    gpuErrchk( cudaMalloc((void**)&deta_ap_ON, source->nf * sizeof(float)) );
    gpuErrchk( cudaMalloc((void**)&deta_ap_OFF, source->nf * sizeof(float)) );
    gpuErrchk( cudaMemcpy(deta_ap_ON, telescope->eta_ap_ON, source->nf * sizeof(float), cudaMemcpyHostToDevice) );
    gpuErrchk( cudaMemcpy(deta_ap_OFF, telescope->eta_ap_OFF, source->nf * sizeof(float), cudaMemcpyHostToDevice) );
    
    gpuErrchk( cudaBindTexture((size_t)0, tex_eta_ap_ON, deta_ap_ON, source->nf * sizeof(float)) );
    gpuErrchk( cudaBindTexture((size_t)0, tex_eta_ap_OFF, deta_ap_OFF, source->nf * sizeof(float)) );

    // Allocate and copy atmosphere arrays
    float *deta_atm;
    int neta_atm = atmosphere->nf * atmosphere->nPWV;
    
    gpuErrchk( cudaMalloc((void**)&deta_atm, neta_atm * sizeof(float)) );
    gpuErrchk( cudaMemcpy(deta_atm, atmosphere->eta_atm, neta_atm * sizeof(float), cudaMemcpyHostToDevice) );

    // Allocate and copy instrument arrays
    float *dfilterbank;
    int nfilterbank = source->nf * instrument->nfreqs_filt;
    gpuErrchk( cudaMalloc((void**)&dfilterbank, nfilterbank * sizeof(float)) );
    
    gpuErrchk( cudaMemcpy(dfilterbank, instrument->filterbank, nfilterbank * sizeof(float), cudaMemcpyHostToDevice) );

    gpuErrchk( cudaBindTexture((size_t)0, tex_filterbank, dfilterbank, nfilterbank * sizeof(float)) );
    
    timer.stop();

    output->t_diag[0] = timer.get();

    gpuErrchk( cudaMalloc((void**)&d_I_nu, source->nI_nu * sizeof(float)) );
    gpuErrchk( cudaMemcpy(d_I_nu, source->I_nu, source->nI_nu * sizeof(float), cudaMemcpyHostToDevice) );

    // CALL TO MAIN SIMULATION KERNEL
    timer.start();

    if(telescope->scantype == 0) {
        
        nShared = 3 * nThreads2D; // Three pointings
        calcPowerNEP_PS<<<gridSize2D, blockSize2D, nShared*sizeof(float)>>>(dI_atm, dI_gnd, dI_tel, dI_CMB,
            d_sigout, d_nepout, d_flagout,
            PWV_out, deta_atm, d_I_nu);
    }

    else {
        runSimulation<<<gridSize2D, blockSize2D>>>(dI_atm, dI_gnd, dI_tel, dI_CMB,
                d_sigout, d_nepout, d_azout, d_elout, d_flagout,
                PWV_out, deta_atm, d_I_nu);
    }
    
    gpuErrchk( cudaDeviceSynchronize() );
    timer.stop();
    //std::string name = "pwv";
    //writeArray<float>(nepout, instrument->nfreqs_filt, name);
    
    gpuErrchk( cudaFree(dI_atm) );
    gpuErrchk( cudaFree(dI_gnd) );
    gpuErrchk( cudaFree(dI_tel) );
    gpuErrchk( cudaFree(dI_CMB) );
    gpuErrchk( cudaFree(d_I_nu) );

    
    getSignal<<<gridSize1D, blockSize1D>>>(d_sigout, d_nepout, devStates);

    gpuErrchk( cudaDeviceSynchronize() );

    output->t_diag[1] = timer.get();
    
    timer.start();

    gpuErrchk( cudaMemcpy(output->signal, d_sigout, nffnt * sizeof(float), cudaMemcpyDeviceToHost) );

    gpuErrchk( cudaMemcpy(output->Az, d_azout, simparams->nTimes * sizeof(int), cudaMemcpyDeviceToHost) );
    gpuErrchk( cudaMemcpy(output->El, d_elout, simparams->nTimes * sizeof(int), cudaMemcpyDeviceToHost) );
    
    gpuErrchk( cudaMemcpy(output->flag, d_flagout, simparams->nTimes * sizeof(int), cudaMemcpyDeviceToHost) );

    gpuErrchk( cudaDeviceReset() );

    delete[] I_atm;
    delete[] I_gnd;
    delete[] I_tel;
    delete[] I_CMB;
    
    timer.stop();
    output->t_diag[2] = timer.get();
}

