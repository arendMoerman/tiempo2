#include "InterfaceCUDA.h"

/*! \file Kernels.cu
    \brief Definitions of CUDA kernels for TiEMPO2.

    author: Arend Moerman
*/

// OBSERVATION-INSTRUMENT PARAMETERS
__constant__ float const_effs[CEFFSSIZE];   // Contains constant efficiencies:chain, gnd, mir, pb 
__constant__ float cdt;                     // Timestep
__constant__ float cfreq_chop;              // Chopping frequency
__constant__ float cfreq_nod;               // Nodding frequency
__constant__ float cf_sample;            // Sampling frequency of readout
__constant__ float cdAz_chop;               // Chopping throw
__constant__ float cdelta;                  // Bandgap energy of MKID
__constant__ int cnt;                       // Number of time evals
__constant__ int cnf_ch;                    // Number of filter freqs
__constant__ int cchop_mode;                // What chopping scheme to use

// ATMOSPHERE PARAMETERS
__constant__ float ch_column;               // Column height
__constant__ float cv_wind;                 // Windspeed

// SCAN PARAMETERS
__constant__ int cscantype;
__constant__ float ccscEl0;
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
texture<float, cudaTextureType1D, cudaReadModeElementType> tex_I_atm;
texture<float, cudaTextureType1D, cudaReadModeElementType> tex_I_gnd;
texture<float, cudaTextureType1D, cudaReadModeElementType> tex_I_tel;
texture<float, cudaTextureType1D, cudaReadModeElementType> tex_I_CMB;

#define KB              1.380649E-23f
#define CL              2.9979246E8f
#define HP              6.62607015E-34f

#define NTHREADS1D      256
#define NTHREADS2DX     32
#define NTHREADS2DY     16

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
  Write a CUDA device array to a file, for debugging.

  @param array Pointer to device array of type T.
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
  Convert angle in arcseconds to degrees.

  @param ang angle in arcseconds.
 */
__host__ __device__ __inline__ float as2deg(float ang) {
    return ang / 3600;
}

/**
  Convert angle in degrees to radian.

  @param ang Angle in degrees.
 */
__host__ __device__ __inline__ float deg2rad(float ang) {
    return ang * 0.017453292;
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
    
    out->Az = center->Az + offset + cAx*sinf(cwx*t)*cosf(cwx*t + deg2rad(cphix)) + cAxmin*sinf(cwxmin*t)*cosf(cwxmin*t + deg2rad(cphix));
    out->El = center->El + cAy*sinf(cwy*t)*sinf(cwy*t + deg2rad(cphiy)) + cAymin*sinf(cwymin*t)*sinf(cwymin*t + deg2rad(cphiy)) - cAy;
}

/**
  Convert an Az-El co-ordinate to a projected x-y co-ordinate on the atmosphere.

  @param angles Az-El co-ordinate to convert.
  @param out Container for storing the calculated x-y point.
 */
__device__ __inline__ void convertAnglesToSpatialAtm(AzEl* angles, xy_atm* out) {
    
    float coord = tanf(deg2rad(angles->Az)) * ch_column;
    
    out->xAz = coord;
    coord = tanf(deg2rad(angles->El)) * ch_column;
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
  @param source CuSource object containing source definitions.
  @param atmosphere CuAtmosphere object containing atmosphere parameters.
  @param nTimes number of time evaluations in simulation.

  @return BT Array of two dim3 objects, containing number of blocks per grid and number of threads per block.
 */
__host__ void initCUDA(Instrument<float> *instrument, Telescope<float> *telescope, Source<float> *source, Atmosphere<float> *atmosphere, int nTimes) {
    // Pack constant array
    float _con[CEFFSSIZE] = {instrument->eta_inst * instrument->eta_misc * telescope->eta_fwd * telescope->eta_mir * 0.5,
        instrument->eta_inst * instrument->eta_misc * (1 - telescope->eta_fwd) * telescope->eta_mir * 0.5,
        instrument->eta_inst * instrument->eta_misc * (1 - telescope->eta_mir) * 0.5, 
        instrument->eta_pb};

    float dt = 1. / instrument->f_sample;
     
    // OBSERVATION-INSTRUMENT PARAMETERS
    gpuErrchk( cudaMemcpyToSymbol(const_effs, &_con, CEFFSSIZE * sizeof(float)) );
    gpuErrchk( cudaMemcpyToSymbol(cdt, &dt, sizeof(float)) );
    gpuErrchk( cudaMemcpyToSymbol(cfreq_chop, &(telescope->freq_chop), sizeof(float)) );
    gpuErrchk( cudaMemcpyToSymbol(cfreq_nod, &(telescope->freq_nod), sizeof(float)) );
    gpuErrchk( cudaMemcpyToSymbol(cf_sample, &(instrument->f_sample), sizeof(float)) );
    gpuErrchk( cudaMemcpyToSymbol(cdAz_chop, &(telescope->dAz_chop), sizeof(float)) );
    gpuErrchk( cudaMemcpyToSymbol(cdelta, &(instrument->delta), sizeof(float)) );
    gpuErrchk( cudaMemcpyToSymbol(cnt, &nTimes, sizeof(int)) );
    gpuErrchk( cudaMemcpyToSymbol(cnf_ch, &(instrument->nf_ch), sizeof(int)) );
    gpuErrchk( cudaMemcpyToSymbol(cchop_mode, &(telescope->chop_mode), sizeof(int)) );
    
    // ATMOSPHERE PARAMETERS
    gpuErrchk( cudaMemcpyToSymbol(ch_column, &(atmosphere->h_column), sizeof(float)) );
    gpuErrchk( cudaMemcpyToSymbol(cv_wind, &(atmosphere->v_wind), sizeof(float)) );

    // SCAN PARAMETERS
    float cscEl0 = 1. / sinf(deg2rad(telescope->El0));

    gpuErrchk( cudaMemcpyToSymbol(cscantype, &(telescope->scantype), sizeof(int)) );
    gpuErrchk( cudaMemcpyToSymbol(ccscEl0, &cscEl0, sizeof(float)) );
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
__global__ void get_chop_pwv_rng(ArrSpec<float> Az_spec, ArrSpec<float> El_spec, 
                                 ArrSpec<float> x_atm, ArrSpec<float> y_atm, 
                                 float *center, float *PWV_screen, float *PWV_out, 
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

        // INTEGERS
        int flag;           // Flag for storing chop position.

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
    
        interpValue(point_atm.xAz, point_atm.yEl,
                    &x_atm, &y_atm, PWV_screen, 0, _PWV_out);            
    
        __syncthreads();
        curand_init(seed, idx, 0, &state[idx]);
        azout[idx] = pointing.Az;
        elout[idx] = pointing.El;
        flagout[idx] = flag;
        
        PWV_out[idx] = _PWV_out;
    }
}

__device__ void commonJob(ArrSpec<float> *f_src, ArrSpec<float> *f_atm, ArrSpec<float> *PWV_atm, int idx, int idy, float *sigout, float *nepout,
        float *PWV_trace, float *eta_atm, float I_nu, int flag) {
        
    // FLOATS
    float eta_atm_interp;   // Interpolated eta_atm, over frequency and PWV
    float freq;             // Bin frequency
    float PSD_nu;           // Local variable for storing PSD.
    float eta_kj;           // Filter efficiency for bin j, at channel k.
    float PWV_tr;           // Local variable for storing PWV value at time.
    float eta_ap;           // Local variable for storing aperture efficiency
    float sigfactor;        // Factor for calculating power. Perform outside of channel loop for speed.
    float nepfactor1;       // Factor 1 for calculating NEP. Perform outside of channel loop for speed.
    float nepfactor2;       // Factor 2 for calculating NEP. Perform outside of channel loop for speed.

    // Reusable symbols for interpolation stuff - listed separately for readability
    PWV_tr = PWV_trace[idx];

    freq = f_src->start + f_src->step * idy;

    interpValue(PWV_tr, freq,
                PWV_atm, f_atm,
                eta_atm, 0, eta_atm_interp);

    if(flag == 0 or flag == -1) {
        eta_ap = tex1Dfetch(tex_eta_ap_ON, idy); 
    }

    else {
        eta_ap = tex1Dfetch(tex_eta_ap_OFF, idy);
    }

    eta_atm_interp = powf(eta_atm_interp, ccscEl0);

    PSD_nu = eta_ap * eta_atm_interp * const_effs[0] * I_nu
        + ( const_effs[0] * (1 - eta_atm_interp) * tex1Dfetch(tex_I_atm, idy)
        + const_effs[1] * tex1Dfetch(tex_I_gnd, idy)
        + const_effs[2] * tex1Dfetch(tex_I_tel, idy)) 
        * CL*CL / (freq*freq);

    sigfactor = PSD_nu * f_src->step;
    nepfactor1 = sigfactor * (HP * freq + 2 * cdelta / const_effs[3]);
    nepfactor2 = sigfactor * PSD_nu;

    #pragma unroll 
    for(int k=0; k<cnf_ch; k++) {
        eta_kj = tex1Dfetch( tex_filterbank, k*f_src->num + idy);
        atomicAdd(&sigout[k*cnt + idx], __fmul_rn(eta_kj, sigfactor)); 
        atomicAdd(&nepout[k*cnt + idx], __fmul_rn(eta_kj, __fmaf_rn(nepfactor2, eta_kj, nepfactor1))); 
    }
}

/**
  Calculate power and NEP in each channel.
  This kernel is optimised for strict single-pointing (i.e., no sky-chopping) observations.

  @param sigout Array for storing output power, for each channel, for each time, in SI units.
  @param nepout Array for storing output NEP, for each channel, for each time, in SI units.
  @param flagout Array for storing wether beam is in chop A or B, in nod AB or BA.
  @param PWV_trace Array containing PWV value of atmosphere as seen by telescope over observation, in millimeters.
  @param eta_atm Array with transmission parameters as function of PWV and frequency.
  @param source Array containing source intensity at three pointings, as function of frequency, in SI units.
 */
__global__ void calcPowerNEP_1(ArrSpec<float> f_src, ArrSpec<float> f_atm, ArrSpec<float> PWV_atm, float *sigout, float *nepout, int *flagout,
        float *PWV_trace, float *eta_atm, float *source) {
    
    int idx = blockIdx.x * blockDim.x + threadIdx.x; 
    int idy = blockIdx.y * blockDim.y + threadIdx.y; 

    if (idx < cnt && idy < f_src.num) {
        float I_nu;             // Specific intensity of source.

        // INTEGERS
        int flag = flagout[idx];

        I_nu = source[idy];

        commonJob(&f_src, &f_atm, &PWV_atm, idx, idy, sigout, nepout, PWV_trace, eta_atm, I_nu, flag);
    }
}

/**
  Calculate power and NEP in each channel.
  This kernel is optimised for observations involving two pointings, encountered in single-pointing ON-OFF chopping.
        
  @param sigout Array for storing output power, for each channel, for each time, in SI units.
  @param nepout Array for storing output NEP, for each channel, for each time, in SI units.
  @param flagout Array for storing wether beam is in chop A or B, in nod AB or BA.
  @param PWV_trace Array containing PWV value of atmosphere as seen by telescope over observation, in millimeters.
  @param eta_atm Array with transmission parameters as function of PWV and frequency.
  @param source Array containing source intensity at three pointings, as function of frequency, in SI units.
 */
__global__ void calcPowerNEP_2(ArrSpec<float> f_src, ArrSpec<float> f_atm, ArrSpec<float> PWV_atm, float *sigout, float *nepout, int *flagout,
        float *PWV_trace, float *eta_atm, float *source) {
    
    int idx = blockIdx.x * blockDim.x + threadIdx.x; 
    int idy = blockIdx.y * blockDim.y + threadIdx.y; 

    if (idx < cnt && idy < f_src.num) {
        float I_nu;             // Specific intensity of source.
        int idx_point;          // Index for pointing of source (0 is OFF-B, 1 is ON-AB, 2 is OFF-A)

        // INTEGERS
        int flag = flagout[idx];
        
        // Determine idx in source from chopping flag
        if(flag==0) {idx_point = 0;}
        else {idx_point = 1;}
    
        I_nu = source[idx_point*f_src.num + idy];

        commonJob(&f_src, &f_atm, &PWV_atm, idx, idy, sigout, nepout, PWV_trace, eta_atm, I_nu, flag);
    }
}

/**
  Calculate power and NEP in each channel.
  This kernel is optimised for observations involving three pointings, encountered in single-pointing ABBA chopping.
        
  @param sigout Array for storing output power, for each channel, for each time, in SI units.
  @param nepout Array for storing output NEP, for each channel, for each time, in SI units.
  @param flagout Array for storing wether beam is in chop A or B, in nod AB or BA.
  @param PWV_trace Array containing PWV value of atmosphere as seen by telescope over observation, in millimeters.
  @param eta_atm Array with transmission parameters as function of PWV and frequency.
  @param source Array containing source intensity at three pointings, as function of frequency, in SI units.
 */
__global__ void calcPowerNEP_3(ArrSpec<float> f_src, ArrSpec<float> f_atm, ArrSpec<float> PWV_atm, float *sigout, float *nepout, int *flagout,
        float *PWV_trace, float *eta_atm, float *source) {
    
    int idx = blockIdx.x * blockDim.x + threadIdx.x; 
    int idy = blockIdx.y * blockDim.y + threadIdx.y; 

    if (idx < cnt && idy < f_src.num) {
        float I_nu;             // Specific intensity of source.
        int idx_point;          // Index for pointing of source (0 is OFF-B, 1 is ON-AB, 2 is OFF-A)

        //printf("%d\n", idx);
        // INTEGERS
        int flag = flagout[idx];
        
        // Determine idx in source from chopping flag
        if(flag==0 or flag==2) {idx_point = 1;}
        else if(flag==1) {idx_point = 2;}
        else {idx_point = 0;}
    
        I_nu = source[idx_point*f_src.num + idy];

        commonJob(&f_src, &f_atm, &PWV_atm, idx, idy, sigout, nepout, PWV_trace, eta_atm, I_nu, flag);
    }
}

/**
  Main simulation kernel. This is where the magic happens.

  @param sigout Array for storing output power, for each channel, for each time, in SI units.
  @param nepout Array for storing output NEP, for each channel, for each time, in SI units.
  @param azout Array containing Azimuth coordinates as function of time.
  @param elout Array containing Elevation coordinates as function of time.
  @param flagout Array for storing wether beam is in chop A or B, in nod AB or BA.
  @param PWV_trace Array containing PWV value of atmosphere as seen by telescope over observation, in millimeters.
  @param eta_atm Array with transmission parameters as function of freqs_atm and PWV_atm.
  @param source Array containing source intensity, as function of azsrc, elsrc and freqs_src, in SI units.
 */
__global__ void calcPowerNEP(ArrSpec<float> f_src, ArrSpec<float> f_atm, ArrSpec<float> PWV_atm, ArrSpec<float> Az_src, ArrSpec<float> El_src, float *sigout, float *nepout, float *azout, float *elout, int *flagout,
        float *PWV_trace, float *eta_atm, float *source) {
    
    int idx = blockIdx.x * blockDim.x + threadIdx.x; 
    int idy = blockIdx.y * blockDim.y + threadIdx.y; 

    if (idx < cnt && idy < f_src.num) {
        float I_nu;             // Specific intensity of source.
        // Reusable symbols for interpolation stuff
        int x0y0, x1y0, x0y1, x1y1;
        float t, u;
            
        AzEl pointing;

        pointing.Az = azout[idx];
        pointing.El = elout[idx];

        int iAz = floorf((pointing.Az - Az_src.start) / Az_src.step);
        int iEl = floorf((pointing.El - El_src.start) / El_src.step);

        int flag = flagout[idx];

        float Az_src_max = Az_src.start + Az_src.step * (Az_src.num - 1);
        float El_src_max = El_src.start + El_src.step * (El_src.num - 1);

        bool offsource = ((pointing.Az < Az_src.start) or (pointing.Az > Az_src_max)) or 
                         ((pointing.El < El_src.start) or (pointing.El > El_src_max));

        if(offsource) {I_nu = tex1Dfetch(tex_I_CMB, idy);}
        
        else {
            x0y0 = f_src.num * (iAz + iEl * Az_src.num);
            x1y0 = f_src.num * (iAz + 1 + iEl * Az_src.num);
            x0y1 = f_src.num * (iAz + (iEl+1) * Az_src.num);
            x1y1 = f_src.num * (iAz + 1 + (iEl+1) * Az_src.num);
            
            t = (pointing.Az - (Az_src.start + Az_src.step*iAz)) / Az_src.step;
            u = (pointing.El - (El_src.start + El_src.step*iEl)) / El_src.step;
            
            I_nu = (1-t)*(1-u) * source[x0y0 + idy];
            I_nu += t*(1-u) * source[x1y0 + idy];
            I_nu += (1-t)*u * source[x0y1 + idy];
            I_nu += t*u * source[x1y1 + idy];
        }

        commonJob(&f_src, &f_atm, &PWV_atm, idx, idy, sigout, nepout, PWV_trace, eta_atm, I_nu, flag);
    }
}

/**
  Calculate the total photon noise in a filter channel.

  After calculating the noise std from the NEP, a random number from a Gaussian is drawn and added to the total power in a channel.
  Note that, because we do not need the NEP after this step, we replace the value with a random Gaussian. 
  This is necessary for the TLS noise calculation, which comes after.

  @param sigout Array for storing output power, for each channel, for each time, in SI units.
  @param nepout Array for storing output NEP, for each channel, for each time, in SI units.
  @param state Array with states for drawing random Gaussian values for noise calculations.
 */
__global__ void calcPhotonNoise(float *sigout, float *nepout, curandState *state) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x; 
    if (idx < cnt) {
        curandState localState = state[idx];
        float sqrt_samp = sqrtf(0.5 * cf_sample); // Constant term needed for noise calculation
        float sigma_k, P_k;


        for(int k=0; k<cnf_ch; k++) {
            sigma_k = sqrtf(2 * nepout[k*cnt + idx]) * sqrt_samp;
            P_k = sigma_k * curand_normal(&localState);

            state[idx] = localState;

            sigout[k*cnt + idx] += P_k;

            nepout[k*cnt + idx] = curand_normal(&localState);
            state[idx] = localState;
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
  @param output CuOutput object for storing simulation output.
  @param nTimes Number of time evaluations in simulation.
 */
void runTiEMPO2_CUDA(Instrument<float> *instrument, Telescope<float> *telescope, Atmosphere<float> *atmosphere, Source<float> *source, Output<float> *output, int nTimesTotal, char *outpath) {
    // FLOATS
    float *d_sigout;        // Device pointer for output power array
    float *d_nepout;        // Device pointer for output NEP array
    float *d_azout;         // Device pointer for output Azimuth array 
    float *d_elout;         // Device pointer for output Elevation array
    float *d_I_nu;          // Device pointer for source intensities
    
    // INTEGERS
    int *d_flagout;         // Device pointer for output chopping flags
    int nffnt;              // Number of filter frequencies times number of time evaluations
    int nf_src;             // Number of frequency points in source.
    int numSMs;             // Number of streaming multiprocessors on GPU
    int nBlocks1D;          // Number of 1D blocks, in terms of number of SMs
    int nBlocks2Dx;         // Number of 2D blocks, in terms of number of SMs, along the x-dimension
    int nBlocks2Dy;         // Number of 2D blocks, in terms of number of SMs, along the y-dimension

    // OTHER DECLARATIONS
    dim3 blockSize1D;       // Size of 1D block (same as nThreads1D, but dim3 type)
    dim3 gridSize1D;        // Number of 1D blocks per grid
    dim3 blockSize2D;       // Size of 2D block, along x and y dimension
    dim3 gridSize2D;        // Number of 2D blocks per grid
    Timer timer;            // Timer class for timing kernel invocations

    // ALLOCATE ARRAY SPECIFICATION COPIES
    struct ArrSpec<float> _f_spec = instrument->f_spec;
    struct ArrSpec<float> _Az_src = source->Az_spec;
    struct ArrSpec<float> _El_src = source->El_spec;
    
    struct ArrSpec<float> _f_atm;
    struct ArrSpec<float> _PWV_atm;
    float *eta_atm;


    readEtaATM<float, ArrSpec<float>>(&eta_atm, &_PWV_atm, &_f_atm);
    
    std::string str_path(atmosphere->path);
    std::string str_outpath(outpath);

    int *meta;
    readAtmMeta(&meta, str_path);

    // Calculate lengths of x and y of single screen
    float lx = meta[1]*atmosphere->dx;
    float ly = meta[2]*atmosphere->dy;
    float lx_av = lx - ly;
    float t_obs_av = lx_av / atmosphere->v_wind; // Max available time per screen

    float timeTotal = nTimesTotal / instrument->f_sample;

    int nJobs = ceil(timeTotal / t_obs_av);
    int nTimesScreen = floor(t_obs_av * instrument->f_sample); // If error, change ceil to floor

    struct ArrSpec<float> _x_atm;
    struct ArrSpec<float> _y_atm;

    _x_atm.start = -ly/2;
    _x_atm.step = atmosphere->dx;
    _x_atm.num = meta[1];
    
    _y_atm.start = -ly/2;
    _y_atm.step = atmosphere->dy;
    _y_atm.num = meta[2];

    // Initialize constant memory
    initCUDA(instrument, telescope, source, atmosphere, nTimesScreen); 

    
    nf_src = _f_spec.num; // Number of spectral points in source
    
    gpuErrchk( cudaDeviceGetAttribute(&numSMs, cudaDevAttrMultiProcessorCount, 0) );

    // TiEMPO2 prefers larger L1 cache over shared memory.
    gpuErrchk( cudaDeviceSetCacheConfig(cudaFuncCachePreferL1) );

    timer.start();
    

    float freq;    // Frequency, used for initialising background sources.

    // Allocate and copy blackbodies
    std::vector<float> I_atm(nf_src);
    std::vector<float> I_gnd(nf_src);
    std::vector<float> I_tel(nf_src);
    std::vector<float> I_CMB(nf_src);

    for(int j=0; j<nf_src; j++)
    {
        freq = _f_spec.start + _f_spec.step * j;
        
        I_atm[j] = getPlanck(atmosphere->Tatm, freq); 
        I_gnd[j] = getPlanck(telescope->Tgnd, freq); 
        I_tel[j] = getPlanck(telescope->Ttel, freq);
        I_CMB[j] = getPlanck(2.725, freq);
    }
    
    float *dI_atm, *dI_gnd, *dI_tel, *dI_CMB;
    
    gpuErrchk( cudaMalloc((void**)&dI_atm, nf_src * sizeof(float)) );
    gpuErrchk( cudaMalloc((void**)&dI_gnd, nf_src * sizeof(float)) );
    gpuErrchk( cudaMalloc((void**)&dI_tel, nf_src * sizeof(float)) );
    gpuErrchk( cudaMalloc((void**)&dI_CMB, nf_src * sizeof(float)) );

    gpuErrchk( cudaMemcpy(dI_atm, I_atm.data(), nf_src * sizeof(float), cudaMemcpyHostToDevice) );
    gpuErrchk( cudaMemcpy(dI_gnd, I_gnd.data(), nf_src * sizeof(float), cudaMemcpyHostToDevice) );
    gpuErrchk( cudaMemcpy(dI_tel, I_tel.data(), nf_src * sizeof(float), cudaMemcpyHostToDevice) );
    gpuErrchk( cudaMemcpy(dI_CMB, I_CMB.data(), nf_src * sizeof(float), cudaMemcpyHostToDevice) );
    
    gpuErrchk( cudaBindTexture((size_t)0, tex_I_atm, dI_atm, nf_src * sizeof(float)) );
    gpuErrchk( cudaBindTexture((size_t)0, tex_I_gnd, dI_gnd, nf_src * sizeof(float)) );
    gpuErrchk( cudaBindTexture((size_t)0, tex_I_tel, dI_tel, nf_src * sizeof(float)) );
    gpuErrchk( cudaBindTexture((size_t)0, tex_I_CMB, dI_CMB, nf_src * sizeof(float)) );
    
    
    // Allocate and copy telescope arrays
    float *deta_ap_ON, *deta_ap_OFF;
    gpuErrchk( cudaMalloc((void**)&deta_ap_ON, nf_src * sizeof(float)) );
    gpuErrchk( cudaMalloc((void**)&deta_ap_OFF, nf_src * sizeof(float)) );
    gpuErrchk( cudaMemcpy(deta_ap_ON, telescope->eta_ap_ON, nf_src * sizeof(float), cudaMemcpyHostToDevice) );
    gpuErrchk( cudaMemcpy(deta_ap_OFF, telescope->eta_ap_OFF, nf_src * sizeof(float), cudaMemcpyHostToDevice) );
    
    gpuErrchk( cudaBindTexture((size_t)0, tex_eta_ap_ON, deta_ap_ON, nf_src * sizeof(float)) );
    gpuErrchk( cudaBindTexture((size_t)0, tex_eta_ap_OFF, deta_ap_OFF, nf_src * sizeof(float)) );

    // Allocate and copy atmosphere arrays
    float *deta_atm;
    int neta_atm = _f_atm.num * _PWV_atm.num;
    
    gpuErrchk( cudaMalloc((void**)&deta_atm, neta_atm * sizeof(float)) );
    gpuErrchk( cudaMemcpy(deta_atm, eta_atm, neta_atm * sizeof(float), cudaMemcpyHostToDevice) );
    delete[] eta_atm;

    // Allocate and copy instrument arrays
    float *dfilterbank;
    int nfilterbank = nf_src * instrument->nf_ch;
    gpuErrchk( cudaMalloc((void**)&dfilterbank, nfilterbank * sizeof(float)) );
    
    gpuErrchk( cudaMemcpy(dfilterbank, instrument->filterbank, nfilterbank * sizeof(float), cudaMemcpyHostToDevice) );

    gpuErrchk( cudaBindTexture((size_t)0, tex_filterbank, dfilterbank, nfilterbank * sizeof(float)) );
    
    timer.stop();

    //output->t_diag[0] = timer.get();

    gpuErrchk( cudaMalloc((void**)&d_I_nu, source->nI_nu * sizeof(float)) );
    gpuErrchk( cudaMemcpy(d_I_nu, source->I_nu, source->nI_nu * sizeof(float), cudaMemcpyHostToDevice) );

    std::string datp;

    // Loop starts here
    printf("\033[92m");
    int idx_wrap = 0;
    int time_counter = 0;
    for(int idx=0; idx<nJobs; idx++) {
        if (idx_wrap == meta[0]) {
            idx_wrap = 0;
        }

        if (idx == (nJobs - 1)) {
            nTimesScreen = nTimesTotal - nTimesScreen*(nJobs-1);
        }
        time_counter += nTimesScreen;

        printf("*** Progress: %d / 100 ***\r", time_counter*100 / nTimesTotal);
        fflush(stdout);

        nffnt = instrument->nf_ch * nTimesScreen; // Number of elements in single-screen output.
        gpuErrchk( cudaMemcpyToSymbol(cnt, &nTimesScreen, sizeof(int)) );
        
        nBlocks1D = ceilf((float)nTimesScreen / NTHREADS1D / numSMs);
        blockSize1D = NTHREADS1D;
        gridSize1D = nBlocks1D*numSMs;
        nBlocks2Dx = ceilf((float)nTimesScreen / NTHREADS2DX / numSMs);
        nBlocks2Dy = ceilf((float)nf_src / NTHREADS2DY / numSMs);

        blockSize2D = dim3(NTHREADS2DX, NTHREADS2DY);
        gridSize2D = dim3(nBlocks2Dx * numSMs, nBlocks2Dy * numSMs);
    
        // Allocate output arrays
        gpuErrchk( cudaMalloc((void**)&d_sigout, nffnt * sizeof(float)) );
        gpuErrchk( cudaMalloc((void**)&d_nepout, nffnt * sizeof(float)) );
        gpuErrchk( cudaMalloc((void**)&d_azout, nTimesScreen * sizeof(float)) );
        gpuErrchk( cudaMalloc((void**)&d_elout, nTimesScreen * sizeof(float)) );
        gpuErrchk( cudaMalloc((void**)&d_flagout, nTimesScreen * sizeof(int)) );

        // Allocate PWV screen now, delete CUDA allocation after first kernel call
        float *PWV_screen;
        float *dPWV_screen;
        
        int nPWV_screen = _x_atm.num * _y_atm.num;
        
        float *PWV_out;
        gpuErrchk( cudaMalloc((void**)&PWV_out, nTimesScreen * sizeof(float)) );
        
        float pointing_center[2] = {0., 0.};
        float *dpointing_center;
        gpuErrchk( cudaMalloc((void**)&dpointing_center, 2 * sizeof(float)) );
        gpuErrchk( cudaMemcpy(dpointing_center, pointing_center, 2 * sizeof(float), cudaMemcpyHostToDevice) );
        
        curandState *devStates;
        gpuErrchk( cudaMalloc((void **)&devStates, nTimesScreen * sizeof(curandState)) );

        datp = std::to_string(idx_wrap) + ".datp";
        readAtmScreen<float, ArrSpec<float>>(&PWV_screen, &_x_atm, &_y_atm, str_path, datp);
        
        gpuErrchk( cudaMalloc((void**)&dPWV_screen, nPWV_screen * sizeof(float)) );
        gpuErrchk( cudaMemcpy(dPWV_screen, PWV_screen, nPWV_screen * sizeof(float), cudaMemcpyHostToDevice) );
        

        get_chop_pwv_rng<<<gridSize1D, blockSize1D>>>(_Az_src, _El_src, _x_atm, _y_atm, dpointing_center, dPWV_screen, PWV_out, d_flagout, d_azout, d_elout, devStates);
       
        gpuErrchk( cudaFree(dpointing_center) );
        gpuErrchk( cudaFree(dPWV_screen) );

        // CALL TO MAIN SIMULATION KERNEL
        timer.start();
        
        // SINGLE POINTING, NO CHOP
        if(telescope->scantype == 0 && telescope->chop_mode == 0) {
            calcPowerNEP_1<<<gridSize2D, blockSize2D>>>(_f_spec, _f_atm, _PWV_atm, d_sigout, d_nepout, d_flagout,
                PWV_out, deta_atm, d_I_nu);
        }
        
        // SINGLE POINTING, STRICT ON-OFF
        else if(telescope->scantype == 0 && telescope->chop_mode == 1) {
            calcPowerNEP_2<<<gridSize2D, blockSize2D>>>(_f_spec, _f_atm, _PWV_atm, d_sigout, d_nepout, d_flagout,
                PWV_out, deta_atm, d_I_nu);
        }


        // SINGLE POINTING, ABBA CHOPPING
        else if(telescope->scantype == 0 && telescope->chop_mode == 2) {
            calcPowerNEP_3<<<gridSize2D, blockSize2D>>>(_f_spec, _f_atm, _PWV_atm, d_sigout, d_nepout, d_flagout,
                PWV_out, deta_atm, d_I_nu);
        }

        else {
            calcPowerNEP<<<gridSize2D, blockSize2D>>>(_f_spec, _f_atm, _PWV_atm, _Az_src, _El_src, d_sigout, d_nepout, d_azout, d_elout, d_flagout,
                    PWV_out, deta_atm, d_I_nu);
        }
        
        gpuErrchk( cudaDeviceSynchronize() );

        gpuErrchk( cudaFree(PWV_out) );
        
        calcPhotonNoise<<<gridSize1D, blockSize1D>>>(d_sigout, d_nepout, devStates);

        gpuErrchk( cudaDeviceSynchronize() );
        timer.stop();
        
        gpuErrchk( cudaFree(devStates) );
        gpuErrchk( cudaFree(d_nepout) );

        //output->t_diag[1] = timer.get();
        

        //timer.start();
        // ALLOCATE STRINGS FOR WRITING OUTPUT
        std::string signame = std::to_string(idx) + "signal.out";
        std::string azname = std::to_string(idx) + "az.out";
        std::string elname = std::to_string(idx) + "el.out";
        std::string flagname = std::to_string(idx) + "flag.out";

        std::vector<float> sigout(nffnt);
        std::vector<float> azout(nTimesScreen);
        std::vector<float> elout(nTimesScreen);
        std::vector<int> flagout(nTimesScreen);

        gpuErrchk( cudaMemcpy(sigout.data(), d_sigout, nffnt * sizeof(float), cudaMemcpyDeviceToHost) );
        gpuErrchk( cudaMemcpy(azout.data(), d_azout, nTimesScreen * sizeof(float), cudaMemcpyDeviceToHost) );
        gpuErrchk( cudaMemcpy(elout.data(), d_elout, nTimesScreen * sizeof(float), cudaMemcpyDeviceToHost) );
        gpuErrchk( cudaMemcpy(flagout.data(), d_flagout, nTimesScreen * sizeof(int), cudaMemcpyDeviceToHost) );

        timer.start();

        write1DArray<float>(sigout, str_outpath, signame);
        write1DArray<float>(azout, str_outpath, azname);
        write1DArray<float>(elout, str_outpath, elname);
        write1DArray<int>(flagout, str_outpath, flagname);
        
        timer.stop();
        
        gpuErrchk( cudaFree(d_sigout) );
        gpuErrchk( cudaFree(d_azout) );
        gpuErrchk( cudaFree(d_elout) );
        gpuErrchk( cudaFree(d_flagout) );

        idx_wrap++;
    }
    gpuErrchk( cudaDeviceReset() );
    timer.stop();
    output->t_diag[2] = timer.get();
    printf("\033[0m\n");
}

