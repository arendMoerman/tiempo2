#include "InterfaceCUDA.h"

/*! \file Kernels.cu
    \brief Definitions of CUDA kernels for TiEMPO2.

    author: Arend Moerman
*/

__constant__ float const_effs[CEFFSSIZE];   // Contains constant efficiencies:chain, gnd, mir, pb 
__constant__ float cPI;                     // Pi
__constant__ float cCL;                     // speed of light
__constant__ float cHP;                     // Planck constant
__constant__ float cKB;                     // Boltzmann constant
__constant__ float cdt;                     // Timestep

__constant__ float cfreq_chop;              // Chopping frequency
__constant__ float cfreq_nod;               // Nodding frequency
__constant__ float cdAz_chop;               // Chopping throw
__constant__ float ct0;                     // Starting time

__constant__ float ch_column;               // Column height
__constant__ float cv_wind;                 // Windspeed
__constant__ int cnx;                       // Number of x points in screen
__constant__ int cny;                       // Number of y points in screen

__constant__ float cdelta;                  // Bandgap energy of MKID
__constant__ float cfreq_sample;            // Sampling frequency of readout

__constant__ int cnt;                       // Number of time evals
__constant__ int cnf_filt;                  // Number of filter freqs
__constant__ int cnf_src;                   // Number of source frequencies
__constant__ int cnf_atm;                   // Number of atmosphere frequencies
__constant__ int cnPWV_atm;                 // Number of atmosphere PWV values

__constant__ int cnAz;                      // Number of az points per freq slice
__constant__ int cnEl;                      // Number of el points per freq slice

__constant__ int cOFF_empty;                // Interpolate on source or no source in OFF position

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
 __host__ std::array<dim3, 2> initCUDA(CuInstrument *instrument, CuTelescope *telescope, CuSimParams *simparams, CuSource *source, CuAtmosphere *atmosphere, int nThreads)
 {
     int nBlocks = ceil(simparams->nTimes / nThreads);

     // Calculate nr of blocks per grid and nr of threads per block
     dim3 nrb(nBlocks); dim3 nrt(nThreads);

     float PI = 3.1415926; /* pi */
     float CL = 2.9979246e8; // m s^-1
     float HP = 6.62607015e-34;
     float KB = 1.380649e-23;

     // Pack constant array
     float _con[CEFFSSIZE] = {instrument->eta_inst * telescope->eta_fwd * telescope->eta_mir * 0.5,
         instrument->eta_inst * (1 - telescope->eta_fwd) * telescope->eta_mir * 0.5,
         instrument->eta_inst * (1 - telescope->eta_mir) * 0.5, instrument->eta_pb};

     float dt = 1. / instrument->freq_sample;
     
         // Copy constant array to Device constant memory
     gpuErrchk( cudaMemcpyToSymbol(const_effs, &_con, CEFFSSIZE * sizeof(float)) );
     gpuErrchk( cudaMemcpyToSymbol(cPI, &PI, sizeof(float)) );
     gpuErrchk( cudaMemcpyToSymbol(cCL, &CL, sizeof(float)) );
     gpuErrchk( cudaMemcpyToSymbol(cHP, &HP, sizeof(float)) );
     gpuErrchk( cudaMemcpyToSymbol(cKB, &KB, sizeof(float)) );
     gpuErrchk( cudaMemcpyToSymbol(cdt, &dt, sizeof(float)) );
     gpuErrchk( cudaMemcpyToSymbol(cfreq_chop, &(telescope->freq_chop), sizeof(float)) );
     gpuErrchk( cudaMemcpyToSymbol(cfreq_nod, &(telescope->freq_nod), sizeof(float)) );
     gpuErrchk( cudaMemcpyToSymbol(cfreq_sample, &(instrument->freq_sample), sizeof(float)) );
     gpuErrchk( cudaMemcpyToSymbol(cdAz_chop, &(telescope->dAz_chop), sizeof(float)) );
     gpuErrchk( cudaMemcpyToSymbol(cdelta, &(instrument->delta), sizeof(float)) );
     gpuErrchk( cudaMemcpyToSymbol(ct0, &(simparams->t0), sizeof(float)) );
     gpuErrchk( cudaMemcpyToSymbol(ch_column, &(atmosphere->h_column), sizeof(float)) );
     gpuErrchk( cudaMemcpyToSymbol(cv_wind, &(atmosphere->v_wind), sizeof(float)) );
    

     gpuErrchk( cudaMemcpyToSymbol(cnt, &(simparams->nTimes), sizeof(int)) );
     gpuErrchk( cudaMemcpyToSymbol(cnf_filt, &(instrument->nfreqs_filt), sizeof(int)) );
     gpuErrchk( cudaMemcpyToSymbol(cnf_src, &(source->nfreqs_src), sizeof(int)) );
     gpuErrchk( cudaMemcpyToSymbol(cnf_atm, &(atmosphere->nfreqs_atm), sizeof(int)) );
     gpuErrchk( cudaMemcpyToSymbol(cnAz, &(source->nAz), sizeof(int)) );
     gpuErrchk( cudaMemcpyToSymbol(cnEl, &(source->nEl), sizeof(int)) );
     gpuErrchk( cudaMemcpyToSymbol(cnPWV_atm, &(atmosphere->nPWV_atm), sizeof(int)) );
     gpuErrchk( cudaMemcpyToSymbol(cnx, &(atmosphere->nx), sizeof(int)) );
     gpuErrchk( cudaMemcpyToSymbol(cny, &(atmosphere->ny), sizeof(int)) );
     
     gpuErrchk( cudaMemcpyToSymbol(cOFF_empty, &(simparams->OFF_empty), sizeof(int)) );

     std::array<dim3, 2> BT;
     BT[0] = nrb;
     BT[1] = nrt;

     return BT;
}

/**
  Initialise an array of random states for pseudorandom number generation.
  Given a seed, each thread is assigned a sequence based on its index.
  Snippet taken from the CUDA handbook.

  @param state Array of curand states. Should be initialised and sized to total number of threads in grid.
  @param seed Integer describing the seed of the generator.
 */
__global__ void setup_kernel(curandState *state, unsigned long long int seed = 0) {
    if (!seed) {
        seed = clock64();
    }
    
    int idx = blockDim.x*blockIdx.x + threadIdx.x;
    curand_init(seed, idx, 0, &state[idx]);
}

/**
  Main simulation kernel. This is where the magic happens.

  @param I_atm Array containing blackbody intensity of atmosphere, in SI units.
  @param I_gnd Array containing blackbody intensity of ground, in SI units.
  @param I_tel Array containing blackbody intensity of telescope, in SI units.
  @param sigout Array for storing output power, for each channel, for each time, in SI units.
  @param flagout Array for storing wether beam is in chop A or B, in nod AB or BA.
  @param freqs_src Array containing bin frequencies, in Hertz.
  @param azsrc Array containing Azimuth co-ordinates of source, in degrees.
  @param elsrc Array containing Elevation co-ordinates of source, in degrees.
  @param eta_ap Array containing aperture efficiencies, for each bin frequency.
  @param x_atm Array containing x co-ordinates of the atmosphere screen, in meters.
  @param y_atm Array containing y co-ordinates of the atmosphere screen, in meters.
  @param PWV_screen Array containing PWV value of atmosphere, over the range described by x_atm and y_atm, in millimeters.
  @param freqs_atm Array containing frequencies over which to interpolate atmospheric transmission, in Hertz.
  @param PWV_atm Array containing PWV values over which to interpolate atmospheric transmission, in millimeters.
  @param eta_atm Array with transmission parameters as fuiunction of freqs_atm and PWV_atm.
  @param filterbank Array containing filterbank of instrument.
  @param source Array containing source intensity, as function of azsrc, elsrc and freqs_src, in SI units.
  @param state Array with states for drawing random Gaussian values for noise calculations.
 */
__global__ void runSimulation(float *I_atm, float *I_gnd, float *I_tel,
        float *sigout, float *azout, float *elout, int *flagout,
        float *freqs_src, float *azsrc, float *elsrc, float *eta_ap,
        float *x_atm, float *y_atm, float *PWV_screen, float *freqs_atm,
        float *PWV_atm, float *eta_atm, float *filterbank, float *source, curandState *state) {
    
    int idx = blockDim.x*blockIdx.x + threadIdx.x;

    if (idx < cnt)
    {
        float t_real; // Real time.
        float t_start; // Time from start of observation.
        float PWV_Gauss_interp; // Interpolated PWV of Gaussian smoothed screen.
        float eta_atm_interp; // Interpolated eta_atm, over frequency and PWV
        float freq; // Bin frequency
        float I_nu; // Specific intensity of source.
        float sigma_k; // Noise per channel.
        float eta_kj; // Filter efficiency for bin j, at channel k
        float sqrt_samp = sqrtf(0.5 * cfreq_sample); // Constant term needed for noise calculation
        int n_chop, n_nod, start_slice;
        int position;

        float dfreq = freqs_src[1] - freqs_src[0];
    
        curandState localState = state[idx];
        bool chop_flag;

        float is_in_lower_half;
        int nod_flag;

        AzEl center;

        center.Az = 0;
        center.El = 0;

        AzEl pointing;
        xy_atm point_atm;

        t_real = idx * cdt + ct0;
        t_start = idx * cdt;

        n_chop = floor(t_start * cfreq_chop);
        n_nod = floor(t_start * cfreq_nod);
        
        chop_flag = (n_chop % 2 != 0); // If even (false), ON. Odd (true), OFF.
        nod_flag = -1 + 2 * (n_nod % 2 != 0); // If even (false), AB. Odd (true), BA.
        
        is_in_lower_half = (t_start - n_nod / cfreq_nod) - (1 / cfreq_nod / 2);
        sgn(is_in_lower_half, position);
        position *= nod_flag;
        
        scanPoint(&center, &pointing, chop_flag, position * cdAz_chop);
        flagout[idx] = chop_flag * position + (1 - chop_flag) * (1 - position); 
        
        // STORAGE: Add current pointing to output array
        azout[idx] = pointing.Az;
        elout[idx] = pointing.El;

        convertAnglesToSpatialAtm(&pointing, &point_atm, ch_column);

        // Add wind to this - currently only along x-axis and pretty manual
        point_atm.xAz = point_atm.xAz + cv_wind * t_real;

        // Interpolate on PWV_Gauss
        PWV_Gauss_interp = interpValue(point_atm.xAz, point_atm.yEl, 
                x_atm, y_atm, cnx, cny, PWV_screen, 0);

        float* PSD_nu = new float[cnf_src];
        int chop_flag_sgn = abs(chop_flag);
        
        for(int j=0; j<cnf_src; j++)
        {   
            freq = freqs_src[j];
            eta_atm_interp = interpValue(freq, PWV_Gauss_interp, 
                    freqs_atm, PWV_atm, cnf_atm, cnPWV_atm, eta_atm, 0);

            start_slice = cnAz * cnEl * j;
            I_nu = interpValue(pointing.Az, pointing.El, 
                azsrc, elsrc, cnAz, cnEl, source, start_slice);
            
            //if (cOFF_empty && chop_flag_sgn) {

            //    I_nu = 0.;
            //}

            //else {
            //    I_nu = interpValue(pointing.Az, pointing.El, 
            //        azsrc, elsrc, cnAz, cnEl, source, start_slice);
            //}
            //if (idx == 0){
            //    printf("%d %d %.12e\n", cOFF_empty, chop_flag_sgn, I_nu);
            //}
            PSD_nu[j] = (eta_ap[j] * eta_atm_interp * const_effs[0] * I_nu
                + const_effs[0] * (1 - eta_atm_interp) * I_atm[j] 
                + const_effs[1] * I_gnd[j] 
                + const_effs[2] * I_tel[j]) * cCL*cCL / (freq*freq);
        }
        
        // In this loop, calculate P_k, NEP_k and noise
        for(int k=0; k<cnf_filt; k++) {
            float P_k = 0.; // Initialise each channel to zero, for each timestep
            float NEP_accum = 0.;

            // Can loop over bins again, cheap operations this time
            for(int j=0; j<cnf_src; j++) {   
                freq = freqs_src[j];
                eta_kj = filterbank[k*cnf_src + j];
                
                //printf("%.12e %d hallooooo\n", PSD_nu[j], j); 
                NEP_accum += PSD_nu[j] * eta_kj * (cHP * freq + PSD_nu[j] * eta_kj + 2 * cdelta / const_effs[3]);
                P_k += PSD_nu[j] * eta_kj;
            }

            sigma_k = sqrt(2 * NEP_accum * dfreq) * sqrt_samp;
            P_k *= dfreq;

            P_k += sigma_k * curand_normal(&localState);
            state[idx] = localState;
           
            // STORAGE: Add signal to signal array in output
            sigout[idx * cnf_filt + k] = P_k; 

        }
        delete[] PSD_nu;
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

    Timer timer;

    int nThreads = 256;
    int totalThreads = ceil(simparams->nTimes / nThreads) * nThreads;

    timer.start();

    
    std::array<dim3, 2> BT = initCUDA(instrument, telescope, simparams, source, atmosphere, nThreads);
    float freq;    // Frequency, used for initialising background sources.

    // Allocate and copy blackbodies
    float *I_atm = new float[source->nfreqs_src];
    float *I_gnd = new float[source->nfreqs_src];
    float *I_tel = new float[source->nfreqs_src];

    for(int j=0; j<source->nfreqs_src; j++)
    {
        freq = source->freqs_src[j];
        
        I_atm[j] = getPlanck(atmosphere->Tatm, freq); 
        I_gnd[j] = getPlanck(telescope->Tgnd, freq); 
        I_tel[j] = getPlanck(telescope->Ttel, freq);
    }
   
    float *dI_atm, *dI_gnd, *dI_tel;
    gpuErrchk( cudaMalloc((void**)&dI_atm, source->nfreqs_src * sizeof(float)) );
    gpuErrchk( cudaMalloc((void**)&dI_gnd, source->nfreqs_src * sizeof(float)) );
    gpuErrchk( cudaMalloc((void**)&dI_tel, source->nfreqs_src * sizeof(float)) );

    gpuErrchk( cudaMemcpy(dI_atm, I_atm, source->nfreqs_src * sizeof(float), cudaMemcpyHostToDevice) );
    gpuErrchk( cudaMemcpy(dI_gnd, I_gnd, source->nfreqs_src * sizeof(float), cudaMemcpyHostToDevice) );
    gpuErrchk( cudaMemcpy(dI_tel, I_tel, source->nfreqs_src * sizeof(float), cudaMemcpyHostToDevice) );
    
    // Allocate and copy source arrays: freqs, az, el and I_nu
    float *dAz_src, *dEl_src, *dfreqs_src, *dI_nu;
    int nI_nu = source->nAz * source->nEl * source->nfreqs_src;

    gpuErrchk( cudaMalloc((void**)&dAz_src, source->nAz * sizeof(float)) );
    gpuErrchk( cudaMalloc((void**)&dEl_src, source->nEl * sizeof(float)) );
    gpuErrchk( cudaMalloc((void**)&dfreqs_src, source->nfreqs_src * sizeof(float)) );
    gpuErrchk( cudaMalloc((void**)&dI_nu, nI_nu * sizeof(float)) );
    
    gpuErrchk( cudaMemcpy(dAz_src, source->Az, source->nAz * sizeof(float), cudaMemcpyHostToDevice) );
    gpuErrchk( cudaMemcpy(dEl_src, source->El, source->nEl * sizeof(float), cudaMemcpyHostToDevice) );
    gpuErrchk( cudaMemcpy(dfreqs_src, source->freqs_src, source->nfreqs_src * sizeof(float), cudaMemcpyHostToDevice) );
    gpuErrchk( cudaMemcpy(dI_nu, source->I_nu, nI_nu * sizeof(float), cudaMemcpyHostToDevice) );
    
    // Allocate and copy telescope arrays
    float *deta_ap;
    gpuErrchk( cudaMalloc((void**)&deta_ap, source->nfreqs_src * sizeof(float)) );
    gpuErrchk( cudaMemcpy(deta_ap, telescope->eta_ap, source->nfreqs_src * sizeof(float), cudaMemcpyHostToDevice) );

    // Allocate and copy atmosphere arrays
    float *dx_atm, *dy_atm, *dPWV_screen, *dfreqs_atm, *dPWV_atm, *deta_atm;
    int nPWV_screen = atmosphere->nx * atmosphere->ny;
    int neta_atm = atmosphere->nfreqs_atm * atmosphere->nPWV_atm;
    
    gpuErrchk( cudaMalloc((void**)&dx_atm, atmosphere->nx * sizeof(float)) );
    gpuErrchk( cudaMalloc((void**)&dy_atm, atmosphere->ny * sizeof(float)) );
    gpuErrchk( cudaMalloc((void**)&dPWV_screen, nPWV_screen * sizeof(float)) );
    gpuErrchk( cudaMalloc((void**)&dfreqs_atm, atmosphere->nfreqs_atm * sizeof(float)) );
    gpuErrchk( cudaMalloc((void**)&dPWV_atm, atmosphere->nPWV_atm * sizeof(float)) );
    gpuErrchk( cudaMalloc((void**)&deta_atm, neta_atm * sizeof(float)) );

    gpuErrchk( cudaMemcpy(dx_atm, atmosphere->x_atm, atmosphere->nx * sizeof(float), cudaMemcpyHostToDevice) );
    gpuErrchk( cudaMemcpy(dy_atm, atmosphere->y_atm, atmosphere->ny * sizeof(float), cudaMemcpyHostToDevice) );
    gpuErrchk( cudaMemcpy(dPWV_screen, atmosphere->PWV, nPWV_screen * sizeof(float), cudaMemcpyHostToDevice) );
    gpuErrchk( cudaMemcpy(dfreqs_atm, atmosphere->freqs_atm, atmosphere->nfreqs_atm * sizeof(float), cudaMemcpyHostToDevice) );
    gpuErrchk( cudaMemcpy(dPWV_atm, atmosphere->PWV_atm, atmosphere->nPWV_atm * sizeof(float), cudaMemcpyHostToDevice) );
    gpuErrchk( cudaMemcpy(deta_atm, atmosphere->eta_atm, neta_atm * sizeof(float), cudaMemcpyHostToDevice) );

    // Allocate and copy instrument arrays
    float *dfilterbank;
    int nfilterbank = source->nfreqs_src * instrument->nfreqs_filt;
    gpuErrchk( cudaMalloc((void**)&dfilterbank, nfilterbank * sizeof(float)) );
    
    gpuErrchk( cudaMemcpy(dfilterbank, instrument->filterbank, nfilterbank * sizeof(float), cudaMemcpyHostToDevice) );

    // Allocate output arrays
    float *sigout;
    gpuErrchk( cudaMalloc((void**)&sigout, (source->nfreqs_src * simparams->nTimes) * sizeof(float)) );
    float *azout, *elout;
    gpuErrchk( cudaMalloc((void**)&azout, simparams->nTimes * sizeof(float)) );
    gpuErrchk( cudaMalloc((void**)&elout, simparams->nTimes * sizeof(float)) );
    
    int *flagout;
    gpuErrchk( cudaMalloc((void**)&flagout, simparams->nTimes * sizeof(int)) );

    // Setup rng for noise calcs
    curandState *devStates;
    gpuErrchk( cudaMalloc((void **)&devStates, totalThreads * sizeof(curandState)) );
    setup_kernel<<<BT[0], BT[1]>>>(devStates);
    
    gpuErrchk( cudaDeviceSynchronize() );
    
    // Set total heap reservation to 128 mb
    gpuErrchk( cudaDeviceSetLimit(cudaLimitMallocHeapSize, 128*1024*1024) );
    
    timer.stop();
    output->t_diag[0] = timer.get();
    
    // CALL TO MAIN SIMULATION KERNEL
    timer.start();

    runSimulation<<<BT[0], BT[1]>>>(dI_atm, dI_gnd, dI_tel,
            sigout, azout, elout, flagout,
            dfreqs_src, dAz_src, dEl_src, deta_ap, dx_atm, dy_atm, dPWV_screen, dfreqs_atm,
            dPWV_atm, deta_atm, dfilterbank, dI_nu, devStates);

    gpuErrchk( cudaDeviceSynchronize() );
    
    timer.stop();
    output->t_diag[1] = timer.get();
    
    timer.start();

    gpuErrchk( cudaMemcpy(output->signal, sigout, (instrument->nfreqs_filt * simparams->nTimes) * sizeof(float), cudaMemcpyDeviceToHost) );

    gpuErrchk( cudaMemcpy(output->Az, azout, simparams->nTimes * sizeof(int), cudaMemcpyDeviceToHost) );
    gpuErrchk( cudaMemcpy(output->El, elout, simparams->nTimes * sizeof(int), cudaMemcpyDeviceToHost) );
    
    gpuErrchk( cudaMemcpy(output->flag, flagout, simparams->nTimes * sizeof(int), cudaMemcpyDeviceToHost) );

    gpuErrchk( cudaDeviceReset() );

    delete[] I_atm;
    delete[] I_gnd;
    delete[] I_tel;
    
    timer.stop();
    output->t_diag[2] = timer.get();
}

