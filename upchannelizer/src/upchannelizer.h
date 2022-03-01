#include <stdio.h>
#include <stdlib.h>
//#include <complex.h>
#include <math.h>


#define N_POL 2 //2                     // Number of polarizations
#define N_TIME 1024 // 8192 //16384 //1 // 8                   // Number of time samples
#define N_TIME_STI 8
#define N_STI N_TIME/N_TIME_STI
#define N_STI_BLOC 32
#define N_STREAMS 1                     // Number of CUDA streams
//#define N_COARSE_FREQ 64 //32               // Number of coarse channels processed at a time
#define MAX_COARSE_FREQ 512                 // Max number of coarse channels is the number of channels in 32k mode
#define N_FINE_FREQ 1 //16384               // Number of fine channels per coarse channel 2^14 = 16384
//#define N_FREQ (N_COARSE_FREQ*N_FINE_FREQ) // Number of frequency bins after second FFT.  Should actually be 2^14, but due to limited memory on my laptop, arbitrarily 10
#define N_FREQ (MAX_COARSE_FREQ*N_FINE_FREQ) // Number of frequency bins after second FFT.  Should actually be 2^14, but due to limited memory on my laptop, arbitrarily 10
//#define N_FREQ_STREAM (N_FREQ/N_STREAMS) // (N_COARSE_FREQ*N_FINE_FREQ)/N_STREAMS // Number of frequency bins after second FFT.  Should actually be 2^14, but due to limited memory on my laptop, arbitrarily 10
//#define N_FREQ_STREAM (N_COARSE_FREQ/N_STREAMS) // (N_COARSE_FREQ*N_FINE_FREQ)/N_STREAMS // Number of frequency bins after second FFT.  Should actually be 2^14, but due to limited memory on my laptop, arbitrarily 10
#define N_ANT 64 // 64                  // Number of possible antennas (64 is also computationally efficient since it is a multiple of 32 which is the size of a warp)
#define N_REAL_ANT 58                   // Number of antennas transmitting data downstream
#define N_BEAM 64 // 64                 // Number of beams

// "2" for inphase and quadrature
#define N_INPUT       (unsigned long int)(2*N_POL*N_TIME*N_FREQ*N_ANT)       // Size of input. Currently, same size as output
#define N_REAL_INPUT  (unsigned long int)(2*N_POL*N_TIME*N_FREQ*N_REAL_ANT)  // Size of input. Currently, same size as output
#define N_COEFF       (unsigned long int)(2*N_POL*N_ANT*N_BEAM*N_FREQ)       // Size of beamformer coefficients
//#define N_PHASE       (unsigned long int)(2*N_POL*N_ANT)                     // Size of telstate phase solution array - Include frequency dimension from meta data
#define DELAY_POLYS   (unsigned long int)(2)                                 // Number of coefficients in polynomial
#define N_DELAYS      (unsigned long int)(DELAY_POLYS*N_ANT*N_BEAM)          // Size of first order polynomial delay array
#define N_OUTPUT      (unsigned long int)(2*N_POL*N_BEAM*N_FREQ*N_TIME)      // Size of beamformer output
//#define N_BF_POW      (unsigned long int)(N_BEAM*N_FREQ*N_TIME)              // Size of beamformer output after abs()^2
#define N_BF_POW      (unsigned long int)(N_BEAM*N_FREQ*N_STI)               // Size of beamformer output after abs()^2 and short time integration
// For cuFFT
#define BATCHES(Nf)     Nf*N_POL*N_ANT
//#define N_POINT(Nt)     Nt/N_TIME_STI
//#define ISTRIDE(Nf,Nw)  Nf*Nw*N_POL*N_ANT
//#define IDIST           1
//#define OSTRIDE(Nw)     Nw*N_POL*N_ANT
//#define ODIST           1   

#ifndef min
#define min(a,b) ((a < b) ? a : b)
#endif
#ifndef max
#define max(a,b) ((a > b) ? a : b)
#endif

#define PI 3.14159265

// p - polarization index
// t - time sample index
// w - time window index
// c - coarse channel frequency index
// f - fine channel frequency index
// a - antenna index
// b - beam index

#define data_in_idx(p, w, t, c, a, Np, Nw, Nt, Nc)     (p + Np*t + Nt*Np*w + Nw*Nt*Np*c + Nc*Nw*Nt*Np*a)
#define data_tr_idx(a, p, w, c, t, Np, Nw, Nc)         (a + N_ANT*p + Np*N_ANT*w + Nw*Np*N_ANT*c + Nc*Nw*Np*N_ANT*t)
// The "Nf" below is equal in value to "Nt*Nc" that is the dimension of "t" since this is the number of FFT points muliplied by the number of coarse channels
#define data_fft_idx(a, p, w, f, Np, Nw)               (a + N_ANT*p + Np*N_ANT*w + Nw*Np*N_ANT*f)

#ifdef __cplusplus
extern "C" {
#endif
void init_FFT(); // Allocate memory to all arrays 
void set_to_zero(); // Set arrays to zero after a block is processed
signed char* simulate_data(int n_chan, int nt);
void input_data_pin(signed char * data_in_pin);
void output_data_pin(float * data_out_pin);
void unregister_data(void * data_unregister);
void Cleanup_FFT();
cuComplex* run_FFT(signed char* data_in, int n_pol, int n_chan, int nt); // Run FFT
#ifdef __cplusplus
}
#endif
