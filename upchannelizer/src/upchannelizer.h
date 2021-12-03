#include <stdio.h>
#include <stdlib.h>
//#include <complex.h>
#include <math.h>


#define N_POL 2 //2                     // Number of polarizations
#define N_TIME 8192 //16384 //1 // 8                   // Number of time samples
#define N_STREAMS 1                     // Number of CUDA streams
#define N_COARSE_FREQ 64 //32               // Number of coarse channels processed at a time
#define N_FINE_FREQ 1 //16384               // Number of fine channels per coarse channel 2^14 = 16384
#define N_FREQ (N_COARSE_FREQ*N_FINE_FREQ) // Number of frequency bins after second FFT.  Should actually be 2^14, but due to limited memory on my laptop, arbitrarily 10
#define N_FREQ_STREAM N_FREQ/N_STREAMS // (N_COARSE_FREQ*N_FINE_FREQ)/N_STREAMS // Number of frequency bins after second FFT.  Should actually be 2^14, but due to limited memory on my laptop, arbitrarily 10
#define N_ANT 64 // 64                  // Number of possible antennas (64 is also computationally efficient since it is a multiple of 32 which is the size of a warp)
#define N_REAL_ANT 58                   // Number of antennas transmitting data downstream
#define N_BEAM 64 // 64                 // Number of beams
//#define N_POL_OUT 4 //2    // Number of output polarizations 

// "2" for inphase and quadrature
#define N_INPUT       (unsigned long int)(2*N_POL*N_TIME*N_FREQ*N_ANT)       // Size of input. Currently, same size as output
#define N_REAL_INPUT  (unsigned long int)(2*N_POL*N_TIME*N_FREQ*N_REAL_ANT)       // Size of input. Currently, same size as output
//#define N_COEFF  (unsigned long int)(2*N_ANT*N_POL*N_BEAM*N_FREQ)     // Size of beamformer coefficients
#define N_COEFF       (unsigned long int)(2*N_ANT*N_BEAM)                    // Size of beamformer coefficients
#define DELAY_POLYS   (unsigned long int)(2)                                 // Number of coefficients in polynomial
#define N_DELAYS      (unsigned long int)(DELAY_POLYS*N_ANT*N_BEAM)          // Size of first order polynomial delay array
#define N_OUTPUT      (unsigned long int)(2*N_POL*N_BEAM*N_FREQ*N_TIME)      // Size of beamformer output

#ifndef min
#define min(a,b) ((a < b) ? a : b)
#endif
#ifndef max
#define max(a,b) ((a > b) ? a : b)
#endif

#define PI 3.14159265

#ifdef __cplusplus
extern "C" {
#endif
void init_fft(); // Allocate memory to all arrays 
signed char* simulate_data();

void input_data_pin(signed char * data_in_pin);
void unregister_data(void * data_unregister);
void fftCleanup();

float* performFFT(signed char* data_in); // Run beamformer
#ifdef __cplusplus
}
#endif
