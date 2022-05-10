#include <stdio.h>
#include <stdlib.h>
//#include <complex.h>
#include <math.h>

#define MAX_THREADS (1024)  // Maximum number of threads in a block
#define N_POL (2) //2                     // Number of polarizations
#define N_TIME (131072) // (1024) // 8192 //16384 //1 // 8                   // Number of time samples
#define N_TIME_STI (8)
#define N_STI (N_TIME/N_TIME_STI)
#define N_STI_BLOC (32)
#define N_STREAMS (1)                     // Number of CUDA streams
//#define N_COARSE_FREQ 64 //32               // Number of coarse channels processed at a time
#define MAX_COARSE_FREQ (32) // (512)                 // Max number of coarse channels is the number of channels in 32k mode
#define N_FINE_FREQ (1) //16384               // Number of fine channels per coarse channel 2^14 = 16384
//#define N_FREQ (N_COARSE_FREQ*N_FINE_FREQ) // Number of frequency bins after second FFT.  Should actually be 2^14, but due to limited memory on my laptop, arbitrarily 10
#define N_FREQ (MAX_COARSE_FREQ*N_FINE_FREQ) // Number of frequency bins after second FFT.  Should actually be 2^14, but due to limited memory on my laptop, arbitrarily 10
//#define N_FREQ_STREAM (N_FREQ/N_STREAMS) // (N_COARSE_FREQ*N_FINE_FREQ)/N_STREAMS // Number of frequency bins after second FFT.  Should actually be 2^14, but due to limited memory on my laptop, arbitrarily 10
//#define N_FREQ_STREAM (N_COARSE_FREQ/N_STREAMS) // (N_COARSE_FREQ*N_FINE_FREQ)/N_STREAMS // Number of frequency bins after second FFT.  Should actually be 2^14, but due to limited memory on my laptop, arbitrarily 10
#define N_ANT (64) // 64                  // Number of possible antennas (64 is also computationally efficient since it is a multiple of 32 which is the size of a warp)
#define N_REAL_ANT (58)                   // Number of antennas transmitting data downstream
#define N_BEAM (64) // 64                 // Number of beams

// "2" for inphase and quadrature
#define N_INPUT       (unsigned long int)(2*N_POL*N_TIME*N_FREQ*N_ANT)                  // Size of input. Currently, same size as output
#define N_REAL_INPUT  (unsigned long int)(2*N_POL*N_TIME*N_FREQ*N_REAL_ANT)             // Size of input. Currently, same size as output
#define N_COEFF       (unsigned long int)(2*N_POL*N_ANT*N_BEAM*N_FREQ)                  // Size of beamformer coefficients
#define DELAY_POLYS   (unsigned long int)(2)                                            // Number of coefficients in polynomial
#define N_DELAYS      (unsigned long int)(DELAY_POLYS*N_ANT*N_BEAM)                     // Size of first order polynomial delay array
#define N_OUTPUT      (unsigned long int)(2*N_POL*N_BEAM*N_FREQ*N_TIME)                 // Size of beamformer output
#define N_FFT         (16384)                                                           // N-points of FFT in 32k mode
#define N_BF_POW      (unsigned long int)(N_BEAM*N_FREQ*N_FFT)                          // Size of beamformer output after abs()^2 and short time integration
// For cuFFT
#define RANK                (1)
//#define BATCH(Np,Nw,Nf)     (N_ANT)*(Np)*(Nw)*(Nf)
//#define BATCH(Np,Nf)        (N_ANT)*(Np)*(Nf)
#define BATCH(Np)           (N_ANT)*(Np)
#define ISTRIDE             (1)
#define IDIST(Nt)           (Nt)
#define OSTRIDE             (1)
#define ODIST(Nt)           (Nt)

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

//#define data_in_idx(p, t, w, c, a, Np, Nt, Nw, Nc)           ((p) + (Np)*(t) + (Nt)*(Np)*(w) + (Nw)*(Nt)*(Np)*(c) + (Nc)*(Nw)*(Nt)*(Np)*(a))
#define data_in_idx(p, c, a, t, w, Np, Nc, Na, Nt)           ((p) + (Np)*(c) + (Nc)*(Np)*(a) + (Na)*(Nc)*(Np)*(t) + (Nt)*(Na)*(Nc)*(Np)*(w))
#define data_tr_idx(t, a, p, c, w, Nt, Np, Nc)               ((t) + (Nt)*(a) + (N_ANT)*(Nt)*(p) + (Np)*(N_ANT)*(Nt)*(c) + (Nc)*(Np)*(N_ANT)*(Nt)*(w))
#define data_fft_out_idx(f, a, p, c, w, Nf, Np, Nc)          ((f) + (Nf)*(a) + (N_ANT)*(Nf)*(p) + (Np)*(N_ANT)*(Nf)*(c) + (Nc)*(Np)*(N_ANT)*(Nf)*(w))
// The "Nf" below is equal in value to "Nt*Nc" that is the dimension of "t" since this is the number of FFT points muliplied by the number of coarse channels
#define data_fftshift_idx(a, p, f, c, w, Np, Nf, Nc)         ((a) + (N_ANT)*(p) + (Np)*(N_ANT)*(f) + (Nf)*(Np)*(N_ANT)*(c) + (Nc)*(Nf)*(Np)*(N_ANT)*(w))
#define coeff_idx(a, p, b, f, Np, Nb)                        ((a) + (N_ANT)*(p) + (Np)*(N_ANT)*(b) + (Nb)*(Np)*(N_ANT)*(f))
//#define phase_idx(a, p, f, Np)                               ((a) + (N_ANT)*(p) + (Np)*(N_ANT)*(f))
//#define delay_idx(d, a, b, Na)                               ((d) + DELAY_POLYS*(a) + DELAY_POLYS*(Na)*(b)) // Should be correct indexing
#define cal_all_idx(a, p, f, Na, Np)                         ((a) + (Na)*(p) + (Np)*(Na)*(f))
#define delay_idx(a, b, t, Na, Nb)                           ((a) + (Na)*(b) + (Nb)*(Na)*(t))
#define coh_bf_idx(p, b, f, c, w, Np, Nb, Nc, Nf)            ((p) + (Np)*(b) + (Nb)*(Np)*(f) + (Nf)*(Nb)*(Np)*(c) + (Nc)*(Nf)*(Nb)*(Np)*(w))
#define pow_bf_idx(f, c, s, b, Nf, Nc, Ns)                   ((f) + (Nf)*(c) + (Nc)*(Nf)*(s) + (Ns)*(Nc)*(Nf)*(b)) // Changed to efficiently write each beam to a filterbank file

typedef struct complex_t{
	float re;
	float im;
}complex_t;

#ifdef __cplusplus
extern "C" {
#endif
void init_upchan_beamformer(); // Allocate memory to all arrays 
void set_to_zero_ubf(); // Set arrays to zero after a block is processed
signed char* simulate_data_ubf(int n_pol, int n_chan, int nt);
float* simulate_coefficients_ubf(int n_pol, int n_beam, int n_chan);
float* generate_coefficients_ubf(complex_t* phase_up, double* delay, int n, double* coarse_chan, int n_pol, int n_beam, int schan, int n_coarse, int subband_idx, uint64_t n_real_ant);
//void input_data_pin(signed char * data_in_pin);
//void output_data_pin(float * data_out_pin);
//void coeff_pin(float * data_coeff_pin);
//void unregister_data(void * data_unregister);
void Cleanup_beamformer();
void upchannelize(complex_t* data_tra, int n_pol, int n_chan, int n_samp); // Upchannelization
float* run_upchannelizer_beamformer(signed char* data_in, float* h_coefficient, int n_pol, int n_ant, int n_beam, int n_chan, int n_win, int n_time_int, int n_samp); // Run upchannelizer and beamformer
#ifdef __cplusplus
}
#endif
