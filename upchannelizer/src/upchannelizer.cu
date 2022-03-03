#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <curand.h>
#include <assert.h>
#include <cublas_v2.h>
#include <time.h>
#include <sys/time.h>
#include <iostream>
#include <string.h>
#include <math.h>
#include <cuComplex.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <cufft.h>
#include <device_launch_parameters.h>
#include "upchannelizer.h"

using namespace std;

// Perform transpose on the data and convert to floats
__global__
void data_transpose(signed char* data_in, cuComplex* data_tra, int offset, int n_pol, int n_chan, int n_win, int n_samp);

// Convenience function for checking CUDA runtime API results
// can be wrapped around any runtime API call. No-op in release builds.
inline
cudaError_t checkCuda(cudaError_t result)
{
  if (result != cudaSuccess) {
    fprintf(stderr, "CUDA Runtime Error: %s\n", 
            cudaGetErrorString(result));
    assert(result == cudaSuccess);
  }
  return result;
}

signed char* d_data_char = NULL;
cuComplex* d_data_comp = NULL;
float* h_fft = NULL;

// Allocate memory to all arrays 
void init_FFT() {
	printf("Here In init_FFT()! \n");

	// Allocate memory for input data float type
	checkCuda(cudaMalloc((void **)&d_data_char, (N_INPUT) * sizeof(signed char)));
	printf("Here 1st cudaMalloc! \n");

	// Allocate memory for input data cuComplex type
	checkCuda(cudaMalloc((void **)&d_data_comp, (N_INPUT) * sizeof(cuComplex) / 2));
	printf("Here 2nd cudaMalloc! \n");

	checkCuda(cudaMallocHost((void **)&h_fft, (N_INPUT) * sizeof(float)));

	return;
}

// Perform transpose on the data and convert to floats
__global__
void data_transpose(signed char* data_in, cuComplex* data_tra, int offset, int n_pol, int n_chan, int n_win, int n_samp) {
	int a = threadIdx.x; // Antenna index
	int p = threadIdx.y; // Polarization index
	int c = blockIdx.y;  // Coarse channel index
	int w = blockIdx.x;  // Time window index
        int t = blockIdx.z;  // Time sample index

	// If the input data is not float e.g. signed char, just multiply it by '1.0f' to convert it to a float
	int h_in = data_in_idx(p, w, t, (c + offset), a, n_pol, n_win, n_samp, n_chan); // data_in_idx(p, t, (f + offset), a, nt, n_chan);
	int h_tr = data_tr_idx(a, p, w, (c + offset), t, n_pol, n_win, n_chan); // data_tr_idx(a, p, (f + offset), t, n_chan);

	data_tra[h_tr].x = data_in[2*h_in]*1.0f;
	data_tra[h_tr].y = data_in[2*h_in + 1]*1.0f;
	

	return;
}


// Perform FFT
void upchannelize(cuComplex* data_tra, int n_pol, int n_chan, int n_win, int n_samp){
        cufftHandle plan;

        // Number of branches to perform FFT on
        int n_branches = N_ANT*n_pol*n_chan*n_win;

	// Setup the cuFFT plan
    	cufftPlan1d(&plan, n_samp, CUFFT_C2C, n_branches);
    	
    	// Execute a complex-to-complex 1D FFT
    	cufftExecC2C(plan, (cufftComplex *)data_tra, (cufftComplex *)data_tra, CUFFT_FORWARD);
}

// Run FFT
float* run_FFT(signed char* data_in, int n_pol, int n_chan, int n_win, int n_samp) {

	cudaError_t err_code;

        // Total number of time samples
        int nt = n_win*n_samp;

	// Transpose kernel: Specify grid and block dimensions
	dim3 dimBlock_transpose(N_ANT, n_pol, 1);
	dim3 dimGrid_transpose(n_samp, n_chan, n_win);

	signed char* d_data_in = d_data_char;
	cuComplex* d_data_tra = d_data_comp;
	float* data_out = h_fft;

	//printf("Before cudaMemcpy(HtoD) coefficients! \n");
	// Copy input data from host to device
	checkCuda(cudaMemcpy(d_data_in, data_in, 2*N_ANT*n_pol*nt*n_chan*sizeof(signed char), cudaMemcpyHostToDevice));

        // Perform transpose on the data and convert to floats  
        data_transpose<<<dimGrid_transpose, dimBlock_transpose>>>(d_data_in, d_data_tra, 0, n_pol, n_chan, n_win, n_samp);
        err_code = cudaGetLastError();
	if (err_code != cudaSuccess) {
		printf("FFT: data_transpose() kernel Failed: %s\n", cudaGetErrorString(err_code));
	}

        // Upchannelize the data
        upchannelize(d_data_tra, n_pol, n_chan, n_win, n_samp);

        // Copy input data from device to host
        checkCuda(cudaMemcpy(data_out, (float *)d_data_tra, 2*N_ANT*n_pol*nt*n_chan*sizeof(float), cudaMemcpyDeviceToHost));

        return data_out;
}

// Generate simulated data
signed char* simulate_data(int n_pol, int n_chan, int nt) {
	signed char* data_sim;
	data_sim = (signed char*)calloc(N_INPUT, sizeof(signed char));

	/*
	'sim_flag' is a flag that indicates the kind of data that is simulated.
	sim_flag = 0 -> Ones
	sim_flag = 1 -> Repeating sequence of 1 to 64
	sim_flag = 2 -> Sequence of 1 to 64 placed in a particular bin (bin 6 for now)
	sim flag = 3 -> Simulated sine wave
	*/
	int sim_flag = 3;
	if (sim_flag == 0) {
		for (int i = 0; i < (N_INPUT / 2); i++) {
			if(i < (N_REAL_INPUT/2)){
				data_sim[2 * i] = 1;
			}else{
				data_sim[2 * i] = 0;
			}
		}
	}
	if (sim_flag == 1) {
		int tmp = 0;
		for (int p = 0; p < n_pol; p++) {
			for (int t = 0; t < nt; t++) {
				for (int f = 0; f < n_chan; f++) {
					for (int a = 0; a < N_ANT; a++) {
						if (tmp >= N_ANT) {
							tmp = 0;
						}
						tmp = (tmp + 1) % (N_ANT+1);
						if(a < N_REAL_ANT){
							data_sim[2 * data_in_idx(p, 0, t, f, a, n_pol, 1, nt, n_chan)] = tmp;
						}else{
							data_sim[2 * data_in_idx(p, 0, t, f, a, n_pol, 1, nt, n_chan)] = 0;
						}
					}
				}
			}
		}
	}
	if (sim_flag == 2) {
		int tmp = 0;
		for (int p = 0; p < n_pol; p++) {
			for (int t = 0; t < nt; t++) {
				for (int a = 0; a < N_ANT; a++) {
					if (tmp >= N_ANT) {
						tmp = 0;
					}
					tmp = (tmp + 1) % (N_ANT+1);
					if(a < N_REAL_ANT){
						data_sim[2 * data_in_idx(p, 0, t, 5, a, n_pol, 1, nt, n_chan)] = tmp;
						data_sim[2 * data_in_idx(p, 0, t, 2, a, n_pol, 1, nt, n_chan)] = tmp;
					}else{
						data_sim[2 * data_in_idx(p, 0, t, 5, a, n_pol, 1, nt, n_chan)] = 0;
						data_sim[2 * data_in_idx(p, 0, t, 2, a, n_pol, 1, nt, n_chan)] = 0;
					}
				}
			}
		}
	}
	if (sim_flag == 3) {
		float freq = 1e3; // Resonant frequency

                float tmp_max = 1.0;
		float tmp_min = -1.0;

		for (int t = 0; t < nt; t++) {
			for (int f = 0; f < n_chan; f++) {
				for (int a = 0; a < N_ANT; a++) {
					if(a < N_REAL_ANT){
						// Requantize from doubles/floats to signed chars with a range from -128 to 127
						// X polarization
						data_sim[2 * data_in_idx(0, 0, t, f, a, n_pol, 1, nt, n_chan)] = (signed char)((((cos(2 * PI * freq * t*0.001) - tmp_min)/(tmp_max-tmp_min)) - 0.5)*256);
						data_sim[2 * data_in_idx(0, 0, t, f, a, n_pol, 1, nt, n_chan) + 1] = 0;
						// Y polarization
						data_sim[2 * data_in_idx(1, 0, t, f, a, n_pol, 1, nt, n_chan)] = (signed char)((((2*cos(2 * PI * freq * t*0.001) - tmp_min)/(tmp_max-tmp_min)) - 0.5)*256);
						data_sim[2 * data_in_idx(1, 0, t, f, a, n_pol, 1, nt, n_chan) + 1] = 0;
					}else{
						// X polarization
						data_sim[2 * data_in_idx(0, 0, t, f, a, n_pol, 1, nt, n_chan)] = 0;
						data_sim[2 * data_in_idx(0, 0, t, f, a, n_pol, 1, nt, n_chan) + 1] = 0;
						// Y polarization
						data_sim[2 * data_in_idx(1, 0, t, f, a, n_pol, 1, nt, n_chan)] = 0;
						data_sim[2 * data_in_idx(1, 0, t, f, a, n_pol, 1, nt, n_chan) + 1] = 0; // Make this negative if a different polarization is tested
					}
				}
			}
		}
	}
	return data_sim;
}

// Free memory
void Cleanup_FFT() {
	// Free up GPU memory at the end of a program
	if (d_data_char != NULL) {
		cudaFree(d_data_char);
	}
	if (d_data_comp != NULL) {
		cudaFree(d_data_comp);
	}
}

//Comment out main() function when compiling for hpguppi
// <----Uncomment here if testing standalone code
// Test all of the kernels and functions, and write the output to
// a text file for analysis
int main() {
	printf("Here!\n");
        int n_pol = 2;
	// 1k mode
	//int n_chan = 16; 
        //int nt = 32768;
	// 4k mode
    	int n_chan = 64;
        int nt = 8192;
	// 32k mode
    	//int n_chan = 512;
        //int nt = 1024;

        int n_win = N_TIME_STI;
        int n_samp = nt/n_win;

	// Allocate memory to all arrays used by run_FFT() 
	init_FFT();

        printf("After init_FFT() \n");

	// Generate simulated data
	signed char* sim_data = simulate_data(n_pol, n_chan, nt);

        printf("After simulate_data() \n");

	// Allocate memory for output array
	float* output_data;

	printf("Here5!\n");

	float time_taken = 0;
	float fft_time = 0;
	int num_runs = 10;

	// Start timing beamformer computation //
	struct timespec tval_before, tval_after;

	for(int ii = 0; ii < num_runs; ii++){
		// Start timing beamformer computation //
		clock_gettime(CLOCK_MONOTONIC, &tval_before);

		// Run beamformer 
                output_data = run_FFT(sim_data, n_pol, n_chan, n_win, n_samp);
		//output_data = run_beamformer(sim_data, sim_coefficients, n_chan, nt);
		//run_beamformer(h_data, h_coeff, output_data);

		// Stop timing beamforming computation //
		clock_gettime(CLOCK_MONOTONIC, &tval_after);
		time_taken = (float)(tval_after.tv_sec - tval_before.tv_sec); //*1e6; // Time in seconds since epoch
		time_taken = time_taken + (float)(tval_after.tv_nsec - tval_before.tv_nsec)*1e-9; // Time in nanoseconds since 'tv_sec - start and end'
		fft_time += time_taken;
		//printf("Time taken: %f s\n", time_taken);
	}
	printf("Average FFT processing time: %f s\n", fft_time/num_runs);

	printf("Here6, FFT output: %f \n", output_data[0]);
	
	// Write data to text file for analysis
	char output_filename[128];

	printf("Here7!\n");

	strcpy(output_filename, "output_d_cufft.txt");

	printf("Here8!\n");

	FILE* output_file;

	printf("Here9!\n");

	output_file = fopen(output_filename, "w");

	printf("Here10!\n");

	for (int ii = 0; ii < ((N_INPUT*n_pol*n_chan*nt)/(N_POL*N_FREQ*N_TIME)); ii++) { // Write up to the size of the data corresponding to 1k, 4k or 32k mode
		//fprintf(output_file, "%c\n", output_data[ii]);
		fprintf(output_file, "%g\n", output_data[ii]);
	}

	printf("Here11!\n");

	fclose(output_file);

	printf("Closed output file.\n");

	//free(sim_data);
	printf("After freeing coefficients.\n");
	//free(output_data);	

	printf("Freed output array and unregistered arrays in pinned memory.\n");

	// Free up device memory
	//cudaFreeHost(h_data);
	//cudaFreeHost(h_coeff);
	Cleanup_FFT();

	printf("Here11!\n");

	return 0;
}
// <----Uncomment here if testing standalone code
