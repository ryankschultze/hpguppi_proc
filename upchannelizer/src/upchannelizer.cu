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
#include <device_launch_parameters.h>
#include "coherent_beamformer_char_in.h"

using namespace std;

// Perform transpose on the data and convert to floats
__global__
void data_transpose(signed char* data_in, cuComplex* data_tra, int offset, int n_chan, int nt);

// Perform beamforming operation
__global__
void coherent_beamformer(cuComplex* input_data, float* coeff, float* output_data, int offset, int n_chan);

// Compute power of beamformer output (abs()^2)
__global__
void beamformer_power(float* bf_volt, float* bf_power, int offset, int n_chan, int nt);

// Compute power of beamformer output with STI (abs()^2)
__global__
void beamformer_power_sti(float* bf_volt, float* bf_power, int offset, int n_chan, int ns);

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

float* d_data_float = NULL;
signed char* d_data_char = NULL;
cuComplex* d_data_comp = NULL;
float* d_coeff = NULL;
float* d_coh_bf_pow = NULL;
float* h_bf_pow = NULL;

// Allocate memory to all arrays 
void init_beamformer() {
	printf("Here In init_beamformer()! \n");

	// Allocate memory for input data float type
	checkCuda(cudaMalloc((void **)&d_data_float, (N_INPUT) * sizeof(float)));
	checkCuda(cudaMalloc((void **)&d_data_char, (N_INPUT) * sizeof(signed char)));
	printf("Here 1st cudaMalloc! \n");

	// Allocate memory for input data cuComplex type
	checkCuda(cudaMalloc((void **)&d_data_comp, (N_INPUT) * sizeof(cuComplex) / 2));
	printf("Here 2nd cudaMalloc! \n");

	// Allocate memory for coefficients float type
	checkCuda(cudaMalloc((void **)&d_coeff, N_COEFF * sizeof(float)));
	printf("Here 3rd cudaMalloc! \n");

	// Allocate memory for output power of coherent beamformer
        checkCuda(cudaMalloc((void **)&d_coh_bf_pow, (N_BF_POW) * sizeof(float)));
	printf("Here 4th cudaMalloc! \n");

	checkCuda(cudaMallocHost((void **)&h_bf_pow, (N_BF_POW) * sizeof(float)));

	return;
}

// Set arrays to zero after a block is processed
void set_to_zero(){
	checkCuda(cudaMemset(d_data_float, 0, (N_INPUT) * sizeof(float)));
}

// Perform transpose on the data and convert to floats
__global__
void data_transpose(signed char* data_in, cuComplex* data_tra, int offset, int n_chan, int nt) {
	int a = threadIdx.x; // Antenna index
	int p = threadIdx.y; // Polarization index
	int f = blockIdx.y;  // Frequency index
	int t = blockIdx.x;  // Time sample index

	int n_freq_streams = n_chan/N_STREAMS;

	// If the input data is not float e.g. signed char, just multiply it by '1.0f' to convert it to a float
	if(f < n_freq_streams){
		int h_in = data_in_idx(p, t, (f + offset), a, nt, n_chan);
		int h_tr = data_tr_idx(a, p, (f + offset), t, n_chan);

		data_tra[h_tr].x = data_in[2*h_in]*1.0f;
		data_tra[h_tr].y = data_in[2*h_in + 1]*1.0f;
	}

	return;
}

// Perform beamforming operation
__global__
void coherent_beamformer(cuComplex* input_data, float* coeff, float* output_data, int offset, int n_chan) {
	int a = threadIdx.x; // Antenna index
	int f = blockIdx.y;  // Frequency index
	int t = blockIdx.x;  // Time sample index
	int b = blockIdx.z;  // Beam index

	__shared__ cuFloatComplex reduced_mul[N_ANT];

	int n_freq_streams = n_chan/N_STREAMS;

	if(f < n_freq_streams){
		for (int p = 0; p < N_POL; p++) { // Polarization index

	
			int i = data_tr_idx(a, p, (f + offset), t, n_chan);
                        int w = coeff_idx(a, p, b, f);

			if (a < N_ANT) {
				reduced_mul[a].x = input_data[i].x * coeff[2*w] + input_data[i].y * coeff[2*w + 1];
				reduced_mul[a].y = input_data[i].y * coeff[2*w] - input_data[i].x * coeff[2*w + 1];
			}
			else {
				reduced_mul[a].x = 0;
				reduced_mul[a].y = 0;
			}
			__syncthreads();

			for (int k = blockDim.x / 2; k > 0; k >>= 1) {
				if (a < k) {
					reduced_mul[a].x += reduced_mul[a + k].x;
					reduced_mul[a].y += reduced_mul[a + k].y;
				}
				__syncthreads();
			}
			if (a == 0) {
				int h1 = coh_bf_idx(p, b, (f + offset), t, n_chan);
				output_data[2 * h1] += reduced_mul[0].x;
				output_data[2 * h1 + 1] += reduced_mul[0].y;
			}
		
		}
	}
	return;
}


// Compute power of beamformer output (abs()^2)
__global__
void beamformer_power(float* bf_volt, float* bf_power, int offset, int n_chan, int nt) {
	int b = threadIdx.x; // Beam index
	int f = blockIdx.y;  // Frequency bin index
	int t = blockIdx.x;  // Time sample index

	int n_freq_streams = n_chan/N_STREAMS;

	if(f < n_freq_streams){	
		// Power = Absolute value squared of output -> r^2 + i^2
		int xp = coh_bf_idx(0, b, (f + offset), t, n_chan); // X polarization
		int yp = coh_bf_idx(1, b, (f + offset), t, n_chan); // Y polarization
		
		float x_pol_pow = (bf_volt[2 * xp] * bf_volt[2 * xp]) + (bf_volt[2 * xp + 1] * bf_volt[2 * xp + 1]); // XX*
		float y_pol_pow = (bf_volt[2 * yp] * bf_volt[2 * yp]) + (bf_volt[2 * yp + 1] * bf_volt[2 * yp + 1]); // YY*

		int h = pow_bf_nosti_idx(b, (f + offset), t, n_chan, nt);
		bf_power[h] = x_pol_pow + y_pol_pow; // XX* + YY*
	}
	
	return;
}

// Compute power of beamformer output (abs()^2)
__global__
void beamformer_power_sti(float* bf_volt, float* bf_power, int offset, int n_chan, int ns) {
	int f = blockIdx.x;  // Frequency bin index
	int b = blockIdx.y;  // Beam index
	int s = blockIdx.z;  // STI window index
	int t = threadIdx.x; // Time sample index

	int n_freq_streams = n_chan/N_STREAMS;
	
	int xp; // X polarization
	int yp; // Y polarization
	float x_pol_pow; // XX*
	float y_pol_pow; // YY*
	float beam_power;
	float scale = 1.0/N_TIME_STI;

	__shared__ float reduced_array[N_STI_BLOC];

	if(f < n_freq_streams){	
		if (t < N_TIME_STI) {
			// Power = Absolute value squared of output -> r^2 + i^2
			xp = coh_bf_idx(0, b, (f + offset), s*N_TIME_STI + t, n_chan); // X polarization
			yp = coh_bf_idx(1, b, (f + offset), s*N_TIME_STI + t, n_chan); // Y polarization

			x_pol_pow = (bf_volt[2 * xp] * bf_volt[2 * xp]) + (bf_volt[2 * xp + 1] * bf_volt[2 * xp + 1]); // XX*
			y_pol_pow = (bf_volt[2 * yp] * bf_volt[2 * yp]) + (bf_volt[2 * yp + 1] * bf_volt[2 * yp + 1]); // YY*
			beam_power = x_pol_pow + y_pol_pow; // XX* + YY*
			reduced_array[t] = beam_power;
		}
		else{
			reduced_array[t] = 0.0;
		}
		__syncthreads();

		// Reduction is performed by splitting up the threads in each block and summing them all up.
		// The number of threads in each block needs to be a power of two in order for the reduction to work. (No left over threads).
		for(int k = blockDim.x/2; k>0; k>>=1){
			if(t<k){
				reduced_array[t] += reduced_array[t+k];
			}
			__syncthreads();
		}
	
		// After reduction is complete, assign each reduced value to appropriate position in output array.
		if(t == 0){
			int h = pow_bf_idx(b, (f + offset), s, n_chan, ns);
			bf_power[h] = reduced_array[0]*scale;
		}
	}

	return;
}

// Run beamformer
float* run_beamformer(signed char* data_in, float* h_coefficient, int n_chan, int nt) {

	cudaError_t err_code;

	//const int freq_chans = N_FREQ_STREAM;
	int n_freq_streams = n_chan/N_STREAMS;
        int ns = nt/N_TIME_STI;

	// Transpose kernel: Specify grid and block dimensions
	dim3 dimBlock_transpose(N_ANT, N_POL, 1);
	dim3 dimGrid_transpose(nt, n_chan, 1);

	// Coherent beamformer kernel: Specify grid and block dimensions
	dim3 dimBlock_coh_bf(N_ANT, 1, 1);
	dim3 dimGrid_coh_bf(nt, n_chan, N_BEAM);

	// Output power of beamformer kernel: Specify grid and block dimensions
	//dim3 dimBlock_bf_pow(N_BEAM, 1, 1);
	//dim3 dimGrid_bf_pow(nt, n_chan, 1);

	// Output power of beamformer kernel with STI: Specify grid and block dimensions
	dim3 dimBlock_bf_pow(N_STI_BLOC, 1, 1);
	dim3 dimGrid_bf_pow(n_chan, N_BEAM, ns);

	float* d_data_bf = d_data_float;
	signed char* d_data_in = d_data_char;
	cuComplex* d_data_tra = d_data_comp;
	float* d_coefficient = d_coeff;
	float* d_bf_pow = d_coh_bf_pow;
	float* data_out = h_bf_pow;

	//printf("Before cudaMemcpy(HtoD) coefficients! \n");
	// Copy beamformer coefficients from host to device
	checkCuda(cudaMemcpy(d_coefficient, h_coefficient, N_COEFF/(MAX_COARSE_FREQ/n_chan) * sizeof(float), cudaMemcpyHostToDevice));

	// CUDA streams and events applied for optimization to possibly eliminate stalls.
	// cudaMemcpy(HtoD) and data_restructure kernel	
	const int streamSizeIn = (2*N_ANT*N_POL*nt*n_freq_streams);
	const unsigned long int streamBytesIn = (2*N_ANT*N_POL*nt*n_freq_streams*sizeof(signed char));

	// beamformer_power kernel and cudaMemcpy(DtoH)
	//const int streamSizePow = (N_BEAM*nt*n_freq_streams);
	//const unsigned long int streamBytesPow = (N_BEAM*nt*n_freq_streams*sizeof(float));
	const int streamSizePow = (N_BEAM*ns*n_freq_streams);
	const unsigned long int streamBytesPow = (N_BEAM*ns*n_freq_streams*sizeof(float));

	// Create events and streams
	// Events ////////////////////////////////////
	cudaEvent_t startEvent, stopEvent;
	checkCuda(cudaEventCreate(&startEvent));
	checkCuda(cudaEventCreate(&stopEvent));		
	checkCuda(cudaEventRecord(startEvent, 0));
	/////////////////////////////////////////////

	cudaStream_t stream[N_STREAMS];

	for (int i = 0; i < N_STREAMS; ++i) {
		checkCuda(cudaStreamCreate(&stream[i]));
	}

	for (int i = 0; i < N_STREAMS; ++i){

		//int offset = i * n_freq_streams;
		int offset_in = i * streamSizeIn;
		// Copy input data from host to device
		checkCuda(cudaMemcpyAsync(&d_data_in[offset_in], &data_in[offset_in], streamBytesIn, cudaMemcpyHostToDevice, stream[i]));
	}
	for (int i = 0; i < N_STREAMS; ++i){
		int offset = i * n_freq_streams;
		// Perform transpose on the data and convert to floats  
		data_transpose<<<dimGrid_transpose, dimBlock_transpose, 0, stream[i]>>>(d_data_in, d_data_tra, offset, n_chan, nt);
		err_code = cudaGetLastError();
		if (err_code != cudaSuccess) {
			printf("BF: data_transpose() kernel Failed: %s\n", cudaGetErrorString(err_code));
		}
	}
	for (int i = 0; i < N_STREAMS; ++i){
		int offset = i * n_freq_streams;
		// Perform beamforming operation
		// Use d_data_in for output since it is no longer being utilized,
		// and it is the same size as the output (4 GiB).
		coherent_beamformer<<<dimGrid_coh_bf, dimBlock_coh_bf, 0, stream[i]>>>(d_data_tra, d_coefficient, d_data_bf, offset, n_chan);
		err_code = cudaGetLastError();
		if (err_code != cudaSuccess) {
			printf("BF: coherent_beamformer() kernel Failed: %s\n", cudaGetErrorString(err_code));
		}
		//printf("Here coherent_beamformer, offset = %d \n", offset);	
	}
	for (int i = 0; i < N_STREAMS; ++i){
		int offset = i * n_freq_streams;
		// Compute power of beamformer output (abs()^2)
		//beamformer_power<<<dimGrid_bf_pow, dimBlock_bf_pow, 0, stream[i]>>>(d_data_bf, d_bf_pow, offset, n_chan);
		beamformer_power_sti<<<dimGrid_bf_pow, dimBlock_bf_pow, 0, stream[i]>>>(d_data_bf, d_bf_pow, offset, n_chan, ns);
		err_code = cudaGetLastError();
		if (err_code != cudaSuccess) {
			printf("BF: beamformer_power() kernel Failed: %s\n", cudaGetErrorString(err_code));
		}
	}
	for (int i = 0; i < N_STREAMS; ++i){
		int offset_pow = i * streamSizePow;
		// Copy output power from device to host
		checkCuda(cudaMemcpyAsync(&data_out[offset_pow], &d_bf_pow[offset_pow], streamBytesPow, cudaMemcpyDeviceToHost, stream[i]));
	}

	for (int i = 0; i < N_STREAMS; ++i) {
		checkCuda(cudaStreamSynchronize(stream[i]));
	}

	// Events ////////////////////////////////////
	checkCuda(cudaEventRecord(stopEvent, 0));
	checkCuda(cudaEventSynchronize(stopEvent));
	/////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////
	// Clean up streams
	checkCuda(cudaEventDestroy(startEvent));
	checkCuda(cudaEventDestroy(stopEvent));
	for (int i = 0; i < N_STREAMS; ++i) {
		checkCuda(cudaStreamDestroy(stream[i]));
	}

	return data_out;
}

// Generate simulated data
signed char* simulate_data(int n_chan, int nt) {
	signed char* data_sim;
	data_sim = (signed char*)calloc(N_INPUT, sizeof(signed char));

	/*
	'sim_flag' is a flag that indicates the kind of data that is simulated.
	sim_flag = 0 -> Ones
	sim_flag = 1 -> Repeating sequence of 1 to 64
	sim_flag = 2 -> Sequence of 1 to 64 placed in a particular bin (bin 6 for now)
	sim flag = 3 -> Simulated radio source in center beam assuming ULA
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
		for (int p = 0; p < N_POL; p++) {
			for (int t = 0; t < nt; t++) {
				for (int f = 0; f < n_chan; f++) {
					for (int a = 0; a < N_ANT; a++) {
						if (tmp >= N_ANT) {
							tmp = 0;
						}
						tmp = (tmp + 1) % (N_ANT+1);
						if(a < N_REAL_ANT){
							data_sim[2 * data_in_idx(p, t, f, a, nt, n_chan)] = tmp;
						}else{
							data_sim[2 * data_in_idx(p, t, f, a, nt, n_chan)] = 0;
						}
					}
				}
			}
		}
	}
	if (sim_flag == 2) {
		int tmp = 0;
		for (int p = 0; p < N_POL; p++) {
			for (int t = 0; t < nt; t++) {
				for (int a = 0; a < N_ANT; a++) {
					if (tmp >= N_ANT) {
						tmp = 0;
					}
					tmp = (tmp + 1) % (N_ANT+1);
					if(a < N_REAL_ANT){
						data_sim[2 * data_in_idx(p, t, 5, a, nt, n_chan)] = tmp;
						data_sim[2 * data_in_idx(p, t, 2, a, nt, n_chan)] = tmp;
					}else{
						data_sim[2 * data_in_idx(p, t, 5, a, nt, n_chan)] = 0;
						data_sim[2 * data_in_idx(p, t, 2, a, nt, n_chan)] = 0;
					}
				}
			}
		}
	}
	if (sim_flag == 3) {
		float c = 3e8; // Speed of light
		float c_freq = 1.25e9; // Center frequency
		float lambda = c / c_freq; // Wavelength
		float d = lambda / 2; // Distance between antennas
		float chan_band = 1; // Fine channel bandwidth in Hz
		double theta = 0; // SOI direction/angle of arrival
		double tau = 0; // Delay
		double rf_freqs = 0;
		double cb = 90; // Center beam in degrees

		float tmp_max = 1.0;
		float tmp_min = -1.0;

		for (int t = 0; t < nt; t++) {
			// Reduce the range of angles in order to prevent wrap around - That's what the 100 and 200 are for.
			theta = ((t/50 - (nt / 100)) + cb)*PI/180; // SOI direction/angle of arrival -> Moving across array over time i.e. angle changes each time sample
			tau = d * cos(theta) / c; // Delay
			for (int f = 0; f < n_chan; f++) {
				rf_freqs = chan_band * f + c_freq;
				for (int a = 0; a < N_ANT; a++) {
					if(a < N_REAL_ANT){
						// Requantize from doubles/floats to signed chars with a range from -128 to 127
						// X polarization
						data_sim[2 * data_in_idx(0, t, f, a, nt, n_chan)] = (signed char)((((cos(2 * PI * rf_freqs * a * tau) - tmp_min)/(tmp_max-tmp_min)) - 0.5)*256);
						data_sim[2 * data_in_idx(0, t, f, a, nt, n_chan) + 1] = (signed char)((((sin(2 * PI * rf_freqs * a * tau) - tmp_min)/(tmp_max-tmp_min)) - 0.5)*256);
						// Y polarization
						data_sim[2 * data_in_idx(1, t, f, a, nt, n_chan)] = (signed char)((((cos(2 * PI * rf_freqs * a * tau) - tmp_min)/(tmp_max-tmp_min)) - 0.5)*256);
						data_sim[2 * data_in_idx(1, t, f, a, nt, n_chan) + 1] = (signed char)((((sin(2 * PI * rf_freqs * a * tau) - tmp_min)/(tmp_max-tmp_min)) - 0.5)*256); // Make this negative if a different polarization is tested
					}else{
						// X polarization
						data_sim[2 * data_in_idx(0, t, f, a, nt, n_chan)] = 0;
						data_sim[2 * data_in_idx(0, t, f, a, nt, n_chan) + 1] = 0;
						// Y polarization
						data_sim[2 * data_in_idx(1, t, f, a, nt, n_chan)] = 0;
						data_sim[2 * data_in_idx(1, t, f, a, nt, n_chan) + 1] = 0; // Make this negative if a different polarization is tested
					}
				}
			}
		}
	}
	return data_sim;
}

// Generate simulated weights or coefficients
float* simulate_coefficients(int n_chan) {
	float* coeff_sim;
	coeff_sim = (float*)calloc(N_COEFF, sizeof(float));

	/*
	'sim_flag' is a flag that indicates the kind of data that is simulated.
	sim_flag = 0 -> Ones
	sim_flag = 1 -> Scale each beam by incrementing value i.e. beam 1 = 1, beam 2 = 2, ..., beam 64 = 64
	sim_flag = 2 -> Scale each beam by incrementing value in a particular bin (bin 3 and 6 for now). Match simulated data sim_flag = 2
	sim flag = 3 -> Simulated beams from 58 to 122 degrees. Assuming a ULA.
	*/
	int sim_flag = 3;
	if (sim_flag == 0) {
		for (int i = 0; i < ((N_COEFF/(MAX_COARSE_FREQ/n_chan)) / 2); i++) {
			coeff_sim[2 * i] = 1;
		}
	}
	if (sim_flag == 1) {
		int tmp = 0;
		
                for (int p = 0; p < N_POL; p++) {
                	for (int f = 0; f < n_chan; f++) {
				for (int a = 0; a < N_ANT; a++) {
					for (int b = 0; b < N_BEAM; b++) {
						if (tmp >= N_BEAM) {
							tmp = 0;
						}
						tmp = (tmp + 1) % (N_BEAM + 1);
                                		coeff_sim[2 * coeff_idx(a, p, b, f)] = tmp;
					}
				}
                	}
                }
		
	}
	if (sim_flag == 2) {
		int tmp = 0;
		for (int p = 0; p < N_POL; p++) {
                	for (int f = 0; f < n_chan; f++) {
				for (int a = 0; a < N_ANT; a++) {
					for (int b = 0; b < N_BEAM; b++) {
						if (tmp >= N_BEAM) {
							tmp = 0;
						}
						tmp = (tmp + 1) % (N_BEAM + 1);
                                		coeff_sim[2 * coeff_idx(a, p, b, f)] = tmp;
					}
				}
                	}
                }
	}
	if (sim_flag == 3) {
		float c = 3e8; // Speed of light
		float c_freq = 1.25e9; // Center frequency
		float lambda = c / c_freq; // Wavelength
		float d = lambda / 2; // Distance between antennas

		float chan_band = 1; // Fine channel bandwidth in Hz
		double rf_freqs = 0;

		float theta = 0; // Beam angle from 58 to 122 degrees
		float tau_beam = 0; // Delay


		for (int f = 0; f < n_chan; f++) {
			rf_freqs = chan_band * f + c_freq;
			for (int b = 0; b < N_BEAM; b++) {
				theta = ((b - (N_BEAM / 2)) + 90)*PI/180; // Beam angle from 58 to 122 degrees - Given SOI at 90 deg or moving across array, the beam with the most power is beamm 33
				tau_beam = d * cos(theta) / c; // Delay
				for (int a = 0; a < N_ANT; a++) {
					// X polarization
					coeff_sim[2 * coeff_idx(a, 0, b, f)] = cos(2 * PI * rf_freqs * a * tau_beam);
					coeff_sim[2 * coeff_idx(a, 0, b, f) + 1] = sin(2 * PI * rf_freqs * a * tau_beam);
					// Y polarization
					coeff_sim[2 * coeff_idx(a, 1, b, f)] = cos(2 * PI * rf_freqs * a * tau_beam);
					coeff_sim[2 * coeff_idx(a, 1, b, f) + 1] = sin(2 * PI * rf_freqs * a * tau_beam);
				}
			}
		}
	}

	return coeff_sim;
}

// Generate weights or coefficients with calculated delays (with delay polynomials (tau), coarse frequency channel (coarse_chan), and epoch (t))
//float* generate_coefficients(float* tau, float* telstate_phase, double* coarse_chan, int n_chan, uint64_t n_real_ant, float t) {
float* generate_coefficients(float* tau, double* coarse_chan, int n_chan, uint64_t n_real_ant, float t) {
	float* coefficients;
	coefficients = (float*)calloc(N_COEFF, sizeof(float));
	float delay_rate = 0;
	float delay_offset = 0;

	for (int p = 0; p < N_POL; p++) {
		for (int f = 0; f < n_chan; f++) {
			for (int b = 0; b < N_BEAM; b++) {
				for (int a = 0; a < N_ANT; a++) {
					if(a < n_real_ant){
						delay_rate = tau[delay_idx(1,a,b,n_real_ant)];
						delay_offset = tau[delay_idx(0,a,b,n_real_ant)];
						coefficients[2 * coeff_idx(a, p, b, f)] = cos(2 * PI * coarse_chan[f] * (t*delay_rate + delay_offset));
						coefficients[2 * coeff_idx(a, p, b, f) + 1] = sin(2 * PI * coarse_chan[f] * (t*delay_rate + delay_offset));
						//coefficients[2 * coeff_idx(a, p, b, f)] = telstate_phase[2*phase_idx(a, p, f)]*cos(2 * PI * coarse_chan[f] * (t*delay_rate + delay_offset)) - telstate_phase[2*phase_idx(a, p, f) + 1]*sin(2 * PI * coarse_chan[f] * (t*delay_rate + delay_offset));
						//coefficients[2 * coeff_idx(a, p, b, f) + 1] = telstate_phase[2*phase_idx(a, p, f)]*sin(2 * PI * coarse_chan[f] * (t*delay_rate + delay_offset)) + telstate_phase[2*phase_idx(a, p, f) + 1]*cos(2 * PI * coarse_chan[f] * (t*delay_rate + delay_offset));
					}else{
						coefficients[2 * coeff_idx(a, p, b, f)] = 0;
						coefficients[2 * coeff_idx(a, p, b, f) + 1] = 0;
					}
				}
			}
		}
	}

	return coefficients;
}

// The input_data_pin() function uses cudaHostRegister() to allocate the input host
// array in pinned memory.
// This speeds up the cudaMemcpy() and enables implementation into HASHPIPE/RTOS.
//void input_data_pin(float * data_in_pin) {
void input_data_pin(signed char * data_in_pin) {
	checkCuda(cudaHostRegister(data_in_pin, N_INPUT*sizeof(signed char), cudaHostRegisterPortable));
}

// The data_coeff_pin() function uses cudaHostRegister() to allocate the input host
// array in pinned memory.
// This speeds up the cudaMemcpy() and enables implementation into HASHPIPE/RTOS.
void coeff_pin(float * data_coeff_pin) {
	checkCuda(cudaHostRegister(data_coeff_pin, N_COEFF*sizeof(float), cudaHostRegisterPortable));
}

// The output_data_pin() function uses cudaHostRegister() to allocate the output host
// array in pinned memory.
// This speeds up the cudaMemcpy() and enables implementation into HASHPIPE/RTOS.
void output_data_pin(float * data_out_pin) {
	checkCuda(cudaHostRegister(data_out_pin, N_BF_POW*sizeof(float), cudaHostRegisterPortable));
}

// Unregister host arrays from pinned memory
void unregister_data(void * data_unregister){
	checkCuda(cudaHostUnregister(data_unregister));
}

// Free memory
void cohbfCleanup() {
	// Free up GPU memory at the end of a program
	if (d_data_char != NULL) {
		cudaFree(d_data_char);
	}
	if (d_data_float != NULL) {
		cudaFree(d_data_float);
	}
	if (d_data_comp != NULL) {
		cudaFree(d_data_comp);
	}
	if (d_coeff != NULL) {
		cudaFree(d_coeff);
	}
	if (d_coh_bf_pow != NULL) {
		cudaFree(d_coh_bf_pow);
	}
}

//Comment out main() function when compiling for hpguppi
/* // <----Uncomment here if testing standalone code
// Test all of the kernels and functions, and write the output to
// a text file for analysis
int main() {
	printf("Here!\n");
	// 1k mode
	//int n_chan = 16; 
        //int nt = 32768;
	// 4k mode
    	int n_chan = 64;
        int nt = 8192;
	// 32k mode
    	//int n_chan = 512;
        //int nt = 1024;

	// Allocate memory to all arrays used by run_beamformer() 
	init_beamformer();

        printf("After init_beamformer() \n");

	// Generate simulated data
	signed char* sim_data = simulate_data(n_chan, nt);

        printf("After simulate_data() \n");

	// Register the array in pinned memory to speed HtoD mem copy
	input_data_pin(sim_data);

	printf("After input_data_pin(sim_data); \n");

	// Generate simulated weights or coefficients
	float* sim_coefficients = simulate_coefficients(n_chan);

	printf("After simulate_coefficients(); \n");

	//printf("Here3!\n");
	// Register the array in pinned memory to speed HtoD mem copy
	coeff_pin(sim_coefficients);

	printf("After coeff_pin(sim_coefficients); \n");

	//printf("real sim_data: %d and imag sim_data: %d\n", sim_data[10485768], sim_data[10485769]);
	//printf("real sim_coef: %f and imag sim_coef: %f\n", sim_coefficients[104], sim_coefficients[105]);

	//printf("real sim_data2: %d and imag sim_data2: %d\n", sim_data[8388616], sim_data[8388617]);
	//printf("real sim_coef2: %f and imag sim_coef2: %f\n", sim_coefficients[106], sim_coefficients[107]);
	// Allocate memory for output array
	float* output_data;
	//output_data = (float*)calloc(N_BF_POW, sizeof(float));
	//output_data = (float*)calloc(N_OUTPUT, sizeof(float));
	//output_data_pin(output_data);

	printf("Here5!\n");

	float time_taken = 0;
	float bf_time = 0;
	int num_runs = 10;

	// Start timing beamformer computation //
	struct timespec tval_before, tval_after;

	for(int ii = 0; ii < num_runs; ii++){
		// Start timing beamformer computation //
		clock_gettime(CLOCK_MONOTONIC, &tval_before);

		// Run beamformer 
		output_data = run_beamformer(sim_data, sim_coefficients, n_chan, nt);
		//run_beamformer(h_data, h_coeff, output_data);

		// Stop timing beamforming computation //
		clock_gettime(CLOCK_MONOTONIC, &tval_after);
		time_taken = (float)(tval_after.tv_sec - tval_before.tv_sec); //*1e6; // Time in seconds since epoch
		time_taken = time_taken + (float)(tval_after.tv_nsec - tval_before.tv_nsec)*1e-9; // Time in nanoseconds since 'tv_sec - start and end'
		bf_time += time_taken;
		//printf("Time taken: %f s\n", time_taken);
	}
	printf("Average beamformer processing time: %f s\n", bf_time/num_runs);

	printf("Here6, Beamformer output: %f \n", output_data[0]);
	
	// Write data to text file for analysis
	char output_filename[128];

	printf("Here7!\n");

	//strcpy(output_filename, "C:\Users\ruzie\OneDrive\Desktop\Work\CUDA_code\output_d.txt");
	strcpy(output_filename, "output_d_cuda.txt");

	printf("Here8!\n");

	FILE* output_file;

	printf("Here9!\n");

	output_file = fopen(output_filename, "w");

	printf("Here10!\n");

	for (int ii = 0; ii < ((N_BF_POW*n_chan*nt)/(MAX_COARSE_FREQ*N_TIME)); ii++) { // Write up to the size of the data corresponding to 1k, 4k or 32k mode
		//fprintf(output_file, "%c\n", output_data[ii]);
		fprintf(output_file, "%g\n", output_data[ii]);
	}

	printf("Here11!\n");

	fclose(output_file);

	printf("Closed output file.\n");

	//unregister_data((float *)sim_data);
	//free(sim_data);
	printf("After unregister.\n");	
	//free(sim_coefficients);
	printf("After freeing coefficients.\n");
	//free(output_data);	

	printf("Freed output array and unregistered arrays in pinned memory.\n");

	// Free up device memory
	//cudaFreeHost(h_data);
	//cudaFreeHost(h_coeff);
	cohbfCleanup();

	printf("Here11!\n");

	return 0;
}
*/ // <----Uncomment here if testing standalone code
