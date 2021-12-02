#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <curand.h>
#include <assert.h>
//#include <unistd.h>
#include <cublas_v2.h>
#include <time.h>
#include <sys/time.h>
#include <iostream>
#include <string.h>
//#include <complex.h>
#include <math.h>
#include <cuComplex.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <device_launch_parameters.h>
#include "coherent_beamformer_char_in.h"

using namespace std;

// Generate simulated data
//signed char* simulate_data();

// Generate simulated weights or coefficients
//float* simulate_coefficients();

// Perform transpose on the data and convert to floats
__global__
void data_transpose(signed char* data_in, cuComplex* data_tra, int offset);

// Convert weights from float to cuComplex
//__global__
//void beamformer_coefficient(float* coeff_float, cuComplex* coeff_complex);

// Perform beamforming operation
__global__
void coherent_beamformer(cuComplex* input_data, float* coeff, float* output_data, int offset);

// Compute power of beamformer output (abs()^2)
__global__
void beamformer_power(float* bf_volt, float* bf_power, int offset);
//void beamformer_power(float* bf_volt, signed char* bf_power, int offset);


// Convenience function for checking CUDA runtime API results
// can be wrapped around any runtime API call. No-op in release builds.
inline
cudaError_t checkCuda(cudaError_t result)
{
//#if defined(DEBUG) || defined(_DEBUG)
  if (result != cudaSuccess) {
    fprintf(stderr, "CUDA Runtime Error: %s\n", 
            cudaGetErrorString(result));
    assert(result == cudaSuccess);
  }
//#endif
  return result;
}

//float* h_data = NULL;
//float* h_coeff = NULL;
float* d_data_float = NULL;
signed char* d_data_char = NULL;
cuComplex* d_data_comp = NULL;
float* d_coeff = NULL;
//cuComplex* d_coeff_comp = NULL;
//float* d_coh_bf_out = NULL;
float* d_coh_bf_pow = NULL;
float* h_bf_pow = NULL;
//signed char* d_coh_bf_pow = NULL;
//signed char* h_bf_pow = NULL;
// Allocate memory to all arrays 
void init_beamformer() {
	printf("Here In init_beamformer()! \n");
	// Allocate pinned memory for input data
	//checkCuda(cudaMallocHost((void **)&h_data, N_INPUT * sizeof(float)));
	//printf("Here 1st cudaMallocHost! \n");
	
	// Allocate pinnted memery for beamformer coefficients
	//checkCuda(cudaMallocHost((void **)&h_coeff, N_COEFF * sizeof(float)));
	//printf("Here 2nd cudaMallocHost! \n");
	//cudaHostAlloc((void **)&h_data, N_INPUT * sizeof(float));

	// Allocate memory for input data float type
	checkCuda(cudaMalloc((void **)&d_data_float, (N_INPUT) * sizeof(float)));
	checkCuda(cudaMalloc((void **)&d_data_char, (N_INPUT) * sizeof(signed char)));
	printf("Here 1st cudaMalloc! \n");

	// Allocate memory for input data cuComplex type
	checkCuda(cudaMalloc((void **)&d_data_comp, (N_INPUT) * sizeof(cuComplex) / 2));
	printf("Here 2nd cudaMalloc! \n");

	/*
	size_t f, t;
    	cudaSetDevice(0);
    	cudaMemGetInfo(&f, &t);
    	fprintf(stdout,"Free: %zu bytes, Available: %zu bytes \n",f,t);
	*/

	// Allocate memory for coefficients float type
	checkCuda(cudaMalloc((void **)&d_coeff, N_COEFF * sizeof(float)));
	printf("Here 3rd cudaMalloc! \n");

	// Allocate memory for coefficients cuComplex type
	//checkCuda(cudaMalloc((void **)&d_coeff_comp, N_COEFF * sizeof(cuComplex) / 2));
	//printf("Here 4th cudaMalloc! \n");

	// Allocate memory for coherent beamformer output
	//checkCuda(cudaMalloc((void **)&d_coh_bf_out, N_OUTPUT * sizeof(float)));
	//printf("Here 5th cudaMalloc! \n");

	// Allocate memory for output power of coherent beamformer
        checkCuda(cudaMalloc((void **)&d_coh_bf_pow, (N_BF_POW) * sizeof(float)));
	printf("Here 4th cudaMalloc! \n");

	checkCuda(cudaMallocHost((void **)&h_bf_pow, (N_BF_POW) * sizeof(float)));

	//checkCuda(cudaMalloc((void **)&d_coh_bf_pow, (N_BF_POW) * sizeof(signed char)));
	//printf("Here 4th cudaMalloc! \n");

	//checkCuda(cudaMallocHost((void **)&h_bf_pow, (N_BF_POW) * sizeof(signed char)));

	return;
}

// Set arrays to zero after a block is processed
void set_to_zero(){
	checkCuda(cudaMemset(d_data_float, 0, (N_INPUT) * sizeof(float)));
}

// Perform transpose on the data and convert to floats
__global__
void data_transpose(signed char* data_in, cuComplex* data_tra, int offset) {
	int a = threadIdx.x; // Antenna index
	int p = threadIdx.y; // Polarization index
	int f = blockIdx.y;  // Frequency index
	int t = blockIdx.x;  // Time sample index

	// If the input data is not float e.g. signed char, just multiply it by '1.0f' to convert it to a float
	if(f < N_FREQ_STREAM){
		int h_in = data_in_idx(a, p, (f + offset), t);
		int h_tr = data_tr_idx(a, p, (f + offset), t);
		//if(a < N_REAL_ANT){
		data_tra[h_tr].x = data_in[2*h_in]*1.0f;
		data_tra[h_tr].y = data_in[2*h_in + 1]*1.0f;
		//}else{
		//	data_tra[h_tr].x = 0*1.0f;
		//	data_tra[h_tr].y = 0*1.0f;
		//}
	}

	return;
}

/*
// Convert weights from float to cuComplex
__global__
void beamformer_coefficient(float* coeff_float, cuComplex* coeff_complex) {
	// Product of antenna and beam dimensions exceeds 1024 so beams
	// are blocks rather than threads to allow for increase in numbers  
	int a = threadIdx.x; // Antenna index
	int p = threadIdx.y; // Polarization index
	int b = blockIdx.y;  // Beam index
	int f = blockIdx.x;  // Frequency bin index 
	coeff_complex[coeff_idx(a, p, b, f)].x = coeff_float[2*coeff_idx(a, p, b, f)];
	coeff_complex[coeff_idx(a, p, b, f)].y = coeff_float[2*coeff_idx(a, p, b, f) + 1];
	
	return;
}
*/

// Perform beamforming operation
__global__
void coherent_beamformer(cuComplex* input_data, float* coeff, float* output_data, int offset) {
	/*
	int p = threadIdx.x; // Polarization index
	int f = blockIdx.x;  // Frequency index
	int t = blockIdx.y;  // Time sample index
	int b = blockIdx.z;  // Beam index
	cuComplex bf_product;
	cuComplex bf_in_data;
	cuComplex bf_coeff;
	for (int a = 0; a < N_ANT; a++) { // Antenna index
		bf_in_data.x = input_data[data_tr_idx(a, p, f, t)].x;
		bf_in_data.y = input_data[data_tr_idx(a, p, f, t)].y;
		bf_coeff.x = coeff[coeff_idx(a, p, b, f)].x;
		bf_coeff.y = coeff[coeff_idx(a, p, b, f)].y;
		// Complex multiplication of data and coefficients
		bf_product.x = (bf_in_data.x * bf_coeff.x) - (bf_in_data.y * bf_coeff.y);
		bf_product.y = (bf_in_data.x * bf_coeff.y) + (bf_in_data.y * bf_coeff.x);
		// Beamform (Sum all antennas)
		output_data[2*coh_bf_idx(p, b, f, t)] += bf_product.x;
		output_data[2*coh_bf_idx(p, b, f, t) + 1] += bf_product.y;
	}
	*/
	int a = threadIdx.x; // Antenna index
	int f = blockIdx.y;  // Frequency index
	int t = blockIdx.x;  // Time sample index
	int b = blockIdx.z;  // Beam index

	__shared__ cuFloatComplex reduced_mul[N_ANT];

	if(f < N_FREQ_STREAM){
		for (int p = 0; p < N_POL; p++) { // Polarization index
			// Reinitialize output_data since we are using the input data array to be more efficient
			//int h = coh_bf_idx(p, b, (f + offset), t);
			//output_data[2 * h] = 0;
			//output_data[2 * h + 1] = 0;

	
			int i = data_tr_idx(a, p, (f + offset), t);
			int w = coeff_idx(a, b);

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
				int h1 = coh_bf_idx(p, b, (f + offset), t);
				output_data[2 * h1] += reduced_mul[0].x;
				output_data[2 * h1 + 1] += reduced_mul[0].y;
			}
		
		}
	}
	return;
}


// Compute power of beamformer output (abs()^2)
__global__
void beamformer_power(float* bf_volt, float* bf_power, int offset) {
	int b = threadIdx.x; // Beam index
	int f = blockIdx.y;  // Frequency bin index
	int t = blockIdx.x;  // Time sample index

	if(f < N_FREQ_STREAM){	
		// Power = Absolute value squared of output -> r^2 + i^2
		int xp = coh_bf_idx(0, b, (f + offset), t); // X polarization
		int yp = coh_bf_idx(1, b, (f + offset), t); // Y polarization
		
		float x_pol_pow = (bf_volt[2 * xp] * bf_volt[2 * xp]) + (bf_volt[2 * xp + 1] * bf_volt[2 * xp + 1]); // XX*
		float y_pol_pow = (bf_volt[2 * yp] * bf_volt[2 * yp]) + (bf_volt[2 * yp + 1] * bf_volt[2 * yp + 1]); // YY*

		int h = pow_bf_idx(b, (f + offset), t);
		bf_power[h] = x_pol_pow + y_pol_pow; // XX* + YY*
	}
	
	//bf_power[pow_bf_idx(0, b, f, t)] = (bf_volt[2*xp]*bf_volt[2*xp]) + (bf_volt[2*xp + 1]*bf_volt[2*xp + 1]); // XX*
	//bf_power[pow_bf_idx(1, b, f, t)] = (bf_volt[2*yp]*bf_volt[2*yp]) + (bf_volt[2*yp + 1]*bf_volt[2*yp + 1]); // YY*
	//bf_power[pow_bf_idx(2, b, f, t)] = (bf_volt[2*xp]*bf_volt[2*yp]) + (bf_volt[2*xp + 1]*bf_volt[2*yp + 1]); // XY* real
	//bf_power[pow_bf_idx(3, b, f, t)] = (bf_volt[2*xp + 1]*bf_volt[2*yp]) - (bf_volt[2*xp]*bf_volt[2*yp + 1]); // XY* imag
	
	return;
}

/*
// Compute power of beamformer output (abs()^2)
__global__
void beamformer_power(float* bf_volt, signed char* bf_power, int offset) {
//void beamformer_power(float* bf_volt, float* bf_power, int offset) {
	int b = threadIdx.x; // Beam index
	int f = blockIdx.x;  // Frequency bin index
	int t = blockIdx.y;  // Time sample index

	if(f < N_FREQ_STREAM){	
		// Power = Absolute value squared of output -> r^2 + i^2
		int xp = coh_bf_idx(0, b, (f + offset), t); // X polarization
		int yp = coh_bf_idx(1, b, (f + offset), t); // Y polarization
		
		float x_pol_pow = (bf_volt[2 * xp] * bf_volt[2 * xp]) + (bf_volt[2 * xp + 1] * bf_volt[2 * xp + 1]); // XX*
		float y_pol_pow = (bf_volt[2 * yp] * bf_volt[2 * yp]) + (bf_volt[2 * yp + 1] * bf_volt[2 * yp + 1]); // YY*

		float tmp_pow = x_pol_pow + y_pol_pow; // XX* + YY*
		// int h = pow_bf_idx(b, (f + offset), t);
		// bf_power[h] = x_pol_pow + y_pol_pow; // XX* + YY*

		int xi = coh_bf_idx(0, 0, 0, 0);
		int yi = coh_bf_idx(1, 0, 0, 0);
		float tmp_x_pol_pow = (bf_volt[2 * xi] * bf_volt[2 * xi]) + (bf_volt[2 * xi + 1] * bf_volt[2 * xi + 1]); // XX*
		float tmp_y_pol_pow = (bf_volt[2 * yi] * bf_volt[2 * yi]) + (bf_volt[2 * yi + 1] * bf_volt[2 * yi + 1]); // YY*
		float tmp_pow_init = tmp_x_pol_pow + tmp_y_pol_pow;


		// Find minimum value of total power
		float min_pow = tmp_pow_init;
		min_pow = min(min_pow, tmp_pow);

		// Find maximum value of total power
		float max_pow = tmp_pow_init;
		max_pow = max(max_pow, tmp_pow);

		// Requantize output from 32 bit floats to 8 bit signed chars
		int h = pow_bf_idx(b, (f + offset), t);
		bf_power[h] = (signed char)((((tmp_pow - min_pow)/(max_pow - min_pow)) - 0.5)*256); // XX* + YY*
	}
	
	
	//bf_power[pow_bf_idx(0, b, f, t)] = (bf_volt[2*xp]*bf_volt[2*xp]) + (bf_volt[2*xp + 1]*bf_volt[2*xp + 1]); // XX*
	//bf_power[pow_bf_idx(1, b, f, t)] = (bf_volt[2*yp]*bf_volt[2*yp]) + (bf_volt[2*yp + 1]*bf_volt[2*yp + 1]); // YY*
	//bf_power[pow_bf_idx(2, b, f, t)] = (bf_volt[2*xp]*bf_volt[2*yp]) + (bf_volt[2*xp + 1]*bf_volt[2*yp + 1]); // XY* real
	//bf_power[pow_bf_idx(3, b, f, t)] = (bf_volt[2*xp + 1]*bf_volt[2*yp]) - (bf_volt[2*xp]*bf_volt[2*yp + 1]); // XY* imag
	

	return;
}
*/
// Run beamformer
//void run_beamformer(float* data_in, float* h_coefficient, float* data_out) {
float* run_beamformer(signed char* data_in, float* h_coefficient) {
//signed char* run_beamformer(signed char* data_in, float* h_coefficient) {
	/*
	// Allocate input data in pinned memory
	// (This may take longer than it's worth to implement pinned memory)
	*h_data = *data_in;
	*/

	cudaError_t err_code;

	//const int freq_chans = N_FREQ_STREAM;
	const int nStreams = N_STREAMS; // Number of streams that make up all of the data. Split in frequency blocks (largest dimension)
	//printf("Total frequency channels: %d , num streams: %d \n", freq_chans, nStreams);
	//printf("Total frequency channels: %d , num streams: %d \n", N_FREQ_STREAM, nStreams);

	// Transpose kernel: Specify grid and block dimensions
	dim3 dimBlock_transpose(N_ANT, N_POL, 1);
	dim3 dimGrid_transpose(N_TIME, N_FREQ, 1);
	//dim3 dimGrid_transpose(N_FREQ, N_TIME, 1);

	// Beamformer coefficient kernel (float to complex): Specify grid and block dimensions
	//dim3 dimBlock_bf_coeff(N_ANT, N_POL, 1);
	//dim3 dimGrid_bf_coeff(N_FREQ, N_BEAM, 1);

	// Coherent beamformer kernel: Specify grid and block dimensions
	//dim3 dimBlock_coh_bf(N_ANT, N_POL, 1);
	dim3 dimBlock_coh_bf(N_ANT, 1, 1);
	dim3 dimGrid_coh_bf(N_TIME, N_FREQ, N_BEAM);

	// Output power of beamformer kernel: Specify grid and block dimensions
	dim3 dimBlock_bf_pow(N_BEAM, 1, 1);
	dim3 dimGrid_bf_pow(N_TIME, N_FREQ, 1);

	float* d_data_bf = d_data_float;
	signed char* d_data_in = d_data_char;
	cuComplex* d_data_tra = d_data_comp;
	float* d_coefficient = d_coeff;
	//cuComplex* d_coeff_c = d_coeff_comp;
	//float* d_bf_output = d_coh_bf_out;
	float* d_bf_pow = d_coh_bf_pow;
	float* data_out = h_bf_pow;

	//signed char* d_bf_pow = d_coh_bf_pow;
	//signed char* data_out = h_bf_pow;

	//printf("Before cudaMemcpy(HtoD) coefficients! \n");
	// Copy beamformer coefficients from host to device
	checkCuda(cudaMemcpy(d_coefficient, h_coefficient, N_COEFF * sizeof(float), cudaMemcpyHostToDevice));
	//printf("Here cudaMemcpy(HtoD) coefficients! \n");

	// CUDA streams and events applied for optimization to possibly eliminate stalls.
	// cudaMemcpy(HtoD) and data_restructure kernel	
	const int streamSizeIn = (2*N_ANT*N_POL*N_TIME*N_FREQ_STREAM);
	//const int streamSizeIn = freq_chans;
	//const unsigned long int streamBytesIn = (2*N_ANT*N_POL*N_TIME*N_FREQ_STREAM*sizeof(float));
	const unsigned long int streamBytesIn = (2*N_ANT*N_POL*N_TIME*N_FREQ_STREAM*sizeof(signed char));
	//printf("Input size: %d in bytes: %lu \n", streamSizeIn, streamBytesIn);

	// cudaMemcpy(HtoD) coefficients
	//const int streamSizeCo = N_ANT*N_BEAM*N_TIME/nStreams;
	//const int streamBytesCo = streamSizeCo * sizeof(float);

	// coherent_beamformer kernel  
	//const unsigned int streamSizeBF = (unsigned int)(2*N_BEAM*N_POL*N_TIME*N_FREQ_STREAM);
	//const unsigned int streamBytesBF = streamSizeBF * sizeof(float);
	//printf("BF output size: %u \n", streamSizeBF);

	// beamformer_power kernel and cudaMemcpy(DtoH)
	const int streamSizePow = (N_BEAM*N_TIME*N_FREQ_STREAM);
	//const int streamSizePow = freq_chans;
	const unsigned long int streamBytesPow = (N_BEAM*N_TIME*N_FREQ_STREAM*sizeof(float));
	//const unsigned long int streamBytesPow = (N_BEAM*N_TIME*N_FREQ_STREAM*sizeof(signed char));	
	//printf("BF power output size: %d in bytes: %lu \n", streamSizePow, streamBytesPow);

	// Create events and streams
	// Events ////////////////////////////////////
	cudaEvent_t startEvent, stopEvent;
	checkCuda(cudaEventCreate(&startEvent));
	checkCuda(cudaEventCreate(&stopEvent));		
	checkCuda(cudaEventRecord(startEvent, 0));
	/////////////////////////////////////////////

	cudaStream_t stream[nStreams];

	for (int i = 0; i < nStreams; ++i) {
		checkCuda(cudaStreamCreate(&stream[i]));
	}

	for (int i = 0; i < nStreams; ++i){

		//int offset = i * N_FREQ_STREAM;
		int offset_in = i * streamSizeIn;
		// Copy input data from host to device
		checkCuda(cudaMemcpyAsync(&d_data_in[offset_in], &data_in[offset_in], streamBytesIn, cudaMemcpyHostToDevice, stream[i]));
		//checkCuda(cudaMemcpyAsync(d_data_in, &data_in[offset_in], streamBytesIn, cudaMemcpyHostToDevice, stream[i]));
		//printf("First cudaMemcpyAsync(HtoD) in run_beamformer(), offset_in = %d and offset = %d \n", offset_in, offset);
	}
	for (int i = 0; i < nStreams; ++i){
		int offset = i * N_FREQ_STREAM;
		// Perform transpose on the data and convert to floats  
		data_transpose<<<dimGrid_transpose, dimBlock_transpose, 0, stream[i]>>>(d_data_in, d_data_tra, offset);
		err_code = cudaGetLastError();
		if (err_code != cudaSuccess) {
			printf("BF: data_transpose() kernel Failed: %s\n", cudaGetErrorString(err_code));
		}
		//printf("Here data_transpose! \n");
	}
	for (int i = 0; i < nStreams; ++i){
		int offset = i * N_FREQ_STREAM;
		// Perform beamforming operation
		// Use d_data_in for output since it is no longer being utilized,
		// and it is the same size as the output (4 GiB).
		//unsigned long int offset_bf = i * streamSizeBF;
		coherent_beamformer<<<dimGrid_coh_bf, dimBlock_coh_bf, 0, stream[i]>>>(d_data_tra, d_coefficient, d_data_bf, offset);
		err_code = cudaGetLastError();
		if (err_code != cudaSuccess) {
			printf("BF: coherent_beamformer() kernel Failed: %s\n", cudaGetErrorString(err_code));
		}
		//printf("Here coherent_beamformer, offset = %d \n", offset);	
	}
	for (int i = 0; i < nStreams; ++i){
		int offset = i * N_FREQ_STREAM;
		// Compute power of beamformer output (abs()^2)
		//int offset_pow = i * streamSizePow;
		beamformer_power<<<dimGrid_bf_pow, dimBlock_bf_pow, 0, stream[i]>>>(d_data_bf, d_bf_pow, offset);
		err_code = cudaGetLastError();
		if (err_code != cudaSuccess) {
			printf("BF: beamformer_power() kernel Failed: %s\n", cudaGetErrorString(err_code));
		}
		//printf("Here beamformer_power, offset_pow = %d and offset = %d \n", offset_pow, offset);
	}
	for (int i = 0; i < nStreams; ++i){
		int offset_pow = i * streamSizePow;
		// Copy output power from device to host
		checkCuda(cudaMemcpyAsync(&data_out[offset_pow], &d_bf_pow[offset_pow], streamBytesPow, cudaMemcpyDeviceToHost, stream[i]));
		//checkCuda(cudaMemcpyAsync(&data_out[offset_pow], d_bf_pow, streamBytesPow, cudaMemcpyDeviceToHost, stream[i]));
		//printf("Here cudaMemcpyAsync(DtoH)! \n");

	}

	for (int i = 0; i < nStreams; ++i) {
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
	for (int i = 0; i < nStreams; ++i) {
		checkCuda(cudaStreamDestroy(stream[i]));
	}

	return data_out;
}

// Generate simulated data
signed char* simulate_data() {
	signed char* data_sim;
	data_sim = (signed char*)calloc(N_INPUT, sizeof(signed char));
	//checkCuda(cudaMallocHost((void **)&data_sim, N_INPUT * sizeof(float)));

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
			for (int t = 0; t < N_TIME; t++) {
				for (int f = 0; f < N_FREQ; f++) {
					for (int a = 0; a < N_ANT; a++) {
						if (tmp >= N_ANT) {
							tmp = 0;
						}
						tmp = (tmp + 1) % (N_ANT+1);
						if(a < N_REAL_ANT){
							data_sim[2 * data_in_idx(a, p, f, t)] = tmp;
						}else{
							data_sim[2 * data_in_idx(a, p, f, t)] = 0;
						}
					}
				}
			}
		}
	}
	if (sim_flag == 2) {
		int tmp = 0;
		for (int p = 0; p < N_POL; p++) {
			for (int t = 0; t < N_TIME; t++) {
				for (int a = 0; a < N_ANT; a++) {
					if (tmp >= N_ANT) {
						tmp = 0;
					}
					tmp = (tmp + 1) % (N_ANT+1);
					if(a < N_REAL_ANT){
						data_sim[2 * data_in_idx(a, p, 5, t)] = tmp;
						data_sim[2 * data_in_idx(a, p, 2, t)] = tmp;
					}else{
						data_sim[2 * data_in_idx(a, p, 5, t)] = 0;
						data_sim[2 * data_in_idx(a, p, 2, t)] = 0;
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

		//float* rf_freqs = (float*)calloc(N_FREQ, sizeof(float));
		//for (int i = 0; i < N_FREQ; i++) {
		//	rf_freqs[i] = chan_band * i + c_freq;
		//}

		//float* theta = (float*)calloc(N_TIME, sizeof(float)); // SOI direction/angle of arrival
		//float* tau = (float*)calloc(N_TIME, sizeof(float)); // Delay

		double theta = 0; // SOI direction/angle of arrival
		double tau = 0; // Delay
		double rf_freqs = 0;
		double cb = 90; // Center beam in degrees

		float tmp_max = 1.0;
		float tmp_min = -1.0;

		for (int t = 0; t < N_TIME; t++) {
			// Reduce the range of angles in order to prevent wrap around - That's what the 100 and 200 are for.
			theta = ((t/50 - (N_TIME / 100)) + cb)*PI/180; // SOI direction/angle of arrival -> Moving across array over time i.e. angle changes each time sample
			tau = d * cos(theta) / c; // Delay
			for (int f = 0; f < N_FREQ; f++) {
				rf_freqs = chan_band * f + c_freq;
				for (int a = 0; a < N_ANT; a++) {
					if(a < N_REAL_ANT){
						// Requantize from doubles/floats to signed chars with a range from -128 to 127
						// X polarization
						data_sim[2 * data_in_idx(a, 0, f, t)] = (signed char)((((cos(2 * PI * rf_freqs * a * tau) - tmp_min)/(tmp_max-tmp_min)) - 0.5)*256);
						data_sim[2 * data_in_idx(a, 0, f, t) + 1] = (signed char)((((sin(2 * PI * rf_freqs * a * tau) - tmp_min)/(tmp_max-tmp_min)) - 0.5)*256);
						// Y polarization
						data_sim[2 * data_in_idx(a, 1, f, t)] = (signed char)((((cos(2 * PI * rf_freqs * a * tau) - tmp_min)/(tmp_max-tmp_min)) - 0.5)*256);
						data_sim[2 * data_in_idx(a, 1, f, t) + 1] = (signed char)((((sin(2 * PI * rf_freqs * a * tau) - tmp_min)/(tmp_max-tmp_min)) - 0.5)*256); // Make this negative if a different polarization is tested
					}else{
						// X polarization
						data_sim[2 * data_in_idx(a, 0, f, t)] = 0;
						data_sim[2 * data_in_idx(a, 0, f, t) + 1] = 0;
						// Y polarization
						data_sim[2 * data_in_idx(a, 1, f, t)] = 0;
						data_sim[2 * data_in_idx(a, 1, f, t) + 1] = 0; // Make this negative if a different polarization is tested
					}
				}
			}
		}
	}
	//for (int a = 0; a < N_ANT; a++){
	//	printf("Antenna %d = %f\n", (a+1), data_sim[2 * data_in_idx(a, 0, 0, 0)]);
	//}
	return data_sim;
}

// Generate simulated weights or coefficients
float* simulate_coefficients() {
	float* coeff_sim;
	coeff_sim = (float*)calloc(N_COEFF, sizeof(float));
	//checkCuda(cudaMallocHost((void **)&coeff_sim, N_COEFF * sizeof(float)));

	/*
	'sim_flag' is a flag that indicates the kind of data that is simulated.
	sim_flag = 0 -> Ones
	sim_flag = 1 -> Scale each beam by incrementing value i.e. beam 1 = 1, beam 2 = 2, ..., beam 64 = 64
	sim_flag = 2 -> Scale each beam by incrementing value in a particular bin (bin 3 and 6 for now). Match simulated data sim_flag = 2
	sim flag = 3 -> Simulated beams from 58 to 122 degrees. Assuming a ULA.
	*/
	int sim_flag = 0;
	if (sim_flag == 0) {
		for (int i = 0; i < (N_COEFF / 2); i++) {
			coeff_sim[2 * i] = 1;
		}
	}
	if (sim_flag == 1) {
		int tmp = 0;
		
		for (int a = 0; a < N_ANT; a++) {
			for (int b = 0; b < N_BEAM; b++) {
				if (tmp >= N_BEAM) {
					tmp = 0;
				}
				tmp = (tmp + 1) % (N_BEAM + 1);
				coeff_sim[2 * coeff_idx(a, b)] = tmp;
			}
		}
		
	}
	if (sim_flag == 2) {
		int tmp = 0;
		for (int a = 0; a < N_ANT; a++) {
			for (int b = 0; b < N_BEAM; b++) {
				if (tmp >= N_BEAM) {
					tmp = 0;
				}
				tmp = (tmp + 1) % (N_BEAM + 1);
				coeff_sim[2 * coeff_idx(a, b)] = tmp;
			}
		}
	}
	if (sim_flag == 3) {
		float c = 3e8; // Speed of light
		float c_freq = 1.25e9; // Center frequency
		float lambda = c / c_freq; // Wavelength
		float d = lambda / 2; // Distance between antennas
		//float chan_band = 1.59; // Fine channel bandwidth in Hz

		//float* rf_freqs = (float*)calloc(N_FREQ, sizeof(float));
		//for (int i = 0; i < N_FREQ; i++) {
		//	rf_freqs[i] = chan_band * i + c_freq;
		//}

		//float* theta = (float*)calloc(N_TIME, sizeof(float)); // Beam angle from 58 to 122 degrees
		//float* tau_beam = (float*)calloc(N_BEAM, sizeof(float)); // Delay

		float theta = 0; // Beam angle from 58 to 122 degrees
		float tau_beam = 0; // Delay

		for (int b = 0; b < N_BEAM; b++) {
			theta = ((b - (N_BEAM / 2)) + 90)*PI/180; // Beam angle from 58 to 122 degrees - Given SOI at 90 deg or moving across array, the beam with the most power is beamm 33
			tau_beam = d * cos(theta) / c; // Delay
			for (int a = 0; a < N_ANT; a++) {
				coeff_sim[2 * coeff_idx(a, b)] = cos(2 * PI * c_freq * a * tau_beam);
				coeff_sim[2 * coeff_idx(a, b) + 1] = sin(2 * PI * c_freq * a * tau_beam);
			}
		}
	}
	//for (int b = 0; b < N_BEAM; b++){
	//	printf("Beam %d = %f\n", (b+1), coeff_sim[2 * coeff_idx(0, b)]);
	//}

	return coeff_sim;
}

// Generate weights or coefficients with calculated delays (with delay polynomials (tau), coarse frequency channel (coarse_chan), and epoch (t))
float* generate_coefficients(float* tau, float coarse_chan, float t) {
	float* coefficients;
	coefficients = (float*)calloc(N_COEFF, sizeof(float));
	float delay_rate = 0;
	float delay_offset = 0;

	for (int b = 0; b < N_BEAM; b++) {
		for (int a = 0; a < N_ANT; a++) {
			delay_rate = tau[delay_idx(1,a,b)];
			delay_offset = tau[delay_idx(0,a,b)];
			coefficients[2 * coeff_idx(a, b)] = cos(2 * PI * coarse_chan * (t*delay_rate + delay_offset));
			coefficients[2 * coeff_idx(a, b) + 1] = sin(2 * PI * coarse_chan * (t*delay_rate + delay_offset));
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

// The coefficient_pin() function uses cudaHostRegister() to allocate the input host
// array in pinned memory.
// This speeds up the cudaMemcpy() and enables implementation into HASHPIPE/RTOS.
//void coefficient_pin(float * coeff_pin) {
//	checkCuda(cudaHostRegister(coeff_pin, N_COEFF*sizeof(float), cudaHostRegisterPortable));
//}

// Unregister host arrays from pinned memory
void unregister_data(void * data_unregister){
	checkCuda(cudaHostUnregister(data_unregister));
}

// Free memory
void cohbfCleanup() {
	// Free up GPU memory at the end of a program
	//if (h_data != NULL) {
	//	cudaFreeHost(h_data);
	//}
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
	//if (d_coeff_comp != NULL) {
	//	cudaFree(d_coeff_comp);
	//}
	//if (d_coh_bf_out != NULL) {
	//	cudaFree(d_coh_bf_out);
	//}
	if (d_coh_bf_pow != NULL) {
		cudaFree(d_coh_bf_pow);
	}
}

//Comment out main() function when compiling for hpguppi
/*// <----Uncomment here if testing standalone code
// Test all of the kernels and functions, and write the output to
// a text file for analysis
int main() {
	printf("Here!\n");

	// Allocate memory to all arrays used by run_beamformer() 
	init_beamformer();

	// Generate simulated data
	signed char* sim_data = simulate_data();
	// Register the array in pinned memory to speed HtoD mem copy
	input_data_pin(sim_data);

	// Generate simulated weights or coefficients
	float* sim_coefficients = simulate_coefficients();
	//printf("Here3!\n");
	// Register the array in pinned memory to speed HtoD mem copy
	coeff_pin(sim_coefficients);

	printf("real sim_data: %d and imag sim_data: %d\n", sim_data[10485768], sim_data[10485769]);
	//printf("real sim_coef: %f and imag sim_coef: %f\n", sim_coefficients[104], sim_coefficients[105]);

	printf("real sim_data2: %d and imag sim_data2: %d\n", sim_data[8388616], sim_data[8388617]);
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
		output_data = run_beamformer(sim_data, sim_coefficients);
		//run_beamformer(h_data, h_coeff, output_data);

		// Stop timing beamforming computation //
		clock_gettime(CLOCK_MONOTONIC, &tval_after);
		time_taken = (float)(tval_after.tv_sec - tval_before.tv_sec); //*1e6; // Time in seconds since epoch
		time_taken = time_taken + (float)(tval_after.tv_nsec - tval_before.tv_nsec)*1e-9; // Time in nanoseconds since 'tv_sec - start and end'
		bf_time += time_taken;
		//printf("Time taken: %f s\n", time_taken);
	}
	printf("Average delay calculation time: %f s\n", bf_time/num_runs);

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

	for (int ii = 0; ii < N_BF_POW; ii++) {
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
