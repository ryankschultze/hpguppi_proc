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
#include "upchannelizer_beamformer.h"

using namespace std;

// Perform transpose on the data and convert to floats
__global__
void data_transpose_ubf(signed char* data_in, cuComplex* data_tra, int offset, int n_pol, int n_chan, int n_ant, int n_ant_config, int n_samp, int n_win);

// Perform transpose on the output of the FFT
__global__
void fft_shift(cuComplex* data_in, cuComplex* data_tra, int offset, int n_ant_config, int n_pol, int n_coarse, int n_win, int n_fine);

// Perform beamforming operation
__global__
void coherent_beamformer_ubf(cuComplex* input_data, float* coeff, float* output_data, int offset, int n_ant_config, int n_pol, int n_fine, int n_coarse, int n_beam, int n_win);

// Compute power of beamformer output with STI (abs()^2)
__global__
void beamformer_power_sti_ubf(cuComplex* bf_volt, float* bf_power, int offset, int n_pol, int n_beam, int n_coarse, int n_fine, int n_time_int, int n_sti);

// Convenience function for checking CUDA runtime API results
// can be wrapped around any runtime API call. No-op in release builds.
inline
cudaError_t checkCuda_ubf(cudaError_t result)
{
  if (result != cudaSuccess) {
    fprintf(stderr, "CUDA Runtime Error: %s\n", 
            cudaGetErrorString(result));
    assert(result == cudaSuccess);
  }
  return result;
}

inline
cufftResult checkCufft(cufftResult result)
{
  if (result != CUFFT_SUCCESS) {
    //fprintf(stderr, "CUDA Runtime Error: %d\n", 
    //        result);
    switch(result){
      case CUFFT_INVALID_PLAN:
        fprintf(stderr, "CUDA Runtime Error: CUFFT_INVALID_PLAN \n");
	break;
      case CUFFT_ALLOC_FAILED:
        fprintf(stderr, "CUDA Runtime Error: CUFFT_ALLOC_FAILED \n");
	break;
      case CUFFT_INVALID_TYPE:
        fprintf(stderr, "CUDA Runtime Error: CUFFT_INVALID_TYPE \n");
	break;
      case CUFFT_INVALID_VALUE:
        fprintf(stderr, "CUDA Runtime Error: CUFFT_INVALID_VALUE \n");
	break;
      case CUFFT_INTERNAL_ERROR:
        fprintf(stderr, "CUDA Runtime Error: CUFFT_INTERNAL_ERROR \n");
	break;
      case CUFFT_EXEC_FAILED:
        fprintf(stderr, "CUDA Runtime Error: CUFFT_EXEC_FAILED \n");
	break;
      case CUFFT_SETUP_FAILED:
        fprintf(stderr, "CUDA Runtime Error: CUFFT_SETUP_FAILED \n");
	break;
      case CUFFT_INVALID_SIZE:
        fprintf(stderr, "CUDA Runtime Error: CUFFT_INVALID_SIZE \n");
	break;
      case CUFFT_UNALIGNED_DATA:
        fprintf(stderr, "CUDA Runtime Error: CUFFT_UNALIGNED_DATA \n");
	break;
    }
    assert(result == CUFFT_SUCCESS);
  }
  return result;
}

signed char* d_data_char = NULL;
cuComplex* d_data_comp = NULL;
cuComplex* d_data_shift = NULL;
float* d_coeff = NULL;
//float* d_coh_bf_pow = NULL;
float* h_bf_pow = NULL;

// Allocate memory to all arrays 
void init_upchan_beamformer(int telescope_flag) {
	printf("Here In init_upchan_beamformer()! \n");

	if(telescope_flag == 0){
		// Allocate memory for input data float type
		checkCuda_ubf(cudaMalloc((void **)&d_data_char, (N_INPUT) * sizeof(signed char)));
		printf("Here 1st cudaMalloc! \n");

		// Allocate memory for input data cuComplex type
		checkCuda_ubf(cudaMalloc((void **)&d_data_comp, (N_INPUT) * sizeof(cuComplex) / 2));
		printf("Here 2nd cudaMalloc! \n");

		// Allocate memory for data with FFT shift cuComplex type
		checkCuda_ubf(cudaMalloc((void **)&d_data_shift, (N_INPUT) * sizeof(cuComplex) / 2));
		printf("Here 3rd cudaMalloc! \n");

		// Allocate memory for coefficients float type
		checkCuda_ubf(cudaMalloc((void **)&d_coeff, N_COEFF * sizeof(float)));
		printf("Here 4th cudaMalloc! \n");

		// Allocate memory for output power of coherent beamformer
	        //checkCuda_ubf(cudaMalloc((void **)&d_coh_bf_pow, (N_BF_POW) * sizeof(float)));
		//printf("Here 5th cudaMalloc! \n");

		checkCuda_ubf(cudaMallocHost((void **)&h_bf_pow, (N_BF_POW) * sizeof(float)));
	}else if(telescope_flag == 1){
		// Allocate memory for input data float type
		checkCuda_ubf(cudaMalloc((void **)&d_data_char, (VLASS_N_INPUT) * sizeof(signed char)));
		printf("Here 1st cudaMalloc! \n");

		// Allocate memory for input data cuComplex type
		checkCuda_ubf(cudaMalloc((void **)&d_data_comp, (VLASS_N_INPUT) * sizeof(cuComplex) / 2));
		printf("Here 2nd cudaMalloc! \n");

		// Allocate memory for data with FFT shift cuComplex type
		checkCuda_ubf(cudaMalloc((void **)&d_data_shift, (VLASS_N_INPUT) * sizeof(cuComplex) / 2));
		printf("Here 3rd cudaMalloc! \n");

		// Allocate memory for coefficients float type
		checkCuda_ubf(cudaMalloc((void **)&d_coeff, VLASS_N_COEFF * sizeof(float)));
		printf("Here 4th cudaMalloc! \n");

		// Allocate memory for output power of coherent beamformer
	        //checkCuda_ubf(cudaMalloc((void **)&d_coh_bf_pow, (VLASS_N_BF_POW) * sizeof(float)));
		//printf("Here 5th cudaMalloc! \n");

		checkCuda_ubf(cudaMallocHost((void **)&h_bf_pow, (VLASS_N_BF_POW) * sizeof(float)));
	}


	return;
}

// Set arrays to zero after a block is processed
void set_to_zero_ubf(int telescope_flag){
	if(telescope_flag == 0){
		checkCuda_ubf(cudaMemset(d_data_comp, 0, (N_INPUT) * sizeof(cuComplex)/2));
	}else if(telescope_flag == 1){
		checkCuda_ubf(cudaMemset(d_data_comp, 0, (VLASS_N_INPUT) * sizeof(cuComplex)/2));
	}
}

// Set arrays to zero after a block is processed
void set_second_to_zero(int telescope_flag){
	if(telescope_flag == 0){
		checkCuda_ubf(cudaMemset(d_data_shift, 0, (N_INPUT) * sizeof(cuComplex)/2));
	}else if(telescope_flag == 1){
		checkCuda_ubf(cudaMemset(d_data_shift, 0, (VLASS_N_INPUT) * sizeof(cuComplex)/2));
	}
}
/*
// Perform transpose on the data and convert to floats
__global__
void data_transpose_ubf(signed char* data_in, cuComplex* data_tra, int offset, int n_pol, int n_chan, int n_ant, int n_ant_config, int n_samp) {
	int t = threadIdx.x; // Time sample index
	int a = blockIdx.x;  // Antenna index
	int w = blockIdx.y;  // Time window index
	int c = blockIdx.z;  // Coarse channel index
        int p = 0;           // Polarization index

	int tb = 0; // Index for block of time samples to compensate max number of threads
	int TS = n_samp/MAX_THREADS; // Number of blocks of time samples to process

	for(p=0; p<n_pol; p++){
		for(tb = 0; tb < TS; tb++){
		// If the input data is not float e.g. signed char, just multiply it by '1.0f' to convert it to a float
                        int h_in = data_in_idx(p, (c + offset), a, t + tb*MAX_THREADS, w, n_pol, n_chan, n_ant, n_samp); 
			int h_tr = data_tr_idx(t + tb*MAX_THREADS, a, p, (c + offset), w, n_samp, n_ant_config, n_pol, n_chan);

			data_tra[h_tr].x = data_in[2*h_in]*1.0f;
			data_tra[h_tr].y = data_in[2*h_in + 1]*1.0f;
		}
	}

	return;
}
*/
// Perform transpose on the data and convert to floats
__global__
void data_transpose_ubf(signed char* data_in, cuComplex* data_tra, int offset, int n_pol, int n_chan, int n_ant, int n_ant_config, int n_samp, int n_win) {
	int t = threadIdx.x; // Time sample index
	int a = blockIdx.x;  // Antenna index
        int w = blockIdx.y;  // Time window index
	int c = blockIdx.z;  // Coarse channel index
	int p = 0;           // Polarization index

	int tb = 0; // Index for block of time samples to compensate max number of threads
	int TS = n_samp/MAX_THREADS; // Number of blocks of time samples to process

	for(p=0; p<n_pol; p++){
		for(tb = 0; tb < TS; tb++){
			// If the input data is not float e.g. signed char, just multiply it by '1.0f' to convert it to a float
			//int h_in = data_in_idx(p, (c + offset), a, t + tb*MAX_THREADS, w, n_pol, n_chan, n_ant, n_samp);
			int h_in = data_in_idx(p, t + tb*MAX_THREADS, w, a, (c + offset), n_pol, n_samp, n_win, n_ant);
			int h_tr = data_tr_idx(t + tb*MAX_THREADS, a, p, (c + offset), w, n_samp, n_ant_config, n_pol, n_chan);

			data_tra[h_tr].x = data_in[2*h_in]*1.0f;
			data_tra[h_tr].y = data_in[2*h_in + 1]*1.0f;
		}
	}

	return;
}

// Perform FFT
void upchannelize(complex_t* data_tra, int n_ant_config, int n_pol, int n_chan, int n_win, int n_samp){
        cufftHandle plan;

	/*
	// Determine worksize then create/make plane
	size_t worksize;
	checkCufft(cufftEstimate1d(n_samp, CUFFT_C2C, BATCH(n_pol), &worksize));
	checkCufft(cufftCreate(&plan));
	checkCufft(cufftMakePlan1d(plan, n_samp, CUFFT_C2C, BATCH(n_pol), &worksize));
	*/

	// Setup the cuFFT plan
        checkCufft(cufftPlan1d(&plan, n_samp, CUFFT_C2C, BATCH(n_ant_config,n_pol)));
    	// Execute a complex-to-complex 1D FFT
	// Reduce the amount of memory utilized by cufft by iterating through coarese channels (c) and time windows (w)
	int h = 0;
	for(int w = 0; w < n_win; w++){
		for(int c = 0; c < n_chan; c++){
			//printf("w = %d, c = %d\n", w,c);
			h = data_tr_idx(0, 0, 0, c, w, n_samp, n_ant_config, n_pol, n_chan);
                	checkCufft(cufftExecC2C(plan, (cufftComplex*)&data_tra[h], (cufftComplex*)&data_tra[h], CUFFT_FORWARD));
		}
	}
	checkCufft(cufftDestroy(plan));
}


// This kernel should only be used in the standalone code to test the FFT shift.
// If used in the upchannelized beamformer code, the amount of memory allocated will be too large
// The FFT shift will be performed during the beamforming computation to reduce the amount of memory allocated on the device (GPU)
// Perform transpose on the output of the FFT
__global__
void fft_shift(cuComplex* data_in, cuComplex* data_tra, int offset, int n_ant_config, int n_pol, int n_coarse, int n_win, int n_fine) {
        // 'f' is the largest dimension and is sometimes larger than 1024 which is the max number of threads
        // So 'f' should be the blockIdx.x which has the largest max value (over 2e9 elements)
        // Since 'f' has to be the fasted moving index for the cufftExecC2C(), but 'a' needs to be the fastest moving index
        // for beamforming, then threadIdx.x can't be used for 'a' in this case, and threadIdx.y must be used instead.
        // blockIdx.x can't be a faster moving dimension than threadIdx.x
	int a = threadIdx.y; // Antenna index
	int p = threadIdx.z; // Polarization index
	int f = blockIdx.x;  // Fine channel index
	int c = blockIdx.y;  // Coarse channel index
	int w = blockIdx.z;  // Time window index

        int h_in = 0;
	int h_sh = 0;

	if(f < (n_fine/2)){
		h_in = data_fft_out_idx(f, a, p, (c + offset), w, n_fine, n_ant_config, n_pol, n_coarse); 
		h_sh = data_fftshift_idx(a, p, (f+(n_fine/2)), (c + offset), w, n_ant_config, n_pol, n_fine, n_coarse);

		data_tra[h_sh].x = data_in[h_in].x;
		data_tra[h_sh].y = data_in[h_in].y;
	}else if((f >= (n_fine/2)) && (f < n_fine)){
		h_in = data_fft_out_idx(f, a, p, (c + offset), w, n_fine, n_ant_config, n_pol, n_coarse);
		h_sh = data_fftshift_idx(a, p, (f-(n_fine/2)), (c + offset), w, n_ant_config, n_pol, n_fine, n_coarse);

		data_tra[h_sh].x = data_in[h_in].x;
		data_tra[h_sh].y = data_in[h_in].y;
	}

	return;
}

// Perform beamforming operation
__global__
void coherent_beamformer_ubf(cuComplex* input_data, float* coeff, cuComplex* output_data, int offset, int n_ant_config, int n_pol, int n_fine, int n_coarse, int n_beam, int n_win) {
	int a = threadIdx.x; // Antenna index
	int f = blockIdx.x;  // Fine channel index
	int c = blockIdx.y;  // Coarse channel index
	int b = blockIdx.z;  // Beam index

	__shared__ cuFloatComplex reduced_mul[N_ANT];

	int n_freq_streams = n_coarse/N_STREAMS;

	if(c < n_freq_streams){
		for (int p = 0; p < n_pol; p++) { // Polarization index
			for (int s = 0; s < n_win; s++) { // STI window index

				int i = data_fftshift_idx(a, p, f, (c + offset), s, n_ant_config, n_pol, n_fine, n_coarse);
	                        int w = coeff_idx(a, p, b, (c + offset), n_ant_config, n_pol, n_beam);

				if (a < n_ant_config) {
	       	                         // Conjugated coefficients places the minus operation in the imaginary component
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
					int h1 = coh_bf_idx(p, b, f, (c + offset), s, n_pol, n_beam, n_coarse, n_fine);
					output_data[h1].x = reduced_mul[0].x;
					output_data[h1].y = reduced_mul[0].y;
				}
			}
		}
	}
	return;
}

// Compute power of beamformer output (abs()^2)
__global__
void beamformer_power_sti_ubf(cuComplex* bf_volt, float* bf_power, int offset, int n_pol, int n_beam, int n_coarse, int n_fine, int n_time_int, int n_sti) {
	int t = threadIdx.x; // Time sample index
	int f = blockIdx.x;  // Fine channel index
	int c = blockIdx.y;  // Coarse channel index
	int b = blockIdx.z;  // Beam index
	int s = 0;           // STI window index
	int h = 0;
	int xp = 0; // X polarization
	int yp = 0; // Y polarization

	int n_freq_streams = n_coarse/N_STREAMS;
	
	float x_pol_pow; // XX*
	float y_pol_pow; // YY*
	float beam_power;
	//float scale = 1.0; // /n_sti;

	__shared__ float reduced_array[N_STI_BLOC];

	for(s = 0; s<n_sti; s++){
		xp = coh_bf_idx(0, b, f, (c + offset), (s*n_time_int + t), n_pol, n_beam, n_coarse, n_fine); // X polarization
		yp = coh_bf_idx(1, b, f, (c + offset), (s*n_time_int + t), n_pol, n_beam, n_coarse, n_fine); // Y polarization
		h = pow_bf_idx(f, (c + offset), s, b, n_fine, n_coarse, n_sti);
		if(c < n_freq_streams){	
			if(n_pol == 1){
				if (t < n_time_int) {
					// Power = Absolute value squared of output -> r^2 + i^2
					x_pol_pow = (bf_volt[xp].x * bf_volt[xp].x) + (bf_volt[xp].y * bf_volt[xp].y); // XX*

					beam_power = x_pol_pow; // XX*
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
					bf_power[h] = reduced_array[0];
				}
			}else if(n_pol == 2){
				if (t < n_time_int) {
					// Power = Absolute value squared of output -> r^2 + i^2
					x_pol_pow = (bf_volt[xp].x * bf_volt[xp].x) + (bf_volt[xp].y * bf_volt[xp].y); // XX*
					y_pol_pow = (bf_volt[yp].x * bf_volt[yp].x) + (bf_volt[yp].y * bf_volt[yp].y); // YY*
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
					bf_power[h] = reduced_array[0];
				}
			}
		}
	}

	return;
}

// Run upchannelizer and beamformer
float* run_upchannelizer_beamformer(signed char* data_in, float* h_coefficient, int n_pol, int n_ant, int n_beam, int n_chan, int n_win, int n_time_int, int n_samp, int telescope_flag) {

	cudaError_t err_code;

        // Total number of time samples
        int nt = n_win*n_samp;

	int n_sti = n_win/n_time_int; // Number of STI windows

        int n_ant_config = 0; // Number of antennas depending on configuration

	// If the number of antennas is half of N_ANT (subarray configuration),
	// then the antenna dimension needs to be half of N_ANT
	if(n_ant <= N_ANT/2){
		n_ant_config = N_ANT/2;
	}else{
		n_ant_config = N_ANT;
	}

	//dim3 dimBlock_transpose(n_samp, 1, 1);
	dim3 dimBlock_transpose(MAX_THREADS, 1, 1);
	dim3 dimGrid_transpose(n_ant_config, n_win, n_chan);

	// FFT shift kernel: Specify grid and block dimensions
	dim3 dimBlock_fftshift(1, n_ant_config, n_pol);
	dim3 dimGrid_fftshift(n_samp, n_chan, n_win);

	// Coherent beamformer kernel: Specify grid and block dimensions
	dim3 dimBlock_coh_bf(n_ant_config, 1, 1);
	dim3 dimGrid_coh_bf(n_samp, n_chan, n_beam);

	// Output power of beamformer kernel with STI: Specify grid and block dimensions
	dim3 dimBlock_bf_pow(N_STI_BLOC, 1, 1);
	dim3 dimGrid_bf_pow(n_samp, n_chan, n_beam);

	signed char* d_data_in = d_data_char;
	cuComplex* d_data_tra = d_data_comp;
	cuComplex* d_data_tra2 = d_data_shift;
	float* d_coefficient = d_coeff;
	//float* d_bf_pow = d_coh_bf_pow;
	float* data_out = h_bf_pow;

	// Copy beamformer coefficients from host to device
	checkCuda_ubf(cudaMemcpy(d_coefficient, h_coefficient, (2*n_pol*n_ant_config*n_beam*n_chan) * sizeof(float), cudaMemcpyHostToDevice));

	//printf("Before cudaMemcpy(HtoD) coefficients! \n");
	// Copy input data from host to device
	checkCuda_ubf(cudaMemcpy(d_data_in, data_in, 2*n_ant_config*n_pol*nt*n_chan*sizeof(signed char), cudaMemcpyHostToDevice));

        // Perform transpose on the data and convert to floats  
        data_transpose_ubf<<<dimGrid_transpose, dimBlock_transpose>>>(d_data_in, d_data_tra, 0, n_pol, n_chan, n_ant, n_ant_config, n_samp, n_win);
        err_code = cudaGetLastError();
	if (err_code != cudaSuccess) {
		printf("FFT: data_transpose() kernel Failed: %s\n", cudaGetErrorString(err_code));
	}

        // Upchannelize the data
        upchannelize((complex_t*)d_data_tra, n_ant_config, n_pol, n_chan, n_win, n_samp);

	// FFT shift and transpose
	fft_shift<<<dimGrid_fftshift, dimBlock_fftshift>>>(d_data_tra, d_data_tra2, 0, n_ant_config, n_pol, n_chan, n_win, n_samp);
	err_code = cudaGetLastError();
	if (err_code != cudaSuccess) {
		printf("FFT: fft_shift() kernel Failed: %s\n", cudaGetErrorString(err_code));
	}

	// Set input of FFT shift to zero so it can be used as the output of the coherent_beamformer
	set_to_zero_ubf(telescope_flag);

	// Coherent beamformer
	coherent_beamformer_ubf<<<dimGrid_coh_bf, dimBlock_coh_bf>>>(d_data_tra2, d_coefficient, d_data_tra, 0, n_ant_config, n_pol, n_samp, n_chan, n_beam, n_win);
	err_code = cudaGetLastError();
	if (err_code != cudaSuccess) {
		printf("FFT: coherent_beamformer_ubf() kernel Failed: %s\n", cudaGetErrorString(err_code));
	}

	// Set input of coherent_beamformer_ubf() to zero so it can be used as the output of the beamformer_power_sti_ubf() kernel
	set_second_to_zero(telescope_flag);

	// Short time integration after beamforming
	beamformer_power_sti_ubf<<<dimGrid_bf_pow, dimBlock_bf_pow>>>(d_data_tra, (float*)d_data_tra2, 0, n_pol, n_beam, n_chan, n_samp, n_time_int, n_sti);
	err_code = cudaGetLastError();
	if (err_code != cudaSuccess) {
		printf("FFT: beamformer_power_sti() kernel Failed: %s\n", cudaGetErrorString(err_code));
	}

	// Copy input data from device to host
	checkCuda_ubf(cudaMemcpy(data_out, (float*)d_data_tra2, n_beam*n_sti*n_samp*n_chan*sizeof(float), cudaMemcpyDeviceToHost));

        return data_out;
}

// Generate simulated data
signed char* simulate_data_ubf(int n_sim_ant, int nants, int n_pol, int n_chan, int nt, int n_win, int telescope_flag) {
	int n_input = 0;
	int n_ant_config = 0;
	if(telescope_flag == 0){
		n_input = N_INPUT;
		n_ant_config = N_ANT;
	}else if(telescope_flag == 1){
		n_input = VLASS_N_INPUT;
		n_ant_config = N_ANT/2;
	}

	signed char* data_sim;
	data_sim = (signed char*)calloc(n_input, sizeof(signed char));

	/*
	'sim_flag' is a flag that indicates the kind of data that is simulated.
	sim_flag = 0 -> Ones
	sim_flag = 1 -> Ones placed in a particular bin (bin 3 for now)
	sim_flag = 2 -> Ones placed in a particular bin at a particular antenna (bin 3 and antenna 3 for now)
	sim_flag = 3 -> Rect placed in a particular bin at a particular antenna (bin 3 and antenna 3 for now)
	sim flag = 4 -> Simulated cosine wave
	sim flag = 5 -> Simulated complex exponential i.e. exp(j*2*pi*f0*t)
	*/
	int sim_flag = 5;
	if (sim_flag == 0) {
		for (int i = 0; i < (n_input / 2); i++) {
			if(i < ((n_sim_ant*n_input)/(2*n_ant_config))){
				data_sim[2 * i] = 1;
			}else{
				data_sim[2 * i] = 0;
			}
		}
	}
	if (sim_flag == 1) {
		for (int p = 0; p < n_pol; p++) {
			for (int t = 0; t < nt; t++) {
				for (int a = 0; a < nants; a++) {
					if(a < n_sim_ant){ 
						data_sim[2 * data_in_idx(p, t, 0, a, 2, n_pol, nt, n_win, nants)] = 1;
					}
				}
			}
		}
	}
	if (sim_flag == 2) {
		for (int p = 0; p < n_pol; p++) {
			for (int t = 0; t < nt; t++) {
				data_sim[2 * data_in_idx(p, t, 0, 2, 2, n_pol, nt, n_win, nants)] = 1;
			}
		}
	}
	if (sim_flag == 3) {
		// data_in_idx(p, t, w, c, a, Np, Nt, Nw, Nc)
		for (int p = 0; p < n_pol; p++) {
			for (int t = (1024*10); t < (nt-(1024*10)); t++) {
				//data_sim[2 * data_in_idx(p, 0, 2, t, 0, n_pol, n_chan, nants, nt)] = 1;
				//data_sim[2 * data_in_idx(p, 1, 2, t, 0, n_pol, n_chan, nants, nt)] = 1;
				data_sim[2 * data_in_idx(p, t, 0, 2, 2, n_pol, nt, n_win, nants)] = 1;
				//data_sim[2 * data_in_idx(p, 3, 2, t, 0, n_pol, n_chan, nants, nt)] = 1;

			}
		}
	}
	if (sim_flag == 4) {
		float freq = 1e3; // Resonant frequency

                float tmp_max = 1.0;
		float tmp_min = -1.0;

		for (int t = 0; t < nt; t++) {
			for (int f = 0; f < n_chan; f++) {
				for (int a = 0; a < nants; a++) {
					if(a < n_sim_ant){
						// Requantize from doubles/floats to signed chars with a range from -128 to 127 
						// X polarization
						data_sim[2 * data_in_idx(0, t, 0, a, f, n_pol, nt, n_win, nants)] = (signed char)((((cos(2 * PI * freq * t*0.000001) - tmp_min)/(tmp_max-tmp_min)) - 0.5)*256);
						//data_sim[2 * data_in_idx(0, f, a, t, 0, n_pol, n_chan, nants, nt) + 1] = 0;
						// Y polarization
						data_sim[2 * data_in_idx(1, t, 0, a, f, n_pol, nt, n_win, nants)] = (signed char)((((2*cos(2 * PI * freq * t*0.000001) - tmp_min)/(tmp_max-tmp_min)) - 0.5)*256);
						//data_sim[2 * data_in_idx(1, f, a, t, 0, n_pol, n_chan, nants, nt) + 1] = 0;
					}
				}
			}
		}
	}
	if (sim_flag == 5) {
		float freq = 1e3; // Resonant frequency

                float tmp_max = 1.0;
		float tmp_min = -1.0;

		for (int w = 0; w < n_win; w++) {
		for (int t = 0; t < nt; t++) {
			for (int f = 0; f < n_chan; f++) {
				for (int a = 0; a < nants; a++) {
					if(a < n_sim_ant){
						// Requantize from doubles/floats to signed chars with a range from -128 to 127 
						// X polarization
						data_sim[2 * data_in_idx(0, t, w, a, f, n_pol, nt, n_win, nants)] = (signed char)((((cos(2 * PI * freq * t*0.000001) - tmp_min)/(tmp_max-tmp_min)) - 0.5)*256);
						data_sim[2 * data_in_idx(0, t, w, a, f, n_pol, nt, n_win, nants) + 1] = (signed char)((((sin(2 * PI * freq * t*0.000001) - tmp_min)/(tmp_max-tmp_min)) - 0.5)*256);
						//data_sim[2 * data_in_idx(0, f, a, t, 0, n_pol, n_chan, nants, nt) + 1] = 0;
						// Y polarization
						data_sim[2 * data_in_idx(1, t, w, a, f, n_pol, nt, n_win, nants)] = (signed char)((((2*cos(2 * PI * freq * t*0.000001) - tmp_min)/(tmp_max-tmp_min)) - 0.5)*256);
						data_sim[2 * data_in_idx(1, t, w, a, f, n_pol, nt, n_win, nants) + 1] = (signed char)((((sin(2 * PI * freq * t*0.000001) - tmp_min)/(tmp_max-tmp_min)) - 0.5)*256);
						//data_sim[2 * data_in_idx(1, f, a, t, 0, n_pol, n_chan, nants, nt) + 1] = 0;
					}
				}
			}
		}
		}
	}
	return data_sim;
}

// Generate simulated weights or coefficients
float* simulate_coefficients_ubf(int n_sim_ant, int nants, int n_pol, int n_beam, int n_chan, int telescope_flag) {
	int n_coeff = 0;
	int n_ant_config = 0;
	if(telescope_flag == 0){
		n_coeff = N_COEFF;
		n_ant_config = N_ANT;
	}else if(telescope_flag == 1){
		n_coeff = VLASS_N_COEFF;
		n_ant_config = N_ANT/2;
	}

	float* coeff_sim;
	coeff_sim = (float*)calloc(n_coeff, sizeof(float));
	/*
	'sim_flag' is a flag that indicates the kind of data that is simulated.
	sim_flag = 0 -> Ones
	sim_flag = 1 -> Scale each beam by incrementing value i.e. beam 1 = 1, beam 2 = 2, ..., beam 64 = 64
	sim_flag = 2 -> Scale each beam by incrementing value in a particular bin (bin 3 and 6 for now). Match simulated data sim_flag = 2
	sim flag = 3 -> Simulated beams from 58 to 122 degrees. Assuming a ULA.
        sim_flag = 4 -> One value at one polarization, one element, one beam, and one frequency bin
	*/
	int sim_flag = 0;
	if (sim_flag == 0) {
		for (int i = 0; i < (n_pol*n_ant_config*n_beam*n_chan); i++) {
			coeff_sim[2*i] = 1;
			//coeff_sim[2*i + 1] = 1;
		}
	}
	if (sim_flag == 1) {
		int tmp = 0;
		
                for (int p = 0; p < n_pol; p++) {
                	for (int f = 0; f < n_chan; f++) {
				for (int a = 0; a < nants; a++) {
					for (int b = 0; b < n_beam; b++) {
						if (tmp >= n_beam) {
							tmp = 0;
						}
                                		coeff_sim[2 * coeff_idx(a, p, b, f, nants, n_pol, n_beam)] = tmp*0.01;
						tmp = (tmp + 1) % (n_beam + 1);
					}
				}
                	}
                }
		
	}
	if (sim_flag == 2) {
		int tmp = 0;
		for (int p = 0; p < n_pol; p++) {
                	for (int f = 0; f < n_chan; f++) {
				for (int a = 0; a < nants; a++) {
					for (int b = 0; b < n_beam; b++) {
						if (tmp >= n_beam) {
							tmp = 0;
						}
						tmp = (tmp + 1) % (n_beam + 1);
                                		coeff_sim[2 * coeff_idx(a, p, b, f, nants, n_pol, n_beam)] = tmp;
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

		int b = 2;
		for (int f = 0; f < n_chan; f++) {
			rf_freqs = chan_band * f + c_freq;
			theta = ((b - (n_beam / 2)) + 90)*PI/180; // Beam angle from 58 to 122 degrees - Given SOI at 90 deg or moving across array, the beam with the most power is beamm 33
			tau_beam = d * cos(theta) / c; // Delay
			for (int a = 0; a < nants; a++) {
                        	if(n_pol == 1){
					// X polarization
					coeff_sim[2 * coeff_idx(a, 0, b, f, nants, n_pol, n_beam)] = cos(2 * PI * rf_freqs * a * tau_beam);
					coeff_sim[2 * coeff_idx(a, 0, b, f, nants, n_pol, n_beam) + 1] = sin(2 * PI * rf_freqs * a * tau_beam);
				}else if(n_pol == 2){
					// X polarization
					coeff_sim[2 * coeff_idx(a, 0, b, f, nants, n_pol, n_beam)] = cos(2 * PI * rf_freqs * a * tau_beam);
					coeff_sim[2 * coeff_idx(a, 0, b, f, nants, n_pol, n_beam) + 1] = sin(2 * PI * rf_freqs * a * tau_beam);
					// Y polarization
					coeff_sim[2 * coeff_idx(a, 1, b, f, nants, n_pol, n_beam)] = cos(2 * PI * rf_freqs * a * tau_beam);
					coeff_sim[2 * coeff_idx(a, 1, b, f, nants, n_pol, n_beam) + 1] = sin(2 * PI * rf_freqs * a * tau_beam);
				}
			}
		}
	}
        if (sim_flag == 4) {
		for(int a = 0; a < nants; a++){
			for(int f = 0; f < n_chan; f++){
		        	coeff_sim[2 * coeff_idx(a, 0, 2, f, nants, n_pol, n_beam)] = 1;
			}
		}
        }

	return coeff_sim;
}

// Generate coefficients with delays and phase up solutions from HDF5 file
float* generate_coefficients_ubf(complex_t* phase_up, double* delay, int n, double* coarse_chan, int n_ant_config, int n_pol, int n_beam, int actual_n_beam, int schan, int n_coarse, int subband_idx, uint64_t n_real_ant, int telescope_flag) {
	int n_coeff = 0;
	if(telescope_flag == 0){
		n_coeff = N_COEFF;
	}else if(telescope_flag == 1){
		n_coeff = VLASS_N_COEFF;
	}

	float* coefficients;
	coefficients = (float*)calloc(n_coeff, sizeof(float));
	double tau = 0;
        int fc = 0; // First frequency channel in the RAW file and on this node

	for (int p = 0; p < n_pol; p++) {
                // 'schan' is the absolute channel index of the first channel in the RAW file.
                // The beamformer recipe files contain all of the channels in the band
                // So 'schan' offsets to start processing with the correct section/range of frequency channels
		//for (int f = subband_idx*n_coarse; f < ((subband_idx*n_coarse) + n_coarse); f++) {
		for (int f = 0; f < n_coarse; f++) {
			for (int b = 0; b < n_beam; b++) {
				for (int a = 0; a < n_ant_config; a++) {
					if(a < n_real_ant){
						tau = delay[delay_idx(a, b, n, n_real_ant, actual_n_beam)];
                                                fc = (f + (subband_idx*n_coarse))+schan; // First frequency channel in the RAW file and on this node

						coefficients[2 * coeff_idx(a, p, b, f, n_ant_config, n_pol, n_beam)] = (float)(phase_up[cal_all_idx(a, p, fc, n_real_ant, n_pol)].re*cos(2 * PI * coarse_chan[f] * tau) - phase_up[cal_all_idx(a, p, fc, n_real_ant, n_pol)].im*sin(2 * PI * coarse_chan[f] * tau));
						coefficients[2 * coeff_idx(a, p, b, f, n_ant_config, n_pol, n_beam) + 1] = (float)(phase_up[cal_all_idx(a, p, fc, n_real_ant, n_pol)].re*sin(2 * PI * coarse_chan[f] * tau) + phase_up[cal_all_idx(a, p, fc, n_real_ant, n_pol)].im*cos(2 * PI * coarse_chan[f] * tau));
					}else{
						coefficients[2 * coeff_idx(a, p, b, f, n_ant_config, n_pol, n_beam)] = 0;
						coefficients[2 * coeff_idx(a, p, b, f, n_ant_config, n_pol, n_beam) + 1] = 0;
					}
				}
			}
		}
	}

	return coefficients;
}

// Free memory
void Cleanup_beamformer() {
	// Free up GPU memory at the end of a program
	if (d_data_char != NULL) {
		cudaFree(d_data_char);
	}
	if (d_data_comp != NULL) {
		cudaFree(d_data_comp);
	}
	if (d_data_shift != NULL) {
		cudaFree(d_data_shift);
	}
	if (d_coeff != NULL) {
		cudaFree(d_coeff);
	}
	//if (d_coh_bf_pow != NULL) {
	//	cudaFree(d_coh_bf_pow);
	//}
}


// Testing to see whether input data is as I expect it to be
float* data_test(signed char *sim_data){
	float* data_float;
	data_float = (float*)calloc(N_INPUT, sizeof(float));
	for (int ii = 0; ii < N_INPUT; ii++) { // Write up to the size of the data corresponding to 1k, 4k or 32k mode
		data_float[ii] = sim_data[ii]*1.0f;
	}
	return data_float;
}

//Comment out main() function when compiling for hpguppi
// <----Uncomment here if testing standalone code
// Test all of the kernels and functions, and write the output to
// a text file for analysis
int main() {
	printf("Here!\n");
	// ---------------------------- //
	// To run in regular array configuration, enter values between 33 and 64 in n_beam and n_ant
	// To run in subarray configuration, enter values 32 or less (and greater than 1 otherwise, beamforming can't be done)
	// ---------------------------- //
	int n_beam = 0;
        int n_pol = 0;
	int n_sim_ant = 0;
	int n_ant_config = 0;
	int n_chan = 0;
	int nt = 0;
        int n_win = 0;
	int n_time_int = 0;
	int n_input = 0;

	int telescope_flag = 0;
	int spec_flag = 0;
        // ---------------- MeerKAT specs --------------- //
	if(telescope_flag == 0){
		n_input = N_INPUT;
		n_beam = 61;
		n_pol = 2;
		n_sim_ant = 58;
		if(n_sim_ant <= N_ANT/2){ // Subarray configuration
			n_ant_config = N_ANT/2;
			// 5 seconds worth of processing at a time
			// 1k mode
			//n_chan = 1; 
		        //nt = 2*4096*1024; // 4194304; // 2^22
			// 4k mode
		    	n_chan = 4; // 64
		        nt = 2*1024*1024; // 1048576; // 2^20
			// 32k mode
			//n_chan = 32;
			//nt = 2*128*1024; // 131072; // 2^17

			n_win = 16;
			n_time_int = 16;
		}else{ // Regular array configuration
			n_ant_config = N_ANT;
			// 5 seconds worth of processing at a time
			// 1k mode
			//n_chan = 1; 
		        //nt = 4096*1024; // 4194304; // 2^22
			// 4k mode
			n_chan = 4; // 64
		        nt = 1024*1024; // 1048576; // 2^20
			// 32k mode
			//n_chan = 32;
			//nt = 128*1024; // 131072; // 2^17

			n_win = 8;
			n_time_int = 8;
		}
	}
	// -----------------------------------------------//
	// ------------------ VLASS specs ----------------//
	else if(telescope_flag == 1){
		n_input = VLASS_N_INPUT;
                // Required Specification
		if(spec_flag == 0){
			n_beam = 5;
			n_pol = 2;
			n_sim_ant = 27;
			n_ant_config = N_ANT/2;
			n_chan = 1;
			nt = 5120000;
			n_win = 40;
			n_time_int = 1;
		}// Desired Specification
		else if(spec_flag == 1){
			n_beam = 31;
			n_pol = 2;
			n_sim_ant = 27;
			n_ant_config = N_ANT/2;
			n_chan = 1;
			nt = 10240000; // 5120000
			n_win = 10000;
			n_time_int = 1;
		}
	}
	// -----------------------------------------------//

        int n_samp = nt/n_win;
	int n_sti = n_win/n_time_int;

	// Allocate memory to all arrays used by run_FFT() 
	init_upchan_beamformer(telescope_flag);

        printf("After init_upchan_beamformer() \n");

	// Generate simulated data
	signed char* sim_data = simulate_data_ubf(n_sim_ant, n_ant_config, n_pol, n_chan, n_samp, n_win, telescope_flag);

        printf("After simulate_data() \n");

	// Generate simulated weights or coefficients
	float* sim_coefficients = simulate_coefficients_ubf(n_sim_ant, n_ant_config, n_pol, n_beam, n_chan, telescope_flag);

	printf("After simulate_coefficients() \n");

	// --------------------- Input data test --------------------- //
	int input_write = 0; // If input_write is set to 1, the simulated data will be written to a binary file for testing/verification

	if(input_write == 1){
		float* input_test = data_test(sim_data);

		// Write data to binary file for analysis
		char input_filename[128];

		printf("Here1!\n");

		strcpy(input_filename, "/mydatag/Unknown/GUPPI/input_h_fft_bf.bin");

		printf("Here2!\n");

		FILE* input_file;

		printf("Here3!\n");

		input_file = fopen(input_filename, "w");

		printf("Here4!\n");

		fwrite(input_test, sizeof(float), n_input, input_file);

		printf("Here5!\n");

		fclose(input_file);

		printf("Closed input file.\n");
	}
	// --------------------- Input data test end ------------------- //


	// Allocate memory for output array
	float* output_data;

	printf("Here6!\n");

	float time_taken = 0;
	float bf_time = 0;
	int num_runs = 10;

	// Start timing FFT computation //
	struct timespec tval_before, tval_after;

	for(int ii = 0; ii < num_runs; ii++){
		// Start timing beamformer computation //
		clock_gettime(CLOCK_MONOTONIC, &tval_before);

		// Run FFT 
                // Things to keep in mind about FFT output:
		// - FFT shift possibly required after FFT if too much memory is allocated
		// - Output may need to be divided number of FFT points
                output_data = run_upchannelizer_beamformer(sim_data, sim_coefficients, n_pol, n_sim_ant, n_beam, n_chan, n_win, n_time_int, n_samp, telescope_flag);

		// Stop timing FFT computation //
		clock_gettime(CLOCK_MONOTONIC, &tval_after);
		time_taken = (float)(tval_after.tv_sec - tval_before.tv_sec); //*1e6; // Time in seconds since epoch
		time_taken = time_taken + (float)(tval_after.tv_nsec - tval_before.tv_nsec)*1e-9; // Time in nanoseconds since 'tv_sec - start and end'
		bf_time += time_taken;
	}
	printf("Average FFT and beamform processing time: %f s\n", bf_time/num_runs);

	printf("Here7, beamformer output: %f \n", output_data[0]);
	
	// Write data to binary file for analysis
	char output_filename[128];

	printf("Here8!\n");

	strcpy(output_filename, "/datag/users/mruzinda/o/output_d_fft_bf.bin");
	//strcpy(output_filename, "/mydatag/Unknown/GUPPI/output_d_fft_bf.bin");
	//strcpy(output_filename, "/home/mruzinda/tmp_output/output_d_fft_bf.bin");

	printf("Here9!\n");

	FILE* output_file;

	printf("Here10!\n");

	output_file = fopen(output_filename, "wb");

	printf("Here11!\n");

	fwrite(output_data, sizeof(float), n_beam*n_chan*n_samp*n_sti, output_file);

	printf("Here12!\n");

	fclose(output_file);

	printf("Closed output file.\n");

	//free(sim_data);
	//free(output_data);	

	printf("Freed output array and unregistered arrays in pinned memory.\n");

	// Free up device memory
	//cudaFreeHost(h_data);
	//cudaFreeHost(h_coeff);
	Cleanup_beamformer();

	printf("Here13!\n");

	return 0;
}
// <----Uncomment here if testing standalone code
