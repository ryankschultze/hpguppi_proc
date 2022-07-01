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

// Testing to see whether input data is as I expect it to be
float* data_test(signed char *sim_data);
// Print out how to use standalone mode to the user
void printUsage();

int main() {
	if(argc > 2) {
		printUsage();
		return -1;
	}

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
