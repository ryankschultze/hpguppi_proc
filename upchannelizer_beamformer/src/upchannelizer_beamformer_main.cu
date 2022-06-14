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

// Testing to see whether input data is as I expect it to be
float* data_test(signed char *sim_data){
	float* data_float;
	data_float = (float*)calloc(N_INPUT, sizeof(float));
	for (int ii = 0; ii < N_INPUT; ii++) { // Write up to the size of the data corresponding to 1k, 4k or 32k mode
		data_float[ii] = sim_data[ii]*1.0f;
	}
	return data_float;
}

int main(int argc, char **argv) {
	if((argc > 6) || (argc < 2)) {
		printf("The minimum requirement is to set the output file as the first argument. Enter -h option for help. \n");
		return -1;
	}

	if(strcmp(argv[1], "-h") == 0){
		printf("To execute this program enter the following command:\n");
		printf("    ./upchannelizer_beamformer_main.exe <output file> <simulated data flag> <simulated coefficients flag> <telescope flag> <mode flag or VLASS specifications depending on telescope flag>\n");
		printf("    <> are not used in the command, but are just used to indicate arguments in this description\n");
		printf("Descriptions of the arguments in the command above:\n");
		printf("    <output file> - Enter the binary file along with it's path e.g. /datag/users/mruzinda/o/output_d_fft_bf.bin \n");
		printf("    The minimum requirement is to set the output file as the first argument. \n");
		printf("    <simulated data flag> - Enter the flag for the kind of simulated data that you would like to use. Default is 0. The following are the options:\n");
		printf("        sim_data = 0 -> Ones (Default)\n");
		printf("        sim_data = 1 -> Ones placed in a particular bin (bin 3 for now)\n");
		printf("        sim_data = 2 -> Ones placed in a particular bin at a particular antenna (bin 3 and antenna 3 for now)\n");
		printf("        sim_data = 3 -> Rect placed in a particular bin at a particular antenna (bin 3 and antenna 3 for now)\n");
		printf("        sim_data = 4 -> Simulated cosine wave\n");
		printf("        sim_data = 5 -> Simulated complex exponential i.e. exp(j*2*pi*f0*t)\n");
		printf("    <simulated coefficients flag> - Enter the flag for the kind of simulated coefficients that you would like to use. Default is 0. The following are the options:\n");
		printf("        sim_coef = 0 -> Ones (Default)\n");
		printf("        sim_coef = 1 -> Scale each beam by incrementing value i.e. beam 1 = 1, beam 2 = 2, ..., beam 64 = 64\n");
		printf("        sim_coef = 2 -> Scale each beam by incrementing value in a particular bin (bin 3 and 6 for now). Match simulated data sim_flag = 2\n");
		printf("        sim_coef = 3 -> Simulated beams from 58 to 122 degrees. Assuming a ULA.\n");
		printf("        sim_coef = 4 -> One value at one polarization, one element, one beam, and one frequency bin\n");
		printf("    <telescope flag> - Indicate the observatory specifications that you would like to use:\n");
		printf("        MK  -> MeeKAT specifications \n");
		printf("        VLA -> VLA specifications \n");
		printf("    <mode flag> - If MK is selected, then the next argument is the mode, which one of the 3 should be entered as 1k, 4k or 32k. \n");
		printf("    <VLASS specification> - Then if VLA is specified, indicate whether the specifications are the required or desired ones. The required are the default\n");
		printf("    If VLA is specified, the next argument should be input as:\n");
		printf("        required -> Required specifications \n");
		printf("        desired  -> Desired specifications \n");

		return 0;
	}

	int sim_data_flag = 0;
	int sim_coef_flag = 0;

	// If only the output file is entered
	if(argc == 2){
		sim_data_flag = 5;
		sim_coef_flag = 4;
	}// If the output file and another argument are entered
	else if(argc == 3){
		sim_data_flag = atoi(argv[2]);
		if(sim_data_flag < 0 || sim_data_flag > 5){
			printf("sim_data_flag is out of bounds i.e. this option doesn't exist. The flag has been set to 5, the default. \n");
			sim_data_flag = 5;
		}
		sim_coef_flag = 4;
	}// If the output file and simulated flags are entered
	else if(argc == 4){
		sim_data_flag = atoi(argv[2]);
		sim_coef_flag = atoi(argv[3]);
		if(sim_data_flag < 0 || sim_data_flag > 5){
			printf("sim_data_flag is out of bounds i.e. this option doesn't exist. The flag has been set to 5, the default. \n");
			sim_data_flag = 5;
		}
		if(sim_coef_flag < 0 || sim_coef_flag > 4){
			printf("sim_coef_flag is out of bounds i.e. this option doesn't exist. The flag has been set to 0, the default. \n");
			sim_coef_flag = 4;
		}
	}

	int telescope_flag = 0;
	// Check for telescope flag
	if(argc > 4){
		if(strcmp(argv[4], "MK") == 0){
			telescope_flag = 0;
		}else if(strcmp(argv[4], "VLA") == 0){
			telescope_flag = 1;
		}
	}// Default telescope
	else{
		printf("The observatory was not entered. The default is MK -> MeerKAT.\n");
		printf("Enter -h as argument for help.\n");
		telescope_flag = 0;
	}

	char mode_flag[5]; // Flag for operational mode for MeerKAT
	int spec_flag = 0; // Specification flag for VLASS
	// If MK is chosen, also get the mode
	if((argc > 5) && (strcmp(argv[4], "MK") == 0)){
		strcpy(mode_flag, argv[5]);
	}// If VLA is chosen, specify required or desired spcs
	else if((argc > 5) && (strcmp(argv[4], "VLA") == 0)){
		if(strcmp(argv[5], "required") == 0){
			spec_flag = 0;
		}else if(strcmp(argv[5], "desired") == 0){
			spec_flag = 1;
		}else{
			printf("Incorrect option, enter <required> or <desired>. The default is with <required> specifications for VLASS.\n");
			printf("Enter -h as argument for help.\n");
			spec_flag = 0;
		}
	}// Default mode
	else{
		printf("4k mode is the default with MeerKAT chosen or if the telescope is not specified. \n");
		strcpy(mode_flag, "4k");
	}


	printf("Mode = %s \n", mode_flag);
	char output_filename[128];

	strcpy(output_filename, argv[1]);
	printf("Output filename = %s\n", output_filename);
	//strcpy(output_filename, "/datag/users/mruzinda/o/output_d_fft_bf.bin");
	//strcpy(output_filename, "/mydatag/Unknown/GUPPI/output_d_fft_bf.bin");
	//strcpy(output_filename, "/home/mruzinda/tmp_output/output_d_fft_bf.bin");

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
			if(strcmp(mode_flag, "1k") == 0){
				n_chan = 1; 
		        	nt = 2*4096*1024; // 4194304; // 2^22
			}// 4k mode
			else if(strcmp(mode_flag, "4k") == 0){
			    	n_chan = 4; // 64
			        nt = 2*1024*1024; // 1048576; // 2^20
			}// 32k mode
			else if(strcmp(mode_flag, "32k") == 0){
				n_chan = 32;
				nt = 2*128*1024; // 131072; // 2^17
			}

			n_win = 16;
			n_time_int = 16;
		}else{ // Regular array configuration
			n_ant_config = N_ANT;
			// 5 seconds worth of processing at a time
			// 1k mode
			if(strcmp(mode_flag, "1k") == 0){
				n_chan = 1; 
			        nt = 4096*1024; // 4194304; // 2^22
			}// 4k mode
			else if(strcmp(mode_flag, "4k") == 0){
				n_chan = 4; // 64
		        	nt = 1024*1024; // 1048576; // 2^20
			}// 32k mode
			if(strcmp(mode_flag, "32k") == 0){
				n_chan = 32;
				nt = 128*1024; // 131072; // 2^17
			}
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

	printf("Specifications are as follows:\n");
	printf("n_beam = %d\n", n_beam);
	printf("n_pol  = %d\n", n_pol);
	printf("n_ant  = %d\n", n_sim_ant);
	printf("n_freq = %d (Number of coarse channels)\n", n_chan);
	printf("n_time = %d (Number of time samples)\n", nt);
	printf("n_fft  = %d (Number of points in FFT)\n", n_samp);
	printf("n_win  = %d (Number of spectral windowns after FFT)\n", n_win);
	printf("n_int  = %d (Number of integrated windows)\n", n_time_int);

	// Allocate memory to all arrays used by run_FFT() 
	init_upchan_beamformer(telescope_flag);
	printf("Allocated memory \n");

	// Generate simulated data
	signed char* sim_data = simulate_data_ubf(n_sim_ant, n_ant_config, n_pol, n_chan, n_samp, n_win, sim_data_flag, telescope_flag);
        printf("Simulated data \n");

	// Generate simulated weights or coefficients
	float* sim_coefficients = simulate_coefficients_ubf(n_sim_ant, n_ant_config, n_pol, n_beam, n_chan, sim_coef_flag, telescope_flag);
	printf("Simulated coefficients \n");

	// --------------------- Input data test --------------------- //
	int input_write = 0; // If input_write is set to 1, the simulated data will be written to a binary file for testing/verification

	if(input_write == 1){
		float* input_test = data_test(sim_data);

		// Write data to binary file for analysis
		char input_filename[128];

		strcpy(input_filename, "/mydatag/Unknown/GUPPI/input_h_fft_bf.bin");

		FILE* input_file;

		input_file = fopen(input_filename, "w");

		fwrite(input_test, sizeof(float), n_input, input_file);

		fclose(input_file);
	}
	// --------------------- Input data test end ------------------- //


	// Allocate memory for output array
	float* output_data;

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

	printf("First element of beamformer output: %f \n", output_data[0]);
	
	// Write data to binary file for analysis
	FILE* output_file;

	output_file = fopen(output_filename, "wb");

	fwrite(output_data, sizeof(float), n_beam*n_chan*n_samp*n_sti, output_file);

	fclose(output_file);

	printf("Closed output file.\n");

	//free(sim_data);
	//free(output_data);	

	// Free up device memory
	//cudaFreeHost(h_data);
	//cudaFreeHost(h_coeff);
	Cleanup_beamformer();

	return 0;
}
