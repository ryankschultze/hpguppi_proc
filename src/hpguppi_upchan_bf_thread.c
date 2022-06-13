/* hpguppi_upchan_bf_thread_fb.c
 *
 * Reads HDF5 files containing phase solutions for phasing up, delays and rates for forming
 * beams, and other useful metadata.
 * Perform upchannelization, coherent beamforming and write databuf blocks out to filterbank files on disk.
 */

#define _GNU_SOURCE 1
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <assert.h>
#include <stdint.h>
#include <endian.h>
#include <math.h>

#include <unistd.h>
#include <string.h>
#include <pthread.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>
#include <sys/time.h>
#include <fcntl.h>
#include <errno.h>

#include "hdf5.h"
#include "ioprio.h"

//#include "hpguppi_time.h"
#include "upchannelizer_beamformer.h"

#include "hashpipe.h"

// Use rawpsec_fbutils because it contains all of the necessary functions to write headers to filterbank files
// This might change in the near future to make this library completely separate from rawspec
#include "rawspec_fbutils.h"
// Use rawpsec_rawutils because it contains all of the necessary functions to parse raw headers orginally from GUPPI RAW files
// This might change in the near future to make this library completely separate from rawspec
#include "rawspec_rawutils.h"

#include "hpguppi_databuf.h"
#include "hpguppi_params.h"
//#include "hpguppi_pksuwl.h"
#include "hpguppi_util.h"

// 80 character string for the BACKEND header record.
static const char BACKEND_RECORD[] =
// 0000000000111111111122222222223333333333
// 0123456789012345678901234567890123456789
  "BACKEND = 'GUPPI   '                    " \
  "                                        ";

static int safe_close_all(int *pfd) {
  if (pfd==NULL) return 0;
  int pfd_size = sizeof(pfd)/sizeof(int);
  int fd_val = -1;
  for(int i = 0; i < pfd_size; i++){
     if(pfd[i] == -1){
       fsync(pfd[i]);
       fd_val = close(pfd[i]);
       if(fd_val == -1){
         printf("A file was not successfully closed! \n");
       }
     }
  }
  return 0;
}

static void
update_fb_hdrs_from_raw_hdr_ubf(fb_hdr_t fb_hdr, const char *p_rawhdr)
{
  rawspec_raw_hdr_t raw_hdr;

  rawspec_raw_parse_header(p_rawhdr, &raw_hdr);
  hashpipe_info(__FUNCTION__,
      "beam_id = %d/%d", raw_hdr.beam_id, raw_hdr.nbeam);

  // Update filterbank headers based on raw params and Nts etc.
  // Same for all products
  fb_hdr.src_raj = raw_hdr.ra;
  fb_hdr.src_dej = raw_hdr.dec;
  fb_hdr.ibeam = raw_hdr.beam_id;
  fb_hdr.nbeams = raw_hdr.nbeam;
  strncpy(fb_hdr.source_name, raw_hdr.src_name, 80);
  fb_hdr.source_name[80] = '\0';
  //fb_hdr.nchans = 16; // 1k mode
  fb_hdr.nchans = 64; // 4k mode
  //fb_hdr.nchans = 512; // 32k mode
  // TODO az_start, za_start
}

static void *run(hashpipe_thread_args_t * args)
{
  // Local aliases to shorten access to args fields
  // Our output buffer happens to be a hpguppi_input_databuf
  hpguppi_input_databuf_t *db = (hpguppi_input_databuf_t *)args->ibuf;
  hashpipe_status_t *st = &args->st;
  const char * thread_name = args->thread_desc->name;
  const char * status_key = args->thread_desc->skey;

  // Read in general parameters
  struct hpguppi_params gp;
  struct psrfits pf;
  pf.sub.dat_freqs = NULL;
  pf.sub.dat_weights = NULL;
  pf.sub.dat_offsets = NULL;
  pf.sub.dat_scales = NULL;
  pthread_cleanup_push((void *)hpguppi_free_psrfits, &pf);

  // Init output file descriptor (-1 means no file open)
  static int fdraw[N_BEAM];
  memset(fdraw, -1, N_BEAM*sizeof(int));
  pthread_cleanup_push((void *)safe_close_all, fdraw);

  // Set I/O priority class for this thread to "real time" 
  if(ioprio_set(IOPRIO_WHO_PROCESS, 0, IOPRIO_PRIO_VALUE(IOPRIO_CLASS_RT, 7))) {
    hashpipe_error(thread_name, "ioprio_set IOPRIO_CLASS_RT");
  }

  // Loop
  int64_t pktidx=0, pktstart=0, pktstop=0;
  int blocksize=0; // Size of beamformer output (output block size)
  int curblock=0;
  int block_count=0; //filenum=0;
  int got_packet_0=0, first=1;
  char *ptr;
  int rv = 0;
  double obsbw;
  double tbin;
  double obsfreq;
  char fb_basefilename[200];
  char hdf5_basefilename[200];
  char coeff_basefilename[200];
  char coeff_fname[256];
  char outdir[200];
  char bfrdir[200];
  char src_name[200];
  uint64_t stt_imjd = 0;
  uint64_t stt_smjd = 0;
  char character = '_';
  char *char_offset;
  long int last_underscore_pos;
  char base_src[200];
  char base_no_src[200];
  char raw_basefilename[200];
  char raw_filename[200];
  char prev_basefilename[200];
  strcpy(prev_basefilename, "prev_basefilename"); // Initialize as different string that raw_basefilename
  char raw_obsid[200];
  
  // Filterbank file header initialization
  fb_hdr_t fb_hdr;
  fb_hdr.machine_id = 20;
  fb_hdr.telescope_id = -1; // Unknown (updated later)
  fb_hdr.data_type = 1;
  fb_hdr.nbeams =  1;
  fb_hdr.ibeam  =  -1; //Not used
  fb_hdr.nbits  = 32;
  fb_hdr.nifs   = 1;

  // Timing variables
  float bf_time = 0;
  float write_time = 0;
  float time_taken = 0;
  float time_taken_w = 0;
  
  int telescope_flag = 0;
  uint64_t synctime = 0;
  uint64_t hclocks = 1;
  uint32_t fenchan = 1;
  uint32_t schan = 0;
  double chan_bw = 1.0;
  uint64_t nants = 0;
  uint64_t obsnchan = 0;
  uint64_t piperblk = 0;
  uint64_t raw_blocsize = 0; // Raw file block size
  int n_fft = 0; // Number of time samples used for FFT (Number of FFT points)
  int n_samp = 0; // Number of time samples per block in a shared memory buffer block
  int n_samp_spec = 0; // Specification of the number of time samples
  int n_win = 0; // Number of spectra windows in output
  int n_win_spec = 8; // Specification of the number of spectra windows in output
  int n_time_int = 8; // Number of time samples to integrate
  int n_sti = 0; // Number of STI windows in output
  int subband_idx = 0; // Index of subband that ranges from 0 to 15
  int prev_subband_idx = -1; // Previous subband index
  double pktidx_time = 0;
  double block_midtime = 0;
  double time_array_midpoint = 0;

  hid_t file_id, npol_id, nbeams_id, obs_id, src_id, cal_all_id, delays_id, time_array_id, ra_id, dec_id, sid1, sid2, sid4, sid5, sid6, src_dspace_id, obs_type, native_obs_type, native_src_type; // identifiers //
  herr_t status, cal_all_elements, delays_elements, time_array_elements, ra_elements, dec_elements, src_elements;

  hid_t reim_tid;
  reim_tid = H5Tcreate(H5T_COMPOUND, sizeof(complex_t));
  H5Tinsert(reim_tid, "r", HOFFSET(complex_t, re), H5T_IEEE_F32LE);
  H5Tinsert(reim_tid, "i", HOFFSET(complex_t, im), H5T_IEEE_F32LE);

  complex_t *cal_all_data = NULL;
  double *delays_data;
  double *time_array_data;
  double *ra_data;
  double *dec_data;
  uint64_t nbeams;
  int actual_nbeams = 0;
  uint64_t npol;
  hvl_t *src_names_str;
  int hdf5_obsidsize;

  int time_array_idx = 0; // time_array index
  FILE* coeffptr;
  int coeff_flag = 0 ;

  float* bf_coefficients; // Beamformer coefficients
  float* tmp_coefficients; // Temporary coefficients

  // Frequency parameter initializations
  double coarse_chan_freq[N_FREQ]; // Coarse channel center frequencies in the band 
  int n_chan_per_node = 0; // Number of coarse channels per compute node
  int n_coarse_proc = 0; // Number of coarse channels processed at a time
  int half_n_coarse_proc = 0; // Half of the number of coarse channels processed at a time
  int n_subband = 16; // Number of subbands
  double center_freq = 0; // Center frequency of the subband and currently the filterbank file
  int n_samp_per_rawblk = 0; // Number of time samples in a RAW block used to calculate pktidx_time which is the approximate unix time at a given PKTIDX
  int n_ant_config = 0;

  int sim_flag = 0; // Flag to use simulated coefficients (set to 1) or calculated beamformer coefficients (set to 0)
  // Add if statement for generate_coefficients() function option which has 3 arguments - tau, coarse frequency channel, and epoch
  if(sim_flag == 1){
    int n_sim_ant = 58;
    if(n_sim_ant <= N_ANT/2){
      n_ant_config = N_ANT/2;
    }else{
      n_ant_config = N_ANT;
    }
    // Generate weights or coefficients (using simulate_coefficients() for now)
    // Used with simulated data when no RAW file is read
    //int n_chan = 16;  // 1k mode
    int n_chan = 64;  // 4k mode
    //int n_chan = 512; // 32k mode
    int n_pol = 2;
    int n_beam = 1;
    tmp_coefficients = simulate_coefficients_ubf(n_sim_ant, n_ant_config, n_pol, n_beam, n_chan, telescope_flag);
    // Register the array in pinned memory to speed HtoD mem copy
    //coeff_pin(tmp_coefficients);
  }
  if(sim_flag == 0){
    bf_coefficients = (float*)calloc(N_COEFF, sizeof(float)); // Beamformer coefficients
    //coeff_pin(bf_coefficients);
  }

  // -----------------Get phase solutions----------------- //

  // Make all initializations before while loop
  // Initialize beamformer (allocate all memory on the device)
  init_upchan_beamformer(telescope_flag);

  // Initialize output data array
  float* output_data;

  /*
  printf("UBF: Using host arrays allocated in pinned memory\n\r");
  for (int i = 0; i < N_INPUT_BLOCKS; i++){
    ///////////////////////////////////////////////////////////////////////////////
    //>>>>   Register host array in pinned memory <<<<
    ///////////////////////////////////////////////////////////////////////////////
    input_data_pin((signed char *)&db->block[i].data);
  }
  */

  char fname[256];    // Filterbank file name
  char datadir[1024]; // Output directory
  char *last_slash;   // Helps find last slash of path
  int rec_stop = 0; // Flag to free memory and show that processing stopped
  while (run_threads()) {

    // Note waiting status
    hashpipe_status_lock_safe(st);
    hputs(st->buf, status_key, "waiting");
    hashpipe_status_unlock_safe(st);

    // Wait for buf to have data
    rv = hpguppi_input_databuf_wait_filled(db, curblock);
    if(rv!=0)continue;
    
    // Get subband index
    hashpipe_status_lock_safe(st);
    hgeti4(st->buf, "SUBBAND", &subband_idx); // Get current index of subband being processed
    hashpipe_status_unlock_safe(st);

    // Read param struct for this block
    ptr = hpguppi_databuf_header(db, curblock);
    if (first) {
      hpguppi_read_obs_params(ptr, &gp, &pf);
      first = 0;
    } else {
      hpguppi_read_subint_params(ptr, &gp, &pf);
    }
    
    // Read pktidx, pktstart, pktstop from header
    hgeti8(ptr, "PKTIDX", &pktidx);
    hgeti8(ptr, "PKTSTART", &pktstart);
    hgeti8(ptr, "PKTSTOP", &pktstop);

    if((pktidx >= pktstop) && (subband_idx == (n_subband-1) || subband_idx == 0)){
      coeff_flag = 0;
      // Possibly free memory here so it can be reallocated at the beginning of a scan to compensate for a change in size
      if(cal_all_data != NULL){
        free(cal_all_data);
        cal_all_data = NULL;
        free(delays_data);
        delays_data = NULL;
        free(time_array_data);
        time_array_data = NULL;
        free(ra_data);
        ra_data = NULL;
        free(dec_data);
        dec_data = NULL;
        status = H5Dvlen_reclaim(native_src_type, src_dspace_id, H5P_DEFAULT, src_names_str);
        // Mark as free
        hpguppi_input_databuf_set_free(db, curblock);
        // Go to next block
        curblock = (curblock + 1) % db->header.n_block;
      }
      for(int b = 0; b < nbeams; b++){
        // If file open, close it
        if(fdraw[b] != -1) {
          // Close file
          close(fdraw[b]);
          // Reset fdraw, got_packet_0, filenum, block_count
          fdraw[b] = -1;
        }
      }
      got_packet_0 = 0;
      block_count=0;
      // Print end of recording conditions only once
      if(rec_stop == 0){
        rec_stop = 1;
        hashpipe_info(thread_name, "reached end of scan: "
        "pktstart %ld pktstop %ld pktidx %ld",
        pktstart, pktstop, pktidx);
        // Inform status buffer of processing status
        hashpipe_status_lock_safe(st);
        hgets(st->buf, "PROCSTAT", sizeof("END"), "END"); // Inform status buffer that the scan processing has ended
        hashpipe_status_unlock_safe(st);
      }
      continue;
    }
    // Reset rec_stop flag to notify the user of the end of a scan
    rec_stop = 0;

    // Inform status buffer of processing status
    hashpipe_status_lock_safe(st);
    hgets(st->buf, "PROCSTAT", sizeof("START"), "START"); // Inform status buffer that the scan processing has started
    hashpipe_status_unlock_safe(st);

    // Get values for calculations at varying points in processing
    hgetu8(ptr, "SYNCTIME", &synctime);
    hgetu8(ptr, "HCLOCKS", &hclocks);
    hgetu4(ptr, "FENCHAN", &fenchan);
    hgetr8(ptr, "CHAN_BW", &chan_bw); // In MHz
    hgetr8(ptr, "OBSFREQ", &obsfreq);
    hgetu8(ptr, "NANTS", &nants);
    hgetu8(ptr, "OBSNCHAN", &obsnchan);
    hgetu4(ptr, "SCHAN", &schan);
    hgetu8(ptr, "PIPERBLK", &piperblk);
    hgetu8(ptr, "STT_SMJD", &stt_smjd);
    hgetu8(ptr, "STT_IMJD", &stt_imjd);
    hgetr8(ptr,"TBIN", &tbin);
    hgets(ptr, "OBSID", sizeof(raw_obsid), raw_obsid);
    hgets(ptr, "SRC_NAME", sizeof(src_name), src_name);
    hgetr8(ptr, "OBSBW", &obsbw);
    hgetu8(ptr, "BLOCSIZE", &raw_blocsize); // Raw file block size
    
    hashpipe_status_lock_safe(st);
    hgets(st->buf, "BASEFILE", sizeof(raw_basefilename), raw_basefilename);
    hgets(st->buf, "OUTDIR", sizeof(outdir), outdir);
    hashpipe_status_unlock_safe(st);

    // If the subarray configuration is in use (half the number of antennas)
    if(nants <= N_ANT/2){
      n_ant_config = N_ANT/2;
    }else{
      n_ant_config = N_ANT;
    }

    // Get the appropriate basefile name from the stride_input_thread 
    // Get HDF5 file data at the beginning of the processing
    if(strcmp(prev_basefilename, raw_basefilename) != 0 || subband_idx != prev_subband_idx){
      strcpy(prev_basefilename, raw_basefilename);
      prev_subband_idx = subband_idx;

      printf("UBF: RAW file base filename:          %s \n", raw_basefilename);
      printf("UBF: Previous RAW file base filename: %s \n", prev_basefilename);

      // Calculate coarse channel center frequencies depending on the mode and center frequency that spans the RAW file
      n_chan_per_node = (int)obsnchan/nants; //((int)fenchan)/n_nodes;
      n_coarse_proc = n_chan_per_node/n_subband;

      // Skip zeroth index since the number of coarse channels is even and the center frequency is between the 2 middle channels
      for(int i=0; i<n_coarse_proc; i++){
        //coarse_chan_freq[i] = (i-((n_chan_per_node-1)/2))*coarse_chan_band + (obsfreq*1e6);
        coarse_chan_freq[i] = ((i + subband_idx*n_coarse_proc)-((n_chan_per_node-1)/2))*(chan_bw*1e-3) + (obsfreq*1e-3);
        // Equivalent equation //
        //coarse_chan_freq[i] = ((i + subband_idx*n_coarse_proc)-(n_chan_per_node/2))*(chan_bw*1e-3) + (obsfreq*1e-3) + ((chan_bw/2)*1e-3);
      }

      if((sim_flag == 0) && (block_count == 0)){
        hashpipe_status_lock_safe(st);
        hgets(st->buf, "BFRDIR", sizeof(bfrdir), bfrdir);
        hashpipe_status_unlock_safe(st);

        // Set specified path to read from HDF5 files
        strcpy(hdf5_basefilename, bfrdir);
        strcat(hdf5_basefilename, "/");
        strcat(hdf5_basefilename, raw_obsid);
        strcat(hdf5_basefilename, ".bfr5");

        // Read HDF5 file and get all necessary parameters (obsid, cal_all, delays, rates, time_array, ras, decs)
        // Open an existing file. //
        file_id = H5Fopen(hdf5_basefilename, H5F_ACC_RDONLY, H5P_DEFAULT);

        // -------------Read obsid first----------------- //
        // Open an existing dataset. //
        obs_id = H5Dopen(file_id, "/obsinfo/obsid", H5P_DEFAULT);
        // Get obsid data type //
        obs_type = H5Dget_type(obs_id);
        native_obs_type = H5Tget_native_type(obs_type, H5T_DIR_DEFAULT);
        hdf5_obsidsize = (int)H5Tget_size(native_obs_type);
        // Allocate memory to string array
        char hdf5_obsid[hdf5_obsidsize+1];
        hdf5_obsid[hdf5_obsidsize] = '\0'; // Terminate string
        // Read the dataset. // 
        status = H5Dread(obs_id, native_obs_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, hdf5_obsid);
        // Close the dataset. //
        status = H5Dclose(obs_id);
        // -----------------------------------------------//

        // -------------Read source names----------------- //
        // Open an existing dataset. //
        src_id = H5Dopen(file_id, "/beaminfo/src_names", H5P_DEFAULT);
        // Get dataspace ID //
        src_dspace_id = H5Dget_space(src_id);
        // Gets the number of elements in the data set //
        src_elements=H5Sget_simple_extent_npoints(src_dspace_id);
        //printf("Number of elements in the src_names dataset is : %d\n", src_elements);
        // Create src_names data type //
        native_src_type = H5Tvlen_create(H5T_NATIVE_CHAR);
        // Allocate memory to string array
        src_names_str = malloc((int)src_elements*sizeof(hvl_t));
        // Read the dataset. //
        status = H5Dread(src_id, native_src_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, src_names_str);
        for(int i=0; i<src_elements; i++) {
          printf("%d: len: %d, str is: %s\n", i, (int)src_names_str[i].len, (char *)src_names_str[i].p);
        }

        // Close the dataset //
        status = H5Dclose(src_id);
        // -----------------------------------------------//

        // Need to initialize these obsid variables
        // Figure out how to read them first
        if(raw_obsid == hdf5_obsid){
          printf("UBF: OBSID in RAW file and HDF5 file match!\n");
          printf("UBF: raw_obsid  = %s \n", raw_obsid);
          printf("UBF: hdf5_obsid = %s \n", hdf5_obsid);
        }else{
          printf("UBF (Warning): OBSID in RAW file and HDF5 file DO NOT match!\n");
          printf("UBF: raw_obsid  = %s \n", raw_obsid);
          printf("UBF: hdf5_obsid = %s \n", hdf5_obsid);
        }
        // Read cal_all once per HDF5 file. It doesn't change through out the entire recording.
        // Read delayinfo once, get all values, and update when required
        // Open an existing datasets //
        cal_all_id = H5Dopen(file_id, "/calinfo/cal_all", H5P_DEFAULT);
        delays_id = H5Dopen(file_id, "/delayinfo/delays", H5P_DEFAULT);
        time_array_id = H5Dopen(file_id, "/delayinfo/time_array", H5P_DEFAULT);
        ra_id = H5Dopen(file_id, "/beaminfo/ras", H5P_DEFAULT);
        dec_id = H5Dopen(file_id, "/beaminfo/decs", H5P_DEFAULT);
        npol_id = H5Dopen(file_id, "/diminfo/npol", H5P_DEFAULT);
        nbeams_id = H5Dopen(file_id, "/diminfo/nbeams", H5P_DEFAULT);

        // Get dataspace ID //
        sid1 = H5Dget_space(cal_all_id);
        sid2 = H5Dget_space(delays_id);
        sid4 = H5Dget_space(time_array_id);
        sid5 = H5Dget_space(ra_id);
        sid6 = H5Dget_space(dec_id);
  
        // Gets the number of elements in the data set //
        cal_all_elements=H5Sget_simple_extent_npoints(sid1);
        delays_elements=H5Sget_simple_extent_npoints(sid2);
        time_array_elements=H5Sget_simple_extent_npoints(sid4);
        ra_elements=H5Sget_simple_extent_npoints(sid5);
        dec_elements=H5Sget_simple_extent_npoints(sid6);

        // Allocate memory for array
        cal_all_data = malloc((int)cal_all_elements*sizeof(complex_t));
        delays_data = malloc((int)delays_elements*sizeof(double));
        time_array_data = malloc((int)time_array_elements*sizeof(double));
        ra_data = malloc((int)ra_elements*sizeof(double));
        dec_data = malloc((int)dec_elements*sizeof(double));

        // Read the dataset. //
        status = H5Dread(npol_id, H5T_STD_I64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &npol);

        status = H5Dread(nbeams_id, H5T_STD_I64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nbeams);

        status = H5Dread(cal_all_id, reim_tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, cal_all_data);

        status = H5Dread(delays_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, delays_data);

        status = H5Dread(time_array_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, time_array_data);

        status = H5Dread(ra_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ra_data);

        status = H5Dread(dec_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dec_data);

        // Close the dataset. //
        status = H5Dclose(cal_all_id);
        status = H5Dclose(delays_id);
        status = H5Dclose(time_array_id);
        status = H5Dclose(ra_id);
        status = H5Dclose(dec_id);
        status = H5Dclose(npol_id);
        status = H5Dclose(nbeams_id);

        // Close the file. //
        status = H5Fclose(file_id);
        printf("status of file close: %d\n", status);

        // The actual number of beams set to compensate for the beams being more than 64 in the BFR5 file
        actual_nbeams = (int)nbeams;
        // Set the number of beams to 64 (number set in upchannelizer_beamformer.h) if they are greater than 64 in the BFR5 file
        if(nbeams > N_BEAM){
          nbeams = N_BEAM;
        }
        // In the subarray configuration the number of beams needs to be half due to mempry constraints with the 2080 Tis
        if(nants <= N_ANT/2){
          nbeams = N_BEAM/2;
        }

        // Reset indexing of time array since new HDF5 file is in use
        time_array_idx = 0;

        // Assign values to tmp variable then copy values from it to pinned memory pointer (bf_coefficients)
        tmp_coefficients = generate_coefficients_ubf(cal_all_data, delays_data, time_array_idx, coarse_chan_freq, n_ant_config, (int)npol, (int)nbeams, actual_nbeams,(int)schan, n_coarse_proc, subband_idx, nants, telescope_flag);
        memcpy(bf_coefficients, tmp_coefficients, N_COEFF*sizeof(float));

        // Get basefilename with no source name
        // strrchr() finds the last occurence of the specified character
        char_offset = strrchr(raw_basefilename, character);
        last_underscore_pos = char_offset-raw_basefilename;

        // Get file name with source name
        memcpy(base_src, &raw_basefilename[0], last_underscore_pos);
        base_src[last_underscore_pos] = '\0';

        // Get basefilename with source name
        // strrchr() finds the last occurence of the specified character
        char_offset = strrchr(base_src, character);
        last_underscore_pos = char_offset-base_src;

        // Get file name with no source name
        memcpy(base_no_src, &base_src[0], last_underscore_pos+1);
        base_no_src[last_underscore_pos+1] = '\0';

        // Number of FFT points depending on mode
        if(n_coarse_proc == 1){        // 1k mode
          n_fft = 524288; // 2^19 point FFT
        }else if(n_coarse_proc == 4){  // 4k mode
          n_fft = 131072; // 2^17 point FFT
        }else if(n_coarse_proc == 32){ // 32k mode
          n_fft = 16384;  // 2^14 point FFT
        }
        n_samp_spec = n_fft*n_win_spec;
      }

      if(sim_flag == 1){
        nbeams = 1; // Generate one filterbank file with the boresight beam
        npol = 2; // Number of polarizations needs to be specified
        // Set specified path to write filterbank files
        strcpy(fb_basefilename, outdir);
        strcat(fb_basefilename, "/");
        strcat(fb_basefilename, raw_basefilename);
      }

      // If this is a new scan, write headers and so on to new filterbank files
      got_packet_0 = 0;
    }

    // Number of time samples in a RAW file
    hashpipe_status_lock_safe(st);
    hgeti4(st->buf, "NSAMP", &n_samp); 
    hashpipe_status_unlock_safe(st);

    printf("UBF: n_samp = %d\n", n_samp);

    // Zero padded if the number of time samples is less than the specification
    // If the number of time samples is greater by 2, then the number of antennas is smaller by half (subarray configuration)
    // So there is a situation where one of the RAW files may have less time samples than 2 times the 
    // specified value here in the subarray configuration. Try to compensate for that with the else if()
    if(n_samp < n_samp_spec){
      n_samp = n_samp_spec;
    }else if(n_samp > (3*n_samp_spec/2) && n_samp < 2*n_samp_spec){
      n_samp = 2*n_samp_spec;
    }

    // Number of time samples used for FFT (Number of FFT points)
    n_time_int = n_samp/n_fft;

    // Number of spectra windows
    n_win = n_time_int;

    // Number of STI windows 
    n_sti = n_win/n_time_int;

    // Size of beamformer output
    blocksize= n_coarse_proc*n_fft*n_sti*sizeof(float); 

    if(sim_flag == 0){
      // Write coefficients to binary file for analysis with CASA
      if((time_array_idx == 0)  && (coeff_flag == 0)){
        strcpy(coeff_basefilename, outdir);
        strcat(coeff_basefilename, "/");
        strcat(coeff_basefilename, raw_basefilename);
        printf("UBF: Beamformer coefficient file name with path and no extension yet: %s \n", coeff_basefilename);

        sprintf(coeff_fname, "%s.coeff.start.bin",  coeff_basefilename);

        coeffptr = fopen(coeff_fname,"wb");

        fwrite(bf_coefficients, sizeof(float), N_COEFF, coeffptr);
         
        fclose(coeffptr);
      } 

      // Number of time samples in a RAW block used to calculate pktidx_time which is the approximate unix time at a given PKTIDX
      n_samp_per_rawblk = (int)(raw_blocsize/(2*obsnchan*npol));

      // Compute unix time in the middle of the block and avg. value current and next time_array value for comparison
      // Update coefficient if unix time is in the middle of the shared memory buffer block (not the RAW block so use 'n_samp' 
      // rather than 'n_samp_per_rawblk') is greater than avg. of the two time array values
      pktidx_time = synctime + (pktidx*tbin*n_samp_per_rawblk/piperblk);
      block_midtime = pktidx_time + tbin*n_samp/2;
      if(time_array_idx < (time_array_elements-1)){
        time_array_midpoint = (time_array_data[time_array_idx] + time_array_data[time_array_idx+1])/2;
      }else if(time_array_idx >= (time_array_elements-1)){
        time_array_midpoint = time_array_data[time_array_idx];
        time_array_idx = time_array_elements-2;
      }

      // Update coefficients every specified number of blocks
      if(block_midtime >= time_array_midpoint){
        // Then update with delays[n+1]
        time_array_idx += 1;

        // Update coefficients with new delay
        // Assign values to tmp variable then copy values from it to pinned memory pointer (bf_coefficients)
        tmp_coefficients = generate_coefficients_ubf(cal_all_data, delays_data, time_array_idx, coarse_chan_freq, n_ant_config, (int)npol, (int)nbeams, actual_nbeams, (int)schan, n_coarse_proc, subband_idx, nants, telescope_flag);

        memcpy(bf_coefficients, tmp_coefficients, N_COEFF*sizeof(float));

        // Write coefficients to binary file for analysis with CASA
        if((time_array_idx == 15)  && (coeff_flag == 0)){
          coeff_flag = 1; // Only necessary to do once in scan

          sprintf(coeff_fname, "%s.coeff.middle.bin",  coeff_basefilename);

          coeffptr = fopen(coeff_fname,"wb");

          fwrite(bf_coefficients, sizeof(float), N_COEFF, coeffptr);
          
          fclose(coeffptr);
        } 
       
      }
    }

    /* Set up data ptr for quant routines */
    pf.sub.data = (unsigned char *)hpguppi_databuf_data(db, curblock);

    // Wait for packet 0 before starting write
    // "packet 0" is the first packet/block of the new recording,
    // it is not necessarily pktidx == 0.
    if ((got_packet_0==0 && gp.stt_valid==1)) {
      got_packet_0 = 1;

      // Create the output directory if needed
      strncpy(datadir, fb_basefilename, 1023);
      last_slash = strrchr(datadir, '/');
      if (last_slash!=NULL && last_slash!=datadir) {
	*last_slash = '\0';
	hashpipe_info(thread_name,
	  "Using directory '%s' for output", datadir);
        if(mkdir_p(datadir, 0755) == -1) {
	  hashpipe_error(thread_name, "mkdir_p(%s)", datadir);
          pthread_exit(NULL);
        }
      }

      // Update filterbank headers based on raw params, Nts, and BFR5 params etc.
      // Technically unnecessary for now, but I might move all of the filterbank header info below to this function so leaving it here for now
      update_fb_hdrs_from_raw_hdr_ubf(fb_hdr, ptr);

      // Open nbeams filterbank files to save a beam per file i.e. N_BIN*n_fft*sizeof(float) per file.
      //printf("UBF: Opening filterbank files \n");
      for(int b = 0; b < nbeams; b++){
        if(b >= 0 && b < 10) {
          if(sim_flag == 0){
            // Set specified path to write filterbank files
            strcpy(fb_basefilename, outdir);
            strcat(fb_basefilename, "/");
            strcat(fb_basefilename, base_no_src);
            strcat(fb_basefilename, (char *)src_names_str[b].p);

            if(subband_idx >= 0 && subband_idx < 10) {
              sprintf(fname, "%s.SB0%d.B0%d.fil",  fb_basefilename, subband_idx, b);
            }else{
              sprintf(fname, "%s.SB%d.B0%d.fil",  fb_basefilename, subband_idx, b);
            }
          }else if(sim_flag == 1){
            sprintf(fname, "%s.B0.fil",  fb_basefilename);
          }
        }else{
          // Set specified path to write filterbank files
          strcpy(fb_basefilename, outdir);
          strcat(fb_basefilename, "/");
          strcat(fb_basefilename, base_no_src);
          strcat(fb_basefilename, (char *)src_names_str[b].p);

          if(subband_idx >= 0 && subband_idx < 10) {
            sprintf(fname, "%s.SB0%d.B%d.fil",  fb_basefilename, subband_idx, b);
          }else{
            sprintf(fname, "%s.SB%d.B%d.fil",  fb_basefilename, subband_idx, b);
          }
        }
        hashpipe_info(thread_name, "Opening fil file '%s'", fname);

                
        fdraw[b] = open(fname, O_CREAT|O_WRONLY|O_TRUNC, 0644);
        //fdraw[b] = open(fname, O_CREAT|O_WRONLY|O_APPEND, 0644);
        if(fdraw[b] == -1) {
	  // If we can't open this output file, we probably won't be able to
          // open any more output files, so print message and bail out.
          hashpipe_error(thread_name,
	    "cannot open filterbank output file, giving up");
            pthread_exit(NULL);
        }
    	posix_fadvise(fdraw[b], 0, 0, POSIX_FADV_DONTNEED);
      }

      // Get center frequency depending on mode (1k, 4k or 32k)
      if(n_coarse_proc == 1){
        center_freq = coarse_chan_freq[0]*1e3;
      }else{
        half_n_coarse_proc = (n_coarse_proc/2);
        center_freq = ((coarse_chan_freq[half_n_coarse_proc-1]+coarse_chan_freq[half_n_coarse_proc])/2)*1e3;
      }

      // Filterbank header values for now
      fb_hdr.telescope_id = 64; // MeerKAT ID (Don't know why it's 64, maybe associated with the number of antennas)
      fb_hdr.foff = chan_bw/n_fft; // Filterbank channel bandwidth in MHz
      fb_hdr.nchans = n_coarse_proc*n_fft; // Number of channels in a filterbank file
      fb_hdr.fch1 = center_freq; // Center frequency in MHz
      fb_hdr.tsamp = tbin*n_fft*n_time_int; // Time interval between output samples
      fb_hdr.tstart = stt_imjd + stt_smjd/86400.0; // tstart for now
      // Write RAW filename to filterbank header
      sprintf(raw_filename, "%s.0000.raw",  raw_basefilename);
      strncpy(fb_hdr.rawdatafile, raw_filename, 80);
      fb_hdr.rawdatafile[80] = '\0';

      // Write filterbank header to output file
      for(int b = 0; b < nbeams; b++){
        fb_hdr.ibeam =  b;
        if(sim_flag == 0){
          fb_hdr.src_raj = ra_data[b];
          fb_hdr.src_dej = dec_data[b];
          strncpy(fb_hdr.source_name, (char *)src_names_str[b].p, 80);
          fb_hdr.source_name[80] = '\0';
        }
        fb_fd_write_header(fdraw[b], &fb_hdr);
      }
    }

    // If we got packet 0, process and write data to disk
    if (got_packet_0) {

      // Note writing status
      hashpipe_status_lock_safe(st);
      hputs(st->buf, status_key, "writing");
      hashpipe_status_unlock_safe(st);

      // Start timing beamforming computation
      struct timespec tval_before, tval_after;
      clock_gettime(CLOCK_MONOTONIC, &tval_before);

      if(sim_flag == 0){
        output_data = run_upchannelizer_beamformer((signed char *)&db->block[curblock].data, bf_coefficients, (int)npol, (int)nants, (int)nbeams, (int)n_coarse_proc, n_win, n_time_int, n_fft, telescope_flag);
      }else if(sim_flag == 1){
        output_data = run_upchannelizer_beamformer((signed char *)&db->block[curblock].data, tmp_coefficients, (int)npol, (int)nants, (int)nbeams, (int)n_coarse_proc, n_win, n_time_int, n_fft, telescope_flag);
      }

      // Set beamformer output (CUDA kernel before conversion to power), that is summing, to zero before moving on to next block
      //set_to_zero_ubf();

      // Stop timing beamforming computation
      clock_gettime(CLOCK_MONOTONIC, &tval_after);
      time_taken = (float)(tval_after.tv_sec - tval_before.tv_sec); //*1e6; // Time in seconds since epoch
      time_taken = time_taken + (float)(tval_after.tv_nsec - tval_before.tv_nsec)*1e-9; // Time in nanoseconds since 'tv_sec - start and end'
      bf_time = time_taken;

      printf("UBF: run_beamformer() plus set_to_zero() time: %f s\n", bf_time);

      // Start timing write
      struct timespec tval_before_w, tval_after_w;
      clock_gettime(CLOCK_MONOTONIC, &tval_before_w);

      for(int b = 0; b < nbeams; b++){
        rv = write(fdraw[b], &output_data[b*n_coarse_proc*n_fft*n_sti], (size_t)blocksize);
        if(rv != blocksize){
          char msg[100];
          perror(thread_name);
          sprintf(msg, "Error writing data (output_data=%p, blocksize=%d, rv=%d)", output_data, blocksize, rv);
          hashpipe_error(thread_name, msg);
        }

        // flush output
        fsync(fdraw[b]);
      }
      // Stop timing write
      clock_gettime(CLOCK_MONOTONIC, &tval_after_w);
      time_taken_w = (float)(tval_after_w.tv_sec - tval_before_w.tv_sec); //*1e6; // Time in seconds since epoch
      time_taken_w = time_taken_w + (float)(tval_after_w.tv_nsec - tval_before_w.tv_nsec)*1e-9; // Time in nanoseconds since 'tv_sec - start and end'
      write_time = time_taken_w;
      printf("UBF: Time taken to write block of size, %d bytes, to disk = %f s \n", (int)nbeams*blocksize, write_time);
      
      printf("UBF: After write() function! Block index = %d \n", block_count);

      // Increment counter
      block_count++;
    }

    // Mark as free
    hpguppi_input_databuf_set_free(db, curblock);

    // Go to next block
    curblock = (curblock + 1) % db->header.n_block;

    // Check for cancel
    pthread_testcancel();

  }

  pthread_cleanup_pop(0); // Closes safe_close 

  pthread_cleanup_pop(0); // Closes hpguppi_free_psrfits

  // Free up device memory
  Cleanup_beamformer();

  hashpipe_info(thread_name, "exiting!");
  pthread_exit(NULL);
}

static hashpipe_thread_desc_t hpguppi_upchan_bf_thread = {
  name: "hpguppi_upchan_bf_thread",
  skey: "DISKSTAT",
  init: NULL,
  run:  run,
  ibuf_desc: {hpguppi_input_databuf_create},
  obuf_desc: {NULL}
};

static __attribute__((constructor)) void ctor()
{
  register_hashpipe_thread(&hpguppi_upchan_bf_thread);
}

// vi: set ts=8 sw=2 noet :
