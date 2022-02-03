/* hpguppi_coherent_bf_thread_fb.c
 *
 * Reads HDF5 files containing phase solutions for phasing up, delays and rates for forming
 * beams, and other useful metadata.
 * Perform coherent beamforming and write databuf blocks out to filterbank files on disk.
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
#include "coherent_beamformer_char_in.h"

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

void
update_fb_hdrs_from_raw_hdr_cbf(fb_hdr_t fb_hdr, const char *p_rawhdr)
{
  rawspec_raw_hdr_t raw_hdr;

  rawspec_raw_parse_header(p_rawhdr, &raw_hdr);
  hashpipe_info(__FUNCTION__,
      "beam_id = %d/%d", raw_hdr.beam_id, raw_hdr.nbeam);

  // Update filterbank headers based on raw params and Nts etc.
  // Same for all products
  fb_hdr.telescope_id = fb_telescope_id(raw_hdr.telescop);
  fb_hdr.src_raj = raw_hdr.ra;
  fb_hdr.src_dej = raw_hdr.dec;
  fb_hdr.tstart = raw_hdr.mjd;
  fb_hdr.ibeam = raw_hdr.beam_id;
  fb_hdr.nbeams = raw_hdr.nbeam;
  strncpy(fb_hdr.source_name, raw_hdr.src_name, 80);
  fb_hdr.source_name[80] = '\0';
  // Output product dependent
  fb_hdr.foff = raw_hdr.obsbw/raw_hdr.obsnchan/N_TIME;
  // This computes correct fch1 for odd or even number of fine channels
  fb_hdr.fch1 = raw_hdr.obsfreq
    - raw_hdr.obsbw*(raw_hdr.obsnchan-1)/(2*raw_hdr.obsnchan)
    - (N_TIME/2) * fb_hdr.foff
    ;//TODO + schan * raw_hdr.obsbw / raw_hdr.obsnchan; // Adjust for schan
  //fb_hdr.nchans = 16; // 1k mode
  fb_hdr.nchans = 64; // 4k mode
  //fb_hdr.nchans = 512; // 32k mode
  fb_hdr.tsamp = raw_hdr.tbin * N_TIME * 256*1024;
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

  /* Read in general parameters */
  struct hpguppi_params gp;
  struct psrfits pf;
  pf.sub.dat_freqs = NULL;
  pf.sub.dat_weights = NULL;
  pf.sub.dat_offsets = NULL;
  pf.sub.dat_scales = NULL;
  pthread_cleanup_push((void *)hpguppi_free_psrfits, &pf);

  /* Init output file descriptor (-1 means no file open) */
  static int fdraw[N_BEAM];
  memset(fdraw, -1, N_BEAM*sizeof(int));
  pthread_cleanup_push((void *)safe_close_all, fdraw);

  // Set I/O priority class for this thread to "real time" 
  if(ioprio_set(IOPRIO_WHO_PROCESS, 0, IOPRIO_PRIO_VALUE(IOPRIO_CLASS_RT, 7))) {
    hashpipe_error(thread_name, "ioprio_set IOPRIO_CLASS_RT");
  }

  printf("CBF: After ioprio_set...\n");

  /* Loop */
  int64_t pktidx=0, pktstart=0, pktstop=0;
  int blocksize=0; // Size of beamformer output (output block size)
  int curblock=0;
  int block_count=0, filenum=0;
  int got_packet_0=0, first=1;
  char *ptr;
  int rv = 0;
  double obsbw;
  double tbin;
  double obsfreq;
  char fb_basefilename[200];
  char hdf5_basefilename[200];
  char outdir[200];
  char bfrdir[200];
  //char character = '/';
  //char *char_offset; 
  //long int slash_pos;
  char raw_basefilename[200];
  char prev_basefilename[200];
  strcpy(prev_basefilename, "prev_basefilename"); // Initialize as different string that raw_basefilename
  //char new_base[200];
  char raw_obsid[200];
  //int header_size = 0;

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
  
  uint64_t synctime = 0;
  uint64_t hclocks = 1;
  uint32_t fenchan = 1;
  double chan_bw = 1.0;
  uint64_t nants = 0;
  uint64_t obsnchan = 0;
  uint64_t piperblk = 0;
  uint64_t raw_blocsize = 0; // Raw file block size
  int n_samp = 0; // Number of time samples per block in a RAW file
  int n_win = 0; // Number of STI windows in output
  double pktidx_time = 0;
  double block_midtime = 0;
  double time_array_midpoint = 0;

  hid_t file_id, npol_id, nbeams_id, obs_id, cal_all_id, delays_id, time_array_id, sid1, sid2, sid4, obs_type, native_obs_type; // identifiers //
  herr_t status, cal_all_elements, delays_elements, time_array_elements;

  //typedef struct complex_t{
  //  float re;
  //  float im;
  //}complex_t;

  hid_t reim_tid;
  reim_tid = H5Tcreate(H5T_COMPOUND, sizeof(complex_t));
  H5Tinsert(reim_tid, "r", HOFFSET(complex_t, re), H5T_IEEE_F32LE);
  H5Tinsert(reim_tid, "i", HOFFSET(complex_t, im), H5T_IEEE_F32LE);

  complex_t *cal_all_data;
  double *delays_data;
  double *time_array_data;
  uint64_t nbeams;
  uint64_t npol;
  int hdf5_obsidsize;

  int time_array_idx = 0; // time_array index

  int a = 34; // Antenna index
  int b = 1;  // Beam index
  int t = 1;  // Time stamp index
  int p = 1;  // Polarization index
  int c = 223;// Coarse channel index

  float* bf_coefficients; // Beamformer coefficients
  float* tmp_coefficients; // Temporary coefficients

  // Frequency parameter initializations
  float full_bw = 856e6; // Hz
  float coarse_chan_band = 0; // Coarse channel bandwidth
  double coarse_chan_freq[N_FREQ]; // Coarse channel center frequencies in the band 
  int n_nodes = 64; // Number of compute nodes
  int n_chan_per_node = 0; // Number of coarse channels per compute node

  int sim_flag = 1; // Flag to use simulated coefficients (set to 1) or calculated beamformer coefficients (set to 0)
  // Add if statement for generate_coefficients() function option which has 3 arguments - tau, coarse frequency channel, and epoch
  if(sim_flag == 1){
    // Generate weights or coefficients (using simulate_coefficients() for now)
    // Used with simulated data when no RAW file is read
    //int n_chan = 16;  // 1k mode
    //int n_chan = 64;  // 4k mode
    int n_chan = 512; // 32k mode
    int n_pol = 2;
    int n_beam = 1;
    tmp_coefficients = simulate_coefficients(n_pol, n_beam, n_chan);
    // Register the array in pinned memory to speed HtoD mem copy
    coeff_pin(tmp_coefficients);
  }
  if(sim_flag == 0){
    bf_coefficients = (float*)calloc(N_COEFF, sizeof(float)); // Beamformer coefficients
    coeff_pin(bf_coefficients);
  }

  // -----------------Get phase solutions----------------- //

  // Make all initializations before while loop
  // Initialize beamformer (allocate all memory on the device)
  printf("CBF: Initializing beamformer...\n");
  init_beamformer();

  // Initialize output data array
  float* output_data;

  printf("CBF: Using host arrays allocated in pinned memory\n\r");
  for (int i = 0; i < N_INPUT_BLOCKS; i++){
    ///////////////////////////////////////////////////////////////////////////////
    //>>>>   Register host array in pinned memory <<<<
    ///////////////////////////////////////////////////////////////////////////////
    input_data_pin((signed char *)&db->block[i].data);
  }

  while (run_threads()) {

    /* Note waiting status */
    hashpipe_status_lock_safe(st);
    hputs(st->buf, status_key, "waiting");
    hashpipe_status_unlock_safe(st);

    /* Wait for buf to have data */
    rv = hpguppi_input_databuf_wait_filled(db, curblock);
    if (rv!=0) continue;

    /* Read param struct for this block */
    ptr = hpguppi_databuf_header(db, curblock);
    if (first) {
      hpguppi_read_obs_params(ptr, &gp, &pf);
      first = 0;
    } else {
      hpguppi_read_subint_params(ptr, &gp, &pf);
    }

    /* Read pktidx, pktstart, pktstop from header */
    hgeti8(ptr, "PKTIDX", &pktidx);
    hgeti8(ptr, "PKTSTART", &pktstart);
    hgeti8(ptr, "PKTSTOP", &pktstop);

    /* Get values for calculations at varying points in processing */
    hgetu8(ptr, "SYNCTIME", &synctime);
    hgetu8(ptr, "HCLOCKS", &hclocks);
    hgetu4(ptr, "FENCHAN", &fenchan);
    hgetr8(ptr, "CHAN_BW", &chan_bw); // In MHz
    hgetr8(ptr, "OBSFREQ", &obsfreq);
    hgetu8(ptr, "NANTS", &nants);
    hgetu8(ptr, "BLOCSIZE", &raw_blocsize); // Raw file block size
    hgetu8(ptr, "OBSNCHAN", &obsnchan);
    //hgetu8(ptr, "NPOL", &npol);
    hgetu8(ptr, "PIPERBLK", &piperblk);
    hgetr8(ptr,"TBIN", &tbin);
    hgets(ptr, "OBSID", sizeof(raw_obsid), raw_obsid);
    hgetr8(ptr, "OBSBW", &obsbw);

    //printf("CBF: Number of polarizations (from RAW file) = %d\n", (int)npol);
    
    hashpipe_status_lock_safe(st);
    hgets(st->buf, "BASEFILE", sizeof(raw_basefilename), raw_basefilename);
    hgets(st->buf, "OUTDIR", sizeof(outdir), outdir);
    hashpipe_status_unlock_safe(st);
    printf("CBF: RAW file base filename from command: %s and outdir is: %s \n", raw_basefilename, outdir);
      
    // Get the appropriate basefile name from the rawfile_input_thread 
    // Get HDF5 file data at the beginning of the processing
    if(strcmp(prev_basefilename, raw_basefilename) != 0){
      strcpy(prev_basefilename, raw_basefilename);

      printf("CBF: RAW file base filename:          %s \n", raw_basefilename);
      printf("CBF: Previous RAW file base filename: %s \n", prev_basefilename);

      // Calculate coarse channel center frequencies depending on the mode and center frequency that spans the RAW file
      n_chan_per_node = (int)obsnchan/nants; //((int)fenchan)/n_nodes;
      coarse_chan_band = obsbw/n_chan_per_node; // full_bw/fenchan;

      //printf("CBF: Number of channels per node (from RAW file) = %d\n", n_chan_per_node);
      //printf("CBF: Number of antennas (from RAW file) = %d\n", (int)nants);

      // Skip zeroth index since the number of coarse channels is even and the center frequency is between the 2 middle channels
      for(int i=0; i<n_chan_per_node; i++){
        //coarse_chan_freq[i] = (i-((n_chan_per_node-1)/2))*coarse_chan_band + (obsfreq*1e6);
        coarse_chan_freq[i] = (i-((n_chan_per_node-1)/2))*(chan_bw*1e6) + (obsfreq*1e6);
        // Equivalent equation //
        //coarse_chan_freq[i] = (i-(n_chan_per_node/2))*(chan_bw*1e6) + (obsfreq*1e6) + ((chan_bw/2)*1e6);
      }

      printf("CBF: coarse_chan_freq[%d] = %lf Hz \n", 0, coarse_chan_freq[0]);
      printf("CBF: coarse_chan_freq[%d] = %lf Hz \n", n_chan_per_node/2, coarse_chan_freq[n_chan_per_node/2]);
      printf("CBF: coarse_chan_freq[%d] = %lf Hz \n", n_chan_per_node-1, coarse_chan_freq[n_chan_per_node-1]);

      printf("CBF: coarse_chan_band = %lf Hz \n", coarse_chan_band);
      printf("CBF: chan_bw = %lf MHz \n", chan_bw);

/*
      // Get filterbank file path
      // strrchr() finds the last occurence of the specified character
      char_offset = strrchr(raw_basefilename, character);
      slash_pos = char_offset-raw_basefilename;
      printf("CBF: The last position of %c is %ld \n", character, slash_pos);

      // Get file name with no path
      memcpy(new_base, &raw_basefilename[slash_pos+1], sizeof(raw_basefilename)-slash_pos);
      printf("CBF: File name with no path: %s \n", new_base);
*/
      // Set specified path to write filterbank files
      strcpy(fb_basefilename, outdir);
      strcat(fb_basefilename, "/");
      strcat(fb_basefilename, raw_basefilename);
      printf("CBF: Filterbank file name with new path and no file number or extension yet: %s \n", fb_basefilename);

      if(sim_flag == 0){
        hashpipe_status_lock_safe(st);
        hgets(st->buf, "BFRDIR", sizeof(bfrdir), bfrdir);
        hashpipe_status_unlock_safe(st);
        printf("CBF: bfrdir is: %s \n", bfrdir);

        // Set specified path to read from HDF5 files
        strcpy(hdf5_basefilename, bfrdir);
        strcat(hdf5_basefilename, "/");
        strcat(hdf5_basefilename, raw_basefilename);
        strcat(hdf5_basefilename, ".bfr5");
        printf("CBF: HDF5 file name with path: %s \n", hdf5_basefilename);

        // Read HDF5 file and get all necessary parameters (obsid, cal_all, delays, rates, time_array)
        // Open an existing file. //
        file_id = H5Fopen(hdf5_basefilename, H5F_ACC_RDONLY, H5P_DEFAULT);

        // -------------Read obsid first----------------- //
        // Open an existing dataset. //
        obs_id = H5Dopen(file_id, "/obsinfo/obsid", H5P_DEFAULT);
        // Get obsid data type //
        obs_type = H5Dget_type(obs_id);
        native_obs_type = H5Tget_native_type(obs_type, H5T_DIR_DEFAULT);
        hdf5_obsidsize = (int)H5Tget_size(native_obs_type);
        printf("obsid string size = %d\n", hdf5_obsidsize);
        // Allocate memory to string array
        char hdf5_obsid[hdf5_obsidsize+1];
        hdf5_obsid[hdf5_obsidsize] = '\0'; // Terminate string
        // Read the dataset. // 
        status = H5Dread(obs_id, native_obs_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, hdf5_obsid);
        printf("CBF: obsid = %s \n", hdf5_obsid);
        // Close the dataset. //
        status = H5Dclose(obs_id);
        // -----------------------------------------------//

        // Need to initialize these obsid variables
        // Figure out how to read them first
        if(raw_obsid == hdf5_obsid){
          printf("CBF: OBSID in RAW file and HDF5 file match!\n");
          printf("CBF: raw_obsid  = %s \n", raw_obsid);
          printf("CBF: hdf5_obsid = %s \n", hdf5_obsid);
        }else{
          printf("CBF (Warning): OBSID in RAW file and HDF5 file DO NOT match!\n");
          printf("CBF: raw_obsid  = %s \n", raw_obsid);
          printf("CBF: hdf5_obsid = %s \n", hdf5_obsid);
        }
        // Read cal_all once per HDF5 file. It doesn't change through out the entire recording.
        // Read delayinfo once, get all values, and update when required
        // Open an existing datasets //
        cal_all_id = H5Dopen(file_id, "/calinfo/cal_all", H5P_DEFAULT);
        delays_id = H5Dopen(file_id, "/delayinfo/delays", H5P_DEFAULT);
        time_array_id = H5Dopen(file_id, "/delayinfo/time_array", H5P_DEFAULT);
        npol_id = H5Dopen(file_id, "/diminfo/npol", H5P_DEFAULT);
        nbeams_id = H5Dopen(file_id, "/diminfo/nbeams", H5P_DEFAULT);

        // Get dataspace ID //
        sid1 =  H5Dget_space(cal_all_id);
        sid2 =  H5Dget_space(delays_id);
        sid4 =  H5Dget_space(time_array_id);
  
        // Gets the number of elements in the data set //
        cal_all_elements=H5Sget_simple_extent_npoints(sid1);
        delays_elements=H5Sget_simple_extent_npoints(sid2);
        time_array_elements=H5Sget_simple_extent_npoints(sid4);
        printf("CBF: Number of elements in the cal_all dataset is : %d\n", cal_all_elements);
        printf("CBF: Number of elements in the delays dataset is : %d\n", delays_elements);
        printf("CBF: Number of elements in the time_array dataset is : %d\n", time_array_elements);

        // Allocate memory for array
        cal_all_data = malloc((int)cal_all_elements*sizeof(complex_t));
        delays_data = malloc((int)delays_elements*sizeof(double));
        time_array_data = malloc((int)time_array_elements*sizeof(double));

        // Read the dataset. //
        status = H5Dread(npol_id, H5T_STD_I64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &npol);
        printf("npol = %lu \n", npol);

        status = H5Dread(nbeams_id, H5T_STD_I64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nbeams);
        printf("nbeams = %lu \n", nbeams);

        status = H5Dread(cal_all_id, reim_tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, cal_all_data);
        printf("CBF: cal_all_data[%d].re = %f \n", cal_all_idx(a, p, c, (int)nants, (int)npol), cal_all_data[cal_all_idx(a, p, c, (int)nants, (int)npol)].re);

        status = H5Dread(delays_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, delays_data);
        printf("CBF: delays_data[%d] = %lf \n", delay_idx(a, b, t, (int)nants, (int)nbeams), delays_data[delay_idx(a, b, t, (int)nants, (int)nbeams)]);

        status = H5Dread(time_array_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, time_array_data);
        printf("CBF: time_array_data[0] = %lf \n", time_array_data[0]);

        // Close the dataset. //
        status = H5Dclose(cal_all_id);
        status = H5Dclose(delays_id);
        status = H5Dclose(time_array_id);
        status = H5Dclose(npol_id);
        status = H5Dclose(nbeams_id);

        // Close the file. //
        status = H5Fclose(file_id);
        printf("status of file close: %d\n", status);

        // Reset indexing of time array since new HDF5 file is in use
        time_array_idx = 0;

        // Assign values to tmp variable then copy values from it to pinned memory pointer (bf_coefficients)
        tmp_coefficients = generate_coefficients(cal_all_data, delays_data, time_array_idx, coarse_chan_freq, (int)npol, (int)nbeams, n_chan_per_node, nants);
        memcpy(bf_coefficients, tmp_coefficients, N_COEFF*sizeof(float));
      }

      if(sim_flag == 1){
        nbeams = 1; // Generate one filterbank file with the boresight beam
        npol = 2; // Number of polarizations needs to be specified
        printf("Number of beams = %lu, just the boresight beam \n", nbeams);
      }

      // Number of time samples per block in a RAW file
      n_samp = (int)(raw_blocsize/(2*obsnchan*npol));

      // Number of STI windows
      n_win = n_samp/N_TIME_STI;

      printf("CBF: Got center frequency, obsfreq = %lf MHz, n. time samples = %d \n", obsfreq, n_samp);

      // Size of beamformer output
      blocksize=((N_BF_POW*n_chan_per_node*n_win)/(N_BEAM*MAX_COARSE_FREQ*N_STI))*sizeof(float); 

      // If this is a new scan, write headers and so on to new filterbank files
      got_packet_0 = 0;
    }

    if(sim_flag == 0){
      // Compute unix time in the middle of the block and avg. value current and next time_array value for comparison
      // Updata coefficient if unix time in the middle of the block is greater than avg. of the two time array values
      pktidx_time = synctime + (pktidx*tbin*n_samp/piperblk);
      block_midtime = pktidx_time + tbin*n_samp/2;
      if(time_array_idx < (time_array_elements-1)){
        time_array_midpoint = (time_array_data[time_array_idx] + time_array_data[time_array_idx+1])/2;
      }else if(time_array_idx >= (time_array_elements-1)){
        time_array_midpoint = time_array_data[time_array_idx];
        time_array_idx = time_array_elements-2;
      }

      printf("Unix time in the middle of the block (pktidx_time + tbin*n_samp/2) = %lf\n", block_midtime);
      printf("(time_array_data[%d] + time_array_data[%d])/2 = %lf\n", time_array_idx, time_array_idx+1, time_array_midpoint);

      // Update coefficients every specified number of blocks
      if(block_midtime >= time_array_midpoint){
        printf("Updating coefficients since block_midtime >= time_array_midpoint! \n");
        // Then update with delays[n+1]
        time_array_idx += 1;

        // Update coefficients with new delay
        // Assign values to tmp variable then copy values from it to pinned memory pointer (bf_coefficients)
        tmp_coefficients = generate_coefficients(cal_all_data, delays_data, time_array_idx, coarse_chan_freq, (int)npol, (int)nbeams, n_chan_per_node, nants);
        memcpy(bf_coefficients, tmp_coefficients, N_COEFF*sizeof(float));
       
      }
    }

    // If packet idx is NOT within start/stop range
    if(pktidx < pktstart || pktstop <= pktidx) {
      printf("CBF: Before checking whether files are open \n");
      for(int b = 0; b < nbeams; b++){
	// If file open, close it
	if(fdraw[b] != -1) {
	  // Close file
	  close(fdraw[b]);
	  // Reset fdraw, got_packet_0, filenum, block_count
	  fdraw[b] = -1;
	  if(b == 0){ // These variables only need to be set to zero once
	    got_packet_0 = 0;
	    //filenum = 0;
	    block_count=0;
	    // Print end of recording conditions
	    hashpipe_info(thread_name, "recording stopped: "
	      "pktstart %lu pktstop %lu pktidx %lu",
	      pktstart, pktstop, pktidx);
	  }
	}
      }
      printf("CBF: Before marking as free \n");
      /* Mark as free */
      hpguppi_input_databuf_set_free(db, curblock);

      /* Go to next block */
      curblock = (curblock + 1) % db->header.n_block;

      continue;
    }

    /* Set up data ptr for quant routines */
    pf.sub.data = (unsigned char *)hpguppi_databuf_data(db, curblock);

    // Wait for packet 0 before starting write
    // "packet 0" is the first packet/block of the new recording,
    // it is not necessarily pktidx == 0.
    if (got_packet_0==0 && gp.stt_valid==1) {
      got_packet_0 = 1;

      char fname[256];
      // Create the output directory if needed
      char datadir[1024];
      strncpy(datadir, fb_basefilename, 1023);
      char *last_slash = strrchr(datadir, '/');
      if (last_slash!=NULL && last_slash!=datadir) {
	*last_slash = '\0';
	hashpipe_info(thread_name,
	  "Using directory '%s' for output", datadir);
        if(mkdir_p(datadir, 0755) == -1) {
	  hashpipe_error(thread_name, "mkdir_p(%s)", datadir);
          pthread_exit(NULL);
        }
      }

      // Update filterbank headers based on raw params and Nts etc.
      // Possibly here
      printf("CBF: update_fb_hdrs_from_raw_hdr_cbf(fb_hdr, ptr) \n");
      update_fb_hdrs_from_raw_hdr_cbf(fb_hdr, ptr);

      // Open nbeams filterbank files to save a beam per file i.e. N_BIN*n_samp*sizeof(float) per file.
      printf("CBF: Opening filterbank files \n");
      for(int b = 0; b < nbeams; b++){
        if(b >= 0 && b < 10) {
          if(sim_flag == 0){
	    sprintf(fname, "%s_B0%d.%04d.fil",  fb_basefilename, b, filenum);
          }else if(sim_flag == 1){
            sprintf(fname, "%s_B0.%04d.fil",  fb_basefilename, filenum);
          }
        }else{
          sprintf(fname, "%s_B%d.%04d.fil", fb_basefilename, b, filenum);
        }
        hashpipe_info(thread_name, "Opening fil file '%s'", fname);
	  last_slash = strrchr(fname, '/');
        if(last_slash) {
	  strncpy(fb_hdr.rawdatafile, last_slash+1, 80);
        } else {
	  strncpy(fb_hdr.rawdatafile, fname, 80);
        }
        fb_hdr.rawdatafile[80] = '\0';
                
                
        fdraw[b] = open(fname, O_CREAT|O_WRONLY|O_TRUNC, 0644);
        if(fdraw[b] == -1) {
	  // If we can't open this output file, we probably won't be able to
          // open any more output files, so print message and bail out.
          hashpipe_error(thread_name,
	    "cannot open filterbank output file, giving up");
            pthread_exit(NULL);
        }
    	posix_fadvise(fdraw[b], 0, 0, POSIX_FADV_DONTNEED);
      }
      printf("CBF: Opened filterbank files after for() loop \n");

      //Fix some header stuff here due to multi-antennas
      // Need to understand and appropriately modify these values if necessary
      // Default values for now
      fb_hdr.telescope_id = 64; // MeerKAT ID (Don't know why it's 64, maybe associated with the number of antennas)
      fb_hdr.foff = chan_bw; // Filterbank channel bandwidth
      fb_hdr.nchans = n_chan_per_node; // Number of channels in a filterbank file
      fb_hdr.fch1 = coarse_chan_freq[0]*1e-6; // Center frequency in MHz
      fb_hdr.tsamp = tbin*N_TIME_STI; // Time interval between output samples
      
      /* Write filterbank header to output file */
      printf("CBF: Writing headers to filterbank files! \n");
      for(int b = 0; b < nbeams; b++){
        fb_hdr.ibeam =  b;
        //fb_hdr.source_name = ;
        fb_fd_write_header(fdraw[b], &fb_hdr);
      }
    }

    /* If we got packet 0, process and write data to disk */
    if (got_packet_0) {

      /* Note writing status */
      hashpipe_status_lock_safe(st);
      hputs(st->buf, status_key, "writing");
      hashpipe_status_unlock_safe(st);

      printf("CBF: Before run_beamformer! \n");
      //printf("CBF: Number of antennas (from RAW file) = %d\n", (int)nants);
      //printf("CBF: Number of channels per node (from RAW file) = %d\n", n_chan_per_node);
      //printf("CBF: Number of time samples = %d\n", n_samp);
      //printf("CBF: Number of beams (from HDF5 file) = %d\n", (int)nbeams);
      //printf("CBF: Number of polarizations (from HDF5 file) = %d\n", (int)npol);
      /* Write data */
      // gpu processing function here, I think...

      // Start timing beamforming computation
      struct timespec tval_before, tval_after;
      clock_gettime(CLOCK_MONOTONIC, &tval_before);

      if(sim_flag == 0){
        output_data = run_beamformer((signed char *)&db->block[curblock].data, bf_coefficients, (int)npol, (int)nbeams, n_chan_per_node, n_samp);
      }else if(sim_flag == 1){
        output_data = run_beamformer((signed char *)&db->block[curblock].data, tmp_coefficients, (int)npol, (int)nbeams, n_chan_per_node, n_samp);
      }

      /* Set beamformer output (CUDA kernel before conversion to power), that is summing, to zero before moving on to next block*/
      set_to_zero();

      // Stop timing beamforming computation
      clock_gettime(CLOCK_MONOTONIC, &tval_after);
      time_taken = (float)(tval_after.tv_sec - tval_before.tv_sec); //*1e6; // Time in seconds since epoch
      time_taken = time_taken + (float)(tval_after.tv_nsec - tval_before.tv_nsec)*1e-9; // Time in nanoseconds since 'tv_sec - start and end'
      bf_time = time_taken;

      printf("CBF: run_beamformer() plus set_to_zero() time: %f s\n", bf_time);

      printf("CBF: First element of output data: %f and writing to filterbank files \n", output_data[0]);

      // Start timing write
      struct timespec tval_before_w, tval_after_w;
      clock_gettime(CLOCK_MONOTONIC, &tval_before_w);

      // This may be okay to write to filterbank files, but I'm not entirely confident
      for(int b = 0; b < nbeams; b++){
	//rv = write(fdraw[b], &output_data[b*n_samp*n_chan_per_node], (size_t)blocksize);
        rv = write(fdraw[b], &output_data[b*n_win*n_chan_per_node], (size_t)blocksize);
	if(rv != blocksize){
	  char msg[100];
          perror(thread_name);
	  sprintf(msg, "Error writing data (output_data=%p, blocksize=%d, rv=%d)", output_data, blocksize, rv);
          hashpipe_error(thread_name, msg);
        }

	/* flush output */
	fsync(fdraw[b]);
      }

      // Stop timing write
      clock_gettime(CLOCK_MONOTONIC, &tval_after_w);
      time_taken_w = (float)(tval_after_w.tv_sec - tval_before_w.tv_sec); //*1e6; // Time in seconds since epoch
      time_taken_w = time_taken_w + (float)(tval_after_w.tv_nsec - tval_before_w.tv_nsec)*1e-9; // Time in nanoseconds since 'tv_sec - start and end'
      write_time = time_taken_w;
      printf("CBF: Time taken to write block of size, %d bytes, to disk = %f s \n", (int)nbeams*blocksize, write_time);
      
      printf("CBF: After write() function! Block index = %d \n", block_count);

      /* Increment counter */
      block_count++;
    }

    /* Mark as free */
    hpguppi_input_databuf_set_free(db, curblock);

    /* Go to next block */
    curblock = (curblock + 1) % db->header.n_block;

    /* Check for cancel */
    pthread_testcancel();

  }

  pthread_cleanup_pop(0); // Closes safe_close 

  pthread_cleanup_pop(0); /* Closes hpguppi_free_psrfits */

  // Free up device memory
  cohbfCleanup();

  hashpipe_info(thread_name, "exiting!");
  pthread_exit(NULL);
}

static hashpipe_thread_desc_t rawdisk_thread = {
  name: "hpguppi_coherent_bf_thread",
  skey: "DISKSTAT",
  init: NULL,
  run:  run,
  ibuf_desc: {hpguppi_input_databuf_create},
  obuf_desc: {NULL}
};

static __attribute__((constructor)) void ctor()
{
  register_hashpipe_thread(&rawdisk_thread);
}

// vi: set ts=8 sw=2 noet :
