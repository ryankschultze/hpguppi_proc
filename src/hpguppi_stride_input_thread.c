/* hpguppi_stride_input_thread.c
 *
 * Routine to read GUPPI RAW files and put them 
 * into shared memory blocks. 
 * Can specify output dir if want it different 
 * from input.
 * Author: Mark R. and Cherry Ng
 */
#define MAX_HDR_SIZE (256000)

#define _GNU_SOURCE 1
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <fcntl.h>
#include <string.h>
#include <pthread.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <unistd.h>
#include "hashpipe.h"
#include "hpguppi_databuf.h"
#include "hpguppi_params.h"

#include "upchannelizer_beamformer.h"

// The multiplication by 2 is for real/inphase and imaginary/quadrature components
#define Niq    (2)
#define data_blk_in_idx(p, t, c, a, Np, Nt, Nc)                  (Niq*((p) + (Np)*(t) + (Nt)*(Np)*(c) + (Nc)*(Nt)*(Np)*(a)))
#define data_shm_blk_idx(p, t, b, a, c, Np, Nt, Nb, Na)          (Niq*((p) + (Np)*(t) + (Nt)*(Np)*(b) + (Nb)*(Nt)*(Np)*(a) + (Na)*(Nb)*(Nt)*(Np)*(c)))
#define data_sim_in_idx(p, t, w, a, c, Np, Nt, Nw, Na)           (Niq*((p) + (Np)*(t) + (Nt)*(Np)*(w) + (Nw)*(Nt)*(Np)*(a) + (Na)*(Nw)*(Nt)*(Np)*(c)))
//#define data_shm_blk_idx(p, c, a, t, b, Np, Nc, Na, Nt)          (Niq*((p) + (Np)*(c) + (Nc)*(Np)*(a) + (Na)*(Nc)*(Np)*(t) + (Nt)*(Na)*(Nc)*(Np)*(b)))
//#define data_sim_in_idx(p, c, a, t, w, Np, Nc, Na, Nt)           (Niq*((p) + (Np)*(c) + (Nc)*(Np)*(a) + (Na)*(Nc)*(Np)*(t) + (Nt)*(Na)*(Nc)*(Np)*(w)))
//#define data_sim_in_idx(p, c, a, t, Np, Nc, Na)                  (Niq*((p) + (Np)*(c) + (Nc)*(Np)*(a) + (Na)*(Nc)*(Np)*(t)))
#define VERBOSE 0
#define VERBOSE2 0
#define TIMING 0

// -------------------------------------------------------------- //
// Function that gets the header size of a block in a RAW file
// -------------------------------------------------------------- //
static int get_header_size(int fdin, char * header_buf, size_t len)
{
    int rv;
    int i=0;
    char header_tmp[MAX_HDR_SIZE];
    // Store header_buf info for dummy block at the end of a scan
    memcpy(header_tmp, &header_buf, MAX_HDR_SIZE);
    rv = read(fdin, header_buf, MAX_HDR_SIZE);
    if (rv == -1) {
        hashpipe_error("hpguppi_stride_input_thread", "error reading file");
    } else if (rv > 0) {
        //Read header loop over the 80-byte records
        for (i=0; i<len; i += 80) {
            // If we found the "END " record
            if(!strncmp(header_buf+i, "END ", 4)) {
                // Move to just after END record
                i += 80;
                break;
            }
        }
    }else if(rv == 0){ // End of file has been reached
        // Reinitialize header buffer so it retains the same header info as the previous block for the dummy block
        memcpy(header_buf, &header_tmp, MAX_HDR_SIZE);
        printf("STRIDE INPUT: End of file in get_header_size \n");
    }
    return i;
}

// -------------------------------------------------------------- //
// Get file size to determine number of blocks in the file
// -------------------------------------------------------------- //
static long int get_file_size(int fdin){
    off_t cur_pos = lseek(fdin, (size_t)0, SEEK_CUR);
    off_t file_size = lseek(fdin, (size_t)0, SEEK_END);
    lseek(fdin, cur_pos, SEEK_SET);
    return file_size;
}

static void *run(hashpipe_thread_args_t * args)
{
    hpguppi_input_databuf_t *db  = (hpguppi_input_databuf_t *)args->obuf;
    hashpipe_status_t st = args->st;
    const char * status_key = args->thread_desc->skey;

    // -------------------------------------------------------------- //
    // Reinitialize subband index to 0
    // -------------------------------------------------------------- //
    hputi4(st.buf, "SUBBAND", 0);

    // -------------------------------------------------------------- //
    // Variables and pointers used in main loop
    // -------------------------------------------------------------- //
    int rv;
    int block_idx = 0;
    int block_count=0, filenum=0;
    long int raw_file_size = 0;
    long int cur_pos = 0;
    long int payload_start = 0;
    int end_of_scan = 0; // End of file flag
    int blocsize;
    int nblocks = 0;
    int nants = 0;            // Number of antennas stated in RAW file
    int npol = 0;             // Number of polarizations stated in RAW file
    int obsnchan = 0;         // Number of coarse channels X antennas stated in RAW file
    int n_coarse = 0;         // Number of coarse channels stated in RAW file
    int n_subband = 16;       // Number of sub bands processed serially
    int n_coarse_proc = 0;    // Number of coarse channels processed at one time
    int n_samp_per_block = 0; // Number of time sample in one block
    int piperblk;
    int64_t pktstart;
    int64_t pktstop;
    int64_t cur_pktidx;
    int64_t prev_pktidx;
    int64_t zero_blk_pktidx;
    int n_missed_blks = 0;
    int telescope_flag = 0;
    char *zero_blk;
    zero_blk = (char*)calloc(N_INPUT, sizeof(char));
    char *ptr;

    //Filenames and paths
    char basefilename[200];
    char fname[256];
    
    hgets(st.buf, "BASEFILE", sizeof(basefilename), basefilename);

    char cur_fname[200] = {0};
    char prev_fname[200] = {0};
    char indir[200] = {0};
    strcpy(prev_fname, "prev_fname"); // Initialize as different string that cur_fname
    char *base_pos;
    long int period_pos;
    char character = '/';
    char *char_offset; 
    long int slash_pos;
    char new_base[200];

    char outdir[256];
    hgets(st.buf, "OUTDIR", sizeof(outdir), outdir);
    static int fdin = -1; // Init output file descriptor (-1 means no file open)
    char *header;
    char header_buf[MAX_HDR_SIZE];
    int headersize; // Possibly padded depending on directio flag
    int open_flags = O_RDONLY;
    int directio = 0;

    // -------------------------------------------------------------- //
    // Simulated data set up. Set sim_flag to 1 if simulated data is wanted.
    // -------------------------------------------------------------- //
    int sim_flag = 0; // Set to 1 if you'd like to use simulated data rather than the payload from the RAW file
    char * sim_data; // Initialize simulated data array
    int n_ant_config = 0;
    int n_chan = 0;
    int n_samp = 0;
    int n_pol = 0;
    int n_sim_ant = 0;
    int n_win = 0;
    if(sim_flag == 1){
      n_sim_ant = 58;
      if(n_sim_ant <= N_ANT/2){
        n_ant_config = N_ANT/2;
        n_win = 16;
        // 5 seconds worth of processing at a time
        // 1k mode
        //n_chan = 1; 
        //n_samp = (2*4096*1024)/n_win; // 4194304; // 2^22
        // 4k mode
        n_chan = 4; // 64
        n_samp = (2*1024*1024)/n_win; // 1048576; // 2^20
        // 32k mode
        //n_chan = 32;
        //n_samp = (2*128*1024)/n_win; // 131072; // 2^17
      }else{
        n_ant_config = N_ANT;
        n_win = 8;
        // 5 seconds worth of processing at a time
        // 1k mode
        //n_chan = 1; 
        //n_samp = (4096*1024)/n_win; // 4194304; // 2^22
        // 4k mode
        n_chan = 4; // 64
        n_samp = (1024*1024)/n_win; // 1048576; // 2^20
        // 32k mode
        //n_chan = 32;
        //n_samp = (128*1024)/n_win; // 131072; // 2^17
      }
      n_pol = 2; 
      sim_data = (char *)simulate_data_ubf(n_sim_ant, n_ant_config, n_pol, n_chan, n_samp, n_win, telescope_flag); // Generate block of simulated data
    }
    ssize_t read_blocsize;
#if TIMING
    float read_time = 0;
    float time_taken_r = 0;
#endif

    int wait_filename = 0; // Flag to print "waiting for new RAW file name" only once
    int a = 0; // Antenna index
    int c = 0; // Coarse channel index
    //int t = 0; // Time sample index
    while (run_threads()) {
        hashpipe_status_lock_safe(&st);
        hputs(st.buf, status_key, "waiting");
        hputi4(st.buf, "NETBKOUT", block_idx);
        hashpipe_status_unlock_safe(&st);
        // -------------------------------------------------------------- //
        // Wait for data
        // Wait for new block to be free, then clear it
        // if necessary and fill its header with new values.
        // -------------------------------------------------------------- //
#if VERBOSE
	printf("STRIDE INPUT: while ((rv=hpguppi_input_databuf_wait_free(db, block_idx)) \n");
#endif
        while ((rv=hpguppi_input_databuf_wait_free(db, block_idx)) 
                != HASHPIPE_OK) {
            if (rv==HASHPIPE_TIMEOUT) {
                hashpipe_status_lock_safe(&st);
                hputs(st.buf, status_key, "blocked");
                hashpipe_status_unlock_safe(&st);
                continue;
            } else {
            hashpipe_error(__FUNCTION__, "error waiting for free databuf");
                pthread_exit(NULL);
                break;
            }
        }

#if VERBOSE
	printf("STRIDE INPUT: Before file open if{} \n");
#endif

        // -------------------------------------------------------------- //
        // Stride through RAW file to get all time samples for specific 
        // subband in the RAW file and place in buffer accordingly
        // -------------------------------------------------------------- //
        for(int s = 0; s<n_subband; s++){ // Shift to next subband of 16
            // -------------------------------------------------------------- //
            // Write index of subband to status buffer
            // -------------------------------------------------------------- //
            hputi4(st.buf, "SUBBAND", s);

            // -------------------------------------------------------------- //
            // At the beginning of a file so set end_of_scan == 0
            // -------------------------------------------------------------- //
            end_of_scan = 0;

            // -------------------------------------------------------------- //
            // Iterate through blocks until the end of a sequence of RAW files
            // -------------------------------------------------------------- //
            while(!end_of_scan){
                // -------------------------------------------------------------- //
                // Read raw files
                // -------------------------------------------------------------- //
                if (fdin == -1) { //no file opened
                    // If there is no file ready at the beginning of processing then wait for it to be written to the buffer.
                    while(strlen(cur_fname) == 0){
                        hashpipe_status_lock_safe(&st);
                        hgets(st.buf, "RAWFILE", sizeof(cur_fname), cur_fname);
                        hashpipe_status_unlock_safe(&st);
                    }
                    // Check to see if RAWFILE is an absolute path
                    if(cur_fname[0] != '/'){
                        // If it's not an absolute path, make it so
                        // So wait for the directory of the RAW files. Get it from shared memory
                        while(strlen(indir) == 0){
                            hashpipe_status_lock_safe(&st);
                            hgets(st.buf, "INPUTDIR", sizeof(indir), indir); // Don't have a name for this keyword yet, just going with 'INPUTDIR' for now
                            hashpipe_status_unlock_safe(&st);
                            // Ensure there's a slash at the end of the path
                            if((strlen(indir) != 0) && (indir[(strlen(indir)-1)] != '/')){
                                strcat(indir, "/");
                            }
                            if(strlen(indir) != 0){
                                strcat(indir, cur_fname); // Concatenate the directory and filename
                                strcpy(cur_fname, indir); // Use cur_fname as the current file name variable moving forward
                            }
                        }
                    }
                    // Now create the basefilename
                    // If a '...0000.raw' file exists, that is different from the previous '...0000.raw' file
                    if (strcmp(prev_fname, cur_fname) != 0){
                        strcpy(prev_fname, cur_fname); // Save this file name for comparison on the next iteration

                        base_pos = strchr(cur_fname, '.'); // Finds the first occurence of a period in the filename
                        period_pos = base_pos-cur_fname;

                        memcpy(basefilename, cur_fname, period_pos); // Copy base filename portion of file name to tmp_basefilename variable

                        // Get basefilename with no path and place in status buffer
                        // strrchr() finds the last occurence of the specified character
                        char_offset = strrchr(basefilename, character);
                        slash_pos = char_offset-basefilename;

                        // Get file name with no path
                        memcpy(new_base, &basefilename[slash_pos+1], sizeof(basefilename)-slash_pos);

#if VERBOSE2
                        printf("STRIDE INPUT: Got current RAW file: %s\n", cur_fname);
                        printf("STRIDE INPUT: The last position of . is %ld \n", period_pos);
                        printf("STRIDE INPUT: Base filename from command: %s \n", basefilename);
                        printf("STRIDE INPUT: The last position of %c is %ld \n", character, slash_pos);
                        printf("STRIDE INPUT: File name with no path: %s \n", new_base);
#endif

                        hputs(st.buf, "BASEFILE", new_base);

                    }
                    else{
                        // The RAW file hasn't changed so wait for new file to show up in the buffer
                        if(wait_filename == 0){ // Print "waiting for new RAW file name" only once
                            wait_filename = 1;
                            printf("STRIDE INPUT: Waiting for new RAW file name! \n");
                        }

                        // Will exit if thread has been cancelled
                        pthread_testcancel();

                        continue;
                    }
                    wait_filename = 0; // Print "waiting for new RAW file name" only once
                    sprintf(fname, "%s.%04d.raw", basefilename, filenum);
            
                    printf("STRIDE INPUT: Opening first raw file '%s'\n", fname);
                    fdin = open(fname, open_flags, 0644);
                    if (fdin==-1) {
                        hashpipe_error(__FUNCTION__,"Error opening file.");
                        pthread_exit(NULL);
                    }

                    // Get raw file size in order to calculate the number of blocks in the file
                    raw_file_size = get_file_size(fdin);

                    headersize= get_header_size(fdin, header_buf, MAX_HDR_SIZE);

                    hgeti4(header_buf, "BLOCSIZE", &blocsize);

                    nblocks = raw_file_size/(headersize+blocsize); // Assume this value is always an integer value for now.
                    printf("STRIDE INPUT: Number of RAW blocks = %d \n", nblocks);
#if VERBOSE
                    printf("STRIDE INPUT: raw_file_size = %ld \n", raw_file_size);
                    printf("STRIDE INPUT: headersize = %d \n", headersize);
                    printf("STRIDE INPUT: blocsize = %d \n", blocsize);
#endif

                    // Reset file position back to the beginning of the file
                    lseek(fdin, 0, SEEK_SET);
                }

                // -------------------------------------------------------------- //
                // Get the current position of the file
                // -------------------------------------------------------------- //
                cur_pos = lseek(fdin, (size_t)0, SEEK_CUR);

#if VERBOSE
                printf("STRIDE INPUT: current position = %ld and file number = %d \n", cur_pos, filenum);
#endif

                // -------------------------------------------------------------- //
                //Handling header - size, output path, directio
                // -------------------------------------------------------------- //
                headersize= get_header_size(fdin, header_buf, MAX_HDR_SIZE);
                printf("STRIDE INPUT: current position = %ld, raw_file_size = %ld and headersize = %d \n", cur_pos, raw_file_size, headersize);

                // -------------------------------------------------------------- //
                // If we are not at the end of the file, read blocks and transfer 
                // to shared mem buffer else if we are at the end, move on to the 
                // next file
                // -------------------------------------------------------------- //
                if(((raw_file_size-cur_pos)>=(blocsize+headersize))){
#if VERBOSE
                    printf("STRIDE INPUT: In if((raw_file_size-cur_pos)>=(blocsize+headersize)){ \n");
#endif
                    // If a block is missing, copy a block of zeros to the buffer in it's place
                    // Otherwise, write the data from a block to the buffer
                    if(n_missed_blks == 0){
#if VERBOSE
                        printf("STRIDE INPUT: In if(n_missed_blks == 0){ \n");
#endif

                        if(block_count == 0){
#if VERBOSE
                        printf("STRIDE INPUT: In if(block_count == 0){ \n");
#endif
                            //set_output_path(header_buf, outdir, MAX_HDR_SIZE);
                            hputs(header_buf, "DATADIR", outdir);

                            // Initialize header of block in buffer
                            header = hpguppi_databuf_header(db, block_idx);

                            // Initialize block in buffer
                            ptr = hpguppi_databuf_data(db, block_idx);
                        }

                        hashpipe_status_lock_safe(&st);
                        hputs(st.buf, status_key, "receiving");
                        memcpy(header, &header_buf, headersize);
                        hashpipe_status_unlock_safe(&st);

                        directio = hpguppi_read_directio_mode(header);

                        // Adjust length for any padding required for DirectIO
                        if(directio) {
                            // Round up to next multiple of 512
                            headersize = (headersize+511) & ~511;
                        }
                        payload_start = lseek(fdin, headersize-MAX_HDR_SIZE, SEEK_CUR);
                        hgeti4(header_buf, "BLOCSIZE", &blocsize);
                        hgeti8(header_buf, "PKTIDX", &cur_pktidx);
                        hgeti8(header_buf, "PKTSTART", &pktstart);
                        hgeti8(header_buf, "PKTSTOP", &pktstop);
                        hgeti4(header_buf, "PIPERBLK", &piperblk);
                        // Descriptions of the variables below are at the initializations
                        hgeti4(header_buf, "NANTS", &nants);
                        hgeti4(header_buf, "NPOL", &npol);
                        if(npol > 1){
                            npol = 2;
                        }
                        hgeti4(header_buf, "OBSNCHAN", &obsnchan);
                        n_coarse = (int)obsnchan/nants;
                        n_coarse_proc = n_coarse/n_subband;
                        n_samp_per_block = (int)(blocsize/(2*obsnchan*npol));

#if VERBOSE2
                        printf("STRIDE INPUT: headersize = %d, directio = %d\n", headersize, directio);
                        printf("STRIDE INPUT: payload_start = %ld \n", payload_start);
                        printf("STRIDE INPUT: blocsize = %d \n", blocsize);
                        printf("STRIDE INPUT: cur_pktidx = %ld \n", cur_pktidx);
                        printf("STRIDE INPUT: pktstart = %ld \n", pktstart);
                        printf("STRIDE INPUT: pktstop = %ld \n", pktstop);
                        printf("STRIDE INPUT: piperblk = %d \n", piperblk);
                        printf("STRIDE INPUT: nants = %d \n", nants);
                        printf("STRIDE INPUT: npol = %d \n", npol);
                        printf("STRIDE INPUT: obsnchan = %d \n", obsnchan);
                        printf("STRIDE INPUT: n_coarse = %d \n", n_coarse);
                        printf("STRIDE INPUT: n_coarse_proc = %d \n", n_coarse_proc);
                        printf("STRIDE INPUT: n_samp_per_block = %d \n", n_samp_per_block);
                        printf("STRIDE INPUT: Current pktidx = %ld, pktstart = %ld, and pktstop = %ld \n", cur_pktidx, pktstart, pktstop);
                        printf("STRIDE INPUT: piperblk = %d \n", piperblk);
#endif

                        // If current packet index is greater than packet start,
                        // check to see whether there are missed blocks
                        if(cur_pktidx > pktstart){
                            // At the beginning of a RAW file, prev_pktidx is set to 0 
                            // in order to verify whether there are missing blocks at 
                            // the beginning of the file
                            if(prev_pktidx == 0){
                                prev_pktidx = pktstart;
                            }
                            // If (cur_pktidx - prev_pktidx) > piperblk, calculate the number of missed blocks
                            // and start over to write zero blocks to shared memory buffer
                            if((cur_pktidx - prev_pktidx) > piperblk){
                                n_missed_blks = ((cur_pktidx - prev_pktidx)/piperblk)-1; // Minus 1 because the block with the prev_pktidx has data
                                zero_blk_pktidx = cur_pktidx;
                                continue;
                            }
                            // Set previous PKTIDX
                            prev_pktidx = cur_pktidx;
                        }

                        //Read data--------------------------------------------------
                        // Start timing read
#if TIMING
                        struct timespec tval_before, tval_after;
                        clock_gettime(CLOCK_MONOTONIC, &tval_before);
#endif

#if VERBOSE
                        printf("STRIDE INPUT: RAW block size: %d, and  BLOCK_DATA_SIZE: %d \n", blocsize, BLOCK_DATA_SIZE);
                        printf("STRIDE INPUT: header size: %d, and  MAX_HDR_SIZE: %d \n", headersize, MAX_HDR_SIZE);
#endif
                        for(a = 0; a<nants; a++){
                            for(c = 0; c<n_coarse_proc; c++){
                                // Reset to beginning of block after headert (start position of the payload)
                                lseek(fdin, payload_start, SEEK_SET);
                                // Offset by antenna and coarse channel index taking the subband index into account as well
                                lseek(fdin, data_blk_in_idx(0, 0, (c + n_coarse_proc*s), a, npol, n_samp_per_block, n_coarse), SEEK_CUR);
                                // Place input data in shared memory buffer (from RAW file or simulated)
                                if(sim_flag == 0){
                                    // Read the remaining dimensions in the RAW block (npol and n_samp_per_block) and place in the shared memory buffer block
                                    //for(t = 0; t<n_samp_per_block; t++){
                                    read_blocsize = read(fdin, &ptr[data_shm_blk_idx(0,0,block_count,a,c,npol,n_samp_per_block,nblocks,nants)], Niq*npol*n_samp_per_block);
                                    //}
                                } else{
                                    // Copy simulated data to shared memory buffer block - data_sim_in_idx(p, t, w, a, c, Np, Nt, Nw, Na)
                                    //for(int t = 0; t<n_samp; t++){
                                    memcpy(&ptr[data_sim_in_idx(0,0,0,a,c,n_pol,n_samp,n_win,N_ANT)], &sim_data[data_sim_in_idx(0,0,0,a,c,n_pol,n_samp,n_win,N_ANT)], Niq*n_pol*n_samp*n_win);
                                    //}
                                }
                            }
                        }

                        // Offset to end of a block (move through all subbands to get to the next block)
                        cur_pos = lseek(fdin, data_blk_in_idx(0, 0, (n_coarse-(n_coarse_proc*(s+1))), 0, npol, n_samp_per_block, n_coarse), SEEK_CUR);

#if VERBOSE
                        printf("STRIDE INPUT: After read current position = %ld \n", cur_pos);
#endif

#if TIMING
                        // Stop timing read
                        clock_gettime(CLOCK_MONOTONIC, &tval_after);
                        time_taken_r = (float)(tval_after.tv_sec - tval_before.tv_sec); //*1e6; // Time in seconds since epoch
                        time_taken_r = time_taken_r + (float)(tval_after.tv_nsec - tval_before.tv_nsec)*1e-9; //*1e-6; // Time in nanoseconds since 'tv_sec - start and end'
                        read_time = time_taken_r;

                        printf("STRIDE INPUT: Time taken to read from RAW file = %f ms \n", read_time);
#endif

                        printf("STRIDE INPUT: Subband index = %d and Block count = %d \n", s, block_count);
                        // Iterate through blocks in RAW files for a scan corresponding to a subband
                        block_count++;
                    } else if(n_missed_blks > 0){
                        if(block_count == 0){
                            headersize= get_header_size(fdin, header_buf, MAX_HDR_SIZE);
                            hgeti4(header_buf, "PIPERBLK", &piperblk);

                            // Copy the same header info from previous block to zero block in buffer
                            header = hpguppi_databuf_header(db, block_idx);
                            hashpipe_status_lock_safe(&st);
                            hputs(st.buf, status_key, "receiving");
                            memcpy(header, &header_buf, headersize);
                            hashpipe_status_unlock_safe(&st);

                            //printf("STRIDE INPUT: directio = hpguppi_read_directio_mode(header); \n");
                            directio = hpguppi_read_directio_mode(header);

#if VERBOSE2
                            printf("STRIDE INPUT: headersize = %d, directio = %d\n", headersize, directio);
#endif

                            // Adjust length for any padding required for DirectIO
                            if(directio) {
                                // Round up to next multiple of 512
                                headersize = (headersize+511) & ~511;
                            }

                            ptr = hpguppi_databuf_data(db, block_idx);
                        }

                        // Set pktidx in the header of the shared memory buffer block (header_buf)
                        zero_blk_pktidx -= n_missed_blks*piperblk;
                        hputi8(header_buf, "PKTIDX", zero_blk_pktidx);

                        // Copy block of zeros to block in buffer
                        memcpy(&ptr[block_count*Niq*npol*n_coarse_proc*nants*n_samp_per_block], zero_blk, Niq*npol*n_coarse_proc*nants*n_samp_per_block);
            
                        // Decrement n_missed_blks by 1
                        n_missed_blks -= 1;

                        // If n_missed_blks == 0, then process the data at the packet index after blocks of zeros (cur_pktidx)
                        if(n_missed_blks == 0){

                            for(a = 0; a<nants; a++){
                                for(c = 0; c<n_coarse_proc; c++){
                                    // Reset to beginning of block after headert (start position of the payload)
                                    lseek(fdin, payload_start, SEEK_SET);
                                    // Offset by antenna and coarse channel index taking the subband index into account as well
                                    lseek(fdin, data_blk_in_idx(0, 0, (c + n_coarse_proc*s), a, npol, n_samp_per_block, n_coarse), SEEK_CUR);
                                    // Place input data in shared memory buffer from RAW file
                                    // Read the remaining dimensions in the RAW block (npol and n_samp_per_block) and place in the shared memory buffer block
                                    //for(int t = 0; t<n_samp_per_block; t++){
                                    read_blocsize = read(fdin, &ptr[data_shm_blk_idx(0,0,block_count,a,c,npol,n_samp_per_block,nblocks,nants)], Niq*npol*n_samp_per_block);
                                    //}
                                }
                            }

                            cur_pos = lseek(fdin, data_blk_in_idx(0, 0, (n_coarse-(n_coarse_proc*(s+1))), 0, npol, n_samp_per_block, n_coarse), SEEK_CUR);
                            printf("STRIDE INPUT: Missed blocks replaced with zeros. After read current position = %ld \n", cur_pos);

                            if(block_count >= 0 && (block_count <= 5)){
                                printf("STRIDE INPUT: Number of bytes read in read(): %zd \n", read_blocsize);
                            }
#if VERBOSE2
                            printf("STRIDE INPUT: First element of buffer: %d \n", ptr[0]);
#endif
                            // Set previous PKTIDX
                            prev_pktidx = cur_pktidx;
                        }

                        // Increment the block count even with zero blocks
                        block_count++;
                    }
                }// If (cur_pos > (raw_file_size-blocsize)) && (cur_pos <= raw_file_size),
                // then we have reached the end of a file so move on to next file or wait for new file and set PKTIDX == PKTSTOP if necessary
                //else if((cur_pos > (raw_file_size-(blocsize+headersize))) && (cur_pos < raw_file_size)){
                else if(((raw_file_size-cur_pos)<(blocsize+headersize))){
                    printf("STRIDE INPUT: Reached end of file! \n");
                    close(fdin);
                    printf("STRIDE INPUT: Closed RAW file! \n");
                    filenum++;

                    // Inform downstream thread about the number of time samples in a RAW file
                    hputi4(st.buf, "NSAMP", block_count*n_samp_per_block);

                    sprintf(fname, "%s.%4.4d.raw", basefilename, filenum);
                    printf("STRIDE INPUT: Opening next raw file '%s'\n", fname);
                    fdin = open(fname, open_flags, 0644);

                    if(fdin != -1){
                        // Get raw file size in order to calculate the number of blocks in the file
                        raw_file_size = get_file_size(fdin);
                    }
                    printf("STRIDE INPUT: fdin = %d, n_samp = %d and raw_file_size = %ld \n", fdin, block_count*n_samp_per_block, raw_file_size);
                    if ((fdin==-1) && (s < (n_subband-1))) { // End of a sequence of RAW files corresponding to a scan and less than the no. of subbands
                        // Start the next scan with file number 0
                        filenum=0;

                        sprintf(fname, "%s.%4.4d.raw", basefilename, filenum);
                        printf("STRIDE INPUT: Opening first raw file '%s' of scan for subband = %d of %d \n", fname, (s+1), n_subband);
                        fdin = open(fname, open_flags, 0644);

                        // Get raw file size in order to calculate the number of blocks in the file
                        raw_file_size = get_file_size(fdin);

                        // Reset previous pktidx to 0 for the next scan
                        prev_pktidx = 0;

                        // End of scan so move on to the next subband
                        end_of_scan = 1;
                    }else if ((fdin==-1) && (s == (n_subband-1))) { // End of a sequence of RAW files corresponding to a scan and end of subbands
                        // Start the next scan with file number 0
                        filenum=0;

                        printf("STRIDE INPUT: cur_pktidx = %ld and pktstop = %ld \n", cur_pktidx, pktstop);

                        // Inform the downstream thread that we have reached the end of a scan
                        if(cur_pktidx < pktstop){
                            hgeti8(header_buf, "PKTIDX", &cur_pktidx);

                            // Create a header for a dummy block
                            header = hpguppi_databuf_header(db, block_idx);
                            hashpipe_status_lock_safe(&st);
                            hputs(st.buf, status_key, "receiving");
                            memcpy(header, &header_buf, headersize);
                            hashpipe_status_unlock_safe(&st);

                            // Send dummy block with PKTIDX set to PKTSTOP (Make sure that PKTIDX is set to PKTSTOP)
                            hputi8(header, "PKTIDX", pktstop);

                            // Initialize block
                            ptr = hpguppi_databuf_data(db, block_idx);

                            // Copy block of zeros to block in buffer
                            memcpy(ptr, zero_blk, N_INPUT);
                            printf("STRIDE INPUT: After memcpy() N_INPUT = %ld \n", N_INPUT);
                        }

                        // Reset previous pktidx to 0 for the next scan
                        prev_pktidx = 0;

                        // End of scan so move on to the next subband
                        end_of_scan = 1;
                    }

                    // Reset block_count to 0 (the index of blocks in the RAW files of a subband)
                    block_count=0;

                    // Mark block as full
                    hpguppi_input_databuf_set_filled(db, block_idx);
                    printf("STRIDE INPUT: After hpguppi_input_databuf_set_filled() block_idx = %d \n", block_idx);

                    // Setup for next block
                    block_idx = (block_idx + 1) % N_INPUT_BLOCKS;

                    // Wait to allow the processing thread to open and read bfr5 file, and open and write to filterbank files
                    sleep(2);
                }// Sometimes there might be headers with no data blocks remaining in the RAW file. So copy zero blocks to buffer blocks with the appropriate header
                else if((cur_pos+(blocsize+headersize)) > raw_file_size){
                    // If we haven't reached the last header, then increment block_count (the index of blocks in the RAW files of a subband)
                    if(cur_pos < (raw_file_size-headersize)){
                        //Handling header - size, output path, directio----------------
#if VERBOSE
                        printf("STRIDE INPUT: int headersize= get_header_size(fdin, header_buf, MAX_HDR_SIZE); \n");
#endif
                        headersize= get_header_size(fdin, header_buf, MAX_HDR_SIZE);

                        if(block_count == 0){
                            //set_output_path(header_buf, outdir, MAX_HDR_SIZE);
                            hputs(header_buf, "DATADIR", outdir);
                            printf("STRIDE INPUT: hputs(header_buf, DATADIR, outdir); \n");

                            // Initialize header of block in buffer
                            header = hpguppi_databuf_header(db, block_idx);

                            // Initialize block in buffer
                            ptr = hpguppi_databuf_data(db, block_idx);
                        }

                        hashpipe_status_lock_safe(&st);
                        hputs(st.buf, status_key, "receiving");
                        memcpy(header, &header_buf, headersize);
                        hashpipe_status_unlock_safe(&st);

                        directio = hpguppi_read_directio_mode(header);

                        // Adjust length for any padding required for DirectIO
                        if(directio) {
                            // Round up to next multiple of 512
                            headersize = (headersize+511) & ~511;
                        }

#if VERBOSE2
                        printf("STRIDE INPUT: headersize = %d, directio = %d\n", headersize, directio);
#endif

                        payload_start = lseek(fdin, headersize-MAX_HDR_SIZE, SEEK_CUR);
                        hgeti4(header_buf, "BLOCSIZE", &blocsize);
                        hgeti8(header_buf, "PKTIDX", &cur_pktidx);
                        hgeti8(header_buf, "PKTSTART", &pktstart);
                        hgeti8(header_buf, "PKTSTOP", &pktstop);
                        hgeti4(header_buf, "PIPERBLK", &piperblk);
                        // Descriptions of the variables below are at the initializations
                        hgeti4(header_buf, "NANTS", &nants);
                        hgeti4(header_buf, "NPOL", &npol);
                        if(npol > 1){
                            npol = 2;
                        }
                        hgeti4(header_buf, "OBSNCHAN", &obsnchan);
                        n_coarse = (int)obsnchan/nants;
                        n_coarse_proc = n_coarse/n_subband;
                        n_samp_per_block = (int)(blocsize/(2*obsnchan*npol));


#if VERBOSE2
                        printf("STRIDE INPUT: payload_start = %ld \n", payload_start);
                        printf("STRIDE INPUT: blocsize = %d \n", blocsize);
                        printf("STRIDE INPUT: cur_pktidx = %ld \n", cur_pktidx);
                        printf("STRIDE INPUT: pktstart = %ld \n", pktstart);
                        printf("STRIDE INPUT: pktstop = %ld \n", pktstop);
                        printf("STRIDE INPUT: piperblk = %d \n", piperblk);
                        printf("STRIDE INPUT: nants = %d \n", nants);
                        printf("STRIDE INPUT: npol = %d \n", npol);
                        printf("STRIDE INPUT: obsnchan = %d \n", obsnchan);
                        printf("STRIDE INPUT: n_coarse = %d \n", n_coarse);
                        printf("STRIDE INPUT: n_coarse_proc = %d \n", n_coarse_proc);
                        printf("STRIDE INPUT: n_samp_per_block = %d \n", n_samp_per_block);
                        printf("STRIDE INPUT: Current pktidx = %ld, pktstart = %ld, and pktstop = %ld \n", cur_pktidx, pktstart, pktstop);
                        printf("STRIDE INPUT: piperblk = %d \n", piperblk);
#endif

                        // Write blocks of zeros to the shared memory buffer blocks
                        // Copy block of zeros to block in buffer
                        memcpy(&ptr[block_count*Niq*npol*n_coarse_proc*nants*n_samp_per_block], zero_blk, Niq*npol*n_coarse_proc*nants*n_samp_per_block);

                        block_count++;
                    }
                    // If we reach the last header with no data block,
                    // then move on to the next RAW file or wait for one
                    if(cur_pos >= (raw_file_size-headersize)){
                        close(fdin);
                        filenum++;

                        // Inform downstream thread about the number of time samples in a RAW file
                        hputi4(st.buf, "NSAMP", block_count*n_samp_per_block);

                        sprintf(fname, "%s.%4.4d.raw", basefilename, filenum);
                        printf("STRIDE INPUT: Opening next raw file '%s'\n", fname);
                        fdin = open(fname, open_flags, 0644);
                        if(fdin != -1){
                            // Get raw file size in order to calculate the number of blocks in the file
                            raw_file_size = get_file_size(fdin);
                        }
                        printf("STRIDE INPUT: fdin = %d, n_samp = %d and raw_file_size = %ld \n", fdin, block_count*n_samp_per_block, raw_file_size);
                        if ((fdin==-1) && (s < (n_subband-1))) { // End of a sequence of RAW files corresponding to a scan and less than the no. of subbands
                            // Start the next scan with file number 0
                            filenum=0;

                            sprintf(fname, "%s.%4.4d.raw", basefilename, filenum);
                            printf("STRIDE INPUT: Opening first raw file '%s' of scan for subband = %d of %d \n", fname, s, n_subband);
                            fdin = open(fname, open_flags, 0644);

                            // Get raw file size in order to calculate the number of blocks in the file
                            raw_file_size = get_file_size(fdin);

                            // Reset previous pktidx to 0 for the next scan
                            prev_pktidx = 0;

                            // End of scan so move on to the next subband
                            end_of_scan = 1;
                        }else if ((fdin==-1) && (s == (n_subband-1))) { // End of a sequence of RAW files corresponding to a scan and end of subbands
                            // Start the next scan with file number 0
                            filenum=0;

                            // Inform the downstream thread that we have reached the end of a scan
                            if(cur_pktidx < pktstop){
                                // Create a header for a dummy block
                                header = hpguppi_databuf_header(db, block_idx);
                                hashpipe_status_lock_safe(&st);
                                hputs(st.buf, status_key, "receiving");
                                memcpy(header, &header_buf, headersize);
                                hashpipe_status_unlock_safe(&st);

                                hputi8(header, "PKTIDX", pktstop);

                                // Initialize block
                                ptr = hpguppi_databuf_data(db, block_idx);

                                // Copy block of zeros to block in buffer
                                memcpy(ptr, zero_blk, N_INPUT);
                            }

                            // Reset previous pktidx to 0 for the next scan
                            prev_pktidx = 0;

                            // End of scan so move on to the next subband
                            end_of_scan = 1;
                        }
                    }

                    // Mark block as full
                    hpguppi_input_databuf_set_filled(db, block_idx);

                    // Setup for next block
                    block_idx = (block_idx + 1) % N_INPUT_BLOCKS;

                    // Wait to allow the beamformer to open and read bfr5 file, and open and write to filterbank files
                    sleep(2);
                }
                // Will exit if thread has been cancelled
                pthread_testcancel();
            }
            // Will exit if thread has been cancelled
            pthread_testcancel();
        }
        // Will exit if thread has been cancelled
        pthread_testcancel();
    }

    // Thread success!
    if (fdin!=-1) {
	printf("STRIDE INPUT: Closing file! \n");
        close(fdin);
    }
    return THREAD_OK;

}

static hashpipe_thread_desc_t hpguppi_stride_input_thread = {
    name: "hpguppi_stride_input_thread",
    skey: "NETSTAT",
    init: NULL,
    run:  run,
    ibuf_desc: {NULL},
    obuf_desc: {hpguppi_input_databuf_create}
};

static __attribute__((constructor)) void ctor()
{
    register_hashpipe_thread(&hpguppi_stride_input_thread);
}
