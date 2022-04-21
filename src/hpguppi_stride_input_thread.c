/* hpguppi_rawfile_input_thread.c
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

#define data_blk_in_idx(p, t, c, a, Np, Nt, Nc)              ((p) + (Np)*(t) + (Nt)*(Np)*(c) + (Nc)*(Nt)*(Np)*(a))
//#define data_shm_blk_idx(p, t, b, f, a, Np, Nt, Nb, Nf)      ((p) + (Np)*(t) + (Nt)*(Np)*(b) + (Nb)*(Nt)*(Np)*(f) + (Nf)*(Nb)*(Nt)*(Np)*(a))
//#define data_in_idx(p, t, w, c, a, Np, Nt, Nw, Nc)           ((p) + (Np)*(t) + (Nt)*(Np)*(w) + (Nw)*(Nt)*(Np)*(c) + (Nc)*(Nw)*(Nt)*(Np)*(a))
#define data_shm_blk_idx(p, c, a, t, b, Np, Nc, Na, Nt)      ((p) + (Np)*(c) + (Nc)*(Np)*(a) + (Na)*(Nc)*(Np)*(t) + (Nt)*(Na)*(Nc)*(Np)*(b))
#define data_in_idx(p, c, a, t, w, Np, Nc, Na, Nt)           ((p) + (Np)*(c) + (Nc)*(Np)*(a) + (Na)*(Nc)*(Np)*(t) + (Nt)*(Na)*(Nc)*(Np)*(w))
#define VERBOSE 0
#define TIMING 0

static int get_header_size(int fdin, char * header_buf, size_t len)
{
    int rv;
    int i=0;
    char header_reset[MAX_HDR_SIZE];
    rv = read(fdin, header_buf, MAX_HDR_SIZE);
    if (rv == -1) {
        hashpipe_error("hpguppi_rawfile_input_thread", "error reading file");
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
        // Reinitialize header buffer so it doesn't retain the same header info as the block that has been read
        memcpy(header_buf, &header_reset, MAX_HDR_SIZE);
    }
    return i;
}

// Get file size to determine number of blocks in the file
static long int get_file_size(int fdin){
    off_t cur_pos = lseek(fdin, (size_t)0, SEEK_CUR);
    off_t file_size = lseek(fdin, (size_t)0, SEEK_END);
    lseek(fdin, cur_pos, SEEK_SET);
    return file_size;
}

static int get_block_size(char * header_buf, size_t len)
{
    int i;
    char bs_str[32];
    int blocsize = 0;
    //Read header loop over the 80-byte records
    for (i=0; i<len; i += 80) {
        if(!strncmp(header_buf+i, "BLOCSIZE", 8)) {
            strncpy(bs_str,header_buf+i+16, 32);
            blocsize = strtoul(bs_str,NULL,0);
            break;
        }
    }
    return blocsize;
}

static int get_nants(char * header_buf, size_t len)
{
    int i;
    char ants_str[32];
    int nants = 0;
    //Read header loop over the 80-byte records
    for (i=0; i<len; i += 80) {
        if(!strncmp(header_buf+i, "NANTS", 8)) {
            strncpy(ants_str,header_buf+i+16, 32);
            nants = strtoul(ants_str,NULL,0);
            break;
        }
    }
    return nants;
}

static int get_npol(char * header_buf, size_t len)
{
    int i;
    char npol_str[32];
    int npol = 0;
    //Read header loop over the 80-byte records
    for (i=0; i<len; i += 80) {
        if(!strncmp(header_buf+i, "NPOL", 8)) {
            strncpy(npol_str,header_buf+i+16, 32);
            npol = strtoul(npol_str,NULL,0);
            if(npol>1){
                npol = 2;
            }
            break;
        }
    }
    return npol;
}

static int get_obsnchan(char * header_buf, size_t len)
{
    int i;
    char obs_str[32];
    int obsnchan = 0;
    //Read header loop over the 80-byte records
    for (i=0; i<len; i += 80) {
        if(!strncmp(header_buf+i, "OBSNCHAN", 8)) {
            strncpy(obs_str,header_buf+i+16, 32);
            obsnchan = strtoul(obs_str,NULL,0);
            break;
        }
    }
    return obsnchan;
}

static int64_t get_cur_pktidx(char * header_buf, size_t len)
{
    int i;
    char bs_str[32];
    int64_t pktidx = 0;
    //Read header loop over the 80-byte records
    for (i=0; i<len; i += 80) {
        if(!strncmp(header_buf+i, "PKTIDX", 6)) {
            strncpy(bs_str,header_buf+i+16, 32);
            pktidx = strtoul(bs_str,NULL,0);
            break;
        }
    }
    return pktidx;
}

static int64_t get_pktstart(char * header_buf, size_t len)
{
    int i;
    char bs_str[32];
    int64_t pktstart = 0;
    //Read header loop over the 80-byte records
    for (i=0; i<len; i += 80) {
        if(!strncmp(header_buf+i, "PKTSTART", 6)) {
            strncpy(bs_str,header_buf+i+16, 32);
            pktstart = strtoul(bs_str,NULL,0);
            break;
        }
    }
    return pktstart;
}

static int64_t get_pktstop(char * header_buf, size_t len)
{
    int i;
    char bs_str[32];
    int64_t pktstop = 0;
    //Read header loop over the 80-byte records
    for (i=0; i<len; i += 80) {
        if(!strncmp(header_buf+i, "PKTSTOP", 6)) {
            strncpy(bs_str,header_buf+i+16, 32);
            pktstop = strtoul(bs_str,NULL,0);
            break;
        }
    }
    return pktstop;
}
/*
static int64_t get_nxt_pktidx(int fdin, int blocsize, char * header_buf, size_t len)
{
    int i;
    int rv;
    char bs_str[32];
    int64_t pktidx = 0;
    // Go to start of next block
    lseek(fdin, blocsize, SEEK_CUR);

    // I don't know what the next block's header size is and I shouldn't make assumptions so read the max header size and offset back based off that
    rv = read(fdin, header_buf, MAX_HDR_SIZE);
    if (rv == -1) {
        hashpipe_error("hpguppi_rawfile_input_thread", "error reading file");
    } else if (rv > 0) {
        //Read header loop over the 80-byte records
        for (i=0; i<len; i += 80) {
            if(!strncmp(header_buf+i, "PKTIDX", 6)) {
                strncpy(bs_str,header_buf+i+16, 32);
                pktidx = strtoul(bs_str,NULL,0);
                break;
            }
        }
        // Reset to the previous block
        lseek(fdin,-1*(blocsize+MAX_HDR_SIZE),SEEK_CUR);
    }else if(rv == 0){ // End of file has been reached
        // Reset to the current packet's payload
        lseek(fdin, -1*blocsize, SEEK_CUR);
    }
    return pktidx;
}
*/
static int get_piperblk(char * header_buf, size_t len)
{
    int i;
    char bs_str[32];
    int piperblk = 0;
    //Read header loop over the 80-byte records
    for (i=0; i<len; i += 80) {
        if(!strncmp(header_buf+i, "PIPERBLK", 6)) {
            strncpy(bs_str,header_buf+i+16, 32);
            piperblk = strtoul(bs_str,NULL,0);
            break;
        }
    }
    return piperblk;
}

static void set_output_path(char * header_buf, char * outdir, size_t len)
{
    int i;
    //Read header loop over the 80-byte records
    for (i=0; i<len; i += 80) {
        if(!strncmp(header_buf+i, "DATADIR", 7)) {
            hputs(header_buf, "DATADIR", outdir);
            break;
        }
    }
}

static void set_pktidx(char * header_buf, int pktidx, size_t len)
{
    int i;
    //Read header loop over the 80-byte records
    for (i=0; i<len; i += 80) {
        if(!strncmp(header_buf+i, "PKTIDX", 8)) {
            hputi4(header_buf, "PKTIDX", pktidx);
            break;
        }
    }
}

/*
static void set_blocksize(char * header_buf, int blocsize, size_t len)
{
    int i;
    //Read header loop over the 80-byte records
    for (i=0; i<len; i += 80) {
        if(!strncmp(header_buf+i, "BLOCSIZE", 8)) {
            hputi4(header_buf, "BLOCSIZE", blocsize);
            break;
        }
    }
}

static ssize_t read_fully(int fd, void * buf, size_t bytes_to_read)
{
    ssize_t bytes_read;
    ssize_t total_bytes_read = 0;

    while(bytes_to_read > 0) {
        bytes_read = read(fd, buf, bytes_to_read);
        if(bytes_read <= 0) {
            if(bytes_read == 0) {
                break;
            } else {
                return -1;
            }
        }
        buf += bytes_read;
        bytes_to_read -= bytes_read;
        total_bytes_read += bytes_read;
    }
    return total_bytes_read;
}
*/ 

static void *run(hashpipe_thread_args_t * args)
{
    hpguppi_input_databuf_t *db  = (hpguppi_input_databuf_t *)args->obuf;
    hashpipe_status_t st = args->st;
    const char * status_key = args->thread_desc->skey;

    /* Main loop */
    int rv;
    int block_idx = 0;
    int block_count=0, filenum=0;
    long int raw_file_size = 0;
    long int cur_pos = 0;
    int end_of_scan = 0; // End of file flag
    int blocsize;
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
    /* Init output file descriptor (-1 means no file open) */
    static int fdin = -1;
    char *header;
    char header_buf[MAX_HDR_SIZE];
    int headersize;
    int open_flags = O_RDONLY;
    int directio = 0;
    int sim_flag = 0; // Set to 1 if you'd like to use simulated data rather than the payload from the RAW file
    char * sim_data; // Initialize simulated data array
    // 5 seconds worth of processing at a time
    // 1k mode
    //int n_chan = 1; 
    //int n_samp = 4096*1024; // 4194304; // 2^22
    // 4k mode
    int n_chan = 4; // 64
    int n_samp = 1024*1024; // 1048576; // 2^20
    // 32k mode
    //int n_chan = 32;
    //int n_samp = 128*1024; // 131072; // 2^17
    int n_pol = 2; 
    sim_data = (char *)simulate_data(n_pol, n_chan, n_samp); // Generate block of simulated data
    ssize_t read_blocsize;
#if TIMING
    float read_time = 0;
    float time_taken_r = 0;
#endif

    while (run_threads()) {
        hashpipe_status_lock_safe(&st);
        hputs(st.buf, status_key, "waiting");
        hputi4(st.buf, "NETBKOUT", block_idx);
        hashpipe_status_unlock_safe(&st);
        // Wait for data
        /* Wait for new block to be free, then clear it
         * if necessary and fill its header with new values.
         */
#if VERBOSE
	printf("RAW INPUT: while ((rv=hpguppi_input_databuf_wait_free(db, block_idx)) \n");
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
	printf("RAW INPUT: Before file open if{} \n");
#endif

        // Stride through RAW file to get all time samples for specific subband in the RAW file and place in buffer accordingly
        for(int s = 0; s<n_subband; s++){ // Shift to next subband of 16
            // Write index of subband to status buffer
            hputi4(st.buf, "SUBBAND", s);

            // At the beginning of a file so set end_of_scan == 0
            end_of_scan = 0;

            // Iterate through blocks until the end of a sequence of RAW files
            while(!end_of_scan){
                // Read raw files
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
                        //printf("RAW INPUT: Got current RAW file: %s\n", cur_fname);
                    }
                    // Now create the basefilename
                    // If a '...0000.raw' file exists, that is different from the previous '...0000.raw' file
                    if (strcmp(prev_fname, cur_fname) != 0){
                        strcpy(prev_fname, cur_fname); // Save this file name for comparison on the next iteration
                        printf("RAW INPUT: Got current RAW file: %s\n", cur_fname);
                        base_pos = strchr(cur_fname, '.'); // Finds the first occurence of a period in the filename
                        period_pos = base_pos-cur_fname;
                        printf("RAW INPUT: The last position of . is %ld \n", period_pos);
                        memcpy(basefilename, cur_fname, period_pos); // Copy base filename portion of file name to tmp_basefilename variable
                
                        printf("RAW INPUT: Base filename from command: %s \n", basefilename);

                        // Get basefilename with no path and place in status buffer
                        // strrchr() finds the last occurence of the specified character
                        char_offset = strrchr(basefilename, character);
                        slash_pos = char_offset-basefilename;
                        printf("RAW INPUT: The last position of %c is %ld \n", character, slash_pos);

                        // Get file name with no path
                        memcpy(new_base, &basefilename[slash_pos+1], sizeof(basefilename)-slash_pos);
                        printf("RAW INPUT: File name with no path: %s \n", new_base);

                        hputs(st.buf, "BASEFILE", new_base);

                    }
                    else{
                        // The RAW file hasn't changed so wait for new file to show up in the buffer
                        //printf("RAW INPUT: Waiting for new RAW file name! \n");
                        continue;
                    }
                    sprintf(fname, "%s.%04d.raw", basefilename, filenum);
            
                    printf("RAW INPUT: Opening first raw file '%s'\n", fname);
                    fdin = open(fname, open_flags, 0644);
                    if (fdin==-1) {
                        hashpipe_error(__FUNCTION__,"Error opening file.");
                        pthread_exit(NULL);
                    }

                    // Get raw file size in order to calculate the number of blocks in the file
                    raw_file_size = get_file_size(fdin);

                    headersize= get_header_size(fdin, header_buf, MAX_HDR_SIZE);

                    blocsize = get_block_size(header_buf, MAX_HDR_SIZE);
                    lseek(fdin, 0, SEEK_SET);
                }

                // Get the current position of the file
                cur_pos = lseek(fdin, (size_t)0, SEEK_CUR);

                // If we are not at the end of the file, read blocks and transfer to shared mem buffer
                // else if we are at the end, move on to the next file
                if((raw_file_size-cur_pos)>=blocsize){
                    // If a block is missing, copy a block of zeros to the buffer in it's place
                    // Otherwise, write the data from a block to the buffer
                    if(n_missed_blks == 0){
                        //Handling header - size, output path, directio----------------
#if VERBOSE
                        printf("RAW INPUT: int headersize= get_header_size(fdin, header_buf, MAX_HDR_SIZE); \n");
#endif
                        headersize= get_header_size(fdin, header_buf, MAX_HDR_SIZE);
                        set_output_path(header_buf, outdir, MAX_HDR_SIZE);
#if VERBOSE
                        printf("RAW INPUT: char *header = hpguppi_databuf_header(db, block_idx); \n");
#endif
                        if(block_count == 0){
                            header = hpguppi_databuf_header(db, block_idx);
                            hashpipe_status_lock_safe(&st);
                            hputs(st.buf, status_key, "receiving");
                            memcpy(header, &header_buf, headersize);
                            hashpipe_status_unlock_safe(&st);

                            //printf("RAW INPUT: directio = hpguppi_read_directio_mode(header); \n");
                            directio = hpguppi_read_directio_mode(header);

#if VERBOSE
                            printf("RAW INPUT: headersize = %d, directio = %d\n", headersize, directio);
#endif

                            // Adjust length for any padding required for DirectIO
                            if(directio) {
                                // Round up to next multiple of 512
                                headersize = (headersize+511) & ~511;
                            }

                            // Initialize block in buffer
                            ptr = hpguppi_databuf_data(db, block_idx);
                        }
                        lseek(fdin, headersize-MAX_HDR_SIZE, SEEK_CUR);
                        blocsize = get_block_size(header_buf, MAX_HDR_SIZE);
                        cur_pktidx = get_cur_pktidx(header_buf, MAX_HDR_SIZE);
                        pktstart = get_pktstart(header_buf, MAX_HDR_SIZE);
                        pktstop = get_pktstop(header_buf, MAX_HDR_SIZE);
                        piperblk = get_piperblk(header_buf, MAX_HDR_SIZE);
                        // Descriptions of the variables below are at the initializations
                        nants = get_nants(header_buf, MAX_HDR_SIZE);
                        npol = get_npol(header_buf, MAX_HDR_SIZE);
                        obsnchan = get_obsnchan(header_buf, MAX_HDR_SIZE);
                        n_coarse = (int)obsnchan/nants;
                        n_coarse_proc = n_coarse/n_subband;
                        n_samp_per_block = (int)(blocsize/(2*obsnchan*npol));

                        printf("RAW INPUT: Current pktidx = %ld, pktstart = %ld, and pktstop = %ld \n", cur_pktidx, pktstart, pktstop);
                        printf("RAW INPUT: piperblk = %d \n", piperblk);

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
                        printf("RAW INPUT: Block size: %d, and  BLOCK_DATA_SIZE: %d \n", blocsize, BLOCK_DATA_SIZE);
                        printf("RAW INPUT: header size: %d, and  MAX_HDR_SIZE: %d \n", headersize, MAX_HDR_SIZE);
#endif
                        for(int a = 0; a<nants; a++){
                            for(int c = 0; c<n_coarse_proc; c++){
                                // Offset by antenna and coarse channel index taking the subband index into account as well
                                lseek(fdin, data_blk_in_idx(0, 0, (c + n_coarse_proc*s), a, npol, n_samp_per_block, n_coarse_proc), SEEK_CUR);
                                // Place input data in shared memory buffer (from RAW file or simulated)
                                if(sim_flag == 0){
                                    // Read the remaining dimensions in the RAW block (npol and n_samp_per_block) and place in the shared memory buffer block
                                    //read_blocsize = read(fdin, &ptr[data_shm_blk_idx(0,0,block_count,c,a,npol,n_samp_per_block,n_blocks,n_coarse_proc)], npol*n_samp_per_block);
                                    read_blocsize = read(fdin, &ptr[data_shm_blk_idx(0,c,a,0,block_count,npol,n_coarse_proc,nants,n_samp_per_block)], npol*n_samp_per_block);
                                } else{
                                    // Copy simulated data to shared memory buffer block - data_in_idx(p, c, a, t, w, Np, Nc, Na, Nt)
                                    //memcpy(&ptr[data_in_idx(0,0,0,c,a,n_pol,n_samp,1,n_chan)], &sim_data[data_in_idx(0,0,0,c,a,n_pol,n_samp,1,n_chan)], n_pol*n_samp);
                                    memcpy(&ptr[data_in_idx(0,c,a,0,0,n_pol,n_chan,N_ANT,1)], &sim_data[data_in_idx(0,c,a,0,0,n_pol,n_chan,N_ANT,1)], n_pol*n_samp);
                                }
                            }
                        }
#if TIMING
                        // Stop timing read
                        clock_gettime(CLOCK_MONOTONIC, &tval_after);
                        time_taken_r = (float)(tval_after.tv_sec - tval_before.tv_sec); //*1e6; // Time in seconds since epoch
                        time_taken_r = time_taken_r + (float)(tval_after.tv_nsec - tval_before.tv_nsec)*1e-9; //*1e-6; // Time in nanoseconds since 'tv_sec - start and end'
                        read_time = time_taken_r;

                        printf("RAW INPUT: Time taken to read from RAW file = %f ms \n", read_time);
#endif
                        // Iterate through blocks in RAW files for a scan corresponding to a subband
                        block_count++;
                    } else if(n_missed_blks > 0){
                        if(block_count == 0){
                            headersize= get_header_size(fdin, header_buf, MAX_HDR_SIZE);
                            piperblk = get_piperblk(header_buf, MAX_HDR_SIZE);

                            // Copy the same header info from previous block to zero block in buffer
                            header = hpguppi_databuf_header(db, block_idx);
                            hashpipe_status_lock_safe(&st);
                            hputs(st.buf, status_key, "receiving");
                            memcpy(header, &header_buf, headersize);
                            hashpipe_status_unlock_safe(&st);

                            //printf("RAW INPUT: directio = hpguppi_read_directio_mode(header); \n");
                            directio = hpguppi_read_directio_mode(header);

#if VERBOSE
                            printf("RAW INPUT: headersize = %d, directio = %d\n", headersize, directio);
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
                        set_pktidx(header_buf, zero_blk_pktidx, MAX_HDR_SIZE);

                        // Copy block of zeros to block in buffer
                        memcpy(&ptr[block_count*npol*n_coarse_proc*nants*n_samp_per_block], zero_blk, npol*n_coarse_proc*nants*n_samp_per_block);
            
                        // Decrement n_missed_blks by 1
                        n_missed_blks -= 1;

                        // If n_missed_blks == 0, then process the data at the packet index after blocks of zeros (cur_pktidx)
                        if(n_missed_blks == 0){

                            for(int a = 0; a<nants; a++){
                                for(int c = 0; c<n_coarse_proc; c++){
                                    // Offset by antenna and coarse channel index taking the subband index into account as well
                                    lseek(fdin, data_blk_in_idx(0, 0, (c + n_coarse_proc*s), a, npol, n_samp_per_block, n_coarse_proc), SEEK_CUR);
                                    // Place input data in shared memory buffer from RAW file
                                    // Read the remaining dimensions in the RAW block (npol and n_samp_per_block) and place in the shared memory buffer block
                                    read_blocsize = read(fdin, &ptr[data_shm_blk_idx(0,c,a,0,block_count,npol,n_coarse_proc,nants,n_samp_per_block)], npol*n_samp_per_block);
                                }
                            }
                            if(block_count >= 0 && (block_count <= 5)){
                                printf("RAW INPUT: Number of bytes read in read(): %zd \n", read_blocsize);
                            }
#if VERBOSE
                            printf("RAW INPUT: First element of buffer: %d \n", ptr[0]);
#endif
                            // Set previous PKTIDX
                            prev_pktidx = cur_pktidx;
                        }

                        // Increment the block count even with zero blocks
                        block_count++;
                    }
                }// Sometimes there might be headers with no data blocks remaining in the RAW file. So copy zero blocks to buffer blocks with the appropriate header
                else if((cur_pos+blocsize) > raw_file_size){
                    // If we reach the last header with no data block,
                    // then move on to the next RAW file or wait for one
                    if(cur_pos >= (raw_file_size-headersize)){
                        close(fdin);
                        filenum++;
            
                        sprintf(fname, "%s.%4.4d.raw", basefilename, filenum);
                        printf("RAW INPUT: Opening next raw file '%s'\n", fname);
                        fdin = open(fname, open_flags, 0644);
                        if (fdin==-1) { // End of a sequence of RAW files corresponding to a scan
                            filenum=0;

                            // If PKTIDX < PKTSTOP and we're at the end of the scan then set PKTIDX == PKTSTOP 
                            // and copy a dummy block to the buffer
                            if(cur_pktidx < pktstop){
                                // Send dummy block with PKTIDX set to PKTSTOP
                                set_pktidx(header_buf, pktstop, MAX_HDR_SIZE);
                            }
                            // Reset previous pktidx to 0 for the next scan
                            prev_pktidx = 0;

                            // End of scan so move on to the next subband
                            end_of_scan = 1;

                            // Reset block_count to 0 (the block index of the RAW file)
                            block_count=0;

                            // Mark block as full
                            hpguppi_input_databuf_set_filled(db, block_idx);

                            // Setup for next block
                            block_idx = (block_idx + 1) % N_INPUT_BLOCKS;
                        }
                    }

                    // Write blocks of zeros to the shared memory buffer blocks
                    // Copy block of zeros to block in buffer
                    memcpy(&ptr[block_count*npol*n_coarse_proc*nants*n_samp_per_block], zero_blk, npol*n_coarse_proc*nants*n_samp_per_block);

                    // If we haven't reached the last header, then increment block_count (the index of blocks in the RAW files of a subband)
                    if(cur_pos < (raw_file_size-headersize)){
                        block_count++;
                    }
                }// If (cur_pos > (raw_file_size-blocsize)) && (cur_pos <= raw_file_size),
                // then we have reached the end of a file so move on to next file or wait for new file and set PKTIDX == PKTSTOP if necessary
                else if((cur_pos > (raw_file_size-blocsize)) && (cur_pos <= raw_file_size)){
                    close(fdin);
                    filenum++;
            
                    sprintf(fname, "%s.%4.4d.raw", basefilename, filenum);
                    printf("RAW INPUT: Opening next raw file '%s'\n", fname);
                    fdin = open(fname, open_flags, 0644);
                    if ((fdin==-1) && (s == (n_subband-1))) { // End of a sequence of RAW files corresponding to a scan and end of subbands
                        filenum=0;

                        // Inform the downstream thread that we have reached the end of a scan
                        if(cur_pktidx < pktstop){
                            // Mark block as full
                            hpguppi_input_databuf_set_filled(db, block_idx);

                            // Setup for next block
                            block_idx = (block_idx + 1) % N_INPUT_BLOCKS;

                            // Send dummy block with PKTIDX set to PKTSTOP
                            set_pktidx(header_buf, pktstop, MAX_HDR_SIZE);

                            // Create a header for a dummy block
                            header = hpguppi_databuf_header(db, block_idx);
                            hashpipe_status_lock_safe(&st);
                            hputs(st.buf, status_key, "receiving");
                            memcpy(header, &header_buf, headersize);
                            hashpipe_status_unlock_safe(&st);

                            // Initialize block
                            ptr = hpguppi_databuf_data(db, block_idx);

                            // Copy block of zeros to block in buffer
                            memcpy(ptr, zero_blk, N_INPUT);
                        }

                        // Mark block as full
                        hpguppi_input_databuf_set_filled(db, block_idx);

                        // Setup for next block
                        block_idx = (block_idx + 1) % N_INPUT_BLOCKS;

                        // Reset block_count to 0 (the index of blocks in the RAW files of a subband)
                        block_count=0;

                        // Reset previous pktidx to 0 for the next scan
                        prev_pktidx = 0;

                        // End of scan so move on to the next subband
                        end_of_scan = 1;
                    }else if ((fdin==-1) && (s < n_subband)) { // End of a sequence of RAW files corresponding to a scan and less than the no. of subbands
                        // File number is reset to 0
                        filenum=0;

                        sprintf(fname, "%s.%4.4d.raw", basefilename, filenum);
                        printf("RAW INPUT: Opening first raw file '%s' of scan for subband = %d of %d \n", fname, s, n_subband);
                        fdin = open(fname, open_flags, 0644);

                        // Start from the beginning of scan (position 0 of filenum=0), but process next subband
                        lseek(fdin, 0, SEEK_SET); // Probably unnecessary, but doesn't hurt

                        // Mark block as full
                        hpguppi_input_databuf_set_filled(db, block_idx);

                        // Setup for next block
                        block_idx = (block_idx + 1) % N_INPUT_BLOCKS;

                        // Reset block_count to 0 (the index of blocks in the RAW files of a subband)
                        block_count=0;

                        // Reset previous pktidx to 0 for the next scan
                        prev_pktidx = 0;

                        // End of scan so move on to the next subband
                        end_of_scan = 1;
                    }
                }

                /* Will exit if thread has been cancelled */
                pthread_testcancel();
            }
        }
    }

    // Thread success!
    if (fdin!=-1) {
	printf("RAW INPUT: Closing file! \n");
        close(fdin);
    }
    return THREAD_OK;

}

static hashpipe_thread_desc_t hpguppi_rawfile_input_thread = {
    name: "hpguppi_rawfile_input_thread",
    skey: "NETSTAT",
    init: NULL,
    run:  run,
    ibuf_desc: {NULL},
    obuf_desc: {hpguppi_input_databuf_create}
};

static __attribute__((constructor)) void ctor()
{
    register_hashpipe_thread(&hpguppi_rawfile_input_thread);
}
