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

#include "coherent_beamformer_char_in.h"

#define data_blk_in_idx(p, t, f, a, Np, Nt, Nf)              ((p) + (Np)*(t) + (Nt)*(Np)*(f) + (Nf)*(Nt)*(Np)*(a))
#define data_shm_blk_idx(p, t, b, f, a, Np, Nt, Nb, Nf)      ((p) + (Np)*(t) + (Nt)*(Np)*(b) + (Nb)*(Nt)*(Np)*(f) + (Nf)*(Nb)*(Nt)*(Np)*(a))
#define data_in_idx(p, t, w, c, a, Np, Nt, Nw, Nc)           ((p) + (Np)*(t) + (Nt)*(Np)*(w) + (Nw)*(Nt)*(Np)*(c) + (Nc)*(Nw)*(Nt)*(Np)*(a))
#define VERBOSE 0
#define TIMING 0

int get_header_size(int fdin, char * header_buf, size_t len)
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
long int get_file_size(int fdin){
    off_t cur_pos = lseek(fdin, (size_t)0, SEEK_CUR);
    off_t file_size = lseek(fdin, (size_t)0, SEEK_END);
    lseek(fdin, cur_pos, SEEK_SET);
    return file_size;
}

int get_block_size(char * header_buf, size_t len)
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

int get_nants(char * header_buf, size_t len)
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

int get_npol(char * header_buf, size_t len)
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

int get_obsnchan(char * header_buf, size_t len)
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

int64_t get_cur_pktidx(char * header_buf, size_t len)
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

int64_t get_nxt_pktidx(int fdin, int blocsize, char * header_buf, size_t len)
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

int get_piperblk(char * header_buf, size_t len)
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

void set_output_path(char * header_buf, char * outdir, size_t len)
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

void set_blocksize(char * header_buf, int blocsize, size_t len)
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

ssize_t read_fully(int fd, void * buf, size_t bytes_to_read)
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

static void *run(hashpipe_thread_args_t * args)
{
    hpguppi_input_databuf_t *db  = (hpguppi_input_databuf_t *)args->obuf;
    hashpipe_status_t st = args->st;
    const char * status_key = args->thread_desc->skey;

    /* Main loop */
    int rv;
    int block_idx = 0;
    int block_count=0, filenum=0;
    int blocsize;
    int nants = 0;            // Number of antennas stated in RAW file
    int npol = 0;             // Number of polarizations stated in RAW file
    int obsnchan = 0;         // Number of coarse channels X antennas stated in RAW file
    int n_blocks = 0;         // Number of blocks in a RAW file
    int n_coarse = 0;         // Number of coarse channels stated in RAW file
    int n_subband = 16;       // Number of sub bands processed serially
    int n_coarse_proc = 0;    // Number of coarse channels processed at one time
    int n_samp_per_block = 0; // Number of time sample in one block
    int n_samp = 0;           // Number of time samples in a RAW file
    int n_fft = 0;            // Number of FFT points depending on 1k, 4k, or 32k mode
    int n_windows = 0;        // Number of short time integration windows (time samples after FFT)
    int piperblk;
    int64_t cur_pktidx;
    int64_t nxt_pktidx;
    int n_missed_blks = 0;
    int zero_blk_idx = 0;
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
    long int raw_file_size = 0;
    char *header;
    char header_buf[MAX_HDR_SIZE];
    int headersize;
    int open_flags = O_RDONLY;
    int directio = 0;
    int sim_flag = 0; // Set to 1 if you'd like to use simulated data rather than the payload from the RAW file
    char * sim_data; // Initialize simulated data array
    int n_chan = 64; // Value set in coherent_beamformer_char_in.h and used with simulated data when no RAW file is read
    int n_samp = 8192;
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

        //Read raw files
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
		sleep(0.1);
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
        }

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

            //Read data--------------------------------------------------
	    // Start timing read
#if TIMING
            struct timespec tval_before, tval_after;
            clock_gettime(CLOCK_MONOTONIC, &tval_before);
#endif
            // Initialize block in buffer
            ptr = hpguppi_databuf_data(db, block_idx);
            lseek(fdin, headersize-MAX_HDR_SIZE, SEEK_CUR);
            blocsize = get_block_size(header_buf, MAX_HDR_SIZE);
            // Descriptions of the variables below are at the initializations
            nants = get_nants(header_buf, MAX_HDR_SIZE);
            npol = get_npol(header_buf, MAX_HDR_SIZE);
            obsnchan = get_obsnchan(header_buf, MAX_HDR_SIZE);
            n_blocks = raw_file_size/(headersize + blocsize);
            n_coarse = (int)obsnchan/nants;
            // Number of FFT points in 1k, 4k, or 32k mode
            if(n_coarse == 16){ // 1k mode
                n_fft = 524288;
            }
            if(n_coarse == 64){ // 4k mode
                n_fft = 131072;
            }
            if(n_coarse == 512){ // 32k mode
                n_fft = 16384;
            }
            n_coarse_proc = n_coarse/n_subband;
            n_samp_per_block = (int)(blocsize/(2*obsnchan*npol));
            n_samp = n_samp_per_block*n_blocks;
            n_windows = n_samp/n_fft;

            // Set block size in header of the block in the buffer
            // If this is not done, it will set the block size of in the header of the RAW file which is never 0
            //set_blocksize(header, blocsize, MAX_HDR_SIZE);

            //hashpipe_status_lock_safe(&st);
            //hputu8(st.buf, "RBLKSIZE", (uint64_t)blocsize);
            //hashpipe_status_unlock_safe(&st);

            // Can only check next block if there IS one so make sure one exists
            if(blocsize > 0){
                piperblk = get_piperblk(header_buf, MAX_HDR_SIZE);
                cur_pktidx = get_cur_pktidx(header_buf, MAX_HDR_SIZE);
                nxt_pktidx = get_nxt_pktidx(fdin, blocsize, header_buf, MAX_HDR_SIZE);

                printf("RAW INPUT: Current pktidx = %ld and next pktidx = %ld \n", cur_pktidx, nxt_pktidx);
                
                // Check to see whether there will be any following blocks that are missed
                // If the next packet index is greater than the current pkt idx + piperblk, then N blocks were missed
                // So replace with zero blocks
                if((nxt_pktidx > 0) && (nxt_pktidx > (cur_pktidx+piperblk))){
                    // Assuming piperblk is the number of packets per block
                    // And PKTIDX is the index of the first packet in a block
                    n_missed_blks = (nxt_pktidx - cur_pktidx)/piperblk;
                }
            }// This (else if) means that we are at the end of the RAW file so advance to the next one to put the appropriate block of data in the buffer
            else if(blocsize <= 0){
                close(fdin);
                filenum++;
            
                sprintf(fname, "%s.%4.4d.raw", basefilename, filenum);
                printf("RAW INPUT: Opening next raw file '%s'\n", fname);
                fdin = open(fname, open_flags, 0644);
                if (fdin==-1) {
                    filenum=0;
                }
            
                block_count=0;
                continue;
            }

#if VERBOSE
            printf("RAW INPUT: Block size: %d, and  BLOCK_DATA_SIZE: %d \n", blocsize, BLOCK_DATA_SIZE);
            printf("RAW INPUT: header size: %d, and  MAX_HDR_SIZE: %d \n", headersize, MAX_HDR_SIZE);
#endif
	    if(sim_flag == 0){
                // Stride through RAW file to get all time samples for specific subband in the RAW file and place in buffer accordingly
                for(int s = 0; s<n_subband; s++){ // Shift to next subband of 16
                    for(int b = 0; b<n_blocks; b++){
                        // Offset to next block - Make sure header size is taken into account
                        // The position is already offset by the header size of the first block above
                        lseek(fdin, b*(blocsize + headersize), SEEK_CUR);
                        for(int a = 0; a<nants; a++){
                            for(int c = 0; c<n_coarse_proc; c++){
                                // Offset by antenna and coarse channel index taking the subband index into account as well
                                lseek(fdin, data_blk_in_idx(0, 0, (c + n_coarse_proc*s), a, npol, n_samp_per_block, n_coarse_proc), SEEK_CUR);
                                // Read the remaining dimensions in the RAW block (npol and n_samp) and place in the shared memory buffer block
                                read_blocsize = read(fdin, ptr + data_shm_blk_idx(0, 0, b, c, a, npol, n_samp_per_block, n_blocks, n_coarse_proc), npol*n_samp_per_block);
                            }
                        }
                        // Set position in file back to beginning after first block's header
                        lseek(fdin, headersize, SEEK_SET);
                    }
                    // Mark shared memory buffer block as full
                    hpguppi_input_databuf_set_filled(db, block_idx);

                    // Setup for next shared memory buffer block
                    block_idx = (block_idx + 1) % N_INPUT_BLOCKS;
                    block_count++;
                    if(block_count == 0){
                        printf("RAW INPUT: Number of bytes read in read(): %zd \n", read_blocsize);
                    }
                }
#if VERBOSE
                printf("RAW INPUT: First element of buffer: %d \n", ptr[0]);
#endif
            } else{
	        memcpy(ptr, sim_data, N_INPUT);
	        //ptr += N_INPUT;
            }
#if TIMING
	    // Stop timing read
            clock_gettime(CLOCK_MONOTONIC, &tval_after);
            time_taken_r = (float)(tval_after.tv_sec - tval_before.tv_sec); //*1e6; // Time in seconds since epoch
            time_taken_r = time_taken_r + (float)(tval_after.tv_nsec - tval_before.tv_nsec)*1e-9; //*1e-6; // Time in nanoseconds since 'tv_sec - start and end'
            read_time = time_taken_r;

            printf("RAW INPUT: Time taken to read from RAW file = %f ms \n", read_time);
#endif
        } else if(n_missed_blks > 0){
            // Copy the same header info from previous block to zero block in buffer
            header = hpguppi_databuf_header(db, block_idx);
            hashpipe_status_lock_safe(&st);
            hputs(st.buf, status_key, "receiving");
            memcpy(header, &header_buf, headersize);
            hashpipe_status_unlock_safe(&st);

            // Copy block of zeros to block in buffer
            ptr = hpguppi_databuf_data(db, block_idx);
            memcpy(ptr, zero_blk, N_INPUT);
            
            // Increment the blocks by one each time
            // And check to see that the number of missed blocks has been reached
            zero_blk_idx += 1;
            if(zero_blk_idx == n_missed_blks){
                // Reset the number of missed blocks and block index to 0
                n_missed_blks = 0;
                zero_blk_idx = 0;
            }
        }

        // Mark block as full
        hpguppi_input_databuf_set_filled(db, block_idx);

        // Setup for next block
        block_idx = (block_idx + 1) % N_INPUT_BLOCKS;
        block_count++;

        /* Will exit if thread has been cancelled */
        pthread_testcancel();

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
