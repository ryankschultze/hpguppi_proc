# hpguppi_proc

There are a few different sub-directories in this repo. Here is a brief description of each sub-directory:
- **coherent_beamformer** - Contains the code used to generate and test the beamformer library without upchannelization. 
- **include** - Contains the header files for the different libraries in use by the hpguppi\_proc threads.
- **init** - Contains the shell script used to run the hpguppi\_proc plugins/modes e.g. the beamformer with and without the upchannelization.
- **post_processing** - Contains scripts used for general analysis of the output from the hpguppi\_proc threads.
- **src** - Contains the hpguppi\_proc thread code and scripts used for reading, processing, and writing.
- **upchannelizer_beamformer** - Contains the code used to generate and test the beamformer library with upchannelization.

## Configuring and compiling hpguppi_proc threads

To configure and compile the first time:
```
 $ cd src
 $ libtoolize
 $ autoconf
 $ autoreconf -is
 $ ./configure --with-libsla=/usr/local/listen/lib --with-libcoherent\_beamformer=../lib --with-libupchannelizer\_beamformer=../lib
 $ make
 ```

To configure and compile any other time after:
```
 $ cd src
 $ autoreconf -is
 $ ./configure --with-libsla=/usr/local/listen/lib --with-libcoherent\_beamformer=../lib --with-libupchannelizer\_beamformer=../lib
 $ make
 ```

If everything is configured as it should be, you can just:
```
 $ make
 ```

NOTE: libsla may be in a different directory depending on the machine.

If python library is required:

```
 $ cd src
 $ autoreconf -is
 $ ./configure --with-libsla=/usr/local/listen/lib --with-libcoherent\_beamformer=../lib --with-libupchan\_beamformer=../lib --with-libpython3.7m=/opt/conda/lib
 $ make
 ```

## Compiling and Running the Beamformer with and without Upchannelization

The coherent beamformer can be compiled to generate an executable that is standalone (for testing) as well as a shared object library for a processing thread to link (e.g. in the hpguppi\_proc repo - [link](https://github.com/UCBerkeleySETI/hpguppi_proc)). 

Before doing any of the following steps, ensure that you have the correct paths setup for **nvcc**, **hashpipe**, **rawspec**, and any other functions or libraries that you may need.

### Standalone version:

To compile the standalone code, uncomment the main() in the .cu script and use the command below. Or the user can write their own script containing a main() function that calls the coherent beamformer functions and kernels. Also be sure to use the SM architecture that corresponds to the GPU in use e.g sm\_86 for NVIDIA A5000.

This link helps show which SM architectures correspond to a particular GPU: [link](https://arnon.dk/matching-sm-architectures-arch-and-gencode-for-various-nvidia-cards/)
There are two beamformers, one without upchannelization, and one with upchannelization. The beamformer without upchannelization can be found in the coherent\_beamformer directory in hpguppi\_proc, and can be compiled with the following command:
```
 $ nvcc -o coherent\_beamformer\_char\_in.exe -arch=sm\_86 coherent\_beamformer\_char\_in.cu
```
And to run the executable:
```
 $ ./coherent\_beamformer\_char\_in.exe
```
 
The beamformer with upchannelization can be found in the upchannelizer\_beamformer directory in hpguppi_proc, and can be compiled with the following command:
```
 $ nvcc -o upchannelizer\_beamformer.exe -arch=sm\_86 upchannelizer\_beamformer.cu -lcufft
```

And to run the executable:
```
 $ ./upchannelizer\_beamformer.exe
```

### Shared Object library:


To compile, create shared objects, and install the libraries, and header files:
```
 $ make
 $ make install
```

Ensure that the SM architecture is correct in the Makefile, and comment out any main() functions that may be in the CUDA script being compiled.

## Beamformer without Upchannelization

Processing with HASHPIPE:

Currently, the coherent beamformer is a plug-in to hashpipe (an application that helps setup threads and shared memory buffers between them) with 2 threads; the hpguppi\_rawfile\_input\_thread, and the hpguppi\_coherent\_bf\_thread. More information about hashpipe can be found here: [link](https://casper.astro.berkeley.edu/wiki/HASHPIPE)

hpguppi\_rawfile\_input\_thread reads GUPPI RAW files, and places blocks of data into shared memory buffers that provide data to the processing thread which in this case is hpguppi\_coherent\_bf\_thread.

hpguppi\_coherent\_bf\_thread calls the functions in the coherent beamformer library and processes the data given phase solutions and delays from an HDF5 file, called a beamformer recipe file, that corresponds to the GUPPI RAW files from a recording.

Before compiling this code, ensure that the coherent beamformer library has been generated, and installed in hpguppi\_proc/lib and the header in hpguppi\_proc/include. This code is compiled in hpguppi\_proc/src by running the make command. 

**READ THE NEXT FEW PARAGRAPHS BEFORE RUNNING THE NEXT COMMAND:** 
After the code is compiled, a shell script called readraw\_init.sh is run with the following command to process data from the GUPPI raw files of a recording:
```
 $ ./readraw\_init.sh cbf
```

where cbf is the mode that will run, and stands for coherent beamformer.

When this script is run, it will wait until there is a GUPPI RAW file name along with it's appropriate directory in the status buffer. It will also be expecting the filterbank (output) and beamformer recipe file paths in order to generate appropriate coefficients.

However, for the coherent beamformer to run smoothly, before the GUPPI RAW file name is placed in the status buffer, all the other required parameters (GUPPI RAW file path, filterbank, and beamformer recipe file paths) should be placed there first.

**NOTE**: The coherent beamformer performs computation with the dimensions set in the GUPPI RAW, and beamformer recipe files.

The paths and file names are placed in the buffer with the following command:
```
hashpipe\_check\_status -k \<keyword\> -s \<keyword string\>
```
The keywords are INPUTDIR, RAWFILE, BFR5DIR, and OUTDIR which are the GUPPI RAW file path, GUPPI RAW file name, beamformer recipe file path, and the filterbank file path respectively.

As previously mentioned, INPUTDIR, BFR5DIR, and OUTDIR must be set before RAWFILE since the coherent beamformer is waiting for that name in the status buffer. These keywords can also all be set before running the readraw\_init.sh script when performing tests. 

The following are examples of how to set the keywords given paths that were being used for testing at the time. 

To set the GUPPI RAW file path:
```
hashpipe\_check\_status -k INPUTDIR -s /mydatag/20220120/0024/Unknown/GUPPI
```

To set the beamformer recipe file path:
```
hashpipe\_check\_status -k BFRDIR -s /home/obs/20220120/0024
```

To set the filterbank (output) file path:
```
hashpipe\_check\_status -k OUTDIR -s /mydatag/20220120/0024
```

And to set the GUPPI RAW file name:
```
hashpipe\_check\_status -k RAWFILE -s guppi...0000.raw
```

Adding the -I \<instance number\> option, e.g. -I 1 for instance number 1, should set up these keywords in the status buffer of the specified instance number.

If the user is doing some testing, the readraw_init.sh script can also setup the parameters above with the following command:
```
./readraw\_init.sh cbf INPUTDIR BFRDIR OUTDIR RAWFILE 
```

For example (the additional spaces are to make things clearer for the reader, but are unnecessary in practice):
```
./readraw\_init.sh cbf /mydatag/20220120/0024/Unknown/GUPPI /home/obs/20220120/0024 /mydatag/20220120/0024 guppi...0000.raw
```

Or with the hashpipe command (this command is one line, but was separated to make it a little more clear, the "\" character is continuing the command on a new line):
```
hashpipe -p hpguppi\_proc hpguppi\_rawfile\_input\_thread hpguppi\_coherent\_bf\_thread \
-o INPUTDIR=/mydatag/20220120/0024/Unknown/GUPPI \
-o BFRDIR=/home/obs/20220120/0024 \
-o OUTDIR=/mydatag/20220120/0024 \
-o RAWFILE=guppi...0000.raw
```

After setting the file name and paths, the program should run until all of the files of the same basefilename have been read, then it will wait for a new RAW file name to be placed in the status buffer. The program will run until it receives a signal interrupt or an error occurs. This is designed to work for the ping pong approach discussed on page 3 here: [link](https://docs.google.com/document/d/1uj7vAF1FXq7kQcGdi2lr7K2eg98MFW3d3eqsAB2Z3LQ/edit#heading=h.twuqnlahbx18)

This coherent beamformer currently generates filterbank files with power rather than raw voltage, summed polarizations, and 8 time sample integration windows. The order of the dimensions from fastest to slowest are coarse frequency channel, and short time integration window. After the program has finished processing the data and has written it to the filterbank files, the user can enter a signal interrupt, Ctrl+C, to stop the program.

## Beamformer WITH Upchannelization

The process is pretty much exactly the same as it is above except with different script names and a different mode name for the shell script. The inner workings of the code are also different, and a brief summary of the threads will be provided.

The upchannelizer beamformer is also a plug-in to hashpipe (an application that helps setup threads and shared memory buffers between them) with 2 threads; the hpguppi\_stride\_input\_thread, and the hpguppi\_upchan\_bf\_thread. 

hpguppi\_stride\_input\_thread reads GUPPI RAW files, strides through each block in the file and takes a chunk corresponding to a 16th or the bandwidth at the compute node and places that data into shared memory buffers that provide it to the processing thread which in this case is hpguppi\_upchan\_bf\_thread.

hpguppi\_upchan\_bf\_thread calls the functions in the upchannelizer beamformer library and processes the data given phase solutions and delays from an HDF5 file, called a beamformer recipe file, that corresponds to the GUPPI RAW files from a recording.

Before compiling this code, ensure that the coherent beamformer library has been generated, and installed in hpguppi\_proc/lib and the header in hpguppi\_proc/include. This code is compiled in hpguppi\_proc/src by running the make command. 

The setup for running everything is pretty much the same as it is above except the threads are the ones discussed in this section and rather than cbf, the mode placed in the argument of readraw\_init.sh is ubf as shown below:
```
./readraw\_init.sh ubf
```

The output is written to filterbank files which for now can be analyzed with plot\_waterfall\_blimpy.py in the post\_processing sub-directory which uses blimpy to generate a waterfall plot. The script is run with the filterbank file as the first argument.

In the near future, the beamformer will be optimized to increase the computation speed. It's still a little slow.

