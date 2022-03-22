Upchannelizer library

To compile the standalone code, make sure the main() function in the .cu script is uncommented and use this command:
- nvcc -o upchannelizer.exe -arch=sm_86 upchannelizer.cu -lcufft

And to run the executable:
- ./upchannelizer.exe

The post-processing scripts in the post_processing directory are used for analysis and verification.

The plot_cufft_output_binary.py reads a raw binary file with the output of the data that doesn't have an FFT shift and it applies one to the data.

The plot_cufft_output_with_fftshift.py reads a raw binary file with the output of the data that has an FFT shift applied and generates plots assuming an FFT shift was applied.

And the fft_verification.py script reads a raw binary files with the simulated input data that is also going to processed by the cuFFT and performs an FFT with numpy.fft for comparison with the output from the upchannelizer code.
