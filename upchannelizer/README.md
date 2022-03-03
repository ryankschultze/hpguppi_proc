Upchannelizer library

To compile the standalone code, make sure the main() function in the .cu script is uncommented and use this command:
- nvcc -o upchannelizer.exe -arch=sm_86 upchannelizer.cu -lcufft

And to run the executable:
- ./upchannelizer.exe
