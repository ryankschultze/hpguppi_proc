Coherent beamformer library

To compile, create shared objects, and install the libraries, and header files:
- In the src/ directory:
  - make
  - make install

To compile the standalone code, make sure the main() function in the .cu script is uncommented and use this command:
- nvcc -o coherent_beamformer_char_in.exe -arch=sm_86 coherent_beamformer_char_in.cu

And to run the executable:
- ./coherent_beamformer_char_in.exe
