-include Makefile.local

prefix ?= /usr/local
libdir ?= $(prefix)/lib
includedir ?= $(prefix)/include
bindir ?= $(prefix)/bin
CUDA_DIR ?= $(prefix)/NVIDIA/cuda
ifneq ($(strip $(OSTYPE)),osx)
CUDA_LIBDIR ?= $(CUDA_DIR)/lib64
else
CUDA_LIBDIR ?= $(CUDA_DIR)/lib
endif

#Don't install in /usr/local yet
bflibdir ?= /home/mruzinda/hpguppi_proc/coherent_beamformer/lib
bfincludedir ?= /home/mruzinda/hpguppi_proc/coherent_beamformer/include

#On blc03 with NVIDIA A5000
CUDA_ARCH ?= "-gencode=arch=compute_80,code=sm_80"
#On blc02 with NVIDIA Xp
#CUDA_ARCH ?= "-gencode=arch=compute_53,code=sm_53"

NVCC = $(CUDA_DIR)/bin/nvcc

LINKER = -Xlinker -export-dynamic

libcoherent_beamformer.so: coherent_beamformer.o
	$(NVCC) $(CUDA_ARCH) -Xcompiler -fPIC -Xcompiler -D_REENTRANT -O3 coherent_beamformer.o -o libcoherent_beamformer.so -L$(CUDA_LIBDIR) $(LINKER) --shared

coherent_beamformer.o: coherent_beamformer_char_in.cu coherent_beamformer_char_in.h
	$(NVCC) $(CUDA_ARCH) -c -o coherent_beamformer.o -Xcompiler -fPIC -Xcompiler -D_REENTRANT -O3 coherent_beamformer_char_in.cu

#libflagbeamformer.so: cublas_beamformer.o
#	$(NVCC) -Xcompiler -fPIC -Xcompiler -D_REENTRANT -O3 cublas_beamformer.o -o libflagbeamformer.so -L$(CUDA_LIBDIR) $(LINKER) --shared

#cublas_main.o: cublas_main.cu cublas_beamformer.h
#	$(NVCC) $(CUDA_ARCH) -I $(CUDA_DIR)/include -c -o cublas_main.o cublas_main.cu
 
clean:
	rm -f *.o
	rm -f *.so
	rm -f coherent_beamformer

install: coherent_beamformer_char_in.h libcoherent_beamformer.so
	mkdir -p $(bfincludedir)
	cp coherent_beamformer_char_in.h $(bfincludedir)
	mkdir -p $(bflibdir)
	cp libcoherent_beamformer.so $(bflibdir)

