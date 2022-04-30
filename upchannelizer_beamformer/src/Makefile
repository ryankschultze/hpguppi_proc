-include Makefile.local

prefix ?= /usr/local
libdir ?= $(prefix)/lib
includedir ?= $(prefix)/include
bindir ?= $(prefix)/bin
CUDA_DIR ?= $(prefix)/cuda
ifneq ($(strip $(OSTYPE)),osx)
CUDA_LIBDIR ?= $(CUDA_DIR)/lib64
else
CUDA_LIBDIR ?= $(CUDA_DIR)/lib
endif

#Don't install in /usr/local yet
#Path to repo lib and include
repo_prefix ?= ../..
bflibdir ?= $(repo_prefix)/lib
bfincludedir ?= $(repo_prefix)/include

#On blc03 with NVIDIA RTX 2080
CUDA_ARCH ?= "-gencode=arch=compute_75,code=sm_75"
#On blc03 with NVIDIA A5000
#CUDA_ARCH ?= "-gencode=arch=compute_80,code=sm_80"
#On blc02 with NVIDIA Xp
#CUDA_ARCH ?= "-gencode=arch=compute_53,code=sm_53"

NVCC = $(CUDA_DIR)/bin/nvcc

CU_LIB = -lcufft -lcudart

LINKER = -Xlinker -export-dynamic

libupchannelizer_beamformer.so: upchannelizer_beamformer.o
	$(NVCC) $(CUDA_ARCH) -Xcompiler -fPIC -Xcompiler -D_REENTRANT -O3 upchannelizer_beamformer.o -o libupchannelizer_beamformer.so -L$(CUDA_LIBDIR) $(CU_LIB) $(LINKER) --shared

upchannelizer_beamformer.o: upchannelizer_beamformer.cu upchannelizer_beamformer.h
	$(NVCC) $(CUDA_ARCH) $(CU_LIB) -c -o upchannelizer_beamformer.o -Xcompiler -fPIC -Xcompiler -D_REENTRANT -O3 upchannelizer_beamformer.cu
 
clean:
	rm -f *.o
	rm -f *.so

install: upchannelizer_beamformer.h libupchannelizer_beamformer.so
	mkdir -p $(bfincludedir)
	cp upchannelizer_beamformer.h $(bfincludedir)
	mkdir -p $(bflibdir)
	cp libupchannelizer_beamformer.so $(bflibdir)

