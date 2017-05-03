
-include ../common.mk
-include ../flags.mk
-include ../user.mk

SRCDIR=.

CFLAGS+=-mpclmul
CPPFLAGS+=-mpclmul
INCFLAGS=-I$(FFTINC)
TARGET=$(BLDDIR)/libFFT.a

CUDA_LDFLAGS	=
CUDA_SOURCES	=
CPP_SOURCES	= Basis.cpp Definitions.cpp Element.cpp Polynomials.cpp FFT.cpp
ifeq ($(BUILD_BIN_TYPE),cuda)
CUDA_SOURCES	= Chunk.cu GPU_FFT.cu FFT_cuda.cu
CPP_SOURCES	= Basis.cpp Definitions.cpp Element.cpp Polynomials.cpp FFT.cpp
CUDA_LDFLAGS	+= -L$(CUDA_TOOLKIT_PATH)/../lib64 --cudart static -Xlinker -lgomp --relocatable-device-code=true -gencode arch=compute_35,code=compute_35 -gencode arch=compute_35,code=sm_35
EXTRA_NVCCFLAGS	+= -D__GPU
CPPFLAGS	+= -D__GPU
INCFLAGS	+= -I$(CUDA_TOOLKIT_PATH)/include
endif	#ifeq ($(BUILD_BIN_TYPE),cuda)

SOURCES=$(CUDA_SOURCES) $(CPP_SOURCES)

SRCS=$(addprefix src/, $(SOURCES))

-include ../shared.mk

$(TARGET): $(OBJS)
	@echo 'Building target: $@'
#	@echo 'Invoking: GCC Archiver'
#	@ar -r  "$@" $(OBJS) $(LIBS) 2> /dev/null
ifeq ($(BUILD_BIN_TYPE),cuda)
	$(NVCC) $(OBJS) $(LIBS) $(CUDA_LDFLAGS) -lib -o "$@"  
else   #ifeq ($(BUILD_BIN_TYPE),cuda)
	ar -r  "$@" $(OBJS) $(LIBS)
endif   #ifeq ($(BUILD_BIN_TYPE),cuda)
#	@echo 'Finished building target: $@'
#	@echo ' '
