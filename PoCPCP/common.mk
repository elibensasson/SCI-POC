WHICH				:= $(shell which which)
ECHO				:= $(shell $(WHICH) echo)
MKDIR				:= $(shell $(WHICH) mkdir)
CP					:= $(shell $(WHICH) cp)
MV					:= $(shell $(WHICH) mv)
RM					:= $(shell $(WHICH) rm)
PWD					:= $(shell $(WHICH) pwd)
DIRNAME				:= $(shell $(WHICH) dirname)
GREP				:= $(shell $(WHICH) grep)
LSPCI				:= $(shell $(WHICH) lspci)
DPKG_DEB			:= $(shell $(WHICH) dpkg-deb)

PLATFORM			:= $(shell if [ `uname -a | grep -ci cygwin` != 0 ]; then echo windows; else echo linux; fi )
CUDA_INSTALLED		:= $(shell nvcc -V 2> /dev/null | $(GREP) -c 'Cuda compilation tools')
NVIDIA_HARDWARE		:= $(shell lspci | $(GREP) -ic nvidia)
BUILD_BIN_TYPE		:= cpu
ifeq ($(PLATFORM),linux)
ifeq ($(FFTLIB),cuda)
ifeq ($(CUDA_INSTALLED), 1)
# Location of the CUDA Toolkit
CUDA_TOOLKIT_PATH	:= $(shell $(WHICH) nvcc | xargs $(DIRNAME))/..
BUILD_BIN_TYPE		:= cuda
else	#ifeq ($(CUDA_INSTALLED), 1)
$(error "nVidia Toolkit not installed in the system")
endif	#ifeq ($(CUDA_INSTALLED), 1)
else	#ifeq ($(FFTLIB),cuda)
endif	#ifeq ($(FFTLIB),cuda)
else	#ifeq ($(PLATFORM),linux)
CYGPATH				:= $(shell $(WHICH) cygpath)
#$(error "Building GPU supported FFT library is only supported on linux")
endif	#ifeq ($(PLATFORM),linux)

