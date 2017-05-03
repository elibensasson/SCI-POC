
OBJS=$(addprefix $(BLDDIR)/, $(SRCS:.cpp=.o))
DEPS=$(addprefix $(BLDDIR)/, $(SRCS:.cpp=.d))

ifeq ($(BUILD_BIN_TYPE),cuda)

ifeq ($(x86_64),1)
	OS_SIZE = 64
	OS_ARCH = x86_64
endif

NVCC	:= $(CUDA_TOOLKIT_PATH)/bin/nvcc -ccbin $(CC)
CC 		:= $(NVCC)
CFLAGS  += -D__GPU
INCFLAGS        += -I$(CUDA_TOOLKIT_PATH)/include
# internal flags
NVCCFLAGS   := -m64 ${ARCH_FLAGS}
CCFLAGS     :=
LDFLAGS     :=
# Extra user flags
EXTRA_NVCCFLAGS   ?=
EXTRA_LDFLAGS     ?=
EXTRA_CCFLAGS     ?=

CUDA_CCFLAGS :=
CUDA_CCFLAGS += $(NVCCFLAGS)
CUDA_CCFLAGS += $(EXTRA_NVCCFLAGS)
CUDA_CCFLAGS += $(addprefix -Xcompiler ,$(CCFLAGS))
CUDA_CCFLAGS += $(addprefix -Xcompiler ,$(EXTRA_CCFLAGS))

CUDA_LDFLAGS :=
CUDA_LDFLAGS += $(CUDA_CCFLAGS)
CUDA_LDFLAGS += $(addprefix -Xlinker ,$(LDFLAGS))
CUDA_LDFLAGS += $(addprefix -Xlinker ,$(EXTRA_LDFLAGS))
CUDA_LDFLAGS += --cudart static --relocatable-device-code=true -lcudart -lcudadevrt 

# Common includes and paths for CUDA
CUDA_INCLUDES  := -I../../common/inc
CUDA_LIBRARIES :=

CFLAGS:=$(addprefix -Xcompiler ,$(CFLAGS))
CPPFLAGS:=$(addprefix -Xcompiler ,$(CPPFLAGS))
INCFLAGS:=$(addprefix -Xcompiler ,$(INCFLAGS))
LIBFLAGS:=$(addprefix -Xlinker ,$(LIBFLAGS))
LNKFLAGS:=$(addprefix -Xlinker ,$(LNKFLAGS))
OTHER_CFLAGS=$(addprefix -Xcompiler , -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" )

# Gencode arguments
ifeq ($(OS_ARCH),armv7l)
#SMS ?= 30 32 35 37 50
SMS ?= 35 #37 50
else
#SMS ?= 30 35 37 50
SMS ?= 35 #37 50
endif

ifeq ($(SMS),)
$(info >>> WARNING - no SM architectures have been specified - waiving sample <<<)
SAMPLE_ENABLED := 0
endif

ifeq ($(GENCODE_FLAGS),)
# Generate SASS code for each SM architecture listed in $(SMS)
$(foreach sm,$(SMS),$(eval GENCODE_FLAGS += -gencode arch=compute_$(sm),code=sm_$(sm)))

# Generate PTX code from the highest SM architecture in $(SMS) to guarantee forward-compatibility
HIGHEST_SM := $(lastword $(sort $(SMS)))
ifneq ($(HIGHEST_SM),)
GENCODE_FLAGS += -gencode arch=compute_$(HIGHEST_SM),code=compute_$(HIGHEST_SM)
endif
endif

CUDA_LDFLAGS += $(GENCODE_FLAGS)

ifeq ($(SAMPLE_ENABLED),0)
EXEC ?= @echo "[@]"
endif

OBJS:=$(OBJS:.cu=.o)
DEPS:=$(DEPS:.cu=.d)

$(BLDDIR)/%.o: %.cu
	@echo 'Building file: $@ ($<)'
#	@echo 'Invoking: GCC C++ Compiler'
	@$(MKDIR) -p $(shell $(DIRNAME) $@)
#	@$(EXEC) $(NVCC) $(CUDA_INCLUDES) $(CUDA_CCFLAGS) $(GENCODE_FLAGS) -c -o "$@" "$<" 2> /dev/null
	$(EXEC) $(NVCC) $(CUDA_INCLUDES) $(CUDA_CCFLAGS) $(GENCODE_FLAGS) -Xcompiler -march=native -Xcompiler -fopenmp -Xptxas -v --compile --relocatable-device-code=true -c -o "$@" "$<"
#	@echo 'Finished building: $<'

endif	#ifeq ($(BUILD_BIN_TYPE),cuda)

#$(error OBJS=$(OBJS))

# Each subdirectory must supply rules for building sources it contributes
$(BLDDIR)/%.o: %.cpp
	@echo 'Building file: $@ ($<)'
#	@echo 'Invoking: GCC C++ Compiler'
	@$(MKDIR) -p $(shell $(DIRNAME) $@)
#	@$(CC) $(CFLAGS) $(CPPFLAGS) $(INCFLAGS) -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<" 2> /dev/null
	$(CC) $(CFLAGS) $(CPPFLAGS) $(INCFLAGS) $(OTHER_CFLAGS) $(PCPFLAGS) -c -o "$@" "$<"
#	@echo 'Finished building: $<'

all: $(TARGET)

clean:
	$(RM) -f $(TARGET) $(OBJS) $(DEPS)

