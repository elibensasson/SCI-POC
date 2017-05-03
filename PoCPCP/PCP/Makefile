
-include ../common.mk
-include ../flags.mk
-include ../user.mk

CFLAGS+=
CPPFLAGS+=
INCFLAGS=-I$(PCPINC) -I$(GADGET3INC) -I$(ALGEBRAINC) -I$(ELFIOINC) -I$(FFTINC) -I$(FFTINC)/.. -I$(FFTINC)/src -isystem$(NTLINC)
LIBS=						\
	boost_serialization		\
	gadgetlib				\
	boost_program_options	\
	boost_regex				\
	gmp						\
	gtest					\
	rt						\
	boost_system			\
	boost_filesystem		\
	crypto					\
	algebraLib				\
	ntl
LIBFLAGS=$(addprefix -l, $(LIBS))
LNKFLAGS=-L"$(GADGET3LNKDIR)" -L"$(NTLLIBLNKDIR)" -L"$(ALGEBRALIBLINKDIR)" -lgomp #-rdynamic -gp
ifeq ($(BUILD_BIN_TYPE),cuda)
else   #ifeq ($(BUILD_BIN_TYPE),cuda)
LNKFLAGS	+= -L"$(FFTLIBLNKDIR)"
LIBS		+= FFT
endif   #ifeq ($(BUILD_BIN_TYPE),cuda)

TARGET=$(BLDDIR)/PCP

# M$ is insensitive, file names and dir names need to be processed
SOURCES:=$(shell cat PCP.vcxproj | grep 'ClCompile Include' | sed -e 's/.*ClCompile Include="\(.*\)".*/\1/' -e 's/\\/\//g')
SRCS:=$(foreach file,$(SOURCES),$(shell find . -iregex  ".*/$(file)"))

-include ../shared.mk

$(TARGET): $(OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC Linker'
ifeq ($(BUILD_BIN_TYPE),cuda)
	$(CC) -o "$@" $(OBJS) $(USER_OBJS) $(CUDA_LDFLAGS) "$(FFTLIBLNKDIR)/libFFT.a" $(LNKFLAGS) $(LIBFLAGS) -link 
else   #ifeq ($(BUILD_BIN_TYPE),cuda)
	$(CC) -o "$@" $(OBJS) $(USER_OBJS) -fopenmp $(LNKFLAGS) $(LIBFLAGS) 
endif   #ifeq ($(BUILD_BIN_TYPE),cuda)
	@echo 'Finished building target: $@'
