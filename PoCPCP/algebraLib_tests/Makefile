
-include ../common.mk
-include ../flags.mk

SRCDIR=.

CFLAGS+=
CPPFLAGS+=
INCFLAGS=-I$(ALGEBRAINC) -I$(FFTINC) -isystem$(NTLINC)
LIBS=						\
	algebraLib				\
	gtest					\
	ntl						

LIBFLAGS=$(addprefix -l, $(LIBS))
LNKFLAGS=-L"$(ALGEBRALNKDIR)" -L"$(NTLLIBLNKDIR)" -lgomp
ifeq ($(BUILD_BIN_TYPE),cuda)
else   #ifeq ($(BUILD_BIN_TYPE),cuda)
LNKFLAGS	+= -L"$(FFTLIBLNKDIR)"
LIBS		+= FFT
endif   #ifeq ($(BUILD_BIN_TYPE),cuda)

TARGET=$(BLDDIR)/algebralib_tests

SRCS:=	$(shell ls *.cpp)

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
