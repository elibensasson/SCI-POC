
-include ../../common.mk
-include ../../flags.mk
-include ../../user.mk

SRCDIR=.

CFLAGS+=
CPPFLAGS+=
INCFLAGS=-I$(GADGETINC) -I$(GADGETINC)/..  -I$(ALGEBRAINC) -I$(FFTINC) -I$(FFTINC)/.. -I$(FFTINC)/src -isystem$(NTLINC)
TARGET=$(BLDDIR)/libgadgetlib.a

SRCS=					\
	common_use.cpp  	\
	constraint.cpp  	\
	gadget.cpp  		\
	infrastructure.cpp  \
	protoboard.cpp

-include ../../shared.mk

$(TARGET): $(OBJS)
	@echo 'Building target: $@'
#	@echo 'Invoking: GCC Archiver'
	ar -r  "$@" $(OBJS) $(LIBS)
#	@echo 'Finished building target: $@'
#	@echo ' '
