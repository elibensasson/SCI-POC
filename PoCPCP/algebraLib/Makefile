
-include ../common.mk
-include ../flags.mk

CFLAGS+=
CPPFLAGS+=
INCFLAGS=-I$(NTLINC) -I$(ALGEBRAINC) -I$(FFTINC)
TARGET=$(BLDDIR)/libalgebraLib.a

SRCS:= $(shell ls src/*.cpp)

-include ../shared.mk

$(TARGET): $(OBJS)
	@echo 'Building target: $@'
	ar -r  "$@" $(OBJS) $(LIBS)
