# required: sudo apt-get install flex gcc g++ make texinfo texi2html eclipse

# tools
-include common.mk

# 
WD						:= $(shell $(PWD))
BUILD_DIR				:= SCIPR-build
BLDDIR					:= $(WD)/$(BUILD_DIR)
ifeq ($(PLATFORM),linux)
BINTOOLS_MAKE_EXTRA		= MAKEINFO=makeinfo
ASSEMBLY_FILES_DIR		= $(WD)/AssemblyFiles
else
BINTOOLS_MAKE_EXTRA		= 
ASSEMBLY_FILES_DIR		= $(shell $(CYGPATH) -w $(WD)/AssemblyFiles)
endif
ifeq ($(PLATFORM),linux)
MSBUILD					:= 
else
MSBUILD					:= $(shell find /cygdrive/c/Windows/Microsoft.NET/Framework64 -maxdepth 1 -type d | sort | tail -n 1)/MSBuild.exe
endif
DIETLIBC_SRC_DIR		:= libc/dietlibc
COMPILER_SRC_DIR		:= $(WD)/gcc/gcc-3.4.6
BINTOOLS_SRC_DIR		:= $(WD)/binutils-2.24
GADGETLIB2_SRC_DIR		:= constraints/gadgetlib2
GADGETLIB3_SRC_DIR		:= gadgetlib/gadgetlib
FFTLIB_SRC_DIR			:= FFT/src
NTL_SRC_DIR				:= ../ntllib
TMP_DPKG_GEN_VER		:= 1.0-0
TMP_DPKG_GEN_ARC		:= amd64
TMP_DPKG_GEN_NAME		:= tinyram
TMP_DPKG_GEN_DIR_PFIX	:= /tmp/$(TMP_DPKG_GEN_NAME)

## DEV_TOOLS_INSTALL_DIR	:= /usr/local
DEV_TOOLS_INSTALL_DIR	:= $(WD)/../inst/release/$(PLATFORM)
COMPILER_BUILD_DIR		:= $(COMPILER_SRC_DIR)/$(BUILD_DIR)
BINUTILS_BUILD_DIR		:= $(BINTOOLS_SRC_DIR)/$(BUILD_DIR)

COMPILER_BUILD_DIR_16	:= $(COMPILER_BUILD_DIR)-16
BINUTILS_BUILD_DIR_16	:= $(BINUTILS_BUILD_DIR)-16
COMPILER_BUILD_DIR_32	:= $(COMPILER_BUILD_DIR)-32
BINUTILS_BUILD_DIR_32	:= $(BINUTILS_BUILD_DIR)-32

ALGEBRALIB_DIR			:= $(WD)/algebraLib
GADGETLIB2_DIR			:= $(WD)/$(GADGETLIB2_SRC_DIR)
GADGETLIB3_DIR			:= $(WD)/$(GADGETLIB3_SRC_DIR)
PCP_DIR				:= $(WD)/PCP
NTLLIB_DIR				:= $(WD)/../ntllib
FFTLIB_DIR				:= $(WD)/FFT
SCRIPT_EXT				:= $(WD)/scipr-ext.git

ALGEBRALIB_TESTS_DIR	:= $(WD)/algebraLib_tests
GADGETLIB2_TESTS_DIR	:= $(WD)/gadgetlib2_tests
GADGETLIB3_TESTS_DIR	:= $(WD)/gadgetlib_tests
PCP_TESTS_DIR			:= $(WD)/PCP
FFT_TESTS_DIR			:= $(WD)/FFT_tests

# VC target path should point to the msvc version used by the msvc IDE, 
# currently our code only supports msvc2013
VC_TARGETS_PATH			:= "C:\Program Files (x86)\MSBuild\Microsoft.Cpp\v4.0\V120"


FFTLIB					:= cpu

.PHONY: all compiler16 binutil16 dietlibc16 compiler32 binutil32 dietlibc32 install \
	compiler16-install binutil16-install dietlibc16-install compiler32-install binutil32-install dietlibc32-install \
	compiler16-clean binutil16-clean dietlibc16-clean compiler32-clean binutil32-clean dietlibc32-clean clean help \
	gadgetlib ntl pcp gadgetlib-clean ntl-clean pcp-clean
	gadgetlib-install ntl-install pcp-install
	tinyram16-deb

help:
	@$(ECHO) "make [options] [DEV_TOOLS_INSTALL_DIR=<installation directory>] [FFTLIB=<cpu|cuda>] [PCPFLAGS=\"override default flags\"]"
	@$(ECHO) "	Where options include:"
	@$(ECHO) "	all                 : You know what this is :)"
	@$(ECHO) "	clean               : Cleanup all build files"
	@$(ECHO) "	all16               : Build all tinyram16 tools"
	@$(ECHO) "	install             : Build and install all the tools"
	@$(ECHO) "	install16           : Build and install all tinyram16 the tools"
	@$(ECHO) "	compiler16          : Build tinyram16 gcc"
	@$(ECHO) "	binutil16           : Build tinyram16 binutil"
	@$(ECHO) "	dietlibc16          : Build tinyram16 dietlibc"
	@$(ECHO) "	compiler16-install  : Install tinyram16 gcc"
	@$(ECHO) "	binutil16-install   : Install tinyram16 binutil"
	@$(ECHO) "	dietlibc16-install  : Install tinyram16 dietlibc"
	@$(ECHO) "	compiler16-clean    : Cleanup tinyram16 compiler build files"
	@$(ECHO) "	binutil16-clean     : Cleanup tinyram16 binutil build files"
	@$(ECHO) "	dietlibc16-clean    : Cleanup tinyram16 dietlibc build files"
	@$(ECHO) "	all32               : Build all tinyram32 tools"
	@$(ECHO) "	install32           : Build and install all tinyram32 the tools"
	@$(ECHO) "	compiler32          : Build tinyram32 gcc"
	@$(ECHO) "	binutil32           : Build tinyram32 binutil"
	@$(ECHO) "	dietlibc32          : Build tinyram32 dietlibc"
	@$(ECHO) "	compiler32-install  : Install tinyram32 gcc"
	@$(ECHO) "	binutil32-install   : Install tinyram32 binutil"
	@$(ECHO) "	dietlibc32-install  : Install tinyram32 dietlibc"
	@$(ECHO) "	compiler32-clean    : Cleanup tinyram32 compiler build files"
	@$(ECHO) "	binutil32-clean     : Cleanup tinyram32 binutil build files"
	@$(ECHO) "	dietlibc32-clean    : Cleanup tinyram32 dietlibc build files"
	@$(ECHO) "	pcp               : Build pcp simulator"
	@$(ECHO) "	pcp-install       : Install pcp simulator"
	@$(ECHO) "	pcp-clean         : Cleanup pcp simulator"
	@$(ECHO) "	pcp-clean-all     : Cleanup pcp simulator, constaints, fftlib and ntllib"
	@$(ECHO) "	gadgetlib           : Build gadgetlib library"
	@$(ECHO) "	gadgetlib-clean     : Cleanup gadgetlib library"
	@$(ECHO) "	gadgetlib-tests     : Build gadgetlib test utility"
	@$(ECHO) "	gadgetlib-tests-clean: Cleanup gadgetlib test utility"
	@$(ECHO) "	algebralib			: build algebralib library"
	@$(ECHO) "	algebralib-tests     : Build algebralib test utility"
	@$(ECHO) "	algebralib-tests-clean: Cleanup algebralib test utility"
	@$(ECHO) "	algebralib-clean	: cleanup algebralib library"
	@$(ECHO) "	ntl                 : Build NTL library"
	@$(ECHO) "	ntl-clean           : Cleanup NTL library"
	@$(ECHO) "	fft                 : Build FFT library"
	@$(ECHO) "	fft-clean           : Cleanup FFT library"
	@$(ECHO) "	fft-tests           : Build FFT test utility"
	@$(ECHO) "	fft-tests-clean     : Cleanup FFT test utility"
	@$(ECHO) "	run-fft-tests     	: Run fft unit tests"
	@$(ECHO) "	run-gadgetlib-tests : Run gadgetlib unit tests"
	@$(ECHO) "	run-algebralib-tests : Run algebralib unit tests"
	@$(ECHO) "	run-pcp-tests     : Run pcp unit tests"
	@$(ECHO) "  tinyram16-deb     : Build debian package for tiny16 devtools (compiuler, assmebler, linker, libc)"
	@$(ECHO) "  tinyram32-deb     : Build debian package for tiny32 devtools (compiuler, assmebler, linker, libc)"
	@$(ECHO) "	help                : Print this text"
	@$(ECHO) "In order to override the default switches add PCPFLAGS=\"my gcc switches\" to the file user.mk in the same directory as this Makefile. The command line switch PCPFLAGS will supersede that of user.mk file."

all: all16 all32

all16: binutil16 compiler16 dietlibc16 pcp

all32: binutil32 compiler32 dietlibc32 pcp

install: install16 install32
	@$(ECHO)
	@$(ECHO) '===================================================================================================================================='
	@$(ECHO) 'Installation complete, configure you path to include "$(DEV_TOOLS_INSTALL_DIR)/bin"'
	@$(ECHO) 'Or run this command line: export PATH=$(DEV_TOOLS_INSTALL_DIR)/bin:$${PATH}'
	@$(ECHO) '===================================================================================================================================='

install16: $(DEV_TOOLS_INSTALL_DIR)/bin binutil16-install pcp-install compiler16-install dietlibc16-install tinyram16-headers-install tinyram16-crt0-install
	@$(ECHO)
	@$(ECHO) '===================================================================================================================================='
	@$(ECHO) 'Installation complete, configure you path to include "$(DEV_TOOLS_INSTALL_DIR)/bin"'
	@$(ECHO) 'Or run this command line: export PATH=$(DEV_TOOLS_INSTALL_DIR)/bin:$${PATH}'
	@$(ECHO) '===================================================================================================================================='

install32: $(DEV_TOOLS_INSTALL_DIR)/bin binutil32-install pcp-install compiler32-install dietlibc32-install tinyram32-headers-install tinyram32-crt0-install
	@$(ECHO)
	@$(ECHO) '===================================================================================================================================='
	@$(ECHO) 'Installation complete, configure you path to include "$(DEV_TOOLS_INSTALL_DIR)/bin"'
	@$(ECHO) 'Or run this command line: export PATH=$(DEV_TOOLS_INSTALL_DIR)/bin:$${PATH}'
	@$(ECHO) '===================================================================================================================================='

$(DEV_TOOLS_INSTALL_DIR)/bin:
	@$(MKDIR) -p $(DEV_TOOLS_INSTALL_DIR)/bin
	
	
tinyram16-headers-install:
	@$(ECHO) 'Copying include files'
	@$(CP) ./libtinyram/include/tinyram.h $(DEV_TOOLS_INSTALL_DIR)/lib/gcc/tinyram16-elf/3.4.6/include/

tinyram16-crt0-install: tinyram16-headers-install
	@$(ECHO) 'Building crt0.o'
	@PATH="$(DEV_TOOLS_INSTALL_DIR)/bin:${PATH}" $(DEV_TOOLS_INSTALL_DIR)/bin/tinyram16-elf-gcc ./libtinyram/crt0.c -c -o $(DEV_TOOLS_INSTALL_DIR)/lib/gcc/tinyram16-elf/3.4.6/crt0.o

tinyram16-deb:
	TMP_DIR=$(TMP_DPKG_GEN_DIR_PFIX)16_$(TMP_DPKG_GEN_VER)_$(TMP_DPKG_GEN_ARC);	\
	@if [ `bison --version | grep -c '2.7.12'` != 1 ];							\
	then																		\
		$(ECHO) "bison version missmatch, to fix try: 'wget http://launchpadlibrarian.net/140087283/libbison-dev_2.7.1.dfsg-1_amd64.deb http://launchpadlibrarian.net/140087282/bison_2.7.1.dfsg-1_amd64.deb' and 'sudo dpkg -i libbison-dev_2.7.1.dfsg-1_amd64.deb bison_2.7.1.dfsg-1_amd64.deb'";	\
	else																		\
		$(RM) -rf $$TMP_DIR;													\
		#$(MAKE) -C . clean compiler16-install tinyram16-headers-install tinyram16-crt0-install dietlibc16-install DEV_TOOLS_INSTALL_DIR=$$TMP_DIR;	\
		$(MAKE) -C .  compiler16-install tinyram16-headers-install tinyram16-crt0-install dietlibc16-install DEV_TOOLS_INSTALL_DIR=$$TMP_DIR;	\
		$(MKDIR) $$TMP_DIR/DEBIAN;												\
		$(ECHO) "Package: tinyram16-dev-tools" > $$TMP_DIR/DEBIAN/control;		\
		$(ECHO) "Version: 1.0-0" >> $$TMP_DIR/DEBIAN/control;					\
		$(ECHO) "Section: contrib" >> $$TMP_DIR/DEBIAN/control;					\
		$(ECHO) "Architecture: amd64" >> $$TMP_DIR/DEBIAN/control;				\
		$(ECHO) "Maintainer: lior.greenblatt@gmail.com" >> $$TMP_DIR/DEBIAN/control;							\
		$(ECHO) "Description: 16bit Tinyram development tools" >> $$TMP_DIR/DEBIAN/control;						\
		$(ECHO) " 16bit Tinyram development tools which include" >> $$TMP_DIR/DEBIAN/control;					\
		$(ECHO) " Compiler, assembler, linker, and ANSI complient libc " >> $$TMP_DIR/DEBIAN/control;			\
		$(DPKG_DEB) --build $$TMP_DIR ./$(TMP_DPKG_GEN_NAME)16_$(TMP_DPKG_GEN_VER)_$(TMP_DPKG_GEN_ARC).deb;	\
	fi

tinyram32-deb:
	TMP_DIR=$(TMP_DPKG_GEN_DIR_PFIX)32_$(TMP_DPKG_GEN_VER)_$(TMP_DPKG_GEN_ARC);	\
	@if [ `bison --version | grep -c '2.7.12'` != 1 ];							\
	then																		\
		$(ECHO) "bison version missmatch, to fix try: 'wget http://launchpadlibrarian.net/140087283/libbison-dev_2.7.1.dfsg-1_amd64.deb http://launchpadlibrarian.net/140087282/bison_2.7.1.dfsg-1_amd64.deb' and 'sudo dpkg -i libbison-dev_2.7.1.dfsg-1_amd64.deb bison_2.7.1.dfsg-1_amd64.deb'";	\
	else																		\
		$(RM) -rf $$TMP_DIR;													\
		#$(MAKE) -C . clean compiler32-install tinyram32-headers-install tinyram32-crt0-install dietlibc32-install DEV_TOOLS_INSTALL_DIR=$$TMP_DIR;	\
		$(MAKE) -C .  compiler32-install tinyram32-headers-install tinyram32-crt0-install dietlibc32-install DEV_TOOLS_INSTALL_DIR=$$TMP_DIR;	\
		$(MKDIR) $$TMP_DIR/DEBIAN;												\
		$(ECHO) "Package: tinyram32-dev-tools" > $$TMP_DIR/DEBIAN/control;		\
		$(ECHO) "Version: 1.0-0" >> $$TMP_DIR/DEBIAN/control;					\
		$(ECHO) "Section: contrib" >> $$TMP_DIR/DEBIAN/control;					\
		$(ECHO) "Architecture: amd64" >> $$TMP_DIR/DEBIAN/control;				\
		$(ECHO) "Maintainer: lior.greenblatt@gmail.com" >> $$TMP_DIR/DEBIAN/control;							\
		$(ECHO) "Description: 32bit Tinyram development tools" >> $$TMP_DIR/DEBIAN/control;						\
		$(ECHO) " 32bit Tinyram development tools which include" >> $$TMP_DIR/DEBIAN/control;					\
		$(ECHO) " Compiler, assembler, linker, and ANSI complient libc " >> $$TMP_DIR/DEBIAN/control;			\
		$(DPKG_DEB) --build $$TMP_DIR ./$(TMP_DPKG_GEN_NAME)32_$(TMP_DPKG_GEN_VER)_$(TMP_DPKG_GEN_ARC).deb;	\
	fi

$(BINUTILS_BUILD_DIR_16):
	$(MKDIR) -p "$@"
	( cd "$@"; ../configure --target=tinyram16-elf --prefix="$(DEV_TOOLS_INSTALL_DIR)"	)

binutil16: $(BINUTILS_BUILD_DIR_16)
	$(MAKE) -C "$(BINUTILS_BUILD_DIR_16)" all $(BINTOOLS_MAKE_EXTRA);

binutil16-install: $(DEV_TOOLS_INSTALL_DIR)/bin binutil16
	$(MAKE) -C "$(BINUTILS_BUILD_DIR_16)" install;

binutil16-clean:
	$(RM) -rf "$(BINUTILS_BUILD_DIR_16)" 

tinyram32-headers-install:
	@$(ECHO) 'Copying include files'
	@$(CP) ./libtinyram/include/tinyram.h $(DEV_TOOLS_INSTALL_DIR)/lib/gcc/tinyram32-elf/3.4.6/include/

tinyram32-crt0-install: tinyram32-headers-install
	@$(ECHO) 'Building crt0.o'
	@PATH="$(DEV_TOOLS_INSTALL_DIR)/bin:${PATH}" $(DEV_TOOLS_INSTALL_DIR)/bin/tinyram32-elf-gcc ./libtinyram/crt0.c -c -o $(DEV_TOOLS_INSTALL_DIR)/lib/gcc/tinyram32-elf/3.4.6/crt0.o

$(BINUTILS_BUILD_DIR_32):
	$(MKDIR) -p "$@"
	( cd "$@"; ../configure --target=tinyram32-elf --prefix="$(DEV_TOOLS_INSTALL_DIR)"	)

binutil32: $(BINUTILS_BUILD_DIR_32)
	$(MAKE) -C "$(BINUTILS_BUILD_DIR_32)" all $(BINTOOLS_MAKE_EXTRA);

binutil32-install: $(DEV_TOOLS_INSTALL_DIR)/bin binutil32
	$(MAKE) -C "$(BINUTILS_BUILD_DIR_32)" install;

binutil32-clean:
	$(RM) -rf "$(BINUTILS_BUILD_DIR_32)" 

$(COMPILER_BUILD_DIR_16):
	$(MKDIR) -p "$@";
	(																				\
		cd "$@"; PATH="$(DEV_TOOLS_INSTALL_DIR)/bin:${PATH}"						\
		../configure																\
		--prefix="$(DEV_TOOLS_INSTALL_DIR)" 										\
		--target=tinyram16-elf 														\
		--disable-multilib 															\
		--with-build-time-tools="$(DEV_TOOLS_INSTALL_DIR)/bin"						\
		--with-as="$(DEV_TOOLS_INSTALL_DIR)/bin/tinyram16-elf-as"					\
		--with-ld="$(DEV_TOOLS_INSTALL_DIR)/bin/tinyram16-elf-ld"					\
		--with-ar="$(DEV_TOOLS_INSTALL_DIR)/bin/tinyram16-elf-ar"					\
	)																				\
	
compiler16: binutil16-install $(COMPILER_BUILD_DIR_16)
	PATH="$(DEV_TOOLS_INSTALL_DIR)/bin:${PATH}" $(MAKE) -C "$(COMPILER_BUILD_DIR_16)"

compiler16-install: $(DEV_TOOLS_INSTALL_DIR)/bin compiler16
	PATH="$(DEV_TOOLS_INSTALL_DIR)/bin:${PATH}" $(MAKE) -C "$(COMPILER_BUILD_DIR_16)" install

compiler16-clean:
	$(RM) -rf "$(COMPILER_BUILD_DIR_16)"

$(COMPILER_BUILD_DIR_32):
	$(MKDIR) -p "$@";
	(																				\
		cd "$@"; PATH="$(DEV_TOOLS_INSTALL_DIR)/bin:${PATH}"						\
		../configure																\
		--prefix="$(DEV_TOOLS_INSTALL_DIR)" 										\
		--target=tinyram32-elf 														\
		--disable-multilib 															\
		--with-build-time-tools="$(DEV_TOOLS_INSTALL_DIR)/bin"						\
		--with-as="$(DEV_TOOLS_INSTALL_DIR)/bin/tinyram32-elf-as"					\
		--with-ld="$(DEV_TOOLS_INSTALL_DIR)/bin/tinyram32-elf-ld"					\
		--with-ar="$(DEV_TOOLS_INSTALL_DIR)/bin/tinyram32-elf-ar"					\
	)																				\
	
compiler32: binutil32-install $(COMPILER_BUILD_DIR_32)
	PATH="$(DEV_TOOLS_INSTALL_DIR)/bin:${PATH}" $(MAKE) -C "$(COMPILER_BUILD_DIR_32)"

compiler32-install: $(DEV_TOOLS_INSTALL_DIR)/bin compiler32
	PATH="$(DEV_TOOLS_INSTALL_DIR)/bin:${PATH}" $(MAKE) -C "$(COMPILER_BUILD_DIR_32)" install

compiler32-clean:
	$(RM) -rf "$(COMPILER_BUILD_DIR_32)"

dietlibc16: binutil16-install compiler16-install
	PATH="$(DEV_TOOLS_INSTALL_DIR)/bin:${PATH}" $(MAKE) -C "$(DIETLIBC_SRC_DIR)" MYARCH=tinyram16 DESTDIR="$(DEV_TOOLS_INSTALL_DIR)" all;

dietlibc16-install: $(DEV_TOOLS_INSTALL_DIR)/bin dietlibc16
	PATH="$(DEV_TOOLS_INSTALL_DIR)/bin:${PATH}" $(MAKE) -C "$(DIETLIBC_SRC_DIR)" MYARCH=tinyram16 DESTDIR="$(DEV_TOOLS_INSTALL_DIR)" install;

dietlibc16-clean:
	$(MAKE) -C "$(DIETLIBC_SRC_DIR)" clean

dietlibc32: binutil32-install compiler32-install
	PATH="$(DEV_TOOLS_INSTALL_DIR)/bin:${PATH}" $(MAKE) -C "$(DIETLIBC_SRC_DIR)" MYARCH=tinyram32 DESTDIR="$(DEV_TOOLS_INSTALL_DIR)" all;

dietlibc32-install: $(DEV_TOOLS_INSTALL_DIR)/bin dietlibc32
	PATH="$(DEV_TOOLS_INSTALL_DIR)/bin:${PATH}" $(MAKE) -C "$(DIETLIBC_SRC_DIR)" MYARCH=tinyram32 DESTDIR="$(DEV_TOOLS_INSTALL_DIR)" install;

dietlibc32-clean:
	$(MAKE) -C "$(DIETLIBC_SRC_DIR)" clean

fft:
ifeq ($(PLATFORM),linux)
	$(MAKE) -C $(FFTLIB_DIR) FFTLIB=$(FFTLIB) FFTINC=$(FFTLIB_DIR)/src/. BLDDIR=$(BLDDIR)/fft
else
	$(MSBUILD) Code.sln /p:Platform=x64 /p:Configuration=Release /p:VCTargetsPath=$(VC_TARGETS_PATH) /t:FFT
endif

fft-clean:
ifeq ($(PLATFORM),linux)
	$(MAKE) -C $(FFTLIB_DIR) FFTLIB=$(FFTLIB) BLDDIR=$(BLDDIR)/fft clean
else
	$(MSBUILD) Code.sln /p:Platform=x64 /p:Configuration=Release /p:VCTargetsPath=$(VC_TARGETS_PATH) /t:FFT:clean
endif

fft-tests: fft
ifeq ($(PLATFORM),linux)
	$(MAKE) -C $(FFTLIB_DIR) FFTLIB=$(FFTLIB) FFTINC=$(FFTLIB_DIR)/src/. BLDDIR=$(BLDDIR)/fft
else
	$(MSBUILD) Code.sln /p:Platform=x64 /p:Configuration=Release /p:VCTargetsPath=$(VC_TARGETS_PATH) /t:FFT_tests
endif

fft-tests-clean:
ifeq ($(PLATFORM),linux)
	$(MAKE) -C $(FFTLIB_DIR) FFTLIB=$(FFTLIB) BLDDIR=$(BLDDIR)/fft clean
else
	$(MSBUILD) Code.sln /p:Platform=x64 /p:Configuration=Release /p:VCTargetsPath=$(VC_TARGETS_PATH) /t:FFT_tests:clean
endif

gadgetlib2: ntl
ifeq ($(PLATFORM),linux)
	$(MAKE) -C $(GADGETLIB2_DIR) GADGETINC=$(GADGETLIB2_DIR)/. NTLINC=$(NTLLIB_DIR) BLDDIR=$(BLDDIR)/gadgetlib2
else
	$(MSBUILD) Code.sln /p:Platform=x64 /p:Configuration=Release /p:VCTargetsPath=$(VC_TARGETS_PATH) /t:gadgetlib2
endif

gadgetlib2-clean:
ifeq ($(PLATFORM),linux)
	$(MAKE) -C $(GADGETLIB2_DIR) BLDDIR=$(BLDDIR)/gadgetlib2 clean
else
	$(MSBUILD) Code.sln /p:Platform=x64 /p:Configuration=Release /p:VCTargetsPath=$(VC_TARGETS_PATH) /t:gadgetlib2:clean
endif

gadgetlib2-tests: gadgetlib2
ifeq ($(PLATFORM),linux)
	$(MAKE) -C $(GADGETLIB2_TESTS_DIR)		\
		NTLINC=$(NTLLIB_DIR)				\
		GADGETINC=$(GADGETLIB2_DIR)/. 		\
		GADGETLNKDIR=$(BLDDIR)/gadgetlib2	\
		NTLLIBLNKDIR=$(NTLLIB_DIR)/src		\
		FFTLIBLNKDIR=$(BLDDIR)/fft			\
		BLDDIR=$(BLDDIR)/gadgetlib2_tests
else
	$(MSBUILD) Code.sln /p:Platform=x64 /p:Configuration=Release /p:VCTargetsPath=$(VC_TARGETS_PATH) /t:gadgetlib2_tests
endif

gadgetlib2-tests-clean:
ifeq ($(PLATFORM),linux)
	$(MAKE) -C $(GADGETLIB2_TESTS_DIR)		\
		GADGETINC=$(GADGETLIB2_DIR)/. 		\
		GADGETLNKDIR=$(BLDDIR)/gadgetlib2	\
		NTLLIBLNKDIR=$(NTLLIB_DIR)/src		\
		FFTLIBLNKDIR=$(BLDDIR)/fft			\
		BLDDIR=$(BLDDIR)/gadgetlib2_tests clean
else
	$(MSBUILD) Code.sln /p:Platform=x64 /p:Configuration=Release /p:VCTargetsPath=$(VC_TARGETS_PATH) /t:gadgetlib2_tests:clean
endif

gadgetlib: ntl
ifeq ($(PLATFORM),linux)
	$(MAKE) -C $(GADGETLIB3_DIR)						\
	GADGETINC=$(GADGETLIB3_DIR)/. NTLINC=$(NTLLIB_DIR)	\
	ALGEBRAINC=$(ALGEBRALIB_DIR)/headers 				\
	ALGEBRALIBLINKDIR=$(BLDDIR)/algebralib				\
	FFTLIBLNKDIR=$(BLDDIR)/fft	FFTLIB=$(FFTLIB)		\
	NTLLIBLNKDIR=$(NTLLIB_DIR)/src						\
	FFTINC=$(FFTLIB_DIR) FFTLIBLNKDIR=$(BLDDIR)/fft		\
	BLDDIR=$(BLDDIR)/gadgetlib


#	$(MAKE) -C $(PCP_DIR) FFTLIB=$(FFTLIB) GADGET2INC=$(GADGETLIB2_DIR)/../.  GADGET3INC=$(GADGETLIB3_DIR)/../.			\
#		NTLINC=$(NTLLIB_DIR) ELFIOINC=$(SCRIPT_EXT)/elfio/elfio-2.2 	\
#		PCPINC=$(PCP_DIR)/src GADGET2LNKDIR=$(BLDDIR)/gadgetlib2 GADGET3LNKDIR=$(BLDDIR)/gadgetlib		\
#		NTLLIBLNKDIR=$(NTLLIB_DIR)/src BLDDIR=$(BLDDIR)/scipr_pcp			\
#		FFTINC=$(FFTLIB_DIR) FFTLIBLNKDIR=$(BLDDIR)/fft					\
#		ALGEBRAINC=$(ALGEBRALIB_DIR)/headers ALGEBRALIBLINKDIR=$(BLDDIR)/algebralib
else
	$(MSBUILD) Code.sln /p:Platform=x64 /p:Configuration=Release /p:VCTargetsPath=$(VC_TARGETS_PATH) /t:gadgetlib
endif

gadgetlib-clean:
ifeq ($(PLATFORM),linux)
	$(MAKE) -C $(GADGETLIB3_DIR) BLDDIR=$(BLDDIR)/gadgetlib clean
else
	$(MSBUILD) Code.sln /p:Platform=x64 /p:Configuration=Release /p:VCTargetsPath=$(VC_TARGETS_PATH) /t:gadgetlib:clean
endif

gadgetlib-tests: gadgetlib
ifeq ($(PLATFORM),linux)
	$(MAKE) -C $(GADGETLIB3_TESTS_DIR)		\
		NTLINC=$(NTLLIB_DIR)				\
		GADGETINC=$(GADGETLIB3_DIR)/. 		\
		GADGETLNKDIR=$(BLDDIR)/gadgetlib	\
		NTLLIBLNKDIR=$(NTLLIB_DIR)/src		\
		FFTLIBLNKDIR=$(BLDDIR)/fft			\
		ALGEBRAINC=$(ALGEBRALIB_DIR)/headers \
		ALGEBRALIBLINKDIR=$(BLDDIR)/algebralib	\
		BLDDIR=$(BLDDIR)/gadgetlib_tests
else
	$(MSBUILD) Code.sln /p:Platform=x64 /p:Configuration=Release /p:VCTargetsPath=$(VC_TARGETS_PATH) /t:gadgetlib_tests
endif

gadgetlib-tests-clean:
ifeq ($(PLATFORM),linux)
	$(MAKE) -C $(GADGETLIB3_TESTS_DIR)		\
		GADGETINC=$(GADGETLIB3_DIR)/. 		\
		GADGETLNKDIR=$(BLDDIR)/gadgetlib	\
		NTLLIBLNKDIR=$(NTLLIB_DIR)/src		\
		FFTLIBLNKDIR=$(BLDDIR)/fft			\
		ALGEBRAINC=$(ALGEBRALIB_DIR)/headers \
		ALGEBRALIBLINKDIR=$(BLDDIR)/algebralib	\
		BLDDIR=$(BLDDIR)/gadgetlib_tests clean
else
	$(MSBUILD) Code.sln /p:Platform=x64 /p:Configuration=Release /p:VCTargetsPath=$(VC_TARGETS_PATH) /t:gadgetlib_tests:clean
endif

ntl:
ifeq ($(PLATFORM),linux)
	$(MAKE) -C $(NTLLIB_DIR)/src BLDDIR=$(BUILD_DIR)
else
	$(MSBUILD) Code.sln /p:Platform=x64 /p:Configuration=Release /p:VCTargetsPath=$(VC_TARGETS_PATH) /t:ntllib
endif

ntl-clean:
ifeq ($(PLATFORM),linux)
	$(MAKE) -C $(NTLLIB_DIR)/src clean
else
	$(MSBUILD) Code.sln /p:Platform=x64 /p:Configuration=Release /p:VCTargetsPath=$(VC_TARGETS_PATH) /t:ntllib:clean
endif

algebralib: ntl fft
ifeq ($(PLATFORM),linux)
	$(MAKE) -C $(ALGEBRALIB_DIR) ALGEBRAINC=$(ALGEBRALIB_DIR)/headers \
	NTLINC=$(NTLLIB_DIR) BLDDIR=$(BLDDIR)/algebralib FFTINC=$(FFTLIB_DIR)/src
else
	$(MSBUILD) Code.sln /p:Platform=x64 /p:Configuration=Release /p:VCTargetsPath=$(VC_TARGETS_PATH) /t:algebraLib
endif

algebralib-clean:
ifeq ($(PLATFORM),linux)
	$(MAKE) -C $(ALGEBRALIB_DIR) BLDDIR=$(BLDDIR)/algebralib clean
else
	$(MSBUILD) Code.sln /p:Platform=x64 /p:Configuration=Release /p:VCTargetsPath=$(VC_TARGETS_PATH) /t:algebraLib:clean
endif

algebralib-tests: algebralib
ifeq ($(PLATFORM),linux)
	$(MAKE) -C $(ALGEBRALIB_TESTS_DIR)		\
		NTLINC=$(NTLLIB_DIR)				\
		ALGEBRAINC=$(ALGEBRALIB_DIR)/headers\
		FFTINC=$(FFTLIB_DIR)/src			\
		ALGEBRALNKDIR=$(BLDDIR)/algebralib	\
		NTLLIBLNKDIR=$(NTLLIB_DIR)/src		\
		FFTLIBLNKDIR=$(BLDDIR)/fft			\
		BLDDIR=$(BLDDIR)/algebralib_tests
else
	$(MSBUILD) Code.sln /p:Platform=x64 /p:Configuration=Release /p:VCTargetsPath=$(VC_TARGETS_PATH) /t:algebraLib_tests
endif

algebralib-tests-clean:
ifeq ($(PLATFORM),linux)
	$(MAKE) -C $(ALGEBRALIB_TESTS_DIR)		\
		ALGEBRAINC=$(ALGEBRALIB_DIR)/headers\
		ALGEBRALNKDIR=$(BLDDIR)/algebralib	\
		NTLLIBLNKDIR=$(NTLLIB_DIR)/src		\
		BLDDIR=$(BLDDIR)/algebralib_tests clean
else
	$(MSBUILD) Code.sln /p:Platform=x64 /p:Configuration=Release /p:VCTargetsPath=$(VC_TARGETS_PATH) /t:algebraLib_tests:clean
endif

pcp: ntl gadgetlib fft algebralib
ifeq ($(PLATFORM),linux)
	$(MAKE) -C $(PCP_DIR) FFTLIB=$(FFTLIB) GADGET2INC=$(GADGETLIB2_DIR)/../.  GADGET3INC=$(GADGETLIB3_DIR)/../.			\
		NTLINC=$(NTLLIB_DIR) ELFIOINC=$(SCRIPT_EXT)/elfio/elfio-2.2 	\
		PCPINC=$(PCP_DIR)/src GADGET3LNKDIR=$(BLDDIR)/gadgetlib		\
		NTLLIBLNKDIR=$(NTLLIB_DIR)/src BLDDIR=$(BLDDIR)/scipr_pcp			\
		FFTINC=$(FFTLIB_DIR) FFTLIBLNKDIR=$(BLDDIR)/fft					\
		ALGEBRAINC=$(ALGEBRALIB_DIR)/headers ALGEBRALIBLINKDIR=$(BLDDIR)/algebralib
else
	$(MSBUILD) Code.sln /p:Platform=x64 /p:Configuration=Release /p:VCTargetsPath=$(VC_TARGETS_PATH) /t:PCP
endif

pcp-clean-all: fft-clean gadgetlib-clean ntl-clean pcp-clean

pcp-clean:
ifeq ($(PLATFORM),linux)
	$(MAKE) -C $(PCP_DIR) BLDDIR=$(BLDDIR)/scipr_pcp clean
else
	$(MSBUILD) Code.sln /p:Platform=x64 /p:Configuration=Release /p:VCTargetsPath=$(VC_TARGETS_PATH) /t:PCP:clean
endif

pcp-install: $(DEV_TOOLS_INSTALL_DIR)/bin pcp
ifeq ($(PLATFORM),linux)
	$(CP) $(BUILD_DIR)/scipr_pcp/PCP $(DEV_TOOLS_INSTALL_DIR)/bin/.
else
	$(CP) $(PCP_DIR)/Release/x64/PCP.exe $(DEV_TOOLS_INSTALL_DIR)/bin/.
endif

run-fft-tests: fft-tests
ifeq ($(PLATFORM),linux)
	echo "======================= not supported yet =========================="
	exit 1
else
	echo "======================= not supported yet =========================="
	exit 1
#	$(FFT_TESTS_DIR)/Release/x64/FFT_tests.exe
endif

run-algebralib-tests: algebralib-tests
ifeq ($(PLATFORM),linux)
	$(BLDDIR)/algebralib_tests/algebralib_tests
else
	$(ALGEBRALIB_TESTS_DIR)/Release/x64/algebraLib_tests.exe
endif

run-gadgetlib2-tests: gadgetlib2-tests
ifeq ($(PLATFORM),linux)
	$(BLDDIR)/gadgetlib2_tests/gadgetlib2_tests
else
	$(GADGETLIB2_TESTS_DIR)/Release/x64/gadgetlib2_tests.exe
endif

run-gadgetlib-tests: gadgetlib-tests
ifeq ($(PLATFORM),linux)
	$(BLDDIR)/gadgetlib_tests/gadgetlib_tests
else
	$(GADGETLIB3_TESTS_DIR)/Release/x64/gadgetlib_tests.exe
endif

run-pcp-tests: pcp
ifeq ($(PLATFORM),linux)
	echo "======================= not supported yet =========================="
	$(BLDDIR)/scipr_pcp/PCP --gtest "" --files-path "$(ASSEMBLY_FILES_DIR)"
else
#	$(PCP_TESTS_DIR)/Release/x64/PCP.exe --gtest "--gtest_filter=-*verify_full_stack_degree_LEAKS_SLOW" --files-path "$(ASSEMBLY_FILES_DIR)"
	$(PCP_TESTS_DIR)/Release/x64/PCP.exe --gtest "" --files-path "$(ASSEMBLY_FILES_DIR)"
endif

clean: gadgetlib-clean ntl-clean pcp-clean fft-clean

clean16: binutil16-clean compiler16-clean dietlibc16-clean gadgetlib-clean ntl-clean pcp-clean

clean32: binutil32-clean compiler32-clean dietlibc32-clean gadgetlib-clean ntl-clean pcp-clean

