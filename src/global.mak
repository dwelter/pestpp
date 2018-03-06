# Global options for GNU Make to build PEST++ libraries and programs
#
# Choose one of each:
# SYSTEM ?= linux (default)
#        ?= mac
#        ?= win
# COMPILER ?= gcc (default)
#          ?= intel
# These can be kept in local.mak
-include $(top_builddir)/local.mak

# Autodetect SYSTEM
#$(info $$OS is $(OS))
ifeq ($(OS),Windows_NT)
SYSTEM ?= win
else
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S), Linux)
SYSTEM ?= linux
endif
ifeq ($(UNAME_S), Darwin)
SYSTEM ?= mac
endif
endif

# Defaults (if unset, or missing local.mak)
#SYSTEM ?= mac
COMPILER ?= intel

ifeq ($(SYSTEM),mac)
# macOS
bindir ?= $(top_builddir)/../exe/mac/
else ifeq ($(SYSTEM),linux)
# GNU Linux
bindir ?= $(top_builddir)/../exe/linux/
else ifeq ($(SYSTEM),win)
# Microsoft Windows
bindir ?= $(top_builddir)/../exe/windows/
else
$(error SYSTEM not understood: $(SYSTEM). Use one of mac, linux or win.)
endif

ifeq ($(SYSTEM),win)
# Microsoft Windows
EXE_EXT	= .exe
OBJ_EXT	= .obj
LIB_EXT	= .lib
CP = copy
RM = del /Q
MKDIR = md
else
# POSIX (mac, linux)
OBJ_EXT	= .o
LIB_PRE	= lib
LIB_EXT	= .a
CP = cp
MKDIR = mkdir -p
endif

ifeq ($(COMPILER),intel)
# Intel compilers
ifeq ($(SYSTEM),win)
# Warning: this build method is not well tested
CXX	?= icl
OPT_FLAGS	= /nologo /O2
CXXFLAGS	= $(OPT_FLAGS) /Qstd=c++11 /EHsc
FFLAGS	= $(OPT_FLAGS) /fpp
FFREE   = /free
else # mac,linux
CXX	?= icpc
OPT_FLAGS	= -O2
CXXFLAGS	= $(OPT_FLAGS) -std=c++11
FFLAGS	= $(OPT_FLAGS) -fpp
FFREE	= -free
ifeq ($(SYSTEM),mac)
MKLROOT = /opt/intel/compilers_and_libraries_2018.1.126/mac/mkl
endif
endif
FC	?= ifort

ifeq ($(SYSTEM),win)
EXT_INCLUDES = -I"$(MKLROOT)"\include
EXT_LIBS = \
    mkl_blas95_lp64.lib \
    mkl_lapack95_lp64.lib \
    mkl_intel_lp64.lib \
    mkl_sequential.lib \
    mkl_core.lib
else ifeq ($(SYSTEM),linux)
EXT_INCLUDES = -I${MKLROOT}/include/intel64/ilp64 -I${MKLROOT}/include
EXT_LIBS = \
    ${MKLROOT}/lib/intel64/libmkl_blas95_lp64.a \
    ${MKLROOT}/lib/intel64/libmkl_lapack95_lp64.a \
    -Wl,--start-group \
        ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a \
        ${MKLROOT}/lib/intel64/libmkl_sequential.a \
        ${MKLROOT}/lib/intel64/libmkl_core.a \
    -Wl,--end-group \
    -lifport -lifcore -lpthread -lm -ldl
else ifeq ($(SYSTEM),mac)
EXTRADIR = /opt/intel/compilers_and_libraries_2018.1.126/mac/compiler/lib
EXT_INCLUDES = -I${MKLROOT}/include/intel64/lp64 -I${MKLROOT}/include
EXT_LIBS = \
    ${MKLROOT}/lib/libmkl_lapack95_ilp64.a \
    ${MKLROOT}/lib/libmkl_blas95_ilp64.a \
    ${MKLROOT}/lib/libmkl_intel_ilp64.a \
    ${MKLROOT}/lib/libmkl_sequential.a \
    ${MKLROOT}/lib/libmkl_core.a \
    ${EXTRADIR}/libifcore.a
endif
else ifeq ($(COMPILER),gcc)
# GNU Compiler Collection
CXX	?= g++
FC	?= gfortran
OPT_FLAGS	= -O2 -march=native
CXXFLAGS	= $(OPT_FLAGS) -std=c++11
FFLAGS	= $(OPT_FLAGS) -cpp
FFREE	= -ffree-form
EXT_LIBS	= -lpthread -llapack -lblas -lgfortran -lquadmath
else
$(error COMPILER not understood: $(COMPILER). Use one of intel or gcc.)
endif

# Assume linker is the C++ compiler
LD = $(CXX)
LDFLAGS += -pthread
#LDFLAGS += -static

# r=insert with replacement; c=create archive; s=add index
ARFLAGS := rcs

LIBS_DIR :=	$(top_builddir)/libs

PESTPP_INCLUDES := \
    -I $(LIBS_DIR)/Eigen \
    -I $(LIBS_DIR)/common \
    -I $(LIBS_DIR)/run_managers/abstract_base \
    -I $(LIBS_DIR)/run_managers/yamr \
    -I $(LIBS_DIR)/run_managers/serial \
    -I $(LIBS_DIR)/run_managers/genie_wrapper \
    -I $(LIBS_DIR)/run_managers/external \
    -I $(LIBS_DIR)/run_managers/wrappers \
    -I $(LIBS_DIR)/pestpp_common \
    -I $(LIBS_DIR)/opt \
    -I $(LIBS_DIR)/linear_analysis $(EXT_INCLUDES)

# Be careful with the order of library dependencies
PESTPP_LIBS := \
    -L$(LIBS_DIR)/pestpp_common -lpestpp_com \
    -L$(LIBS_DIR)/run_managers/wrappers -lrm_wrappers \
    -L$(LIBS_DIR)/run_managers/yamr -lrm_yamr \
    -L$(LIBS_DIR)/run_managers/serial -lrm_serial \
    -L$(LIBS_DIR)/run_managers/external -lrm_external \
    -L$(LIBS_DIR)/run_managers/genie_wrapper -lrm_genie_wrapper \
    -L$(LIBS_DIR)/run_managers/abstract_base -lrm_abstract \
    -L$(LIBS_DIR)/mio -lmio \
    -L$(LIBS_DIR)/common -lcommon \
    -L$(LIBS_DIR)/propack -lpropack \
    -L$(LIBS_DIR)/linear_analysis -llinear_analysis \
    -L$(LIBS_DIR)/pest_routines -lpest_routines \
    -L$(LIBS_DIR)/opt -lopt \
     $(EXT_LIBS)


# Generic pattern rules

%$(OBJ_EXT):	%.cpp
	$(CXX) -c $(PESTPP_INCLUDES) $(CXXFLAGS) $(CPPFLAGS) -o $@ $<

%$(OBJ_EXT):	%.f
	$(FC) -c $(FFLAGS) -o $@ $<

%$(OBJ_EXT):	%.for
	$(FC) -c $(FFLAGS) -o $@ $<

%$(OBJ_EXT):	%.f90
	$(FC) -c $(FFLAGS) -o $@ $<

%$(OBJ_EXT):	%.F
	$(FC) -c $(FFLAGS) $(CPPFLAGS) -o $@ $<

%$(OBJ_EXT):	%.FOR
	$(FC) -c $(FFLAGS) $(CPPFLAGS) -o $@ $<

.SUFFIXES: .c .h .cpp .hpp .f .for .F .FOR .mod $(OBJ_EXT) $(LIB_EXT)

# Default rules for handling subdirectories

%-target:
	$(MAKE) -C $*

%-install:
	$(MAKE) -C $* install

%-clean:
	$(MAKE) -C $* clean
