SOURCES = ../sources

# include the original Makefile
include ${SOURCES}/Makefile

#=======================================================================
#                   define the compiler names
#=======================================================================

CC       = gcc
F90      = gfortran
#F90      = ifort
#F90      =  /opt/intel/composer_xe_2015.3.187/bin/intel64/ifort
PYTHON   = python

#=======================================================================
#                     additional flags
#=======================================================================

ifeq ($(F90),gfortran)
	FPP      = $(FC) -E
	FPP_F90FLAGS = -x f95-cpp-input -fPIC
	F90FLAGS = $(RFLAGS)
    FCOMP    = $(FC) #gfortran
    LIBS     =
endif

ifeq ($(F90),ifort)

	FPP      = gfortran -E # gfortran f90wrap temp files only. not compilation
	FPP_F90FLAGS = -x f95-cpp-input -fPIC
	F90FLAGS = -fpscomp logicals -fPIC # use 1 and 0 for True and False
    FCOMP    = intelem # for f2py
    LIBS =
endif

CFLAGS = -fPIC #     ==> universal for ifort, gfortran, pgi

#=======================================================================
#=======================================================================

UNAME = $(shell uname)

ifeq (${UNAME}, Darwin)
  LIBTOOL = libtool -static -o
else
  LIBTOOL = ar src
endif

# ======================================================================
# PROJECT CONFIG, do not put spaced behind the variables
# ======================================================================
# Python module name
PYTHON_MODN = focus
# mapping between Fortran and C types
KIND_MAP = kind_map
# static library name
LIB_NAME = lib$(PYTHON_MODN).a

#=======================================================================
#       List all source files that require a Python interface
#=======================================================================

# names (without suffix), f90 sources
TMP_FILES := $(filter-out focus,$(ALLFILES))
LIBSRC_WRAP_SOURCES := $(TMP_FILES) focuspy

# file names
LIBSRC_WRAP_FILES = $(addsuffix .f90,${LIBSRC_WRAP_SOURCES})

# macro_explained files
LIBSRC_WRAP_MFILES = $(addsuffix _m.F90,${LIBSRC_WRAP_SOURCES})

# object files
LIBSRC_WRAP_OBJECTS = $(addsuffix _r.o,${LIBSRC_WRAP_SOURCES})

# fpp files
LIBSRC_WRAP_FPP_FILES = $(addsuffix .fpp,${LIBSRC_WRAP_SOURCES})

#=======================================================================
#                 Relevant suffixes
#=======================================================================

.SUFFIXES: .f90 .fpp

#=======================================================================
#
#=======================================================================

.PHONY: all all_clean


all: _${PYTHON_MODN}.so  _${PYTHON_MODN}_pkg.so test

all_clean: f90wrap_clean
	make -C ${SOURCES} clean ALLFILES="${ALLFILES} focuspy"
	rm -f ${SOURCES}/focuspy.f90

f90wrap_clean:
	-rm -f ${LIBSRC_OBJECTS} ${LIBSRC_FPP_FILES} $(LIB_NAME) _${PYTHON_MODN}*.so \
	_${PYTHON_MODN}_pkg.so *.mod *.fpp f90wrap*.f90 f90wrap*.o *.o ${PYTHON_MODN}.py
	-rm -rf ${PYTHON_MODN}_pkg
	-rm -rf src.*/ .f2py_f2cmap .libs/ __pycache__/
	-rm -rf $(LIB_NAME)

$(LIBSRC_WRAP_FPP_FILES): %.fpp: ${SOURCES}/%_m.F90
	${FPP} ${FPP_F90FLAGS} $<  -o $@

$(LIB_NAME): $(addprefix, ${SOURCES}/, $(HFILES)) focuspy.f90
	cp focuspy.f90 ${SOURCES}/
	make -C ${SOURCES} xfocus CC=gfortran PFLAGS=${CFLAGS} ALLFILES="${ALLFILES} focuspy"
	cd ${SOURCES}; ${LIBTOOL} ../python/$@ $(ROBJS) $(NUMOBJ) focuspy_r.o
	ranlib $@


_${PYTHON_MODN}.so: $(LIB_NAME) ${LIBSRC_WRAP_FPP_FILES} 
	f90wrap -m ${PYTHON_MODN} ${LIBSRC_WRAP_FPP_FILES} -k ${KIND_MAP} -v
	f2py-f90wrap --f90exec=${FC} --build-dir . -c -m _${PYTHON_MODN} -L. -l$(PYTHON_MODN) -I${SOURCES} $(HDF5) f90wrap*.f90
	sed -i '' 's/""", globals/""", Globals/g' focus.py


_${PYTHON_MODN}_pkg.so: $(LIB_NAME) ${LIBSRC_WRAP_FPP_FILES}
	f90wrap -m ${PYTHON_MODN}_pkg ${LIBSRC_WRAP_FPP_FILES} -k ${KIND_MAP} -v -P
	f2py-f90wrap --f90exec=${FC} --build-dir . -c -m _${PYTHON_MODN}_pkg -L. -l$(PYTHON_MODN) -I${SOURCES} $(HDF5) f90wrap*.f90
	sed -i '' 's/""", globals/""", Globals/g' focus.py

test:
	@python -c "import focus; print('FOCUS version: {:}'.format(focus.globals.version))"
