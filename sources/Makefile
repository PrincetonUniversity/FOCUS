#!/bin/sh

############################################################################################################

 ALLFILES= globals ncsx_ports_mod ncsx_ports_load ncsx_ports_eval wtmgrid poinplot initial  \
           rdsurf rdknot rdboozer rdcoils packdof bfield bmnharm bnormal fdcheck \
           torflux length surfsep minvol bindip totalpm datalloc \
		   solvers descent congrad qnewton lmalg saving diagnos \
           specinp boozer focus            
 HFILES= $(ALLFILES:=.f90)   # raw source files
 FFILES= $(ALLFILES:=_m.F90) # Fortran 90 files 
 PFILES= $(ALLFILES:=.pdf)   # documentations
 ROBJS=$(ALLFILES:=_r.o)     # release version objectives
 DOBJS=$(ALLFILES:=_d.o)     # debug version objectives
 NUMOBJ= ode.o lmder1.o hybrj.o lbfgsb.o cg_descent.o simann.o # numerical libraries

############################################################################################################

 MACROS=macros

#### Default Intel+OpenMPI+HDF5###############
# available env: /p/focus/modules/focus/basis
 CC=intel # if want to use gfortran; make CC=gfortran or other options
 FC=mpif90
 PFLAGS= -D TOPO # for pre-processing compiler flags, like -D dposition
 RFLAGS=-mcmodel=large -O3 -m64 -unroll0 -fno-alias -ip -traceback  $(PFLAGS) #-vec_report0 #-ipo -xhost
 DFLAGS=-O0 -g -traceback $(PFLAGS) -check all -check bounds -check noarg_temp_created -check uninit \
     -ftrapuv -init=snan,arrays -debug all -D DEBUG
 LIBS=-I$(HDF5_HOME)/include
 LINKS=-L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -lpthread -lz -lm

 LMGRID = F
 ifeq ($(LMGRID),T)
   PFLAGS += -D mgrid_poincare
   LIBS += -I$(STELLOPT_PATH)/LIBSTELL/Release
   LINKS += -L$(STELLOPT_PATH)/LIBSTELL/Release -lstell
 endif

#### GFORTRAN+OpenMPI+HDF5###############
# available env: /p/focus/modules/focus/gfortran
ifeq ($(CC),gfortran)
# RFLAGS=-O3 -w -fdefault-real-8 -ffree-line-length-none -march=native -ffast-math
 # PFLAGS+=-fallow-argument-mismatch # only used for GCC-10
 RFLAGS=-O3 -w -ffree-line-length-none -ffast-math $(PFLAGS)
 DFLAGS=-Og $(PFLAGS) -w -ffree-line-length-none -Wextra -Wtarget-lifetime -fbacktrace -fbounds-check  \
     -ffpe-trap=zero -fcheck=all -DDEBUG
endif

#### Stellar Intel+IntelMPI+HDF5###############
# available env: /home/caoxiang/module/focus/intel
ifneq (,$(findstring stellar,$(HOSTNAME)))
 FC=mpiifort  # this is optional after recent updates
endif

#### Traverse + GNU Fortran ################
ifneq (,$(findstring traverse,$(HOSTNAME)))
# available env: /home/caoxiang/Modules/focus/
 RFLAGS=-O3 -w -ffree-line-length-none -mcpu=power9 -ffast-math $(PFLAGS)
 DFLAGS=-Og $(PFLAGS) -w -ffree-line-length-none -Wextra -Wtarget-lifetime -fbacktrace -fbounds-check  \
     -ffpe-trap=zero -fcheck=all -DDEBUG
endif

#### IPP draco Intel+IntelMPI+HDF5###############
ifneq (,$(findstring draco,$(HOSTNAME)))
 FC=mpiifort  # this is optional after recent updates
 RFLAGS+=-xCORE-AVX2
 LINKS=-L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -lpthread -lz -lm -Wl,-rpath -Wl,$(HDF5_HOME)/lib
endif

#### IPP cobra Intel+IntelMPI+HDF5###############
ifneq (,$(findstring cobra,$(HOSTNAME)))
 FC=mpiifort  # this is optional after recent updates
 RFLAGS+=-xCORE-AVX512 -qopt-zmm-usage=high
 LINKS=-L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -lpthread -lz -lm -Wl,-rpath -Wl,$(HDF5_HOME)/lib
endif


#### IPP raven Intel+IntelMPI+HDF5###############
ifneq (,$(findstring raven,$(HOSTNAME)))
 FC=mpiifort  # this is optional after recent updates
 RFLAGS+=-xCORE-AVX512 -qopt-zmm-usage=high
 LINKS=-L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -lpthread -lz -lm -Wl,-rpath -Wl,$(HDF5_HOME)/lib
endif

############################################################################################################

 WEBDIR=$(HOME)/w3_html

 date:=$(shell date)
 text:=$(shell date +%F)

############################################################################################################

xfocus: $(ROBJS) $(NUMOBJ)
	$(FC) -o xfocus $(ROBJS) $(NUMOBJ) $(LINKS)
	@echo "Compiling xfocus finished."
	mkdir -p ../bin ; mv xfocus ../bin/

dfocus: $(DOBJS) $(NUMOBJ)
	$(FC) -o dfocus $(DOBJS) $(NUMOBJ) $(LINKS)
	@echo "Compiling dfocus finished."
	mkdir -p ../bin ; mv dfocus ../bin/

############################################################################################################
ode.o : ode.f90 
	$(FC) -c $(RFLAGS) -o $@ $< $(LIBS)
lmder1.o : lmder1.f
	$(FC) -c $(RFLAGS) -o $@ $< $(LIBS)

hybrj.o: hybrj.f
	$(FC) -c $(RFLAGS) -o $@ $< $(LIBS)

lbfgsb.o: lbfgsb.f
	$(FC) -c $(RFLAGS) -o $@ $< $(LIBS)

cg_descent.o : cg_descent.f
	$(FC) -c $(RFLAGS) -o $@ $< $(LIBS)

simann.o : simann.f90
	$(FC) -c $(RFLAGS) -o $@ $< $(LIBS)

$(ROBJS): %_r.o: %_m.F90 globals_r.o
	$(FC) -c $(RFLAGS) -o $@ $<  $(LIBS)

$(DOBJS): %_d.o: %_m.F90 globals_d.o
	$(FC) -c $(DFLAGS) -o $@ $<  $(LIBS)

$(FFILES): %_m.F90: %.f90
	m4 -P $(MACROS) $< > $@


############################################################################################################

clean:
	 rm -f $(FFILES) ; rm -f $(ROBJS) ;  rm -f $(DOBJS) ; rm -f $(NUMOBJ) ; rm -f *.mod ; rm -f *.pdf

############################################################################################################

$(PFILES): %.pdf: %.f90 head.tex end.tex
#	@ls -lT $*.f90 | cut -c 35-55 > .$*.date
	@ls --full-time $*.f90 | cut -c 32-50 > .$*.date
	@awk -v file=$* -v date=.$*.date 'BEGIN{getline cdate < date ; FS="!latex" ; print "\\input{head} \\code{"file"}"} \
	{if(NF>1) print $$2} \
	END{print "\\vspace{1mm} \\hrule \\vspace{1mm} \\footnotesize $*.f90 last modified on "cdate";" ; print "\\input{end}"}' $*.f90 > $*.tex
	@echo $*.pdf
	@pdflatex -shell-escape -interaction=nonstopmode -file-line-error $*.tex | grep ".*:[0-9]*:.*" ||:
	@pdflatex -shell-escape -interaction=nonstopmode -file-line-error $*.tex | grep ".*:[0-9]*:.*" ||: 
	@pdflatex -shell-escape -interaction=nonstopmode -file-line-error $*.tex | grep ".*:[0-9]*:.*" ||:
	@rm -f $*.tex $*.aux $*.blg $*.log $*.ps .$*.date $*.toc $*.out
	
############################################################################################################

pdfs: $(PFILES)
	@echo "Please read pdfs in this directory!"

############################################################################################################
