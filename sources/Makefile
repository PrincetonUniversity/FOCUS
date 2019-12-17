#!/bin/sh

############################################################################################################

 ALLFILES= globals ncsx_ports_mod ncsx_ports_load ncsx_ports_eval initial rdsurf rdknot rdcoils packdof \
           bfield bmnharm bnormal fdcheck \
           torflux length surfsep minvol datalloc solvers descent congrad qnewton lmalg saving diagnos \
           specinp poinplot boozer wtmgrid focus            
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

#### GFORTRAN+OpenMPI+HDF5###############
# available env: /p/focus/modules/focus/gfortran
ifeq ($(CC),gfortran)
# RFLAGS=-O3 -w -fdefault-real-8 -ffree-line-length-none -march=native -ffast-math
 RFLAGS=-O3 -w -ffree-line-length-none -march=native -ffast-math $(PFLAGS)
 DFLAGS=-Og $(PFLAGS) -w -ffree-line-length-none -Wextra -Wtarget-lifetime -fbacktrace -fbounds-check  \
     -ffpe-trap=zero -fcheck=all -DDEBUG
endif

#### EDDY Intel+IntelMPI+HDF5###############
# available env: /home/caoxiang/Modules/focus/develop
ifneq (,$(findstring eddy,$(HOSTNAME)))
 FC=mpiifort  # this is optional after recent updates
endif

#### Traverse + GNU Fortran ################
ifneq (,$(findstring traverse,$(HOSTNAME)))
# available env: /home/caoxiang/Modules/focus/
 RFLAGS=-O3 -w -ffree-line-length-none -mcpu=power9 -ffast-math $(PFLAGS)
 DFLAGS=-Og $(PFLAGS) -w -ffree-line-length-none -Wextra -Wtarget-lifetime -fbacktrace -fbounds-check  \
     -ffpe-trap=zero -fcheck=all -DDEBUG
endif

############################################################################################################

 HDF5=-I$(HDF5_HOME)/include -L$(HDF5_HOME)/lib -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran \
-lhdf5 -lpthread -lz -lm

 WEBDIR=$(HOME)/w3_html

 date:=$(shell date)
 text:=$(shell date +%F)

############################################################################################################

xfocus: $(ROBJS) $(NUMOBJ)
	$(FC) -o xfocus $(ROBJS) $(NUMOBJ) $(HDF5)
	@echo "Compiling xfocus finished."
	mkdir -p ../bin ; mv xfocus ../bin/

dfocus: $(DOBJS) $(NUMOBJ)
	$(FC) -o dfocus $(DOBJS) $(NUMOBJ) $(HDF5)
	@echo "Compiling dfocus finished."
	mkdir -p ../bin ; mv dfocus ../bin/

############################################################################################################
ode.o : ode.f90
	$(FC) -c $(RFLAGS) -o $@ $<

lmder1.o : lmder1.f
	$(FC) -c $(RFLAGS) -o $@ $<

hybrj.o: hybrj.f
	$(FC) -c $(RFLAGS) -o $@ $<

lbfgsb.o: lbfgsb.f
	$(FC) -c $(RFLAGS) -o $@ $<

cg_descent.o : cg_descent.f
	$(FC) -c $(RFLAGS) -o $@ $<

simann.o : simann.f90
	$(FC) -c $(RFLAGS) -o $@ $<

$(ROBJS): %_r.o: %_m.F90
	$(FC) -c $(RFLAGS) -o $@ $<  $(HDF5)

$(DOBJS): %_d.o: %_m.F90
	$(FC) -c $(DFLAGS) -o $@ $<  $(HDF5)

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
