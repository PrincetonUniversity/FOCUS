#!/bin/sh

############################################################################################################

 ALLFILES= globals initial datalloc rdsurf rdknot rdcoils packdof bfield bnormal bmnharm fdcheck \
           torflux length surfsep solvers descent congrad lmalg saving diagnos specinp poinplot focus
 HFILES= $(ALLFILES:=.h)
 FFILES= $(ALLFILES:=.F90)
 PFILES= $(ALLFILES:=.pdf)
 ROBJS=$(ALLFILES:=_r.o)
 DOBJS=$(ALLFILES:=_d.o)
 NUMOBJ= ode.o lmder1.o

############################################################################################################

 MACROS=macros
 CC=intel # if want to use gfortran; make CC=gfortran xfocus; otherwise using Intel
 FC=mpif90
ifeq ($(CC),gfortran)
# RFLAGS=-O3 -w -fdefault-real-8 -ffree-line-length-none -march=native -ffast-math
 RFLAGS=-O3 -w -ffree-line-length-none -march=native -ffast-math
 DFLAGS=-g3 -Wextra -Wtarget-lifetime -fbacktrace -fbounds-check -ffpe-trap=zero -fcheck=all -DDEBUG
else
# RFLAGS=-r8 -mcmodel=large -O3 -m64 -unroll0 -fno-alias -ip -traceback #-vec_report0 #-ipo -xhost
 RFLAGS=-mcmodel=large -O3 -m64 -unroll0 -fno-alias -ip -traceback #-vec_report0 #-ipo -xhost
 DFLAGS=-check all -check noarg_temp_created -debug full -D DEBUG
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

$(ROBJS): %_r.o: %.F90
	$(FC) -c $(RFLAGS) -o $@ $<  $(HDF5)

$(DOBJS): %_d.o: %.F90
	$(FC) -c $(RFLAGS) $(DFLAGS) -o $@ $<  $(HDF5)

$(FFILES): %.F90: %.h
	m4 -P $(MACROS) $< > $@

############################################################################################################

clean:
	 rm -f $(FFILES) ; rm -f $(ROBJS) ;  rm -f $(DOBJS) ; rm -f $(NUMOBJ) ; rm -f *.mod ; rm -f *.pdf

############################################################################################################

$(PFILES): %.pdf: %.h head.tex end.tex
#	@ls -lT $*.h | cut -c 35-55 > .$*.date
	@ls --full-time $*.h | cut -c 32-50 > .$*.date
	@awk -v file=$* -v date=.$*.date 'BEGIN{getline cdate < date ; FS="!latex" ; print "\\input{head} \\code{"file"}"} \
	{if(NF>1) print $$2} \
	END{print "\\vspace{1mm} \\hrule \\vspace{1mm} \\footnotesize $*.h last modified on "cdate";" ; print "\\input{end}"}' $*.h > $*.tex
	@echo $*.pdf
	@pdflatex -shell-escape -interaction=nonstopmode -file-line-error $*.tex | grep ".*:[0-9]*:.*" ||:
	@pdflatex -shell-escape -interaction=nonstopmode -file-line-error $*.tex | grep ".*:[0-9]*:.*" ||: 
	@pdflatex -shell-escape -interaction=nonstopmode -file-line-error $*.tex | grep ".*:[0-9]*:.*" ||:
	

############################################################################################################

pdfs: $(PFILES)
	@echo "Please read pdfs in this directory!"
	@rm -f $*.tex $*.aux $*.blg $*.log $*.ps .$*.date $*.toc $*.out

############################################################################################################
