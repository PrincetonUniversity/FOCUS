#!/bin/sh

############################################################################################################

 ALLFILES= globals initial datalloc rdsurf rdknot rdcoils packdof bfield bnormal bmnharm fdcheck \
           solvers descent saving focus
 HFILES= $(ALLFILES:=.h)
 FFILES= $(ALLFILES:=.F90)
 PFILES= $(ALLFILES:=.pdf)
 ROBJS=$(ALLFILES:=_r.o)
 DOBJS=$(ALLFILES:=_d.o)
 NUMOBJ= ode.o

############################################################################################################

 MACROS=macros
 CC=intel # if want to use gfortran; make CC=gfortran
 FC=mpif90
 ifeq ($(CC), intel)
	RFLAGS=-r8 -mcmodel=large -O3 -m64 -unroll0 -fno-alias -ip -traceback -vec_report0
	DFLAGS=-check all -check noarg_temp_created -debug full -D DEBUG
 else
	RFLAGS=-O3 -Wall -fdefault-real-8 -ffixed-line-length-none
	GFLAGS=-g3 -Wextra -Wtarget-lifetime -fbacktrace -fbounds-check -ffpe-trap=zero -fcheck=all -DDEBUG
 endif

# FLAGS=-r8 -mcmodel=large -O3 -m64 -unroll0 -fno-alias -ip -traceback -vec_report0
# DFLAGS=-check all -check noarg_temp_created -debug full -D DEBUG
# GFLAGS=-O3 -Wall -fdefault-real-8 -ffixed-line-length-none
# DGFLAGS=-g3 -Wextra -Wtarget-lifetime -fbacktrace -fbounds-check -ffpe-trap=zero -fcheck=all -DDEBUG

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

dfocus: $(DOBJS) $(NUMOBJ)
	$(FC) -o dfocus $(DOBJS) $(NUMOBJ) $(HDF5)
	@echo "Compiling dfocus finished."


############################################################################################################
ode.o : ode.f90
	$(FC) -c $(RFLAGS) -o $@ $<

$(ROBJS): %_r.o: %.F90
	$(FC) -c $(RFLAGS) -o $@ $<  $(HDF5)

$(DOBJS): %_d.o: %.F90
	$(FC) -c $(RFLAGS) $(DFLAGS) -o $@ $<  $(HDF5)

$(FFILES): %.F90: %.h
	m4 -P $(MACROS) $< > $@

############################################################################################################

clean:
	 rm -f $(FFILES) ; rm -f $(ROBJS) ;  rm -f $(DOBJS) ; rm -f *.mod ; rm -f *.pdf ; rm -f *.toc ; rm -f *.dvi ; rm -f *.out ; rm -f .*.date ; rm -f *.tex *.aux *.blg *.log *.ps

############################################################################################################

$(PFILES): %.pdf: %.h head.tex end.tex
	@ls --full-time $*.h | cut -c 32-50 > .$*.date
	@awk -v file=$* -v date=.$*.date 'BEGIN{getline cdate < date ; FS="!latex" ; print "\\input{head} \\code{"file"}"} \
	{if(NF>1) print $$2} \
	END{print "\\hrule \\vspace{1mm} \\footnotesize $*.h last modified on "cdate";" ; print "\\input{end}"}' $*.h > $*.tex
	@echo $*.pdf
	@pdflatex -shell-escape -interaction=nonstopmode -file-line-error $*.tex | grep ".*:[0-9]*:.*" ||:
	@pdflatex -shell-escape -interaction=nonstopmode -file-line-error $*.tex | grep ".*:[0-9]*:.*" ||: 
	@pdflatex -shell-escape -interaction=nonstopmode -file-line-error $*.tex | grep ".*:[0-9]*:.*" ||: 
	

############################################################################################################

pdfs: $(PFILES)
	@echo "Please read pdfs in this directory!"

############################################################################################################
