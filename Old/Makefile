#!/bin/sh

############################################################################################################

 ALLFILES= globals numrec initial surface rdknot rdcoils knotxx iccoil bfield bnormal torflux tlength equarcl \
     specwid coilsep nlinear fdcheck denergy descent restart diagnos exinput identfy bnftran congrad focus

############################################################################################################

 HFILES= $(ALLFILES:=.h)
 FFILES= $(ALLFILES:=.F90)
 PFILES= $(ALLFILES:=.pdf)
 OBJS=$(ALLFILES:=_r.o)
 DOBJS=$(ALLFILES:=_d.o)

############################################################################################################

 MACROS=macros
 CC=intel # if want to use gfortran; make CC=gfortran
 FC=mpif90
 FLAGS=-r8 -mcmodel=large -O3 -m64 -unroll0 -fno-alias -ip -traceback -vec_report0
 DFLAGS=-check all -check noarg_temp_created -debug full -D DEBUG

############################################################################################################

 NAG=-L$(NAG_ROOT)/lib -lnag_nag

 HDF5=-I$(HDF5_HOME)/include -L$(HDF5_HOME)/lib -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran \
-lhdf5 -lpthread -lz -lm

 OCULUSVERSION=17
 OCULUSFLAGS=-L$(OCULUS) -loculus.$(OCULUSVERSION) -I$(OCULUS)

 WEBDIR=$(HOME)/w3_html

 date:=$(shell date)
 text:=$(shell date +%F)

############################################################################################################

xfocus: $(OBJS)
	$(FC) -o xfocus $(OBJS) $(OCULUSFLAGS) $(NAG) $(HDF5)
	@echo "Compiling xfocus finished."

dfocus: $(DOBJS)
	$(FC) -o dfocus $(DOBJS) $(OCULUSFLAGS) $(NAG) $(HDF5)
	@echo "Compiling dfocus finished."

############################################################################################################
$(OBJS): %_r.o: %.F90
	$(FC) -c $(FLAGS) -o $@ $<  $(OCULUSFLAGS) $(NAG) $(HDF5)

$(DOBJS): %_d.o: %.F90
	$(FC) -c $(FLAGS) $(DFLAGS) -o $@ $< $(OCULUSFLAGS) $(NAG) $(HDF5)

$(FFILES): %.F90: %.h
	m4 -P $(MACROS) $< > $@

 date:=$(shell date)
 text:=$(shell date +%F)

clean:
	rm -f *.o ; rm -f *.mod ; rm -f $(FFILES) ; rm -f *.pdf ; rm -f *.toc ; rm -f *.dvi ; rm -f *.out ; rm -f .*.date

############################################################################################################

$(PFILES): %.pdf: %.h head.tex end.tex
ifeq ($(USER),czhu)
	@ls --full-time $*.h | cut -c 32-50 > .$*.date
	@awk -v file=$* -v date=.$*.date 'BEGIN{getline cdate < date ; FS="!latex" ; print "\\input{head} \\code{"file"}"} \
	{if(NF>1) print $$2} \
	END{print "\\hrule \\vspace{1mm} \\footnotesize $*.h last modified on "cdate";" ; print "\\input{end}"}' $*.h > $*.tex
	@echo $*.pdf
	@pdflatex -shell-escape -interaction=nonstopmode -file-line-error $*.tex | grep ".*:[0-9]*:.*" ||:
	@pdflatex -shell-escape -interaction=nonstopmode -file-line-error $*.tex | grep ".*:[0-9]*:.*" ||: 
	@pdflatex -shell-escape -interaction=nonstopmode -file-line-error $*.tex | grep ".*:[0-9]*:.*" ||: 
	@rm -f $*.tex $*.aux $*.blg $*.log $*.ps .$*.date $*.toc $*.out
	@mv $*.pdf ../docs
endif

############################################################################################################

pdfs: $(PFILES)
ifneq ($(USER),czhu)
	@echo "Please read pdfs in this directory!"
else
	@echo "All pdfs are moved to ../docs"
endif
