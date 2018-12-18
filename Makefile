#!/bin/sh

############################################################################################################
# comment out on 20180228 for NAG incompative
# ALLFILES= globals numrec initial surface rdknot rdcoils knotxx iccoil bfield bnormal torflux tlength equarcl \
#     specwid coilsep nlinear fdcheck denergy descent congrad truncnt restart diagnos exinput identfy \
#     bnftran hessian focus

 ALLFILES= globals numrec initial surface boozsurf rdcoils  iccoil bfield bnormal torflux tlength equarcl \
     specwid coilsep nlinear fdcheck denergy descent congrad truncnt restart diagnos exinput identfy \
     bnftran hessian focus

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
 FLAGS=-r8 -mcmodel=large -O3 -m64 -unroll0 -fno-alias -ip -traceback -D NORM #-vec_report0  # normalized to |B|
 DFLAGS=-check all -check noarg_temp_created -debug full -D DEBUG

############################################################################################################

 #NAG=-L$(NAG_ROOT)/lib -lnag
 NAG= -I$(NAG_ROOT)/nag_interface_blocks $(NAG_ROOT)/lib/libnag_mkl.a \
	-Wl,--start-group $(NAG_ROOT)/mkl_intel64_11.3.3/lib/libmkl_intel_lp64.a \
                        $(NAG_ROOT)/mkl_intel64_11.3.3/lib/libmkl_intel_thread.a \
                        $(NAG_ROOT)/mkl_intel64_11.3.3/lib/libmkl_core.a -Wl,--end-group \
        -liomp5 -lpthread -lm -ldl -lstdc++

 HDF5=-I$(HDF5_HOME)/include -L$(HDF5_HOME)/lib -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran \
-lhdf5 -lpthread -lz -lm

# OCULUSVERSION=17
# OCULUSFLAGS=-L$(OCULUS) -loculus.$(OCULUSVERSION) -I$(OCULUS)
 OCULUSFLAGS=

 WEBDIR=$(HOME)/w3_html

 date:=$(shell date)
 text:=$(shell date +%F)

############################################################################################################

xfocus: $(OBJS) tnpack.o hybrj.o modchl.o
	$(FC) -o xfocus $(OBJS) tnpack.o  hybrj.o modchl.o $(NAG) $(HDF5) $(MKL)
	@echo "Compiling xfocus finished."
	mkdir -p bin ; mv xfocus ./bin/

dfocus: $(DOBJS) tnpack.o hybrj.o modchl.o
	$(FC) -o dfocus $(DOBJS) tnpack.o hybrj.o modchl.o $(NAG) $(HDF5) $(MKL)
	@echo "Compiling dfocus finished."
	mkdir -p bin ; mv dfocus ./bin/

############################################################################################################
hybrj.o: hybrj.f
	$(FC) -c $(FLAGS) $(DFLAGS) -o $@ $<

tnpack.o: tnpack.f
	$(FC) -c $(FLAGS) $(DFLAGS) -o $@ $<

modchl.o: modchl.f
	$(FC) -c $(FLAGS) $(DFLAGS) -o $@ $<

$(OBJS): %_r.o: %.F90
	$(FC) -c $(FLAGS) -o $@ $<  $(NAG) $(HDF5) $(MKL)

$(DOBJS): %_d.o: %.F90
	$(FC) -c $(FLAGS) $(DFLAGS) -o $@ $< $(NAG) $(HDF5) $(MKL)

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
