#!/bin/sh

###############################################################################################################################################################

 hfiles=manual globals initial surface rdknot rdcoils knotxx iccoil bfield bnormal torflux tlength equarcl \
     specwid coilsep nlinear fdcheck denergy descent restart diagnos exinput identfy bnftran focus

 allfiles=$(hfiles)

###############################################################################################################################################################

 MACROS=macros

 FC=mpif90

 FLAGS=-r8 -mcmodel=large -O2 -m64 -unroll0 -fno-alias -ip -traceback # -check all -debug full -D DEBUG # man ifort for information;

 DFLAGS=

###############################################################################################################################################################

 NAG=-L$(NAG_ROOT)/lib -lnag_nag

 HDF5=-I$(HDF5_HOME)/include -L$(HDF5_HOME)/lib -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lpthread -lz -lm

 OCULUSVERSION=16
 OCULUSFLAGS=-L$(OCULUS) -loculus.$(OCULUSVERSION) -I$(OCULUS)

 WEBDIR=$(HOME)/w3_html

 date:=$(shell date)
 text:=$(shell date +%F)

###############################################################################################################################################################

xfocus: $(addsuffix .o,$(allfiles)) $(MACROS) Makefile
	$(FC) -o xfocus $(addsuffix .o,$(allfiles)) $(OCULUSFLAGS) $(NAG) $(HDF5)
	date
	/bin/echo -e "\a"

dfocus: $(addsuffix .o,$(allfiles)) $(MACROS) Makefile
	$(FC) -o dfocus $(addsuffix .o,$(allfiles)) $(OCULUSFLAGS) $(NAG) $(HDF5)
	date
	/bin/echo -e "\a"

###############################################################################################################################################################

globals.o: globals.h $(MACROS) Makefile
	m4 -P $(MACROS) globals.h > globals.F90
	$(FC) $(FLAGS) $(DFLAGS) -o globals.o -c globals.F90 $(OCULUSFLAGS)

###############################################################################################################################################################

restart.o: restart.h globals.o $(MACROS) Makefile
	m4 -P $(MACROS) restart.h > $*.F90
	$(FC) $(FLAGS) $(DFLAGS) -o restart.o -c $*.F90 $(HDF5)

bfield.o: bfield.h globals.o $(MACROS) Makefile
	m4 -P $(MACROS) bfield.h > $*.F90
	$(FC) $(FLAGS) $(DFLAGS) -o bfield.o -c $*.F90 $(OCULUSFLAGS) $(HDF5)

nlinear.o: nlinear.h globals.o $(MACROS) Makefile
	m4 -P $(MACROS) nlinear.h > $*.F90
	$(FC) $(FLAGS) $(DFLAGS) -o nlinear.o -c $*.F90 $(HDF5)

###############################################################################################################################################################

%.o: %.h globals.o $(MACROS) Makefile
	m4 -P $(MACROS) $*.h > $*.F90
	$(FC) $(FLAGS) $(DFLAGS) -o $*.o -c $*.F90 $(OCULUSFLAGS)

%.o: %.f90 Makefile
	$(FC) $(FLAGS) $(DFLAGS) -o $*.o -c $*

%.o: %.f Makefile
	$(FC) $(FLAGS) $(DFLAGS) -o $*.o -c $*.f

###############################################################################################################################################################

clean:
	rm -f *.o ; rm -f *.mod ; rm -f $(addsuffix .F90,$(allfiles)) ; rm -f *.pdf ; rm -f *.toc ; rm -f *.dvi ; rm -f *.out ; rm -f .*.date

###############################################################################################################################################################

%.pdf: %.h head.tex end.tex Makefile
ifeq ($(USER),shudson)
	#emacs -r -fn 7x14 -g 160x80+280 $*.h
	@ls --full-time $*.h | cut -c 32-50 > .$*.date
	awk -v file=$* -v date=.$*.date 'BEGIN{getline cdate < date ; FS="!latex" ; print "\\input{head} \\code{"file"}"} \
	{if(NF>1) print $$2} \
	END{print "\\hrule \\vspace{1mm} \\footnotesize $*.h last modified on "cdate";" ; print "\\input{end}"}' $*.h > $*.tex
	@echo "-------------------------------------------------------------------------------------------------------------------------------"
	@echo $*
	@echo "-------------------------------------------------------------------------------------------------------------------------------"
	latex $* ; latex $* ; latex $*
	dvips -P pdf -o $*.ps $*.dvi ; ps2pdf $*.ps
	rm -f $*.tex $*.aux $*.blg $*.log $*.ps

else
	#emacs -r -fn 7x14 -g 160x80+280 $*.h
	@ls --full-time $*.h | cut -c 32-50 > .$*.date
	awk -v file=$* -v date=.$*.date 'BEGIN{getline cdate < date ; FS="!latex" ; print "\\input{head_czhu} \\code{"file"}"} \
	{if(NF>1) print $$2} \
	END{print "\\hrule \\vspace{1mm} \\footnotesize $*.h last modified on "cdate";" ; print "\\input{end_czhu}"}' $*.h > $*.tex
	@echo "-------------------------------------------------------------------------------------------------------------------------------"
	@echo $*
	@echo "-------------------------------------------------------------------------------------------------------------------------------"
	latex $* ; latex $* ; latex $*
	dvips -P pdf -o $*.ps $*.dvi ; ps2pdf $*.ps
	rm -f $*.tex $*.aux $*.blg $*.log $*.ps
endif

###############################################################################################################################################################

pdfs: $(addsuffix .pdf,$(allfiles)) head.html
ifeq ($(USER),shudson)

	cat head.html > $(WEBDIR)/Focus/subroutines.html

	for file in $(allfiles) ; do cp $${file}.pdf $(WEBDIR)/Focus/. ; grep "!title" $${file}.h | cut -c 7- | \
	                           awk -v file=$${file} -F!\
	                            '{print "<tr><td><a href="file".pdf\">"file"</a></td><td>"$$1"</td><td>"$$2"</td></tr>"}' \
	                            >> $(WEBDIR)/Focus/subroutines.html ; \
	                          done

	echo "</table><a href="http://w3.pppl.gov/~shudson/">return to Dr. S.R. Hudson</a></body></html>" >> $(WEBDIR)/Focus/subroutines.html
else
       
	cat head_czhu.html > $(WEBDIR)/Focus/subroutines.html

	for file in $(allfiles) ; do cp $${file}.pdf $(WEBDIR)/Focus/. ; grep "!title" $${file}.h | cut -c 7- | \
	                           awk -v file=$${file} -F!\
	                            '{print "<tr><td><a href="file".pdf\">"file"</a></td><td>"$$1"</td><td>"$$2"</td></tr>"}' \
	                            >> $(WEBDIR)/Focus/subroutines.html ; \
	                          done

	echo "</table><a href="http://w3.pppl.gov/~$(USER)/">return to upper directory</a></body></html>" >> $(WEBDIR)/Focus/subroutines.html
endif

###############################################################################################################################################################
