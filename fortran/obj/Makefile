# Makefile created by mkmf $Id: mkmf,v 18.0 2010/03/02 23:26:08 fms Exp $ 

SRCROOT = /home/i/innes/matrixrt/fortran/

CPPDEFS = -DSHORT_CHAR


include /home/i/innes/matrixrt/fortran/mkmf.template.oxford_ubuntu_1804


.DEFAULT:
	-echo $@ does not exist.
all: main.exe
band_grey_rad.o: $(SRCROOT)src/band_grey_rad.f90 io.o utils.o params.o twostr_module.o nc_wrap.o phys.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/band_grey_rad.f90
condense.o: $(SRCROOT)src/condense.f90 params.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/condense.f90
convection.o: $(SRCROOT)src/convection.f90 params.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/convection.f90
flux.o: $(SRCROOT)src/flux.F90 params.o condense.o toon.o ross.o short_char_ross.o rad_ross.o band_grey_rad.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/flux.F90
io.o: $(SRCROOT)src/io.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/io.f90
main.o: $(SRCROOT)src/main.F90 io.o utils.o params.o timestep.o matrix.o rad_ross.o band_grey_rad.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/main.F90
matrix.o: $(SRCROOT)src/matrix.f90 params.o flux.o utils.o condense.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/matrix.f90
nc_wrap.o: $(SRCROOT)src/nc_wrap.f90 io.o params.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/nc_wrap.f90
params.o: $(SRCROOT)src/params.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/params.f90
phys.o: $(SRCROOT)src/phys.f90 params.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/phys.f90
rad_ross.o: $(SRCROOT)src/rad_ross.f90 params.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/rad_ross.f90
ross.o: $(SRCROOT)src/ross.f90 params.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/ross.f90
short_char_ross.o: $(SRCROOT)src/short_char_ross.f90 params.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/short_char_ross.f90
tau_mod.o: $(SRCROOT)src/tau_mod.f90 params.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/tau_mod.f90
timestep.o: $(SRCROOT)src/timestep.f90 params.o flux.o utils.o convection.o condense.o io.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/timestep.f90
toon.o: $(SRCROOT)src/toon.F90 params.o tau_mod.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/toon.F90
twostr_module.o: $(SRCROOT)src/twostr_module.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/twostr_module.f90
utils.o: $(SRCROOT)src/utils.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/utils.f90
./ross.f90: $(SRCROOT)src/ross.f90
	cp $(SRCROOT)src/ross.f90 .
./twostr_module.f90: $(SRCROOT)src/twostr_module.f90
	cp $(SRCROOT)src/twostr_module.f90 .
./params.f90: $(SRCROOT)src/params.f90
	cp $(SRCROOT)src/params.f90 .
./tau_mod.f90: $(SRCROOT)src/tau_mod.f90
	cp $(SRCROOT)src/tau_mod.f90 .
./toon.F90: $(SRCROOT)src/toon.F90
	cp $(SRCROOT)src/toon.F90 .
./timestep.f90: $(SRCROOT)src/timestep.f90
	cp $(SRCROOT)src/timestep.f90 .
./flux.F90: $(SRCROOT)src/flux.F90
	cp $(SRCROOT)src/flux.F90 .
./condense.f90: $(SRCROOT)src/condense.f90
	cp $(SRCROOT)src/condense.f90 .
./main.F90: $(SRCROOT)src/main.F90
	cp $(SRCROOT)src/main.F90 .
./utils.f90: $(SRCROOT)src/utils.f90
	cp $(SRCROOT)src/utils.f90 .
./convection.f90: $(SRCROOT)src/convection.f90
	cp $(SRCROOT)src/convection.f90 .
./io.f90: $(SRCROOT)src/io.f90
	cp $(SRCROOT)src/io.f90 .
./band_grey_rad.f90: $(SRCROOT)src/band_grey_rad.f90
	cp $(SRCROOT)src/band_grey_rad.f90 .
./matrix.f90: $(SRCROOT)src/matrix.f90
	cp $(SRCROOT)src/matrix.f90 .
./phys.f90: $(SRCROOT)src/phys.f90
	cp $(SRCROOT)src/phys.f90 .
./short_char_ross.f90: $(SRCROOT)src/short_char_ross.f90
	cp $(SRCROOT)src/short_char_ross.f90 .
./nc_wrap.f90: $(SRCROOT)src/nc_wrap.f90
	cp $(SRCROOT)src/nc_wrap.f90 .
./rad_ross.f90: $(SRCROOT)src/rad_ross.f90
	cp $(SRCROOT)src/rad_ross.f90 .
SRC = $(SRCROOT)src/toon.F90 $(SRCROOT)src/phys.f90 $(SRCROOT)src/rad_ross.f90 $(SRCROOT)src/params.f90 $(SRCROOT)src/utils.f90 $(SRCROOT)src/twostr_module.f90 $(SRCROOT)src/band_grey_rad.f90 $(SRCROOT)src/short_char_ross.f90 $(SRCROOT)src/timestep.f90 $(SRCROOT)src/tau_mod.f90 $(SRCROOT)src/io.f90 $(SRCROOT)src/matrix.f90 $(SRCROOT)src/nc_wrap.f90 $(SRCROOT)src/flux.F90 $(SRCROOT)src/convection.f90 $(SRCROOT)src/condense.f90 $(SRCROOT)src/main.F90 $(SRCROOT)src/ross.f90
OBJ = toon.o phys.o rad_ross.o params.o utils.o twostr_module.o band_grey_rad.o short_char_ross.o timestep.o tau_mod.o io.o matrix.o nc_wrap.o flux.o convection.o condense.o main.o ross.o
OFF = $(SRCROOT)src/ross.f90 $(SRCROOT)src/twostr_module.f90 $(SRCROOT)src/params.f90 $(SRCROOT)src/tau_mod.f90 $(SRCROOT)src/toon.F90 $(SRCROOT)src/timestep.f90 $(SRCROOT)src/flux.F90 $(SRCROOT)src/condense.f90 $(SRCROOT)src/main.F90 $(SRCROOT)src/utils.f90 $(SRCROOT)src/convection.f90 $(SRCROOT)src/io.f90 $(SRCROOT)src/band_grey_rad.f90 $(SRCROOT)src/matrix.f90 $(SRCROOT)src/phys.f90 $(SRCROOT)src/short_char_ross.f90 $(SRCROOT)src/nc_wrap.f90 $(SRCROOT)src/rad_ross.f90
clean: neat
	-rm -f .main.exe.cppdefs $(OBJ) main.exe
neat:
	-rm -f $(TMPFILES)
localize: $(OFF)
	cp $(OFF) .
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
main.exe: $(OBJ) 
	$(LD) $(OBJ) -o main.exe  $(LDFLAGS)
