# Makefile created by mkmf $Id: mkmf,v 18.0 2010/03/02 23:26:08 fms Exp $ 

SRCROOT = /home/innes/matrixrt/fortran/

CPPDEFS = -DPICKET


include /home/innes/matrixrt/fortran/mkmf.template.oxford_ubuntu_1804


.DEFAULT:
	-echo $@ does not exist.
all: main.exe
flux.o: $(SRCROOT)src/flux.F90 params.o tsrad_short_char.o ross.o short_char_ross.o rad_ross.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/flux.F90
io.o: $(SRCROOT)src/io.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/io.f90
main.o: $(SRCROOT)src/main.F90 io.o utils.o params.o timestep.o matrix.o rad_ross.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/main.F90
matrix.o: $(SRCROOT)src/matrix.f90 params.o flux.o utils.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/matrix.f90
params.o: $(SRCROOT)src/params.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/params.f90
rad_ross.o: $(SRCROOT)src/rad_ross.f90 params.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/rad_ross.f90
ross.o: $(SRCROOT)src/ross.f90 params.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/ross.f90
short_char_ross.o: $(SRCROOT)src/short_char_ross.f90 params.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/short_char_ross.f90
timestep.o: $(SRCROOT)src/timestep.f90 params.o flux.o utils.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/timestep.f90
tsrad_short_char.o: $(SRCROOT)src/tsrad_short_char.f90 params.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/tsrad_short_char.f90
utils.o: $(SRCROOT)src/utils.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/utils.f90
./rad_ross.f90: $(SRCROOT)src/rad_ross.f90
	cp $(SRCROOT)src/rad_ross.f90 .
./short_char_ross.f90: $(SRCROOT)src/short_char_ross.f90
	cp $(SRCROOT)src/short_char_ross.f90 .
./flux.F90: $(SRCROOT)src/flux.F90
	cp $(SRCROOT)src/flux.F90 .
./matrix.f90: $(SRCROOT)src/matrix.f90
	cp $(SRCROOT)src/matrix.f90 .
./tsrad_short_char.f90: $(SRCROOT)src/tsrad_short_char.f90
	cp $(SRCROOT)src/tsrad_short_char.f90 .
./utils.f90: $(SRCROOT)src/utils.f90
	cp $(SRCROOT)src/utils.f90 .
./timestep.f90: $(SRCROOT)src/timestep.f90
	cp $(SRCROOT)src/timestep.f90 .
./params.f90: $(SRCROOT)src/params.f90
	cp $(SRCROOT)src/params.f90 .
./io.f90: $(SRCROOT)src/io.f90
	cp $(SRCROOT)src/io.f90 .
./ross.f90: $(SRCROOT)src/ross.f90
	cp $(SRCROOT)src/ross.f90 .
./main.F90: $(SRCROOT)src/main.F90
	cp $(SRCROOT)src/main.F90 .
SRC = $(SRCROOT)src/ross.f90 $(SRCROOT)src/matrix.f90 $(SRCROOT)src/tsrad_short_char.f90 $(SRCROOT)src/flux.F90 $(SRCROOT)src/io.f90 $(SRCROOT)src/timestep.f90 $(SRCROOT)src/params.f90 $(SRCROOT)src/rad_ross.f90 $(SRCROOT)src/short_char_ross.f90 $(SRCROOT)src/main.F90 $(SRCROOT)src/utils.f90
OBJ = ross.o matrix.o tsrad_short_char.o flux.o io.o timestep.o params.o rad_ross.o short_char_ross.o main.o utils.o
OFF = $(SRCROOT)src/rad_ross.f90 $(SRCROOT)src/short_char_ross.f90 $(SRCROOT)src/flux.F90 $(SRCROOT)src/matrix.f90 $(SRCROOT)src/tsrad_short_char.f90 $(SRCROOT)src/utils.f90 $(SRCROOT)src/timestep.f90 $(SRCROOT)src/params.f90 $(SRCROOT)src/io.f90 $(SRCROOT)src/ross.f90 $(SRCROOT)src/main.F90
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
