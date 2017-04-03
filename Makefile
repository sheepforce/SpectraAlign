# the installation prefix where it will be installed
PREFIX		= ~/.local/bin

FC		= gfortran

FFLAGS		=  -O2
FFLAGS		+= -fopenmp
FFLAGS		+= -Wall


.PHONY: spectraaligner

spectraaligner: modules
	$(FC) spectraaligner.f95 spectra.o grid_functions.o -o spectraaligner.exe

modules: 
	$(FC) -c spectra.f95 grid_functions.f95

clean:
	rm -f *.mod *.o *.exe

install:
	cp spectraaligner.exe $(PREFIX)/spectraaligner
