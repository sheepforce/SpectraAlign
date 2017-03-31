FC		= gfortran
FFLAGS		=  -O2
FFLAGS		+= -fopenmp
FFLAGS		+= -Wall

modules: 
	$(FC) -c spectra.f95 grid_functions.f95

spectraaligner: modules
	$(FC) spectraaligner.f95 spectra.o grid_functions.o -o spectraaligner.exe

clean:
	rm -f *.mod *.o *.exe
