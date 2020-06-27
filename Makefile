FFLAGS =  -O  #-fopenmp 
FF = -c
FC = gfortran

LIBBLAS = -lblas
LIBLAPACK = -llapack
LIBFFTW3 = -lfftw3

LIBS = $(LIBFFTW3) $(LIBLAPACK) $(LIBBLAS)


OBJ = declarations.o  used_files.o  present_version.o \
      algebra.o declfft.o io.o   \
      electronic.o convolution.o selfenergies.o \
      initialize.o green.o WNncaVIB.o

WNncaVIB.out: $(OBJ)
	$(FC) $(FFLAGS) $(OBJ) $(LIBS) -o WNncaVIB.out

declarations.o: declarations.f90
	$(FC) $(FFLAGS) $(FF)  declarations.f90

declfft.o: declfft.f90
	$(FC) $(FFLAGS) $(FF) declfft.f90 

initialize.o: initialize.f90
	$(FC) $(FFLAGS) $(FF) initialize.f90

io.o: io.f90
	$(FC) $(FFLAGS) $(FF) io.f90

convolution.o: convolution.f90
	$(FC) $(FFLAGS) $(FF) convolution.f90

selfenergies.o: selfenergies.f90
	$(FC) $(FFLAGS) $(FF) selfenergies.f90

algebra.o: algebra.f90
	$(FC) $(FFLAGS) $(FF)  algebra.f90

electronic.o: electronic.f90
	$(FC) $(FFLAGS) $(FF) electronic.f90

green.o: green.f90
	$(FC) $(FFLAGS) $(FF) green.f90

present_version.o: present_version.f90
	$(FC) $(FFLAGS) $(FF) present_version.f90

used_files.o: used_files.f90
	$(FC) $(FFLAGS) $(FF) used_files.f90

WNncaVIB.o: WNncaVIB.f90
	$(FC) $(FFLAGS) $(FF) WNncaVIB.f90

clean:
	@echo "Cleaning the directory"
	rm -f *.o *.mod

