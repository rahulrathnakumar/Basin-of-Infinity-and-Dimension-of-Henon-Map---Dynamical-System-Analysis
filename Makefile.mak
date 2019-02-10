SOURCES = precision.f90 lsq.f90 main.f90 henondim.f90
OBJECTS = precision.o lsq.o henondim.o main.o
MODULES = precision.mod lsq.mod henondim.mod





FC = ifort
OPT = 
LIBS = -mkl
PARALLEL = -qopenmp
DEBUG = -g
OTHER = -c
FFLAGS = $(OPT) $(OTHER) $(PARALLEL)

%.o : %.f90
		  $(FC) $(FFLAGS) $<
%.o :.%f9
		  $(FC) $(FFLAGS) $<
all:  	 test_linfit
test_linfit: $(OBJECTS)
		$(FC) $(PARALLEL) $(OBJECTS) $(LIBS) -o main
precision.o:
lsq.o: precision.o
henondim.o: precision.o
main.o: precision.o
main.o: henondim.o
main.o: lsq.o

clean:
	    rm -rf $(OBJECTS) $(MODULES) a.out core