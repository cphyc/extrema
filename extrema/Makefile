.SUFFIXES:
.SUFFIXES:  .out .o .f90 .for .f
# Edit this to be correct on your system
F90 = gfortran

F90OPTS = -O3 -fopenmp -fPIC -llapack -lblas -lfftw3

#-----------------------------------------------------

OBJ = extrema_types.o extrema_storage.o extrema_mod.o

objects: $(OBJ)
#----------------------------------------------------

%.o: %.f90
	$(F90) $(F90OPTS) $(LFLAGS) -c $^

#----------------------------------------------------------

clean:
	/bin/rm -f *.o

cleanall:
	/bin/rm -f *.o $(MODULEDIR)/[!n]*.mod
