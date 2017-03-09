F2PY=f2py
CCNAME=gfortran
FCNAME=gfortran
FCC=gfortran
SKIPFUN=parse_params read_info_headers read_brick_header read_brick_data read_particle read_mergertree_headers_1 read_mergertree_headers_2 read_mergertree_parent_of read_1 read_1_dynamics read_2_dummy read_2 read_region
SOURCES=
LIBS=$(shell pwd)/../third_party/FLAP/shared/libflap.so

CFLAGS=-O3 -fopenmp -llapack -lblas -lfftw3 -fPIC -g -fcheck=all -fbacktrace -ffpe-trap=invalid,underflow,inexact,denormal

all: types.o extrema/extrema_types.o extrema/extrema_storage.o extrema/extrema_mod.o

extrema.pyf: extrema.f90
	$(F2PY) -h $@  --overwrite-signature -m $(basename $@) $^ skip: extrema_compute_ext : \
		> /dev/null

extrema.so: types.o extrema.pyf extrema/extrema_types.o extrema/extrema_storage.o extrema/extrema_mod.o extrema.f90
	cd extrema && make objects
	$(F2PY) -m $(basename $@) -c $^ skip: extrema_compute_ext_ext : \
		-I./extrema \
		--fcompiler=$(FCNAME) --debug-capi \
		--f90flags="$(CFLAGS)" --link-lapack_opt -lgomp \
		-L./
		> /dev/null

# extrema.so:
# 	cd extrema; make objects
# 	gfortran -Iextrema $(CFLAGS) -shared -fPIC -o libextrema.so

%.o: %.f90
	$(FCC) $(CFLAGS) -Iextrema -c $^ -o $@


clean:
	cd extrema && make cleanall
	rm -f *.so *.pyc *.pyf *.o
