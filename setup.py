from numpy.distutils.core import setup, Extension

debug = False
F90CFLAGS = ["-fopenmp"]
if debug:
    F90CFLAGS += ["-g", "-fcheck=all", "-fbacktrace",
                  "-ffpe-trap=invalid,underflow,inexact,denormal"]

F90LFLAGS = ["-llapack", "-lblas", "-lfftw3", "-fopenmp"]
F90SRC = ['extrema/extrema_mod.f90',
          'extrema/extrema_types.f90',
          'extrema/extrema_storage.f90']

F90OBJ = [e.replace(".f90", ".o") for e in F90SRC]


F2PYOPTS = ['only:', 'extrema_compute', 'extrema_set', ":"]

SOURCES = ['extrema.f90']


def compile_ext():
    import sh
    sh.make()


def compile():
    wrapper = Extension('extrema',
                        sources=SOURCES,
                        extra_f90_compile_args=F90CFLAGS,
                        f2py_options=F2PYOPTS,
                        extra_link_args=F90LFLAGS+F90OBJ)
    setup(
        ext_modules=[wrapper]
    )


if __name__ == '__main__':
    compile_ext()
    compile()
