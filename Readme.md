Extrema finder
--------------

This package computes the extrema of a 3D field.

Requirements
------------
Requires:
* BLAS
* LAPACK
* FFTW3 (will be removed eventually)
* numpy
* python (of course!)
* gfortran (for other compilers, see below)

If you're using ifort or other, you have to clone this repository and modify the Makefile accordingly.


Install
-------
Once you've met all dependencies, just run:
```sh
pip install -e git+https://github.com/cphyc/extrema.git#egg=extrema
```

Example
-------
```python
from extrema import extrema
import numpy as np

# Generate a random field
field = np.random.rand(32**3).reshape(32, 32, 32)

# Mandatory: set parameters of the run
npeaks = 32**3  # max number of peaks
nbins = 32**3   # number of bins of the field
nthreads = 2    # number of threads to use
extrema.set(npeaks, nbins, nthreads)

# Find the peaks
extrema.compute(field)

# Get the data back
data = dict(
     pos=extrema.peak_pos,     # pos of the critical point
     val=extrema.peak_values,  # value at the point
     eigval=extrema.eigval,    # eigen values of the hessian
     eigvec=extrema.eigvect,   # eigen vectors of the hessian
     type=extrema.peak_types   # type of the extrema (4: maximum, 3,2:saddle point, 1:minimum)
     index=extrema.peak_index  # unravelled index of the cell containing the critical point
)

# Unravel indices
indices = np.unraval_index(data['index'], field.shape)
```


Copyright
---------
The Fortran Routines have been borrowed from Dmitri Pogosyan. Any work using this package should name him explicitly.

License
-------
This work is licensed under the CC-BY-SA license. You are allowed to copy, modify and distribute it as long as you keed the license. See more in the LICENSE file.
