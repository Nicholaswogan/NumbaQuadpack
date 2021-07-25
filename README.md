# NumbaQuadpack

`NumbaQuadpack` is a python wrapper to [cquadpack](https://github.com/ESSS/cquadpack), which is for solving integrals. Right now, this package only wraps the `dqags` method, which computes definite integrals.

For most definite integrals, with non-infinite bounds `NumbaQuadpack.dqags` does the same thing as `scipy.integrate.quad`, but NumbaQuadpack is much faster than scipy, because scipy uses the python interpreter to set up the problem. Also, `NumbaQuadpack.dqags` can be called from within numba-jited functions, and scipy's functions can not be.

## Installation
`NumbaQuadpack` will probably only work on MacOS or Linux. You must have `CMake` and a C compiler. You must also have python >3.6.0 with `numpy` and `numba`.

After satisfying the dependencies, install with the pip command below

```
python -m pip install git+git://github.com/Nicholaswogan/NumbaQuadpack.git
```

## Basic usage

```python
from NumbaQuadpack import quadpack_sig, dqags
import numpy as np
import numba as nb

@nb.cfunc(quadpack_sig)
def f(x, data):
    return x**2
funcptr = f.address

a = 0
b = 1

# computes definite integral of function f 
# from a to b
sol, abserr, success = dqags(funcptr, a, b)
```

You can pass an array data to the your integrand in the following way

```python
@nb.cfunc(quadpack_sig)
def f(x, data_):
    data = nb.carray(data_, (2,))
    return np.sum(data)*x**2
funcptr = f.address
a = 0
b = 1
data = np.array([1.0,2.0],np.float64)
sol, abserr, success = dqags(funcptr, a, b, data = data)
```

Above, `data_`, in the integrand is a pointer to an array of floats. Therefore, operations like `np.sum(data_)` or `len(data_)` will not work. However, you can use the function `nb.carray()` to make a numpy array out of the pointers. However, you can use the function `nb.carray()` to make a numpy array out of the pointers (which is what is done above). 


Also, `dqags` can be used from within jitted functions, which is not the case with `scipy.integrate.quad`.

```python
@nb.njit()
def test():
    sol, abserr, success = dqags(funcptr, a, b)
    return sol
sol = test() # this works!
```
