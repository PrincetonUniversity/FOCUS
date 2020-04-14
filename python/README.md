This is the python wrapper for FOCUS.

# Dependencies
You need to install [f90wrap](https://github.com/jameskermode/f90wrap) and [mpi4py](https://mpi4py.readthedocs.io/en/stable/install.html).
You can do it via `pip`
```
pip install f90wrap mpi4py
```
Please note before install mpi4py, you might hvae to specify `mpicc` first.
More directions are referred to the official website.

# Compilation
To compile, use `make all`. Right now, I have only implemented `GCC` compiler. 

# Examples
A simple example can be found in ../examples/rotating_ellipse/test.py.
You can test it by `make test`.

# Usage
More documentation will be available in the future.
