# Download and compile

## Get the source
FOCUS/FAMUS is now not open-sourced.
If you would like to use the code, please contact Dr. Caoxiang Zhu (czhu[at]pppl.gov) with a short justification. 

Once you have access, you can find the repository at https://github.com/PrincetonUniversity/FOCUS.
You can download the sources using the following command
```
git clone https://github.com/PrincetonUniversity/FOCUS.git
```

If you are fresh to GitHub, you can visit this [page](Get_the_code.md) to learn how to get a copy of the code.

## Compile
### Prerequisites

  The current version of FOCUS uses the following Fortran compilers/libraries: 
  
  - **[Intel](https://software.intel.com/en-us/fortran-compilers)/[GCC](https://gcc.gnu.org/wiki/GFortran) Fortran compiler** 
  - **[OpenMPI](https://www.open-mpi.org/) (for parallel computation)**
  - **[HDF5](https://support.hdfgroup.org/HDF5/) (for output)**
    
  If the versions you downloaded are not compatible, please raise an issue.

### Compiling
  
  All the Fortran90 sources are in \*.f90 files. When **make**, \*.f90 file will produce \*_m.F90 with extracted macros (*seen in [macros](https://github.com/PrincetonUniversity/FOCUS/tree/master/sources/macros)*).

  There are two options vailable in [Makefile](https://github.com/PrincetonUniversity/FOCUS/tree/master/sources/Makefile):
  
  **Optimized concise version (recommended)**
  ```
  make xfocus
  ```
  **Debugging version (with more stricted checkings and more error informations)**
  ```
  make dfocus
  ```
  *If you want to use GCC compiler, try `make CC=gfortran xfocus`.*
  
### Direct use at PPPL

  If you are using FOCUS on the PPPL cluster, you can directly use FOCUS by typing the following lines:
  ```
  module load focus
  ```
  You can load FAMUS if you type `module load focus/dipole`.

  There are several different versions availble, please view more informations by typing
  ```
  module avail focus
  module whatis focus/develop
  ```

**After a rencent update to CentOS 7, the new way to use FOCUS/FAMUS module is now*
  ````
  module load mod_focus
  module load focus
  ````

## Compile python wrapper

FOCUS also comes with a python wrapper using [f90wrap](https://github.com/jameskermode/f90wrap).
The `Makefile` and related settings can be found in `./python`.