# FOCUS (Flexible Optimized Coils Using Space curves)

## Synopsis
The FOCUS is mainly written by Caoxiang Zhu and Dr. Stuart Hudson. 
The code uses 3D curves to represent coils in fusion devices.
For a given plamsa configuration (with or without target Bn distributions), FOCUS could find optimal coils, meeting the user-specified physics and engineering constraints.

## Getting Started

### Prerequisites
The current version of the FOCUS is using intel compiler with libraries like NAG, OPENMPI, hdf5-serial and OCULUS.
If you are at PPPL cluster, you should add these terms in your environment variables file (like *~/.cshrc*):
```
module load intel nag/mark22 openmpi/1.8.4 hdf5-serial/1.8.5
setenv OCULUS  /u/shudson/Oculus
```
*The OCULUS library is mainly for post-proceedings. If you cannot access OCULUS, try to comment out the related parts.*

### Source code
To get the source code, please type:
```
git clone https://github.com/PrincetonUniversity/FOCUS/
```
The directory [Old](https://github.com/PrincetonUniversity/FOCUS/tree/master/Old) contains all the sources for existing version of FOCUS, [New](https://github.com/PrincetonUniversity/FOCUS/tree/master/New) has the under-developping new version.
When finishing the new code developping, the *Old* directory would probably be deleted.

There also some python utilities in [pyfocus](https://github.com/PrincetonUniversity/FOCUS/tree/master/pyfocus) and [examples](https://github.com/PrincetonUniversity/FOCUS/tree/master/examples) contains some example cases.

**The below introduction is for the old version FOCUS.**

### Make
All the Fortran90 sources are in *.h files. When *make*, *.h file will produce *.F90 with extracted macros (*seen in [macros](https://github.com/PrincetonUniversity/FOCUS/tree/master/Old/macros)*).

There are several *make* options vailable in [Makefile](https://github.com/PrincetonUniversity/FOCUS/tree/master/Old/Makefile):
* Optimized concise version (recommended)
  ```
  make xfocus [2>&1 | tee make.log]
  ```
* make dfocus (debugging version, with more stricted checkings)
  ```
  make dfocus [2>&1 | tee make.log]
  ```
* produce documentated pdfs
  ```
  make pdfs
  ```
    *This will produce pdf documentations.*
    
### Run the code
Once successfully compiled the code, you will get the executable *xfocus (or dfocus)*. 
You can type
```
mpirun -np 32 xfocus <suffix>
```
There are three basic inputs needed for running:

* **<suffix>.input**

  The input namelist file, details can been seen in [globals](https://github.com/PrincetonUniversity/FOCUS/tree/master/Old/globals.h)
  
* **plasma.boundary**

  The plasma boundary files including Fourier harmonics for the plasma surface and Bnormal distribution.
  See [surface](https://github.com/PrincetonUniversity/FOCUS/tree/master/Old/surface.h).
  
* **intial coils**

  Linitialize = -1 : read *<suffix>.coils* and fit with Fourier series;
  
  Linitialize =  0 : read all the *<suffix>.focus* files in the directory;
  
  Linitialize =  N : N>1, intialize N circular coils toroidally surrounding the plasma;
  
  See [rdcoils](https://github.com/PrincetonUniversity/FOCUS/tree/master/Old/rdcoils.h).
  
### Resulting files
For outputs, the code will provide some basic screen outputings during the running. And it also calls [restart](https://github.com/PrincetonUniversity/FOCUS/tree/master/Old/restart.h) to save a hdf5 file for each step. The *coils* and *.fo.coil.* files are also updated for each step.
