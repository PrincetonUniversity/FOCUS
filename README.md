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
The directory [src](https://github.com/PrincetonUniversity/FOCUS/tree/master/src) contains all the source codes and
[Tools](https://github.com/PrincetonUniversity/FOCUS/tree/master/Tools) has some useful utilities.

### Make
All the Fortran90 sources are in *.h files. When *make*, *.h file will produce *.F90 with extracted macros (*seen in [macros](https://github.com/PrincetonUniversity/FOCUS/blob/src/master/macros)*).

There are several *make* options vailable in [Makefile](https://github.com/PrincetonUniversity/FOCUS/blob/src/master/Makefile):
* Optimized concise version (recommended)
  ```
  make xfocus DFLAGS="-D BNORM" 2>&1 | tee make.log
  ```
    *If you don't have target Bn distribution (e.g. stellarator vacuum field), you can just use "make xfocus".*
* make dfocus (debugging version)
  ```
  make dfocus DFLAGS="-check all -debug full -D DEBUG" 2>&1 | tee make.log
  ```
* produce documentated pdfs
  ```
  make pdfs
  ```
    *This will also move all pdfs to ~/w3_html/Focus/ for online viewing.*

### Run the code
Once successfully compiled the code, you will get the executable *xfocus (or dfocus)*. 
You can type
```
mpirun -np 32 xfocus suffix
```
There are three basic files needed for a run (this example can be seen in [Examples](https://github.com/PrincetonUniversity/FOCUS/tree/master/src/Examples)).

* **suffix.fo**

  The input namelist file, details can been seen in [globals](https://github.com/PrincetonUniversity/FOCUS/blob/master/src/globals.h)
  
* **plasma.boundary**

  The plasma boundary files including Fourier harmonics for the plasma surface and Bnormal distribution.
  See [surface](https://github.com/PrincetonUniversity/FOCUS/blob/master/src/surface.h).
  
* **intial coils**

  Linitialize = -1 : read *coils.suffix* and fit with Fourier series;
  
  Linitialize =  0 : read all the *.fo.coil.xxx* files in the directory;
  
  Linitialize =  N : N>1, intialize N circular coils toroidally surrounding the plasma;
  
  See [rdcoils](https://github.com/PrincetonUniversity/FOCUS/blob/master/src/rdcoils.h).
  
### Resulting files
The code will output some basic status updates to the screen when running. And it also calls [restart](https://github.com/PrincetonUniversity/FOCUS/blob/master/src/restart.h) to save a hdf5 file for each step. The *coils* and *.fo.coil.* files are also updated for each step.
  
The [coilpy](https://github.com/PrincetonUniversity/FOCUS/blob/master/Tools/coilpy.py) contains severay python functions for plotting.
