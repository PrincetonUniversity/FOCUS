# Synopsis

The FOCUS uses 3D curves to represent coils in fusion devices.
For a given plamsa configuration (with or without target Bn distributions), FOCUS could find optimal coils, meeting the user-specified physics and engineering constraints.

&nbsp;

# Table of contents

- [Before running](#before-running)
  * [Get the source](#get-the-source)
  * [Compile](#compile)
- [Input files](#input-files)
- [Running](#running)
- [Output files](#output-files)
- [Plotting](#plotting)
- [Documentation](#documentation)

# Before running

## Get the source
If you are fresh to GitHub, you can visit this [page](https://princetonuniversity.github.io/FOCUS/Get_the_code) to learn how to get a copy of the code.

## Compile
* Prerequisites

  The current version of the FOCUS is using intel compiler with libraries OpenMPI, hdf5-serial.
  If you are at PPPL cluster, you should add these terms in your environment variables file (like *~/.cshrc*):
  ```
  module load intel openmpi/1.8.4 hdf5-serial/1.8.5
  ```

* Makefile

  Go to the *./New/* file.
  
  All the Fortran90 sources are in \*.h files. When **make**, \*.h file will produce \*.F90 with extracted macros (*seen in [macros](https://github.com/PrincetonUniversity/FOCUS/tree/master/New/macros)*).

  There are two options vailable in [Makefile](https://github.com/PrincetonUniversity/FOCUS/tree/master/New/Makefile):
  
  **Optimized concise version (recommended)**
  ```
  make xfocus
  ```
  **Debugging version (with more stricted checkings and more error informations)**
  ```
  make dfocus
  ```
  *If you want to use GCC compiler, try `make CC=gfortran xfocus`.*
  
# Input files

FOCUS needs three input files for running, input namelist, target plasma boundary and initial coils.
Here is a brief description for these three files. 
For instance, you want to run the code with a case name of "example".

* input namelist
  
  The file **example.input** contains all the input variables for FOCUS.
  Here is an example input file for the rotating ellipse case [test.input](https://github.com/PrincetonUniversity/FOCUS/tree/master/examples/New_rotating_ellipse/test.input)
  A detailed description for these variables can be found in [initial.pdf](https://princetonuniversity.github.io/FOCUS/initial.pdf).
  These files should be placed in the same directory as the executable.
  
* target plasma boundary

  There are different options for reading the target plasma boundary.
  It's controled by *case_surface* in the input namelist. 
  
  - *case_surface = 0* : **plasma.boundary**
  
    This is for general unknotted cases, like stellarator and tokamaks. It takes VMEC-like format. 
    Detalis about the format can be seen in [rdsurf.pdf](https://princetonuniversity.github.io/FOCUS/rdsurf.pdf)
	Here is an example for the rotating ellipse case [plasma.boundary](https://github.com/PrincetonUniversity/FOCUS/tree/master/examples/New_rotating_ellipse/plasma.boundary)
    
* initial coils

  FOCUS also requires the user to provide an initial guess for the coils. This is controled by *case_init*.
  
  - *case_init = -1* : **coils.example**
  
    Read the coils data from **coils.example** and fit the coils with Fourier coefficients. 
    The format ofcoils.\* file can be seen in [VMECwiki](http://vmecwiki.pppl.wikispaces.net/MAKEGRID).
	Here is an example for the rotating ellipse case [coils.test](https://github.com/PrincetonUniversity/FOCUS/tree/master/examples/New_rotating_ellipse/coils.test)
    
  - *case_init =  0* : **example.focus**
  
    Read the coils data from **example.focus**. This file contains all the Fourier harmonics and control labels for the coils.
	Here is an example for the rotating ellipse case [test.focus](https://github.com/PrincetonUniversity/FOCUS/tree/master/examples/New_rotating_ellipse/test.focus)
    
  - *case_init =  1*
  
    Initialize *Ncoils* circular coils (r=*init_radius*, I=*init_current*) surrounding the plasma boundary.

If you want to optimize individual Bn spectrum (*weight_bharm>0* in the namelist), you may also need to provide an input file named *target.harmonics*.
Detalis about the format can be seen in [bmnharm.pdf](https://princetonuniversity.github.io/FOCUS/bmnharm.pdf)
Here is an example for the DIIID RMP coils case [target.harmonics](https://github.com/PrincetonUniversity/FOCUS/tree/master/examples/New_d3d_RMP/target.harmonics)

&nbsp;

# Running

Once you finish preparing the input files and successfully compile the code, move to your working directory, which contains the input files and the executable, and then type
```
mpirun -np 32 xfocus example
```
You may need to allocate computating cores first, e.g. try `salloc -p dawson -n 32 -t 24:00:00` at PPPL.

The code shoul print some information on the screen (or in stdout file for *sbatch*).

Use the variable *IsQuiet* in the namelist to control the details you want.

&nbsp;

# Output files

When calling [saving](https://github.com/PrincetonUniversity/FOCUS/tree/master/New/saving.h), current calculating results will be saved to local files.
The frequence of writing outputs is controlled by the variable *save_freq*.
The basic output file is written in hdf5 format.

* hdf5 file
  
  All the results can be seen in *focus_example.h5*. 
  Any hdf5 file reading functions can be used for reading the results. 
  Recommending the class *hdf5* defined in [pyfocus](https://github.com/PrincetonUniversity/FOCUS/blob/master/pyfocus/coil.py) (details can seen in [Plotting](#Plotting)).
  
If you want to see other output files, please turn on the flags in the namelist.

* save_coils

  Optimized coils will be stored in *example.focus* and *example.coils*. 
  They are FOCUS format coils and the standard coils file mentioned above, respectively. 
  So they can be used as the initial coils of next optimization.
  
* save_harmonics

  The Bn spectrum produced by current coils are stored in *example.harmonics*.
  
* save_filaments

Intermediate coil data (XYZ points in space) are store in binary file *.example.filaments.000001*. 
  So you can use them to plot coil evolution movie. (Better way is to use the data in hdf5 file.)

&nbsp;

# Plotting
There are several tools for processing the data, like a [python package](https://github.com/PrincetonUniversity/FOCUS/blob/master/pyfocus/coil.py), 
a powerful GUI interface Echidna in IDL written by Dr. Hudson, and some MATLAB scripts by Dr. Lazerson.

Here are some basic instructions for using the python package (later I will write a detailed one):

* import python package
  ```
  from pyfocus import *
  ```

* plot the plasma boundary
  ```
  plot_plasma_boundary('/your_path/plasma.boundary', 'surface3d')
  ```

* plot coils in coils.xxx file
  ```
  coil = read_coils('/your_path/coils.example')
  plot_coils(coil)
  ```

* load hdf5 file
  ```
  test = hdf5('/your_path/focus_test.h5')
  ```

* plot chi-square converging curve
  ```
  chievolve(test, 'chi')
  ```
  
* plot coil evolution movie
  ```
  coilevolve(test, delay=1000)
  ```
  

&nbsp;

# Documentations
In the source files, the comments starting with "!latex " can be exported into a tex file and generate pdf documentations.
You can find some of them in [Subroutins](https://princetonuniversity.github.io/FOCUS/subroutines).

&nbsp;

# Contact
If you have any questions, please contact Caoxiang Zhu (czhu@pppl.gov or zcxiang@mail.ustc.edu.cn).

Please note that there are no warranties! :octocat:

-----------
