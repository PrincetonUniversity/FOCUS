# Synopsis

FOCUS uses 3D curves to represent coils in fusion devices.
For a given plamsa configuration (with or without target Bn distributions), FOCUS could find optimal coils, meeting the user-specified physics and engineering constraints.

<img alt="FOCUS" src="http://joshburns.net/blog/wp-content/uploads/2013/08/focus.jpg" height="180">
*Photo credit: iStockphoto.com*

&nbsp;

# Table of contents

- [Before running](#before-running)
- [Input files](#input-files)
- [Running](#running)
- [Output files](#output-files)
- [Plotting](#plotting)
- [Documentations](#documentations)
- [Publications](#publications)

# Before running

## Get the source
If you are fresh to GitHub, you can visit this [page](https://princetonuniversity.github.io/FOCUS/Get_the_code) to learn how to get a copy of the code.

## Compile
* Prerequisites

  The current version of FOCUS uses the following compilers/libraries:   
  - **[Intel](https://software.intel.com/en-us/fortran-compilers)/[gfotran](https://gcc.gnu.org/wiki/GFortran) compiler** 
  - **[OpenMPI](https://www.open-mpi.org/)**
  - **[HDF5](https://support.hdfgroup.org/HDF5/)**
    
  If the versions you downloaded are not compatible, please raise an issue.

* Compiling
  
  All the Fortran90 sources are in \*.h files. When **make**, \*.h file will produce \*.F90 with extracted macros (*seen in [macros](https://github.com/PrincetonUniversity/FOCUS/tree/master/macros)*).

  There are two options vailable in [Makefile](https://github.com/PrincetonUniversity/FOCUS/tree/master/Makefile):
  
  **Optimized concise version (recommended)**
  ```
  make xfocus
  ```
  **Debugging version (with more stricted checkings and more error informations)**
  ```
  make dfocus
  ```
  *If you want to use GCC compiler, try `make CC=gfortran xfocus`.*
  
* Usage at PPPL

  If you are using PPPL cluster, you should use FOFUCS by typing the following lines:
  ```
  module use /p/focus/modules/
  module load focus
  ```
  There are several different versions availble, please view more informations by typing
  ```
  module avail focus
  module whatis focus/old
  ```
  
# Input files

FOCUS needs three input files for running, input namelist, target plasma boundary and initial coils.
Here is a brief description for these three files. 
For instance, you want to run the code with a case name of "example".

* input namelist
  
  The file **example.input** contains all the input variables for FOCUS.
  Here is an example input file for the rotating ellipse case [ellipse.input](https://github.com/PrincetonUniversity/FOCUS/tree/master/examples/rotating_ellipse/ellipse.input)
  A detailed description for these variables can be found in [initial.pdf](https://princetonuniversity.github.io/FOCUS/initial.pdf).
  These files should be placed in the same directory as the executable.
  
* target plasma boundary

  There are different options for reading the target plasma boundary.
  It's controled by *case_surface* in the input namelist. 
  
  - *case_surface = 0* : **plasma.boundary**
  
    This is for general unknotted cases, like stellarator and tokamaks. It takes VMEC-like format. 
    Detalis about the format can be seen in [rdsurf.pdf](https://princetonuniversity.github.io/FOCUS/rdsurf.pdf)
	Here is an example for the rotating ellipse case [plasma.boundary](https://github.com/PrincetonUniversity/FOCUS/tree/master/examples/rotating_ellipse/plasma.boundary)
    
* initial coils

  FOCUS also requires the user to provide an initial guess for the coils. This is controled by *case_init*.
  
  - *case_init = -1* : **coils.example**
  
    Read the coils data from **coils.example** and fit the coils with Fourier coefficients. 
    The format ofcoils.\* file can be seen in [VMECwiki](http://vmecwiki.pppl.wikispaces.net/MAKEGRID).
	Here is an example for the rotating ellipse case [coils.ellipse](https://github.com/PrincetonUniversity/FOCUS/tree/master/examples/rotating_ellipse/coils.ellipse)
    
  - *case_init =  0* : **example.focus**
  
    Read the coils data from **example.focus**. This file contains all the Fourier harmonics and control labels for the coils.
	Here is an example for the rotating ellipse case [ellipse.focus](https://github.com/PrincetonUniversity/FOCUS/tree/master/examples/rotating_ellipse/elipse.focus)
    
  - *case_init =  1*
  
    Initialize *Ncoils* circular coils (r=*init_radius*, I=*init_current*) surrounding the plasma boundary.

If you want to optimize individual Bn spectrum (*weight_bharm>0* in the namelist), you may also need to provide an input file named *target.harmonics*.
Detalis about the format can be seen in [bmnharm.pdf](https://princetonuniversity.github.io/FOCUS/bmnharm.pdf)
Here is an example for the DIIID RMP coils case [target.harmonics](https://github.com/PrincetonUniversity/FOCUS/tree/master/examples/d3d_RMP/target.harmonics)

&nbsp;

# Running

Once you finish preparing the input files and successfully compile the code, move to your working directory, which contains the input files and the executable, and then type
```
mpirun -np 32 xfocus example
```
You may need to allocate computating cores first, e.g. try `salloc -p dawson -n 32 -t 24:00:00` at PPPL, or use a *sbatch* command.

The code shoul print some information on the screen (or in stdout file for *sbatch*).

Use the variable *IsQuiet* in the namelist to control the details you want.

&nbsp;

# Output files

Once calling [saving](https://github.com/PrincetonUniversity/FOCUS/tree/master/saving.h), present calculating results will be saved to local files.
The frequency of writing outputs is controlled by the variable *save_freq*.
The basic output file is written in hdf5 format.

* hdf5 file
  
  All the results can be seen in *focus_example.h5*. 
  Any hdf5 file reading functions can be used for reading the results. 
  
If you want to view other output files, please turn on the flags in the namelist.

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
There are several tools for processing the data.
  - **python:** There is an python package written by Dr. Caoxiang Zhu, using [Matplotlib](https://matplotlib.org/) and [Mayavi](http://docs.enthought.com/mayavi/mayavi/).
  - **OMFIT:** In the [OMFIT](https://gafusion.github.io/OMFIT-source/) frame, there is a powerful module [*focus*](https://docs.google.com/document/d/1aGpRUMpYxBJmQXfkOK2OFMZ4P0Mn5H1_OSAjJVvxpHU/edit?ts=5ad10d28#) managed by Dr. Nikolas Logan.
  - **IDL:** There is a GUI interface Echidna in IDL written by Dr. Stuart Hudson.
  - **Matlab:** There are also some MATLAB scripts by Dr. Lazerson.  

&nbsp;

# Documentations
In the source files, the comments starting with "!latex " can be exported into a tex file and generate pdf documentations.
You can find some of them in [Subroutines](https://princetonuniversity.github.io/FOCUS/subroutines).

&nbsp;

# Publications
The first paper introducing FOCUS is [C. Zhu, S.R. Hudson, Y. Song, and Y. Wan, Nuclear Fusion 58, 16008 (2018)](http://iopscience.iop.org/article/10.1088/1741-4326/aa8e0a/). For a full list of publications and presentations, please view [FOCUS Publications](https://princetonuniversity.github.io/FOCUS/Get_the_code).

# Contact
If you have any questions, please contact Dr. Caoxiang Zhu (czhu@pppl.gov or caoxiangzhu@gmail.com).

-----------