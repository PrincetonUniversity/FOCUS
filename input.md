# Input files

FOCUS needs three input files for running, input namelist, target plasma boundary and initial coils.
Here is a brief description for these three files. 
For instance, you want to run the code with a case name of "example".

* input namelist
  
  The file **example.input** contains all the input variables for FOCUS.
  Here is an example input file for the rotating ellipse case [ellipse.input](https://github.com/PrincetonUniversity/FOCUS/tree/master/examples/rotating_ellipse/ellipse.input)
  A detailed description for these variables can be found in [initial.pdf](https://princetonuniversity.github.io/FOCUS/initial.pdf).
  These files should be placed in the same directory as the
  executable.
  To obtain an template of the latest input namelist, you can just run
  `xfocus -i`.
  
* target plasma boundary

  There are different options for reading the target plasma boundary.
  It's controled by *case_surface* in the input namelist. 
  
  - *case_surface = 0* : **input_surf** (default: 'plasma.boundary')
  
    This is for general unknotted cases, like stellarator and tokamaks. It takes VMEC-like format. 
    Detalis about the format can be seen in [rdsurf.pdf](https://princetonuniversity.github.io/FOCUS/rdsurf.pdf)
    Here is an example for the rotating ellipse case [plasma.boundary](misc/plasma.boundary)
    
    For more information about preparing FOCUS boundary from VMEC and BNORM, please view [here](notes/Coil_design_codes_benchmark.html).
    
* initial coils

  FOCUS also requires the user to provide an initial guess for the coils. This is controled by *case_init*.
  
  - *case_init = -1* :  **input_coils** (default: 'coils.example') 
  
    Read the coils data from **coils.example** and fit the coils with Fourier coefficients. 
    The format ofcoils.\* file can be seen in [VMECwiki](https://princetonuniversity.github.io/STELLOPT/MAKEGRID).
	Here is an example for the rotating ellipse case [coils.ellipse](misc/ellipse.coils)
    
  - *case_init =  0* : **input_coils** (default: 'example.focus')
  
    Read the coils data from **example.focus**. This file contains all the Fourier harmonics and control labels for the coils.
	Here is an example for the rotating ellipse case [ellipse.focus](misc/ellipse.focus)
    
  - *case_init =  1*
  
    Initialize *Ncoils* circular coils (r=*init_radius*, I=*init_current*) surrounding the plasma boundary.

  - *case_init =  2*
  
    Initialize *Ncoils-1* magnetic dipoles (r=*init_radius*, I=*init_current*) surrounding the plasma boundary plus one central current.

If you want to optimize individual Bn spectrum (*weight_bharm>0* in the namelist), you may also need to provide an input file named by **input_harm** (default: 'target.harmonics').
Detalis about the format can be seen in [bmnharm.pdf](bmnharm.pdf)
Here is an example for the DIIID RMP coils case [target.harmonics](misc/target.harmonics)

&nbsp;

# Objective functions

There are multiple objective functions available in FOCUS.
If the weight of each function is greater than 0, then it is turned on.
The overall target function will be the weighted summation of each individual one.
Here is a list of input variable related to objective functions that should be stored in `*.input`. 
For more details of each function, users are suggested to check [subroutines](subroutines.md).

```fortran
  ! Normal field error
  REAL                 :: weight_bnorm   =   0.000D+00
  INTEGER              :: case_bnormal   =   0
  ! Bmn resonant harmonics
  REAL                 :: weight_bharm   =   0.000D+00
  INTEGER              :: bharm_jsurf    =   0
  CHARACTER(100)       :: input_harm     = 'target.harmonics' ! input target harmonics file
  ! toroidal flux
  REAL                 :: weight_tflux   =   0.000D+00
  REAL                 :: target_tflux   =   0.000D+00
  ! coil length
  REAL                 :: weight_ttlen   =   0.000D+00
  INTEGER              :: case_length    =   1
  REAL                 :: target_length  =   0.000D+00
  REAL                 :: length_delta   =   0.000D+00
  ! coil curvature
  REAL                 :: weight_curv    =   0.000D+00
  INTEGER              :: case_curv      =   4 
  INTEGER              :: penfun_curv    =   1
  REAL                 :: curv_alpha     =   0.000D+00
  REAL                 :: curv_beta      =   2.000D+00
  REAL                 :: curv_gamma     =   2.000D+00
  REAL                 :: curv_sigma     =   0.000D+00
  REAL                 :: curv_k0        =   1.000D+01
  REAL                 :: curv_k1        =   0.000D+00
  INTEGER              :: curv_k1len     =   0
  ! coil torsion
  REAL                 :: weight_tors    =   0.000D+00
  INTEGER              :: case_tors      =   1
  INTEGER              :: penfun_tors    =   1
  REAL                 :: tors0          =   0.000D+00
  REAL                 :: tors_alpha     =   1.000D+00
  REAL                 :: tors_beta      =   2.000D+00
  REAL                 :: tors_gamma     =   1.000D+00
  ! coil nissin complexity
  REAL                 :: weight_nissin  =   0.000D+00
  INTEGER              :: penfun_nissin  =   1
  REAL                 :: nissin0        =   0.000D+00
  REAL                 :: nissin_alpha   =   1.000D+00
  REAL                 :: nissin_beta    =   2.000D+00
  REAL                 :: nissin_gamma   =   2.000D+00
  REAL                 :: nissin_sigma   =   0.000D+00
  ! coil-surface separation
  REAL                 :: weight_cssep   =   0.000D+00
  REAL                 :: cssep_factor   =   4.000D+00 
  CHARACTER(100)       :: limiter_surf   = 'none'             ! limiter surface
  ! coil-coil separation
  REAL                 :: weight_ccsep   =   0.000D+00
  INTEGER              :: penfun_ccsep   =   1
  INTEGER              :: ccsep_skip     =   0
  REAL                 :: r_delta        =   1.000D-01
  REAL                 :: ccsep_alpha    =   1.000D+01
  REAL                 :: ccsep_beta     =   2.000D+00
```

# Running

Once you finish preparing the input files and successfully compile the code, move to your working directory, which contains the input files and the executable, and then type
```
mpirun -np 32 xfocus example
```
You may need to allocate computating cores first, e.g. try `srun -n 32 -t 12:00:00 --mem 2Gb xfocus example` at PPPL, or use a *sbatch* command.

The code should print some information on the screen (or in stdout file for *sbatch*).

Use the variable *IsQuiet* in the namelist to control the details you want.

To get a brief help message, please type
```
xfocus --help  (or xfocus -h)
```

&nbsp;
