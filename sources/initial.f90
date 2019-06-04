
!title (initialize) ! Read input file, initialize global variables.

!latex \briefly{Reads input files broadcasts, and allocates/initializes global variables.}

!latex \calledby{\link{focus}, \link{globals}}
!latex \calls{\link{}}

!latex \section{Input namelist}
!latex  The \emph{focusin} namelist is the only input namelist needed for FOCUS running. It should be
!latex  written to the file \emph{example.input}, where `example' is the argument passed by command line. 
!latex  Here are the details for the variables. \\
!latex  \bi
!latex  \item \inputvar{IsQuiet = -1} \\
!latex    \textit{Information displayed to the user} \\
!latex    \bi \vspace{-5mm}
!latex    \item[-2:] more details \& update unconstrained cost functions;
!latex    \item[-1:] more details;
!latex    \item[0:]  essential;
!latex    \item[1:] concise.
!latex    \ei
!latex 
!latex  \item \inputvar{IsSymmetric = 0} \\
!latex    \textit{Enforce stellarator symmetry or not} \\
!latex    \bi \vspace{-5mm}
!latex    \item[0:] no stellarator symmetry enforced;
!latex    \item[1:] plasma periodicty enforced;
!latex    \item[2:] coil and plasma periodicity enforced. 
!latex    \ei
!latex 
!latex  \par \begin{tikzpicture} \draw[dashed] (0,1) -- (10,1); \end{tikzpicture}
!latex 
!latex  \item \inputvar{input\_surf = `plasma.boundary'} \\
!latex    \textit{Input file containing plasma boundary information.} 
!latex 
!latex  \item \inputvar{input\_coils = `none'} \\
!latex    \textit{Input file containing initial guess for coils (in either format).}
!latex    If it is 'none' by default, it will be updated to 'coils.example' (case\_init=-1) 
!latex       or 'example.focus' (case\_init=0).
!latex 
!latex  \item \inputvar{input\_harm = `target.harmonics'} \\
!latex    \textit{Input file containing the target harmonics for Bmn optimization.} 
!latex 
!latex  \par \begin{tikzpicture} \draw[dashed] (0,1) -- (10,1); \end{tikzpicture}
!latex 
!latex  \item \inputvar{case\_surface = 0} \\
!latex    \textit{Specify the input plasma boundary format} \\
!latex    \bi \vspace{-5mm}
!latex    \item[0:] general VMEC-like format (Rbc, Rbs, Zbc, Zbs), seen in \link{rdsurf};
!latex    \item[1:] read axis for knots, seen in \link{rdknot}; (not ready)
!latex    %\item[2:] read wout file from VMEC outputs, seen in \link{rdvmec}.
!latex    \ei
!latex 
!latex  \item \inputvar{knotsurf = 0.2} \\
!latex    \textit{Minor plasma radius for knototrans, only valid for \inputvar{case\_surface = 1}} 
!latex 
!latex  \item \inputvar{ellipticity = 0.0} \\
!latex    \textit{Ellipticity of plasma for knototrans, only valid for \inputvar{case\_surface = 1}} 
!latex 
!latex  \item \inputvar{Nteta = 64} \\
!latex    \textit{Poloidal resolution for discretizing the plasma}
!latex 
!latex  \item \inputvar{Nzeta = 64} \\
!latex    \textit{Toroidal resolution for discretizing the plasma}
!latex 
!latex  \par \begin{tikzpicture} \draw[dashed] (0,1) -- (10,1); \end{tikzpicture}
!latex 
!latex  \item \inputvar{case\_init = 0} \\
!latex    \textit{Specify the initializing method for coils, seen in \link{rdcoils}} \\
!latex    \bi \vspace{-5mm}
!latex    \item[-1:] read the standard MAKEGRID format coils from \inputvar{input_coils};
!latex    \item[0:] read FOCUS format data from \inputvar{input_coils};
!latex    \item[1:] toroidally spaced \inputvar{Ncoils} circular coils with radius of \inputvar{init\_radius};
!latex    \item[2:] toroidally spaced \inputvar{Ncoils}-1 magnetic dipoles pointing poloidallly on the toroidal surface 
!latex               with radius of \inputvar{init\_radius} and a central infinitely long current. 
!latex               Dipole magnetizations anc the central current are all set to  \inputvar{init\_current}.
!latex    \ei
!latex 
!latex  \item \inputvar{case\_coils = 1} \\
!latex    \textit{Specify representation used for the initial coils, seen in \link{rdcoils}} \\
!latex    \bi \vspace{-5mm}
!latex    \item[0:] using piecewise linear representation; (not ready)
!latex    \item[1:] using Fourier series representation;
!latex    \ei
!latex 
!latex  \item \inputvar{Ncoils = 0} \\
!latex    \textit{Number of coils initilized, only valid for \inputvar{case\_init = 1}}
!latex 
!latex  \item \inputvar{init\_current = 1.0E6} \\
!latex    \textit{Initial coil current (A), only valid for \inputvar{case\_init = 1}} 
!latex 
!latex  \item \inputvar{init\_radius = 1.0} \\
!latex    \textit{Initial coil radius (m), only valid for \inputvar{case\_init = 1}} 
!latex 
!latex  \item \inputvar{IsVaryCurrent = 1} \\
!latex    \textit{Keep coil currents fixed or not, overriden by \emph{example.focus}} \\
!latex    \bi \vspace{-5mm}
!latex    \item[0:] coil currents are fixed;
!latex    \item[1:] coil currents are free;
!latex    \ei
!latex 
!latex  \item \inputvar{IsVaryGeometry = 1} \\
!latex    \textit{Keep coil geometries fixed or not, overriden by \emph{example.focus}} \\
!latex    \bi \vspace{-5mm}
!latex    \item[0:] coil geometries are fixed;
!latex    \item[1:] coil geometries are free;
!latex    \ei
!latex 
!latex  \item \inputvar{NFcoil = 4} \\
!latex    \textit{Number of Fourier coefficients for coils, only valid for \inputvar{case\_coils = 1},
!latex           overriden by \emph{example.focus}} 
!latex 
!latex  \item \inputvar{Nseg = 128} \\
!latex    \textit{Number of segments for discritizing coils, only valid for \inputvar{case\_coils = 1},
!latex           overriden by \emph{example.focus}} 
!latex 
!latex  \par \begin{tikzpicture} \draw[dashed] (0,1) -- (10,1); \end{tikzpicture}
!latex 
!latex  \item \inputvar{IsNormalize = 1} \\
!latex    \textit{Normalizing coil parameters or not} \\
!latex    \bi \vspace{-5mm}
!latex    \item[0:] keep raw data (normalized to 1.0);
!latex    \item[1:] currents being normalized to averaged absolute current,
!latex              coil geometry parameters being normalized to major radius;
!latex    \ei
!latex 
!latex  \item \inputvar{IsNormWeight = 1} \\
!latex    \textit{each constraints normalized to initial value or not} \\
!latex    \bi \vspace{-5mm}
!latex    \item[0:] keep raw value for constraints;
!latex    \item[1:] $w = w/f_0$ weights normalized to the initial values of each constraints;
!latex    \ei
!latex 
!latex  \item \inputvar{case\_bnormal = 0} \\
!latex    \textit{ Bn error normalized to $|B|$ or not} \\
!latex    \bi \vspace{-5mm}
!latex    \item[0:] keep raw Bn error;
!latex    \item[1:] Bn residue normalized to local $|B|$;
!latex    \ei
!latex 
!latex  \item \inputvar{case\_length = 0} \\
!latex    \textit{ options for constructing coil length constraint, seen in \link{length}} \\
!latex    \bi \vspace{-5mm}
!latex    \item[1:] quadratic format, converging the target\_length;
!latex    \item[2:] exponential format, as short as possible;
!latex    \ei
!latex 
!latex  \item \inputvar{weight\_bnorm = 1.0} \\
!latex    \textit{weight for Bn error, if zero, turned off; seen in \link{bnormal}} 
!latex 
!latex  \item \inputvar{weight\_bharm = 0.0} \\
!latex    \textit{weight for Bn Fourier harmonics error, if zero, turned off; seen in \link{bmnharm}} 
!latex 
!latex  \item \inputvar{weight\_tflux = 0.0} \\
!latex    \textit{weight for toroidal flux error, if zero, turned off; seen in \link{torflux}} 
!latex 
!latex  \item \inputvar{target\_tflux = 0.0} \\
!latex    \textit{target value for the toroidal flux, if zero, automatically set to initial $\Psi_{ave}$; 
!latex    seen in \link{solvers}} 
!latex 
!latex  \item \inputvar{weight\_ttlen = 0.0} \\
!latex    \textit{weight for coils length error, if zero, turned off; seen in \link{length}} 
!latex 
!latex  \item \inputvar{target\_length = 0.0} \\
!latex    \textit{target value (or for normalization) of the coils length, if zero, automatically 
!latex    set to initial actual length; seen in \link{rdcoils}} 
!latex 
!latex  \item \inputvar{weight\_specw = 0.0} \\
!latex    \textit{weight for spectral condensation error, if zero, turned off; seen in \link{specwid}}; (not ready) 
!latex 
!latex  \item \inputvar{weight\_ccsep = 0.0} \\
!latex    \textit{weight for coil-coil separation constraint, if zero, turned off; seen in \link{coilsep}}; (not ready) 
!latex
!latex  \item \inputvar{weight\_Inorm = 1.0} \\
!latex    \textit{additional factor for normalizing currents; the larger, the more optimized for currents;
!latex     seen in \link{rdcoils}} 
!latex 
!latex  \item \inputvar{weight\_Gnorm = 1.0} \\
!latex    \textit{additional factor for normalizing geometric variables; the larger, the more optimized for
!latex     coil shapes; seen in \link{rdcoils}} 
!latex 
!latex  \par \begin{tikzpicture} \draw[dashed] (0,1) -- (10,1); \end{tikzpicture}
!latex 
!latex  \item \inputvar{case\_optimize = 1} \\
!latex    \textit{specify optimizing options.} \\
!latex    \bi \vspace{-5mm}
!latex    \item[-2:] check the 2nd derivatives; seen in\link{fdcheck}; (not ready)
!latex    \item[-1:] check the 1st derivatives; seen in\link{fdcheck};
!latex    \item[ 0:] no optimizations performed; 
!latex    \item[ 1:] optimizing with algorithms using the gradient (DF and/or CG); seen in \link{solvers};
!latex    \item[ 2:] optimizing with algorithms using the Hessian (HT and/or NT); seen in \link{solvers}; (not ready)
!latex    \ei
!latex 
!latex  \item \inputvar{exit\_tol = 1.000D-04} \\
!latex    \textit{additional creteria to judge if the cost function decreases significantly; 
!latex     if $\frac{|\chi^2_i - \chi^2_{i-5}|}{\chi^2_i} < exit\_tol$, send an exit signal; seen in \link{solvers}}
!latex 
!latex  \item \inputvar{DF\_maxiter = 0} \\
!latex    \textit{maximum iterations allowed for using Differential Flow (DF); if zero, turned of; seen in \link{descent}}
!latex 
!latex  \item \inputvar{DF\_xtol = 1.000D-08} \\
!latex    \textit{relative error for ODE solver; seen in \link{descent}}
!latex 
!latex  \item \inputvar{DF\_tausta = 0.000D+00} \\
!latex    \textit{starting value of $\tau$; usually $0.0$ is a good idea; seen in \link{descent}}
!latex 
!latex  \item \inputvar{DF\_tauend = 0.000D+00} \\
!latex    \textit{ending value of $\tau$; the larger value of $\tau_{end} - \tau_{sta}$, the more optimized; seen in \link{descent}}
!latex 
!latex  \item \inputvar{CG\_maxiter = 0} \\
!latex    \textit{maximum iterations allowed for using Conjugate Gradient (CG); if zero, turned of; seen in \link{congrad}}
!latex 
!latex  \item \inputvar{CG\_xtol = 1.000D-08} \\
!latex    \textit{the stopping criteria of finding minimum; if $|\dd{\chi^2} / \dd{\vect{X}}|< CG\_xtol$, exit the optimization; seen in  \link{congrad}};
!latex 
!latex  \item \inputvar{CG\_wolfe\_c1 = 1.000D-04} \\
!latex    \textit{c1 value in the strong wolfe condition for line search; usually $1.0\times 10^{-4}$; seen in  \link{congrad}};
!latex 
!latex  \item \inputvar{CG\_wolfe\_c2 = 0.1} \\
!latex    \textit{c2 value in the strong wolfe condition for line search; if one CG step takes too long, try to increase c2, but remember $0<c1<c2<1$; seen in \link{congrad}};
!latex 
!latex  \item \inputvar{LM\_maxiter = 0} \\
!latex    \textit{maximum iterations allowed for using Levenberg-Marquard (LM); if zero, turned of; seen in \link{lmalg}}
!latex 
!latex  \item \inputvar{LM\_xtol = 1.000D-08} \\
!latex    \textit{the stopping criteria of finding minimum; if the relative error between two consecutivec iterates is at most xtol, the optimization terminates; seen in  \link{lmalg}};
!latex 
!latex  \item \inputvar{LM\_ftol = 1.000D-08} \\
!latex    \textit{the stopping criteria of finding minimum; if both the actual and predicted relative reductions in the sum of squares are at most ftol, the optimization terminates; seen in  \link{lmalg}};
!latex 
!latex  \item \inputvar{LM\_factor = 1.000D+02} \\
!latex    \textit{factor is a positive input variable used in determining the initial step bound. this bound is set to the product of factor and the euclidean norm of diag*x if nonzero, or else to factor itself. in most cases factor should lie in the interval (0.1,100.0). 100 is a generally recommended value. seen in  \link{lmalg}};
!latex 
!latex  \par \begin{tikzpicture} \draw[dashed] (0,1) -- (10,1); \end{tikzpicture}
!latex 
!latex  \item \inputvar{case\_postproc = 1} \\
!latex    \textit{specify post-processing options.} \\
!latex    \bi \vspace{-5mm}
!latex    \item[ 0:] no extra post-processing;
!latex    \item[ 1:] evaluate the present coils for each cost functions, coil curvature, coil-coil separation, and coil-plasma separation, Bn harmonics overlap, coil importance;
!latex    \item[ 2:] diagnos; write SPEC input file;
!latex    \item[ 3:] diagnos; Field-line tracing, axis locating and iota calculating;
!latex    \item[ 4:] diagnos; Field-line tracing; Calculates Bmn coefficients in Boozer coordinates;
!latex    \ei
!latex 
!latex  \item \inputvar{update\_plasma = 0} \\
!latex    \textit{if euqals 1, write example.plasma file with updated Bn coefficients};
!latex 
!latex  \item \inputvar{pp\_phi = 0.0} \\
!latex    \textit{Toroidal angle $\phi = pp\_phi * \pi$ for filed-line tracing, axis locating, etc.}
!latex 
!latex  \item \inputvar{pp\_raxis = 0.0} \\
!latex        \inputvar{pp\_zaxis = 0.0} \\
!latex    \textit{Initial guess for axis positions (raxis, zaxis). 
!latex            If both zero, will be overide to ($\frac{r_1+r_2}{2}, \frac{z_1+z_2}{2}$), 
!latex            where $r_1 = R(0, \phi)$ , $r_2=R(\pi, \phi)$ (likewise for $z_1, z_2$.)}
!latex 
!latex  \item \inputvar{pp\_rmax = 0.0} \\
!latex         \inputvar{pp\_zmax = 0.0} \\
!latex    \textit{Upper bounds for field-line tracing. 
!latex            If both zero, will be overide to ($r_1, z_1$).}
!latex 
!latex  \item \inputvar{pp\_ns = 10} \\
!latex    \textit{Number of surfaces for filed-line tracing, axis locating, etc.
!latex            Starting points on $\phi$ will be interpolated between 
!latex            ($r_{axis}, z_{axis}$) and ($r_{max}, z_{max}$).}
!latex 
!latex  \item \inputvar{pp\_maxiter = 1000} \\
!latex    \textit{Cycles for tracing the field-line, representing the dots for each field-line in Poincare plots.}
!latex 
!latex  \item \inputvar{pp\_tol = 1.0E-6} \\
!latex    \textit{Tolerance of ODE solver used for tracing field-lines.}
!latex 
!latex  \item \inputvar{save\_freq = 1} \\
!latex    \textit{frequency for writing output files; should be positive; seen in  \link{solvers}};
!latex 
!latex  \item \inputvar{save\_coils = 0} \\
!latex    \textit{flag for indicating whether write \emph{example.focus} and \emph{example.coils}; seen in  \link{saving}};
!latex 
!latex  \item \inputvar{save\_harmonics = 0} \\
!latex    \textit{flag for indicating whether write \emph{example.harmonics}; seen in  \link{saving}};
!latex 
!latex  \item \inputvar{save\_filaments = 0} \\
!latex    \textit{flag for indicating whether write \emph{.example.filaments.xxxxxx}; seen in  \link{saving}};
!latex  \ei

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine initial

  use globals
  implicit none
  include "mpif.h"

  LOGICAL :: exist
  INTEGER :: icpu, index_dot

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  machprec = epsilon(pi)         ! get the machine precision
  sqrtmachprec = sqrt(machprec)  ! sqrt of machine precision
  vsmall = ten * machprec        ! very small number
  small = thousand * machprec    ! small number

  !-------------read input namelist----------------------------------------------------------------------
  if( myid == 0 ) then ! only the master node reads the input; 25 Mar 15;
     call getarg(1,ext) ! get argument from command line

     select case(trim(ext))
     case ( '-h', '--help' )
        write(ounit,*)'-------HELP INFORMATION--------------------------'
        write(ounit,*)' Usage: xfocus <options> input_file'
        write(ounit,*)'    <options>'
        write(ounit,*)'     --init / -i  :  Write an example input file'
        write(ounit,*)'     --help / -h  :  Output help message'
        write(ounit,*)'-------------------------------------------------'
        call MPI_ABORT( MPI_COMM_WORLD, 1, ierr )
     case ( '-i', '--init' )
        call write_focus_namelist ! in initial.h
     case default
        index_dot = INDEX(ext,'.input')
        IF (index_dot .gt. 0)  ext = ext(1:index_dot-1)
        write(ounit, '("initial : machine_prec   = ", ES12.5, " ; sqrtmachprec   = ", ES12.5   &
             & )') machprec, sqrtmachprec
#ifdef DEBUG
        write(ounit, '("DEBUG info: extension from command line is "A)') trim(ext)
#endif
     end select
  endif

  ClBCAST( ext,  100,  0 )
  inputfile = trim(ext)//".input"

  !-------------read the namelist-----------------------------------------------------------------------
  if( myid == 0 ) then
     inquire(file=trim(inputfile), EXIST=exist) ! inquire if inputfile existed;
     FATAL( initial, .not.exist, input file ext.input not provided )
  endif

  do icpu = 1, ncpu
     call MPI_BARRIER( MPI_COMM_WORLD, ierr )
     if (myid == icpu-1) then                              ! each cpu read the namelist in turn;
        open(runit, file=trim(inputfile), status="old", action='read')
        read(runit, focusin)
        close(runit)
     endif ! end of if( myid == 0 )
  enddo

  !-------------output files name ---------------------------------------------------------------------------

  hdf5file   = "focus_"//trim(ext)//".h5"
  out_focus  = trim(ext)//".focus"
  out_coils  = trim(ext)//".coils"
  out_harm   = trim(ext)//".harmonics"
  out_plasma = trim(ext)//".plasma"

  !-------------show the namelist for checking----------------------------------------------------------

  if (myid == 0) then ! Not quiet to output more informations;

     write(ounit, *) "-----------INPUT NAMELIST------------------------------------"
     write(ounit, '("initial : Read namelist focusin from : ", A)') trim(inputfile)
     write(ounit, '("        : Read plasma boundary  from : ", A)') trim(input_surf)
     if (weight_bharm > machprec) then
        write(ounit, '("        : Read Bmn harmonics    from : ", A)') trim(input_harm)
     endif

     select case( case_init )
     case(-1 )
        if (trim(input_coils) == 'none') input_coils = "coils."//trim(ext)
        inquire( file=trim(input_coils), exist=exist )
        FATAL( initial, .not.exist, coils file coils.ext not provided )
        FATAL( initial, NFcoil <= 0    , no enough harmonics )
        FATAL( initial, Nseg   <= 0    , no enough segments  )
        FATAL( initial, target_length  < zero, illegal )        
        write(ounit, '("        : Read initial coils    from : ", A, A)') trim(input_coils), '(MAKEGRID format)'
     case( 0 )
        if (trim(input_coils) == 'none') input_coils = trim(ext)//".focus"
        inquire( file=trim(input_coils), exist=exist )
        FATAL( initial, .not.exist, FOCUS coil file ext.focus not provided )
        write(ounit, '("        : Read initial dipole  from : ", A, A)') trim(input_coils), '(Parameters only)'
        inquire( file=trim(fixed_coils), exist=exist )
        FATAL( initial, .not.exist, fixed coil file ext.focus not provided )
        write(ounit, '("        : Read fixed coils    from : ", A, A)') trim(input_coils), '(FOCUS format)'
     case( 1 )
        FATAL( initial, Ncoils < 1, should provide the No. of coils)
        FATAL( initial, init_current == zero, invalid coil current)
        FATAL( initial, init_radius < zero, invalid coil radius)
        FATAL( initial, NFcoil <= 0    , no enough harmonics )
        FATAL( initial, Nseg   <= 0    , no enough segments  )
        FATAL( initial, target_length  < zero, illegal )
        if (IsQuiet < 1) write(ounit, 1000) 'case_init', case_init, 'Initialize circular coils.'
     case( 2 )
        FATAL( initial, Ncoils < 1, should provide the No. of coils)
        FATAL( initial, init_current == zero, invalid coil current)
        FATAL( initial, init_radius < zero, invalid coil radius)
        FATAL( initial, target_length  < zero, illegal )
        if (IsQuiet < 1) write(ounit, 1000) 'case_init', case_init, 'Initialize magnetic dipoles.'
     case default
        FATAL( initial, .true., selected case_init is not supported )
     end select

     select case (IsQuiet)
     case (:-2)
        write(ounit, 1000) 'IsQuiet', IsQuiet, 'Output all information.'
     case (-1 )
        write(ounit, 1000) 'IsQuiet', IsQuiet, 'Output detailed information.'
     case ( 0 )
        write(ounit, 1000) 'IsQuiet', IsQuiet, 'Output essential information.'
     case ( 1:)
        write(ounit, 1000) 'IsQuiet', IsQuiet, 'Output concise information.'
     case default
        FATAL( initial, .true., IsQuiet /= integer unspported option)
     end select
     
1000 format(8X, ": ", A15, " = ", I6, " ; ", A)        

     select case (IsSymmetric)
     case (0)
        if (IsQuiet < 0) write(ounit, 1000) 'IsSymmetric', IsSymmetric, & 
             &  'No stellarator symmetry or periodicity enforced.'
     case (1)
        if (IsQuiet < 0) write(ounit, 1000) 'IsSymmetric', IsSymmetric, &
             &  'Plasma boundary periodicity is enforced.'
        FATAL( initial, .true., This would cause unbalanced coils please use IsSymmetric=0 instead)
     case (2)
        if (IsQuiet < 0) write(ounit, 1000) 'IsSymmetric', IsSymmetric, &
             &  'Plasma boundary and coil periodicity are both enforced.'
     case default
        FATAL( initial, .true., IsSymmetric /= 0 or 2 unspported option)
     end select

     FATAL( initial, Nteta   <=   0, illegal surface resolution )
     FATAL( initial, Nzeta   <=   0, illegal surface resolution )     

     FATAL( initial, case_coils /= 1, only fourier representation is valid )
     if (IsQuiet < 0) write(ounit, 1000) 'case_coils', case_coils, 'Using Fourier series as the basic representation.'

     select case (case_optimize)
     case ( -2 )
        write(ounit, 1000) 'case_optimize', case_optimize, 'Test the 2nd derivatives.'
        FATAL( initial, .true., 2nd deriavtives are not ready.)
     case ( -1 )
        write(ounit, 1000) 'case_optimize', case_optimize, 'Test the 1st derivatives.'
     case (  0 )
        write(ounit, 1000) 'case_optimize', case_optimize, 'No optimization will be performed.'
     case (  1 )
        write(ounit, 1000) 'case_optimize', case_optimize, 'Gradient-based optimizations will be performed.'

        if (DF_maxiter > 0) then
           if (IsQuiet < 1) write(ounit, '(26X,": Differential flow will be used, maxiter=", I6)') DF_maxiter
           if (IsQuiet < 0) write(ounit, '(26X,": DF_tausta = "ES12.5" ; DF_tauend = "ES12.5 &
                &  " ; DF_xtol = "ES12.5)') DF_tausta, DF_tauend, DF_xtol
        endif

        if (CG_maxiter > 0) then        
           FATAL( Initial, CG_wolfe_c1 <= zero .or. CG_wolfe_c1 >= one, should be 0<c1<1 )
           FATAL( Initial, CG_wolfe_c2 <= zero .or. CG_wolfe_c2 >= one, should be 0<c2<1 )
           FATAL( Initial, CG_wolfe_c1 >= CG_wolfe_c2, should be c1<c2)
           if (IsQuiet < 1) write(ounit, '(26X,": Nonlinear Conjugate Gradient method will be used, maxiter="&
                &  , I6)') CG_maxiter
           if (IsQuiet < 0) write(ounit, '(26X,": CG_wolfe_c1 = "ES12.5" ; CG_wolfe_c2 = "ES12.5 &
                &  " ; CG_xtol = " ES12.5)') CG_wolfe_c1, CG_wolfe_c2, CG_xtol
        endif

     case (  2 )
        write(ounit, 1000) 'case_optimize', case_optimize, 'Hessian-based optimizations will be performed.'
        FATAL( initial, .true., Hessian-based optimizations are not ready.)

        if (DF_maxiter > 0) then
           if (IsQuiet < 1) write(ounit, '(26X,": Differential flow will be used, maxiter=", I6)') DF_maxiter
           if (IsQuiet < 0) write(ounit,'(26X,": DF_tausta = "ES12.5" ; DF_tauend = "ES12.5 &
                &  " ; DF_xtol = "ES12.5)') DF_tausta, DF_tauend, DF_xtol
        endif

        if (CG_maxiter > 0) then        
           FATAL( Initial, CG_wolfe_c1 <= zero .or. CG_wolfe_c1 >= one, should be 0<c1<1 )
           FATAL( Initial, CG_wolfe_c2 <= zero .or. CG_wolfe_c2 >= one, should be 0<c2<1 )
           FATAL( Initial, CG_wolfe_c1 >= CG_wolfe_c2, should be c1<c2)
           if (IsQuiet < 1) write(ounit, '(26X,": Nonlinear Conjugate Gradient method will be used, maxiter="&
                &  , I6)') CG_maxiter
           if (IsQuiet < 0) write(ounit,'(26X,": CG_wolfe_c1 = "ES12.5" ; CG_wolfe_c2 = "ES12.5 &
                &  " ; CG_xtol = "ES12.5)') CG_wolfe_c1, CG_wolfe_c2, CG_xtol
        endif

        if (HN_maxiter > 0) then
           if (IsQuiet < 1) write(ounit, '(26X,": Hybrid Newton method will be used, maxiter=", I6)') HN_maxiter
           if (IsQuiet < 0) write(ounit,'(26X,": HN_factor = "ES12.5" ; HN_xtol = " &
                &  ES12.5)') HN_factor, HN_xtol
        endif

        if (TN_maxiter > 0) then        
           FATAL( Initial, TN_cr <= zero .or. TN_cr > one, should be 0<cr<=1 )
           if (IsQuiet < 1) write(ounit, '(26X,": Truncated Newton method will be used, maxiter=", I6)') TN_maxiter
           if (IsQuiet < 0) write(ounit,'(26X,": TN_cr = "ES12.5" ; TN_reorder = "I1" ; TN_xtol = " &
                &  ES12.5)') TN_cr, TN_reorder, TN_xtol
        endif

     case default
        FATAL( initial, .true., selected case_optimize is not supported )
     end select

     if (case_optimize > 0) then
        FATAL( initial, DF_maxiter < 0, must be non-negative )
        FATAL( initial, CG_maxiter < 0, must be non-negative )
        FATAL( initial, HN_maxiter < 0, must be non-negative )
        FATAL( initial, TN_maxiter < 0, must be non-negative )
     endif

     select case (IsNormalize)
     case (0)
        write(ounit, 1000) 'IsNormalize', IsNormalize, 'No normalization on coil parameters.'
     case default
        write(ounit, 1000) 'IsNormalize', IsNormalize, 'Normalization on coil parameters will be performed.'
     end select

     select case (IsNormWeight)
     case (0)
        write(ounit, 1000) 'IsNormWeight', IsNormWeight, 'No normalization on weights for constraints.'
     case default
        write(ounit, 1000) 'IsNormWeight', IsNormWeight, 'Normalization on weights for constraints will be performed.'
     end select

     select case ( case_bnormal )
     case ( 0 )
        if (IsQuiet < 1) write(ounit, 1000) 'case_bnormal', case_bnormal, 'No normalization on Bn.'
     case ( 1 )
        if (IsQuiet < 1) write(ounit, 1000) 'case_bnormal', case_bnormal, 'Bn normalized to |B|.'
     case default
        FATAL( initial, .true., selected case_bnormal is not supported )
     end select

     select case ( case_length )
     case ( 1 )
        if (IsQuiet < 1) write(ounit, 1000) 'case_length', case_length, 'Quadratic format of length penalty.'
     case ( 2 )
        if (IsQuiet < 1) write(ounit, 1000) 'case_length', case_length, 'Exponential format of length penalty.'
     case default
        FATAL( initial, .true., selected case_length is not supported )
     end select

     FATAL( initial, weight_bnorm  < zero, illegal )
     FATAL( initial, weight_bharm  < zero, illegal )
     FATAL( initial, weight_tflux  < zero, illegal )
     FATAL( initial, weight_ttlen  < zero, illegal )
     FATAL( initial, weight_specw  < zero, illegal )
     FATAL( initial, weight_ccsep  < zero, illegal )
     FATAL( initial, weight_cssep  < zero, illegal )

     select case ( case_postproc )
     case ( 0 )
        if (IsQuiet < 1) write(ounit, 1000) 'case_postproc', case_postproc, 'No extra post-processings.'
     case ( 1 )
        if (IsQuiet < 1) write(ounit, 1000) 'case_postproc', case_postproc, 'Coil evaluations will be performed.'
     case ( 2 )
        if (IsQuiet < 1) write(ounit, 1000) 'case_postproc', case_postproc, & 
             &  'Coil evaluations and writing SPEC input will be performed.'
     case ( 3 )
        if (IsQuiet < 1) write(ounit, 1000) 'case_postproc', case_postproc, & 
             &  'Coil evaluations and field-line tracing will be performed.'
     case ( 4 )
        if (IsQuiet < 1) write(ounit, 1000) 'case_postproc', case_postproc, & 
             &  'Vacuum Boozer coordinates decompostion will be performed.'
     case ( 5 )
        if (IsQuiet < 1) write(ounit, 1000) 'case_postproc', case_postproc, & 
             &  'A binary mgrid file will be saved.'
     case default
        FATAL( initial, .true., selected case_postproc is not supported )
     end select

     FATAL( initial, save_freq <= 0, should not be negative )
     write(ounit, '("outputs : HDF5 outputs           are saved in : ", A)') trim(hdf5file)
     if (save_coils /= 0) then
        write(ounit, '("outputs : Optimizated coils      are saved in : ", A, " ; ", A)') &
             trim(out_focus), trim(out_coils)
     endif
     if (weight_bharm > machprec) then
        write(ounit, '("outputs : Realized Bn harmonics  are saved in : ", A)') trim(out_harm)
     endif
     if (update_plasma/=0) then
        write(ounit, '("outputs : Updated plasma boundary is saved in : ", A)') trim(out_plasma)
     endif

  endif

  FATAL( initial, ncpu >= 1000 , too macy cpus, modify nodelabel)
  write(nodelabel,'(i3.3)') myid ! nodelabel is global; 30 Oct 15;

  ! initialize iteration and total iterations;
  iout = 1 ; Nouts = 1
  if (case_optimize >0) Nouts = DF_maxiter + CG_maxiter + LM_maxiter + HN_maxiter + TN_maxiter

  discretefactor = (pi2/Nteta) * (pi2/Nzeta)

  !save weights before normalized
  tmpw_bnorm = weight_bnorm
  tmpw_bharm = weight_bharm
  tmpw_tflux = weight_tflux
  tmpt_tflux = target_tflux
  tmpw_ttlen = weight_ttlen
 !tmpw_specw = weight_specw
  tmpw_ccsep = weight_ccsep

  call MPI_BARRIER( MPI_COMM_WORLD, ierr )

  return

end subroutine initial

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

SUBROUTINE write_focus_namelist
  use globals
  implicit none
  include "mpif.h"

  LOGICAL :: exist
  CHARACTER(LEN=100) :: example = 'example.input'

  if( myid == 0 ) then
     inquire(file=trim(example), EXIST=exist) ! inquire if inputfile existed;
     FATAL( initial, exist, example input file example.input already existed )
     write(ounit, *) 'Writing an template input file in ', trim(example)
     open(wunit, file=trim(example), status='unknown', action='write')
     write(wunit, focusin)
     close(wunit)
  endif

  call MPI_BARRIER( MPI_COMM_WORLD, ierr )
  call MPI_FINALIZE( ierr )
  stop

  return
END SUBROUTINE write_focus_namelist
