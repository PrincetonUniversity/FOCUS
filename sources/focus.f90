
!title (main) ! Main program.

!latex \briefly{Main program.}

!latex \calledby{\link{}}
!latex \calls{\link{initial}, \link{rdsurf}, \link{rdcoils}, 
!latex \link{saving}, \link{solvers}}.

!latex \tableofcontents

!latex  \subsection{Brief introduction}
!latex  This code is used for designing 3-D coils in stellarators, knotrans and tokamaks. 
!latex  FAMUS (Finding/Flexible Optimized Coils Using Space curves) gets rid of winding surface 
!latex  by representing 
!latex  coils using space curves (either Fundamental theorem of space curves or Fourier series or other 
!latex  representations.). And for the first time, the derivatives (both the first and the second ones)
!latex  are analytically calculated. \par
!latex  For more information, please visti \href{https://princetonuniversity.github.io/FAMUS/}
!latex  If you have any questions, please send a email to czhu@pppl.gov (or zcxiang@mail.ustc.edu.cn).

!latex  \subsection{How to execute}
!latex  A brief help message will be printed if you just type `xFAMUS --help`

!latex \subsection{Misc}
!latex \bi
!latex \item Details about the macros used are listed in \emph{Macros};
!latex \item The unit used in FAMUS is the International System of Units (SI);
!latex \item All the shared variables are defined in \link{globals};
!latex \item Variables starting with {\bf case\_xxx} are used for controling the flow direction;
!latex \item Variables starting with {\bf Isxxx} are used for logical judgement.
!latex \ei

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

PROGRAM FAMUS
  
  use globals, only: dp, ncpu, myid, ounit, ierr, astat, eunit, case_surface, case_coils, case_optimize, &
       case_postproc, xdof, time_initialize, time_optimize, time_postproc, &
       version, MPI_COMM_FAMUS
  use mpi  
  implicit none

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  INTEGER :: secs, mins, hrs
  REAL    :: tstart, tfinish ! local variables

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  call initial ! read input namelist and broadcast;

  call check_input

  select case( case_surface )

  case( 0 ) ; call fousurf   ! general format (VMEC-like) plasma boundary;
  case( 1 ) ; call rdknot    ! knototran-like plasma boundary;
 !case( 2 ) ; call readwout  ! read vmec output for plasma boundary and Boozer coordinates; for future;

  end select

  select case( case_coils )

 !case( 0 )   ; call coilpwl ! piece-wise linear; for future;
  case( 1 )   ; call rdcoils

  end select


  call packdof(xdof)  ! packdof in xdof array;

  call MPI_BARRIER( MPI_COMM_FAMUS, ierr )

  tfinish = MPI_Wtime()
  time_initialize = tfinish - tstart
  if( myid  ==  0 ) write(ounit, '(A, ES12.5, A3)') "FAMUS   : Initialization took ", time_initialize," S."

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  ! if (lsqmin) call leastsq   ! least-square minimization

  if (case_optimize /= 0) call solvers       ! call different solvers;

  call MPI_BARRIER( MPI_COMM_FAMUS, ierr )

  tstart = MPI_Wtime()
  time_optimize = tstart - tfinish
  if( myid  ==  0 ) then
     secs = int(time_optimize)
     hrs = secs/(60*60)
     mins = (secs-hrs*60*60)/60
     secs = secs-hrs*60*60-mins*60
     if(hrs>0)then
         write(ounit, '(A, 3(I6, A3))') "FAMUS   : Optimization took ",hrs," H ", mins," M ",secs," S."
     elseif(mins>0)then
         write(ounit, '(A, 2(I6, A3))') "FAMUS   : Optimization took ", mins," M ",secs," S."
     else
         write(ounit, '(A, ES12.5, A3)') "FAMUS   : Optimization took ", time_optimize," S."
     endif
  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  call unpacking(xdof)  ! unpack the optimized xdof array;

  if (myid == 0) write(ounit, *) "-----------POST-PROCESSING-----------------------------------"

  select case( case_postproc )

  case( 0 ) 
  case( 1 ) ; call diagnos ; 
  case( 2 ) ; call diagnos ; call specinp !; call saving 
 !case( 2 ) ; call saving  ; call diagnos ; call wtmgrid  ! write mgrid file;
  case( 3 ) ; call diagnos ; call poinplot ! Poincare plots; for future; 
 ! case( 3 ) ;  call poinplot ! Poincare plots; for future; 
  case( 4 ) ; call diagnos ; call boozmn ; call poinplot ! Last closed surface
  case( 5 ) ; call diagnos ; call wtmgrid  ! write mgrid file
  case(6)
     call poinplot ; call wtmgrid 
 !case( 4 ) ; call saving  ; call diagnos ; call resonant ! resonant harmonics analysis; for future; 
 case( 7 ) ; call diagnos ; call wtmgrid ; call poinplot

  end select

  call saving ! save all the outputs

  call MPI_BARRIER( MPI_COMM_FAMUS, ierr )

  tfinish = MPI_Wtime()
  time_postproc = tfinish - tstart
  if( myid  ==  0 )then
     secs = int(time_postproc)
     hrs = secs/(60*60)
     mins = (secs-hrs*60*60)/60
     secs = secs-hrs*60*60-mins*60
     if(hrs>0)then
         write(ounit, '(A, 3(I6, A3))') "FAMUS   : Post-processing took ",hrs," H ", mins," M ",secs," S."
     elseif(mins>0)then
         write(ounit, '(A, 2(I6, A3))') "FAMUS   : Post-processing took ", mins," M ",secs," S."
     else
         write(ounit, '(A, ES12.5, A3)') "FAMUS   : Post-processing took ", time_postproc," S."
     endif
  endif

  if(myid == 0) write(ounit, *) "-------------------------------------------------------------"

  call MPI_FINALIZE( ierr )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  !call cleanup
  
END PROGRAM FAMUS

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
