
!title (main) ! Main program.

!latex \briefly{Main program.}

!latex \calledby{\link{}}
!latex \calls{\link{initial}, \link{rdsurf}, \link{rdcoils}, 
!latex \link{saving}, \link{solvers}}.

!latex \tableofcontents

!latex  \subsection{Brief introduction}
!latex  This code is used for designing 3-D coils in stellarators, knotrans and tokamaks. 
!latex  FOCUS (Finding/Flexible Optimized Coils Using Space curves) gets rid of winding surface 
!latex  by representing 
!latex  coils using space curves (either Fundamental theorem of space curves or Fourier series or other 
!latex  representations.). And for the first time, the derivatives (both the first and the second ones)
!latex  are analytically calculated. \par
!latex  For more information, please visti \href{https://princetonuniversity.github.io/FOCUS/}
!latex  If you have any questions, please send a email to czhu@pppl.gov (or zcxiang@mail.ustc.edu.cn).

!latex  \subsection{How to execute}
!latex  A brief help message will be printed if you just type `xfocus --help`

!latex \subsection{Misc}
!latex \bi
!latex \item Details about the macros used are listed in \emph{Macros};
!latex \item The unit used in FOCUS is the International System of Units (SI);
!latex \item All the shared variables are defined in \link{globals};
!latex \item Variables starting with {\bf case\_xxx} are used for controling the flow direction;
!latex \item Variables starting with {\bf Isxxx} are used for logical judgement.
!latex \ei

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

PROGRAM focus

  use globals, only: dp, ncpu, myid, ounit, ierr, astat, eunit, case_surface, case_coils, case_optimize, &
       case_postproc, xdof, time_initialize, time_optimize, time_postproc, &
       version, MPI_COMM_FOCUS
  use mpi  !to enable gfortran mpi_wtime bugs; 07/20/2017
  implicit none

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  INTEGER :: secs, mins, hrs
  REAL    :: tstart, tfinish, x, y, z, B(3)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  call initial

  tstart =  MPI_WTIME()

  ! check input namelist 
  call check_input

  select case( case_coils )

 !case( 0 )   ; call coilpwl ! piece-wise linear; for future;
  case( 1 )   ; call rdcoils

  end select

  x = 5.83444249; y = 0.33921676; z = 0.06219241
  call bfield(1, x, y, z, B(1), B(2), B(3))
  if (myid==0) print *,B

  tfinish = MPI_Wtime()
  time_initialize = tfinish - tstart
  if( myid  ==  0 ) write(ounit, '(A, ES12.5, A3)') "focus   : Initialization took ", time_initialize," S."

  if (myid == 0) write(ounit, *) "-----------POST-PROCESSING-----------------------------------"

  select case( case_postproc )

  case( 0 ) ; call poinplot
  case( 1 ) ; call wtmgrid 
  case( 2 ) ; call boozmn ; call poinplot 

  end select

  call saving ! save all the outputs

  call MPI_BARRIER( MPI_COMM_FOCUS, ierr )

  tfinish = MPI_Wtime()
  time_postproc = tfinish - tstart
  if( myid  ==  0 )then
     secs = int(time_postproc)
     hrs = secs/(60*60)
     mins = (secs-hrs*60*60)/60
     secs = secs-hrs*60*60-mins*60
     if(hrs>0)then
         write(ounit, '(A, 3(I6, A3))') "focus   : Post-processing took ",hrs," H ", mins," M ",secs," S."
     elseif(mins>0)then
         write(ounit, '(A, 2(I6, A3))') "focus   : Post-processing took ", mins," M ",secs," S."
     else
         write(ounit, '(A, ES12.5, A3)') "focus   : Post-processing took ", time_postproc," S."
     endif
  endif

  if(myid == 0) write(ounit, *) "-------------------------------------------------------------"

  call MPI_FINALIZE( ierr )
  
END PROGRAM focus

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
