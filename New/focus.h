
!title (main) ! Main program.

!latex \briefly{Main program.}

!latex \calledby{\link{}}
!latex \calls{\link{initial}, \link{surface}, \link{rdknot}, \link{rdcoils}, 
!latex \link{restart}, \link{optimizer}, \link{cleanup}}.

!latex \tableofcontents

!latex  \subsection{Brief introduction}
!latex  This code is used for designing 3-D coils in stellarators, knotrans and tokamaks. 
!latex  FOCUS (Finding/Flexible Optimized Coils Using Space curves) gets rid of winding surface 
!latex  by representing 
!latex  coils using space curves (either Fundamental theorem of space curves or Fourier series or other 
!latex  representations.). And for the first time, the derivatives (both the first and the second ones)
!latex  are analytically calculated. \par

!latex  Parts of the code were first written by 
!latex  \href{http://w3.pppl.gov/~shudson/}{\blu{Dr. Stuart R. Hudson}} in April 2016.
!latex  Then Caoxiang Zhu (CZHU) took over the whole project and it's currently under developping. 
!latex  If you have any questions, please send a email to czhu@pppl.gov (or zcxiang@mail.ustc.edu.cn).

!latex \subsection{Update Diary}
!latex 2015/10/30: Dr. Stuart Hudson wrote the code {\bf OPTIM} using 
!latex \nag{http://www.nag.co.uk/numeric/FL/manual19/pdf/E04/e04jyf_fl19.pdf}{E04JYF} to find the 
!latex optimal Fourier series for coils on a given winding surface. \par
!latex 2016/04/xx: New code {\bf KNOTOPT} that represents coils using 3D Fourier series was written.\par
!latex 2016/04/21: CZHU began to join the project and mainly took over the project.  \par
!latex 2016/xx/xx: A lot of new stuffs were added in but not well documented. \par
!latex 2016/11/01: The code was renmaed to FOCUS and a poster was presented by CZHU at the APS-DPP 
!latex meeting in San Jose, CA. \par
!latex 2017/02/15: A re-writing for debugging and better structure began by CZHU. \par
!latex 2017/04/04: The code repository was tranported to Princeton University @ GitHub/ \par
!latex 2017/05/15: Nonlinear Conjugate Gradient method was implemented. \par
!latex 2017/05/20: Truncated Newton Method with Preconditioning CG method was implemented.\par
!latex 2017/06/04: The first paper introducing FOCUS was submitted to Nuclear Fusion. \par
!latex 2017/06/07: Hybrid Newton method was implemented. \par
 
!latex \subsection{Structure of the code}
!latex \begin{tikzpicture}[node distance=2cm, auto]
!latex \node [block] (start) {Main program in \link{focus}};
!latex \node [io, below of=start] (input) {read input in \link{initial} \& allocate data in \link{dataalloc}};
!latex \node [io, below of=input] (surface) {read \& discretize surface data in \link{surface}};
!latex \node [io, below of=surface] (coils) {initialize coils data in \link{fourier}};
!latex \node [cloud, left of=coils, xshift=-4cm, yshift=1.0cm] (diagnos) 
!latex               {coils evaluation in \link{diagnos}};
!latex \node [block, below of=coils] (pack) {Packing degrees of freedom in \link{packdof}};
!latex \node [decision, below of=pack,]  (optimizer) {Optimizing in \link{solvers}};
!latex \node [block, right of=optimizer, xshift=2.5cm] (unpack) 
!latex               {unpack DOF to coils in \link{packdof}};
!latex \node [block, right of=unpack, xshift=2cm] (costfun) 
!latex               {calculate the cost functions in \link{solvers}};
!latex \node [block, below of=optimizer, yshift=-1.5cm] (postproc) {post proceeding in \link{wtmgrid}};
!latex \node [io, below of=postproc] (output) {saving all the data in \link{savedata}};
!latex \node [block, below of=output] (clean) {clean and finish in \link{cleanup}};

!latex \path [line] (start) -- (input);
!latex \path [line] (input) -- (surface);
!latex \path [line] (surface) -- (coils);
!latex \path [line, dashed] (surface) -| (diagnos);
!latex \path [line, dashed] (coils) -| (diagnos);
!latex \path [line] (coils) -- (pack);
!latex \path [line] (pack) -- (optimizer);
!latex \path [line] (optimizer) -- node {iterations} (unpack);
!latex \path [line] (unpack) -- (costfun);
!latex \path [line] (costfun) |- (pack);
!latex \path [line] (optimizer) -- node {is over} (postproc);
!latex \path [line] (postproc) -- (output);
!latex \path [line, dashed] (postproc) -| (diagnos);
!latex \path [line] (output) -- (clean);
!latex \end{tikzpicture}

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
  
  use globals, only: ncpu, myid, ounit, ierr, astat, eunit, case_surface, case_coils, &
       case_postproc, xdof, tstart, tfinish, time_initialize, time_optimize, time_postproc
  
  implicit none
  
  include "mpif.h"

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  INTEGER :: secs, mins, hrs

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  myid = 0 ; ncpu = 1

  ! MPI initialize
  call MPI_INIT( ierr )
  call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
  call MPI_COMM_SIZE( MPI_COMM_WORLD, ncpu, ierr )

  tstart =  MPI_WTIME()
  if(myid == 0) write(ounit, *) "-------------------------------------------------------------"
  if(myid == 0) write(ounit,'("focus   : Begin execution with ncpu =",i5)') ncpu
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  call initial ! read input namelist and broadcast;
  
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

  call MPI_BARRIER( MPI_COMM_WORLD, ierr )

  tfinish = MPI_Wtime()
  time_initialize = tfinish - tstart
  if( myid  ==  0 ) write(ounit, '(A30, ES12.5, A3)') "focus   : Initialization took ", time_initialize," S."

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  call solvers       ! call different solvers;

  call MPI_BARRIER( MPI_COMM_WORLD, ierr )

  tstart = MPI_Wtime()
  time_optimize = tstart - tfinish
  if( myid  ==  0 ) then
     secs = int(time_optimize)
     hrs = secs/(60*60)
     mins = (secs-hrs*60*60)/60
     secs = secs-hrs*60*60-mins*60
     if(hrs>0)then
         write(ounit, '(A30, 3(I6, A3))') "focus   : Optimization took ",hrs," H ", mins," M ",secs," S."
     elseif(mins>0)then
         write(ounit, '(A30, 2(I6, A3))') "focus   : Optimization took ", mins," M ",secs," S."
     else
         write(ounit, '(A30, ES12.5, A3)') "focus   : Optimization took ", time_optimize," S."
     endif
  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  call unpacking(xdof)  ! unpack the optimized xdof array;

  select case( case_postproc )

  case( 0 ) ; call saving
 !case( 1 ) ; call saving ; call diagnos
 !case( 2 ) ; call saving ; call diagnos ; call wtmgrid  ! write mgrid file;
 !case( 3 ) ; call saving ; call diagnos ; call poinplot ! Poincare plots; for future; 
 !case( 4 ) ; call saving ; call diagnos ; call resonant ! resonant harmonics analysis; for future; 

  end select

  call MPI_BARRIER( MPI_COMM_WORLD, ierr )

  tfinish = MPI_Wtime()
  time_postproc = tfinish - tstart
  if( myid  ==  0 )then
     secs = int(time_postproc)
     hrs = secs/(60*60)
     mins = (secs-hrs*60*60)/60
     secs = secs-hrs*60*60-mins*60
     if(hrs>0)then
         write(ounit, '(A30, 3(I6, A3))') "focus   : Post-processing took ",hrs," H ", mins," M ",secs," S."
     elseif(mins>0)then
         write(ounit, '(A30, 2(I6, A3))') "focus   : Post-processing took ", mins," M ",secs," S."
     else
         write(ounit, '(A30, ES12.5, A3)') "focus   : Post-processing took ", time_postproc," S."
     endif
  endif

  if(myid == 0) write(ounit, *) "-------------------------------------------------------------"

  call MPI_FINALIZE( ierr )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  !call cleanup
  
END PROGRAM focus

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
