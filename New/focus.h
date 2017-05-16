
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
 
!latex \subsection{Structure of the code}
!latex \begin{tikzpicture}[node distance=2cm, auto]
!latex \node [block] (start) {start from \link{focus}};
!latex \node [io, below of=start] (input) {read input namelist in \link{initial}};
!latex \node [io, below of=input] (surface) {read \& process surface data in \link{surface}};
!latex \node [io, below of=surface] (coils) {initialize coils data in \link{rdcoils}};
!latex \node [cloud, left of=coils, xshift=-4cm, yshift=1.0cm] (diagnos) 
!latex               {coils evaluation in \link{diagnos}};
!latex \node [block, below of=coils] (pack) {Packing degrees of freedom in \link{packdof}};
!latex \node [decision, below of=pack,]  (optimizer) {Optimizing in \link{optimizer}};
!latex \node [block, right of=optimizer, xshift=2.5cm] (unpack) 
!latex               {unpack DOF to coils in \link{unpacking}};
!latex \node [block, right of=unpack, xshift=2cm] (costfun) 
!latex               {calculate the cost functions in \link{costfun}};
!latex \node [block, below of=optimizer, yshift=-1.5cm] (postproc) {post proceeding in \link{postproc}};
!latex \node [io, below of=postproc] (output) {saving all the data in \link{restart}};
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
!latex \item Details about the macros are listed in Macros;
!latex \item All the shared variables are defined in \link{globals};
!latex \item Variables starting with case\_xxx are used for controling the flow direction;
!latex \item Most the data are saved in hdf5 format in \link{restart};
!latex \item Some routines are from the \oculus{subroutines}.
!latex \item The text width is set to 108, so that two windows can be dispalayed at the same time.
!latex \ei

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

PROGRAM focus
  
  use globals, only: ncpu, myid, ounit, case_surface, case_coils, case_optimizer, IsNormWeight, &
       case_postproc, xdof, tstart, tfinish, time_initialize, time_optimize, time_postproc
  
  use oculus, only : pp00aa
  
  implicit none
  
  include "mpif.h"

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  INTEGER :: ierr, astat, irestart, itmp  ! for error indicators; 2017/02/16
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
  
  call initial ! read input and broadcast;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  select case( case_surface )
     case( 0 ) ; call generic
     case( 1 ) ; call rdknot
  end select
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  select case( case_coils )
    !case( 0 )   ; call coilpwl
     case( 1 )   ; call initfou
  end select

  call packdof(xdof)
  
  tfinish = MPI_Wtime()
  time_initialize = tfinish - tstart
  if( myid  ==  0 ) write(ounit,'("focus   : Initializations take ", es23.15," seconds;")') time_initialize
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!!$  
!!$  if( Lnormalize .ne. 0 ) call weightnormalize
!!$  
!!$!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  select case( case_optimizer )
  case ( -1)      !finite-difference checking the first  derivatives;
     itmp = 1
     call AllocData(itmp)
     irestart = 0
     call fdcheck(itmp)
  case(  1)       ! differential flow;
     itmp = 1
     call AllocData(itmp)
     irestart = 1
     call descent
  case default
!!$     itmp = 1
!!$     call AllocData(itmp)
!!$     irestart = 1
!!$     call costfun(itmp)
  end select

  tstart = MPI_Wtime()
  time_optimize = tstart - tfinish
  if( myid  ==  0 ) then
     secs = int(time_optimize)
     hrs = secs/(60*60)
     mins = (secs-hrs*60*60)/60
     secs = secs-hrs*60*60-mins*60
     if(hrs>0)then
         write(ounit, *) "focus   : Optimization took ",hrs," hours, ", mins," minutes, ",secs," seconds"
     elseif(mins>0)then
         write(ounit, *) "focus   : Optimization took ", mins," minutes, ",secs," seconds"
     else
         write(ounit, *) "focus   : Optimization took ", secs," seconds;"
     endif
  endif

  call restart(irestart)
  !call identfy

!!$!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!!$
!!$  ! binary selection;
!!$  ! 0 : off
!!$  ! 1 : diagnos; 2: exinput; 4: wtmgrid
!!$  select case( Lpoincare )
!!$  case( 0  ) ; 
!!$  case( 1  ) ; call diagnos ! normal poicare plot from coils file 
!!$  case( 2  ) ; call exinput ! call Oculus:bn00aa to prepare input for free-boundary SPEC vacuum verfication; 04 Aug 16;
!!$  case( 3  ) ; call exinput
!!$             ; call diagnos
!!$  case( 4  ) ; call wtmgrid !write mgrid
!!$  case( 5  ) ; call wtmgrid
!!$             ; call diagnos
!!$  case( 6  ) ; call exinput
!!$             ; call wtmgrid
!!$  case( 7  ) ; call exinput
!!$             ; call wtmgrid
!!$             ; call diagnos
!!$  case( 8: ) ; call wtmgrid ! 8 only write mgrid; > 8 write both
!!$  end select
!!$
!!$ !call identfy
!!$
!!$  call cpu_time(tstart)
!!$  
!!$  if( myid  ==  0 ) write(ounit,'("focus   : ",10x," : post-proceed finished ; totally takes ", es23.15, " seconds;")') tstart - tfinish

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  tfinish = MPI_Wtime()
  time_postproc = tfinish - tstart
  if( myid  ==  0 )then
     secs = int(time_postproc)
     hrs = secs/(60*60)
     mins = (secs-hrs*60*60)/60
     secs = secs-hrs*60*60-mins*60
     if(hrs>0)then
        write(ounit, *) "focus   : Post-processing took ",hrs," hours, ", mins," minutes, ",secs," seconds"
     elseif(mins>0)then
        write(ounit, *) "focus   : Post-processing took ", mins," minutes, ",secs," seconds"
     elseif(secs > 0)then
        write(ounit, *) "focus   : Post-processing took ", secs," seconds;"
     else
        write(ounit,'("focus   : Post-processing took ", es10.3," seconds;")') time_postproc
     endif
  endif

  call MPI_FINALIZE( ierr )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
END PROGRAM focus

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
