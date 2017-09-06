!title (solvers) ! Integrated interface for all solvers

!latex \briefly{After initializing the surface and coils, the core part is calling minization algorithms 
!latex to optimize coil parameters.}

!latex \calledby{\link{focus}}
!latex \calls{\link{bnormal}, \link{bmnharm}, \link{torflux}, \link{length}, \link{coilsep}, 
!latex         \link{descent}, \link{congrad}, \link{hybridnt}, \link{truncnt}}

!latex  \section{Cost functions}
!latex  The chi-squared optimization method is used here. The single target function is composed of 
!latex  several chosen object functions with user-supplied weights. The general formula is
!latex  \be
!latex  \ds \chi^2(\vect{X}) = \sum_j w_j \left( \frac{f_j({\vect X}) - f_{j,o}}{f_{j,o}} \right)^2 .
!latex  \ee
!latex  
!latex  Currently, we have implemented constraints on Bnormal, Bmn harmonics, toroidal fulx, coil length, 
!latex  coil-coil separation. For details, please view the docuentation of each constraint.
!latex 
!latex  \section{Normalization}
!latex  Besides the normalization terms in each constraint, like $|{\bf B}|$ in Bnormal, there is also an 
!latex  option to normalize the object function values to its initial value.
!latex  
!latex  When \inputvar{IsNormWeight = 1}, all the nonzero weights will be divided by the current object
!latex  function values. For example, in the beginning, the Bnormal error is $f_{B_0} = 0.1$ and 
!latex  input $w_B = 1.0$. Then the updated $w'_B = w_B/f_{B_0} = 10.0$, such that at every step
!latex  \be
!latex  \ds w'_B f_B = w_B \frac{f_B}{f_{B_0}} \ .
!latex  \ee
!latex  
!latex  \emph{* Please note that when writing the output file, the original weights (as same as input) 
!latex          and {\bf IsNormWeight=1} are stored. So when you restart, the updated weights could be
!latex          different.}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine solvers

  use globals, only : ierr, iout, myid, ounit, IsQuiet, IsNormWeight, Ndof, Nouts, xdof, &
       & case_optimize, DF_maxiter, CG_maxiter, HN_maxiter, TN_maxiter, coil, DoF

  implicit none

  include "mpif.h"

  REAL :: start, finish, ii, lxdof(1:Ndof), ff, gg(1:Ndof)

  if ( myid == 0 ) write(ounit, *) "-----------OPTIMIZATIONS-------------------------------------"
  
  if ( myid == 0 .and. IsQuiet < 1 ) write(ounit, '("solvers : #degrees-of-freedom =",i6," ;")') Ndof

  if ( abs(case_optimize) >= 1 ) call AllocData(1)
  if ( abs(case_optimize) >= 2 ) call AllocData(2)
  
  if ( case_optimize < 0 ) then ; call fdcheck(case_optimize) ; return ! finite difference checking derivatives;
  endif
  
  if( IsNormWeight /= 0 ) call normweight
  
  if (myid == 0 .and. IsQuiet < 0) write(ounit, *) "---------------------------------------------------"
  if (myid == 0 .and. IsQuiet < 0) write(ounit, *) " Initial status"
  if (myid == 0) write(ounit, '("output  : "A6" : "9(A12," ; "))') "iout", "mark", "chi", "dE_norm", &
       "Bnormal", "Bmn harmonics", "tor. flux", "coil length", "spectral", "c-c sep." 

  call costfun(1)

  call saveBmn    ! in bmnharm.h;

  iout = 0 ! reset output counter;

  call output(0.0)
  
  !--------------------------------DF--------------------------------------------------------------------
  
  if( DF_maxiter > 0 ) then
   
   if (myid == 0 .and. IsQuiet < 0) write(ounit, *) "---------------------------------------------------"
   if (myid == 0 .and. IsQuiet < 0) write(ounit, *) " Optimizing with Differential Flow (DF)"
   
   call MPI_BARRIER( MPI_COMM_WORLD, ierr ) ! wait all cpus;
   
   start = MPI_Wtime()
   
   call unpacking(xdof)
   call descent
   call packdof(xdof)
   
   finish = MPI_Wtime()
   
   if (myid  ==  0) write(ounit,'("solvers : DF takes ", es23.15," seconds;")') finish - start
   
  endif
  
  !--------------------------------CG--------------------------------------------------------------------
  if (CG_maxiter > 0)  then
     if (myid == 0 .and. IsQuiet < 0) write(ounit, *) "---------------------------------------------------"
     if (myid == 0 .and. IsQuiet < 0) write(ounit, *) " Optimizing with Nonlinear Conjugate Gradient (CG)"
     call MPI_BARRIER( MPI_COMM_WORLD, ierr ) ! wait all cpus;
     start = MPI_Wtime()
     call unpacking(xdof)
     call congrad
     call packdof(xdof)
     finish = MPI_Wtime()
     if (myid  ==  0) write(ounit,'("solvers : CG takes ", es23.15," seconds;")') finish - start
  endif
  
  !--------------------------------HN--------------------------------------------------------------------
  if (HN_maxiter > 0)  then
     if (myid == 0 .and. IsQuiet < 0) write(ounit, *) "---------------------------------------------------"
     if (myid == 0 .and. IsQuiet < 0) write(ounit, *) " Optimizing with Hybrid Newton Method (HN) "
     call MPI_BARRIER( MPI_COMM_WORLD, ierr ) ! wait all cpus;
     start = MPI_Wtime()
     call unpacking(xdof)
     !call hybridnt
     finish = MPI_Wtime()
     if (myid  ==  0) write(ounit,'("solvers : HN takes ", es23.15," seconds;")') finish - start
  endif
  
  !--------------------------------TN--------------------------------------------------------------------
  if (TN_maxiter > 0)  then
     if (myid == 0 .and. IsQuiet < 0) write(ounit, *) "---------------------------------------------------"
     if (myid == 0 .and. IsQuiet < 0) write(ounit, *) " Optimizing with Truncated Newton Method (TN) "
     call MPI_BARRIER( MPI_COMM_WORLD, ierr ) ! wait all cpus;
     start = MPI_Wtime()
     call unpacking(xdof)
     !call truncnt
     call packdof(xdof)
     finish = MPI_Wtime()
     if (myid  ==  0) write(ounit,'("solvers : TN takes ", es23.15," seconds;")') finish - start
  endif
  
  !------------------------------------------------------------------------------------------------------
  
  return
end subroutine solvers

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
