
!title (descent) ! Minimize the target function via a differential flow.

!latex \briefly{The minimization problem is solved by integrating a system of ODEs.}

!latex \calledby{\link{solvers}}
!latex \calls{\link{solvers}}

!latex \section{overview}
!latex \bi \vspace{5mm}
!latex \item[1.] The evolution is described mathematically as a system of coupled, first-order equations:
!latex           \be \frac{\partial {\bf x}}{\partial \tau} = - \frac{\partial E}{\partial {\bf x}},
!latex           \ee
!latex           where $\tau$ is an artifical time, ${\bf x}$ is the free coil parameters.
!latex \item[2.] The integration was originally performed using 
!latex            \nag{http://www.nag.co.uk/numeric/FL/manual19/pdf/D02/d02bjf_fl19.pdf}{D02BJF}. 
!latex \item[3.] It's now implemented with 
!latex           \href{http://people.sc.fsu.edu/~jburkardt/f77_src/odepack/odepack.html}{\blu{ODEPAC}}.
!latex            And the source code is in \code{ode.F90}.
!latex \item[4.] The minimization is controlled by \inputvar{DF\_tausta}, \inputvar{DF\_tauend}, 
!latex           \inputvar{DF\_xtol} and \inputvar{DF\_maxiter}.
!latex \ei


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine descent
  !---------------------------------------------------------------------------------------------
  ! Using the differential flow to optimize the coils;
  ! DATE: 2017/04/05
  !---------------------------------------------------------------------------------------------    
  use globals, only : zero, half, myid, ncpu, ounit, astat, ierr, sqrtmachprec, &
        Ndof, iout, DF_tausta, DF_tauend, DF_xtol, DF_maxiter

  implicit none  
  include "mpif.h"

  !---------------------------------------------------------------------------------------------     
  INTEGER              :: itau, iflag, iwork(5)
  REAL                 :: tau, relerr, abserr, lxdof(1:Ndof), dE(1:Ndof)
  REAL   , allocatable :: work(:)

  external             :: denergy
  !---------------------------------------------------------------------------------------------    

  call packdof(lxdof) !copy to local;

  iflag = 1 ; tau = DF_tausta
  relerr = DF_xtol ; abserr = sqrtmachprec
  SALLOCATE( work, (1:100+21*Ndof), zero )

  call denergy(tau, lxdof, dE)
  if (myid == 0) write(ounit, '("output  : "A6" : "9(A12," ; "))') "iout", "tau", "chi", "dE_norm", &
       "Bnormal", "Bmn harmonics", "toroidal flux", "coil length", "spectral", "c-c separation" 
  call output(DF_tausta)

  do itau = 1, DF_maxiter

    tau = DF_tausta + itau * (DF_tauend - DF_tausta) / DF_maxiter
    call ode ( denergy, Ndof, dE, tau, DF_tauend, relerr, abserr, iflag, work, iwork )
    call output(tau)

  end do  

  call unpacking(lxdof)

  call mpi_barrier(MPI_COMM_WORLD, ierr)

  DALLOCATE( work )

  return

end subroutine descent

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine denergy( tau, lxdof, dE )
  
  use globals, only :  Ndof, myid, ounit, t1E
  implicit none
  include "mpif.h"
  !---------------------------------------------------------------------------------------------      
  REAL                 :: tau, lxdof(1:Ndof), dE(1:Ndof)
  
  INTEGER              :: iorder
  external             :: unpacking, costfun
  !---------------------------------------------------------------------------------------------    

  call unpacking(lxdof(1:Ndof))

  iorder = 1
  call costfun(iorder)
 
  dE(1:Ndof) = - t1E(1:Ndof)

  return
end subroutine denergy


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
