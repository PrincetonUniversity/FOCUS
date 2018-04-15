
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
  use globals, only : zero, half, myid, ncpu, ounit, IsQuiet, astat, ierr, sqrtmachprec, &
        Ndof, iout, DF_tausta, DF_tauend, DF_xtol, DF_maxiter, exit_signal

  implicit none  
  include "mpif.h"

  !---------------------------------------------------------------------------------------------     
  INTEGER              :: itau, iflag, iwork(5)
  REAL                 :: t0, tau, relerr, abserr, lxdof(1:Ndof), dE(1:Ndof)
  REAL   , allocatable :: work(:)

  external             :: denergy
  !---------------------------------------------------------------------------------------------    

  call packdof(lxdof) !copy to local;

  iflag = 1 ; t0 = DF_tausta ; tau = t0
  relerr = DF_xtol ; abserr = sqrtmachprec
  SALLOCATE( work, (1:100+21*Ndof), zero )

  call denergy(tau, lxdof, dE)
  if (myid == 0) write(ounit, '("output  : "A6" : "8(A12," ; "))') "iout", "tau", "chi", "dE_norm", &
       "Bnormal", "Bmn harmonics", "tor. flux", "coil length", "c-s sep." 
  !call output(t0)

  do itau = 1, DF_maxiter
     
     tau = DF_tausta + itau * (DF_tauend - DF_tausta) / DF_maxiter

     call mpi_barrier(MPI_COMM_WORLD, ierr)
     call ode ( denergy, Ndof, lxdof, t0, tau, relerr, abserr, iflag, work, iwork )

     if ( iflag /= 2 .and. myid == 0) then
        write ( ounit, '(A,I3)' ) 'descent : ODE solver ERROR; returned IFLAG = ', iflag
        if ( IsQuiet < 0 ) then
           select case ( iflag )
           case ( 3 )
              write(ounit, '("descent : DF_xtol or abserr too small.")')
           case ( 4 )
              write(ounit, '("descent : tau not reached after 500 steps.")')
           case ( 5 )
              write(ounit, '("descent : tau not reached because equation to be stiff.")')
           case ( 6 )
              write(ounit, '("descent : INVALID input parameters.")')
           end select
        end if
        call MPI_ABORT( MPI_COMM_WORLD, 1, ierr )
     end if

     call unpacking(lxdof)
     call costfun(1)
     call output(t0)

     if ( exit_signal ) then
        if(myid .eq. 0) write(ounit, '("descent : EXITING-------No obvious change in last 5 outputs!")')
        exit         ! no obvious changes in past 5 iterations; 07/20/2017
     endif

  end do  

  DALLOCATE( work )

  return

end subroutine descent

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine denergy( tau, lxdof, dE )
  
  use globals, only :  Ndof, myid, ounit, t1E
  implicit none
  include "mpif.h"
  !---------------------------------------------------------------------------------------------      
  REAL                 :: tau, lxdof(*), dE(*)
  
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
