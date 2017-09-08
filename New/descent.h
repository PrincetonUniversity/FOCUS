!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

!title (descent) ! Minimize the target function via a differential flow.

!latex \briefly{The minimization problem is solved by integrating a system of ODEs.}

!latex \calledby{\link{solvers}}
!latex \calls{\link{solvers}}

!latex \section{overview}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

subroutine descent( Ndof, lxdof )
  
  use globals, only : zero, half, sqrtmachprec, myid, ncpu, ounit, &
                      DF_tausta, DF_xtol

  implicit none  
  
  include "mpif.h"
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  INTEGER              :: Ndof
  REAL                 :: lxdof(1:Ndof)
  
  INTEGER              :: ierr, astat, iflag, iwork(1:5)
  REAL                 :: tstart, tend, relerr, abserr
  REAL   , allocatable :: work(:)
  
  external             :: descentodes
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
!FATAL( descent, .true., under reconstruction )
  
!packorunpack = 'U' ; call packdof( packorunpack )
  
  FATAL( descent, DF_tausta.le.zero, illegal integration range )
  FATAL( descent, DF_xtol  .le.zero, illegal integration tolerance )
  
  iflag = 1 ; tstart = zero ; tend = DF_tausta ; relerr = DF_xtol ; abserr = sqrtmachprec
  
  SALLOCATE( work, (1:100+21*Ndof), zero )
  
!  call denergy( tau, lxdof, dE )
  
! if (myid == 0) write(ounit, '("output  : "A6" : "9(A12," ; "))') "iout", "tau", "chi", "dE_norm", &
!      "Bnormal", "Bmn harmonics", "tor. flux", "coil length", "spectral", "c-c sep." 
!call output(t0)
  
! do itau = 1, DF_maxiter
  
!  tau = DF_tausta + itau * (DF_tauend - DF_tausta) / DF_maxiter
  
!  call mpi_barrier(MPI_COMM_WORLD, ierr)

  write(ounit,'("descent : " 10x " : calling ode:descentodes ;")')
  
  call ode( descentodes, Ndof, lxdof(1:Ndof), tstart, tend, relerr, abserr, iflag, work(1:100+21*Ndof), iwork(1:5) )
  
  write(ounit,'("descent : " 10x " : called  ode/descentodes ; iflag ="i3" ;")') iflag
  
!   if ( iflag /= 2 .and. myid == 0) then
!    write ( ounit, '(A,I3)' ) 'descent : ODE solver ERROR; returned IFLAG = ', iflag
!    if ( IsQuiet < 0 ) then
!     select case ( iflag )
!     case ( 3 )
!      write(ounit, '("descent : DF_xtol or abserr too small.")')
!     case ( 4 )
!      write(ounit, '("descent : tau not reached after 500 steps.")')
!     case ( 5 )
!      write(ounit, '("descent : tau not reached because equation to be stiff.")')
!     case ( 6 )
!      write(ounit, '("descent : INVALID input parameters.")')
!     end select
!    end if
!    call MPI_ABORT( MPI_COMM_WORLD, 1, ierr )
!   endif
  
!  call unpacking(lxdof)
!  call costfun(1)
!  call output(t0)
  
!   if ( exit_signal ) then
!    if(myid .eq. 0) write(ounit, '("descent : EXITING-------No obvious change in last 5 outputs!")')
!    exit         ! no obvious changes in past 5 iterations; 07/20/2017
!   endif
  
!  enddo ! end of do itau; 
  
  DALLOCATE( work )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  return
  
end subroutine descent

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

subroutine descentodes( tau, lxdof, lfdof )
  
  use globals, only : myid, ounit, Ndof, xdof, myid, Ltflux, Llength
  
  implicit none
  
  include "mpif.h"
  
  INTEGER              :: ierr, ideriv
  REAL                 :: tau, lxdof(*), lfdof(*)
  CHARACTER            :: packorunpack
  
  xdof(1:Ndof) = lxdof(1:Ndof)
  
  packorunpack = 'U' ; call packdof( packorunpack )
  
  ideriv = 1
  
  ;             call bnormal( ideriv )
  if( Ltflux  ) call torflux( ideriv )
  if( Llength ) call  length( ideriv )
  
  FATAL( descent, .true., under reconstruction )
  
  call mpi_barrier(MPI_COMM_WORLD, ierr)
  
  return
  
end subroutine descentodes

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

