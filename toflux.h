
!title (brief) ! Description.

!latex \briefly{Extended description.}

!latex \calledby{\link{}}
!latex \calls{\link{}}

!latex \tableofcontents

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine toflux( zeta, lflux )
  
  use knotopt, only : zero, pi2, myid, ounit, &
                      bsnlimit, bstol, &
                      tfzeta
  
  implicit none
  
  include "mpif.h"
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  REAL                 :: zeta, lflux
  
  INTEGER              :: astat, id01ajf, LRwork, LIwork
  INTEGER, allocatable :: Iwork(:)
  REAL                 :: lowlimit, upplimit, EPSABS, EPSREL, ABSERR, dflux
  REAL   , allocatable :: Rwork(:)
  
  external             :: dflux
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LRwork = bsnlimit
  SALLOCATE(Rwork,(1:LRwork),zero)
  
  LIwork = LRwork / 4
  SALLOCATE(Iwork,(1:LIwork),zero)
  
  lowlimit = zero ; upplimit = pi2 ; EPSABS = bstol ; EPSREL = bstol
  
  tfzeta = zeta ! passed through to dflux; 04 Nov 15;
  
  lowlimit = zero ; upplimit = pi2 ; EPSABS = bstol ; EPSREL = bstol
  
  id01ajf = 1 ; call D01AJF( dflux, lowlimit, upplimit, EPSABS, EPSREL, lflux, ABSERR, Rwork(1:LRwork), LRwork, Iwork(1:LIwork), LIwork, id01ajf )
  
  DALLOCATE(Iwork)
  DALLOCATE(Rwork)

 !if( myid.eq.0 ) write(ounit,'("toflux : " 10x " : zeta ="f9.4" ; lflux ="es12.5" ;")') zeta, lflux
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  return
  
end subroutine toflux

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
REAL function dflux( tt )
  
  use knotopt, only : zero, one, myid, ncpu, &
                      Ncoils, icoil, &
                      tfzeta, &
                      bsfield
  
  use oculus, only : bs00aa
  
  implicit none
  
  include "mpif.h"
  
  REAL                 :: tt
  
  INTEGER              :: ibs00aa, ierr
  REAL                 :: teta, zeta, xx(1:3), xt(1:3), xz(1:3), lAx, lAy, lAz, Ax, Ay, Az
  
!  teta = tt ; zeta = tfzeta ; call plassf( teta, zeta, xx(1:3), xt(1:3), xz(1:3) )
!  
!  bsfield%x = xx(1)
!  bsfield%y = xx(2)
!  bsfield%z = xx(3)
!  
!  lAx = zero
!  lAy = zero
!  lAz = zero
!  
!  do icoil = 1, Ncoils ! icoil is a global variable which is passed through to iccoil.h; 11 Oct 15;
!   
!   if( myid.ne.modulo(icoil-1,ncpu) ) cycle ! parallelization loop; 11 Oct 15;
!   
!   ibs00aa = 0 ; bsfield%LB = .false. ; bsfield%LA = .true. ; call bs00aa( bsfield, ibs00aa )
!   
!   lAx = lAx + bsfield%Ax * xcu(icoil)
!   lAy = lAy + bsfield%Ay * xcu(icoil)
!   lAz = lAz + bsfield%Az * xcu(icoil)
!   
!  enddo ! end of do icoil; 11 Oct 15;
!  
!  call MPI_REDUCE( lAx, Ax, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
!  call MPI_REDUCE( lAy, Ay, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
!  call MPI_REDUCE( lAz, Az, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
!  
!  RlBCAST( Ax, 1, 0 )
!  RlBCAST( Ay, 1, 0 )
!  RlBCAST( Az, 1, 0 )
!  
!  dflux = Ax*xt(1) + Ay*xt(2) + Az*xt(3)

  dflux = zero
  
  return
  
end function dflux

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
