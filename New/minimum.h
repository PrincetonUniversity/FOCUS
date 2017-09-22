!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

!title (minimum) ! Minimize the target function via a differential flow.

!latex \briefly{The minimization problem is solved by integrating a system of ODEs.}

!latex \calledby{\link{solvers}}
!latex \calls{\link{solvers}}

!latex \section{overview}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

subroutine minimum( Ndof, xdof )
  
  use globals          , only : ounit, myid
  
  implicit none  
  
  include "mpif.h"
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  INTEGER              :: Ndof
  REAL                 :: xdof(1:Ndof)

  INTEGER              :: Liwork, Lrwork, iwork(1:Ndof+2), ibound, ie04kzf, iuser(1:1), ierr
  REAL                 :: ff, fdof(1:Ndof), rwork(1:Ndof*(Ndof+7)), BL(1:Ndof), BU(1:Ndof), ruser(1:1)
  
  external             :: gradobj
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  Liwork = Ndof + 2
  
  Lrwork = Ndof * ( Ndof + 7 )  
  
  FATAL( minimum, Lrwork.lt.10, see NAG documentation for E04KZF )
  
  ibound = 1 ; ie04kzf = 1
  
  call E04KZF( Ndof, ibound, gradobj, BL(1:Ndof), BU(1:Ndof), xdof(1:Ndof), ff, fdof(1:Ndof), iwork(1:Liwork), Liwork, rwork(1:Lrwork), Lrwork, &
               iuser(1:1), ruser(1:1), ie04kzf )
   
  select case( ie04kzf )
  case( 0 )    ; write(ounit,'("focus   : " 10x " : ie04kzf =",i3," ; success ;           ")') ie04kzf
  case( 1 )    ; write(ounit,'("focus   : " 10x " : ie04kzf =",i3," ; input error ;       ")') ie04kzf
  case( 2 )    ; write(ounit,'("focus   : " 10x " : ie04kzf =",i3," ; large #calc. ;      ")') ie04kzf
  case( 3 )    ; write(ounit,'("focus   : " 10x " : ie04kzf =",i3," ; quasi minimum ;     ")') ie04kzf
! case( 4 )    ; write(ounit,'("focus   : " 10x " : ie04kzf =",i3," ; not used ;          ")') ie04kzf
  case( 5 )    ; write(ounit,'("focus   : " 10x " : ie04kzf =",i3," ; some doubt ;        ")') ie04kzf
  case( 6 )    ; write(ounit,'("focus   : " 10x " : ie04kzf =",i3," ; some doubt ;        ")') ie04kzf
  case( 7 )    ; write(ounit,'("focus   : " 10x " : ie04kzf =",i3," ; some doubt ;        ")') ie04kzf
  case( 8 )    ; write(ounit,'("focus   : " 10x " : ie04kzf =",i3," ; some doubt ;        ")') ie04kzf
  case( 9 )    ; write(ounit,'("focus   : " 10x " : ie04kzf =",i3," ; unbounded ;         ")') ie04kzf
  case(10 )    ; write(ounit,'("focus   : " 10x " : ie04kzf =",i3," ; gradient error ;    ")') ie04kzf
  case default ; write(ounit,'("focus   : " 10x " : ie04kzf =",i3," ; unrecognized ifail ;")') ie04kzf
  end select
  
  return
  
end subroutine minimum

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

subroutine gradobj( Ndof, xdof, ff, fdof, iuser, ruser )
  
  use globals, only : myid, ounit, tstart, totlengt, Tfluxave, target_tflux, Bdotnsqd, ffbest
  
  implicit none
  
  include "mpif.h"
  
  INTEGER   :: Ndof, iuser(1:1)
  REAL      :: xdof(1:Ndof), ff, fdof(1:Ndof ), ruser(1:1)
  
  REAL      :: ferr, tnow
  
  call dforce( Ndof, xdof(1:Ndof), ff, fdof(1:Ndof) )
  
  if( ff.lt.ffbest ) then ; ffbest = ff ; call archive
  endif
  
  ferr = sqrt( sum(fdof(1:Ndof)*fdof(1:Ndof)) / Ndof )
  
  tnow = MPI_WTIME()
  
  if( myid.eq.0 ) write(ounit,1010) tnow-tstart, "       ", totlengt(0), Tfluxave(0)-target_tflux, Bdotnsqd(0), ff, ferr
  
1010 format("gradobj : ",f10.1," : ",a7," : L =",es17.10," ; F =",es18.10," ; ":"B =",es17.10," ; O =",es17.10," ; |dO| =",es12.05," ; ":"time =",f9.2,"s ;")
  
  return
  
end subroutine gradobj

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
