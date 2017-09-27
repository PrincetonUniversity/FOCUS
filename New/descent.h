!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

!title (descent) ! Minimize the target function via a differential flow.

!latex \briefly{The minimization problem is solved by integrating a system of ODEs.}

!latex \calledby{\link{solvers}}
!latex \calls{\link{solvers}}

!latex \section{overview}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

module descentparameters
  
  INTEGER :: Ndegreeoffreedom ! do you know a better way of passing this through to subroutine:odes ;
  
end module descentparameters

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

subroutine descent( Ndof, xdof )
  
  use globals          , only : zero, one, six, half, sqrtmachprec, small, myid, ncpu, ounit, tstart, &
                               !tauend, Ntauout, tautol, &
                                NRKsave, NRKstep, RKstep, &
                                totlengt, Tfluxave, Bdotnsqd, &
                                target_length, target_tflux, converged, &
                                fforig, ffbest
  use descentparameters

  implicit none
  
  include "mpif.h"
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  INTEGER              :: Ndof
  REAL                 :: xdof(1:Ndof)

  INTEGER              :: irksave, irkstep
  REAL                 :: ffold, ydof(1:Ndof), fk(1:Ndof,1:4) ! Runge-Kutta ;

  REAL                 :: fjac(1:Ndof,1:Ndof), RR(1:Ndof*(Ndof+1)/2), QTF(1:Ndof), Wk(1:Ndof,1:4), xtol, epsfcn, factor
  INTEGER              :: ML, MU, MODE, Ldfjac, irevcm, ic05ndf, LR

! INTEGER              :: Mdof
! REAL                 :: qdof(1:2*Ndof)
  
  INTEGER, parameter   :: Liwork = 5
  INTEGER              :: ierr, astat, iflag, Lrwork, iwork(1:Liwork), iter
  REAL                 :: tnow, told, taustart, ltauend, relerr, abserr, ff, fdof(1:Ndof), ferr, fold
  REAL   , allocatable :: rwork(:)
  
  external             :: xodes! qodes
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
! Ndegreeoffreedom = Ndof ! this needs to be passed through to subroutines:odes ;
  
! Mdof = 2 * Ndof
! qdof(     1:Ndof) =   xdof(1:Ndof) ! position ;
! qdof(Ndof+1:Mdof) = - fdof(1:Ndof) ! velocity ;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
! Lrwork = 100 + ( 21 * Ndof )
! Lrwork = 100 + ( 21 * Mdof )
  
! SALLOCATE( rwork, (1:Lrwork), zero )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  tnow = MPI_WTIME()
  
  if( myid.eq.0 ) write(ounit,1000) tnow-tstart, NRKsave, NRKstep, RKstep, converged
1000 format("descent : ",f10.1," : Runge-Kutta ; NRKsave (outer loop) =",i9," ; NRKstep (inner loop) =",i5," ; RKstep =",es12.5," ; converged =",es13.5," ;")
  
  told = tnow
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  call dforce( Ndof, xdof(1:Ndof), ff, fdof(1:Ndof) ) ; fk(1:Ndof,1) = - fdof(1:Ndof)
  
  ferr = sqrt( sum(fdof(1:Ndof)*fdof(1:Ndof)) / Ndof )

  ffold = ff ! initialize;
  
  do irksave = 1, NRKsave ! allows intermediate output & archive;
   
   do irkstep = 1, NRKstep ! allows intermediate output & archive;
    
    ydof(1:Ndof) = xdof(1:Ndof) + fk(1:Ndof,1) * RKstep * half ; call dforce( Ndof, ydof(1:Ndof), ff, fdof(1:Ndof) ) ; fk(1:Ndof,2) = - fdof(1:Ndof)
    ydof(1:Ndof) = xdof(1:Ndof) + fk(1:Ndof,2) * RKstep * half ; call dforce( Ndof, ydof(1:Ndof), ff, fdof(1:Ndof) ) ; fk(1:Ndof,3) = - fdof(1:Ndof)
    ydof(1:Ndof) = xdof(1:Ndof) + fk(1:Ndof,3) * RKstep        ; call dforce( Ndof, ydof(1:Ndof), ff, fdof(1:Ndof) ) ; fk(1:Ndof,4) = - fdof(1:Ndof)
    
    xdof(1:Ndof) = xdof(1:Ndof) + ( fk(1:Ndof,1) + 2 * fk(1:Ndof,2) + 2 * fk(1:Ndof,3) + fk(1:Ndof,4) ) * RKstep / six ! Runge-Kutta step;
    
    call dforce( Ndof, xdof(1:Ndof), ff, fdof(1:Ndof) ) ; fk(1:Ndof,1) = - fdof(1:Ndof)

    ferr = sqrt( sum(fdof(1:Ndof)*fdof(1:Ndof)) / Ndof )
    
    if( ff.gt.ffold ) then
     if( myid.eq.0 ) write(ounit,1010) tnow-tstart, "difference", totlengt(0)-target_length, Tfluxave(0)-target_tflux, Bdotnsqd(0), ff, ferr, tnow-told
     return
    endif
    
    ffold = ff

   enddo ! end of do irkstep ;
      
   call archive( Ndof, xdof(1:Ndof), ferr )
   
   tnow = MPI_WTIME()
   
   if( myid.eq.0 ) write(ounit,1010) tnow-tstart, "difference", totlengt(0)-target_length, Tfluxave(0)-target_tflux, Bdotnsqd(0), ff, ferr, tnow-told
   
   told = tnow

   if( ferr.lt.converged ) exit
    
  enddo ! end of do irksave ;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
!  if( myid.eq.0 ) write(ounit,'("descent : ",f10.1," : calling C05NDF ;")') tnow-tstart
!  
!  xtol = sqrtmachprec ; Ml = Ndof-1  ; Mu = Ndof-1 ; epsfcn = sqrtmachprec ; mode = 1 ; factor = 1.0E-01 ; Ldfjac = Ndof ; LR = Ndof * ( Ndof + 1 ) / 2
!  
!  irevcm = 0 ; ic05ndf = 1
!  
!  do
!   
!   call C05NDF( irevcm, Ndof, xdof(1:Ndof), fdof(1:Ndof), xtol, Ml, Mu, epsfcn, fk(1:Ndof,1), mode, factor, fjac, Ldfjac, RR, LR, QTF, Wk, ic05ndf )
!   
!   select case( irevcm )
!    
!   case( 0 ) ! final exit;
!    
!    write(ounit,'("descent : " 10x " : finished ; ic05ndf ="i3" ;")') 
!    
!    exit
!    
!   case( 1 ) ! no action required;
!    
!    ferr = sqrt( sum(fdof(1:Ndof)*fdof(1:Ndof)) / Ndof )
!    
!    call archive( Ndof, xdof(1:Ndof), ferr )
!    
!    tnow = MPI_WTIME()
!    
!    if( myid.eq.0 ) write(ounit,1010) tnow-tstart, "intmediate", totlengt(0)-target_length, Tfluxave(0)-target_tflux, Bdotnsqd(0), ff, ferr, tnow-told
!    
!    told = tnow
!    
!   case( 2 ) ! compute fvec;
!    
!    call dforce( Ndof, xdof(1:Ndof), ff, fdof(1:Ndof) )
!    
!   end select ! end select case( irevcm );
!   
!  enddo ! end of do;
!  
!  call dforce( Ndof, xdof(1:Ndof), ff, fdof(1:Ndof) )
!  
!  ferr = sqrt( sum(fdof(1:Ndof)*fdof(1:Ndof)) / Ndof )
!  
!  call archive( Ndof, xdof(1:Ndof), ferr )
!  
!  tnow = MPI_WTIME()
!  
!  if( myid.eq.0 ) write(ounit,1010) tnow-tstart, " finished ", totlengt(0)-target_length, Tfluxave(0)-target_tflux, Bdotnsqd(0), ff, ferr, tnow-told
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

!  call ode( xodes, Ndof, xdof(1:Ndof), taustart, ltauend, relerr, abserr, iflag, rwork(1:Lrwork), iwork(1:Liwork) )
!  call ode( qodes, Mdof, qdof(1:Mdof), taustart, ltauend, relerr, abserr, iflag, rwork(1:Lrwork), iwork(1:Liwork) )
   
  !tnow = MPI_WTIME()

!  xdof(1:Ndof) = qdof(1:Ndof)
   
  !call dforce( Ndof, xdof(1:Ndof), ff, fdof(1:Ndof) )
   
  !if( myid.eq. 0) then
    
  ! select case( iflag )
  ! case( 2 )
  !  write(ounit,1010) tnow-tstart, "difference", totlengt(0)-target_length, Tfluxave(0)-target_tflux, Bdotnsqd(0), ff, ferr, tnow-told  
  ! case( 3 )
  !  write(ounit,'("descent : ",f10.1," ; called  ode/odes ; iflag =",i3," ; tautol or abserr too small ;")') tnow-tstart, iflag
  ! case ( 4 )
  !  write(ounit,'("descent : ",f10.1," ; called  ode/odes ; iflag =",i3," ; tauend not reached ; steps > 500 ;")') tnow-tstart, iflag
  ! case ( 5 )
  !  write(ounit,'("descent : ",f10.1," ; called  ode/odes ; iflag =",i3," ; tauend not reached ; equation is stiff ;")') tnow-tstart, iflag
  ! case ( 6 )
  !  write(ounit,'("descent : ",f10.1," ; called  ode/odes ; iflag =",i3," ; invalid input ;")') tnow-tstart, iflag
  ! case default
  !  write(ounit,'("descent : ",f10.1," ; called  ode/odes ; iflag =",i3," ; unrecognized ;")') tnow-tstart, iflag
  ! end select
    
  !endif ! end of if( myid.eq.0 ) ;
   
!  if( sum(qdof(Ndof+1:Mdof)*fdof(1:Ndof)).gt.zero ) qdof(Ndof+1:Mdof) = -fdof(1:Ndof) ! going `uphill';

!  if( ff.gt.fold ) then
! ! qdof(Ndof+1:Mdof) = -fdof(1:Ndof) ! velocity ;
!   relerr = relerr * half
!   if( myid.eq.0 ) write(ounit,'("descent : ", 10x ," : decreasing xtol =",es9.2," ;")') relerr
!  endif

!  if( ff.lt.ffbest ) ffbest = ff
   
  !if( ff.gt.( 1 * fforig + 4 * ffbest ) / 5 ) exit ! exit iter loop;

!  fold = ff

!  enddo ! end of do iter;
  
1010 format("descent : ",f10.1," : ",a10," : L =",es18.10," ; F =",es18.10," ; ":"B =",es17.10," ; O =",es17.10," ; |dO| =",es12.05," ; time =",f9.2,"s ;")

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
! DALLOCATE( rwork )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  return
  
end subroutine descent

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

subroutine xodes( tau, xdof, fdof )
  
  use globals          , only :
  use descentparameters
  
  implicit none
  
  include "mpif.h"
  
  REAL                 :: tau, xdof(*), fdof(*)
  
  INTEGER              :: Ndof, ierr
  REAL                 :: ff
  
  Ndof = Ndegreeoffreedom
  
  call dforce( Ndof, xdof(1:Ndof), ff, fdof(1:Ndof) )
  
  fdof(1:Ndof) = - fdof(1:Ndof)

  return
  
end subroutine xodes

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

!subroutine qodes( tau, qdof, fdof )
!  
!  use globals          , only : friction
!  use descentparameters
!  
!  implicit none
!  
!  include "mpif.h"
!  
!  REAL                 :: tau, qdof(*), fdof(*)
!  
!  INTEGER              :: Ndof, ierr
!  REAL                 :: ff
!  
!  Ndof = Ndegreeoffreedom
!  
!  call dforce( Ndof, qdof(1:Ndof), ff, fdof(1:Ndof) )
!  
!  fdof(     1:  Ndof) =                             qdof(Ndof+1:2*Ndof)
!  fdof(Ndof+1:2*Ndof) = - fdof(1:Ndof) - friction * qdof(Ndof+1:2*Ndof)
!
!  return
!  
!end subroutine qodes

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
