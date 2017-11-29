!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

!title (descent) ! Minimize the target function via a differential flow.

!latex \briefly{The minimization problem is solved by integrating a system of ODEs.}

!latex \calledby{\link{solvers}}
!latex \calls{\link{solvers}}

!latex \section{overview}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

subroutine descent( Ndof, xdof, ferr )
  
  use globals          , only : zero, one, six, half, pi2, sqrtmachprec, small, myid, ncpu, ounit, tstart, &
                                Ncoils, &
                                NRKsave, NRKstep, RKstep, &
                                totlengt, Tfluxave, Bdotnsqd, &
                                target_tflux, converged, &
                                fforig, ffbest

  implicit none
  
  include "mpif.h"
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  INTEGER              :: Ndof
  REAL                 :: xdof(1:Ndof), ferr

  INTEGER              :: irksave, irkstep, isurf
  REAL                 :: ffold, ydof(1:Ndof), fk(1:Ndof,1:4), tnow, told, ff, fdof(1:Ndof), fold

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
    
  tnow = MPI_WTIME()
  
  if( myid.eq.0 ) write(ounit,1000) tnow-tstart, NRKsave, NRKstep, RKstep, converged
1000 format("descent : ",f10.1," : Runge-Kutta ; NRKsave (outer loop) =",i9," ; NRKstep (inner loop) =",i5," ; RKstep =",es12.5," ; converged =",es13.5," ;")
  
  told = tnow
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  isurf = 1
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

  call dforce( isurf, Ndof, xdof(1:Ndof), ff, fdof(1:Ndof) ) ; fk(1:Ndof,1) = - fdof(1:Ndof)
  
  ferr = sqrt( sum(fdof(1:Ndof)*fdof(1:Ndof)) / Ndof ) / ff
  
  if( ferr.lt.converged ) return
  
  ffold = ff ! initialize;
  
  do irksave = 1, NRKsave ! allows intermediate output & archive;
   
   do irkstep = 1, NRKstep ! allows intermediate output & archive;
    
    ydof(1:Ndof) = xdof(1:Ndof) + fk(1:Ndof,1) * RKstep * half ; call dforce( isurf, Ndof, ydof(1:Ndof), ff, fdof(1:Ndof) ) ; fk(1:Ndof,2) = - fdof(1:Ndof)
    ydof(1:Ndof) = xdof(1:Ndof) + fk(1:Ndof,2) * RKstep * half ; call dforce( isurf, Ndof, ydof(1:Ndof), ff, fdof(1:Ndof) ) ; fk(1:Ndof,3) = - fdof(1:Ndof)
    ydof(1:Ndof) = xdof(1:Ndof) + fk(1:Ndof,3) * RKstep        ; call dforce( isurf, Ndof, ydof(1:Ndof), ff, fdof(1:Ndof) ) ; fk(1:Ndof,4) = - fdof(1:Ndof)
    
    xdof(1:Ndof) = xdof(1:Ndof) + ( fk(1:Ndof,1) + 2 * fk(1:Ndof,2) + 2 * fk(1:Ndof,3) + fk(1:Ndof,4) ) * RKstep / six ! Runge-Kutta step;
    
    call dforce( isurf, Ndof, xdof(1:Ndof), ff, fdof(1:Ndof) ) ; fk(1:Ndof,1) = - fdof(1:Ndof)

    ferr = sqrt( sum(fdof(1:Ndof)*fdof(1:Ndof)) / Ndof ) / ff
    
    if( ff.gt.ffold ) then
     if( myid.eq.0 ) write(ounit,1010) tnow-tstart, "INCREASING", totlengt(0)/pi2, Tfluxave(0)-target_tflux, Bdotnsqd(0), ff, ferr, tnow-told
     RKstep = RKstep * half
    !return
    endif
    
    ffold = ff

   enddo ! end of do irkstep ;
      
   call archive( Ndof, xdof(1:Ndof), ferr )
   
   tnow = MPI_WTIME()
   
   if( myid.eq.0 ) write(ounit,1010) tnow-tstart, "difference", totlengt(0)/pi2, Tfluxave(0)-target_tflux, Bdotnsqd(0), ff, ferr, tnow-told
   
   told = tnow

   if( ferr.lt.converged ) exit
    
  enddo ! end of do irksave ;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
1010 format("descent : ",f10.1," : ",a10," : L =",es18.10," ; F =",es18.10," ; ":"B =",es17.10," ; O =",es17.10," ; |dO| =",es12.05," ; time =",f9.2,"s ;")

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  return
  
end subroutine descent

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
