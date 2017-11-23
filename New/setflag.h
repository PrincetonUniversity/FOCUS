!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

!title (main) ! Main program.

!latex \briefly{Main program.}

!latex \calledby{\link{}}
!latex \calls{\link{initial}, \link{rdsurf}, \link{rdcoils}, 
!latex \link{saving}, \link{solvers}}.

!latex \tableofcontents

!latex \subsection{Brief introduction}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

subroutine setflag( Ndof, xdof )
  
  use globals, only : zero, one, two, half, pi2, sqrtmachprec, myid, ncpu, ounit, tstart, &
                      Icheck, Ns, Ncoils, coil, &
                      totlengt, Tfluxave, Bdotnsqd, Tfluxerr, &
                      weight_tflux , target_tflux, &
                      NRKsave, NRKstep, RKstep, Ldescent
  
  implicit none
  
  include "mpif.h"
  
  INTEGER :: Ndof
  REAL    ::  xdof(1:Ndof)
  REAL    :: oxdof(1:Ndof), off(-1:1), ofdof(1:Ndof,-1:1)
  
  INTEGER :: icoil, mm, idof, ifd, ierr, isurf
  REAL    :: fd, est, tnow

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  FATAL( setflag,         Ndof.le.   0, illegal )
  FATAL( setflag, sqrtmachprec.le.zero, illegal )
  FATAL( setflag,         half.le.zero, illegal )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  if( Icheck.eq.1 ) then
   
   oxdof(1:Ndof) = xdof(1:Ndof) ! original;
   
   fd = sqrt(sqrtmachprec)
   
   do idof = 1, Ndof
    
    do ifd = -1, 1
     
     xdof(1:Ndof) = oxdof(1:Ndof) ; xdof(idof) = xdof(idof) + ifd * fd * half
     
     isurf = 1 ; call dforce( isurf, Ndof, xdof(1:Ndof), off(ifd), ofdof(1:Ndof,ifd) )
     
    enddo ! end of do ifd;
    
    CHECK( setflag, abs(fd).lt. sqrtmachprec, divide by zero )
    
    est = ( off(+1) - off(-1) ) / fd
    
    if( myid.eq.0 ) then
     if( abs(ofdof(idof,0)).gt.sqrtmachprec ) write(ounit,1000) idof, ofdof(idof,0), est, ( ofdof(idof,0) - est ) / ofdof(idof,0)
     if( abs(ofdof(idof,0)).le.sqrtmachprec ) write(ounit,1000) idof, ofdof(idof,0), est
    endif
    
   enddo ! end of do idof;
   
1000 format("setflag : ", 10x ," : check derivatives : ",i3," ; analytic =",es23.15," =",es23.15," = finite diff. ; ",:," relative err. =",es10.02," ;")
   

   ifd = 0
   
   xdof(1:Ndof) = oxdof(1:Ndof)
   
   isurf = 1 ; call dforce( isurf, Ndof, xdof(1:Ndof), off(ifd), ofdof(1:Ndof,ifd) ) ! just to recover initial values of force;


  endif ! end of if( Icheck.eq.1 ) ;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

!  select case( IsNormalize )
!   
!  case( 0 )
!   
!   Gnorm = one
!   Inorm = one
!   
!  case( 1 )
!   
!   Gnorm = zero ! initialize summation;
!   Inorm = zero
!   
!   do icoil = 1, Ncoils
!    
!    do mm = 0, coil(icoil)%NF
!     Gnorm = Gnorm + coil(icoil)%xc(mm)**2 + coil(icoil)%xs(mm)**2
!     Gnorm = Gnorm + coil(icoil)%yc(mm)**2 + coil(icoil)%ys(mm)**2
!     Gnorm = Gnorm + coil(icoil)%zc(mm)**2 + coil(icoil)%zs(mm)**2
!    enddo
!    
!    Inorm = Inorm + coil(icoil)%I**2
!    
!   enddo
!   
!   Gnorm = sqrt(Gnorm) * weight_gnorm
!   Inorm = sqrt(Inorm) * weight_inorm
!   
!   FATAL( rdcoils, weight_gnorm.lt.zero, invalid weight )
!   FATAL( rdcoils, weight_inorm.lt.zero, invalid weight )
!   
!   if( myid.eq.0 ) write(ounit, '("setflag : " 10x " : currents normalized by ",es23.15," ; geometry normalized by ",es23.15," ;")') Inorm, Gnorm
!   
!  case default
!   
!   FATAL( rdcoils, .true., selected IsNormalize is not supported )
!   
!  end select
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  tnow = MPI_WTIME()
  
  if( NRKsave.gt.0 .and. NRKstep.gt.0 ) then
   
   Ldescent = .true.
   
   if( myid.eq.0 ) write(ounit,'("setflag : "f10.1" : {{ NRKsave > 0 & NRKstep > 0 }} => {{ Ldescent ="L2" }} ; ")') tnow-tstart, Ldescent
   
   if( RKstep.lt.zero ) RKstep = 1.0E-02
   
  else
   
   Ldescent = .false.
   
   if( myid.eq.0 ) write(ounit,'("setflag : "f10.1" : {{ NRKsave < 0 || NRKstep < 0 }} => {{ Ldescent ="L2" }} ; ")') tnow-tstart, Ldescent
   
  endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  return
  
end subroutine setflag

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
