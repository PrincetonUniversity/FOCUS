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
                      weight_ttlen , target_length, weight_tflux , target_tflux, &
                      tauend, Ntauout, Ldescent
  
  implicit none
  
  include "mpif.h"
  
  INTEGER :: Ndof
  REAL    ::  xdof(1:Ndof)
  REAL    :: oxdof(1:Ndof), off(-1:1), ofdof(1:Ndof,-1:1)
  
  INTEGER :: icoil, mm, idof, ifd, ierr
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
     
     call dforce( Ndof, xdof(1:Ndof), off(ifd), ofdof(1:Ndof,ifd) )
     
    enddo ! end of do ifd;
    
    CHECK( setflag, abs(fd).lt. sqrtmachprec, divide by zero )
    
    est = ( off(+1) - off(-1) ) / fd
    
    if( myid.eq.0 ) then
     if( abs(ofdof(idof,0)).gt.sqrtmachprec ) write(ounit,1000) idof, ofdof(idof,0), est, ( ofdof(idof,0) - est ) / ofdof(idof,0)
     if( abs(ofdof(idof,0)).le.sqrtmachprec ) write(ounit,1000) idof, ofdof(idof,0), est
    endif
    
   enddo ! end of do idof;
   
1000 format("setflag : ", 10x ," : ",i3," ; analytic =",es23.15," =",es23.15," = finite diff. ; ",:,"err =",es10.02," ;")
   
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
  
!  if    ( weight_ttlen .lt.-sqrtmachprec                                      ) then ! target coil length   is initial value   ;
!   
!   Llength = .true. ; target_length = totlengt(0)
!   
!   write(ounit,'("setflag : " 10x " : L_t =",es12.5," ; ":"L_i ="99(es12.5,","))') totlengt(0)!, ( coil(icoil)%dL(0), icoil = 1, Ncoils )
!   
!  elseif( target_length.ge.-sqrtmachprec .and. target_length.le.+sqrtmachprec ) then !        coil length   is not a constraint;
!   
!   Llength = .false.
!   
!  elseif(                                      target_length.gt.+sqrtmachprec ) then ! target coil length   is input           ;
!   
!   Llength = .true.
!   
!  endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

!  if    ( weight_tflux.lt.-sqrtmachprec                                     ) then ! target toroidal flux is initial value   ;
!   
!   Ltflux = .true. ; target_tflux = Tfluxave(0)
!   
!   Tfluxerr = sum( (surf%dT(0:Nz-1,0)-target_tflux)**2 ) * half
!   
!   write(ounit,'("setflag : " 10x " : Tfluxave =",es23.15" ;":" Tfluxerr ="es12.5" ;")') Tfluxave(0), Tfluxerr
!   
!  elseif( weight_tflux.ge.-sqrtmachprec .and. weight_tflux.le.+sqrtmachprec ) then !        toroidal flux is not a constraint;
!   
!   Ltflux = .false. ; weight_tflux = zero
!   
!  elseif(                                     weight_tflux.gt.+sqrtmachprec ) then ! target toroidal flux is input           ;
!   
!   Ltflux = .true.
!   
!  endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

!  if    ( weight_bnorm.lt.-sqrtmachprec                                     ) then ! target toroidal flux is initial value   ;
!   
!   Lbnorm = .true. ; target_bnorm = zero
!   
!  !Tfluxerr = sum( (surf%dT(0:Nz-1,0)-target_tflux)**2 ) * half
!   
!   write(ounit,'("setflag : " 10x " : Bdotnsqd =",es12.5" ;")') Bdotnsqd
!   
!  elseif( weight_bnorm.ge.-sqrtmachprec .and. weight_bnorm.le.+sqrtmachprec ) then !        toroidal flux is not a constraint;
!   
!   FATAL( setflag, .true., nonsense weight_bnorm )
!   
!   Lbnorm = .false. ; weight_bnorm = zero
!   
!  elseif(                                     weight_bnorm.gt.+sqrtmachprec ) then ! target toroidal flux is input           ;
!   
!   Ltflux = .true.
!   
!  endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

  if( tauend.gt.zero ) then
   
   Ldescent = .true.
   
   tnow = MPI_WTIME()
   
   if( myid.eq.0 ) write(ounit,'("setflag : "f10.1" : {{ tauend > 0.0 }} => {{ Ldescent ="L2" }} ; ")') tnow-tstart, Ldescent
   
   if( Ntauout.le.0 ) Ntauout = 1
   
  endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
!  if( CG_maxiter.gt.0 ) then
!
!   Lcgradient = .true.
!
!   if( myid.eq.0 ) write(ounit,'("setflag : " 10x " : Lcgradient ="L2" ; ")') Lcgradient
!
!  endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

  return
  
end subroutine setflag

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
