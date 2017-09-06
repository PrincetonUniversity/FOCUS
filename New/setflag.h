!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

!title (main) ! Main program.

!latex \briefly{Main program.}

!latex \calledby{\link{}}
!latex \calls{\link{initial}, \link{rdsurf}, \link{rdcoils}, 
!latex \link{saving}, \link{solvers}}.

!latex \tableofcontents

!latex \subsection{Brief introduction}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

subroutine setflag
  
  use globals, only : zero, one, pi2, sqrtmachprec, myid, ounit, &
                      Nseg, Ncoils, coil, cmt, smt, &
                      Ndof, xdof, &
                      IsNormalize, weight_gnorm, weight_inorm, Gnorm, Inorm, &
                      bsconstant, Bdotnsquared, &
                      weight_tflux , target_tflux , Ltflux , toroidalfluxaverage, toroidalfluxerror, &
                      weight_length, target_length, Llength, &
                      DF_maxiter, CG_maxiter, Ldescent, Lcgradient, &
                      iout, Nouts
  
  implicit none
  
  include "mpif.h"

  INTEGER   :: icoil, imn, ierr, ideriv, astat, NF, ii, mm
  REAL      :: tt
  CHARACTER :: packorunpack
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  NF = 0
  do icoil = 1, Ncoils
   NF = max( coil(icoil)%NF, NF )
  enddo
  write(ounit,'("setflag : " 10x " : NF =",i3," ;")') NF

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  SALLOCATE( cmt, (0:Nseg, 0:NF), zero )
  SALLOCATE( smt, (0:Nseg, 0:NF), zero )
  
  do ii = 0, Nseg ; tt = ii * pi2 / Nseg
   do mm = 0, NF
    cmt(ii,mm) = cos( mm * tt )
    smt(ii,mm) = sin( mm * tt )
   enddo
  enddo

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

  do icoil = 1, Ncoils

   call coilxyz( icoil )
   
  enddo
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  Ndof = 0 

  do icoil = 1, Ncoils
   
   if( coil(icoil)%Lc.ne.0 ) Ndof = Ndof + 6*coil(icoil)%NF+3 ! allow geometry to vary; 04 Sep 17;
   if( coil(icoil)%Ic.ne.0 ) Ndof = Ndof + 1                  ! allow current to vary; 04 Sep 17;
   
  enddo

  write(ounit,'("setflag : " 10x " : Ndof ="i9" ;")') Ndof

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

  SALLOCATE( xdof, (1:Ndof), zero )

  packorunpack = 'P' ; call packdof( packorunpack )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  select case( IsNormalize )
   
  case( 0 )
   
   Gnorm = one
   Inorm = one
   
  case( 1 )
   
   Gnorm = zero
   Inorm = zero ! initialize summation; SRH; 29 Sep 17;
   
   do icoil = 1, Ncoils
    
    do imn = 0, coil(icoil)%NF
     Gnorm = Gnorm + coil(icoil)%xc(imn)**2 + coil(icoil)%xs(imn)**2
     Gnorm = Gnorm + coil(icoil)%yc(imn)**2 + coil(icoil)%ys(imn)**2
     Gnorm = Gnorm + coil(icoil)%zc(imn)**2 + coil(icoil)%zs(imn)**2
    enddo
    
    Inorm = Inorm + coil(icoil)%I**2
    
   enddo
   
   Gnorm = sqrt(Gnorm) * weight_gnorm
   Inorm = sqrt(Inorm) * weight_inorm
   
   FATAL( rdcoils, weight_gnorm.lt.zero, invalid weight )
   FATAL( rdcoils, weight_inorm.lt.zero, invalid weight )
   
   if( myid == 0 ) write(ounit, '("setflag : " 10x " : currents normalized by ",es23.15," ; geometry normalized by ",es23.15," ;")') Inorm, Gnorm
   
  case default
   
   FATAL( rdcoils, .true., selected IsNormalize is not supported )
   
  end select

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

  packorunpack = 'U' ; call packdof( packorunpack )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  ideriv = 0 ; call bnormal( ideriv )
  
  write(ounit,'("setflag : " 10x " : bsconstant ="es12.5" ; Bdotnsquared ="es12.5" ;")') bsconstant, Bdotnsquared

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  if    ( weight_tflux.lt.-sqrtmachprec                                     ) then ! target toroidal flux is initial value   ; 04 Sep 17;
   
   Ltflux = .true.
   
   target_tflux = zero
   
   ideriv = 0 ; call torflux( ideriv ) ! first calculate toroidal flux; 04 Sep 17;
   
   target_tflux = toroidalfluxaverage ! set target toroidal flux;04 Sep 17;
   
   ideriv = 0 ; call torflux( ideriv ) ! then calculated toroidal flux error; 04 Sep 17;
   
   write(ounit,'("setflag : " 10x " : toroidalfluxaverage =",es13.5," ; toroidalfluxerror =",es12.5" ;")') toroidalfluxaverage, toroidalfluxerror
   
  elseif( weight_tflux.ge.-sqrtmachprec .and. weight_tflux.le.+sqrtmachprec ) then !        toroidal flux is not a constraint; 04 Sep 17;
   
   Ltflux = .false.
   
  elseif(                                     weight_tflux.gt.+sqrtmachprec ) then ! target toroidal flux is input           ; 04 Sep 17;
   
   Ltflux = .true.
   
  endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  if    ( weight_length.lt.-sqrtmachprec                                      ) then ! target coil length   is initial value   ; 04 Sep 17;
   
   Llength = .true.
   
   do icoil = 1, Ncoils
    ideriv = 0 ; call length( icoil, ideriv )
   enddo
   
   write(ounit,'("setflag : " 10x " : coil length =",99(es12.5,","))') ( coil(icoil)%L, icoil = 1, Ncoils )
   
  elseif( target_length.ge.-sqrtmachprec .and. target_length.le.+sqrtmachprec ) then !        coil length   is not a constraint; 04 Sep 17;
   
   Llength = .false.
   
  elseif(                                      target_length.gt.+sqrtmachprec ) then ! target coil length   is input           ; 04 Sep 17;
   
   Llength = .true.
   
  endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  if( DF_maxiter.gt.0 ) then
   
   Ldescent = .true.

   write(ounit,'("setflag : " 10x " : Ldescent ="L2" ; ")') Ldescent

  endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  if( CG_maxiter.gt.0 ) then

   Lcgradient = .true.

   write(ounit,'("setflag : " 10x " : Lcgradient ="L2" ; ")') Lcgradient

  endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

! iout = 1 ; Nouts = DF_maxiter + CG_maxiter !+ HN_maxiter + TN_maxiter

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  return
  
end subroutine setflag

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
