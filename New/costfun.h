!title (costfun) ! cost funtion

!latex \briefly{cost function}

!latex \calledby{\link{focus}}
!latex \calls{\link{bnormal}, \link{bmnharm}, \link{torflux}, \link{length}, \link{coilsep}, 
!latex         \link{descent}, \link{congrad}, \link{hybridnt}, \link{truncnt}}

subroutine costfun(ideriv)
  
  use globals, only : zero, one, sqrtmachprec, myid, ounit, astat, ierr, IsQuiet, &
                      Ncoils, deriv, Ndof, xdof, dofnorm, coil, &
                      chi, t1E, t2E, &
                      Bdotnsquared      , t1B, t2B, weight_bnorm,  &
                      bharm      , t1H, t2H, weight_bharm,  &
                      toroidalfluxerror      , t1F, t2F, weight_tflux, target_tflux, toroidalfluxaverage, &
                      ttlen      , t1L, t2L, weight_ttlen, case_length, &
                      specw      , t1S, t2S, weight_specw, &
                      ccsep      , t1C, t2C, weight_ccsep
  
  implicit none
  
  include "mpif.h"
  
  INTEGER, INTENT(in) :: ideriv

  REAL                :: start, finish
  
  call mpi_barrier(MPI_COMM_WORLD, ierr)
  
  if (IsQuiet <= -2) then
   
     call bnormal(0)

     if( weight_bharm > sqrtmachprec ) call bmnharm(0)

     if( abs(target_tflux) < sqrtmachprec ) then
        call torflux(0)
        target_tflux = toroidalfluxaverage
        if (myid == 0) write(ounit,'("costfun : reset target toroidal flux to ",es23.15," ;")') target_tflux
     endif

     call torflux(0)

     if ( (case_length == 1) .and. (sum(coil(1:Ncoils)%Lo) < sqrtmachprec) ) then
        coil(1:Ncoils)%Lo = one
        call length(0)
        coil(1:Ncoils)%Lo = coil(1:Ncoils)%L
        if(myid .eq. 0) write(ounit,'("solvers : reset target coil length to the current actual length. ")')
     endif

     call length(0)

  endif

  chi = zero

  if ( ideriv == 1 ) then
     t1E = zero
  elseif ( ideriv == 2 ) then
     t1E = zero; t2E = zero
  endif

  !call unpacking(xdof)
  ! Bnormal surface intergration;
  if (weight_bnorm > sqrtmachprec) then
 
     call bnormal(ideriv)
     chi = chi + weight_bnorm * Bdotnsquared
     if     ( ideriv == 1 ) then
        t1E = t1E +  weight_bnorm * t1B
     elseif ( ideriv == 2 ) then
        t1E = t1E +  weight_bnorm * t1B
        t2E = t2E +  weight_bnorm * t2B
     endif
  endif

  ! individual Bn harmonics;
  if (weight_bharm > sqrtmachprec) then

     call bmnharm(ideriv)
     chi = chi + weight_bharm * bharm
     if     ( ideriv == 1 ) then
        t1E = t1E +  weight_bharm * t1H
     elseif ( ideriv == 2 ) then
        t1E = t1E +  weight_bharm * t1H
        t2E = t2E +  weight_bharm * t2H
     endif
  endif

  ! toroidal flux;
  if (weight_tflux > sqrtmachprec) then
     if ( abs(target_tflux) < sqrtmachprec ) then
        call torflux(0)
        target_tflux = toroidalfluxaverage        
        if (myid==0) write(ounit,'("solvers : Reset target toroidal flux to "ES23.15)') target_tflux
     endif

     call torflux(ideriv)
     chi = chi + weight_tflux * toroidalfluxerror / target_tflux**2 ! normalization
     if     ( ideriv == 1 ) then
        t1E = t1E +  weight_tflux * t1F / target_tflux**2
     elseif ( ideriv == 2 ) then
        t1E = t1E +  weight_tflux * t1F / target_tflux**2
        t2E = t2E +  weight_tflux * t2F / target_tflux**2
     endif
  endif

  ! coil length;
  if (weight_ttlen > sqrtmachprec) then

     if ( (case_length == 1) .and. (sum(coil(1:Ncoils)%Lo) < sqrtmachprec) ) then
        coil(1:Ncoils)%Lo = one
        call length(0)
        coil(1:Ncoils)%Lo = coil(1:Ncoils)%L
        if(myid .eq. 0) write(ounit,'("solvers : reset target coil length to the current actual length. ")')
     endif

     call length(ideriv)

     chi = chi + weight_ttlen * ttlen
     if     ( ideriv == 1 ) then
        t1E = t1E +  weight_ttlen * t1L
     elseif ( ideriv == 2 ) then
        t1E = t1E +  weight_ttlen * t1L
        t2E = t2E +  weight_ttlen * t2L
     endif

  endif
!!$
!!$  ! if (myid == 0) write(ounit,'("calling tlength used",f10.5,"seconds.")') finish-start
!!$
!!$  if (weight_eqarc .ge. sqrtmachprec) then
!!$
!!$  ! call equarcl(ideriv)
!!$  ! call specwid_df(ideriv)
!!$   if ( Loptimize == 3) then
!!$    call langrange(ideriv)
!!$   else
!!$    !call specwid(ideriv)
!!$    call specwid_df(ideriv)
!!$   endif
!!$   
!!$   chi = chi + weight_eqarc * eqarc
!!$   if     ( ideriv == 1 ) then
!!$    t1E = t1E + weight_eqarc * t1A
!!$   elseif ( ideriv == 2 ) then
!!$    t1E = t1E + weight_eqarc * t1A
!!$    t2E = t2E + weight_eqarc * t2A
!!$   endif
!!$
!!$  endif
!!$
!!$  if (weight_ccsep .ge. sqrtmachprec) then
!!$
!!$   call coilsep(ideriv)
!!$   chi = chi + weight_ccsep * ccsep
!!$   if     ( ideriv == 1 ) then
!!$    t1E = t1E + weight_ccsep * t1C
!!$   elseif ( ideriv == 2 ) then
!!$    t1E = t1E + weight_ccsep * t1C
!!$    t2E = t2E + weight_ccsep * t2C
!!$   endif
!!$
!!$  endif
!!$
!!$  if (allocated(deriv)) then
!!$     deriv = zero
!!$     do ii = 1, Ndof
!!$        call DoFconvert(ii, icl, inf)
!!$        if (allocated(t1E)) deriv(ii, 0) = t1E(icl, inf)
!!$        if (allocated(t1B)) deriv(ii, 1) = t1B(icl, inf)
!!$        if (allocated(t1F)) deriv(ii, 2) = t1F(icl, inf)
!!$        if (allocated(t1L)) deriv(ii, 3) = t1L(icl, inf)
!!$        if (allocated(t1A)) deriv(ii, 4) = t1A(icl, inf)
!!$        if (allocated(t1C)) deriv(ii, 5) = t1C(icl, inf)
!!$     enddo
!!$  endif

  if ( ideriv == 1 ) then                ! multiply t1E & t2E with normalized terms; 06/09/2017
     t1E = t1E * dofnorm
!!$   elseif ( ideriv == 2 ) then
!!$      t1E = t1E * dofnorm
!!$      t2E = t2E * hesnorm
  endif

  call mpi_barrier(MPI_COMM_WORLD, ierr)

  return
end subroutine costfun

