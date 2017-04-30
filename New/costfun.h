!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine costfun(ideriv)
  use globals, only: zero, sqrtmachprec, myid, ounit, IsQuiet, iter, &
       Ncoils, deriv, Ndof, xdof, &
       totalenergy, t1E, t2E, &
       bnorm      , t1B, t2B, weight_bnorm,  &
       bharm      , t1H, t2H, weight_bharm,  &
       tflux      , t1F, t2F, weight_tflux, target_tflux, isign, &
       ttlen      , t1L, t2L, weight_ttlen, &
       specw      , t1S, t2S, weight_specw, &
       ccsep      , t1C, t2C, weight_ccsep

  implicit none
  include "mpif.h"

  INTEGER, INTENT(in) :: ideriv

  INTEGER      :: astat, ierr

  REAL         :: start, finish
  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  iter = iter + 1
  ! if(IsQuiet .eq. -1 .and. myid .eq. 0) write(ounit,'("Costfun :"12X": begin the "I6"th run.")') iter
  if(IsQuiet <= -1) then
     call bnormal(0)
     if (weight_bharm .gt. sqrtmachprec) call bmnharm(0)
!!$
!!$   if ( target_tflux .eq. 0.0 ) then
!!$    call torflux(0)
!!$    target_tflux = isign*sqrt(2.0*tflux)             
!!$    if(myid .eq. 0) write(ounit,'("Costfun :"11X" : Reset target toroidal flux to"ES23.15)') target_tflux
!!$   endif
!!$
!!$   call torflux(0)
!!$
!!$   if (lc .eq. 1) then
!!$    call tlength(0)
!!$   else
!!$    call tlengthExp(0)
!!$   endif
!!$   
!!$   call specwid(0)
!!$   call coilsep(0)
  endif
  
  totalenergy = zero
  
  !if ( ideriv .ge. 1 .and. .not. allocated(t1E) ) ALLOCATE(t1E(1:Ndof))
  !if ( ideriv .eq. 2 .and. .not. allocated(t2E) ) ALLOCATE(t2E(1:Ncoils, 0:Cdof, 1:Ncoils, 0:Cdof))

  if     ( ideriv .eq. 1 ) then
   t1E = zero
  elseif ( ideriv .eq. 2 ) then
   t1E = zero; t2E = zero
  endif

  !call unpacking(xdof)

  ! Bnormal surface intergration;
  if (weight_bnorm .gt. sqrtmachprec) then

   call bnormal(ideriv)
   totalenergy = totalenergy + weight_bnorm * bnorm
   if     ( ideriv .eq. 1 ) then
    t1E = t1E +  weight_bnorm * t1B
   elseif ( ideriv .eq. 2 ) then
    t1E = t1E +  weight_bnorm * t1B
    t2E = t2E +  weight_bnorm * t2B
   endif
  endif

  ! individual Bn harmonics;
  if (weight_bharm .gt. sqrtmachprec) then

   call bmnharm(ideriv)
   totalenergy = totalenergy + weight_bharm * bharm
   if     ( ideriv .eq. 1 ) then
    t1E = t1E +  weight_bharm * t1H
   elseif ( ideriv .eq. 2 ) then
    t1E = t1E +  weight_bharm * t1H
    t2E = t2E +  weight_bharm * t2H
   endif
  endif

!!$  ! if(myid .eq. 0) write(ounit,'("calling bnormal used",f10.5,"seconds.")') finish-start
!!$
!!$  if (weight_tflux .gt. sqrtmachprec) then
!!$
!!$   if ( target_tflux .eq. 0.0 ) then
!!$    call torflux(0)
!!$    target_tflux = isign*sqrt(2.0*tflux)             
!!$    if(myid .eq. 0) write(ounit,'("Costfun :"11X" : Reset target toroidal flux to"ES23.15)') target_tflux
!!$   endif
!!$
!!$   call torflux(ideriv)
!!$   totalenergy = totalenergy + weight_tflux * tflux / target_tflux**2 ! normalization
!!$   if     ( ideriv .eq. 1 ) then
!!$    t1E = t1E +  weight_tflux * t1F / target_tflux**2
!!$   elseif ( ideriv .eq. 2 ) then
!!$    t1E = t1E +  weight_tflux * t1F / target_tflux**2
!!$    t2E = t2E +  weight_tflux * t2F / target_tflux**2
!!$   endif
!!$
!!$  endif
!!$
!!$  ! if(myid .eq. 0) write(ounit,'("calling torflux used",f10.5,"seconds.")') start-finish
!!$
!!$  if (weight_ttlen .gt. sqrtmachprec) then
!!$
!!$   if( lc .eq. 1 ) then 
!!$    call tlength(ideriv)
!!$   elseif (lc .eq. 2) then
!!$    call tlengthExp(ideriv)
!!$   else
!!$    FATAL( dnergy, .true., Conflicts between lc and weight_ttlen )
!!$   !call MPI_ABORT( MPI_COMM_WORLD, 1, ierr )
!!$   !stop "COSTFUN: Conflicts between lc and weight_ttlen"
!!$   endif
!!$   
!!$   totalenergy = totalenergy + weight_ttlen * ttlen
!!$   if     ( ideriv .eq. 1 ) then
!!$    t1E = t1E +  weight_ttlen * t1L
!!$   elseif ( ideriv .eq. 2 ) then
!!$    t1E = t1E +  weight_ttlen * t1L
!!$    t2E = t2E +  weight_ttlen * t2L
!!$   endif
!!$
!!$  endif
!!$
!!$  ! if(myid .eq. 0) write(ounit,'("calling tlength used",f10.5,"seconds.")') finish-start
!!$
!!$  if (weight_eqarc .ge. sqrtmachprec) then
!!$
!!$  ! call equarcl(ideriv)
!!$  ! call specwid_df(ideriv)
!!$   if ( Loptimize .eq. 3) then
!!$    call langrange(ideriv)
!!$   else
!!$    !call specwid(ideriv)
!!$    call specwid_df(ideriv)
!!$   endif
!!$   
!!$   totalenergy = totalenergy + weight_eqarc * eqarc
!!$   if     ( ideriv .eq. 1 ) then
!!$    t1E = t1E + weight_eqarc * t1A
!!$   elseif ( ideriv .eq. 2 ) then
!!$    t1E = t1E + weight_eqarc * t1A
!!$    t2E = t2E + weight_eqarc * t2A
!!$   endif
!!$
!!$  endif
!!$
!!$  if (weight_ccsep .ge. sqrtmachprec) then
!!$
!!$   call coilsep(ideriv)
!!$   totalenergy = totalenergy + weight_ccsep * ccsep
!!$   if     ( ideriv .eq. 1 ) then
!!$    t1E = t1E + weight_ccsep * t1C
!!$   elseif ( ideriv .eq. 2 ) then
!!$    t1E = t1E + weight_ccsep * t1C
!!$    t2E = t2E + weight_ccsep * t2C
!!$   endif
!!$
!!$  endif
!!$
!!$  if(allocated(deriv)) then
!!$     deriv = zero
!!$     do ii = 1, Ndof
!!$        call DoFconvert(ii, icl, inf)
!!$        if(allocated(t1E)) deriv(ii, 0) = t1E(icl, inf)
!!$        if(allocated(t1B)) deriv(ii, 1) = t1B(icl, inf)
!!$        if(allocated(t1F)) deriv(ii, 2) = t1F(icl, inf)
!!$        if(allocated(t1L)) deriv(ii, 3) = t1L(icl, inf)
!!$        if(allocated(t1A)) deriv(ii, 4) = t1A(icl, inf)
!!$        if(allocated(t1C)) deriv(ii, 5) = t1C(icl, inf)
!!$     enddo
!!$  endif


  FATAL( costfun, iter > 100000, too many iterations for one single call)
 !if(iter .ge. 1E5) call MPI_ABORT( MPI_COMM_WORLD, 10, ierr )

  return
end subroutine costfun

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
