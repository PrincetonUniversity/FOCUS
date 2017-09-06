
subroutine normweight
  use globals, only : zero, one, machprec, ounit, myid, xdof, Bdotnsquared, bharm, toroidalfluxerror, ttlen, specw, ccsep, &
       weight_bnorm, weight_bharm, weight_tflux, weight_ttlen, weight_specw, weight_ccsep, target_tflux, &
       toroidalfluxaverage, coil, Ncoils, case_length

  implicit none  
  include "mpif.h"

  INTEGER    :: ierr, icoil
  REAL       :: tmp, cur_tflux

  !----------------------------------------------------------------------------------------------------

  call unpacking(xdof)

  !-!-!-!-!-!-!-!-!-!-tflux-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( weight_tflux .ge. machprec ) then

     call torflux(0)
     if ( abs(target_tflux) < machprec ) then
        target_tflux = toroidalfluxaverage
        if(myid .eq. 0) write(ounit,'("solvers : Reset target toroidal flux to "ES12.5)') target_tflux
     else if (sum(abs(coil(1:Ncoils)%Ic)) == Ncoils) then !only valid when all currents are free;
           do icoil = 1, Ncoils
              coil(icoil)%I = coil(icoil)%I * target_tflux / toroidalfluxaverage
           enddo
           if(myid .eq. 0) write(ounit,'("solvers : rescale coil currents with a factor of "ES12.5)') &
                target_tflux / toroidalfluxaverage
     endif

     call torflux(0)
     if (abs(toroidalfluxerror) > machprec) weight_tflux = weight_tflux / toroidalfluxerror * target_tflux**2
     if( myid .eq. 0 ) write(ounit, 1000) "weight_tflux", weight_tflux

  endif

  !-!-!-!-!-!-!-!-!-!-bnorm-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( weight_bnorm >= machprec ) then

     call bnormal(0)   
     if (abs(Bdotnsquared) > machprec) weight_bnorm = weight_bnorm / Bdotnsquared
     if( myid == 0 ) write(ounit, 1000) "weight_bnorm", weight_bnorm

  endif

  !-!-!-!-!-!-!-!-!-!-bnorm-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( weight_bharm >= machprec ) then

     call bmnharm(0)   
     if (abs(bharm) > machprec) weight_bharm = weight_bharm / bharm
     if( myid == 0 ) write(ounit, 1000) "weight_bharm", weight_bharm

  endif

  !-!-!-!-!-!-!-!-!-!-ttlen-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( weight_ttlen .ge. machprec ) then

     if ( (case_length == 1) .and. (sum(coil(1:Ncoils)%Lo) < machprec) ) then
        coil(1:Ncoils)%Lo = one
        call length(0)
        coil(1:Ncoils)%Lo = coil(1:Ncoils)%L
        if(myid .eq. 0) write(ounit,'("solvers : reset target coil length to the current actual length. ")')
     endif

     call length(0)

     if (abs(ttlen) .gt. machprec) weight_ttlen = weight_ttlen / ttlen
     if( myid .eq. 0 ) write(ounit, 1000) "weight_ttlen", weight_ttlen

  endif

!!$!-!-!-!-!-!-!-!-!-!-eqarc-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!!$
!!$  if( weight_eqarc .ge. machprec ) then
!!$
!!$   call specwid(0)
!!$   if (abs(eqarc) .gt. machprec) weight_eqarc = weight_eqarc / eqarc
!!$   if( myid .eq. 0 ) write(ounit, 1000) "weight_eqarc", weight_eqarc
!!$   
!!$  endif 
!!$
!!$!-!-!-!-!-!-!-!-!-!-ccsep-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!!$
!!$  if( weight_ccsep .ge. machprec ) then
!!$
!!$   call coilsep(0)
!!$   if (abs(ccsep) .gt. machprec) weight_ccsep = weight_ccsep / ccsep
!!$   if( myid .eq. 0 ) write(ounit, 1000) "weight_ccsep", weight_ccsep
!!$   
!!$  endif

1000 format("solvers : "A12" is normalized to" ES23.15)

  call packdof(xdof)

  return

end subroutine normweight

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

