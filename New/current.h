!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! not parallelized; communications may take more time;
subroutine current(ideriv)
  use globals, only : zero, half, two, pi2, machprec, ncpu, myid, ounit, &
       coil, DoF, Ncoils, Nfixcur, Ndof, icost, t1I, target_current

  implicit none
  include "mpif.h"
  INTEGER, INTENT(in) :: ideriv

  INTEGER             :: astat, ierr, icoil, idof, ND

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  icost = zero

  !FATAL( current , Nfixcur .eq. Ncoils, No free currents )

  if( ideriv >= 0 ) then

        do icoil = 1, Ncoils
           if ( coil(icoil)%Ic /= 0 ) icost = icost + exp( (coil(icoil)%I / target_current)**2 )
        enddo

     icost = icost / (Ncoils - Nfixcur + machprec)

  endif

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  FATAL( current, .not. allocated(t1I), Please allocate t1I first )

  if ( ideriv >= 1 ) then

     t1I = zero; idof = 0

     do icoil = 1, Ncoils

        ND = DoF(icoil)%ND

        if ( coil(icoil)%Ic /= 0 ) then !if current is free;
           idof = idof +1
           t1I(idof) = coil(icoil)%I * exp( (coil(icoil)%I / target_current)**2 )
        endif

        if ( coil(icoil)%Lc /= 0 ) then !if geometry is free;
           idof = idof + ND
        endif

     enddo !end icoil;

     t1I = two * t1I / (Ncoils - Nfixcur + machprec) / target_current**2

  endif

  return
end subroutine current

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
