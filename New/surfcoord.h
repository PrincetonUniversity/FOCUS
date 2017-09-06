subroutine surfcoord( theta, zeta, r, z)
  use globals, only: zero, Nfou, bim, bin, Rbc, Rbs, Zbc, Zbs
  implicit none
  include "mpif.h"

  REAL, INTENT(in ) :: theta, zeta
  REAL, INTENT(out) :: r, z

  INTEGER           :: imn
  REAL              :: arg
  !-------------calculate r, z coodinates for theta, zeta------------------------------------------------  
  if( .not. allocated(bim) ) STOP  "please allocate surface data first!"

  r = zero; z = zero  
  do imn = 1, Nfou
     arg = bim(imn) * theta - bin(imn) * zeta
     R =  R + Rbc(imn) * cos(arg) + Rbs(imn) * sin(arg)
     Z =  Z + Zbc(imn) * cos(arg) + Zbs(imn) * sin(arg)
  enddo

  return
end subroutine surfcoord

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
