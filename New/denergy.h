subroutine denergy( tau, lxdof, dE )
  
  use globals, only :  Ndof, myid, ounit, t1E
  
  implicit none

  include "mpif.h"
  
  REAL                 :: tau, lxdof(*), dE(*)  
  INTEGER              :: iorder
  
  external             :: unpacking, costfun
  
  call unpacking(lxdof(1:Ndof))
  
  iorder = 1 ; call costfun(iorder)
 
  dE(1:Ndof) = - t1E(1:Ndof)

  return

end subroutine denergy

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

