subroutine hessian
  use kmodule, only : zero, half, ounit, myid, Ndof, totalenergy, t1E, t2E
  implicit none
  include "mpif.h"

  INTEGER :: N, ii, jj, c1, c2, n1, n2, idof, ierr, astat, nstep
  REAL    :: gradnorm, df, php, f0
  REAL, allocatable :: xdof(:), bdof(:), px(:), hess(:,:)
  REAL, parameter :: dx = 0.01
  REAL, parameter :: small = 1.0D-06
  REAL    :: pxhxp

  N = Ndof
  FATAL( hessian, N <= 0, invalid # of dof )

  SALLOCATE(xdof, (1:N), zero)
  call pack(xdof)

  ! test if gradient small enough;
  call costfun(1)
  gradnorm = sqrt(sum(t1E**2))
  if (gradnorm > small) then
     write(ounit, '("hessian : "10X" : |gradient| = "ES23.15" ; not small enough!")') gradnorm
     return
  endif

  ! construct hessian matrix;
  SALLOCATE( hess, (n,n), zero )
  call costfun(2)
  do jj = 1, N
     call DoFconvert(jj,c2,n2)
     do ii = 1,N
        call DoFconvert(ii,c1,n1)
        hess(ii,jj) = t2E(c1,n1,c2,n2)
     enddo
  enddo

  ! perturb the dof vector;
  SALLOCATE(bdof, (1:N), zero)
  SALLOCATE(px  , (1:N), zero)

  bdof = xdof ! backup vector;
  f0 = totalenergy ! original value
  if (myid == 0) write(ounit, '("hessian :idof / Ndof : ", 3(A23" ; "))') 'actual value', 'hessian approxi', 'abs diff'

  do idof = 1, N
     px = zero; px(idof) = dx 
     xdof = bdof + px
     call unpack(xdof) ; call costfun(0)
     df = totalenergy - f0
     php = half * pxhxp(px, hess, N)
     if (myid == 0) write(ounit, '("hessian :"I4" / "I4" : ", 3(ES23.15" ; "))') idof, N, df, php, abs(df-php)
  enddo

  xdof = bdof 

  DALLOCATE( xdof )
  DALLOCATE( bdof )
  DALLOCATE(  px  )
  DALLOCATE( hess )

  return

end subroutine hessian

!--------------------------------------------------

real function pxhxp(x, h, n)
  implicit none
  include "mpif.h"

  INTEGER, INTENT(IN) :: n
  REAL   , INTENT(IN) :: x(1:n, 1:1), h(1:n, 1:n)

  REAL                :: xp(1:1, 1:n), result(1:1, 1:1)

  xp(1,1:n) = x(1:n,1)

  result = matmul( matmul(xp, h), x )
  pxhxp = result(1, 1)

end function pxhxp
