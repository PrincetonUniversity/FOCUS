
!title (coil geometry) ! Specification of coil geometry.

!latex \briefly{The coil geometry, in real space, is determined.}

!latex \calledby{\oculus{bs00aa}}

!latex \tableofcontents

!latex \subsection{overview}
!latex \bi
!latex \item[2.] The coils are assumed to be closed (i.e. periodic), one-dimensional loops embedded in three-dimensional space, 
!latex           with position described by ${\bf x}(t) \equiv x(t) {\bf i} + y(t) {\bf j} + z(t) {\bf k}$, 
!latex           with the arbitrary curve parameter $t\in[0,2\pi]$, and ${\bf x}(t+2\pi) = {\bf x}(t)$,
!latex \item[2.] Presently, a Fourier representation is assumed, e.g.
!latex           \be x & = & \sum_{n=0}^{N} x_{n,c} \cos(m t) + \sum_{n=1}^{N} x_{n,s} \sin(m t),
!l tex               y & = & \sum_{n=0}^{N} y_{n,c} \cos(m t) + \sum_{n=1}^{N} y_{n,s} \sin(m t), \\
!l tex               z & = & \sum_{n=0}^{N} z_{n,c} \cos(m t) + \sum_{n=1}^{N} z_{n,s} \sin(m t).
!latex           \ee
!latex \ei

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine iccoil( tt, xx, yy, zz, ifail )
  
  use kmodule, only : zero, myid, coil, icoil
  
  implicit none
  
  include "mpif.h"
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  REAL                 :: tt, xx(0:1), yy(0:1), zz(0:1)
  INTEGER              :: ifail
  
  INTEGER              :: ierr, mm
  REAL                 :: carg, sarg
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  xx(0:1) = zero ; yy(0:1) = zero ; zz(0:1) = zero

  do mm = 0, coil(icoil)%N ; carg = cos(mm*tt) ; sarg = sin(mm*tt)
   xx(0) = xx(0) +     coil(icoil)%xc(mm) * carg + coil(icoil)%xs(mm) * sarg
   xx(1) = xx(1) + ( - coil(icoil)%xc(mm) * sarg + coil(icoil)%xs(mm) * carg ) * mm
   yy(0) = yy(0) +     coil(icoil)%yc(mm) * carg + coil(icoil)%ys(mm) * sarg
   yy(1) = yy(1) + ( - coil(icoil)%yc(mm) * sarg + coil(icoil)%ys(mm) * carg ) * mm
   zz(0) = zz(0) +     coil(icoil)%zc(mm) * carg + coil(icoil)%zs(mm) * sarg
   zz(1) = zz(1) + ( - coil(icoil)%zc(mm) * sarg + coil(icoil)%zs(mm) * carg ) * mm
  enddo
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
    
  return
  
end subroutine iccoil

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
