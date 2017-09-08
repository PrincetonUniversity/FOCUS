!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

!title (field) ! Computes magnetic vector potential and field.

!latex \briefly{Computes magnetic field given coil geometry.}

!latex \calledby{\link{bnormal}}
!latex \calls{}

!latex \tableofcontents

!latex \subsection{magnetic field}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

subroutine abfield( icoil, iteta, jzeta, NS, dAdx, dBdx )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  use globals, only : zero, one, two, three, sqrtmachprec, myid, ounit, &
                      surf, Nteta, Nzeta, &
                      Ncoils, coil, Nseg, deltacurveparameter, &
                      dBdx
  
  implicit none
  
  include "mpif.h"
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  INTEGER, intent(in ) :: icoil, iteta, jzeta, NS
  REAL   , intent(out) :: dAdx(0:NS,1:3,0:3), dBdx(0:NS,1:3,0:3)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  INTEGER :: ierr, astat, kk
  REAL    :: dx, dy, dz, rr(1:5), dd(1:5), tx, ty, tz, rdotdl, yz, zx, xy, xt, yt, zt
  REAL    :: dxtx, dxty, dxtz, dytx, dyty, dytz, dztx, dzty, dztz
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  CHECK( abfield , icoil.lt.1 .or. icoil.gt.Ncoils, icoil not in range )
  CHECK( abfield , iteta.lt.0 .or. iteta.ge.Nteta , iteta not in range )
  CHECK( abfield , jzeta.lt.0 .or. jzeta.ge.Nzeta , jzeta not in range )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  do kk = 1, Nseg
   
   dx = surf%xx(iteta,jzeta) - coil(icoil)%xx(kk) ! distance from evaluation point to curve; 04 Sep 17;
   dy = surf%yy(iteta,jzeta) - coil(icoil)%yy(kk)
   dz = surf%zz(iteta,jzeta) - coil(icoil)%zz(kk)
   
   rr(2) = dx**2 + dy**2 + dz**2 ; rr(1) = sqrt(rr(2)) ; rr(3) = rr(1) * rr(2) ; rr(5) = rr(3) * rr(2)
   
   CHECK( abfield , rr(1).lt.sqrtmachprec, divide by zero : coils intersect plasma boundary )
   CHECK( abfield , rr(3).lt.sqrtmachprec, divide by zero : coils intersect plasma boundary )
   CHECK( abfield , rr(5).lt.sqrtmachprec, divide by zero : coils intersect plasma boundary )
   
   dd(1) = one / rr(1) ; dd(3) = one / rr(3) ; dd(5) = one / rr(5)
   
   tx = coil(icoil)%xt(kk) ! shorthand tangent; 04 Sep 17;
   ty = coil(icoil)%yt(kk)
   tz = coil(icoil)%zt(kk)
   
   xt = surf%xt(iteta,jzeta)
   yt = surf%yt(iteta,jzeta)
   zt = surf%zt(iteta,jzeta)

   dxtx = dx*tx ; dxty = dx*ty ; dxtz = dx*tz
   dytx = dy*tx ; dyty = dy*ty ; dytz = dy*tz
   dztx = dz*tx ; dzty = dz*ty ; dztz = dz*tz
   
   xy = dxty - dytx ! cross products; 04 Sep 17;
   yz = dytz - dzty
   zx = dztx - dxtz
   
   dAdx(kk,1,0) = tx * dd(1) ! Ax; 04 Sep 17;
   dAdx(kk,2,0) = ty * dd(1) ! Ay; 04 Sep 17;
   dAdx(kk,3,0) = tz * dd(1) ! Az; 04 Sep 17;

   dBdx(kk,1,0) = yz * dd(3) ! Bx; 04 Sep 17;
   dBdx(kk,2,0) = zx * dd(3) ! By; 04 Sep 17;
   dBdx(kk,3,0) = xy * dd(3) ! Bz; 04 Sep 17;
   
  !dAdx(kk,1,1) = - ( dyty + dztz ) * dd(3) ! dAx/dx; 04 Sep 17;
  !dAdx(kk,1,2) =     dytx          * dd(3) ! dAx/dy; 04 Sep 17;
  !dAdx(kk,1,3) =     dztx          * dd(3) ! dAx/dz; 04 Sep 17;
   
  !dAdx(kk,2,1) =     dxty          * dd(3) ! dAy/dx; 04 Sep 17;
  !dAdx(kk,2,2) = - ( dxtx + dztz ) * dd(3) ! dAy/dy; 04 Sep 17;
  !dAdx(kk,2,3) =     dzty          * dd(3) ! dAy/dz; 04 Sep 17;
   
  !dAdx(kk,3,1) =     dxtz          * dd(3) ! dAz/dx; 04 Sep 17;
  !dAdx(kk,3,2) =     dytz          * dd(3) ! dAz/dy; 04 Sep 17;
  !dAdx(kk,3,3) = - ( dxtx + dyty ) * dd(3) ! dAz/dz; 04 Sep 17;
   
   rdotdl = dxtx + dyty + dztz ! r dot dl; 04 Sep 17;
   
   dBdx(kk,1,1) = three * yz * dx * dd(5)                                                  ! dBx/dx; 04 Sep 17;
   dBdx(kk,1,2) = three * yz * dy * dd(5) - three * dz * rdotdl * dd(5) + two * tz * dd(3) ! dBx/dy; 04 Sep 17;
   dBdx(kk,1,3) = three * yz * dz * dd(5) + three * dy * rdotdl * dd(5) - two * ty * dd(3) ! dBx/dz; 04 Sep 17;
   
   dBdx(kk,2,1) = three * zx * dx * dd(5) + three * dz * rdotdl * dd(5) - two * tz * dd(3) ! dBy/dx; 04 Sep 17;
   dBdx(kk,2,2) = three * zx * dy * dd(5)                                                  ! dBy/dy; 04 Sep 17;
   dBdx(kk,2,3) = three * zx * dz * dd(5) - three * dx * rdotdl * dd(5) + two * tx * dd(3) ! dBy/dz; 04 Sep 17;
   
   dBdx(kk,3,1) = three * xy * dx * dd(5) - three * dy * rdotdl * dd(5) + two * ty * dd(3) ! dBz/dx; 04 Sep 17;
   dBdx(kk,3,2) = three * xy * dy * dd(5) + three * dx * rdotdl * dd(5) - two * tx * dd(3) ! dBz/dy; 04 Sep 17;
   dBdx(kk,3,3) = three * xy * dz * dd(5)                                                  ! dBz/dz; 04 Sep 17;

  enddo ! end of do kk;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

 !dAdx(0,1,0) = sum( dAdx(1:NS,1,0) ) ! * deltacurveparameter; 04 Sep 17;
 !dAdx(0,2,0) = sum( dAdx(1:NS,2,0) )
 !dAdx(0,3,0) = sum( dAdx(1:NS,3,0) )

  dBdx(0,1,0) = sum( dBdx(1:NS,1,0) ) ! * deltacurveparameter; 04 Sep 17;
  dBdx(0,2,0) = sum( dBdx(1:NS,2,0) )
  dBdx(0,3,0) = sum( dBdx(1:NS,3,0) )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  return
  
end subroutine abfield

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
