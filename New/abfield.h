!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

!title (field) ! Computes magnetic vector potential and field.

!latex \briefly{Computes magnetic field given coil geometry.}

!latex \calledby{\link{bnormal}}
!latex \calls{}

!latex \tableofcontents

!latex \subsection{magnetic field}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

subroutine abfield( icoil, ii, jj, Ns, dFdx, dBdx )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  use globals, only : zero, one, two, three, sqrtmachprec, myid, ounit, &
                      surf, Nt, Nz, &
                      Ncoils, coil
  
  implicit none
  
  include "mpif.h"
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  INTEGER, intent(in ) :: icoil, ii, jj, Ns
  REAL   , intent(out) :: dFdx(0:Ns,0:3), dBdx(0:Ns,0:3)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  INTEGER :: ierr, astat, kk
  
  REAL    :: dx, dy, dz, rr(1:5), dd(1:5), tx, ty, tz, rdl, yz, zx, xy, xt, yt, zt, xn, yn, zn
  REAL    :: dxtx, dxty, dxtz, dytx, dyty, dytz, dztx, dzty, dztz, td3
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  CHECK( abfield, icoil.lt.1 .or. icoil.gt.Ncoils, icoil not in range )
  CHECK( abfield, ii   .lt.0 .or. ii   .ge.Nt    , ii    not in range )
  CHECK( abfield, jj   .lt.0 .or. jj   .ge.Nz    , jj    not in range )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  xt = surf%xt(ii,jj) ! poloidal tangent to surface (shorthand);
  yt = surf%yt(ii,jj)
  zt = surf%zt(ii,jj)

  xn = surf%nx(ii,jj) !          normal  to surface (shorthand);
  yn = surf%ny(ii,jj)
  zn = surf%nz(ii,jj)

  do kk = 1, Ns
   
   dx = surf%xx(ii,jj) - coil(icoil)%xx(kk) ! distance from evaluation point to curve;
   dy = surf%yy(ii,jj) - coil(icoil)%yy(kk)
   dz = surf%zz(ii,jj) - coil(icoil)%zz(kk)
   
   rr(2) = dx**2 + dy**2 + dz**2 ; rr(1) = sqrt(rr(2)) ; rr(3) = rr(1) * rr(2) ; rr(5) = rr(3) * rr(2) !         distance;
   
   FATAL( abfield, rr(1).lt.sqrtmachprec, divide by zero : coils intersect plasma boundary )
   CHECK( abfield, rr(3).lt.sqrtmachprec, divide by zero : coils intersect plasma boundary )
   CHECK( abfield, rr(5).lt.sqrtmachprec, divide by zero : coils intersect plasma boundary )
   
   ;                             ; dd(1) = one / rr(1) ; dd(3) = one / rr(3)   ; dd(5) = three / rr(5) ! inverse distance;
   
   tx = coil(icoil)%xt(kk) ! tangent to coil (shorthand);
   ty = coil(icoil)%yt(kk)
   tz = coil(icoil)%zt(kk)
   
   dxtx = dx*tx ; dxty = dx*ty ; dxtz = dx*tz
   dytx = dy*tx ; dyty = dy*ty ; dytz = dy*tz
   dztx = dz*tx ; dzty = dz*ty ; dztz = dz*tz
   
   rdl = dxtx + dyty + dztz ! r dot dl;
   
   dFdx(kk,0) = (   tx*xt             +   ty*yt             +   tz*zt             ) * dd(1) !  flux   ;
   
   dFdx(kk,1) = ( ( dxtx - rdl ) * xt + ( dxty       ) * yt + ( dxtz       ) * zt ) * dd(3) ! dflux/dx;
   dFdx(kk,2) = ( ( dytx       ) * xt + ( dyty - rdl ) * yt + ( dytz       ) * zt ) * dd(3) ! dflux/dy;
   dFdx(kk,3) = ( ( dztx       ) * xt + ( dzty       ) * yt + ( dztz - rdl ) * zt ) * dd(3) ! dflux/dz;
   
   yz = dzty - dytz
   zx = dxtz - dztx
   xy = dytx - dxty

   td3 = two * dd(3)
   
   dBdx(kk,0) = ( xn*yz + yn*zx + zn*xy ) * td3 ! Bn;

   dBdx(kk,1) = xn * ( dd(5) * ( yz * dx            )            ) &
              + yn * ( dd(5) * ( zx * dx + dz * rdl ) - tz * td3 ) &
              + zn * ( dd(5) * ( xy * dx - dy * rdl ) + ty * td3 ) ! dBn/dx;

   dBdx(kk,2) = xn * ( dd(5) * ( yz * dy - dz * rdl ) + tz * td3 ) &
              + yn * ( dd(5) * ( zx * dy            )            ) &
              + zn * ( dd(5) * ( xy * dy + dx * rdl ) - tx * td3 ) ! dBn/dy;

   dBdx(kk,3) = xn * ( dd(5) * ( yz * dz + dy * rdl ) - ty * td3 ) &
              + yn * ( dd(5) * ( zx * dz - dx * rdl ) + tx * td3 ) &
              + zn * ( dd(5) * ( xy * dz            )            ) ! dBn/dz;


   dBdx(kk,1:3) = dBdx(kk,1:3) * two ! WHERE DID THIS FACTOR OF TWO COME FROM;

  enddo ! end of do kk;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

  return
  
end subroutine abfield

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

!dAdx(kk,1,0) = tx * dd(1) ! Ax;
!dAdx(kk,2,0) = ty * dd(1) ! Ay;
!dAdx(kk,3,0) = tz * dd(1) ! Az;

!dAdx(kk,1,1) = - ( dyty + dztz ) * dd(3) ! dAx/dx;
!dAdx(kk,1,2) =     dytx          * dd(3) ! dAx/dy;
!dAdx(kk,1,3) =     dztx          * dd(3) ! dAx/dz;

!dAdx(kk,2,1) =     dxty          * dd(3) ! dAy/dx;
!dAdx(kk,2,2) = - ( dxtx + dztz ) * dd(3) ! dAy/dy;
!dAdx(kk,2,3) =     dzty          * dd(3) ! dAy/dz;

!dAdx(kk,3,1) =     dxtz          * dd(3) ! dAz/dx;
!dAdx(kk,3,2) =     dytz          * dd(3) ! dAz/dy;
!dAdx(kk,3,3) = - ( dxtx + dyty ) * dd(3) ! dAz/dz;

!dAdx(kk,1,1) = ( - rdl + dxtx               ) * dd(3) ! dAx/dx;
!dAdx(kk,1,2) = (                   dytx        ) * dd(3) ! dAx/dy;
!dAdx(kk,1,3) = (                          dztx ) * dd(3) ! dAx/dz;

!dAdx(kk,2,1) = (            dxty               ) * dd(3) ! dAy/dx;
!dAdx(kk,2,2) = ( - rdl        + dyty        ) * dd(3) ! dAy/dy;
!dAdx(kk,2,3) = (                          dzty ) * dd(3) ! dAy/dz;

!dAdx(kk,3,1) = (            dxtz               ) * dd(3) ! dAz/dx;
!dAdx(kk,3,2) = (                   dytz        ) * dd(3) ! dAz/dy;
!dAdx(kk,3,3) = ( - rdl               + dztz ) * dd(3) ! dAz/dz;

!dBdx(kk,1,0) = yz * dd(3) ! Bx;
!dBdx(kk,2,0) = zx * dd(3) ! By;
!dBdx(kk,3,0) = xy * dd(3) ! Bz;

!dBdx(kk,1,1) = three * yz * dx * dd(5)                                                  ! dBx/dx;
!dBdx(kk,1,2) = three * yz * dy * dd(5) - three * dz * rdl * dd(5) + two * tz * dd(3) ! dBx/dy;
!dBdx(kk,1,3) = three * yz * dz * dd(5) + three * dy * rdl * dd(5) - two * ty * dd(3) ! dBx/dz;

!dBdx(kk,2,1) = three * zx * dx * dd(5) + three * dz * rdl * dd(5) - two * tz * dd(3) ! dBy/dx;
!dBdx(kk,2,2) = three * zx * dy * dd(5)                                                  ! dBy/dy;
!dBdx(kk,2,3) = three * zx * dz * dd(5) - three * dx * rdl * dd(5) + two * tx * dd(3) ! dBy/dz;

!dBdx(kk,3,1) = three * xy * dx * dd(5) - three * dy * rdl * dd(5) + two * ty * dd(3) ! dBz/dx;
!dBdx(kk,3,2) = three * xy * dy * dd(5) + three * dx * rdl * dd(5) - two * tx * dd(3) ! dBz/dy;
!dBdx(kk,3,3) = three * xy * dz * dd(5)                                                  ! dBz/dz;
