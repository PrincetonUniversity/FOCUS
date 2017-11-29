!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

!title (field) ! Computes magnetic vector potential and field.

!latex \briefly{Computes magnetic field given coil geometry.}

!latex \calledby{\link{bnormal}}
!latex \calls{}

!latex \tableofcontents

!latex \subsection{magnetic field}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

subroutine abfield( isurf, icoil, ii, jj, Ns, dFdx, dBdx )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  use globals, only : zero, one, two, three, sqrtmachprec, myid, ounit, &
                      surf, Nt, Nz, &
                      Ncoils, coil
  
  implicit none
  
  include "mpif.h"
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  INTEGER, intent(in ) :: icoil, ii, jj, Ns, isurf
  REAL   , intent(out) :: dFdx(0:Ns-1,0:3), dBdx(0:Ns-1,0:3)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  INTEGER :: ierr, astat, kk
  
  REAL    :: dx, dy, dz, rr(1:5), dd(1:5), tx, ty, tz, rdl, yz, zx, xy, xt, yt, zt, xn, yn, zn
  REAL    :: dxtx, dxty, dxtz, dytx, dyty, dytz, dztx, dzty, dztz, td3
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  CHECK( abfield, icoil.lt.1 .or. icoil.gt.Ncoils, icoil not in range )
  CHECK( abfield, ii   .lt.0 .or. ii   .ge.Nt    , ii    not in range )
  CHECK( abfield, jj   .lt.0 .or. jj   .ge.Nz    , jj    not in range )
  CHECK( abfield, isurf.lt.1 .or. isurf.gt.2     , isurf not in range ) ! 28 Nov 17;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

  xt = surf(isurf)%xu(1,ii,jj) ! poloidal tangent to surface (shorthand);
  yt = surf(isurf)%xu(2,ii,jj)
  zt = surf(isurf)%xu(3,ii,jj)
  
  xn = surf(isurf)%nn(1,ii,jj) !          normal  to surface (shorthand);
  yn = surf(isurf)%nn(2,ii,jj)
  zn = surf(isurf)%nn(3,ii,jj)
  
  do kk = 0, Ns-1 ! 12 Nov 17;
   
   dx = surf(isurf)%xx(1,ii,jj) - coil(icoil)%xx(1,kk) ! distance from evaluation point to curve;
   dy = surf(isurf)%xx(2,ii,jj) - coil(icoil)%xx(2,kk)
   dz = surf(isurf)%xx(3,ii,jj) - coil(icoil)%xx(3,kk)
   
   rr(2) = dx**2 + dy**2 + dz**2 ; rr(1) = sqrt(rr(2)) ; rr(3) = rr(1) * rr(2) ; rr(5) = rr(3) * rr(2) !         distance;
   
   FATAL( abfield, rr(1).lt.sqrtmachprec, divide by zero : coils intersect plasma boundary )
   CHECK( abfield, rr(3).lt.sqrtmachprec, divide by zero : coils intersect plasma boundary )
   CHECK( abfield, rr(5).lt.sqrtmachprec, divide by zero : coils intersect plasma boundary )
   
   ;                             ; dd(1) = one / rr(1) ; dd(3) = one / rr(3)   ; dd(5) = three / rr(5) ! inverse distance;
   
   tx = coil(icoil)%xt(1,kk) ! tangent to coil (shorthand);
   ty = coil(icoil)%xt(2,kk)
   tz = coil(icoil)%xt(3,kk)
   
   dxtx = dx * tx ; dxty = dx * ty ; dxtz = dx * tz
   dytx = dy * tx ; dyty = dy * ty ; dytz = dy * tz
   dztx = dz * tx ; dzty = dz * ty ; dztz = dz * tz
   
   rdl = dxtx + dyty + dztz ! r dot dl;
   
   dFdx(kk,0) = (   tx*xt             +   ty*yt             +   tz*zt             ) * dd(1) !  flux   ;
   
   dFdx(kk,1) = ( ( dxtx - rdl ) * xt + ( dxty       ) * yt + ( dxtz       ) * zt ) * dd(3) ! dflux/dx;
   dFdx(kk,2) = ( ( dytx       ) * xt + ( dyty - rdl ) * yt + ( dytz       ) * zt ) * dd(3) ! dflux/dy;
   dFdx(kk,3) = ( ( dztx       ) * xt + ( dzty       ) * yt + ( dztz - rdl ) * zt ) * dd(3) ! dflux/dz;
   
   yz = dzty - dytz
   zx = dxtz - dztx
   xy = dytx - dxty
   
   td3 =       dd(3)

   dBdx(kk,0  ) = ( xn*yz + yn*zx + zn*xy ) * td3 ! Bn;
 
   td3 = two * dd(3)
  
   dBdx(kk,1  ) = xn * ( dd(5) * ( yz * dx            )            ) &
                + yn * ( dd(5) * ( zx * dx + dz * rdl ) - tz * td3 ) &
                + zn * ( dd(5) * ( xy * dx - dy * rdl ) + ty * td3 ) ! dBn/dx;

   dBdx(kk,2  ) = xn * ( dd(5) * ( yz * dy - dz * rdl ) + tz * td3 ) &
                + yn * ( dd(5) * ( zx * dy            )            ) &
                + zn * ( dd(5) * ( xy * dy + dx * rdl ) - tx * td3 ) ! dBn/dy;

   dBdx(kk,3  ) = xn * ( dd(5) * ( yz * dz + dy * rdl ) - ty * td3 ) &
                + yn * ( dd(5) * ( zx * dz - dx * rdl ) + tx * td3 ) &
                + zn * ( dd(5) * ( xy * dz            )            ) ! dBn/dz;

!  dBdx(kk,1:3) = dBdx(kk,1:3) * two ! WHERE DID THIS FACTOR OF TWO COME FROM; 12 Nov 17;

!  td3 = two * dd(3) ! WHERE DID THE FACTOR OF TWO COME FROM; 12 Nov 17;
   
!  dBdx(kk,0  ) = ( xn*yz + yn*zx + zn*xy ) * td3 ! Bn;
   
!  dBdx(kk,1  ) = xn * ( dd(5) * ( yz * dx            )            ) &
!               + yn * ( dd(5) * ( zx * dx + dz * rdl ) - tz * td3 ) &
!               + zn * ( dd(5) * ( xy * dx - dy * rdl ) + ty * td3 ) ! dBn/dx;

!  dBdx(kk,2  ) = xn * ( dd(5) * ( yz * dy - dz * rdl ) + tz * td3 ) &
!               + yn * ( dd(5) * ( zx * dy            )            ) &
!               + zn * ( dd(5) * ( xy * dy + dx * rdl ) - tx * td3 ) ! dBn/dy;

!  dBdx(kk,3  ) = xn * ( dd(5) * ( yz * dz + dy * rdl ) - ty * td3 ) &
!               + yn * ( dd(5) * ( zx * dz - dx * rdl ) + tx * td3 ) &
!               + zn * ( dd(5) * ( xy * dz            )            ) ! dBn/dz;

!  dBdx(kk,1:3) = dBdx(kk,1:3) * two ! WHERE DID THIS FACTOR OF TWO COME FROM; 12 Nov 17;

  enddo ! end of do kk;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

  return
  
end subroutine abfield

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
