!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

!title (field) ! Computes magnetic vector potential and field.

!latex \briefly{Computes magnetic field given coil geometry.}

!latex \calledby{\link{bnormal}}
!latex \calls{}

!latex \tableofcontents

!latex \subsection{magnetic field}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

subroutine brfield( icoil, ii, jj, Ns, BB )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  use globals, only : zero, one, two, three, sqrtmachprec, myid, ounit, &
                      surf, Nt, Nz, &
                      Ncoils, coil
  
  implicit none
  
  include "mpif.h"
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  INTEGER, intent(in ) :: icoil, ii, jj, Ns
  REAL   , intent(out) :: BB(0:Ns-1,1:3)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  INTEGER :: ierr, astat, kk, isurf
  
  REAL    :: dx, dy, dz, rr(1:5), dd(1:5), tx, ty, tz, rdl, yz, zx, xy, xt, yt, zt, xn, yn, zn
  REAL    :: dxtx, dxty, dxtz, dytx, dyty, dytz, dztx, dzty, dztz, td3
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  CHECK( brfield, icoil.lt.1 .or. icoil.gt.Ncoils, icoil not in range )
  CHECK( brfield, ii   .lt.0 .or. ii   .ge.Nt    , ii    not in range )
  CHECK( brfield, jj   .lt.0 .or. jj   .ge.Nz    , jj    not in range )
  
  isurf = 1

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
 !xt = surf(isurf)%xt(ii,jj) ! poloidal tangent to surface (shorthand); required for flux; 16 Nov 17;
 !yt = surf(isurf)%yt(ii,jj)
 !zt = surf(isurf)%zt(ii,jj)
  
 !xn = surf(isurf)%nx(ii,jj) !          normal  to surface (shorthand);
 !yn = surf(isurf)%ny(ii,jj)
 !zn = surf(isurf)%nz(ii,jj)
  
  do kk = 0, Ns-1 ! 12 Nov 17;
   
   dx = surf(isurf)%xx(1,ii,jj) - coil(icoil)%xx(1,kk) ! distance from evaluation point to curve;
   dy = surf(isurf)%xx(2,ii,jj) - coil(icoil)%xx(2,kk)
   dz = surf(isurf)%xx(3,ii,jj) - coil(icoil)%xx(3,kk)
   
   rr(2) = dx**2 + dy**2 + dz**2 ; rr(1) = sqrt(rr(2)) ; rr(3) = rr(1) * rr(2) ! rr(5) = rr(3) * rr(2) !         distance;
   
   FATAL( brfield, rr(1).lt.sqrtmachprec, divide by zero : coils intersect plasma boundary )
   CHECK( brfield, rr(3).lt.sqrtmachprec, divide by zero : coils intersect plasma boundary )
  !CHECK( brfield, rr(5).lt.sqrtmachprec, divide by zero : coils intersect plasma boundary )
   
   ;                             ; dd(1) = one / rr(1) ; dd(3) = one / rr(3)   ! dd(5) = three / rr(5) ! inverse distance;
   
   tx = coil(icoil)%xt(1,kk) ! tangent to coil (shorthand);
   ty = coil(icoil)%xt(2,kk)
   tz = coil(icoil)%xt(3,kk)
   
   dxtx = dx * tx ; dxty = dx * ty ; dxtz = dx * tz
   dytx = dy * tx ; dyty = dy * ty ; dytz = dy * tz
   dztx = dz * tx ; dzty = dz * ty ; dztz = dz * tz
   
  !rdl = dxtx + dyty + dztz ! r dot dl; required for flux; 16 Nov 17;
   
  !dFdx(kk,0) = (   tx*xt             +   ty*yt             +   tz*zt             ) * dd(1) !  flux   ;
   
  !dFdx(kk,1) = ( ( dxtx - rdl ) * xt + ( dxty       ) * yt + ( dxtz       ) * zt ) * dd(3) ! dflux/dx;
  !dFdx(kk,2) = ( ( dytx       ) * xt + ( dyty - rdl ) * yt + ( dytz       ) * zt ) * dd(3) ! dflux/dy;
  !dFdx(kk,3) = ( ( dztx       ) * xt + ( dzty       ) * yt + ( dztz - rdl ) * zt ) * dd(3) ! dflux/dz;
   
   yz = dzty - dytz
   zx = dxtz - dztx
   xy = dytx - dxty

!  td3 = two * dd(3) ! seems unneccessary; 16 Nov 17;
   td3 =       dd(3) !                   ;          ; 28 Nov 17;
   
   BB(kk,1) =      yz                   * td3 ! Bx;
   BB(kk,2) =              zx           * td3 ! By;
   BB(kk,3) =                      xy   * td3 ! Bz;
   
  enddo ! end of do kk;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

  return
  
end subroutine brfield

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
