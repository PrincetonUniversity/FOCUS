!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

!title (field) ! Computes magnetic field.

!latex \briefly{Computes magnetic field given coil geometry.}

!latex \calledby{\link{bnormal}}
!latex \calls{}

!latex \tableofcontents

!latex \subsection{magnetic field}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

subroutine bfield( icoil, iteta, jzeta, Bx, By, Bz )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  use globals, only : zero, one, two, three, sqrtmachprec, myid, ounit, &
                      surf, Nteta, Nzeta, &
                      Ncoils, coil, deltacurveparameter

  implicit none
  
  include "mpif.h"
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  INTEGER, intent(in ) :: icoil, iteta, jzeta
  REAL   , intent(out) :: Bx, By, Bz
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  INTEGER                           :: ierr, astat, ii
  REAL                              :: dx, dy, dz, rr(1:5), tx, ty, tz, invr3, invr5, rdotdl, yz, zx, xy
  REAL, dimension(1:coil(icoil)%NS) :: dBxx, dBxy, dBxz, dByx, dByy, dByz, dBzx, dBzy, dBzz
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  CHECK( bfield , icoil .lt. 1 .or. icoil .gt. Ncoils, icoil not in right range )
  CHECK( bfield , iteta .lt. 0 .or. iteta .ge. Nteta , iteta not in right range )
  CHECK( bfield , jzeta .lt. 0 .or. jzeta .ge. Nzeta , jzeta not in right range )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  Bx = zero ; By = zero ; Bz = zero
  
  do ii = 1, coil(icoil)%NS
   
   dx = surf(1)%xx(iteta,jzeta) - coil(icoil)%xx(ii) ! distance from evaluation point to curve; 04 Sep 17;
   dy = surf(1)%yy(iteta,jzeta) - coil(icoil)%yy(ii)
   dz = surf(1)%zz(iteta,jzeta) - coil(icoil)%zz(ii)
   
   rr(2) = dx**2 + dy**2 + dz**2 ; rr(1) = sqrt(rr(2)) ; rr(3) = rr(1) * rr(2) ; rr(5) = rr(3) * rr(2)
   
   CHECK( bfield , rr(3).lt.sqrtmachprec, divide by zero : coils intersect plasma boundary )
   CHECK( bfield , rr(5).lt.sqrtmachprec, divide by zero : coils intersect plasma boundary )
   
   invr3 = one / rr(3) ; invr5 = one / rr(5)
   
   tx = coil(icoil)%xt(ii) ! shorthand tangent; 04 Sep 17;
   ty = coil(icoil)%yt(ii)
   tz = coil(icoil)%zt(ii)
   
   yz = dz * ty - dy * tz ! cross products; 04 Sep 17;
   zx = dx * tz - dz * tx
   xy = dy * tx - dx * ty
   
   Bx = Bx + yz * invr3
   By = By + zx * invr3
   Bz = Bz + xy * invr3
   
   rdotdl = dx * tx  +  dy * ty  +  dz * tz ! r dot x'
   
   dBxx(ii) = ( three * yz * dx * invr5                                                  ) * deltacurveparameter !Bx/x
   dBxy(ii) = ( three * yz * dy * invr5 - three * dz * rdotdl * invr5 + two * tz * invr3 ) * deltacurveparameter !Bx/y
   dBxz(ii) = ( three * yz * dz * invr5 + three * dy * rdotdl * invr5 - two * ty * invr3 ) * deltacurveparameter !Bx/z
   
   dByx(ii) = ( three * zx * dx * invr5 + three * dz * rdotdl * invr5 - two * tz * invr3 ) * deltacurveparameter !By/x
   dByy(ii) = ( three * zx * dy * invr5                                                  ) * deltacurveparameter !By/y
   dByz(ii) = ( three * zx * dz * invr5 - three * dx * rdotdl * invr5 + two * tx * invr3 ) * deltacurveparameter !By/z
   
   dBzx(ii) = ( three * xy * dx * invr5 - three * dy * rdotdl * invr5 + two * ty * invr3 ) * deltacurveparameter !Bz/x
   dBzy(ii) = ( three * xy * dy * invr5 + three * dx * rdotdl * invr5 - two * tx * invr3 ) * deltacurveparameter !Bz/y
   dBzz(ii) = ( three * xy * dz * invr5                                                  ) * deltacurveparameter !Bz/z
   
  enddo ! end of do ii;
  
  Bx = Bx * deltacurveparameter
  By = By * deltacurveparameter
  Bz = Bz * deltacurveparameter
  
  FATAL( bfield , .true., DoF has not been defined )
  
  FATAL( bfield , .true., Bx By and Bz need to be given correct dimensions )
  
! Bx(1:ND) = matmul(dBxx, DoF(icoil)%xof) + matmul(dBxy, DoF(icoil)%yof) + matmul(dBxz, DoF(icoil)%zof)
! By(1:ND) = matmul(dByx, DoF(icoil)%xof) + matmul(dByy, DoF(icoil)%yof) + matmul(dByz, DoF(icoil)%zof)
! Bz(1:ND) = matmul(dBzx, DoF(icoil)%xof) + matmul(dBzy, DoF(icoil)%yof) + matmul(dBzz, DoF(icoil)%zof)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  return
  
end subroutine bfield

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

!subroutine bfield1(icoil, iteta, jzeta, Bx, By, Bz, ND)
!
!  use globals, only: coil, DoF, surf, NFcoil, Ncoils, Nteta, Nzeta, &
!                     zero, myid, ounit
!  implicit none
!  include "mpif.h"
!
!  INTEGER, intent(in ) :: icoil, iteta, jzeta, ND
!  REAL, dimension(1:1, 1:ND), intent(inout) :: Bx, By, Bz
!
!  INTEGER              :: ierr, astat, kseg, NS
!  REAL                 :: dlx, dly, dlz, r, rm3, rm5, ltx, lty, ltz, rxp
!  REAL, dimension(1:1, 1:coil(icoil)%NS)   :: dBxx, dBxy, dBxz, dByx, dByy, dByz, dBzx, dBzy, dBzz
!
!  FATAL( bfield1, icoil .lt. 1 .or. icoil .gt. Ncoils, icoil not in right range )
!  FATAL( bfield1, iteta .lt. 0 .or. iteta .gt. Nteta     , iteta not in right range )
!  FATAL( bfield1, jzeta .lt. 0 .or. jzeta .gt. Nzeta     , jzeta not in right range )
!  FATAL( bfield1, ND <= 0, wrong inout dimension of ND )
!  
!  NS = coil(icoil)%NS
!
!  dlx = zero; ltx = zero; Bx = zero
!  dly = zero; lty = zero; By = zero
!  dlz = zero; ltz = zero; Bz = zero
!
!  do kseg = 1, NS
!     
!     dlx = surf(1)%xx(iteta,jzeta) - coil(icoil)%xx(kseg)
!     dly = surf(1)%yy(iteta,jzeta) - coil(icoil)%yy(kseg)
!     dlz = surf(1)%zz(iteta,jzeta) - coil(icoil)%zz(kseg)
!
!     r = sqrt(dlx**2 + dly**2 + dlz**2); rm3 = r**(-3);  rm5 = r**(-5)
!
!     ltx = coil(icoil)%xt(kseg)
!     lty = coil(icoil)%yt(kseg)
!     ltz = coil(icoil)%zt(kseg)
!
!     rxp = dlx*ltx + dly*lty + dlz*ltz !r dot x'
!
!     dBxx(1,kseg) = ( 3*(dlz*lty-dly*ltz)*dlx*rm5                             ) * coil(icoil)%dd(kseg) !Bx/x
!     dBxy(1,kseg) = ( 3*(dlz*lty-dly*ltz)*dly*rm5 - 3*dlz*rxp*rm5 + 2*ltz*rm3 ) * coil(icoil)%dd(kseg) !Bx/y
!     dBxz(1,kseg) = ( 3*(dlz*lty-dly*ltz)*dlz*rm5 + 3*dly*rxp*rm5 - 2*lty*rm3 ) * coil(icoil)%dd(kseg) !Bx/z
!                                                                               
!     dByx(1,kseg) = ( 3*(dlx*ltz-dlz*ltx)*dlx*rm5 + 3*dlz*rxp*rm5 - 2*ltz*rm3 ) * coil(icoil)%dd(kseg) !By/x
!     dByy(1,kseg) = ( 3*(dlx*ltz-dlz*ltx)*dly*rm5                             ) * coil(icoil)%dd(kseg) !By/y
!     dByz(1,kseg) = ( 3*(dlx*ltz-dlz*ltx)*dlz*rm5 - 3*dlx*rxp*rm5 + 2*ltx*rm3 ) * coil(icoil)%dd(kseg) !By/z
!                                                                               
!     dBzx(1,kseg) = ( 3*(dly*ltx-dlx*lty)*dlx*rm5 - 3*dly*rxp*rm5 + 2*lty*rm3 ) * coil(icoil)%dd(kseg) !Bz/x
!     dBzy(1,kseg) = ( 3*(dly*ltx-dlx*lty)*dly*rm5 + 3*dlx*rxp*rm5 - 2*ltx*rm3 ) * coil(icoil)%dd(kseg) !Bz/y
!     dBzz(1,kseg) = ( 3*(dly*ltx-dlx*lty)*dlz*rm5                             ) * coil(icoil)%dd(kseg) !Bz/z
!
!  enddo    ! enddo kseg
!
!  Bx(1:1, 1:ND) = matmul(dBxx, DoF(icoil)%xof) + matmul(dBxy, DoF(icoil)%yof) + matmul(dBxz, DoF(icoil)%zof)
!  By(1:1, 1:ND) = matmul(dByx, DoF(icoil)%xof) + matmul(dByy, DoF(icoil)%yof) + matmul(dByz, DoF(icoil)%zof)
!  Bz(1:1, 1:ND) = matmul(dBzx, DoF(icoil)%xof) + matmul(dBzy, DoF(icoil)%yof) + matmul(dBzz, DoF(icoil)%zof)
!
!  return
!
!end subroutine bfield1
