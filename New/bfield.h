
!title (field) ! Computes magnetic field.

!latex \briefly{Computes magnetic field given coil geometry.}

!latex \calledby{\link{bnormal}}
!latex \calls{}

!latex \tableofcontents

!latex \subsection{magnetic field}
!latex \bi
!latex \item The magnetic field of filamentary coils is calculated bt Biot-Savart Law, involving a line integral.
!latex J. Hanson and S. Hirshman had a better representation for straight segments to avoid unnecessary sigularities
!latex and improve numerical error at points neary the coil.
!latex \item But currently, we use the normal expression of Biot-Savart Law and derivatives of B with repsect to 
!latex x, y, z is also calculated.
!latex \item Later, error analysis and comparison to Hanson's method should be carried out. 
!latex \ei

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine bfield0(icoil, iteta, jzeta, Bx, By, Bz)
!------------------------------------------------------------------------------------------------------ 
! DATE:  06/15/2016; 03/26/2017
! calculate the magnetic field of icoil using manually discretized coils. 
! Biot-Savart constant and currents are not included for later simplication. 
! Be careful if coils have different resolutions.
!------------------------------------------------------------------------------------------------------   
  use globals, only: coil, surf, Ncoils, Nteta, Nzeta, &
                     zero, myid, ounit
  implicit none
  include "mpif.h"

  INTEGER, intent(in ) :: icoil, iteta, jzeta
  REAL   , intent(out) :: Bx, By, Bz

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  INTEGER              :: ierr, astat, kseg
  REAL                 :: dlx, dly, dlz, rm3, ltx, lty, ltz

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  FATAL( bfield0, icoil .lt. 1 .or. icoil .gt. Ncoils, icoil not in right range )
  FATAL( bfield0, iteta .lt. 0 .or. iteta .gt. Nteta , iteta not in right range )
  FATAL( bfield0, jzeta .lt. 0 .or. jzeta .gt. Nzeta , jzeta not in right range )
  
  dlx = zero; ltx = zero; Bx = zero
  dly = zero; lty = zero; By = zero
  dlz = zero; ltz = zero; Bz = zero

  do kseg = 1, coil(icoil)%NS
        
   dlx = surf(1)%xx(iteta,jzeta) - coil(icoil)%xx(kseg)
   dly = surf(1)%yy(iteta,jzeta) - coil(icoil)%yy(kseg)
   dlz = surf(1)%zz(iteta,jzeta) - coil(icoil)%zz(kseg)
   rm3 = (sqrt(dlx**2 + dly**2 + dlz**2))**(-3)

   ltx = coil(icoil)%xt(kseg)
   lty = coil(icoil)%yt(kseg)
   ltz = coil(icoil)%zt(kseg)
   
   Bx = Bx + ( dlz*lty - dly*ltz ) * rm3 * coil(icoil)%dd(kseg)
   By = By + ( dlx*ltz - dlz*ltx ) * rm3 * coil(icoil)%dd(kseg)
   Bz = Bz + ( dly*ltx - dlx*lty ) * rm3 * coil(icoil)%dd(kseg)

  enddo    ! enddo kseg

  return

end subroutine bfield0

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine bfield1(icoil, iteta, jzeta, Bx, By, Bz, ND)
!------------------------------------------------------------------------------------------------------ 
! DATE:  06/15/2016; 03/26/2017
! calculate the magnetic field and the first dirivatives of icoil using manually discretized coils;
! Biot-Savart constant and currents are not included for later simplication;
! Discretizing factor is includeed; coil(icoil)%dd(kseg)
!------------------------------------------------------------------------------------------------------   
  use globals, only: coil, DoF, fouCoil, surf, NFcoil, Ncoils, Cdof, Nteta, Nzeta, &
                     zero, myid, ounit
  implicit none
  include "mpif.h"

  INTEGER, intent(in ) :: icoil, iteta, jzeta, ND
  REAL, dimension(1:1, 1:ND), intent(inout) :: Bx, By, Bz

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  INTEGER              :: ierr, astat, kseg, NS
  REAL                 :: dlx, dly, dlz, r, rm3, rm5, ltx, lty, ltz, rxp
  REAL, dimension(1:1, 1:coil(icoil)%NS)   :: dBxx, dBxy, dBxz, dByx, dByy, dByz, dBzx, dBzy, dBzz

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  FATAL( bfield1, icoil .lt. 1 .or. icoil .gt. Ncoils, icoil not in right range )
  FATAL( bfield1, iteta .lt. 0 .or. iteta .gt. Nteta , iteta not in right range )
  FATAL( bfield1, jzeta .lt. 0 .or. jzeta .gt. Nzeta , jzeta not in right range )
  FATAL( bfield1, ND <= 0, wrong inout dimension of ND )
  
  NS = coil(icoil)%NS

  dlx = zero; ltx = zero; Bx = zero
  dly = zero; lty = zero; By = zero
  dlz = zero; ltz = zero; Bz = zero

  do kseg = 1, NS
     
     dlx = surf(1)%xx(iteta,jzeta) - coil(icoil)%xx(kseg)
     dly = surf(1)%yy(iteta,jzeta) - coil(icoil)%yy(kseg)
     dlz = surf(1)%zz(iteta,jzeta) - coil(icoil)%zz(kseg)

     r = sqrt(dlx**2 + dly**2 + dlz**2); rm3 = r**(-3);  rm5 = r**(-5)

     ltx = coil(icoil)%xt(kseg)
     lty = coil(icoil)%yt(kseg)
     ltz = coil(icoil)%zt(kseg)

     rxp = dlx*ltx + dly*lty + dlz*ltz !r dot x'

     dBxx(1,kseg) = ( 3*(dlz*lty-dly*ltz)*dlx*rm5                             ) * coil(icoil)%dd(kseg) !Bx/x
     dBxy(1,kseg) = ( 3*(dlz*lty-dly*ltz)*dly*rm5 - 3*dlz*rxp*rm5 + 2*ltz*rm3 ) * coil(icoil)%dd(kseg) !Bx/y
     dBxz(1,kseg) = ( 3*(dlz*lty-dly*ltz)*dlz*rm5 + 3*dly*rxp*rm5 - 2*lty*rm3 ) * coil(icoil)%dd(kseg) !Bx/z
                                                                               
     dByx(1,kseg) = ( 3*(dlx*ltz-dlz*ltx)*dlx*rm5 + 3*dlz*rxp*rm5 - 2*ltz*rm3 ) * coil(icoil)%dd(kseg) !By/x
     dByy(1,kseg) = ( 3*(dlx*ltz-dlz*ltx)*dly*rm5                             ) * coil(icoil)%dd(kseg) !By/y
     dByz(1,kseg) = ( 3*(dlx*ltz-dlz*ltx)*dlz*rm5 - 3*dlx*rxp*rm5 + 2*ltx*rm3 ) * coil(icoil)%dd(kseg) !By/z
                                                                               
     dBzx(1,kseg) = ( 3*(dly*ltx-dlx*lty)*dlx*rm5 - 3*dly*rxp*rm5 + 2*lty*rm3 ) * coil(icoil)%dd(kseg) !Bz/x
     dBzy(1,kseg) = ( 3*(dly*ltx-dlx*lty)*dly*rm5 + 3*dlx*rxp*rm5 - 2*ltx*rm3 ) * coil(icoil)%dd(kseg) !Bz/y
     dBzz(1,kseg) = ( 3*(dly*ltx-dlx*lty)*dlz*rm5                             ) * coil(icoil)%dd(kseg) !Bz/z

  enddo    ! enddo kseg

  Bx(1:1, 1:ND) = matmul(dBxx, DoF(icoil)%xof) + matmul(dBxy, DoF(icoil)%yof) + matmul(dBxz, DoF(icoil)%zof)
  By(1:1, 1:ND) = matmul(dByx, DoF(icoil)%xof) + matmul(dByy, DoF(icoil)%yof) + matmul(dByz, DoF(icoil)%zof)
  Bz(1:1, 1:ND) = matmul(dBzx, DoF(icoil)%xof) + matmul(dBzy, DoF(icoil)%yof) + matmul(dBzz, DoF(icoil)%zof)


  return

end subroutine bfield1

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

SUBROUTINE test
  
  use globals, only: coil, DoF, surf, NFcoil, Ncoils, Nteta, Nzeta, Ndof, bnorm, t1B, xdof, &
                     zero, pi2, myid, ounit
  implicit none
  include "mpif.h"

  INTEGER :: icoil, iteta, jzeta, ierr, astat, Cdof
  
  Cdof = 6*NFcoil + 3
  icoil = 1 ; iteta = 0; jzeta = 0

  call bnormal(0)
  if(myid == 0) write(ounit,*) "bnorm = ", bnorm

  if(.not. allocated(coil(icoil)%Bx)) then 
     SALLOCATE( coil(icoil)%Bx, (0:Cdof, 0:Cdof), zero )
     SALLOCATE( coil(icoil)%By, (0:Cdof, 0:Cdof), zero )
     SALLOCATE( coil(icoil)%Bz, (0:Cdof, 0:Cdof), zero )
  endif

  call bfield0(icoil, iteta, jzeta, coil(icoil)%Bx(0,0), coil(icoil)%By(0,0), &
          coil(icoil)%Bz(0,0))
  
  call bfield1(icoil, iteta, jzeta, coil(icoil)%Bx(1:Cdof,0), coil(icoil)%By(1:Cdof,0), &
          coil(icoil)%Bz(1:Cdof,0), DoF(icoil)%ND)

  if(myid == 0) write(ounit,*) "Bx/x = ", coil(icoil)%Bx(0:Cdof,0)
  if(myid == 0) write(ounit,*) "By/x = ", coil(icoil)%By(0:Cdof,0)
  if(myid == 0) write(ounit,*) "Bz/x = ", coil(icoil)%Bz(0:Cdof,0)

  if(myid == 0) write(ounit,*) coil(1)%I, xdof(1)
  !write(ounit,*) "xof(1) = ", DoF(icoil)%xof(1, 1:Cdof)
  !write(ounit,*) "xof(32) = ", DoF(icoil)%xof(32, 1:Cdof) 
  !call bnormal(1)
  !write(ounit,*) "OB/x = ", t1B(1:Ndof)

end SUBROUTINE test
