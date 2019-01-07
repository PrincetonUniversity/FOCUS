
!title (bfield) ! Computes magnetic field.

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
  use globals, only: dp, coil, surf, Ncoils, Nteta, Nzeta, &
                     zero, myid, ounit, Npc, Nfp, pi2, half, two, one, bsconstant
  implicit none
  include "mpif.h"

  INTEGER, intent(in ) :: icoil, iteta, jzeta
  REAL   , intent(out) :: Bx, By, Bz

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  INTEGER              :: ierr, astat, kseg
  REAL                 :: dlx, dly, dlz, rm3, ltx, lty, ltz, rr, r2, m_dot_r, phi, mx, my, mz

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  FATAL( bfield0, icoil .lt. 1 .or. icoil .gt. Ncoils*Npc, icoil not in right range )
  FATAL( bfield0, iteta .lt. 0 .or. iteta .gt. Nteta     , iteta not in right range )
  FATAL( bfield0, jzeta .lt. 0 .or. jzeta .gt. Nzeta     , jzeta not in right range )

  Bx = zero; By = zero; Bz = zero
  
  select case (coil(icoil)%itype)
  !--------------------------------------------------------------------------------------------- 
  case(1)

     dlx = zero; ltx = zero
     dly = zero; lty = zero
     dlz = zero; ltz = zero

     do kseg = 0, coil(icoil)%NS-1

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

     Bx = Bx * coil(icoil)%I * bsconstant
     By = By * coil(icoil)%I * bsconstant
     Bz = Bz * coil(icoil)%I * bsconstant

  !--------------------------------------------------------------------------------------------- 
  case(2)

     dlx = surf(1)%xx(iteta,jzeta) - coil(icoil)%ox
     dly = surf(1)%yy(iteta,jzeta) - coil(icoil)%oy
     dlz = surf(1)%zz(iteta,jzeta) - coil(icoil)%oz
     r2  = dlx**2 + dly**2 + dlz**2
     rm3 = one/(sqrt(r2)*r2)
     mx = sin(coil(icoil)%mt) * cos(coil(icoil)%mp)
     my = sin(coil(icoil)%mt) * sin(coil(icoil)%mp)
     mz = cos(coil(icoil)%mt)
     m_dot_r = mx * dlx + my * dly + mz * dlz

     Bx = 3.0_dp * m_dot_r * rm3 / r2 * dlx - mx * rm3
     By = 3.0_dp * m_dot_r * rm3 / r2 * dly - my * rm3
     Bz = 3.0_dp * m_dot_r * rm3 / r2 * dlz - mz * rm3

     Bx = Bx * coil(icoil)%I * bsconstant
     By = By * coil(icoil)%I * bsconstant
     Bz = Bz * coil(icoil)%I * bsconstant

  !--------------------------------------------------------------------------------------------- 
  case(3)
     ! might be only valid for cylindrical coordinates
     ! Bt = u0*I/(2 pi R)
     phi = ( jzeta + half ) * pi2 / ( Nzeta*Nfp )
     rr = sqrt( surf(1)%xx(iteta,jzeta)**2 + surf(1)%yy(iteta,jzeta)**2 )
     coil(icoil)%Bt = two/rr * coil(icoil)%I * bsconstant

     Bx = - coil(icoil)%Bt * sin(phi)
     By =   coil(icoil)%Bt * cos(phi)
     Bz =   coil(icoil)%Bz 

  !---------------------------------------------------------------------------------------------
  case default
     FATAL(bfield0, .true., not supported coil types)
  end select

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
  use globals, only: dp, coil, DoF, surf, NFcoil, Ncoils, Nteta, Nzeta, &
                     zero, myid, ounit, Npc, one, bsconstant
  implicit none
  include "mpif.h"

  INTEGER, intent(in ) :: icoil, iteta, jzeta, ND
  REAL, dimension(1:1, 1:ND), intent(inout) :: Bx, By, Bz

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  INTEGER              :: ierr, astat, kseg, NS
  REAL                 :: dlx, dly, dlz, r2, rm3, rm5, rm7, m_dot_r, ltx, lty, ltz, rxp, &
                          sinp, sint, cosp, cost, mx, my, mz
  REAL, dimension(1:1, 0:coil(icoil)%NS-1)   :: dBxx, dBxy, dBxz, dByx, dByy, dByz, dBzx, dBzy, dBzz

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  FATAL( bfield1, icoil .lt. 1 .or. icoil .gt. Ncoils*Npc, icoil not in right range )
  FATAL( bfield1, iteta .lt. 0 .or. iteta .gt. Nteta     , iteta not in right range )
  FATAL( bfield1, jzeta .lt. 0 .or. jzeta .gt. Nzeta     , jzeta not in right range )
  FATAL( bfield1, ND <= 0, wrong inout dimension of ND )

  Bx = zero; By = zero; Bz = zero

  select case (coil(icoil)%itype)
  !--------------------------------------------------------------------------------------------- 
  case(1)
     
     NS = coil(icoil)%NS

     dlx = zero; ltx = zero
     dly = zero; lty = zero
     dlz = zero; ltz = zero

     do kseg = 0, NS-1

        dlx = surf(1)%xx(iteta,jzeta) - coil(icoil)%xx(kseg)
        dly = surf(1)%yy(iteta,jzeta) - coil(icoil)%yy(kseg)
        dlz = surf(1)%zz(iteta,jzeta) - coil(icoil)%zz(kseg)

        r2 = dlx**2 + dly**2 + dlz**2; rm3 = one/(sqrt(r2)*r2); rm5 = rm3/r2;

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

     Bx = Bx * coil(icoil)%I * bsconstant
     By = By * coil(icoil)%I * bsconstant
     Bz = Bz * coil(icoil)%I * bsconstant
  !--------------------------------------------------------------------------------------------- 
  case(2)  ! permanent dipoles

     dlx = surf(1)%xx(iteta,jzeta) - coil(icoil)%ox
     dly = surf(1)%yy(iteta,jzeta) - coil(icoil)%oy
     dlz = surf(1)%zz(iteta,jzeta) - coil(icoil)%oz
     r2  = dlx**2 + dly**2 + dlz**2
     rm3 = one/(sqrt(r2)*r2)
     rm5 = rm3/r2
     rm7 = rm5/r2

     cost = cos(coil(icoil)%mt) ; sint = sin(coil(icoil)%mt)
     cosp = cos(coil(icoil)%mp) ; sinp = sin(coil(icoil)%mp)    
     mx = sint*cosp ; my = sint*sinp ; mz = cost
     m_dot_r = mx*dlx + my*dly + mz*dlz

     Bx(1, 1) = 15.0_dp*m_dot_r*dlx*dlx*rm7 - 3.0_dp*mx*dlx*rm5 - 3.0_dp*mx*dlx*rm5 - 3.0_dp*m_dot_r*rm5
     By(1, 1) = 15.0_dp*m_dot_r*dlx*dly*rm7 - 3.0_dp*mx*dly*rm5 - 3.0_dp*my*dlx*rm5
     Bz(1, 1) = 15.0_dp*m_dot_r*dlx*dlz*rm7 - 3.0_dp*mx*dlz*rm5 - 3.0_dp*mz*dlx*rm5

     Bx(1, 2) = 15.0_dp*m_dot_r*dly*dlx*rm7 - 3.0_dp*my*dlx*rm5 - 3.0_dp*mx*dly*rm5
     By(1, 2) = 15.0_dp*m_dot_r*dly*dly*rm7 - 3.0_dp*my*dly*rm5 - 3.0_dp*my*dly*rm5 - 3.0_dp*m_dot_r*rm5
     Bz(1, 2) = 15.0_dp*m_dot_r*dly*dlz*rm7 - 3.0_dp*my*dlz*rm5 - 3.0_dp*mz*dly*rm5

     Bx(1, 3) = 15.0_dp*m_dot_r*dlz*dlx*rm7 - 3.0_dp*mz*dlx*rm5 - 3.0_dp*mx*dlz*rm5
     By(1, 3) = 15.0_dp*m_dot_r*dlz*dly*rm7 - 3.0_dp*mz*dly*rm5 - 3.0_dp*my*dlz*rm5
     Bz(1, 3) = 15.0_dp*m_dot_r*dlz*dlz*rm7 - 3.0_dp*mz*dlz*rm5 - 3.0_dp*mz*dlz*rm5 - 3.0_dp*m_dot_r*rm5 


!!$     Bx(1, 4) = 3.0_dp*dlx*dlx*rm5 - rm3
!!$     By(1, 4) = 3.0_dp*dlx*dly*rm5
!!$     Bz(1, 4) = 3.0_dp*dlx*dlz*rm5
!!$
!!$     Bx(1, 5) = 3.0_dp*dly*dlx*rm5
!!$     By(1, 5) = 3.0_dp*dly*dly*rm5 - rm3
!!$     Bz(1, 5) = 3.0_dp*dly*dlz*rm5
!!$
!!$     Bx(1, 6) = 3.0_dp*dlz*dlx*rm5
!!$     By(1, 6) = 3.0_dp*dlz*dly*rm5
!!$     Bz(1, 6) = 3.0_dp*dlz*dlz*rm5 - rm3
!!$
!!$     Bx = Bx * bsconstant
!!$     By = By * bsconstant
!!$     Bz = Bz * bsconstant

     Bx(1, 4) = 3.0_dp*dlx*( cost*cosp*dlx + cost*sinp*dly - sint*dlz)*rm5 - cost*cosp*rm3
     By(1, 4) = 3.0_dp*dly*( cost*cosp*dlx + cost*sinp*dly - sint*dlz)*rm5 - cost*sinp*rm3 
     Bz(1, 4) = 3.0_dp*dlz*( cost*cosp*dlx + cost*sinp*dly - sint*dlz)*rm5 + sint     *rm3 

     Bx(1, 5) = 3.0_dp*dlx*(-sint*sinp*dlx + sint*cosp*dly           )*rm5 + sint*sinp*rm3
     By(1, 5) = 3.0_dp*dly*(-sint*sinp*dlx + sint*cosp*dly           )*rm5 - sint*cosp*rm3
     Bz(1, 5) = 3.0_dp*dlz*(-sint*sinp*dlx + sint*cosp*dly           )*rm5 

     Bx = Bx * coil(icoil)%I * bsconstant
     By = By * coil(icoil)%I * bsconstant
     Bz = Bz * coil(icoil)%I * bsconstant

  !--------------------------------------------------------------------------------------------- 
  case(3)  ! only for Bz
     
     Bx = zero
     By = zero
     Bz = one

  end select

  return

end subroutine bfield1

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
