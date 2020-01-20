
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

subroutine bfield0(icoil, x, y, z, tBx, tBy, tBz)
!------------------------------------------------------------------------------------------------------ 
! DATE:  06/15/2016; 03/26/2017
! calculate the magnetic field of icoil using manually discretized coils. 
! Biot-Savart constant and currents are not included for later simplication. 
! Be careful if coils have different resolutions.
!------------------------------------------------------------------------------------------------------   
  use globals, only: dp, coil, surf, Ncoils, Nteta, Nzeta, cosnfp, sinnfp, &
                     zero, myid, ounit, Nfp, pi2, half, two, one, bsconstant
  use mpi
  implicit none

  INTEGER, intent(in ) :: icoil
  REAL,    intent(in ) :: x, y, z
  REAL   , intent(out) :: tBx, tBy, tBz

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  INTEGER              :: ierr, astat, kseg, ip, is, cs, Npc
  REAL                 :: dlx, dly, dlz, rm3, ltx, lty, ltz, rr, r2, m_dot_r, &
       &                  mx, my, mz, xx, yy, zz, Bx, By, Bz

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  FATAL( bfield0, icoil .lt. 1 .or. icoil .gt. Ncoils, icoil not in right range )
  ! initialization
  Npc = 1 ; cs = 0
  tBx = zero ; tBy = zero ; tBz = zero
  dlx = zero ; dly = zero ; dlz = zero
  ltx = zero ; lty = zero ; ltz = zero
  ! check if the coil is stellarator symmetric
  select case (coil(icoil)%symm) 
  case ( 0 )
     cs  = 0
     Npc = 1
  case ( 1 )
     cs  = 0
     Npc = Nfp
  case ( 2) 
     cs  = 1
     Npc = Nfp
  end select
  ! periodicity and stellarator symmetry
  do ip = 1, Npc
     do is = 0, cs
        ! find the point on plasma by rotating in reverse direction. + symmetric
        xx = ( x*cosnfp(ip) + y*sinnfp(ip) )
        yy = (-x*sinnfp(ip) + y*cosnfp(ip) ) * (-1)**is
        zz =  z * (-1)**is
        Bx = zero; By = zero; Bz = zero
        select case (coil(icoil)%type)
        ! Fourier coils
        case(1)
           ! Biot-Savart law
           do kseg = 0, coil(icoil)%NS-1
              dlx = xx - coil(icoil)%xx(kseg)
              dly = yy - coil(icoil)%yy(kseg)
              dlz = zz - coil(icoil)%zz(kseg)
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
        ! magnetic dipoles
        case(2)
           ! Biot-Savart law
           dlx = xx - coil(icoil)%ox
           dly = yy - coil(icoil)%oy
           dlz = zz - coil(icoil)%oz
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
        ! toroidal field and verticle field
        case(3)
           ! might be only valid for cylindrical coordinates
           ! Bt = u0*I/(2 pi R)
           rr = sqrt( xx**2 + yy**2 )
           coil(icoil)%Bt = two/rr * coil(icoil)%I * bsconstant
           Bx = - coil(icoil)%Bt * yy/rr
           By =   coil(icoil)%Bt * xx/rr
           Bz =   coil(icoil)%Bz 
        case default
           FATAL(bfield0, .true., not supported coil types)
        end select     
        ! sum all the contributions
        tBx = tBx + (Bx*cosnfp(ip) - By*sinnfp(ip))*(-1)**is
        tBy = tBy + (By*cosnfp(ip) + Bx*sinnfp(ip))
        tBz = tBz +  Bz
     enddo
  enddo

  return

end subroutine bfield0

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine bfield1(icoil, xx, yy, zz, Bx, By, Bz, ND)
!------------------------------------------------------------------------------------------------------ 
! DATE:  06/15/2016; 03/26/2017
! calculate the magnetic field and the first dirivatives of icoil using manually discretized coils;
! Biot-Savart constant and currents are not included for later simplication;
! Discretizing factor is includeed; coil(icoil)%dd(kseg)
!------------------------------------------------------------------------------------------------------   
  use globals, only: dp, coil, DoF, surf, NFcoil, Ncoils, Nteta, Nzeta, &
                     zero, myid, ounit, Nfp, one, bsconstant
  implicit none
  include "mpif.h"

  INTEGER, intent(in ) :: icoil, ND
  REAL,    intent(in ) :: xx, yy, zz
  REAL, dimension(1:1, 1:ND), intent(inout) :: Bx, By, Bz

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  INTEGER              :: ierr, astat, kseg, NS
  REAL                 :: dlx, dly, dlz, r2, rm3, rm5, rm7, m_dot_r, ltx, lty, ltz, rxp, &
                          sinp, sint, cosp, cost, mx, my, mz
  REAL, dimension(1:1, 0:coil(icoil)%NS-1)   :: dBxx, dBxy, dBxz, dByx, dByy, dByz, dBzx, dBzy, dBzz

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  FATAL( bfield1, icoil .lt. 1 .or. icoil .gt. Ncoils, icoil not in right range )
  FATAL( bfield1, ND <= 0, wrong inout dimension of ND )

  Bx = zero; By = zero; Bz = zero
  dlx = zero ; dly = zero ; dlz = zero
  ltx = zero ; lty = zero ; ltz = zero

  select case (coil(icoil)%type)
  !--------------------------------------------------------------------------------------------- 
  case(1)
     
     NS = coil(icoil)%NS

     do kseg = 0, NS-1

        dlx = xx - coil(icoil)%xx(kseg)
        dly = yy - coil(icoil)%yy(kseg)
        dlz = zz - coil(icoil)%zz(kseg)

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

     dlx = xx - coil(icoil)%ox
     dly = yy - coil(icoil)%oy
     dlz = zz - coil(icoil)%oz
     r2  = dlx**2 + dly**2 + dlz**2
     rm3 = one/(sqrt(r2)*r2)
     rm5 = rm3/r2
     rm7 = rm5/r2

     cost = cos(coil(icoil)%mt) ; sint = sin(coil(icoil)%mt)
     cosp = cos(coil(icoil)%mp) ; sinp = sin(coil(icoil)%mp)    
     mx = sint*cosp ; my = sint*sinp ; mz = cost
     m_dot_r = mx*dlx + my*dly + mz*dlz

#ifdef dposition
     ! dipole position is variable
     Bx(1, 1) = 15.0_dp*m_dot_r*dlx*dlx*rm7 - 3.0_dp*mx*dlx*rm5 - 3.0_dp*mx*dlx*rm5 - 3.0_dp*m_dot_r*rm5
     By(1, 1) = 15.0_dp*m_dot_r*dlx*dly*rm7 - 3.0_dp*mx*dly*rm5 - 3.0_dp*my*dlx*rm5
     Bz(1, 1) = 15.0_dp*m_dot_r*dlx*dlz*rm7 - 3.0_dp*mx*dlz*rm5 - 3.0_dp*mz*dlx*rm5

     Bx(1, 2) = 15.0_dp*m_dot_r*dly*dlx*rm7 - 3.0_dp*my*dlx*rm5 - 3.0_dp*mx*dly*rm5
     By(1, 2) = 15.0_dp*m_dot_r*dly*dly*rm7 - 3.0_dp*my*dly*rm5 - 3.0_dp*my*dly*rm5 - 3.0_dp*m_dot_r*rm5
     Bz(1, 2) = 15.0_dp*m_dot_r*dly*dlz*rm7 - 3.0_dp*my*dlz*rm5 - 3.0_dp*mz*dly*rm5

     Bx(1, 3) = 15.0_dp*m_dot_r*dlz*dlx*rm7 - 3.0_dp*mz*dlx*rm5 - 3.0_dp*mx*dlz*rm5
     By(1, 3) = 15.0_dp*m_dot_r*dlz*dly*rm7 - 3.0_dp*mz*dly*rm5 - 3.0_dp*my*dlz*rm5
     Bz(1, 3) = 15.0_dp*m_dot_r*dlz*dlz*rm7 - 3.0_dp*mz*dlz*rm5 - 3.0_dp*mz*dlz*rm5 - 3.0_dp*m_dot_r*rm5 

     Bx(1, 4) = 3.0_dp*dlx*( cost*cosp*dlx + cost*sinp*dly - sint*dlz)*rm5 - cost*cosp*rm3
     By(1, 4) = 3.0_dp*dly*( cost*cosp*dlx + cost*sinp*dly - sint*dlz)*rm5 - cost*sinp*rm3 
     Bz(1, 4) = 3.0_dp*dlz*( cost*cosp*dlx + cost*sinp*dly - sint*dlz)*rm5 + sint     *rm3 

     Bx(1, 5) = 3.0_dp*dlx*(-sint*sinp*dlx + sint*cosp*dly           )*rm5 + sint*sinp*rm3
     By(1, 5) = 3.0_dp*dly*(-sint*sinp*dlx + sint*cosp*dly           )*rm5 - sint*cosp*rm3
     Bz(1, 5) = 3.0_dp*dlz*(-sint*sinp*dlx + sint*cosp*dly           )*rm5 
#else
     ! dipole origins are fixed
     Bx(1, 1) = 3.0_dp*dlx*( cost*cosp*dlx + cost*sinp*dly - sint*dlz)*rm5 - cost*cosp*rm3
     By(1, 1) = 3.0_dp*dly*( cost*cosp*dlx + cost*sinp*dly - sint*dlz)*rm5 - cost*sinp*rm3 
     Bz(1, 1) = 3.0_dp*dlz*( cost*cosp*dlx + cost*sinp*dly - sint*dlz)*rm5 + sint     *rm3 

     Bx(1, 2) = 3.0_dp*dlx*(-sint*sinp*dlx + sint*cosp*dly           )*rm5 + sint*sinp*rm3
     By(1, 2) = 3.0_dp*dly*(-sint*sinp*dlx + sint*cosp*dly           )*rm5 - sint*cosp*rm3
     Bz(1, 2) = 3.0_dp*dlz*(-sint*sinp*dlx + sint*cosp*dly           )*rm5 
#endif

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

subroutine coils_bfield(B,x,y,z)
  use globals, only: dp, coil, surf, Ncoils, Nteta, Nzeta, &
       zero, myid, ounit, Nfp, bsconstant, one, two, ncpu, &
       master, nworker, myworkid, MPI_COMM_MASTERS, MPI_COMM_MYWORLD, MPI_COMM_WORKERS
  use mpi
  implicit none

  REAL  , intent( in)     :: x, y, z
  REAL  , intent(inout)   :: B(3)
  !INTEGER, INTENT(in)     :: comm ! MPI communicator

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  INTEGER              :: ierr, astat
  REAL                 :: Bx, By, Bz
  INTEGER              :: icoil

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  call MPI_BARRIER(MPI_COMM_MYWORLD, ierr ) ! wait all cpus;

  B = zero
  do icoil = 1, Ncoils
     if ( myworkid /= modulo(icoil-1, nworker) ) cycle ! MPI
     ! Bx = zero; By = zero; Bz = zero
     call bfield0( icoil, x, y, z, Bx, By, Bz )
     B(1) = B(1) + Bx
     B(2) = B(2) + By
     B(3) = B(3) + Bz
  enddo

  call MPI_ALLREDUCE(MPI_IN_PLACE, B, 3, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_MYWORLD, ierr )

  return

end subroutine coils_bfield

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
