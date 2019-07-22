
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
  use globals, only: dp, coil, surf, Ncoils, Nteta, Nzeta, cosnfp, sinnfp, Npc, &
                     zero, myid, ounit, Npc, Nfp, pi2, half, two, one, bsconstant, momentq, machprec
  implicit none
  include "mpif.h"

  INTEGER, intent(in ) :: icoil
  REAL,    intent(in ) :: x, y, z
  REAL   , intent(out) :: tBx, tBy, tBz

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  INTEGER              :: ierr, astat, kseg, ip, is, symmetry
  REAL                 :: dlx, dly, dlz, rm3, ltx, lty, ltz, rr, r2, m_dot_r, &
                        & mx, my, mz, xx, yy, zz, Bx, By, Bz, sBx, sBy, sBz

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  FATAL( bfield0, icoil .lt. 1 .or. icoil .gt. Ncoils, icoil not in right range )

  tBx = zero ; tBy = zero ; tBz = zero
  dlx = zero ; dly = zero ; dlz = zero
  ltx = zero ; lty = zero ; ltz = zero

  ! check if stellarator symmetric
  if (coil(icoil)%symmetry == 2) then
     symmetry = 1
  else
     symmetry = 0
  endif
  
  do ip = 1, Npc
     xx =  x*cosnfp(ip) + y*sinnfp(ip) ! find the point on plasma
     yy = -x*sinnfp(ip) + y*cosnfp(ip) ! by rotating in reverse direction.
     zz =  z
     sBx = zero; sBy = zero; sBz = zero
     do is = 0, symmetry         ! stellarator symmetry
        Bx = zero; By = zero; Bz = zero
        ! calculate magnetic field at (xx,yy,zz)
        select case (coil(icoil)%itype)     
        case(1) ! Fourier coil
           do kseg = 0, coil(icoil)%NS-1
              ! position vector
              dlx = xx - coil(icoil)%xx(kseg)
              dly = yy - coil(icoil)%yy(kseg)
              dlz = zz - coil(icoil)%zz(kseg)
              rm3 = (sqrt(dlx**2 + dly**2 + dlz**2))**(-3)
              ! tangential coil
              ltx = coil(icoil)%xt(kseg)
              lty = coil(icoil)%yt(kseg)
              ltz = coil(icoil)%zt(kseg)
              ! Biot-Savart Law
              Bx = Bx + ( dlz*lty - dly*ltz ) * rm3 * coil(icoil)%dd(kseg)
              By = By + ( dlx*ltz - dlz*ltx ) * rm3 * coil(icoil)%dd(kseg)
              Bz = Bz + ( dly*ltx - dlx*lty ) * rm3 * coil(icoil)%dd(kseg)
           enddo    ! enddo kseg
           ! normalize currents and mu0/4pi
!!$        Bx = Bx * coil(icoil)%I * bsconstant
!!$        By = By * coil(icoil)%I * bsconstant
!!$        Bz = Bz * coil(icoil)%I * bsconstant
        case(2) ! magnetic dipole
           dlx = xx -            coil(icoil)%ox
           dly = yy - (-1)**is * coil(icoil)%oy
           dlz = zz - (-1)**is * coil(icoil)%oz
           r2  = dlx**2 + dly**2 + dlz**2
           rm3 = one/(sqrt(r2)*r2)
           mx = sin(coil(icoil)%mt) * cos(coil(icoil)%mp) * (-1)**is
           my = sin(coil(icoil)%mt) * sin(coil(icoil)%mp)
           mz = cos(coil(icoil)%mt)
           m_dot_r = mx * dlx + my * dly + mz * dlz
           ! Magnetic field
           Bx = 3.0_dp * m_dot_r * rm3 / r2 * dlx - mx * rm3 
           By = 3.0_dp * m_dot_r * rm3 / r2 * dly - my * rm3
           Bz = 3.0_dp * m_dot_r * rm3 / r2 * dlz - mz * rm3

           !coil(icoil)%I = coil(icoil)%moment*sin(coil(icoil)%pho)**momentq
!!$        Bx = Bx * coil(icoil)%I * bsconstant
!!$        By = By * coil(icoil)%I * bsconstant
!!$        Bz = Bz * coil(icoil)%I * bsconstant
        case(3)
           ! might be only valid for cylindrical coordinates
           ! Bt = u0*I/(2 pi R)
           rr = sqrt( xx**2 + yy**2 )
           coil(icoil)%Bt = two/rr !* coil(icoil)%I * bsconstant

           Bx = - coil(icoil)%Bt * yy/rr
           By =   coil(icoil)%Bt * xx/rr
           Bz =   coil(icoil)%Bz / (coil(icoil)%I * bsconstant + machprec) ! add machprec to avoid divide by zero
        case default
           FATAL(bfield0, .true., not supported coil types)
        end select
        sBx = sBx + Bx
        sBy = sBy + By
        sBz = sBz + Bz
     end do ! end do is
     
     ! exit for those have no symmetries
     if (coil(icoil)%symmetry == 0) then
        tBx = sBx
        tBy = sBy
        tBz = sBz
        exit
     endif
     ! otherwise, sum the symmetric points
     tBx = tBx + sBx*cosnfp(ip) - sBy*sinnfp(ip)
     tBy = tBy + sBy*cosnfp(ip) + sBx*sinnfp(ip)
     tBz = tBz + sBz
  enddo ! end do ip

  return

end subroutine bfield0

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine bfield1(icoil, x, y, z, tBx, tBy, tBz, ND)
!------------------------------------------------------------------------------------------------------ 
! DATE:  06/15/2016; 03/26/2017
! calculate the magnetic field and the first dirivatives of icoil using manually discretized coils;
! Biot-Savart constant and currents are not included for later simplication;
! Discretizing factor is includeed; coil(icoil)%dd(kseg)
!------------------------------------------------------------------------------------------------------   
  use globals, only: dp, coil, DoF, surf, NFcoil, Ncoils, Nteta, Nzeta, &
                     zero, myid, ounit, Npc, one, bsconstant, cosnfp, sinnfp, Npc, momentq
  implicit none
  include "mpif.h"

  INTEGER, intent(in ) :: icoil, ND
  REAL,    intent(in ) :: x, y, z
  REAL, dimension(1:1, 1:ND), intent(inout) :: tBx, tBy, tBz

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  INTEGER              :: ierr, astat, kseg, NS, ip, is, symmetry
  REAL                 :: dlx, dly, dlz, r2, rm3, rm5, rm7, m_dot_r, ltx, lty, ltz, rxp, &
                          sinp, sint, cosp, cost, mx, my, mz, xx, yy, zz
  REAL, dimension(1:1, 1:ND) :: Bx, By, Bz, sBx, sBy, sBz
  REAL, dimension(1:1, 0:coil(icoil)%NS-1) :: dBxx, dBxy, dBxz, dByx, dByy, dByz, dBzx, dBzy, dBzz

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  FATAL( bfield1, icoil .lt. 1 .or. icoil .gt. Ncoils, icoil not in right range )
  FATAL( bfield1, ND <= 0, wrong inout dimension of ND )

  ! check if stellarator symmetric
  if (coil(icoil)%symmetry == 2) then
     symmetry = 1
  else
     symmetry = 0
  endif

  tBx = zero ; tBy = zero ; tBz = zero
  dlx = zero ; dly = zero ; dlz = zero
  ltx = zero ; lty = zero ; ltz = zero

  do ip = 1, Npc
     xx =  x*cosnfp(ip) + y*sinnfp(ip) ! find the point on plasma
     yy = -x*sinnfp(ip) + y*cosnfp(ip) ! by rotating in reverse direction.
     zz = z
     sBx = zero; sBy = zero; sBz = zero
     do is = 0, symmetry         ! stellarator symmetry
        Bx = zero ; By = zero ; Bz = zero
        select case (coil(icoil)%itype)

        case(1) ! Fourier coils
           NS = coil(icoil)%NS
           do kseg = 0, NS-1

              dlx = xx - coil(icoil)%xx(kseg)
              dly = yy - coil(icoil)%yy(kseg)
              dlz = zz - coil(icoil)%zz(kseg)

              r2 = dlx**2 + dly**2 + dlz**2
              rm3 = one/(sqrt(r2)*r2)
              rm5 = rm3/r2;

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
           dly = yy - coil(icoil)%oy * (-1)**is
           dlz = zz - coil(icoil)%oz * (-1)**is
           r2  = dlx**2 + dly**2 + dlz**2
           rm3 = one/(sqrt(r2)*r2)
           rm5 = rm3/r2
           rm7 = rm5/r2

           cost = cos(coil(icoil)%mt) ; sint = sin(coil(icoil)%mt)
           cosp = cos(coil(icoil)%mp) ; sinp = sin(coil(icoil)%mp)    
           mx = sint*cosp * (-1)**is 
           my = sint*sinp
           mz = cost
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

           !coil(icoil)%I = coil(icoil)%moment*sin(coil(icoil)%pho)**momentq

           Bx = Bx * coil(icoil)%I * bsconstant
           By = By * coil(icoil)%I * bsconstant
           Bz = Bz * coil(icoil)%I * bsconstant

           !--------------------------------------------------------------------------------------------- 
        case(3)  ! only for Bz

           Bx = zero
           By = zero
           Bz = one

        end select
        sBx = sBx + Bx
        sBy = sBy + By
        sBz = sBz + Bz
     end do ! end do is
     
     ! exit for those have no symmetries
     if (coil(icoil)%symmetry == 0) then
        tBx = Bx
        tBy = By
        tBz = Bz
        exit
     endif
     ! otherwise, sum the symmetric points
     tBx = tBx + sBx*cosnfp(ip) - By*sinnfp(ip)
     tBy = tBy + sBy*cosnfp(ip) + Bx*sinnfp(ip)
     tBz = tBz + sBz
  enddo ! end do ip

  return

end subroutine bfield1

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine coils_bfield(B,x,y,z)
  use globals, only: dp, coil, surf, Ncoils, Nteta, Nzeta, &
       zero, myid, ounit, Npc, bsconstant, one, two, ncpu, Npc, cosnfp, sinnfp, &
       master, nworker, myworkid, MPI_COMM_MASTERS, MPI_COMM_MYWORLD, MPI_COMM_WORKERS
  use mpi
  implicit none

  REAL  , intent( in)     :: x, y, z
  REAL  , intent(inout)   :: B(3)
  !INTEGER, INTENT(in)     :: comm ! MPI communicator

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  INTEGER              :: ierr, astat, ip
  REAL                 :: Bx, By, Bz, tBx, tBy, tBz
  INTEGER              :: icoil

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  call MPI_BARRIER(MPI_COMM_WORLD, ierr ) ! wait all cpus;
  B = zero
  do icoil = 1, Ncoils
     ! if ( myworkid /= modulo(icoil-1, nworker) ) cycle ! MPI
     ! Bx = zero; By = zero; Bz = zero
     call bfield0( icoil, x, y, z, Bx, By, Bz )
     B(1) = B(1) + Bx*coil(icoil)%I*bsconstant
     B(2) = B(2) + By*coil(icoil)%I*bsconstant
     B(3) = B(3) + Bz*coil(icoil)%I*bsconstant
  enddo
!!$  do ip = 1, Npc
!!$     tBx = zero ; tBy = zero ; tBz = zero
!!$     do icoil = 1, Ncoils
!!$        ! if ( myworkid /= modulo(icoil-1, nworker) ) cycle ! MPI
!!$        ! Bx = zero; By = zero; Bz = zero
!!$        call bfield0( icoil, x, y, z, Bx, By, Bz )
!!$        tBx = tBx + Bx
!!$        tBy = tBy + By
!!$        tBz = tBz + Bz        
!!$     enddo
!!$     B(1) = B(1) + tBx*cosnfp(ip) - tBy*sinnfp(ip)
!!$     B(2) = B(2) + tBy*cosnfp(ip) + tBx*sinnfp(ip)
!!$     B(3) = B(3) + tBz
!!$  enddo
  call MPI_ALLREDUCE(MPI_IN_PLACE, B, 3, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )

  return

end subroutine coils_bfield

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
