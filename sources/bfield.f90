
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
  use globals, only: dp, coil, surf, Ncoils, Nteta, Nzeta, cosnfp, sinnfp, Nfp_raw, MPI_COMM_FAMUS, &
                     zero, myid, ounit, Nfp, pi2, half, two, one, bsconstant, momentq, machprec
  implicit none
  include "mpif.h"

  INTEGER, intent(in ) :: icoil
  REAL,    intent(in ) :: x, y, z
  REAL   , intent(out) :: tBx, tBy, tBz

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  INTEGER              :: ierr, astat, kseg, ip, is, symmetry, npc ! local symmetry
  REAL                 :: dlx, dly, dlz, rm3, ltx, lty, ltz, rr, r2, m_dot_r, &
                        & mx, my, mz, xx, yy, zz, Bx, By, Bz, sBx, sBy, sBz

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  FATAL( bfield0, icoil .lt. 1 .or. icoil .gt. Ncoils, icoil not in right range )

  tBx = zero ; tBy = zero ; tBz = zero
  dlx = zero ; dly = zero ; dlz = zero
  ltx = zero ; lty = zero ; ltz = zero
  
  ! check if stellarator symmetric
  npc = Nfp_raw
  symmetry = 0 
  if (coil(icoil)%symmetry == 0) Npc = 1
  if (coil(icoil)%symmetry == 2) symmetry = 1
  
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
           ! if(r2<1E-4) write(ounit,*) myid, icoil, r2
           rm3 = one/(sqrt(r2)*r2)
           mx = coil(icoil)%mx * (-1)**is / coil(icoil)%I
           my = coil(icoil)%my / coil(icoil)%I
           mz = coil(icoil)%mz / coil(icoil)%I
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

subroutine fixed_bfield0(icoil, x, y, z, tBx, tBy, tBz)
   !------------------------------------------------------------------------------------------------------ 
   ! DATE:  06/15/2016; 03/26/2017
   ! calculate the magnetic field of icoil using manually discretized coils. 
   ! Biot-Savart constant and currents are not included for later simplication. 
   ! Be careful if coils have different resolutions.
   !------------------------------------------------------------------------------------------------------   
     use globals, only: dp, fixed_coil, surf, fixed_Ncoils, Nteta, Nzeta, cosnfp, sinnfp, Nfp_raw, MPI_COMM_FAMUS, &
                        zero, myid, ounit, Nfp, pi2, half, two, one, bsconstant, momentq, machprec
     implicit none
     include "mpif.h"
   
     INTEGER, intent(in ) :: icoil
     REAL,    intent(in ) :: x, y, z
     REAL   , intent(out) :: tBx, tBy, tBz
   
   !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
     INTEGER              :: ierr, astat, kseg, ip, is, symmetry, npc ! local symmetry
     REAL                 :: dlx, dly, dlz, rm3, ltx, lty, ltz, rr, r2, m_dot_r, &
                           & mx, my, mz, xx, yy, zz, Bx, By, Bz, sBx, sBy, sBz
   
   !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
     FATAL( bfield0, icoil .lt. 1 .or. icoil .gt. fixed_Ncoils, icoil not in right range )
   
     tBx = zero ; tBy = zero ; tBz = zero
     dlx = zero ; dly = zero ; dlz = zero
     ltx = zero ; lty = zero ; ltz = zero
     
     ! check if stellarator symmetric
     npc = Nfp_raw
     symmetry = 0 
     if (fixed_coil(icoil)%symmetry == 0) Npc = 1
     if (fixed_coil(icoil)%symmetry == 2) symmetry = 1
     
     do ip = 1, Npc
        xx =  x*cosnfp(ip) + y*sinnfp(ip) ! find the point on plasma
        yy = -x*sinnfp(ip) + y*cosnfp(ip) ! by rotating in reverse direction.
        zz =  z
        sBx = zero; sBy = zero; sBz = zero
        do is = 0, symmetry         ! stellarator symmetry
           Bx = zero; By = zero; Bz = zero
           ! calculate magnetic field at (xx,yy,zz)
           select case (fixed_coil(icoil)%itype)     
           case(1) ! Fourier coil
              do kseg = 0, fixed_coil(icoil)%NS-1
                 ! position vector
                 dlx = xx - fixed_coil(icoil)%xx(kseg)
                 dly = yy - fixed_coil(icoil)%yy(kseg)
                 dlz = zz - fixed_coil(icoil)%zz(kseg)
                 rm3 = (sqrt(dlx**2 + dly**2 + dlz**2))**(-3)
                 ! tangential coil
                 ltx = fixed_coil(icoil)%xt(kseg)
                 lty = fixed_coil(icoil)%yt(kseg)
                 ltz = fixed_coil(icoil)%zt(kseg)
                 ! Biot-Savart Law
                 Bx = Bx + ( dlz*lty - dly*ltz ) * rm3 * fixed_coil(icoil)%dd(kseg)
                 By = By + ( dlx*ltz - dlz*ltx ) * rm3 * fixed_coil(icoil)%dd(kseg)
                 Bz = Bz + ( dly*ltx - dlx*lty ) * rm3 * fixed_coil(icoil)%dd(kseg)
              enddo    ! enddo kseg
           case(2) ! magnetic dipole
              dlx = xx -            fixed_coil(icoil)%ox
              dly = yy - (-1)**is * fixed_coil(icoil)%oy
              dlz = zz - (-1)**is * fixed_coil(icoil)%oz
              r2  = dlx**2 + dly**2 + dlz**2
              ! if(r2<1E-4) write(ounit,*) myid, icoil, r2
              rm3 = one/(sqrt(r2)*r2)
              mx = fixed_coil(icoil)%mx * (-1)**is / fixed_coil(icoil)%I
              my = fixed_coil(icoil)%my / fixed_coil(icoil)%I
              mz = fixed_coil(icoil)%mz / fixed_coil(icoil)%I
              m_dot_r = mx * dlx + my * dly + mz * dlz
              ! Magnetic field
              Bx = 3.0_dp * m_dot_r * rm3 / r2 * dlx - mx * rm3 
              By = 3.0_dp * m_dot_r * rm3 / r2 * dly - my * rm3
              Bz = 3.0_dp * m_dot_r * rm3 / r2 * dlz - mz * rm3
   
           case(3)
              ! might be only valid for cylindrical coordinates
              ! Bt = u0*I/(2 pi R)
              rr = sqrt( xx**2 + yy**2 )
              fixed_coil(icoil)%Bt = two/rr !* coil(icoil)%I * bsconstant
   
              Bx = - fixed_coil(icoil)%Bt * yy/rr
              By =   fixed_coil(icoil)%Bt * xx/rr
              Bz =   fixed_coil(icoil)%Bz / (fixed_coil(icoil)%I * bsconstant + machprec) ! add machprec to avoid divide by zero
           case default
              FATAL(bfield0, .true., not supported coil types)
           end select
           sBx = sBx + Bx
           sBy = sBy + By
           sBz = sBz + Bz
        end do ! end do is
        
        ! exit for those have no symmetries
        if (fixed_coil(icoil)%symmetry == 0) then
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
   
   end subroutine fixed_bfield0   

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine coils_bfield(B,x,y,z)
  use globals, only: dp, coil, surf, Ncoils, Nteta, Nzeta, &
       zero, myid, ounit, bsconstant, one, two, ncpu, cosnfp, sinnfp, &
       master, nworker, myworkid, MPI_COMM_MASTERS, MPI_COMM_MYWORLD, MPI_COMM_WORKERS, MPI_COMM_FAMUS
  use mpi
  implicit none

  REAL  , intent( in)     :: x, y, z
  REAL  , intent(inout)   :: B(3)
  !INTEGER, INTENT(in)     :: comm ! MPI communicator

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  INTEGER              :: ierr, astat, ip
  REAL                 :: Bx, By, Bz
  INTEGER              :: icoil

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  call MPI_BARRIER(MPI_COMM_FAMUS, ierr ) ! wait all cpus;
  B = zero
  do icoil = 1, Ncoils
     ! if ( myworkid /= modulo(icoil-1, nworker) ) cycle ! MPI
     ! Bx = zero; By = zero; Bz = zero
     call bfield0( icoil, x, y, z, Bx, By, Bz )
     B(1) = B(1) + Bx*coil(icoil)%I*bsconstant
     B(2) = B(2) + By*coil(icoil)%I*bsconstant
     B(3) = B(3) + Bz*coil(icoil)%I*bsconstant
  enddo
  call MPI_ALLREDUCE(MPI_IN_PLACE, B, 3, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FAMUS, ierr )

  return

end subroutine coils_bfield

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
