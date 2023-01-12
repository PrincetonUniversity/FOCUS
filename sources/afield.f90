
!title (afield) ! Computes magnetic vector potential field.

!latex \briefly{Computes magnetic vector potential field given coil geometry.}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine afield0(icoil, x, y, z, tAx, tAy, tAz)

  use globals, only: dp, coil, surf, Ncoils, Nteta, Nzeta, cosnfp, sinnfp, &
                     zero, myid, ounit, Nfp, pi2, half, two, one, bsconstant, &
                     coil_type_multi, coil_type_spline, MPI_COMM_FOCUS
  use mpi
  implicit none

  INTEGER, intent(in)  :: icoil
  REAL,    intent(in)  :: x, y, z
  REAL   , intent(out) :: tAx, tAy, tAz

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  INTEGER              :: ierr, astat, kseg, ip, is, cs, Npc
  REAL                 :: dlx, dly, dlz, rm1, rr, r2, m_dot_r, &
       &                  mx, my, mz, xx, yy, zz, Ax, Ay, Az, absrp

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  FATAL( afield0, icoil .lt. 1 .or. icoil .gt. Ncoils, icoil not in right range )
  ! initialization
  tAx = zero ; tAy = zero ; tAz = zero
  ! check if the coil is stellarator symmetric
  select case (coil(icoil)%symm) 
  case ( 0 )
     cs  = 0
     Npc = 1
  case ( 1 )
     cs  = 0
     Npc = Nfp
  case ( 2 ) 
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
        Ax = zero; Ay = zero; Az = zero
        select case (coil(icoil)%type)
        ! Fourier coils
        case(1,coil_type_multi,coil_type_spline)
           ! Biot-Savart law
           do kseg = 0, coil(icoil)%NS-1
              dlx = xx - coil(icoil)%xx(kseg)                                   ! r-r_c
              dly = yy - coil(icoil)%yy(kseg)
              dlz = zz - coil(icoil)%zz(kseg)
              rm1 = (dlx**2 + dly**2 + dlz**2)**(-0.5)                      ! |r-r_c|^-1
              Ax = Ax + coil(icoil)%xt(kseg)*rm1
              Ay = Ay + coil(icoil)%yt(kseg)*rm1
              Az = Az + coil(icoil)%zt(kseg)*rm1
           enddo    ! enddo kseg
           Ax = Ax * coil(icoil)%I * bsconstant
           Ay = Ay * coil(icoil)%I * bsconstant
           Az = Az * coil(icoil)%I * bsconstant
        ! magnetic dipoles
        case(2)
           FATAL(afield0, .true., not supported coil types)
        ! toroidal field and verticle field
        case(3)
           FATAL(afield0, .true., not supported coil types)
        case default
           FATAL(afield0, .true., not supported coil types)
        end select     
        ! sum all the contributions
        Ax = Ax*(-1)**is
        tAx = tAx + (Ax*cosnfp(ip) - Ay*sinnfp(ip))
        tAy = tAy + (Ay*cosnfp(ip) + Ax*sinnfp(ip))
        tAz = tAz +  Az
     enddo
  enddo
  tAx = pi2*tAx/coil(icoil)%NS
  tAy = pi2*tAy/coil(icoil)%NS
  tAz = pi2*tAz/coil(icoil)%NS

  return

end subroutine afield0

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine afield1(icoil, x, y, z, tAx, tAy, tAz, ND)

  use globals, only: dp, coil, DoF, surf, NFcoil, Ncoils, Nteta, Nzeta, &
                     zero, myid, ounit, Nfp, one, bsconstant, cosnfp, sinnfp, &
                     pi2, coil_type_multi, coil_type_spline, MPI_COMM_FOCUS
  use mpi
  implicit none

  INTEGER, intent(in) :: icoil, ND
  REAL,    intent(in) :: x, y, z
  REAL, dimension(1:1, 1:ND), intent(inout) :: tAx, tAy, tAz

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  INTEGER              :: ierr, astat, kseg, NS, ip, is, cs, Npc
  REAL                 :: dlx, dly, dlz, r2, rm3, m_dot_r, dots, &
                          sinp, sint, cosp, cost, mx, my, mz, xx, yy, zz
  REAL, dimension(1:1, 1:ND) :: Ax, Ay, Az
  REAL, dimension(1:1, 0:coil(icoil)%NS-1)   :: dAxx, dAxy, dAxz, dAyx, dAyy, dAyz, dAzx, dAzy, dAzz

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  FATAL( afield1, icoil .lt. 1 .or. icoil .gt. Ncoils, icoil not in right range )
  FATAL( afield1, ND <= 0, wrong inout dimension of ND )

  tAx = zero ; tAy = zero ; tAz = zero

  ! check if the coil is stellarator symmetric
  select case (coil(icoil)%symm) 
  case ( 0 )
     cs  = 0
     Npc = 1
  case ( 1 )
     cs  = 0
     Npc = Nfp
  case ( 2 ) 
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
        Ax = zero; Ay = zero; Az = zero
        select case (coil(icoil)%type)
        ! Fourier coils
        case(1,coil_type_multi,coil_type_spline)
           do kseg = 0, coil(icoil)%NS-1
              dlx = xx - coil(icoil)%xx(kseg)          ! r-r_c
              dly = yy - coil(icoil)%yy(kseg)
              dlz = zz - coil(icoil)%zz(kseg)
              r2 = dlx**2 + dly**2 + dlz**2 
              rm3 = one/(sqrt(r2)*r2)                  ! |r-r_c|^-3
              dots = dlx*coil(icoil)%xt(kseg) + dly*coil(icoil)%yt(kseg) + dlz*coil(icoil)%zt(kseg)  ! (r-r_c) dot r_c'
              dAxx(1,kseg) = rm3*(dlx*coil(icoil)%xt(kseg)-dots)
              dAyx(1,kseg) = rm3* dly*coil(icoil)%xt(kseg)
              dAzx(1,kseg) = rm3* dlz*coil(icoil)%xt(kseg)
              dAxy(1,kseg) = rm3* dlx*coil(icoil)%yt(kseg)
              dAyy(1,kseg) = rm3*(dly*coil(icoil)%yt(kseg)-dots)
              dAzy(1,kseg) = rm3* dlz*coil(icoil)%yt(kseg) 
              dAxz(1,kseg) = rm3* dlx*coil(icoil)%zt(kseg) 
              dAyz(1,kseg) = rm3* dly*coil(icoil)%zt(kseg)
              dAzz(1,kseg) = rm3*(dlz*coil(icoil)%zt(kseg)-dots)
           enddo
           Ax(1:1, 1:ND) = matmul(dAxx, DoF(icoil)%xof) + matmul(dAyx, DoF(icoil)%yof) + matmul(dAzx, DoF(icoil)%zof)
           Ay(1:1, 1:ND) = matmul(dAxy, DoF(icoil)%xof) + matmul(dAyy, DoF(icoil)%yof) + matmul(dAzy, DoF(icoil)%zof)
           Az(1:1, 1:ND) = matmul(dAxz, DoF(icoil)%xof) + matmul(dAyz, DoF(icoil)%yof) + matmul(dAzz, DoF(icoil)%zof)
           Ax = Ax * coil(icoil)%I * bsconstant
           Ay = Ay * coil(icoil)%I * bsconstant
           Az = Az * coil(icoil)%I * bsconstant
        ! Magnetic dipoles 
        case(2)
           FATAL(afield1, .true., not supported coil types)
        ! toroidal field and verticle field
        case(3)
           FATAL(afield1, .true., not supported coil types)
        case default
           FATAL(afield1, .true., not supported coil types)
        end select
        ! sum all the contributions
        Ax = Ax*(-1)**is
        tAx = tAx + (Ax*cosnfp(ip) - Ay*sinnfp(ip))
        tAy = tAy + (Ay*cosnfp(ip) + Ax*sinnfp(ip))
        tAz = tAz +  Az
     enddo
  enddo
  tAx = pi2*tAx/coil(icoil)%NS
  tAy = pi2*tAy/coil(icoil)%NS
  tAz = pi2*tAz/coil(icoil)%NS

  return

end subroutine afield1

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine deltaafield(icoil, x, y, z, dAxx, dAxy, dAxz, dAyx, dAyy, dAyz, dAzx, dAzy, dAzz)

  use globals, only: dp, coil, DoF, surf, NFcoil, Ncoils, Nteta, Nzeta, &
                     zero, myid, ounit, Nfp, one, bsconstant, cosnfp, sinnfp, coil_type_multi, &
                     coil_type_spline, MPI_COMM_FOCUS
  use mpi
  implicit none

  INTEGER, intent(in) :: icoil
  REAL,    intent(in) :: x, y, z
  REAL, dimension(1:coil(icoil)%NS), intent(out) :: dAxx, dAxy, dAxz, dAyx, dAyy, dAyz, dAzx, dAzy, dAzz

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  INTEGER              :: ierr, astat, kseg, NS, ip, is, cs, Npc
  REAL                 :: dlx, dly, dlz, r2, rm3, m_dot_r, ltx, lty, ltz, dots, &
                          sinp, sint, cosp, cost, mx, my, mz, xx, yy, zz

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  FATAL( deltaafield, icoil .lt. 1 .or. icoil .gt. Ncoils, icoil not in right range )

  dlx = zero ; dly = zero ; dlz = zero
  ltx = zero ; lty = zero ; ltz = zero

  xx = x
  yy = y
  zz = z

  select case (coil(icoil)%type)
  ! Fourier coils
  case(1,coil_type_multi,coil_type_spline)
     NS = coil(icoil)%NS
     ! Should be vectorized
     do kseg = 1, NS
        dlx = xx - coil(icoil)%xx(kseg)                                              ! r-r_c
        dly = yy - coil(icoil)%yy(kseg)
        dlz = zz - coil(icoil)%zz(kseg)
        r2 = dlx**2 + dly**2 + dlz**2                                                ! !r_c-r|^2
        rm3 = one/(sqrt(r2)*r2)                                                      ! |r_c-r|^-1
        ltx = coil(icoil)%xt(kseg)                                                   ! r_c'
        lty = coil(icoil)%yt(kseg)
        ltz = coil(icoil)%zt(kseg)
        dots = dlx*ltx + dly*lty + dlz*ltz                                           ! (r-r_c) dot r_c'

        dAxx(kseg) = coil(icoil)%I*bsconstant*( rm3*dlx*ltx - rm3*dots ) !Ax/x
        dAxy(kseg) = coil(icoil)%I*bsconstant*( rm3*dly*ltx )            !Ax/y
        dAxz(kseg) = coil(icoil)%I*bsconstant*( rm3*dlz*ltx )            !Ax/z
        dAyx(kseg) = coil(icoil)%I*bsconstant*( rm3*dlx*lty )            !Ay/x
        dAyy(kseg) = coil(icoil)%I*bsconstant*( rm3*dly*lty - rm3*dots ) !Ay/y
        dAyz(kseg) = coil(icoil)%I*bsconstant*( rm3*dlz*lty )            !Ay/z
        dAzx(kseg) = coil(icoil)%I*bsconstant*( rm3*dlx*ltz )            !Az/x
        dAzy(kseg) = coil(icoil)%I*bsconstant*( rm3*dly*ltz )            !Az/y
        dAzz(kseg) = coil(icoil)%I*bsconstant*( rm3*dlz*ltz - rm3*dots ) !Az/z
     enddo
  ! Magnetic dipoles 
  case(2)
     FATAL(deltaafield, .true., not supported coil types)
  ! toroidal field and verticle field
  case(3)
     FATAL(deltaafield, .true., not supported coil types)
  case default
     FATAL(deltaafield, .true., not supported coil types)
  end select

  return

end subroutine deltaafield
