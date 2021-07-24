
!title (afield) ! Computes magnetic vector potential field.

!latex \briefly{Computes magnetic vector potential field given coil geometry.}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine afield0(icoil, x, y, z, tAx, tAy, tAz)

  use globals, only: dp, coil, surf, Ncoils, Nteta, Nzeta, cosnfp, sinnfp, &
                     zero, myid, ounit, Nfp, pi2, half, two, one, bsconstant, MPI_COMM_FOCUS
  use mpi
  implicit none

  INTEGER, intent(in ) :: icoil
  REAL,    intent(in ) :: x, y, z
  REAL   , intent(out) :: tAx, tAy, tAz

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  INTEGER              :: ierr, astat, kseg, ip, is, cs, Npc
  REAL                 :: dlx, dly, dlz, rm1, ltx, lty, ltz, rr, r2, m_dot_r, &
       &                  mx, my, mz, xx, yy, zz, Ax, Ay, Az, absrp

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  FATAL( afield0, icoil .lt. 1 .or. icoil .gt. Ncoils, icoil not in right range )
  ! initialization
  Npc = 1 ; cs = 0 ; ip = 1
  tAx = zero ; tAy = zero ; tAz = zero
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
        Ax = zero; Ay = zero; Az = zero
        select case (coil(icoil)%type)
        ! Fourier coils
        case(1)
           ! Biot-Savart law
           do kseg = 0, coil(icoil)%NS-1
              dlx = xx - coil(icoil)%xx(kseg)                                   ! r-r_c
              dly = yy - coil(icoil)%yy(kseg)
              dlz = zz - coil(icoil)%zz(kseg)
              rm1 = (sqrt(dlx**2 + dly**2 + dlz**2))**(-1)                      ! |r_c-r|^-1
              ltx = coil(icoil)%xt(kseg)                                        ! r_c'
              lty = coil(icoil)%yt(kseg)
              ltz = coil(icoil)%zt(kseg)
              Ax = Ax + ltx*rm1
              Ay = Ay + lty*rm1
              Az = Az + ltz*rm1
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

  return

end subroutine afield0

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine afield1(icoil, x, y, z, tAx, tAy, tAz, ND)

  use globals, only: dp, coil, DoF, surf, NFcoil, Ncoils, Nteta, Nzeta, &
                     zero, myid, ounit, Nfp, one, bsconstant, cosnfp, sinnfp, MPI_COMM_FOCUS
  use mpi
  implicit none

  INTEGER, intent(in ) :: icoil, ND
  REAL,    intent(in ) :: x, y, z
  REAL, dimension(1:1, 1:ND), intent(inout) :: tAx, tAy, tAz

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  INTEGER              :: ierr, astat, kseg, NS, ip, is, cs, Npc
  REAL                 :: dlx, dly, dlz, r2, rm3, m_dot_r, ltx, lty, ltz, dots, &
                          sinp, sint, cosp, cost, mx, my, mz, xx, yy, zz
  REAL, dimension(1:1, 1:ND) :: Ax, Ay, Az
  REAL, dimension(1:1, 0:coil(icoil)%NS-1)   :: dAxx, dAxy, dAxz, dAyx, dAyy, dAyz, dAzx, dAzy, dAzz

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  FATAL( afield1, icoil .lt. 1 .or. icoil .gt. Ncoils, icoil not in right range )
  FATAL( afield1, ND <= 0, wrong inout dimension of ND )

  ! initialization
  Npc = 1 ; cs = 0 ; ip = 1
  tAx = zero ; tAy = zero ; tAz = zero
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
        Ax = zero; Ay = zero; Az = zero

        select case (coil(icoil)%type)
        ! Fourier coils
        case(1)
           NS = coil(icoil)%NS
           do kseg = 0, NS-1
              dlx = xx - coil(icoil)%xx(kseg)                                              ! r-r_c
              dly = yy - coil(icoil)%yy(kseg)
              dlz = zz - coil(icoil)%zz(kseg)
              r2 = dlx**2 + dly**2 + dlz**2 
              rm3 = one/(sqrt(r2)*r2)                                                      ! |r_c-r|^-1
              !rm5 = rm3/r2
              ltx = coil(icoil)%xt(kseg)                                                   ! r_c'
              lty = coil(icoil)%yt(kseg)
              ltz = coil(icoil)%zt(kseg)
              dots = -1.0*dlx*ltx + -1.0*dly*lty + -1.0*dlz*ltz                            ! (r_c-r) dot r_c'

              dAxx(1,kseg) = ( rm3*dlx*ltx + rm3*dots ) !Ax/x
              dAxy(1,kseg) = ( rm3*dly*ltx )            !Ax/y
              dAxz(1,kseg) = ( rm3*dlz*ltx )            !Ax/z
              dAyx(1,kseg) = ( rm3*dlx*lty )            !Ay/x
              dAyy(1,kseg) = ( rm3*dly*lty + rm3*dots ) !Ay/y
              dAyz(1,kseg) = ( rm3*dlz*lty )            !Ay/z
              dAzx(1,kseg) = ( rm3*dlx*ltz )            !Az/x
              dAzy(1,kseg) = ( rm3*dly*ltz )            !Az/y
              dAzz(1,kseg) = ( rm3*dlz*ltz + rm3*dots ) !Az/z
           enddo    ! enddo kseg
           ! db/dv = dB/dx * dx/dv  v->variables
           Ax(1:1, 1:ND) = matmul(dAxx, DoF(icoil)%xof) + matmul(dAxy, DoF(icoil)%yof) + matmul(dAxz, DoF(icoil)%zof)
           Ay(1:1, 1:ND) = matmul(dAyx, DoF(icoil)%xof) + matmul(dAyy, DoF(icoil)%yof) + matmul(dAyz, DoF(icoil)%zof)
           Az(1:1, 1:ND) = matmul(dAzx, DoF(icoil)%xof) + matmul(dAzy, DoF(icoil)%yof) + matmul(dAzz, DoF(icoil)%zof)
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

  return

end subroutine afield1

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine deltaafield(icoil, x, y, z, dAxx, dAxy, dAxz, dAyx, dAyy, dAyz, dAzx, dAzy, dAzz)

  use globals, only: dp, coil, DoF, surf, NFcoil, Ncoils, Nteta, Nzeta, &
                     zero, myid, ounit, Nfp, one, bsconstant, cosnfp, sinnfp, MPI_COMM_FOCUS
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
  case(1)
     NS = coil(icoil)%NS
     ! Should be vectorized
     do kseg = 1, NS
        dlx = xx - coil(icoil)%xx(kseg)                                              ! r-r_c
        dly = yy - coil(icoil)%yy(kseg)
        dlz = zz - coil(icoil)%zz(kseg)
        r2 = dlx**2 + dly**2 + dlz**2 
        rm3 = one/(sqrt(r2)*r2)                                                      ! |r_c-r|^-1
        ltx = coil(icoil)%xt(kseg)                                                   ! r_c'
        lty = coil(icoil)%yt(kseg)
        ltz = coil(icoil)%zt(kseg)
        dots = -1.0*dlx*ltx + -1.0*dly*lty + -1.0*dlz*ltz                            ! (r_c-r) dot r_c'

        dAxx(kseg) = coil(icoil)%I*bsconstant*( rm3*dlx*ltx + rm3*dots ) !Ax/x
        dAxy(kseg) = coil(icoil)%I*bsconstant*( rm3*dly*ltx )            !Ax/y
        dAxz(kseg) = coil(icoil)%I*bsconstant*( rm3*dlz*ltx )            !Ax/z
        dAyx(kseg) = coil(icoil)%I*bsconstant*( rm3*dlx*lty )            !Ay/x
        dAyy(kseg) = coil(icoil)%I*bsconstant*( rm3*dly*lty + rm3*dots ) !Ay/y
        dAyz(kseg) = coil(icoil)%I*bsconstant*( rm3*dlz*lty )            !Ay/z
        dAzx(kseg) = coil(icoil)%I*bsconstant*( rm3*dlx*ltz )            !Az/x
        dAzy(kseg) = coil(icoil)%I*bsconstant*( rm3*dly*ltz )            !Az/y
        dAzz(kseg) = coil(icoil)%I*bsconstant*( rm3*dlz*ltz + rm3*dots ) !Az/z
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
