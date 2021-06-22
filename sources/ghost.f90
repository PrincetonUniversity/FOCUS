!title (ghost) ! Calculate ghost/quadratic flux minimizing surfaces. (tkruger)

!latex \briefly{This function calculates a ghost surface to be used in the 
!latex         resonant Fourier harmonic calculation. It is important that we use
!latex         these ghost surfaces instead of a rational surface from VMEC. Rational
!latex         surfaces from VMEC will inherently have some error associated with 
!latex         their construction and will cause errors in the resonant Fouer 
!latex         harmonic calculation and moreover its sensitivity.}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! not parallelized; communications may take more time;
subroutine ghost(index)
  use globals, only: dp, zero, half, pi2, machprec, ncpu, myid, ounit, &
       coil, DoF, Ncoils, Nfixgeo, Ndof, resbn_m, gsurf, MPI_COMM_FOCUS
  
  use mpi
  implicit none
  INTEGER, INTENT(in) :: index

  INTEGER             :: astat, ierr

  ! Calculate f and g functions for o and x
  call calcfg(index)

  ! Calculate derivatives of f and g
  call calcfg_deriv(index)

  ! Optimize F = of^2 + xf^2 + og^2 + xg^2 functional
  call congrad_stable(index)

  ! Construct ghost surface

  return

end subroutine ghost

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine calcfg(index)
  use globals, only: dp, zero, half, pi2, ncpu, myid, resbn_m, resbn_n, gsurf, Ncoils, ounit, MPI_COMM_FOCUS
  
  use mpi
  implicit none
  INTEGER, INTENT(IN) :: index

  INTEGER             :: astat, ierr, i, Nseg_stable, p, q, icoil
  REAL                :: Bxhold, Byhold, Bzhold
  REAL, ALLOCATABLE   :: zeta(:), gradosx(:), gradosy(:), gradosz(:), gradxsx(:), gradxsy(:), &
                         gradxsz(:), gradozetax(:), gradozetay(:), gradozetaz(:), gradxzetax(:), &
                         gradxzetay(:), gradxzetaz(:), gradothetax(:), gradothetay(:), &
                         gradothetaz(:), gradxthetax(:), gradxthetay(:), gradxthetaz(:), oBx(:), &
                         oBy(:), oBz(:), xBx(:), xBy(:), xBz(:), odRadx(:), odRady(:), odZadx(:), &
                         odZady(:), xdRadx(:), xdRady(:), xdZadx(:), xdZady(:), Radot(:), Zadot(:)

  Nseg_stable = gsurf(index)%Nseg_stable
  SALLOCATE( zeta, (1:Nseg_stable), 0.0)
  zeta(1:Nseg_stable) = gsurf(index)%zeta(1:Nseg_stable)

  ! Calculate s, theta, sdot, thetadot for o and x
  p = resbn_n
  q = resbn_m

  gsurf(index)%os(1:Nseg_stable) = 0.0
  gsurf(index)%xs(1:Nseg_stable) = 0.0
  gsurf(index)%otheta(1:Nseg_stable) = 0.0
  gsurf(index)%xtheta(1:Nseg_stable) = 0.0
  gsurf(index)%osdot(1:Nseg_stable) = 0.0
  gsurf(index)%xsdot(1:Nseg_stable) = 0.0
  gsurf(index)%othetadot(1:Nseg_stable) = 0.0
  gsurf(index)%xthetadot(1:Nseg_stable) = 0.0

  do i = 1, gsurf(index)%NF_stable
     gsurf(index)%os(1:Nseg_stable) = gsurf(index)%os(1:Nseg_stable) + &
             gsurf(index)%osnc(i)*cos(gsurf(index)%on(i)*zeta(1:Nseg_stable)/q) + &
             gsurf(index)%osns(i)*sin(gsurf(index)%on(i)*zeta(1:Nseg_stable)/q)
     gsurf(index)%xs(1:Nseg_stable) = gsurf(index)%xs(1:Nseg_stable) + &
             gsurf(index)%xsnc(i)*cos(gsurf(index)%xn(i)*zeta(1:Nseg_stable)/q) + &
             gsurf(index)%xsns(i)*sin(gsurf(index)%xn(i)*zeta(1:Nseg_stable)/q)
     gsurf(index)%otheta(1:Nseg_stable) = gsurf(index)%otheta(1:Nseg_stable) + &
             gsurf(index)%othetanc(i)*cos(gsurf(index)%on(i)*zeta(1:Nseg_stable)/q) + &
             gsurf(index)%othetans(i)*sin(gsurf(index)%on(i)*zeta(1:Nseg_stable)/q)
     gsurf(index)%xtheta(1:Nseg_stable) = gsurf(index)%xtheta(1:Nseg_stable) + &
             gsurf(index)%xthetanc(i)*cos(gsurf(index)%xn(i)*zeta(1:Nseg_stable)/q) + &
             gsurf(index)%xthetans(i)*sin(gsurf(index)%xn(i)*zeta(1:Nseg_stable)/q)

     gsurf(index)%osdot(1:Nseg_stable) = gsurf(index)%osdot(1:Nseg_stable) - &
             gsurf(index)%on(i)*gsurf(index)%osnc(i)*sin(gsurf(index)%on(i)*zeta(1:Nseg_stable)/q)/q + &
             gsurf(index)%on(i)*gsurf(index)%osns(i)*cos(gsurf(index)%on(i)*zeta(1:Nseg_stable)/q)/q
     gsurf(index)%xsdot(1:Nseg_stable) = gsurf(index)%xsdot(1:Nseg_stable) - &
             gsurf(index)%xn(i)*gsurf(index)%xsnc(i)*sin(gsurf(index)%xn(i)*zeta(1:Nseg_stable)/q)/q + &
             gsurf(index)%xn(i)*gsurf(index)%xsns(i)*cos(gsurf(index)%xn(i)*zeta(1:Nseg_stable)/q)/q
     gsurf(index)%othetadot(1:Nseg_stable) = gsurf(index)%othetadot(1:Nseg_stable) - &
             gsurf(index)%on(i)*gsurf(index)%othetanc(i)*sin(gsurf(index)%on(i)*zeta(1:Nseg_stable)/q)/q + &
             gsurf(index)%on(i)*gsurf(index)%othetans(i)*cos(gsurf(index)%on(i)*zeta(1:Nseg_stable)/q)/q
     gsurf(index)%xthetadot(1:Nseg_stable) = gsurf(index)%xthetadot(1:Nseg_stable) - &
             gsurf(index)%xn(i)*gsurf(index)%xthetanc(i)*sin(gsurf(index)%xn(i)*zeta(1:Nseg_stable)/q)/q + &
             gsurf(index)%xn(i)*gsurf(index)%xthetans(i)*cos(gsurf(index)%xn(i)*zeta(1:Nseg_stable)/q)/q
  enddo
  gsurf(index)%otheta(1:Nseg_stable) = gsurf(index)%otheta(1:Nseg_stable) + p*zeta(1:Nseg_stable)/q
  gsurf(index)%xtheta(1:Nseg_stable) = gsurf(index)%xtheta(1:Nseg_stable) + p*zeta(1:Nseg_stable)/q

  gsurf(index)%othetadot(1:Nseg_stable) = gsurf(index)%othetadot(1:Nseg_stable) + p/q
  gsurf(index)%xthetadot(1:Nseg_stable) = gsurf(index)%xthetadot(1:Nseg_stable) + p/q

  ! Calculate axis, and derivative of axis w.r.t zeta
  gsurf(index)%Ra(1:Nseg_stable) = 0.0
  gsurf(index)%Za(1:Nseg_stable) = 0.0
  SALLOCATE( Radot, (1:Nseg_stable), 0.0 )
  SALLOCATE( Zadot, (1:Nseg_stable), 0.0 )
  do i = 1, gsurf(index)%NF_axis
     gsurf(index)%Ra(1:Nseg_stable) = gsurf(index)%Ra(1:Nseg_stable) + &
             gsurf(index)%axisrnc(i)*cos(gsurf(index)%axisn(i)*zeta(1:Nseg_stable))
     gsurf(index)%Za(1:Nseg_stable) = gsurf(index)%Za(1:Nseg_stable) + &
             gsurf(index)%axiszns(i)*sin(gsurf(index)%axisn(i)*zeta(1:Nseg_stable))
     Radot(1:Nseg_stable) = Radot(1:Nseg_stable) - &
             gsurf(index)%axisn(i)*gsurf(index)%axisrnc(i)*sin(gsurf(index)%axisn(i)*zeta(1:Nseg_stable))
     Zadot(1:Nseg_stable) = Zadot(1:Nseg_stable) + &
             gsurf(index)%axisn(i)*gsurf(index)%axiszns(i)*cos(gsurf(index)%axisn(i)*zeta(1:Nseg_stable))
  enddo

  ! Transfrom field line from s,zeta,theta to x,y,z
  gsurf(index)%ox(1:Nseg_stable) = cos(zeta(1:Nseg_stable))*(gsurf(index)%Ra(1:Nseg_stable) + &
          gsurf(index)%os(1:Nseg_stable)*cos(gsurf(index)%otheta(1:Nseg_stable)))
  gsurf(index)%xx(1:Nseg_stable) = cos(zeta(1:Nseg_stable))*(gsurf(index)%Ra(1:Nseg_stable) + &
          gsurf(index)%xs(1:Nseg_stable)*cos(gsurf(index)%xtheta(1:Nseg_stable)))
  gsurf(index)%oy(1:Nseg_stable) = sin(zeta(1:Nseg_stable))*(gsurf(index)%Ra(1:Nseg_stable) + &
          gsurf(index)%os(1:Nseg_stable)*cos(gsurf(index)%otheta(1:Nseg_stable)))
  gsurf(index)%xy(1:Nseg_stable) = sin(zeta(1:Nseg_stable))*(gsurf(index)%Ra(1:Nseg_stable) + &
          gsurf(index)%xs(1:Nseg_stable)*cos(gsurf(index)%xtheta(1:Nseg_stable)))
  gsurf(index)%oz(1:Nseg_stable) = gsurf(index)%os(1:Nseg_stable)*sin(gsurf(index)%otheta(1:Nseg_stable)) + &
          gsurf(index)%Za(1:Nseg_stable)
  gsurf(index)%xz(1:Nseg_stable) = gsurf(index)%xs(1:Nseg_stable)*sin(gsurf(index)%xtheta(1:Nseg_stable)) + &
          gsurf(index)%Za(1:Nseg_stable)

  ! Calculate derivative of field line w.r.t zeta
  gsurf(index)%oxdot(1:Nseg_stable) = -1.0*sin(zeta(1:Nseg_stable))*(gsurf(index)%Ra(1:Nseg_stable)&
         + gsurf(index)%os(1:Nseg_stable)*cos(gsurf(index)%otheta(1:Nseg_stable))) &
         + cos(zeta(1:Nseg_stable))*( Radot(1:Nseg_stable) + gsurf(index)%osdot(1:Nseg_stable)&
         *cos(gsurf(index)%otheta(1:Nseg_stable)) - gsurf(index)%os(1:Nseg_stable)*&
         sin(gsurf(index)%otheta(1:Nseg_stable))*gsurf(index)%othetadot(1:Nseg_stable) )
  gsurf(index)%xxdot(1:Nseg_stable) = -1.0*sin(zeta(1:Nseg_stable))*(gsurf(index)%Ra(1:Nseg_stable)& 
         + gsurf(index)%xs(1:Nseg_stable)*cos(gsurf(index)%xtheta(1:Nseg_stable))) &
         + cos(zeta(1:Nseg_stable))*( Radot(1:Nseg_stable) + gsurf(index)%xsdot(1:Nseg_stable)&
         *cos(gsurf(index)%xtheta(1:Nseg_stable)) - gsurf(index)%xs(1:Nseg_stable)*&
         sin(gsurf(index)%xtheta(1:Nseg_stable))*gsurf(index)%xthetadot(1:Nseg_stable) )

  gsurf(index)%oydot(1:Nseg_stable) =      cos(zeta(1:Nseg_stable))*(gsurf(index)%Ra(1:Nseg_stable)& 
         + gsurf(index)%os(1:Nseg_stable)*cos(gsurf(index)%otheta(1:Nseg_stable))) &
         + sin(zeta(1:Nseg_stable))*( Radot(1:Nseg_stable) + gsurf(index)%osdot(1:Nseg_stable)&
         *cos(gsurf(index)%otheta(1:Nseg_stable)) - gsurf(index)%os(1:Nseg_stable)*&
         sin(gsurf(index)%otheta(1:Nseg_stable))*gsurf(index)%othetadot(1:Nseg_stable) )
  gsurf(index)%xydot(1:Nseg_stable) =      cos(zeta(1:Nseg_stable))*(gsurf(index)%Ra(1:Nseg_stable)&
         + gsurf(index)%xs(1:Nseg_stable)*cos(gsurf(index)%xtheta(1:Nseg_stable))) &
         + sin(zeta(1:Nseg_stable))*( Radot(1:Nseg_stable) + gsurf(index)%xsdot(1:Nseg_stable)&
         *cos(gsurf(index)%xtheta(1:Nseg_stable)) - gsurf(index)%xs(1:Nseg_stable)*&
         sin(gsurf(index)%xtheta(1:Nseg_stable))*gsurf(index)%xthetadot(1:Nseg_stable) )

  gsurf(index)%ozdot(1:Nseg_stable) = gsurf(index)%osdot(1:Nseg_stable)*&
          sin(gsurf(index)%otheta(1:Nseg_stable)) + gsurf(index)%os(1:Nseg_stable)&
          *cos(gsurf(index)%otheta(1:Nseg_stable))*gsurf(index)%othetadot(1:Nseg_stable) + Zadot(1:Nseg_stable)
  gsurf(index)%xzdot(1:Nseg_stable) = gsurf(index)%xsdot(1:Nseg_stable)*&
          sin(gsurf(index)%xtheta(1:Nseg_stable)) + gsurf(index)%xs(1:Nseg_stable)&
          *cos(gsurf(index)%xtheta(1:Nseg_stable))*gsurf(index)%xthetadot(1:Nseg_stable) + Zadot(1:Nseg_stable)

  ! Calculate derivative of axis w.r.t x and y, variables are local
  SALLOCATE( odRadx, (1:Nseg_stable), 0.0 )
  SALLOCATE( odRady, (1:Nseg_stable), 0.0 )
  SALLOCATE( odZadx, (1:Nseg_stable), 0.0 )
  SALLOCATE( odZady, (1:Nseg_stable), 0.0 )
  SALLOCATE( xdRadx, (1:Nseg_stable), 0.0 )
  SALLOCATE( xdRady, (1:Nseg_stable), 0.0 )
  SALLOCATE( xdZadx, (1:Nseg_stable), 0.0 )
  SALLOCATE( xdZady, (1:Nseg_stable), 0.0 )
  do i = 1, gsurf(index)%NF_axis
     odRadx(1:Nseg_stable) = odRadx(1:Nseg_stable) + &
             gsurf(index)%axisn(i)*gsurf(index)%axisrnc(i)*sin(gsurf(index)%axisn(i)*zeta(1:Nseg_stable))*gsurf(index)%oy(1:Nseg_stable) / &
             ( gsurf(index)%ox(1:Nseg_stable)**2 * ( gsurf(index)%oy(1:Nseg_stable)**2/gsurf(index)%ox(1:Nseg_stable)**2 + 1 ) )
     odRady(1:Nseg_stable) = odRady(1:Nseg_stable) - &
             gsurf(index)%axisn(i)*gsurf(index)%axisrnc(i)*sin(gsurf(index)%axisn(i)*zeta(1:Nseg_stable)) / &
             ( gsurf(index)%ox(1:Nseg_stable) * ( gsurf(index)%oy(1:Nseg_stable)**2/gsurf(index)%ox(1:Nseg_stable)**2 + 1 ) )
     odZadx(1:Nseg_stable) = odZadx(1:Nseg_stable) - &
             gsurf(index)%axisn(i)*gsurf(index)%axiszns(i)*cos(gsurf(index)%axisn(i)*zeta(1:Nseg_stable))*gsurf(index)%oy(1:Nseg_stable) / &
             ( gsurf(index)%ox(1:Nseg_stable)**2 * ( gsurf(index)%oy(1:Nseg_stable)**2/gsurf(index)%ox(1:Nseg_stable)**2 + 1 ) )
     odZady(1:Nseg_stable) = odZady(1:Nseg_stable) + &
             gsurf(index)%axisn(i)*gsurf(index)%axiszns(i)*cos(gsurf(index)%axisn(i)*zeta(1:Nseg_stable)) / &
             ( gsurf(index)%ox(1:Nseg_stable) * ( gsurf(index)%oy(1:Nseg_stable)**2/gsurf(index)%ox(1:Nseg_stable)**2 + 1 ) )

     xdRadx(1:Nseg_stable) = xdRadx(1:Nseg_stable) + &
             gsurf(index)%axisn(i)*gsurf(index)%axisrnc(i)*sin(gsurf(index)%axisn(i)*zeta(1:Nseg_stable))*gsurf(index)%xy(1:Nseg_stable) / &
             ( gsurf(index)%xx(1:Nseg_stable)**2 * ( gsurf(index)%xy(1:Nseg_stable)**2/gsurf(index)%xx(1:Nseg_stable)**2 + 1 ) )
     xdRady(1:Nseg_stable) = xdRady(1:Nseg_stable) - &
             gsurf(index)%axisn(i)*gsurf(index)%axisrnc(i)*sin(gsurf(index)%axisn(i)*zeta(1:Nseg_stable)) / &
             ( gsurf(index)%xx(1:Nseg_stable) * ( gsurf(index)%xy(1:Nseg_stable)**2/gsurf(index)%xx(1:Nseg_stable)**2 + 1 ) )
     xdZadx(1:Nseg_stable) = xdZadx(1:Nseg_stable) - &
             gsurf(index)%axisn(i)*gsurf(index)%axiszns(i)*cos(gsurf(index)%axisn(i)*zeta(1:Nseg_stable))*gsurf(index)%xy(1:Nseg_stable) / &
             ( gsurf(index)%xx(1:Nseg_stable)**2 * ( gsurf(index)%xy(1:Nseg_stable)**2/gsurf(index)%xx(1:Nseg_stable)**2 + 1 ) )
     xdZady(1:Nseg_stable) = xdZady(1:Nseg_stable) + &
             gsurf(index)%axisn(i)*gsurf(index)%axiszns(i)*cos(gsurf(index)%axisn(i)*zeta(1:Nseg_stable)) / &
             ( gsurf(index)%xx(1:Nseg_stable) * ( gsurf(index)%xy(1:Nseg_stable)**2/gsurf(index)%xx(1:Nseg_stable)**2 + 1 ) )
  enddo 

  ! Calculate contravariant basis, variables are local
  SALLOCATE( gradosx, (1:Nseg_stable), 0.0 )
  SALLOCATE( gradosy, (1:Nseg_stable), 0.0 )
  SALLOCATE( gradosz, (1:Nseg_stable), 0.0 )
  SALLOCATE( gradxsx, (1:Nseg_stable), 0.0 )
  SALLOCATE( gradxsy, (1:Nseg_stable), 0.0 )
  SALLOCATE( gradxsz, (1:Nseg_stable), 0.0 )
  
  gradosx(1:Nseg_stable) = ( sqrt( gsurf(index)%ox(1:Nseg_stable)**2 + gsurf(index)%oy(1:Nseg_stable)**2 ) - &
          gsurf(index)%Ra(1:Nseg_stable) )*( (gsurf(index)%ox(1:Nseg_stable)**2 + &
          gsurf(index)%oy(1:Nseg_stable)**2)**(-.5)*gsurf(index)%ox(1:Nseg_stable) - &
          odRadx(1:Nseg_stable) )/gsurf(index)%os(1:Nseg_stable) - &
          (gsurf(index)%oz(1:Nseg_stable) - gsurf(index)%Za(1:Nseg_stable))*odZadx(1:Nseg_stable)/gsurf(index)%os(1:Nseg_stable)
  gradxsx(1:Nseg_stable) = ( sqrt( gsurf(index)%xx(1:Nseg_stable)**2 + gsurf(index)%xy(1:Nseg_stable)**2 ) - &
          gsurf(index)%Ra(1:Nseg_stable) )*( (gsurf(index)%xx(1:Nseg_stable)**2 + &
          gsurf(index)%xy(1:Nseg_stable)**2)**(-.5)*gsurf(index)%xx(1:Nseg_stable) - & 
          xdRadx(1:Nseg_stable) )/gsurf(index)%xs(1:Nseg_stable) - &
          (gsurf(index)%xz(1:Nseg_stable) - gsurf(index)%Za(1:Nseg_stable))*xdZadx(1:Nseg_stable)/gsurf(index)%xs(1:Nseg_stable)
  
  gradosy(1:Nseg_stable) = ( sqrt( gsurf(index)%ox(1:Nseg_stable)**2 + gsurf(index)%oy(1:Nseg_stable)**2 ) - &
          gsurf(index)%Ra(1:Nseg_stable) )*( (gsurf(index)%ox(1:Nseg_stable)**2 + &
          gsurf(index)%oy(1:Nseg_stable)**2)**(-.5)*gsurf(index)%oy(1:Nseg_stable) - &
          odRady(1:Nseg_stable) )/gsurf(index)%os(1:Nseg_stable) - &
          (gsurf(index)%oz(1:Nseg_stable) - gsurf(index)%Za(1:Nseg_stable))*odZady(1:Nseg_stable)/gsurf(index)%os(1:Nseg_stable)
  gradxsy(1:Nseg_stable) = ( sqrt( gsurf(index)%xx(1:Nseg_stable)**2 + gsurf(index)%xy(1:Nseg_stable)**2 ) - &
          gsurf(index)%Ra(1:Nseg_stable) )*( (gsurf(index)%xx(1:Nseg_stable)**2 + &
          gsurf(index)%xy(1:Nseg_stable)**2)**(-.5)*gsurf(index)%xy(1:Nseg_stable) - &
          xdRady(1:Nseg_stable) )/gsurf(index)%xs(1:Nseg_stable) - &
          (gsurf(index)%xz(1:Nseg_stable) - gsurf(index)%Za(1:Nseg_stable))*xdZady(1:Nseg_stable)/gsurf(index)%xs(1:Nseg_stable)
 
  gradosz(1:Nseg_stable) = ( gsurf(index)%oz(1:Nseg_stable) - gsurf(index)%Za(1:Nseg_stable) )/gsurf(index)%os(1:Nseg_stable)
  gradxsz(1:Nseg_stable) = ( gsurf(index)%xz(1:Nseg_stable) - gsurf(index)%Za(1:Nseg_stable) )/gsurf(index)%xs(1:Nseg_stable)

  SALLOCATE( gradozetax, (1:Nseg_stable), 0.0 )
  SALLOCATE( gradozetay, (1:Nseg_stable), 0.0 )
  SALLOCATE( gradozetaz, (1:Nseg_stable), 0.0 )
  SALLOCATE( gradxzetax, (1:Nseg_stable), 0.0 )
  SALLOCATE( gradxzetay, (1:Nseg_stable), 0.0 )
  SALLOCATE( gradxzetaz, (1:Nseg_stable), 0.0 )
  
  gradozetax(1:Nseg_stable) = -1.0*(gsurf(index)%oy(1:Nseg_stable)**2/gsurf(index)%ox(1:Nseg_stable)**2+1)**-1.0 * &
          gsurf(index)%oy(1:Nseg_stable)/gsurf(index)%ox(1:Nseg_stable)**2
  gradxzetax(1:Nseg_stable) = -1.0*(gsurf(index)%xy(1:Nseg_stable)**2/gsurf(index)%xx(1:Nseg_stable)**2+1)**-1.0 * &
          gsurf(index)%xy(1:Nseg_stable)/gsurf(index)%xx(1:Nseg_stable)**2
  
  gradozetay(1:Nseg_stable) = (gsurf(index)%oy(1:Nseg_stable)**2/gsurf(index)%ox(1:Nseg_stable)**2+1)**-1.0 * &
          1.0/gsurf(index)%ox(1:Nseg_stable)
  gradxzetay(1:Nseg_stable) = (gsurf(index)%xy(1:Nseg_stable)**2/gsurf(index)%xx(1:Nseg_stable)**2+1)**-1.0 * &
          1.0/gsurf(index)%xx(1:Nseg_stable)

  !gradozetaz(1:Nseg_stable) = 0.0
  !gradxzetaz(1:Nseg_stable) = 0.0

  SALLOCATE( gradothetax, (1:Nseg_stable), 0.0 )
  SALLOCATE( gradothetay, (1:Nseg_stable), 0.0 )
  SALLOCATE( gradothetaz, (1:Nseg_stable), 0.0 )
  SALLOCATE( gradxthetax, (1:Nseg_stable), 0.0 )
  SALLOCATE( gradxthetay, (1:Nseg_stable), 0.0 )
  SALLOCATE( gradxthetaz, (1:Nseg_stable), 0.0 )

  gradothetax(1:Nseg_stable) = -1.0*( (gsurf(index)%oz(1:Nseg_stable)-gsurf(index)%Za(1:Nseg_stable)) * &
          (sqrt(gsurf(index)%ox(1:Nseg_stable)**2+gsurf(index)%oy(1:Nseg_stable)**2) - gsurf(index)%Ra(1:Nseg_stable))**-2.0 * &
          ( (gsurf(index)%ox(1:Nseg_stable)**2+gsurf(index)%oy(1:Nseg_stable)**2)**-.5 * gsurf(index)%ox(1:Nseg_stable) - odRadx(1:Nseg_stable) ) &
          + odZadx(1:Nseg_stable) / (sqrt(gsurf(index)%ox(1:Nseg_stable)**2+gsurf(index)%oy(1:Nseg_stable)**2) &
          - gsurf(index)%Ra(1:Nseg_stable)) ) / (tan(gsurf(index)%otheta(1:Nseg_stable))**2+1)
  gradxthetax(1:Nseg_stable) = -1.0*( (gsurf(index)%xz(1:Nseg_stable)-gsurf(index)%Za(1:Nseg_stable)) * &
          (sqrt(gsurf(index)%xx(1:Nseg_stable)**2+gsurf(index)%xy(1:Nseg_stable)**2) - gsurf(index)%Ra(1:Nseg_stable))**-2.0 * &
          ( (gsurf(index)%xx(1:Nseg_stable)**2+gsurf(index)%xy(1:Nseg_stable)**2)**-.5 * gsurf(index)%xx(1:Nseg_stable) - xdRadx(1:Nseg_stable) ) &
          + xdZadx(1:Nseg_stable) / (sqrt(gsurf(index)%xx(1:Nseg_stable)**2+gsurf(index)%xy(1:Nseg_stable)**2) &
          - gsurf(index)%Ra(1:Nseg_stable)) ) / (tan(gsurf(index)%xtheta(1:Nseg_stable))**2+1)

  gradothetay(1:Nseg_stable) = -1.0*( (gsurf(index)%oz(1:Nseg_stable)-gsurf(index)%Za(1:Nseg_stable)) * &
          (sqrt(gsurf(index)%ox(1:Nseg_stable)**2+gsurf(index)%oy(1:Nseg_stable)**2) - gsurf(index)%Ra(1:Nseg_stable))**-2.0 * &
          ( (gsurf(index)%ox(1:Nseg_stable)**2+gsurf(index)%oy(1:Nseg_stable)**2)**-.5 * gsurf(index)%oy(1:Nseg_stable) - odRady(1:Nseg_stable) ) &
          + odZady(1:Nseg_stable) / (sqrt(gsurf(index)%ox(1:Nseg_stable)**2+gsurf(index)%oy(1:Nseg_stable)**2) &
          - gsurf(index)%Ra(1:Nseg_stable)) ) / (tan(gsurf(index)%otheta(1:Nseg_stable))**2+1)
  gradxthetay(1:Nseg_stable) = -1.0*( (gsurf(index)%xz(1:Nseg_stable)-gsurf(index)%Za(1:Nseg_stable)) * &
          (sqrt(gsurf(index)%xx(1:Nseg_stable)**2+gsurf(index)%xy(1:Nseg_stable)**2) - gsurf(index)%Ra(1:Nseg_stable))**-2.0 * &
          ( (gsurf(index)%xx(1:Nseg_stable)**2+gsurf(index)%xy(1:Nseg_stable)**2)**-.5 * gsurf(index)%xy(1:Nseg_stable) - xdRady(1:Nseg_stable) ) &
          + xdZady(1:Nseg_stable) / (sqrt(gsurf(index)%xx(1:Nseg_stable)**2+gsurf(index)%xy(1:Nseg_stable)**2) &
          - gsurf(index)%Ra(1:Nseg_stable)) ) / (tan(gsurf(index)%xtheta(1:Nseg_stable))**2+1)

  gradothetaz(1:Nseg_stable) = ( (tan(gsurf(index)%otheta(1:Nseg_stable))**2+1) * &
          ( sqrt(gsurf(index)%ox(1:Nseg_stable)**2+gsurf(index)%oy(1:Nseg_stable)**2) - gsurf(index)%Ra(1:Nseg_stable) ) )**-1.0
  gradxthetaz(1:Nseg_stable) = ( (tan(gsurf(index)%xtheta(1:Nseg_stable))**2+1) * &
          ( sqrt(gsurf(index)%xx(1:Nseg_stable)**2+gsurf(index)%xy(1:Nseg_stable)**2) - gsurf(index)%Ra(1:Nseg_stable) ) )**-1.0

  ! Dont need to calculate Jacobian since not calculating covariant basis
  !J =     gradsx*( gradthetay*gradzetaz - gradthetaz*gradzetay )
  !J = J + gradsy*( gradthetaz*gradzetax - gradthetax*gradzetaz )
  !J = J + gradsz*( gradthetax*gradzetay - gradthetay*gradzetax )
  !J = -1.0*J**-1.0 ! Made negative since switched angles 

  ! Calculate field on stable field lines
  SALLOCATE( oBx, (1:Nseg_stable), 0.0 )
  SALLOCATE( oBy, (1:Nseg_stable), 0.0 )
  SALLOCATE( oBz, (1:Nseg_stable), 0.0 )
  SALLOCATE( xBx, (1:Nseg_stable), 0.0 )
  SALLOCATE( xBy, (1:Nseg_stable), 0.0 )
  SALLOCATE( xBz, (1:Nseg_stable), 0.0 )

  ! Parallelized, maybe change to be parallelized elsewhere
  do i = 1, Nseg_stable
     if( myid.ne.modulo(i,ncpu) ) cycle ! parallelization loop;
     do icoil = 1, Ncoils
        call bfield0(icoil, gsurf(index)%ox(i), gsurf(index)%oy(i), gsurf(index)%oz(i), Bxhold, Byhold, Bzhold)
        oBx(i) = oBx(i) + Bxhold
        oBy(i) = oBy(i) + Byhold
        oBz(i) = oBz(i) + Bzhold

        call bfield0(icoil, gsurf(index)%xx(i), gsurf(index)%xy(i), gsurf(index)%xz(i), Bxhold, Byhold, Bzhold)
        xBx(i) = xBx(i) + Bxhold
        xBy(i) = xBy(i) + Byhold
        xBz(i) = xBz(i) + Bzhold
     enddo
  enddo
  call MPI_BARRIER( MPI_COMM_FOCUS, ierr )
  call MPI_ALLREDUCE( MPI_IN_PLACE, oBx, Nseg_stable, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FOCUS, ierr )
  call MPI_ALLREDUCE( MPI_IN_PLACE, oBy, Nseg_stable, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FOCUS, ierr )
  call MPI_ALLREDUCE( MPI_IN_PLACE, oBz, Nseg_stable, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FOCUS, ierr )
  call MPI_ALLREDUCE( MPI_IN_PLACE, xBx, Nseg_stable, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FOCUS, ierr )
  call MPI_ALLREDUCE( MPI_IN_PLACE, xBy, Nseg_stable, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FOCUS, ierr )
  call MPI_ALLREDUCE( MPI_IN_PLACE, xBz, Nseg_stable, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FOCUS, ierr )

  ! Calculate B^s, B^zeta, B^theta on stable field lines
  gsurf(index)%obsups(1:Nseg_stable) = oBx(1:Nseg_stable)*gradosx(1:Nseg_stable) + &
          oBy(1:Nseg_stable)*gradosy(1:Nseg_stable) + oBz(1:Nseg_stable)*gradosz(1:Nseg_stable)
  gsurf(index)%xbsups(1:Nseg_stable) = xBx(1:Nseg_stable)*gradxsx(1:Nseg_stable) + &
          xBy(1:Nseg_stable)*gradxsy(1:Nseg_stable) + xBz(1:Nseg_stable)*gradxsz(1:Nseg_stable)

  gsurf(index)%obsupzeta(1:Nseg_stable) = oBx(1:Nseg_stable)*gradozetax(1:Nseg_stable) + &
          oBy(1:Nseg_stable)*gradozetay(1:Nseg_stable) + oBz(1:Nseg_stable)*gradozetaz(1:Nseg_stable)
  gsurf(index)%xbsupzeta(1:Nseg_stable) = xBx(1:Nseg_stable)*gradxzetax(1:Nseg_stable) + &
          xBy(1:Nseg_stable)*gradxzetay(1:Nseg_stable) + xBz(1:Nseg_stable)*gradxzetaz(1:Nseg_stable)

  gsurf(index)%obsuptheta(1:Nseg_stable) = oBx(1:Nseg_stable)*gradothetax(1:Nseg_stable) + &
          oBy(1:Nseg_stable)*gradothetay(1:Nseg_stable) + oBz(1:Nseg_stable)*gradothetaz(1:Nseg_stable)
  gsurf(index)%xbsuptheta(1:Nseg_stable) = xBx(1:Nseg_stable)*gradxthetax(1:Nseg_stable) + &
          xBy(1:Nseg_stable)*gradxthetay(1:Nseg_stable) + xBz(1:Nseg_stable)*gradxthetaz(1:Nseg_stable)

  ! Calculate f,g and F
  gsurf(index)%of(1:Nseg_stable) = gsurf(index)%obsups(1:Nseg_stable)/gsurf(index)%obsupzeta(1:Nseg_stable) - &
          gsurf(index)%osdot(1:Nseg_stable)
  gsurf(index)%xf(1:Nseg_stable) = gsurf(index)%xbsups(1:Nseg_stable)/gsurf(index)%xbsupzeta(1:Nseg_stable) - &
          gsurf(index)%xsdot(1:Nseg_stable)

  gsurf(index)%og(1:Nseg_stable) = gsurf(index)%obsuptheta(1:Nseg_stable)/gsurf(index)%obsupzeta(1:Nseg_stable) - &
          gsurf(index)%othetadot(1:Nseg_stable)
  gsurf(index)%xg(1:Nseg_stable) = gsurf(index)%xbsuptheta(1:Nseg_stable)/gsurf(index)%xbsupzeta(1:Nseg_stable) - &
          gsurf(index)%xthetadot(1:Nseg_stable)

  gsurf(index)%F = 0.0
  do i = 1, Nseg_stable-1
     gsurf(index)%F = gsurf(index)%F + gsurf(index)%of(i)**2 + gsurf(index)%og(i)**2 + &
             gsurf(index)%xf(i)**2 + gsurf(index)%xg(i)**2
  enddo
  gsurf(index)%F = pi2*q*gsurf(index)%F/(Nseg_stable-1)

  ! DALLOCATE
  DALLOCATE( zeta )
  DALLOCATE( odRadx )
  DALLOCATE( odRady )
  DALLOCATE( odZadx )
  DALLOCATE( odZady )
  DALLOCATE( xdRadx )
  DALLOCATE( xdRady )
  DALLOCATE( xdZadx )
  DALLOCATE( xdZady )
  DALLOCATE( gradosx )
  DALLOCATE( gradosy )
  DALLOCATE( gradosz )
  DALLOCATE( gradxsx )
  DALLOCATE( gradxsy )
  DALLOCATE( gradxsz )
  DALLOCATE( gradozetax )
  DALLOCATE( gradozetay )
  DALLOCATE( gradozetaz )
  DALLOCATE( gradxzetax )
  DALLOCATE( gradxzetay )
  DALLOCATE( gradxzetaz )
  DALLOCATE( gradothetax )
  DALLOCATE( gradothetay )
  DALLOCATE( gradothetaz )
  DALLOCATE( gradxthetax )
  DALLOCATE( gradxthetay )
  DALLOCATE( gradxthetaz )
  DALLOCATE( oBx )
  DALLOCATE( oBy )
  DALLOCATE( oBz )
  DALLOCATE( xBx )
  DALLOCATE( xBy )
  DALLOCATE( xBz )

  return

end subroutine calcfg

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine calcfg_deriv(index)
  use globals, only: dp, zero, half, machprec, sqrtmachprec, ncpu, myid, ounit, MPI_COMM_FOCUS, &
                               psmall, gsurf
  use mpi
  implicit none

  INTEGER, INTENT(in)  :: index
  
  INTEGER              :: astat, ierr, idof, NF_stable
  REAL                 :: tmp_xdof(1:gsurf(index)%Ndof_stable), fd, negvalue, posvalue
  REAL                 :: start, finish, small
  
  do idof = 1, gsurf(index)%Ndof_stable
     ! perturbation will be relative.
     small = gsurf(index)%xdof_stable(idof) * psmall
     if (abs(small)<machprec) small = psmall
     !backward pertubation;
     tmp_xdof = gsurf(index)%xdof_stable
     tmp_xdof(idof) = tmp_xdof(idof) - half * small
     call unpacking_stable(tmp_xdof,index)
     call calcfg(index)
     negvalue = gsurf(index)%F
     !forward pertubation;
     tmp_xdof = gsurf(index)%xdof_stable
     tmp_xdof(idof) = tmp_xdof(idof) + half * small
     call unpacking_stable(tmp_xdof,index)
     call calcfg(index)
     posvalue = gsurf(index)%F
     !finite difference;
     gsurf(index)%dFdxdof_stable(idof) = (posvalue - negvalue) / small
  enddo

  !NF_stable = gsurf(index)%NF_stable
  !
  !gsurf(index)%dFdxdof_stable(            1:  NF_stable) = dFdosnc(1:NF_stable)
  !gsurf(index)%dFdxdof_stable(  NF_stabel+1:2*NF_stable) = dFdosns(1:NF_stable)
  !gsurf(index)%dFdxdof_stable(2*NF_stable+1:3*NF_stable) = dFdothetanc(1:NF_stable)
  !gsurf(index)%dFdxdof_stable(3*NF_stable+1:4*NF_stable) = dFdothetans(1:NF_stable)
  ! inlude x terms also 

  return

end subroutine calcfg_deriv

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine congrad_stable(index)
  use globals, only: dp, sqrtmachprec, myid, ounit, CG_maxiter, CG_xtol, &
       exit_signal, tstart, tfinish, gsurf, MPI_COMM_FOCUS
  
  use mpi
  implicit none

  INTEGER, INTENT(in)     :: index
  INTEGER                 :: ierr, astat, iter, n, nfunc, ngrad, status, CG_maxiter_hold, NF_stable
  REAL                    :: f, gnorm, CG_xtol_hold
  REAL, dimension(1:gsurf(index)%Ndof_stable) :: x, g, d, xtemp, gtemp
  EXTERNAL                :: myvalue_stable, mygrad_stable

  iter = 0
  n = gsurf(index)%Ndof_stable
  x(1:n) = gsurf(index)%xdof_stable(1:n)

  CG_maxiter_hold = CG_maxiter
  CG_maxiter = 10
  CG_xtol_hold = CG_xtol
  CG_xtol = 1.0E-6

  call cg_descent (CG_xtol, x, n, myvalue_stable, mygrad_stable, status, gnorm, f, iter, nfunc, ngrad, d, g, xtemp, gtemp)

  if (myid == 0) then
     select case (status)
     case (0)
        write(ounit, '("congrad : status="I1": convergence tolerance satisfied.")')  status
     case (1)
        write(ounit, '("congrad : status="I1": change in func <= feps*|f|.")')  status
     case (2)
        write(ounit, '("congrad : status="I1": total iterations exceeded maxit.")')  status
     case (3)
        write(ounit, '("congrad : status="I1": slope always negative in line search.")')  status
     case (4)
        write(ounit, '("congrad : status="I1": number secant iterations exceed nsecant.")')  status
     case (5)
        write(ounit, '("congrad : status="I1": search direction not a descent direction.")')  status
     case (6)
        write(ounit, '("congrad : status="I1": line search fails in initial interval.")')  status
     case (7)
        write(ounit, '("congrad : status="I1": line search fails during bisection.")')  status
     case (8)
        write(ounit, '("congrad : status="I1": line search fails during interval update.")')  status
     case default
        write(ounit, '("congrad : status="I1": unknow options!")')  status
     end select
  end if

  if(myid .eq. 0) write(ounit, '("congrad : Stable conjugate gradient finished.")')

  CG_maxiter = CG_maxiter_hold
  CG_xtol = CG_xtol_hold

  NF_stable = gsurf(index)%NF_stable
  gsurf(index)%xdof_stable(            1:  NF_stable) = gsurf(index)%osnc(1:NF_stable)
  gsurf(index)%xdof_stable(  NF_stable+1:2*NF_stable) = gsurf(index)%osns(1:NF_stable)
  gsurf(index)%xdof_stable(2*NF_stable+1:3*NF_stable) = gsurf(index)%othetanc(1:NF_stable)
  gsurf(index)%xdof_stable(3*NF_stable+1:4*NF_stable) = gsurf(index)%othetans(1:NF_stable)
  gsurf(index)%xdof_stable(4*NF_stable+1:5*NF_stable) = gsurf(index)%xsnc(1:NF_stable)
  gsurf(index)%xdof_stable(5*NF_stable+1:6*NF_stable) = gsurf(index)%xsns(1:NF_stable)
  gsurf(index)%xdof_stable(6*NF_stable+1:7*NF_stable) = gsurf(index)%xthetanc(1:NF_stable)
  gsurf(index)%xdof_stable(7*NF_stable+1:8*NF_stable) = gsurf(index)%xthetans(1:NF_stable)

  return

end subroutine congrad_stable

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

subroutine myvalue_stable(f, x, n)
  use globals, only: dp, myid, ounit, ierr, gsurf, MPI_COMM_FOCUS
  use mpi
  implicit none

  INTEGER, INTENT(in) :: n
  REAL, INTENT(in)    :: x(n)
  REAL, INTENT(out)   :: f

  INTEGER             :: index

  call MPI_BARRIER( MPI_COMM_FOCUS, ierr ) ! wait all cpus;

  index = 1
  
  call unpacking_stable(x,index)
  call calcfg(index)
  f = gsurf(index)%F

  return

end subroutine myvalue_stable

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

subroutine mygrad_stable(g, x, n)
  use globals, only: dp, myid, ounit, ierr, gsurf, MPI_COMM_FOCUS
  use mpi
  implicit none

  INTEGER, INTENT(in) :: n
  REAL, INTENT(in)    :: x(n)
  REAL, INTENT(out)   :: g(n)

  INTEGER             :: index

  call MPI_BARRIER( MPI_COMM_FOCUS, ierr ) ! wait all cpus;
  
  index = 1

  call unpacking_stable(x,index)
  call calcfg_deriv(index)
  g = gsurf(index)%dFdxdof_stable

  return

end subroutine mygrad_stable

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

subroutine unpacking_stable(x, index)
  use globals, only: dp, myid, ounit, ierr, gsurf, MPI_COMM_FOCUS
  use mpi
  implicit none

  INTEGER, INTENT(in) :: index
  REAL, INTENT(in)    :: x(gsurf(index)%Ndof_stable)

  INTEGER             :: NF_stable

  call MPI_BARRIER( MPI_COMM_FOCUS, ierr ) ! wait all cpus;

  NF_stable = gsurf(index)%NF_stable
  gsurf(index)%osnc(1:NF_stable)     = x(            1:  NF_stable)
  gsurf(index)%osns(1:NF_stable)     = x(  NF_stable+1:2*NF_stable)
  gsurf(index)%othetanc(1:NF_stable) = x(2*NF_stable+1:3*NF_stable)
  gsurf(index)%othetans(1:NF_stable) = x(3*NF_stable+1:4*NF_stable)
  gsurf(index)%xsnc(1:NF_stable)     = x(4*NF_stable+1:5*NF_stable)
  gsurf(index)%xsns(1:NF_stable)     = x(5*NF_stable+1:6*NF_stable)
  gsurf(index)%xthetanc(1:NF_stable) = x(6*NF_stable+1:7*NF_stable)
  gsurf(index)%xthetans(1:NF_stable) = x(7*NF_stable+1:8*NF_stable)

  return

end subroutine unpacking_stable
