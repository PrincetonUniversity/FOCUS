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
       coil, Ncoils, Nfixgeo, resbn_m, gsurf, MPI_COMM_FOCUS
  
  use mpi
  implicit none
  INTEGER, INTENT(in) :: index

  INTEGER             :: astat, ierr

  ! Set hold variables if first call of iteration
  !if ( gsurf(index)%iter_track .eq. iter ) then
  !   gsurf(index)%xdof_stable_hold(1:gsurf(index)%Ndof_stable) = gsurf(index)%xdof_stable(1:gsurf(index)%Ndof_stable)
  !   gsurf(index)%iter_track = gsurf(index)%iter_track + 1
  !endif

  ! Reset variables if > first call of iteration
  !if ( gsurf(index)%iter_track > iter ) then
  !   gsurf(index)%xdof_stable(1:gsurf(index)%Ndof_stable) = gsurf(index)%xdof_stable_hold(1:gsurf(index)%Ndof_stable)
  !   call unpacking_stable(gsurf(index)%xdof_stable,index)
  !endif

  ! Reset variables to be initial 
  gsurf(index)%xdof_stable(1:gsurf(index)%Ndof_stable) = gsurf(index)%xdof_stable_hold(1:gsurf(index)%Ndof_stable)
  call unpacking_stable(gsurf(index)%xdof_stable,index)

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

  INTEGER             :: astat, ierr, i, j, Nseg_stable, icoil
  REAL                :: Bxhold, Byhold, Bzhold, p, q
  REAL, ALLOCATABLE   :: zeta(:), gradosx(:), gradosy(:), gradosz(:), gradxsx(:), gradxsy(:), &
                         gradxsz(:), gradozetax(:), gradozetay(:), gradozetaz(:), gradxzetax(:), &
                         gradxzetay(:), gradxzetaz(:), gradothetax(:), gradothetay(:), &
                         gradothetaz(:), gradxthetax(:), gradxthetay(:), gradxthetaz(:), oBx(:), &
                         oBy(:), oBz(:), xBx(:), xBy(:), xBz(:), odRadx(:), odRady(:), odZadx(:), &
                         odZady(:), xdRadx(:), xdRady(:), xdZadx(:), xdZady(:), Radot(:), Zadot(:)

  Nseg_stable = gsurf(index)%Nseg_stable
  SALLOCATE( zeta, (1:Nseg_stable), 0.0)
  zeta(1:Nseg_stable) = gsurf(index)%zeta(1:Nseg_stable)

  p = real(resbn_n)
  q = real(resbn_m)

  ! Allocate everything
  SALLOCATE( Radot, (1:Nseg_stable), 0.0 )
  SALLOCATE( Zadot, (1:Nseg_stable), 0.0 ) 
  SALLOCATE( odRadx, (1:Nseg_stable), 0.0 )
  SALLOCATE( odRady, (1:Nseg_stable), 0.0 )
  SALLOCATE( odZadx, (1:Nseg_stable), 0.0 )
  SALLOCATE( odZady, (1:Nseg_stable), 0.0 )
  SALLOCATE( xdRadx, (1:Nseg_stable), 0.0 )
  SALLOCATE( xdRady, (1:Nseg_stable), 0.0 )
  SALLOCATE( xdZadx, (1:Nseg_stable), 0.0 )
  SALLOCATE( xdZady, (1:Nseg_stable), 0.0 )
  SALLOCATE( gradosx, (1:Nseg_stable), 0.0 )
  SALLOCATE( gradosy, (1:Nseg_stable), 0.0 )
  SALLOCATE( gradosz, (1:Nseg_stable), 0.0 )
  SALLOCATE( gradxsx, (1:Nseg_stable), 0.0 )
  SALLOCATE( gradxsy, (1:Nseg_stable), 0.0 )
  SALLOCATE( gradxsz, (1:Nseg_stable), 0.0 )
  SALLOCATE( gradozetax, (1:Nseg_stable), 0.0 )
  SALLOCATE( gradozetay, (1:Nseg_stable), 0.0 )
  SALLOCATE( gradozetaz, (1:Nseg_stable), 0.0 )
  SALLOCATE( gradxzetax, (1:Nseg_stable), 0.0 )
  SALLOCATE( gradxzetay, (1:Nseg_stable), 0.0 )
  SALLOCATE( gradxzetaz, (1:Nseg_stable), 0.0 )
  SALLOCATE( gradothetax, (1:Nseg_stable), 0.0 )
  SALLOCATE( gradothetay, (1:Nseg_stable), 0.0 )
  SALLOCATE( gradothetaz, (1:Nseg_stable), 0.0 )
  SALLOCATE( gradxthetax, (1:Nseg_stable), 0.0 )
  SALLOCATE( gradxthetay, (1:Nseg_stable), 0.0 )
  SALLOCATE( gradxthetaz, (1:Nseg_stable), 0.0 )
  SALLOCATE( oBx, (1:Nseg_stable), 0.0 )
  SALLOCATE( oBy, (1:Nseg_stable), 0.0 )
  SALLOCATE( oBz, (1:Nseg_stable), 0.0 )
  SALLOCATE( xBx, (1:Nseg_stable), 0.0 )
  SALLOCATE( xBy, (1:Nseg_stable), 0.0 )
  SALLOCATE( xBz, (1:Nseg_stable), 0.0 )

  gsurf(index)%of(1:Nseg_stable) = 0.0
  gsurf(index)%xf(1:Nseg_stable) = 0.0
  gsurf(index)%og(1:Nseg_stable) = 0.0
  gsurf(index)%xg(1:Nseg_stable) = 0.0
  
  do j = 1, Nseg_stable ! start of unindented big loop
  if( myid.ne.modulo(j,ncpu) ) cycle ! parallelization 

  gsurf(index)%os(j) = 0.0
  gsurf(index)%xs(j) = 0.0
  gsurf(index)%otheta(j) = 0.0
  gsurf(index)%xtheta(j) = 0.0
  gsurf(index)%osdot(j) = 0.0
  gsurf(index)%xsdot(j) = 0.0
  gsurf(index)%othetadot(j) = 0.0
  gsurf(index)%xthetadot(j) = 0.0
  do i = 1, gsurf(index)%NF_stable
     gsurf(index)%os(j) = gsurf(index)%os(j) + &
             gsurf(index)%osnc(i)*cos(gsurf(index)%on(i)*zeta(j)/q) + &
             gsurf(index)%osns(i)*sin(gsurf(index)%on(i)*zeta(j)/q)
     gsurf(index)%xs(j) = gsurf(index)%xs(j) + &
             gsurf(index)%xsnc(i)*cos(gsurf(index)%xn(i)*zeta(j)/q) + &
             gsurf(index)%xsns(i)*sin(gsurf(index)%xn(i)*zeta(j)/q)
     gsurf(index)%otheta(j) = gsurf(index)%otheta(j) + &
             gsurf(index)%othetanc(i)*cos(gsurf(index)%on(i)*zeta(j)/q) + &
             gsurf(index)%othetans(i)*sin(gsurf(index)%on(i)*zeta(j)/q)
     gsurf(index)%xtheta(j) = gsurf(index)%xtheta(j) + &
             gsurf(index)%xthetanc(i)*cos(gsurf(index)%xn(i)*zeta(j)/q) + &
             gsurf(index)%xthetans(i)*sin(gsurf(index)%xn(i)*zeta(j)/q)
     gsurf(index)%osdot(j) = gsurf(index)%osdot(j) - &
             gsurf(index)%on(i)*gsurf(index)%osnc(i)*sin(gsurf(index)%on(i)*zeta(j)/q)/q + &
             gsurf(index)%on(i)*gsurf(index)%osns(i)*cos(gsurf(index)%on(i)*zeta(j)/q)/q
     gsurf(index)%xsdot(j) = gsurf(index)%xsdot(j) - &
             gsurf(index)%xn(i)*gsurf(index)%xsnc(i)*sin(gsurf(index)%xn(i)*zeta(j)/q)/q + &
             gsurf(index)%xn(i)*gsurf(index)%xsns(i)*cos(gsurf(index)%xn(i)*zeta(j)/q)/q
     gsurf(index)%othetadot(j) = gsurf(index)%othetadot(j) - &
             gsurf(index)%on(i)*gsurf(index)%othetanc(i)*sin(gsurf(index)%on(i)*zeta(j)/q)/q + &
             gsurf(index)%on(i)*gsurf(index)%othetans(i)*cos(gsurf(index)%on(i)*zeta(j)/q)/q
     gsurf(index)%xthetadot(j) = gsurf(index)%xthetadot(j) - &
             gsurf(index)%xn(i)*gsurf(index)%xthetanc(i)*sin(gsurf(index)%xn(i)*zeta(j)/q)/q + &
             gsurf(index)%xn(i)*gsurf(index)%xthetans(i)*cos(gsurf(index)%xn(i)*zeta(j)/q)/q
  enddo
  gsurf(index)%otheta(j) = gsurf(index)%otheta(j) + p*zeta(j)/q
  gsurf(index)%xtheta(j) = gsurf(index)%xtheta(j) + p*zeta(j)/q
  gsurf(index)%othetadot(j) = gsurf(index)%othetadot(j) + p/q
  gsurf(index)%xthetadot(j) = gsurf(index)%xthetadot(j) + p/q
  ! Calculate axis, and derivative of axis w.r.t zeta
  gsurf(index)%Ra(j) = 0.0
  gsurf(index)%Za(j) = 0.0
  do i = 1, gsurf(index)%NF_axis
     gsurf(index)%Ra(j) = gsurf(index)%Ra(j) + &
             gsurf(index)%axisrnc(i)*cos(gsurf(index)%axisn(i)*zeta(j))
     gsurf(index)%Za(j) = gsurf(index)%Za(j) + &
             gsurf(index)%axiszns(i)*sin(gsurf(index)%axisn(i)*zeta(j))
     Radot(j) = Radot(j) - &
             gsurf(index)%axisn(i)*gsurf(index)%axisrnc(i)*sin(gsurf(index)%axisn(i)*zeta(j))
     Zadot(j) = Zadot(j) + &
             gsurf(index)%axisn(i)*gsurf(index)%axiszns(i)*cos(gsurf(index)%axisn(i)*zeta(j))
  enddo
  ! Transfrom field line from s,zeta,theta to x,y,z
  gsurf(index)%ox(j) = cos(zeta(j))*(gsurf(index)%Ra(j) + &
          gsurf(index)%os(j)*cos(gsurf(index)%otheta(j)))
  gsurf(index)%xx(j) = cos(zeta(j))*(gsurf(index)%Ra(j) + &
          gsurf(index)%xs(j)*cos(gsurf(index)%xtheta(j)))
  gsurf(index)%oy(j) = sin(zeta(j))*(gsurf(index)%Ra(j) + &
          gsurf(index)%os(j)*cos(gsurf(index)%otheta(j)))
  gsurf(index)%xy(j) = sin(zeta(j))*(gsurf(index)%Ra(j) + &
          gsurf(index)%xs(j)*cos(gsurf(index)%xtheta(j)))
  gsurf(index)%oz(j) = gsurf(index)%os(j)*sin(gsurf(index)%otheta(j)) + &
          gsurf(index)%Za(j)
  gsurf(index)%xz(j) = gsurf(index)%xs(j)*sin(gsurf(index)%xtheta(j)) + &
          gsurf(index)%Za(j)
  ! Calculate derivative of field line w.r.t zeta
  gsurf(index)%oxdot(j) = -1.0*sin(zeta(j))*(gsurf(index)%Ra(j)&
         + gsurf(index)%os(j)*cos(gsurf(index)%otheta(j))) &
         + cos(zeta(j))*( Radot(j) + gsurf(index)%osdot(j)&
         *cos(gsurf(index)%otheta(j)) - gsurf(index)%os(j)*&
         sin(gsurf(index)%otheta(j))*gsurf(index)%othetadot(j) )
  gsurf(index)%xxdot(j) = -1.0*sin(zeta(j))*(gsurf(index)%Ra(j)& 
         + gsurf(index)%xs(j)*cos(gsurf(index)%xtheta(j))) &
         + cos(zeta(j))*( Radot(j) + gsurf(index)%xsdot(j)&
         *cos(gsurf(index)%xtheta(j)) - gsurf(index)%xs(j)*&
         sin(gsurf(index)%xtheta(j))*gsurf(index)%xthetadot(j) )
  gsurf(index)%oydot(j) =      cos(zeta(j))*(gsurf(index)%Ra(j)& 
         + gsurf(index)%os(j)*cos(gsurf(index)%otheta(j))) &
         + sin(zeta(j))*( Radot(j) + gsurf(index)%osdot(j)&
         *cos(gsurf(index)%otheta(j)) - gsurf(index)%os(j)*&
         sin(gsurf(index)%otheta(j))*gsurf(index)%othetadot(j) )
  gsurf(index)%xydot(j) =      cos(zeta(j))*(gsurf(index)%Ra(j)&
         + gsurf(index)%xs(j)*cos(gsurf(index)%xtheta(j))) &
         + sin(zeta(j))*( Radot(j) + gsurf(index)%xsdot(j)&
         *cos(gsurf(index)%xtheta(j)) - gsurf(index)%xs(j)*&
         sin(gsurf(index)%xtheta(j))*gsurf(index)%xthetadot(j) )
  gsurf(index)%ozdot(j) = gsurf(index)%osdot(j)*&
          sin(gsurf(index)%otheta(j)) + gsurf(index)%os(j)&
          *cos(gsurf(index)%otheta(j))*gsurf(index)%othetadot(j) + Zadot(j)
  gsurf(index)%xzdot(j) = gsurf(index)%xsdot(j)*&
          sin(gsurf(index)%xtheta(j)) + gsurf(index)%xs(j)&
          *cos(gsurf(index)%xtheta(j))*gsurf(index)%xthetadot(j) + Zadot(j)
  ! Calculate derivative of axis w.r.t x and y, variables are local
  do i = 1, gsurf(index)%NF_axis
     odRadx(j) = odRadx(j) + &
             gsurf(index)%axisn(i)*gsurf(index)%axisrnc(i)*sin(gsurf(index)%axisn(i)*zeta(j))*gsurf(index)%oy(j) / &
             ( gsurf(index)%ox(j)**2 * ( gsurf(index)%oy(j)**2/gsurf(index)%ox(j)**2 + 1 ) )
     odRady(j) = odRady(j) - &
             gsurf(index)%axisn(i)*gsurf(index)%axisrnc(i)*sin(gsurf(index)%axisn(i)*zeta(j)) / &
             ( gsurf(index)%ox(j) * ( gsurf(index)%oy(j)**2/gsurf(index)%ox(j)**2 + 1 ) )
     odZadx(j) = odZadx(j) - &
             gsurf(index)%axisn(i)*gsurf(index)%axiszns(i)*cos(gsurf(index)%axisn(i)*zeta(j))*gsurf(index)%oy(j) / &
             ( gsurf(index)%ox(j)**2 * ( gsurf(index)%oy(j)**2/gsurf(index)%ox(j)**2 + 1 ) )
     odZady(j) = odZady(j) + &
             gsurf(index)%axisn(i)*gsurf(index)%axiszns(i)*cos(gsurf(index)%axisn(i)*zeta(j)) / &
             ( gsurf(index)%ox(j) * ( gsurf(index)%oy(j)**2/gsurf(index)%ox(j)**2 + 1 ) )
     xdRadx(j) = xdRadx(j) + &
             gsurf(index)%axisn(i)*gsurf(index)%axisrnc(i)*sin(gsurf(index)%axisn(i)*zeta(j))*gsurf(index)%xy(j) / &
             ( gsurf(index)%xx(j)**2 * ( gsurf(index)%xy(j)**2/gsurf(index)%xx(j)**2 + 1 ) )
     xdRady(j) = xdRady(j) - &
             gsurf(index)%axisn(i)*gsurf(index)%axisrnc(i)*sin(gsurf(index)%axisn(i)*zeta(j)) / &
             ( gsurf(index)%xx(j) * ( gsurf(index)%xy(j)**2/gsurf(index)%xx(j)**2 + 1 ) )
     xdZadx(j) = xdZadx(j) - &
             gsurf(index)%axisn(i)*gsurf(index)%axiszns(i)*cos(gsurf(index)%axisn(i)*zeta(j))*gsurf(index)%xy(j) / &
             ( gsurf(index)%xx(j)**2 * ( gsurf(index)%xy(j)**2/gsurf(index)%xx(j)**2 + 1 ) )
     xdZady(j) = xdZady(j) + &
             gsurf(index)%axisn(i)*gsurf(index)%axiszns(i)*cos(gsurf(index)%axisn(i)*zeta(j)) / &
             ( gsurf(index)%xx(j) * ( gsurf(index)%xy(j)**2/gsurf(index)%xx(j)**2 + 1 ) )
  enddo 
  ! Calculate contravariant basis, variables are local
  gradosx(j) = ( sqrt( gsurf(index)%ox(j)**2 + gsurf(index)%oy(j)**2 ) - &
          gsurf(index)%Ra(j) )*( (gsurf(index)%ox(j)**2 + &
          gsurf(index)%oy(j)**2)**(-.5)*gsurf(index)%ox(j) - &
          odRadx(j) )/gsurf(index)%os(j) - &
          (gsurf(index)%oz(j) - gsurf(index)%Za(j))*odZadx(j)/gsurf(index)%os(j)
  gradxsx(j) = ( sqrt( gsurf(index)%xx(j)**2 + gsurf(index)%xy(j)**2 ) - &
          gsurf(index)%Ra(j) )*( (gsurf(index)%xx(j)**2 + &
          gsurf(index)%xy(j)**2)**(-.5)*gsurf(index)%xx(j) - & 
          xdRadx(j) )/gsurf(index)%xs(j) - &
          (gsurf(index)%xz(j) - gsurf(index)%Za(j))*xdZadx(j)/gsurf(index)%xs(j)
  gradosy(j) = ( sqrt( gsurf(index)%ox(j)**2 + gsurf(index)%oy(j)**2 ) - &
          gsurf(index)%Ra(j) )*( (gsurf(index)%ox(j)**2 + &
          gsurf(index)%oy(j)**2)**(-.5)*gsurf(index)%oy(j) - &
          odRady(j) )/gsurf(index)%os(j) - &
          (gsurf(index)%oz(j) - gsurf(index)%Za(j))*odZady(j)/gsurf(index)%os(j)
  gradxsy(j) = ( sqrt( gsurf(index)%xx(j)**2 + gsurf(index)%xy(j)**2 ) - &
          gsurf(index)%Ra(j) )*( (gsurf(index)%xx(j)**2 + &
          gsurf(index)%xy(j)**2)**(-.5)*gsurf(index)%xy(j) - &
          xdRady(j) )/gsurf(index)%xs(j) - &
          (gsurf(index)%xz(j) - gsurf(index)%Za(j))*xdZady(j)/gsurf(index)%xs(j)
  gradosz(j) = ( gsurf(index)%oz(j) - gsurf(index)%Za(j) )/gsurf(index)%os(j)
  gradxsz(j) = ( gsurf(index)%xz(j) - gsurf(index)%Za(j) )/gsurf(index)%xs(j)
  gradozetax(j) = -1.0*(gsurf(index)%oy(j)**2/gsurf(index)%ox(j)**2+1)**-1.0 * &
          gsurf(index)%oy(j)/gsurf(index)%ox(j)**2
  gradxzetax(j) = -1.0*(gsurf(index)%xy(j)**2/gsurf(index)%xx(j)**2+1)**-1.0 * &
          gsurf(index)%xy(j)/gsurf(index)%xx(j)**2
  gradozetay(j) = (gsurf(index)%oy(j)**2/gsurf(index)%ox(j)**2+1)**-1.0 * &
          1.0/gsurf(index)%ox(j)
  gradxzetay(j) = (gsurf(index)%xy(j)**2/gsurf(index)%xx(j)**2+1)**-1.0 * &
          1.0/gsurf(index)%xx(j)
  !gradozetaz(j) = 0.0
  !gradxzetaz(j) = 0.0
  gradothetax(j) = -1.0*( (gsurf(index)%oz(j)-gsurf(index)%Za(j)) * &
          (sqrt(gsurf(index)%ox(j)**2+gsurf(index)%oy(j)**2) - gsurf(index)%Ra(j))**-2.0 * &
          ( (gsurf(index)%ox(j)**2+gsurf(index)%oy(j)**2)**-.5 * gsurf(index)%ox(j) - odRadx(j) ) &
          + odZadx(j) / (sqrt(gsurf(index)%ox(j)**2+gsurf(index)%oy(j)**2) &
          - gsurf(index)%Ra(j)) ) / (tan(gsurf(index)%otheta(j))**2+1)
  gradxthetax(j) = -1.0*( (gsurf(index)%xz(j)-gsurf(index)%Za(j)) * &
          (sqrt(gsurf(index)%xx(j)**2+gsurf(index)%xy(j)**2) - gsurf(index)%Ra(j))**-2.0 * &
          ( (gsurf(index)%xx(j)**2+gsurf(index)%xy(j)**2)**-.5 * gsurf(index)%xx(j) - xdRadx(j) ) &
          + xdZadx(j) / (sqrt(gsurf(index)%xx(j)**2+gsurf(index)%xy(j)**2) &
          - gsurf(index)%Ra(j)) ) / (tan(gsurf(index)%xtheta(j))**2+1)
  gradothetay(j) = -1.0*( (gsurf(index)%oz(j)-gsurf(index)%Za(j)) * &
          (sqrt(gsurf(index)%ox(j)**2+gsurf(index)%oy(j)**2) - gsurf(index)%Ra(j))**-2.0 * &
          ( (gsurf(index)%ox(j)**2+gsurf(index)%oy(j)**2)**-.5 * gsurf(index)%oy(j) - odRady(j) ) &
          + odZady(j) / (sqrt(gsurf(index)%ox(j)**2+gsurf(index)%oy(j)**2) &
          - gsurf(index)%Ra(j)) ) / (tan(gsurf(index)%otheta(j))**2+1)
  gradxthetay(j) = -1.0*( (gsurf(index)%xz(j)-gsurf(index)%Za(j)) * &
          (sqrt(gsurf(index)%xx(j)**2+gsurf(index)%xy(j)**2) - gsurf(index)%Ra(j))**-2.0 * &
          ( (gsurf(index)%xx(j)**2+gsurf(index)%xy(j)**2)**-.5 * gsurf(index)%xy(j) - xdRady(j) ) &
          + xdZady(j) / (sqrt(gsurf(index)%xx(j)**2+gsurf(index)%xy(j)**2) &
          - gsurf(index)%Ra(j)) ) / (tan(gsurf(index)%xtheta(j))**2+1)
  gradothetaz(j) = ( (tan(gsurf(index)%otheta(j))**2+1) * &
          ( sqrt(gsurf(index)%ox(j)**2+gsurf(index)%oy(j)**2) - gsurf(index)%Ra(j) ) )**-1.0
  gradxthetaz(j) = ( (tan(gsurf(index)%xtheta(j))**2+1) * &
          ( sqrt(gsurf(index)%xx(j)**2+gsurf(index)%xy(j)**2) - gsurf(index)%Ra(j) ) )**-1.0
  ! Calculate field on stable field lines
  do icoil = 1, Ncoils
     call bfield0(icoil, gsurf(index)%ox(j), gsurf(index)%oy(j), gsurf(index)%oz(j), Bxhold, Byhold, Bzhold)
     oBx(j) = oBx(j) + Bxhold
     oBy(j) = oBy(j) + Byhold
     oBz(j) = oBz(j) + Bzhold
     call bfield0(icoil, gsurf(index)%xx(i), gsurf(index)%xy(i), gsurf(index)%xz(i), Bxhold, Byhold, Bzhold)
     xBx(j) = xBx(j) + Bxhold
     xBy(j) = xBy(j) + Byhold
     xBz(j) = xBz(j) + Bzhold
  enddo
  ! Calculate B^s, B^zeta, B^theta on stable field lines
  gsurf(index)%obsups(j) = oBx(j)*gradosx(j) + &
          oBy(j)*gradosy(j) + oBz(j)*gradosz(j)
  gsurf(index)%xbsups(j) = xBx(j)*gradxsx(j) + &
          xBy(j)*gradxsy(j) + xBz(j)*gradxsz(j)
  gsurf(index)%obsupzeta(j) = oBx(j)*gradozetax(j) + &
          oBy(j)*gradozetay(j) + oBz(j)*gradozetaz(j)
  gsurf(index)%xbsupzeta(j) = xBx(j)*gradxzetax(j) + &
          xBy(j)*gradxzetay(j) + xBz(j)*gradxzetaz(j)
  gsurf(index)%obsuptheta(j) = oBx(j)*gradothetax(j) + &
          oBy(j)*gradothetay(j) + oBz(j)*gradothetaz(j)
  gsurf(index)%xbsuptheta(j) = xBx(j)*gradxthetax(j) + &
          xBy(j)*gradxthetay(j) + xBz(j)*gradxthetaz(j)
  ! Calculate f,g and F
  gsurf(index)%of(j) = gsurf(index)%obsups(j)/gsurf(index)%obsupzeta(j) - &
          gsurf(index)%osdot(j)
  gsurf(index)%xf(j) = gsurf(index)%xbsups(j)/gsurf(index)%xbsupzeta(j) - &
          gsurf(index)%xsdot(j)
  gsurf(index)%og(j) = gsurf(index)%obsuptheta(j)/gsurf(index)%obsupzeta(j) - &
          gsurf(index)%othetadot(j)
  gsurf(index)%xg(j) = gsurf(index)%xbsuptheta(j)/gsurf(index)%xbsupzeta(j) - &
          gsurf(index)%xthetadot(j)

  enddo ! end of unindented big loop

  call MPI_BARRIER( MPI_COMM_FOCUS, ierr )
  call MPI_ALLREDUCE( MPI_IN_PLACE, gsurf(index)%of, Nseg_stable, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FOCUS, ierr )
  call MPI_ALLREDUCE( MPI_IN_PLACE, gsurf(index)%xf, Nseg_stable, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FOCUS, ierr )
  call MPI_ALLREDUCE( MPI_IN_PLACE, gsurf(index)%og, Nseg_stable, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FOCUS, ierr )
  call MPI_ALLREDUCE( MPI_IN_PLACE, gsurf(index)%xg, Nseg_stable, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FOCUS, ierr )

  gsurf(index)%F = 0.0
  do i = j, Nseg_stable-1
     gsurf(index)%F = gsurf(index)%F + gsurf(index)%of(j)**2 + gsurf(index)%og(j)**2 + &
             gsurf(index)%xf(j)**2 + gsurf(index)%xg(j)**2
  enddo
  gsurf(index)%F = pi2*q*gsurf(index)%F/(Nseg_stable-1)

  ! DALLOCATE
  DALLOCATE( zeta )
  DALLOCATE( Radot )
  DALLOCATE( Zadot )
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
                               resbn_m, pi2, psmall, gsurf
  use mpi
  implicit none

  INTEGER, INTENT(in)  :: index
  
  INTEGER              :: astat, ierr, idof, NF_stable, Nseg_stable, i, j
  REAL                 :: tmp_xdof(1:gsurf(index)%Ndof_stable), fd, negvalue, posvalue
  REAL                 :: start, finish, small, q
  REAL, ALLOCATABLE    :: zeta(:), gradobsupsx(:), gradobsupsy(:), gradobsupsz(:), gradobsupzetax(:), &
  gradobsupzetay(:), gradobsupzetaz(:), gradobsupthetax(:), gradobsupthetay(:), gradobsupthetaz(:), &
  gradxbsupsx(:), gradxbsupsy(:), gradxbsupsz(:), gradxbsupzetax(:), gradxbsupzetay(:), gradxbsupzetaz(:), &
  gradxbsupthetax(:), gradxbsupthetay(:), gradxbsupthetaz(:), oxneg(:), oyneg(:), ozneg(:), oxpos(:), &
  oypos(:), ozpos(:), xxneg(:), xyneg(:), xzneg(:), xxpos(:), xypos(:), xzpos(:), negbsups(:), negbsupzeta(:), &
  negbsuptheta(:), posbsups(:), posbsupzeta(:), posbsuptheta(:), dosdosnc(:,:), dosdosns(:,:), &
  dothetadothetanc(:,:), dothetadothetans(:,:), dxsdxsnc(:,:), dxsdxsns(:,:), dxthetadxthetanc(:,:), &
  dxthetadxthetans(:,:), dosdotdosnc(:,:), dosdotdosns(:,:), dothetadotdothetanc(:,:), dothetadotdothetans(:,:), &
  dxsdotdxsnc(:,:), dxsdotdxsns(:,:), dxthetadotdxthetanc(:,:), dxthetadotdxthetans(:,:), doxdosnc(:,:), &
  doxdosns(:,:), doydosnc(:,:), doydosns(:,:), dozdosnc(:,:), dozdosns(:,:), dxxdxsnc(:,:), dxxdxsns(:,:), &
  dxydxsnc(:,:), dxydxsns(:,:), dxzdxsnc(:,:), dxzdxsns(:,:), doxdothetanc(:,:), doxdothetans(:,:), & 
  doydothetanc(:,:), doydothetans(:,:), dozdothetanc(:,:), dozdothetans(:,:), dxxdxthetanc(:,:), &
  dxxdxthetans(:,:), dxydxthetanc(:,:), dxydxthetans(:,:), dxzdxthetanc(:,:), dxzdxthetans(:,:), &
  dobsupsdosnc(:,:), dobsupsdosns(:,:), dobsupzetadosnc(:,:), dobsupzetadosns(:,:), dobsupthetadosnc(:,:), &
  dobsupthetadosns(:,:), dxbsupsdxsnc(:,:), dxbsupsdxsns(:,:), dxbsupzetadxsnc(:,:), dxbsupzetadxsns(:,:), &
  dxbsupthetadxsnc(:,:), dxbsupthetadxsns(:,:), dobsupsdothetanc(:,:), dobsupsdothetans(:,:), &
  dobsupzetadothetanc(:,:), dobsupzetadothetans(:,:), dobsupthetadothetanc(:,:), dobsupthetadothetans(:,:), &
  dxbsupsdxthetanc(:,:), dxbsupsdxthetans(:,:), dxbsupzetadxthetanc(:,:), dxbsupzetadxthetans(:,:), &
  dxbsupthetadxthetanc(:,:), dxbsupthetadxthetans(:,:), dofdosnc(:,:), dofdosns(:,:), dofdothetanc(:,:), &
  dofdothetans(:,:), dogdosnc(:,:), dogdosns(:,:), dogdothetanc(:,:), dogdothetans(:,:), dxfdxsnc(:,:), &
  dxfdxsns(:,:), dxfdxthetanc(:,:), dxfdxthetans(:,:), dxgdxsnc(:,:), dxgdxsns(:,:), dxgdxthetanc(:,:), &
  dxgdxthetans(:,:), dFdosnc(:), dFdosns(:), dFdothetanc(:), dFdothetans(:), dFdxsnc(:), dFdxsns(:), &
  dFdxthetanc(:), dFdxthetans(:), oxsmall(:), oysmall(:), ozsmall(:), xxsmall(:), xysmall(:), xzsmall(:)

  ! figure out what dimensions are needed

  !do idof = 1, gsurf(index)%Ndof_stable
  !   ! perturbation will be relative.
  !   small = gsurf(index)%xdof_stable(idof) * psmall
  !   if (abs(small)<machprec) small = psmall
  !   !backward pertubation;
  !   tmp_xdof = gsurf(index)%xdof_stable
  !   tmp_xdof(idof) = tmp_xdof(idof) - half * small
  !   call unpacking_stable(tmp_xdof,index)
  !   call calcfg(index)
  !   negvalue = gsurf(index)%F
  !   !forward pertubation;
  !   tmp_xdof = gsurf(index)%xdof_stable
  !   tmp_xdof(idof) = tmp_xdof(idof) + half * small
  !   call unpacking_stable(tmp_xdof,index)
  !   call calcfg(index)
  !   posvalue = gsurf(index)%F
  !   !finite difference;
  !   gsurf(index)%dFdxdof_stable(idof) = (posvalue - negvalue) / small
  !enddo
  
  q = real(resbn_m)

  NF_stable = gsurf(index)%NF_stable
  Nseg_stable = gsurf(index)%Nseg_stable

  SALLOCATE( zeta, (1:Nseg_stable), 0.0 )
  zeta(1:Nseg_stable) = gsurf(index)%zeta(1:Nseg_stable)

  ! Consider deallocating arrays sooner
  
  SALLOCATE( gradobsupsx, (1:Nseg_stable), 0.0 )
  SALLOCATE( gradobsupsy, (1:Nseg_stable), 0.0 )
  SALLOCATE( gradobsupsz, (1:Nseg_stable), 0.0 )
  SALLOCATE( gradobsupzetax, (1:Nseg_stable), 0.0 )
  SALLOCATE( gradobsupzetay, (1:Nseg_stable), 0.0 )
  SALLOCATE( gradobsupzetaz, (1:Nseg_stable), 0.0 )
  SALLOCATE( gradobsupthetax, (1:Nseg_stable), 0.0 )
  SALLOCATE( gradobsupthetay, (1:Nseg_stable), 0.0 )
  SALLOCATE( gradobsupthetaz, (1:Nseg_stable), 0.0 )

  SALLOCATE( gradxbsupsx, (1:Nseg_stable), 0.0 )
  SALLOCATE( gradxbsupsy, (1:Nseg_stable), 0.0 )
  SALLOCATE( gradxbsupsz, (1:Nseg_stable), 0.0 )
  SALLOCATE( gradxbsupzetax, (1:Nseg_stable), 0.0 )
  SALLOCATE( gradxbsupzetay, (1:Nseg_stable), 0.0 )
  SALLOCATE( gradxbsupzetaz, (1:Nseg_stable), 0.0 )
  SALLOCATE( gradxbsupthetax, (1:Nseg_stable), 0.0 )
  SALLOCATE( gradxbsupthetay, (1:Nseg_stable), 0.0 )
  SALLOCATE( gradxbsupthetaz, (1:Nseg_stable), 0.0 )

  SALLOCATE( oxneg, (1:Nseg_stable), 0.0 )
  SALLOCATE( oyneg, (1:Nseg_stable), 0.0 )
  SALLOCATE( ozneg, (1:Nseg_stable), 0.0 )
  SALLOCATE( oxpos, (1:Nseg_stable), 0.0 )
  SALLOCATE( oypos, (1:Nseg_stable), 0.0 )
  SALLOCATE( ozpos, (1:Nseg_stable), 0.0 )

  SALLOCATE( xxneg, (1:Nseg_stable), 0.0 )
  SALLOCATE( xyneg, (1:Nseg_stable), 0.0 )
  SALLOCATE( xzneg, (1:Nseg_stable), 0.0 )
  SALLOCATE( xxpos, (1:Nseg_stable), 0.0 )
  SALLOCATE( xypos, (1:Nseg_stable), 0.0 )
  SALLOCATE( xzpos, (1:Nseg_stable), 0.0 )

  SALLOCATE( negbsups, (1:Nseg_stable), 0.0 )
  SALLOCATE( negbsupzeta, (1:Nseg_stable), 0.0 )
  SALLOCATE( negbsuptheta, (1:Nseg_stable), 0.0 )
  SALLOCATE( posbsups, (1:Nseg_stable), 0.0 )
  SALLOCATE( posbsupzeta, (1:Nseg_stable), 0.0 )
  SALLOCATE( posbsuptheta, (1:Nseg_stable), 0.0 )

  SALLOCATE( oxsmall, (1:Nseg_stable), 0.0 )
  SALLOCATE( oysmall, (1:Nseg_stable), 0.0 )
  SALLOCATE( ozsmall, (1:Nseg_stable), 0.0 )
  SALLOCATE( xxsmall, (1:Nseg_stable), 0.0 )
  SALLOCATE( xysmall, (1:Nseg_stable), 0.0 )
  SALLOCATE( xzsmall, (1:Nseg_stable), 0.0 )

  ! Functions to take FD for gradient of contravariant field 
  oxsmall(1:Nseg_stable) = gsurf(index)%ox(1:Nseg_stable) * psmall
  oysmall(1:Nseg_stable) = gsurf(index)%oy(1:Nseg_stable) * psmall
  ozsmall(1:Nseg_stable) = gsurf(index)%oz(1:Nseg_stable) * psmall
  do j = 1, Nseg_stable
     if (abs(oxsmall(j))<machprec) oxsmall(j) = psmall
     if (abs(oysmall(j))<machprec) oysmall(j) = psmall
     if (abs(ozsmall(j))<machprec) ozsmall(j) = psmall
  enddo

  oxneg(1:Nseg_stable) = gsurf(index)%ox(1:Nseg_stable) - half*oxsmall(1:Nseg_stable)
  oyneg(1:Nseg_stable) = gsurf(index)%oy(1:Nseg_stable) - half*oysmall(1:Nseg_stable)
  ozneg(1:Nseg_stable) = gsurf(index)%oz(1:Nseg_stable) - half*ozsmall(1:Nseg_stable)
  oxpos(1:Nseg_stable) = gsurf(index)%ox(1:Nseg_stable) + half*oxsmall(1:Nseg_stable)
  oypos(1:Nseg_stable) = gsurf(index)%oy(1:Nseg_stable) + half*oysmall(1:Nseg_stable)
  ozpos(1:Nseg_stable) = gsurf(index)%oz(1:Nseg_stable) + half*ozsmall(1:Nseg_stable)
     
  call calcbsup(index,oxneg,gsurf(index)%oy,gsurf(index)%oz,negbsups,negbsupzeta,negbsuptheta)
  call calcbsup(index,oxpos,gsurf(index)%oy,gsurf(index)%oz,posbsups,posbsupzeta,posbsuptheta)
  gradobsupsx(1:Nseg_stable) = (posbsups(1:Nseg_stable) - negbsups(1:Nseg_stable)) / oxsmall(1:Nseg_stable)
  gradobsupzetax(1:Nseg_stable) = (posbsupzeta(1:Nseg_stable) - negbsupzeta(1:Nseg_stable)) / oxsmall(1:Nseg_stable)
  gradobsupthetax(1:Nseg_stable) = (posbsuptheta(1:Nseg_stable) - negbsuptheta(1:Nseg_stable)) / oxsmall(1:Nseg_stable)

  call calcbsup(index,gsurf(index)%ox,oyneg,gsurf(index)%oz,negbsups,negbsupzeta,negbsuptheta)
  call calcbsup(index,gsurf(index)%ox,oypos,gsurf(index)%oz,posbsups,posbsupzeta,posbsuptheta)
  gradobsupsy(1:Nseg_stable) = (posbsups(1:Nseg_stable) - negbsups(1:Nseg_stable)) / oysmall(1:Nseg_stable)
  gradobsupzetay(1:Nseg_stable) = (posbsupzeta(1:Nseg_stable) - negbsupzeta(1:Nseg_stable)) / oysmall(1:Nseg_stable)
  gradobsupthetay(1:Nseg_stable) = (posbsuptheta(1:Nseg_stable) - negbsuptheta(1:Nseg_stable)) / oysmall(1:Nseg_stable) 
    
  call calcbsup(index,gsurf(index)%ox,gsurf(index)%oy,ozneg,negbsups,negbsupzeta,negbsuptheta)
  call calcbsup(index,gsurf(index)%ox,gsurf(index)%oy,ozpos,posbsups,posbsupzeta,posbsuptheta)
  gradobsupsz(1:Nseg_stable) = (posbsups(1:Nseg_stable) - negbsups(1:Nseg_stable)) / ozsmall(1:Nseg_stable)
  gradobsupzetaz(1:Nseg_stable) = (posbsupzeta(1:Nseg_stable) - negbsupzeta(1:Nseg_stable)) / ozsmall(1:Nseg_stable)
  gradobsupthetaz(1:Nseg_stable) = (posbsuptheta(1:Nseg_stable) - negbsuptheta(1:Nseg_stable)) / ozsmall(1:Nseg_stable)

  ! Everything for x

  xxsmall(1:Nseg_stable) = gsurf(index)%xx(1:Nseg_stable) * psmall
  xysmall(1:Nseg_stable) = gsurf(index)%xy(1:Nseg_stable) * psmall
  xzsmall(1:Nseg_stable) = gsurf(index)%xz(1:Nseg_stable) * psmall
  do j = 1, Nseg_stable
     if (abs(xxsmall(j))<machprec) xxsmall(j) = psmall
     if (abs(xysmall(j))<machprec) xysmall(j) = psmall
     if (abs(xzsmall(j))<machprec) xzsmall(j) = psmall
  enddo

  xxneg(1:Nseg_stable) = gsurf(index)%xx(1:Nseg_stable) - half*xxsmall(1:Nseg_stable)
  xyneg(1:Nseg_stable) = gsurf(index)%xy(1:Nseg_stable) - half*xysmall(1:Nseg_stable)
  xzneg(1:Nseg_stable) = gsurf(index)%xz(1:Nseg_stable) - half*xzsmall(1:Nseg_stable)
  xxpos(1:Nseg_stable) = gsurf(index)%xx(1:Nseg_stable) + half*xxsmall(1:Nseg_stable)
  xypos(1:Nseg_stable) = gsurf(index)%xy(1:Nseg_stable) + half*xysmall(1:Nseg_stable)
  xzpos(1:Nseg_stable) = gsurf(index)%xz(1:Nseg_stable) + half*xzsmall(1:Nseg_stable)

  call calcbsup(index,xxneg,gsurf(index)%xy,gsurf(index)%xz,negbsups,negbsupzeta,negbsuptheta)
  call calcbsup(index,xxpos,gsurf(index)%xy,gsurf(index)%xz,posbsups,posbsupzeta,posbsuptheta)
  gradxbsupsx(1:Nseg_stable) = (posbsups(1:Nseg_stable) - negbsups(1:Nseg_stable)) / xxsmall(1:Nseg_stable)
  gradxbsupzetax(1:Nseg_stable) = (posbsupzeta(1:Nseg_stable) - negbsupzeta(1:Nseg_stable)) / xxsmall(1:Nseg_stable)
  gradxbsupthetax(1:Nseg_stable) = (posbsuptheta(1:Nseg_stable) - negbsuptheta(1:Nseg_stable)) / xxsmall(1:Nseg_stable)

  call calcbsup(index,gsurf(index)%xx,xyneg,gsurf(index)%xz,negbsups,negbsupzeta,negbsuptheta)
  call calcbsup(index,gsurf(index)%xx,xypos,gsurf(index)%xz,posbsups,posbsupzeta,posbsuptheta)
  gradxbsupsy(1:Nseg_stable) = (posbsups(1:Nseg_stable) - negbsups(1:Nseg_stable)) / xysmall(1:Nseg_stable)
  gradxbsupzetay(1:Nseg_stable) = (posbsupzeta(1:Nseg_stable) - negbsupzeta(1:Nseg_stable)) / xysmall(1:Nseg_stable)
  gradxbsupthetay(1:Nseg_stable) = (posbsuptheta(1:Nseg_stable) - negbsuptheta(1:Nseg_stable)) / xysmall(1:Nseg_stable)   

  call calcbsup(index,gsurf(index)%xx,gsurf(index)%xy,xzneg,negbsups,negbsupzeta,negbsuptheta)
  call calcbsup(index,gsurf(index)%xx,gsurf(index)%xy,xzpos,posbsups,posbsupzeta,posbsuptheta)
  gradxbsupsz(1:Nseg_stable) = (posbsups(1:Nseg_stable) - negbsups(1:Nseg_stable)) / xzsmall(1:Nseg_stable)
  gradxbsupzetaz(1:Nseg_stable) = (posbsupzeta(1:Nseg_stable) - negbsupzeta(1:Nseg_stable)) / xzsmall(1:Nseg_stable)
  gradxbsupthetaz(1:Nseg_stable) = (posbsuptheta(1:Nseg_stable) - negbsuptheta(1:Nseg_stable)) / xzsmall(1:Nseg_stable) 

  SALLOCATE( dosdosnc, (1:Nseg_stable,1:NF_stable), 0.0 )
  SALLOCATE( dosdosns, (1:Nseg_stable,1:NF_stable), 0.0 )
  SALLOCATE( dothetadothetanc, (1:Nseg_stable,1:NF_stable), 0.0 )
  SALLOCATE( dothetadothetans, (1:Nseg_stable,1:NF_stable), 0.0 )
  SALLOCATE( dxsdxsnc, (1:Nseg_stable,1:NF_stable), 0.0 )
  SALLOCATE( dxsdxsns, (1:Nseg_stable,1:NF_stable), 0.0 )
  SALLOCATE( dxthetadxthetanc, (1:Nseg_stable,1:NF_stable), 0.0 )
  SALLOCATE( dxthetadxthetans, (1:Nseg_stable,1:NF_stable), 0.0 )
  
  SALLOCATE( dosdotdosnc, (1:Nseg_stable,1:NF_stable), 0.0 )
  SALLOCATE( dosdotdosns, (1:Nseg_stable,1:NF_stable), 0.0 )
  SALLOCATE( dothetadotdothetanc, (1:Nseg_stable,1:NF_stable), 0.0 )
  SALLOCATE( dothetadotdothetans, (1:Nseg_stable,1:NF_stable), 0.0 )
  SALLOCATE( dxsdotdxsnc, (1:Nseg_stable,1:NF_stable), 0.0 )
  SALLOCATE( dxsdotdxsns, (1:Nseg_stable,1:NF_stable), 0.0 )
  SALLOCATE( dxthetadotdxthetanc, (1:Nseg_stable,1:NF_stable), 0.0 )
  SALLOCATE( dxthetadotdxthetans, (1:Nseg_stable,1:NF_stable), 0.0 )
  
  do i = 1, NF_stable
     dosdosnc(1:Nseg_stable,i) = cos(gsurf(index)%on(i)*zeta(1:Nseg_stable)/q)
     dosdosns(1:Nseg_stable,i) = sin(gsurf(index)%on(i)*zeta(1:Nseg_stable)/q)
     dothetadothetanc(1:Nseg_stable,i) = cos(gsurf(index)%on(i)*zeta(1:Nseg_stable)/q)
     dothetadothetans(1:Nseg_stable,i) = sin(gsurf(index)%on(i)*zeta(1:Nseg_stable)/q)
     dxsdxsnc(1:Nseg_stable,i) = cos(gsurf(index)%xn(i)*zeta(1:Nseg_stable)/q)
     dxsdxsns(1:Nseg_stable,i) = sin(gsurf(index)%xn(i)*zeta(1:Nseg_stable)/q)
     dxthetadxthetanc(1:Nseg_stable,i) = cos(gsurf(index)%xn(i)*zeta(1:Nseg_stable)/q)
     dxthetadxthetans(1:Nseg_stable,i) = sin(gsurf(index)%xn(i)*zeta(1:Nseg_stable)/q)

     dosdotdosnc(1:Nseg_stable,i)         = -1.0*gsurf(index)%on(i)*sin(gsurf(index)%on(i)*zeta(1:Nseg_stable)/q)/q
     dosdotdosns(1:Nseg_stable,i)         =      gsurf(index)%on(i)*cos(gsurf(index)%on(i)*zeta(1:Nseg_stable)/q)/q
     dothetadotdothetanc(1:Nseg_stable,i) = -1.0*gsurf(index)%on(i)*sin(gsurf(index)%on(i)*zeta(1:Nseg_stable)/q)/q
     dothetadotdothetans(1:Nseg_stable,i) =      gsurf(index)%on(i)*cos(gsurf(index)%on(i)*zeta(1:Nseg_stable)/q)/q
     dxsdotdxsnc(1:Nseg_stable,i)         = -1.0*gsurf(index)%xn(i)*sin(gsurf(index)%xn(i)*zeta(1:Nseg_stable)/q)/q
     dxsdotdxsns(1:Nseg_stable,i)         =      gsurf(index)%xn(i)*cos(gsurf(index)%xn(i)*zeta(1:Nseg_stable)/q)/q
     dxthetadotdxthetanc(1:Nseg_stable,i) = -1.0*gsurf(index)%xn(i)*sin(gsurf(index)%xn(i)*zeta(1:Nseg_stable)/q)/q 
     dxthetadotdxthetans(1:Nseg_stable,i) =      gsurf(index)%xn(i)*cos(gsurf(index)%xn(i)*zeta(1:Nseg_stable)/q)/q
  enddo

  SALLOCATE( doxdosnc, (1:Nseg_stable,1:NF_stable), 0.0 )
  SALLOCATE( doxdosns, (1:Nseg_stable,1:NF_stable), 0.0 )
  SALLOCATE( doydosnc, (1:Nseg_stable,1:NF_stable), 0.0 )
  SALLOCATE( doydosns, (1:Nseg_stable,1:NF_stable), 0.0 )
  SALLOCATE( dozdosnc, (1:Nseg_stable,1:NF_stable), 0.0 )
  SALLOCATE( dozdosns, (1:Nseg_stable,1:NF_stable), 0.0 )

  SALLOCATE( dxxdxsnc, (1:Nseg_stable,1:NF_stable), 0.0 )
  SALLOCATE( dxxdxsns, (1:Nseg_stable,1:NF_stable), 0.0 )
  SALLOCATE( dxydxsnc, (1:Nseg_stable,1:NF_stable), 0.0 )
  SALLOCATE( dxydxsns, (1:Nseg_stable,1:NF_stable), 0.0 )
  SALLOCATE( dxzdxsnc, (1:Nseg_stable,1:NF_stable), 0.0 )
  SALLOCATE( dxzdxsns, (1:Nseg_stable,1:NF_stable), 0.0 )

  SALLOCATE( doxdothetanc, (1:Nseg_stable,1:NF_stable), 0.0 )
  SALLOCATE( doxdothetans, (1:Nseg_stable,1:NF_stable), 0.0 )
  SALLOCATE( doydothetanc, (1:Nseg_stable,1:NF_stable), 0.0 )
  SALLOCATE( doydothetans, (1:Nseg_stable,1:NF_stable), 0.0 )
  SALLOCATE( dozdothetanc, (1:Nseg_stable,1:NF_stable), 0.0 )
  SALLOCATE( dozdothetans, (1:Nseg_stable,1:NF_stable), 0.0 )

  SALLOCATE( dxxdxthetanc, (1:Nseg_stable,1:NF_stable), 0.0 )
  SALLOCATE( dxxdxthetans, (1:Nseg_stable,1:NF_stable), 0.0 )
  SALLOCATE( dxydxthetanc, (1:Nseg_stable,1:NF_stable), 0.0 )
  SALLOCATE( dxydxthetans, (1:Nseg_stable,1:NF_stable), 0.0 )
  SALLOCATE( dxzdxthetanc, (1:Nseg_stable,1:NF_stable), 0.0 )
  SALLOCATE( dxzdxthetans, (1:Nseg_stable,1:NF_stable), 0.0 )

  do i = 1, NF_stable
     doxdosnc(1:Nseg_stable,i) = cos(zeta(1:Nseg_stable))*dosdosnc(1:Nseg_stable,i)*cos(gsurf(index)%otheta(1:Nseg_stable))
     doxdosns(1:Nseg_stable,i) = cos(zeta(1:Nseg_stable))*dosdosns(1:Nseg_stable,i)*cos(gsurf(index)%otheta(1:Nseg_stable))
     doydosnc(1:Nseg_stable,i) = sin(zeta(1:Nseg_stable))*dosdosnc(1:Nseg_stable,i)*cos(gsurf(index)%otheta(1:Nseg_stable))
     doydosns(1:Nseg_stable,i) = sin(zeta(1:Nseg_stable))*dosdosns(1:Nseg_stable,i)*cos(gsurf(index)%otheta(1:Nseg_stable))
     dozdosnc(1:Nseg_stable,i) =                          dosdosnc(1:Nseg_stable,i)*sin(gsurf(index)%otheta(1:Nseg_stable))
     dozdosns(1:Nseg_stable,i) =                          dosdosns(1:Nseg_stable,i)*sin(gsurf(index)%otheta(1:Nseg_stable))

     dxxdxsnc(1:Nseg_stable,i) = cos(zeta(1:Nseg_stable))*dxsdxsnc(1:Nseg_stable,i)*cos(gsurf(index)%xtheta(1:Nseg_stable))
     dxxdxsns(1:Nseg_stable,i) = cos(zeta(1:Nseg_stable))*dxsdxsns(1:Nseg_stable,i)*cos(gsurf(index)%xtheta(1:Nseg_stable))
     dxydxsnc(1:Nseg_stable,i) = sin(zeta(1:Nseg_stable))*dxsdxsnc(1:Nseg_stable,i)*cos(gsurf(index)%xtheta(1:Nseg_stable))
     dxydxsns(1:Nseg_stable,i) = sin(zeta(1:Nseg_stable))*dxsdxsns(1:Nseg_stable,i)*cos(gsurf(index)%xtheta(1:Nseg_stable))
     dxzdxsnc(1:Nseg_stable,i) =                          dxsdxsnc(1:Nseg_stable,i)*sin(gsurf(index)%xtheta(1:Nseg_stable))
     dxzdxsns(1:Nseg_stable,i) =                          dxsdxsns(1:Nseg_stable,i)*sin(gsurf(index)%xtheta(1:Nseg_stable))

     doxdothetanc(1:Nseg_stable,i) = -1.0*cos(zeta(1:Nseg_stable))*gsurf(index)%os(1:Nseg_stable)&
             *sin(gsurf(index)%otheta(1:Nseg_stable))*dothetadothetanc(1:Nseg_stable,i)
     doxdothetans(1:Nseg_stable,i) = -1.0*cos(zeta(1:Nseg_stable))*gsurf(index)%os(1:Nseg_stable)&
             *sin(gsurf(index)%otheta(1:Nseg_stable))*dothetadothetans(1:Nseg_stable,i)
     doydothetanc(1:Nseg_stable,i) = -1.0*sin(zeta(1:Nseg_stable))*gsurf(index)%os(1:Nseg_stable)&
             *sin(gsurf(index)%otheta(1:Nseg_stable))*dothetadothetanc(1:Nseg_stable,i)
     doydothetans(1:Nseg_stable,i) = -1.0*sin(zeta(1:Nseg_stable))*gsurf(index)%os(1:Nseg_stable)&
             *sin(gsurf(index)%otheta(1:Nseg_stable))*dothetadothetans(1:Nseg_stable,i)
     dozdothetanc(1:Nseg_stable,i) =                               gsurf(index)%os(1:Nseg_stable)&
             *cos(gsurf(index)%otheta(1:Nseg_stable))*dothetadothetanc(1:Nseg_stable,i)
     dozdothetans(1:Nseg_stable,i) =                               gsurf(index)%os(1:Nseg_stable)&
             *cos(gsurf(index)%otheta(1:Nseg_stable))*dothetadothetans(1:Nseg_stable,i)

     dxxdxthetanc(1:Nseg_stable,i) = -1.0*cos(zeta(1:Nseg_stable))*gsurf(index)%xs(1:Nseg_stable)&
             *sin(gsurf(index)%xtheta(1:Nseg_stable))*dxthetadxthetanc(1:Nseg_stable,i)
     dxxdxthetans(1:Nseg_stable,i) = -1.0*cos(zeta(1:Nseg_stable))*gsurf(index)%xs(1:Nseg_stable)&
             *sin(gsurf(index)%xtheta(1:Nseg_stable))*dxthetadxthetans(1:Nseg_stable,i)
     dxydxthetanc(1:Nseg_stable,i) = -1.0*sin(zeta(1:Nseg_stable))*gsurf(index)%xs(1:Nseg_stable)&
             *sin(gsurf(index)%xtheta(1:Nseg_stable))*dxthetadxthetanc(1:Nseg_stable,i)
     dxydxthetans(1:Nseg_stable,i) = -1.0*sin(zeta(1:Nseg_stable))*gsurf(index)%xs(1:Nseg_stable)&
             *sin(gsurf(index)%xtheta(1:Nseg_stable))*dxthetadxthetans(1:Nseg_stable,i)
     dxzdxthetanc(1:Nseg_stable,i) =                               gsurf(index)%xs(1:Nseg_stable)&
             *cos(gsurf(index)%xtheta(1:Nseg_stable))*dxthetadxthetanc(1:Nseg_stable,i)
     dxzdxthetans(1:Nseg_stable,i) =                               gsurf(index)%xs(1:Nseg_stable)&
             *cos(gsurf(index)%xtheta(1:Nseg_stable))*dxthetadxthetans(1:Nseg_stable,i)
  enddo

  SALLOCATE( dobsupsdosnc, (1:Nseg_stable,1:NF_stable), 0.0 )
  SALLOCATE( dobsupsdosns, (1:Nseg_stable,1:NF_stable), 0.0 )
  SALLOCATE( dobsupzetadosnc, (1:Nseg_stable,1:NF_stable), 0.0 )
  SALLOCATE( dobsupzetadosns, (1:Nseg_stable,1:NF_stable), 0.0 )
  SALLOCATE( dobsupthetadosnc, (1:Nseg_stable,1:NF_stable), 0.0 )
  SALLOCATE( dobsupthetadosns, (1:Nseg_stable,1:NF_stable), 0.0 )
  
  SALLOCATE( dxbsupsdxsnc, (1:Nseg_stable,1:NF_stable), 0.0 )
  SALLOCATE( dxbsupsdxsns, (1:Nseg_stable,1:NF_stable), 0.0 )
  SALLOCATE( dxbsupzetadxsnc, (1:Nseg_stable,1:NF_stable), 0.0 )
  SALLOCATE( dxbsupzetadxsns, (1:Nseg_stable,1:NF_stable), 0.0 )
  SALLOCATE( dxbsupthetadxsnc, (1:Nseg_stable,1:NF_stable), 0.0 )
  SALLOCATE( dxbsupthetadxsns, (1:Nseg_stable,1:NF_stable), 0.0 )

  SALLOCATE( dobsupsdothetanc, (1:Nseg_stable,1:NF_stable), 0.0 )
  SALLOCATE( dobsupsdothetans, (1:Nseg_stable,1:NF_stable), 0.0 )
  SALLOCATE( dobsupzetadothetanc, (1:Nseg_stable,1:NF_stable), 0.0 )
  SALLOCATE( dobsupzetadothetans, (1:Nseg_stable,1:NF_stable), 0.0 )
  SALLOCATE( dobsupthetadothetanc, (1:Nseg_stable,1:NF_stable), 0.0 )
  SALLOCATE( dobsupthetadothetans, (1:Nseg_stable,1:NF_stable), 0.0 )

  SALLOCATE( dxbsupsdxthetanc, (1:Nseg_stable,1:NF_stable), 0.0 )
  SALLOCATE( dxbsupsdxthetans, (1:Nseg_stable,1:NF_stable), 0.0 )
  SALLOCATE( dxbsupzetadxthetanc, (1:Nseg_stable,1:NF_stable), 0.0 )
  SALLOCATE( dxbsupzetadxthetans, (1:Nseg_stable,1:NF_stable), 0.0 )
  SALLOCATE( dxbsupthetadxthetanc, (1:Nseg_stable,1:NF_stable), 0.0 )
  SALLOCATE( dxbsupthetadxthetans, (1:Nseg_stable,1:NF_stable), 0.0 )

  do i = 1, NF_stable
     dobsupsdosnc(1:Nseg_stable,i) = doxdosnc(1:Nseg_stable,i)*gradobsupsx(1:Nseg_stable) + & 
             doydosnc(1:Nseg_stable,i)*gradobsupsy(1:Nseg_stable) + dozdosnc(1:Nseg_stable,i)*gradobsupsz(1:Nseg_stable)
     dobsupsdosns(1:Nseg_stable,i) = doxdosns(1:Nseg_stable,i)*gradobsupsx(1:Nseg_stable) + &
             doydosns(1:Nseg_stable,i)*gradobsupsy(1:Nseg_stable) + dozdosns(1:Nseg_stable,i)*gradobsupsz(1:Nseg_stable)
     dobsupzetadosnc(1:Nseg_stable,i) = doxdosnc(1:Nseg_stable,i)*gradobsupzetax(1:Nseg_stable) + &
             doydosnc(1:Nseg_stable,i)*gradobsupzetay(1:Nseg_stable) + dozdosnc(1:Nseg_stable,i)*gradobsupzetaz(1:Nseg_stable)
     dobsupzetadosns(1:Nseg_stable,i) = doxdosns(1:Nseg_stable,i)*gradobsupzetax(1:Nseg_stable) + &
             doydosns(1:Nseg_stable,i)*gradobsupzetay(1:Nseg_stable) + dozdosns(1:Nseg_stable,i)*gradobsupzetaz(1:Nseg_stable)
     dobsupthetadosnc(1:Nseg_stable,i) = doxdosnc(1:Nseg_stable,i)*gradobsupthetax(1:Nseg_stable) + &
             doydosnc(1:Nseg_stable,i)*gradobsupthetay(1:Nseg_stable) + dozdosnc(1:Nseg_stable,i)*gradobsupthetaz(1:Nseg_stable)
     dobsupthetadosns(1:Nseg_stable,i) = doxdosns(1:Nseg_stable,i)*gradobsupthetax(1:Nseg_stable) + &
             doydosns(1:Nseg_stable,i)*gradobsupthetay(1:Nseg_stable) + dozdosns(1:Nseg_stable,i)*gradobsupthetaz(1:Nseg_stable)

     dxbsupsdxsnc(1:Nseg_stable,i) = dxxdxsnc(1:Nseg_stable,i)*gradxbsupsx(1:Nseg_stable) + &
             dxydxsnc(1:Nseg_stable,i)*gradxbsupsy(1:Nseg_stable) + dxzdxsnc(1:Nseg_stable,i)*gradxbsupsz(1:Nseg_stable)
     dxbsupsdxsns(1:Nseg_stable,i) = dxxdxsns(1:Nseg_stable,i)*gradxbsupsx(1:Nseg_stable) + &
             dxydxsns(1:Nseg_stable,i)*gradxbsupsy(1:Nseg_stable) + dxzdxsns(1:Nseg_stable,i)*gradxbsupsz(1:Nseg_stable)
     dxbsupzetadxsnc(1:Nseg_stable,i) = dxxdxsnc(1:Nseg_stable,i)*gradxbsupzetax(1:Nseg_stable) + &
             dxydxsnc(1:Nseg_stable,i)*gradxbsupzetay(1:Nseg_stable) + dxzdxsnc(1:Nseg_stable,i)*gradxbsupzetaz(1:Nseg_stable)
     dxbsupzetadxsns(1:Nseg_stable,i) = dxxdxsns(1:Nseg_stable,i)*gradxbsupzetax(1:Nseg_stable) + &
             dxydxsns(1:Nseg_stable,i)*gradxbsupzetay(1:Nseg_stable) + dxzdxsns(1:Nseg_stable,i)*gradxbsupzetaz(1:Nseg_stable)
     dxbsupthetadxsnc(1:Nseg_stable,i) = dxxdxsnc(1:Nseg_stable,i)*gradxbsupthetax(1:Nseg_stable) + &
             dxydxsnc(1:Nseg_stable,i)*gradxbsupthetay(1:Nseg_stable) + dxzdxsnc(1:Nseg_stable,i)*gradxbsupthetaz(1:Nseg_stable)
     dxbsupthetadxsns(1:Nseg_stable,i) = dxxdxsns(1:Nseg_stable,i)*gradxbsupthetax(1:Nseg_stable) + &
             dxydxsns(1:Nseg_stable,i)*gradxbsupthetay(1:Nseg_stable) + dxzdxsns(1:Nseg_stable,i)*gradxbsupthetaz(1:Nseg_stable)

     dobsupsdothetanc(1:Nseg_stable,i) = doxdothetanc(1:Nseg_stable,i)*gradobsupsx(1:Nseg_stable) + &
             doydothetanc(1:Nseg_stable,i)*gradobsupsy(1:Nseg_stable) + dozdothetanc(1:Nseg_stable,i)*gradobsupsz(1:Nseg_stable)
     dobsupsdothetans(1:Nseg_stable,i) = doxdothetans(1:Nseg_stable,i)*gradobsupsx(1:Nseg_stable) + &
             doydothetans(1:Nseg_stable,i)*gradobsupsy(1:Nseg_stable) + dozdothetans(1:Nseg_stable,i)*gradobsupsz(1:Nseg_stable)
     dobsupzetadothetanc(1:Nseg_stable,i) = doxdothetanc(1:Nseg_stable,i)*gradobsupzetax(1:Nseg_stable) + &
             doydothetanc(1:Nseg_stable,i)*gradobsupzetay(1:Nseg_stable) + dozdothetanc(1:Nseg_stable,i)*gradobsupzetaz(1:Nseg_stable)
     dobsupzetadothetans(1:Nseg_stable,i) = doxdothetans(1:Nseg_stable,i)*gradobsupzetax(1:Nseg_stable) + &
             doydothetans(1:Nseg_stable,i)*gradobsupzetay(1:Nseg_stable) + dozdothetans(1:Nseg_stable,i)*gradobsupzetaz(1:Nseg_stable)
     dobsupthetadothetanc(1:Nseg_stable,i) = doxdothetanc(1:Nseg_stable,i)*gradobsupthetax(1:Nseg_stable) + &
             doydothetanc(1:Nseg_stable,i)*gradobsupthetay(1:Nseg_stable) + dozdothetanc(1:Nseg_stable,i)*gradobsupthetaz(1:Nseg_stable)
     dobsupthetadothetans(1:Nseg_stable,i) = doxdothetans(1:Nseg_stable,i)*gradobsupthetax(1:Nseg_stable) + &
             doydothetans(1:Nseg_stable,i)*gradobsupthetay(1:Nseg_stable) + dozdothetans(1:Nseg_stable,i)*gradobsupthetaz(1:Nseg_stable) 
     
     dxbsupsdxthetanc(1:Nseg_stable,i) = dxxdxthetanc(1:Nseg_stable,i)*gradxbsupsx(1:Nseg_stable) + &
             dxydxthetanc(1:Nseg_stable,i)*gradxbsupsy(1:Nseg_stable) + dxzdxthetanc(1:Nseg_stable,i)*gradxbsupsz(1:Nseg_stable)
     dxbsupsdxthetans(1:Nseg_stable,i) = dxxdxthetans(1:Nseg_stable,i)*gradxbsupsx(1:Nseg_stable) + &
             dxydxthetans(1:Nseg_stable,i)*gradxbsupsy(1:Nseg_stable) + dxzdxthetans(1:Nseg_stable,i)*gradxbsupsz(1:Nseg_stable)
     dxbsupzetadxthetanc(1:Nseg_stable,i) = dxxdxthetanc(1:Nseg_stable,i)*gradxbsupzetax(1:Nseg_stable) + &
             dxydxthetanc(1:Nseg_stable,i)*gradxbsupzetay(1:Nseg_stable) + dxzdxthetanc(1:Nseg_stable,i)*gradxbsupzetaz(1:Nseg_stable)
     dxbsupzetadxthetans(1:Nseg_stable,i) = dxxdxthetans(1:Nseg_stable,i)*gradxbsupzetax(1:Nseg_stable) + &
             dxydxthetans(1:Nseg_stable,i)*gradxbsupzetay(1:Nseg_stable) + dxzdxthetans(1:Nseg_stable,i)*gradxbsupzetaz(1:Nseg_stable)
     dxbsupthetadxthetanc(1:Nseg_stable,i) = dxxdxthetanc(1:Nseg_stable,i)*gradxbsupthetax(1:Nseg_stable) + &
             dxydxthetanc(1:Nseg_stable,i)*gradxbsupthetay(1:Nseg_stable) + dxzdxthetanc(1:Nseg_stable,i)*gradxbsupthetaz(1:Nseg_stable)
     dxbsupthetadxthetans(1:Nseg_stable,i) = dxxdxthetans(1:Nseg_stable,i)*gradxbsupthetax(1:Nseg_stable) + &
             dxydxthetans(1:Nseg_stable,i)*gradxbsupthetay(1:Nseg_stable) + dxzdxthetans(1:Nseg_stable,i)*gradxbsupthetaz(1:Nseg_stable)
  enddo

  SALLOCATE( dofdosnc, (1:Nseg_stable,1:NF_stable), 0.0 ) 
  SALLOCATE( dofdosns, (1:Nseg_stable,1:NF_stable), 0.0 )
  SALLOCATE( dofdothetanc, (1:Nseg_stable,1:NF_stable), 0.0 )
  SALLOCATE( dofdothetans, (1:Nseg_stable,1:NF_stable), 0.0 )

  SALLOCATE( dogdosnc, (1:Nseg_stable,1:NF_stable), 0.0 )
  SALLOCATE( dogdosns, (1:Nseg_stable,1:NF_stable), 0.0 )
  SALLOCATE( dogdothetanc, (1:Nseg_stable,1:NF_stable), 0.0 )
  SALLOCATE( dogdothetans, (1:Nseg_stable,1:NF_stable), 0.0 )

  SALLOCATE( dxfdxsnc, (1:Nseg_stable,1:NF_stable), 0.0 )
  SALLOCATE( dxfdxsns, (1:Nseg_stable,1:NF_stable), 0.0 )
  SALLOCATE( dxfdxthetanc, (1:Nseg_stable,1:NF_stable), 0.0 )
  SALLOCATE( dxfdxthetans, (1:Nseg_stable,1:NF_stable), 0.0 )

  SALLOCATE( dxgdxsnc, (1:Nseg_stable,1:NF_stable), 0.0 )
  SALLOCATE( dxgdxsns, (1:Nseg_stable,1:NF_stable), 0.0 )
  SALLOCATE( dxgdxthetanc, (1:Nseg_stable,1:NF_stable), 0.0 )
  SALLOCATE( dxgdxthetans, (1:Nseg_stable,1:NF_stable), 0.0 )
  
  do i = 1, NF_stable
     dofdosnc(1:Nseg_stable,i) = dobsupsdosnc(1:Nseg_stable,i) / gsurf(index)%obsupzeta(1:Nseg_stable) &
             - gsurf(index)%obsups(1:Nseg_stable)*dobsupzetadosnc(1:Nseg_stable,i) & 
             / gsurf(index)%obsupzeta(1:Nseg_stable)**2 - dosdotdosnc(1:Nseg_stable,i)
     dofdosns(1:Nseg_stable,i) = dobsupsdosns(1:Nseg_stable,i) / gsurf(index)%obsupzeta(1:Nseg_stable) &
             - gsurf(index)%obsups(1:Nseg_stable)*dobsupzetadosns(1:Nseg_stable,i) & 
             / gsurf(index)%obsupzeta(1:Nseg_stable)**2 - dosdotdosns(1:Nseg_stable,i)
     dofdothetanc(1:Nseg_stable,i) = dobsupsdothetanc(1:Nseg_stable,i) / gsurf(index)%obsupzeta(1:Nseg_stable) &
             - gsurf(index)%obsups(1:Nseg_stable)*dobsupzetadothetanc(1:Nseg_stable,i) &
             / gsurf(index)%obsupzeta(1:Nseg_stable)**2
     dofdothetans(1:Nseg_stable,i) = dobsupsdothetans(1:Nseg_stable,i) / gsurf(index)%obsupzeta(1:Nseg_stable) &
             - gsurf(index)%obsups(1:Nseg_stable)*dobsupzetadothetans(1:Nseg_stable,i) &
             / gsurf(index)%obsupzeta(1:Nseg_stable)**2

     dogdosnc(1:Nseg_stable,i) = dobsupthetadosnc(1:Nseg_stable,i) / gsurf(index)%obsupzeta(1:Nseg_stable) &
             - gsurf(index)%obsuptheta(1:Nseg_stable)*dobsupzetadosnc(1:Nseg_stable,i) &
             / gsurf(index)%obsupzeta(1:Nseg_stable)**2
     dogdosns(1:Nseg_stable,i) = dobsupthetadosns(1:Nseg_stable,i) / gsurf(index)%obsupzeta(1:Nseg_stable) &
             - gsurf(index)%obsuptheta(1:Nseg_stable)*dobsupzetadosns(1:Nseg_stable,i) &
             / gsurf(index)%obsupzeta(1:Nseg_stable)**2
     dogdothetanc(1:Nseg_stable,i) = dobsupthetadothetanc(1:Nseg_stable,i) / gsurf(index)%obsupzeta(1:Nseg_stable) &
             - gsurf(index)%obsuptheta(1:Nseg_stable)*dobsupzetadothetanc(1:Nseg_stable,i) &
             / gsurf(index)%obsupzeta(1:Nseg_stable)**2 - dothetadotdothetanc(1:Nseg_stable,i)
     dogdothetans(1:Nseg_stable,i) = dobsupthetadothetans(1:Nseg_stable,i) / gsurf(index)%obsupzeta(1:Nseg_stable) &
             - gsurf(index)%obsuptheta(1:Nseg_stable)*dobsupzetadothetans(1:Nseg_stable,i) &
             / gsurf(index)%obsupzeta(1:Nseg_stable)**2 - dothetadotdothetans(1:Nseg_stable,i)
     
     dxfdxsnc(1:Nseg_stable,i) = dxbsupsdxsnc(1:Nseg_stable,i) / gsurf(index)%xbsupzeta(1:Nseg_stable) &
             - gsurf(index)%xbsups(1:Nseg_stable)*dxbsupzetadxsnc(1:Nseg_stable,i) &
             / gsurf(index)%xbsupzeta(1:Nseg_stable)**2 - dxsdotdxsnc(1:Nseg_stable,i)
     dxfdxsns(1:Nseg_stable,i) = dxbsupsdxsns(1:Nseg_stable,i) / gsurf(index)%xbsupzeta(1:Nseg_stable) &
             - gsurf(index)%xbsups(1:Nseg_stable)*dxbsupzetadxsns(1:Nseg_stable,i) &
             / gsurf(index)%xbsupzeta(1:Nseg_stable)**2 - dxsdotdxsns(1:Nseg_stable,i)
     dxfdxthetanc(1:Nseg_stable,i) = dxbsupsdxthetanc(1:Nseg_stable,i) / gsurf(index)%xbsupzeta(1:Nseg_stable) &
             - gsurf(index)%xbsups(1:Nseg_stable)*dxbsupzetadxthetanc(1:Nseg_stable,i) &
             / gsurf(index)%xbsupzeta(1:Nseg_stable)**2 
     dxfdxthetans(1:Nseg_stable,i) = dxbsupsdxthetans(1:Nseg_stable,i) / gsurf(index)%xbsupzeta(1:Nseg_stable) &
             - gsurf(index)%xbsups(1:Nseg_stable)*dxbsupzetadxthetans(1:Nseg_stable,i) &
             / gsurf(index)%xbsupzeta(1:Nseg_stable)**2
     
     dxgdxsnc(1:Nseg_stable,i) = dxbsupthetadxsnc(1:Nseg_stable,i) / gsurf(index)%xbsupzeta(1:Nseg_stable) &
             - gsurf(index)%xbsuptheta(1:Nseg_stable)*dxbsupzetadxsnc(1:Nseg_stable,i) &
             / gsurf(index)%xbsupzeta(1:Nseg_stable)**2
     dxgdxsns(1:Nseg_stable,i) = dxbsupthetadxsns(1:Nseg_stable,i) / gsurf(index)%xbsupzeta(1:Nseg_stable) &
             - gsurf(index)%xbsuptheta(1:Nseg_stable)*dxbsupzetadxsns(1:Nseg_stable,i) &
             / gsurf(index)%xbsupzeta(1:Nseg_stable)**2
     dxgdxthetanc(1:Nseg_stable,i) = dxbsupthetadxthetanc(1:Nseg_stable,i) / gsurf(index)%xbsupzeta(1:Nseg_stable) &
             - gsurf(index)%xbsuptheta(1:Nseg_stable)*dxbsupzetadxthetanc(1:Nseg_stable,i) &
             / gsurf(index)%xbsupzeta(1:Nseg_stable)**2 - dxthetadotdxthetanc(1:Nseg_stable,i)
     dxgdxthetans(1:Nseg_stable,i) = dxbsupthetadxthetans(1:Nseg_stable,i) / gsurf(index)%xbsupzeta(1:Nseg_stable) &
             - gsurf(index)%xbsuptheta(1:Nseg_stable)*dxbsupzetadxthetans(1:Nseg_stable,i) &
             / gsurf(index)%xbsupzeta(1:Nseg_stable)**2 - dxthetadotdxthetans(1:Nseg_stable,i)
  enddo

  SALLOCATE( dFdosnc, (1:NF_stable), 0.0 )
  SALLOCATE( dFdosns, (1:NF_stable), 0.0 )
  SALLOCATE( dFdothetanc, (1:NF_stable), 0.0 )
  SALLOCATE( dFdothetans, (1:NF_stable), 0.0 )
  
  SALLOCATE( dFdxsnc, (1:NF_stable), 0.0 )
  SALLOCATE( dFdxsns, (1:NF_stable), 0.0 )
  SALLOCATE( dFdxthetanc, (1:NF_stable), 0.0 )
  SALLOCATE( dFdxthetans, (1:NF_stable), 0.0 )

  do i = 1, NF_stable
     dFdosnc(i) = sum( gsurf(index)%of(1:Nseg_stable-1)*dofdosnc(1:Nseg_stable-1,i) + &
             gsurf(index)%og(1:Nseg_stable-1)*dogdosnc(1:Nseg_stable-1,i) )
     dFdosns(i) = sum( gsurf(index)%of(1:Nseg_stable-1)*dofdosns(1:Nseg_stable-1,i) + &
             gsurf(index)%og(1:Nseg_stable-1)*dogdosns(1:Nseg_stable-1,i) )
     dFdothetanc(i) = sum( gsurf(index)%of(1:Nseg_stable-1)*dofdothetanc(1:Nseg_stable-1,i) + &
             gsurf(index)%og(1:Nseg_stable-1)*dogdothetanc(1:Nseg_stable-1,i) )
     dFdothetans(i) = sum( gsurf(index)%of(1:Nseg_stable-1)*dofdothetans(1:Nseg_stable-1,i) + &
             gsurf(index)%og(1:Nseg_stable-1)*dogdothetans(1:Nseg_stable-1,i) )
     dFdxsnc(i) = sum( gsurf(index)%xf(1:Nseg_stable-1)*dxfdxsnc(1:Nseg_stable-1,i) + &
             gsurf(index)%xg(1:Nseg_stable-1)*dxgdxsnc(1:Nseg_stable-1,i) )
     dFdxsns(i) = sum( gsurf(index)%xf(1:Nseg_stable-1)*dxfdxsns(1:Nseg_stable-1,i) + &
             gsurf(index)%xg(1:Nseg_stable-1)*dxgdxsns(1:Nseg_stable-1,i) )
     dFdxthetanc(i) = sum( gsurf(index)%xf(1:Nseg_stable-1)*dxfdxthetanc(1:Nseg_stable-1,i) + &
             gsurf(index)%xg(1:Nseg_stable-1)*dxgdxthetanc(1:Nseg_stable-1,i) )
     dFdxthetans(i) = sum( gsurf(index)%xf(1:Nseg_stable-1)*dxfdxthetans(1:Nseg_stable-1,i) + &
             gsurf(index)%xg(1:Nseg_stable-1)*dxgdxthetans(1:Nseg_stable-1,i) )
  enddo
  dFdosnc(1:NF_stable) = 2.0*pi2*q*dFdosnc(1:NF_stable)/(Nseg_stable-1)
  dFdosns(1:NF_stable) = 2.0*pi2*q*dFdosns(1:NF_stable)/(Nseg_stable-1)
  dFdothetanc(1:NF_stable) = 2.0*pi2*q*dFdothetanc(1:NF_stable)/(Nseg_stable-1)
  dFdothetans(1:NF_stable) = 2.0*pi2*q*dFdothetans(1:NF_stable)/(Nseg_stable-1)
  dFdxsnc(1:NF_stable) = 2.0*pi2*q*dFdxsnc(1:NF_stable)/(Nseg_stable-1)
  dFdxsns(1:NF_stable) = 2.0*pi2*q*dFdxsns(1:NF_stable)/(Nseg_stable-1)
  dFdxthetanc(1:NF_stable) = 2.0*pi2*q*dFdxthetanc(1:NF_stable)/(Nseg_stable-1)
  dFdxthetans(1:NF_stable) = 2.0*pi2*q*dFdxthetans(1:NF_stable)/(Nseg_stable-1)

  gsurf(index)%dFdxdof_stable(            1:  NF_stable) = dFdosnc(1:NF_stable)
  gsurf(index)%dFdxdof_stable(  NF_stable+1:2*NF_stable) = dFdosns(1:NF_stable)
  gsurf(index)%dFdxdof_stable(2*NF_stable+1:3*NF_stable) = dFdothetanc(1:NF_stable)
  gsurf(index)%dFdxdof_stable(3*NF_stable+1:4*NF_stable) = dFdothetans(1:NF_stable)
  gsurf(index)%dFdxdof_stable(4*NF_stable+1:5*NF_stable) = dFdxsnc(1:NF_stable)
  gsurf(index)%dFdxdof_stable(5*NF_stable+1:6*NF_stable) = dFdxsns(1:NF_stable)
  gsurf(index)%dFdxdof_stable(6*NF_stable+1:7*NF_stable) = dFdxthetanc(1:NF_stable)
  gsurf(index)%dFdxdof_stable(7*NF_stable+1:8*NF_stable) = dFdxthetans(1:NF_stable)

  ! Deallocation statements, consider moving
  DALLOCATE( gradobsupsx )
  DALLOCATE( gradobsupsy )
  DALLOCATE( gradobsupsz )
  DALLOCATE( gradobsupzetax )
  DALLOCATE( gradobsupzetay )
  DALLOCATE( gradobsupzetaz )
  DALLOCATE( gradobsupthetax )
  DALLOCATE( gradobsupthetay )
  DALLOCATE( gradobsupthetaz )
  DALLOCATE( gradxbsupsx )
  DALLOCATE( gradxbsupsy )
  DALLOCATE( gradxbsupsz )
  DALLOCATE( gradxbsupzetax )
  DALLOCATE( gradxbsupzetay )
  DALLOCATE( gradxbsupzetaz )
  DALLOCATE( gradxbsupthetax )
  DALLOCATE( gradxbsupthetay )
  DALLOCATE( gradxbsupthetaz )
  DALLOCATE( oxneg )
  DALLOCATE( oyneg )
  DALLOCATE( ozneg )
  DALLOCATE( oxpos )
  DALLOCATE( oypos )
  DALLOCATE( ozpos )
  DALLOCATE( xxneg )
  DALLOCATE( xyneg )
  DALLOCATE( xzneg )
  DALLOCATE( xxpos )
  DALLOCATE( xypos )
  DALLOCATE( xzpos )
  DALLOCATE( negbsups )
  DALLOCATE( negbsupzeta )
  DALLOCATE( negbsuptheta )
  DALLOCATE( posbsups )
  DALLOCATE( posbsupzeta )
  DALLOCATE( posbsuptheta )
  DALLOCATE( dosdosnc )
  DALLOCATE( dosdosns )
  DALLOCATE( dothetadothetanc )
  DALLOCATE( dothetadothetans )
  DALLOCATE( dxsdxsnc )
  DALLOCATE( dxsdxsns )
  DALLOCATE( dxthetadxthetanc )
  DALLOCATE( dxthetadxthetans )
  DALLOCATE( dosdotdosnc )
  DALLOCATE( dosdotdosns )
  DALLOCATE( dothetadotdothetanc )
  DALLOCATE( dothetadotdothetans )
  DALLOCATE( dxsdotdxsnc )
  DALLOCATE( dxsdotdxsns )
  DALLOCATE( dxthetadotdxthetanc )
  DALLOCATE( dxthetadotdxthetans )
  DALLOCATE( doxdosnc )
  DALLOCATE( doxdosns )
  DALLOCATE( doydosnc )
  DALLOCATE( doydosns )
  DALLOCATE( dozdosnc )
  DALLOCATE( dozdosns )
  DALLOCATE( dxxdxsnc )
  DALLOCATE( dxxdxsns )
  DALLOCATE( dxydxsnc )
  DALLOCATE( dxydxsns )
  DALLOCATE( dxzdxsnc )
  DALLOCATE( dxzdxsns )
  DALLOCATE( doxdothetanc )
  DALLOCATE( doxdothetans )
  DALLOCATE( doydothetanc )
  DALLOCATE( doydothetans )
  DALLOCATE( dozdothetanc )
  DALLOCATE( dozdothetans )
  DALLOCATE( dxxdxthetanc )
  DALLOCATE( dxxdxthetans )
  DALLOCATE( dxydxthetanc )
  DALLOCATE( dxydxthetans )
  DALLOCATE( dxzdxthetanc )
  DALLOCATE( dxzdxthetans )
  DALLOCATE( dobsupsdosnc )
  DALLOCATE( dobsupsdosns )
  DALLOCATE( dobsupzetadosnc )
  DALLOCATE( dobsupzetadosns )
  DALLOCATE( dobsupthetadosnc )
  DALLOCATE( dobsupthetadosns )
  DALLOCATE( dxbsupsdxsnc )
  DALLOCATE( dxbsupsdxsns )
  DALLOCATE( dxbsupzetadxsnc )
  DALLOCATE( dxbsupzetadxsns )
  DALLOCATE( dxbsupthetadxsnc )
  DALLOCATE( dxbsupthetadxsns )
  DALLOCATE( dobsupsdothetanc )
  DALLOCATE( dobsupsdothetans )
  DALLOCATE( dobsupzetadothetanc )
  DALLOCATE( dobsupzetadothetans )
  DALLOCATE( dobsupthetadothetanc )
  DALLOCATE( dobsupthetadothetans )
  DALLOCATE( dxbsupsdxthetanc )
  DALLOCATE( dxbsupsdxthetans )
  DALLOCATE( dxbsupzetadxthetanc )
  DALLOCATE( dxbsupzetadxthetans )
  DALLOCATE( dxbsupthetadxthetanc )
  DALLOCATE( dxbsupthetadxthetans )  
  DALLOCATE( dofdosnc )
  DALLOCATE( dofdosns )
  DALLOCATE( dofdothetanc )
  DALLOCATE( dofdothetans )
  DALLOCATE( dogdosnc )
  DALLOCATE( dogdosns )
  DALLOCATE( dogdothetanc )
  DALLOCATE( dogdothetans )
  DALLOCATE( dxfdxsnc )
  DALLOCATE( dxfdxsns )
  DALLOCATE( dxfdxthetanc )
  DALLOCATE( dxfdxthetans )
  DALLOCATE( dxgdxsnc )
  DALLOCATE( dxgdxsns )
  DALLOCATE( dxgdxthetanc )
  DALLOCATE( dxgdxthetans ) 
  DALLOCATE( dFdosnc )
  DALLOCATE( dFdosns )
  DALLOCATE( dFdothetanc )
  DALLOCATE( dFdothetans )
  DALLOCATE( dFdxsnc )
  DALLOCATE( dFdxsns )
  DALLOCATE( dFdxthetanc )
  DALLOCATE( dFdxthetans )
  DALLOCATE( oxsmall )
  DALLOCATE( oysmall )
  DALLOCATE( ozsmall )
  DALLOCATE( xxsmall )
  DALLOCATE( xysmall )
  DALLOCATE( xzsmall )

  return

end subroutine calcfg_deriv

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine congrad_stable(index)
  use globals, only: dp, sqrtmachprec, myid, ounit, CG_maxiter, CG_xtol, &
       exit_signal, tstart, tfinish, gsurf, output_use, CG_maxiter_s, &
       CG_xtol_s, MPI_COMM_FOCUS
  
  use mpi
  implicit none

  INTEGER, INTENT(in)     :: index
  INTEGER                 :: ierr, astat, iter_stable, n_stable, nfunc_stable, ngrad_stable, status, &
                             CG_maxiter_hold, NF_stable, j, cont
  REAL                    :: f_stable, f_hold, gnorm_stable, CG_xtol_hold, step, initial_step
  REAL, dimension(1:gsurf(index)%Ndof_stable) :: x_stable, g_stable, d_stable, xtemp_stable, gtemp_stable, &
                             x_hold
  EXTERNAL                :: myvalue_stable, mygrad_stable

  iter_stable = 0
  n_stable = gsurf(index)%Ndof_stable
  x_stable(1:n_stable) = gsurf(index)%xdof_stable(1:n_stable)

  CG_maxiter_hold = CG_maxiter
  CG_maxiter = CG_maxiter_s
  CG_xtol_hold = CG_xtol
  CG_xtol = CG_xtol_s

  output_use = 0

  !call cg_descent_stable (CG_xtol, x_stable, n_stable, myvalue_stable, mygrad_stable, status, &
  !       gnorm_stable, f_stable, iter_stable, nfunc_stable, ngrad_stable, d_stable, g_stable, xtemp_stable, gtemp_stable)

  !call cg_descent (CG_xtol, x_stable, n_stable, myvalue_stable, mygrad_stable, status, &
  !       gnorm_stable, f_stable, iter_stable, nfunc_stable, ngrad_stable, d_stable, g_stable, xtemp_stable, gtemp_stable)

  !if (myid == 0) then
  !   select case (status)
  !   case (0)
  !      write(ounit, '("congrad : status="I1": convergence tolerance satisfied.")')  status
  !   case (1)
  !      write(ounit, '("congrad : status="I1": change in func <= feps*|f|.")')  status
  !   case (2)
  !      write(ounit, '("congrad : status="I1": total iterations exceeded maxit.")')  status
  !   case (3)
  !      write(ounit, '("congrad : status="I1": slope always negative in line search.")')  status
  !   case (4)
  !      write(ounit, '("congrad : status="I1": number secant iterations exceed nsecant.")')  status
  !   case (5)
  !      write(ounit, '("congrad : status="I1": search direction not a descent direction.")')  status
  !   case (6)
  !      write(ounit, '("congrad : status="I1": line search fails in initial interval.")')  status
  !   case (7)
  !      write(ounit, '("congrad : status="I1": line search fails during bisection.")')  status
  !   case (8)
  !      write(ounit, '("congrad : status="I1": line search fails during interval update.")')  status
  !   case default
  !      write(ounit, '("congrad : status="I1": unknow options!")')  status
  !   end select
  !end if

  !if(myid .eq. 0) write(ounit, '("congrad : Stable conjugate gradient finished.")')

  initial_step = 1.0E-14
  do j = 1, CG_maxiter
     call myvalue_stable(f_stable, x_stable, n_stable)
     call mygrad_stable(g_stable, x_stable, n_stable)
     step = initial_step
     x_hold(1:n_stable) = x_stable(1:n_stable)
     f_hold = f_stable
     cont = 1
     if (myid == 0 ) write(ounit, '("Value of F:   "ES12.5)') f_stable
     do while (cont .eq. 1)
        x_stable(1:n_stable) = x_hold(1:n_stable) - step*g_stable(1:n_stable)
        call myvalue_stable(f_stable, x_stable, n_stable)
        !if (myid == 0 ) write(ounit, '("Value of F:   "ES12.5)') f_stable
        if ( f_stable > f_hold ) then
           cont = 0
           x_stable(1:n_stable) = x_hold(1:n_stable) - step*g_stable(1:n_stable)/2.0
           call myvalue_stable(f_stable, x_stable, n_stable)
        endif
        f_hold = f_stable
        step = step*2.0
     enddo
  enddo
  
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

  output_use = 1

  return

end subroutine congrad_stable

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

subroutine myvalue_stable(f_stable, x_stable, n_stable)
  use globals, only: dp, myid, ounit, ierr, gsurf, MPI_COMM_FOCUS
  use mpi
  implicit none

  INTEGER, INTENT(in) :: n_stable
  REAL, INTENT(in)    :: x_stable(n_stable)
  REAL, INTENT(out)   :: f_stable

  INTEGER             :: index

  call MPI_BARRIER( MPI_COMM_FOCUS, ierr ) ! wait all cpus;

  index = 1
  
  call unpacking_stable(x_stable,index)
  call calcfg(index)
  f_stable = gsurf(index)%F

  return

end subroutine myvalue_stable

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

subroutine mygrad_stable(g_stable, x_stable, n_stable)
  use globals, only: dp, myid, ounit, ierr, gsurf, MPI_COMM_FOCUS
  use mpi
  implicit none

  INTEGER, INTENT(in) :: n_stable
  REAL, INTENT(in)    :: x_stable(n_stable)
  REAL, INTENT(out)   :: g_stable(n_stable)

  INTEGER             :: index

  call MPI_BARRIER( MPI_COMM_FOCUS, ierr ) ! wait all cpus;
  
  index = 1

  call unpacking_stable(x_stable,index)
  call calcfg_deriv(index)
  g_stable = gsurf(index)%dFdxdof_stable

  return

end subroutine mygrad_stable

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

subroutine unpacking_stable(x_stable, index)
  use globals, only: dp, myid, ounit, ierr, gsurf, MPI_COMM_FOCUS
  use mpi
  implicit none

  INTEGER, INTENT(in) :: index
  REAL, INTENT(in)    :: x_stable(gsurf(index)%Ndof_stable)

  INTEGER             :: NF_stable

  call MPI_BARRIER( MPI_COMM_FOCUS, ierr ) ! wait all cpus;

  NF_stable = gsurf(index)%NF_stable
  gsurf(index)%osnc(1:NF_stable)     = x_stable(            1:  NF_stable)
  gsurf(index)%osns(1:NF_stable)     = x_stable(  NF_stable+1:2*NF_stable)
  gsurf(index)%othetanc(1:NF_stable) = x_stable(2*NF_stable+1:3*NF_stable)
  gsurf(index)%othetans(1:NF_stable) = x_stable(3*NF_stable+1:4*NF_stable)
  gsurf(index)%xsnc(1:NF_stable)     = x_stable(4*NF_stable+1:5*NF_stable)
  gsurf(index)%xsns(1:NF_stable)     = x_stable(5*NF_stable+1:6*NF_stable)
  gsurf(index)%xthetanc(1:NF_stable) = x_stable(6*NF_stable+1:7*NF_stable)
  gsurf(index)%xthetans(1:NF_stable) = x_stable(7*NF_stable+1:8*NF_stable)

  return

end subroutine unpacking_stable

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

subroutine calcbsup(index,ox,oy,oz,obsups,obsupzeta,obsuptheta)
  use globals, only: dp, zero, half, pi2, ncpu, myid, resbn_m, gsurf, Ncoils, ounit, &
                     MPI_COMM_FOCUS
  
  use mpi
  implicit none
  INTEGER, INTENT(IN) :: index
  REAL, INTENT(IN)    :: ox(1:gsurf(index)%Nseg_stable), oy(1:gsurf(index)%Nseg_stable), &
                         oz(1:gsurf(index)%Nseg_stable)
  REAL, INTENT(OUT)   :: obsups(1:gsurf(index)%Nseg_stable), obsupzeta(1:gsurf(index)%Nseg_stable), &
                         obsuptheta(1:gsurf(index)%Nseg_stable)

  INTEGER             :: astat, ierr, i, Nseg_stable, icoil
  REAL                :: Bxhold, Byhold, Bzhold, q
  REAL, ALLOCATABLE   :: zeta(:), gradosx(:), gradosy(:), gradosz(:), os(:), gradozetax(:), &
                         gradozetay(:), gradozetaz(:), gradothetax(:), gradothetay(:), &
                         gradothetaz(:), tanotheta(:), oBx(:), oBy(:), oBz(:), odRadx(:), odRady(:), &
                         odZadx(:), odZady(:)

  Nseg_stable = gsurf(index)%Nseg_stable
  SALLOCATE( zeta, (1:Nseg_stable), 0.0)
  zeta(1:Nseg_stable) = gsurf(index)%zeta(1:Nseg_stable)

  q = real(resbn_m)

  ! Calculate derivative of axis w.r.t x and y, variables are local
  SALLOCATE( odRadx, (1:Nseg_stable), 0.0 )
  SALLOCATE( odRady, (1:Nseg_stable), 0.0 )
  SALLOCATE( odZadx, (1:Nseg_stable), 0.0 )
  SALLOCATE( odZady, (1:Nseg_stable), 0.0 )
  do i = 1, gsurf(index)%NF_axis
     odRadx(1:Nseg_stable) = odRadx(1:Nseg_stable) + &
             gsurf(index)%axisn(i)*gsurf(index)%axisrnc(i)*sin(gsurf(index)%axisn(i)*zeta(1:Nseg_stable))*oy(1:Nseg_stable) / &
             ( ox(1:Nseg_stable)**2 * ( oy(1:Nseg_stable)**2/ox(1:Nseg_stable)**2 + 1 ) )
     odRady(1:Nseg_stable) = odRady(1:Nseg_stable) - &
             gsurf(index)%axisn(i)*gsurf(index)%axisrnc(i)*sin(gsurf(index)%axisn(i)*zeta(1:Nseg_stable)) / &
             ( ox(1:Nseg_stable) * ( oy(1:Nseg_stable)**2/ox(1:Nseg_stable)**2 + 1 ) )
     odZadx(1:Nseg_stable) = odZadx(1:Nseg_stable) - &
             gsurf(index)%axisn(i)*gsurf(index)%axiszns(i)*cos(gsurf(index)%axisn(i)*zeta(1:Nseg_stable))*oy(1:Nseg_stable) / &
             ( ox(1:Nseg_stable)**2 * ( oy(1:Nseg_stable)**2/ox(1:Nseg_stable)**2 + 1 ) )
     odZady(1:Nseg_stable) = odZady(1:Nseg_stable) + &
             gsurf(index)%axisn(i)*gsurf(index)%axiszns(i)*cos(gsurf(index)%axisn(i)*zeta(1:Nseg_stable)) / &
             ( ox(1:Nseg_stable) * ( oy(1:Nseg_stable)**2/ox(1:Nseg_stable)**2 + 1 ) )
  enddo 

  ! Calculate contravariant basis, variables are local
  SALLOCATE( gradosx, (1:Nseg_stable), 0.0 )
  SALLOCATE( gradosy, (1:Nseg_stable), 0.0 )
  SALLOCATE( gradosz, (1:Nseg_stable), 0.0 )
  SALLOCATE( os, (1:Nseg_stable), 0.0 )
  
  os(1:Nseg_stable) = sqrt( (sqrt(ox(1:Nseg_stable)**2+oy(1:Nseg_stable)**2)-gsurf(index)%Ra(1:Nseg_stable))**2 + &
          (oz(1:Nseg_stable)-gsurf(index)%Za(1:Nseg_stable))**2 )
 
  gradosx(1:Nseg_stable) = ( sqrt( ox(1:Nseg_stable)**2 + oy(1:Nseg_stable)**2 ) - &
          gsurf(index)%Ra(1:Nseg_stable) )*( (ox(1:Nseg_stable)**2 + &
          oy(1:Nseg_stable)**2)**(-.5)*ox(1:Nseg_stable) - &
          odRadx(1:Nseg_stable) )/os(1:Nseg_stable) - &
          (oz(1:Nseg_stable) - gsurf(index)%Za(1:Nseg_stable))*odZadx(1:Nseg_stable)/os(1:Nseg_stable)
  
  gradosy(1:Nseg_stable) = ( sqrt( ox(1:Nseg_stable)**2 + oy(1:Nseg_stable)**2 ) - &
          gsurf(index)%Ra(1:Nseg_stable) )*( (ox(1:Nseg_stable)**2 + &
          oy(1:Nseg_stable)**2)**(-.5)*oy(1:Nseg_stable) - &
          odRady(1:Nseg_stable) )/os(1:Nseg_stable) - &
          (oz(1:Nseg_stable) - gsurf(index)%Za(1:Nseg_stable))*odZady(1:Nseg_stable)/os(1:Nseg_stable)
 
  gradosz(1:Nseg_stable) = ( oz(1:Nseg_stable) - gsurf(index)%Za(1:Nseg_stable) )/os(1:Nseg_stable)

  SALLOCATE( gradozetax, (1:Nseg_stable), 0.0 )
  SALLOCATE( gradozetay, (1:Nseg_stable), 0.0 )
  SALLOCATE( gradozetaz, (1:Nseg_stable), 0.0 )
  
  gradozetax(1:Nseg_stable) = -1.0*(oy(1:Nseg_stable)**2/ox(1:Nseg_stable)**2+1)**-1.0 * &
          oy(1:Nseg_stable)/ox(1:Nseg_stable)**2
  
  gradozetay(1:Nseg_stable) = (oy(1:Nseg_stable)**2/ox(1:Nseg_stable)**2+1)**-1.0 * &
          1.0/ox(1:Nseg_stable)

  !gradozetaz(1:Nseg_stable) = 0.0

  SALLOCATE( gradothetax, (1:Nseg_stable), 0.0 )
  SALLOCATE( gradothetay, (1:Nseg_stable), 0.0 )
  SALLOCATE( gradothetaz, (1:Nseg_stable), 0.0 )
  SALLOCATE( tanotheta, (1:Nseg_stable), 0.0 )

  tanotheta(1:Nseg_stable) = (oz(1:Nseg_stable)-gsurf(index)%Za(1:Nseg_stable)) / &
          (sqrt(ox(1:Nseg_stable)**2+oy(1:Nseg_stable)**2)-gsurf(index)%Ra(1:Nseg_stable))

  gradothetax(1:Nseg_stable) = -1.0*( (oz(1:Nseg_stable)-gsurf(index)%Za(1:Nseg_stable)) * &
          (sqrt(ox(1:Nseg_stable)**2+oy(1:Nseg_stable)**2) - gsurf(index)%Ra(1:Nseg_stable))**-2.0 * &
          ( (ox(1:Nseg_stable)**2+oy(1:Nseg_stable)**2)**-.5 * ox(1:Nseg_stable) - odRadx(1:Nseg_stable) ) &
          + odZadx(1:Nseg_stable) / (sqrt(ox(1:Nseg_stable)**2+oy(1:Nseg_stable)**2) &
          - gsurf(index)%Ra(1:Nseg_stable)) ) / (tanotheta(1:Nseg_stable)**2+1)

  gradothetay(1:Nseg_stable) = -1.0*( (oz(1:Nseg_stable)-gsurf(index)%Za(1:Nseg_stable)) * &
          (sqrt(ox(1:Nseg_stable)**2+oy(1:Nseg_stable)**2) - gsurf(index)%Ra(1:Nseg_stable))**-2.0 * &
          ( (ox(1:Nseg_stable)**2+oy(1:Nseg_stable)**2)**-.5 * oy(1:Nseg_stable) - odRady(1:Nseg_stable) ) &
          + odZady(1:Nseg_stable) / (sqrt(ox(1:Nseg_stable)**2+oy(1:Nseg_stable)**2) &
          - gsurf(index)%Ra(1:Nseg_stable)) ) / (tanotheta(1:Nseg_stable)**2+1)

  gradothetaz(1:Nseg_stable) = ( (tanotheta(1:Nseg_stable)**2+1) * &
          ( sqrt(ox(1:Nseg_stable)**2+oy(1:Nseg_stable)**2) - gsurf(index)%Ra(1:Nseg_stable) ) )**-1.0

  ! Calculate field on stable field lines
  SALLOCATE( oBx, (1:Nseg_stable), 0.0 )
  SALLOCATE( oBy, (1:Nseg_stable), 0.0 )
  SALLOCATE( oBz, (1:Nseg_stable), 0.0 )

  ! Parallelized, maybe change to be parallelized elsewhere
  do i = 1, Nseg_stable
     if( myid.ne.modulo(i,ncpu) ) cycle ! parallelization loop;
     do icoil = 1, Ncoils
        call bfield0(icoil, ox(i), oy(i), oz(i), Bxhold, Byhold, Bzhold)
        oBx(i) = oBx(i) + Bxhold
        oBy(i) = oBy(i) + Byhold
        oBz(i) = oBz(i) + Bzhold
     enddo
  enddo
  call MPI_BARRIER( MPI_COMM_FOCUS, ierr )
  call MPI_ALLREDUCE( MPI_IN_PLACE, oBx, Nseg_stable, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FOCUS, ierr )
  call MPI_ALLREDUCE( MPI_IN_PLACE, oBy, Nseg_stable, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FOCUS, ierr )
  call MPI_ALLREDUCE( MPI_IN_PLACE, oBz, Nseg_stable, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FOCUS, ierr )

  ! Calculate B^s, B^zeta, B^theta on stable field lines
  obsups(1:Nseg_stable) = oBx(1:Nseg_stable)*gradosx(1:Nseg_stable) + &
          oBy(1:Nseg_stable)*gradosy(1:Nseg_stable) + oBz(1:Nseg_stable)*gradosz(1:Nseg_stable)

  obsupzeta(1:Nseg_stable) = oBx(1:Nseg_stable)*gradozetax(1:Nseg_stable) + &
          oBy(1:Nseg_stable)*gradozetay(1:Nseg_stable) + oBz(1:Nseg_stable)*gradozetaz(1:Nseg_stable)

  obsuptheta(1:Nseg_stable) = oBx(1:Nseg_stable)*gradothetax(1:Nseg_stable) + &
          oBy(1:Nseg_stable)*gradothetay(1:Nseg_stable) + oBz(1:Nseg_stable)*gradothetaz(1:Nseg_stable)

  ! DALLOCATE
  DALLOCATE( zeta )
  DALLOCATE( odRadx )
  DALLOCATE( odRady )
  DALLOCATE( odZadx )
  DALLOCATE( odZady )
  DALLOCATE( gradosx )
  DALLOCATE( gradosy )
  DALLOCATE( gradosz )
  DALLOCATE( os )
  DALLOCATE( gradozetax )
  DALLOCATE( gradozetay )
  DALLOCATE( gradozetaz )
  DALLOCATE( gradothetax )
  DALLOCATE( gradothetay )
  DALLOCATE( gradothetaz )
  DALLOCATE( tanotheta )
  DALLOCATE( oBx )
  DALLOCATE( oBy )
  DALLOCATE( oBz )

  return

end subroutine calcbsup
