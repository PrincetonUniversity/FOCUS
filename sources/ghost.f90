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

  implicit none
  include "mpif.h"
  INTEGER, INTENT(in)    :: index

  INTEGER             :: astat, ierr, i, Nseg_stable

  Nseg_stable = gsurf(index)%Nseg_stable
  do i = 1,Nseg_stable
     gsurf(index)%zeta(i) = (i-1)*pi2*resbn_m/(Nseg_stable-1)
  enddo

  call make_axis(index) ! Keep fixed, future could solve every iteration 

  ! Calculate f and g functions for o and x
  call calcfg(index)

  ! Optimize F = of^2 + xf^2 + og^2 + xg^2 functional

  ! Construct ghost surface

  return

end subroutine ghost

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine make_axis(index)
  use globals, only: dp, zero, half, pi2, ncpu, myid, resbn_m, gsurf, ounit, MPI_COMM_FOCUS

  implicit none
  include "mpif.h"
  INTEGER, INTENT(IN) :: index

  INTEGER             :: astat, ierr, i, Nseg_stable

  Nseg_stable = gsurf(index)%Nseg_stable

  gsurf(index)%Ra(1:Nseg_stable) = 0.0
  gsurf(index)%Za(1:Nseg_stable) = 0.0

  do i = 1, gsurf(index)%NF_axis
     gsurf(index)%Ra(1:Nseg_stable) = gsurf(index)%Ra(1:Nseg_stable) + &
             gsurf(index)%axisrnc(i)*cos(gsurf(index)%axisn(i)*gsurf(index)%zeta(1:Nseg_stable))
     gsurf(index)%Za(1:Nseg_stable) = gsurf(index)%Za(1:Nseg_stable) + &
             gsurf(index)%axiszns(i)*sin(gsurf(index)%axisn(i)*gsurf(index)%zeta(1:Nseg_stable))
  enddo

end subroutine make_axis

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine calcfg(index)
  use globals, only: dp, zero, half, pi2, ncpu, myid, resbn_m, resbn_n, gsurf, Ncoils, ounit, MPI_COMM_FOCUS

  implicit none
  include "mpif.h"
  INTEGER, INTENT(IN) :: index

  INTEGER             :: astat, ierr, i, Nseg_stable, p, q, icoil
  REAL                :: Bxhold, Byhold, Bzhold
  REAL, ALLOCATABLE   :: zeta(:), gradosx(:), gradosy(:), gradosz(:), gradxsx(:), gradxsy(:), &
                         gradxsz(:), gradozetax(:), gradozetay(:), gradozetaz(:), gradxzetax(:), &
                         gradxzetay(:), gradxzetaz(:), gradothetax(:), gradothetay(:), gradothetaz(:), &
                         gradxthetax(:), gradxthetay(:), gradxthetaz(:), oblah(:), xblah(:), oBx(:), &
                         oBy(:), oBz(:), xBx(:), xBy(:), xBz(:), obsups(:), obsupzeta(:), &
                         obsuptheta(:), xbsups(:), xbsupzeta(:), xbsuptheta(:)

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

  ! Calculate xdot, etc for dl. Use fd for now
  do i = 1, Nseg_stable-1
     gsurf(index)%oxdot(i) = gsurf(index)%ox(i+1) - gsurf(index)%ox(i)
     gsurf(index)%xxdot(i) = gsurf(index)%xx(i+1) - gsurf(index)%xx(i)
     gsurf(index)%oydot(i) = gsurf(index)%oy(i+1) - gsurf(index)%oy(i)
     gsurf(index)%xydot(i) = gsurf(index)%xy(i+1) - gsurf(index)%xy(i)
     gsurf(index)%ozdot(i) = gsurf(index)%oz(i+1) - gsurf(index)%oz(i)
     gsurf(index)%xzdot(i) = gsurf(index)%xz(i+1) - gsurf(index)%xz(i)
  enddo
  gsurf(index)%oxdot(Nseg_stable) = gsurf(index)%oxdot(1)
  gsurf(index)%xxdot(Nseg_stable) = gsurf(index)%xxdot(1)
  gsurf(index)%oydot(Nseg_stable) = gsurf(index)%oydot(1)
  gsurf(index)%xydot(Nseg_stable) = gsurf(index)%xydot(1)
  gsurf(index)%ozdot(Nseg_stable) = gsurf(index)%ozdot(1)
  gsurf(index)%xzdot(Nseg_stable) = gsurf(index)%xzdot(1)
  
  gsurf(index)%oxdot(1:Nseg_stable) = (Nseg_stable-1)*gsurf(index)%oxdot(1:Nseg_stable)/(pi2*resbn_m)
  gsurf(index)%xxdot(1:Nseg_stable) = (Nseg_stable-1)*gsurf(index)%xxdot(1:Nseg_stable)/(pi2*resbn_m)
  gsurf(index)%oydot(1:Nseg_stable) = (Nseg_stable-1)*gsurf(index)%oydot(1:Nseg_stable)/(pi2*resbn_m)
  gsurf(index)%xydot(1:Nseg_stable) = (Nseg_stable-1)*gsurf(index)%xydot(1:Nseg_stable)/(pi2*resbn_m)
  gsurf(index)%ozdot(1:Nseg_stable) = (Nseg_stable-1)*gsurf(index)%ozdot(1:Nseg_stable)/(pi2*resbn_m)
  gsurf(index)%xzdot(1:Nseg_stable) = (Nseg_stable-1)*gsurf(index)%xzdot(1:Nseg_stable)/(pi2*resbn_m)
 

  ! Calculate contravariant basis, make variables local, MATH IS INCORRECT
  SALLOCATE( gradosx, (1:Nseg_stable), 0.0 )
  SALLOCATE( gradosy, (1:Nseg_stable), 0.0 )
  SALLOCATE( gradosz, (1:Nseg_stable), 0.0 )
  SALLOCATE( gradxsx, (1:Nseg_stable), 0.0 )
  SALLOCATE( gradxsy, (1:Nseg_stable), 0.0 )
  SALLOCATE( gradxsz, (1:Nseg_stable), 0.0 )
  
  gradosx(1:Nseg_stable) = ( sqrt( gsurf(index)%ox(1:Nseg_stable)**2 + gsurf(index)%oy(1:Nseg_stable)**2 ) - &
          gsurf(index)%Ra(1:Nseg_stable) )*(gsurf(index)%ox(1:Nseg_stable)**2 + &
          gsurf(index)%oy(1:Nseg_stable)**2)**(-.5)*gsurf(index)%ox(1:Nseg_stable)/gsurf(index)%os(1:Nseg_stable)
  gradxsx(1:Nseg_stable) = ( sqrt( gsurf(index)%xx(1:Nseg_stable)**2 + gsurf(index)%xy(1:Nseg_stable)**2 ) - &
          gsurf(index)%Ra(1:Nseg_stable) )*(gsurf(index)%xx(1:Nseg_stable)**2 + &
          gsurf(index)%xy(1:Nseg_stable)**2)**(-.5)*gsurf(index)%xx(1:Nseg_stable)/gsurf(index)%xs(1:Nseg_stable)
  
  gradosy(1:Nseg_stable) = ( sqrt( gsurf(index)%ox(1:Nseg_stable)**2 + gsurf(index)%oy(1:Nseg_stable)**2 ) - &
          gsurf(index)%Ra(1:Nseg_stable) )*(gsurf(index)%ox(1:Nseg_stable)**2 + &
          gsurf(index)%oy(1:Nseg_stable)**2)**(-.5)*gsurf(index)%oy(1:Nseg_stable)/gsurf(index)%os(1:Nseg_stable)
  gradxsy(1:Nseg_stable) = ( sqrt( gsurf(index)%xx(1:Nseg_stable)**2 + gsurf(index)%xy(1:Nseg_stable)**2 ) - &
          gsurf(index)%Ra(1:Nseg_stable) )*(gsurf(index)%xx(1:Nseg_stable)**2 + &
          gsurf(index)%xy(1:Nseg_stable)**2)**(-.5)*gsurf(index)%xy(1:Nseg_stable)/gsurf(index)%xs(1:Nseg_stable)
  
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
  SALLOCATE( oblah, (1:Nseg_stable), 0.0 )
  SALLOCATE( xblah, (1:Nseg_stable), 0.0 )

  oblah(1:Nseg_stable) = (gsurf(index)%oz(1:Nseg_stable) - gsurf(index)%Za(1:Nseg_stable)) & 
          /( sqrt(gsurf(index)%ox(1:Nseg_stable)**2+gsurf(index)%oy(1:Nseg_stable)**2) - gsurf(index)%Ra(1:Nseg_stable) )
  xblah(1:Nseg_stable) = (gsurf(index)%xz(1:Nseg_stable) - gsurf(index)%Za(1:Nseg_stable)) &
          /( sqrt(gsurf(index)%xx(1:Nseg_stable)**2+gsurf(index)%xy(1:Nseg_stable)**2) - gsurf(index)%Ra(1:Nseg_stable) )

  gradothetax(1:Nseg_stable) = -1.0*( oblah(1:Nseg_stable)**2 + 1 )**-1.0 &
          * oblah(1:Nseg_stable)**2 * (gsurf(index)%ox(1:Nseg_stable)**2+gsurf(index)%oy(1:Nseg_stable)**2)**-.5 &
          * gsurf(index)%ox(1:Nseg_stable)/(gsurf(index)%oz(1:Nseg_stable)-gsurf(index)%Za(1:Nseg_stable))
  gradxthetax(1:Nseg_stable) = -1.0*( xblah(1:Nseg_stable)**2 + 1 )**-1.0 &
          * xblah(1:Nseg_stable)**2 * (gsurf(index)%xx(1:Nseg_stable)**2+gsurf(index)%xy(1:Nseg_stable)**2)**-.5 &
          * gsurf(index)%xx(1:Nseg_stable)/(gsurf(index)%xz(1:Nseg_stable)-gsurf(index)%Za(1:Nseg_stable))

  gradothetay(1:Nseg_stable) = -1.0*( oblah(1:Nseg_stable)**2 + 1 )**-1.0 &
          * oblah(1:Nseg_stable)**2 * (gsurf(index)%ox(1:Nseg_stable)**2+gsurf(index)%oy(1:Nseg_stable)**2)**-.5 &
          * gsurf(index)%oy(1:Nseg_stable)/(gsurf(index)%oz(1:Nseg_stable)-gsurf(index)%Za(1:Nseg_stable))
  gradxthetay(1:Nseg_stable) = -1.0*( xblah(1:Nseg_stable)**2 + 1 )**-1.0 &
          * xblah(1:Nseg_stable)**2 * (gsurf(index)%xx(1:Nseg_stable)**2+gsurf(index)%xy(1:Nseg_stable)**2)**-.5 &
          * gsurf(index)%xy(1:Nseg_stable)/(gsurf(index)%xz(1:Nseg_stable)-gsurf(index)%Za(1:Nseg_stable))

  gradothetaz(1:Nseg_stable) = ( oblah(1:Nseg_stable)**2 + 1 )**-1.0 &
          * oblah(1:Nseg_stable) / (gsurf(index)%oz(1:Nseg_stable)-gsurf(index)%Za(1:Nseg_stable))
  gradxthetaz(1:Nseg_stable) = ( xblah(1:Nseg_stable)**2 + 1 )**-1.0 &
          * xblah(1:Nseg_stable) / (gsurf(index)%xz(1:Nseg_stable)-gsurf(index)%Za(1:Nseg_stable))

  ! Dont need to calculate Jacobian since not calculating covariant basis? 
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

  do i = 1, Nseg_stable
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

   ! Calculate B^s, B^zeta, B^theta on stable field lines
   SALLOCATE( obsups, (1:Nseg_stable), 0.0 )
   SALLOCATE( obsupzeta, (1:Nseg_stable), 0.0 )
   SALLOCATE( obsuptheta, (1:Nseg_stable), 0.0 )
   SALLOCATE( xbsups, (1:Nseg_stable), 0.0 )
   SALLOCATE( xbsupzeta, (1:Nseg_stable), 0.0 )
   SALLOCATE( xbsuptheta, (1:Nseg_stable), 0.0 )

   obsups(1:Nseg_stable) = oBx(1:Nseg_stable)*gradosx(1:Nseg_stable) + &
           oBy(1:Nseg_stable)*gradosy(1:Nseg_stable) + oBz(1:Nseg_stable)*gradosz(1:Nseg_stable)
   xbsups(1:Nseg_stable) = xBx(1:Nseg_stable)*gradxsx(1:Nseg_stable) + &
           xBy(1:Nseg_stable)*gradxsy(1:Nseg_stable) + xBz(1:Nseg_stable)*gradxsz(1:Nseg_stable)

   obsupzeta(1:Nseg_stable) = oBx(1:Nseg_stable)*gradozetax(1:Nseg_stable) + &
           oBy(1:Nseg_stable)*gradozetay(1:Nseg_stable) + oBz(1:Nseg_stable)*gradozetaz(1:Nseg_stable)
   xbsupzeta(1:Nseg_stable) = xBx(1:Nseg_stable)*gradxzetax(1:Nseg_stable) + &
           xBy(1:Nseg_stable)*gradxzetay(1:Nseg_stable) + xBz(1:Nseg_stable)*gradxzetaz(1:Nseg_stable)

   obsuptheta(1:Nseg_stable) = oBx(1:Nseg_stable)*gradothetax(1:Nseg_stable) + &
           oBy(1:Nseg_stable)*gradothetay(1:Nseg_stable) + oBz(1:Nseg_stable)*gradothetaz(1:Nseg_stable)
   xbsuptheta(1:Nseg_stable) = xBx(1:Nseg_stable)*gradxthetax(1:Nseg_stable) + &
           xBy(1:Nseg_stable)*gradxthetay(1:Nseg_stable) + xBz(1:Nseg_stable)*gradxthetaz(1:Nseg_stable)

  ! Calculate f,g and F
  gsurf(index)%of(1:Nseg_stable) = obsups(1:Nseg_stable)/obsupzeta(1:Nseg_stable) - &
          gsurf(index)%osdot(1:Nseg_stable)
  gsurf(index)%xf(1:Nseg_stable) = xbsups(1:Nseg_stable)/xbsupzeta(1:Nseg_stable) - &
          gsurf(index)%xsdot(1:Nseg_stable)

  gsurf(index)%og(1:Nseg_stable) = obsuptheta(1:Nseg_stable)/obsupzeta(1:Nseg_stable) - &
          gsurf(index)%othetadot(1:Nseg_stable)
  gsurf(index)%xg(1:Nseg_stable) = xbsuptheta(1:Nseg_stable)/xbsupzeta(1:Nseg_stable) - &
          gsurf(index)%xthetadot(1:Nseg_stable)

  gsurf(index)%F = 0.0
  do i = 1, Nseg_stable-1
     gsurf(index)%F = gsurf(index)%F + gsurf(index)%of(i)**2 + gsurf(index)%og(i)**2 + &
             gsurf(index)%xf(i)**2 + gsurf(index)%xg(i)**2
  enddo

  ! Ignore multiplying by 2 pi q 
  gsurf(index)%F = gsurf(index)%F/(Nseg_stable-1)

  ! DALLOCATE
  DALLOCATE( zeta )
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
  DALLOCATE( gradothetax )
  DALLOCATE( gradothetay )
  DALLOCATE( gradothetaz )
  DALLOCATE( gradxthetax )
  DALLOCATE( gradxthetay )
  DALLOCATE( gradxthetaz )
  DALLOCATE( oblah )
  DALLOCATE( xblah )
  DALLOCATE( oBx )
  DALLOCATE( oBy )
  DALLOCATE( oBz )
  DALLOCATE( xBx )
  DALLOCATE( xBy )
  DALLOCATE( xBz )
  DALLOCATE( obsups )
  DALLOCATE( obsupzeta )
  DALLOCATE( obsuptheta )
  DALLOCATE( xbsups )
  DALLOCATE( xbsupzeta )
  DALLOCATE( xbsuptheta )

  return

end subroutine calcfg







