!title (ghost) ! Calculate ghost/quadratic flux minimizing surfaces. (tkruger)

!latex \briefly{This function calculates a ghost surface to be used in the 
!latex         resonant Fourier harmonic calculation. It is important that we use
!latex         these ghost surfaces instead of a rational surface from VMEC. Rational
!latex         surfaces from VMEC will inherently have some error associated with 
!latex         their construction and will cause errors in the resonant Fouer 
!latex         harmonic calculation and moreover its sensitivity.}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! not parallelized; communications may take more time;
subroutine ghost(s,theta,zeta,bsups,bsuptheta,bsupzeta)
  use globals, only: dp, zero, half, pi2, machprec, ncpu, myid, ounit, &
       coil, DoF, Ncoils, Nfixgeo, Ndof, LM_fvec, LM_fjac, R0, &
       MPI_COMM_FOCUS

  implicit none
  include "mpif.h"
  REAL, INTENT(in)    :: s, theta, zeta
  REAL, INTENT(out)   :: bsups, bsuptheta, bsupzeta

  INTEGER             :: icoil
  REAL                :: x, y, z, gradsx, gradsy, gradsz, gradthetax, & 
                         gradthetay, gradthetaz, gradzetax, gradzetay, &
                         gradzetaz, blah, J, Bx, By, Bz, Bxicoil, Byicoil, &
                         Bzicoil

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  x =    cos(zeta)*(R0 + s*cos(theta))
  y = -1*sin(zeta)*(R0 + s*cos(theta))
  z = s*sin(theta)

  gradsx = 2*( sqrt(x**2+y**2) - R0 )*(x**2+y**2)**(-.5)*x
  gradsy = 2*( sqrt(x**2+y**2) - R0 )*(x**2+y**2)**(-.5)*y
  gradsz = 2*z

  blah = z/( sqrt(x**2+y**2) - R0 )
  gradthetax = -1.0*( blah**2 + 1 )**-1.0 * blah**2 * (x**2+y**2)**-.5 * x/z
  gradthetay = -1.0*( blah**2 + 1 )**-1.0 * blah**2 * (x**2+y**2)**-.5 * y/z
  gradthetaz =      ( blah**2 + 1 )**-1.0 * blah/z

  gradzetax = (y**2/x**2 + 1)**-1.0 * y/x**2
  gradzetay = (y**2/x**2 + 1)**-1.0 * -1.0 / x
  gradzetaz = 0.0

  J =     gradsx*( gradthetay*gradzetaz - gradthetaz*gradzetay )
  J = J + gradsy*( gradthetaz*gradzetax - gradthetax*gradzetaz )
  J = J + gradsz*( gradthetax*gradzetay - gradthetay*gradzetax )
  J = J**-1.0

  Bx = 0.0
  By = 0.0
  Bz = 0.0

  do icoil = 1, Ncoils
     call bfield0(icoil, x, y, z, Bxicoil, Byicoil, Bzicoil)
     Bx = Bx + Bxicoil
     By = By + Byicoil
     Bz = Bz + Bzicoil
  enddo

  bsups     = Bx*gradsx     + By*gradsy     + Bz*gradsz
  bsuptheta = Bx*gradthetax + By*gradthetay + Bz*gradthetaz
  bsupzeta  = Bx*gradzetax  + By*gradzetay  + Bz*gradzetaz

  return

end subroutine ghost


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
