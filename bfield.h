
!title (field) ! Computes magnetic field.

!latex \briefly{Computes magnetic field given coil geometry and currents.}

!latex \calledby{\link{}}
!latex \calls{\oculus{bs00aa}}

!latex \tableofcontents

!latex \subsection{magnetic field}
!latex \bi
!latex \item See \oculus{bs00aa}.
!latex \ei


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!!$subroutine bfield( RpZ, itangent, dBRpZ, ibfield ) 
!!$
!!$! DATE : April 2016
!!$! using bs00aa in oculus ( actually NAG adaptive integration routines ) calculating magnetic field;
!!$! only be used at the beginning stage; may be deleted later.
!!$
!!$!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!!$  
!!$  use kmodule, only : zero, one, pi2, vsmall, myid, ounit, Ncoils, coil, icoil, bsfield
!!$  
!!$  use oculus , only : bs00aa
!!$  
!!$!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!!$  
!!$  implicit none
!!$  
!!$  include "mpif.h"
!!$  
!!$  INTEGER            :: itangent, ibfield
!!$  REAL               :: RpZ(1:3), dBRpZ(1:3,0:3)
!!$  
!!$  INTEGER            :: ierr, ibs00aa
!!$  REAL               :: zeta, Bx, By, Bz, czeta, szeta
!!$  
!!$!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!!$  
!!$  dBRpZ(1:3,0:3) = zero ! set default intent(out) ; 11 Oct 15;
!!$  
!!$!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!!$  
!!$  zeta = modulo( RpZ(2), pi2)
!!$  
!!$  czeta = cos(zeta)
!!$  szeta = sin(zeta)
!!$  
!!$  bsfield%x = RpZ(1) * czeta
!!$  bsfield%y = RpZ(1) * szeta
!!$  bsfield%z = RpZ(3)
!!$  
!!$  Bx = zero
!!$  By = zero
!!$  Bz = zero
!!$  
!!$  do icoil = 1, Ncoils ! icoil is a global variable which is passed through to auxiliary.h; 11 Oct 15;
!!$   
!!$   ibs00aa =  0 ; bsfield%LB = .true. ; bsfield%LA = .false. ; bsfield%LL = .false. ; call bs00aa( bsfield, ibs00aa )
!!$   
!!$   Bx = Bx + bsfield%Bx * coil(icoil)%I
!!$   By = By + bsfield%By * coil(icoil)%I
!!$   Bz = Bz + bsfield%Bz * coil(icoil)%I
!!$   
!!$  enddo ! end of do icoil; 11 Oct 15;
!!$  
!!$  dBRpZ(1,0) = (   Bx * czeta + By * szeta )
!!$  dBRpZ(2,0) = ( - Bx * szeta + By * czeta ) / RpZ(1)
!!$  dBRpZ(3,0) =     Bz
!!$  
!!$  ibfield = 0
!!$  
!!$  return
!!$  
!!$end subroutine bfield

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
! DATE:  06/15/2016
! calculate the magnetic field of icoil using manually discretized coils. 
! Biot-Savart constant and integral factor are not included for later simplication. Be careful if coils have different resolutions.

subroutine bfield0(icoil, iteta, jzeta, Bx, By, Bz)
  
  use kmodule, only: coil, surf, NDcoil, Ncoils, Nteta, Nzeta, &
                     zero, myid, ounit
  implicit none
  include "mpif.h"

  INTEGER, intent(in ) :: icoil, iteta, jzeta
  REAL   , intent(out) :: Bx, By, Bz

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  INTEGER              :: ierr, astat, kseg
  REAL                 :: dlx, dly, dlz, rm3, ltx, lty, ltz

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  FATAL( bfield0, icoil .lt. 1 .or. icoil .gt. Ncoils, icoil not in right range )
  FATAL( bfield0, iteta .lt. 0 .or. iteta .gt. Nteta , iteta not in right range )
  FATAL( bfield0, jzeta .lt. 0 .or. jzeta .gt. Nzeta , jzeta not in right range )
  
  dlx = zero; ltx = zero; Bx = zero
  dly = zero; lty = zero; By = zero
  dlz = zero; ltz = zero; Bz = zero

  do kseg = 0, coil(icoil)%D - 1       ! NDcoil can be changed to coil(icoil)%D; 06/15/2016
        
   dlx = coil(icoil)%xx(kseg) - surf(1)%xx(iteta,jzeta)
   dly = coil(icoil)%yy(kseg) - surf(1)%yy(iteta,jzeta)
   dlz = coil(icoil)%zz(kseg) - surf(1)%zz(iteta,jzeta)
   rm3 = (sqrt(dlx**2 + dly**2 + dlz**2))**(-3)

   ltx = coil(icoil)%xt(kseg)
   lty = coil(icoil)%yt(kseg)
   ltz = coil(icoil)%zt(kseg)
   
   Bx = Bx + ( dly*ltz - dlz*lty ) * rm3
   By = By + ( dlz*ltx - dlx*ltz ) * rm3
   Bz = Bz + ( dlx*lty - dly*ltx ) * rm3

  enddo    ! enddo kseg

  return

end subroutine bfield0

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
! DATE:  06/15/2016
! calculate the magnetic field and its first dirivatives of icoil using manually discretized coils data. 
! Biot-Savart constant and integral factor are not included for later simplication. Be careful if coils have different resolutions.

subroutine bfield1(icoil, iteta, jzeta, Bx, By, Bz)
  
  use kmodule, only: coil, surf, NFcoil, NDcoil, Ncoils, Cdof, Nteta, Nzeta, cmt, smt, &
                     zero, myid, ounit
  implicit none
  include "mpif.h"

  INTEGER, intent(in ) :: icoil, iteta, jzeta
  REAL   , intent(out) :: Bx(0:Cdof), By(0:Cdof), Bz(0:Cdof)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  INTEGER              :: ierr, astat, kseg, NN, ll
  REAL                 :: dlx, dly, dlz, r, rm3, rm5, ltx, lty, ltz

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  FATAL( bfield1, icoil .lt. 1 .or. icoil .gt. Ncoils, icoil not in right range )
  FATAL( bfield1, iteta .lt. 0 .or. iteta .gt. Nteta , iteta not in right range )
  FATAL( bfield1, jzeta .lt. 0 .or. jzeta .gt. Nzeta , jzeta not in right range )
  
  NN = NFcoil
  dlx = zero; ltx = zero; Bx = zero
  dly = zero; lty = zero; By = zero
  dlz = zero; ltz = zero; Bz = zero

  do kseg = 0, NDcoil - 1       ! NDcoil can be changed to coil(icoil)%D; 06/15/2016
        
   dlx = coil(icoil)%xx(kseg) - surf(1)%xx(iteta,jzeta)
   dly = coil(icoil)%yy(kseg) - surf(1)%yy(iteta,jzeta)
   dlz = coil(icoil)%zz(kseg) - surf(1)%zz(iteta,jzeta)
   r = sqrt(dlx**2 + dly**2 + dlz**2); rm3 = r**(-3);  rm5 = r**(-5)

   ltx = coil(icoil)%xt(kseg)
   lty = coil(icoil)%yt(kseg)
   ltz = coil(icoil)%zt(kseg)
   
   Bx(0) = Bx(0) + ( dly*ltz - dlz*lty ) * rm3
   By(0) = By(0) + ( dlz*ltx - dlx*ltz ) * rm3
   Bz(0) = Bz(0) + ( dlx*lty - dly*ltx ) * rm3

   do ll = 0, NN

    Bx(ll     +1) = Bx(ll     +1)                                                 - 3*(dly*ltz-dlz*lty)*dlx*cmt(kseg,ll)*rm5 !Bx/xc
    Bx(ll+  NN+2) = Bx(ll+  NN+2)                                                 - 3*(dly*ltz-dlz*lty)*dlx*smt(kseg,ll)*rm5 !Bx/xs
    Bx(ll+2*NN+3) = Bx(ll+2*NN+3) + ( cmt(kseg,ll)*ltz + ll*smt(kseg,ll)*dlz)*rm3 - 3*(dly*ltz-dlz*lty)*dly*cmt(kseg,ll)*rm5 !Bx/yc
    Bx(ll+3*NN+4) = Bx(ll+3*NN+4) + ( smt(kseg,ll)*ltz - ll*cmt(kseg,ll)*dlz)*rm3 - 3*(dly*ltz-dlz*lty)*dly*smt(kseg,ll)*rm5 !Bx/ys
    Bx(ll+4*NN+5) = Bx(ll+4*NN+5) + (-cmt(kseg,ll)*lty - ll*smt(kseg,ll)*dly)*rm3 - 3*(dly*ltz-dlz*lty)*dlz*cmt(kseg,ll)*rm5 !Bx/zc
    Bx(ll+5*NN+6) = Bx(ll+5*NN+6) + (-smt(kseg,ll)*lty + ll*cmt(kseg,ll)*dly)*rm3 - 3*(dly*ltz-dlz*lty)*dlz*smt(kseg,ll)*rm5 !Bx/zs

    By(ll     +1) = By(ll     +1) + (-cmt(kseg,ll)*ltz - ll*smt(kseg,ll)*dlz)*rm3 - 3*(dlz*ltx-dlx*ltz)*dlx*cmt(kseg,ll)*rm5 !By/xc
    By(ll+  NN+2) = By(ll+  NN+2) + (-smt(kseg,ll)*ltz + ll*cmt(kseg,ll)*dlz)*rm3 - 3*(dlz*ltx-dlx*ltz)*dlx*smt(kseg,ll)*rm5 !By/xs
    By(ll+2*NN+3) = By(ll+2*NN+3)                                                 - 3*(dlz*ltx-dlx*ltz)*dly*cmt(kseg,ll)*rm5 !By/yc
    By(ll+3*NN+4) = By(ll+3*NN+4)                                                 - 3*(dlz*ltx-dlx*ltz)*dly*smt(kseg,ll)*rm5 !By/ys
    By(ll+4*NN+5) = By(ll+4*NN+5) + ( cmt(kseg,ll)*ltx + ll*smt(kseg,ll)*dlx)*rm3 - 3*(dlz*ltx-dlx*ltz)*dlz*cmt(kseg,ll)*rm5 !By/zc
    By(ll+5*NN+6) = By(ll+5*NN+6) + ( smt(kseg,ll)*ltx - ll*cmt(kseg,ll)*dlx)*rm3 - 3*(dlz*ltx-dlx*ltz)*dlz*smt(kseg,ll)*rm5 !By/zs

    Bz(ll     +1) = Bz(ll     +1) + ( cmt(kseg,ll)*lty + ll*smt(kseg,ll)*dly)*rm3 - 3*(dlx*lty-dly*ltx)*dlx*cmt(kseg,ll)*rm5 !Bz/xc
    Bz(ll+  NN+2) = Bz(ll+  NN+2) + ( smt(kseg,ll)*lty - ll*cmt(kseg,ll)*dly)*rm3 - 3*(dlx*lty-dly*ltx)*dlx*smt(kseg,ll)*rm5 !Bz/xs
    Bz(ll+2*NN+3) = Bz(ll+2*NN+3) + (-cmt(kseg,ll)*ltx - ll*smt(kseg,ll)*dlx)*rm3 - 3*(dlx*lty-dly*ltx)*dly*cmt(kseg,ll)*rm5 !Bz/yc
    Bz(ll+3*NN+4) = Bz(ll+3*NN+4) + (-smt(kseg,ll)*ltx + ll*cmt(kseg,ll)*dlx)*rm3 - 3*(dlx*lty-dly*ltx)*dly*smt(kseg,ll)*rm5 !Bz/ys
    Bz(ll+4*NN+5) = Bz(ll+4*NN+5)                                                 - 3*(dlx*lty-dly*ltx)*dlz*cmt(kseg,ll)*rm5 !Bz/zc
    Bz(ll+5*NN+6) = Bz(ll+5*NN+6)                                                 - 3*(dlx*lty-dly*ltx)*dlz*smt(kseg,ll)*rm5 !Bz/zs

   enddo ! enddo ll
   
  enddo    ! enddo kseg

  return

end subroutine bfield1

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
! DATE:  06/15/2016
! calculate the magnetic field, its first dirivatives and its second dirivatives of icoil using manually discretized coils data. 
! Biot-Savart constant and integral factor are not included for later simplication. Be careful if coils have different resolutions.

subroutine bfield2(icoil, iteta, jzeta, Bx, By, Bz)
  
  use kmodule, only: coil, surf, NFcoil, NDcoil, Ncoils, Cdof, Nteta, Nzeta, cmt, smt, &
                     zero, myid, ounit
  implicit none
  include "mpif.h"

  INTEGER, intent(in ) :: icoil, iteta, jzeta
  REAL   , intent(out) :: Bx(0:Cdof,0:Cdof), By(0:Cdof,0:Cdof), Bz(0:Cdof,0:Cdof)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  INTEGER              :: ierr, astat, kseg, NN, ll, mm
  REAL                 :: dlx, dly, dlz, r, rm3, rm5, rm7, ltx, lty, ltz, lBx, lBy, lBz

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  FATAL( bfield2, icoil .lt. 1 .or. icoil .gt. Ncoils, icoil not in right range )
  FATAL( bfield2, iteta .lt. 0 .or. iteta .gt. Nteta , iteta not in right range )
  FATAL( bfield2, jzeta .lt. 0 .or. jzeta .gt. Nzeta , jzeta not in right range )
  
  NN = NFcoil
  dlx = zero; ltx = zero; lBx = zero; Bx = zero
  dly = zero; lty = zero; lBy = zero; By = zero
  dlz = zero; ltz = zero; lBz = zero; Bz = zero

  do kseg = 0, NDcoil - 1       ! NDcoil can be changed to coil(icoil)%D; 06/15/2016
        
   dlx = coil(icoil)%xx(kseg) - surf(1)%xx(iteta,jzeta)
   dly = coil(icoil)%yy(kseg) - surf(1)%yy(iteta,jzeta)
   dlz = coil(icoil)%zz(kseg) - surf(1)%zz(iteta,jzeta)
   r = sqrt(dlx**2 + dly**2 + dlz**2); rm3 = r**(-3);  rm5 = r**(-5); rm7 = r**(-7)

   ltx = coil(icoil)%xt(kseg)
   lty = coil(icoil)%yt(kseg)
   ltz = coil(icoil)%zt(kseg)

   lBx = dly*ltz - dlz*lty
   lBy = dlz*ltx - dlx*ltz
   lBz = dlx*lty - dly*ltx
   
   Bx(0,0) = Bx(0,0) + lBx * rm3
   By(0,0) = By(0,0) + lBy * rm3
   Bz(0,0) = Bz(0,0) + lBz * rm3

   do ll = 0, NN

    Bx(ll     +1,0) = Bx(ll     +1,0)                                                 - 3*lBx*dlx*cmt(kseg,ll)*rm5 !Bx/xc
    Bx(ll+  NN+2,0) = Bx(ll+  NN+2,0)                                                 - 3*lBx*dlx*smt(kseg,ll)*rm5 !Bx/xs
    Bx(ll+2*NN+3,0) = Bx(ll+2*NN+3,0) + ( cmt(kseg,ll)*ltz + ll*smt(kseg,ll)*dlz)*rm3 - 3*lBx*dly*cmt(kseg,ll)*rm5 !Bx/yc
    Bx(ll+3*NN+4,0) = Bx(ll+3*NN+4,0) + ( smt(kseg,ll)*ltz - ll*cmt(kseg,ll)*dlz)*rm3 - 3*lBx*dly*smt(kseg,ll)*rm5 !Bx/ys
    Bx(ll+4*NN+5,0) = Bx(ll+4*NN+5,0) + (-cmt(kseg,ll)*lty - ll*smt(kseg,ll)*dly)*rm3 - 3*lBx*dlz*cmt(kseg,ll)*rm5 !Bx/zc
    Bx(ll+5*NN+6,0) = Bx(ll+5*NN+6,0) + (-smt(kseg,ll)*lty + ll*cmt(kseg,ll)*dly)*rm3 - 3*lBx*dlz*smt(kseg,ll)*rm5 !Bx/zs

    By(ll     +1,0) = By(ll     +1,0) + (-cmt(kseg,ll)*ltz - ll*smt(kseg,ll)*dlz)*rm3 - 3*lBy*dlx*cmt(kseg,ll)*rm5 !By/xc
    By(ll+  NN+2,0) = By(ll+  NN+2,0) + (-smt(kseg,ll)*ltz + ll*cmt(kseg,ll)*dlz)*rm3 - 3*lBy*dlx*smt(kseg,ll)*rm5 !By/xs
    By(ll+2*NN+3,0) = By(ll+2*NN+3,0)                                                 - 3*lBy*dly*cmt(kseg,ll)*rm5 !By/yc
    By(ll+3*NN+4,0) = By(ll+3*NN+4,0)                                                 - 3*lBy*dly*smt(kseg,ll)*rm5 !By/ys
    By(ll+4*NN+5,0) = By(ll+4*NN+5,0) + ( cmt(kseg,ll)*ltx + ll*smt(kseg,ll)*dlx)*rm3 - 3*lBy*dlz*cmt(kseg,ll)*rm5 !By/zc
    By(ll+5*NN+6,0) = By(ll+5*NN+6,0) + ( smt(kseg,ll)*ltx - ll*cmt(kseg,ll)*dlx)*rm3 - 3*lBy*dlz*smt(kseg,ll)*rm5 !By/zs

    Bz(ll     +1,0) = Bz(ll     +1,0) + ( cmt(kseg,ll)*lty + ll*smt(kseg,ll)*dly)*rm3 - 3*lBz*dlx*cmt(kseg,ll)*rm5 !Bz/xc
    Bz(ll+  NN+2,0) = Bz(ll+  NN+2,0) + ( smt(kseg,ll)*lty - ll*cmt(kseg,ll)*dly)*rm3 - 3*lBz*dlx*smt(kseg,ll)*rm5 !Bz/xs
    Bz(ll+2*NN+3,0) = Bz(ll+2*NN+3,0) + (-cmt(kseg,ll)*ltx - ll*smt(kseg,ll)*dlx)*rm3 - 3*lBz*dly*cmt(kseg,ll)*rm5 !Bz/yc
    Bz(ll+3*NN+4,0) = Bz(ll+3*NN+4,0) + (-smt(kseg,ll)*ltx + ll*cmt(kseg,ll)*dlx)*rm3 - 3*lBz*dly*smt(kseg,ll)*rm5 !Bz/ys
    Bz(ll+4*NN+5,0) = Bz(ll+4*NN+5,0)                                                 - 3*lBz*dlz*cmt(kseg,ll)*rm5 !Bz/zc
    Bz(ll+5*NN+6,0) = Bz(ll+5*NN+6,0)                                                 - 3*lBz*dlz*smt(kseg,ll)*rm5 !Bz/zs

   enddo ! enddo ll

   ll = 0; mm = 0

   do mm = 0, NN
    do ll = 0, NN

     !------------------------------------------------------------------Bx on xc--------------------------------------------------------------------------

     Bx(ll     +1,mm     +1) = Bx(ll     +1,mm     +1) - 3*lBx*cmt(kseg,ll)*cmt(kseg,mm)*rm5 + 15*lBx*dlx*dlx*cmt(kseg,ll)*cmt(kseg,mm)*rm7    !Bx/xc / xc

     Bx(ll+  NN+2,mm     +1) = Bx(ll+  NN+2,mm     +1) - 3*lBx*smt(kseg,ll)*cmt(kseg,mm)*rm5 + 15*lBx*dlx*dlx*smt(kseg,ll)*cmt(kseg,mm)*rm7    !Bx/xs / xc

     Bx(ll+2*NN+3,mm     +1) = Bx(ll+2*NN+3,mm     +1) -                  0                  + 15*lBx*dlx*dly*cmt(kseg,ll)*cmt(kseg,mm)*rm7  & 
                                                       - 3*( cmt(kseg,ll)*ltz + ll*smt(kseg,ll)*dlz)*dlx*cmt(kseg,mm)*rm5                      !Bx/yc / xc

     Bx(ll+3*NN+4,mm     +1) = Bx(ll+3*NN+4,mm     +1) -                  0                  + 15*lBx*dlx*dly*smt(kseg,ll)*cmt(kseg,mm)*rm7  & 
                                                       - 3*( smt(kseg,ll)*ltz - ll*cmt(kseg,ll)*dlz)*dlx*cmt(kseg,mm)*rm5                      !Bx/ys / xc

     Bx(ll+4*NN+5,mm     +1) = Bx(ll+4*NN+5,mm     +1) -                  0                  + 15*lBx*dlx*dlz*cmt(kseg,ll)*cmt(kseg,mm)*rm7  &
                                                       - 3*(-cmt(kseg,ll)*lty - ll*smt(kseg,ll)*dly)*dlx*cmt(kseg,mm)*rm5                      !Bx/zc / xc

     Bx(ll+5*NN+6,mm     +1) = Bx(ll+5*NN+6,mm     +1) -                  0                  + 15*lBx*dlx*dlz*smt(kseg,ll)*cmt(kseg,mm)*rm7  &
                                                       - 3*(-smt(kseg,ll)*lty + ll*cmt(kseg,ll)*dly)*dlx*cmt(kseg,mm)*rm5                      !Bx/zs / xc
     !------------------------------------------------------------------Bx on xs--------------------------------------------------------------------------

     Bx(ll+  NN+2,mm+  NN+2) = Bx(ll+  NN+2,mm+  NN+2) - 3*lBx*smt(kseg,ll)*smt(kseg,mm)*rm5 + 15*lBx*dlx*dlx*smt(kseg,ll)*smt(kseg,mm)*rm7    !Bx/xs / xs

     Bx(ll+2*NN+3,mm+  NN+2) = Bx(ll+2*NN+3,mm+  NN+2) -                  0                  + 15*lBx*dlx*dly*cmt(kseg,ll)*smt(kseg,mm)*rm7  & 
                                                       - 3*( cmt(kseg,ll)*ltz + ll*smt(kseg,ll)*dlz)*dlx*smt(kseg,mm)*rm5                      !Bx/yc / xs

     Bx(ll+3*NN+4,mm+  NN+2) = Bx(ll+3*NN+4,mm+  NN+2) -                  0                  + 15*lBx*dlx*dly*smt(kseg,ll)*smt(kseg,mm)*rm7  & 
                                                       - 3*( smt(kseg,ll)*ltz - ll*cmt(kseg,ll)*dlz)*dlx*smt(kseg,mm)*rm5                      !Bx/ys / xs

     Bx(ll+4*NN+5,mm+  NN+2) = Bx(ll+4*NN+5,mm+  NN+2) -                  0                  + 15*lBx*dlx*dlz*cmt(kseg,ll)*smt(kseg,mm)*rm7  &
                                                       - 3*(-cmt(kseg,ll)*lty - ll*smt(kseg,ll)*dly)*dlx*smt(kseg,mm)*rm5                      !Bx/zc / xs

     Bx(ll+5*NN+6,mm+  NN+2) = Bx(ll+5*NN+6,mm+  NN+2) -                  0                  + 15*lBx*dlx*dlz*smt(kseg,ll)*smt(kseg,mm)*rm7  &
                                                       - 3*(-smt(kseg,ll)*lty + ll*cmt(kseg,ll)*dly)*dlx*smt(kseg,mm)*rm5                      !Bx/zs / xs
     !------------------------------------------------------------------Bx on yc--------------------------------------------------------------------------

     Bx(ll+2*NN+3,mm+2*NN+3) = Bx(ll+2*NN+3,mm+2*NN+3) - 3*lBx*cmt(kseg,ll)*cmt(kseg,mm)*rm5 + 15*lBx*dly*dly*cmt(kseg,ll)*cmt(kseg,mm)*rm7  & 
                                                       - 3*( cmt(kseg,ll)*ltz + ll*smt(kseg,ll)*dlz)*dly*cmt(kseg,mm)*rm5                    &
                                                       - 3*( cmt(kseg,mm)*ltz + mm*smt(kseg,mm)*dlz)*dly*cmt(kseg,ll)*rm5                      !Bx/yc / yc

     Bx(ll+3*NN+4,mm+2*NN+3) = Bx(ll+3*NN+4,mm+2*NN+3) - 3*lBx*smt(kseg,ll)*cmt(kseg,mm)*rm5 + 15*lBx*dly*dly*smt(kseg,ll)*cmt(kseg,mm)*rm7  & 
                                                       - 3*( smt(kseg,ll)*ltz - ll*cmt(kseg,ll)*dlz)*dly*cmt(kseg,mm)*rm5                    &
                                                       - 3*( cmt(kseg,mm)*ltz + mm*smt(kseg,mm)*dlz)*dly*smt(kseg,ll)*rm5                      !Bx/ys / yc

     Bx(ll+4*NN+5,mm+2*NN+3) = Bx(ll+4*NN+5,mm+2*NN+3) -                  0                  + 15*lBx*dly*dlz*cmt(kseg,ll)*cmt(kseg,mm)*rm7  &
                                                       - 3*(-cmt(kseg,ll)*lty - ll*smt(kseg,ll)*dly)*dly*cmt(kseg,mm)*rm5                    &
                                                       - 3*( cmt(kseg,mm)*ltz + mm*smt(kseg,mm)*dlz)*dlz*cmt(kseg,ll)*rm5                    &
                                                       +   ( mm*cmt(kseg,ll)*smt(kseg,mm) - ll*smt(kseg,ll)*cmt(kseg,mm) )*rm3                 !Bx/zc / yc

     Bx(ll+5*NN+6,mm+2*NN+3) = Bx(ll+5*NN+6,mm+2*NN+3) -                  0                  + 15*lBx*dly*dlz*smt(kseg,ll)*cmt(kseg,mm)*rm7  &
                                                       - 3*(-smt(kseg,ll)*lty + ll*cmt(kseg,ll)*dly)*dly*cmt(kseg,mm)*rm5                    &
                                                       - 3*( cmt(kseg,mm)*ltz + mm*smt(kseg,mm)*dlz)*dlz*smt(kseg,ll)*rm5                    &
                                                       +   ( mm*smt(kseg,ll)*smt(kseg,mm) + ll*cmt(kseg,ll)*cmt(kseg,mm) )*rm3                 !Bx/zs / yc    
     !------------------------------------------------------------------Bx on ys-------------------------------------------------------------------------- 

     Bx(ll+3*NN+4,mm+3*NN+4) = Bx(ll+3*NN+4,mm+3*NN+4) - 3*lBx*smt(kseg,ll)*smt(kseg,mm)*rm5 + 15*lBx*dly*dly*smt(kseg,ll)*smt(kseg,mm)*rm7  & 
                                                       - 3*( smt(kseg,ll)*ltz - ll*cmt(kseg,ll)*dlz)*dly*smt(kseg,mm)*rm5                    &
                                                       - 3*( smt(kseg,mm)*ltz - mm*cmt(kseg,mm)*dlz)*dly*smt(kseg,ll)*rm5                      !Bx/ys / ys

     Bx(ll+4*NN+5,mm+3*NN+4) = Bx(ll+4*NN+5,mm+3*NN+4) -                  0                  + 15*lBx*dly*dlz*cmt(kseg,ll)*smt(kseg,mm)*rm7  &
                                                       - 3*(-cmt(kseg,ll)*lty - ll*smt(kseg,ll)*dly)*dly*smt(kseg,mm)*rm5                    &
                                                       - 3*( smt(kseg,mm)*ltz - mm*cmt(kseg,mm)*dlz)*dlz*cmt(kseg,ll)*rm5                    &
                                                       +   (-mm*cmt(kseg,ll)*cmt(kseg,mm) - ll*smt(kseg,ll)*smt(kseg,mm) )*rm3                 !Bx/zc / ys

     Bx(ll+5*NN+6,mm+3*NN+4) = Bx(ll+5*NN+6,mm+3*NN+4) -                  0                  + 15*lBx*dly*dlz*smt(kseg,ll)*smt(kseg,mm)*rm7  &
                                                       - 3*(-smt(kseg,ll)*lty + ll*cmt(kseg,ll)*dly)*dly*smt(kseg,mm)*rm5                    &
                                                       - 3*( smt(kseg,mm)*ltz - mm*cmt(kseg,mm)*dlz)*dlz*smt(kseg,ll)*rm5                    &
                                                       +   (-mm*smt(kseg,ll)*cmt(kseg,mm) + ll*cmt(kseg,ll)*smt(kseg,mm) )*rm3                 !Bx/zs / ys  
     !------------------------------------------------------------------Bx on zc-------------------------------------------------------------------------- 

     Bx(ll+4*NN+5,mm+4*NN+5) = Bx(ll+4*NN+5,mm+4*NN+5) - 3*lBx*cmt(kseg,ll)*cmt(kseg,mm)*rm5 + 15*lBx*dlz*dlz*cmt(kseg,ll)*cmt(kseg,mm)*rm7  &
                                                       - 3*(-cmt(kseg,ll)*lty - ll*smt(kseg,ll)*dly)*dlz*cmt(kseg,mm)*rm5                    &
                                                       - 3*(-cmt(kseg,mm)*lty - mm*smt(kseg,mm)*dly)*dlz*cmt(kseg,ll)*rm5                      !Bx/zc / zc

     Bx(ll+5*NN+6,mm+4*NN+5) = Bx(ll+5*NN+6,mm+4*NN+5) - 3*lBx*smt(kseg,ll)*cmt(kseg,mm)*rm5 + 15*lBx*dlz*dlz*smt(kseg,ll)*cmt(kseg,mm)*rm7  &
                                                       - 3*(-smt(kseg,ll)*lty + ll*cmt(kseg,ll)*dly)*dlz*cmt(kseg,mm)*rm5                    &
                                                       - 3*(-cmt(kseg,mm)*lty - mm*smt(kseg,mm)*dly)*dlz*smt(kseg,ll)*rm5                      !Bx/zs / zc
     !------------------------------------------------------------------Bx on zs-------------------------------------------------------------------------- 

     Bx(ll+5*NN+6,mm+5*NN+6) = Bx(ll+5*NN+6,mm+5*NN+6) - 3*lBx*smt(kseg,ll)*smt(kseg,mm)*rm5 + 15*lBx*dlz*dlz*smt(kseg,ll)*smt(kseg,mm)*rm7  &
                                                       - 3*(-smt(kseg,ll)*lty + ll*cmt(kseg,ll)*dly)*dlz*smt(kseg,mm)*rm5                    &
                                                       - 3*(-smt(kseg,mm)*lty + mm*cmt(kseg,mm)*dly)*dlz*smt(kseg,ll)*rm5                      !Bx/zs / zs

!------------------------------------------------------------------By set-------------------------------------------------------------------------------------

     !------------------------------------------------------------------By on xc-------------------------------------------------------------------------- !similar as
     By(ll     +1,mm     +1) = By(ll     +1,mm     +1) - 3*lBy*cmt(kseg,ll)*cmt(kseg,mm)*rm5 + 15*lBy*dlx*dlx*cmt(kseg,ll)*cmt(kseg,mm)*rm7  &
                                                       - 3*(-cmt(kseg,ll)*ltz - ll*smt(kseg,ll)*dlz)*dlx*cmt(kseg,mm)*rm5                    &
                                                       - 3*(-cmt(kseg,mm)*ltz - mm*smt(kseg,mm)*dlz)*dlx*cmt(kseg,ll)*rm5                      !By/xc / xc !Bx/zc/zc

     By(ll+  NN+2,mm     +1) = By(ll+  NN+2,mm     +1) - 3*lBy*smt(kseg,ll)*cmt(kseg,mm)*rm5 + 15*lBy*dlx*dlx*smt(kseg,ll)*cmt(kseg,mm)*rm7  &
                                                       - 3*(-smt(kseg,ll)*ltz + ll*cmt(kseg,ll)*dlz)*dlx*cmt(kseg,mm)*rm5                    &
                                                       - 3*(-cmt(kseg,mm)*ltz - mm*smt(kseg,mm)*dlz)*dlx*smt(kseg,ll)*rm5                      !By/xs / xc !Bx/zs/zc

     By(ll+2*NN+3,mm     +1) = By(ll+2*NN+3,mm     +1) -                  0                  + 15*lBy*dlx*dly*cmt(kseg,mm)*cmt(kseg,ll)*rm7  &
                                                       - 3*(-cmt(kseg,mm)*ltz - mm*smt(kseg,mm)*dlz)*dly*cmt(kseg,ll)*rm5                      !By/yc / xc !Bx/xc/zc

     By(ll+3*NN+4,mm     +1) = By(ll+3*NN+4,mm     +1) -                  0                  + 15*lBy*dlx*dly*cmt(kseg,mm)*smt(kseg,ll)*rm7  &
                                                       - 3*(-cmt(kseg,mm)*ltz - mm*smt(kseg,mm)*dlz)*dly*smt(kseg,ll)*rm5                      !By/ys / xc !Bx/xs/zc

     By(ll+4*NN+5,mm     +1) = By(ll+4*NN+5,mm     +1) -                  0                  + 15*lBy*dlx*dlz*cmt(kseg,mm)*cmt(kseg,ll)*rm7  &
                                                       - 3*(-cmt(kseg,mm)*ltz - mm*smt(kseg,mm)*dlz)*dlz*cmt(kseg,ll)*rm5                    &
                                                       - 3*( cmt(kseg,ll)*ltx + ll*smt(kseg,ll)*dlx)*dlx*cmt(kseg,mm)*rm5                    &
                                                       +   ( ll*cmt(kseg,mm)*smt(kseg,ll) - mm*smt(kseg,mm)*cmt(kseg,ll) )*rm3                 !By/zc / xc !Bx/yc/zc  

     By(ll+5*NN+6,mm     +1) = By(ll+5*NN+6,mm     +1) -                  0                  + 15*lBy*dlx*dlz*cmt(kseg,mm)*smt(kseg,ll)*rm7  &
                                                       - 3*(-cmt(kseg,mm)*ltz - mm*smt(kseg,mm)*dlz)*dlz*smt(kseg,ll)*rm5                    &
                                                       - 3*( smt(kseg,ll)*ltx - ll*cmt(kseg,ll)*dlx)*dlx*cmt(kseg,mm)*rm5                    &
                                                       +   (-ll*cmt(kseg,mm)*cmt(kseg,ll) - mm*smt(kseg,mm)*smt(kseg,ll) )*rm3                 !By/zs / xc !Bx/ys/zc  
     !------------------------------------------------------------------By on xs--------------------------------------------------------------------------

     By(ll+  NN+2,mm+  NN+2) = By(ll+  NN+2,mm+  NN+2) - 3*lBy*smt(kseg,ll)*smt(kseg,mm)*rm5 + 15*lBy*dlx*dlx*smt(kseg,ll)*smt(kseg,mm)*rm7  &
                                                       - 3*(-smt(kseg,ll)*ltz + ll*cmt(kseg,ll)*dlz)*dlx*smt(kseg,mm)*rm5                    &
                                                       - 3*(-smt(kseg,mm)*ltz + mm*cmt(kseg,mm)*dlz)*dlx*smt(kseg,ll)*rm5                      !By/xs / xs !Bx/zs/zs

     By(ll+2*NN+3,mm+  NN+2) = By(ll+2*NN+3,mm+  NN+2) -                  0                  + 15*lBy*dlx*dly*smt(kseg,mm)*cmt(kseg,ll)*rm7  &
                                                       - 3*(-smt(kseg,mm)*ltz + mm*cmt(kseg,mm)*dlz)*dly*cmt(kseg,ll)*rm5                      !By/yc / xs !Bx/xc/zs

     By(ll+3*NN+4,mm+  NN+2) = By(ll+3*NN+4,mm+  NN+2) -                  0                  + 15*lBy*dlx*dly*smt(kseg,mm)*smt(kseg,ll)*rm7  &
                                                       - 3*(-smt(kseg,mm)*ltz + mm*cmt(kseg,mm)*dlz)*dly*smt(kseg,ll)*rm5                      !By/ys / xs !Bx/xs/zs

     By(ll+4*NN+5,mm+  NN+2) = By(ll+4*NN+5,mm+  NN+2) -                  0                  + 15*lBy*dlx*dlz*smt(kseg,mm)*cmt(kseg,ll)*rm7  &
                                                       - 3*(-smt(kseg,mm)*ltz + mm*cmt(kseg,mm)*dlz)*dlz*cmt(kseg,ll)*rm5                    &
                                                       - 3*( cmt(kseg,ll)*ltx + ll*smt(kseg,ll)*dlx)*dlx*smt(kseg,mm)*rm5                    &
                                                       +   ( ll*smt(kseg,mm)*smt(kseg,ll) + mm*cmt(kseg,mm)*cmt(kseg,ll) )*rm3                 !By/zc / xs !Bx/yc/zs  

     By(ll+5*NN+6,mm+  NN+2) = By(ll+5*NN+6,mm+  NN+2) -                  0                  + 15*lBy*dlx*dlz*smt(kseg,mm)*smt(kseg,ll)*rm7  &
                                                       - 3*(-smt(kseg,mm)*ltz + mm*cmt(kseg,mm)*dlz)*dlz*smt(kseg,ll)*rm5                    &
                                                       - 3*( smt(kseg,ll)*ltx - ll*cmt(kseg,ll)*dlx)*dlx*smt(kseg,mm)*rm5                    &
                                                       +   (-ll*smt(kseg,mm)*cmt(kseg,ll) + mm*cmt(kseg,mm)*smt(kseg,ll) )*rm3                 !By/zs / xs !Bx/ys/zs  
     !------------------------------------------------------------------By on yc--------------------------------------------------------------------------

     By(ll+2*NN+3,mm+2*NN+3) = By(ll+2*NN+3,mm+2*NN+3) - 3*lBy*cmt(kseg,ll)*cmt(kseg,mm)*rm5 + 15*lBy*dly*dly*cmt(kseg,ll)*cmt(kseg,mm)*rm7    !By/yc / yc !Bx/xc/xc

     By(ll+3*NN+4,mm+2*NN+3) = By(ll+3*NN+4,mm+2*NN+3) - 3*lBy*smt(kseg,ll)*cmt(kseg,mm)*rm5 + 15*lBy*dly*dly*smt(kseg,ll)*cmt(kseg,mm)*rm7    !By/ys / yc !Bx/xs/xc

     By(ll+4*NN+5,mm+2*NN+3) = By(ll+4*NN+5,mm+2*NN+3) -                  0                  + 15*lBy*dlz*dly*cmt(kseg,ll)*cmt(kseg,mm)*rm7  & 
                                                       - 3*( cmt(kseg,ll)*ltx + ll*smt(kseg,ll)*dlx)*dly*cmt(kseg,mm)*rm5                      !By/zc / yc !Bx/yc/xc

     By(ll+5*NN+6,mm+2*NN+3) = By(ll+5*NN+6,mm+2*NN+3) -                  0                  + 15*lBy*dlz*dly*smt(kseg,ll)*cmt(kseg,mm)*rm7  & 
                                                       - 3*( smt(kseg,ll)*ltx - ll*cmt(kseg,ll)*dlx)*dly*cmt(kseg,mm)*rm5                      !By/zs / yc !Bx/ys/xc   
     !------------------------------------------------------------------By on ys-------------------------------------------------------------------------- 
     
     By(ll+3*NN+4,mm+3*NN+4) = By(ll+3*NN+4,mm+3*NN+4) - 3*lBy*smt(kseg,ll)*smt(kseg,mm)*rm5 + 15*lBy*dly*dly*smt(kseg,ll)*smt(kseg,mm)*rm7    !By/ys / ys !Bx/xs/xs


     By(ll+4*NN+5,mm+3*NN+4) = By(ll+4*NN+5,mm+3*NN+4) -                  0                  + 15*lBy*dlz*dly*cmt(kseg,ll)*smt(kseg,mm)*rm7  & 
                                                       - 3*( cmt(kseg,ll)*ltx + ll*smt(kseg,ll)*dlx)*dly*smt(kseg,mm)*rm5                      !By/zc / ys !Bx/yc/xs

     By(ll+5*NN+6,mm+3*NN+4) = By(ll+5*NN+6,mm+3*NN+4) -                  0                  + 15*lBy*dlz*dly*smt(kseg,ll)*smt(kseg,mm)*rm7  & 
                                                       - 3*( smt(kseg,ll)*ltx - ll*cmt(kseg,ll)*dlx)*dly*smt(kseg,mm)*rm5                      !By/zs / ys !Bx/ys/xs  
     !------------------------------------------------------------------By on zc-------------------------------------------------------------------------- 

     By(ll+4*NN+5,mm+4*NN+5) = By(ll+4*NN+5,mm+4*NN+5) - 3*lBy*cmt(kseg,ll)*cmt(kseg,mm)*rm5 + 15*lBy*dlz*dlz*cmt(kseg,ll)*cmt(kseg,mm)*rm7  & 
                                                       - 3*( cmt(kseg,ll)*ltx + ll*smt(kseg,ll)*dlx)*dlz*cmt(kseg,mm)*rm5                    &
                                                       - 3*( cmt(kseg,mm)*ltx + mm*smt(kseg,mm)*dlx)*dlz*cmt(kseg,ll)*rm5                      !By/zc / zc !Bx/yc/yc

     By(ll+5*NN+6,mm+4*NN+5) = By(ll+5*NN+6,mm+4*NN+5) - 3*lBy*smt(kseg,ll)*cmt(kseg,mm)*rm5 + 15*lBy*dlz*dlz*smt(kseg,ll)*cmt(kseg,mm)*rm7  & 
                                                       - 3*( smt(kseg,ll)*ltx - ll*cmt(kseg,ll)*dlx)*dlz*cmt(kseg,mm)*rm5                    &
                                                       - 3*( cmt(kseg,mm)*ltx + mm*smt(kseg,mm)*dlx)*dlz*smt(kseg,ll)*rm5                      !By/zs / zc !Bx/ys/yc
     !------------------------------------------------------------------By on zs-------------------------------------------------------------------------- 

     By(ll+5*NN+6,mm+5*NN+6) = By(ll+5*NN+6,mm+5*NN+6) - 3*lBy*smt(kseg,ll)*smt(kseg,mm)*rm5 + 15*lBy*dlz*dlz*smt(kseg,ll)*smt(kseg,mm)*rm7  & 
                                                       - 3*( smt(kseg,ll)*ltx - ll*cmt(kseg,ll)*dlx)*dlz*smt(kseg,mm)*rm5                    &
                                                       - 3*( smt(kseg,mm)*ltx - mm*cmt(kseg,mm)*dlx)*dlz*smt(kseg,ll)*rm5                      !By/zs / zs !Bx/ys/ys


!------------------------------------------------------------------Bz set-------------------------------------------------------------------------------------
     !------------------------------------------------------------------Bz on xc--------------------------------------------------------------------------

     Bz(ll     +1,mm     +1) = Bz(ll     +1,mm     +1) - 3*lBz*cmt(kseg,ll)*cmt(kseg,mm)*rm5 + 15*lBz*dlx*dlx*cmt(kseg,ll)*cmt(kseg,mm)*rm7  & 
                                                       - 3*( cmt(kseg,ll)*lty + ll*smt(kseg,ll)*dly)*dlx*cmt(kseg,mm)*rm5                    &
                                                       - 3*( cmt(kseg,mm)*lty + mm*smt(kseg,mm)*dly)*dlx*cmt(kseg,ll)*rm5                      !Bz/xc / xc !Bx/yc/yc

     Bz(ll+  NN+2,mm     +1) = Bz(ll+  NN+2,mm     +1) - 3*lBz*smt(kseg,ll)*cmt(kseg,mm)*rm5 + 15*lBz*dlx*dlx*smt(kseg,ll)*cmt(kseg,mm)*rm7  & 
                                                       - 3*( smt(kseg,ll)*lty - ll*cmt(kseg,ll)*dly)*dlx*cmt(kseg,mm)*rm5                    &
                                                       - 3*( cmt(kseg,mm)*lty + mm*smt(kseg,mm)*dly)*dlx*smt(kseg,ll)*rm5                      !Bz/xs / xc !Bx/ys/yc

     Bz(ll+2*NN+3,mm     +1) = Bz(ll+2*NN+3,mm     +1) -                  0                  + 15*lBz*dly*dlx*cmt(kseg,ll)*cmt(kseg,mm)*rm7  &
                                                       - 3*(-cmt(kseg,ll)*ltx - ll*smt(kseg,ll)*dlx)*dlx*cmt(kseg,mm)*rm5                    &
                                                       - 3*( cmt(kseg,mm)*lty + mm*smt(kseg,mm)*dly)*dly*cmt(kseg,ll)*rm5                    &
                                                       +   ( mm*cmt(kseg,ll)*smt(kseg,mm) - ll*smt(kseg,ll)*cmt(kseg,mm) )*rm3                 !Bz/yc / xc !Bx/zc/yc

     Bz(ll+3*NN+4,mm     +1) = Bz(ll+3*NN+4,mm     +1) -                  0                  + 15*lBz*dly*dlx*smt(kseg,ll)*cmt(kseg,mm)*rm7  &
                                                       - 3*(-smt(kseg,ll)*ltx + ll*cmt(kseg,ll)*dlx)*dlx*cmt(kseg,mm)*rm5                    &
                                                       - 3*( cmt(kseg,mm)*lty + mm*smt(kseg,mm)*dly)*dly*smt(kseg,ll)*rm5                    &
                                                       +   ( mm*smt(kseg,ll)*smt(kseg,mm) + ll*cmt(kseg,ll)*cmt(kseg,mm) )*rm3                 !Bz/ys / xc !Bx/zs/yc

     Bz(ll+4*NN+5,mm     +1) = Bz(ll+4*NN+5,mm     +1) -                  0                  + 15*lBz*dlz*dlx*cmt(kseg,mm)*cmt(kseg,ll)*rm7  & 
                                                       - 3*( cmt(kseg,mm)*lty + mm*smt(kseg,mm)*dly)*dlz*cmt(kseg,ll)*rm5                      !Bz/zc / xc !Bx/xc/yc

     Bz(ll+5*NN+6,mm     +1) = Bz(ll+5*NN+6,mm     +1) -                  0                  + 15*lBz*dlz*dlx*cmt(kseg,mm)*smt(kseg,ll)*rm7  & 
                                                       - 3*( cmt(kseg,mm)*lty + mm*smt(kseg,mm)*dly)*dlz*smt(kseg,ll)*rm5                      !Bz/zs / xc !Bx/xs/yc
     !------------------------------------------------------------------Bz on xs--------------------------------------------------------------------------

     Bz(ll+  NN+2,mm+  NN+2) = Bz(ll+  NN+2,mm+  NN+2) - 3*lBz*smt(kseg,ll)*smt(kseg,mm)*rm5 + 15*lBz*dlx*dlx*smt(kseg,ll)*smt(kseg,mm)*rm7  & 
                                                       - 3*( smt(kseg,ll)*lty - ll*cmt(kseg,ll)*dly)*dlx*smt(kseg,mm)*rm5                    &
                                                       - 3*( smt(kseg,mm)*lty - mm*cmt(kseg,mm)*dly)*dlx*smt(kseg,ll)*rm5                      !Bz/xs / xs !Bx/ys/ys

     Bz(ll+2*NN+3,mm+  NN+2) = Bz(ll+2*NN+3,mm+  NN+2) -                  0                  + 15*lBz*dly*dlx*cmt(kseg,ll)*smt(kseg,mm)*rm7  &
                                                       - 3*(-cmt(kseg,ll)*ltx - ll*smt(kseg,ll)*dlx)*dlx*smt(kseg,mm)*rm5                    &
                                                       - 3*( smt(kseg,mm)*lty - mm*cmt(kseg,mm)*dly)*dly*cmt(kseg,ll)*rm5                    &
                                                       +   (-mm*cmt(kseg,ll)*cmt(kseg,mm) - ll*smt(kseg,ll)*smt(kseg,mm) )*rm3                 !Bz/yc / xs !Bx/zc/ys

     Bz(ll+3*NN+4,mm+  NN+2) = Bz(ll+3*NN+4,mm+  NN+2) -                  0                  + 15*lBz*dly*dlx*smt(kseg,ll)*smt(kseg,mm)*rm7  &
                                                       - 3*(-smt(kseg,ll)*ltx + ll*cmt(kseg,ll)*dlx)*dlx*smt(kseg,mm)*rm5                    &
                                                       - 3*( smt(kseg,mm)*lty - mm*cmt(kseg,mm)*dly)*dly*smt(kseg,ll)*rm5                    &
                                                       +   (-mm*smt(kseg,ll)*cmt(kseg,mm) + ll*cmt(kseg,ll)*smt(kseg,mm) )*rm3                 !Bz/ys / xs !Bx/zs/ys

     Bz(ll+4*NN+5,mm+  NN+2) = Bz(ll+4*NN+5,mm+  NN+2) -                  0                  + 15*lBz*dlz*dlx*smt(kseg,mm)*cmt(kseg,ll)*rm7  & 
                                                       - 3*( smt(kseg,mm)*lty - mm*cmt(kseg,mm)*dly)*dlz*cmt(kseg,ll)*rm5                      !Bz/zc / xs !Bx/xc/ys

     Bz(ll+5*NN+6,mm+  NN+2) = Bz(ll+5*NN+6,mm+  NN+2) -                  0                  + 15*lBz*dlz*dlx*smt(kseg,mm)*smt(kseg,ll)*rm7  & 
                                                       - 3*( smt(kseg,mm)*lty - mm*cmt(kseg,mm)*dly)*dlz*smt(kseg,ll)*rm5                      !Bz/zs / xs !Bx/xs/ys
     !------------------------------------------------------------------Bz on yc--------------------------------------------------------------------------

     Bz(ll+2*NN+3,mm+2*NN+3) = Bz(ll+2*NN+3,mm+2*NN+3) - 3*lBz*cmt(kseg,ll)*cmt(kseg,mm)*rm5 + 15*lBz*dly*dly*cmt(kseg,ll)*cmt(kseg,mm)*rm7  &
                                                       - 3*(-cmt(kseg,ll)*ltx - ll*smt(kseg,ll)*dlx)*dly*cmt(kseg,mm)*rm5                    &
                                                       - 3*(-cmt(kseg,mm)*ltx - mm*smt(kseg,mm)*dlx)*dly*cmt(kseg,ll)*rm5                      !Bz/yc / yc !Bx/zc/zc

     Bz(ll+3*NN+4,mm+2*NN+3) = Bz(ll+3*NN+4,mm+2*NN+3) - 3*lBz*smt(kseg,ll)*cmt(kseg,mm)*rm5 + 15*lBz*dly*dly*smt(kseg,ll)*cmt(kseg,mm)*rm7  &
                                                       - 3*(-smt(kseg,ll)*ltx + ll*cmt(kseg,ll)*dlx)*dly*cmt(kseg,mm)*rm5                    &
                                                       - 3*(-cmt(kseg,mm)*ltx - mm*smt(kseg,mm)*dlx)*dly*smt(kseg,ll)*rm5                      !Bz/ys / yc !Bx/zs/zc

     Bz(ll+4*NN+5,mm+2*NN+3) = Bz(ll+4*NN+5,mm+2*NN+3) -                  0                  + 15*lBz*dlz*dly*cmt(kseg,mm)*cmt(kseg,ll)*rm7  &
                                                       - 3*(-cmt(kseg,mm)*ltx - mm*smt(kseg,mm)*dlx)*dlz*cmt(kseg,ll)*rm5                      !Bz/zc / yc !Bx/xc/zc


     Bz(ll+5*NN+6,mm+2*NN+3) = Bz(ll+5*NN+6,mm+2*NN+3) -                  0                  + 15*lBz*dlz*dly*cmt(kseg,mm)*smt(kseg,ll)*rm7  &
                                                       - 3*(-cmt(kseg,mm)*ltx - mm*smt(kseg,mm)*dlx)*dlz*smt(kseg,ll)*rm5                      !Bz/zs / yc !Bx/xs/zc   
     !------------------------------------------------------------------Bz on ys-------------------------------------------------------------------------- 

     Bz(ll+3*NN+4,mm+3*NN+4) = Bz(ll+3*NN+4,mm+3*NN+4) - 3*lBz*smt(kseg,ll)*smt(kseg,mm)*rm5 + 15*lBz*dly*dly*smt(kseg,ll)*smt(kseg,mm)*rm7  &
                                                       - 3*(-smt(kseg,ll)*ltx + ll*cmt(kseg,ll)*dlx)*dly*smt(kseg,mm)*rm5                    &
                                                       - 3*(-smt(kseg,mm)*ltx + mm*cmt(kseg,mm)*dlx)*dly*smt(kseg,ll)*rm5                      !Bz/ys / ys !Bx/zs/zs

     Bz(ll+4*NN+5,mm+3*NN+4) = Bz(ll+4*NN+5,mm+3*NN+4) -                  0                  + 15*lBz*dlz*dly*smt(kseg,mm)*cmt(kseg,ll)*rm7  &
                                                       - 3*(-smt(kseg,mm)*ltx + mm*cmt(kseg,mm)*dlx)*dlz*cmt(kseg,ll)*rm5                      !Bz/zc / ys !Bx/xc/zs


     Bz(ll+5*NN+6,mm+3*NN+4) = Bz(ll+5*NN+6,mm+3*NN+4) -                  0                  + 15*lBz*dlz*dly*smt(kseg,mm)*smt(kseg,ll)*rm7  &
                                                       - 3*(-smt(kseg,mm)*ltx + mm*cmt(kseg,mm)*dlx)*dlz*smt(kseg,ll)*rm5                      !Bz/zs / ys !Bx/xs/zs    
     !------------------------------------------------------------------Bz on zc--------------------------------------------------------------------------
 
     Bz(ll+4*NN+5,mm+4*NN+5) = Bz(ll+4*NN+5,mm+4*NN+5) - 3*lBz*cmt(kseg,ll)*cmt(kseg,mm)*rm5 + 15*lBz*dlz*dlz*cmt(kseg,ll)*cmt(kseg,mm)*rm7    !Bz/zc / zc !Bx/xc/xc

     Bz(ll+5*NN+6,mm+4*NN+5) = Bz(ll+5*NN+6,mm+4*NN+5) - 3*lBz*smt(kseg,ll)*cmt(kseg,mm)*rm5 + 15*lBz*dlz*dlz*smt(kseg,ll)*cmt(kseg,mm)*rm7    !Bz/zs / zc !Bx/xs/xc
     !------------------------------------------------------------------Bz on zs-------------------------------------------------------------------------- 

     Bz(ll+5*NN+6,mm+5*NN+6) = Bz(ll+5*NN+6,mm+5*NN+6) - 3*lBz*smt(kseg,ll)*smt(kseg,mm)*rm5 + 15*lBz*dlz*dlz*smt(kseg,ll)*smt(kseg,mm)*rm7    !Bz/zs / zs !Bx/xs/xs

    enddo !enddo ll
   enddo !enddo mm
   
  enddo    ! enddo kseg

  do mm = 2, Cdof
   do ll = 1, mm-1

    Bx(ll,mm) = Bx(mm,ll)
    By(ll,mm) = By(mm,ll)
    Bz(ll,mm) = Bz(mm,ll)

   enddo
  enddo

  return

end subroutine bfield2
