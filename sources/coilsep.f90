!title (coilsep) ! Calculate coil separation objective functon and its derivatives. (tkruger)

!latex \briefly{The constraint on coil to coil separation solves for coils that have
!latex         reasonable finite builds and can be assembled without interferance.

!latex \calledby{\link{solvers}}

!latex  \section{General}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! ccsep is total penalty
! chi = chi + weight_ccsep*ccsep
! t1C is total derivative of penalty
! LM implemented
subroutine coilsep(ideriv)
  use globals, only: dp, zero, half, pi2, machprec, ncpu, myid, ounit, &
       coil, DoF, Ncoils, Nfixgeo, Ndof, ccsep, t1C, t2C, Nfp, &
       iccsep, mccsep, LM_fvec, LM_fjac, weight_ccsep, &
       MPI_COMM_FOCUS

  implicit none
  include "mpif.h"
  INTEGER, INTENT(in) :: ideriv

  INTEGER             :: astat, ierr, icoil, idof, ND, ivec, numCoil
  REAL                :: d1C(1:Ndof), ccsep0

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  ccsep = zero
  ivec = 1
  numCoil = zero
  ccsep0 = zero
 
  ! Do calculation for free and fixed coils 
  if( ideriv >= 0 ) then
     if ( Ncoils .eq. 1 ) return
     do icoil = 1, Ncoils
        call CoilSepDeriv0(icoil, ccsep0)
        if ( coil(icoil)%symm == 0 ) numCoil = numCoil + 1
        if ( coil(icoil)%symm == 1 ) numCoil = numCoil + Nfp
        if ( coil(icoil)%symm == 2 ) numCoil = numCoil + 2*Nfp
        ccsep = ccsep + ccsep0
        if ( coil(icoil)%Lc /= 0 ) then
           if (mccsep > 0) then ! L-M format of targets
              LM_fvec(iccsep+ivec) = weight_ccsep * ccsep0
              ivec = ivec + 1
           endif
        endif
     enddo
     if (mccsep > 0) then ! L-M format of targets
        FATAL( coilsep, ivec == mccsep, Errors in counting ivec for L-M )
     endif
     ccsep = 2.0*ccsep / ( numCoil * (numCoil - 1) + machprec) ! Triangle number 
  endif

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  ! Free coil derivatives summed over fixed and free 
  if ( ideriv >= 1 ) then
     t1C = zero ; d1C = zero 
     if ( Ncoils .eq. 1 ) return
     idof = 0 ; ivec = 1
     do icoil = 1, Ncoils
        ND = DoF(icoil)%ND
        if ( coil(icoil)%Ic /= 0 ) then !if current is free
           idof = idof +1
        endif
        if ( coil(icoil)%Lc /= 0 ) then !if geometry is free
           call CoilSepDeriv1( icoil, t1C(idof+1:idof+ND), ND )
           if (mccsep > 0) then ! L-M format of targets
              LM_fjac(iccsep+ivec, idof+1:idof+ND) = weight_ccsep * t1C(idof+1:idof+ND)
              ivec = ivec + 1
           endif
           idof = idof + ND
        endif
     enddo
     FATAL( coilsep , idof .ne. Ndof, counting error in packing )
     if (mccsep > 0) then ! L-M format of targets
        FATAL( coilsep, ivec == mccsep, Errors in counting ivec for L-M )
     endif
     t1C = 2.0*t1C / ( numCoil * (numCoil - 1) + machprec )
  endif
  
  return
end subroutine coilsep

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! calculate coil length

subroutine CoilSepDeriv0(icoil, coilsep0)

  use globals, only: dp, zero, pi2, coil, myid, ncpu, ounit, Ncoils, Nfp, machprec, r_delta, ccsep_alpha, penfun_ccsep, ccsep_beta, ccsep_skip, MPI_COMM_FOCUS 
  implicit none
  include "mpif.h"

  INTEGER, intent(in)  :: icoil
  REAL   , intent(out) :: coilsep0
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  INTEGER              :: kseg0, kseg1, astat, ierr, i, per0, per1, ss0, ss1, j0, j00, j1, l0, l00, l1, NS, NS1
  REAL                 :: rdifff, hypc, coilsepHold, xc0, yc0, zc0, xc1, yc1, zc1, coilseph
  REAL, allocatable    :: x0(:), y0(:), z0(:), x1(:), y1(:), z1(:), absrp0(:), absrp1(:)

  NS = coil(icoil)%NS

  SALLOCATE(x0, (0:NS-1), zero)
  SALLOCATE(y0, (0:NS-1), zero)
  SALLOCATE(z0, (0:NS-1), zero)
  SALLOCATE(absrp0, (0:NS-1), zero)

  FATAL( CoilSepDeriv0, icoil .lt. 1 .or. icoil .gt. Ncoils, icoil not in right range )
  FATAL( CoilSepDeriv0, penfun_ccsep .ne. 1 .and. penfun_ccsep .ne. 2, invalid choice of penfun_ccsep, pick 1 or 2 )
  FATAL( CurvDeriv0, r_delta .le. 0.0 , r_delta needs to be positive )
  FATAL( CurvDeriv0, ccsep_alpha < 0.0 , ccsep_alpha cannot be negative )
  FATAL( CurvDeriv0, ccsep_beta < 2.0 , ccsep_beta needs to be >= 2 )

  coilsep0 = zero
  coilseph = zero
  coilsepHold = zero

  if ( coil(icoil)%symm == 0 ) then 
     per0 = 1
     ss0 = 0
  elseif ( coil(icoil)%symm == 1 ) then
     per0 = Nfp
     ss0 = 0
  elseif ( coil(icoil)%symm == 2 ) then
     per0 = Nfp
     ss0 = 1
  else
     FATAL( CoilSepDeriv0, .true. , Errors in coil symmetry )
  endif
  do j0 = 1, per0
     do l0 = 0, ss0
        x0(0:NS-1) = (coil(icoil)%xx(0:NS-1))*cos(pi2*(j0-1)/Nfp) - (coil(icoil)%yy(0:NS-1))*sin(pi2*(j0-1)/Nfp)
        y0(0:NS-1) = ((-1.0)**l0)*((coil(icoil)%yy(0:NS-1))*cos(pi2*(j0-1)/Nfp) + (coil(icoil)%xx(0:NS-1))*sin(pi2*(j0-1)/Nfp))
        z0(0:NS-1) = (coil(icoil)%zz(0:NS-1))*((-1.0)**l0)
        if ( ccsep_skip .eq. 1 ) then
           absrp0(0:NS-1) = sqrt( (coil(icoil)%xt(0:NS-1))**2 + (coil(icoil)%yt(0:NS-1))**2 + (coil(icoil)%zt(0:NS-1))**2 )
           xc0 = sum( x0(0:NS-1)*absrp0(0:NS-1) ) / sum( absrp0(0:NS-1) )
           yc0 = sum( y0(0:NS-1)*absrp0(0:NS-1) ) / sum( absrp0(0:NS-1) )
           zc0 = sum( z0(0:NS-1)*absrp0(0:NS-1) ) / sum( absrp0(0:NS-1) )
        endif
        if ( coil(icoil)%symm /= 0 ) then
           SALLOCATE(x1, (0:NS-1), zero)
           SALLOCATE(y1, (0:NS-1), zero)
           SALLOCATE(z1, (0:NS-1), zero)
           do j00 = j0, per0
              do l00 = l0, ss0
                 if ( j0 .eq. j00 .and. l0 .eq. l00 ) cycle
                 x1(0:NS-1) = (coil(icoil)%xx(0:NS-1))*cos(pi2*(j00-1)/Nfp) - (coil(icoil)%yy(0:NS-1))*sin(pi2*(j00-1)/Nfp)
                 y1(0:NS-1) = ((-1.0)**l00)*((coil(icoil)%yy(0:NS-1))*cos(pi2*(j00-1)/Nfp) + (coil(icoil)%xx(0:NS-1))*sin(pi2*(j00-1)/Nfp))
                 z1(0:NS-1) = (coil(icoil)%zz(0:NS-1))*((-1.0)**l00)
                 if ( ccsep_skip .eq. 1 ) then
                    xc1 = sum( x1(0:NS-1)*absrp0(0:NS-1) ) / sum( absrp0(0:NS-1) )
                    yc1 = sum( y1(0:NS-1)*absrp0(0:NS-1) ) / sum( absrp0(0:NS-1) )
                    zc1 = sum( z1(0:NS-1)*absrp0(0:NS-1) ) / sum( absrp0(0:NS-1) )
                    if ( sqrt((xc0-xc1)**2+(yc0-yc1)**2+(zc0-zc1)**2) .ge. (coil(icoil)%Lo+coil(icoil)%Lo)*.25 ) cycle
                 endif
                 do kseg0 = 0, NS-1
                    if( myid.ne.modulo(kseg0,ncpu) ) cycle ! parallelization loop;
                    do kseg1 = 0, NS-1
                       rdifff = sqrt( ( x0(kseg0) - x1(kseg1) )**2 + ( y0(kseg0) - y1(kseg1) )**2 + ( z0(kseg0) - z1(kseg1) )**2 )
                       if ( rdifff < r_delta ) then
                          if ( penfun_ccsep .eq. 1 ) then
                             hypc = 0.5*exp(ccsep_alpha*( r_delta - rdifff )) + 0.5*exp(-1.0*ccsep_alpha*( r_delta - rdifff ))
                             coilseph = coilseph + (hypc - 1.0)**2
                          else 
                             coilseph = coilseph + (ccsep_alpha*( r_delta - rdifff ))**ccsep_beta
                          endif
                       endif
                    enddo
                 enddo
                 call MPI_BARRIER( MPI_COMM_FOCUS, ierr )
                 call MPI_ALLREDUCE( MPI_IN_PLACE, coilseph, 1 , MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FOCUS, ierr )
                 coilsepHold = coilsepHold + coilseph
                 coilseph = zero
              enddo
           enddo
           coilsep0 = coilsep0 + pi2*pi2*coilsepHold/((NS)*(NS))
           coilsepHold = 0.0
           DALLOCATE(x1)
           DALLOCATE(y1)
           DALLOCATE(z1) 
        endif
        do i = icoil+1, Ncoils 
           if (coil(icoil)%Lc == 0 .and. coil(i)%Lc == 0) cycle
           NS1 = coil(i)%NS
           SALLOCATE(x1, (0:NS1-1), zero)
           SALLOCATE(y1, (0:NS1-1), zero)
           SALLOCATE(z1, (0:NS1-1), zero)
           SALLOCATE(absrp1, (0:NS1-1), zero)
           if ( coil(i)%symm == 0 ) then
              per1 = 1
              ss1 = 0
           elseif ( coil(i)%symm == 1 ) then
              per1 = Nfp
              ss1 = 0
           elseif ( coil(i)%symm == 2 ) then
              per1 = Nfp
              ss1 = 1
           else
              FATAL( CoilSepDeriv0, .true. , Errors in coil symmetry )
           endif
           do j1 = 1, per1
              do l1 = 0, ss1
                 x1(0:NS1-1) = (coil(i)%xx(0:NS1-1))*cos(pi2*(j1-1)/Nfp) - (coil(i)%yy(0:NS1-1))*sin(pi2*(j1-1)/Nfp)
                 y1(0:NS1-1) = ((-1.0)**l1)*((coil(i)%yy(0:NS1-1))*cos(pi2*(j1-1)/Nfp) + (coil(i)%xx(0:NS1-1))*sin(pi2*(j1-1)/Nfp))
                 z1(0:NS1-1) = (coil(i)%zz(0:NS1-1))*((-1.0)**l1)
                 if ( ccsep_skip .eq. 1 ) then
                    absrp1(0:NS-1) = sqrt( (coil(i)%xt(0:NS-1))**2 + (coil(i)%yt(0:NS-1))**2 + (coil(i)%zt(0:NS-1))**2 )
                    xc1 = sum( x1(0:NS-1)*absrp1(0:NS-1) ) / sum( absrp1(0:NS-1) )
                    yc1 = sum( y1(0:NS-1)*absrp1(0:NS-1) ) / sum( absrp1(0:NS-1) )
                    zc1 = sum( z1(0:NS-1)*absrp1(0:NS-1) ) / sum( absrp1(0:NS-1) )
                    if ( sqrt((xc0-xc1)**2+(yc0-yc1)**2+(zc0-zc1)**2 ) .ge. (coil(i)%Lo+coil(icoil)%Lo)*.25 ) cycle
                 endif
                 do kseg0 = 0, NS-1
                    if( myid.ne.modulo(kseg0,ncpu) ) cycle ! parallelization loop;
                    do kseg1 = 0, NS1-1
                       rdifff = sqrt( ( x0(kseg0) - x1(kseg1) )**2 + ( y0(kseg0) - y1(kseg1) )**2 + ( z0(kseg0) - z1(kseg1) )**2 )
                       if ( rdifff < r_delta ) then
                          if ( penfun_ccsep .eq. 1 ) then
                             hypc = 0.5*exp(ccsep_alpha*( r_delta - rdifff )) + 0.5*exp(-1.0*ccsep_alpha*( r_delta - rdifff ))
                             coilseph = coilseph + (hypc - 1.0)**2
                          else 
                             coilseph = coilseph + (ccsep_alpha*( r_delta - rdifff ))**ccsep_beta
                          endif
                       endif
                    enddo
                 enddo
                 call MPI_BARRIER( MPI_COMM_FOCUS, ierr )
                 call MPI_ALLREDUCE( MPI_IN_PLACE, coilseph, 1 , MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FOCUS, ierr )
                 coilsepHold = coilsepHold + coilseph
                 coilseph = zero
              enddo
           enddo
           coilsep0 = coilsep0 + pi2*pi2*coilsepHold/(NS*NS1)
           coilsepHold = 0.0
           DALLOCATE(x1)
           DALLOCATE(y1)
           DALLOCATE(z1)
           DALLOCATE(absrp1)
        enddo
     enddo
  enddo
  DALLOCATE(x0)
  DALLOCATE(y0)
  DALLOCATE(z0)
  DALLOCATE(absrp0)

  return

end subroutine CoilSepDeriv0

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!calculate coil length in derivs(0) and first derivatives in derivs(1:Cdof)

subroutine CoilSepDeriv1(icoil, derivs, ND)

  use globals, only: dp, zero, pi2, coil, DoF, myid, ncpu, ounit, machprec, Ncoils, Nfp, r_delta, ccsep_alpha, penfun_ccsep, ccsep_beta, ccsep_skip, MPI_COMM_FOCUS
  implicit none
  include "mpif.h"

  INTEGER, intent(in)  :: icoil, ND
  REAL   , intent(out) :: derivs(1:1, 1:ND)
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  INTEGER              :: kseg0, kseg1, astat, ierr, per0, per1, ss0, ss1, j0, j00, j1, l0, l00, l1, i, NS, NS1
  REAL                 :: dl3, hypc, hyps, rdifff, int_hold, xc0, yc0, zc0, xc1, yc1, zc1
  REAL, allocatable    :: x0(:), y0(:), z0(:), x1(:), y1(:), z1(:), derivs_hold(:,:), absrp0(:), absrp1(:), derivsh(:)

  FATAL( CoilSepDeriv1, icoil .lt. 1 .or. icoil .gt. Ncoils, icoil not in right range )
  
  derivs = zero

  SALLOCATE(derivs_hold, (1:1,1:ND), zero)
  SALLOCATE(derivsh, (1:ND), zero)
  
  NS = coil(icoil)%NS

  SALLOCATE(x0, (0:NS-1), zero)
  SALLOCATE(y0, (0:NS-1), zero)
  SALLOCATE(z0, (0:NS-1), zero)
  SALLOCATE(absrp0, (0:NS-1), zero)

  if ( coil(icoil)%symm == 0 ) then
     per0 = 1
     ss0 = 0
  elseif ( coil(icoil)%symm == 1 ) then
     per0 = Nfp
     ss0 = 0
  elseif ( coil(icoil)%symm == 2 ) then
     per0 = Nfp
     ss0 = 1
  else
     FATAL( CoilSepDeriv1, .true. , Errors in coil symmetry )
  endif 
  do j0 = 1, per0
     do l0 = 0, ss0
        x0(0:NS-1) = (coil(icoil)%xx(0:NS-1))*cos(pi2*(j0-1)/Nfp) - (coil(icoil)%yy(0:NS-1))*sin(pi2*(j0-1)/Nfp)
        y0(0:NS-1) = ((-1.0)**l0)*((coil(icoil)%yy(0:NS-1))*cos(pi2*(j0-1)/Nfp) + (coil(icoil)%xx(0:NS-1))*sin(pi2*(j0-1)/Nfp))
        z0(0:NS-1) = (coil(icoil)%zz(0:NS-1))*((-1.0)**l0)
        if ( ccsep_skip .eq. 1 ) then
           absrp0(0:NS-1) = sqrt( (coil(icoil)%xt(0:NS-1))**2 + (coil(icoil)%yt(0:NS-1))**2 + (coil(icoil)%zt(0:NS-1))**2 )
           xc0 = sum( x0(0:NS-1)*absrp0(0:NS-1) ) / sum( absrp0(0:NS-1) )
           yc0 = sum( y0(0:NS-1)*absrp0(0:NS-1) ) / sum( absrp0(0:NS-1) )
           zc0 = sum( z0(0:NS-1)*absrp0(0:NS-1) ) / sum( absrp0(0:NS-1) )
        endif
        if ( coil(icoil)%symm /= 0 ) then
           SALLOCATE(x1, (0:NS-1), zero)
           SALLOCATE(y1, (0:NS-1), zero)
           SALLOCATE(z1, (0:NS-1), zero)
           do j00 = j0, per0
              do l00 = l0, ss0
                 if ( j0 .eq. j00 .and. l0 .eq. l00 ) cycle
                 x1(0:NS-1) = (coil(icoil)%xx(0:NS-1))*cos(pi2*(j00-1)/Nfp) - (coil(icoil)%yy(0:NS-1))*sin(pi2*(j00-1)/Nfp)
                 y1(0:NS-1) = ((-1.0)**l00)*((coil(icoil)%yy(0:NS-1))*cos(pi2*(j00-1)/Nfp) + (coil(icoil)%xx(0:NS-1))*sin(pi2*(j00-1)/Nfp))
                 z1(0:NS-1) = (coil(icoil)%zz(0:NS-1))*((-1.0)**l00)
                 if ( ccsep_skip .eq. 1 ) then
                    xc1 = sum( x1(0:NS-1)*absrp0(0:NS-1) ) / sum( absrp0(0:NS-1) )
                    yc1 = sum( y1(0:NS-1)*absrp0(0:NS-1) ) / sum( absrp0(0:NS-1) )
                    zc1 = sum( z1(0:NS-1)*absrp0(0:NS-1) ) / sum( absrp0(0:NS-1) )
                    if ( sqrt((xc0-xc1)**2+(yc0-yc1)**2+(zc0-zc1)**2) .ge. (coil(icoil)%Lo+coil(icoil)%Lo)*.25 ) cycle
                 endif
                 do kseg0 = 0, NS-1
                    if( myid.ne.modulo(kseg0,ncpu) ) cycle ! parallelization loop;
                    do kseg1 = 0, NS-1
                       rdifff = sqrt( ( x0(kseg0) - x1(kseg1) )**2 + ( y0(kseg0) - y1(kseg1) )**2 + ( z0(kseg0) - z1(kseg1) )**2 )
                       if ( rdifff < r_delta ) then
                          if ( penfun_ccsep .eq. 1 ) then
                             hypc = 0.5*exp(ccsep_alpha*( r_delta - rdifff )) + 0.5*exp(-1.0*ccsep_alpha*( r_delta - rdifff ))
                             hyps = 0.5*exp(ccsep_alpha*( r_delta - rdifff )) - 0.5*exp(-1.0*ccsep_alpha*( r_delta - rdifff ))
                             int_hold = -2.0*ccsep_alpha*(hypc-1.0)*hyps/rdifff
                          else
                             int_hold = -1.0*ccsep_beta*ccsep_alpha*(ccsep_alpha*(r_delta-rdifff))**(ccsep_beta-1.0)/rdifff
                          endif
                          derivsh(:) = derivsh(:) + int_hold*( &
                                  (x0(kseg0)-x1(kseg1))*(DoF(icoil)%xof(kseg0,:)*cos(pi2*(j0-1)/Nfp)-Dof(icoil)%yof(kseg0,:)*sin(pi2*(j0-1)/Nfp)) + &
                                  (y0(kseg0)-y1(kseg1))*((-1.0)**l0)*(DoF(icoil)%yof(kseg0,:)*cos(pi2*(j0-1)/Nfp)+Dof(icoil)%xof(kseg0,:)*sin(pi2*(j0-1)/Nfp)) + &
                                  (z0(kseg0)-z1(kseg1))*Dof(icoil)%zof(kseg0,:)*((-1.0)**l0) - &
                                  (x0(kseg0)-x1(kseg1))*(DoF(icoil)%xof(kseg1,:)*cos(pi2*(j00-1)/Nfp)-Dof(icoil)%yof(kseg1,:)*sin(pi2*(j00-1)/Nfp)) - &
                                  (y0(kseg0)-y1(kseg1))*((-1.0)**l00)*(DoF(icoil)%yof(kseg1,:)*cos(pi2*(j00-1)/Nfp)+Dof(icoil)%xof(kseg1,:)*sin(pi2*(j00-1)/Nfp)) - &
                                  (z0(kseg0)-z1(kseg1))*Dof(icoil)%zof(kseg1,:)*((-1.0)**l00) )
                       endif
                    enddo
                 enddo
                 call MPI_BARRIER( MPI_COMM_FOCUS, ierr )
                 call MPI_ALLREDUCE( MPI_IN_PLACE, derivsh, ND , MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FOCUS, ierr )
                 derivs_hold(1,:) = derivs_hold(1,:) + derivsh(:)
                 derivsh(:) = zero
              enddo
           enddo
           derivs = derivs + pi2*pi2*derivs_hold/(NS*NS)
           derivs_hold = zero
           DALLOCATE(x1)
           DALLOCATE(y1)
           DALLOCATE(z1)
        endif
        do i = 1, Ncoils
           if ( i == icoil ) cycle
           if (coil(icoil)%Lc == 0 .and. coil(i)%Lc == 0) cycle
           if ( coil(i)%symm == 0 ) then
              per1 = 1
              ss1 = 0
           elseif ( coil(i)%symm == 1 ) then
              per1 = Nfp
              ss1 = 0
           elseif ( coil(i)%symm == 2 ) then
              per1 = Nfp
              ss1 = 1
           else
              FATAL( CoilSepDeriv1, .true. , Errors in coil symmetry )
           endif
           NS1 = coil(i)%NS
           SALLOCATE(x1, (0:NS1-1), zero)
           SALLOCATE(y1, (0:NS1-1), zero)
           SALLOCATE(z1, (0:NS1-1), zero)
           SALLOCATE(absrp1, (0:NS1-1), zero)
           do j1 = 1, per1
              do l1 = 0, ss1
                 x1(0:NS1-1)=(coil(i)%xx(0:NS1-1))*cos(pi2*(j1-1)/Nfp)-(coil(i)%yy(0:NS1-1))*sin(pi2*(j1-1)/Nfp)
                 y1(0:NS1-1)=((-1.0)**l1)*((coil(i)%yy(0:NS1-1))*cos(pi2*(j1-1)/Nfp)+(coil(i)%xx(0:NS1-1))*sin(pi2*(j1-1)/Nfp))
                 z1(0:NS1-1)=(coil(i)%zz(0:NS1-1))*((-1.0)**l1)
                 if ( ccsep_skip .eq. 1 ) then
                    absrp1(0:NS-1) = sqrt( (coil(i)%xt(0:NS-1))**2 + (coil(i)%yt(0:NS-1))**2 + (coil(i)%zt(0:NS-1))**2 )
                    xc1 = sum( x1(0:NS-1)*absrp1(0:NS-1) ) / sum( absrp1(0:NS-1) )
                    yc1 = sum( y1(0:NS-1)*absrp1(0:NS-1) ) / sum( absrp1(0:NS-1) )
                    zc1 = sum( z1(0:NS-1)*absrp1(0:NS-1) ) / sum( absrp1(0:NS-1) )
                    if ( sqrt((xc0-xc1)**2+(yc0-yc1)**2+(zc0-zc1)**2 ) .ge. (coil(i)%Lo+coil(icoil)%Lo)*.25 ) cycle
                 endif
                 do kseg0 = 0, NS-1
                    if( myid.ne.modulo(kseg0,ncpu) ) cycle ! parallelization loop;
                    do kseg1 = 0, NS1-1
                       rdifff = sqrt( ( x0(kseg0) - x1(kseg1) )**2 + ( y0(kseg0) - y1(kseg1) )**2 + ( z0(kseg0) - z1(kseg1) )**2 )
                       if ( rdifff < r_delta ) then
                          if ( penfun_ccsep .eq. 1 ) then
                             hypc = 0.5*exp(ccsep_alpha*( r_delta - rdifff )) + 0.5*exp(-1.0*ccsep_alpha*( r_delta - rdifff ))
                             hyps = 0.5*exp(ccsep_alpha*( r_delta - rdifff )) - 0.5*exp(-1.0*ccsep_alpha*( r_delta - rdifff ))
                             int_hold = -2.0*ccsep_alpha*(hypc-1.0)*hyps/rdifff
                          else
                             int_hold = -1.0*ccsep_beta*ccsep_alpha*(ccsep_alpha*(r_delta-rdifff))**(ccsep_beta-1.0)/rdifff
                          endif
                          derivsh(:) = derivsh(:) + int_hold*( &
                                  (x0(kseg0)-x1(kseg1))*(DoF(icoil)%xof(kseg0,:)*cos(pi2*(j0-1)/Nfp)-Dof(icoil)%yof(kseg0,:)*sin(pi2*(j0-1)/Nfp)) + &
                                  (y0(kseg0)-y1(kseg1))*((-1.0)**l0)*(DoF(icoil)%yof(kseg0,:)*cos(pi2*(j0-1)/Nfp)+Dof(icoil)%xof(kseg0,:)*sin(pi2*(j0-1)/Nfp)) + &
                                  (z0(kseg0)-z1(kseg1))*Dof(icoil)%zof(kseg0,:)*((-1.0)**l0) )
                       endif
                    enddo
                 enddo
                 call MPI_BARRIER( MPI_COMM_FOCUS, ierr )
                 call MPI_ALLREDUCE( MPI_IN_PLACE, derivsh, ND , MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FOCUS, ierr )
                 derivs_hold(1,:) = derivs_hold(1,:) + derivsh(:)
                 derivsh(:) = zero
              enddo
           enddo
           derivs = derivs + pi2*pi2*derivs_hold/(NS*NS1)
           derivs_hold = zero
           DALLOCATE(x1)
           DALLOCATE(y1)
           DALLOCATE(z1)
           DALLOCATE(absrp1)
        enddo
     enddo
  enddo

  DALLOCATE(derivs_hold)
  DALLOCATE(derivsh)
  DALLOCATE(x0)
  DALLOCATE(y0)
  DALLOCATE(z0)
  DALLOCATE(absrp0)

  return

end subroutine CoilSepDeriv1


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
