!title (coilsep) ! Calculate coil separation objective functon and its derivatives. (tkruger)

!latex \briefly{The constraint on coil to coil separation solves for coils that have
!latex         reasonable finite builds and can be assembled without interferance.
!latex         This function is still under development. emph{targt\_length}.}

!latex \calledby{\link{solvers}}

!latex  \section{General}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! ccsep is total penalty
! chi = chi + weight_ccsep*ccsep
! t1CU is total derivative of penalty
! LM implemented
! not parallelized, at some point check to see how long takes to run
subroutine coilsep(ideriv)
  use globals, only: dp, zero, half, pi2, machprec, ncpu, myid, ounit, &
       coil, DoF, Ncoils, Nfixgeo, Ndof, ccsep, t1C, t2C, Nfp, &
       iccsep, mccsep, LM_fvec, LM_fjac, weight_ccsep, FouCoil, &
       MPI_COMM_FOCUS

  implicit none
  include "mpif.h"
  INTEGER, INTENT(in) :: ideriv

  INTEGER             :: astat, ierr, icoil, idof, ND, NF, ivec, numCoil
  REAL                :: d1C(1:Ndof), ccsep0

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  ccsep = zero
  ivec = 1
  numCoil = zero
  ccsep0 = zero
 
  ! Do calculation for free and fixed coils 
  if( ideriv >= 0 ) then
     if ( Ncoils .eq. 1 ) return
     do icoil = 1, Ncoils     !only care about unique coils;
        if ( coil(icoil)%type == 1 ) then  ! only for Fourier, Probably should delete this 
           !if( myid.ne.modulo(icoil-1,ncpu) ) cycle ! parallelization loop;
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
     t1C = zero ; d1C = zero !; norm = zero
     !if(myid .eq. 0) write(ounit, '("HIIIIIIIIII")')
     if ( Ncoils .eq. 1 ) return
     idof = 0 ; ivec = 1
     do icoil = 1, Ncoils
        ND = DoF(icoil)%ND
        if ( coil(icoil)%Ic /= 0 ) then !if current is free;
           idof = idof +1
        endif
        if ( coil(icoil)%Lc /= 0 ) then !if geometry is free;
           if(coil(icoil)%type .eq. 1) then ! only for Fourier
              NF = FouCoil(icoil)%NF
              ! calculate normalization
              !norm(icoil) = (coil(icoil)%L - coil(icoil)%Lo) / coil(icoil)%Lo**2  ! quadratic;
              ! call lederiv1 to calculate the 1st derivatives
              call CoilSepDeriv1( icoil, d1C(idof+1:idof+ND), ND , NF )
              !t1C(idof+1:idof+ND) = d1C(idof+1:idof+ND) * norm(icoil)
              t1C(idof+1:idof+ND) = d1C(idof+1:idof+ND) * 1.0
              ! I THINK LINE ABOVE IS WRONG, DO NOT THINK IT IS NEEDED, UNLESS NORM IS COEFFICENT
              if (mccsep > 0) then ! L-M format of targets
                 LM_fjac(iccsep+ivec, idof+1:idof+ND) = weight_ccsep * d1C(idof+1:idof+ND)
                 !if (case_length == 2) &
                 !     & LM_fjac(ivec, idof+1:idof+ND) = LM_fjac(ivec, idof+1:idof+ND) &
                 !     & * exp(coil(icoil)%L) / exp(coil(icoil)%Lo)
                 ivec = ivec + 1
              endif
           endif 
           idof = idof + ND
        endif

     enddo !end icoil;
     FATAL( coilsep , idof .ne. Ndof, counting error in packing )

     if (mccsep > 0) then ! L-M format of targets
        FATAL( coilsep, ivec == mccsep, Errors in counting ivec for L-M )
     endif
     !t1C = t1C / (Ncoils - Nfixgeo + machprec)
     t1C = 2.0*t1C / ( numCoil * (numCoil - 1) + machprec )
  endif

  return
end subroutine coilsep

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! calculate coil length

subroutine CoilSepDeriv0(icoil, coilsep0)

  use globals, only: dp, zero, pi2, coil, myid, ounit, Ncoils, Nfp, machprec, r_delta, ccsep_alpha, MPI_COMM_FOCUS 
  implicit none
  include "mpif.h"

  INTEGER, intent(in)  :: icoil
  REAL   , intent(out) :: coilsep0
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  INTEGER              :: kseg0, kseg1, astat, ierr, i, per0, per1, ss0, ss1, j0, j1, l0, l1
  REAL                 :: rdiff, H, hypc, coilsepHold 
  REAL, allocatable    :: x0(:), y0(:), z0(:), x1(:), y1(:), z1(:)

  !SALLOCATE(x0, (1:coil(icoil)%NS), zero)
  !SALLOCATE(y0, (1:coil(icoil)%NS), zero)
  !SALLOCATE(z0, (1:coil(icoil)%NS), zero)
  SALLOCATE(x0, (0:coil(icoil)%NS-1), zero)
  SALLOCATE(y0, (0:coil(icoil)%NS-1), zero)
  SALLOCATE(z0, (0:coil(icoil)%NS-1), zero)

  FATAL( CoilSepDeriv0, icoil .lt. 1 .or. icoil .gt. Ncoils, icoil not in right range )

  coilsep0 = zero
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
     FATAL( CoilSepDeriv0, coil(icoil)%symm == 3, Errors in coil symmetry )
  endif
  do j0 = 1, per0
     do l0 = 0, ss0
        !x0 = (coil(icoil)%xx)*cos(pi2*(j0-1)/Nfp) - (coil(icoil)%yy)*sin(pi2*(j0-1)/Nfp)
        !y0 = (coil(icoil)%yy)*cos(pi2*(j0-1)/Nfp) + ((-1.0)**ss0)*(coil(icoil)%xx)*sin(pi2*(j0-1)/Nfp)
        !z0 = (coil(icoil)%zz)*(-1.0)**ss0
        x0(0:coil(icoil)%NS-1) = (coil(icoil)%xx(0:coil(icoil)%NS-1))*cos(pi2*(j0-1)/Nfp) - (coil(icoil)%yy(0:coil(icoil)%NS-1))*sin(pi2*(j0-1)/Nfp)
        y0(0:coil(icoil)%NS-1) = ((-1.0)**l0)*((coil(icoil)%yy(0:coil(icoil)%NS-1))*cos(pi2*(j0-1)/Nfp) + (coil(icoil)%xx(0:coil(icoil)%NS-1))*sin(pi2*(j0-1)/Nfp))
        z0(0:coil(icoil)%NS-1) = (coil(icoil)%zz(0:coil(icoil)%NS-1))*((-1.0)**l0)
        ! I THINK THIS IS RIGHT
        do i = icoil+1, Ncoils 
           if (coil(icoil)%Lc == 0 .and. coil(i)%Lc == 0) then
              ! Do nothing
              !if(myid .eq. 0) write(ounit, '(8X": The minimum BLAH distance is "4X" :" I3" ; at coils")') ss0
           else
              !SALLOCATE(x1, (1:coil(i)%NS), zero)
              !SALLOCATE(y1, (1:coil(i)%NS), zero)
              !SALLOCATE(z1, (1:coil(i)%NS), zero)
              SALLOCATE(x1, (0:coil(i)%NS-1), zero)
              SALLOCATE(y1, (0:coil(i)%NS-1), zero)
              SALLOCATE(z1, (0:coil(i)%NS-1), zero)
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
                 FATAL( CoilSepDeriv0, coil(i)%symm == 3, Errors in coil symmetry )
              endif
              do j1 = 1, per1
                 do l1 = 0, ss1
                    x1(0:coil(i)%NS-1) = (coil(i)%xx(0:coil(i)%NS-1))*cos(pi2*(j1-1)/Nfp) - (coil(i)%yy(0:coil(i)%NS-1))*sin(pi2*(j1-1)/Nfp)
                    y1(0:coil(i)%NS-1) = ((-1.0)**l1)*((coil(i)%yy(0:coil(i)%NS-1))*cos(pi2*(j1-1)/Nfp) + (coil(i)%xx(0:coil(i)%NS-1))*sin(pi2*(j1-1)/Nfp))
                    z1(0:coil(i)%NS-1) = (coil(i)%zz(0:coil(i)%NS-1))*((-1.0)**l1)
                    do kseg0 = 1, coil(icoil)%NS
                       do kseg1 = 1, coil(i)%NS
                          rdiff = sqrt( ( x0(kseg0) - x1(kseg1) )**2 + ( y0(kseg0) - y1(kseg1) )**2 + ( z0(kseg0) - z1(kseg1) )**2 )
                          if ( rdiff < r_delta ) H = 1.0
                          if ( rdiff .ge. r_delta ) H = 0.0
                          hypc = 0.5*exp(ccsep_alpha*( r_delta - rdiff )) + 0.5*exp(-1.0*ccsep_alpha*( r_delta - rdiff ))
                          coilsepHold = coilsepHold + H*((hypc - 1.0)**2)
                       enddo
                    enddo
                 enddo
              enddo
              !coilsep0 = pi2*pi2*coilsepHold/((coil(icoil)%NS)*(coil(i)%NS)+machprec)
              coilsep0 = coilsep0 + pi2*pi2*coilsepHold/((coil(icoil)%NS)*(coil(i)%NS)+machprec)
              coilsepHold = 0.0
           endif
        enddo
     enddo
  enddo

  return

end subroutine CoilSepDeriv0

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!calculate coil length in derivs(0) and first derivatives in derivs(1:Cdof)

subroutine CoilSepDeriv1(icoil, derivs, ND, NF)

  use globals, only: dp, zero, pi2, coil, DoF, myid, ounit, machprec, Ncoils, Nfp, r_delta, ccsep_alpha, FouCoil, MPI_COMM_FOCUS
  implicit none
  include "mpif.h"

  INTEGER, intent(in)  :: icoil, ND
  REAL   , intent(out) :: derivs(1:1, 1:ND)
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  INTEGER              :: kseg0, kseg1, astat, ierr, per0, per1, ss0, ss1, j0, j1, l0, l1, i, nff, NF
  REAL                 :: dl3, xt, yt, zt, H, hypc, hyps, rdiff, int_hold
  REAL, allocatable    :: x0(:), y0(:), z0(:), x1(:), y1(:), z1(:), derivs_hold(:,:)

  FATAL( CoilSepDeriv1, icoil .lt. 1 .or. icoil .gt. Ncoils, icoil not in right range )
  
  derivs = zero
  !derivs_hold = derivs

  SALLOCATE(derivs_hold, (1:1,1:3+6*NF), zero)

  !SALLOCATE(x0, (1:coil(icoil)%NS), zero)
  !SALLOCATE(y0, (1:coil(icoil)%NS), zero)
  !SALLOCATE(z0, (1:coil(icoil)%NS), zero)
  SALLOCATE(x0, (0:coil(icoil)%NS-1), zero)
  SALLOCATE(y0, (0:coil(icoil)%NS-1), zero)
  SALLOCATE(z0, (0:coil(icoil)%NS-1), zero)

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
     FATAL( CoilSepDeriv1, coil(icoil)%symm == 3, Errors in coil symmetry )
  endif 
  do j0 = 1, per0
     do l0 = 0, ss0
        !x0 = (coil(icoil)%xx)*cos(pi2*(j0-1)/Nfp) - (coil(icoil)%yy)*sin(pi2*(j0-1)/Nfp)
        !x0(1:coil(icoil)%NS) = (coil(icoil)%xx(1:coil(icoil)%NS))*cos(pi2*(j0-1)/Nfp) - (coil(icoil)%yy(1:coil(icoil)%NS))*sin(pi2*(j0-1)/Nfp)
        x0(0:coil(icoil)%NS-1) = (coil(icoil)%xx(0:coil(icoil)%NS)-1)*cos(pi2*(j0-1)/Nfp) - (coil(icoil)%yy(0:coil(icoil)%NS-1))*sin(pi2*(j0-1)/Nfp)
        !y0 = (coil(icoil)%yy)*cos(pi2*(j0-1)/Nfp) + ((-1.0)**ss0)*(coil(icoil)%xx)*sin(pi2*(j0-1)/Nfp)
        !y0(1:coil(icoil)%NS) = (coil(icoil)%yy(1:coil(icoil)%NS))*cos(pi2*(j0-1)/Nfp) + ((-1.0)**ss0)*(coil(icoil)%xx(1:coil(icoil)%NS))*sin(pi2*(j0-1)/Nfp)
        y0(0:coil(icoil)%NS-1) = ((-1.0)**l0)*((coil(icoil)%yy(0:coil(icoil)%NS-1))*cos(pi2*(j0-1)/Nfp) + (coil(icoil)%xx(0:coil(icoil)%NS-1))*sin(pi2*(j0-1)/Nfp))
        !z0 = coil(icoil)%zz*(-1.0)**ss0
        !z0(1:coil(icoil)%NS) = coil(icoil)%zz(1:coil(icoil)%NS)*(-1.0)**ss0
        z0(0:coil(icoil)%NS-1) = (coil(icoil)%zz(0:coil(icoil)%NS-1))*((-1.0)**l0)
        ! I THINK THIS IS RIGHT
        do i = 1, Ncoils
           if ( i == icoil ) cycle
           if (coil(icoil)%Lc == 0 .and. coil(i)%Lc == 0) then
              ! Do nothing
           else
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
                 FATAL( CoilSepDeriv1, coil(i)%symm == 3, Errors in coil symmetry )
              endif
              !SALLOCATE(x1, (1:coil(i)%NS), zero)
              !SALLOCATE(y1, (1:coil(i)%NS), zero)
              !SALLOCATE(z1, (1:coil(i)%NS), zero)
              SALLOCATE(x1, (0:coil(i)%NS-1), zero)
              SALLOCATE(y1, (0:coil(i)%NS-1), zero)
              SALLOCATE(z1, (0:coil(i)%NS-1), zero)
              do j1 = 1, per1
                 do l1 = 0, ss1
                    !x1 = (coil(i)%xx)*cos(pi2*(j1-1)/Nfp) - (coil(i)%yy)*sin(pi2*(j1-1)/Nfp)
                    !x1(1:coil(i)%NS)=(coil(i)%xx(1:coil(i)%NS))*cos(pi2*(j1-1)/Nfp)-(coil(i)%yy(1:coil(i)%NS))*sin(pi2*(j1-1)/Nfp)
                    x1(0:coil(i)%NS-1)=(coil(i)%xx(0:coil(i)%NS-1))*cos(pi2*(j1-1)/Nfp)-(coil(i)%yy(0:coil(i)%NS-1))*sin(pi2*(j1-1)/Nfp)
                    !y1 = (coil(i)%yy)*cos(pi2*(j1-1)/Nfp) + ((-1.0)**ss1)*(coil(i)%xx)*sin(pi2*(j1-1)/Nfp)
                    !y1(1:coil(i)%NS)=(coil(i)%yy(1:coil(i)%NS))*cos(pi2*(j1-1)/Nfp)+((-1.0)**ss1)*(coil(i)%xx(1:coil(i)%NS))*sin(pi2*(j1-1)/Nfp)
                    y1(0:coil(i)%NS-1)=((-1.0)**l1)*((coil(i)%yy(0:coil(i)%NS-1))*cos(pi2*(j1-1)/Nfp)+(coil(i)%xx(0:coil(i)%NS-1))*sin(pi2*(j1-1)/Nfp))
                    !z1 = coil(i)%zz*(-1.0)**ss1
                    !z1(1:coil(i)%NS)=coil(i)%zz(1:coil(i)%NS)*(-1.0)**ss1
                    z1(0:coil(i)%NS-1)=(coil(i)%zz(0:coil(i)%NS-1))*((-1.0)**l1)
                    do kseg0 = 1, coil(icoil)%NS
                       do kseg1 = 1, coil(i)%NS
                          rdiff = sqrt( ( x0(kseg0) - x1(kseg1) )**2 + ( y0(kseg0) - y1(kseg1) )**2 + ( z0(kseg0) - z1(kseg1) )**2 )
                          !if(myid .eq. 0) write(ounit, '(8X": The minimum BLAH distance is "4X" :" ES23.15" ; at coils")') rdiff
                          if ( rdiff < r_delta ) then
                             H = 1.0
                             !if(myid .eq. 0) write(ounit, '(8X": The minimum BLAH distance is "4X" :" ES23.15" ; at coils")') rdiff
                          else
                             H = 0.0
                          endif
                          hypc = 0.5*exp(ccsep_alpha*( r_delta - rdiff )) + 0.5*exp(-1.0*ccsep_alpha*( r_delta - rdiff ))
                          hyps = 0.5*exp(ccsep_alpha*( r_delta - rdiff )) - 0.5*exp(-1.0*ccsep_alpha*( r_delta - rdiff ))
                          int_hold = -2.0*ccsep_alpha*H*(hypc-1.0)*hyps/rdiff
                          !if(myid .eq. 0) write(ounit, '(8X": The minimum BLAH distance is "4X" :" ES23.15" ; at coils")'),z0(kseg0)-coil(icoil)%zz(kseg0-1)
                          !if(myid .eq. 0) write(ounit, '(8X": The minimum BLAH distance is "4X" :" ES23.15" ; at coils")'),x1(kseg1)-coil(i)%xx(kseg1)
                          !do nff = 0,NF
                             !ncc = -1.0*nff    *FouCoil(icoil)%smt(kseg,nff)
                             !ncp = -1.0*nff*nff*FouCoil(icoil)%cmt(kseg,nff)
                             !nss =      nff*FouCoil(icoil)%cmt(kseg,nff)
                             !nsp = -1.0*nff*nff*FouCoil(icoil)%smt(kseg,nff)
                             ! Xc
                             !derivs_hold(1,1+nff)      = derivs_hold(1,1+nff)      + int_hold*(x0(kseg0)-x1(kseg1))*FouCoil(icoil)%cmt(kseg0,nff)
                             ! Xs
                             !derivs_hold(1,1+nff+NF)   = derivs_hold(1,1+nff+NF)   + int_hold*(x0(kseg0)-x1(kseg1))*FouCoil(icoil)%smt(kseg0,nff)
                             ! Yc  
                             !derivs_hold(1,2+nff+2*NF) = derivs_hold(1,2+nff+2*NF) + int_hold*(y0(kseg0)-y1(kseg1))*FouCoil(icoil)%cmt(kseg0,nff)
                             ! Ys
                             !derivs_hold(1,2+nff+3*NF) = derivs_hold(1,2+nff+3*NF) + int_hold*(y0(kseg0)-y1(kseg1))*FouCoil(icoil)%smt(kseg0,nff)
                             ! Zc
                             !derivs_hold(1,3+nff+4*NF) = derivs_hold(1,3+nff+4*NF) + int_hold*(z0(kseg0)-z1(kseg1))*FouCoil(icoil)%cmt(kseg0,nff)
                             ! Zs
                             !derivs_hold(1,3+nff+5*NF) = derivs_hold(1,3+nff+5*NF) + int_hold*(z0(kseg0)-z1(kseg1))*FouCoil(icoil)%smt(kseg0,nff)
                             ! TK: NEED TO FIX THIS FOR ANY SYMMETRIC CASE
                          !enddo
                          do nff = 0, NF
                             ! Xc
                             derivs_hold(1,1+nff)=derivs_hold(1,1+nff)+int_hold*((x0(kseg0)-x1(kseg1))*FouCoil(icoil)%cmt(kseg0,nff)*cos(pi2*(j0-1)/Nfp)-(y0(kseg0)-y1(kseg1))*FouCoil(icoil)%cmt(kseg0,nff)*sin(pi2*(j0-1)/Nfp))
                             ! Yc
                             derivs_hold(1,2+nff+2*NF)=derivs_hold(1,2+nff+2*NF)+int_hold*((y0(kseg0)-y1(kseg1))*FouCoil(icoil)%cmt(kseg0,nff)*cos(pi2*(j0-1)/Nfp)+(x0(kseg0)-x1(kseg1))*FouCoil(icoil)%cmt(kseg0,nff)*((-1.0)**ss0)*sin(pi2*(j0-1)/Nfp))
                             ! Zc
                             derivs_hold(1,3+nff+4*NF)=derivs_hold(1,3+nff+4*NF)+int_hold*(z0(kseg0)-z1(kseg1))*FouCoil(icoil)%cmt(kseg0,nff)*(-1.0)**ss1
                          enddo
                          do nff = 1, NF
                             ! Xs
                             derivs_hold(1,1+nff+NF)=derivs_hold(1,1+nff+NF)+int_hold*((x0(kseg0)-x1(kseg1))*FouCoil(icoil)%smt(kseg0,nff)*cos(pi2*(j0-1)/Nfp)-(y0(kseg0)-y1(kseg1))*FouCoil(icoil)%smt(kseg0,nff)*sin(pi2*(j0-1)/Nfp))
                             ! Ys
                             derivs_hold(1,2+nff+3*NF)=derivs_hold(1,2+nff+3*NF)+int_hold*((y0(kseg0)-y1(kseg1))*FouCoil(icoil)%smt(kseg0,nff)*cos(pi2*(j0-1)/Nfp)+(x0(kseg0)-x1(kseg1))*FouCoil(icoil)%smt(kseg0,nff)*((-1.0)**ss0)*sin(pi2*(j0-1)/Nfp))
                             ! Zs
                             derivs_hold(1,3+nff+5*NF)=derivs_hold(1,3+nff+5*NF)+int_hold*(z0(kseg0)-z1(kseg1))*FouCoil(icoil)%smt(kseg0,nff)*(-1.0)**ss1
                             ! OTHER FIELD PERIODS MIGHT HAVE NEGATIVE ERROR AND CANCEL 
                          enddo
                       enddo
                    enddo
                 enddo
              enddo
              derivs = derivs + pi2*pi2*derivs_hold/((coil(icoil)%NS)*(coil(i)%NS)+machprec)
              derivs_hold = zero
           endif
        enddo
     enddo
  enddo

  return

end subroutine CoilSepDeriv1


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
