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
       iccsep, mccsep, LM_fvec, LM_fjac, weight_ccsep, &
       MPI_COMM_FOCUS

  implicit none
  include "mpif.h"
  INTEGER, INTENT(in) :: ideriv

  INTEGER             :: astat, ierr, icoil, idof, ND, ivec, numCoil
  REAL                :: d1C(1:Ndof), norm(1:Ncoils), coilsep0

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  ccsep = zero
  !coilsep0 = zero
  ivec = 1
  numCoil = zero
 
  ! What if Ncoils is equal to 1???
  if(Ncoils .eq. 1) return ! Maybe change 

  ! Do calculation for free and fixed coils 
  if( ideriv >= 0 ) then
     do icoil = 1, Ncoils     !only care about unique coils;
        if(coil(icoil)%type == 1) then  ! only for Fourier, Probably should delete this 
           !if( myid.ne.modulo(icoil-1,ncpu) ) cycle ! parallelization loop;
           call CoilSepDeriv0(icoil, coilsep0)
           if ( coil(icoil)%symm == 0 ) numCoil = numCoil + 1
           if ( coil(icoil)%symm == 1 ) numCoil = numCoil + Nfp
           if ( coil(icoil)%symm == 2 ) numCoil = numCoil + 2*Nfp
           ccsep = ccsep + coilsep0
           if ( coil(icoil)%Lc /= 0 ) then
              if (mccsep > 0) then ! L-M format of targets
                 LM_fvec(iccsep+ivec) = weight_ccsep * coilsep0
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
     t1C = zero ; d1C = zero ; norm = zero
     idof = 0 ; ivec = 1
     do icoil = 1, Ncoils
        ND = DoF(icoil)%ND
        if ( coil(icoil)%Ic /= 0 ) then !if current is free;
           idof = idof +1
        endif
        if ( coil(icoil)%Lc /= 0 ) then !if geometry is free;
           if(coil(icoil)%type .eq. 1) then ! only for Fourier
              ! calculate normalization
              !norm(icoil) = (coil(icoil)%L - coil(icoil)%Lo) / coil(icoil)%Lo**2  ! quadratic;
              ! call lederiv1 to calculate the 1st derivatives
              call CoilSepDeriv1( icoil, d1C(idof+1:idof+ND), ND )
              t1C(idof+1:idof+ND) = d1C(idof+1:idof+ND) * norm(icoil)
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

     t1C = t1C / (Ncoils - Nfixgeo + machprec)
  endif

  return
end subroutine coilsep

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! calculate coil length

subroutine CoilSepDeriv0(icoil, coilsep0)

  use globals, only: dp, zero, coil, myid, ounit, Ncoils, MPI_COMM_FOCUS 
  implicit none
  include "mpif.h"

  INTEGER, intent(in)  :: icoil
  REAL   , intent(out) :: coilsep0
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  INTEGER              :: kseg, astat, ierr, iicoil
  REAL                 :: dlength

  FATAL( LenDeriv0, icoil .lt. 1 .or. icoil .gt. Ncoils, icoil not in right range )
  
  coilsep0 = zero 

  ! Loop over symmetric coils
  ! Loop over nonsymmetric coils for each symmetric coil 

  do iicoil = icoil+1, Ncoils 
     do kseg = 0, coil(icoil)%NS-1

     enddo
  enddo 

  return

end subroutine CoilSepDeriv0

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!calculate coil length in derivs(0) and first derivatives in derivs(1:Cdof)

subroutine CoilSepDeriv1(icoil, derivs, ND)

  use globals, only: dp, zero, pi2, coil, DoF, myid, ounit, Ncoils, MPI_COMM_FOCUS
  implicit none
  include "mpif.h"

  INTEGER, intent(in)  :: icoil, ND
  REAL   , intent(out) :: derivs(1:1, 1:ND)
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  INTEGER              :: kseg, astat, ierr
  REAL                 :: dl3, xt, yt, zt, xa, ya, za
  REAL, dimension(1:1, 0:coil(icoil)%NS-1) :: dLx, dLy, dLz

  FATAL( LenDeriv1, icoil .lt. 1 .or. icoil .gt. Ncoils, icoil not in right range )
  
  derivs = zero
  
  ! Loop through all symmetric cases separations first 
  ! Loop through all coils for each symmetric coil

  do kseg = 0, coil(icoil)%NS-1
     
     xt = coil(icoil)%xt(kseg) ; yt = coil(icoil)%yt(kseg) ; zt = coil(icoil)%zt(kseg)
     
     !dl3 = sqrt(xt*xt + yt*yt + zt*zt)**3

     !dLx(1,kseg) = ( yt*ya*xt + zt*za*xt - yt*yt*xa - zt*zt*xa ) / dl3 * coil(icoil)%dd(kseg)
     !dLy(1,kseg) = ( xt*xa*yt + zt*za*yt - xt*xt*ya - zt*zt*ya ) / dl3 * coil(icoil)%dd(kseg)
     !dLz(1,kseg) = ( xt*xa*zt + yt*ya*zt - xt*xt*za - yt*yt*za ) / dl3 * coil(icoil)%dd(kseg)

  enddo ! end kseg

  !derivs(1:1, 1:ND) = matmul(dLx, DoF(icoil)%xof) + matmul(dLy, DoF(icoil)%yof) + matmul(dLz, DoF(icoil)%zof)

  return

end subroutine CoilSepDeriv1


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
