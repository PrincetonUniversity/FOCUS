!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

!title (tlength) ! Calculate total length cost functon and its derivatives.

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

!latex \briefly{The constraint on coil length . . . }

!latex \calledby{\link{solvers}}

!latex  \section{General}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

subroutine length( icoil, ideriv )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  use globals, only : zero, half, myid, Ncoils, Nseg, coil, discretecurve, case_length, coillengtherror
!, half, pi2, machprec, ncpu, myid, ounit, coil, DoF, Ncoils, Nfixgeo, Ndof, coillength, t1L, t2L, case_length
  
  implicit none
  
  include "mpif.h"
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  INTEGER, INTENT(in) :: icoil, ideriv
  
  INTEGER             :: astat, ierr, ii
  REAL                :: dlength !d1L(1:Ndof), norm(1:Ncoils)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  coil(icoil)%L = zero
  
  do ii = 0, Nseg-1
   
   coil(icoil)%L = coil(icoil)%L + sqrt( coil(icoil)%xt(ii)**2 + coil(icoil)%yt(ii)**2 + coil(icoil)%zt(ii)**2 )
   
  enddo ! end ii
  
  coil(icoil)%L = coil(icoil)%L ! * discretecurve
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  coil(icoil)%Le = zero ! length error; 04 Sep 17;
  
  select case( case_length )
   
  case( 1 ) ! quadratic; 04 Sep 17;
   
   if ( coil(icoil)%Lf /= 0 ) coil(icoil)%Le = coil(icoil)%Le + half * (coil(icoil)%L - coil(icoil)%Lo)**2 / coil(icoil)%Lo**2
   
  case( 2 ) ! exponential; 04 Sep 17;
   
   if ( coil(icoil)%Lf /= 0 ) coil(icoil)%Le = coil(icoil)%Le + exp(coil(icoil)%L) / exp(coil(icoil)%Lo)
   
  case default
   
   FATAL( length, .true., selected case_length not supported )
   
  end select
  
!  coil(icoil)%L = coil(icoil)%L / ( Ncoils - Nfixgeo + machprec )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  if ( ideriv >= 1 ) then
   
   FATAL( length, .true., under reconstruction )

!     t1L = zero ; d1L = zero ; norm = zero
!
!     idof = 0
!     do icoil = 1, Ncoils
!
!        ND = DoF(icoil)%ND
!
!        if (case_length == 1) then
!           norm(icoil) = (coil(icoil)%L - coil(icoil)%Lo) / coil(icoil)%Lo**2  ! quadratic;
!        elseif (case_length == 2) then
!           norm(icoil) = exp(coil(icoil)%L) / exp(coil(icoil)%Lo)       ! exponential;
!        else
!           FATAL( length, .true. , invalid case_length option )
!        end if
!
!        if ( coil(icoil)%Ic /= 0 ) then !if current is free;
!           idof = idof +1
!        endif
!
!        if ( coil(icoil)%Lf /= 0 ) then !if geometry is free;
!           call lenDeriv1( icoil, d1L(idof+1:idof+ND), ND )
!           t1L(idof+1:idof+ND) = d1L(idof+1:idof+ND) * norm(icoil)
!           idof = idof + ND
!        endif
!
!     enddo !end icoil;
!     FATAL( torflux , idof .ne. Ndof, counting error in packing )
!
!     t1L = t1L / (Ncoils - Nfixgeo + machprec)
!
  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

  return

end subroutine length

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

!subroutine dlength(icoil, locallength)
!
!  use globals, only: zero, coil, myid, ounit, Ncoils
!  implicit none
!  include "mpif.h"
!
!  INTEGER, intent(in)  :: icoil
!  REAL   , intent(out) :: locallength
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!  INTEGER              :: kseg, astat, ierr
!  REAL                 :: dlength
!
!  FATAL( dlength, icoil .lt. 1 .or. icoil .gt. Ncoils, icoil not in right range )
!  
!  locallength = zero
!  
!  do kseg = 0, Nseg-1
!
!     dlength = sqrt(coil(icoil)%xt(kseg)**2 + coil(icoil)%yt(kseg)**2 + coil(icoil)%zt(kseg)**2)
!     locallength  = locallength + dlength * coil(icoil)%dd(kseg)
!
!  enddo ! end kseg
!
!  return
!
!end subroutine dlength

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!subroutine LenDeriv1(icoil, derivs, ND)
!
!  use globals, only: zero, pi2, coil, DoF, myid, ounit, Ncoils
!  implicit none
!  include "mpif.h"
!
!  INTEGER, intent(in)  :: icoil, ND
!  REAL   , intent(out) :: derivs(1:1, 1:ND)
!
!  INTEGER              :: kseg, astat, ierr
!  REAL                 :: dl3, xt, yt, zt, xa, ya, za
!  REAL, dimension(1:1, 1:Nseg) :: dLx, dLy, dLz
!
!  FATAL( LenDeriv1, icoil .lt. 1 .or. icoil .gt. Ncoils, icoil not in right range )
!  
!  derivs = zero
!  
!  do kseg = 1, Nseg
!     
!     xt = coil(icoil)%xt(kseg) ; yt = coil(icoil)%yt(kseg) ; zt = coil(icoil)%zt(kseg)
!     xa = coil(icoil)%xa(kseg) ; ya = coil(icoil)%ya(kseg) ; za = coil(icoil)%za(kseg)
!
!     dl3 = sqrt(xt*xt + yt*yt + zt*zt)**3
!
!     dLx(1,kseg) = ( yt*ya*xt + zt*za*xt - yt*yt*xa - zt*zt*xa ) / dl3 * coil(icoil)%dd(kseg)
!     dLy(1,kseg) = ( xt*xa*yt + zt*za*yt - xt*xt*ya - zt*zt*ya ) / dl3 * coil(icoil)%dd(kseg)
!     dLz(1,kseg) = ( xt*xa*zt + yt*ya*zt - xt*xt*za - yt*yt*za ) / dl3 * coil(icoil)%dd(kseg)
!
!  enddo ! end kseg
!
!  derivs(1:1, 1:ND) = matmul(dLx, DoF(icoil)%xof) + matmul(dLy, DoF(icoil)%yof) + matmul(dLz, DoF(icoil)%zof)
!
!  return
!
!end subroutine LenDeriv1


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
