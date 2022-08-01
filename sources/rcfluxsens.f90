!title (rcfluxsens) ! Calculate denominator of coil tolerance functional (tkruger)
subroutine rcfluxsens(ideriv)
  use globals, only: dp, zero, half, pi2, machprec, ncpu, myid, ounit, &
       coil, DoF, Ncoils, Nfixgeo, Ndof, LM_fvec, LM_fjac, Nfp, MPI_COMM_FOCUS
  use mpi
  implicit none
  INTEGER, INTENT(in) :: ideriv
  INTEGER             :: astat, ierr, icoil, idof, ND, ivec, mult
  REAL                :: d1Z(1:Ndof), norm(1:Ncoils), t1Z(1:Ndof) ! Change this
  REAL                :: absdpsihold, absdpsi ! Change this

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  ! Solve for periodic field lines?

  absdpsi = zero

  if( ideriv >= 0 ) then
     do icoil = 1, Ncoils     !only care about unique coils;
        if(coil(icoil)%type == 1) then  ! only for Fourier
           !if( myid.ne.modulo(icoil-1,ncpu) ) cycle ! parallelization loop;
           call dpsi0(icoil, absdpsihold)
           !RlBCAST( coil(icoil)%L, 1, modulo(icoil-1,ncpu) ) !broadcast each coil's length
           ! Multiply absdpsihold by coil symmetry
           select case (coil(icoil)%symm)
           case ( 0 )
              mult = 1
           case ( 1 )
              mult = Nfp
           case ( 2)
              mult = Nfp*2
           end select
           absdpsi = absdpsi + absdpsihold*real(mult)
        endif
     enddo
  endif

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if ( ideriv >= 1 ) then
     t1Z = zero
     idof = 0
     do icoil = 1, Ncoils
        ND = DoF(icoil)%ND
        if ( coil(icoil)%Ic /= 0 ) then !if current is free;
           ! Include this later
           idof = idof + 1
        endif
        if ( coil(icoil)%Lc /= 0 ) then !if geometry is free;
           select case (coil(icoil)%symm)
           case ( 0 )
              mult = 1
           case ( 1 )
              mult = Nfp
           case ( 2)
              mult = Nfp*2
           end select
           if(coil(icoil)%type .eq. 1) then ! only for Fourier
              call dpsi1( icoil, d1Z(idof+1:idof+ND), ND )
              !t1Z(idof+1:idof+ND) = d1Z(idof+1:idof+ND) * real(mult)
           endif 
           idof = idof + ND
        endif
     enddo !end icoil;
     FATAL( length , idof .ne. Ndof, counting error in packing )
  endif
  return
end subroutine rcfluxsens

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine dpsi0(icoil, absdpsi)

  use globals, only: dp, zero, coil, gsurf, myid, ounit, Ncoils, resbn_m, pi2, MPI_COMM_FOCUS 
  implicit none
  include "mpif.h"

  INTEGER, intent(in)  :: icoil
  REAL   , intent(out) :: absdpsi
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  INTEGER              :: i, NS, astat, ierr
  REAL, dimension(0:coil(icoil)%NS) :: dAxx, dAxy, dAxz, dAyx, dAyy, dAyz, dAzx, dAzy, dAzz

  FATAL( dpsi0, icoil .lt. 1 .or. icoil .gt. Ncoils, icoil not in right range )
  
  coil(icoil)%dpsidx = zero
  coil(icoil)%dpsidy = zero
  coil(icoil)%dpsidz = zero
  NS = coil(icoil)%NS
  do i = 1, gsurf(1)%Nseg_stable-1
     call deltaafield(icoil, gsurf(1)%ox(i), gsurf(1)%oy(i), gsurf(1)%oz(i), dAxx, dAxy, dAxz, dAyx, dAyy, dAyz, dAzx, dAzy, dAzz)
     coil(icoil)%dpsidx(0:NS) = coil(icoil)%dpsidx(0:NS) + dAxx(0:NS)*gsurf(1)%oxdot(i) + dAyx(0:NS)*gsurf(1)%oydot(i) + dAzx(0:NS)*gsurf(1)%ozdot(i)
     coil(icoil)%dpsidy(0:NS) = coil(icoil)%dpsidy(0:NS) + dAxy(0:NS)*gsurf(1)%oxdot(i) + dAyy(0:NS)*gsurf(1)%oydot(i) + dAzy(0:NS)*gsurf(1)%ozdot(i)
     coil(icoil)%dpsidz(0:NS) = coil(icoil)%dpsidz(0:NS) + dAxz(0:NS)*gsurf(1)%oxdot(i) + dAyz(0:NS)*gsurf(1)%oydot(i) + dAzz(0:NS)*gsurf(1)%ozdot(i)
     call deltaafield(icoil, gsurf(1)%xx(i), gsurf(1)%xy(i), gsurf(1)%xz(i), dAxx, dAxy, dAxz, dAyx, dAyy, dAyz, dAzx, dAzy, dAzz)
     coil(icoil)%dpsidx(0:NS) = coil(icoil)%dpsidx(0:NS) - dAxx(0:NS)*gsurf(1)%xxdot(i) - dAyx(0:NS)*gsurf(1)%xydot(i) - dAzx(0:NS)*gsurf(1)%xzdot(i)
     coil(icoil)%dpsidy(0:NS) = coil(icoil)%dpsidy(0:NS) - dAxy(0:NS)*gsurf(1)%xxdot(i) - dAyy(0:NS)*gsurf(1)%xydot(i) - dAzy(0:NS)*gsurf(1)%xzdot(i)
     coil(icoil)%dpsidz(0:NS) = coil(icoil)%dpsidz(0:NS) - dAxz(0:NS)*gsurf(1)%xxdot(i) - dAyz(0:NS)*gsurf(1)%xydot(i) - dAzz(0:NS)*gsurf(1)%xzdot(i)
  enddo
  coil(icoil)%dpsidx(0:NS) = pi2*resbn_m*coil(icoil)%dpsidx(0:NS)/real(gsurf(1)%Nseg_stable-1)
  coil(icoil)%dpsidy(0:NS) = pi2*resbn_m*coil(icoil)%dpsidy(0:NS)/real(gsurf(1)%Nseg_stable-1)
  coil(icoil)%dpsidz(0:NS) = pi2*resbn_m*coil(icoil)%dpsidz(0:NS)/real(gsurf(1)%Nseg_stable-1)
  absdpsi = sum( sqrt( coil(icoil)%dpsidx(1:NS)**2 + coil(icoil)%dpsidy(1:NS)**2 + coil(icoil)%dpsidz(1:NS)**2 ) )/real(NS)
  
  return

end subroutine dpsi0

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine dpsi1(icoil, derivs, ND)

  use globals, only: dp, zero, pi2, coil, DoF, myid, ounit, Ncoils, MPI_COMM_FOCUS
  implicit none
  include "mpif.h"

  INTEGER, intent(in)  :: icoil, ND
  REAL   , intent(out) :: derivs(1:1, 1:ND)
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  INTEGER              :: kseg, astat, ierr
  REAL                 :: xt, yt, zt, xa, ya, za
!  REAL, dimension(1:1, 0:coil(icoil)%NS-1) :: dLx, dLy, dLz

  FATAL( dpsi1, icoil .lt. 1 .or. icoil .gt. Ncoils, icoil not in right range )
  
  derivs = zero
  
  do kseg = 0, coil(icoil)%NS-1
     
     xt = coil(icoil)%xt(kseg) ; yt = coil(icoil)%yt(kseg) ; zt = coil(icoil)%zt(kseg)
     xa = coil(icoil)%xa(kseg) ; ya = coil(icoil)%ya(kseg) ; za = coil(icoil)%za(kseg)

!     dl3 = sqrt(xt*xt + yt*yt + zt*zt)**3

!     dLx(1,kseg) = ( yt*ya*xt + zt*za*xt - yt*yt*xa - zt*zt*xa ) / dl3 * coil(icoil)%dd(kseg)
!     dLy(1,kseg) = ( xt*xa*yt + zt*za*yt - xt*xt*ya - zt*zt*ya ) / dl3 * coil(icoil)%dd(kseg)
!     dLz(1,kseg) = ( xt*xa*zt + yt*ya*zt - xt*xt*za - yt*yt*za ) / dl3 * coil(icoil)%dd(kseg)

  enddo ! end kseg

!  derivs(1:1, 1:ND) = matmul(dLx, DoF(icoil)%xof) + matmul(dLy, DoF(icoil)%yof) + matmul(dLz, DoF(icoil)%zof)

  return

end subroutine dpsi1


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
