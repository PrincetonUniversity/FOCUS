!title (rcfluxsens) ! Calculate denominator of coil tolerance functional (tkruger)
subroutine rcfluxsens(ideriv)
  use globals, only: dp, zero, half, pi2, machprec, ncpu, myid, ounit, &
       coil, DoF, Ncoils, Nfixgeo, Ndof, LM_fvec, LM_fjac, Nfp, dpsidr, &
       ghost_call, ghost_use, ghost_once, weight_resbn, t1Z, pflsuc, MPI_COMM_FOCUS
  use mpi
  implicit none
  INTEGER, INTENT(in) :: ideriv
  INTEGER             :: astat, ierr, icoil, idof, ND, ivec, mult
  REAL                :: d1Z(1:Ndof)
  REAL                :: absdpsihold

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if ( weight_resbn .eq. 0.0 ) then
     if (ghost_use .eq. 1 .and. ghost_call .eq. 1) then
        call ghost(1, pflsuc)
        !FATAL( rcfluxsens , pflsuc .eq. 0, Periodic field line solve failed )
        if ( pflsuc .eq. 0 ) then
           dpsidr = 1000.0
           d1Z = zero
           return
        endif
        if ( ghost_once .eq. 1 ) then
           ghost_call = 0
        endif
     endif
  endif

  dpsidr = zero

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
           dpsidr = dpsidr + absdpsihold*real(mult)
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
              t1Z(idof+1:idof+ND) = t1Z(idof+1:idof+ND) + d1Z(idof+1:idof+ND)*real(mult)
           endif 
           idof = idof + ND
        endif
     enddo !end icoil;
     FATAL( rcfluxsens , idof .ne. Ndof, counting error in packing )
  endif
  return
end subroutine rcfluxsens

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine dpsi0(icoil, absdpsi)

  use globals, only: dp, zero, coil, gsurf, myid, ounit, Ncoils, resbn_m, pi2, Npert, dpsi_linear, MPI_COMM_FOCUS 
  implicit none
  include "mpif.h"

  INTEGER, intent(in)  :: icoil
  REAL   , intent(out) :: absdpsi
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  INTEGER              :: i, NS, astat, ierr, k
  REAL, dimension(1:coil(icoil)%NS) :: dAxx, dAxy, dAxz, dAyx, dAyy, dAyz, dAzx, dAzy, dAzz

  FATAL( dpsi0, icoil .lt. 1 .or. icoil .gt. Ncoils, icoil not in right range )
  
  coil(icoil)%dpsidx = zero
  coil(icoil)%dpsidy = zero
  coil(icoil)%dpsidz = zero
  NS = coil(icoil)%NS
  do i = 1, gsurf(1)%Nseg_stable-1
     call deltaafield(icoil, gsurf(1)%ox(i), gsurf(1)%oy(i), gsurf(1)%oz(i), dAxx, dAxy, dAxz, dAyx, dAyy, dAyz, dAzx, dAzy, dAzz)
     coil(icoil)%dpsidx(1:NS) = coil(icoil)%dpsidx(1:NS) + dAxx(1:NS)*gsurf(1)%oxdot(i) + dAyx(1:NS)*gsurf(1)%oydot(i) + dAzx(1:NS)*gsurf(1)%ozdot(i)
     coil(icoil)%dpsidy(1:NS) = coil(icoil)%dpsidy(1:NS) + dAxy(1:NS)*gsurf(1)%oxdot(i) + dAyy(1:NS)*gsurf(1)%oydot(i) + dAzy(1:NS)*gsurf(1)%ozdot(i)
     coil(icoil)%dpsidz(1:NS) = coil(icoil)%dpsidz(1:NS) + dAxz(1:NS)*gsurf(1)%oxdot(i) + dAyz(1:NS)*gsurf(1)%oydot(i) + dAzz(1:NS)*gsurf(1)%ozdot(i)
     call deltaafield(icoil, gsurf(1)%xx(i), gsurf(1)%xy(i), gsurf(1)%xz(i), dAxx, dAxy, dAxz, dAyx, dAyy, dAyz, dAzx, dAzy, dAzz)
     coil(icoil)%dpsidx(1:NS) = coil(icoil)%dpsidx(1:NS) - dAxx(1:NS)*gsurf(1)%xxdot(i) - dAyx(1:NS)*gsurf(1)%xydot(i) - dAzx(1:NS)*gsurf(1)%xzdot(i)
     coil(icoil)%dpsidy(1:NS) = coil(icoil)%dpsidy(1:NS) - dAxy(1:NS)*gsurf(1)%xxdot(i) - dAyy(1:NS)*gsurf(1)%xydot(i) - dAzy(1:NS)*gsurf(1)%xzdot(i)
     coil(icoil)%dpsidz(1:NS) = coil(icoil)%dpsidz(1:NS) - dAxz(1:NS)*gsurf(1)%xxdot(i) - dAyz(1:NS)*gsurf(1)%xydot(i) - dAzz(1:NS)*gsurf(1)%xzdot(i)
  enddo
  coil(icoil)%dpsidx(1:NS) = pi2*real(resbn_m)*coil(icoil)%dpsidx(1:NS)/real(gsurf(1)%Nseg_stable-1)
  coil(icoil)%dpsidy(1:NS) = pi2*real(resbn_m)*coil(icoil)%dpsidy(1:NS)/real(gsurf(1)%Nseg_stable-1)
  coil(icoil)%dpsidz(1:NS) = pi2*real(resbn_m)*coil(icoil)%dpsidz(1:NS)/real(gsurf(1)%Nseg_stable-1)
  if ( dpsi_linear .eq. 0 ) then
    absdpsi = 0.5*pi2*sum( coil(icoil)%dpsidx(1:NS)**2 + coil(icoil)%dpsidy(1:NS)**2 + coil(icoil)%dpsidz(1:NS)**2 )/real(NS)
  else
    absdpsi = 0.0
    ! Can be parallelized
    do k = 1, Npert
       absdpsi = absdpsi + ( pi2*sum( coil(icoil)%dpsidx(1:NS)*coil(icoil)%pertx(1:NS,k) + coil(icoil)%dpsidy(1:NS)*coil(icoil)%perty(1:NS,k) + &
            coil(icoil)%dpsidz(1:NS)*coil(icoil)%pertz(1:NS,k) ) / real(NS) )**2.0
    enddo
    absdpsi = absdpsi / real(Npert)
  endif
  
  return

end subroutine dpsi0

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine dpsi1(icoil, derivs, ND)

  use globals, only: dp, zero, coil, gsurf, myid, ounit, Ncoils, resbn_m, pi2, DoF, bsconstant, MPI_COMM_FOCUS
  implicit none
  include "mpif.h"

  INTEGER, intent(in)  :: icoil, ND
  REAL   , intent(out) :: derivs(1:1, 1:ND)
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  INTEGER              :: i, j, k, NS, astat, ierr
  REAL, dimension(1:1, 0:coil(icoil)%NS-1) :: dLLx, dLy, dLz
  REAL, dimension(0:coil(icoil)%NS) :: xx, yy, zz, xt, yt, zt, xa, ya, za, absdr
  REAL, dimension(0:coil(icoil)%NS,1:3) :: dpsidr, dpsidrp, dr, rp, term1, term2, term3, term4
  REAL, dimension(0:coil(icoil)%NS,1:3,1:3) :: drdpsidr, drpdpsidr, drpdpsidrp
  REAL, dimension(1:3,1:gsurf(1)%Nseg_stable-1) :: dlo, dlx

  FATAL( dpsi1, icoil .lt. 1 .or. icoil .gt. Ncoils, icoil not in right range )
  NS = coil(icoil)%NS
  
  dpsidr = zero
  dpsidrp = zero
  drdpsidr = zero
  drpdpsidr = zero
  drpdpsidrp = zero
  
  dpsidr(1:NS,1) = coil(icoil)%dpsidx(1:NS)
  dpsidr(1:NS,2) = coil(icoil)%dpsidy(1:NS)
  dpsidr(1:NS,3) = coil(icoil)%dpsidz(1:NS)
  dpsidr(0,1) = dpsidr(NS,1)
  dpsidr(0,2) = dpsidr(NS,2)
  dpsidr(0,3) = dpsidr(NS,3)

  dlo(1,1:gsurf(1)%Nseg_stable-1) = gsurf(1)%oxdot(1:gsurf(1)%Nseg_stable-1)
  dlo(2,1:gsurf(1)%Nseg_stable-1) = gsurf(1)%oydot(1:gsurf(1)%Nseg_stable-1)
  dlo(3,1:gsurf(1)%Nseg_stable-1) = gsurf(1)%ozdot(1:gsurf(1)%Nseg_stable-1)
  dlx(1,1:gsurf(1)%Nseg_stable-1) = gsurf(1)%xxdot(1:gsurf(1)%Nseg_stable-1)
  dlx(2,1:gsurf(1)%Nseg_stable-1) = gsurf(1)%xydot(1:gsurf(1)%Nseg_stable-1)
  dlx(3,1:gsurf(1)%Nseg_stable-1) = gsurf(1)%xzdot(1:gsurf(1)%Nseg_stable-1)
  
  xx(0:NS) = coil(icoil)%xx(0:NS)
  yy(0:NS) = coil(icoil)%yy(0:NS)
  zz(0:NS) = coil(icoil)%zz(0:NS)
  xt(0:NS) = coil(icoil)%xt(0:NS)
  yt(0:NS) = coil(icoil)%yt(0:NS)
  zt(0:NS) = coil(icoil)%zt(0:NS)
  xa(0:NS) = coil(icoil)%xa(0:NS)
  ya(0:NS) = coil(icoil)%ya(0:NS)
  za(0:NS) = coil(icoil)%za(0:NS)
  rp(0:NS,1) = xt(0:NS)
  rp(0:NS,2) = yt(0:NS)
  rp(0:NS,3) = zt(0:NS)

  do i = 1, gsurf(1)%Nseg_stable-1
     dr(0:NS,1) = gsurf(1)%ox(i) - xx(0:NS)
     dr(0:NS,2) = gsurf(1)%oy(i) - yy(0:NS)
     dr(0:NS,3) = gsurf(1)%oz(i) - zz(0:NS)
     absdr(0:NS) = sqrt( ( gsurf(1)%ox(i) - xx(0:NS) )**2.0 + ( gsurf(1)%oy(i) - yy(0:NS) )**2.0 + ( gsurf(1)%oz(i) - zz(0:NS) )**2.0 )
     do j = 1, 3
        dpsidrp(0:NS,j) = dpsidrp(0:NS,j) + 3.0*(xt(0:NS)*dr(0:NS,1)+yt(0:NS)*dr(0:NS,2)+zt(0:NS)*dr(0:NS,3))*( (xt(0:NS)*dlo(1,i)+yt(0:NS)*dlo(2,i)+zt(0:NS)*dlo(3,i))*dr(0:NS,j) + &
             (xt(0:NS)*dr(0:NS,1)+yt(0:NS)*dr(0:NS,2)+zt(0:NS)*dr(0:NS,3))*dlo(j,i) )/absdr(0:NS)**5.0
        dpsidrp(0:NS,j) = dpsidrp(0:NS,j) + (xa(0:NS)*dlo(1,i)+ya(0:NS)*dlo(2,i)+za(0:NS)*dlo(3,i))*dr(0:NS,j)/absdr(0:NS)**3.0
        dpsidrp(0:NS,j) = dpsidrp(0:NS,j) - (xt(0:NS)*xt(0:NS)+yt(0:NS)*yt(0:NS)+zt(0:NS)*zt(0:NS))*dlo(j,i)/absdr(0:NS)**3.0
        dpsidrp(0:NS,j) = dpsidrp(0:NS,j) - (xt(0:NS)*dlo(1,i)+yt(0:NS)*dlo(2,i)+zt(0:NS)*dlo(3,i))*rp(0:NS,j)/absdr(0:NS)**3.0
        dpsidrp(0:NS,j) = dpsidrp(0:NS,j) + (xa(0:NS)*dr(0:NS,1)+ya(0:NS)*dr(0:NS,2)+za(0:NS)*dr(0:NS,3))*dlo(j,i)/absdr(0:NS)**3.0
        do k = 1, 3
           drdpsidr(0:NS,j,k) = drdpsidr(0:NS,j,k) + ( 3.0*dr(0:NS,j)*( (xt(0:NS)*dlo(1,i)+yt(0:NS)*dlo(2,i)+zt(0:NS)*dlo(3,i))*dr(0:NS,k) + &
                (xt(0:NS)*dr(0:NS,1)+yt(0:NS)*dr(0:NS,2)+zt(0:NS)*dr(0:NS,3))*dlo(k,i) )/absdr(0:NS)**2.0 )/absdr(0:NS)**3.0
           if( j .eq. k ) drdpsidr(0:NS,j,k) = drdpsidr(0:NS,j,k) - (xt(0:NS)*dlo(1,i)+yt(0:NS)*dlo(2,i)+zt(0:NS)*dlo(3,i))/absdr(0:NS)**3.0
           drdpsidr(0:NS,j,k) = drdpsidr(0:NS,j,k) - rp(0:NS,j)*dlo(k,i)/absdr(0:NS)**3.0
           drpdpsidr(0:NS,j,k) = drpdpsidr(0:NS,j,k) + ( dlo(j,i)*dr(0:NS,k) + dr(0:NS,j)*dlo(k,i) )/absdr(0:NS)**3.0
           drpdpsidrp(0:NS,j,k) = drpdpsidrp(0:NS,j,k) + 3.0*(xt(0:NS)*dr(0:NS,1)+yt(0:NS)*dr(0:NS,2)+zt(0:NS)*dr(0:NS,3))*(dlo(j,i)*dr(0:NS,k)+dr(0:NS,j)*dlo(k,i))/absdr(0:NS)**5.0
           drpdpsidrp(0:NS,j,k) = drpdpsidrp(0:NS,j,k) - ( dlo(j,i)*rp(0:NS,k) + rp(0:NS,j)*dlo(k,i) )/absdr(0:NS)**3.0
        enddo
     enddo
     dr(0:NS,1) = gsurf(1)%xx(i) - xx(0:NS)
     dr(0:NS,2) = gsurf(1)%xy(i) - yy(0:NS)
     dr(0:NS,3) = gsurf(1)%xz(i) - zz(0:NS)
     absdr(0:NS) = sqrt( ( gsurf(1)%xx(i) - xx(0:NS) )**2.0 + ( gsurf(1)%xy(i) - yy(0:NS) )**2.0 + ( gsurf(1)%xz(i) - zz(0:NS) )**2.0 )
     do j = 1, 3
        dpsidrp(0:NS,j) = dpsidrp(0:NS,j) - 3.0*(xt(0:NS)*dr(0:NS,1)+yt(0:NS)*dr(0:NS,2)+zt(0:NS)*dr(0:NS,3))*( (xt(0:NS)*dlx(1,i)+yt(0:NS)*dlx(2,i)+zt(0:NS)*dlx(3,i))*dr(0:NS,j) + &
             (xt(0:NS)*dr(0:NS,1)+yt(0:NS)*dr(0:NS,2)+zt(0:NS)*dr(0:NS,3))*dlx(j,i) )/absdr(0:NS)**5.0
        dpsidrp(0:NS,j) = dpsidrp(0:NS,j) - (xa(0:NS)*dlx(1,i)+ya(0:NS)*dlx(2,i)+za(0:NS)*dlx(3,i))*dr(0:NS,j)/absdr(0:NS)**3.0
        dpsidrp(0:NS,j) = dpsidrp(0:NS,j) + (xt(0:NS)*xt(0:NS)+yt(0:NS)*yt(0:NS)+zt(0:NS)*zt(0:NS))*dlx(j,i)/absdr(0:NS)**3.0
        dpsidrp(0:NS,j) = dpsidrp(0:NS,j) + (xt(0:NS)*dlx(1,i)+yt(0:NS)*dlx(2,i)+zt(0:NS)*dlx(3,i))*rp(0:NS,j)/absdr(0:NS)**3.0
        dpsidrp(0:NS,j) = dpsidrp(0:NS,j) - (xa(0:NS)*dr(0:NS,1)+ya(0:NS)*dr(0:NS,2)+za(0:NS)*dr(0:NS,3))*dlx(j,i)/absdr(0:NS)**3.0
        do k = 1, 3
           drdpsidr(0:NS,j,k) = drdpsidr(0:NS,j,k) - ( 3.0*dr(0:NS,j)*( (xt(0:NS)*dlx(1,i)+yt(0:NS)*dlx(2,i)+zt(0:NS)*dlx(3,i))*dr(0:NS,k) + &
                (xt(0:NS)*dr(0:NS,1)+yt(0:NS)*dr(0:NS,2)+zt(0:NS)*dr(0:NS,3))*dlx(k,i) )/absdr(0:NS)**2.0 )/absdr(0:NS)**3.0
           if( j .eq. k ) drdpsidr(0:NS,j,k) = drdpsidr(0:NS,j,k) + (xt(0:NS)*dlx(1,i)+yt(0:NS)*dlx(2,i)+zt(0:NS)*dlx(3,i))/absdr(0:NS)**3.0
           drdpsidr(0:NS,j,k) = drdpsidr(0:NS,j,k) + rp(0:NS,j)*dlx(k,i)/absdr(0:NS)**3.0
           drpdpsidr(0:NS,j,k) = drpdpsidr(0:NS,j,k) - ( dlx(j,i)*dr(0:NS,k) + dr(0:NS,j)*dlx(k,i) )/absdr(0:NS)**3.0
           drpdpsidrp(0:NS,j,k) = drpdpsidrp(0:NS,j,k) - 3.0*(xt(0:NS)*dr(0:NS,1)+yt(0:NS)*dr(0:NS,2)+zt(0:NS)*dr(0:NS,3))*(dlx(j,i)*dr(0:NS,k)+dr(0:NS,j)*dlx(k,i))/absdr(0:NS)**5.0
           drpdpsidrp(0:NS,j,k) = drpdpsidrp(0:NS,j,k) + ( dlx(j,i)*rp(0:NS,k) + rp(0:NS,j)*dlx(k,i) )/absdr(0:NS)**3.0
        enddo
     enddo
  enddo
  dpsidrp(0:NS,1:3)        = coil(icoil)%I*bsconstant*pi2*real(resbn_m)*dpsidrp(   0:NS,1:3)     / real(gsurf(1)%Nseg_stable-1)
  drdpsidr(0:NS,1:3,1:3)   = coil(icoil)%I*bsconstant*pi2*real(resbn_m)*drdpsidr(  0:NS,1:3,1:3) / real(gsurf(1)%Nseg_stable-1)
  drpdpsidr(0:NS,1:3,1:3)  = coil(icoil)%I*bsconstant*pi2*real(resbn_m)*drpdpsidr( 0:NS,1:3,1:3) / real(gsurf(1)%Nseg_stable-1)
  drpdpsidrp(0:NS,1:3,1:3) = coil(icoil)%I*bsconstant*pi2*real(resbn_m)*drpdpsidrp(0:NS,1:3,1:3) / real(gsurf(1)%Nseg_stable-1)

  do j = 1, 3
     !term1(0:NS,j) = drdpsidr(0:NS,j,1)*dpsidr(0:NS,1) + drdpsidr(0:NS,j,2)*dpsidr(0:NS,2) + drdpsidr(0:NS,j,3)*dpsidr(0:NS,3)
     !term2(0:NS,j) = (dpsidr(0:NS,1)*dpsidrp(0:NS,1)+dpsidr(0:NS,2)*dpsidrp(0:NS,2)+dpsidr(0:NS,3)*dpsidrp(0:NS,3))* &
     !                (drpdpsidr(0:NS,j,1)*dpsidr(0:NS,1)+drpdpsidr(0:NS,j,2)*dpsidr(0:NS,2)+drpdpsidr(0:NS,j,3)*dpsidr(0:NS,3)) &
     !                / (dpsidr(0:NS,1)*dpsidr(0:NS,1)+dpsidr(0:NS,2)*dpsidr(0:NS,2)+dpsidr(0:NS,3)*dpsidr(0:NS,3))
     !term3(0:NS,j) = -1.0*( drpdpsidrp(0:NS,j,1)*dpsidr(0:NS,1) + drpdpsidrp(0:NS,j,2)*dpsidr(0:NS,2) + drpdpsidrp(0:NS,j,3)*dpsidr(0:NS,3) )
     !term4(0:NS,j) = -1.0*( drpdpsidr(0:NS,j,1)*dpsidrp(0:NS,1) + drpdpsidr(0:NS,j,2)*dpsidrp(0:NS,2) + drpdpsidr(0:NS,j,3)*dpsidrp(0:NS,3) )

     term1(0:NS,j) = drdpsidr(0:NS,j,1)*dpsidr(0:NS,1) + drdpsidr(0:NS,j,2)*dpsidr(0:NS,2) + drdpsidr(0:NS,j,3)*dpsidr(0:NS,3)
     term2(0:NS,j) = -1.0*( drpdpsidrp(0:NS,j,1)*dpsidr(0:NS,1) + drpdpsidrp(0:NS,j,2)*dpsidr(0:NS,2) + drpdpsidrp(0:NS,j,3)*dpsidr(0:NS,3) )
     term3(0:NS,j) = -1.0*( drpdpsidr(0:NS,j,1)*dpsidrp(0:NS,1) + drpdpsidr(0:NS,j,2)*dpsidrp(0:NS,2) + drpdpsidr(0:NS,j,3)*dpsidrp(0:NS,3) )
  enddo
  
!  do j = 1, 3
!     term1(0:NS,j) = term1(0:NS,j) / sqrt(dpsidr(0:NS,1)*dpsidr(0:NS,1)+dpsidr(0:NS,2)*dpsidr(0:NS,2)+dpsidr(0:NS,3)*dpsidr(0:NS,3))
!     term2(0:NS,j) = term2(0:NS,j) / sqrt(dpsidr(0:NS,1)*dpsidr(0:NS,1)+dpsidr(0:NS,2)*dpsidr(0:NS,2)+dpsidr(0:NS,3)*dpsidr(0:NS,3)) 
!     term3(0:NS,j) = term3(0:NS,j) / sqrt(dpsidr(0:NS,1)*dpsidr(0:NS,1)+dpsidr(0:NS,2)*dpsidr(0:NS,2)+dpsidr(0:NS,3)*dpsidr(0:NS,3))
!     term4(0:NS,j) = term4(0:NS,j) / sqrt(dpsidr(0:NS,1)*dpsidr(0:NS,1)+dpsidr(0:NS,2)*dpsidr(0:NS,2)+dpsidr(0:NS,3)*dpsidr(0:NS,3))
!  enddo
  
  dLLx(1,0:NS-1) = term1(0:NS-1,1) + term2(0:NS-1,1) + term3(0:NS-1,1)! + term4(0:NS-1,1)
  dLy(1,0:NS-1)  = term1(0:NS-1,2) + term2(0:NS-1,2) + term3(0:NS-1,2)! + term4(0:NS-1,2)
  dLz(1,0:NS-1)  = term1(0:NS-1,3) + term2(0:NS-1,3) + term3(0:NS-1,3)! + term4(0:NS-1,3)
  derivs(1:1, 1:ND) = matmul(dLLx, DoF(icoil)%xof) + matmul(dLy, DoF(icoil)%yof) + matmul(dLz, DoF(icoil)%zof)
  derivs(1:1, 1:ND) = pi2 * derivs(1:1, 1:ND) / real(NS)

  return

end subroutine dpsi1

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
