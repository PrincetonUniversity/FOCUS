!title (nissin) ! Calculate nissin complexity objective functon and its derivatives. (tkruger)

!latex \briefly{The objective function minimizes the coil's derivative of curvature and 
!latex         the coil's curvature*torsion. This objective function was designed to 
!latex         make coils easier to fabricate using the Nissin Freeform Bender. 
!latex         This function is still under development. 

!latex \calledby{\link{solvers}}

!latex  \section{General}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! nis is total penalty 
! chi = chi + weight_nis*nis
! t1N is total derivative of penalty
! LM implemented
! not parallelized, at some point check to see how long takes to run
subroutine nissin(ideriv)
  use globals, only: dp, zero, half, pi2, machprec, ncpu, myid, ounit, MPI_COMM_FOCUS, &
       coil, DoF, Ncoils, Nfixgeo, Ndof, nis, t1N, t2N, weight_nis, FouCoil, &
       mnis, inis, LM_fvec, LM_fjac, nis_alpha, case_nis, nis0

  implicit none
  include "mpif.h"
  INTEGER, INTENT(in) :: ideriv

  INTEGER             :: astat, ierr, icoil, idof, ND, NF, ivec
  REAL                :: nisAdd, nisHold

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  nis = zero
  nisAdd = zero
  nisHold = zero
  ivec = 1

  FATAL( nissin , case_nis .ne. 1 , invalid choice of case_nis )

  if( ideriv >= 0 ) then

     do icoil = 1, Ncoils
        if( coil(icoil)%type .ne. 1 ) exit ! only for Fourier
        if( coil(icoil)%Lc     /=  0 ) then ! if geometry is free
           call NisDeriv0(icoil,nisAdd)
           if( case_nis .eq. 1 ) then
              !FATAL( nissin, nis_alpha < 1 , nis_alpha >= 1 for case_nis == 1 ) ! Change 
              !nisHold = (nisAdd**2 + 1)**nis_alpha - 1
              !nis = nis + nisHold
              nis = nis + nisAdd
              if (mnis > 0) then ! L-M format of targets
                 !LM_fvec(inis+ivec) = weight_nis*nisHold
                 LM_fvec(inis+ivec) = weight_nis*nisAdd
                 ivec = ivec + 1
              endif
           else
              !FATAL( nissin, nis_alpha < 0 , nis_alpha >= 0 for case_nis == 2 ) ! Change
              !if( nisAdd**2 > nis0**2 ) then 
              !   nisHold = ( cosh(nis_alpha*(nisAdd**2 - nis0**2)) - 1 )**2
              !   nis = nis + nisHold
              !endif
              !if (mnis > 0) then ! L-M format of targets
              !   LM_fvec(inis+ivec) = weight_nis*nisHold
              !   ivec = ivec + 1
              !endif
           endif
        endif 
     enddo

     nis = nis / (Ncoils - Nfixgeo + machprec)

  endif

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if ( ideriv >= 1 ) then

     t1N = zero
     ivec = 1
     idof = 0
     do icoil = 1, Ncoils

        if(coil(icoil)%type .ne. 1) exit ! only for Fourier

        ND = DoF(icoil)%ND
        NF = FouCoil(icoil)%NF

        if ( coil(icoil)%Ic /= 0 ) then !if current is free;
           idof = idof + 1
        endif

        if ( coil(icoil)%Lc /= 0 ) then !if geometry is free;
           call NisDeriv1( icoil, t1N(idof+1:idof+ND), ND, NF )
           if (mnis > 0) then ! L-M format of targets
              LM_fjac(inis+ivec, idof+1:idof+ND) = weight_nis * t1N(idof+1:idof+ND)
              ivec = ivec + 1
           endif
           idof = idof + ND
        endif

     enddo
     FATAL( nissin , idof .ne. Ndof, counting error in packing )

     t1N = t1N / (Ncoils - Nfixgeo + machprec)

  endif

  return
end subroutine nissin

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! calculate Nissin Complexity

subroutine NisDeriv0(icoil,nisRet)

  use globals, only: dp, zero, pi2, ncpu, astat, ierr, myid, ounit, coil, NFcoil, Nseg, Ncoils, &
          nis_alpha, nis0, MPI_COMM_FOCUS

  implicit none
  include "mpif.h"

  INTEGER, intent(in)  :: icoil
  REAL   , intent(out) :: nisRet 
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  INTEGER              :: kseg, NS
  REAL                 :: nis_hold
  REAL,allocatable     :: niss(:), xt(:), yt(:), zt(:), xa(:), ya(:), za(:), xb(:), yb(:), zb(:), &
     lambdax(:), lambday(:), lambdaz(:), lambdanorm(:), lambdaunitx(:), lambdaunity(:), lambdaunitz(:), &
     lambdapx(:), lambdapy(:), lambdapz(:), S(:), S1(:), S2(:), absrp(:)

  NS = coil(icoil)%NS 
  SALLOCATE(niss, (0:NS), zero)
  SALLOCATE(S, (0:NS), zero)
  SALLOCATE(S1, (0:NS), zero)
  SALLOCATE(S2, (0:NS), zero)
  SALLOCATE(lambdax, (0:NS), zero)
  SALLOCATE(lambday, (0:NS), zero)
  SALLOCATE(lambdaz, (0:NS), zero)
  SALLOCATE(lambdapx, (0:NS), zero)
  SALLOCATE(lambdapy, (0:NS), zero)
  SALLOCATE(lambdapz, (0:NS), zero)
  SALLOCATE(lambdanorm, (0:NS), zero)
  SALLOCATE(lambdaunitx, (0:NS), zero)
  SALLOCATE(lambdaunity, (0:NS), zero)
  SALLOCATE(lambdaunitz, (0:NS), zero)
  SALLOCATE(absrp, (0:NS), zero )
  SALLOCATE(xt, (0:NS), zero)
  SALLOCATE(yt, (0:NS), zero)
  SALLOCATE(zt, (0:NS), zero)
  SALLOCATE(xa, (0:NS), zero)
  SALLOCATE(ya, (0:NS), zero)
  SALLOCATE(za, (0:NS), zero)
  SALLOCATE(xb, (0:NS), zero)
  SALLOCATE(yb, (0:NS), zero)
  SALLOCATE(zb, (0:NS), zero)
  xt(0:coil(icoil)%NS) = coil(icoil)%xt(0:coil(icoil)%NS)
  yt(0:coil(icoil)%NS) = coil(icoil)%yt(0:coil(icoil)%NS)
  zt(0:coil(icoil)%NS) = coil(icoil)%zt(0:coil(icoil)%NS)
  xa(0:coil(icoil)%NS) = coil(icoil)%xa(0:coil(icoil)%NS)
  ya(0:coil(icoil)%NS) = coil(icoil)%ya(0:coil(icoil)%NS)
  za(0:coil(icoil)%NS) = coil(icoil)%za(0:coil(icoil)%NS)
  xb(0:coil(icoil)%NS) = coil(icoil)%xb(0:coil(icoil)%NS)
  yb(0:coil(icoil)%NS) = coil(icoil)%yb(0:coil(icoil)%NS)
  zb(0:coil(icoil)%NS) = coil(icoil)%zb(0:coil(icoil)%NS)

  FATAL( NisDeriv0, icoil .lt. 1 .or. icoil .gt. Ncoils, icoil not in right range )

  niss = zero
  nisRet = zero

  lambdax(0:NS) = yt(0:NS)*za(0:NS) - zt(0:NS)*ya(0:NS)
  lambday(0:NS) = zt(0:NS)*xa(0:NS) - xt(0:NS)*za(0:NS)
  lambdaz(0:NS) = xt(0:NS)*ya(0:NS) - yt(0:NS)*xa(0:NS)
  lambdapx(0:NS) = yt(0:NS)*zb(0:NS) - zt(0:NS)*yb(0:NS)
  lambdapy(0:NS) = zt(0:NS)*xb(0:NS) - xt(0:NS)*zb(0:NS)
  lambdapz(0:NS) = xt(0:NS)*yb(0:NS) - yt(0:NS)*xb(0:NS)
  lambdanorm(0:NS) = sqrt( lambdax(0:NS)**2 + lambday(0:NS)**2 + lambdaz(0:NS)**2 ) 
  lambdaunitx(0:NS) = lambdax(0:NS) / lambdanorm(0:NS)
  lambdaunity(0:NS) = lambday(0:NS) / lambdanorm(0:NS)
  lambdaunitz(0:NS) = lambdaz(0:NS) / lambdanorm(0:NS)
  absrp(0:NS) = sqrt( xt(0:NS)**2 + yt(0:NS)**2 + zt(0:NS)**2 )

  S1(0:NS) = (lambdaunitx(0:NS)*lambdapx(0:NS) + lambdaunity(0:NS)*lambdapy(0:NS) + lambdaunitz(0:NS)*lambdapz(0:NS)) / absrp(0:NS)**3
  S1(0:NS) = S1(0:NS) - 3*lambdanorm(0:NS)*(xt(0:NS)*xa(0:NS) + yt(0:NS)*ya(0:NS) + zt(0:NS)*za(0:NS)) / absrp(0:NS)**5
  S1(0:NS) = S1(0:NS) / absrp(0:NS)

  S2(0:NS) = (lambdaunitx(0:NS)*xb(0:NS) + lambdaunity(0:NS)*yb(0:NS) + lambdaunitz(0:NS)*zb(0:NS)) / absrp(0:NS)**3

  S(0:NS) = sqrt(S1(0:NS)**2 + S2(0:NS)**2)

  coil(icoil)%maxs = maxval(S)

  do kseg = 0, NS
     if( S(kseg) .ge. nis0) then
        !niss(kseg) = (cosh(nis_alpha*(S(kseg)-nis0)) - 1)**2
        niss(kseg) = (nis_alpha*(S(kseg)-nis0))**2
     endif
  enddo
  
  !coil(icoil)%maxcurv = maxval(curvv)

  niss(0:NS) = niss(0:NS)*absrp(0:NS)

  nisRet = sum(niss)-niss(0)
  call lenDeriv0( icoil, coil(icoil)%L )
  nisRet = pi2*nisRet/(NS*coil(icoil)%L)

  DALLOCATE(niss)
  DALLOCATE(S)
  DALLOCATE(S1)
  DALLOCATE(S2)
  DALLOCATE(lambdax)
  DALLOCATE(lambday)
  DALLOCATE(lambdaz)
  DALLOCATE(lambdapx)
  DALLOCATE(lambdapy)
  DALLOCATE(lambdapz)
  DALLOCATE(lambdanorm)
  DALLOCATE(lambdaunitx)
  DALLOCATE(lambdaunity)
  DALLOCATE(lambdaunitz)
  DALLOCATE(absrp)
  DALLOCATE(xt)
  DALLOCATE(yt)
  DALLOCATE(zt)
  DALLOCATE(xa)
  DALLOCATE(ya)
  DALLOCATE(za)
  DALLOCATE(xb)
  DALLOCATE(yb)
  DALLOCATE(zb)

  return

end subroutine NisDeriv0

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine NisDeriv1(icoil, derivs, ND, NF) !Calculate all derivatives for a coil

  use globals, only: dp, zero, pi2, coil, DoF, myid, ounit, Ncoils, &
          case_nis, nis_alpha, nis0, FouCoil, MPI_COMM_FOCUS
  implicit none
  include "mpif.h"

  INTEGER, intent(in)  :: icoil, ND , NF
  REAL   , intent(out) :: derivs(1:1, 1:ND)
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  INTEGER              :: kseg, astat, ierr, doff, n, i, NS
  REAL                 :: L, nisRet
  REAL                 :: d1L(1:1, 1:ND)
  REAL,allocatable     :: absrp(:), xt(:), yt(:), zt(:), xa(:), ya(:), za(:), xb(:), yb(:), zb(:), &
     lambdax(:), lambday(:), lambdaz(:), lambdanorm(:), lambdaunitx(:), lambdaunity(:), lambdaunitz(:), &
     dxtdDoF(:,:), dytdDoF(:,:), dztdDoF(:,:), dxadDoF(:,:), &
     dyadDoF(:,:), dzadDoF(:,:), dxbdDoF(:,:), dybdDoF(:,:), dzbdDoF(:,:), dlambdaxdDof(:,:), &
     dlambdaydDof(:,:), dlambdazdDof(:,:), dtaudDof(:,:), dabsrpdDof(:,:), lambdapx(:), &
     lambdapy(:), lambdapz(:), S(:), S1(:), S2(:), dSdDof(:,:), dS1dDof(:,:), dS2dDof(:,:), dnissdDof(:,:), &
     dlambdapxdDof(:,:), dlambdapydDof(:,:), dlambdapzdDof(:,:), dlambdaunitxdDof(:,:), dlambdaunitydDof(:,:), &
     dlambdaunitzdDof(:,:), niss(:)

  FATAL( NisDeriv1, icoil .lt. 1 .or. icoil .gt. Ncoils, icoil not in right range )
  FATAL( NisDeriv1, ND .ne. 3+6*NF, Degrees of freedom are incorrect )

  NS = coil(icoil)%NS

  SALLOCATE(niss, (0:NS), zero)
  SALLOCATE(S, (0:NS), zero)
  SALLOCATE(S1, (0:NS), zero)
  SALLOCATE(S2, (0:NS), zero)
  SALLOCATE(absrp, (0:NS), zero)
  SALLOCATE(lambdax, (0:NS), zero)
  SALLOCATE(lambday, (0:NS), zero)
  SALLOCATE(lambdaz, (0:NS), zero)
  SALLOCATE(lambdapx, (0:NS), zero)
  SALLOCATE(lambdapy, (0:NS), zero)
  SALLOCATE(lambdapz, (0:NS), zero)
  SALLOCATE(lambdanorm, (0:NS), zero)
  SALLOCATE(lambdaunitx, (0:NS), zero)
  SALLOCATE(lambdaunity, (0:NS), zero)
  SALLOCATE(lambdaunitz, (0:NS), zero)
  SALLOCATE(xt, (0:NS), zero)
  SALLOCATE(yt, (0:NS), zero)
  SALLOCATE(zt, (0:NS), zero)
  SALLOCATE(xa, (0:NS), zero)
  SALLOCATE(ya, (0:NS), zero)
  SALLOCATE(za, (0:NS), zero)
  SALLOCATE(xb, (0:NS), zero)
  SALLOCATE(yb, (0:NS), zero)
  SALLOCATE(zb, (0:NS), zero)

  SALLOCATE(dxtdDoF, (0:NS,1:ND), zero)
  SALLOCATE(dytdDoF, (0:NS,1:ND), zero)
  SALLOCATE(dztdDoF, (0:NS,1:ND), zero)
  SALLOCATE(dxadDoF, (0:NS,1:ND), zero)
  SALLOCATE(dyadDoF, (0:NS,1:ND), zero)
  SALLOCATE(dzadDoF, (0:NS,1:ND), zero)
  SALLOCATE(dxbdDoF, (0:NS,1:ND), zero)
  SALLOCATE(dybdDoF, (0:NS,1:ND), zero)
  SALLOCATE(dzbdDoF, (0:NS,1:ND), zero)
  SALLOCATE(dlambdaxdDof, (0:NS,1:ND), zero)
  SALLOCATE(dlambdaydDof, (0:NS,1:ND), zero)
  SALLOCATE(dlambdazdDof, (0:NS,1:ND), zero)
  SALLOCATE(dtaudDof, (0:NS,1:ND), zero)
  SALLOCATE(dabsrpdDof, (0:NS,1:ND), zero)
  SALLOCATE(dSdDof, (0:NS,1:ND), zero)
  SALLOCATE(dS1dDof, (0:NS,1:ND), zero)
  SALLOCATE(dS2dDof, (0:NS,1:ND), zero)
  SALLOCATE(dnissdDof, (0:NS,1:ND), zero)
  SALLOCATE(dlambdapxdDof, (0:NS,1:ND), zero)
  SALLOCATE(dlambdapydDof, (0:NS,1:ND), zero)
  SALLOCATE(dlambdapzdDof, (0:NS,1:ND), zero)
  SALLOCATE(dlambdaunitxdDof, (0:NS,1:ND), zero)
  SALLOCATE(dlambdaunitydDof, (0:NS,1:ND), zero)
  SALLOCATE(dlambdaunitzdDof, (0:NS,1:ND), zero)

  xt(0:coil(icoil)%NS) = coil(icoil)%xt(0:coil(icoil)%NS)
  yt(0:coil(icoil)%NS) = coil(icoil)%yt(0:coil(icoil)%NS)
  zt(0:coil(icoil)%NS) = coil(icoil)%zt(0:coil(icoil)%NS)
  xa(0:coil(icoil)%NS) = coil(icoil)%xa(0:coil(icoil)%NS)
  ya(0:coil(icoil)%NS) = coil(icoil)%ya(0:coil(icoil)%NS)
  za(0:coil(icoil)%NS) = coil(icoil)%za(0:coil(icoil)%NS)
  xb(0:coil(icoil)%NS) = coil(icoil)%xb(0:coil(icoil)%NS)
  yb(0:coil(icoil)%NS) = coil(icoil)%yb(0:coil(icoil)%NS)
  zb(0:coil(icoil)%NS) = coil(icoil)%zb(0:coil(icoil)%NS)
  
  do n = 1,NF
     dxtdDof(0:NS,n+1)      = -1*FouCoil(icoil)%smt(0:NS,n) * n
     dxtdDof(0:NS,n+NF+1)   =    FouCoil(icoil)%cmt(0:NS,n) * n
     dytdDof(0:NS,n+2*NF+2) = -1*FouCoil(icoil)%smt(0:NS,n) * n
     dytdDof(0:NS,n+3*NF+2) =    FouCoil(icoil)%cmt(0:NS,n) * n
     dztdDof(0:NS,n+4*NF+3) = -1*FouCoil(icoil)%smt(0:NS,n) * n
     dztdDof(0:NS,n+5*NF+3) =    FouCoil(icoil)%cmt(0:NS,n) * n

     dxadDof(0:NS,n+1)      = -1*FouCoil(icoil)%cmt(0:NS,n) * n*n
     dxadDof(0:NS,n+NF+1)   = -1*FouCoil(icoil)%smt(0:NS,n) * n*n
     dyadDof(0:NS,n+2*NF+2) = -1*FouCoil(icoil)%cmt(0:NS,n) * n*n
     dyadDof(0:NS,n+3*NF+2) = -1*FouCoil(icoil)%smt(0:NS,n) * n*n
     dzadDof(0:NS,n+4*NF+3) = -1*FouCoil(icoil)%cmt(0:NS,n) * n*n
     dzadDof(0:NS,n+5*NF+3) = -1*FouCoil(icoil)%smt(0:NS,n) * n*n

     dxbdDof(0:NS,n+1)      =    FouCoil(icoil)%smt(0:NS,n) * n*n*n
     dxbdDof(0:NS,n+NF+1)   = -1*FouCoil(icoil)%cmt(0:NS,n) * n*n*n
     dybdDof(0:NS,n+2*NF+2) =    FouCoil(icoil)%smt(0:NS,n) * n*n*n
     dybdDof(0:NS,n+3*NF+2) = -1*FouCoil(icoil)%cmt(0:NS,n) * n*n*n
     dzbdDof(0:NS,n+4*NF+3) =    FouCoil(icoil)%smt(0:NS,n) * n*n*n
     dzbdDof(0:NS,n+5*NF+3) = -1*FouCoil(icoil)%cmt(0:NS,n) * n*n*n
  enddo

  derivs = zero
  d1L = zero
  call lenDeriv0( icoil, coil(icoil)%L )
  L = coil(icoil)%L
  call lenDeriv1( icoil, d1L(1:1,1:ND), ND )

  lambdax(0:NS) = yt(0:NS)*za(0:NS) - zt(0:NS)*ya(0:NS) 
  lambday(0:NS) = zt(0:NS)*xa(0:NS) - xt(0:NS)*za(0:NS)
  lambdaz(0:NS) = xt(0:NS)*ya(0:NS) - yt(0:NS)*xa(0:NS)
  lambdapx(0:NS) = yt(0:NS)*zb(0:NS) - zt(0:NS)*yb(0:NS)
  lambdapy(0:NS) = zt(0:NS)*xb(0:NS) - xt(0:NS)*zb(0:NS)
  lambdapz(0:NS) = xt(0:NS)*yb(0:NS) - yt(0:NS)*xb(0:NS)
  lambdanorm(0:NS) = sqrt( lambdax(0:NS)**2 + lambday(0:NS)**2 + lambdaz(0:NS)**2 )
  lambdaunitx(0:NS) = lambdax(0:NS) / lambdanorm(0:NS)
  lambdaunity(0:NS) = lambday(0:NS) / lambdanorm(0:NS)
  lambdaunitz(0:NS) = lambdaz(0:NS) / lambdanorm(0:NS)

  absrp(0:NS) = sqrt( xt(0:NS)**2 + yt(0:NS)**2 + zt(0:NS)**2 )

  do i = 1,ND
     dlambdaxdDof(0:NS,i) = dytdDof(0:NS,i)*za(0:NS) - dztdDof(0:NS,i)*ya(0:NS) &
                             + yt(0:NS)*dzadDof(0:NS,i) - zt(0:NS)*dyadDof(0:NS,i)
     dlambdaydDof(0:NS,i) = dztdDof(0:NS,i)*xa(0:NS) - dxtdDof(0:NS,i)*za(0:NS) &
                             + zt(0:NS)*dxadDof(0:NS,i) - xt(0:NS)*dzadDof(0:NS,i)
     dlambdazdDof(0:NS,i) = dxtdDof(0:NS,i)*ya(0:NS) - dytdDof(0:NS,i)*xa(0:NS) &
                             + xt(0:NS)*dyadDof(0:NS,i) - yt(0:NS)*dxadDof(0:NS,i)

     dlambdapxdDof(0:NS,i) = dytdDof(0:NS,i)*zb(0:NS) - dztdDof(0:NS,i)*yb(0:NS) &
                             + yt(0:NS)*dzbdDof(0:NS,i) - zt(0:NS)*dybdDof(0:NS,i)
     dlambdapydDof(0:NS,i) = dztdDof(0:NS,i)*xb(0:NS) - dxtdDof(0:NS,i)*zb(0:NS) &
                             + zt(0:NS)*dxbdDof(0:NS,i) - xt(0:NS)*dzbdDof(0:NS,i)
     dlambdapzdDof(0:NS,i) = dxtdDof(0:NS,i)*yb(0:NS) - dytdDof(0:NS,i)*xb(0:NS) &
                             + xt(0:NS)*dybdDof(0:NS,i) - yt(0:NS)*dxbdDof(0:NS,i)

     dabsrpdDof(0:NS,i) = ( xt(0:NS)**2 + yt(0:NS)**2 + zt(0:NS)**2 )**(-.5) * ( xt(0:NS)*dxtdDof(0:NS,i) + & 
                          yt(0:NS)*dytdDof(0:NS,i) + zt(0:NS)*dztdDof(0:NS,i) )
     
     dlambdaunitxdDof(0:NS,i) = dlambdaxdDof(0:NS,i) / lambdanorm(0:NS) - lambdax(0:NS) * lambdanorm(0:NS)**(-3) * &
             (lambdax(0:NS)*dlambdaxdDof(0:NS,i) + lambday(0:NS)*dlambdaydDof(0:NS,i) + lambdaz(0:NS)*dlambdazdDof(0:NS,i))
     dlambdaunitydDof(0:NS,i) = dlambdaydDof(0:NS,i) / lambdanorm(0:NS) - lambday(0:NS) * lambdanorm(0:NS)**(-3) * &
             (lambdax(0:NS)*dlambdaxdDof(0:NS,i) + lambday(0:NS)*dlambdaydDof(0:NS,i) + lambdaz(0:NS)*dlambdazdDof(0:NS,i))
     dlambdaunitzdDof(0:NS,i) = dlambdazdDof(0:NS,i) / lambdanorm(0:NS) - lambdaz(0:NS) * lambdanorm(0:NS)**(-3) * &
             (lambdax(0:NS)*dlambdaxdDof(0:NS,i) + lambday(0:NS)*dlambdaydDof(0:NS,i) + lambdaz(0:NS)*dlambdazdDof(0:NS,i))
     
  enddo

  S1(0:NS) = (lambdaunitx(0:NS)*lambdapx(0:NS) + lambdaunity(0:NS)*lambdapy(0:NS) + lambdaunitz(0:NS)*lambdapz(0:NS)) / absrp(0:NS)**3
  S1(0:NS) = S1(0:NS) - 3*lambdanorm(0:NS)*(xt(0:NS)*xa(0:NS) + yt(0:NS)*ya(0:NS) + zt(0:NS)*za(0:NS)) / absrp(0:NS)**5
  S1(0:NS) = S1(0:NS) / absrp(0:NS)

  S2(0:NS) = (lambdaunitx(0:NS)*xb(0:NS) + lambdaunity(0:NS)*yb(0:NS) + lambdaunitz(0:NS)*zb(0:NS)) / absrp(0:NS)**3

  S(0:NS) = sqrt(S1(0:NS)**2 + S2(0:NS)**2)

  do kseg = 0, NS
     if( S(kseg) .ge. nis0) then
        !niss(kseg) = (cosh(nis_alpha*(S(kseg)-nis0)) - 1)**2
        niss(kseg) = (nis_alpha*(S(kseg)-nis0))**2
     endif
  enddo

  niss(0:NS) = niss(0:NS)*absrp(0:NS)

  nisRet = sum(niss)-niss(0)
  nisRet = pi2*nisRet/(NS*coil(icoil)%L)
  
  do i = 1,ND
     
     dS1dDof(0:NS,i) = ( dlambdaunitxdDof(0:NS,i)*lambdapx(0:NS) + dlambdaunitydDof(0:NS,i)*lambdapy(0:NS) + &
             dlambdaunitzdDof(0:NS,i)*lambdapz(0:NS) )*absrp(0:NS)**(-3)
     dS1dDof(0:NS,i) = dS1dDof(0:NS,i) + (lambdaunitx(0:NS)*dlambdapxdDof(0:NS,i) + lambdaunity(0:NS)*dlambdapydDof(0:NS,i) + &
             lambdaunitz(0:NS)*dlambdapzdDof(0:NS,i)) * absrp(0:NS)**(-3)
     dS1dDof(0:NS,i) = dS1dDof(0:NS,i) - 3*(lambdaunitx(0:NS)*lambdapx(0:NS) + lambdaunity(0:NS)*lambdapy(0:NS) + &
             lambdaunitz(0:NS)*lambdapz(0:NS) )*absrp(0:NS)**(-5)*(xt(0:NS)*dxtdDof(0:NS,i)+yt(0:NS)*dytdDof(0:NS,i)+zt(0:NS)*dztdDof(0:NS,i))
     dS1dDof(0:NS,i) = dS1dDof(0:NS,i) - 3*(lambdaunitx(0:NS)*dlambdaxdDof(0:NS,i)+lambdaunity(0:NS)*dlambdaydDof(0:NS,i)+lambdaunitz(0:NS)*dlambdazdDof(0:NS,i)) &
             * absrp(0:NS)**(-5)*(xt(0:NS)*xa(0:NS)+yt(0:NS)*ya(0:NS)+zt(0:NS)*za(0:NS))
     dS1dDof(0:NS,i) = dS1dDof(0:NS,i) + 15*lambdanorm(0:NS)*absrp(0:NS)**(-7)*(xt(0:NS)*dxtdDof(0:NS,i) + yt(0:NS)*dytdDof(0:NS,i) + &
             zt(0:NS)*dztdDof(0:NS,i))*(xt(0:NS)*xa(0:NS)+yt(0:NS)*ya(0:NS)+zt(0:NS)*za(0:NS))
     dS1dDof(0:NS,i) = dS1dDof(0:NS,i) - 3*lambdanorm(0:NS)*absrp(0:NS)**(-5)*(dxtdDof(0:NS,i)*xa(0:NS) + dytdDof(0:NS,i)*ya(0:NS) + &
             dztdDof(0:NS,i)*za(0:NS) + xt(0:NS)*dxadDof(0:NS,i) + yt(0:NS)*dyadDof(0:NS,i) + zt(0:NS)*dzadDof(0:NS,i))
     
     dS1dDof(0:NS,i) = dS1dDof(0:NS,i) / absrp(0:NS)
     !dS1dDof(0:NS,i) = dS1dDof(0:NS,i) - S1(0:NS) * absrp(0:NS)**(-3) * (xt(0:NS)*dxtdDof(0:NS,i) + yt(0:NS)*dytdDof(0:NS,i) + zt(0:NS)*dztdDof(0:NS,i))
     dS1dDof(0:NS,i) = dS1dDof(0:NS,i) - S1(0:NS) * absrp(0:NS)**(-2) * (xt(0:NS)*dxtdDof(0:NS,i) + yt(0:NS)*dytdDof(0:NS,i) + zt(0:NS)*dztdDof(0:NS,i))
     
     dS2dDof(0:NS,i) = (dlambdaunitxdDof(0:NS,i)*xb(0:NS) + dlambdaunitydDof(0:NS,i)*yb(0:NS) + dlambdaunitzdDof(0:NS,i)*zb(0:NS)) &
             * absrp(0:NS)**(-3)
     dS2dDof(0:NS,i) = dS2dDof(0:NS,i) + (lambdaunitx(0:NS)*dxbdDof(0:NS,i) + lambdaunity(0:NS)*dybdDof(0:NS,i) + &
             lambdaunitz(0:NS)*dzbdDof(0:NS,i) ) * absrp(0:NS)**(-3)
     dS2dDof(0:NS,i) = dS2dDof(0:NS,i) - 3*(lambdaunitx(0:NS)*xb(0:NS) + lambdaunity(0:NS)*yb(0:NS) + lambdaunitz(0:NS)*zb(0:NS) ) &
              * absrp(0:NS)**(-5) * (xt(0:NS)*dxtdDof(0:NS,i) + yt(0:NS)*dytdDof(0:NS,i) + zt(0:NS)*dztdDof(0:NS,i))
  
     dSdDof(0:NS,i) = ( S1(0:NS)**2 + S2(0:NS)**2 )**(-.5) * ( S1(0:NS)*dS1dDof(0:NS,i) + S2(0:NS)*dS2dDof(0:NS,i) )

     do kseg = 0, NS
        if( S(kseg) .ge. nis0 ) then
           !dnissdDof(kseg,i) = 2*((cosh(nis_alpha*(S(kseg)-nis0))-1))*sinh(nis_alpha*(S(kseg)-nis0))*nis_alpha*dSdDof(kseg,i)*absrp(kseg)
           !dnissdDof(kseg,i) = dnissdDof(kseg,i) + (cosh(nis_alpha*(S(kseg)-nis0)) - 1)**2 * dabsrpdDof(kseg,i)

           dnissdDof(kseg,i) = 2*nis_alpha*(nis_alpha*(S(kseg)-nis0))*dSdDof(kseg,i)*absrp(kseg)
           dnissdDof(kseg,i) = dnissdDof(kseg,i) + (nis_alpha*(S(kseg)-nis0))**2 * dabsrpdDof(kseg,i)
           derivs(1,i) = derivs(1,i) + dnissdDof(kseg,i)
        endif
     enddo
     derivs(1,i) = derivs(1,i) - dnissdDof(0,i)
  enddo

  derivs(1,1:ND) = derivs(1,1:ND)*pi2/(NS*coil(icoil)%L)
  derivs(1,1:ND) = derivs(1,1:ND) - d1L(1,1:ND)*nisRet / (coil(icoil)%L)

  DALLOCATE(niss)  
  DALLOCATE(S)
  DALLOCATE(S1)
  DALLOCATE(S2)
  DALLOCATE(absrp)
  DALLOCATE(lambdax)
  DALLOCATE(lambday)
  DALLOCATE(lambdaz)
  DALLOCATE(lambdapx)
  DALLOCATE(lambdapy)
  DALLOCATE(lambdapz)
  DALLOCATE(lambdanorm)
  DALLOCATE(lambdaunitx)
  DALLOCATE(lambdaunity)
  DALLOCATE(lambdaunitz)
  DALLOCATE(xt)
  DALLOCATE(yt)
  DALLOCATE(zt)
  DALLOCATE(xa)
  DALLOCATE(ya)
  DALLOCATE(za)
  DALLOCATE(xb)
  DALLOCATE(yb)
  DALLOCATE(zb)

  DALLOCATE(dxtdDof)
  DALLOCATE(dytdDof)
  DALLOCATE(dztdDof)
  DALLOCATE(dxadDof)
  DALLOCATE(dyadDof)
  DALLOCATE(dzadDof)
  DALLOCATE(dxbdDof)
  DALLOCATE(dybdDof)
  DALLOCATE(dzbdDof)
  DALLOCATE(dlambdaxdDof)
  DALLOCATE(dlambdaydDof)
  DALLOCATE(dlambdazdDof)
  DALLOCATE(dtaudDof)
  DALLOCATE(dabsrpdDof)
  DALLOCATE(dSdDof)
  DALLOCATE(dS1dDof)
  DALLOCATE(dS2dDof)
  DALLOCATE(dnissdDof)
  DALLOCATE(dlambdapxdDof)
  DALLOCATE(dlambdapydDof)
  DALLOCATE(dlambdapzdDof)
  DALLOCATE(dlambdaunitxdDof)
  DALLOCATE(dlambdaunitydDof)
  DALLOCATE(dlambdaunitzdDof)

  return

end subroutine NisDeriv1


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
