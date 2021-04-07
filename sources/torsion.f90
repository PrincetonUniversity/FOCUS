!title (torsion) ! Calculate torsion objective functon and its derivatives. (tkruger)

!latex \briefly{The objective function optimizes coils to have zero average torsion
!latex         which makes the multi-filaments models less complex.
!latex         This function is still under development. emph{targt\_length}.}

!latex \calledby{\link{solvers}}

!latex  \section{General}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! tors is total penalty 
! chi = chi + weight_tors*tors
! t1T is total derivative of penalty
! LM implemented
! not parallelized, does not take long to run
subroutine torsion(ideriv)
  use globals, only: dp, zero, half, pi2, machprec, ncpu, myid, ounit, MPI_COMM_FOCUS, &
       coil, DoF, Ncoils, Nfixgeo, Ndof, tors, t1T, t2T, weight_tors, FouCoil, &
       mtors, itors, LM_fvec, LM_fjac, tors_alpha, tors_beta, tors_gamma, case_tors, &
       penfun_tors, tors0

  implicit none
  include "mpif.h"
  INTEGER, INTENT(in) :: ideriv

  INTEGER             :: astat, ierr, icoil, idof, ND, NF, ivec
  REAL                :: torsAdd, torsHold

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  tors = zero
  torsAdd = zero
  torsHold = zero
  ivec = 1

  FATAL( torsion , case_tors .ne. 1 .and. case_tors .ne. 2 , invalid choice of case_tors )
  ! Put in other statements 

  if( ideriv >= 0 ) then

     do icoil = 1, Ncoils
        if( coil(icoil)%type .ne. 1 ) exit ! only for Fourier
        if( coil(icoil)%Lc /= 0 ) then ! if geometry is free
           call TorsDeriv0(icoil,torsAdd)
           if( case_tors .eq. 1 ) then
              if( torsAdd**2 .ge. tors0**2 ) then
                 if( penfun_tors .eq. 1 ) then
                    torsHold = ( cosh(tors_alpha*(torsAdd**2 - tors0**2)) - 1 )**2
                 else
                    torsHold = (tors_alpha*(torsAdd**2-tors0**2))**tors_beta
                 endif
              endif
           else
              torsHold = (torsAdd**2 + 1)**tors_gamma - 1
           endif
           tors = tors + torsHold
           if (mtors > 0) then ! L-M format of targets 
              LM_fvec(itors+ivec) = weight_tors*torsHold
              ivec = ivec + 1
           endif
        endif 
        torsHold = 0
     enddo

     tors = tors / (Ncoils - Nfixgeo + machprec)

  endif

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if ( ideriv >= 1 ) then

     t1T = zero
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
           call TorsDeriv1( icoil, t1T(idof+1:idof+ND), ND, NF )
           if (mtors > 0) then ! L-M format of targets
              LM_fjac(itors+ivec, idof+1:idof+ND) = weight_tors * t1T(idof+1:idof+ND)
              ivec = ivec + 1
           endif
           idof = idof + ND
        endif

     enddo
     FATAL( torsion , idof .ne. Ndof, counting error in packing )

     t1T = t1T / (Ncoils - Nfixgeo + machprec)

  endif

  return
end subroutine torsion

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! calculate average coil torsion

subroutine TorsDeriv0(icoil,torsRet)

  use globals, only: dp, zero, pi2, ncpu, astat, ierr, myid, ounit, coil, NFcoil, Nseg, Ncoils, &
          case_tau, MPI_COMM_FOCUS

  implicit none
  include "mpif.h"

  INTEGER, intent(in)  :: icoil
  REAL   , intent(out) :: torsRet
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  INTEGER              :: kseg, NS
  REAL                 :: tors_hold
  REAL,allocatable     :: torss(:), xt(:), yt(:), zt(:), xa(:), ya(:), za(:), xb(:), yb(:), zb(:), &
     lambdax(:), lambday(:), lambdaz(:), lambdanorm(:), lambdaunitx(:), lambdaunity(:), lambdaunitz(:), &
     absrp(:)

  NS = coil(icoil)%NS 
  SALLOCATE(torss, (0:NS), zero)
  SALLOCATE(lambdax, (0:NS), zero)
  SALLOCATE(lambday, (0:NS), zero)
  SALLOCATE(lambdaz, (0:NS), zero)
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
  SALLOCATE(absrp, (0:NS), zero)
  xt(0:coil(icoil)%NS) = coil(icoil)%xt(0:coil(icoil)%NS)
  yt(0:coil(icoil)%NS) = coil(icoil)%yt(0:coil(icoil)%NS)
  zt(0:coil(icoil)%NS) = coil(icoil)%zt(0:coil(icoil)%NS)
  xa(0:coil(icoil)%NS) = coil(icoil)%xa(0:coil(icoil)%NS)
  ya(0:coil(icoil)%NS) = coil(icoil)%ya(0:coil(icoil)%NS)
  za(0:coil(icoil)%NS) = coil(icoil)%za(0:coil(icoil)%NS)
  xb(0:coil(icoil)%NS) = coil(icoil)%xb(0:coil(icoil)%NS)
  yb(0:coil(icoil)%NS) = coil(icoil)%yb(0:coil(icoil)%NS)
  zb(0:coil(icoil)%NS) = coil(icoil)%zb(0:coil(icoil)%NS)

  FATAL( TorsDeriv0, icoil .lt. 1 .or. icoil .gt. Ncoils, icoil not in right range )

  torss = zero
  torsRet = zero

  lambdax(0:NS) = yt(0:NS)*za(0:NS) - zt(0:NS)*ya(0:NS)
  lambday(0:NS) = zt(0:NS)*xa(0:NS) - xt(0:NS)*za(0:NS)
  lambdaz(0:NS) = xt(0:NS)*ya(0:NS) - yt(0:NS)*xa(0:NS)
  lambdanorm(0:NS) = sqrt( lambdax(0:NS)**2 + lambday(0:NS)**2 + lambdaz(0:NS)**2 ) 
  lambdaunitx(0:NS) = lambdax(0:NS) / lambdanorm(0:NS)
  lambdaunity(0:NS) = lambday(0:NS) / lambdanorm(0:NS)
  lambdaunitz(0:NS) = lambdaz(0:NS) / lambdanorm(0:NS)
 
  ! Put in conditional for case_tau

  torss(0:NS) = ( lambdaunitx(0:NS)*xb(0:NS) + lambdaunity(0:NS)*yb(0:NS) + lambdaunitz(0:NS)*zb(0:NS) ) / lambdanorm(0:NS)
  
  coil(icoil)%minlambda = minval(lambdanorm)

  torss(0:NS) = torss(0:NS)*sqrt( xt(0:NS)**2 + yt(0:NS)**2 + zt(0:NS)**2 )

  torsRet = sum(torss)-torss(0)
  call lenDeriv0( icoil, coil(icoil)%L )
  torsRet = pi2*torsRet/(NS*coil(icoil)%L)

  DALLOCATE(torss)
  DALLOCATE(lambdax)
  DALLOCATE(lambday)
  DALLOCATE(lambdaz)
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
  DALLOCATE(absrp)

  return

end subroutine TorsDeriv0

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine TorsDeriv1(icoil, derivs, ND, NF) !Calculate all derivatives for a coil

  use globals, only: dp, zero, pi2, coil, DoF, myid, ounit, Ncoils, &
          case_tau, case_tors, tors_alpha, tors_beta, tors_gamma, tors0, &
          penfun_tors, FouCoil, MPI_COMM_FOCUS
  implicit none
  include "mpif.h"

  INTEGER, intent(in)  :: icoil, ND , NF
  REAL   , intent(out) :: derivs(1:1, 1:ND)
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  INTEGER              :: kseg, astat, ierr, doff, n, i, NS
  REAL                 :: L, tauavg
  REAL                 :: d1L(1:1, 1:ND)
  REAL,allocatable     :: tau(:), absrp(:), xt(:), yt(:), zt(:), xa(:), ya(:), za(:), xb(:), yb(:), zb(:), &
     lambdax(:), lambday(:), lambdaz(:), lambdanorm(:), lambdaunitx(:), lambdaunity(:), lambdaunitz(:), &
     dxdDoF(:,:), dydDoF(:,:), dzdDoF(:,:), dxtdDoF(:,:), dytdDoF(:,:), dztdDoF(:,:), dxadDoF(:,:), &
     dyadDoF(:,:), dzadDoF(:,:), dxbdDoF(:,:), dybdDoF(:,:), dzbdDoF(:,:), dlambdaxdDof(:,:), &
     dlambdaydDof(:,:), dlambdazdDof(:,:), dtaudDof(:,:), dtauavgdDof(:), dabsrpdDof(:,:), &
     dlambdanormdDof(:,:)

  FATAL( TorsDeriv1, icoil .lt. 1 .or. icoil .gt. Ncoils, icoil not in right range )
  FATAL( TorsDeriv1, ND .ne. 3+6*NF, Degrees of freedom are incorrect )

  NS = coil(icoil)%NS

  SALLOCATE(tau, (0:NS), zero)
  SALLOCATE(absrp, (0:NS), zero)
  SALLOCATE(lambdax, (0:NS), zero)
  SALLOCATE(lambday, (0:NS), zero)
  SALLOCATE(lambdaz, (0:NS), zero)
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
  xt(0:coil(icoil)%NS) = coil(icoil)%xt(0:coil(icoil)%NS)
  yt(0:coil(icoil)%NS) = coil(icoil)%yt(0:coil(icoil)%NS)
  zt(0:coil(icoil)%NS) = coil(icoil)%zt(0:coil(icoil)%NS)
  xa(0:coil(icoil)%NS) = coil(icoil)%xa(0:coil(icoil)%NS)
  ya(0:coil(icoil)%NS) = coil(icoil)%ya(0:coil(icoil)%NS)
  za(0:coil(icoil)%NS) = coil(icoil)%za(0:coil(icoil)%NS)
  xb(0:coil(icoil)%NS) = coil(icoil)%xb(0:coil(icoil)%NS)
  yb(0:coil(icoil)%NS) = coil(icoil)%yb(0:coil(icoil)%NS)
  zb(0:coil(icoil)%NS) = coil(icoil)%zb(0:coil(icoil)%NS)

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
  SALLOCATE(dtauavgdDof, (1:ND), zero)
  SALLOCATE(dabsrpdDof, (0:NS,1:ND), zero)
  SALLOCATE(dlambdanormdDof, (0:NS,1:ND), zero)
  
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
  lambdanorm(0:NS) = sqrt( lambdax(0:NS)**2 + lambday(0:NS)**2 + lambdaz(0:NS)**2 )
  lambdaunitx(0:NS) = lambdax(0:NS) / lambdanorm(0:NS)
  lambdaunity(0:NS) = lambday(0:NS) / lambdanorm(0:NS)
  lambdaunitz(0:NS) = lambdaz(0:NS) / lambdanorm(0:NS)

  ! Put in conditional for case tau 

  tau(0:NS) = ( lambdaunitx(0:NS)*xb(0:NS) + lambdaunity(0:NS)*yb(0:NS) + lambdaunitz(0:NS)*zb(0:NS) ) / lambdanorm(0:NS)
  absrp(0:NS) = sqrt( xt(0:NS)**2 + yt(0:NS)**2 + zt(0:NS)**2 )
  tauavg = pi2*( sum(tau(0:NS)*absrp(0:NS)) - tau(0)*absrp(0) ) / ( L * NS )

  do i = 1,ND
     dlambdaxdDof(0:NS,i) = dytdDof(0:NS,i)*za(0:NS) - dztdDof(0:NS,i)*ya(0:NS) &
                             + yt(0:NS)*dzadDof(0:NS,i) - zt(0:NS)*dyadDof(0:NS,i)
     dlambdaydDof(0:NS,i) = dztdDof(0:NS,i)*xa(0:NS) - dxtdDof(0:NS,i)*za(0:NS) &
                             + zt(0:NS)*dxadDof(0:NS,i) - xt(0:NS)*dzadDof(0:NS,i)
     dlambdazdDof(0:NS,i) = dxtdDof(0:NS,i)*ya(0:NS) - dytdDof(0:NS,i)*xa(0:NS) &
                             + xt(0:NS)*dyadDof(0:NS,i) - yt(0:NS)*dxadDof(0:NS,i)

     dtaudDof(0:NS,i) = -2*(lambdaunitx(0:NS)*dlambdaxdDof(0:NS,i) + lambdaunity(0:NS)*dlambdaydDof(0:NS,i) + &
     lambdaunitz(0:NS)*dlambdazdDof(0:NS,i) )*(lambdaunitx(0:NS)*xb(0:NS) + lambdaunity(0:NS)*yb(0:NS) + & 
     lambdaunitz(0:NS)*zb(0:NS)) / lambdanorm(0:NS)**2 + (dlambdaxdDof(0:NS,i)*xb(0:NS) + dlambdaydDof(0:NS,i)*yb(0:NS) + &
     dlambdazdDof(0:NS,i)*zb(0:NS) + lambdax(0:NS)*dxbdDof(0:NS,i) + lambday(0:NS)*dybdDof(0:NS,i) + &
     lambdaz(0:NS)*dzbdDof(0:NS,i)) / lambdanorm(0:NS)**2

     dabsrpdDof(0:NS,i) = ( xt(0:NS)**2 + yt(0:NS)**2 + zt(0:NS)**2 )**(-.5) * ( xt(0:NS)*dxtdDof(0:NS,i) + & 
                          yt(0:NS)*dytdDof(0:NS,i) + zt(0:NS)*dztdDof(0:NS,i) )

     dlambdanormdDof(0:NS,i) = lambdanorm(0:NS)**(-1)*(lambdax(0:NS)*dlambdaxdDof(0:NS,i) + &
             lambday(0:NS)*dlambdaydDof(0:NS,i) + lambdaz(0:NS)*dlambdazdDof(0:NS,i))
  enddo

  do i = 1,ND
     dtauavgdDof(i) = -1*d1L(1,i)*( sum(tau(0:NS)*absrp(0:NS)) - tau(0)*absrp(0) ) / L + &
                           sum(dtaudDof(0:NS,i)*absrp(0:NS)) - dtaudDof(0,i)*absrp(0) + &
                           sum(tau(0:NS)*dabsrpdDof(0:NS,i)) - tau(0)*dabsrpdDof(0,i)
  enddo
  dtauavgdDof(1:ND) = dtauavgdDof(1:ND)*pi2/(L*NS)

  if( case_tors .eq. 1 ) then
     if( tauavg**2 .ge. tors0**2 ) then
        if( penfun_tors .eq. 1 ) then
           derivs(1,1:ND) = 4*tors_alpha*( cosh(tors_alpha*(tauavg**2 - tors0**2)) - 1 )* &
                   sinh(tors_alpha*(tauavg**2 - tors0**2))*tauavg*dtauavgdDof(1:ND)
        else
           derivs(1,1:ND) = 2*tors_beta*tors_alpha*(tors_alpha*(tauavg**2-tors0**2))**(tors_beta-1)*tauavg*dtauavgdDof(1:ND)
        endif
     endif
  else
     derivs(1,1:ND) = 2*tors_gamma*(tauavg**2 + 1)**(tors_gamma-1)*tauavg*dtauavgdDof(1:ND)
  endif

  DALLOCATE(tau)
  DALLOCATE(absrp)
  DALLOCATE(lambdax)
  DALLOCATE(lambday)
  DALLOCATE(lambdaz)
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
  DALLOCATE(dtauavgdDof)
  DALLOCATE(dabsrpdDof)
  DALLOCATE(dlambdanormdDof)

  return

end subroutine TorsDeriv1


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
