!title (curvature) ! Calculate curvature objective functon and its derivatives. (tkruger)

!latex \briefly{The constraint on curvature prevents the coil from becoming too complex and
!latex         is an important engineering constraint for the realization of feasible coils.
!latex         This function is still under development. emph{targt\_length}.}

!latex \calledby{\link{solvers}}

!latex  \section{General}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! curv is total penalty 
! chi = chi + weight_curv*curv
! t1CU is total derivative of penalty
! LM implemented
! not parallelized, at some point check to see how long takes to run
subroutine curvature(ideriv)
  use globals, only: dp, zero, half, pi2, machprec, ncpu, myid, ounit, MPI_COMM_FOCUS, &
       coil, DoF, Ncoils, Nfixgeo, Ndof, curv, t1CU, t2CU, weight_curv, FouCoil, &
       mcurv, icurv, LM_fvec, LM_fjac

  implicit none
  include "mpif.h"
  INTEGER, INTENT(in) :: ideriv

  INTEGER             :: astat, ierr, icoil, idof, ND, NF, ivec
  REAL                :: curvAdd

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  curv = zero
  curvAdd = zero
  ivec = 1

  if( ideriv >= 0 ) then

     do icoil = 1, Ncoils
        if( coil(icoil)%type .ne. 1 ) exit ! only for Fourier
        if( coil(icoil)%Lc     /=  0 ) then ! if geometry is free
           call CurvDeriv0(icoil,curvAdd)
           curv = curv + curvAdd
           if (mcurv > 0) then ! L-M format of targets
              LM_fvec(icurv+ivec) = weight_curv*curvAdd
              ivec = ivec + 1
           endif
        endif 
     enddo

     curv = curv / (Ncoils - Nfixgeo + machprec)

  endif

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if ( ideriv >= 1 ) then

     t1CU = zero
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
           call CurvDeriv1( icoil, t1CU(idof+1:idof+ND), ND, NF )
           if (mcurv > 0) then ! L-M format of targets
              LM_fjac(icurv+ivec, idof+1:idof+ND) = weight_curv * t1CU(idof+1:idof+ND)
              ivec = ivec + 1
           endif
           idof = idof + ND
        endif

     enddo
     FATAL( curvature , idof .ne. Ndof, counting error in packing )

     t1CU = t1CU / (Ncoils - Nfixgeo + machprec)

  endif

  return
end subroutine curvature

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! calculate coil curvature

subroutine CurvDeriv0(icoil,curvRet)

  use globals, only: dp, zero, pi2, ncpu, astat, ierr, myid, ounit, coil, NFcoil, Nseg, Ncoils, &
          case_curv, curv_alpha, curv_c, k0, k1, curv_beta, curv_gamma, curv_sigma, penfun_curv, k1_len, MPI_COMM_FOCUS

  implicit none
  include "mpif.h"

  INTEGER, intent(in)  :: icoil
  REAL   , intent(out) :: curvRet 
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  INTEGER              :: kseg, NS
  REAL                 :: hypc, curv_hold, k1_use
  REAL,allocatable     :: curvv(:), xt(:), yt(:), zt(:), xa(:), ya(:), za(:)

  NS = coil(icoil)%NS 
  SALLOCATE(curvv, (0:NS), zero)
  SALLOCATE(xt, (0:NS), zero)
  SALLOCATE(yt, (0:NS), zero)
  SALLOCATE(zt, (0:NS), zero)
  SALLOCATE(xa, (0:NS), zero)
  SALLOCATE(ya, (0:NS), zero)
  SALLOCATE(za, (0:NS), zero)
  xt(0:coil(icoil)%NS) = coil(icoil)%xt(0:coil(icoil)%NS)
  yt(0:coil(icoil)%NS) = coil(icoil)%yt(0:coil(icoil)%NS)
  zt(0:coil(icoil)%NS) = coil(icoil)%zt(0:coil(icoil)%NS)
  xa(0:coil(icoil)%NS) = coil(icoil)%xa(0:coil(icoil)%NS)
  ya(0:coil(icoil)%NS) = coil(icoil)%ya(0:coil(icoil)%NS)
  za(0:coil(icoil)%NS) = coil(icoil)%za(0:coil(icoil)%NS)

  FATAL( CurvDeriv0, icoil .lt. 1 .or. icoil .gt. Ncoils, icoil not in right range )

  if ( penfun_curv .ne. 1 .and. penfun_curv .ne. 2 ) then
     FATAL( CurvDeriv0, .true. , invalid choice of penfun_curv, pick 1 or 2 )
  endif
  if ( k0 .le. 0.0 ) then
     FATAL( CurvDeriv0, .true. , k0 needs to be positive )
  endif
  if ( k1 < 0.0 ) then
     FATAL( CurvDeriv0, .true. , k1 needs to be greater than or equal to 0 )
  endif
  if ( curv_alpha < 0.0 ) then
     FATAL( CurvDeriv0, .true. , curv_alpha needs to be greater than or equal to 0 )
  endif
  if ( curv_beta .le. 1.0 ) then
     FATAL( CurvDeriv0, .true. , curv_beta needs to be greater than 1 )
  endif
  if ( curv_gamma < 1.0 ) then
     FATAL( CurvDeriv0, .true. , curv_gamma needs to be greater than 1 )
  endif
  if ( curv_sigma < 0.0 ) then
     FATAL( CurvDeriv0, .true. , curv_sigma needs to be greater than or equal to 0 )
  endif
 
  curvv = zero
  curvRet = zero

  curvv = sqrt( (za*yt-zt*ya)**2 + (xa*zt-xt*za)**2  + (ya*xt-yt*xa)**2 ) / (xt**2+yt**2+zt**2)**(1.5)
  coil(icoil)%maxcurv = maxval(curvv)

  if( case_curv == 1 ) then ! linear
     curvRet = sum(curvv)-curvv(0)
     curvRet = pi2*curvRet/NS
  elseif( case_curv == 2 ) then ! quadratic
     curvv = curvv*curvv
     curvRet = sum(curvv)-curvv(0)
     curvRet = pi2*curvRet/NS
  elseif( case_curv == 3 ) then ! penalty method 
     do kseg = 0,NS-1
        if( curvv(kseg) > k0 ) then
            hypc = 0.5*exp( curv_alpha*( curvv(kseg) - k0 ) ) + 0.5*exp( -1.0*curv_alpha*( curvv(kseg) - k0 ) )
            curvRet = curvRet + (hypc - 1.0)**2
        endif
     enddo
     curvRet = pi2*curvRet/NS
  elseif( case_curv == 4 ) then ! penalty plus complexity 
     if ( k1_len .eq. 1 ) then
        k1_use = pi2/coil(icoil)%Lo
     else
        k1_use = k1
     endif
     do kseg = 0,NS-1
        curv_hold = 0.0
        if ( curvv(kseg) > k0 ) then
           if ( penfun_curv .eq. 1 ) then
              hypc = 0.5*exp( curv_alpha*( curvv(kseg) - k0 ) ) + 0.5*exp( -1.0*curv_alpha*( curvv(kseg) - k0 ) )
              curv_hold = (hypc - 1.0)**2
           else
              curv_hold = ( curv_alpha*(curvv(kseg)-k0) )**curv_beta
           endif
        endif
        if ( curvv(kseg) > k1_use ) then
           curv_hold = curv_hold + curv_sigma*( ( curvv(kseg) - k1_use )**curv_gamma )
        endif
        curvRet = curvRet + curv_hold*sqrt(xt(kseg)**2+yt(kseg)**2+zt(kseg)**2)
     enddo
     call lenDeriv0( icoil, coil(icoil)%L )
     curvRet = pi2*curvRet/(NS*coil(icoil)%L)
  else   
     FATAL( CurvDeriv0, .true. , invalid case_curv option ) 
  endif

  DALLOCATE(curvv)
  DALLOCATE(xt)
  DALLOCATE(yt)
  DALLOCATE(zt)
  DALLOCATE(xa)
  DALLOCATE(ya)
  DALLOCATE(za)

  return

end subroutine CurvDeriv0

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine CurvDeriv1(icoil, derivs, ND, NF) !Calculate all derivatives for a coil

  use globals, only: dp, zero, pi2, coil, DoF, myid, ounit, Ncoils, &
                case_curv, curv_alpha, k0, curv_c, k1, curv_beta, curv_gamma, curv_sigma, penfun_curv, k1_len, FouCoil, MPI_COMM_FOCUS
  implicit none
  include "mpif.h"

  INTEGER, intent(in)  :: icoil, ND , NF
  REAL   , intent(out) :: derivs(1:1, 1:ND)
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  INTEGER              :: kseg, astat, ierr, doff, nff, NS
  REAL                 :: dl3, xt, yt, zt, xa, ya, za, f1, f2, ncc, nss, ncp, nsp, curvHold, penCurv, curvH, H0, H1, leng, hypc, hyps, curv_deriv, k1_use
  REAL                 :: d1L(1:1, 1:ND)

  FATAL( CurvDeriv1, icoil .lt. 1 .or. icoil .gt. Ncoils, icoil not in right range )

  NS = coil(icoil)%NS
  derivs = zero
  d1L = zero
  call lenDeriv0( icoil, coil(icoil)%L )
  leng = coil(icoil)%L
  call lenDeriv1( icoil, d1L(1:1,1:ND), ND )
  
  do kseg = 0,NS-1

     xt = coil(icoil)%xt(kseg) ; yt = coil(icoil)%yt(kseg) ; zt = coil(icoil)%zt(kseg) ;
     xa = coil(icoil)%xa(kseg) ; ya = coil(icoil)%ya(kseg) ; za = coil(icoil)%za(kseg) ;
     f1 = sqrt( (xt*ya-xa*yt)**2 + (xt*za-xa*zt)**2 + (yt*za-ya*zt)**2 )
     f2 = (xt**2+yt**2+zt**2)**(1.5)

     do nff = 1,NF
        ncc = -1.0*nff    *FouCoil(icoil)%smt(kseg,nff)
        ncp = -1.0*nff*nff*FouCoil(icoil)%cmt(kseg,nff)
        nss =      nff    *FouCoil(icoil)%cmt(kseg,nff)
        nsp = -1.0*nff*nff*FouCoil(icoil)%smt(kseg,nff)

        if( case_curv == 1 ) then
           ! Xc 
           derivs(1,1+nff)      = derivs(1,1+nff)      + -1.0*(f1/(f2**2))*( 3.0*xt*sqrt(xt**2+yt**2+zt**2)*ncc) &
           + ((f1**(-1))*((xt*ya-xa*yt)*(ncc*ya-ncp*yt)+(xt*za-xa*zt)*(ncc*za-ncp*zt)) )/f2 
           ! Xs
           derivs(1,1+nff+NF)   = derivs(1,1+nff+NF)   + -1.0*(f1/(f2**2))*( 3.0*xt*sqrt(xt**2+yt**2+zt**2)*nss) &
           + ((f1**(-1))*((xt*ya-xa*yt)*(nss*ya-nsp*yt)+(xt*za-xa*zt)*(nss*za-nsp*zt)) )/f2
           ! Yc
           derivs(1,2+nff+2*NF) = derivs(1,2+nff+2*NF) + -1.0*(f1/(f2**2))*( 3.0*yt*sqrt(xt**2+yt**2+zt**2)*ncc) &
           + ((f1**(-1))*((xt*ya-xa*yt)*(ncp*xt-ncc*xa)+(yt*za-ya*zt)*(ncc*za-ncp*zt)) )/f2  
           ! Ys
           derivs(1,2+nff+3*NF) = derivs(1,2+nff+3*NF) + -1.0*(f1/(f2**2))*( 3.0*yt*sqrt(xt**2+yt**2+zt**2)*nss) &
           + ((f1**(-1))*((xt*ya-xa*yt)*(nsp*xt-nss*xa)+(yt*za-ya*zt)*(nss*za-nsp*zt)) )/f2
           ! Zc
           derivs(1,3+nff+4*NF) = derivs(1,3+nff+4*NF) + -1.0*(f1/(f2**2))*( 3.0*zt*sqrt(xt**2+yt**2+zt**2)*ncc) &
           + ((f1**(-1))*((xt*za-xa*zt)*(xt*ncp-xa*ncc)+(yt*za-ya*zt)*(yt*ncp-ya*ncc)) )/f2
           ! Zs
           derivs(1,3+nff+5*NF) = derivs(1,3+nff+5*NF) + -1.0*(f1/(f2**2))*( 3.0*zt*sqrt(xt**2+yt**2+zt**2)*nss) &
           + ((f1**(-1))*((xt*za-xa*zt)*(xt*nsp-xa*nss)+(yt*za-ya*zt)*(yt*nsp-ya*nss)) )/f2

        elseif( case_curv == 2 ) then
           curvHold = sqrt( (za*yt-zt*ya)**2 + (xa*zt-xt*za)**2 + (ya*xt-yt*xa)**2 ) / ((xt)**2+(yt)**2+(zt)**2)**(1.5)
           ! Xc 
           derivs(1,1+nff)      = derivs(1,1+nff)      + (-1.0*(f1/(f2**2))*( 3.0*xt*sqrt(xt**2+yt**2+zt**2)*ncc) &
           + ((f1**(-1))*((xt*ya-xa*yt)*(ncc*ya-ncp*yt)+(xt*za-xa*zt)*(ncc*za-ncp*zt)) )/f2)*2.0*curvHold 
           ! Xs
           derivs(1,1+nff+NF)   = derivs(1,1+nff+NF)   + (-1.0*(f1/(f2**2))*( 3.0*xt*sqrt(xt**2+yt**2+zt**2)*nss) &
           + ((f1**(-1))*((xt*ya-xa*yt)*(nss*ya-nsp*yt)+(xt*za-xa*zt)*(nss*za-nsp*zt)) )/f2)*2.0*curvHold
           ! Yc
           derivs(1,2+nff+2*NF) = derivs(1,2+nff+2*NF) + (-1.0*(f1/(f2**2))*( 3.0*yt*sqrt(xt**2+yt**2+zt**2)*ncc) &
           + ((f1**(-1))*((xt*ya-xa*yt)*(ncp*xt-ncc*xa)+(yt*za-ya*zt)*(ncc*za-ncp*zt)) )/f2)*2.0*curvHold  
           ! Ys
           derivs(1,2+nff+3*NF) = derivs(1,2+nff+3*NF) + (-1.0*(f1/(f2**2))*( 3.0*yt*sqrt(xt**2+yt**2+zt**2)*nss) &
           + ((f1**(-1))*((xt*ya-xa*yt)*(nsp*xt-nss*xa)+(yt*za-ya*zt)*(nss*za-nsp*zt)) )/f2)*2.0*curvHold
           ! Zc
           derivs(1,3+nff+4*NF) = derivs(1,3+nff+4*NF) + (-1.0*(f1/(f2**2))*( 3.0*zt*sqrt(xt**2+yt**2+zt**2)*ncc) &
           + ((f1**(-1))*((xt*za-xa*zt)*(xt*ncp-xa*ncc)+(yt*za-ya*zt)*(yt*ncp-ya*ncc)) )/f2)*2.0*curvHold
           ! Zs
           derivs(1,3+nff+5*NF) = derivs(1,3+nff+5*NF) + (-1.0*(f1/(f2**2))*( 3.0*zt*sqrt(xt**2+yt**2+zt**2)*nss) &
           + ((f1**(-1))*((xt*za-xa*zt)*(xt*nsp-xa*nss)+(yt*za-ya*zt)*(yt*nsp-ya*nss)) )/f2)*2.0*curvHold

        elseif( case_curv == 3 ) then
           curvHold = sqrt( (za*yt-zt*ya)**2 + (xa*zt-xt*za)**2 + (ya*xt-yt*xa)**2 ) / ((xt)**2+(yt)**2+(zt)**2)**(1.5)
           hypc = 0.5*exp( curv_alpha*(curvHold-k0) ) + 0.5*exp( -1.0*curv_alpha*(curvHold-k0) )
           hyps = 0.5*exp( curv_alpha*(curvHold-k0) ) - 0.5*exp( -1.0*curv_alpha*(curvHold-k0) )
           penCurv = 2.0*( hypc - 1.0 )*hyps*curv_alpha
           if( curvHold > k0 ) then
              ! Xc
              derivs(1,1+nff)      = derivs(1,1+nff)      + (-1.0*(f1/(f2**2))*( 3.0*xt*sqrt(xt**2+yt**2+zt**2)*ncc) &
                            + ((f1**(-1))*((xt*ya-xa*yt)*(ncc*ya-ncp*yt)+(xt*za-xa*zt)*(ncc*za-ncp*zt)) )/f2)*penCurv
              ! Xs
              derivs(1,1+nff+NF)   = derivs(1,1+nff+NF)   + (-1.0*(f1/(f2**2))*( 3.0*xt*sqrt(xt**2+yt**2+zt**2)*nss) &
                            + ((f1**(-1))*((xt*ya-xa*yt)*(nss*ya-nsp*yt)+(xt*za-xa*zt)*(nss*za-nsp*zt)) )/f2)*penCurv
              ! Yc
              derivs(1,2+nff+2*NF) = derivs(1,2+nff+2*NF) + (-1.0*(f1/(f2**2))*( 3.0*yt*sqrt(xt**2+yt**2+zt**2)*ncc) &
                            + ((f1**(-1))*((xt*ya-xa*yt)*(ncp*xt-ncc*xa)+(yt*za-ya*zt)*(ncc*za-ncp*zt)) )/f2)*penCurv
              ! Ys
              derivs(1,2+nff+3*NF) = derivs(1,2+nff+3*NF) + (-1.0*(f1/(f2**2))*( 3.0*yt*sqrt(xt**2+yt**2+zt**2)*nss) &
                            + ((f1**(-1))*((xt*ya-xa*yt)*(nsp*xt-nss*xa)+(yt*za-ya*zt)*(nss*za-nsp*zt)) )/f2)*penCurv
              ! Zc
              derivs(1,3+nff+4*NF) = derivs(1,3+nff+4*NF) + (-1.0*(f1/(f2**2))*( 3.0*zt*sqrt(xt**2+yt**2+zt**2)*ncc) &
                            + ((f1**(-1))*((xt*za-xa*zt)*(xt*ncp-xa*ncc)+(yt*za-ya*zt)*(yt*ncp-ya*ncc)) )/f2)*penCurv
              ! Zs
              derivs(1,3+nff+5*NF) = derivs(1,3+nff+5*NF) + (-1.0*(f1/(f2**2))*( 3.0*zt*sqrt(xt**2+yt**2+zt**2)*nss) &
                            + ((f1**(-1))*((xt*za-xa*zt)*(xt*nsp-xa*nss)+(yt*za-ya*zt)*(yt*nsp-ya*nss)) )/f2)*penCurv
           else
              derivs(1,1+nff) = derivs(1,1+nff)
              derivs(1,1+nff+NF) = derivs(1,1+nff+NF)
              derivs(1,2+nff+2*NF) = derivs(1,2+nff+2*NF)
              derivs(1,2+nff+3*NF) = derivs(1,2+nff+3*NF)
              derivs(1,3+nff+4*NF) = derivs(1,3+nff+4*NF)
              derivs(1,3+nff+5*NF) = derivs(1,3+nff+5*NF)
           endif

        elseif( case_curv == 4 ) then
           curvHold = f1/f2
           penCurv = 0.0
           curv_deriv = 0.0
           if ( k1_len .eq. 1 ) then
              k1_use = pi2/coil(icoil)%Lo
           else
              k1_use = k1
           endif
           if ( penfun_curv .eq. 1 ) then
              if ( curvHold > k0 ) then
                 hypc = 0.5*exp( curv_alpha*(curvHold-k0) ) + 0.5*exp( -1.0*curv_alpha*(curvHold-k0) )
                 hyps = 0.5*exp( curv_alpha*(curvHold-k0) ) - 0.5*exp( -1.0*curv_alpha*(curvHold-k0) )
                 penCurv = ( hypc - 1.0 )**2
                 curv_deriv = 2.0*curv_alpha*( hypc - 1.0 )*hyps
              endif
           else
              if ( curvHold > k0 ) then
                 penCurv = (curv_alpha*(curvHold-k0))**curv_beta
                 curv_deriv = curv_beta*curv_alpha*( (curv_alpha*(curvHold-k0))**(curv_beta-1.0) )
              endif
           endif
           if ( curvHold > k1_use ) then
              curv_deriv = curv_deriv + curv_sigma*curv_gamma*( (curvHold-k1_use)**(curv_gamma-1.0) )
              penCurv = penCurv + curv_sigma*( (curvHold-k1_use)**curv_gamma )
           endif
           ! Xc
           derivs(1,1+nff)     = derivs(1,1+nff)       + sqrt(xt**2+yt**2+zt**2)*penCurv*-1.0*d1L(   1,1+nff)/leng
           derivs(1,1+nff)     = derivs(1,1+nff)       + xt*ncc*(xt**2+yt**2+zt**2)**(-.5)*penCurv
           derivs(1,1+nff)     = derivs(1,1+nff)       + sqrt(xt**2+yt**2+zt**2) &
           *(-1.0*(f1/(f2**2))*( 3.0*xt*sqrt(xt**2+yt**2+zt**2)*ncc) + ((f1**(-1))*((xt*ya-xa*yt)*(ncc*ya-ncp*yt)+(xt*za-xa*zt)*(ncc*za-ncp*zt)) )/f2)*curv_deriv
           ! Xs
           derivs(1,1+nff+NF)  = derivs(1,1+nff+NF)    + sqrt(xt**2+yt**2+zt**2)*penCurv*-1.0*d1L(1,1+nff+NF)/leng
           derivs(1,1+nff+NF)  = derivs(1,1+nff+NF)    + xt*nss*(xt**2+yt**2+zt**2)**(-.5)*penCurv
           derivs(1,1+nff+NF)  = derivs(1,1+nff+NF)    + sqrt(xt**2+yt**2+zt**2) &
           *(-1.0*(f1/(f2**2))*( 3.0*xt*sqrt(xt**2+yt**2+zt**2)*nss) + ((f1**(-1))*((xt*ya-xa*yt)*(nss*ya-nsp*yt)+(xt*za-xa*zt)*(nss*za-nsp*zt)) )/f2)*curv_deriv
           ! Yc
           derivs(1,2+nff+2*NF) = derivs(1,2+nff+2*NF) + sqrt(xt**2+yt**2+zt**2)*penCurv*-1.0*d1L(1,2+nff+2*NF)/leng
           derivs(1,2+nff+2*NF) = derivs(1,2+nff+2*NF) + yt*ncc*(xt**2+yt**2+zt**2)**(-.5)*penCurv
           derivs(1,2+nff+2*NF) = derivs(1,2+nff+2*NF) + sqrt(xt**2+yt**2+zt**2) &
           *(-1.0*(f1/(f2**2))*( 3.0*yt*sqrt(xt**2+yt**2+zt**2)*ncc) + ((f1**(-1))*((xt*ya-xa*yt)*(ncp*xt-ncc*xa)+(yt*za-ya*zt)*(ncc*za-ncp*zt)) )/f2)*curv_deriv
           ! Ys
           derivs(1,2+nff+3*NF) = derivs(1,2+nff+3*NF) + sqrt(xt**2+yt**2+zt**2)*penCurv*-1.0*d1L(1,2+nff+3*NF)/leng
           derivs(1,2+nff+3*NF) = derivs(1,2+nff+3*NF) + yt*nss*(xt**2+yt**2+zt**2)**(-.5)*penCurv
           derivs(1,2+nff+3*NF) = derivs(1,2+nff+3*NF) + sqrt(xt**2+yt**2+zt**2) &
           *(-1.0*(f1/(f2**2))*( 3.0*yt*sqrt(xt**2+yt**2+zt**2)*nss) + ((f1**(-1))*((xt*ya-xa*yt)*(nsp*xt-nss*xa)+(yt*za-ya*zt)*(nss*za-nsp*zt)) )/f2)*curv_deriv
           ! Zc
           derivs(1,3+nff+4*NF) = derivs(1,3+nff+4*NF) + sqrt(xt**2+yt**2+zt**2)*penCurv*-1.0*d1L(1,3+nff+4*NF)/leng
           derivs(1,3+nff+4*NF) = derivs(1,3+nff+4*NF) + zt*ncc*(xt**2+yt**2+zt**2)**(-.5)*penCurv
           derivs(1,3+nff+4*NF) = derivs(1,3+nff+4*NF) + sqrt(xt**2+yt**2+zt**2) &
           *(-1.0*(f1/(f2**2))*( 3.0*zt*sqrt(xt**2+yt**2+zt**2)*ncc) + ((f1**(-1))*((xt*za-xa*zt)*(xt*ncp-xa*ncc)+(yt*za-ya*zt)*(yt*ncp-ya*ncc)) )/f2)*curv_deriv
           ! Zs
           derivs(1,3+nff+5*NF) = derivs(1,3+nff+5*NF) + sqrt(xt**2+yt**2+zt**2)*penCurv*-1.0*d1L(1,3+nff+5*NF)/leng
           derivs(1,3+nff+5*NF) = derivs(1,3+nff+5*NF) + zt*nss*(xt**2+yt**2+zt**2)**(-.5)*penCurv
           derivs(1,3+nff+5*NF) = derivs(1,3+nff+5*NF) + sqrt(xt**2+yt**2+zt**2) &
           *(-1.0*(f1/(f2**2))*( 3.0*zt*sqrt(xt**2+yt**2+zt**2)*nss) + ((f1**(-1))*((xt*za-xa*zt)*(xt*nsp-xa*nss)+(yt*za-ya*zt)*(yt*nsp-ya*nss)) )/f2)*curv_deriv
        else
           FATAL( CurvDeriv1, .true. , invalid case_curv option )  
        endif
     enddo
  enddo
  
  derivs = pi2*derivs/NS
  if( case_curv == 4 ) derivs = derivs/leng

  return

end subroutine CurvDeriv1


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
