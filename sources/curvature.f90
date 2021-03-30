!title (curvature) ! Calculate curvature objective functon and its derivatives. (tkruger)

!latex \briefly{The curvature objective function can be used to minimize the coils' 
!latex         average curvature or it can constrain the curvature or both. 
!latex         Functional derivatives are not implemented. This implementation only works
!latex         for coils that use the Fourier series parameterization.
!latex         emph{targt\_length}.}

!latex \calledby{\link{solvers}}

!latex  \section{General}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! curv is total curvature penalty 
! chi = chi + weight_curv*curv
! t1K is total derivative of curvature penalty
! LM implemented
! not parallelized, does not take long to calculate 
subroutine curvature(ideriv)
  use globals, only: dp, zero, half, pi2, machprec, ncpu, myid, ounit, MPI_COMM_FOCUS, &
       coil, DoF, Ncoils, Nfixgeo, Ndof, curv, t1K, t2K, weight_curv, FouCoil, &
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

     t1K = zero
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
           call CurvDeriv1( icoil, t1K(idof+1:idof+ND), ND, NF )
           if (mcurv > 0) then ! L-M format of targets
              LM_fjac(icurv+ivec, idof+1:idof+ND) = weight_curv * t1K(idof+1:idof+ND)
              ivec = ivec + 1
           endif
           idof = idof + ND
        endif

     enddo
     FATAL( curvature , idof .ne. Ndof, counting error in packing )

     t1K = t1K / (Ncoils - Nfixgeo + machprec)

  endif

  return
end subroutine curvature

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! calculate coil curvature

subroutine CurvDeriv0(icoil,curvRet)

  use globals, only: dp, zero, pi2, ncpu, astat, ierr, myid, ounit, coil, NFcoil, Nseg, Ncoils, &
          case_curv, curv_alpha, k0, k1, curv_beta, curv_gamma, curv_sigma, penfun_curv, k1_len, MPI_COMM_FOCUS

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
  if ( k0 < 0.0 ) then
     FATAL( CurvDeriv0, .true. , k0 cannot be negative )
  endif
  if ( k1 < 0.0 ) then
     FATAL( CurvDeriv0, .true. , k1 cannot be negative )
  endif
  if ( curv_alpha < 0.0 ) then
     FATAL( CurvDeriv0, .true. , curv_alpha cannot be negative )
  endif
  if ( curv_beta < 2.0 ) then
     FATAL( CurvDeriv0, .true. , curv_beta needs to be >= 2 )
  endif
  if ( curv_gamma < 1.0 ) then
     FATAL( CurvDeriv0, .true. , curv_gamma needs to be >= 1 )
  endif
  if ( curv_sigma < 0.0 ) then
     FATAL( CurvDeriv0, .true. , curv_sigma cannot be negative )
  endif
  if ( curv_gamma .eq. 1.0 .and. k1 .ne. 0.0 ) then
     FATAL( CurvDeriv0, .true. , if curv_gamma = 1, k1 must = 0 )
  endif

  ! Set variable based on case_curv
  if ( case_curv .eq. 1 ) then
     curv_alpha = 0.0
     curv_sigma = 1.0
     curv_gamma = 1.0
     k1 = 0.0
  elseif ( case_curv .eq. 2 ) then
     curv_alpha = 0.0
     curv_sigma = 1.0
     curv_gamma = 2.0
     k1 = 0.0
  elseif (case_curv .eq. 3 ) then
     curv_sigma = 0.0
  else 
     ! Do nothing 
  endif
 
  curvv = zero
  curvRet = zero

  curvv = sqrt( (za*yt-zt*ya)**2 + (xa*zt-xt*za)**2  + (ya*xt-yt*xa)**2 ) / (xt**2+yt**2+zt**2)**(1.5)
  coil(icoil)%maxcurv = maxval(curvv)

  if ( k1_len .eq. 1 ) then
     k1_use = pi2/coil(icoil)%Lo
  else
     k1_use = k1
  endif

  do kseg = 0,NS-1
     curv_hold = 0.0
     if ( curv_alpha .ne. 0.0 ) then
        if ( curvv(kseg) > k0 ) then
           if ( penfun_curv .eq. 1 ) then
              hypc = 0.5*exp( curv_alpha*( curvv(kseg) - k0 ) ) + 0.5*exp( -1.0*curv_alpha*( curvv(kseg) - k0 ) )
              curv_hold = (hypc - 1.0)**2
           else
              curv_hold = ( curv_alpha*(curvv(kseg)-k0) )**curv_beta
           endif
        endif
     endif
     if ( curv_sigma .ne. 0.0 ) then
        if ( curvv(kseg) > k1_use ) then
           curv_hold = curv_hold + curv_sigma*( ( curvv(kseg) - k1_use )**curv_gamma )
        endif
     endif
     curvRet = curvRet + curv_hold*sqrt(xt(kseg)**2+yt(kseg)**2+zt(kseg)**2)
  enddo

  call lenDeriv0( icoil, coil(icoil)%L )
  curvRet = pi2*curvRet/(NS*coil(icoil)%L)

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
                case_curv, curv_alpha, k0, k1, curv_beta, curv_gamma, curv_sigma, penfun_curv, k1_len, FouCoil, MPI_COMM_FOCUS
  implicit none
  include "mpif.h"

  INTEGER, intent(in)  :: icoil, ND , NF
  REAL   , intent(out) :: derivs(1:1, 1:ND)
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  INTEGER              :: kseg, astat, ierr, doff, nff, NS, n
  REAL                 :: xt, yt, zt, xa, ya, za, f1, f2, curvHold, penCurv, leng, hypc, hyps, curv_deriv, &
                          k1_use, rtxrax, rtxray, rtxraz
  REAL                 :: d1L(1:1, 1:ND)
  REAL,allocatable     :: dxtdDoF(:,:), dytdDoF(:,:), dztdDoF(:,:), dxadDoF(:,:), dyadDoF(:,:), dzadDoF(:,:)

  FATAL( CurvDeriv1, icoil .lt. 1 .or. icoil .gt. Ncoils, icoil not in right range )

  NS = coil(icoil)%NS
  
  SALLOCATE(dxtdDoF, (0:NS,1:ND), zero)
  SALLOCATE(dytdDoF, (0:NS,1:ND), zero)
  SALLOCATE(dztdDoF, (0:NS,1:ND), zero)
  SALLOCATE(dxadDoF, (0:NS,1:ND), zero)
  SALLOCATE(dyadDoF, (0:NS,1:ND), zero)
  SALLOCATE(dzadDoF, (0:NS,1:ND), zero)
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
  enddo

  derivs = zero
  d1L = zero
  call lenDeriv0( icoil, coil(icoil)%L )
  leng = coil(icoil)%L
  call lenDeriv1( icoil, d1L(1:1,1:ND), ND )

  if ( case_curv .eq. 1 ) then
     curv_alpha = 0.0
     curv_sigma = 1.0
     curv_gamma = 1.0
     k1 = 0.0
  elseif ( case_curv .eq. 2 ) then
     curv_alpha = 0.0
     curv_sigma = 1.0
     curv_gamma = 2.0
     k1 = 0.0
  elseif (case_curv .eq. 3 ) then
     curv_sigma = 0.0
  else
     ! Do nothing
  endif
  
  do kseg = 0,NS-1
     xt = coil(icoil)%xt(kseg) ; yt = coil(icoil)%yt(kseg) ; zt = coil(icoil)%zt(kseg) ;
     xa = coil(icoil)%xa(kseg) ; ya = coil(icoil)%ya(kseg) ; za = coil(icoil)%za(kseg) ;
     rtxrax = yt*za - zt*ya
     rtxray = zt*xa - xt*za
     rtxraz = xt*ya - yt*xa
     f1 = sqrt( (xt*ya-xa*yt)**2 + (xt*za-xa*zt)**2 + (yt*za-ya*zt)**2 )
     f2 = (xt**2+yt**2+zt**2)**1.5
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

     derivs(1,1:ND) = derivs(1,1:ND) + sqrt(xt**2+yt**2+zt**2)*penCurv*-1.0*d1L(1,1:ND)/leng
     derivs(1,1:ND) = derivs(1,1:ND) + (xt*dxtdDof(kseg,1:ND)+yt*dytdDof(kseg,1:ND)+zt*dztdDof(kseg,1:ND))*(xt**2+yt**2+zt**2)**(-.5)*penCurv
     derivs(1,1:ND) = derivs(1,1:ND) + sqrt(xt**2+yt**2+zt**2) &
     *(-3.0*(f1/f2**2)*sqrt(xt**2+yt**2+zt**2)*(xt*dxtdDof(kseg,1:ND)+yt*dytdDof(kseg,1:ND)+zt*dztdDof(kseg,1:ND)) &
     + ( rtxrax*(dytdDof(kseg,1:ND)*za - dztdDof(kseg,1:ND)*ya + yt*dzadDof(kseg,1:ND) - zt*dyadDof(kseg,1:ND)) &
     +   rtxray*(dztdDof(kseg,1:ND)*xa - dxtdDof(kseg,1:ND)*za + zt*dxadDof(kseg,1:ND) - xt*dzadDof(kseg,1:ND)) &
     +   rtxraz*(dxtdDof(kseg,1:ND)*ya - dytdDof(kseg,1:ND)*xa + xt*dyadDof(kseg,1:ND) - yt*dxadDof(kseg,1:ND)) )/(f1*f2) )*curv_deriv
  enddo
  
  derivs = pi2*derivs/(NS*leng)

  return

end subroutine CurvDeriv1


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
