!title (curvature) ! Calculate total curvature cost functon and its derivatives. (tkruger)

!latex \briefly{The constraint on curvature prevents the coil from becoming too complex and
!latex         is an important engineering constraint for the realization of feasible coils.
!latex         This function is still under development by Thomas Kruger and once the function
!latex         is ready and approved by Caoxiang Zhu, it will be pushed \emph{targt\_length}.}

!latex \calledby{\link{solvers}}

!latex  \section{General}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! curv is total penalty 
! chi = chi + weight_curv*curv
! t1CU is total derivative of penalty

! Not doing any L-M work

! not parallelized; still in development;
subroutine curvature(ideriv)
  use globals, only: dp, zero, half, pi2, machprec, ncpu, myid, ounit, &
       coil, DoF, Ncoils, Nfixgeo, Ndof, curv, t1CU, t2CU, weight_curv, FouCoil

  implicit none
  include "mpif.h"
  INTEGER, INTENT(in) :: ideriv

  INTEGER             :: astat, ierr, icoil, idof, ND, NF
  REAL                :: curvAdd 

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  curv = zero
  curvAdd = zero

  if( ideriv >= 0 ) then

     do icoil = 1, Ncoils
        if(coil(icoil)%itype .ne. 1) exit ! only for Fourier
        call CurvDeriv0(icoil,curvAdd)
        curv = curv + curvAdd; 
     enddo

     curv = curv / (Ncoils - Nfixgeo + machprec)

  endif

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if ( ideriv >= 1 ) then

     t1CU = zero

     idof = 0
     do icoil = 1, Ncoils

        if(coil(icoil)%itype .ne. 1) exit ! only for Fourier

        ND = DoF(icoil)%ND
        NF = FouCoil(icoil)%NF

        if ( coil(icoil)%Ic /= 0 ) then !if current is free;
           idof = idof + 1
        endif

        if ( coil(icoil)%Lc /= 0 ) then !if geometry is free;
           call CurvDeriv1( icoil, t1CU(idof+1:idof+ND), ND, NF )
           idof = idof + ND
        endif

     enddo !end icoil;
     FATAL( curvature , idof .ne. Ndof, counting error in packing )

     t1CU = t1CU / (Ncoils - Nfixgeo + machprec)

  endif

  return
end subroutine curvature

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! calculate coil curvature

subroutine CurvDeriv0(icoil,curvRet)

  use globals, only: dp, zero, pi2, ncpu, astat, ierr, myid, ounit, coil, NFcoil, Nseg, Ncoils, case_curv, k0, curv_alpha 

  implicit none
  include "mpif.h"

  INTEGER, intent(in)  :: icoil
  REAL   , intent(out) :: curvRet 
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  INTEGER              :: kseg, NS
  REAL,allocatable     :: curvv(:)

  SALLOCATE(curvv, (0:coil(icoil)%NS),zero)

  FATAL( CurvDeriv0, icoil .lt. 1 .or. icoil .gt. Ncoils, icoil not in right range )
 
  curvv = zero
  curvRet = zero

  curvv = sqrt( (coil(icoil)%za*coil(icoil)%yt-coil(icoil)%zt*coil(icoil)%ya)**2 &
             + (coil(icoil)%xa*coil(icoil)%zt-coil(icoil)%xt*coil(icoil)%za)**2  & 
             + (coil(icoil)%ya*coil(icoil)%xt-coil(icoil)%yt*coil(icoil)%xa)**2 )& 
             / ((coil(icoil)%xt)**2+(coil(icoil)%yt)**2+(coil(icoil)%zt)**2)**(1.5)
  coil(icoil)%maxcurv = maxval(curvv)

  if( case_curv == 1 ) then
     curvRet = sum(curvv)
     curvRet = curvRet/coil(icoil)%NS
  elseif( case_curv == 2 ) then
     curvv = curvv*curvv
     curvRet = sum(curvv) 
     curvRet = curvRet/coil(icoil)%NS
  elseif( case_curv == 3 ) then
     NS = coil(icoil)%NS
     ! Put in penalty function
     do kseg = 0,NS
        if( curvv(kseg) > k0 ) then
           
        end if
     enddo 
  else   
     FATAL( CurvDeriv0, .true. , invalid case_curv option ) 
  end if

  return

end subroutine CurvDeriv0

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine CurvDeriv1(icoil, derivs, ND, NF) !Calculate all derivatives for a coil

        use globals, only: dp, zero, pi2, coil, DoF, myid, ounit, Ncoils, case_curv
  implicit none
  include "mpif.h"

  INTEGER, intent(in)  :: icoil, ND, NF
  REAL   , intent(out) :: derivs(1:1, 1:ND)
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  INTEGER              :: kseg, astat, ierr, doff, nff, NS
  REAL                 :: dl3, xt, yt, zt, xa, ya, za, f1, f2, ncc, nss, ncp, nsp, curvHold
  REAL, dimension(1:1, 0:coil(icoil)%NS-1) :: dLx, dLy, dLz

  FATAL( CurvDeriv1, icoil .lt. 1 .or. icoil .gt. Ncoils, icoil not in right range )
  
  derivs = zero

  ! NF=0  
  derivs(1,1) = 0.0
  derivs(1,2*NF+2) = 0.0
  derivs(1,4*NF+3) = 0.0
  NS = coil(icoil)%NS

  do kseg = 0,NS

     xt = coil(icoil)%xt(kseg) ; yt = coil(icoil)%yt(kseg) ; zt = coil(icoil)%zt(kseg)
     xa = coil(icoil)%xa(kseg) ; ya = coil(icoil)%ya(kseg) ; za = coil(icoil)%za(kseg)
     f1 = sqrt( (xt*ya-xa*yt)**2 + (xt*za-xa*zt)**2 + (yt*za-ya*zt)**2  );
     f2 = (xt**2+yt**2+zt**2)**(1.5);

     do nff = 1,NF

        ncc =  -1.0*nff*sin(nff*pi2*kseg/NS)     !Use for x y and z
        ncp =  -1.0*nff*nff*cos(nff*pi2*kseg/NS)
        nss =       nff*cos(nff*pi2*kseg/NS)
        nsp = - 1.0*nff*nff*sin(nff*pi2*kseg/NS)
        
        if( case_curv == 1 ) then
           ! Xc 
           derivs(1,1+nff)      = derivs(1,1+nff)      + -1.0*(f1/(f2**2))*( 3.0*xt*sqrt(xt**2+yt**2+zt**2)*ncc) \
           + ((f1**(-1))*((xt*ya-xa*yt)*(ncc*ya-ncp*yt)+(xt*za-xa*zt)*(ncc*za-ncp*zt)) )/f2; 
           ! Xs
           derivs(1,1+nff+NF)   = derivs(1,1+nff+NF)   + -1.0*(f1/(f2**2))*( 3.0*xt*sqrt(xt**2+yt**2+zt**2)*nss) \
           + ((f1**(-1))*((xt*ya-xa*yt)*(nss*ya-nsp*yt)+(xt*za-xa*zt)*(nss*za-nsp*zt)) )/f2;
           ! Yc
           derivs(1,2+nff+2*NF) = derivs(1,2+nff+2*NF) + -1.0*(f1/(f2**2))*( 3.0*yt*sqrt(xt**2+yt**2+zt**2)*ncc) \
           + ((f1**(-1))*((xt*ya-xa*yt)*(ncp*xt-ncc*xa)+(yt*za-ya*zt)*(ncc*za-ncp*zt)) )/f2;  
           ! Ys
           derivs(1,2+nff+3*NF) = derivs(1,2+nff+3*NF) + -1.0*(f1/(f2**2))*( 3.0*yt*sqrt(xt**2+yt**2+zt**2)*nss) \
           + ((f1**(-1))*((xt*ya-xa*yt)*(nsp*xt-nss*xa)+(yt*za-ya*zt)*(nss*za-nsp*zt)) )/f2;
           ! Zc
           derivs(1,3+nff+4*NF) = derivs(1,3+nff+4*NF) + -1.0*(f1/(f2**2))*( 3.0*zt*sqrt(xt**2+yt**2+zt**2)*ncc) \
           + ((f1**(-1))*((xt*za-xa*zt)*(xt*ncp-xa*ncc)+(yt*za-ya*zt)*(yt*ncp-ya*ncc)) )/f2;
           ! Zs
           derivs(1,3+nff+5*NF) = derivs(1,3+nff+5*NF) + -1.0*(f1/(f2**2))*( 3.0*zt*sqrt(xt**2+yt**2+zt**2)*nss) \
           + ((f1**(-1))*((xt*za-xa*zt)*(xt*nsp-xa*nss)+(yt*za-ya*zt)*(yt*nsp-ya*nss)) )/f2;

        else if( case_curv == 2 ) then
           curvHold = sqrt( (za*yt-zt*ya)**2 + (xa*zt-xt*za)**2 + (ya*xt-yt*xa)**2 ) / ((xt)**2+(yt)**2+(zt)**2)**(1.5)
           ! Xc 
           derivs(1,1+nff)      = derivs(1,1+nff)      + (-1.0*(f1/(f2**2))*( 3.0*xt*sqrt(xt**2+yt**2+zt**2)*ncc) \
           + ((f1**(-1))*((xt*ya-xa*yt)*(ncc*ya-ncp*yt)+(xt*za-xa*zt)*(ncc*za-ncp*zt)) )/f2)*2.0*curvHold; 
           ! Xs
           derivs(1,1+nff+NF)   = derivs(1,1+nff+NF)   + (-1.0*(f1/(f2**2))*( 3.0*xt*sqrt(xt**2+yt**2+zt**2)*nss) \
           + ((f1**(-1))*((xt*ya-xa*yt)*(nss*ya-nsp*yt)+(xt*za-xa*zt)*(nss*za-nsp*zt)) )/f2)*2.0*curvHold;
           ! Yc
           derivs(1,2+nff+2*NF) = derivs(1,2+nff+2*NF) + (-1.0*(f1/(f2**2))*( 3.0*yt*sqrt(xt**2+yt**2+zt**2)*ncc) \
           + ((f1**(-1))*((xt*ya-xa*yt)*(ncp*xt-ncc*xa)+(yt*za-ya*zt)*(ncc*za-ncp*zt)) )/f2)*2.0*curvHold;  
           ! Ys
           derivs(1,2+nff+3*NF) = derivs(1,2+nff+3*NF) + (-1.0*(f1/(f2**2))*( 3.0*yt*sqrt(xt**2+yt**2+zt**2)*nss) \
           + ((f1**(-1))*((xt*ya-xa*yt)*(nsp*xt-nss*xa)+(yt*za-ya*zt)*(nss*za-nsp*zt)) )/f2)*2.0*curvHold;
           ! Zc
           derivs(1,3+nff+4*NF) = derivs(1,3+nff+4*NF) + (-1.0*(f1/(f2**2))*( 3.0*zt*sqrt(xt**2+yt**2+zt**2)*ncc) \
           + ((f1**(-1))*((xt*za-xa*zt)*(xt*ncp-xa*ncc)+(yt*za-ya*zt)*(yt*ncp-ya*ncc)) )/f2)*2.0*curvHold;
           ! Zs
           derivs(1,3+nff+5*NF) = derivs(1,3+nff+5*NF) + (-1.0*(f1/(f2**2))*( 3.0*zt*sqrt(xt**2+yt**2+zt**2)*nss) \
           + ((f1**(-1))*((xt*za-xa*zt)*(xt*nsp-xa*nss)+(yt*za-ya*zt)*(yt*nsp-ya*nss)) )/f2)*2.0*curvHold;
 
        else
           FATAL( CurvDeriv1, .true. , invalid case_curv option )  
        end if
     enddo
  enddo

  derivs = derivs/NS
  
  return

end subroutine CurvDeriv1


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
