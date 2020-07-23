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
       mcurv, icurv, LM_fvec, LM_fjac,coil_type_spline,Splines

  implicit none
  include "mpif.h"
  INTEGER, INTENT(in) :: ideriv

  INTEGER             :: astat, ierr, icoil, idof, ND, NF,NCP, ivec
  REAL                :: curvAdd

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  curv = zero
  curvAdd = zero
  ivec = 1

  if( ideriv >= 0 ) then

     do icoil = 1, Ncoils
        if( (coil(icoil)%type .ne. 1) .AND. (coil(icoil)%type .ne. coil_type_spline) ) exit
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

     if( (coil(icoil)%type .ne. 1) .AND. (coil(icoil)%type .ne. coil_type_spline) ) exit

        ND = DoF(icoil)%ND
        if (coil(icoil)%type == 1) NF = FouCoil(icoil)%NF
        if (coil(icoil)%type == coil_type_spline) NCP = Splines(icoil)%NCP

        if ( coil(icoil)%Ic /= 0 ) then !if current is free;
           idof = idof + 1
        endif

        if ( coil(icoil)%Lc /= 0 ) then !if geometry is free;
           if (coil(icoil)%type == coil_type_spline) call CurvDeriv1( icoil, t1CU(idof+1:idof+ND), ND, NCP )
           if (coil(icoil)%type == 1) call CurvDeriv1( icoil, t1CU(idof+1:idof+ND), ND, NF )
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
          case_curv, curv_alpha, curv_c, k0, MPI_COMM_FOCUS,coil_type_spline,Splines

  implicit none
  include "mpif.h"

  INTEGER, intent(in)  :: icoil
  REAL   , intent(out) :: curvRet 
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  INTEGER              :: kseg, NS
  REAL,allocatable     :: curvv(:)

  NS = coil(icoil)%NS 
  SALLOCATE(curvv, (0:NS),zero)

  FATAL( CurvDeriv0, icoil .lt. 1 .or. icoil .gt. Ncoils, icoil not in right range )
 
  curvv = zero
  curvRet = zero

  curvv = sqrt( (coil(icoil)%za*coil(icoil)%yt-coil(icoil)%zt*coil(icoil)%ya)**2 &
             + (coil(icoil)%xa*coil(icoil)%zt-coil(icoil)%xt*coil(icoil)%za)**2  & 
             + (coil(icoil)%ya*coil(icoil)%xt-coil(icoil)%yt*coil(icoil)%xa)**2 )& 
             / ((coil(icoil)%xt)**2+(coil(icoil)%yt)**2+(coil(icoil)%zt)**2)**(1.5)
  coil(icoil)%maxcurv = maxval(curvv)

  if( case_curv == 1 ) then ! linear
     curvRet = sum(curvv)-curvv(0)
     if (coil(icoil)%type==1) curvRet = pi2*curvRet/NS
     if (coil(icoil)%type==coil_type_spline) curvRet = 1.0*curvRet/NS	
  elseif( case_curv == 2 ) then ! quadratic
     curvv = curvv*curvv
     curvRet = sum(curvv)-curvv(0)
     if (coil(icoil)%type==1) curvRet = pi2*curvRet/NS
     if (coil(icoil)%type==coil_type_spline) curvRet = 1.0*curvRet/NS	
  elseif( case_curv == 3 ) then ! penalty method 
     if( curv_alpha < 2.0 ) then
          FATAL( CurvDeriv0, .true. , curv_alpha must be 2 or greater )
     endif
     do kseg = 0,NS-1
        if( curvv(kseg) > k0 ) then
           curvRet = curvRet + ( curvv(kseg) - k0 )**curv_alpha
        endif
     enddo
     if (coil(icoil)%type==1) curvRet = pi2*curvRet/NS
     if (coil(icoil)%type==coil_type_spline) curvRet = 1.0*curvRet/NS	
  elseif( case_curv == 4 ) then ! penalty plus linear 
     if( curv_alpha < 2.0 ) then 
          FATAL( CurvDeriv0, .true. , curv_alpha must be 2 or greater )
     endif
     do kseg = 0,NS-1
        if( curvv(kseg) > k0 ) then
           curvRet = curvRet + sqrt(coil(icoil)%xt(kseg)**2+coil(icoil)%yt(kseg)**2+coil(icoil)%zt(kseg)**2)*( (curvv(kseg) - k0 )**curv_alpha + curv_c*curvv(kseg) )
        else
           curvRet = curvRet + sqrt(coil(icoil)%xt(kseg)**2+coil(icoil)%yt(kseg)**2+coil(icoil)%zt(kseg)**2)*curv_c*curvv(kseg)
        endif
     enddo
     call lenDeriv0( icoil, coil(icoil)%L )
     if (coil(icoil)%type==1) curvRet = pi2*curvRet/(NS*coil(icoil)%L)
     if (coil(icoil)%type==coil_type_spline) curvRet = 1.0*curvRet/(NS*coil(icoil)%L)	
  else   
     FATAL( CurvDeriv0, .true. , invalid case_curv option ) 
  endif

  return

end subroutine CurvDeriv0

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine CurvDeriv1(icoil, derivs, ND, NC) !Calculate all derivatives for a coil

  use globals, only: dp, zero, pi2, coil, DoF, myid, ounit, Ncoils, &
                case_curv, curv_alpha, k0, curv_c, FouCoil, MPI_COMM_FOCUS,Splines,coil_type_spline
  implicit none
  include "mpif.h"

  INTEGER, intent(in)  :: icoil, ND , NC  !NC is actually NCP for splines and NF for FouCoil
  REAL   , intent(out) :: derivs(1:1, 1:ND)
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  INTEGER              :: kseg, astat, ierr, doff, nff, NS
  REAL                 :: dl3, xt, yt, zt, xa, ya, za, f1, f2, ncc, nss, ncp, nsp, curvHold, penCurv, curvH, leng
  REAL, dimension(1:1, 0:coil(icoil)%NS-1) :: dLx, dLy, dLz
  REAL                 :: d1L(1:1, 1:ND)

  FATAL( CurvDeriv1, icoil .lt. 1 .or. icoil .gt. Ncoils, icoil not in right range )
  
  derivs = zero

  derivs(1,1:ND) = 0.0
  NS = coil(icoil)%NS
  d1L = zero
  call lenDeriv0( icoil, coil(icoil)%L )
  leng = coil(icoil)%L
  call lenDeriv1( icoil, d1L(1:1,1:ND), ND )
  
  do kseg = 0,NS-1

     xt = coil(icoil)%xt(kseg) ; yt = coil(icoil)%yt(kseg) ; zt = coil(icoil)%zt(kseg)
     xa = coil(icoil)%xa(kseg) ; ya = coil(icoil)%ya(kseg) ; za = coil(icoil)%za(kseg)
     f1 = sqrt( (xt*ya-xa*yt)**2 + (xt*za-xa*zt)**2 + (yt*za-ya*zt)**2 );
     f2 = (xt**2+yt**2+zt**2)**(1.5);
     select case (coil(icoil)%type)
	case(1)
           do nff = 1,NC
        	ncc = -1.0*nff    *FouCoil(icoil)%smt(kseg,nff)
        	ncp = -1.0*nff*nff*FouCoil(icoil)%cmt(kseg,nff)
        	nss =      nff    *FouCoil(icoil)%cmt(kseg,nff)
        	nsp = -1.0*nff*nff*FouCoil(icoil)%smt(kseg,nff)

  	  	if( case_curv == 1 ) then
			!if (kseg < 10 .or. kseg > 250) write(ounit,'(I7.3 "seg" I7.3 " CP" 6F20.10)')kseg,nff,derivs(1,nff)        	        		
			! Xc 
           		derivs(1,1+nff)      = derivs(1,1+nff)      + -1.0*(f1/(f2**2))*( 3.0*xt*sqrt(xt**2+yt**2+zt**2)*ncc) &
           		+ ((f1**(-1))*((xt*ya-xa*yt)*(ncc*ya-ncp*yt)+(xt*za-xa*zt)*(ncc*za-ncp*zt)) )/f2; 
           		! Xs
           		derivs(1,1+nff+NC)   = derivs(1,1+nff+NC)   + -1.0*(f1/(f2**2))*( 3.0*xt*sqrt(xt**2+yt**2+zt**2)*nss) &
          	 	+ ((f1**(-1))*((xt*ya-xa*yt)*(nss*ya-nsp*yt)+(xt*za-xa*zt)*(nss*za-nsp*zt)) )/f2;
           		! Yc
           		derivs(1,2+nff+2*NC) = derivs(1,2+nff+2*NC) + -1.0*(f1/(f2**2))*( 3.0*yt*sqrt(xt**2+yt**2+zt**2)*ncc) &
           		+ ((f1**(-1))*((xt*ya-xa*yt)*(ncp*xt-ncc*xa)+(yt*za-ya*zt)*(ncc*za-ncp*zt)) )/f2;  
           		! Ys
           		derivs(1,2+nff+3*NC) = derivs(1,2+nff+3*NC) + -1.0*(f1/(f2**2))*( 3.0*yt*sqrt(xt**2+yt**2+zt**2)*nss) &
           		+ ((f1**(-1))*((xt*ya-xa*yt)*(nsp*xt-nss*xa)+(yt*za-ya*zt)*(nss*za-nsp*zt)) )/f2;
           		! Zc
           		derivs(1,3+nff+4*NC) = derivs(1,3+nff+4*NC) + -1.0*(f1/(f2**2))*( 3.0*zt*sqrt(xt**2+yt**2+zt**2)*ncc) &
           		+ ((f1**(-1))*((xt*za-xa*zt)*(xt*ncp-xa*ncc)+(yt*za-ya*zt)*(yt*ncp-ya*ncc)) )/f2;
           		! Zs
           		derivs(1,3+nff+5*NC) = derivs(1,3+nff+5*NC) + -1.0*(f1/(f2**2))*( 3.0*zt*sqrt(xt**2+yt**2+zt**2)*nss) &
           		+ ((f1**(-1))*((xt*za-xa*zt)*(xt*nsp-xa*nss)+(yt*za-ya*zt)*(yt*nsp-ya*nss)) )/f2;

	        elseif( case_curv == 2 ) then
        		curvHold = sqrt( (za*yt-zt*ya)**2 + (xa*zt-xt*za)**2 + (ya*xt-yt*xa)**2 ) / ((xt)**2+(yt)**2+(zt)**2)**(1.5);
           		! Xc 
           		derivs(1,1+nff)      = derivs(1,1+nff)      + (-1.0*(f1/(f2**2))*( 3.0*xt*sqrt(xt**2+yt**2+zt**2)*ncc) &
          		 + ((f1**(-1))*((xt*ya-xa*yt)*(ncc*ya-ncp*yt)+(xt*za-xa*zt)*(ncc*za-ncp*zt)) )/f2)*2.0*curvHold; 
        		! Xs
           		derivs(1,1+nff+NC)   = derivs(1,1+nff+NC)   + (-1.0*(f1/(f2**2))*( 3.0*xt*sqrt(xt**2+yt**2+zt**2)*nss) &
           		+ ((f1**(-1))*((xt*ya-xa*yt)*(nss*ya-nsp*yt)+(xt*za-xa*zt)*(nss*za-nsp*zt)) )/f2)*2.0*curvHold;
           		! Yc
           		derivs(1,2+nff+2*NC) = derivs(1,2+nff+2*NC) + (-1.0*(f1/(f2**2))*( 3.0*yt*sqrt(xt**2+yt**2+zt**2)*ncc) &
           		+ ((f1**(-1))*((xt*ya-xa*yt)*(ncp*xt-ncc*xa)+(yt*za-ya*zt)*(ncc*za-ncp*zt)) )/f2)*2.0*curvHold;  
           		! Ys
           		derivs(1,2+nff+3*NC) = derivs(1,2+nff+3*NC) + (-1.0*(f1/(f2**2))*( 3.0*yt*sqrt(xt**2+yt**2+zt**2)*nss) &
           		+ ((f1**(-1))*((xt*ya-xa*yt)*(nsp*xt-nss*xa)+(yt*za-ya*zt)*(nss*za-nsp*zt)) )/f2)*2.0*curvHold;
           		! Zc
           		derivs(1,3+nff+4*NC) = derivs(1,3+nff+4*NC) + (-1.0*(f1/(f2**2))*( 3.0*zt*sqrt(xt**2+yt**2+zt**2)*ncc) &
           		+ ((f1**(-1))*((xt*za-xa*zt)*(xt*ncp-xa*ncc)+(yt*za-ya*zt)*(yt*ncp-ya*ncc)) )/f2)*2.0*curvHold;
           		! Zs
           		derivs(1,3+nff+5*NC) = derivs(1,3+nff+5*NC) + (-1.0*(f1/(f2**2))*( 3.0*zt*sqrt(xt**2+yt**2+zt**2)*nss) &
           		+ ((f1**(-1))*((xt*za-xa*zt)*(xt*nsp-xa*nss)+(yt*za-ya*zt)*(yt*nsp-ya*nss)) )/f2)*2.0*curvHold;

        	elseif( case_curv == 3 ) then
           		curvHold = sqrt( (za*yt-zt*ya)**2 + (xa*zt-xt*za)**2 + (ya*xt-yt*xa)**2 ) / ((xt)**2+(yt)**2+(zt)**2)**(1.5);
           		penCurv = curv_alpha*((curvHold-k0)**(curv_alpha-1.0));
           		if( curvHold > k0 ) then
              		! Xc 
              		derivs(1,1+nff)      = derivs(1,1+nff)      + (-1.0*(f1/(f2**2))*( 3.0*xt*sqrt(xt**2+yt**2+zt**2)*ncc) &
              		+ ((f1**(-1))*((xt*ya-xa*yt)*(ncc*ya-ncp*yt)+(xt*za-xa*zt)*(ncc*za-ncp*zt)) )/f2)*penCurv; 
              		! Xs
              		derivs(1,1+nff+NC)   = derivs(1,1+nff+NC)   + (-1.0*(f1/(f2**2))*( 3.0*xt*sqrt(xt**2+yt**2+zt**2)*nss) &
              		+ ((f1**(-1))*((xt*ya-xa*yt)*(nss*ya-nsp*yt)+(xt*za-xa*zt)*(nss*za-nsp*zt)) )/f2)*penCurv;
              		! Yc
              		derivs(1,2+nff+2*NC) = derivs(1,2+nff+2*NC) + (-1.0*(f1/(f2**2))*( 3.0*yt*sqrt(xt**2+yt**2+zt**2)*ncc) &
              		+ ((f1**(-1))*((xt*ya-xa*yt)*(ncp*xt-ncc*xa)+(yt*za-ya*zt)*(ncc*za-ncp*zt)) )/f2)*penCurv;  
              		! Ys
              		derivs(1,2+nff+3*NC) = derivs(1,2+nff+3*NC) + (-1.0*(f1/(f2**2))*( 3.0*yt*sqrt(xt**2+yt**2+zt**2)*nss) &
              		+ ((f1**(-1))*((xt*ya-xa*yt)*(nsp*xt-nss*xa)+(yt*za-ya*zt)*(nss*za-nsp*zt)) )/f2)*penCurv;
              		! Zc
              		derivs(1,3+nff+4*NC) = derivs(1,3+nff+4*NC) + (-1.0*(f1/(f2**2))*( 3.0*zt*sqrt(xt**2+yt**2+zt**2)*ncc) &
              		+ ((f1**(-1))*((xt*za-xa*zt)*(xt*ncp-xa*ncc)+(yt*za-ya*zt)*(yt*ncp-ya*ncc)) )/f2)*penCurv;
              		! Zs
              		derivs(1,3+nff+5*NC) = derivs(1,3+nff+5*NC) + (-1.0*(f1/(f2**2))*( 3.0*zt*sqrt(xt**2+yt**2+zt**2)*nss) &
              		+ ((f1**(-1))*((xt*za-xa*zt)*(xt*nsp-xa*nss)+(yt*za-ya*zt)*(yt*nsp-ya*nss)) )/f2)*penCurv;
           	else
              		derivs(1,1+nff) = derivs(1,1+nff);
              		derivs(1,1+nff+NC) = derivs(1,1+nff+NC);
              		derivs(1,2+nff+2*NC) = derivs(1,2+nff+2*NC);
              		derivs(1,2+nff+3*NC) = derivs(1,2+nff+3*NC);
              		derivs(1,3+nff+4*NC) = derivs(1,3+nff+4*NC);
              		derivs(1,3+nff+5*NC) = derivs(1,3+nff+5*NC);
           	endif

        	elseif( case_curv == 4 ) then
           		curvHold = sqrt( (za*yt-zt*ya)**2 + (xa*zt-xt*za)**2 + (ya*xt-yt*xa)**2 ) / ((xt)**2+(yt)**2+(zt)**2)**(1.5);
           		if( curvHold > k0 ) then
              			curvH = 1.0
           		else 
              			curvH = 0.0
           		endif
           		! Xc
           		derivs(1,1+nff)     = derivs(1,1+nff)       + sqrt(xt**2+yt**2+zt**2) &
           		*(curvH*( curvHold - k0 )**curv_alpha + curv_c*curvHold)*-1.0*d1L(   1,1+nff)/leng
           		derivs(1,1+nff)     = derivs(1,1+nff)       + xt*ncc*(xt**2+yt**2+zt**2)**(-.5) &
           		*(curvH*( curvHold - k0 )**curv_alpha + curv_c*curvHold)
           		derivs(1,1+nff)     = derivs(1,1+nff)       + sqrt(xt**2+yt**2+zt**2) &
           		*(-1.0*(f1/(f2**2))*( 3.0*xt*sqrt(xt**2+yt**2+zt**2)*ncc) + ((f1**(-1))*((xt*ya-xa*yt)*(ncc*ya-ncp*yt)+(xt*za-xa*zt)*(ncc*za-ncp*zt)) )/f2) &
           		*(curv_alpha*curvH*( curvHold - k0 )**(curv_alpha - 1.0) + curv_c)
           		! Xs
           		derivs(1,1+nff+NC)  = derivs(1,1+nff+NC)    + sqrt(xt**2+yt**2+zt**2) &
           		*(curvH*( curvHold - k0 )**curv_alpha + curv_c*curvHold)*-1.0*d1L(1,1+nff+NC)/leng
           		derivs(1,1+nff+NC)  = derivs(1,1+nff+NC)    + xt*nss*(xt**2+yt**2+zt**2)**(-.5) &
           		*(curvH*( curvHold - k0 )**curv_alpha + curv_c*curvHold)
           		derivs(1,1+nff+NC)  = derivs(1,1+nff+NC)    + sqrt(xt**2+yt**2+zt**2) &
           		*(-1.0*(f1/(f2**2))*( 3.0*xt*sqrt(xt**2+yt**2+zt**2)*nss) + ((f1**(-1))*((xt*ya-xa*yt)*(nss*ya-nsp*yt)+(xt*za-xa*zt)*(nss*za-nsp*zt)) )/f2) &
           		*(curv_alpha*curvH*( curvHold - k0 )**(curv_alpha - 1.0) + curv_c)
           		! Yc
           		derivs(1,2+nff+2*NC) = derivs(1,2+nff+2*NC) + sqrt(xt**2+yt**2+zt**2) &
           		*(curvH*( curvHold - k0 )**curv_alpha + curv_c*curvHold)*-1.0*d1L(1,2+nff+2*NC)/leng
           		derivs(1,2+nff+2*NC) = derivs(1,2+nff+2*NC) + yt*ncc*(xt**2+yt**2+zt**2)**(-.5) &
           		*(curvH*( curvHold - k0 )**curv_alpha + curv_c*curvHold)
           		derivs(1,2+nff+2*NC) = derivs(1,2+nff+2*NC) + sqrt(xt**2+yt**2+zt**2) &
           		*(-1.0*(f1/(f2**2))*( 3.0*yt*sqrt(xt**2+yt**2+zt**2)*ncc) + ((f1**(-1))*((xt*ya-xa*yt)*(ncp*xt-ncc*xa)+(yt*za-ya*zt)*(ncc*za-ncp*zt)) )/f2) &
           		*(curv_alpha*curvH*( curvHold - k0 )**(curv_alpha - 1.0) + curv_c)
           		! Ys
           		derivs(1,2+nff+3*NC) = derivs(1,2+nff+3*NC) + sqrt(xt**2+yt**2+zt**2) &
           		*(curvH*( curvHold - k0 )**curv_alpha + curv_c*curvHold)*-1.0*d1L(1,2+nff+3*NC)/leng
           		derivs(1,2+nff+3*NC) = derivs(1,2+nff+3*NC) + yt*nss*(xt**2+yt**2+zt**2)**(-.5) &
           		*(curvH*( curvHold - k0 )**curv_alpha + curv_c*curvHold)
           		derivs(1,2+nff+3*NC) = derivs(1,2+nff+3*NC) + sqrt(xt**2+yt**2+zt**2) &
           		*(-1.0*(f1/(f2**2))*( 3.0*yt*sqrt(xt**2+yt**2+zt**2)*nss) + ((f1**(-1))*((xt*ya-xa*yt)*(nsp*xt-nss*xa)+(yt*za-ya*zt)*(nss*za-nsp*zt)) )/f2) &
           		*(curv_alpha*curvH*( curvHold - k0 )**(curv_alpha - 1.0) + curv_c)
           		! Zc
           		derivs(1,3+nff+4*NC) = derivs(1,3+nff+4*NC) + sqrt(xt**2+yt**2+zt**2) &
           		*(curvH*( curvHold - k0 )**curv_alpha + curv_c*curvHold)*-1.0*d1L(1,3+nff+4*NC)/leng
           		derivs(1,3+nff+4*NC) = derivs(1,3+nff+4*NC) + zt*ncc*(xt**2+yt**2+zt**2)**(-.5) &
           		*(curvH*( curvHold - k0 )**curv_alpha + curv_c*curvHold)
           		derivs(1,3+nff+4*NC) = derivs(1,3+nff+4*NC) + sqrt(xt**2+yt**2+zt**2) &
           		*(-1.0*(f1/(f2**2))*( 3.0*zt*sqrt(xt**2+yt**2+zt**2)*ncc) + ((f1**(-1))*((xt*za-xa*zt)*(xt*ncp-xa*ncc)+(yt*za-ya*zt)*(yt*ncp-ya*ncc)) )/f2) &
           		*(curv_alpha*curvH*( curvHold - k0 )**(curv_alpha - 1.0) + curv_c)
           		! Zs
           		derivs(1,3+nff+5*NC) = derivs(1,3+nff+5*NC) + sqrt(xt**2+yt**2+zt**2) &
           		*(curvH*( curvHold - k0 )**curv_alpha + curv_c*curvHold)*-1.0*d1L(1,3+nff+5*NC)/leng
           		derivs(1,3+nff+5*NC) = derivs(1,3+nff+5*NC) + zt*nss*(xt**2+yt**2+zt**2)**(-.5) &                
           		*(curvH*( curvHold - k0 )**curv_alpha + curv_c*curvHold)
           		derivs(1,3+nff+5*NC) = derivs(1,3+nff+5*NC) + sqrt(xt**2+yt**2+zt**2) &
           		*(-1.0*(f1/(f2**2))*( 3.0*zt*sqrt(xt**2+yt**2+zt**2)*nss) + ((f1**(-1))*((xt*za-xa*zt)*(xt*nsp-xa*nss)+(yt*za-ya*zt)*(yt*nsp-ya*nss)) )/f2) &
          		 *(curv_alpha*curvH*( curvHold - k0 )**(curv_alpha - 1.0) + curv_c)
        	else
           		FATAL( CurvDeriv1, .true. , invalid case_curv option )  
        	endif
       	  enddo
	case(coil_type_spline)
	    do nff = 1,NC
        	ncc = Splines(icoil)%db_dt(kseg,nff-1)
        	ncp = Splines(icoil)%db_dt_2(kseg,nff-1)

  	  	if( case_curv == 1 ) then

			! X 
           		derivs(1,nff)      = derivs(1,nff)      + -1.0*(f1/(f2**2))*( 3.0*xt*sqrt(xt**2+yt**2+zt**2)*ncc) &
           		+ ((f1**(-1))*((xt*ya-xa*yt)*(ncc*ya-ncp*yt)+(xt*za-xa*zt)*(ncc*za-ncp*zt)) )/f2; 
			! Y
           		derivs(1,nff+NC) = derivs(1,nff+NC) + -1.0*(f1/(f2**2))*( 3.0*yt*sqrt(xt**2+yt**2+zt**2)*ncc) &
           		+ ((f1**(-1))*((xt*ya-xa*yt)*(ncp*xt-ncc*xa)+(yt*za-ya*zt)*(ncc*za-ncp*zt)) )/f2;
			! Z
           		derivs(1,nff+2*NC) = derivs(1,nff+2*NC) + -1.0*(f1/(f2**2))*( 3.0*zt*sqrt(xt**2+yt**2+zt**2)*ncc) &
           		+ ((f1**(-1))*((xt*za-xa*zt)*(xt*ncp-xa*ncc)+(yt*za-ya*zt)*(yt*ncp-ya*ncc)) )/f2;


	        elseif( case_curv == 2 ) then
        		curvHold = sqrt( (za*yt-zt*ya)**2 + (xa*zt-xt*za)**2 + (ya*xt-yt*xa)**2 ) / ((xt)**2+(yt)**2+(zt)**2)**(1.5);
           		! X 
           		derivs(1,nff)      = derivs(1,nff)      + (-1.0*(f1/(f2**2))*( 3.0*xt*sqrt(xt**2+yt**2+zt**2)*ncc) &
          		 + ((f1**(-1))*((xt*ya-xa*yt)*(ncc*ya-ncp*yt)+(xt*za-xa*zt)*(ncc*za-ncp*zt)) )/f2)*2.0*curvHold; 
			! Y
           		derivs(1,nff+NC) = derivs(1,nff+NC) + (-1.0*(f1/(f2**2))*( 3.0*yt*sqrt(xt**2+yt**2+zt**2)*ncc) &
           		+ ((f1**(-1))*((xt*ya-xa*yt)*(ncp*xt-ncc*xa)+(yt*za-ya*zt)*(ncc*za-ncp*zt)) )/f2)*2.0*curvHold; 
			! Z
           		derivs(1,nff+2*NC) = derivs(1,nff+2*NC) + (-1.0*(f1/(f2**2))*( 3.0*zt*sqrt(xt**2+yt**2+zt**2)*ncc) &
           		+ ((f1**(-1))*((xt*za-xa*zt)*(xt*ncp-xa*ncc)+(yt*za-ya*zt)*(yt*ncp-ya*ncc)) )/f2)*2.0*curvHold; 
        	elseif( case_curv == 3 ) then
           		curvHold = sqrt( (za*yt-zt*ya)**2 + (xa*zt-xt*za)**2 + (ya*xt-yt*xa)**2 ) / ((xt)**2+(yt)**2+(zt)**2)**(1.5);
           		penCurv = curv_alpha*((curvHold-k0)**(curv_alpha-1.0));
           		if( curvHold > k0 ) then
              			! X 
              			derivs(1,nff)      = derivs(1,nff)      + (-1.0*(f1/(f2**2))*( 3.0*xt*sqrt(xt**2+yt**2+zt**2)*ncc) &
              			+ ((f1**(-1))*((xt*ya-xa*yt)*(ncc*ya-ncp*yt)+(xt*za-xa*zt)*(ncc*za-ncp*zt)) )/f2)*penCurv; 
				! Y
				derivs(1,nff + NC) = derivs(1,nff+NC) + (-1.0*(f1/(f2**2))*( 3.0*yt*sqrt(xt**2+yt**2+zt**2)*ncc) &
              			+ ((f1**(-1))*((xt*ya-xa*yt)*(ncp*xt-ncc*xa)+(yt*za-ya*zt)*(ncc*za-ncp*zt)) )/f2)*penCurv; 
				! Z
              			derivs(1,nff+2*NC) = derivs(1,nff+2*NC) + (-1.0*(f1/(f2**2))*( 3.0*zt*sqrt(xt**2+yt**2+zt**2)*ncc) &
              			+ ((f1**(-1))*((xt*za-xa*zt)*(xt*ncp-xa*ncc)+(yt*za-ya*zt)*(yt*ncp-ya*ncc)) )/f2)*penCurv;
           		else
              			derivs(1,nff) = derivs(1,nff);
				derivs(1,nff + NC) = derivs(1,nff + NC);
				derivs(1,nff + 2*NC) = derivs(1,nff + 2*NC);
           		endif

        	elseif( case_curv == 4 ) then
           		curvHold = sqrt( (za*yt-zt*ya)**2 + (xa*zt-xt*za)**2 + (ya*xt-yt*xa)**2 ) / ((xt)**2+(yt)**2+(zt)**2)**(1.5);
           		if( curvHold > k0 ) then
              			curvH = 1.0
           		else 
              			curvH = 0.0
           		endif
           		! X
           		derivs(1,nff)     = derivs(1,nff)       + sqrt(xt**2+yt**2+zt**2) &
           		*(curvH*( curvHold - k0 )**curv_alpha + curv_c*curvHold)*-1.0*d1L(   1,nff)/leng
           		derivs(1,nff)     = derivs(1,nff)       + xt*ncc*(xt**2+yt**2+zt**2)**(-.5) &
           		*(curvH*( curvHold - k0 )**curv_alpha + curv_c*curvHold)
           		derivs(1,nff)     = derivs(1,nff)       + sqrt(xt**2+yt**2+zt**2) &
           		*(-1.0*(f1/(f2**2))*( 3.0*xt*sqrt(xt**2+yt**2+zt**2)*ncc) + ((f1**(-1))*((xt*ya-xa*yt)*(ncc*ya-ncp*yt)+(xt*za-xa*zt)*(ncc*za-ncp*zt)) )/f2) &
           		*(curv_alpha*curvH*( curvHold - k0 )**(curv_alpha - 1.0) + curv_c)
           		! Y
           		derivs(1,nff+NC) = derivs(1,nff+NC) + sqrt(xt**2+yt**2+zt**2) &
           		*(curvH*( curvHold - k0 )**curv_alpha + curv_c*curvHold)*-1.0*d1L(1,nff+NC)/leng
           		derivs(1,nff+NC) = derivs(1,nff+NC) + yt*ncc*(xt**2+yt**2+zt**2)**(-.5) &
           		*(curvH*( curvHold - k0 )**curv_alpha + curv_c*curvHold)
           		derivs(1,nff+NC) = derivs(1,nff+NC) + sqrt(xt**2+yt**2+zt**2) &
           		*(-1.0*(f1/(f2**2))*( 3.0*yt*sqrt(xt**2+yt**2+zt**2)*ncc) + ((f1**(-1))*((xt*ya-xa*yt)*(ncp*xt-ncc*xa)+(yt*za-ya*zt)*(ncc*za-ncp*zt)) )/f2) &
           		*(curv_alpha*curvH*( curvHold - k0 )**(curv_alpha - 1.0) + curv_c)
           		! Z
           		derivs(1,nff+2*NC) = derivs(1,nff+2*NC) + sqrt(xt**2+yt**2+zt**2) &
           		*(curvH*( curvHold - k0 )**curv_alpha + curv_c*curvHold)*-1.0*d1L(1,nff+2*NC)/leng
           		derivs(1,nff+2*NC) = derivs(1,nff+2*NC) + zt*ncc*(xt**2+yt**2+zt**2)**(-.5) &
           		*(curvH*( curvHold - k0 )**curv_alpha + curv_c*curvHold)
           		derivs(1,nff+2*NC) = derivs(1,nff+2*NC) + sqrt(xt**2+yt**2+zt**2) &
           		*(-1.0*(f1/(f2**2))*( 3.0*zt*sqrt(xt**2+yt**2+zt**2)*ncc) + ((f1**(-1))*((xt*za-xa*zt)*(xt*ncp-xa*ncc)+(yt*za-ya*zt)*(yt*ncp-ya*ncc)) )/f2) &
           		*(curv_alpha*curvH*( curvHold - k0 )**(curv_alpha - 1.0) + curv_c)
        	else
           		FATAL( CurvDeriv1, .true. , invalid case_curv option )  
        	endif
       	  enddo
	case default 
		FATAL( CurvDeriv1, .true. , invalid coil_type option )	
	end select
  enddo
  if (coil(icoil)%type == 1)  derivs = pi2*derivs/NS
  if (coil(icoil)%type == coil_type_spline )derivs = derivs/(NS)
  if( case_curv == 4 ) derivs = derivs/leng

  return

end subroutine CurvDeriv1


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
