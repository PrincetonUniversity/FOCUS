
!title (curvature) ! Calculate objective functon for straight-out coils and its derivatives. (tkruger)

!latex \briefly{The cost function will allow for straight-out coils, improving 
!latex         This function is still under development. emph{targt\_length}.}

!latex \calledby{\link{solvers}}

!latex  \section{General}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! curv is total penalty 
! chi = chi + weight_curv*curv
! t1CU is total derivative of penalty
! LM implemented
! not parallelized, at some point check to see how long takes to run
subroutine straight(ideriv)
  use globals, only: dp, zero, half, pi2, machprec, ncpu, myid, ounit, MPI_COMM_FOCUS, &
       coil, DoF, Ncoils, Nfixgeo, Ndof, str, t1Str, t2Str, weight_straight, FouCoil, &
       mstr, istr, LM_fvec, LM_fjac,coil_type_spline,Splines,origin_surface_x,origin_surface_y,origin_surface_z,coeff_disp_straight

  implicit none
  include "mpif.h"
  INTEGER, INTENT(in) :: ideriv

  INTEGER             :: astat, ierr, icoil, idof, ND, NF,NCP, ivec
  REAL                :: strAdd

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  str = zero
  strAdd = zero
  ivec = 1



  if( ideriv >= 0 ) then

     do icoil = 1, Ncoils
        if( (coil(icoil)%type .ne. 1) .AND. (coil(icoil)%type .ne. coil_type_spline) ) exit ! only for Fourier
        if( coil(icoil)%Lc     /=  0 ) then ! if geometry is free

           call StrDeriv0(icoil,strAdd)
           str = str + strAdd

           if (mstr > 0) then ! L-M format of targets
              LM_fvec(istr+ivec) = weight_straight*strAdd
              ivec = ivec + 1
           endif
        endif 
     enddo

     str = str / (Ncoils - Nfixgeo + machprec)

  endif

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if ( ideriv >= 1 ) then

     t1Str = zero
     ivec = 1
     idof = 0
     icoil = 1	
     !call unpacking(xdof)
     !if (myid == 0 .AND. icoil==1) write(ounit,'("R**2  " 7F20.10)')(coil(icoil)%xx-origin_surface_x)**2 + &
     !								(coil(icoil)%yy-origin_surface_y)**2
     do icoil = 1, Ncoils

	if( (coil(icoil)%type .ne. 1) .AND. (coil(icoil)%type .ne. coil_type_spline) ) exit ! only for Fourier       

        ND = DoF(icoil)%ND
        if (coil(icoil)%type == 1) NF = FouCoil(icoil)%NF
        if (coil(icoil)%type == coil_type_spline) NCP = Splines(icoil)%NCP

        if ( coil(icoil)%Ic /= 0 ) then !if current is free;
           idof = idof + 1
        endif

        if ( coil(icoil)%Lc /= 0 ) then !if geometry is free;
           if (coil(icoil)%type == coil_type_spline) call StrDeriv1( icoil, t1Str(idof+1:idof+ND), ND, NCP )
           if (coil(icoil)%type == 1) call StrDeriv1( icoil, t1Str(idof+1:idof+ND), ND, NF )
           if (mstr > 0) then ! L-M format of targets
              LM_fjac(istr+ivec, idof+1:idof+ND) = weight_straight * t1Str(idof+1:idof+ND)
              ivec = ivec + 1
           endif
           idof = idof + ND
        endif

     enddo
     FATAL( straight-out , idof .ne. Ndof, counting error in packing )

     t1Str = t1Str / (Ncoils - Nfixgeo + machprec)

  endif

  return
end subroutine straight

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
! calculate coil curvature

subroutine StrDeriv0(icoil,strRet)

  use globals, only: dp, zero, pi2, ncpu, astat, ierr, myid, ounit, coil, NFcoil, Nseg, Ncoils, &
          case_straight, straight_alpha, str_c, str_k0, MPI_COMM_FOCUS,coil_type_spline,Splines, &
	  origin_surface_x, origin_surface_y, origin_surface_z,coeff_disp_straight


  implicit none
  include "mpif.h"

  INTEGER, intent(in)  :: icoil
  REAL   , intent(out) :: strRet 
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  INTEGER              :: kseg, NS
  REAL,allocatable     :: strv(:)
  REAL                 :: mean_xy_distance, dispersion

  NS = coil(icoil)%NS 

  SALLOCATE(strv,   (0:NS-1),zero)	

  FATAL( StrDeriv0, icoil .lt. 1 .or. icoil .gt. Ncoils, icoil not in right range )

  strv = zero
  strRet = zero
  mean_xy_distance = SUM((coil(icoil)%xx-origin_surface_x)**2 + (coil(icoil)%yy-origin_surface_y)**2)/(coil(icoil)%NS)
  !dispersion = (MAXVAL((coil(icoil)%xx-origin_surface_x)**2 + (coil(icoil)%yy-origin_surface_y)**2)- &
!	       MINVAL((coil(icoil)%xx-origin_surface_x)**2 + (coil(icoil)%yy-origin_surface_y)**2))/2
  dispersion = (MAXVAL((coil(icoil)%xx-origin_surface_x)**2 + (coil(icoil)%yy-origin_surface_y)**2) - mean_xy_distance)/2
      do kseg = 0,coil(icoil)%NS-1
       if ((coil(icoil)%xx(kseg)-origin_surface_x)**2 + (coil(icoil)%yy(kseg)-origin_surface_y)**2 > &
	   (mean_xy_distance + coeff_disp_straight*dispersion)) then
           strv(kseg) = sqrt( (coil(icoil)%za(kseg)*coil(icoil)%yt(kseg)-coil(icoil)%zt(kseg)*coil(icoil)%ya(kseg))**2 &
                            + (coil(icoil)%xa(kseg)*coil(icoil)%zt(kseg)-coil(icoil)%xt(kseg)*coil(icoil)%za(kseg))**2  & 
                            + (coil(icoil)%ya(kseg)*coil(icoil)%xt(kseg)-coil(icoil)%yt(kseg)*coil(icoil)%xa(kseg))**2 )& 
                            / ((coil(icoil)%xt(kseg))**2+(coil(icoil)%yt(kseg))**2+(coil(icoil)%zt(kseg))**2)**(1.5)
	else 
		strv(kseg) = 0
	endif
      enddo
  !if (myid==0 .AND. icoil==1) write(ounit,'("strv "7F20.10)')strv

  if( case_straight == 1 ) then ! linear

     strRet = sum(strv)
     if (coil(icoil)%type==1) strRet = pi2*strRet/NS
     if (coil(icoil)%type==coil_type_spline) strRet = 1.0*strRet/NS	
  elseif( case_straight == 2 ) then ! quadratic
     strv = strv*strv
     strRet = sum(strv)
     if (coil(icoil)%type==1) strRet = pi2*strRet/NS
     if (coil(icoil)%type==coil_type_spline) strRet = 1.0*strRet/NS	
  elseif( case_straight == 3 ) then ! penalty method 
     if( straight_alpha < 2.0 ) then
          FATAL( StrDeriv0, .true. , straight_alpha must be 2 or greater )
     endif
     do kseg = 0,coil(icoil)%NS-1
	if ((coil(icoil)%xx(kseg)-origin_surface_x)**2 + (coil(icoil)%yy(kseg)-origin_surface_y)**2 > &
	   (mean_xy_distance + coeff_disp_straight*dispersion)) then
           if( strv(kseg) > str_k0 ) then
              strRet = strRet + ( strv(kseg) - str_k0 )**straight_alpha
           endif
	endif
     enddo
     if (coil(icoil)%type==1) strRet = pi2*strRet/NS
     if (coil(icoil)%type==coil_type_spline) strRet = 1.0*strRet/NS	
  elseif( case_straight == 4 ) then ! penalty plus linear 
     if( straight_alpha < 2.0 ) then 
          FATAL( StrDeriv0, .true. , straight_alpha must be 2 or greater )
     endif

     do kseg = 0 ,coil(icoil)%NS-1
	if ((coil(icoil)%xx(kseg)-origin_surface_x)**2 + (coil(icoil)%yy(kseg)-origin_surface_y)**2 > &
	   (mean_xy_distance + coeff_disp_straight*dispersion)) then
           if( strv(kseg) > str_k0 ) then
              strRet = strRet + sqrt(coil(icoil)%xt(kseg)**2+coil(icoil)%yt(kseg)**2+coil(icoil)%zt(kseg)**2)*( (strv(kseg) - str_k0 )**straight_alpha &
		       + str_c*strv(kseg) )
           else
              strRet = strRet + sqrt(coil(icoil)%xt(kseg)**2+coil(icoil)%yt(kseg)**2+coil(icoil)%zt(kseg)**2)*str_c*strv(kseg)
           endif
	endif
     enddo

     call lenDeriv0( icoil, coil(icoil)%L )
     if (coil(icoil)%type==1) strRet = pi2*strRet/NS
     if (coil(icoil)%type==coil_type_spline) strRet = 1.0*strRet/NS	

  else   
     FATAL( StrDeriv0, .true. , invalid case_straight option ) 
  endif

  return

end subroutine StrDeriv0

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
subroutine StrDeriv1(icoil, derivs, ND, NC ) !Calculate all derivatives for a coil

  use globals, only: dp, zero, pi2, coil, DoF, myid, ounit, Ncoils, &
                case_straight, straight_alpha, str_k0, str_c, FouCoil, MPI_COMM_FOCUS,Splines,coil_type_spline, &
	        origin_surface_x, origin_surface_y, origin_surface_z,xdof,coeff_disp_straight
  implicit none
  !include "mpif.h"

  INTEGER, intent(in)  :: icoil, ND , NC  !NC is actually NCP for spline coils
  REAL   , intent(out) :: derivs(1:1, 1:ND)
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  INTEGER              :: kseg, astat, ierr, doff, nff, NS
  REAL                 :: dl3, xt, yt, zt, xa, ya, za, f1, f2, ncc, nss, ncp, nsp, strHold, penStr, strH, leng,mean_xy_distance,dispersion
  REAL, dimension(1:1, 0:coil(icoil)%NS-1) :: dLx, dLy, dLz
  REAL                 :: d1L(1:1, 1:ND)

  FATAL( StrDeriv1, icoil .lt. 1 .or. icoil .gt. Ncoils, icoil not in right range )

  derivs = zero
  mean_xy_distance = SUM((coil(icoil)%xx-origin_surface_x)**2 + (coil(icoil)%yy-origin_surface_y)**2)/(coil(icoil)%NS)
  !dispersion = (MAXVAL((coil(icoil)%xx-origin_surface_x)**2 + (coil(icoil)%yy-origin_surface_y)**2)- &
!	       MINVAL((coil(icoil)%xx-origin_surface_x)**2 + (coil(icoil)%yy-origin_surface_y)**2))/2

  dispersion = (MAXVAL((coil(icoil)%xx-origin_surface_x)**2 + (coil(icoil)%yy-origin_surface_y)**2) - mean_xy_distance)/2

 !write(ounit, '( " dist " F20.10 " disp " F20.10)') &
!						     mean_xy_distance,dispersion

  derivs(1,1:ND) = 0.0
  NS = coil(icoil)%NS
  d1L = zero
  call lenDeriv0( icoil, coil(icoil)%L )
  leng = coil(icoil)%L
  call lenDeriv1( icoil, d1L(1:1,1:ND), ND )  
     do kseg = 0,coil(icoil)%NS-1

        if (((coil(icoil)%xx(kseg)-origin_surface_x)**2 + (coil(icoil)%yy(kseg)-origin_surface_y)**2) > &
	   (mean_xy_distance + coeff_disp_straight*dispersion)) then

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

   	     	    if( case_straight == 1 ) then        	        		
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

	            elseif( case_straight == 2 ) then
        		 strHold = sqrt( (za*yt-zt*ya)**2 + (xa*zt-xt*za)**2 + (ya*xt-yt*xa)**2 ) / ((xt)**2+(yt)**2+(zt)**2)**(1.5);
           		 ! Xc 
           		 derivs(1,1+nff)      = derivs(1,1+nff)      + (-1.0*(f1/(f2**2))*( 3.0*xt*sqrt(xt**2+yt**2+zt**2)*ncc) &
          		  + ((f1**(-1))*((xt*ya-xa*yt)*(ncc*ya-ncp*yt)+(xt*za-xa*zt)*(ncc*za-ncp*zt)) )/f2)*2.0*strHold; 
        		 ! Xs
           		 derivs(1,1+nff+NC)   = derivs(1,1+nff+NC)   + (-1.0*(f1/(f2**2))*( 3.0*xt*sqrt(xt**2+yt**2+zt**2)*nss) &
           		 + ((f1**(-1))*((xt*ya-xa*yt)*(nss*ya-nsp*yt)+(xt*za-xa*zt)*(nss*za-nsp*zt)) )/f2)*2.0*strHold;
           		 ! Yc
           		 derivs(1,2+nff+2*NC) = derivs(1,2+nff+2*NC) + (-1.0*(f1/(f2**2))*( 3.0*yt*sqrt(xt**2+yt**2+zt**2)*ncc) &
           		 + ((f1**(-1))*((xt*ya-xa*yt)*(ncp*xt-ncc*xa)+(yt*za-ya*zt)*(ncc*za-ncp*zt)) )/f2)*2.0*strHold;  
           		 ! Ys
           		 derivs(1,2+nff+3*NC) = derivs(1,2+nff+3*NC) + (-1.0*(f1/(f2**2))*( 3.0*yt*sqrt(xt**2+yt**2+zt**2)*nss) &
           		 + ((f1**(-1))*((xt*ya-xa*yt)*(nsp*xt-nss*xa)+(yt*za-ya*zt)*(nss*za-nsp*zt)) )/f2)*2.0*strHold;
           		 ! Zc
           		 derivs(1,3+nff+4*NC) = derivs(1,3+nff+4*NC) + (-1.0*(f1/(f2**2))*( 3.0*zt*sqrt(xt**2+yt**2+zt**2)*ncc) &
           		 + ((f1**(-1))*((xt*za-xa*zt)*(xt*ncp-xa*ncc)+(yt*za-ya*zt)*(yt*ncp-ya*ncc)) )/f2)*2.0*strHold;
           		 ! Zs
           		 derivs(1,3+nff+5*NC) = derivs(1,3+nff+5*NC) + (-1.0*(f1/(f2**2))*( 3.0*zt*sqrt(xt**2+yt**2+zt**2)*nss) &
           		 + ((f1**(-1))*((xt*za-xa*zt)*(xt*nsp-xa*nss)+(yt*za-ya*zt)*(yt*nsp-ya*nss)) )/f2)*2.0*strHold;

        	    elseif( case_straight == 3 ) then
           	       strHold = sqrt( (za*yt-zt*ya)**2 + (xa*zt-xt*za)**2 + (ya*xt-yt*xa)**2 ) / ((xt)**2+(yt)**2+(zt)**2)**(1.5);
           	       penStr = straight_alpha*((strHold-str_k0)**(straight_alpha-1.0));
                       if( strHold > str_k0 ) then
              	         ! Xc 
              		 derivs(1,1+nff)      = derivs(1,1+nff)      + (-1.0*(f1/(f2**2))*( 3.0*xt*sqrt(xt**2+yt**2+zt**2)*ncc) &
              		 + ((f1**(-1))*((xt*ya-xa*yt)*(ncc*ya-ncp*yt)+(xt*za-xa*zt)*(ncc*za-ncp*zt)) )/f2)*penStr; 
              		 ! Xs
              		 derivs(1,1+nff+NC)   = derivs(1,1+nff+NC)   + (-1.0*(f1/(f2**2))*( 3.0*xt*sqrt(xt**2+yt**2+zt**2)*nss) &
              		 + ((f1**(-1))*((xt*ya-xa*yt)*(nss*ya-nsp*yt)+(xt*za-xa*zt)*(nss*za-nsp*zt)) )/f2)*penStr;
              		 ! Yc
              		 derivs(1,2+nff+2*NC) = derivs(1,2+nff+2*NC) + (-1.0*(f1/(f2**2))*( 3.0*yt*sqrt(xt**2+yt**2+zt**2)*ncc) &
              		 + ((f1**(-1))*((xt*ya-xa*yt)*(ncp*xt-ncc*xa)+(yt*za-ya*zt)*(ncc*za-ncp*zt)) )/f2)*penStr;  
              		 ! Ys
              		 derivs(1,2+nff+3*NC) = derivs(1,2+nff+3*NC) + (-1.0*(f1/(f2**2))*( 3.0*yt*sqrt(xt**2+yt**2+zt**2)*nss) &
              		 + ((f1**(-1))*((xt*ya-xa*yt)*(nsp*xt-nss*xa)+(yt*za-ya*zt)*(nss*za-nsp*zt)) )/f2)*penStr;
              		 ! Zc
              		 derivs(1,3+nff+4*NC) = derivs(1,3+nff+4*NC) + (-1.0*(f1/(f2**2))*( 3.0*zt*sqrt(xt**2+yt**2+zt**2)*ncc) &
              		 + ((f1**(-1))*((xt*za-xa*zt)*(xt*ncp-xa*ncc)+(yt*za-ya*zt)*(yt*ncp-ya*ncc)) )/f2)*penStr;
              		 ! Zs
              		 derivs(1,3+nff+5*NC) = derivs(1,3+nff+5*NC) + (-1.0*(f1/(f2**2))*( 3.0*zt*sqrt(xt**2+yt**2+zt**2)*nss) &
              		 + ((f1**(-1))*((xt*za-xa*zt)*(xt*nsp-xa*nss)+(yt*za-ya*zt)*(yt*nsp-ya*nss)) )/f2)*penStr;
           	       else
              		 derivs(1,1+nff) = derivs(1,1+nff);
              		 derivs(1,1+nff+NC) = derivs(1,1+nff+NC);
              		 derivs(1,2+nff+2*NC) = derivs(1,2+nff+2*NC);
              		 derivs(1,2+nff+3*NC) = derivs(1,2+nff+3*NC);
              		 derivs(1,3+nff+4*NC) = derivs(1,3+nff+4*NC);
              		 derivs(1,3+nff+5*NC) = derivs(1,3+nff+5*NC);
           	       endif

        	    elseif( case_straight == 4 ) then
           	        strHold = sqrt( (za*yt-zt*ya)**2 + (xa*zt-xt*za)**2 + (ya*xt-yt*xa)**2 ) / ((xt)**2+(yt)**2+(zt)**2)**(1.5);
           	        if( strHold > str_k0 ) then
              			strH = 1.0
           	        else 
              			strH = 0.0
           	        endif
           		! Xc
           		derivs(1,1+nff)     = derivs(1,1+nff)       + sqrt(xt**2+yt**2+zt**2) &
           		*(strH*( strHold - str_k0 )**straight_alpha + str_c*strHold)*-1.0*d1L(   1,1+nff)/leng
           		derivs(1,1+nff)     = derivs(1,1+nff)       + xt*ncc*(xt**2+yt**2+zt**2)**(-.5) &
           		*(strH*( strHold - str_k0 )**straight_alpha + str_c*strHold)
           		derivs(1,1+nff)     = derivs(1,1+nff)       + sqrt(xt**2+yt**2+zt**2) &
           		*(-1.0*(f1/(f2**2))*( 3.0*xt*sqrt(xt**2+yt**2+zt**2)*ncc) + ((f1**(-1))*((xt*ya-xa*yt)*(ncc*ya-ncp*yt)+(xt*za-xa*zt)* &
			(ncc*za-ncp*zt)) )/f2) &
           		*(straight_alpha*strH*( strHold - str_k0 )**(straight_alpha - 1.0) + str_c)
           		! Xs
           		derivs(1,1+nff+NC)  = derivs(1,1+nff+NC)    + sqrt(xt**2+yt**2+zt**2) &
           		*(strH*( strHold - str_k0 )**straight_alpha + str_c*strHold)*-1.0*d1L(1,1+nff+NC)/leng
           		derivs(1,1+nff+NC)  = derivs(1,1+nff+NC)    + xt*nss*(xt**2+yt**2+zt**2)**(-.5) &
           		*(strH*( strHold - str_k0 )**straight_alpha + str_c*strHold)
           		derivs(1,1+nff+NC)  = derivs(1,1+nff+NC)    + sqrt(xt**2+yt**2+zt**2) &
           		*(-1.0*(f1/(f2**2))*( 3.0*xt*sqrt(xt**2+yt**2+zt**2)*nss) + ((f1**(-1))*((xt*ya-xa*yt)*(nss*ya-nsp*yt)+(xt*za-xa*zt)*(nss*za-nsp*zt)) )/f2) &
           		*(straight_alpha*strH*( strHold - str_k0 )**(straight_alpha - 1.0) + str_c)
           		! Yc
           		derivs(1,2+nff+2*NC) = derivs(1,2+nff+2*NC) + sqrt(xt**2+yt**2+zt**2) &
           		*(strH*( strHold - str_k0 )**straight_alpha + str_c*strHold)*-1.0*d1L(1,2+nff+2*NC)/leng
           		derivs(1,2+nff+2*NC) = derivs(1,2+nff+2*NC) + yt*ncc*(xt**2+yt**2+zt**2)**(-.5) &
           		*(strH*( strHold - str_k0 )**straight_alpha + str_c*strHold)
           		derivs(1,2+nff+2*NC) = derivs(1,2+nff+2*NC) + sqrt(xt**2+yt**2+zt**2) &
           		*(-1.0*(f1/(f2**2))*( 3.0*yt*sqrt(xt**2+yt**2+zt**2)*ncc) + ((f1**(-1))*((xt*ya-xa*yt)*(ncp*xt-ncc*xa)+(yt*za-ya*zt)*(ncc*za-ncp*zt)) )/f2) &
           		*(straight_alpha*strH*( strHold - str_k0 )**(straight_alpha - 1.0) + str_c)
           		! Ys
           		derivs(1,2+nff+3*NC) = derivs(1,2+nff+3*NC) + sqrt(xt**2+yt**2+zt**2) &
           		*(strH*( strHold - str_k0 )**straight_alpha + str_c*strHold)*-1.0*d1L(1,2+nff+3*NC)/leng
           		derivs(1,2+nff+3*NC) = derivs(1,2+nff+3*NC) + yt*nss*(xt**2+yt**2+zt**2)**(-.5) &
           		*(strH*( strHold - str_k0 )**straight_alpha + str_c*strHold)
           		derivs(1,2+nff+3*NC) = derivs(1,2+nff+3*NC) + sqrt(xt**2+yt**2+zt**2) &
           		*(-1.0*(f1/(f2**2))*( 3.0*yt*sqrt(xt**2+yt**2+zt**2)*nss) + ((f1**(-1))*((xt*ya-xa*yt)*(nsp*xt-nss*xa)+(yt*za-ya*zt)*(nss*za-nsp*zt)) )/f2) &
           		*(straight_alpha*strH*( strHold - str_k0 )**(straight_alpha - 1.0) + str_c)
           		! Zc
           		derivs(1,3+nff+4*NC) = derivs(1,3+nff+4*NC) + sqrt(xt**2+yt**2+zt**2) &
           		*(strH*( strHold - str_k0 )**straight_alpha + str_c*strHold)*-1.0*d1L(1,3+nff+4*NC)/leng
           		derivs(1,3+nff+4*NC) = derivs(1,3+nff+4*NC) + zt*ncc*(xt**2+yt**2+zt**2)**(-.5) &
           		*(strH*( strHold - str_k0 )**straight_alpha + str_c*strHold)
           		derivs(1,3+nff+4*NC) = derivs(1,3+nff+4*NC) + sqrt(xt**2+yt**2+zt**2) &
           		*(-1.0*(f1/(f2**2))*( 3.0*zt*sqrt(xt**2+yt**2+zt**2)*ncc) + ((f1**(-1))*((xt*za-xa*zt)*(xt*ncp-xa*ncc)+(yt*za-ya*zt)*(yt*ncp-ya*ncc)) )/f2) &
           		*(straight_alpha*strH*( strHold - str_k0 )**(straight_alpha - 1.0) + str_c)
           		! Zs
           		derivs(1,3+nff+5*NC) = derivs(1,3+nff+5*NC) + sqrt(xt**2+yt**2+zt**2) &
           		*(strH*( strHold - str_k0 )**straight_alpha + str_c*strHold)*-1.0*d1L(1,3+nff+5*NC)/leng
           		derivs(1,3+nff+5*NC) = derivs(1,3+nff+5*NC) + zt*nss*(xt**2+yt**2+zt**2)**(-.5) &                
           		*(strH*( strHold - str_k0 )**straight_alpha + str_c*strHold)
           		derivs(1,3+nff+5*NC) = derivs(1,3+nff+5*NC) + sqrt(xt**2+yt**2+zt**2) &
           		*(-1.0*(f1/(f2**2))*( 3.0*zt*sqrt(xt**2+yt**2+zt**2)*nss) + ((f1**(-1))*((xt*za-xa*zt)*(xt*nsp-xa*nss)+(yt*za-ya*zt)*(yt*nsp-ya*nss)) )/f2) &
          		 *(straight_alpha*strH*( strHold - str_k0 )**(straight_alpha - 1.0) + str_c)
        	    else
           		FATAL( StrDeriv1, .true. , invalid case_straight option )  
        	    endif
       	         enddo
	      case(coil_type_spline)

	         do nff = 1,NC
        	    ncc = Splines(icoil)%db_dt(kseg,nff-1)
        	    ncp = Splines(icoil)%db_dt_2(kseg,nff-1)

  	  	    if( case_straight == 1 ) then
			! X 
           		derivs(1,nff)      = derivs(1,nff)      + -1.0*(f1/(f2**2))*( 3.0*xt*sqrt(xt**2+yt**2+zt**2)*ncc) &
           		+ ((f1**(-1))*((xt*ya-xa*yt)*(ncc*ya-ncp*yt)+(xt*za-xa*zt)*(ncc*za-ncp*zt)) )/f2; 
			! Y
           		derivs(1,nff+NC) = derivs(1,nff+NC) + -1.0*(f1/(f2**2))*( 3.0*yt*sqrt(xt**2+yt**2+zt**2)*ncc) &
           		+ ((f1**(-1))*((xt*ya-xa*yt)*(ncp*xt-ncc*xa)+(yt*za-ya*zt)*(ncc*za-ncp*zt)) )/f2;
			! Z
           		derivs(1,nff+2*NC) = derivs(1,nff+2*NC) + -1.0*(f1/(f2**2))*( 3.0*zt*sqrt(xt**2+yt**2+zt**2)*ncc) &
           		+ ((f1**(-1))*((xt*za-xa*zt)*(xt*ncp-xa*ncc)+(yt*za-ya*zt)*(yt*ncp-ya*ncc)) )/f2;
  			!if (myid==0 .AND. icoil==1) write(ounit,'("derivs iter " I10 "  " 6F20.10)')kseg,derivs
			!if (myid==0 .AND. icoil==1) write(ounit,'("derivs iter " I10.3 )')kseg

	            elseif( case_straight == 2 ) then
        		strHold = sqrt( (za*yt-zt*ya)**2 + (xa*zt-xt*za)**2 + (ya*xt-yt*xa)**2 ) / ((xt)**2+(yt)**2+(zt)**2)**(1.5);
           		! X 
           		derivs(1,nff)      = derivs(1,nff)      + (-1.0*(f1/(f2**2))*( 3.0*xt*sqrt(xt**2+yt**2+zt**2)*ncc) &
          		 + ((f1**(-1))*((xt*ya-xa*yt)*(ncc*ya-ncp*yt)+(xt*za-xa*zt)*(ncc*za-ncp*zt)) )/f2)*2.0*strHold; 
			! Y
           		derivs(1,nff+NC) = derivs(1,nff+NC) + (-1.0*(f1/(f2**2))*( 3.0*yt*sqrt(xt**2+yt**2+zt**2)*ncc) &
           		+ ((f1**(-1))*((xt*ya-xa*yt)*(ncp*xt-ncc*xa)+(yt*za-ya*zt)*(ncc*za-ncp*zt)) )/f2)*2.0*strHold; 
			! Z
           		derivs(1,nff+2*NC) = derivs(1,nff+2*NC) + (-1.0*(f1/(f2**2))*( 3.0*zt*sqrt(xt**2+yt**2+zt**2)*ncc) &
           		+ ((f1**(-1))*((xt*za-xa*zt)*(xt*ncp-xa*ncc)+(yt*za-ya*zt)*(yt*ncp-ya*ncc)) )/f2)*2.0*strHold; 
        	    elseif( case_straight == 3 ) then
           		strHold = sqrt( (za*yt-zt*ya)**2 + (xa*zt-xt*za)**2 + (ya*xt-yt*xa)**2 ) / ((xt)**2+(yt)**2+(zt)**2)**(1.5);
           		penStr = straight_alpha*((strHold-str_k0)**(straight_alpha-1.0));
           		if( strHold > str_k0 ) then
              			! X 
              			derivs(1,nff)      = derivs(1,nff)      + (-1.0*(f1/(f2**2))*( 3.0*xt*sqrt(xt**2+yt**2+zt**2)*ncc) &
              			+ ((f1**(-1))*((xt*ya-xa*yt)*(ncc*ya-ncp*yt)+(xt*za-xa*zt)*(ncc*za-ncp*zt)) )/f2)*penStr; 
				! Y
				derivs(1,nff + NC) = derivs(1,nff+NC) + (-1.0*(f1/(f2**2))*( 3.0*yt*sqrt(xt**2+yt**2+zt**2)*ncc) &
              			+ ((f1**(-1))*((xt*ya-xa*yt)*(ncp*xt-ncc*xa)+(yt*za-ya*zt)*(ncc*za-ncp*zt)) )/f2)*penStr; 
				! Z
              			derivs(1,nff+2*NC) = derivs(1,nff+2*NC) + (-1.0*(f1/(f2**2))*( 3.0*zt*sqrt(xt**2+yt**2+zt**2)*ncc) &
              			+ ((f1**(-1))*((xt*za-xa*zt)*(xt*ncp-xa*ncc)+(yt*za-ya*zt)*(yt*ncp-ya*ncc)) )/f2)*penStr;
           		else
              			derivs(1,nff) = derivs(1,nff);
				derivs(1,nff + NC) = derivs(1,nff + NC);
				derivs(1,nff + 2*NC) = derivs(1,nff + 2*NC);
           		endif

        	    elseif( case_straight == 4 ) then
           		strHold = sqrt( (za*yt-zt*ya)**2 + (xa*zt-xt*za)**2 + (ya*xt-yt*xa)**2 ) / ((xt)**2+(yt)**2+(zt)**2)**(1.5);
           		if( strHold > str_k0 ) then
              			strH = 1.0
           		else 
              			strH = 0.0
           		endif
           		! X
           		derivs(1,nff)     = derivs(1,nff)       + sqrt(xt**2+yt**2+zt**2) &
           		*(strH*( strHold - str_k0 )**straight_alpha + str_c*strHold)*-1.0*d1L(   1,nff)/leng
           		derivs(1,nff)     = derivs(1,nff)       + xt*ncc*(xt**2+yt**2+zt**2)**(-.5) &
           		*(strH*( strHold - str_k0 )**straight_alpha + str_c*strHold)
           		derivs(1,nff)     = derivs(1,nff)       + sqrt(xt**2+yt**2+zt**2) &
           		*(-1.0*(f1/(f2**2))*( 3.0*xt*sqrt(xt**2+yt**2+zt**2)*ncc) + ((f1**(-1))*((xt*ya-xa*yt)*(ncc*ya-ncp*yt)+(xt*za-xa*zt)*(ncc*za-ncp*zt)) )/f2) &
           		*(straight_alpha*strH*( strHold - str_k0 )**(straight_alpha - 1.0) + str_c)
           		! Y
           		derivs(1,nff+NC) = derivs(1,nff+NC) + sqrt(xt**2+yt**2+zt**2) &
           		*(strH*( strHold - str_k0 )**straight_alpha + str_c*strHold)*-1.0*d1L(1,nff+NC)/leng
           		derivs(1,nff+NC) = derivs(1,nff+NC) + yt*ncc*(xt**2+yt**2+zt**2)**(-.5) &
           		*(strH*( strHold - str_k0 )**straight_alpha + str_c*strHold)
           		derivs(1,nff+NC) = derivs(1,nff+NC) + sqrt(xt**2+yt**2+zt**2) &
           		*(-1.0*(f1/(f2**2))*( 3.0*yt*sqrt(xt**2+yt**2+zt**2)*ncc) + ((f1**(-1))*((xt*ya-xa*yt)*(ncp*xt-ncc*xa)+(yt*za-ya*zt)*(ncc*za-ncp*zt)) )/f2) &
           		*(straight_alpha*strH*( strHold - str_k0 )**(straight_alpha - 1.0) + str_c)
           		! Z
           		derivs(1,nff+2*NC) = derivs(1,nff+2*NC) + sqrt(xt**2+yt**2+zt**2) &
           		*(strH*( strHold - str_k0 )**straight_alpha + str_c*strHold)*-1.0*d1L(1,nff+2*NC)/leng
           		derivs(1,nff+2*NC) = derivs(1,nff+2*NC) + zt*ncc*(xt**2+yt**2+zt**2)**(-.5) &
           		*(strH*( strHold - str_k0 )**straight_alpha + str_c*strHold)
           		derivs(1,nff+2*NC) = derivs(1,nff+2*NC) + sqrt(xt**2+yt**2+zt**2) &
           		*(-1.0*(f1/(f2**2))*( 3.0*zt*sqrt(xt**2+yt**2+zt**2)*ncc) + ((f1**(-1))*((xt*za-xa*zt)*(xt*ncp-xa*ncc)+(yt*za-ya*zt)*(yt*ncp-ya*ncc)) )/f2) &
           		*(straight_alpha*strH*( strHold - str_k0 )**(straight_alpha - 1.0) + str_c)
        	    else
           		FATAL( StrDeriv1, .true. , invalid case_straight option )  
           	    endif
       	         enddo
	      case default 
		 FATAL( StrDeriv1, .true. , invalid coil_type option )	
	      end select	 
	endif
     enddo

  if (coil(icoil)%type == 1)  derivs = pi2*derivs/NS
  if (coil(icoil)%type == coil_type_spline )derivs = derivs/(NS)
  if( case_straight == 4 ) derivs = derivs/leng

 !if (myid==0 .AND. icoil==1) write(ounit,'(" str derivs " 7F20.10)')derivs
  !if (myid==0 .AND. icoil==3) write(ounit,'(" 3 derivs " 7F20.10)')derivs

  return

end subroutine StrDeriv1

