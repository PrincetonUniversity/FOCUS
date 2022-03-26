
!title (Straight out) ! Calculate objective functon for straight-out coils and its derivatives. (tkruger)

!latex \briefly{The cost function will allow for straight-out coils.}

!latex \calledby{\link{solvers}}

!latex    \section{General}
!latex    Computes the straight out penalty in one of three possible forms. It is similar to the curvature penality but only applied to the outer part of the coil.
!latex    \item[1.]  A linear objective function 
!latex    \begin{eqnarray}
!latex    f_{curv} = \frac{1}{N_c}\sum_{i=1}^{N_c} \int_0^{2\pi}W(t)\kappa_i dt
!latex    \end{eqnarray}
!latex        \item A quadratic objective function 
!latex    \begin{eqnarray}
!latex    f_{curv} = \frac{1}{N_c}\sum_{i=1}^{N_c} \int_0^{2\pi}W(t)\kappa_i^2 dt
!latex    \end{eqnarray}
!latex           \item A maximum curvature objective function 
!latex    \begin{eqnarray}
!latex    f_{curv} = \frac{1}{N_c}\sum_{i=1}^{N_c} \int_0^{2\pi} W(t) H_{\kappa_o}(\kappa_i) (\kappa_i - \kappa_o)^\alpha dt
!latex    \end{eqnarray}
    
!latex    where $H_{\kappa_o}(\kappa_i)$ is the step function,  $\kappa_o$ and $\alpha$ are user-defined parameters and $\kappa$ is the curvature at a point defined as 
!latex    \begin{equation}
!latex    \kappa = \frac{( (z^\prime^\prime y^\prime - y^\prime^\prime z^\prime)^2 + (x^\prime^\prime z^\prime - z^\prime^\prime x^\prime)^2 + (y^\prime^\prime x^\prime - x^\prime^\prime y^\prime)^2)^\frac{1}{2}}{(x^\prime^2 + y^\prime^2 + z^\prime^2 )^\frac{3}{2}}
!latex    \end{equation}    
!latex     In particular $\kappa_o$ has the meaning of a maximum allowed curvature and the cost function has no effect on points with a curvature smaller than this parameter.

!latex     W(t) is a weight function to enure the penalty is only applied to the outer part of the coil
!latex     \begin{equation}
!latex     W(t) = \begin{cases} 1, \hspace{2em} P_{xy}(t) > P_m + \beta P_d \\ 0, \hspace{2em}  otherwise \end{cases}
!latex     where $P_m$ is the mean squared distance of the coil projection in the x-y plane from a user-defined point and $P_d$
!latex 	   is a measure for half the maximum projection in the x-y plane of the coil radius.
!latex	   Considering the distance of the projection of a point along the coil in the x-y plane from the user-defined point ($x_0,y_0$)
!latex     \begin{equation}
!latex     P_{xy}(t) = (x_i(t) - x_0)^2 + (y_i(t) - y_0)^2 \ ,
!latex     \end{equation}
!latex     they are defined as 
!latex     \begin{align}
!latex      & P_m = \overline{P_{xy}} \ , & P_d = \frac{\max(P_{xy}) - P_m}{2} \ .
!latex      \end{align}
!latex      $\beta$ is a user-defined parameter that can be used to tweak the section of the coil affected by the optimization.
!latex      With a value $\beta = 0$ the section of the coil with a projected distance greater than the mean value for that coil will all be affected 
!latex      while for a value $\beta = 2$ the new cost function will not be applied to any point. 

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! str is total penalty 
! chi = chi + weight_straight*str
! t1Str is total derivative of penalty
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
          case_straight, str_alpha, str_k0, str_k1, str_beta, str_gamma, str_sigma, penfun_str, &
          str_k1len, MPI_COMM_FOCUS,coil_type_spline,Splines, origin_surface_x, origin_surface_y, origin_surface_z,coeff_disp_straight


  implicit none
  include "mpif.h"

  INTEGER, intent(in)  :: icoil
  REAL   , intent(out) :: strRet 
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  INTEGER              :: kseg, NS
  REAL,allocatable     :: strv(:), xt(:), yt(:), zt(:), xa(:), ya(:), za(:)
  REAL                 :: hypc, curv_hold, k1_use,mean_xy_distance, dispersion

  NS = coil(icoil)%NS 

  SALLOCATE(strv,   (0:NS-1),zero)
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

  if ( case_straight .eq. 1 ) then
     str_alpha = 0.0
     str_sigma = 1.0
     str_gamma = 1.0
     str_k1 = 0.0
  elseif ( case_straight .eq. 2 ) then
     str_alpha = 0.0
     str_sigma = 1.0
     str_gamma = 2.0
     str_k1 = 0.0
  elseif (case_straight .eq. 3 ) then
     str_sigma = 0.0
  else 
     ! Do nothing 
  endif

  FATAL( StrDeriv0, icoil .lt. 1 .or. icoil .gt. Ncoils, icoil not in right range )
  FATAL( StrDeriv0, penfun_str .ne. 1 .and. penfun_str .ne. 2 , invalid choice of penfun_str, pick 1 or 2 )
  FATAL( StrDeriv0, str_k0 < 0.0 , str_k0 cannot be negative )
  FATAL( StrDeriv0, str_k1 < 0.0 , str_k1 cannot be negative )
  FATAL( StrDeriv0, str_alpha < 0.0 , str_alpha cannot be negative )
  FATAL( StrDeriv0, str_beta < 2.0 , str_beta needs to be >= 2 )
  FATAL( StrDeriv0, str_gamma < 1.0 , str_gamma needs to be >= 1 )
  FATAL( StrDeriv0, str_sigma < 0.0 , str_sigma cannot be negative )
  FATAL( StrDeriv0, str_gamma .eq. 1.0 .and. str_k1 .ne. 0.0 , if str_gamma = 1, str_k1 must = 0 )

  strv = zero
  strRet = zero
  mean_xy_distance = SUM((coil(icoil)%xx-origin_surface_x)**2 + (coil(icoil)%yy-origin_surface_y)**2)/(coil(icoil)%NS)
  !dispersion = (MAXVAL((coil(icoil)%xx-origin_surface_x)**2 + (coil(icoil)%yy-origin_surface_y)**2)- &
!	       MINVAL((coil(icoil)%xx-origin_surface_x)**2 + (coil(icoil)%yy-origin_surface_y)**2))/2
  dispersion = (MAXVAL((coil(icoil)%xx-origin_surface_x)**2 + (coil(icoil)%yy-origin_surface_y)**2) - mean_xy_distance)/2
      do kseg = 0,coil(icoil)%NS-1
       if ((coil(icoil)%xx(kseg)-origin_surface_x)**2 + (coil(icoil)%yy(kseg)-origin_surface_y)**2 > &
	   (mean_xy_distance + coeff_disp_straight*dispersion)) then
           strv(kseg) = sqrt( (za(kseg)*yt(kseg)-zt(kseg)*ya(kseg))**2 + (xa(kseg)*zt(kseg)-xt(kseg)*za(kseg))**2 + (ya(kseg)*xt(kseg)-yt(kseg)*xa(kseg))**2 )& 
                            / ((xt(kseg))**2+(yt(kseg))**2+(zt(kseg))**2)**(1.5)
	else 
		strv(kseg) = 0
	endif
      enddo
  !if (myid==0 .AND. icoil==1) write(ounit,'("strv "7F20.10)')strv
  coil(icoil)%straight = strv 
  
  
  if ( str_k1len .eq. 1 ) then
     k1_use = pi2/coil(icoil)%Lo
  else
     k1_use = str_k1
  endif

  do kseg = 0,NS-1
     str_hold = 0.0
     if ( str_alpha .ne. 0.0 ) then
        if ( strv(kseg) > str_k0 ) then
           if ( penfun_str .eq. 1 ) then
              hypc = 0.5*exp( str_alpha*( strv(kseg) - str_k0 ) ) + 0.5*exp( -1.0*str_alpha*( strv(kseg) - str_k0 ) )
              str_hold = (hypc - 1.0)**2
           else
              str_hold = ( str_alpha*(strv(kseg)-str_k0) )**str_beta
           endif
        endif
     endif
     if ( str_sigma .ne. 0.0 ) then
        if ( strv(kseg) > k1_use ) then
           str_hold = str_hold + str_sigma*( ( strv(kseg) - k1_use )**str_gamma )
        endif

     endif
     strRet = strRet + str_hold*sqrt(xt(kseg)**2+yt(kseg)**2+zt(kseg)**2)
  enddo

  call lenDeriv0( icoil, coil(icoil)%L )
  
  if (coil(icoil)%type==1) strRet = pi2*strRet/(NS*coil(icoil)%L)
  if (coil(icoil)%type==coil_type_spline) strRet = 1.0*strRet/(NS*coil(icoil)%L)	
  
  DALLOCATE(strv)
  DALLOCATE(xt)
  DALLOCATE(yt)
  DALLOCATE(zt)
  DALLOCATE(xa)
  DALLOCATE(ya)
  DALLOCATE(za)
 
  return

end subroutine StrDeriv0

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
subroutine StrDeriv1(icoil, derivs, ND, NC ) !Calculate all derivatives for a coil

  use globals, only: dp, zero, pi2, coil, DoF, myid, ounit, Ncoils, &
                case_straight, straight_alpha, str_k0, str_k1, str_beta, str_gamma &
		,str_sigma, penfun_str, str_k1len, FouCoil, MPI_COMM_FOCUS,Splines,coil_type_spline, &
	        origin_surface_x, origin_surface_y, origin_surface_z,xdof,coeff_disp_straight
  implicit none
  !include "mpif.h"

  INTEGER, intent(in)  :: icoil, ND , NC  !NC is actually NCP for spline coils
  REAL   , intent(out) :: derivs(1:1, 1:ND)
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  INTEGER              :: kseg, astat, ierr, doff, nff, NS,n
  REAL                 :: dl3, xt, yt, zt, xa, ya, za, f1, f2, ncc, nss, ncp, nsp, strHold, penStr, strH, leng, hypc, hyps, str_deriv, &
                          k1_use, rtxrax, rtxray, rtxraz,mean_xy_distance,dispersion
  REAL, dimension(1:1, 0:coil(icoil)%NS-1) :: dLx, dLy, dLz
  REAL                 :: d1L(1:1, 1:ND)
  REAL,allocatable     :: dxtdDoF(:,:), dytdDoF(:,:), dztdDoF(:,:), dxadDoF(:,:), dyadDoF(:,:), dzadDoF(:,:)

  FATAL( StrDeriv1, icoil .lt. 1 .or. icoil .gt. Ncoils, icoil not in right range )

  mean_xy_distance = SUM((coil(icoil)%xx-origin_surface_x)**2 + (coil(icoil)%yy-origin_surface_y)**2)/(coil(icoil)%NS)
  !dispersion = (MAXVAL((coil(icoil)%xx-origin_surface_x)**2 + (coil(icoil)%yy-origin_surface_y)**2)- &
!	       MINVAL((coil(icoil)%xx-origin_surface_x)**2 + (coil(icoil)%yy-origin_surface_y)**2))/2

  dispersion = (MAXVAL((coil(icoil)%xx-origin_surface_x)**2 + (coil(icoil)%yy-origin_surface_y)**2) - mean_xy_distance)/2

 !write(ounit, '( " dist " F20.10 " disp " F20.10)') &
!						     mean_xy_distance,dispersion
  derivs = zero
  derivs(1,1:ND) = 0.0
  NS = coil(icoil)%NS
  d1L = zero
  
  SALLOCATE(dxtdDoF, (0:NS,1:ND), zero)
  SALLOCATE(dytdDoF, (0:NS,1:ND), zero)
  SALLOCATE(dztdDoF, (0:NS,1:ND), zero)
  SALLOCATE(dxadDoF, (0:NS,1:ND), zero)
  SALLOCATE(dyadDoF, (0:NS,1:ND), zero)
  SALLOCATE(dzadDoF, (0:NS,1:ND), zero)
  
  select case (coil(icoil)%type)    
 	 case(1)
     		do n = 1,NC
            	    	dxtdDof(0:NS,n+1)      = -1*FouCoil(icoil)%smt(0:NS,n) * n
     		    	dxtdDof(0:NS,n+NC+1)   =    FouCoil(icoil)%cmt(0:NS,n) * n
     	    		dytdDof(0:NS,n+2*NC+2) = -1*FouCoil(icoil)%smt(0:NS,n) * n
     	    		dytdDof(0:NS,n+3*NC+2) =    FouCoil(icoil)%cmt(0:NS,n) * n
     	    		dztdDof(0:NS,n+4*NC+3) = -1*FouCoil(icoil)%smt(0:NS,n) * n
     	    		dztdDof(0:NS,n+5*NC+3) =    FouCoil(icoil)%cmt(0:NS,n) * n

     	    		dxadDof(0:NS,n+1)      = -1*FouCoil(icoil)%cmt(0:NS,n) * n*n
     	    		dxadDof(0:NS,n+NC+1)   = -1*FouCoil(icoil)%smt(0:NS,n) * n*n
     	    		dyadDof(0:NS,n+2*NC+2) = -1*FouCoil(icoil)%cmt(0:NS,n) * n*n
     	    		dyadDof(0:NS,n+3*NC+2) = -1*FouCoil(icoil)%smt(0:NS,n) * n*n
     	    		dzadDof(0:NS,n+4*NC+3) = -1*FouCoil(icoil)%cmt(0:NS,n) * n*n
     	    		dzadDof(0:NS,n+5*NC+3) = -1*FouCoil(icoil)%smt(0:NS,n) * n*n
     		enddo
  	case(coil_type_spline) 
	  		dxtdDof(0:NS,0:ND)      =    Splines(icoil)%db_dt(0:NS,0:ND)
	  		dytdDof(0:NS,0:ND)      =    Splines(icoil)%db_dt(0:NS,0:ND)
	  		dztdDof(0:NS,0:ND)      =    Splines(icoil)%db_dt(0:NS,0:ND)
			
     	    		dxadDof(0:NS,0:ND)       =    Splines(icoil)%db_dt_2(0:NS,0:ND)
     	    		dyadDof(0:NS,0:ND)       =    Splines(icoil)%db_dt_2(0:NS,0:ND)
     	    		dzadDof(0:NS,0:ND)       =    Splines(icoil)%db_dt_2(0:NS,0:ND)
  	case default 
		FATAL( StrDeriv1, .true. , invalid coil_type option )	
  end select
  
  
  call lenDeriv0( icoil, coil(icoil)%L )
  leng = coil(icoil)%L
  call lenDeriv1( icoil, d1L(1:1,1:ND), ND )  
  
  if ( case_straight .eq. 1 ) then
     str_alpha = 0.0
     str_sigma = 1.0
     str_gamma = 1.0
     str_k1 = 0.0
  elseif ( case_straight .eq. 2 ) then
     str_alpha = 0.0
     str_sigma = 1.0
     str_gamma = 2.0
     str_k1 = 0.0
  elseif (case_straight .eq. 3 ) then
     str_sigma = 0.0
  else
     ! Do nothing
  endif
  
  do kseg = 0,NS-1
  
     if (((coil(icoil)%xx(kseg)-origin_surface_x)**2 + (coil(icoil)%yy(kseg)-origin_surface_y)**2) > &
	   (mean_xy_distance + coeff_disp_straight*dispersion)) then
	   
	     xt = coil(icoil)%xt(kseg) ; yt = coil(icoil)%yt(kseg) ; zt = coil(icoil)%zt(kseg) ;
	     xa = coil(icoil)%xa(kseg) ; ya = coil(icoil)%ya(kseg) ; za = coil(icoil)%za(kseg) ;
	     rtxrax = yt*za - zt*ya
	     rtxray = zt*xa - xt*za
	     rtxraz = xt*ya - yt*xa
	     f1 = sqrt( (xt*ya-xa*yt)**2 + (xt*za-xa*zt)**2 + (yt*za-ya*zt)**2 );
	     f2 = (xt**2+yt**2+zt**2)**(1.5);
	     strHold = f1/f2
	     penStr = 0.0
	     str_deriv = 0.0

	     if ( str_k1len .eq. 1 ) then
		k1_use = pi2/coil(icoil)%Lo
	     else
		k1_use = str_k1
	     endif
	     if ( penfun_str .eq. 1 ) then
		if ( strHold > str_k0 ) then
		   hypc = 0.5*exp( str_alpha*(strHold-str_k0) ) + 0.5*exp( -1.0*str_alpha*(strHold-str_k0) )
		   hyps = 0.5*exp( str_alpha*(strHold-str_k0) ) - 0.5*exp( -1.0*str_alpha*(strHold-str_k0) )
		   penStr = ( hypc - 1.0 )**2
		   str_deriv = 2.0*str_alpha*( hypc - 1.0 )*hyps
		endif
	     else
		if ( strHold > str_k0 ) then
		   strCurv = (str_alpha*(strHold-str_k0))**str_beta
		   str_deriv = str_beta*str_alpha*( (str_alpha*(strHold-str_k0))**(str_beta-1.0) )
		endif
	     endif
	     if ( strHold > k1_use ) then
		str_deriv = str_deriv + str_sigma*str_gamma*( (strHold-k1_use)**(str_gamma-1.0) )
		penStr = penStr + str_sigma*( (strHold-k1_use)**str_gamma )
	     endif

	     derivs(1,1:ND) = derivs(1,1:ND) + sqrt(xt**2+yt**2+zt**2)*penStr*-1.0*d1L(1,1:ND)/leng
	     derivs(1,1:ND) = derivs(1,1:ND) + (xt*dxtdDof(kseg,1:ND)+yt*dytdDof(kseg,1:ND)+zt*dztdDof(kseg,1:ND))*(xt**2+yt**2+zt**2)**(-.5)*penStr
	     derivs(1,1:ND) = derivs(1,1:ND) + sqrt(xt**2+yt**2+zt**2) &
	     *(-3.0*(f1/f2**2)*sqrt(xt**2+yt**2+zt**2)*(xt*dxtdDof(kseg,1:ND)+yt*dytdDof(kseg,1:ND)+zt*dztdDof(kseg,1:ND)) &
	     + ( rtxrax*(dytdDof(kseg,1:ND)*za - dztdDof(kseg,1:ND)*ya + yt*dzadDof(kseg,1:ND) - zt*dyadDof(kseg,1:ND)) &
	     +   rtxray*(dztdDof(kseg,1:ND)*xa - dxtdDof(kseg,1:ND)*za + zt*dxadDof(kseg,1:ND) - xt*dzadDof(kseg,1:ND)) &
	     +   rtxraz*(dxtdDof(kseg,1:ND)*ya - dytdDof(kseg,1:ND)*xa + xt*dyadDof(kseg,1:ND) - yt*dxadDof(kseg,1:ND)) )/(f1*f2) )*str_deriv
     endif
  enddo
  
  if (coil(icoil)%type == 1)  derivs = pi2*derivs/(NS*leng)
  if (coil(icoil)%type == coil_type_spline )derivs = derivs/(NS*leng)


 !if (myid==0 .AND. icoil==1) write(ounit,'(" str derivs " 7F20.10)')derivs
  !if (myid==0 .AND. icoil==3) write(ounit,'(" 3 derivs " 7F20.10)')derivs

  return

end subroutine StrDeriv1

