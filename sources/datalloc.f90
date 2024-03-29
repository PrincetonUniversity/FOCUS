!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine AllocData(type)
!------------------------------------------------------------------------------------------------------ 
! DATE:  04/05/2017
! Allocate data before using them, especially for those used several times;
! part can be : -1('dof'), 0('costfun0'), 1('costfun1')
!------------------------------------------------------------------------------------------------------
  use globals
  use bnorm_mod
  use mpi
  implicit none

  INTEGER, intent(in) :: type

  INTEGER             :: icoil, idof, ND, NF, NCP, icur, imag, isurf, NS, mm, iseg, i, iter_cp, n
  REAL                :: xtmp, mtmp, tt

  isurf = plasma

  !-------------------------------------------------------------------------------------------
  if (type == -1) then ! dof related data;
     Cdof = 0; Ndof = 0; Tdof = 0
     do icoil = 1, Ncoils 
        select case (coil(icoil)%type)       
        case(1)
           ! get number of DoF for each coil and allocate arrays;
           NS = coil(icoil)%NS
           NF = FouCoil(icoil)%NF
           ND = (6*NF + 3) ! total variables for geometry
           DoF(icoil)%ND = coil(icoil)%Lc * ND !# of DoF for icoil;
           SALLOCATE(DoF(icoil)%xdof, (1:DoF(icoil)%ND), zero)
           SALLOCATE(DoF(icoil)%xof , (0:coil(icoil)%NS-1, 1:ND), zero)
           SALLOCATE(DoF(icoil)%yof , (0:coil(icoil)%NS-1, 1:ND), zero)
           SALLOCATE(DoF(icoil)%zof , (0:coil(icoil)%NS-1, 1:ND), zero)
           SALLOCATE(DoF(icoil)%xtof, (0:coil(icoil)%NS, 1:ND), zero)
           SALLOCATE(DoF(icoil)%ytof, (0:coil(icoil)%NS, 1:ND), zero)
           SALLOCATE(DoF(icoil)%ztof, (0:coil(icoil)%NS, 1:ND), zero)
           ! allocate and calculate trignometric functions for re-use           
           SALLOCATE( FouCoil(icoil)%cmt, (0:NS, 0:NF), zero )
           SALLOCATE( FouCoil(icoil)%smt, (0:NS, 0:NF), zero )           
           do iseg = 0, NS
              tt = iseg * pi2 / NS
              do mm = 0, NF
                 FouCoil(icoil)%cmt(iseg,mm) = cos( mm * tt )
                 FouCoil(icoil)%smt(iseg,mm) = sin( mm * tt )
              enddo
           enddo

!!$           ip = (icoil-1)/Ncoils  ! the integer is the period number;
!!$           DoF(icoil)%xof(0:NS-1,      1:  NF+1) =  cosip(ip) * cmt(0:NS-1, 0:NF)  !x/xc
!!$           DoF(icoil)%xof(0:NS-1,   NF+2:2*NF+1) =  cosip(ip) * smt(0:NS-1, 1:NF)  !x/xs
!!$           DoF(icoil)%xof(0:NS-1, 2*NF+2:3*NF+2) = -sinip(ip) * cmt(0:NS-1, 0:NF)  !x/yc ; valid for ip>0 ;
!!$           DoF(icoil)%xof(0:NS-1, 3*NF+3:4*NF+2) = -sinip(ip) * smt(0:NS-1, 1:NF)  !x/ys ; valid for ip>0 ;
!!$           DoF(icoil)%yof(0:NS-1,      1:  NF+1) =  sinip(ip) * cmt(0:NS-1, 0:NF)  !y/xc ; valid for ip>0 ;
!!$           DoF(icoil)%yof(0:NS-1,   NF+2:2*NF+1) =  sinip(ip) * smt(0:NS-1, 1:NF)  !y/xs ; valid for ip>0 ;
!!$           DoF(icoil)%yof(0:NS-1, 2*NF+2:3*NF+2) =  cosip(ip) * cmt(0:NS-1, 0:NF)  !y/yc
!!$           DoF(icoil)%yof(0:NS-1, 3*NF+3:4*NF+2) =  cosip(ip) * smt(0:NS-1, 1:NF)  !y/ys
!!$           DoF(icoil)%zof(0:NS-1, 4*NF+3:5*NF+3) =              cmt(0:NS-1, 0:NF)  !z/zc
!!$           DoF(icoil)%zof(0:NS-1, 5*NF+4:6*NF+3) =              smt(0:NS-1, 1:NF)  !z/zs

           ! the derivatives of dx/dv 
           DoF(icoil)%xof(0:NS-1,      1:  NF+1) = FouCoil(icoil)%cmt(0:NS-1, 0:NF)  !x/xc
           DoF(icoil)%xof(0:NS-1,   NF+2:2*NF+1) = FouCoil(icoil)%smt(0:NS-1, 1:NF)  !x/xs
           !DoF(icoil)%xof(0:NS-1, 2*NF+2:3*NF+2) = FouCoil(icoil)%cmt(0:NS-1, 0:NF)  !x/yc 
           !DoF(icoil)%xof(0:NS-1, 3*NF+3:4*NF+2) = FouCoil(icoil)%smt(0:NS-1, 1:NF)  !x/ys 
           !DoF(icoil)%yof(0:NS-1,      1:  NF+1) = FouCoil(icoil)%cmt(0:NS-1, 0:NF)  !y/xc 
           !DoF(icoil)%yof(0:NS-1,   NF+2:2*NF+1) = FouCoil(icoil)%smt(0:NS-1, 1:NF)  !y/xs 
           DoF(icoil)%yof(0:NS-1, 2*NF+2:3*NF+2) = FouCoil(icoil)%cmt(0:NS-1, 0:NF)  !y/yc
           DoF(icoil)%yof(0:NS-1, 3*NF+3:4*NF+2) = FouCoil(icoil)%smt(0:NS-1, 1:NF)  !y/ys
           DoF(icoil)%zof(0:NS-1, 4*NF+3:5*NF+3) = FouCoil(icoil)%cmt(0:NS-1, 0:NF)  !z/zc
           DoF(icoil)%zof(0:NS-1, 5*NF+4:6*NF+3) = FouCoil(icoil)%smt(0:NS-1, 1:NF)  !z/zs

           do n = 1, NF
              DoF(icoil)%xtof(0:NS,n+1)      = -1.0*FouCoil(icoil)%smt(0:NS,n) * real(n)  !xt/xc
              DoF(icoil)%xtof(0:NS,n+NF+1)   =      FouCoil(icoil)%cmt(0:NS,n) * real(n)  !xt/xs
              DoF(icoil)%ytof(0:NS,n+2*NF+2) = -1.0*FouCoil(icoil)%smt(0:NS,n) * real(n)  !yt/yc
              DoF(icoil)%ytof(0:NS,n+3*NF+2) =      FouCoil(icoil)%cmt(0:NS,n) * real(n)  !yt/ys
              DoF(icoil)%ztof(0:NS,n+4*NF+3) = -1.0*FouCoil(icoil)%smt(0:NS,n) * real(n)  !zt/zc
              DoF(icoil)%ztof(0:NS,n+5*NF+3) =      FouCoil(icoil)%cmt(0:NS,n) * real(n)  !zt/zs
           enddo
           
           ! allocate xyz data
           SALLOCATE( coil(icoil)%xx, (0:coil(icoil)%NS), zero )
           SALLOCATE( coil(icoil)%yy, (0:coil(icoil)%NS), zero )
           SALLOCATE( coil(icoil)%zz, (0:coil(icoil)%NS), zero )
           SALLOCATE( coil(icoil)%xt, (0:coil(icoil)%NS), zero )
           SALLOCATE( coil(icoil)%yt, (0:coil(icoil)%NS), zero )
           SALLOCATE( coil(icoil)%zt, (0:coil(icoil)%NS), zero )
           SALLOCATE( coil(icoil)%xa, (0:coil(icoil)%NS), zero )
           SALLOCATE( coil(icoil)%ya, (0:coil(icoil)%NS), zero )
           SALLOCATE( coil(icoil)%za, (0:coil(icoil)%NS), zero )
           SALLOCATE( coil(icoil)%xb, (0:coil(icoil)%NS), zero )
           SALLOCATE( coil(icoil)%yb, (0:coil(icoil)%NS), zero )
           SALLOCATE( coil(icoil)%zb, (0:coil(icoil)%NS), zero )
           SALLOCATE( coil(icoil)%dl, (0:coil(icoil)%NS), zero )
           SALLOCATE( coil(icoil)%dd, (0:coil(icoil)%NS), zero )

           coil(icoil)%dd = pi2 / NS  ! discretizing factor also set in rdcoils?
           
           SALLOCATE( coil(icoil)%curvature, (0:coil(icoil)%NS), zero )
           SALLOCATE( coil(icoil)%straight, (0:coil(icoil)%NS), zero )
           
           if ( filforce .eq. 1 ) then
              SALLOCATE( coil(icoil)%Fx, (0:coil(icoil)%NS), 0.0 )
              SALLOCATE( coil(icoil)%Fy, (0:coil(icoil)%NS), 0.0 )
              SALLOCATE( coil(icoil)%Fz, (0:coil(icoil)%NS), 0.0 )
              SALLOCATE( coil(icoil)%Bxx, (0:coil(icoil)%NS), 0.0 )
              SALLOCATE( coil(icoil)%Byy, (0:coil(icoil)%NS), 0.0 )
              SALLOCATE( coil(icoil)%Bzz, (0:coil(icoil)%NS), 0.0 )
           endif
           
           if ( calcfb .eq. 1 ) then
              SALLOCATE( coil(icoil)%nfbx, (0:coil(icoil)%NS), 0.0 )
              SALLOCATE( coil(icoil)%nfby, (0:coil(icoil)%NS), 0.0 )
              SALLOCATE( coil(icoil)%nfbz, (0:coil(icoil)%NS), 0.0 )
              SALLOCATE( coil(icoil)%bfbx, (0:coil(icoil)%NS), 0.0 )
              SALLOCATE( coil(icoil)%bfby, (0:coil(icoil)%NS), 0.0 )
              SALLOCATE( coil(icoil)%bfbz, (0:coil(icoil)%NS), 0.0 )
           endif
           
        case(2)
#ifdef dposition
           DoF(icoil)%ND = coil(icoil)%Lc * 5 ! number of DoF for permanent magnet
#else
           DoF(icoil)%ND = coil(icoil)%Lc * 2 ! number of DoF for permanent magnet
#endif
           SALLOCATE(DoF(icoil)%xdof, (1:DoF(icoil)%ND), zero)
        case(3) 
           DoF(icoil)%ND = coil(icoil)%Lc * 1 ! number of DoF for background Bt, Bz
           SALLOCATE(DoF(icoil)%xdof, (1:DoF(icoil)%ND), zero)
        case(coil_type_spline)
           ! get number of DoF for each coil and allocate arrays;
           NS = coil(icoil)%NS
           NCP = Splines(icoil)%NCP
           ND = (3*NCP) ! total variables for geometry
           DoF(icoil)%ND = coil(icoil)%Lc * ND !# of DoF for icoil;
           SALLOCATE(DoF(icoil)%xdof, (1:DoF(icoil)%ND), zero)
           SALLOCATE(DoF(icoil)%xof , (0:coil(icoil)%NS-1, 1:ND), zero)
           SALLOCATE(DoF(icoil)%yof , (0:coil(icoil)%NS-1, 1:ND), zero)
           SALLOCATE(DoF(icoil)%zof , (0:coil(icoil)%NS-1, 1:ND), zero)
           ! allocate and calculate trignometric functions for re-use           
           SALLOCATE( Splines(icoil)%basis_0, (0:NS-1, 0:NCP+2), zero )
           SALLOCATE( Splines(icoil)%basis_1, (0:NS-1, 0:NCP+1), zero )
           SALLOCATE( Splines(icoil)%basis_2, (0:NS-1, 0:NCP), zero )
           SALLOCATE( Splines(icoil)%basis_3, (0:NS-1, 0:NCP-1),   zero )
           SALLOCATE( Splines(icoil)%db_dt  , (0:NS-1, 0:NCP-1),   zero )
           SALLOCATE( Splines(icoil)%db_dt_2, (0:NS-1, 0:NCP-1),   zero )

	   do i =0, coil(icoil)%NS-1
                    Splines(icoil)%eval_points(i) = 1.0*(i)/(coil(icoil)%NS)
           enddo

           ! the derivatives of dx/dv 

           call eval_spline_basis(icoil)           !Compute spline basis and derivatives using analytic derivation
           call eval_spline_basis1(icoil)
           call eval_spline_basis2(icoil)
	   !call check_eval_basis(icoil)    !Check derivatives of basis functions numerically	
	   call enforce_spline_periodicity(icoil)  !Forces the derivates in respect to the first three control points to be equal to the last three   

           DoF(icoil)%xof(0:coil(icoil)%NS-1,      1: NCP) = Splines(icoil)%basis_3(0:coil(icoil)%NS-1, 0:  NCP-1)  !x/xc
           DoF(icoil)%yof(0:coil(icoil)%NS-1, NCP+1:2*NCP) = Splines(icoil)%basis_3(0:coil(icoil)%NS-1, 0:  NCP-1)  !y/yc
           DoF(icoil)%zof(0:coil(icoil)%NS-1, 2*NCP+1:3*NCP) = Splines(icoil)%basis_3(0:coil(icoil)%NS-1, 0:  NCP-1)  !z/zc

           ! allocate xyz data
           SALLOCATE( coil(icoil)%xx, (0:coil(icoil)%NS), zero )
           SALLOCATE( coil(icoil)%yy, (0:coil(icoil)%NS), zero )
           SALLOCATE( coil(icoil)%zz, (0:coil(icoil)%NS), zero )
           SALLOCATE( coil(icoil)%xt, (0:coil(icoil)%NS), zero )
           SALLOCATE( coil(icoil)%yt, (0:coil(icoil)%NS), zero )
           SALLOCATE( coil(icoil)%zt, (0:coil(icoil)%NS), zero )
           SALLOCATE( coil(icoil)%xa, (0:coil(icoil)%NS), zero )
           SALLOCATE( coil(icoil)%ya, (0:coil(icoil)%NS), zero )
           SALLOCATE( coil(icoil)%za, (0:coil(icoil)%NS), zero )
           SALLOCATE( coil(icoil)%dl, (0:coil(icoil)%NS), zero )
           SALLOCATE( coil(icoil)%dd, (0:coil(icoil)%NS), zero )
           SALLOCATE( coil(icoil)%curvature, (0:coil(icoil)%NS), zero )
           SALLOCATE( coil(icoil)%straight, (0:coil(icoil)%NS), zero )

	   coil(icoil)%dd = 1.0/(coil(icoil)%NS)


        case default
           FATAL(AllocData, .true., not supported coil types)
        end select
        
     enddo

     do icoil = 1, Ncoils

        Ndof = Ndof + coil(icoil)%Ic + DoF(icoil)%ND
        if (coil(icoil)%type==1) then
           Tdof = Tdof + 1              + 6*(FouCoil(icoil)%NF)+3
        else if (coil(icoil)%type==coil_type_spline) then
           Tdof = Tdof + 1              + 3*(Splines(icoil)%NCP)
        else 
           Tdof = Tdof + 1 + DoF(icoil)%ND
        end if
        if (DoF(icoil)%ND >= Cdof) Cdof = DoF(icoil)%ND ! find the largest ND for single coil;

      enddo

     if(Ndof == 0) then ! no DOF;
        Nouts = 0
        if(myid==0) write(ounit, *) "AllocData : No free variables; no optimization will be performed."
     endif

     SALLOCATE(    xdof, (1:Ndof), zero ) ! dof vector;
     SALLOCATE( dofnorm, (1:Ndof), one ) ! dof normalized value vector;
     SALLOCATE( evolution, (1:Nouts+1, 0:14), zero ) !evolution array;
     if (Tdof >= 1) then
        SALLOCATE( coilspace, (1:Nouts+1, 1:Tdof), zero ) ! all the coil parameters;
     endif 
     
     ! determine dofnorm
     if ( IsNormalize > 0 ) then 
        ! calculate Inorm and Gnorm
        Inorm = zero ; Mnorm = zero
        icur = 0 ; imag = 0 ! icur for coil current count, imag for dipole count
        do icoil = 1, Ncoils
           if(coil(icoil)%type == 1 .or. coil(icoil)%type == 3 .or. coil(icoil)%type == coil_type_spline ) then  
              ! Fourier representation or central currents
              Inorm = Inorm + coil(icoil)%I**2
              icur = icur + 1
           else if (coil(icoil)%type == 2) then
              ! permanent dipole
              Mnorm = Mnorm + coil(icoil)%I**2
              imag = imag + 1
           endif
        enddo
        Gnorm = (surf(plasma)%vol/(pi*pi2))**(one/three)  ! Gnorm is a hybrid of major and minor radius
        Gnorm = Gnorm * weight_gnorm 
        
        icur = max(1, icur) ; imag = max(1, imag)    ! avoid dividing zero
        Inorm = sqrt(Inorm/icur) * weight_inorm      ! quadratic mean
        Mnorm = sqrt(Mnorm/imag) * weight_mnorm      ! quadratic mean

        if (abs(Gnorm) < machprec) Gnorm = one
        if (abs(Inorm) < machprec) Inorm = one
        if (abs(Mnorm) < machprec) Mnorm = one

        if (IsQuiet<1) then
           if (myid==0) then
              write(ounit, '(8X": Parameter normalizations : "3(A6, ES12.5, 2X))') &
                   'Inorm=', Inorm, 'Gnorm=', Gnorm, 'Mnorm=', Mnorm
           endif
        endif

        ! construct dofnorm
        idof = 0
        do icoil = 1, Ncoils

           if(coil(icoil)%type == 1 .or. coil(icoil)%type == coil_type_spline ) then  ! Fourier representation
              if(coil(icoil)%Ic /= 0) then
                 dofnorm(idof+1) = Inorm
                 idof = idof + 1
              endif

              ND = DoF(icoil)%ND
              if(coil(icoil)%Lc /= 0) then
                 dofnorm(idof+1:idof+ND) = Gnorm
                 idof = idof + ND
              endif
           else if (coil(icoil)%type == 2) then  ! permanent magnets
              if(coil(icoil)%Ic /= 0) then
                 dofnorm(idof+1) = Mnorm
                 idof = idof + 1
              endif
              if(coil(icoil)%Lc /= 0) then
                 !xtmp = max(one, sqrt( coil(icoil)%ox**2 + coil(icoil)%oy**2 + coil(icoil)%oz**2 ) ) ! origin position
                 !mtmp = max(one, sqrt( coil(icoil)%mp**2 + coil(icoil)%mt**2 ) ) ! moment orentation
                 xtmp = Gnorm ! position normalized to Gnorm
                 mtmp = pi    ! orentation normalized to pi
#ifdef dposition
                 dofnorm(idof+1:idof+3) = xtmp
                 dofnorm(idof+4:idof+5) = mtmp
                 idof = idof + 5
#else
                 dofnorm(idof+1:idof+2) = mtmp
                 idof = idof + 2
#endif
              endif
           else if (coil(icoil)%type == 3) then  ! backgroud toroidal/vertical field
              if(coil(icoil)%Ic /= 0) then
                 dofnorm(idof+1) = Inorm
                 idof = idof + 1
              endif

              if(coil(icoil)%Lc /= 0) then
                 if(abs(coil(icoil)%Bz) > sqrtmachprec) then
                    dofnorm(idof+1) = coil(icoil)%Bz
                 else
                    dofnorm(idof+1) = one
                 endif
                 idof = idof + 1
              endif
           else
              STOP " wrong coil type in rdcoils"
              call MPI_ABORT(MPI_COMM_FOCUS, 1, ierr)
           endif

        enddo !end do icoil;
        FATAL( AllocData , idof .ne. Ndof, counting error in unpacking )

     endif

  endif

  !--------------------------------------------------------------------------------------------- 
  if (type == 0 .or. type == 1) then  ! 0-order cost functions related arrays;

     ! Bnorm and Bharm needed;
     if (weight_bnorm > sqrtmachprec .or. weight_bharm > sqrtmachprec .or. weight_sbnorm > sqrtmachprec .or. IsQuiet <= -2) then
        SALLOCATE(         bn, (0:Nteta-1,0:Nzeta-1), zero ) ! Bn from coils;        
        SALLOCATE( surf(isurf)%bn, (0:Nteta-1,0:Nzeta-1), zero ) ! total Bn;
        SALLOCATE( surf(isurf)%Bx, (0:Nteta-1,0:Nzeta-1), zero ) ! Bx on the surface;
        SALLOCATE( surf(isurf)%By, (0:Nteta-1,0:Nzeta-1), zero ) ! By on the surface;
        SALLOCATE( surf(isurf)%Bz, (0:Nteta-1,0:Nzeta-1), zero ) ! Bz on the surface;
        SALLOCATE(         Bm, (0:Nteta-1,0:Nzeta-1), zero ) ! |B| on the surface;       
        SALLOCATE( dBx, (0:Cdof,0:Cdof), zero ) ! d^2Bx/(dx1,dx2) on each coil; Cdof is the max coil dof
        SALLOCATE( dBy, (0:Cdof,0:Cdof), zero ) ! d^2By/(dx1,dx2) on each coil;
        SALLOCATE( dBz, (0:Cdof,0:Cdof), zero ) ! d^2Bz/(dx1,dx2) on each coil;        
     endif

     ! Bharm needed;
     if (weight_bharm > sqrtmachprec) then
        call readbmn
        SALLOCATE(  Bmnc , (1:NBmn), zero )  ! current Bmn cos values;
        SALLOCATE(  Bmns , (1:NBmn), zero )  ! current Bmn sin values;
        SALLOCATE( iBmnc , (1:NBmn), zero )
        SALLOCATE( iBmns , (1:NBmn), zero )
     endif
            
  endif

  !---------------------------------------------------------------------------------------------
  if (type == 1) then ! 1st-order cost functions related arrays;
     
     FATAL( AllocData, Ndof < 1, INVALID Ndof value )
     SALLOCATE( t1E, (1:Ndof), zero )
     SALLOCATE( deriv, (1:Ndof, 0:11), zero )

     ! Bnorm related;
     if (weight_bnorm > sqrtmachprec .or. weight_bharm > sqrtmachprec .or. weight_sbnorm > sqrtmachprec ) then
        SALLOCATE( t1B, (1:Ndof), zero )  !total d bnorm / d x;
        SALLOCATE( dBn, (1:Ndof), zero )  !total d Bn / d x;
        SALLOCATE( dBm, (1:Ndof), zero )  !total d Bm / d x;
        SALLOCATE( d1B, (1:Ndof,0:Nteta-1,0:Nzeta-1), zero ) ! discretized dBn
     endif

     ! Bharm related;
     if (weight_bharm > sqrtmachprec) then
        SALLOCATE( t1H, (1:Ndof), zero )
!       SALLOCATE( dB , (1:Ndof, 0:Nteta-1, 0:Nzeta-1), zero ) !distribution of dB/dx;
     endif

     ! tflux needed;
     if (weight_tflux > sqrtmachprec) then
        SALLOCATE( t1F,  (1:Ndof), zero )
     endif

     ! isum needed;
     if (weight_isum > sqrtmachprec) then
        SALLOCATE( t1I,  (1:Ndof), zero )
     endif

     ! ttlen needed;
     if (weight_ttlen > sqrtmachprec) then
        SALLOCATE( t1L,  (1:Ndof), zero )
     endif

     ! curv needed;
     if (weight_curv > sqrtmachprec) then
        SALLOCATE( t1K,  (1:Ndof), zero )
     endif

     ! ccsep needed;
     if (weight_ccsep > sqrtmachprec) then
        SALLOCATE( t1C,  (1:Ndof), zero )
     endif

     if (weight_straight > sqrtmachprec) then
        SALLOCATE( t1Str,  (1:Ndof), zero )
     endif

     ! cssep needed;
     if (weight_cssep > sqrtmachprec) then
        SALLOCATE( t1S,  (1:Ndof), zero )
     endif

     ! tors needed;
     if (weight_tors > sqrtmachprec) then
        SALLOCATE( t1T,  (1:Ndof), zero )
     endif 

     ! nissin needed;
     if (weight_nissin > 0.0_dp) then
        SALLOCATE( t1N,  (1:Ndof), zero )
     endif

     ! stochastic bnorm needed;
     if (weight_sbnorm > 0.0_dp) then
        SALLOCATE( t1Bavg,  (1:Ndof), zero )
     endif

     ! L-M algorithn enabled
     if (LM_maxiter > 0) then
        LM_mfvec = 0 ! number of total cost functions

        if (weight_bnorm > sqrtmachprec) then
           ibnorm = LM_mfvec
           mbnorm = Nteta*Nzeta
           LM_mfvec = LM_mfvec + mbnorm
        endif

        if (weight_bharm > sqrtmachprec) then 
           ibharm = LM_mfvec
           mbharm = 2*NBmn
           LM_mfvec = LM_mfvec + mbharm
        endif
        
        if (weight_tflux > sqrtmachprec) then
           itflux = LM_mfvec
           mtflux = Nzeta
           LM_mfvec = LM_mfvec + mtflux
        endif

        if (weight_isum > sqrtmachprec) then
           iisum = LM_mfvec
           misum = 1 ! put together
           LM_mfvec = LM_mfvec + misum
        endif        
        
        if (weight_ttlen > sqrtmachprec) then
           ittlen = LM_mfvec
           mttlen = Ncoils - Nfixgeo
           LM_mfvec = LM_mfvec + mttlen
        endif
       
        if (weight_curv > sqrtmachprec) then
           icurv = LM_mfvec
           mcurv = Ncoils - Nfixgeo
           LM_mfvec = LM_mfvec + mcurv
        endif
 
        if (weight_cssep > sqrtmachprec) then
           icssep = LM_mfvec
           mcssep = Ncoils - Nfixgeo
           LM_mfvec = LM_mfvec + mcssep
        endif

        if (weight_tors > sqrtmachprec) then
           itors = LM_mfvec
           mtors = Ncoils - Nfixgeo
           LM_mfvec = LM_mfvec + mtors
        endif

        if (weight_nissin > 0.0_dp) then
           inissin = LM_mfvec
           mnissin = Ncoils - Nfixgeo
           LM_mfvec = LM_mfvec + mnissin
        endif
        
        FATAL( AllocData, LM_mfvec <= 0, INVALID number of cost functions )
        SALLOCATE( LM_fvec, (1:LM_mfvec), zero )
        SALLOCATE( LM_fjac , (1:LM_mfvec, 1:Ndof), zero )

        if (myid == 0) write(ounit, '("datalloc: total number of cost functions for L-M is "I0)') LM_mfvec
        
     endif

  endif
  !--------------------------------------------------------------------------------------------- 

  return
end subroutine AllocData
