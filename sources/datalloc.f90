!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine AllocData(itype)
  !------------------------------------------------------------------------------------------------------ 
  ! DATE:  04/05/2017
  ! Allocate data before using them, especially for those used several times;
  ! part can be : -1('dof'), 0('costfun0'), 1('costfun1')
  !------------------------------------------------------------------------------------------------------
  use globals
  use bnorm_mod
  implicit none
  include "mpif.h"

  INTEGER, intent(in) :: itype

  INTEGER             :: icoil, idof, ND, NF, icur, imag, icpu, dof_array(0:ncpu-1)
  REAL                :: xtmp, mtmp

  !-------------------------------------------------------------------------------------------
  if (itype == -1) then ! dof related data;

     Cdof = 0; Ndof = 0; Tdof = 0; ldof=0; dof_offset=0
     dof_array = 0 

     do icoil = 1, Ncoils

        select case (coil(icoil)%itype)       
        case(1)
           ! get number of DoF for each coil and allocate arrays;
           NF = FouCoil(icoil)%NF
           ND = (6*NF + 3) ! total variables for geometry
           DoF(icoil)%ND = coil(icoil)%Lc * ND !# of DoF for icoil;
           SALLOCATE(DoF(icoil)%xdof, (1:DoF(icoil)%ND), zero)
           SALLOCATE(DoF(icoil)%xof , (0:coil(icoil)%NS-1, 1:ND), zero)
           SALLOCATE(DoF(icoil)%yof , (0:coil(icoil)%NS-1, 1:ND), zero)
           SALLOCATE(DoF(icoil)%zof , (0:coil(icoil)%NS-1, 1:ND), zero)
        case(2)
#ifdef dposition
           ND = 5
#else
           ND = 2           
#endif
           DoF(icoil)%ND = coil(icoil)%Lc * ND ! number of DoF for permanent magnet
           SALLOCATE(DoF(icoil)%xdof, (1:ND), zero)
        case(3) 
           ND = 1
           DoF(icoil)%ND = coil(icoil)%Lc * ND ! number of DoF for background Bt, Bz
           SALLOCATE(DoF(icoil)%xdof, (1:ND), zero)
        case default
           FATAL(AllocData, .true., not supported coil types)
        end select

     enddo

     do icoil = 1, Ncoils

        ldof = ldof + coil(icoil)%Ic + DoF(icoil)%ND
        if (allocated(FouCoil)) then
           Tdof = Tdof + 1              + 6*(FouCoil(icoil)%NF)+3
        else 
           Tdof = Tdof + coil(icoil)%Ic + DoF(icoil)%ND
        end if
        if (DoF(icoil)%ND >= Cdof) Cdof = DoF(icoil)%ND ! find the largest ND for single coil;

     enddo

     do icpu = 0, ncpu-1
        if (myid == icpu) then
           dof_array(icpu) = ldof
        endif
     enddo

     CALL MPI_ALLREDUCE(ldof, Ndof, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr )
     CALL MPI_ALLREDUCE(MPI_IN_PLACE, dof_array, ncpu, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr )

     FATAL( datalloc, sum(dof_array) .ne. Ndof, error in counting dof number)

     do icpu = 1, ncpu-1   ! calculate the offset in dof arrays
        if (myid == icpu) dof_offset = sum(dof_array(0:icpu-1))
     enddo     

     if(Ndof == 0) then ! no DOF;
        Nouts = 0
        if(myid==0) write(ounit, *) "AllocData : No free variables; no optimization will be performed."
     endif

     SALLOCATE(    xdof, (1:Ndof), zero ) ! dof vector;
     SALLOCATE( dofnorm, (1:Ndof), one ) ! dof normalized value vector;
     SALLOCATE( evolution, (1:Nouts+1, 0:7), zero ) !evolution array;
     ! SALLOCATE( coilspace, (1:Nouts+1, 1:Tdof), zero ) ! all the coil parameters;

     ! determine dofnorm
     if ( IsNormalize > 0 ) then 
        ! calculate Inorm and Gnorm
        Inorm = zero ; Mnorm = zero
        icur = 0 ; imag = 0 ! icur for coil current count, imag for dipole count
        do icoil = 1, Ncoils
           if(coil(icoil)%itype == 1 .or. coil(icoil)%itype == 3 ) then  
              ! Fourier representation or central currents
              Inorm = Inorm + coil(icoil)%I**2
              icur = icur + 1
           else if (coil(icoil)%itype == 2) then
              ! permanent dipole
              Mnorm = Mnorm + coil(icoil)%I**2
              imag = imag + 1
           endif
        enddo
        Gnorm = (surf(1)%vol/(pi*pi2))**(one/three)  ! Gnorm is a hybrid of major and minor radius
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
        idof = dof_offset
        do icoil = 1, Ncoils

           if(coil(icoil)%itype == 1) then  ! Fourier representation
              if(coil(icoil)%Ic /= 0) then
                 dofnorm(idof+1) = Inorm
                 idof = idof + 1
              endif

              ND = DoF(icoil)%ND
              if(coil(icoil)%Lc /= 0) then
                 dofnorm(idof+1:idof+ND) = Gnorm
                 idof = idof + ND
              endif
           else if (coil(icoil)%itype == 2) then  ! permanent magnets
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
           else if (coil(icoil)%itype == 3) then  ! backgroud toroidal/vertical field
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
              call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
           endif

        enddo !end do icoil;
        FATAL( AllocData , idof-dof_offset .ne. ldof, counting error in unpacking )

     endif

     call MPI_ALLREDUCE( MPI_IN_PLACE, dofnorm, Ndof, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )

  endif

  !--------------------------------------------------------------------------------------------- 
  if (itype == 0 .or. itype == 1) then  ! 0-order cost functions related arrays;

     ! Bnorm and Bharm needed;
     if (weight_bnorm > sqrtmachprec .or. weight_bharm > sqrtmachprec .or. IsQuiet <= -2) then
        SALLOCATE(         bn, (0:Nteta-1,0:Nzeta-1), zero ) ! Bn from coils;        
        SALLOCATE( surf(1)%bn, (0:Nteta-1,0:Nzeta-1), zero ) ! total Bn;
        SALLOCATE( surf(1)%Bx, (0:Nteta-1,0:Nzeta-1), zero ) ! Bx on the surface;
        SALLOCATE( surf(1)%By, (0:Nteta-1,0:Nzeta-1), zero ) ! By on the surface;
        SALLOCATE( surf(1)%Bz, (0:Nteta-1,0:Nzeta-1), zero ) ! Bz on the surface;
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
  if (itype == 1) then ! 1st-order cost functions related arrays;

     FATAL( AllocData, Ndof < 1, INVALID Ndof value )
     SALLOCATE( t1E, (1:Ndof), zero )
     SALLOCATE( deriv, (1:Ndof, 0:6), zero )

     ! Bnorm related;
     if (weight_bnorm > sqrtmachprec .or. weight_bharm > sqrtmachprec) then
        SALLOCATE( t1B, (1:Ndof), zero )  !total d bnorm / d x;
        SALLOCATE( dBn, (1:Ndof), zero )  !total d Bn / d x;
        SALLOCATE( dBm, (1:Ndof), zero )  !total d Bm / d x;   
     endif

     if ( weight_bharm > sqrtmachprec .or. LM_maxiter > 0 ) then
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

     ! ttlen needed;
     if (weight_ttlen > sqrtmachprec) then
        SALLOCATE( t1L,  (1:Ndof), zero )
     endif

     ! cssep needed;
     if (weight_cssep > sqrtmachprec) then
        SALLOCATE( t1S,  (1:Ndof), zero )
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

        if (weight_ttlen > sqrtmachprec) then
           ittlen = LM_mfvec
           mttlen = Ncoils - Nfixgeo
           LM_mfvec = LM_mfvec + mttlen
        endif

        if (weight_cssep > sqrtmachprec) then
           icssep = LM_mfvec
           mcssep = Ncoils - Nfixgeo
           LM_mfvec = LM_mfvec + mcssep
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
