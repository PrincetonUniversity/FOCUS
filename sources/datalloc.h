!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine AllocData(itype)
!------------------------------------------------------------------------------------------------------ 
! DATE:  04/05/2017
! Allocate data before using them, especially for those used several times;
! part can be : -1('dof'), 0('costfun0'), 1('costfun1')
!------------------------------------------------------------------------------------------------------
  use globals
  implicit none
  include "mpif.h"

  INTEGER, intent(in) :: itype

  INTEGER             :: icoil, idof, ND, NF

  !-------------------------------------------------------------------------------------------
  if (itype == -1) then ! dof related data;

     Cdof = 0; Ndof = 0; Tdof = 0

     do icoil = 1, Ncoils*Npc
        
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
        case default
           FATAL(AllocData, .true., not supported coil types)
        end select
        
     enddo

     do icoil = 1, Ncoils

        Ndof = Ndof + coil(icoil)%Ic + DoF(icoil)%ND
        Tdof = Tdof + 1              + 6*(FouCoil(icoil)%NF)+3
        if (DoF(icoil)%ND >= Cdof) Cdof = DoF(icoil)%ND ! find the largest ND for single coil;

      enddo

     if(Ndof == 0) then ! no DOF;
        Nouts = 0
        if(myid==0) write(ounit, *) "AllocData : No free variables; no optimization will be performed."
     endif

     SALLOCATE(    xdof, (1:Ndof), zero ) ! dof vector;
     SALLOCATE( dofnorm, (1:Ndof), zero ) ! dof normalized value vector;
     SALLOCATE( evolution, (1:Nouts+1, 0:7), zero ) !evolution array;
     SALLOCATE( coilspace, (1:Nouts+1, 1:Tdof), zero ) ! all the coil parameters;
     
     idof = 0
     do icoil = 1, Ncoils

        if(coil(icoil)%Ic /= 0) then
           dofnorm(idof+1) = Inorm
           idof = idof + 1
        endif

        ND = DoF(icoil)%ND
        if(coil(icoil)%Lc /= 0) then
           dofnorm(idof+1:idof+ND) = Gnorm
           idof = idof + ND
        endif

     enddo !end do icoil;
     FATAL( AllocData , idof .ne. Ndof, counting error in unpacking )

  endif

  !--------------------------------------------------------------------------------------------- 
  if (itype == 0 .or. itype == 1) then  ! 0-order cost functions related arrays;

     ! Bnorm and Bharm needed;
     if (weight_bnorm > sqrtmachprec .or. weight_bharm > sqrtmachprec .or. IsQuiet <= -2) then
        SALLOCATE(         bn, (0:Nteta-1,0:Nzeta-1), zero ) !Bn from coils;        
        SALLOCATE( surf(1)%bn, (0:Nteta-1,0:Nzeta-1), zero ) !total Bn;
        SALLOCATE( surf(1)%Bx, (0:Nteta-1,0:Nzeta-1), zero ) !Bx on the surface;
        SALLOCATE( surf(1)%By, (0:Nteta-1,0:Nzeta-1), zero ) !By on the surface;
        SALLOCATE( surf(1)%Bz, (0:Nteta-1,0:Nzeta-1), zero ) !Bz on the surface;
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
        SALLOCATE( t1B, (1:Ndof), zero )                       !total dB/dx;
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
