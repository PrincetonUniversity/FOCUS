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

     Ndof = 0; Tdof = 0

     do icoil = 1, Ncoils
        select case (coil(icoil)%itype)
        
        case(1)
           ! get number of DoF for each coil and allocate arrays;
           NF = FouCoil(icoil)%NF
           ND = (6*NF + 3) ! total variables for geometry
           DoF(icoil)%ND = coil(icoil)%Lc * ND !# of DoF for icoil;
           SALLOCATE(DoF(icoil)%xdof, (1:DoF(icoil)%ND), zero)
           SALLOCATE(DoF(icoil)%xof , (1:coil(icoil)%NS, 1:ND), zero)
           SALLOCATE(DoF(icoil)%yof , (1:coil(icoil)%NS, 1:ND), zero)
           SALLOCATE(DoF(icoil)%zof , (1:coil(icoil)%NS, 1:ND), zero)
        case default
           FATAL(AllocData, .true., not supported coil types)
        end select

        Ndof = Ndof + coil(icoil)%Ic + DoF(icoil)%ND
        Tdof = Tdof + 1              + ND

     enddo

     if(Ndof == 0) then ! no DOF;
        Nouts = 0
        if(myid==0) write(ounit, *) "AllocData : No free variables; no optimization will be performed."
     endif

     SALLOCATE(    xdof, (1:Ndof), zero ) ! dof vector;
     SALLOCATE( dofnorm, (1:Ndof), zero ) ! dof normalized value vector;
     SALLOCATE( evolution, (1:Nouts+1, 0:8), zero ) !evolution array;
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

        do icoil = 1, Ncoils
           ND = DoF(icoil)%ND
           SALLOCATE( coil(icoil)%Bx, (0:ND, 0:ND), zero )   ! Bx;
           SALLOCATE( coil(icoil)%By, (0:ND, 0:ND), zero )   ! By;  
           SALLOCATE( coil(icoil)%Bz, (0:ND, 0:ND), zero )   ! Bz;
        enddo
     endif

     ! Bharm needed;
     if (weight_bharm > sqrtmachprec) then
        call readbmn
        SALLOCATE(  Bmnc , (1:NBmn), zero )  ! current Bmn cos values;
        SALLOCATE(  Bmns , (1:NBmn), zero )  ! current Bmn sin values;
     endif

     ! tflux needed;
     if (weight_tflux > sqrtmachprec .or. IsQuiet <= -2) then
        do icoil = 1, Ncoils
           ND = DoF(icoil)%ND
           SALLOCATE( coil(icoil)%Ax, (0:ND, 0:ND), zero )   ! Ax;
           SALLOCATE( coil(icoil)%Ay, (0:ND, 0:ND), zero )   ! Ay;  
           SALLOCATE( coil(icoil)%Az, (0:ND, 0:ND), zero )   ! Az;
        enddo
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
        SALLOCATE( dB , (1:Ndof, 0:Nteta-1, 0:Nzeta-1), zero ) !distribution of dB/dx;
     endif

     ! Bharm related;
     if (weight_bharm > sqrtmachprec) then
        SALLOCATE( t1H,  (1:Ndof), zero )
     endif

     ! tflux needed;
     if (weight_tflux > sqrtmachprec) then
        SALLOCATE( t1F,  (1:Ndof), zero )
     endif

     ! ttlen needed;
     if (weight_ttlen > sqrtmachprec) then
        SALLOCATE( t1L,  (1:Ndof), zero )
     endif

  endif
  !--------------------------------------------------------------------------------------------- 

  return
end subroutine AllocData
