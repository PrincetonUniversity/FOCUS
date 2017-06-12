!title (packdof) ! paking degree of freedom (dof) into one vector

!latex  \briefly{Packing all the free coil parameters into a one rank vector.}
!latex 
!latex  \subsection{Overview}
!latex  \bi
!latex  \item The \inputvar{case\_coils} determines the packing and unpacking patern. 
!latex  \item \inputvar{case\_coils} = 1: Coils are represented with Fourier series.
!latex  \item For each coil, the number of DOF is $6N_F+3$ ($\sin 0$ terms are omitted.)
!latex  \be
!latex  \vect{X_i} = \left[ \overbrace{I, \underbrace{X_{c,0}, \cdots, X_{c,N}}_\text{N+1}, 
!latex  \underbrace{X_{s,1}, \cdots, X_{s,N}}_\text{N}, Y_{c,0}, \cdots, Z_{s,N}}^\text{6N+4} \right ]
!latex  \ee
!latex  \item Coil currents/geometry can also be fixed, and they are determined by coil\%Ic and coil\%Lc.
!latex  \item The total number of DOF $Ndof$ should be
!latex  \be
!latex  Ndof = Ncoils \ * \ (6 * NFcoil + 4) \, - \, Nfixcur \, - \, Nfixgeo \ * \ (6 * NFcoil + 3)
!latex  \ee
!latex  \ei
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

SUBROUTINE packdof(lxdof)
  !---------------------------------------------------------------------------------------------
  ! Pack all DOF into one vector;
  ! DATE: 2017/03/19
  !--------------------------------------------------------------------------------------------- 
  use globals, only : zero, myid, ounit, &
                    & case_coils, Ncoils, coil, DoF, Ndof, Inorm, Gnorm
  implicit none
  include "mpif.h"

  REAL    :: lxdof(1:Ndof)
  INTEGER :: idof, icoil, ND, astat, ierr
  !--------------------------------------------------------------------------------------------- 
     ! reset xdof;
     lxdof = zero

     call packcoil !pack coil parameters into DoF;
     ! packing;
     idof = 0
     do icoil = 1, Ncoils

        if(coil(icoil)%Ic /= 0) then 
           lxdof(idof+1) = coil(icoil)%I / Inorm
           idof = idof + 1
        endif

        ND = DoF(icoil)%ND
        if(coil(icoil)%Lc /= 0) then
           lxdof(idof+1:idof+ND) = DoF(icoil)%xdof(1:ND) / Gnorm
           idof = idof + ND
        endif

     enddo !end do icoil;

  !--------------------------------------------------------------------------------------------- 
  FATAL( packdof , idof .ne. Ndof, counting error in packing )

  !write(ounit, *) "pack ", lxdof(1)

  return
END SUBROUTINE packdof

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

SUBROUTINE unpacking(lxdof)
  !---------------------------------------------------------------------------------------------
  ! UnPack all DOF from one vector;
  ! DATE: 2017/04/03
  !--------------------------------------------------------------------------------------------- 
  use globals, only : zero, myid, ounit, &
                    & case_coils, Ncoils, coil, DoF, Ndof, Inorm, Gnorm
  implicit none
  include "mpif.h"

  REAL    :: lxdof(1:Ndof)
  INTEGER :: idof, icoil, ND, astat, ierr, ifirst
  !---------------------------------------------------------------------------------------------
  !FATAL( packdof, .not. allocated(xdof), Please allocate xdof first. )

  idof = 0 ; ifirst = 0
  do icoil = 1, Ncoils

     if(coil(icoil)%Ic /= 0) then
        coil(icoil)%I = lxdof(idof+1) * Inorm
        idof = idof + 1
     endif

     ND = DoF(icoil)%ND
     if(coil(icoil)%Lc /= 0) then
        DoF(icoil)%xdof(1:ND) = lxdof(idof+1:idof+ND) * Gnorm
        idof = idof + ND
     endif

  enddo !end do icoil;

  !--------------------------------------------------------------------------------------------- 
  FATAL( unpacking , idof .ne. Ndof, counting error in unpacking )

  call unpackcoil !unpack DoF to coil parameters;
  call discoil(ifirst)

  return
END SUBROUTINE unpacking

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

SUBROUTINE packcoil
  !---------------------------------------------------------------------------------------------
  ! pack coil representation variables into DoF (only geometries without currents);
  ! DATE: 2017/03/25
  !--------------------------------------------------------------------------------------------- 
  use globals, only : zero, myid, ounit, case_coils, Ncoils, coil, FouCoil, DoF
  implicit none
  include "mpif.h"

  INTEGER  :: icoil, idof, NF, ierr, astat

  FATAL( packcoil, .not. allocated(coil)   , illegal )
  FATAL( packcoil, .not. allocated(FouCoil), illegal )
  FATAL( packcoil, .not. allocated(DoF)    , illegal )

  do icoil = 1, Ncoils

     select case (coil(icoil)%itype)
     !--------------------------------------------------------------------------------------------- 
     case(1)
        ! get number of DoF for each coil and allocate arrays;
        NF = FouCoil(icoil)%NF
        DoF(icoil)%xdof = zero

        !pack Fourier series;
        idof = 0
        if(coil(icoil)%Lc /= 0) then
           DoF(icoil)%xdof(idof+1 : idof+NF+1) = FouCoil(icoil)%xc(0:NF); idof = idof + NF + 1
           DoF(icoil)%xdof(idof+1 : idof+NF  ) = FouCoil(icoil)%xs(1:NF); idof = idof + NF    
           DoF(icoil)%xdof(idof+1 : idof+NF+1) = FouCoil(icoil)%yc(0:NF); idof = idof + NF + 1
           DoF(icoil)%xdof(idof+1 : idof+NF  ) = FouCoil(icoil)%ys(1:NF); idof = idof + NF    
           DoF(icoil)%xdof(idof+1 : idof+NF+1) = FouCoil(icoil)%zc(0:NF); idof = idof + NF + 1
           DoF(icoil)%xdof(idof+1 : idof+NF  ) = FouCoil(icoil)%zs(1:NF); idof = idof + NF    
        endif
        FATAL( packcoil , idof .ne. DoF(icoil)%ND, counting error in packing )
        
     !---------------------------------------------------------------------------------------------
     case default
        FATAL(packcoil, .true., not supported coil types)
     end select     

  enddo ! end do icoil;

END SUBROUTINE packcoil

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

SUBROUTINE unpackcoil
  !---------------------------------------------------------------------------------------------
  ! pack coil representation variables into DoF (only geometries without currents);
  ! DATE: 2017/03/25
  !--------------------------------------------------------------------------------------------- 
  use globals, only : zero, myid, ounit, case_coils, Ncoils, coil, FouCoil, DoF
  implicit none
  include "mpif.h"

  INTEGER  :: icoil, idof, NF, ierr, astat

  FATAL( unpackcoil, .not. allocated(coil)   , illegal )
  FATAL( unpackcoil, .not. allocated(FouCoil), illegal )
  FATAL( unpackcoil, .not. allocated(DoF)    , illegal )

  do icoil = 1, Ncoils

     select case (coil(icoil)%itype)
        !--------------------------------------------------------------------------------------------- 
     case(1)
        ! get number of DoF for each coil and allocate arrays;
        NF = FouCoil(icoil)%NF
        idof = 0
        if (coil(icoil)%Lc /= 0) then
           !unpack Fourier series;
           FouCoil(icoil)%xc(0:NF) = DoF(icoil)%xdof(idof+1 : idof+NF+1) ; idof = idof + NF + 1
           FouCoil(icoil)%xs(1:NF) = DoF(icoil)%xdof(idof+1 : idof+NF  ) ; idof = idof + NF    
           FouCoil(icoil)%yc(0:NF) = DoF(icoil)%xdof(idof+1 : idof+NF+1) ; idof = idof + NF + 1
           FouCoil(icoil)%ys(1:NF) = DoF(icoil)%xdof(idof+1 : idof+NF  ) ; idof = idof + NF    
           FouCoil(icoil)%zc(0:NF) = DoF(icoil)%xdof(idof+1 : idof+NF+1) ; idof = idof + NF + 1
           FouCoil(icoil)%zs(1:NF) = DoF(icoil)%xdof(idof+1 : idof+NF  ) ; idof = idof + NF   
        endif
        FATAL( packcoil , idof .ne. DoF(icoil)%ND, counting error in packing )

        !---------------------------------------------------------------------------------------------
     case default
        FATAL(packcoil, .true., not supported coil types)
     end select

  enddo ! end do icoil;

END SUBROUTINE unpackcoil
