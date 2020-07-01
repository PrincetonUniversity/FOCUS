!title (packdof) ! paking degree of freedom (dof) into one vector

!latex  \briefly{Packing all the free coil parameters into a one rank vector.}
!latex 
!latex  \subsection{Overview}
!latex  \bi
!latex  \item The \inputvar{case\_coils} determines the packing and unpacking patern. 
!latex  \item \inputvar{case\_coils} = 1: Coils are represented with Fourier series.
!latex  \item For each coil, the number of DOF is $6N_F+4$ ($\sin 0$ terms are omitted.)
!latex  \be
!latex  \vect{X_i} = \left [ \overbrace{I, \underbrace{X_{c,0}, \cdots, X_{c,N}}_{\text{N+1}}, 
!latex  \underbrace{X_{s,1}, \cdots, X_{s,N}}_{\mathrm{N}}, Y_{c,0}, \cdots, Z_{s,N}}^{\mathrm{6N+4}} \right ]
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
  use globals, only : dp, zero, myid, ounit, MPI_COMM_FOCUS, &
                    & case_coils, Ncoils, coil, DoF, Ndof, DoFnorm, coil_type_spline, CPCoil
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

     select case (coil(icoil)%type)
     !--------------------------------------------------------------------------------------------- 
     case(1)

        if(coil(icoil)%Ic /= 0) then 
           lxdof(idof+1) = coil(icoil)%I
           idof = idof + 1
        endif

        ND = DoF(icoil)%ND
        if(coil(icoil)%Lc /= 0) then
           lxdof(idof+1:idof+ND) = DoF(icoil)%xdof(1:ND)
           idof = idof + ND
        endif
     !--------------------------------------------------------------------------------------------- 
     case(2) 
        if(coil(icoil)%Ic /= 0) then 
           lxdof(idof+1) = coil(icoil)%I
           idof = idof + 1
        endif
        ND = DoF(icoil)%ND
        if(coil(icoil)%Lc /= 0) then
           lxdof(idof+1:idof+ND) = DoF(icoil)%xdof(1:ND)
           idof = idof + ND
        endif
     !--------------------------------------------------------------------------------------------- 
     case(3)
        if(coil(icoil)%Ic /= 0) then 
           lxdof(idof+1) = coil(icoil)%I
           idof = idof + 1
        endif

        if(coil(icoil)%Lc /= 0) then
           lxdof(idof+1) = DoF(icoil)%xdof(1)
           idof = idof + 1
        endif
      !--------------------------------------------------------------------------------------------- 
     case(coil_type_spline)

        if(coil(icoil)%Ic /= 0) then 
           lxdof(idof+1) = coil(icoil)%I
           idof = idof + 1
        endif

        ND = DoF(icoil)%ND
        if(coil(icoil)%Lc /= 0) then
           lxdof(idof+1:idof+ND) = DoF(icoil)%xdof(1:ND)
           idof = idof + ND
        endif
       
     !---------------------------------------------------------------------------------------------
     case default
        FATAL(packdof01, .true., not supported coil types)
     end select

  enddo !end do icoil;

  !--------------------------------------------------------------------------------------------- 
  FATAL( packdof02 , idof .ne. Ndof, counting error in packing )

  !write(ounit, *) "pack ", lxdof(1)
  lxdof = lxdof / DoFnorm
  call mpi_barrier(MPI_COMM_FOCUS, ierr)

  return
END SUBROUTINE packdof

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

SUBROUTINE unpacking(lxdof)
  !---------------------------------------------------------------------------------------------
  ! UnPack all DOF from one vector;
  ! DATE: 2017/04/03
  !--------------------------------------------------------------------------------------------- 
  use globals, only: dp, zero, myid, ounit, MPI_COMM_FOCUS, &
       & case_coils, Ncoils, coil, DoF, Ndof, DoFnorm,coil_type_spline, CPCoil
  implicit none
  include "mpif.h"

  REAL    :: lxdof(1:Ndof)
  INTEGER :: idof, icoil, ND, astat, ierr, ifirst
  !---------------------------------------------------------------------------------------------
  !FATAL( packdof, .not. allocated(xdof), Please allocate xdof first. )

  idof = 0 ; ifirst = 0
  do icoil = 1, Ncoils

     select case (coil(icoil)%type)
     !--------------------------------------------------------------------------------------------- 
     case(1)

        if(coil(icoil)%Ic /= 0) then
           coil(icoil)%I = lxdof(idof+1) * dofnorm(idof+1)
           idof = idof + 1
        endif

        ND = DoF(icoil)%ND
        if(coil(icoil)%Lc /= 0) then
           DoF(icoil)%xdof(1:ND) = lxdof(idof+1:idof+ND) * dofnorm(idof+1:idof+ND)
           idof = idof + ND
        endif

     !--------------------------------------------------------------------------------------------- 
     case(2) 
        if(coil(icoil)%Ic /= 0) then 
           coil(icoil)%I = lxdof(idof+1) * dofnorm(idof+1)
           idof = idof + 1
        endif
        ND = DoF(icoil)%ND
        if(coil(icoil)%Lc /= 0) then
           DoF(icoil)%xdof(1:ND) = lxdof(idof+1:idof+ND) * dofnorm(idof+1:idof+ND)
           idof = idof + ND
        endif
     !--------------------------------------------------------------------------------------------- 
     case(3)
        if(coil(icoil)%Ic /= 0) then 
           coil(icoil)%I = lxdof(idof+1) * dofnorm(idof+1)
           idof = idof + 1
        endif

        if(coil(icoil)%Lc /= 0) then
           DoF(icoil)%xdof(1) = lxdof(idof+1) * dofnorm(idof+1)
           idof = idof + 1
        endif
      !--------------------------------------------------------------------------------------------- 
     case(coil_type_spline)

        if(coil(icoil)%Ic /= 0) then
           coil(icoil)%I = lxdof(idof+1) * dofnorm(idof+1)
           idof = idof + 1
        endif

        ND = DoF(icoil)%ND
        if(coil(icoil)%Lc /= 0) then
           DoF(icoil)%xdof(1:ND) = lxdof(idof+1:idof+ND) * dofnorm(idof+1:idof+ND)
           idof = idof + ND
        endif

         
     !---------------------------------------------------------------------------------------------
     case default
        FATAL(unpacking01, .true., not supported coil types)
     end select

  enddo !end do icoil;

  !--------------------------------------------------------------------------------------------- 
  FATAL( unpacking02 , idof .ne. Ndof, counting error in unpacking )

  call unpackcoil !unpack DoF to coil parameters;
  call discoil(ifirst)

  call mpi_barrier(MPI_COMM_FOCUS, ierr)

  return
END SUBROUTINE unpacking

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

SUBROUTINE packcoil
  !---------------------------------------------------------------------------------------------
  ! pack coil representation variables into DoF (only geometries without currents);
  ! DATE: 2017/03/25
  !--------------------------------------------------------------------------------------------- 
  use globals, only: dp, zero, myid, ounit, case_coils, Ncoils, coil, FouCoil, DoF, MPI_COMM_FOCUS,coil_type_spline, CPCoil
  implicit none
  include "mpif.h"

  INTEGER  :: icoil, idof, NF, ierr, astat,NCP

  FATAL( packcoil01, .not. allocated(coil)   , illegal )
  ! FATAL( packcoil, .not. allocated(FouCoil), illegal )
  FATAL( packcoil02, .not. allocated(DoF)    , illegal )

  do icoil = 1, Ncoils

     select case (coil(icoil)%type)
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
        FATAL( packcoil03 , idof .ne. DoF(icoil)%ND, counting error in packing )
        
     !---------------------------------------------------------------------------------------------
     case(2)
        idof = 0
        if(coil(icoil)%Lc /= 0) then
#ifdef dposition
           ! dipole position is variable
           DoF(icoil)%xdof(idof+1:idof+5) = (/ coil(icoil)%ox, coil(icoil)%oy, coil(icoil)%oz, &
                                               coil(icoil)%mt, coil(icoil)%mp /)
           idof = idof + 5      
#else
           DoF(icoil)%xdof(idof+1:idof+2) = (/ coil(icoil)%mt, coil(icoil)%mp /)
           idof = idof + 2     
#endif              
        endif
        FATAL( packcoil04 , idof .ne. DoF(icoil)%ND, counting error in packing )
     !---------------------------------------------------------------------------------------------
     case(3)
        idof = 0
        if(coil(icoil)%Lc /= 0) then
           DoF(icoil)%xdof(idof+1) = coil(icoil)%Bz; idof = idof + 1
        endif
        FATAL( packcoil05 , idof .ne. DoF(icoil)%ND, counting error in packing )
     !--------------------------------------------------------------------------------------------- 
     case(coil_type_spline)
        ! get number of DoF for each coil and allocate arrays;
        NCP = CPCoil(icoil)%NCP
        DoF(icoil)%xdof = zero

        !pack Fourier series;
        idof = 0
        if(coil(icoil)%Lc /= 0) then
           DoF(icoil)%xdof(idof+1 : idof+NCP*3) = CPCoil(icoil)%Cpoints(0:3*NCP-1); idof = idof + NCP*3    
        endif
        FATAL( packcoil03 , idof .ne. DoF(icoil)%ND, counting error in packing )
        
     

     !---------------------------------------------------------------------------------------------
     case default
        FATAL(packcoil06, .true., not supported coil types)
     end select     

  enddo ! end do icoil;

END SUBROUTINE packcoil

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

SUBROUTINE unpackcoil
  !---------------------------------------------------------------------------------------------
  ! pack coil representation variables into DoF (only geometries without currents);
  ! DATE: 2017/03/25
  !--------------------------------------------------------------------------------------------- 
  use globals, only: dp, zero, myid, ounit, case_coils, Ncoils, coil, FouCoil, DoF, MPI_COMM_FOCUS,coil_type_spline, CPCoil
  implicit none
  include "mpif.h"

  INTEGER  :: icoil, idof, NF, ierr, astat,NCP

  FATAL( unpackcoil01, .not. allocated(coil)   , illegal )
  ! FATAL( unpackcoil, .not. allocated(FouCoil), illegal )
  FATAL( unpackcoil02, .not. allocated(DoF)    , illegal )

  do icoil = 1, Ncoils

     select case (coil(icoil)%type)
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
        FATAL( unpackcoil03 , idof .ne. DoF(icoil)%ND, counting error in packing )

     !---------------------------------------------------------------------------------------------
     case(2)
        idof = 0
        if(coil(icoil)%Lc /= 0) then
#ifdef dposition
           ! dipole position is variable
           coil(icoil)%ox = DoF(icoil)%xdof(idof+1) ; idof = idof + 1
           coil(icoil)%oy = DoF(icoil)%xdof(idof+1) ; idof = idof + 1
           coil(icoil)%oz = DoF(icoil)%xdof(idof+1) ; idof = idof + 1
#endif
           coil(icoil)%mt = DoF(icoil)%xdof(idof+1) ; idof = idof + 1
           coil(icoil)%mp = DoF(icoil)%xdof(idof+1) ; idof = idof + 1
        endif
        FATAL( unpackcoil04 , idof .ne. DoF(icoil)%ND, counting error in packing )

     !---------------------------------------------------------------------------------------------
     case(3)
        idof = 0        
        
        if(coil(icoil)%Lc /= 0) then
           coil(icoil)%Bz =  DoF(icoil)%xdof(idof+1) ; idof = idof + 1        
        endif      
        FATAL( unpackcoil05 , idof .ne. DoF(icoil)%ND, counting error in packing )
     !--------------------------------------------------------------------------------------------- 
     case(coil_type_spline)
        ! get number of DoF for each coil and allocate arrays;
        NCP = CPCoil(icoil)%NCP
        idof = 0
        if (coil(icoil)%Lc /= 0) then
           !unpack Fourier series;
           CPCoil(icoil)%Cpoints(0:NCP*3-1) = DoF(icoil)%xdof(idof+1 : idof+3*NCP) ; idof = idof + 3*NCP 
        endif
        FATAL( unpackcoil03 , idof .ne. DoF(icoil)%ND, counting error in packing )


     !---------------------------------------------------------------------------------------------
     case default
        FATAL( unpackcoil06 , .true., not supported coil types)
     end select

  enddo ! end do icoil;

END SUBROUTINE unpackcoil
