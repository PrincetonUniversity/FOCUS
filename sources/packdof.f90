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
  use globals, only : dp, zero, myid, ounit, MPI_COMM_FAMUS, &
                    & case_coils, Ncoils, coil, DoF, Ndof, DoFnorm, dof_offset, ldof
  implicit none
  include "mpif.h"

  INTEGER :: lxdof(1:Ndof)
  INTEGER :: idof, icoil, ND, astat, ierr
  !--------------------------------------------------------------------------------------------- 

  ! reset xdof;
  lxdof = 0

  idof = dof_offset
  do icoil = 1, Ncoils
     if(coil(icoil)%Lc /= 0) then
        lxdof(idof+1) = coil(icoil)%ang
        idof = idof + 1
      endif
  enddo !end do icoil;

  !--------------------------------------------------------------------------------------------- 
  FATAL( packdof02 , idof-dof_offset .ne. ldof, counting error in packing )

  call MPI_ALLREDUCE( MPI_IN_PLACE, lxdof, Ndof, MPI_INTEGER, MPI_SUM, MPI_COMM_FAMUS, ierr )
  call mpi_barrier(MPI_COMM_FAMUS, ierr)

  return
END SUBROUTINE packdof

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

SUBROUTINE unpacking(lxdof)
  !---------------------------------------------------------------------------------------------
  ! UnPack all DOF from one vector;
  ! DATE: 2017/04/03
  !--------------------------------------------------------------------------------------------- 
  use globals, only: dp, zero, myid, ounit, MPI_COMM_FAMUS, &
       & case_coils, Ncoils, coil, DoF, Ndof, DoFnorm, dof_offset, ldof, momentq
  implicit none
  include "mpif.h"

  INTEGER :: lxdof(1:Ndof)
  INTEGER :: idof, icoil, ND, astat, ierr, ifirst
  !---------------------------------------------------------------------------------------------
  !FATAL( packdof, .not. allocated(xdof), Please allocate xdof first. )

  idof = dof_offset ; ifirst = 0
  do icoil = 1, Ncoils
   if(coil(icoil)%Lc /= 0) then
      coil(icoil)%ang = lxdof(idof+1)
      idof = idof + 1
    endif
  enddo !end do icoil;

  !--------------------------------------------------------------------------------------------- 
  FATAL( unpacking02 , idof-dof_offset .ne. ldof, counting error in unpacking )

  call discoil(ifirst)

  call mpi_barrier(MPI_COMM_FAMUS, ierr)

  return
END SUBROUTINE unpacking

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

