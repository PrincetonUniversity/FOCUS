!!$module minvol_mod
!!$  ! contains some common variables used in subroutine minvol
!!$  ! allocating once and re-using them will save allocation time
!!$  use globals, only : dp
!!$  implicit none
!!$
!!$  ! 0-order
!!$  REAL, allocatable :: 
!!$  ! 1st-order
!!$  REAL, allocatable :: 
!!$
!!$end module minvol_mod

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! This routine calculates the volume of dipoles (L-1 norm)
SUBROUTINE total_pm(ideriv)
  use globals, only: dp, zero, ncpu, myid, ounit, Nfp, &
       pmvol, t1U, coil, Ndof, Ncoils, DoF, total_moment, dof_offset, ldof, momentq, MPI_COMM_FAMUS
  implicit none
  include "mpif.h"
  INTEGER, INTENT(in) :: ideriv

  INTEGER             :: astat, ierr, icoil, idof, ND, ivec
  REAL                :: d1L(1:Ndof), norm(1:Ncoils), mag, dVdp

  pmvol = zero
  !-------------------------------calculate pmvol-------------------------------------------------- 
  if( ideriv >= 0 ) then
     do icoil = 1, Ncoils
        if ( coil(icoil)%Ic /= 0 ) then !if current is free;
           if (coil(icoil)%itype == 2) then
              mag = coil(icoil)%moment * abs( coil(icoil)%pho ) ** momentq
              if (coil(icoil)%symmetry == 0) then ! no symmetries
                 pmvol = pmvol + mag
              else if (coil(icoil)%symmetry == 1) then ! periodicity
                 pmvol = pmvol + mag*Nfp
              else if (coil(icoil)%symmetry == 2) then ! stellarator symmetry
                 pmvol = pmvol + mag*Nfp*2
              else
                 FATAL( total_pm01, .true., unspoorted symmetry option )
              end if 
           endif
        endif
     enddo
     call MPI_ALLREDUCE( MPI_IN_PLACE, pmvol, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FAMUS, ierr )
  endif
  !-------------------------------calculate d pmvol / d pho-------------------------------------------------- 
  if ( ideriv >= 1 ) then
     t1U = zero
     idof = dof_offset
     do icoil = 1, Ncoils
        ND = DoF(icoil)%ND
        if ( coil(icoil)%Ic /= 0 ) then !if current is free;
           idof = idof +1
           if (coil(icoil)%itype == 2) then
              dVdp = coil(icoil)%moment * SIGN(1.0_dp, coil(icoil)%pho) & 
                     & * momentq * abs(coil(icoil)%pho)**(momentq-1) 
              if (coil(icoil)%symmetry == 0) then ! no symmetries
                 t1U(idof) = dVdp
              else if (coil(icoil)%symmetry == 1) then ! periodicity
                 t1U(idof) = dVdp * Nfp
              else if (coil(icoil)%symmetry == 2) then ! stellarator symmetry
                 t1U(idof) = dVdp * Nfp*2
              else
                 FATAL( total_pm02, .true., unspoorted symmetry option )
              end if 
           else 
              t1U(idof) = zero
           endif
        endif
        if ( coil(icoil)%Lc /= 0 ) then !if geometry is free;
           t1U(idof+1:idof+ND) = zero
           idof = idof + ND
        endif
     enddo !end icoil;
     FATAL( total_pm , idof-dof_offset .ne. ldof, counting error in packing )
     call MPI_ALLREDUCE( MPI_IN_PLACE, t1U, Ndof, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FAMUS, ierr )
  endif

  call MPI_barrier( MPI_COMM_FAMUS, ierr )
  
  return
END SUBROUTINE total_pm
