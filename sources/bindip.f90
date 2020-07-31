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

SUBROUTINE bindip(ideriv)
  ! copied from minvol
  ! dpsum = \sum J = \sum |p|( 1 - |p| )
  ! meant to encourage binary values for dipole strength
  use globals, only: dp, zero, ncpu, myid, ounit, Nfp, &
       pmsum, dpbin, t1V, coil, Ndof, Ncoils, DoF, total_moment, dof_offset, ldof, momentq
  implicit none
  include "mpif.h"
  INTEGER, INTENT(in) :: ideriv

  INTEGER             :: astat, ierr, icoil, idof, ND, ivec
  REAL                :: d1L(1:Ndof), norm(1:Ncoils)

  pmsum = zero !
  dpbin = zero
  !-------------------------------calculate dpbin-------------------------------------------------- 
  if( ideriv >= 0 ) then
     do icoil = 1, Ncoils
        if ( coil(icoil)%Ic /= 0 ) then !if current is free;
           if (coil(icoil)%itype == 2) then
              if (coil(icoil)%symmetry == 0) then ! no symmetries
                 pmsum = pmsum + coil(icoil)%I*coil(icoil)%I
              else if (coil(icoil)%symmetry == 1) then ! periodicity
                 pmsum = pmsum + coil(icoil)%I*coil(icoil)%I*Nfp
              else if (coil(icoil)%symmetry == 2) then ! stellarator symmetry
                 pmsum = pmsum + coil(icoil)%I*coil(icoil)%I*Nfp*2
              else
                 FATAL( bindip01, .true., unspoorted symmetry option )
              end if 
           endif
        endif
     enddo
     call MPI_ALLREDUCE( MPI_IN_PLACE, pmsum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )
     pmsum = pmsum / total_moment
  endif
  !-------------------------------calculate d dpbin / d pho-------------------------------------------------- 
  if ( ideriv >= 1 ) then
     t1V = zero
     idof = dof_offset
     do icoil = 1, Ncoils
        ND = DoF(icoil)%ND
        if ( coil(icoil)%Ic /= 0 ) then !if current is free;
           idof = idof +1
           if (coil(icoil)%itype == 2) then
              if (coil(icoil)%symmetry == 0) then ! no symmetries
                 t1V(idof) = 2*coil(icoil)%I*coil(icoil)%moment*momentq &
                      &       *(coil(icoil)%pho)**(momentq-1)
              else if (coil(icoil)%symmetry == 1) then ! periodicity
                 t1V(idof) = 2*coil(icoil)%I*coil(icoil)%moment*momentq &
                      &       *(coil(icoil)%pho)**(momentq-1)*Nfp
              else if (coil(icoil)%symmetry == 2) then ! stellarator symmetry
                 t1V(idof) = 2*coil(icoil)%I*coil(icoil)%moment*momentq &
                      &       *(coil(icoil)%pho)**(momentq-1)*Nfp*2
              else
                 FATAL( bindip02, .true., unspoorted symmetry option ) 
                  ! Q: What does this do?
              end if 
              !t1V(idof) = coil(icoil)%moment*momentq*sin(coil(icoil)%pho)**(momentq-1)*cos(coil(icoil)%pho)
           else 
              t1V(idof) = zero
           endif
        endif
        if ( coil(icoil)%Lc /= 0 ) then !if geometry is free;
           t1V(idof+1:idof+ND) = zero
           idof = idof + ND
        endif
     enddo !end icoil;
     FATAL( bindip , idof-dof_offset .ne. ldof, counting error in packing )
     call MPI_ALLREDUCE( MPI_IN_PLACE, t1V, Ndof, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )
     t1V = t1V / total_moment
     !TMPOUT(t1V)
  endif

  call MPI_barrier( MPI_COMM_WORLD, ierr )
  
  return
END SUBROUTINE bindip
