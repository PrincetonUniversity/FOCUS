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
  ! purpose to encourage binary values for dipole strength
  use globals, only: dp, zero, ncpu, myid, ounit, Nfp, &
       pmsum, dpbin, t1V, t1D, coil, Ndof, Ncoils, DoF, total_moment, dof_offset, ldof, momentq
  implicit none
  include "mpif.h"
  INTEGER, INTENT(in) :: ideriv

  INTEGER             :: astat, ierr, icoil, idof, ND, ivec
  REAL                :: d1L(1:Ndof), norm(1:Ncoils), rho, chi_d, dchi_d

  dpbin = zero
  !-------------------------------calculate dpbin-------------------------------------------------- 
  if( ideriv >= 0 ) then
     do icoil = 1, Ncoils
        if ( coil(icoil)%Ic /= 0 ) then !if current is free;
           if (coil(icoil)%itype == 2) then
              rho = coil(icoil)%pho ** momentq
              chi_d = ABS(rho) * ( 1 - ABS(rho) )
              if (coil(icoil)%symmetry == 0) then ! no symmetries
                 !pmsum = pmsum + coil(icoil)%I*coil(icoil)%I
                 dpbin = dpbin + chi_d**2
              else if (coil(icoil)%symmetry == 1) then ! periodicity
                 !pmsum = pmsum + coil(icoil)%I*coil(icoil)%I*Nfp
                 dpbin = dpbin + chi_d**2 * Nfp
              else if (coil(icoil)%symmetry == 2) then ! stellarator symmetry
                 !pmsum = pmsum + coil(icoil)%I*coil(icoil)%I*Nfp*2
                 dpbin = dpbin + chi_d**2 * Nfp * 2
              else
                 FATAL( bindip01, .true., unspoorted symmetry option )
              end if 
           endif
        endif
     enddo
     call MPI_ALLREDUCE( MPI_IN_PLACE, dpbin, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )
     !pmsum = pmsum / total_moment
     dpbin = dpbin / Ndof ! assumes no other dof
  endif
  !-------------------------------calculate d dpbin / d pho-------------------------------------------------- 
  if ( ideriv >= 1 ) then
     t1D = zero
     idof = dof_offset
     do icoil = 1, Ncoils
        ND = DoF(icoil)%ND
        if ( coil(icoil)%Ic /= 0 ) then !if current is free;
           idof = idof +1
           if (coil(icoil)%itype == 2) then
              rho = coil(icoil)%pho ** momentq
              dchi_d = 2 * rho * (1 - abs(rho) ) * ( 1 - 2*abs(rho) ) &
                         & * momentq * abs(coil(icoil)%pho) ** (momentq - 1)
              if (coil(icoil)%symmetry == 0) then ! no symmetries
               !  t1V(idof) = 2*coil(icoil)%I*coil(icoil)%moment*momentq &
               !       &       *(coil(icAoil)%pho)**(momentq-1)
                 t1D(idof) = dchi_d
              else if (coil(icoil)%symmetry == 1) then ! periodicity
               !  t1V(idof) = 2*coil(icoil)%I*coil(icoil)%moment*momentq &
               !       &       *(coil(icoil)%pho)**(momentq-1)*Nfp
                 t1D(idof) = dchi_d * Nfp
              else if (coil(icoil)%symmetry == 2) then ! stellarator symmetry
               !  t1V(idof) = 2*coil(icoil)%I*coil(icoil)%moment*momentq &
               !       &       *(coil(icoil)%pho)**(momentq-1)*Nfp*2
                 t1D(idof) = dchi_d * Nfp*2
              else
                 FATAL( bindip02, .true., unspoorted symmetry option ) 
              end if 
           else 
              t1D(idof) = zero
           endif
        endif
        if ( coil(icoil)%Lc /= 0 ) then !if geometry is free;
           t1D(idof+1:idof+ND) = zero
           idof = idof + ND
        endif
     enddo !end icoil;
     FATAL( bindip , idof-dof_offset .ne. ldof, counting error in packing )
     call MPI_ALLREDUCE( MPI_IN_PLACE, t1D, Ndof, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )
     !t1V = t1V / total_moment
     t1D = t1D / Ndof
  endif

  call MPI_barrier( MPI_COMM_WORLD, ierr )
  
  return
END SUBROUTINE bindip
