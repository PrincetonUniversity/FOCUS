subroutine discrete(ideriv)
   use globals, only: dp, zero, pi, ncpu, myid, ounit, Nfp, MPI_COMM_FAMUS, &
                      disor, t1o, coil, Ndof, Ncoils, DoF, dof_offset, ldof, &
                      momentq
   use mpi
   IMPLICIT NONE

   INTEGER, INTENT(in) :: ideriv
   INTEGER :: astat, ierr, icoil, idof, ND
   REAL :: theta1, phi1, phi2, ldis, &
           theta_term, theta_deriv, phi_term, phi_deriv, dtheta, dphi

   disor = zero
   !-------------------------calculate discrete orientation--------------------------------------------
   if (ideriv >= 0) then
      do icoil = 1, Ncoils
         if (coil(icoil)%Lc /= 0) then !if orientation is free;
            if (coil(icoil)%itype == 2) then
               theta1 = pi/2
               phi1 = atan(coil(icoil)%oy/coil(icoil)%ox)
               phi2 = phi1 + pi/2
               ldis = abs(coil(icoil)%mt)*abs(coil(icoil)%mt - theta1) &
                      + abs(coil(icoil)%mt)*(abs(coil(icoil)%mp - phi1)*abs(coil(icoil)%mp - phi2))
               ldis = ldis * coil(icoil)%pho ** momentq
               if (coil(icoil)%symmetry == 0) then ! no symmetries
                  continue
               else if (coil(icoil)%symmetry == 1) then ! periodicity
                  ldis = ldis*Nfp
               else if (coil(icoil)%symmetry == 2) then ! stellarator symmetry
                  ldis = ldis*Nfp*2
               else
                  FATAL(discrete01, .true., unspoorted symmetry option)
               end if
               disor = disor + ldis
            endif
         endif
      enddo
      call MPI_ALLREDUCE(MPI_IN_PLACE, disor, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FAMUS, ierr)
   endif

   !-------------------------calculate first derivatives--------------------------------------------
   if (ideriv >= 1) then
      t1o = zero
      idof = dof_offset
      do icoil = 1, Ncoils
         if (coil(icoil)%Ic /= 0) then !if current is free;
            ldis = momentq * coil(icoil)%pho ** (momentq-1) &
                  * (abs(coil(icoil)%mt)*abs(coil(icoil)%mt - theta1) &
                  + abs(coil(icoil)%mt)*(abs(coil(icoil)%mp - phi1)*abs(coil(icoil)%mp - phi2)))
            if (coil(icoil)%symmetry == 0) then ! no symmetries
               continue
            else if (coil(icoil)%symmetry == 1) then ! periodicity
               ldis = ldis*Nfp
            else if (coil(icoil)%symmetry == 2) then ! stellarator symmetry
               ldis = ldis*Nfp*2
            else
               FATAL(discrete01, .true., unspoorted symmetry option)
            end if
            idof = idof + 1
            t1O(idof) = ldis
         endif
         if (coil(icoil)%Lc /= 0) then !if orientation is free;
            ND = DoF(icoil)%ND
            if (coil(icoil)%itype == 2) then
               theta1 = pi/2
               phi1 = atan(coil(icoil)%oy/coil(icoil)%ox)
               phi2 = pi/2
               theta_term = abs(coil(icoil)%mt)*abs(coil(icoil)%mt - theta1)
               theta_deriv = 2*coil(icoil)%mt - theta1
               phi_term = (abs(coil(icoil)%mp - phi1)*abs(coil(icoil)%mp - phi2))
               phi_deriv = abs(coil(icoil)%mt)*(2*coil(icoil)%mp - (phi1 + phi2))
               ! calculate d f_o / d theta
               if (coil(icoil)%mt >= theta1) then
                  dtheta = theta_deriv + phi_term
               else if (coil(icoil)%mt >= 0 .and. coil(icoil)%mt < theta1) then
                  dtheta = -theta_deriv + phi_term
               else
                  dtheta = theta_deriv - phi_term
               endif
               ! calculate d f_o / d theta
               if (coil(icoil)%mp >= phi1 .and. coil(icoil)%mp < phi2) then
                  dphi = -phi_deriv
               else
                  dphi = phi_deriv
               endif
               if (coil(icoil)%symmetry == 0) then ! no symmetries
                  continue
               else if (coil(icoil)%symmetry == 1) then ! periodicity
                  dtheta = dtheta*Nfp
                  dphi = dphi*Nfp
               else if (coil(icoil)%symmetry == 2) then ! stellarator symmetry
                  dtheta = dtheta*Nfp*2
                  dphi = dphi*Nfp*2
               else
                  FATAL(discrete02, .true., unspoorted symmetry option)
               end if
               idof = idof + 1; t1o(idof) = dtheta * coil(icoil)%pho ** momentq
               idof = idof + 1; t1o(idof) = dphi * coil(icoil)%pho ** momentq
            endif
         endif
      enddo
      FATAL(discrete03, idof - dof_offset .ne. ldof, counting error in packing)
      call MPI_ALLREDUCE(MPI_IN_PLACE, t1o, Ndof, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FAMUS, ierr)
   endif
end subroutine
