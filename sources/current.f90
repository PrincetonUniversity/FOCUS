!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
subroutine current( ideriv )
!------------------------------------------------------------------------------------------------------ 
! Calculate the difference between the total FREE current and the target.
! Formula is isum = ((SUM(I_i) - target_isum) / target_isum)**2
! ideriv = 0 -> only calculate the total current constraint;
! ideriv = 1 -> calculate the toroidal flux constraint and its first derivatives;
! ideriv = 2 -> calculate the toroidal flux constraint and its first & second derivatives;
!------------------------------------------------------------------------------------------------------   
  use globals, only: dp, zero, half, one, pi2, sqrtmachprec, ncpu, myid, ounit, &
       coil, DoF, Ncoils, Cdof, Nfp, &
       isum, t1I, Ndof, weight_isum, target_isum, &
       iisum, misum, LM_fvec, LM_fjac, MPI_COMM_FOCUS
  use mpi
  implicit none

  INTEGER, INTENT(in)                   :: ideriv
  !--------------------------------------------------------------------------------------------
  INTEGER                               :: astat, ierr
  INTEGER                               :: icoil, idof, ND
  REAL                                  :: total_current
  REAL, dimension(1:Ndof)               :: current_array, symmetry
  !--------------------------initialize and allocate arrays------------------------------------- 
  isum = zero
  total_current = zero
  current_array = zero
  symmetry = zero
  !--------------------------function only ------------------------------------------------------ 
  if( ideriv >= 0 ) then
     do icoil = 1, Ncoils     
        if(coil(icoil)%Ic /= 0) then ! exclude fixed currents
           ! check if the coil is stellarator symmetric
            select case (coil(icoil)%symm) 
            case ( 0 )
                symmetry(icoil) = 1
            case ( 1 )
                symmetry(icoil)  = Nfp
            case ( 2) 
                symmetry(icoil)  = 2*Nfp
            end select
           !if( myid.ne.modulo(icoil-1,ncpu) ) cycle ! parallelization loop;
           current_array(icoil) = coil(icoil)%I
        endif
     enddo
     total_current = SUM(current_array * symmetry)
     isum = ((total_current - target_isum) / target_isum)**2

     if (misum > 0) then ! L-M format of targets
        LM_fvec(iisum+1) = weight_isum * ((total_current - target_isum) / target_isum)
     endif
  endif 

  !--------------------------function & 1st deriv. -----------------------------------------------

  if ( ideriv >= 1 ) then
     t1I = zero
     idof = 0 
     do icoil = 1, Ncoils
        ND = DoF(icoil)%ND
        if ( coil(icoil)%Ic /= 0 ) then !if current is free;
           idof = idof +1
           t1I(idof) = 2 * ((total_current - target_isum) / target_isum) * symmetry(icoil) / target_isum

           if (misum > 0) then ! L-M format of targets
              LM_fjac(iisum + 1, idof) = weight_isum * symmetry(icoil) / target_isum
           endif
        endif
        if ( coil(icoil)%Lc /= 0 ) then !if geometry is free;
           idof = idof + ND
        endif

     enddo ! end icoil;
     FATAL( current , idof .ne. Ndof, counting error in packing )
  endif

  !--------------------------------------------------------------------------------------------

  call MPI_barrier( MPI_COMM_FOCUS, ierr )

  return
end subroutine current

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
