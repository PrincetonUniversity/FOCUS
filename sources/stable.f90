! This is the overall function to handle stable field lines
SUBROUTINE stable
  use globals, only : dp, myid, ounit, input_stable, gsurf, MPI_COMM_FOCUS
  use mpi
  implicit none

  LOGICAL :: exist
  INTEGER :: iosta, astat, ierr

  allocate(gsurf(1:1))

  ! read the stable field lines
!  inquire( file=trim(input_stable), exist=exist)
!  FATAL( stable, .not.exist, input_stable does not exist )
!  call rdstable( input_stable, 1 )

  RETURN
END SUBROUTINE stable
