! This is the overall function to handle stable field lines
SUBROUTINE stable
  use globals, only : dp, myid, ounit, input_stable, gsurf, MPI_COMM_FOCUS, period_Nseg, resbn_m
  use mpi
  implicit none

  LOGICAL :: exist
  INTEGER :: iosta, astat, ierr, index, Nseg_stable

  allocate(gsurf(1:1))

  ! read the stable field lines
!  inquire( file=trim(input_stable), exist=exist)
!  FATAL( stable, .not.exist, input_stable does not exist )
!  call rdstable( input_stable, 1 )

  index = 1
  gsurf(index)%Nseg_stable = period_Nseg*resbn_m + 1
  Nseg_stable = gsurf(index)%Nseg_stable
  SALLOCATE( gsurf(index)%ox, (1:Nseg_stable), 0.0 )
  SALLOCATE( gsurf(index)%oy, (1:Nseg_stable), 0.0 )
  SALLOCATE( gsurf(index)%oz, (1:Nseg_stable), 0.0 )
  SALLOCATE( gsurf(index)%xx, (1:Nseg_stable), 0.0 )
  SALLOCATE( gsurf(index)%xy, (1:Nseg_stable), 0.0 )
  SALLOCATE( gsurf(index)%xz, (1:Nseg_stable), 0.0 )
  SALLOCATE( gsurf(index)%oxdot, (1:Nseg_stable), 0.0 )
  SALLOCATE( gsurf(index)%oydot, (1:Nseg_stable), 0.0 )
  SALLOCATE( gsurf(index)%ozdot, (1:Nseg_stable), 0.0 )
  SALLOCATE( gsurf(index)%xxdot, (1:Nseg_stable), 0.0 )
  SALLOCATE( gsurf(index)%xydot, (1:Nseg_stable), 0.0 )
  SALLOCATE( gsurf(index)%xzdot, (1:Nseg_stable), 0.0 )
  gsurf(1)%donee = 0

  RETURN
END SUBROUTINE stable
