! This is the overall function to handle surfaces
SUBROUTINE surface
  use globals, only : dp, myid, ounit, machprec, surf, plasma, limiter, input_surf, limiter_surf, &
       psurf, weight_cssep, MPI_COMM_FOCUS,plasma_surf_fourier, plasma_surf_knot, plasma_surf_boozer, &
       plasma_surf_hdf5, case_surface
  use mpi
  implicit none

  LOGICAL :: exist
  INTEGER :: iosta, astat, ierr

  ! determine the total number of surfaces
  ! if ( weight_cssep > machprec .and. trim(limiter_surf) /= trim(input_surf) ) then 
  if ( weight_cssep > machprec ) then     
     plasma = 1
     limiter = 2
     if ( limiter_surf .eq. input_surf ) limiter = plasma ! use the plasma surface as limiter surface
  else ! use the plasma surface as limiter
     plasma = 1
     limiter = 1
  endif
  allocate(surf(plasma:limiter))
  psurf = limiter

  ! read the plasma surface  
  inquire( file=trim(input_surf), exist=exist)
  FATAL( surface, .not.exist, input_surf does not exist )

  select case (case_surface)
  case(plasma_surf_fourier);  call fousurf( input_surf, plasma )
  case(plasma_surf_knot);     call rdknot
  case(plasma_surf_boozer);   call rdbooz( input_surf, plasma )
  case(plasma_surf_hdf5);     call rdhdf5( input_surf, plasma )
  ! read wout option missed
  end select 

  ! read the limiter surface
  if (limiter /= plasma) then
      FATAL( surface, limiter <= plasma, something goes wrong the surface indexing )
      select case (case_surface)
      case(plasma_surf_fourier);  call fousurf( input_surf, limiter )
      case(plasma_surf_knot);     FATAL(surface, case_surface==plasma_surf_knot, limiter surface not supported in knot option)
      case(plasma_surf_boozer);   call rdbooz( input_surf, limiter )
      case(plasma_surf_hdf5);     call rdhdf5( input_surf, limiter )
      end select   
  endif 

  RETURN
END SUBROUTINE surface
