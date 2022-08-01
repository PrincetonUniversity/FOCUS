! This is the overall function to handle surfaces
SUBROUTINE surface
  use globals, only : dp, myid, ounit, machprec, surf, plasma, limiter, input_surf, limiter_surf, &
       psurf, weight_cssep, MPI_COMM_FOCUS, plasma_surf_boozer, case_surface
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

  if (case_surface == plasma_surf_boozer)   call rdbooz
  if (case_surface .ne. plasma_surf_boozer)   call fousurf( input_surf, plasma )

  ! read the limiter surface
  if (limiter /= plasma .and. case_surface == plasma_surf_boozer) then
     inquire( file=trim(limiter_surf), exist=exist)  
     FATAL( surface, .not.exist, limiter_surf does not exist )
     FATAL( surface, limiter <= plasma, something goes wrong the surface indexing )
     call rdbooz
  elseif (limiter /= plasma .and. case_surface .ne. plasma_surf_boozer) then
     inquire( file=trim(limiter_surf), exist=exist)  
     FATAL( surface, .not.exist, limiter_surf does not exist )
     FATAL( surface, limiter <= plasma, something goes wrong the surface indexing )
     call fousurf( limiter_surf, limiter )
  endif 

  RETURN
END SUBROUTINE surface
