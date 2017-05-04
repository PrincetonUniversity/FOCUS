!title (fd check) ! Read input file, initialize global variables.

!latex \briefly{Check the derivatives using finite difference method}

!latex \calledby{\link{knotopt}}
!latex \calls{\link{}}

!latex \tableofcontents

!latex \subsection{overview}
!latex \bi
!latex \item[1.] 
!latex \ei
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

SUBROUTINE fdcheck( ideriv )
!------------------------------------------------------------------------------------------------------ 
! DATE:  04/03/2017;
! ideriv = 1 -> check the first  derivatives with finite difference;
! ideriv = 2 -> check the second derivatives with finite difference;
!------------------------------------------------------------------------------------------------------

  use globals, only: zero, half, machprec, sqrtmachprec, ncpu, myid, ounit, &
                     coil, xdof, Ndof, t1E, t2E, totalenergy, dofnorm
                     
  implicit none
  include "mpif.h"

  INTEGER, INTENT(in)  :: ideriv
  !--------------------------------------------------------------------------------------------

  INTEGER              :: astat, ierr, idof
  REAL                 :: tmp_xdof(1:Ndof), small=1.0E-6, fd, negvalue, posvalue, diff, rdiff
  REAL                 :: start, finish
  !--------------------------------------------------------------------------------------------

  if(myid == 0) write(ounit, *) "-----------Checking derivatives------------------------------"
  FATAL( fdcheck, Ndof < 1, No enough DOFs )
  tmp_xdof = xdof  !save the current xdof;
  !--------------------------------------------------------------------------------------------

  select case (ideriv)
  !--------------------------------------------------------------------------------------------

  case( 1 )
  if(myid == 0) write(ounit,'("fdcheck : Checking the first derivatives using finite-difference method")')

  call cpu_time(start)
  call costfun(1); t1E = t1E * dofnorm
  call cpu_time(finish)
  if(myid .eq. 0) write(ounit,'("fdcheck : First order derivatives of energy function takes " &
       ES23.15 " seconds.")') finish - start
  if( myid.eq.0 ) write(ounit,'("fdcheck : idof  /  Ndof ;   analytical value"5X" ;   fd-method value"6X &
       " ;   difference"11X" ;   reletive diff")')

  do idof = 1, Ndof
     !backward pertubation;
     xdof(idof) = xdof(idof) - half * small
     call unpacking(xdof)
     call costfun(0)
     negvalue = totalenergy
     xdof = tmp_xdof
     !forward pertubation;
     xdof(idof) = xdof(idof) + half * small
     call unpacking(xdof)
     call costfun(0)
     posvalue = totalenergy
     xdof = tmp_xdof
     !finite difference;
     fd = (posvalue - negvalue) / small
     diff = abs(t1E(idof) - fd)
     !output;
     if( abs(fd) < machprec ) then
         rdiff = 0
     else
         rdiff = diff / fd
     endif

     if( myid.eq.0 ) then 
         write(ounit,'("fdcheck : ", I6, "/", I6, 4(" ; "ES23.15))') idof, Ndof, t1E(idof), fd, diff, rdiff
         if (diff >= sqrtmachprec) write(ounit, *) "----------suspicious unmatching-----------------------"
      endif
      
  enddo
  !--------------------------------------------------------------------------------------------

  !case( 2 )
  case default
     FATAL(fdcheck, .true., not supported ideriv)
  !--------------------------------------------------------------------------------------------

  end select
  !--------------------------------------------------------------------------------------------
  return
END SUBROUTINE fdcheck
