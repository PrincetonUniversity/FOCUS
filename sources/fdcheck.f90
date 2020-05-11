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

  use globals, only: dp, zero, half, machprec, sqrtmachprec, ncpu, myid, ounit, MPI_COMM_FOCUS, &
                     coil, xdof, Ndof, t1E, t2E, chi, LM_maxiter, LM_fvec, LM_fjac
                     
  implicit none
  include "mpif.h"

  INTEGER, INTENT(in)  :: ideriv
  !--------------------------------------------------------------------------------------------

  INTEGER              :: astat, ierr, idof, ivec
  REAL                 :: tmp_xdof(1:Ndof), fd, negvalue, posvalue, diff, rdiff
  REAL                 :: start, finish
  REAL, parameter      :: small=1.0E-6
  !--------------------------------------------------------------------------------------------

  if(myid == 0) write(ounit, *) "-----------Checking derivatives------------------------------"
  FATAL( fdcheck, Ndof < 1, No enough DOFs )
  !--------------------------------------------------------------------------------------------

  select case (ideriv)
  !--------------------------------------------------------------------------------------------

  case( -1 )
  if(myid == 0) write(ounit,'("fdcheck : Checking the first derivatives using finite-difference method")')

  call cpu_time(start)
  call costfun(1)
  call cpu_time(finish)
  if(myid .eq. 0) write(ounit,'("fdcheck : First order derivatives of energy function takes " &
       ES23.15 " seconds.")') finish - start
  if( myid.eq.0 ) write(ounit,'("fdcheck : idof  /  Ndof ;   analytical value"5X" ;   fd-method value"6X &
       " ;   difference"11X" ;   relative diff")')

  do idof = 1, Ndof
     !backward pertubation;
     tmp_xdof = xdof
     tmp_xdof(idof) = tmp_xdof(idof) - half * small
     call unpacking(tmp_xdof)
     call costfun(0)
     negvalue = chi

     !forward pertubation;
     tmp_xdof = xdof
     tmp_xdof(idof) = tmp_xdof(idof) + half * small
     call unpacking(tmp_xdof)
     call costfun(0)
     posvalue = chi
     
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
         if (diff >= small**2) write(ounit, *) "----------suspicious unmatching-----------------------"
      endif
      
  enddo

  ! L-M format
  if (LM_maxiter > 0) then
     ivec = 1 ! the evaluation term
     if(myid .eq. 0) write(ounit,'("fdcheck : check the derivatives of the ", I0, "-th term in L-M format.")') ivec
     if( myid.eq.0 ) write(ounit,'("fdcheck : idof  /  Ndof ;   analytical value"5X" ;   fd-method value"6X &
          " ;   difference"11X" ;   relative diff")')  
  
     do idof = 1, Ndof
        !backward pertubation;
        tmp_xdof = xdof
        tmp_xdof(idof) = tmp_xdof(idof) - half * small
        call unpacking(tmp_xdof)
        call costfun(0)
        negvalue = LM_fvec(ivec)

        !forward pertubation;
        tmp_xdof = xdof
        tmp_xdof(idof) = tmp_xdof(idof) + half * small
        call unpacking(tmp_xdof)
        call costfun(0)
        posvalue = LM_fvec(ivec)

        !finite difference;
        fd = (posvalue - negvalue) / small
        diff = abs(LM_fjac(ivec, idof) - fd)
        !output;
        if( abs(fd) < machprec ) then
           rdiff = 0
        else
           rdiff = diff / fd
        endif

        if( myid.eq.0 ) then 
           write(ounit,'("fdcheck : ", I6, "/", I6, 4(" ; "ES23.15))') idof, Ndof, LM_fjac(ivec, idof), fd, diff, rdiff
           if (diff >= small**2) write(ounit, *) "----------suspicious unmatching-----------------------"
        endif

     enddo

  endif
  
  !--------------------------------------------------------------------------------------------

  !case( -2 )
  case default
     FATAL(fdcheck, .true., not supported ideriv)
  !--------------------------------------------------------------------------------------------

  end select
  !--------------------------------------------------------------------------------------------
  return
END SUBROUTINE fdcheck

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
