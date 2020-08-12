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

  use globals, only: dp, zero, half, machprec, sqrtmachprec, ncpu, myid, ounit, MPI_COMM_FAMUS, &
                     coil, xdof, Ndof, t1E, t2E, chi, LM_maxiter, LM_fvec, LM_fjac
                     
  implicit none
  include "mpif.h"

  INTEGER, INTENT(in)  :: ideriv
  !--------------------------------------------------------------------------------------------
  INTEGER              :: astat, ierr, idof, ivec, imax
  REAL                 :: tmp_xdof(1:Ndof), fd, negvalue, posvalue, diff, rdiff
  REAL                 :: start, finish, maxdiff, maxrdiff, small
  REAL, parameter      :: psmall=1.0E-4
  !--------------------------------------------------------------------------------------------

  maxdiff = zero ; maxrdiff = zero ; imax = 0
  if(myid == 0) write(ounit, *) "-----------Checking derivatives------------------------------"
  FATAL( fdcheck, Ndof < 1, No enough DOFs )
  !--------------------------------------------------------------------------------------------

  select case (ideriv)
  !--------------------------------------------------------------------------------------------

  case( -1 )
  if (myid == 0) then 
     write(ounit,'("fdcheck : Checking the first derivatives using finite-difference method.")')
     write(ounit,'(8X": Relative perturbation magnitude: delta = "ES12.5)') psmall
  end if

  call cpu_time(start)
  call costfun(1)
  call cpu_time(finish)
  if(myid .eq. 0) write(ounit,'("fdcheck : First order derivatives of energy function takes " &
       ES23.15 " seconds.")') finish - start
  if( myid.eq.0 ) write(ounit,'("fdcheck : idof/Ndof", 5(" ; ", A15))') "magnitude", "analytical",  &
   &  "fd-method", "abs diff", "relative diff"

  do idof = 1, Ndof
     ! perturbation will be relative.
     small = xdof(idof) * psmall
     if (small<machprec) small = psmall
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
         rdiff = diff / abs(fd)
     endif

     if( myid.eq.0 ) then 
         write(ounit,'("fdcheck : ", I4, "/", I4, 5(" ; "ES15.7))') idof, Ndof, small, t1E(idof), fd, diff, rdiff
         if (diff >= psmall**2) write(ounit, *) "----------suspicious unmatching-----------------------"
      endif
      
      ! get the maximum difference
      if (diff > maxdiff) then
         imax = idof
         maxdiff = diff
         maxrdiff = rdiff
      end if
  enddo

  if (myid.eq.0) write(ounit, '(8X": Max. difference: ", ES12.5, "; relative diff: ", ES12.5, "; at i="I6," .")') maxdiff, maxrdiff, imax
  ! return errorcode 1 
  if (maxrdiff > psmall) call MPI_ABORT( MPI_COMM_FAMUS, 1, ierr )

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
           if (diff >= sqrtmachprec) write(ounit, *) "----------suspicious unmatching-----------------------"
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
