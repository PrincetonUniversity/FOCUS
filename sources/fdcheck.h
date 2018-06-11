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
                     coil, xdof, Ndof, t1E, t2E, chi
                     
  implicit none
  include "mpif.h"

  INTEGER, INTENT(in)  :: ideriv
  !--------------------------------------------------------------------------------------------

  INTEGER              :: astat, ierr, idof
  REAL                 :: tmp_xdof(1:Ndof), fd, negvalue, posvalue, diff, rdiff
  REAL                 :: start, finish
  REAL, parameter      :: small=1.0E-4
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
         if (diff >= sqrtmachprec) write(ounit, *) "----------suspicious unmatching-----------------------"
      endif
      
  enddo
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

SUBROUTINE fdcheck2( ideriv )
!------------------------------------------------------------------------------------------------------ 
! DATE:  04/03/2017;
! ideriv = 1 -> check the first  derivatives with finite difference;
! ideriv = 2 -> check the second derivatives with finite difference;
!------------------------------------------------------------------------------------------------------

  use globals, only: zero, half, machprec, sqrtmachprec, ncpu, myid, ounit, &
                     coil, xdof, Ndof, t1E, t2E, chi, Bmnc, Bmns, wBmn, tBmnc, tBmns, NBmn, dB, Bmnim, Bmnin, Nzeta, Nteta, t1H, bharm, dofnorm
                     
  implicit none
  include "mpif.h"

  INTEGER, INTENT(in)  :: ideriv
  !--------------------------------------------------------------------------------------------

  INTEGER              :: astat, ierr, idof, imn
  REAL                 :: tmp_xdof(1:Ndof)
  REAL                 :: start, finish
  REAL, parameter      :: small=1.0E-4
  REAL, allocatable    :: dBc(:), dBs(:) 
  REAL                 :: fvec(2*Nbmn),fjac(2*Nbmn,Ndof),tmp(2*Nbmn), fd(2*Nbmn), diff(2*Nbmn)
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
  SALLOCATE( dBc, (1:NBmn), zero )  ! temporary dB_mn_cos
  SALLOCATE( dBs, (1:NBmn), zero )  ! temporary dB_mn_sin
  do idof = 1, Ndof
     call twodft( dB(idof,  0:Nteta-1, 0:Nzeta-1), dBs, dBc, Bmnim, Bmnin, NBmn )
     fjac(1:Nbmn, idof) = wBmn * dBc
     fjac(Nbmn+1:2*Nbmn, idof) = wBmn * dBs
     fjac(1:2*Nbmn, idof) = fjac(1:2*Nbmn, idof) * dofnorm(idof)
     t1H(idof) = sum( wBmn * ( (Bmnc - tBmnc)*dBc + (Bmns - tBmns)*dBs ) )
  enddo
  DALLOCATE( dBc )
  DALLOCATE( dBs )  
  
  if(myid .eq.0) write(ounit,'("fdcheck : bharm = "2(" ; "ES23.15))') chi, half * sum( wBmn * ((Bmnc - tBmnc)**2 + (Bmns - tBmns)**2) )
  
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
     tmp(1:Nbmn) = wBmn * (Bmnc - tBmnc)
     tmp(Nbmn+1:2*Nbmn) = wBmn * (Bmns - tBmns)

     !forward pertubation;
     tmp_xdof = xdof
     tmp_xdof(idof) = tmp_xdof(idof) + half * small
     call unpacking(tmp_xdof)
     call costfun(0)
     fvec(1:Nbmn) = wBmn * (Bmnc - tBmnc)
     fvec(Nbmn+1:2*Nbmn) = wBmn * (Bmns - tBmns)
     
     !finite difference;
     fd = (fvec - tmp) / small
     diff = abs(fjac(1:2*Nbmn, idof) - fd)

     if( myid.eq.0 ) then 
        do imn = 1, 2*Nbmn
           write(ounit,'("fdcheck : ", I6, "/", I6, 4(" ; "ES23.15))') idof, imn, fjac(imn, idof), fd(imn), diff(imn), diff(imn)/fd(imn)
           if (diff(imn) >= sqrtmachprec) write(ounit, *) "----------suspicious unmatching-----------------------"
        enddo
     endif
      
  enddo

  do idof = 1, Ndof
     if( myid.eq.0 ) then 
        write(ounit,'("fdcheck : ", I6, "/", I6, 2(" ; "ES23.15))') idof, Ndof, t1H(idof), t1E(idof)
     endif
  enddo
  !--------------------------------------------------------------------------------------------

  !case( -2 )
  case default
     FATAL(fdcheck, .true., not supported ideriv)
  !--------------------------------------------------------------------------------------------

  end select
  !--------------------------------------------------------------------------------------------
  return
END SUBROUTINE fdcheck2
