
!title (truncnt) ! Truncated Newton method

!latex \briefly{Using a Truncated Newron algorithm \doilink{10.1145/128745.150973}{TNPACK} .}


SUBROUTINE truncnt

  use kmodule, only: zero, myid, ounit, Idisplay, lunit, Ndof, Ntauout, itau
  implicit none
  include "mpif.h"

  INTEGER              :: INFORM, N, NZ, LW, LIW, OPLIST(20)
  REAL                 :: x(1:Ndof), PLIST(20), F, G(1:Ndof)
  INTEGER, allocatable :: IW(:)
  REAL, allocatable    :: W(:)
  INTEGER              :: ierr, astat
  CHARACTER(20), parameter :: nulfilename = 'NUL'
  EXTERNAL             :: calfgh,calpat,calhdp, SETLIS, TNMIN


  if (myid .eq. 0) write(ounit, '("truncnt : "10X" : Begin using Truncated Newton Method to optimize.")')

  itau = 0
  N = Ndof
  NZ = N*(N+1)/2
  LW = 10*N+5*NZ
  LIW = 7*N+5*NZ + 1

  SALLOCATE(IW, (1:LIW), 0)
  SALLOCATE( W, (1: LW), zero)

  ! initial guess
  call costfun(2)
  call pack(x(1:N))
  call output

  CALL SETLIS(N,OPLIST,PLIST,INFORM)

  open(lunit,file=nulfilename)  ! create a null file to eliminate outputs from other codes;

  OPLIST(1) = 1
  OPLIST(2) = 0
  OPLIST(4) = Ntauout
  OPLIST(10) = 0
  OPLIST(11) = lunit   ! output to null file;
  if (myid .eq. 0 .and. Idisplay .le. 0) OPLIST(11) = ounit ! only root cpu output;
  OPLIST(12) = 1
  OPLIST(15) = 2
  OPLIST(16) = 2
  OPLIST(17) = 2

  PLIST(2) = 1.0D-08
  PLIST(3) = 0.1
  !PLIST(5) = 0.1
  PLIST(8) = 10.0D0

  CALL TNMIN(N,X,F,G,OPLIST,PLIST,INFORM,NZ,W,LW,IW,LIW,calfgh,calpat,calhdp)

  if (myid .eq. 0 .and. Idisplay .le. -1) WRITE(ounit, '("truncnt : INFORM = "I3)') INFORM

  close(lunit)

  return

END SUBROUTINE truncnt

!========================================================================

SUBROUTINE calfgh(N, X, F, G, A, IA,JA, NOUT)
  use kmodule, only : zero, myid, ounit, Ndof, totalenergy, t1E, t2E
  implicit none
  include "mpif.h"

  !-----dummy arguments-------------
  INTEGER, INTENT(in )   :: N, NOUT
  REAL,    INTENT(in )   :: X(1:N)
  
  REAL,    INTENT(out)   :: F, G(N)

  ! ONLY INPUT WHEN NOUT = 3 OR 5; OTHERWISE OUTPUT
  INTEGER                :: IA(1:N+1), JA(*)
  REAL                   :: A(*)


  !-----local variables------------
  INTEGER  :: ierr, nn, idof, jdof, c1, n1, c2, n2

  FATAL( truncnt, N /= Ndof, inconsistent dimension)
  call unpack(x(1:N))
  call discretecoil

  ! NOUT = 0; only F
  call costfun(0)
  F = totalenergy
  if (NOUT .eq. 0) return

  ! NOUT = 1; return F,G;
  call costfun(1)
  do idof = 1, N
     call DoFconvert(idof,c1,n1)
     G(idof) = t1E(c1,n1)
  enddo
  if (NOUT .eq. 1) return   
  
  ! NOUT = 2; return F, G, H;
  call costfun(2)

!!$  do idof = 1, n
!!$     call DoFconvert(idof,c1,n1)
!!$     HESD(idof) = t2E(c1,n1,c1,n1)  !diagnal terms
!!$  enddo
!!$
!!$  nn = 0
!!$  do idof = 1, N-1                    ! row;
!!$     call DoFconvert(idof,c1,n1)
!!$     do jdof = idof+1, N              ! column;
!!$        call DoFconvert(jdof,c2,n2)
!!$        nn = nn + 1
!!$        HESU( nn ) = t2E(c1,n1,c2,n2) ! upper tringular half;
!!$     enddo
!!$  enddo
!!$
!!$  FATAL( truncnt, nn .ne. N*(N-1)/2, counting error)
 
  if (NOUT .eq. 2) return

  ! NOUT = others; let M = tridiagnol H;
  nn = 0
  do idof = 1, N-1
     call DoFconvert(idof,c1,n1)
     do jdof = idof, idof + 1
        call DoFconvert(jdof,c2,n2)
        nn = nn + 1
        A(nn) = t2E(c1,n1,c2,n2)
     enddo
  enddo

  idof = N; nn = nn + 1
  call DoFconvert(idof,c1,n1)
  A(nn) = t2E(c1,n1,c1,n1)
     
  FATAL( truncnt, nn .ne. 2*N-1, counting error)

  return
END SUBROUTINE calfgh

!========================================================================

SUBROUTINE calhdp(N,D,HD,X,G)
  use kmodule, only: zero, myid, ounit, t2E
  implicit none
  include "mpif.h"

  !-----dummy arguments-------------
  INTEGER, INTENT(in )   :: N
  REAL,    INTENT(in )   :: D(1:N), X(1:N), G(1:N)
  
  REAL,    INTENT(out)   :: HD(1:N)

  !-----local variables------------
  INTEGER  :: ierr, ii, jj, c1, n1, c2, n2

  do ii = 1, N
     HD(ii) = zero
     call DoFconvert(ii,c1,n1)
     do jj = 1, N
        call DoFconvert(jj,c2,n2)
        HD(ii) = HD(ii) + t2E(c1, n1, c2, n2) * D(jj)
     enddo
  enddo

  return
END SUBROUTINE calhdp

!========================================================================

SUBROUTINE calpat(N,X,A,IA,JA)
  use kmodule, only: zero, myid, ounit
  implicit none
  include "mpif.h"

  !-----dummy arguments-------------
  INTEGER, INTENT(in )   :: N
  REAL,    INTENT(in )   :: X(1:N)

  INTEGER, INTENT(out)   :: IA(1:N+1), JA(*)

  REAL                   :: A(*)

  !-----local variables------------
  INTEGER :: nn, ii, jj, ierr

  nn = 0
  do ii = 1,N-1
     IA(ii) = 2*ii-1
     do jj = 0, 1
        nn = nn + 1
        JA(nn) = ii + jj
     enddo
  enddo

  ii = N; IA(ii) = 2*ii-1

  nn = nn + 1
  JA(nn) = N
  FATAL( truncnt, nn .ne. 2*N-1, counting error)

  IA(N+1) = 2*N
  return
END SUBROUTINE calpat
