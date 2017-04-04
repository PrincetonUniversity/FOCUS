
!title (main) ! Main program.

!latex \briefly{Main program.}

!l tex \calledby{\link{}}
!latex \calls{\link{al00aa}, \link{rwsurf}, \link{surface}, \link{rdknot}, \link{rdcoil}, \link{restart}, \link{setdof}, \link{knotxx}, \oculus{pp00aa}}

!latex \tableofcontents

!latex \subsection{tem}
!latex This is for temperary using. I will rewrite the whole document someday in the future. \\
!latex 2016/08/04;  write down the coordinate used in the code. The $\theta$ direction is clockwise. Be careful with this. I will add a subroutine in \link{rdcoils}
!latex to test if the coils file is using the same coordinate.
!latex \subsection{overview}
!latex \bi
!latex \item[1.] First, \link{al00aa} is called to read the input file.
!latex \item[2a.] If \inputvar{Itopology} $= 0$, construct unknotatrons;
!latex            the winding surface and plasma boundary are determined by \link{rwsurf} and \link{surface}.
!latex \item[2b.] If \inputvar{Itopology} $= 1$, construct arbitrary-knotatrons; \link{rdknot} is called to determine the knot;
!latex \item[3.] Then, the geometry of the coils are read from file, \link{rdcoil}.
!latex \item[4.] If \inputvar{Loptimize} $\ne 0$, then 
!latex \bi
!latex \item[i.] the coil degrees-of-freedom are assigned by \link{setdof};
!latex \item[ii.] the coil optimization proceeds using \nag{http://www.nag.co.uk/numeric/FL/manual19/pdf/E04/e04jyf_fl19.pdf}{E04JYF}, 
!latex           which iteratively calls \link{ofunct}.
!latex \ei
!latex \item[5.] If \inputvar{Lpoincare} $\ne 0$, then 
!latex \bi
!latex \item[i.] a \Poincare plot is constructed using \oculus{pp00aa}.
!latex \item[ii.] if \inputvar{Lpoincare=2}, input information for free-boundary SPEC vacuum verfication calculation is provided by \oculus{bn00aa};
!latex \ei
!latex \ei

!latex \subsection{Oculus}
!latex \bi
!latex \item[1.] Subroutines from the Oculus library are used, see e.g. \oculus{bs00aa} and \oculus{pp00aa}.
!latex \ei

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

PROGRAM focus
  
  use kmodule, only : ncpu, myid, ounit, Itopology, Loptimize, Lnormalize, Lpoincare, tstart, tfinish
  
  use oculus, only : pp00aa
  
  implicit none
  
  include "mpif.h"

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  INTEGER :: ierr, astat, irestart

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  myid = 0 ; ncpu = 1
  
  call MPI_INIT( ierr )
  call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
  call MPI_COMM_SIZE( MPI_COMM_WORLD, ncpu, ierr )
  
  if( myid.eq.0 ) write(ounit,'("focus   : ", 10x ," : ncpu =",i4," ; begin execution ;")') ncpu
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  call initial ! read input and broadcast; 04 Jun 15;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  select case( Itopology )
  case( 0 ) ; call surface
  case( 1 ) ; call rdknot
  end select
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  call rdcoils
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( Lnormalize .ne. 0 ) call weightnormalize
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  call cpu_time(tstart)
  
  select case( Loptimize )
  case( -2 )   ; call testderivs2 ; irestart=0  ! finite-difference check of the second derivatives  ;
  case( -1 )   ; call testderivs1 ; irestart=0  ! finite-difference check of the first  derivatives  ;
! case( -2 )   ; call NAGderiv2   ; irestart=0  ! NAG routine check of the second derivatives        ;
! case( -1 )   ; call NAGderiv1   ; irestart=0  ! NAG routine check of the first  derivatives        ;
  case(  1 )   ; call descent     ; irestart=1  ! differential flow; shudson            ; 14 Apr 16  ;
  case(  2 )   ; call descent2    ; irestart=1  ! differential flow; czhu               ; 14 Apr 16  ;
  case(  3 )   ; call Powell      ; irestart=1  ! Powell nonlinear solver; NAG C05PDF; with    derivs;
  case(  4 )   ; call Newton      ; irestart=1  ! Newton minimum finder  ; NAG E04LYF; with    derivs;
  case(  9 )   ; call SVD         ; irestart=1  ! Analyze current Heissian matrix using SVD; F08KBF  ;
  case default ;                  ; irestart=1
  end select

  !call test

  call cpu_time(tfinish)
  
  if( myid .eq. 0 ) write(ounit,'("focus   : ",10x," : optimization finished ; totally takes ", es23.15, " seconds;")') tfinish - tstart

  call identfy

  call restart( irestart )
  

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  ! binary selection;
  ! 0 : off
  ! 1 : diagnos; 2: exinput; 4: wtmgrid
  select case( Lpoincare )
  case( 0  ) ; 
  case( 1  ) ; call diagnos ! normal poicare plot from coils file 
  case( 2  ) ; call exinput ! call Oculus:bn00aa to prepare input for free-boundary SPEC vacuum verfication; 04 Aug 16;
  case( 3  ) ; call exinput
             ; call diagnos
  case( 4  ) ; call wtmgrid !write mgrid
  case( 5  ) ; call wtmgrid
             ; call diagnos
  case( 6  ) ; call exinput
             ; call wtmgrid
  case( 7  ) ; call exinput
             ; call wtmgrid
             ; call diagnos
  case( 8: ) ; call wtmgrid ! 8 only write mgrid; > 8 write both
  end select

 !call identfy

  call cpu_time(tstart)
  
  if( myid .eq. 0 ) write(ounit,'("focus   : ",10x," : post-proceed finished ; totally takes ", es23.15, " seconds;")') tstart - tfinish

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  call MPI_FINALIZE( ierr )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
END PROGRAM focus

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!


SUBROUTINE test

!just for benchmarking with the new code; 04/04/2017
  
  use kmodule, only: coil, surf, NFcoil, Ncoils, Nteta, Nzeta, Ndof, bnorm, t1B, Cdof, &
                     zero, pi2, myid, ounit, NDcoil
  implicit none
  include "mpif.h"

  INTEGER :: icoil, iteta, jzeta, ierr, astat
  
  !Cdof = 6*NFcoil + 3
  icoil = 1 ; iteta = 0; jzeta = 0

  if(.not. allocated(coil(icoil)%Bx)) then 
     SALLOCATE( coil(icoil)%Bx, (0:Cdof, 0:Cdof), zero )
     SALLOCATE( coil(icoil)%By, (0:Cdof, 0:Cdof), zero )
     SALLOCATE( coil(icoil)%Bz, (0:Cdof, 0:Cdof), zero )
  endif
  
     call bfield1(icoil, iteta, jzeta, coil(icoil)%Bx(0:Cdof,0), coil(icoil)%By(0:Cdof,0), &
          coil(icoil)%Bz(0:Cdof,0))
    
     write(ounit,*) "Bx/x = ", coil(icoil)%Bx(0:Cdof,0) !* pi2/NDcoil
     write(ounit,*) "By/x = ", coil(icoil)%By(0:Cdof,0) !* pi2/NDcoil
     write(ounit,*) "Bz/x = ", coil(icoil)%Bz(0:Cdof,0) !* pi2/NDcoil

  !endif

  call bnormal2(0)
  write(ounit,*) "bnorm = ", bnorm
  !call bnormal(1)
  !write(ounit,*) "OB/x = ", t1B(1:Ncoils, 0:Cdof)

end SUBROUTINE test
