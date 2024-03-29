
!title (main) ! Main program.

!latex \briefly{Main program.}

!latex \calledby{\link{}}
!latex \calls{\link{initial}, \link{rdsurf}, \link{rdcoils}, 
!latex \link{saving}, \link{solvers}}.

!latex \tableofcontents

!latex  \subsection{Brief introduction}
!latex  This code is used for designing 3-D coils in stellarators, knotrans and tokamaks. 
!latex  FOCUS (Finding/Flexible Optimized Coils Using Space curves) gets rid of winding surface 
!latex  by representing 
!latex  coils using space curves (either Fundamental theorem of space curves or Fourier series or other 
!latex  representations.). And for the first time, the derivatives (both the first and the second ones)
!latex  are analytically calculated. \par
!latex  For more information, please visti \href{https://princetonuniversity.github.io/FOCUS/}
!latex  If you have any questions, please send a email to czhu@pppl.gov (or zcxiang@mail.ustc.edu.cn).

!latex  \subsection{How to execute}
!latex  A brief help message will be printed if you just type `xfocus --help`

!latex \subsection{Misc}
!latex \bi
!latex \item Details about the macros used are listed in \emph{Macros};
!latex \item The unit used in FOCUS is the International System of Units (SI);
!latex \item All the shared variables are defined in \link{globals};
!latex \item Variables starting with {\bf case\_xxx} are used for controling the flow direction;
!latex \item Variables starting with {\bf Isxxx} are used for logical judgement.
!latex \ei

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
! module test_mod
!   use globals, only : dp, zero
!   implicit none

!   type testsurface
!   INTEGER              :: Nteta, Nzeta, Nfou=0, Nfp=0, NBnf=0
!   REAL   , allocatable :: Rbc(:), Zbs(:), Rbs(:), Zbc(:), Bnc(:), Bns(:)
!   REAL   , allocatable :: xx(:,:), yy(:,:), zz(:,:), nx(:,:), ny(:,:), nz(:,:), &
!                           xt(:,:), yt(:,:), zt(:,:), xp(:,:), yp(:,:), zp(:,:), &
!                           ds(:,:), bn(:,:), pb(:,:), &
!                           Bx(:,:), By(:,:), Bz(:,:)
!   INTEGER, allocatable :: bim(:), bin(:), Bnim(:), Bnin(:)
!   REAL                 :: vol, area
! end type testsurface


! type testFreedom
!   INTEGER              :: ND
!   REAL   , allocatable :: xdof(:), xof(:,:), yof(:,:), zof(:,:)
! end type testFreedom

!   type(testsurface), target :: test_surf(10)
!   type(testFreedom), target :: test_DoF(10)
  

! end module test_mod

subroutine focus()

  use globals, only: dp, ncpu, myid, ounit, ierr, astat, eunit, case_surface, case_coils, case_optimize, &
       case_postproc, xdof, time_initialize, time_optimize, time_postproc, &
       version, MPI_COMM_FOCUS, pi, machprec, sqrtmachprec, vsmall, small, ten, thousand
  use mpi  !to enable gfortran mpi_wtime bugs; 07/20/2017
  implicit none

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  INTEGER :: secs, mins, hrs
  REAL    :: tstart, tfinish ! local variables

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  !call MPI_COMM_RANK( MPI_COMM_FOCUS, myid, ierr )
  !call MPI_COMM_SIZE( MPI_COMM_FOCUS, ncpu, ierr )

  if(myid == 0) write(ounit, *) "---------------------  FOCUS ", version, "------------------------------"
  if(myid == 0) write(ounit,'("focus   : Begin execution with ncpu =",i5)') ncpu

  tstart =  MPI_WTIME()

  ! check input namelist 
  call check_input
  
  select case( case_surface )

  case( 0 ) ; call surface   ! general format (VMEC-like) plasma boundary;
  case( 1 ) ; call rdknot    ! knototran-like plasma boundary;
 !case( 2 ) ; call readwout  ! read vmec output for plasma boundary and Boozer coordinates; for future;

  end select
    
  select case( case_coils )

 !case( 0 )   ; call coilpwl ! piece-wise linear; for future;
  case( 1 )   ; call rdcoils

  end select

  call packdof(xdof)  ! packdof in xdof array;

  call MPI_BARRIER( MPI_COMM_FOCUS, ierr )

  tfinish = MPI_Wtime()
  time_initialize = tfinish - tstart
  if( myid  ==  0 ) write(ounit, '(A, ES12.5, A3)') "focus   : Initialization took ", time_initialize," S."

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if (case_optimize /= 0) call solvers       ! call different solvers;

  call MPI_BARRIER( MPI_COMM_FOCUS, ierr )

  tstart = MPI_Wtime()
  time_optimize = tstart - tfinish
  if( myid  ==  0 ) then
     secs = int(time_optimize)
     hrs = secs/(60*60)
     mins = (secs-hrs*60*60)/60
     secs = secs-hrs*60*60-mins*60
     if(hrs>0)then
         write(ounit, '(A, 3(I6, A3))') "focus   : Optimization took ",hrs," H ", mins," M ",secs," S."
     elseif(mins>0)then
         write(ounit, '(A, 2(I6, A3))') "focus   : Optimization took ", mins," M ",secs," S."
     else
         write(ounit, '(A, ES12.5, A3)') "focus   : Optimization took ", time_optimize," S."
     endif
  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  call unpacking(xdof)  ! unpack the optimized xdof array;

  if (myid == 0) write(ounit, *) "-----------POST-PROCESSING-----------------------------------"

  select case( case_postproc )

  case( 0 ) 
  case( 1 ) ; call diagnos ; 
  case( 2 ) ; call diagnos ; call specinp !; call saving 
 !case( 2 ) ; call saving  ; call diagnos ; call wtmgrid  ! write mgrid file;
  case( 3 ) ; call diagnos ; call poinplot ! Poincare plots; for future; 
 ! case( 3 ) ;  call poinplot ! Poincare plots; for future; 
  case( 4 ) ; call diagnos ; call boozmn ; call poinplot ! Last closed surface
  case( 5 ) ; call diagnos ; call wtmgrid  ! write mgrid file
 !case( 4 ) ; call saving  ; call diagnos ; call resonant ! resonant harmonics analysis; for future; 

  end select

  call saving ! save all the outputs
  call MPI_BARRIER( MPI_COMM_FOCUS, ierr )

  tfinish = MPI_Wtime()
  time_postproc = tfinish - tstart
  if( myid  ==  0 )then
     secs = int(time_postproc)
     hrs = secs/(60*60)
     mins = (secs-hrs*60*60)/60
     secs = secs-hrs*60*60-mins*60
     if(hrs>0)then
         write(ounit, '(A, 3(I6, A3))') "focus   : Post-processing took ",hrs," H ", mins," M ",secs," S."
     elseif(mins>0)then
         write(ounit, '(A, 2(I6, A3))') "focus   : Post-processing took ", mins," M ",secs," S."
     else
         write(ounit, '(A, ES12.5, A3)') "focus   : Post-processing took ", time_postproc," S."
     endif
  endif

  if(myid == 0) write(ounit, *) "-------------------------------------------------------------"

  return 

end subroutine focus

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine mute(action)
  use globals, only: ounit
  implicit none 

  INTEGER,intent(in) :: action
  INTEGER, parameter :: iopen = 1, iclose = 0
  INTEGER            :: ios
  CHARACTER(LEN=*), parameter  :: tmp_out = 'tmp.focus_output' 

  ! open a tmp file for screen output
  if (action == iopen) then
    ounit = 37 ! something not used
    open(ounit, file=tmp_out, status="unknown", action="write", iostat=ios) ! create a scratch file?
    if (ios.ne.0) print *, "something wrong with open a tmp file in focuspy.mute. IOSTAT=", ios 
  else
    close(ounit)
    ounit = 6 ! recover to screen output
  endif
  return
end subroutine mute

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine fortran_allocate(array, dim, lbound)
  ! safely allocate 1D array called in python
  use globals, only : dp, ounit
  implicit none

  REAL, ALLOCATABLE, INTENT(INOUT) :: array(:)
  INTEGER, INTENT(IN) :: dim
  INTEGER, OPTIONAL, INTENT(IN) :: lbound
 
  INTEGER :: n0, ierr
  INTEGER :: icount = 0

  if ( allocated(array) ) then
    deallocate(array, stat=ierr)
    if ( ierr.ne.0 ) write(ounit, *) "Problems in deallocating the array, ierr=", ierr, ", icount=", icount
  end if

  n0 = 1
  if (present(lbound)) n0 = lbound
  allocate( array(n0:n0+dim-1), stat=ierr)
  if ( ierr.ne.0 ) write(ounit, *) "Problems in allocating the array, ierr=", ierr, ", icount=", icount, &
    & ", n0=", n0, ", dim=", dim
  return
end subroutine fortran_allocate

subroutine allocate_coils()
  ! allocate the allocatable derived types of coil data
  use globals, only: Ncoils, coil, DoF, FouCoil, ounit
  implicit none

  INTEGER :: ierr

  allocate( FouCoil(1:Ncoils), stat=ierr)
  if ( ierr.ne.0 ) write(ounit, *) "Problems in allocating the array FouCoil, ierr=", ierr
  allocate(    coil(1:Ncoils), stat=ierr)
  if ( ierr.ne.0 ) write(ounit, *) "Problems in allocating the array coil, ierr=", ierr
  allocate(     DoF(1:Ncoils), stat=ierr)
  if ( ierr.ne.0 ) write(ounit, *) "Problems in allocating the array DoF, ierr=", ierr

  return
end subroutine allocate_coils

subroutine get_coils(i, n, xx, yy, zz)
  ! access the data of the i-th coil
  !
  ! import focus
  ! xx = np.zeros(focus.globals.nseg)
  ! yy = np.zeros(focus.globals.nseg)
  ! zz = np.zeros(focus.globals.nseg)
  ! focus.get_coils(1, focus.globals.nseg, xx, yy, zz)
  !
  use globals, only : dp, coil, ounit
  implicit none

  INTEGER,intent(in) :: i, n
  REAL,intent(out) ::  xx(n), yy(n), zz(n)

  if (n .ne. (size(coil(i)%xx)-1)) write(ounit, *) "The size of array doesn't match.", n, size(coil(i)%xx)

  xx = coil(i)%xx(0:n-1)
  yy = coil(i)%yy(0:n-1)
  zz = coil(i)%zz(0:n-1)

  return 
end subroutine get_coils