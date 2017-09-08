!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!!-!-!-!-!-!-!-!-!-!-!-!-

!title (main) ! Main program.

!latex \briefly{Main program.}

!latex \calledby{\link{}}
!latex \calls{\link{initial}, \link{rdsurf}, \link{rdcoils}, 
!latex \link{saving}, \link{solvers}}.

!latex \tableofcontents

!latex  \subsection{Brief introduction}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

PROGRAM focus
  
  use globals, only : zero, ncpu, myid, ounit, ierr, astat, ext, &
                      Nseg, &
                      Ndof, xdof, &
                      tstart, &
                      Ldescent, &
             !         case_surface, case_optimize, &
             !         case_postproc, xdof, tstart, tfinish, time_initialize, time_optimize, time_postproc, surf, pi2, &
             !         Nfp, knotsurf, half, Nteta, Nzeta, zero, IsSymmetric, discretefactor, &
                      inputfile, surffile, coilfile, harmfile, hdf5file, inpcoils, outcoils
             !         Ncoils, coil
                      
  
  use mpi
  
  implicit none
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  INTEGER :: icoil
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  myid = 0 ; ncpu = 1
  
  call MPI_INIT( ierr )
  call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
  call MPI_COMM_SIZE( MPI_COMM_WORLD, ncpu, ierr )
  
  tstart =  MPI_WTIME()
  if( myid == 0 ) write(ounit,'("focus   : " 10x " : start ; ncpu =",i5," ;")') ncpu
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( myid == 0 ) call getarg(1,ext) ! read command line input;
  
  ClBCAST( ext, 100, 0 )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  inputfile = trim(ext)//".input"
  surffile  = "plasma.boundary"
  coilfile  = trim(ext)//".focus"
  harmfile  = trim(ext)//".harmonics"
  hdf5file  = "focus_"//trim(ext)//".h5"
  inpcoils  = "coils."//trim(ext)
  outcoils  = trim(ext)//".coils"
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  call initial ! reads input;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  call readsrf ! reads reference surface;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  call rdcoils ! reads/initializes coil geometries etc. ;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!  if( myid.eq.0 ) then
!   
!   do icoil = 1, Ncoils
!    write(ounit,1000), coil(icoil)%name, coil(icoil)%NF, coil(icoil)%NS, coil(icoil)%I, coil(icoil)%Ic
!    write(ounit,'("focus   : " 10x " : x_c^"i3.3" ="999es13.5)') icoil, coil(icoil)%xc(0:coil(icoil)%NF)
!    write(ounit,'("focus   : " 10x " : x_s^"i3.3" ="999es13.5)') icoil, coil(icoil)%xs(0:coil(icoil)%NF)
!    write(ounit,'("focus   : " 10x " : y_c^"i3.3" ="999es13.5)') icoil, coil(icoil)%yc(0:coil(icoil)%NF)
!    write(ounit,'("focus   : " 10x " : y_s^"i3.3" ="999es13.5)') icoil, coil(icoil)%ys(0:coil(icoil)%NF)
!    write(ounit,'("focus   : " 10x " : z_c^"i3.3" ="999es13.5)') icoil, coil(icoil)%zc(0:coil(icoil)%NF)
!    write(ounit,'("focus   : " 10x " : z_s^"i3.3" ="999es13.5)') icoil, coil(icoil)%zs(0:coil(icoil)%NF)
!   enddo
!   
!1000 format("focus   : " 10x " : name = "a7" ; NF ="i3" ; NS ="i4" ; I ="es13.5" ; Ic ="i2" ;")
!   
!  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  call setflag ! this will count/pack degrees-of-freedom, initialize weights, etc.  . . . 04 Sep 17;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( Ldescent ) call descent( Ndof, xdof(1:Ndof) )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! write(ounit,'("focus   : " 10x " : calling drivers ;")') ; pause

! call drivers ! this will call numerical optimization algorithms; 04 Sep 17;

! write(ounit,'("focus   : " 10x " : called  drivers ;")') ; pause

! if ( case_optimize /= 0 ) call solvers ! call different solvers;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! call MPI_BARRIER( MPI_COMM_WORLD, ierr )

! tstart = MPI_Wtime()
! time_optimize = tstart - tfinish
! if( myid  ==  0 ) then
!    secs = int(time_optimize)
!    hrs = secs/(60*60)
!    mins = (secs-hrs*60*60)/60
!    secs = secs-hrs*60*60-mins*60
!    if(hrs>0)then
!        write(ounit, '(A, 3(I6, A3))') "focus   : Optimization took ",hrs," H ", mins," M ",secs," S."
!    elseif(mins>0)then
!        write(ounit, '(A, 2(I6, A3))') "focus   : Optimization took ", mins," M ",secs," S."
!    else
!        write(ounit, '(A, ES12.5, A3)') "focus   : Optimization took ", time_optimize," S."
!    endif
! endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  call archive ! this will save information to output/restart files; 04 Sep 17;

! call unpacking( xdof )  ! unpack the optimized xdof array;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! if (myid == 0) write(ounit, *) "-----------POST-PROCESSING-----------------------------------"

! select case( case_postproc )
! case( 0 ) ;              ; call saving
! case( 1 ) ; call diagnos ; call saving 
! case( 2 ) ; call saving  ; call diagnos ; call wtmgrid  ! write mgrid file;
! case( 3 ) ; call saving  ; call diagnos ; call poinplot ! Poincare plots; for future; 
! case( 4 ) ; call saving  ; call diagnos ; call resonant ! resonant harmonics analysis; for future; 
! case default
!  FATAL( focus, .true., selected case_postproc not supported )
! end select

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! call MPI_BARRIER( MPI_COMM_WORLD, ierr )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! tfinish = MPI_Wtime()
! time_postproc = tfinish - tstart
! if( myid  ==  0 )then
!     secs = int(time_postproc)
!     hrs = secs/(60*60)
!     mins = (secs-hrs*60*60)/60
!     secs = secs-hrs*60*60-mins*60
!     if(hrs>0)then
!         write(ounit, '(A, 3(I6, A3))') "focus   : Post-processing took ",hrs," H ", mins," M ",secs," S."
!     elseif(mins>0)then
!         write(ounit, '(A, 2(I6, A3))') "focus   : Post-processing took ", mins," M ",secs," S."
!     else
!         write(ounit, '(A, ES12.5, A3)') "focus   : Post-processing took ", time_postproc," S."
!     endif
!  endif
!
!  if(myid == 0) write(ounit, *) "-------------------------------------------------------------"

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  call MPI_FINALIZE( ierr )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  !call cleanup

  stop
  
END PROGRAM focus

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
