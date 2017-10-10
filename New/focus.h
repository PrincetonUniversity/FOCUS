!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!!-!-!-!-!-!-!-!-!-!-!-!-

!title (main) ! Main program.

!latex \briefly{Main program.}

!latex \calledby{\link{}}
!latex \calls{\link{initial}, \link{rdsurf}, \link{rdcoils}, 
!latex \link{saving}, \link{solvers}}.

!latex \tableofcontents

!latex  \subsection{Brief introduction}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! to set keyboard shortcut in emacs                                        

! (1) define macro         , e.g. \C-x \C-( . . . \C-x \C-)              
! (2) name macro           , e.g. Esc-x name-last-kbd-macro arbitraryname ;
! (3) set keyboard shortcut, e.g. Esc-x global-set-key F12 arbitraryname 

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

program focus
  
  use globals, only : zero, one, two, pi2, sqrtmachprec, ncpu, myid, ounit, tstart, &
                      ext, inputfile, surffile, axisfile, coilfile, inpcoils, hdf5file, outcoils, outplots, &
                      Nt, Nz, surf, Ncoils, Ns, coil, &
                      cmt, smt, discretecurve, &
                      tdof, Mdof, &
                      totlengt, Tfluxave, Bdotnsqd, &
                      target_length, target_tflux, &
                      Ldescent, Iminimize, &
                      fforig, ffbest, iarchive, &
                      weight_length, weight_tflux, weight_bnorm, wspectral, pspectral, &
                      Ntrj, Npts

  use mpi
  
  implicit none
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  INTEGER                :: icoil, mm, maxNF, Ndof, ii, ierr, astat
  REAL                   :: tnow, told, ff, ferr
  REAL     , allocatable :: xdof(:), fdof(:)
  CHARACTER              :: packorunpack

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  myid = 0 ; ncpu = 1
  
  call MPI_INIT( ierr )
  call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
  call MPI_COMM_SIZE( MPI_COMM_WORLD, ncpu, ierr )
  
  tstart = MPI_WTIME()

  if( myid.eq.0 ) write(ounit,'("focus   :   cpu time : start ; ncpu =",i5," ;")') ncpu
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( myid.eq.0 ) call getarg(1,ext) ! read command line input;
  
  ClBCAST( ext, 100, 0 )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  inputfile =              trim(ext)//".fo.in"    ! namelist input;
  surffile  =              trim(ext)//".fo.bdy"   ! Fourier harmonics of reference boundary;
  axisfile  =              trim(ext)//".fo.axis"  ! Fourier harmonics of magnetic axis;
  coilfile  =              trim(ext)//".fo.coils" ! Fourier harmonics of coils (FOCUS format);
  inpcoils  =    "coils."//trim(ext)              ! xgrid-style coils (input);
 !harmfile  =              trim(ext)//".fo.Bmn"   ! Fourier harmonics of target normal field;
  hdf5file  =              trim(ext)//".fo.h5"    ! output file;
  outcoils  = "coils.fo."//trim(ext)              ! xgrid-style coils (input);
  outplots  =         "."//trim(ext)//".fo.dat"
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  call initial ! reads input;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  call readsrf ! reads reference surface;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
8000 continue
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

  call rdcoils ! reads/initializes coil geometries etc. ;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  Ndof = 0 ! total number of degrees-of-freedom in coil representation;
  
  do icoil = 1, Ncoils
   
   if( coil(icoil)%Lfree.ne.0 ) Ndof = Ndof + 6 * coil(icoil)%NF + 3 ! allow geometry to vary;
   if( coil(icoil)%Ifree.ne.0 ) Ndof = Ndof + 1                      ! allow current  to vary;
   
  enddo
  
  tnow = MPI_WTIME()
  
  if( myid.eq.0 ) write(ounit,'("focus   : ",f10.1," : Ndof =",i8," ;")') tnow-tstart, Ndof
  
  if( Ndof.le.0 ) goto 9999

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  SALLOCATE( xdof, (1:Ndof), zero ) ! position = geometry of coils             ;
  SALLOCATE( fdof, (1:Ndof), zero ) ! force    = gradient of objective function;
  
  SALLOCATE( tdof, (1:Ndof), zero ) ! tangent direction in Fourier space       ;
  SALLOCATE( Mdof, (1:Ndof), zero ) ! spectral gradient in Fourier space       ;

  packorunpack = 'P' ; call packdof( Ndof, xdof(1:Ndof), packorunpack )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
   SALLOCATE(        surf%dL, (              0:Ndof), zero )
   SALLOCATE(        surf%dT, (       0:Nz-1,0:Ndof), zero )
   SALLOCATE(        surf%dB, (0:Nt-1,0:Nz-1,0:Ndof), zero )
  
  do icoil = 1, Ncoils
   
   SALLOCATE( coil(icoil)%dL, (              0:Ndof), zero )
   SALLOCATE( coil(icoil)%dT, (       0:Nz-1,0:Ndof), zero )
   SALLOCATE( coil(icoil)%dB, (0:Nt-1,0:Nz-1,0:Ndof), zero )
   
  enddo
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  SALLOCATE( totlengt, (0:Ndof), zero )
  SALLOCATE( Tfluxave, (0:Ndof), zero )
  SALLOCATE( Bdotnsqd, (0:Ndof), zero )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  if( myid.eq.0 ) write(ounit,1010) tnow-tstart, "targets   ", target_length, target_tflux
  
  told = MPI_WTIME()

  call dforce( Ndof, xdof(1:Ndof), ff, fdof(1:Ndof) )
  
  fforig = ff ; ffbest = ff
  
  ferr = sqrt( sum(fdof(1:Ndof)*fdof(1:Ndof)) / Ndof )
  
  tnow = MPI_WTIME()

  if( myid.eq.0 ) write(ounit,1010) tnow-tstart, "difference", totlengt(0)-target_length, Tfluxave(0)-target_tflux, Bdotnsqd(0), ff, ferr, tnow-told
  
1010 format("focus   : ",f10.1," : ",a10," : L =",es18.10," ; F =",es18.10," ; ":"B =",es17.10," ; O =",es17.10," ; |dO| =",es12.05," ; time =",f9.2,"s ;")  

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  call setflag( Ndof, xdof(1:Ndof) ) ! initialize weights, etc. ;
  
  tnow = MPI_WTIME()
  
  if( myid.eq.0 ) write(ounit,1000) tnow-tstart, weight_length, weight_tflux, weight_bnorm, wspectral, pspectral
1000 format("focus   : ",f10.1," : weight_length =",es12.5" ; weight_tflux =",es12.5," ; weight_bnorm =",es12.5,&
                               " ; wspectral =",es12.5," ; pspectral =",f5.2," ;")
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  call archive( Ndof, xdof(1:Ndof), ferr )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( Ldescent ) then
   
   call descent( Ndof, xdof(1:Ndof), ferr )
   
  !deallocate( coil ) 
  !DALLOCATE( xdof )
  !DALLOCATE( fdof )
  !DALLOCATE( surf%dL )
  !DALLOCATE( surf%dT )
  !DALLOCATE( surf%dB )
  !DALLOCATE( totlengt )
  !DALLOCATE( Tfluxave )
  !DALLOCATE( Bdotnsqd )
  !friction = friction + two
  !goto 8000
   
  endif ! end of if( Ldescent ) ;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  select case( Iminimize )
   
  case( 1 ) 
   
   call minimum( Ndof, xdof(1:Ndof) )
   
  end select
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( Npts.ge.1 ) then

   call pcplot

   call archive( Ndof, xdof(1:Ndof), ferr )

  endif
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  call MPI_BARRIER( MPI_COMM_WORLD, ierr )
  
  tnow = MPI_WTIME()
  
  if( myid.eq.0 ) write(ounit,'("focus   : ",f10.1," : finished ; time =",f9.1," min =",f9.2," hours ;")') ( tnow - tstart ) / (/ 1, 60, 3600 /)
  
  call MPI_FINALIZE( ierr )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

9999 continue
  
  stop
  
end program focus

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
