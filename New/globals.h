!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

!title (input) ! Description of global variables.

!latex \briefly{Defines input namelists and global variables, details of input namelist can be viwed at 
!latex \link{initial}.}

!latex \calledby{\link{focus}}
!latex \calls{\link{initial}}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

module globals
  
  implicit none
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

  REAL, parameter      :: zero       =         0.0
  REAL, parameter      :: one        =         1.0
  REAL, parameter      :: two        =         2.0
  REAL, parameter      :: three      =         3.0
  REAL, parameter      :: four       =         4.0
  REAL, parameter      :: five       =         5.0
  REAL, parameter      :: six        =         6.0
  REAL, parameter      :: seven      =         7.0
  REAL, parameter      :: eight      =         8.0
  REAL, parameter      :: nine       =         9.0
  REAL, parameter      :: ten        =        10.0
  REAL, parameter      :: eleven     =        11.0
  REAL, parameter      :: twelve     =        12.0
  
  REAL, parameter      :: hundred    =       100.0
  REAL, parameter      :: thousand   =      1000.0
  
  REAL, parameter      :: half       =  one / two
  REAL, parameter      :: third      =  one / three 
  REAL, parameter      :: quart      =  one / four
  REAL, parameter      :: fifth      =  one / five
  REAL, parameter      :: sixth      =  one / six
  
  REAL, parameter      :: pi         =  3.141592653589793238462643383279502884197
  REAL, parameter      :: pi2        =  pi * two
  REAL, parameter      :: bsconstant =  1.0E-07   !biot-savart constant
  REAL, parameter      :: antibscont =  1.0E-07 / bsconstant
  REAL, parameter      :: mu0        =  2.0E-07 * pi2
  REAL, parameter      :: goldenmean =  1.618033988749895 ! golden mean = ( one + sqrt(five) ) / two ;    
  
  INTEGER              :: ounit      =  6 ! output to screen
  INTEGER              :: eunit      =  0 ! error unit
  INTEGER              :: wunit      =  7 ! write unit
  INTEGER              :: runit      =  8 ! read  unit
  INTEGER              :: funit      =  9 ! general I/O
  INTEGER              :: bunit      = 10 ! backup unit for I/O
  
  CHARACTER(LEN=100)   :: ext       ! extention;
  CHARACTER(LEN=100)   :: inputfile ! input namelist;
  CHARACTER(LEN=100)   :: surffile  ! surface file;
  CHARACTER(LEN=100)   :: axisfile  ! axis file;
  CHARACTER(LEN=100)   :: coilfile  ! FOCUS coil file ! set in initial; SRH; 29 Sep 17;
  CHARACTER(LEN=100)   :: inpcoils  ! input coils.ext file
 !CHARACTER(LEN=100)   :: harmfile  ! harmonics file
  CHARACTER(LEN=100)   :: hdf5file  ! hdf5 file
  CHARACTER(LEN=100)   :: outcoils  ! output ext.coils file
  
  INTEGER              :: Icheck         =        0         ! checking;

  INTEGER              :: Isurface       =        1         
  INTEGER              :: case_surface   =        0         ! redundant;
  REAL                 :: minorrad       =        0.100D-00
  REAL                 :: knotsurf       =        0.200D-00 ! redundant;
  REAL                 :: ellipticity    =        0.000D+00

  INTEGER              :: Nteta          =       64
  INTEGER              :: Nzeta          =       64

  INTEGER              :: Initialize     =        0
  INTEGER              :: case_init      =        0 ! redundant;
  INTEGER              :: Ncoils         =        0
  REAL                 :: init_current   =        1.000D+06
  REAL                 :: init_radius    =        1.000D+00
  INTEGER              :: IsVaryCurrent  =        1
  INTEGER              :: IsVaryGeometry =        1
  INTEGER              :: NFcoil         =        2

  INTEGER              :: Nsegments      =       64
  INTEGER              :: Nseg           =      128 ! redundant;

  INTEGER              :: case_length    =        1
  REAL                 :: weight_bnorm   =        1.000D+00
  REAL                 :: weight_bharm   =        0.000D+00
  REAL                 :: weight_tflux   =        0.000D+00
  REAL                 :: target_tflux   =        0.000D+00
  REAL                 :: weight_ttlen   =        0.000D+00
  REAL                 :: target_length  =        0.000D+00

  REAL                 :: converged      =        1.000D-03

  REAL                 ::  tauend        =        1.000D-00
  INTEGER              :: Ntauout        =       10
  REAL                 ::  tautol        =        1.000D-03

  INTEGER              :: DF_maxiter     =       10         ! redundant;
  REAL                 :: DF_xtol        =        1.000D-08 ! redundant;
  REAL                 :: DF_tauend      =        1.000D+00 ! redundant;

  REAL                 :: friction       =        0.500D-00

  INTEGER              :: Iminimize      =        0

  namelist / focusin /  Icheck         , &
                        Isurface       , &
                        case_surface   , & ! redundant;
                        minorrad       , &
                        knotsurf       , & ! redundant;
                        ellipticity    , &
                        Nteta          , &
                        Nzeta          , &
                        Initialize     , &
                        case_init      , & ! redundant;
                        Ncoils         , &
                        init_current   , &
                        init_radius    , &
                        IsVaryCurrent  , &
                        IsVaryGeometry , &
                        NFcoil         , &
                        Nsegments      , &
                        Nseg           , & ! redundant;
                        case_length    , &
                        weight_bnorm   , &
                        weight_tflux   , &
                        target_tflux   , &
                        weight_ttlen   , &
                        target_length  , &
                        converged      , &
                         tauend        , &
                        Ntauout        , &
                         tautol        , &
                        DF_xtol        , & ! redundant;
                        DF_maxiter     , & ! redundant;
                        DF_tauend      , & ! redundant;
                        friction       , &
                        Iminimize

  INTEGER              :: myid, ncpu
  REAL                 :: machprec, vsmall, small, sqrtmachprec
  CHARACTER            :: nodelabel*3

  INTEGER              :: Nt, Nz, Ntz, Ns ! shorthand;

  type toroidalsurface
     INTEGER              :: Nteta, Nzeta
     REAL                 :: area ! surface area;
     REAL   , allocatable :: csarea(:) ! cross section area;
     REAL   , allocatable :: xx(:,:), yy(:,:), zz(:,:), nx(:,:), ny(:,:), nz(:,:), ds(:,:), xt(:,:), yt(:,:), zt(:,:)
     REAL   , allocatable :: dL(:), dT(:,:), dB(:,:,:), Bp(:,:) ! total normal magnetic field; plasma normal magnetic field;
  end type toroidalsurface

  type spacecurve
     INTEGER              :: itype, NF, Ifree, Lfree, NS
     REAL                 ::     I    , Lo   , Le, maxcurv
     REAL   , allocatable :: xc(:), xs(:), yc(:), ys(:), zc(:), zs(:)
     REAL   , allocatable :: xx(:), yy(:), zz(:), xt(:), yt(:), zt(:), xa(:), ya(:), za(:)
     REAL   , allocatable :: dL(:), dT(:,:), dB(:,:,:)
     INTEGER              :: gdof, idof
     character(LEN=10)    :: name
  end type spacecurve

  type(spacecurve)     , allocatable :: coil(:)  
  type(toroidalsurface)              :: surf

  INTEGER              :: Nfp
  INTEGER              :: Nfou = 0, NBnf = 0 ! this should be local to readsrf; 04 Sep 17;
  INTEGER, allocatable :: bim(:), bin(:), Bnim(:), Bnin(:)
  REAL   , allocatable :: Rbc(:), Zbs(:), Rbs(:), Zbc(:), Bnc(:), Bns(:), cosip(:), sinip(:)
    
  REAL                 :: Inorm = one, Gnorm = one                !current and geometry normalizations;

  REAL                 :: Tfluxerr
  REAL   , allocatable :: totlengt(:), Tfluxave(:), Bdotnsqd(:)

  INTEGER              :: axisNF
  REAL   , allocatable :: axisxc(:), axisxs(:), axisyc(:), axisys(:), axiszc(:), axiszs(:)
  
  REAL,    allocatable :: coilsI(:) !coilsI stores the true currents, coil%I stores scaled current;?
  INTEGER, allocatable :: coilseg(:)
  character(LEN=20), allocatable :: coilname(:)

  REAL                 :: tstart, ffbest, fforig

  LOGICAL              :: Ldescent

  REAL                 :: discretecurve, deltatheta, discretesurface
  REAL, allocatable    :: cmt(:,:), smt(:,:)

  INTEGER              :: iarchive

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

end module globals

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
