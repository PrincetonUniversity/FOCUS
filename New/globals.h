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

  REAL, parameter    :: zero       =         0.0
  REAL, parameter    :: one        =         1.0
  REAL, parameter    :: two        =         2.0
  REAL, parameter    :: three      =         3.0
  REAL, parameter    :: four       =         4.0
  REAL, parameter    :: five       =         5.0
  REAL, parameter    :: six        =         6.0
  REAL, parameter    :: seven      =         7.0
  REAL, parameter    :: eight      =         8.0
  REAL, parameter    :: nine       =         9.0
  REAL, parameter    :: ten        =        10.0
  REAL, parameter    :: eleven     =        11.0
  REAL, parameter    :: twelve     =        12.0
  
  REAL, parameter    :: hundred    =       100.0
  REAL, parameter    :: thousand   =      1000.0
  
  REAL, parameter    :: half       =  one / two
  REAL, parameter    :: third      =  one / three 
  REAL, parameter    :: quart      =  one / four
  REAL, parameter    :: fifth      =  one / five
  REAL, parameter    :: sixth      =  one / six
  
  REAL, parameter    :: pi         =  3.141592653589793238462643383279502884197
  REAL, parameter    :: pi2        =  pi * two
  REAL, parameter    :: bsconstant =  1.0E-07   !biot-savart constant
  REAL, parameter    :: antibscont =  1.0E-07 / bsconstant
  REAL, parameter    :: mu0        =  2.0E-07 * pi2
  REAL, parameter    :: goldenmean =  1.618033988749895 ! golden mean = ( one + sqrt(five) ) / two ;    
  
  INTEGER            :: ounit      =  6 ! output to screen
  INTEGER            :: eunit      =  0 ! error unit
  INTEGER            :: wunit      =  7 ! write unit
  INTEGER            :: runit      =  8 ! read  unit
  INTEGER            :: funit      =  9 ! general I/O
  INTEGER            :: bunit      = 10 ! backup unit for I/O
  INTEGER            :: punit      = 11 ! Poincare
  
  CHARACTER(LEN=100) :: ext       ! extention;
  CHARACTER(LEN=100) :: inputfile ! input namelist;
  CHARACTER(LEN=100) :: surffile  ! surface file;
  CHARACTER(LEN=100) :: axisfile  ! axis file;
  CHARACTER(LEN=100) :: coilfile  ! FOCUS coil file ! set in initial; SRH; 29 Sep 17;
  CHARACTER(LEN=100) :: inpcoils  ! input coils.ext file
 !CHARACTER(LEN=100) :: harmfile  ! harmonics file
  CHARACTER(LEN=100) :: hdf5file  ! hdf5 file
  CHARACTER(LEN=100) :: outcoils  ! output ext.coils file
  CHARACTER(LEN=100) :: outplots  ! output ext.coils file
  
  INTEGER            :: Icheck         =        0         ! checking;

  INTEGER            :: Isurface       =        1         
  REAL               :: minorrad       =        0.100D-00
  REAL               :: ellipticity    =        1.000D+00
  INTEGER            :: nrotate        =        0
  REAL               :: zetaoff        =        0.000D+00

  INTEGER            :: Nteta          =       64
  INTEGER            :: Nzeta          =       64

  INTEGER            :: Initialize     =        0
  INTEGER            :: Ncoils         =        0
  REAL               :: init_current   =        1.000D+06
  REAL               :: init_radius    =        1.000D+00
  INTEGER            :: IsVaryCurrent  =        1
  INTEGER            :: IsVaryGeometry =        1
  INTEGER            :: NFcoil         =        2

  INTEGER            :: Nsegments      =       64

  INTEGER            :: case_length    =        1
  REAL               :: weight_bnorm   =        1.000D+00
  REAL               :: weight_bharm   =        0.000D+00
  REAL               :: weight_tflux   =        1.000D+00
  REAL               :: target_tflux   =        1.000D+00
  REAL               :: weight_length  =        1.000D-02
  REAL               :: target_length  =        0.000D+00 ! REDUNDANT; 12 Nov 17;

  REAL               :: wspectral      =        1.000D-06
  REAL               :: pspectral      =        2.000D+00

  REAL               :: converged      =        1.000D-03

  REAL               ::  tauend        =        1.000D-00
  INTEGER            :: Ntauout        =       10
  REAL               ::  tautol        =        1.000D-03

  INTEGER            :: NRKsave        =       10
  INTEGER            :: NRKstep        =       10
  REAL               ::  RKstep        =        1.000D-02

  REAL               :: friction       =        0.500D-00

  INTEGER            :: Iminimize      =        0

  INTEGER            :: Ntrj           =        1
  INTEGER            :: Npts           =        0
  REAL               :: odetol         =        1.0D-06

  namelist / focusin /  Icheck         , &
                        Isurface       , &
                        minorrad       , &
                        ellipticity    , &
                        nrotate        , &
                        zetaoff        , &
                        Nteta          , &
                        Nzeta          , &
                        Initialize     , &
                        Ncoils         , &
                        init_current   , &
                        init_radius    , &
                        IsVaryCurrent  , &
                        IsVaryGeometry , &
                        NFcoil         , &
                        Nsegments      , &
                        case_length    , &
                        weight_bnorm   , &
                        weight_tflux   , &
                        target_tflux   , &
                        weight_length  , &
                        target_length  , & ! REDUNDANT; 12 Nov 17;
                        wspectral      , &
                        pspectral      , &
                        converged      , &
                         tauend        , &
                        Ntauout        , &
                         tautol        , &
                        NRKsave        , &
                        NRKstep        , &
                         RKstep        , &
                        friction       , &
                        Iminimize      , &
                        Ntrj           , &
                        Npts           , &
                        odetol

  INTEGER                 :: myid, ncpu
  REAL                    :: machprec, vsmall, small, sqrtmachprec
  CHARACTER               :: nodelabel*3
  
  INTEGER                 :: Nt, Nz, Ntz, Ns ! shorthand;

  type toroidalsurface
     INTEGER              :: Nteta, Nzeta
     REAL                 :: area, vol ! surface area;
     REAL   , allocatable :: csarea(:) ! cross section area;
     REAL   , allocatable :: xx(:,:), yy(:,:), zz(:,:), nx(:,:), ny(:,:), nz(:,:), ds(:,:), xt(:,:), yt(:,:), zt(:,:)
     REAL   , allocatable :: dL(:), dT(:,:), dB(:,:,:), Bp(:,:) ! total normal magnetic field; plasma normal magnetic field;
     REAL   , allocatable :: Bs(:,:), Bt(:,:), Bz(:,:) ! components of magnetic field at surface; 12 Nov 17;
     REAL   , allocatable :: EE(:,:), FF(:,:), GG(:,:) ! coefficients of 1st fundamental form;
     REAL   , allocatable :: LL(:,:), MM(:,:), PP(:,:) ! coefficients of 2nd fundamental form;
     REAL   , allocatable :: HH(:,:)                   ! mean curvature;
  end type toroidalsurface
  type(toroidalsurface)   :: surf

  type spacecurve
     INTEGER              :: itype, NF, Ifree, Lfree, NS
     REAL                 ::     I           , Le, maxcurv
     REAL   , allocatable :: xc(:), xs(:), yc(:), ys(:), zc(:), zs(:)
     REAL   , allocatable :: xx(:), yy(:), zz(:), xt(:), yt(:), zt(:), xa(:), ya(:), za(:)
     REAL   , allocatable :: dL(:), dT(:,:), dB(:,:,:), RR(:,:,:,:)
     INTEGER              :: gdof, idof
     character(LEN=10)    :: name
  end type spacecurve
  type(spacecurve), allocatable :: coil(:)  

  INTEGER                 :: Nfp
  INTEGER                 :: Nfou = 0, NBnf = 0 ! this should be local to readsrf; 04 Sep 17;
  INTEGER, allocatable    :: bim(:), bin(:), Bnim(:), Bnin(:)
  REAL   , allocatable    :: Rbc(:), Zbs(:), Rbs(:), Zbc(:), Bnc(:), Bns(:), cosip(:), sinip(:)
    
  REAL                    :: Tfluxerr
  REAL   , allocatable    :: totlengt(:), Tfluxave(:), Bdotnsqd(:)

  INTEGER                 :: axisNF
  REAL   , allocatable    :: axisxc(:), axisxs(:), axisyc(:), axisys(:), axiszc(:), axiszs(:), axistc(:), axists(:)
  
  REAL,    allocatable    :: coilsI(:) !coilsI stores the true currents, coil%I stores scaled current;?
  INTEGER, allocatable    :: coilseg(:)
  character(LEN=20), allocatable :: coilname(:)

  REAL                    :: tstart, ffbest, fforig

  LOGICAL                 :: Ldescent

  REAL                    :: discretecurve, deltatheta, discretesurface
  REAL, allocatable       :: cmt(:,:), smt(:,:)

  INTEGER                 :: iarchive
  REAL, allocatable       :: iota(:,:)

  REAL, allocatable       :: tdof(:) ! tangent direction in independent degree-of-freedom (Fourier) space;
  REAL, allocatable       :: Mdof(:) ! spectral gradient in independent degree-of-freedom (Fourier) space;

  REAL                    :: xyaxis(1:2)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

end module globals

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
