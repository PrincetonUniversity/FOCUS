
!title (input) ! Description of global variables.

!latex \briefly{Defines input namelists and global variables, details of input namelist can be viwed at 
!latex \link{initial}.}

!latex \calledby{\link{focus}}
!latex \calls{\link{initial}}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

module globals
  
  implicit none
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  CHARACTER(LEN=10), parameter :: version='v0.2.05' ! version number

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsection{Useful parameters}
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
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsection{IO units}
  INTEGER              :: ounit      =  6 ! output to screen
  INTEGER              :: eunit      =  0 ! error unit
  INTEGER              :: wunit      =  7 ! write unit
  INTEGER              :: runit      =  8 ! read  unit
  INTEGER              :: funit      =  9 ! general I/O
  INTEGER              :: bunit      = 10 ! backup unit for I/O
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  CHARACTER(LEN=100)   :: ext       ! extention
  CHARACTER(LEN=100)   :: inputfile ! input namelist
  CHARACTER(LEN=100)   :: surffile  ! surface file
  CHARACTER(LEN=100)   :: knotfile  ! knototran file
  CHARACTER(LEN=100)   :: coilfile  ! FOCUS coil file
  CHARACTER(LEN=100)   :: harmfile  ! harmonics file
  CHARACTER(LEN=100)   :: hdf5file  ! hdf5 file
  CHARACTER(LEN=100)   :: inpcoils  ! input coils.ext file
  CHARACTER(LEN=100)   :: outcoils  ! output ext.coils file
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!latex \subsection{Input namelist: \type{focusin}}
  INTEGER              :: IsQuiet        =       -1        
  INTEGER              :: IsSymmetric    =        0 
        
  INTEGER              :: case_surface   =        0         
  REAL                 :: knotsurf       =        0.200D-00
  REAL                 :: ellipticity    =        0.000D+00
  INTEGER              :: Nteta          =       64           
  INTEGER              :: Nzeta          =       64  

  INTEGER              :: case_init      =        0
  INTEGER              :: case_coils     =        1
  INTEGER              :: Ncoils         =        0
  REAL                 :: init_current   =        1.000D+06        
  REAL                 :: init_radius    =        1.000D+00
  INTEGER              :: IsVaryCurrent  =        1         
  INTEGER              :: IsVaryGeometry =        1         
  INTEGER              :: NFcoil         =        4         
  INTEGER              :: Nseg           =      128 
              
  INTEGER              :: IsNormalize    =        1
  INTEGER              :: IsNormWeight   =        1
  INTEGER              :: case_bnormal   =        0
  INTEGER              :: case_length    =        1         
  REAL                 :: weight_bnorm   =        1.000D+00
  REAL                 :: weight_bharm   =        0.000D+00
  REAL                 :: weight_tflux   =        0.000D+00
  REAL                 :: target_tflux   =        0.000D+00
  REAL                 :: weight_ttlen   =        0.000D+00
  REAL                 :: target_length  =        0.000D+00
  REAL                 :: weight_cssep   =        0.000D+00
  REAL                 :: cssep_factor   =        1.000D+00 
  REAL                 :: weight_specw   =        0.000D+00
  REAL                 :: weight_ccsep   =        0.000D+00
  REAL                 :: weight_inorm   =        1.000D+00
  REAL                 :: weight_gnorm   =        1.000D+00

  INTEGER              :: case_optimize  =        1
  REAL                 :: exit_tol       =        1.000D-04
  INTEGER              :: DF_maxiter     =        0
  REAL                 :: DF_xtol        =        1.000D-08     
  REAL                 :: DF_tausta      =        0.000D+00
  REAL                 :: DF_tauend      =        1.000D+00               
 
  INTEGER              :: CG_maxiter     =        0
  REAL                 :: CG_xtol        =        1.000D-08
  REAL                 :: CG_wolfe_c1    =        1.000D-04
  REAL                 :: CG_wolfe_c2    =        0.1

  INTEGER              :: HN_maxiter     =        0
  REAL                 :: HN_xtol        =        1.000D-08
  REAL                 :: HN_factor      =      100.0

  INTEGER              :: TN_maxiter     =        0
  REAL                 :: TN_xtol        =        1.000D-08
  INTEGER              :: TN_reorder     =        0
  REAL                 :: TN_cr          =        0.1

  INTEGER              :: case_postproc  =        1
  INTEGER              :: save_freq      =        1
  INTEGER              :: save_coils     =        0 
  INTEGER              :: save_harmonics =        0
  INTEGER              :: save_filaments =        0              
                                                         

  

  namelist / focusin /  IsQuiet        , &
                        IsSymmetric    , & 
                        case_surface   , &
                        knotsurf       , &
                        ellipticity    , & 
                        Nteta          , &
                        Nzeta          , & 
                        case_init      , &
                        case_coils     , &  
                        Ncoils         , &
                        init_current   , & 
                        init_radius    , & 
                        IsVaryCurrent  , & 
                        IsVaryGeometry , & 
                        NFcoil         , &
                        Nseg           , &
                        IsNormalize    , &
                        IsNormWeight   , &
                        case_bnormal   , &
                        case_length    , &
                        weight_bnorm   , &
                        weight_bharm   , &
                        weight_tflux   , &
                        target_tflux   , &
                        weight_ttlen   , &
                        target_length  , &
                        weight_cssep   , &
                        cssep_factor   , &
                        weight_specw   , &
                        weight_ccsep   , &
                        weight_inorm   , &
                        weight_gnorm   , &
                        case_optimize  , &
                        exit_tol       , &
                        DF_maxiter     , & 
                        DF_xtol        , & 
                        DF_tausta      , &  
                        DF_tauend      , &       
                        CG_maxiter     , & 
                        CG_xtol        , & 
                        CG_wolfe_c1    , &
                        CG_wolfe_c2    , &
                        HN_maxiter     , &  
                        HN_xtol        , & 
                        HN_factor      , & 
                        TN_maxiter     , &
                        TN_reorder     , &
                        TN_xtol        , &
                        TN_cr          , &  
                        case_postproc  , & 
                        save_freq      , & 
                        save_coils     , &
                        save_harmonics , &
                        save_filaments
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!latex  \subsection{MPI stuffs}
  INTEGER              :: myid, ncpu
  REAL                 :: machprec, vsmall, small, sqrtmachprec
  CHARACTER            :: nodelabel*3

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsection{surface and coils data}
  type toroidalsurface
     INTEGER              :: Nteta, Nzeta
     REAL   , allocatable :: xx(:,:), yy(:,:), zz(:,:), nx(:,:), ny(:,:), nz(:,:), &
                             ds(:,:), xt(:,:), yt(:,:), zt(:,:), bn(:,:), pb(:,:), &
                             Bx(:,:), By(:,:), Bz(:,:)
  end type toroidalsurface

  type arbitrarycoil
     INTEGER              :: NS, Ic, Lc, itype
     REAL                 :: I,  L, Lo, maxcurv
     REAL   , allocatable :: xx(:), yy(:), zz(:), xt(:), yt(:), zt(:), xa(:), ya(:), za(:), &
                             dl(:), dd(:), &
                             Ax(:,:), Ay(:,:), Az(:,:)
     character(LEN=10)    :: name
  end type arbitrarycoil

  type FourierCoil
     INTEGER              :: NF
     REAL   , allocatable :: xc(:), xs(:), yc(:), ys(:), zc(:), zs(:)
  end type FourierCoil

  type DegreeOfFreedom
     INTEGER              :: ND
     REAL   , allocatable :: xdof(:), xof(:,:), yof(:,:), zof(:,:)
  end type DegreeOfFreedom
  
  type(arbitrarycoil)  , allocatable :: coil(:)  
  type(toroidalsurface), allocatable :: surf(:)
  type(FourierCoil)    , allocatable :: FouCoil(:)
  type(DegreeOfFreedom), allocatable :: DoF(:)

  INTEGER              :: Nfou=0, Nfp=0, NBnf=0, Npc = 1
  INTEGER, allocatable :: bim(:), bin(:), Bnim(:), Bnin(:)
  REAL   , allocatable :: Rbc(:), Zbs(:), Rbs(:), Zbc(:), Bnc(:), Bns(:), cosip(:), sinip(:)
    
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsection{Packing and unpacking}
  INTEGER              :: Cdof, Ndof, nfixcur, nfixgeo, Tdof
  REAL                 :: Inorm = one, Gnorm = one                !current and geometry normalizations;
  REAL   , allocatable :: xdof(:), dofnorm(:)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsection{Optimization}
  ! General target functions;
  INTEGER              :: iout, Nouts
  REAL                 :: chi, discretefactor
  REAL   , allocatable :: t1E(:), t2E(:,:), evolution(:,:), coilspace(:,:), deriv(:,:)
  LOGICAL              :: exit_signal = .False.
  ! Bn surface integration;
  REAL                 :: bnorm
  REAL   , allocatable :: t1B(:), t2B(:,:), bn(:,:), dB(:,:,:)
  ! Bn reasonant harmoics;
  INTEGER              :: NBmn
  INTEGER, allocatable :: Bmnin(:), Bmnim(:)
  REAL                 :: bharm
  REAL   , allocatable :: t1H(:), t2H(:,:), Bmnc(:),Bmns(:), wBmn(:), tBmnc(:), tBmns(:), &
                          carg(:,:), sarg(:,:), iBmnc(:), iBmns(:)
  ! Tflux error;
  REAL                 :: tflux, psi_avg
  REAL   , allocatable :: t1F(:), t2F(:,:)
  ! Length constraint
  REAL                 :: ttlen
  REAL   , allocatable :: t1L(:), t2L(:,:)
  ! Coil-surface spearation
  INTEGER              :: psurf = 1 ! the prevent surface label; default 1 is the plasma boundary
  REAL                 :: cssep
  REAL   , allocatable :: t1S(:), t2S(:,:)
  ! Coil-coil spearation
  REAL                 :: ccsep
  REAL   , allocatable :: t1C(:), t2C(:,:)
  ! Spectral condensation;
  REAL                 :: specw
  REAL   , allocatable :: t1P(:), t2P(:,:)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsection{Knotran parameters}
  REAL                 :: knotphase
  REAL   , allocatable :: xkc(:), xks(:), ykc(:), yks(:), zkc(:), zks(:)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsection{Reading coils file}
  REAL,    allocatable :: coilsX(:,:), coilsY(:,:), coilsZ(:,:), coilsI(:) 
                        !coilsI stores the true currents, coil%I stores scaled current;?
  INTEGER, allocatable :: coilseg(:)
  character(LEN=20), allocatable :: coilname(:)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsection{Time counting}
  REAL                 :: tstart, tfinish, time_initialize, time_optimize, time_postproc

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsection{Miscellaneous}
  REAL                 :: tmpw_bnorm, tmpw_tflux ,tmpt_tflux, tmpw_ttlen, tmpw_specw, tmpw_ccsep, tmpw_bharm
                          !tmp weight for saving to restart file
  REAL, allocatable    :: mincc(:,:)!
  INTEGER              :: ierr, astat

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end module globals