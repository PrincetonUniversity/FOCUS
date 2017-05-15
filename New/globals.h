
!title (input) ! Description of global variables.

!latex \briefly{Defines input namelists and global variables, details of input namelist can be viwed at 
!latex \link{initial}.}

!l tex \calledby{\link{focus}}
!l tex \calls{\link{initialize}}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

module globals
  
  use oculus, only : biotsavart, poincaredata
  
  implicit none
  
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
  
  CHARACTER(LEN=100)   :: ext
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!latex \subsection{Input namelist: \type{focusin}}
  INTEGER              :: IsQuiet        =       -1        
  INTEGER              :: IsSymmetric    =        0 
        
  INTEGER              :: case_surface   =        0         
  REAL                 :: knotsurf       =        0.200D-00
  INTEGER              :: Nteta          =       64           
  INTEGER              :: Nzeta          =       64  

  INTEGER              :: case_init      =        0
  INTEGER              :: case_coils     =        1
  INTEGER              :: Ncoils         =        0
  REAL                 :: init_current   =        1.000D+07        
  REAL                 :: init_radius    =        0.500D+00
  INTEGER              :: IsVaryCurrent  =        1         
  INTEGER              :: IsVaryGeometry =        1         
  REAL                 :: target_length  =        1.000D+00 
  INTEGER              :: NFcoil         =        4         
  INTEGER              :: Nseg           =      128 
        
  INTEGER              :: case_optimizer =        0           
  INTEGER              :: IsNormalize    =        0
  INTEGER              :: IsNormBnormal  =        0
  INTEGER              :: IsNormWeight   =        0       
  REAL                 :: weight_bnorm   =        1.000D+00
  REAL                 :: weight_bharm   =        0.000D+00
  REAL                 :: weight_tflux   =        0.000D+00
  REAL                 :: target_tflux   =        0.000D+00
  REAL                 :: weight_ttlen   =        0.000D+00
  REAL                 :: weight_specw   =        0.000D+00
  REAL                 :: weight_ccsep   =        0.000D+00

  REAL                 :: SD_tausta      =        0.000D-00
  REAL                 :: SD_tauend      =        1.000D-00
  REAL                 :: SD_tautol      =        1.000D-04
  INTEGER              :: SD_Nout        =      100        
  INTEGER              :: SD_savefreq    =        1                
 
  REAL                 :: NT_xtol        =        1.000D-04 
  REAL                 :: NT_eta         =        0.900D+00 
  REAL                 :: NT_stepmx      =        1.000D+05 

  INTEGER              :: case_postproc  =        0         
  REAL                 :: PP_odetol      =        1.000D-10 
  INTEGER              :: PP_Ppts        =      100         
  INTEGER              :: PP_Ptrj        =        8         
  REAL                 :: PP_phi         =        0.000D-00
  REAL                 :: PP_Rmin        =        0.000D-00
  REAL                 :: PP_Rmax        =        0.000D-00
  REAL                 :: PP_Zmin        =        0.000D-00
  REAL                 :: PP_Zmax        =        0.000D-00 
  REAL                 :: PP_bstol       =        1.000D-06 
  INTEGER              :: PP_bsnlimit    =   100000         
                                                         

  

  namelist / focusin /  IsQuiet       , &
                        IsSymmetric   , & 
                        case_surface  , & 
                        knotsurf      , &
                        Nteta         , &   
                        Nzeta         , &
                        case_init     , & 
                        case_coils    , &
                        Ncoils        , &
                        init_current  , &  
                        init_radius   , &
                        IsVaryCurrent , & 
                        IsVaryGeometry, & 
                        target_length , & 
                        NFcoil        , & 
                        Nseg          , &
                        case_optimizer, &
                        IsNormBnormal , &
                        IsNormalize   , &
                        IsNormWeight  , &
                        weight_bnorm  , &
                        weight_bharm  , &
                        weight_tflux  , &
                        target_tflux  , &
                        weight_ttlen  , &
                        weight_specw  , &
                        weight_ccsep  , &
                        SD_tausta     , &
                        SD_tauend     , &
                        SD_tautol     , &
                        SD_Nout       , &
                        SD_savefreq   , &        
                        NT_xtol       , & 
                        NT_eta        , & 
                        NT_stepmx     , & 
                        case_postproc , &
                        PP_odetol     , & 
                        PP_Ppts       , & 
                        PP_Ptrj       , & 
                        PP_phi        , & 
                        PP_Rmin       , &
                        PP_Rmax       , &
                        PP_Zmin       , &
                        PP_Zmax       , &
                        PP_bstol      , & 
                        PP_bsnlimit   
  
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
                             ds(:,:), xt(:,:), yt(:,:), zt(:,:), bn(:,:), tn(:,:), pb(:,:), &
                             Bx(:,:), By(:,:), Bz(:,:)
  end type toroidalsurface

  type arbitrarycoil
     INTEGER              :: NS, Ic, Lc, itype
     REAL                 :: I,  L, Lo, maxcurv
     REAL   , allocatable :: xx(:), yy(:), zz(:), xt(:), yt(:), zt(:), xa(:), ya(:), za(:), dd(:), &
                             Bx(:,:), By(:,:), Bz(:,:),  Ax(:,:), Ay(:,:), Az(:,:)
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

  INTEGER              :: Nfou=0, Nfp=0, NBnf=0
  INTEGER, allocatable :: bim(:), bin(:), Bnim(:), Bnin(:)
  REAL   , allocatable :: Rbc(:), Zbs(:), Rbs(:), Zbc(:), Bnc(:), Bns(:)
    
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsection{Packing and unpacking}
  INTEGER              :: Cdof, Ndof, nfixcur, nfixgeo, Tdof, iter
  REAL                 :: Inorm = one, Gnorm = one                !current and geometry normalizations;
  REAL   , allocatable :: xdof(:), dofnorm(:)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsection{Optimization}
  ! General target functions;
  INTEGER              :: itau
  REAL                 :: totalenergy, discretefactor
  REAL   , allocatable :: t1E(:), t2E(:,:), evolution(:,:), coilspace(:,:), deriv(:,:)
  ! Bn surface integration;
  REAL                 :: bnorm
  REAL   , allocatable :: t1B(:), t2B(:,:), bn(:,:), dBx(:,:,:)
  ! Bn reasonant harmoics;
  INTEGER              :: NBmn
  INTEGER, allocatable :: Bmnin(:), Bmnim(:)
  REAL                 :: bharm
  REAL   , allocatable :: t1H(:), t2H(:,:), Bmnc(:),Bmns(:), wBmn(:), tBmnc(:), tBmns(:), carg(:,:), sarg(:,:)
  ! Tflux error;
  INTEGER              :: isign = 1
  REAL                 :: tflux
  REAL   , allocatable :: t1F(:), t2F(:,:)
  ! Length constraint
  REAL                 :: ttlen
  REAL   , allocatable :: t1L(:), t2L(:,:)
  ! Coil-coil spearation
  REAL                 :: ccsep
  REAL   , allocatable :: t1C(:), t2C(:,:)
  ! Spectral condensation;
  REAL                 :: specw
  REAL   , allocatable :: t1S(:), t2S(:,:)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsection{Knotran parameters}
  REAL                 :: knotphase
  REAL   , allocatable :: xkc(:), xks(:), ykc(:), yks(:), zkc(:), zks(:)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsection{Oculus parameters}
  INTEGER              :: icoil

  type(biotsavart)     :: bsfield
  type(poincaredata)   :: poincare
  
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

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end module globals
