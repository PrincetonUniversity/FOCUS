
!title (input) ! Description of input namelists.

!latex \briefly{Defines input namelists and global variables}

!l tex \calledby{\link{}}
!l tex \calls{\link{}}

!latex \tableofcontents

!latex \subsection{overview}
!latex \bi
!latex \item[1.] Here, and elsewhere, input variables are shown in \inputvar{red}. The input list is read from file and broadcast in \link{al00aa}.
!latex \ei

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

module kmodule
  
  use oculus, only : biotsavart, poincaredata
  
  implicit none
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
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
  REAL, parameter      :: bsconstant =  1.0E-04 ! mu0/4pi
  REAL, parameter      :: antibscont =  1.0E-07 / bsconstant
  REAL, parameter      :: mu0        =  2.0E-07 * pi2
  REAL, parameter      :: goldenmean =  1.618033988749895 ! golden mean = ( one + sqrt(five) ) / two ;    
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  INTEGER              :: ounit      =  0

  INTEGER              :: lunit      =  8
  INTEGER              :: funit      =  9
  INTEGER              :: nunit      = 10
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  CHARACTER            :: ext*100
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!latex \subsection{input list: \type{focusin}}
!latex \bi

  INTEGER              :: Idisplay    =       -1         !latex \item \inputvar{Idisplay      =        0        } : silent output; -1 = details ;
  INTEGER              :: Isymmetric  =        1         !latex \item \inputvar{Isymmetric    =        1        } : enforce stellarator symmetry;
  INTEGER              :: Itopology   =        0         !latex \item \inputvar{Itopology     =        0        } : selects knottedness of plasma:
  REAL                 :: knotsurf    =        0.200D-00 !latex \item \inputvar{knotsurf      =        0.200D-00} : radius of knotted plasma boundary;
  INTEGER              :: Linitialize =        0         !latex \item \inputvar{Linitialize   =        0        } : 
  REAL                 :: Rmaj        =        1.000D+00 !latex \item \inputvar{Rmaj          =        1.000D+00} : major radius of coils;
  REAL                 :: rmin        =        0.500D+00 !latex \item \inputvar{rmin          =        0.500D+00} : minor radius of coils;
  INTEGER              :: Ic          =        0         !latex \item \inputvar{Ic            =        0        } : 
  REAL                 :: Io          =        1.000D+00 !latex \item \inputvar{Io            =        1.000D+00} : 
  REAL                 :: Iw          =        1.000D+00 !latex \item \inputvar{Iw            =        1.000D+00} : 
  INTEGER              :: Lc          =        2         !latex \item \inputvar{Lc            =        0        } : logical flag controlling length weight; see \link{tlength};
  REAL                 :: Lo          =        1.000D+00 !latex \item \inputvar{Lo            =        1.000D+00} : 
  REAL                 :: Lw          =        1.000D+00 !latex \item \inputvar{Lw            =        1.000D+00} : 
  INTEGER              :: NFcoil      =        4         !latex \item \inputvar{NFcoil        =        4        } : Fourier harmonics for each coil;
  INTEGER              :: NDcoil      =      128         !latex \item \inputvar{NDcoil        =      128        } : discrete segments per coil;
  INTEGER              :: Loptimize   =        0         !latex \item \inputvar{Loptimize     =    -2/-1/0/1/2/3} : -1 and -2 are for testing the derivatives; 
                                                         !latex       1 old descent; 2 new descent; 3 Powell nonlinear equations slover; 4 Newton;
  INTEGER              :: Lnormalize  =        0         !latex \item \inputvar{Lnormalize    =        0/1      } : turn off/on normalizing weights;         
  REAL                 :: weight_bnorm=        1.000D+00 !latex \item \inputvar{weight\_bnorm =        1.000D+00} : weight for bnormal constraint; \link{bnormal}
  REAL                 :: weight_tflux=        0.500D+00 !latex \item \inputvar{weight\_tflux =        0.500D+00} : weight for toroidal flux constraint; \link{torflux}
  REAL                 :: target_tflux=        0.000D+00 !latex \item \inputvar{target\_tflux =        1.000D+00} : target toroidal flux; \link{torflux}
  REAL                 :: weight_ttlen=        0.000D+00 !latex \item \inputvar{weight\_ttlen =        0.000D+00} : weight for coil length; \link{tlength}
  REAL                 :: weight_eqarc=        0.000D+00 !latex \item \inputvar{weight\_eqarc =        1.000D+00} : weight for equal arc length constraint; \link{equarcl}
  REAL                 :: weight_ccsep=        0.000D+00 !latex \item \inputvar{weight\_ccsep =        0.000D+00} : weight for coil-coil separation       ; \link{coilsep}
  REAL                 :: tauend      =        1.000D-00 !latex \item \inputvar{tauend        =        1.000D-00} : artificial relaxtion ``time'', \link{evolve};
  REAL                 :: tautol      =        1.000D-04 !latex \item \inputvar{tautol        =        1.000D-04} : o.d.e. integration tolerance;
  INTEGER              :: Ntauout     =      100         !latex \item \inputvar{Ntauout       =      100        } : intermediate time steps; \link{evolve};
  INTEGER              :: Savfreq     =        1         !latex \item \inputvar{Savfreq       =        1        } : Saving Frequency; 1 means saving files for each step;
  INTEGER              :: Nteta       =       64         !latex \item \inputvar{Nteta         =       64        } : 
  INTEGER              :: Nzeta       =       64         !latex \item \inputvar{Nzeta         =       64        } : 
  REAL                 :: absacc      =        1.000D-08 !latex \item \inputvar{absacc        =        1.000D-08} :
  REAL                 :: absreq      =        1.000D-12 !latex \item \inputvar{absreq        =        1.000D-12} :
  REAL                 :: relreq      =        1.000D-01 !latex \item \inputvar{relreq        =        1.000D-01} :
  REAL                 :: xtol        =        1.000D-04 !latex \item \inputvar{xtol          =        0.000D+00} : E04LBF tolerance 10*sqrtmachprec;
  REAL                 :: eta         =        0.900D+00 !latex \item \inputvar{eta           =        0.900D+00} : E04LBF accurance rate (step ration);
  REAL                 :: stepmx      =        1.000D+05 !latex \item \inputvar{stepmx        =        1.000D+05} : E04LBF Euclidean distance between solution and starting;
  INTEGER              :: Lpoincare   =        0         !latex \item \inputvar{Lpoincare     =        0        } : to construct \Poincare plot;
                                                         !latex       \bi \item if \inputvar{Lpoincare} $> 0$, then the fieldline parameter 
                                                         !latex                 is the cylindrical toroidal angle, and so $B^\phi$ must not equal zero;
                                                         !latex           \item if \inputvar{Lpoincare} $=-1$, then the fieldline parameter is the length;
                                                         !latex       \ei
  REAL                 :: odetol      =        1.000D-10 !latex \item \inputvar{odetol      =        1.000D-10} : \Poincare plot, \link{pp00aa};
  INTEGER              :: Ppts        =      100         !latex \item \inputvar{Ppts        =      100        } : \Poincare plot, \link{pp00aa};
  INTEGER              :: Ptrj        =        8         !latex \item \inputvar{Ptrj        =        8        } : \Poincare plot, \link{pp00aa};
  REAL                 :: phi         =        0.000D-00 !latex \item \inputvar{phi         =        0.000D-00} : \Poincare plot, \link{pp00aa};
  REAL                 :: bstol       =        1.000D-06 !latex \item \inputvar{bstol       =        1.000D-06} : 
                                                         !latex       tolerance in Biot-Savart integral; passed to \oculus{bs00aa};
  INTEGER              :: bsnlimit    =   100000         !latex \item \inputvar{bsnlimit    =   100000        } : 
                                                         !latex       max. number of iterations used in Biot-Savart integral; passed to \oculus{bs00aa};
#ifdef FASHION
  REAL                 :: cen_cur     =   1.0E7          !central filament;
  REAL                 :: cen_zmin    =   -100.0D00      !lowest z coordinate;
  REAL                 :: cen_zmax    =    100.0D00      !uppest z coordinate;
!latex \ei
  
  namelist / focusin /   Idisplay                      , &
                         Isymmetric                    , &
                         Itopology                     , &
                         knotsurf                      , &
                         Linitialize                   , &
                         Rmaj                          , &
                         rmin                          , &
                         Ic                            , &
                         Io                            , &
                         Iw                            , &
                         Lc                            , &
                         Lo                            , &
                         Lw                            , &
                         NFcoil                        , &
                         NDcoil                        , &
                         cen_cur                       , &
                         cen_zmin                      , &
                         cen_zmax                      , &
                         Loptimize                     , &
                         Lnormalize                    , &
                         weight_bnorm                  , &
                         weight_tflux                  , &
                         target_tflux                  , &
                         weight_ttlen                  , &
                         weight_eqarc                  , &
                         weight_ccsep                  , &
                         tauend                        , &
                         tautol                        , &
                         Ntauout                       , &
                         Savfreq                       , &
                         Nteta                         , &
                         Nzeta                         , &
                         absacc                        , &
                         absreq                        , & ! redundant; 14 Apr 16;
                         relreq                        , & ! redundant; 14 Apr 16;
                         xtol                          , &
                         eta                           , &
                         stepmx                        , &
                         Lpoincare                     , &
                         odetol                        , &
                         Ppts                          , &
                         Ptrj                          , &
                         phi                           , &
                         bstol                         , &
                         bsnlimit                      
#else
    
  namelist / focusin /   Idisplay                      , &
                         Isymmetric                    , &
                         Itopology                     , &
                         knotsurf                      , &
                         Linitialize                   , &
                         Rmaj                          , &
                         rmin                          , &
                         Ic                            , &
                         Io                            , &
                         Iw                            , &
                         Lc                            , &
                         Lo                            , &
                         Lw                            , &
                         NFcoil                        , &
                         NDcoil                        , &
                         Loptimize                     , &
                         Lnormalize                    , &
                         weight_bnorm                  , &
                         weight_tflux                  , &
                         target_tflux                  , &
                         weight_ttlen                  , &
                         weight_eqarc                  , &
                         weight_ccsep                  , &
                         tauend                        , &
                         tautol                        , &
                         Ntauout                       , &
                         Savfreq                       , &
                         Nteta                         , &
                         Nzeta                         , &
                         absacc                        , &
                         absreq                        , & ! redundant; 14 Apr 16;
                         relreq                        , & ! redundant; 14 Apr 16;
                         xtol                          , &
                         eta                           , &
                         stepmx                        , &
                         Lpoincare                     , &
                         odetol                        , &
                         Ppts                          , &
                         Ptrj                          , &
                         phi                           , &
                         bstol                         , &
                         bsnlimit  
#endif
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  INTEGER              :: myid, ncpu
  REAL                 :: machprec, vsmall, small, sqrtmachprec
  CHARACTER            :: nodelabel*3
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  INTEGER              :: bmn, bNfp, nbf
  INTEGER, allocatable :: bim(:), bin(:), bnim(:), bnin(:)
  REAL   , allocatable :: Rbc(:), Zbs(:), Rbs(:), Zbc(:), bnc(:), bns(:)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  INTEGER              :: Ncoils, itime, iteta, jzeta, nrestart, itau, Cdof, Ndof, Tdof, Ndim, iter, nfixcur, nfixgeo
  REAL                 :: totalenergy, discretefactor
  LOGICAL              :: Langrange = .false.           ! flag for whether including langrange multipiler in DoFs; 08/17/2016

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  INTEGER              :: knotNF
  REAL                 :: knotphase
  REAL   , allocatable :: xkc(:), xks(:), ykc(:), yks(:), zkc(:), zks(:)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  REAL   , allocatable :: cmt(:,:), smt(:,:), shudson(:), newton(:)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  REAL   , allocatable :: evolution(:,:), coilspace(:,:), deriv(:,:)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  type arbitrarycoil
     INTEGER                :: N, D, Ic, Lc
     REAL                   :: I, Io, Iw, L, Lo, Lw, maxcurv
     REAL     , allocatable :: xc(:), xs(:), yc(:), ys(:), zc(:), zs(:), xx(:), yy(:), zz(:), xt(:), yt(:), zt(:), Bx(:,:), By(:,:), Bz(:,:),  &
                                             lmdc(:), lmds(:),                                xa(:), ya(:), za(:), Ax(:,:), Ay(:,:), Az(:,:)
     CHARACTER*10           :: name
  end type arbitrarycoil

  type toroidalsurface
     INTEGER              :: Nteta, Nzeta
     REAL   , allocatable :: xx(:,:), yy(:,:), zz(:,:), nx(:,:), ny(:,:), nz(:,:), ds(:,:), xt(:,:), yt(:,:), zt(:,:), bnt(:,:)
  end type toroidalsurface

  type(arbitrarycoil)  , allocatable :: coil(:)  
  type(toroidalsurface), allocatable :: surf(:)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! use Oculus; 04 Nov 1

  INTEGER              :: icoil

  type(biotsavart)     :: bsfield
  type(poincaredata)   :: poincare
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  INTEGER              :: isign = 1  ! sign symbol for flux
  REAL, allocatable    :: t1E(:,:), t2E(:,:,:,:), t1B(:,:), t2B(:,:,:,:), bn(:,:), bm(:,:), t1F(:,:), t2F(:,:,:,:), t1L(:,:), t2L(:,:,:,:), &
                          t1A(:,:), t2A(:,:,:,:), t1C(:,:), t2C(:,:,:,:), dlc(:,:,:), dls(:,:,:), n1E(:,:), n2E(:,:,:,:), tbn(:,:)
  REAL                 :: bnorm, tflux, ttlen, eqarc, ccsep

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! coils file data 
  REAL,    allocatable :: coilsX(:,:), coilsY(:,:), coilsZ(:,:), coilsI(:) !coilsI stores the true currents, coil%I stores scaled current;
  INTEGER, allocatable :: Nseg(:)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! seperation data
  REAL,    allocatable :: mincc(:,:)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!time benchmarking
  REAL                 :: tstart, tfinish

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!


!tmp weight for saving to restart file
  REAL                 :: tmpw_bnorm, tmpw_tflux ,tmpt_tflux, tmpw_ttlen, tmpw_eqarc, tmpw_ccsep

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  REAL,allocatable     :: SaveBx(:,:), SaveBy(:,:), SaveBz(:,:)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end module kmodule
