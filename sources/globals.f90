
!title (input) ! Description of global variables.

!latex \briefly{Defines input namelists and global variables, details of input namelist can be viwed at 
!latex \link{initial}.}

!latex \calledby{\link{focus}}
!latex \calls{\link{initial}}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

module globals
  
  implicit none
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!


  CHARACTER(10), parameter :: version='v0.17.00' ! version number


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  integer, parameter :: sp = selected_real_kind(6 , 37  )
  integer, parameter :: dp = selected_real_kind(15, 307 ) ! double precision
  integer, parameter :: qp = selected_real_kind(33, 4931)

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
  
  REAL, parameter      :: pi         =  acos(-1.0_dp)
  REAL, parameter      :: pi2        =  pi * two
  REAL, parameter      :: bsconstant =  1.0E-07   !biot-savart constant
  REAL, parameter      :: antibscont =  1.0E-07 / bsconstant
  REAL, parameter      :: mu0        =  2.0E-07 * pi2
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsection{IO units}
  INTEGER              :: ounit      =  6 ! output to screen
  INTEGER              :: eunit      =  0 ! error unit
  INTEGER              :: wunit      =  7 ! write unit
  INTEGER              :: runit      =  8 ! read  unit
  INTEGER              :: funit      =  9 ! general I/O
  INTEGER              :: bunit      = 10 ! backup unit for I/O
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  CHARACTER(100)   :: ext       ! extention
  CHARACTER(100)   :: inputfile ! input namelist
  CHARACTER(100)   :: hdf5file  ! hdf5 file
  CHARACTER(100)   :: out_coils ! output ext.coils file
  CHARACTER(100)   :: out_focus ! output ext.focus file
  CHARACTER(100)   :: out_harm  ! output harmonics file
  CHARACTER(100)   :: out_plasma  ! updated plasma boundary
  CHARACTER(100)   :: out_stable  ! output ext.stable file
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!latex \subsection{Input namelist: \type{focusin}}
  INTEGER              :: IsQuiet        =   0
  INTEGER              :: lim_out        =   1
  ! surface related      
  INTEGER              :: IsSymmetric    =   0 
  INTEGER              :: case_surface   =   0
  INTEGER,parameter    :: plasma_surf_boozer =   5
  REAL                 :: knotsurf       =   0.200D-00
  REAL                 :: ellipticity    =   0.000D+00
  INTEGER              :: Nteta          =   64           
  INTEGER              :: Nzeta          =   64  
  ! coil related
  INTEGER              :: case_init      =   0
  INTEGER              :: case_coils     =   1
  INTEGER              :: Ncoils         =   0
  REAL                 :: init_current   =   1.000D+06        
  REAL                 :: init_radius    =   1.000D+00
  INTEGER              :: IsVaryCurrent  =   1         
  INTEGER              :: IsVaryGeometry =   1         
  INTEGER              :: NFcoil         =   4         
  INTEGER              :: Nseg           =   128 
  ! normalizations            
  INTEGER              :: IsNormalize    =   1
  INTEGER              :: IsNormWeight   =   1
  ! Normal field error
  INTEGER              :: case_bnormal   =   0
  REAL                 :: weight_bnorm   =   1.000D+00
  ! Bmn harmonics
  REAL                 :: weight_bharm   =   0.000D+00
  INTEGER              :: bharm_jsurf    =   0
  ! Resonant Bn
  REAL                 :: weight_resbn   =   0.000D+00
  REAL                 :: target_resbn   =   0.000D+00
  INTEGER              :: resbn_m        =   1
  INTEGER              :: resbn_n        =   1
  INTEGER              :: ghost_use      =   1
  INTEGER              :: ghost_call     =   1
  INTEGER              :: ghost_once     =   0
  INTEGER              :: ghost_failreturn = 0
  INTEGER              :: rcflux_use     =   1
  INTEGER              :: period_Nseg    =   1024
  INTEGER              :: pfl_mult       =   1
  REAL                 :: orpfl          =   1.000D+00
  REAL                 :: ozpfl          =   0.000D+00
  REAL                 :: xrpfl          =   1.000D+00
  REAL                 :: xzpfl          =   0.000D+00
  REAL                 :: pfl_xtol       =   1.000D-08
  REAL                 :: rcflux_target  =   0.000D+00
  INTEGER              :: pfldim         =   1
  ! Island sensitivity
  REAL                 :: weight_dpsidr  =   0.000D+00
  INTEGER              :: dpsi_linear    =   1
  ! toroidal flux
  REAL                 :: weight_tflux   =   0.000D+00
  REAL                 :: target_tflux   =   0.000D+00
  ! coil length
  REAL                 :: weight_ttlen   =   0.000D+00
  INTEGER              :: case_length    =   1
  REAL                 :: target_length  =   0.000D+00
  REAL                 :: length_delta   =   0.000D+00
  ! coil curvature
  REAL                 :: weight_curv    =   0.000D+00
  INTEGER              :: case_curv      =   4 
  INTEGER              :: penfun_curv    =   1
  REAL                 :: curv_alpha     =   0.000D+00
  REAL                 :: curv_beta      =   2.000D+00
  REAL                 :: curv_gamma     =   2.000D+00
  REAL                 :: curv_sigma     =   0.000D+00
  REAL                 :: curv_k0        =   1.000D+01
  REAL                 :: curv_k1        =   0.000D+00
  INTEGER              :: curv_k1len     =   0
  ! straight function
  REAL                 :: weight_straight     =   0.000D+00
  REAL                 :: coeff_disp_straight =   0.50D+00
  INTEGER              :: case_straight       =   3 
  INTEGER              :: penfun_str    =   1
  REAL                 :: str_alpha     =   0.000D+00
  REAL                 :: str_beta      =   2.000D+00
  REAL                 :: str_gamma     =   2.000D+00
  REAL                 :: str_sigma     =   0.000D+00
  REAL                 :: str_k0        =   1.000D+01
  REAL                 :: str_k1        =   0.000D+00
  INTEGER              :: str_k1len     =   0
  ! coil torsion
  REAL                 :: weight_tors    =   0.000D+00
  INTEGER              :: case_tors      =   1
  INTEGER              :: penfun_tors    =   1
  REAL                 :: tors0          =   0.000D+00
  REAL                 :: tors_alpha     =   1.000D+00
  REAL                 :: tors_beta      =   2.000D+00
  REAL                 :: tors_gamma     =   1.000D+00
  ! coil nissin complexity
  REAL                 :: weight_nissin  =   0.000D+00
  INTEGER              :: penfun_nissin  =   1
  REAL                 :: nissin0        =   0.000D+00
  REAL                 :: nissin_alpha   =   1.000D+00
  REAL                 :: nissin_beta    =   2.000D+00
  REAL                 :: nissin_gamma   =   2.000D+00
  REAL                 :: nissin_sigma   =   0.000D+00
  ! coil-surface separation
  REAL                 :: weight_cssep   =   0.000D+00
  REAL                 :: cssep_factor   =   4.000D+00 
  INTEGER              :: case_cssep     =   1
  REAL                 :: mincssep       =   7.000D-01
  REAL                 :: cssep_alpha    =   1.000D+00
  REAL                 :: cssep_beta     =   2.000D+00
  REAL                 :: cssep_gamma    =   2.000D+00
  REAL                 :: cssep_sigma    =   0.000D+00
  ! coil-coil separation
  REAL                 :: weight_ccsep   =   0.000D+00
  REAL                 :: weight_inorm   =   1.000D+00
  REAL                 :: weight_gnorm   =   1.000D+00
  REAL                 :: weight_mnorm   =   1.000D+00
  REAL                 :: origin_surface_x =   0.000D+00
  REAL                 :: origin_surface_y  =   0.000D+00
  REAL                 :: origin_surface_z  =   0.000D+00
  INTEGER              :: penfun_ccsep   =   1
  INTEGER              :: ccsep_skip     =   0
  REAL                 :: r_delta        =   1.000D-01
  REAL                 :: ccsep_alpha    =   1.000D+01
  REAL                 :: ccsep_beta     =   2.000D+00
  REAL                 :: weight_specw   =   0.000D+00
  ! bnrom stochastic
  REAL                 :: weight_sbnorm  =   0.000D+00   ! Stochastic bnormal weight
  ! resbn stochastic
  REAL                 :: weight_sresbn  =   0.000D+00   ! Stochastic bnormal weight
  INTEGER              :: Npert          =   0           ! Number of coil perturbations
  INTEGER              :: Nmax           =   3           ! Max frequency of perturbations
  REAL                 :: sdelta         =   1.000D-02   ! Perturbation magnitude
  ! Multi-filament
  INTEGER              :: Nturns         =   4
  INTEGER              :: Npancakes      =   2
  INTEGER              :: AlphaPert      =   0
  INTEGER              :: Alphafd        =   1
  INTEGER              :: Nalpha         =   1000
!  INTEGER              :: NalphaOpt      =   10
  INTEGER              :: ResMult        =   10          ! This will be added in later
  REAL                 :: Nturn_space    =   1.000D-02
  REAL                 :: Npancake_space =   1.000D-02
  REAL                 :: trans_length   =   1.000D-01
  ! optimize controls
  INTEGER              :: case_optimize  =   0
  REAL                 :: fdiff_delta    =   1.000D-04
  REAL                 :: exit_tol       =   1.000D-04
  INTEGER              :: output_use     =   1
  ! differential flow
  INTEGER              :: DF_maxiter     =   0
  REAL                 :: DF_xtol        =   1.000D-08     
  REAL                 :: DF_tausta      =   0.000D+00
  REAL                 :: DF_tauend      =   1.000D+00               
  ! conjugate gradient
  INTEGER              :: CG_maxiter     =   0
  REAL                 :: CG_xtol        =   1.000D-08
  REAL                 :: CG_wolfe_c1    =   0.1
  REAL                 :: CG_wolfe_c2    =   0.9
  INTEGER              :: CG_maxiter_s   =   10
  REAL                 :: CG_xtol_s      =   1.000D-08
  ! levenberg-marquardt
  INTEGER              :: LM_maxiter     =   0
  REAL                 :: LM_xtol        =   1.000D-08
  REAL                 :: LM_ftol        =   1.000D-08
  REAL                 :: LM_factor      =   1.000D+02
  ! not used
  INTEGER              :: HN_maxiter     =   0
  REAL                 :: HN_xtol        =   1.000D-08
  REAL                 :: HN_factor      =   100.0
  ! not used
  INTEGER              :: TN_maxiter     =   0
  REAL                 :: TN_xtol        =   1.000D-08
  INTEGER              :: TN_reorder     =   0
  REAL                 :: TN_cr          =   0.1
  ! post-processing
  INTEGER              :: case_postproc  =   1
  INTEGER              :: save_freq      =   1
  INTEGER              :: save_coils     =   0 
  INTEGER              :: save_harmonics =   0
  INTEGER              :: save_filaments =   0
  INTEGER              :: update_plasma  =   0   
  INTEGER              :: save_stable    =   0 
  INTEGER              :: filforce       =   0   ! Calculate filament body forces
  INTEGER              :: calcfb         =   0   ! Calculate coil finite-build
  ! poincare plots
  REAL                 :: pp_phi         =  0.000D+00
  REAL                 :: pp_raxis       =  0.000D+00       
  REAL                 :: pp_zaxis       =  0.000D+00
  REAL                 :: pp_rmax        =  0.000D+00
  REAL                 :: pp_zmax        =  0.000D+00
  INTEGER              :: pp_ns          =  10
  INTEGER              :: pp_maxiter     =  1000
  INTEGER              :: pp_nsteps      =  1
  INTEGER              :: pp_nfp         =  1
  REAL                 :: pp_xtol        =  1.000D-06

  CHARACTER(100)   :: input_surf     = 'plasma.boundary'  ! surface file
  CHARACTER(100)   :: input_coils    = 'none'             ! input file for coils
  CHARACTER(100)   :: input_stable   = 'none'             ! input file for stable field lines
  CHARACTER(100)   :: input_harm     = 'target.harmonics' ! input target harmonics file
  CHARACTER(100)   :: limiter_surf   = 'none'             ! limiter surface
                                                         
  namelist / focusin / &
  IsQuiet       ,&
  lim_out       ,&
  IsSymmetric   ,&
  lim_out       ,&
  case_surface  ,&
  input_surf    ,&
  knotsurf      ,&
  ellipticity   ,&
  Nteta         ,&
  Nzeta         ,&
  case_init     ,&
  case_coils    ,&
  input_coils   ,&
  Ncoils        ,&
  init_current  ,&
  init_radius   ,&
  IsVaryCurrent ,&
  IsVaryGeometry,&
  NFcoil        ,&
  Nseg          ,&
  IsNormalize   ,&
  IsNormWeight  ,&
  weight_inorm  ,&
  weight_gnorm  ,&
  weight_mnorm  ,&
  weight_bnorm  ,&
  case_bnormal  ,&
  weight_resbn  ,&
  target_resbn  ,&
  resbn_m       ,&
  resbn_n       ,&
  ghost_use     ,&
  ghost_call    ,&
  ghost_once    ,&
  ghost_failreturn,&
  rcflux_use    ,&
  period_Nseg   ,&
  pfl_mult      ,&
  orpfl         ,&
  ozpfl         ,&
  xrpfl         ,&
  xzpfl         ,&
  pfl_xtol      ,&
  rcflux_target ,&
  pfldim        ,&
  weight_dpsidr ,&
  dpsi_linear   ,&
  weight_bharm  ,&
  bharm_jsurf   ,&
  input_harm    ,&
  weight_tflux  ,&
  target_tflux  ,&
  weight_ttlen  ,&
  case_length   ,&
  target_length ,&
  length_delta  ,&
  weight_curv   ,&
  case_curv     ,&
  penfun_curv   ,&
  curv_alpha    ,&
  curv_beta     ,&
  curv_gamma    ,&
  curv_sigma    ,&
  curv_k0       ,&
  curv_k1       ,&
  curv_k1len    ,&
  weight_straight, &
  coeff_disp_straight,&
  case_straight  , &
  penfun_str     ,&
  str_alpha    ,&
  str_beta     ,&
  str_gamma    ,&
  str_sigma    ,&
  str_k0       ,&
  str_k1       ,&
  str_k1len    ,&
  origin_surface_x,&
  origin_surface_y,&
  origin_surface_z,&
  weight_tors   ,&
  case_tors     ,&
  penfun_tors   ,&
  tors0         ,&
  tors_alpha    ,&
  tors_beta     ,&
  tors_gamma    ,&
  weight_nissin ,&
  penfun_nissin ,&
  nissin0       ,&
  nissin_alpha  ,&
  nissin_beta   ,&
  nissin_gamma  ,&
  nissin_sigma  ,&
  weight_cssep  ,&
  cssep_factor  ,&
  limiter_surf  ,&
  case_cssep    ,&
  mincssep      ,&
  cssep_alpha   ,&
  cssep_beta    ,&
  cssep_gamma   ,&
  cssep_sigma   ,&
  weight_ccsep  ,&
  penfun_ccsep  ,&
  ccsep_skip    ,&
  r_delta       ,&
  ccsep_alpha   ,&
  ccsep_beta    ,&
  weight_specw  ,&
  weight_sbnorm ,&
  weight_sresbn ,&
  Npert         ,&
  Nmax          ,&
  sdelta        ,&
  Nturns        ,&
  Npancakes     ,&
  AlphaPert     ,&
  Alphafd       ,&
  Nalpha        ,&
!  NalphaOpt     ,&
  ResMult       ,&
  Nturn_space   ,&
  Npancake_space,&
  trans_length  ,&
  case_optimize ,&
  fdiff_delta   ,&
  exit_tol      ,&
  DF_maxiter    ,&
  DF_xtol       ,&
  DF_tausta     ,&
  DF_tauend     ,&
  CG_maxiter    ,&
  CG_xtol       ,&
  CG_wolfe_c1   ,&
  CG_wolfe_c2   ,&
  CG_maxiter_s  ,&
  CG_xtol_s     ,&
  output_use    ,&
  LM_maxiter    ,&
  LM_xtol       ,&
  LM_ftol       ,&
  LM_factor     ,&
  HN_maxiter    ,&
  HN_xtol       ,&
  HN_factor     ,&
  TN_maxiter    ,&
  TN_xtol       ,&
  TN_reorder    ,&
  TN_cr         ,&
  case_postproc ,&
  save_freq     ,&
  save_coils    ,&
  save_harmonics,&
  save_filaments,&
  save_stable   ,&
  update_plasma ,&
  filforce      ,&
  calcfb        ,&
  pp_phi        ,&
  pp_raxis      ,&
  pp_zaxis      ,&
  pp_rmax       ,&
  pp_zmax       ,&
  pp_ns         ,&
  pp_maxiter    ,&
  pp_nsteps     ,&
  pp_nfp        ,&
  pp_xtol                            

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!latex  \subsection{MPI stuffs}
  INTEGER, PARAMETER   :: master=0
  INTEGER              :: myid, ncpu, myworkid, color, masterid, nmaster, nworker
  INTEGER              :: MPI_COMM_MASTERS, MPI_COMM_MYWORLD, MPI_COMM_WORKERS, MPI_COMM_FOCUS
  REAL                 :: machprec, vsmall, small, sqrtmachprec

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsection{surface and coils data}
  type toroidalsurface
     INTEGER              :: Nteta, Nzeta, Nfou=0, Nfp=0, NBnf=0
     REAL   , allocatable :: Rbc(:), Zbs(:), Rbs(:), Zbc(:), Bnc(:), Bns(:), Pmnc(:), Pmns(:)
     REAL   , allocatable :: xx(:,:), yy(:,:), zz(:,:), nx(:,:), ny(:,:), nz(:,:), &
                             xt(:,:), yt(:,:), zt(:,:), xp(:,:), yp(:,:), zp(:,:), &
                             ds(:,:), bn(:,:), pb(:,:), &
                             Bx(:,:), By(:,:), Bz(:,:)
     INTEGER, allocatable :: bim(:), bin(:), Bnim(:), Bnin(:)
     REAL                 :: vol, area
  end type toroidalsurface

  type arbitrarycoil
     INTEGER              :: NS, Ic=0, Lc=0, type=0, symm=0
     REAL                 :: I=zero,  L=zero, Lo, maxcurv, ox, oy, oz, mt, mp, Bt, Bz, avgcurv, &
                                   minlambda, maxs
     REAL   , allocatable :: xx(:), yy(:), zz(:), xt(:), yt(:), zt(:), xa(:), ya(:), za(:), &
                             xb(:), yb(:), zb(:), dl(:), dd(:), dpsidx(:), dpsidy(:), dpsidz(:), &
                             psx(:), psy(:), psz(:), Bxx(:), Byy(:), Bzz(:), Fx(:), Fy(:), Fz(:), &
                             nfbx(:), nfby(:), nfbz(:), bfbx(:), bfby(:), bfbz(:), curvature(:), straight(:), &
                             alpha(:), alphap(:), alphapp(:), alphadof(:,:), pertx(:,:), perty(:,:), pertz(:,:), &
                             pertxp(:,:), pertyp(:,:), pertzp(:,:)
     INTEGER, allocatable :: nxx(:), nyy(:), nzz(:)
     character(10)        :: name
  end type arbitrarycoil

  type FourierCoil
     INTEGER              :: NF
     REAL   , allocatable :: xc(:), xs(:), yc(:), ys(:), zc(:), zs(:), cmt(:,:), smt(:,:)
  end type FourierCoil

  type SplineCoil
     INTEGER              :: NCP,NT,lwrk,ier
     INTEGER,allocatable  :: iwrk(:)
     REAL                 :: FP
     REAL   , allocatable :: vect (:),basis_0(:,:),basis_1(:,:),basis_2(:,:),basis_3(:,:),Cpoints_fit(:),eval_points(:),wrk(:),weights(:), &
                              Data_points(:),Cpoints(:),db_dt(:,:),db_dt_2(:,:)
  end type SplineCoil

  type DegreeOfFreedom
     INTEGER              :: ND
     REAL   , allocatable :: xdof(:), xof(:,:), yof(:,:), zof(:,:), xtof(:,:), ytof(:,:), ztof(:,:)
  end type DegreeOfFreedom

  type ghostsurface
     INTEGER              :: NF_stable, NF_axis, Nseg_stable, Ndof_stable, iter_track
     INTEGER, allocatable :: on(:), xn(:), axisn(:)
     REAL   , allocatable :: osnc(:), osns(:), othetanc(:), othetans(:), xsnc(:), xsns(:), &
                             xthetanc(:), xthetans(:), axisrnc(:), axiszns(:), os(:), xs(:), &
                             otheta(:), xtheta(:), zeta(:), Ra(:), Za(:), osdot(:), xsdot(:), &
                             othetadot(:), xthetadot(:), ox(:), oy(:), oz(:), xx(:), xy(:), &
                             xz(:), of(:), og(:), xf(:), xg(:), oxdot(:), oydot(:), ozdot(:), &
                             xxdot(:), xydot(:), xzdot(:), obsups(:), obsupzeta(:), obsuptheta(:),&
                             xbsups(:), xbsupzeta(:), xbsuptheta(:), xdof_stable(:), &
                             dFdxdof_stable(:), xdof_stable_hold(:), dpsidx(:,:), dpsidy(:,:), &
                             dpsidz(:,:)
     REAL                 :: F
  end type ghostsurface
  
  type(arbitrarycoil)  , target, allocatable :: coil(:)  
  type(toroidalsurface), target, allocatable :: surf(:)
  type(FourierCoil)    , target, allocatable :: FouCoil(:)
  type(DegreeOfFreedom), target, allocatable :: DoF(:)
  type(ghostsurface)   , target, allocatable :: gsurf(:)
  type(SplineCoil)     ,         allocatable :: Splines(:)

  INTEGER              :: Nfp = 1, symmetry = 0, surf_Nfp = 1
!  INTEGER              :: symmetry = 0, surf_Nfp = 1
  INTEGER              :: plasma = 1, limiter = 1
  REAL   , allocatable :: cosnfp(:), sinnfp(:)
    
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsection{Packing and unpacking}
  INTEGER              :: Cdof, Ndof, nfixcur, nfixgeo, Tdof
  REAL                 :: Inorm = one, Gnorm = one, Mnorm = one   !current, geometry, and moment normalizations;
  REAL   , allocatable :: xdof(:), dofnorm(:)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsection{Optimization}
  ! General target functions;
  INTEGER              :: iout, Nouts, LM_iter, LM_mfvec
  INTEGER              :: ibnorm = 0, ibharm = 0, itflux = 0, ittlen = 0, icssep = 0, icurv = 0, istr=0,  iccsep = 0, itors = 0, inissin = 0 ! starting number
  INTEGER              :: mbnorm = 0, mbharm = 0, mtflux = 0, mttlen = 0, mcssep = 0, mcurv = 0, mstr=0,  mccsep = 0, mtors = 0, mnissin = 0 ! numbers of targets
  REAL                 :: chi, discretefactor, sumDE
  REAL   , allocatable :: t1E(:), t2E(:,:), evolution(:,:), coilspace(:,:), deriv(:,:)
  REAL   , allocatable :: LM_fvec(:), LM_fjac(:,:)
  LOGICAL              :: exit_signal = .False., LM_output = .False.
  ! Bn surface integration;
  REAL                 :: bnorm
  REAL   , allocatable :: t1B(:), t2B(:,:), bn(:,:)
  ! Bn resonant harmoics in Boozer coordinates
  REAL                 :: resbn, resbn_bnc, resbn_bns, psi
  REAL   , allocatable :: t1R(:), b1s(:), b1c(:)
  ! Bn resonant harmoics;
  INTEGER              :: NBmn
  INTEGER, allocatable :: Bmnin(:), Bmnim(:)
  REAL                 :: bharm, bharm_factor
  REAL   , allocatable :: t1H(:), t2H(:,:), Bmnc(:),Bmns(:), wBmn(:), tBmnc(:), tBmns(:), &
                          carg(:,:), sarg(:,:), iBmnc(:), iBmns(:)
  ! Tflux error;
  INTEGER              :: tflux_sign = -1 ! default theta : counter-clockwise
  REAL                 :: tflux, psi_avg
  REAL   , allocatable :: t1F(:), t2F(:,:)
  ! Length constraint
  REAL                 :: ttlen
  REAL   , allocatable :: t1L(:), t2L(:,:)
  ! Curvature constraint
  REAL                 :: curv
  REAL   , allocatable :: t1K(:), t2K(:,:)
  ! Average Torsion 
  REAL                 :: tors
  REAL   , allocatable :: t1T(:), t2T(:,:)
  ! nissin complexity
  REAL                 :: nissin
  REAL   , allocatable :: t1N(:), t2N(:,:)
  REAL                 :: str
  REAL   , allocatable :: t1Str(:), t2Str(:,:)
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
  ! Stochastic bnorm
  REAL                 :: bnormavg, bnormmax
  REAL   , allocatable :: t1Bavg(:)
  INTEGER              :: stoch_done = 0
  ! Stochastic resbn
  REAL                 :: resbnavg, abspsimax
  REAL, allocatable    :: stochpsi(:), stochpsipred(:)
  REAL   , allocatable :: t1Ravg(:)
  ! Island sensitivity
  REAL                 :: dpsidr
  INTEGER              :: pflsuc
  REAL   , allocatable :: t1Z(:)

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
  INTEGER,parameter    :: coil_type_spline = 5
  INTEGER,parameter    :: coil_type_multi = 4
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsection{Time counting}
  REAL                 :: tstart, tfinish, time_initialize, time_optimize, time_postproc

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsection{Miscellaneous}
  REAL                 :: tmpw_bnorm, tmpw_tflux ,tmpt_tflux, tmpw_ttlen, tmpw_specw, tmpw_ccsep, tmpw_bharm, & 
                          tmpw_curv, tmpw_tors, tmpw_nissin, tmpw_str
  REAL                 :: overlap = 0.0
                          !tmp weight for saving to restart file
  REAL, allocatable    :: mincc(:,:), coil_importance(:)
  INTEGER              :: ierr, astat

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  ! fieldline tracing
  REAL, ALLOCATABLE    :: XYZB(:,:,:), ppr(:,:), ppz(:,:), iota(:)
  INTEGER              :: tor_num, total_num, booz_mpol, booz_ntor, booz_mn
  LOGICAL              :: lboozmn = .false.
  INTEGER, ALLOCATABLE :: bmim(:), bmin(:)
  REAL, ALLOCATABLE    :: booz_mnc(:,:), booz_mns(:,:)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end module globals

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

module bnorm_mod
  ! contains some common variables used in subroutine bnormal
  ! allocating once and re-using them will save allocation time
  use globals, only : dp
  implicit none

  ! 0-order
  REAL, allocatable :: dBx(:,:), dBy(:,:), dBz(:,:), Bm(:,:), dAx(:,:), dAy(:,:), dAz(:,:)
  ! 1st-order
  REAL, allocatable :: dBn(:), dBm(:), d1B(:,:,:)

end module bnorm_mod

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
module bharm_mod
  ! contains some common variables used in subroutine bnormal
  ! allocating once and re-using them will save allocation time
  use globals, only : dp
  implicit none

  ! 0-order
  ! none for now; in future, others should be moved to here. 03/30/2019
  ! 1st-order
  REAL, allocatable :: dBc(:), dBs(:)

end module bharm_mod

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
module mgrid_mod
  use globals, only : dp, zero, pi2
  INTEGER :: NR = 101, NZ=101, NP=72, MFP=0
  REAL    :: Rmin=zero, Rmax=zero, Zmin=zero, Zmax=zero, Pmin=zero, Pmax=pi2
  namelist / mgrid / Rmin, Rmax, Zmin, Zmax, Pmin, Pmax, NR, NZ, NP
end module mgrid_mod
