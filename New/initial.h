
!title (initialize) ! Read input file, initialize global variables.

!latex \briefly{Reads input files broadcasts, and allocates/initializes global variables.}

!latex \calledby{\link{focus}, \link{globals}}
!latex \calls{\link{}}

!latex \section{Input namelist}
!latex \bi
!latex \item \inputvar{IsQuiet        =       -1        } : control the output info; \\
!latex       {\bf $<$-1}: more details; \\ {\bf-1}: some details;\\  {\bf0}: moderate;\\ {\bf1}: concise;
!latex \item \inputvar{IsSymmetric    =        0        } : turn off stellarator symmetry (not ready);
!latex \par \begin{tikzpicture} \draw[dashed] (0,1) -- (10,1); \end{tikzpicture}

!latex \item \inputvar{case\_surface  =        0        } : select plasma shape;  \\
!latex       {\bf0}: stellarator; \\ {\bf1}: knotran; \\ {\bf2}: tokmak (not ready); \\
!latex        seen in \link{surface}.
!latex \item \inputvar{knotsurf       =        0.200D-00} : radius of knotted plasma boundary;
!latex \item \inputvar{Nteta          =       64        } : poloidal resolution of the surface;   
!latex \item \inputvar{Nzeta          =       64        } : toroidal resolution of the surface; 
!latex \par \begin{tikzpicture} \draw[dashed] (0,1) -- (10,1); \end{tikzpicture}

!latex \item \inputvar{case\_init     =        0        } : select initial coils source; \\
!latex       {\bf0}: toroidally place circular coils; \\
!latex       {\bf1}: from .ext.fo.coil.xxx ; \\
!latex       {\bf2}: from coils.ext; \\ seen in \link{rdcoils}
!latex \item \inputvar{case\_coils    =        1        } : select coil representation types; \\
!latex       {\bf0}: piece-wise linear (not ready); \\
!latex       {\bf1}:Fourier series seen in \link{coilfou};
!latex \item \inputvar{Ncoils         =        0        } : number of coils (valid for case\_coils 0 and 1)
!latex \item \inputvar{init\_current  =        1.000D+07} : initial coil current (for case\_coils 0);
!latex \item \inputvar{init\_radius   =        0.500D+00} : initial coil radius  (for case\_coils 0);
!latex \item \inputvar{IsVaryCurrent  =        1        } : Let all the coil currents varying (1) 
!latex \item \inputvar{IsVaryGeometry =        1        } : Let all the coil geometry varying (1) 
!latex \item \inputvar{target\_length  =       1.000D+00} : coil target length; if 0.0, auto calculated;
!latex                                                      seen in \link{coillen} 
!latex \item \inputvar{NFcoil         =        4        } : number of Fourier harmonics for coils;
!latex \item \inputvar{Nseg           =      128        } : discretizing number of segments in coils;
!latex \par \begin{tikzpicture} \draw[dashed] (0,1) -- (10,1); \end{tikzpicture}

!latex \item \inputvar{case\_optimizer =        0        } : select optimizing method; \\
!latex       {\bf-2}: check 2nd the derivatives; {\bf-1}: check 1st derivatives;  (\link{fdcheck}) \\
!latex       {\bf0}: no optimization; \\ {\bf 1}: steepest descent method; (\link{descent})\\
!latex       {\bf2}: Newton method (not ready)      
!latex \item \inputvar{IsNormalize    =        1        } : Normalize coil parameters (1); \link{coilfou}
!latex \item \inputvar{IsNormBnormal  =        0        } : Normalize Bnormal or not; \link{costfun} \\
!latex       {\bf 0}: no normalized surface integral; \\
!latex       {\bf 1}: surface ingtegral normalized with mod B; \\
!latex       {]bf 2}: Optimize Bn Fourier spectrum;
!latex \item \inputvar{IsNormWeight   =        0        } : Normalize weights (1); \link{costfun}
!latex \item \inputvar{weight\_bnorm  =        1.000D+00} : weight for Bn rms error; \link{bnormal}
!latex \item \inputvar{weight\_bharm  =        0.000D+00} : weight for Bn harmonics error; \link{bnftran}
!latex \item \inputvar{weight\_tflux  =        0.000D+00} : weight for toroidal flux error; \link{torflux}
!latex \item \inputvar{target\_tflux  =        0.000D+00} : target toroidal flux; if 0.0, auto calculated;
!latex \item \inputvar{weight\_ttlen  =        0.000D+00} : weight for coil length error; \link{coillen}
!latex \item \inputvar{weight\_specw  =        0.000D+00} : weight for spectral condensation;\link{specwid}
!latex \item \inputvar{weight\_ccsep  =        0.000D+00} : weight for coil-coil separation; \link{coilsep}
!latex \par \begin{tikzpicture} \draw[dashed] (0,1) -- (10,1); \end{tikzpicture}

!latex \item \inputvar{SD\_tausta     =        0.000D-00} : start   step in steepest descent method; 
!latex \item \inputvar{SD\_tauend     =        1.000D-00} : stoping step in steepest descent method; 
!latex \item \inputvar{SD\_tautol     =        1.000D-04} : tolerance for solving ODEs in SD method;
!latex \item \inputvar{SD\_Nout       =      100        } : total monitoring output;
!latex \item \inputvar{SD\_savefreq   =        1        } : output files saving frequency;        
!latex \par \begin{tikzpicture} \draw[dashed] (0,1) -- (10,1); \end{tikzpicture}
 
!latex \item \inputvar{NT\_xtol       =        1.000D-04} : tolerance for Newton method;  
!latex \item \inputvar{NT\_eta        =        0.900D+00} : eta in \link{gaussnt}; 
!latex \item \inputvar{NT\_stepmx     =        1.000D+05} : stepmx in \link{gaussnt}; 
!latex \par \begin{tikzpicture} \draw[dashed] (0,1) -- (10,1); \end{tikzpicture}
 
!latex \item \inputvar{case\_postproc =        0        } : select post-processing method;\\
!latex       {\bf0}: skip; \\ {\bf1}: poincare plot; \\{\bf2}: SPEC input;  \\ {\bf4}: write mgrid; \\
!latex       The three can be combined with each other; seen in \link{postpro}
!latex \item \inputvar{PP\_odetol     =        1.000D-10} : odetol for \Poincare plot; \oculus{pp00aa} 
!latex \item \inputvar{PP\_Ppts       =      100        } : number of line tracing iterations; 
!latex \item \inputvar{PP\_Ptrj       =        8        } : number of trajectories;
!latex \item \inputvar{PP\_phi        =        0.000D-00} : toroidal angle for the plot plane;
!latex \item \inputvar{PP\_Rmin       =        0.000D-00} : minimum R in PP / Rmin in mgrid;
!latex \item \inputvar{PP\_Rmax       =        0.000D-00} : maximum R in PP / Rmax in mgrid;
!latex \item \inputvar{PP\_Zmin       =        0.000D-00} : minimum Z in PP / Zmin in mgrid;
!latex \item \inputvar{PP\_Zmax       =        0.000D-00} : maximum Z in PP / Zmax in mgrid; 
!latex \item \inputvar{PP\_bstol      =        1.000D-06} : biot-savart tolerance; \oculus{bs00aa}
!latex \item \inputvar{PP\_bsnlimit   =   100000        } : biot-savart integral iteration limit; 
!latex \ei

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine Initial

  use globals
  implicit none
  include "mpif.h"


  LOGICAL :: exist
  INTEGER :: ierr, astat

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  machprec = epsilon(pi)         ! get the machine precision
  sqrtmachprec = sqrt(machprec)  ! sqrt of machine precision
  vsmall = ten * machprec        ! very small number
  small = thousand * machprec    ! small number

  !-------------read input namelist----------------------------------------------------------------------
  if( myid == 0 ) then ! only the master node reads the input; 25 Mar 15;
     call getarg(1,ext) ! get argument from command line
     write(ounit, '("initial : machine_prec   = ", ES12.5, " ; sqrtmachprec   = ", ES12.5,   &
          & " ; ext            = ", A100)') machprec, sqrtmachprec, ext
     inquire(file=trim(ext)//".input", EXIST=exist) ! inquire if ext.fo existed
  endif

  LlBCAST( exist, 1, 0 )
  FATAL( initial, .not.exist, input file ext.input not provided )

  if( myid == 0 ) then
     open(runit, file=trim(ext)//".input", status="unknown")
     read(runit, focusin)
     close(runit)
  endif ! end of if( myid == 0 ) ; 25 Mar 15;

  !-------------broadcast the namelist-------------------------------------------------------------------

  ClBCAST( ext           ,  100,  0 )
  IlBCAST( IsQuiet       ,    1,  0 )
  IlBCAST( IsSymmetric   ,    1,  0 )
  IlBCAST( case_surface  ,    1,  0 )
  RlBCAST( knotsurf      ,    1,  0 )
  IlBCAST( Nteta         ,    1,  0 )
  IlBCAST( Nzeta         ,    1,  0 )
  IlBCAST( case_init     ,    1,  0 )
  IlBCAST( case_coils    ,    1,  0 )
  IlBCAST( Ncoils        ,    1,  0 )
  RlBCAST( init_current  ,    1,  0 )
  RlBCAST( init_radius   ,    1,  0 )
  IlBCAST( IsVaryCurrent ,    1,  0 )
  IlBCAST( IsVaryGeometry,    1,  0 )
  RlBCAST( target_length ,    1,  0 )
  IlBCAST( NFcoil        ,    1,  0 )
  IlBCAST( Nseg          ,    1,  0 )
  IlBCAST( case_optimizer,    1,  0 )
  IlBCAST( IsNormalize   ,    1,  0 )
  IlBCAST( IsNormBnormal ,    1,  0 )
  IlBCAST( IsNormWeight  ,    1,  0 )
  RlBCAST( weight_bnorm  ,    1,  0 )
  RlBCAST( weight_bharm  ,    1,  0 )
  RlBCAST( weight_tflux  ,    1,  0 )
  RlBCAST( target_tflux  ,    1,  0 )
  RlBCAST( weight_ttlen  ,    1,  0 )
  RlBCAST( weight_specw  ,    1,  0 )
  RlBCAST( weight_ccsep  ,    1,  0 )
  RlBCAST( SD_tausta     ,    1,  0 )
  RlBCAST( SD_tauend     ,    1,  0 )
  RlBCAST( SD_tautol     ,    1,  0 )
  IlBCAST( SD_Nout       ,    1,  0 )
  IlBCAST( SD_savefreq   ,    1,  0 )
  RlBCAST( NT_xtol       ,    1,  0 )
  RlBCAST( NT_eta        ,    1,  0 )
  RlBCAST( NT_stepmx     ,    1,  0 )
  IlBCAST( case_postproc ,    1,  0 )
  RlBCAST( PP_odetol     ,    1,  0 )
  IlBCAST( PP_Ppts       ,    1,  0 )
  IlBCAST( PP_Ptrj       ,    1,  0 )
  RlBCAST( PP_phi        ,    1,  0 )
  RlBCAST( PP_Rmin       ,    1,  0 )
  RlBCAST( PP_Rmax       ,    1,  0 )
  RlBCAST( PP_Zmin       ,    1,  0 )
  RlBCAST( PP_Zmax       ,    1,  0 )
  RlBCAST( PP_bstol      ,    1,  0 )
  IlBCAST( PP_bsnlimit   ,    1,  0 )

  !-------------assign Poincare plot stuffs--------------------------------------------------------------

  bsfield%tol = PP_bstol ; bsfield%N = PP_bsnlimit 
  ! bsfield is global; passed through to Oculus:bs00aa; 30 Oct 15;

  !-------------show the name list for checking----------------------------------------------------------

  if(myid == 0 .and. IsQuiet <= -1) then ! Not quiet to output more informations;

     write(ounit, *) "-----------INPUT NAMELIST------------------------------------"

     write(ounit, '("initial : IsQuiet        = ", I5,7X , " ; IsSymmetric    = ", I5,7X " ; " )') &
          IsQuiet, IsSymmetric

     write(ounit, '(     8X,": case_surface   = ", I5,7X , " ; knotsurf       = ", ES12.5,   &
          &                " ; Nteta          = ", I5,7X , " ; Nzeta          = ", I5,7X " ; " )') &
          case_surface, knotsurf, Nteta, Nzeta

     write(ounit, '(     8X,": case_init      = ", I5,7X , " ; case_coils     = ", I5,7X ,   &
          &                " ; Ncoils         = ", I5,7X , " ; init_current   = ", ES12.5,   &
          &                " ; init_radius    = ", ES12.5" ; " )') &
          case_init, case_coils, Ncoils, init_current, init_radius
     write(ounit, '(     8X,": IsVaryCurrent  = ", I5,7X , " ; IsVaryGeometry = ", I5,7X ,   &
          &                " ; target_lenth   = ", ES12.5,   &
          &                " ; NFcoil         = ", I5,7X , " ; Nseg           = ", I5,7X " ; " )') &
          IsvaryCurrent, IsvaryGeometry, target_length, NFcoil, Nseg  

     write(ounit, '(     8X,": case_optimizer = ", I5,7X , " ; IsNormalize    = ", I5,7X ,   &
          &                " ; IsNormBnormal  = ", I5,7X , " ; IsNormWeight   = ", I5,7X ,   &
          &                " ; weight_bnorm   = ", ES12.5" ; ")') &
          case_optimizer, IsNormalize, IsNormBnormal, IsNormWeight, weight_bnorm
     write(ounit, '(     8X,": weight_tflux   = ", ES12.5, " ; target_tflux   = ", ES12.5,   &
          &                " ; weight_ttlen   = ", ES12.5, " ; weight_specw   = ", ES12.5,   &
          &                " ; weight_ccsep   = ", ES12.5)')    &
          weight_tflux, target_tflux,weight_ttlen, weight_specw, weight_ccsep

     write(ounit, '(     8X,": SD_tausta      = ", ES12.5, " ; SD_tauend      = ", ES12.5,   &
          &                " ; SD_tautol      = ", ES12.5,                                   &
          &                " ; SD_Nout        = ", I5,7X , " ; SD_savefreq    = ", I5,7X " ; " )') &
          SD_tausta, SD_tauend, SD_tautol, SD_Nout, SD_savefreq

     write(ounit, '(     8X,": NT_xtol        = ", ES12.5, " ; NT_eta         = ", ES12.5,   &
          &                " ; NT_stepmx      = ", ES12.5" ; ")') NT_xtol, NT_eta, NT_stepmx       

     write(ounit, '(     8X,": case_postproc  = ", I5,7X , " ; PP_odetol      = ", ES12.5,   &
          &                " ; PP_Ppts        = ", I5,7X , " ; PP_Ptrj        = ", I5,7X ,   &
          &                " ; PP_bstol       = ", ES12.5, " ; PP_bsnlimit    = ", I10,2X" ; ")') &
          case_postproc, PP_odetol, PP_Ppts, PP_ptrj, PP_bstol, PP_bsnlimit
     write(ounit, '(     8X,": PP_phi         = ", ES12.5, " ; PP_Rmin        = ", ES12.5,   &
          &                " ; PP_Rmax        = ", ES12.5, " ; PP_Zmin        = ", ES12.5,   &
          &                " ; PP_Zmax        = ", ES12.5" ; ")')    &
          PP_phi, PP_Rmin, PP_Rmax, PP_Zmin, PP_Zmax                         

  endif

  !-------------Error checking--------------------------------------------------------------------------- 

  FATAL( initial, IsSymmetric /= 0, Stellarator symmetry is not activated yet)

  select case( case_surface )
  case( 0 )
     inquire( file='plasma.boundary', exist=exist )
     FATAL( initial, .not.exist, plasma boundary file plasma.boundary not provided )
  case( 1 )
     FATAL( initial, knotsurf < zero, illegal)
     !case( 2 )
     ! reading Tokamak plasma
  case default
     FATAL( initial, .true., selected surface type is not supported )
  end select

  FATAL( initial, Nteta   <=   0, illegal surface resolution )
  FATAL( initial, Nzeta   <=   0, illegal surface resolution )

  select case( case_init )
  case( 0 )
     FATAL( initial, Ncoils < 1, should provide the No. of coils)
     FATAL( initial, init_current < zero, invalid coil current)
     FATAL( initial, init_radius < zero, invalid coil radius)
  case( 1 )
     !FATAL( initial, Ncoils < 0, invalid negative value)
  case(-1 )
     inquire( file="coils."//trim(ext), exist=exist )
     FATAL( initial, .not.exist, coils file coils.ext not provided )
  case default
     FATAL( initial, .true., selected case_init is not supported )
  end select

  FATAL( initial, case_coils /= 1, only fourier representation is valid )
  FATAL( initial, NFcoil <= 0    , no enough harmonics )
  FATAL( initial, Nseg   <= 0    , no enough segments  )

  select case( case_optimizer )
  case(-2 )
  case(-1 )
  case( 0 )
  case( 1 )
     FATAL( initial, SD_tausta <  zero     , illegal )
     FATAL( initial, SD_tauend <= SD_tausta, illegal )
     FATAL( initial, SD_tautol <= zero     , illegal )
     FATAL( initial, SD_Nout   <=    0     , illegal )
 !case( 2 )
  case default
     FATAL( initial, .true., selected optimizer is not supported )
  end select

  FATAL( initial, weight_bnorm  < zero, illegal )
  FATAL( initial, weight_tflux  < zero, illegal )
  FATAL( initial, weight_ttlen  < zero, illegal )
  FATAL( initial, weight_specw  < zero, illegal )
  FATAL( initial, weight_ccsep  < zero, illegal )

  if( mod(case_postproc, 2) /= 0 ) then
     FATAL( initial, PP_odetol   <=  zero, illegal )
     FATAL( initial, PP_Ppts     <   0   , illegal )
     FATAL( initial, PP_Ptrj     <   0   , illegal )
     FATAL( initial, PP_bstol    <=  zero, illegal )
     FATAL( initial, PP_bsnlimit <=     0, illegal )
  endif

  FATAL( initial, ncpu >= 1000 , too macy cpus, modify nodelabel)
  write(nodelabel,'(i3.3)') myid ! nodelabel is global; 30 Oct 15;

  discretefactor = (pi2/Nteta) * (pi2/Nzeta)
  
  !save weights before normalized
  tmpw_bnorm = weight_bnorm
  tmpw_bharm = weight_bharm
  tmpw_tflux = weight_tflux
  tmpt_tflux = target_tflux
  tmpw_ttlen = weight_ttlen
  tmpw_specw = weight_specw
  tmpw_ccsep = weight_ccsep

  return

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine initial
