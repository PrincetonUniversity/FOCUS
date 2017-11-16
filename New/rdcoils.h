
!title (coils) ! Initialize coils data using Fourier series.

!latex \briefly{Initialize the coils data with Fourier series.}

!latex \calledby{\link{focus}}
!latex \calls{\link{}}

!latex \subsection{overview}
!latex \bi
!latex \item[1.] If \inputvar{case\_coils=1}, then the Fourier series will be used for represent the coils.
!latex \item[2.] The basic equations about the Fourier representation is,
!latex  \be \label{fouriercoil}
!latex  x & = & X_{c,0} + \sum_ {n=1}^{N} \left [X_{c,n} \cos(nt)  +  X_{s,n} \sin(nt) \right ], \\ 
!latex  y & = & Y_{c,0} + \sum_ {n=1}^{N} \left [Y_{c,n} \cos(nt)  +  Y_{s,n} \sin(nt) \right ], \\
!latex  z & = & Z_{c,0} + \sum_ {n=1}^{N} \left [Z_{c,n} \cos(nt)  +  Z_{s,n} \sin(nt) \right ], 
!latex  \ee
!latex  \ei
!latex  \subsection{Initilization}
!latex  There are several ways to initialize the coils data.
!latex  \bi
!latex  \item[1.] \inputvar{case\_init = 0} : Toroidally placing \inputvar{Ncoils} circular coils with a 
!latex             radius of \inputvar{init\_radius} and current of \inputvar{init\_current}. The $i$th coil 
!latex             is placed at $\z = \frac{i-1}{Ncoils} \frac{2\pi}{Nfp}$.
!latex  \item[2.] \inputvar{case\_init = 1} : Read coils data from {\bf .ext.coil.xxx} files. xxx can vary 
!latex             from $001$ to $999$. Each file has such a format. \red{This is the most flexible way, and
!latex             each coil can be different.}            
!latex  \begin{raw}
!latex  #type of coils;   name
!latex      1      "Module 1"
!latex  #  Nseg       I   If  L   Lf  Lo
!latex     128  1.0E+07  0 6.28  1 3.14
!latex  # NFcoil
!latex  1
!latex  # Fourier harmonics for coils ( xc; xs; yc; ys; zc; zs)
!latex  3.00 0.30
!latex  0.00 0.00
!latex  0.00 0.00
!latex  0.00 0.00
!latex  0.00 0.00
!latex  0.00 0.30
!latex  \end{raw}
!latex  \ei
!latex  \bi
!latex  \item[3.] \inputvar{case\_init = -1} : Get coils data from a standard coils.ext file and 
!latex           then Fourier decomposed (normal Fourier tansformation and truncated with $NFcoil$ harmonics)
!latex  \ei
!latex   \subsection{Discretization}
!latex   \bi
!latex   \item[1.] Discretizing the coils data involves massive triangular functions in nested loops.
!latex   As shown in  \Eqn{fouriercoil}, the outside loop is for different discrete points and for each
!latex   point, a loop is needed to get the summation of the harmonics.
!latex   \item[2.] To avoid calling triangular functions every operations, it's a btter idea to allocate
!latex   the public triangular arrays.
!latex   \be
!latex   cmt(iD, iN) = \cos(iN \ \frac{iD}{D_i} 2\pi); iD = 0, coil(icoil)\%D; \ iN = 0,  coil(icoil)\%N \\
!latex   smt(iD, iN) = \sin(iN \ \frac{iD}{D_i} 2\pi); iD = 0, coil(icoil)\%D; \ iN = 0,  coil(icoil)\%N
!latex   \ee
!latex   \item[3.] Using the concept of vectorization, we can also finish this just through matrix 
!latex   operations. This is in \subroutine{fouriermatrix}.
!latex   \begin{raw}
!latex   subroutine fouriermatrix(xc, xs, xx, NF, ND)
!latex   nn(0:NF, 1:1) : matrix for N; iN
!latex   tt(1:1, 0:ND) : matrix for angle; iD/ND*2pi
!latex   nt(0:NF,0:ND) : grid for nt; nt = matmul(nn, tt)
!latex   xc(1:1, 0:NF) : cosin harmonics;
!latex   xs(1:1, 0:NF) : sin harmonics;
!latex   xx(1:1, 0:ND) : returned disrecte points;
!latex   
!latex   xx = xc * cos(nt) + xs * sin(nt)
!latex   \end{raw}
!latex   \item[4.] Actually, in real tests, the new method is not so fast. And parallelizations are actually
!latex   slowing the speed, both for the normal and vectorized method. 
!latex   \ei

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine rdcoils
  
  use globals, only : zero, one, half, pi, pi2, myid, ounit, tstart, &
                      NFcoil, Ncoils, IsVaryCurrent, IsVaryGeometry, &
                      runit, coilfile, discretecurve, &
                      Ns, init_current, init_radius, coil, Initialize, Isurface, &
                      cmt, smt
  
  implicit none
  
  include "mpif.h"
  
  LOGICAL              :: exist
  INTEGER              :: icoil, maxnseg, ifirst, NF, itmp, ip, icoef, jj, kk, ierr, astat, lNF, maxNF, mm
  REAL                 :: zeta, totalcurrent, Ro, Zo, r1, r2, z1, z2, tt
  REAL                 :: ax(1:3), at(1:3), az(1:3), xx(1:3), xs(1:3), xt(1:3), xz(1:3), xtt(1:3), xtz(1:3), xzz(1:3), v1(1:3), v2(1:3), w1(1:3), w2(1:3)
  REAL                 :: rdummy, tnow
  REAL   , allocatable :: coilsX(:,:), coilsY(:,:), coilsZ(:,:)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
    
  select case( Initialize )
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
  case( -1 ) ! Initialize = -1 : 04 Sep 17;
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
   FATAL( rdcoils, .true., this option is under reconstruction -- contact shudson@pppl.gov to reimplement )
   
!  if (myid == 0) then
!   write(ounit,'("rdcoils : reading coils data from "A)') inpcoils  ! different from coilsfile
!   call readcoils(inpcoils, maxnseg)
!   write(ounit,'("rdcoils : Read ",i6," coils in "A" ;")') Ncoils, trim(inpcoils)
!  endif
   
!  IlBCAST( Ncoils   ,      1, 0 )
!  IlBCAST( maxnseg  ,      1, 0 )

!  if( .not. allocated(coilsX) ) then !allocate arrays on other nodes;
!   SALLOCATE( coilsX, (1:maxnseg,1:Ncoils), zero )    
!   SALLOCATE( coilsY, (1:maxnseg,1:Ncoils), zero )
!   SALLOCATE( coilsZ, (1:maxnseg,1:Ncoils), zero )
!   SALLOCATE( coilsI, (          1:Ncoils), zero )
!   SALLOCATE( coilseg,(          1:Ncoils),    0 )
!   SALLOCATE( coilname,(         1:Ncoils),   '' )
!  endif
   
!  RlBCAST( coilsX, maxnseg*Ncoils, 0 )
!  RlBCAST( coilsY, maxnseg*Ncoils, 0 )
!  RlBCAST( coilsZ, maxnseg*Ncoils, 0 )
!  RlBCAST( coilsI,         Ncoils, 0 )   
!  IlBCAST( coilseg,        Ncoils, 0 )
!  ClBCAST( coilname,       Ncoils, 0 )
   
!  allocate(    coil(1:Ncoils) )
!  allocate( FouCoil(1:Ncoils) )
!  allocate(     DoF(1:Ncoils) )
   
!  Ncoils = Ncoils / Npc ! Ncoils changed to unique number of coils;
   
!  icoil = 0
   
!  do icoil = 1, Ncoils
    
!   coil(icoil)%NS =  Nseg  
!   coil(icoil)%I  =  coilsI(icoil)
!   coil(icoil)%If =  IsVaryCurrent
!   coil(icoil)%L  =  target_length ! irrelevant until re-computed;
!   coil(icoil)%Lf =  IsVaryGeometry
!   coil(icoil)%Lo =  target_length
!   coil(icoil)%name = trim(coilname(icoil))
    
!   FATAL( rdcoils, coil(icoil)%If < 0 .or. coil(icoil)%If > 1, illegal )
!   FATAL( rdcoils, coil(icoil)%Lf < 0 .or. coil(icoil)%Lf > 1, illegal )
!   FATAL( rdcoils, coil(icoil)%Lo < zero                     , illegal )
    
!   if(coil(icoil)%If == 0) Nfixcur = Nfixcur + 1
!   if(coil(icoil)%Lf == 0) Nfixgeo = Nfixgeo + 1
    
!   FouCoil(icoil)%NF = NFcoil
    
!   NF = NFcoil  ! alias
    
!   SALLOCATE( FouCoil(icoil)%xc, (0:NF), zero )
!   SALLOCATE( FouCoil(icoil)%xs, (0:NF), zero )
!   SALLOCATE( FouCoil(icoil)%yc, (0:NF), zero )
!   SALLOCATE( FouCoil(icoil)%ys, (0:NF), zero )
!   SALLOCATE( FouCoil(icoil)%zc, (0:NF), zero )
!   SALLOCATE( FouCoil(icoil)%zs, (0:NF), zero )        
    
!   call fourier( coilsX(1:coilseg(icoil),icoil), Foucoil(icoil)%xc, Foucoil(icoil)%xs, coilseg(icoil), NF)
!   call fourier( coilsY(1:coilseg(icoil),icoil), Foucoil(icoil)%yc, Foucoil(icoil)%ys, coilseg(icoil), NF)
!   call fourier( coilsZ(1:coilseg(icoil),icoil), Foucoil(icoil)%zc, Foucoil(icoil)%zs, coilseg(icoil), NF)
    
!  enddo

!  DALLOCATE( coilsX )
!  DALLOCATE( coilsY )
!  DALLOCATE( coilsZ )
!  DALLOCATE( coilsI )
!  DALLOCATE( coilseg)
!  DALLOCATE(coilname)
   
!  coil(1:Ncoils)%itype = case_coils
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
  case( 0 ) ! Initialize = 0 : read fourier harmonics of coils from file; 29 Sep 17;
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
   tnow = MPI_WTIME()
   
   if( myid.eq.0 ) write(ounit,'("rdcoils : "f10.1" : reading ",A," ;")') tnow-tstart, trim(coilfile)
   
   inquire( file=trim(coilfile), exist=exist )
   
   if( exist ) then
    if( myid.eq.0 ) open( runit, file=trim(coilfile), status="old", action='read' )
   else
    inquire( file=".fo.coils", exist=exist ) ! backup coils file;
    FATAL( rdaxis , .not.exist, coilfile or .fo.coils does not exist )
    if( myid.eq.0 ) open( runit, file=".fo.coils", status="old", action='read' )
   endif
   
   if( myid.eq.0 ) then
    read( runit,* )
    read( runit,* ) Ncoils
    write(ounit,'("rdcoils : " 10x " : Ncoils =",i4," ; Nsegments =",i4," ; NFcoil =",i3," ;")') Ncoils, Ns, NFcoil
   endif
   
   IlBCAST( Ncoils, 1, 0 )
   
   allocate( coil(1:Ncoils) ) ! coil is type( spacecurve );
   
   do icoil = 1, Ncoils
    
    if( myid.eq.0 ) then
     read( runit,* )
     read( runit,* )
     read( runit,* ) coil(icoil)%itype, coil(icoil)%name
    endif
    
    IlBCAST( coil(icoil)%itype,  1,  0 )
    ClBCAST( coil(icoil)%name , 10,  0 )
    
    if( myid.eq.0 ) then
     read( runit,* )
     read( runit,* ) coil(icoil)%NS, coil(icoil)%I, coil(icoil)%Ifree, rdummy, coil(icoil)%Lfree! coil(icoil)%Lo
    endif
    
    IlBCAST( coil(icoil)%NS   , 1, 0 )
    RlBCAST( coil(icoil)%I    , 1, 0 )
    IlBCAST( coil(icoil)%Ifree, 1, 0 )
    IlBCAST( coil(icoil)%Lfree, 1, 0 )
   !RlBCAST( coil(icoil)%Lo   , 1, 0 )
    
    FATAL( rdcoils, coil(icoil)%NS    < 0                           , illegal )
    FATAL( rdcoils, coil(icoil)%Ifree < 0 .or. coil(icoil)%Ifree > 1, illegal )
    FATAL( rdcoils, coil(icoil)%Lfree < 0 .or. coil(icoil)%Lfree > 1, illegal )
   !FATAL( rdcoils, coil(icoil)%Lo    < zero                        , illegal )     
    
    if( myid.eq.0 ) then
     read( runit,* )
     read( runit,* ) coil(icoil)%NF
    endif
    
    IlBCAST( coil(icoil)%NF, 1, 0 )
    
    FATAL( rdcoils, coil(icoil)%NF < 0, illegal )
    
    if( NFcoil.gt.0 ) then ; NF = NFcoil         ! override Fourier resolution;
    else                   ; NF = coil(icoil)%NF
    endif

    coil(icoil)%Ifree = IsVaryCurrent  ! override variation of current ;
    coil(icoil)%I     = init_current   ! override              current ;
    coil(icoil)%Lfree = IsVaryGeometry ! override variation of geometry ;
    
    SALLOCATE( coil(icoil)%xc, (0:NF), zero ) ! this is inside if( myid.eq.0 );
    SALLOCATE( coil(icoil)%xs, (0:NF), zero )
    SALLOCATE( coil(icoil)%yc, (0:NF), zero )
    SALLOCATE( coil(icoil)%ys, (0:NF), zero )
    SALLOCATE( coil(icoil)%zc, (0:NF), zero )
    SALLOCATE( coil(icoil)%zs, (0:NF), zero )
    
    lNF = min( NF, coil(icoil)%NF )
    
    if( myid.eq.0 ) then
     read( runit,* )
     read( runit,* ) coil(icoil)%xc(0:lNF)
     read( runit,* ) coil(icoil)%xs(0:lNF)
     read( runit,* ) coil(icoil)%yc(0:lNF)
     read( runit,* ) coil(icoil)%ys(0:lNF)
     read( runit,* ) coil(icoil)%zc(0:lNF)
     read( runit,* ) coil(icoil)%zs(0:lNF)
    endif
    
    RlBCAST( coil(icoil)%xc(0:lNF), 1+lNF,  0 )
    RlBCAST( coil(icoil)%xs(0:lNF), 1+lNF,  0 )
    RlBCAST( coil(icoil)%yc(0:lNF), 1+lNF,  0 )
    RlBCAST( coil(icoil)%ys(0:lNF), 1+lNF,  0 )
    RlBCAST( coil(icoil)%zc(0:lNF), 1+lNF,  0 )
    RlBCAST( coil(icoil)%zs(0:lNF), 1+lNF,  0 )
    
    coil(icoil)%NF = NF ! over-ride Fourier resolution;
    
   enddo ! end of do icoil;
   
   if( myid.eq.0 ) close( runit )
   
   write( ounit,'("rdcoils : " 10x " : I ="99es09.02" ;")') (/ ( coil(icoil)%I, icoil = 1, Ncoils ) /)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
  case( 1 ) ! Initialize = 1 : construct initial coil set ; 29 Sep 17;
      
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
   if( myid.eq.0 ) then
    
    FATAL( rdcoils, NFcoil.lt.0, illegal )
    
    select case( Isurface )
    case( 0 ) ; write(ounit,1000) Ncoils, NFcoil, Ns, init_current
    case( 1 ) ; write(ounit,1000) Ncoils, NFcoil, Ns, init_current, init_radius
    end select
    
1000 format("rdcoils : ", 10x ," : initializing ",i4," coils ; NF =",i3," ; Ns =",i4," ; initial current =",es13.5," ; ":"initial radius =",f6.2," ;")
    
   endif ! end of if( myid.eq.0 );
   
   allocate( coil(1:Ncoils) )
   
   select case( Isurface )
   case( 0 ) ! coilsX,Y,Z not required if fourier harmonics of coils are read from file; 29 Sep 17;
   case( 1 )
    SALLOCATE( coilsX, (1:Ns,1:Ncoils), zero ) ! THIS IS VERY CLUMSY;
    SALLOCATE( coilsY, (1:Ns,1:Ncoils), zero )
    SALLOCATE( coilsZ, (1:Ns,1:Ncoils), zero )
   case default
    FATAL( rdcoils, .true., selected Isurface not supported )
   end select
   
   do icoil = 1, Ncoils
    
    coil(icoil)%NF    = NFcoil    
    coil(icoil)%I     = init_current
    coil(icoil)%Ifree = IsVaryCurrent
    coil(icoil)%Lfree = IsVaryGeometry
   !coil(icoil)%Lo    = zero ! target_length ! 12 Nov 17;
    
    write(coil(icoil)%name,'("mod_"I3.3)') icoil
    
    SALLOCATE( coil(icoil)%xc, (0:NFcoil), zero )
    SALLOCATE( coil(icoil)%xs, (0:NFcoil), zero )
    SALLOCATE( coil(icoil)%yc, (0:NFcoil), zero )
    SALLOCATE( coil(icoil)%ys, (0:NFcoil), zero )
    SALLOCATE( coil(icoil)%zc, (0:NFcoil), zero )
    SALLOCATE( coil(icoil)%zs, (0:NFcoil), zero )
    
    coil(icoil)%gdof = -1
    coil(icoil)%idof = -1
    
    zeta = ( icoil - 1 + zero ) * pi2 / Ncoils
!   zeta = ( icoil - 1 + half ) * pi2 / Ncoils
    
    select case( Isurface )
     
    case( 0 ) ! VMEC-like boundary; 28 Aug 17;
     
     call surfcoord( zero, zeta, r1, z1)
     call surfcoord(   pi, zeta, r2, z2)
     
     Ro = half * ( r1 + r2 )
     Zo = half * ( z1 + z2 )        
     
     coil(icoil)%xc(0:1) = (/ Ro * cos(zeta), init_radius * cos(zeta) /)
!    coil(icoil)%xs(0:1) = (/          zero ,                   zero  /)
     coil(icoil)%yc(0:1) = (/ Ro * sin(zeta), init_radius * sin(zeta) /)
!    coil(icoil)%ys(0:1) = (/          zero ,                   zero  /)
     coil(icoil)%zc(0:1) = (/ Zo            ,                   zero  /)
     coil(icoil)%zs(0:1) = (/          zero , init_radius             /)
     
    case( 1 ) ! axis-style boundary; 28 Aug 17;
     
     FATAL( rdcoils, init_radius.lt.one, definition of init_radius has changed )

     do jj = 1, Ns ! Ns is input; 29 Sep 17;
      
      tt = (jj-1) * pi2 / Ns !poloidal angle;
      
      call knotxx( init_radius, tt, zeta, ax, at, az, xx, xs, xt, xz, xtt, xtz, xzz, v1, v2, w1, w2 )
      
      coilsX(jj,icoil) = xx(1)
      coilsY(jj,icoil) = xx(2)
      coilsZ(jj,icoil) = xx(3)
      
     enddo ! end of do jj ; 28 Aug 17;
     
     call fourier( coilsX(1:Ns,icoil), coil(icoil)%xc(0:NFcoil), coil(icoil)%xs(0:NFcoil), Ns, NFcoil )
     call fourier( coilsY(1:Ns,icoil), coil(icoil)%yc(0:NFcoil), coil(icoil)%ys(0:NFcoil), Ns, NFcoil )
     call fourier( coilsZ(1:Ns,icoil), coil(icoil)%zc(0:NFcoil), coil(icoil)%zs(0:NFcoil), Ns, NFcoil )
     
    end select ! end of select case( Isurface );
    
   enddo ! end of do icoil;
   
   select case( Isurface )
   case( 0 ) ! coilsX,Y,Z not required if Fourier harmonics of coils are read from file; 29 Sep 17;
   case( 1 )
    DALLOCATE( coilsX )
    DALLOCATE( coilsY )
    DALLOCATE( coilsZ )
   case default
    FATAL( rdcoils, .true., selected Isurface not supported )
   end select
   
  case default
   
   FATAL( rdcoils, .true., selected Initialize not supported )
   
  end select ! end of select case( Initialize);
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  do icoil = 1, Ncoils
   
   FATAL( rdcoils, coil(icoil)%Ifree < 0    .or. coil(icoil)%Ifree > 1, illegal )
   FATAL( rdcoils, coil(icoil)%Lfree < 0    .or. coil(icoil)%Lfree > 1, illegal )
  !FATAL( rdcoils, coil(icoil)%Lo    < zero                           , illegal )
   
  enddo
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  do icoil = 1, Ncoils
   
   SALLOCATE( coil(icoil)%xx, (0:Ns-1), zero )
   SALLOCATE( coil(icoil)%yy, (0:Ns-1), zero )
   SALLOCATE( coil(icoil)%zz, (0:Ns-1), zero )
   SALLOCATE( coil(icoil)%xt, (0:Ns-1), zero )
   SALLOCATE( coil(icoil)%yt, (0:Ns-1), zero )
   SALLOCATE( coil(icoil)%zt, (0:Ns-1), zero )
   SALLOCATE( coil(icoil)%xa, (0:Ns-1), zero )
   SALLOCATE( coil(icoil)%ya, (0:Ns-1), zero )
   SALLOCATE( coil(icoil)%za, (0:Ns-1), zero )

  enddo ! end of do icoil; 29 Sep 17;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! SALLOCATE( cosip, (0:Npc),  one )  ! cos( ip * pi / Np ) ; default  one;
! SALLOCATE( sinip, (0:Npc), zero )  ! sin( ip * pi / Np ) ; default zero;
  
! if ( Npc >= 2 ) then
   
!  do ip = 1, Npc-1
    
!   cosip(ip) = cos( ip * pi2 / Npc )
!   sinip(ip) = sin( ip * pi2 / Npc )
    
!   do icoil = 1, Ncoils
     
!    select case( coil(icoil)%itype )
      
!    case( 1 )
      
     !NF = FouCoil(icoil)%NF
      
!     SALLOCATE( coil(icoil+ip*Ncoils)%xc, (0:coil(icoil)%NF), zero )
!     SALLOCATE( coil(icoil+ip*Ncoils)%xs, (0:coil(icoil)%NF), zero )
!     SALLOCATE( coil(icoil+ip*Ncoils)%yc, (0:coil(icoil)%NF), zero )
!     SALLOCATE( coil(icoil+ip*Ncoils)%ys, (0:coil(icoil)%NF), zero )
!     SALLOCATE( coil(icoil+ip*Ncoils)%zc, (0:coil(icoil)%NF), zero )
!     SALLOCATE( coil(icoil+ip*Ncoils)%zs, (0:coil(icoil)%NF), zero )
      
!    case default
      
!     FATAL(discoil, .true., selected coil type not supported )
      
!    end select
     
!   enddo ! end of do icoil; 29 Sep 17;
    
!  enddo ! end of do ip; 29 Sep 17;
   
!  call mapcoil ! map perodic coils;
   
! endif ! end of if( Npc.ge.2 );

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! totalcurrent = zero
  
! do icoil = 1, Ncoils*Npc

!  FATAL( rdcoils, coil(icoil)%I.lt.zero, has negative current been considered )

!  totalcurrent = totalcurrent + coil(icoil)%I

! enddo
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! if(myid == 0 .and. IsQuiet <= 0) then
!    write(ounit,'("rdcoils : "i3" fixed currents ; "i3" fixed geometries.")') &
!         & Nfixcur, Nfixgeo
!    !write( ounit,'("rdcoils : total current G ="es23.15" ; 2 . pi2 . G = "es23.15" ;")') &
!    !     & totalcurrent, totalcurrent * pi2 * two
! endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!  itmp   = -1 ; call AllocData(itmp)
!
!  ifirst =  1 ; call discoil(ifirst) ; ifirst = 0
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  maxNF = 0
  
  do icoil = 1, Ncoils
   
   maxNF = max( coil(icoil)%NF, maxNF ) ! determine max. Fourier resolution;
   
  enddo
  
  SALLOCATE( cmt, (0:Ns-1, 0:maxNF), zero ) ! trigonometric factors mapping coils from `Fourier' to real space;
  SALLOCATE( smt, (0:Ns-1, 0:maxNF), zero )
  
  do kk = 0, Ns-1 ! 12 Nov 17;
   
   tt = kk * discretecurve
   
   do mm = 0, maxNF
    
    cmt(kk,mm) = cos( mm * tt )
    smt(kk,mm) = sin( mm * tt )
    
   enddo ! end of do mm;
   
  enddo ! end of do kk;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  return

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine rdcoils

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

SUBROUTINE readcoils(filename, maxnseg)
  use globals, only : zero, coilseg, coilname, Ncoils, ounit, myid
  implicit none
  include "mpif.h"

  INTEGER, parameter         :: mcoil = 256, mseg = 1024 ! Largest coils and segments number
  INTEGER                    :: icoil, cunit, istat, astat, lstat, ierr, maxnseg, seg(1:mseg)
  REAL, dimension(mseg,mcoil):: x, y, z, I
  CHARACTER*100              :: filename
  CHARACTER*200              :: line
  REAL                       :: tmp
  CHARACTER (LEN=20), dimension(mcoil) :: name
  REAL,    allocatable :: coilsX(:,:), coilsY(:,:), coilsZ(:,:), coilsI(:)

  cunit = 99; I = 1.0; Ncoils= 1; maxnseg = 0; seg = 0;

  open(cunit,FILE=trim(filename),STATUS='old',IOSTAT=istat)
  if ( istat .ne. 0 ) then
     write(ounit,'("rdcoils : Reading coils data error in "A)') trim(filename)
     call MPI_ABORT( MPI_COMM_WORLD, 1, ierr )
  endif
     
  ! read coils and segments data
  read(cunit,* )
  read(cunit,* )
  read(cunit,* )

  do
     read(cunit,'(a)', IOSTAT = istat) line
     if(istat .ne. 0 .or. line(1:3) == 'end') exit !detect EOF or end

     seg(Ncoils) = seg(Ncoils) + 1
     read(line,*, IOSTAT = lstat) x(seg(Ncoils), Ncoils), y(seg(Ncoils), Ncoils), z(seg(Ncoils), Ncoils), &
          I(seg(Ncoils), Ncoils)
     if ( I(seg(Ncoils), Ncoils) == 0.0 ) then
        seg(Ncoils) = seg(Ncoils) - 1  !remove the duplicated last point
        read(line, *, IOSTAT = lstat) tmp, tmp, tmp, tmp, tmp, name(Ncoils)
        Ncoils = Ncoils + 1
     endif
  enddo

  close(cunit)

  Ncoils = Ncoils - 1
  maxnseg = maxval(seg)
  FATAL( readcoils, Ncoils  .le. 0 .or. Ncoils   >  mcoil, illegal )
  FATAL( readcoils, maxnseg .le. 0 .or. maxnseg  >  mseg , illegal )

#ifdef DEBUG
  write(ounit,'("rdcoils : Finding " I4 " coils and maximum segments number is " I6)') &
       Ncoils, maxnseg
  do icoil = 1, Ncoils
     write(ounit, '("rdcoils : Number of segments in coil " I4 " is " I6)') icoil, seg(icoil)
  enddo
#endif

  SALLOCATE( coilsX  , (1:maxnseg, 1:Ncoils), zero )
  SALLOCATE( coilsY  , (1:maxnseg, 1:Ncoils), zero )
  SALLOCATE( coilsZ  , (1:maxnseg, 1:Ncoils), zero )
  SALLOCATE( coilsI  , (1:Ncoils)           , zero )
  SALLOCATE( coilseg , (1:Ncoils)           ,    0 )
  SALLOCATE( coilname, (1:Ncoils)           , ''   )

  coilsX( 1:maxnseg, 1:Ncoils) = x( 1:maxnseg, 1:Ncoils)
  coilsY( 1:maxnseg, 1:Ncoils) = y( 1:maxnseg, 1:Ncoils)
  coilsZ( 1:maxnseg, 1:Ncoils) = z( 1:maxnseg, 1:Ncoils)
  coilsI(            1:Ncoils) = I( 1        , 1:Ncoils)
  coilseg(           1:Ncoils) = seg(          1:Ncoils)
  coilname(          1:Ncoils) = name(         1:Ncoils)

  return

end SUBROUTINE READCOILS

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

SUBROUTINE fourier( X, XFC, XFS, Nsegs, NFcoil)
  use globals, only: ounit, zero, pi2, half, myid
  implicit none
  include "mpif.h"

  REAL    :: X(1:Nsegs), XFC(0:NFcoil), XFS(0:NFcoil)
  INTEGER :: Nsegs, NFcoil, ifou, iseg, funit, ierr
  REAL, allocatable:: A(:), B(:)

  allocate(A(0:Nsegs-1))
  allocate(B(0:Nsegs-1))

  FATAL(fourier, Nsegs < 2*NFcoil, Nsegs too small)
  A = zero; B = zero

  do ifou = 0, Nsegs-1

     do iseg = 1, Nsegs
        A(ifou) = A(ifou) + X(iseg)*cos(ifou*pi2*(iseg-1)/Nsegs)
        B(ifou) = B(ifou) + X(iseg)*sin(ifou*pi2*(iseg-1)/Nsegs)
     enddo

  enddo

  A = 2.0/Nsegs * A; A(0) = half*A(0)
  B = 2.0/Nsegs * B

  XFC(0:NFcoil) = A(0:NFcoil)
  XFS(0:NFcoil) = B(0:NFcoil)

  deallocate(A, B)

  return

end SUBROUTINE fourier

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine surfcoord( theta, zeta, r, z)
  use globals, only: zero, Nfou, bim, bin, Rbc, Rbs, Zbc, Zbs
  implicit none
  include "mpif.h"

  REAL, INTENT(in ) :: theta, zeta
  REAL, INTENT(out) :: r, z

  INTEGER           :: imn
  REAL              :: arg
  !-------------calculate r, z coodinates for theta, zeta------------------------------------------------  
  if( .not. allocated(bim) ) STOP  "please allocate surface data first!"

  r = zero; z = zero  
  do imn = 1, Nfou
     arg = bim(imn) * theta - bin(imn) * zeta
     R =  R + Rbc(imn) * cos(arg) + Rbs(imn) * sin(arg)
     Z =  Z + Zbc(imn) * cos(arg) + Zbs(imn) * sin(arg)
  enddo

  return
end subroutine surfcoord

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

