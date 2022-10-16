
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
!latex  \item[1.] \inputvar{case\_init = 1} : Toroidally placing \inputvar{Ncoils} circular coils with a 
!latex             radius of \inputvar{init\_radius} and current of \inputvar{init\_current}. The $i$th coil 
!latex             is placed at $\z = \frac{i-1}{Ncoils} \frac{2\pi}{Nfp}$.
!latex  \item[2.] \inputvar{case\_init = 0} : Read coils data from {\bf ext.focus} file. The format is as following. 
!latex     \red{This is the most flexible way, and  each coil can be different.}            
!latex  \begin{raw}
!latex   # Total number of coils
!latex              16
!latex   #------------1--------------------------------
!latex   #coil_type  symm   coil_name
!latex       1   0  Mod_001   
!latex   #Nseg        current         Ifree         Length         Lfree  target_length
!latex     128  9.844910899889484E+05     1  5.889288927667147E+00     1  1.000000000000000E+00
!latex   #NFcoil
!latex    4
!latex   #Fourier harmonics for coils ( xc; xs; yc; ys; zc; zs) 
!latex    3.044612087666170E+00  8.531153655332238E-01  4.194525679767678E-02  2.139790853335835E-02  3.243811555342430E-03
!latex    0.000000000000000E+00  3.542408058492299E-16 -9.108712738922674E-16  1.841880477639364E-16 -1.172175996642087E-16
!latex   -4.456021385977147E-15  8.545613874434043E-16 -3.133154295448265E-16  1.764367073160815E-16 -1.187904023667544E-16
!latex    0.000000000000000E+00 -5.425716121023922E-02 -8.986316303345250E-02 -2.946386365076052E-03 -4.487052148209031E-03
!latex   -4.293247278325474E-17 -1.303273952226587E-15  7.710821807870230E-16 -3.156539892466338E-16  9.395672288215928E-17
!latex    0.000000000000000E+00  9.997301975562740E-01  2.929938238054118E-02  2.436889176706748E-02  1.013941937492003E-03
!latex   #-----------2--permanent magnet---------------
!latex   #coil_type  symm   coil_name
!latex       2  0  dipole_01  
!latex   #  Lc  ox   oy   oz  Ic  I  mt  mp
!latex      1   0.0  0.0  0.0  1 1.0E6  0.0  0.0
!latex   #-----------3--backgound Bt Bz----------------
!latex   #coil_type  symm   coil_name
!latex       3  0  bg_BtBz_01
!latex   # Ic     I    Lc  Bz  (Ic control I; Lc control Bz)
!latex     1    1.0E6  0  0.0
!latex      .
!latex      .
!latex      .
!latex#------------coil_type_spline--------------------------------
!latex   #coil_type  symm   coil_name
!latex       5   0  Mod_001   
!latex   #Nseg        current         Ifree         Length         Lfree  target_length
!latex     128  9.844910899889484E+05     1  5.889288927667147E+00     1  1.000000000000000E+00
!latex   #NCP 
!latex    13    
!latex   #Vector of knots
!latex   -0.3 -0.2 -0.1 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3
!latex   #Control Points Coordinates for coils ( x; y; z)  
!latex    0.58   0.4354 0.5464 0.242  0.456546 0.1414 0.465 0.23455 0.3636 0.75746 0.58 0.4354 0.5464
!latex    0.4353 1.3453 2.356  0.3453 0.3453   0.23444 2.454 5.2453 2.5654 0.53453 0.4353 1.3453 2.356 
!latex    -1.243 2.242  0.234  1.25   2.245    0.68876 1.45  4.4543 1.3445 6.3545  -1.243 2.242  0.234
!latex  
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
  use globals
  use mpi
  implicit none

  LOGICAL   :: exist
  INTEGER   :: icoil, maxnseg, ifirst, NF, itmp, ip, icoef, total_coef, num_pm, num_bg, & 
               num_per_array, num_tor, ipol, itor
  REAL      :: Rmaj, zeta, totalcurrent, z0, r1, r2, z1, z2, rtmp, teta
  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  Nfixcur = 0 ! fixed coil current number
  Nfixgeo = 0 ! fixed coil geometry number
  num_pm  = 0 ! number of permanent magnets
  num_bg  = 0 ! number of background field

  if(myid == 0) write(ounit, *) "-----------INITIALIZE COILS----------------------------------"
  select case( case_init )
  !-------------read coils file--------------------------------------------------------------------------
  case(-1 )
     if (myid == 0) then
        write(ounit,'("rdcoils : Reading coils data (MAKEGRID format) from "A)') trim(input_coils)
        call readcoils(input_coils, maxnseg)
        write(ounit,'("        : Read ",i6," coils.")') Ncoils
        if (IsQuiet < 0) write(ounit, '(8X,": NFcoil = "I3" ; IsVaryCurrent = "I1 &
             " ; IsVaryGeometry = "I1)') NFcoil, IsVaryCurrent, IsVaryGeometry
     endif
     IlBCAST( Ncoils   ,      1, 0 )
     IlBCAST( maxnseg  ,      1, 0 )
     ! allocate arrays on other nodes;
     if( .not. allocated(coilsX) ) then 
        SALLOCATE( coilsX, (1:maxnseg,1:Ncoils), zero )    
        SALLOCATE( coilsY, (1:maxnseg,1:Ncoils), zero )
        SALLOCATE( coilsZ, (1:maxnseg,1:Ncoils), zero )
        SALLOCATE( coilsI, (          1:Ncoils), zero )
        SALLOCATE( coilseg,(          1:Ncoils),    0 )
        SALLOCATE( coilname,(         1:Ncoils),   '' )
     endif
     ! broadcast coils data;
     RlBCAST( coilsX, maxnseg*Ncoils, 0 )
     RlBCAST( coilsY, maxnseg*Ncoils, 0 )
     RlBCAST( coilsZ, maxnseg*Ncoils, 0 )
     RlBCAST( coilsI,         Ncoils, 0 )   
     IlBCAST( coilseg,        Ncoils, 0 )
     ClBCAST( coilname,       Ncoils, 0 )
     ! Ncoils are the number of unique coils
     allocate(    coil(1:Ncoils) )
     allocate( FouCoil(1:Ncoils) )
     allocate(     DoF(1:Ncoils) )
     !Ncoils = Ncoils / Npc ! Ncoils changed to unique number of coils;
     icoil = 0
     do icoil = 1, Ncoils
        ! general coil parameters;
        coil(icoil)%symm = 0 ! no symmetry or periodicity
        coil(icoil)%NS =  Nseg  
        coil(icoil)%I  =  coilsI(icoil)
        coil(icoil)%Ic =  IsVaryCurrent
        coil(icoil)%L  =  target_length ! irrelevant until re-computed;
        !coil(icoil)%curv_k0 =  curv_k0
        coil(icoil)%Lc =  IsVaryGeometry
        coil(icoil)%Lo =  target_length
        coil(icoil)%name = trim(coilname(icoil))
        ! check coil current and length
        FATAL( rdcoils01, coil(icoil)%Ic < 0 .or. coil(icoil)%Ic > 1, illegal )
        FATAL( rdcoils02, coil(icoil)%Lc < 0 .or. coil(icoil)%Lc > 1, illegal )
        FATAL( rdcoils03, coil(icoil)%Lo < zero                     , illegal )
        if(coil(icoil)%Ic == 0) Nfixcur = Nfixcur + 1
        if(coil(icoil)%Lc == 0) Nfixgeo = Nfixgeo + 1
        ! Fourier representation related;
        FouCoil(icoil)%NF = NFcoil
        NF = NFcoil  ! alias
        SALLOCATE( FouCoil(icoil)%xc, (0:NF), zero )
        SALLOCATE( FouCoil(icoil)%xs, (0:NF), zero )
        SALLOCATE( FouCoil(icoil)%yc, (0:NF), zero )
        SALLOCATE( FouCoil(icoil)%ys, (0:NF), zero )
        SALLOCATE( FouCoil(icoil)%zc, (0:NF), zero )
        SALLOCATE( FouCoil(icoil)%zs, (0:NF), zero )        
        !if(myid .ne. modulo(icoil-1, ncpu)) cycle
        ! Fourier transformation (FFT might be appied)
        call Fourier( coilsX(1:coilseg(icoil),icoil), Foucoil(icoil)%xc, Foucoil(icoil)%xs, coilseg(icoil), NF)
        call Fourier( coilsY(1:coilseg(icoil),icoil), Foucoil(icoil)%yc, Foucoil(icoil)%ys, coilseg(icoil), NF)
        call Fourier( coilsZ(1:coilseg(icoil),icoil), Foucoil(icoil)%zc, Foucoil(icoil)%zs, coilseg(icoil), NF)
     enddo
     ! clean space
     DALLOCATE( coilsX )
     DALLOCATE( coilsY )
     DALLOCATE( coilsZ )
     DALLOCATE( coilsI )
     DALLOCATE( coilseg)
     DALLOCATE(coilname)
     ! use Fourier representation by default
     coil(1:Ncoils)%type = 1

  !-------------individual coil file---------------------------------------------------------------------
  case( 0 )
     ! get coil number first
     if( myid==0 ) then  
        open( runit, file=trim(input_coils), status="old", action='read')
        read( runit,*)
        read( runit,*) Ncoils
        write(ounit,'("rdcoils : identified "i6" unique coils in "A" ;")') Ncoils, trim(input_coils)
     endif                     
     ! broadcase and allocate data
     IlBCAST( Ncoils        ,    1,  0 )
     allocate( FouCoil(1:Ncoils) )
     allocate(    coil(1:Ncoils) )
     allocate(     DoF(1:Ncoils) )
     allocate(  Splines(1:Ncoils) )
     ! master CPU read the coils
     if( myid==0 ) then
        do icoil = 1, Ncoils
           read( runit,*)
           read( runit,*)
           read( runit,*) coil(icoil)%type, coil(icoil)%symm, coil(icoil)%name
           FATAL( rdcoils04, (coil(icoil)%type < 1 .or. coil(icoil)%type > 3) .and. &
               (coil(icoil)%type .ne. coil_type_spline) .and. (coil(icoil)%type .ne. coil_type_multi) , illegal )
           FATAL( rdcoils05, coil(icoil)%symm < 0 .or. coil(icoil)%symm > 2, illegal )
           if(coil(icoil)%type == 1 .or. coil(icoil)%type == coil_type_multi) then  ! Fourier representation (and temporarily also for multi)
              read( runit,*)
              read( runit,*) coil(icoil)%NS, coil(icoil)%I, coil(icoil)%Ic, &
                   & coil(icoil)%L, coil(icoil)%Lc, coil(icoil)%Lo
              FATAL( rdcoils06, coil(icoil)%NS < 0                        , illegal )
              FATAL( rdcoils07, coil(icoil)%Ic < 0 .or. coil(icoil)%Ic > 1, illegal )
              FATAL( rdcoils08, coil(icoil)%Lc < 0 .or. coil(icoil)%Lc > 1, illegal )
              FATAL( rdcoils09, coil(icoil)%L  < zero                     , illegal )
              FATAL( rdcoils10, coil(icoil)%Lo < zero                     , illegal )
              read( runit,*)
              read( runit,*) FouCoil(icoil)%NF
              FATAL( rdcoils12, Foucoil(icoil)%NF  < 0                    , illegal )
              SALLOCATE( FouCoil(icoil)%xc, (0:FouCoil(icoil)%NF), zero )
              SALLOCATE( FouCoil(icoil)%xs, (0:FouCoil(icoil)%NF), zero )
              SALLOCATE( FouCoil(icoil)%yc, (0:FouCoil(icoil)%NF), zero )
              SALLOCATE( FouCoil(icoil)%ys, (0:FouCoil(icoil)%NF), zero )
              SALLOCATE( FouCoil(icoil)%zc, (0:FouCoil(icoil)%NF), zero )
              SALLOCATE( FouCoil(icoil)%zs, (0:FouCoil(icoil)%NF), zero )
              read( runit,*)
              read( runit,*) FouCoil(icoil)%xc(0:FouCoil(icoil)%NF)
              read( runit,*) FouCoil(icoil)%xs(0:FouCoil(icoil)%NF)
              read( runit,*) FouCoil(icoil)%yc(0:FouCoil(icoil)%NF)
              read( runit,*) FouCoil(icoil)%ys(0:FouCoil(icoil)%NF)
              read( runit,*) FouCoil(icoil)%zc(0:FouCoil(icoil)%NF)
              read( runit,*) FouCoil(icoil)%zs(0:FouCoil(icoil)%NF)
           else if (coil(icoil)%type == 2) then  ! permanent magnets
              read( runit,*)
              read( runit,*) coil(icoil)%Lc, coil(icoil)%ox, coil(icoil)%oy, coil(icoil)%oz, &
                             coil(icoil)%Ic, coil(icoil)%I , coil(icoil)%mt, coil(icoil)%mp      
           else if (coil(icoil)%type == 3) then  ! backgroud toroidal/vertical field
              read( runit,*)
              read( runit,*) coil(icoil)%Ic, coil(icoil)%I, coil(icoil)%Lc, coil(icoil)%Bz 
              coil(icoil)%symm = 0 ! automatic reset to 0; might not be necessary; 2020/01/17
           else if(coil(icoil)%type == coil_type_spline) then  ! Spline representation
              read( runit,*)
              read( runit,*) coil(icoil)%NS, coil(icoil)%I, coil(icoil)%Ic, &
                   & coil(icoil)%L, coil(icoil)%Lc, coil(icoil)%Lo
              FATAL( rdcoils06, coil(icoil)%NS < 0                        , illegal )
              FATAL( rdcoils07, coil(icoil)%Ic < 0 .or. coil(icoil)%Ic > 1, illegal )
              FATAL( rdcoils08, coil(icoil)%Lc < 0 .or. coil(icoil)%Lc > 1, illegal )
              FATAL( rdcoils09, coil(icoil)%L  < zero                     , illegal )
              FATAL( rdcoils10, coil(icoil)%Lo < zero                     , illegal )
              read( runit,*)
              read( runit,*) Splines(icoil)%NCP
              Splines(icoil)%NT = Splines(icoil)%NCP + 4
              FATAL( rdcoils12, Splines(icoil)%NCP  < 0                    , illegal )
              FATAL( rdcoils12_2, Splines(icoil)%NT  < 0                     , illegal )
              FATAL( rdcoils12_3, Splines(icoil)%NT .NE. Splines(icoil)%NCP +4, illegal )
              FATAL( rdcoils12_4, Splines(icoil)%NT < 10, illegal )
              SALLOCATE( Splines(icoil)%vect, (0:Splines(icoil)%NT-1), zero )
              SALLOCATE( Splines(icoil)%eval_points, (0:coil(icoil)%NS-1), zero )
              SALLOCATE( Splines(icoil)%Cpoints, (0:Splines(icoil)%NCP * 3 - 1 ), zero )
              read( runit,*)
              read( runit,*) Splines(icoil)%vect(0:Splines(icoil)%NT-1)
              Splines(icoil)%vect(Splines(icoil)%NT-3) = 1.0 + Splines(icoil)%vect(4) - Splines(icoil)%vect(3)
              Splines(icoil)%vect(Splines(icoil)%NT-2) = 1.0 + Splines(icoil)%vect(5) - Splines(icoil)%vect(3)
              Splines(icoil)%vect(Splines(icoil)%NT-1) = 1.0 + Splines(icoil)%vect(6) - Splines(icoil)%vect(3)
              Splines(icoil)%vect(0) =  Splines(icoil)%vect(Splines(icoil)%NT-7) - Splines(icoil)%vect(Splines(icoil)%NT-4)
              Splines(icoil)%vect(1) =  Splines(icoil)%vect(Splines(icoil)%NT-6) - Splines(icoil)%vect(Splines(icoil)%NT-4)
              Splines(icoil)%vect(2) =  Splines(icoil)%vect(Splines(icoil)%NT-5) - Splines(icoil)%vect(Splines(icoil)%NT-4)
              read( runit,*)
              read( runit,*) Splines(icoil)%Cpoints(0:Splines(icoil)%NCP-1)
              read( runit,*) Splines(icoil)%Cpoints(Splines(icoil)%NCP:Splines(icoil)%NCP*2-1)
              read( runit,*) Splines(icoil)%Cpoints(Splines(icoil)%NCP*2:Splines(icoil)%NCP*3-1)
           else
              STOP " wrong coil type in rdcoils"
              call MPI_ABORT(MPI_COMM_FOCUS, 1, ierr)
           endif              
        enddo !end do icoil;
        close( runit )
     endif ! end of if( myid==0 );
     ! broad cast the data and allocate space on slavers
     do icoil = 1, Ncoils
        IlBCAST( coil(icoil)%type        , 1        ,  0 )
        IlBCAST( coil(icoil)%symm        , 1        ,  0 )
        ClBCAST( coil(icoil)%name        , 10       ,  0 )
        if(coil(icoil)%type == 1 .or. coil(icoil)%type == coil_type_multi) then  ! Fourier representation
           IlBCAST( coil(icoil)%NS           , 1        ,  0 )
           RlBCAST( coil(icoil)%I            , 1        ,  0 )
           IlBCAST( coil(icoil)%Ic           , 1        ,  0 )
           RlBCAST( coil(icoil)%L            , 1        ,  0 )
           IlBCAST( coil(icoil)%Lc           , 1        ,  0 )
           RlBCAST( coil(icoil)%Lo           , 1        ,  0 )
           IlBCAST( FouCoil(icoil)%NF        , 1        ,  0 )
           if (.not. allocated(FouCoil(icoil)%xc) ) then
              SALLOCATE( FouCoil(icoil)%xc, (0:FouCoil(icoil)%NF), zero )
              SALLOCATE( FouCoil(icoil)%xs, (0:FouCoil(icoil)%NF), zero )
              SALLOCATE( FouCoil(icoil)%yc, (0:FouCoil(icoil)%NF), zero )
              SALLOCATE( FouCoil(icoil)%ys, (0:FouCoil(icoil)%NF), zero )
              SALLOCATE( FouCoil(icoil)%zc, (0:FouCoil(icoil)%NF), zero )
              SALLOCATE( FouCoil(icoil)%zs, (0:FouCoil(icoil)%NF), zero ) 
           endif
           RlBCAST( FouCoil(icoil)%xc(0:FouCoil(icoil)%NF) , 1+FouCoil(icoil)%NF ,  0 )
           RlBCAST( FouCoil(icoil)%xs(0:FouCoil(icoil)%NF) , 1+FouCoil(icoil)%NF ,  0 )
           RlBCAST( FouCoil(icoil)%yc(0:FouCoil(icoil)%NF) , 1+FouCoil(icoil)%NF ,  0 )
           RlBCAST( FouCoil(icoil)%ys(0:FouCoil(icoil)%NF) , 1+FouCoil(icoil)%NF ,  0 )
           RlBCAST( FouCoil(icoil)%zc(0:FouCoil(icoil)%NF) , 1+FouCoil(icoil)%NF ,  0 )
           RlBCAST( FouCoil(icoil)%zs(0:FouCoil(icoil)%NF) , 1+FouCoil(icoil)%NF ,  0 )
           if(coil(icoil)%Ic == 0) Nfixcur = Nfixcur + 1
           if(coil(icoil)%Lc == 0) Nfixgeo = Nfixgeo + 1
        else if (coil(icoil)%type == 2) then  ! permanent magnets
           IlBCAST( coil(icoil)%Ic, 1 , 0 )
           RlBCAST( coil(icoil)%I , 1 , 0 )
           IlBCAST( coil(icoil)%Lc, 1 , 0 )
           RlBCAST( coil(icoil)%ox, 1 , 0 )
           RlBCAST( coil(icoil)%oy, 1 , 0 )
           RlBCAST( coil(icoil)%oz, 1 , 0 )
           RlBCAST( coil(icoil)%mt, 1 , 0 )
           RlBCAST( coil(icoil)%mp, 1 , 0 )
           if(coil(icoil)%Ic == 0) Nfixcur = Nfixcur + 1
           ! if(coil(icoil)%Lc == 0) Nfixgeo = Nfixgeo + 1
           Nfixgeo = Nfixgeo + 1 ! always treat as a fixed geometry
        else if (coil(icoil)%type == 3) then  ! backgroud toroidal/vertical field
           IlBCAST( coil(icoil)%Ic, 1 , 0 )
           RlBCAST( coil(icoil)%I , 1 , 0 )
           IlBCAST( coil(icoil)%Lc, 1 , 0 )
           RlBCAST( coil(icoil)%Bz, 1 , 0 )             
           if(coil(icoil)%Ic == 0) Nfixcur = Nfixcur + 1
           ! if(coil(icoil)%Lc == 0) Nfixgeo = Nfixgeo + 1
           Nfixgeo = Nfixgeo + 1 ! always treat as a fixed geometry
        else if (coil(icoil)%type == coil_type_spline) then
           IlBCAST( coil(icoil)%NS           , 1        ,  0 )
           RlBCAST( coil(icoil)%I            , 1        ,  0 )
           IlBCAST( coil(icoil)%Ic           , 1        ,  0 )
           RlBCAST( coil(icoil)%L            , 1        ,  0 )
           IlBCAST( coil(icoil)%Lc           , 1        ,  0 )
           RlBCAST( coil(icoil)%Lo           , 1        ,  0 )
           IlBCAST( Splines(icoil)%NCP        , 1        ,  0 )
           IlBCAST( Splines(icoil)%NT         , 1        ,  0 )
           if (.not. allocated(Splines(icoil)%vect) ) then
              SALLOCATE( Splines(icoil)%vect, (0:Splines(icoil)%NT-1), zero )
              SALLOCATE( Splines(icoil)%eval_points, (0:coil(icoil)%NS-1), zero )
              SALLOCATE( Splines(icoil)%Cpoints, (0:Splines(icoil)%NCP * 3 -1 ), zero )
           endif
           RlBCAST( Splines(icoil)%vect(0:Splines(icoil)%NT -1) , Splines(icoil)%NT ,  0 )
           RlBCAST( Splines(icoil)%eval_points(0:coil(icoil)%NS-1) , coil(icoil)%NS ,  0 )
           RlBCAST( Splines(icoil)%Cpoints(0:Splines(icoil)%NCP*3 - 1) , Splines(icoil)%NCP*3 ,  0 )
           if(coil(icoil)%Ic == 0) Nfixcur = Nfixcur + 1
           if(coil(icoil)%Lc == 0) Nfixgeo = Nfixgeo + 1
        else
           STOP " wrong coil type in rdcoils"
           call MPI_ABORT(MPI_COMM_FOCUS, 1, ierr)
        endif
     enddo
  !-------------toroidally placed circular coils---------------------------------------------------------
  case( 1 ) ! toroidally placed coils; 2017/03/13
     ! allocate data
     allocate( FouCoil(1:Ncoils) )
     allocate(    coil(1:Ncoils) )
     allocate(     DoF(1:Ncoils) )
     ! screen outputs
     if (myid == 0)  then
        write(ounit, '(8X,": Initialize "I4" unique circular coils with r="ES12.5"m ; I="&
             ES12.5" A")') Ncoils, init_radius, init_current
        if (IsQuiet < 0) write(ounit, '(8X,": NFcoil = "I3" ; IsVaryCurrent = "I1 &
             " ; IsVaryGeometry = "I1)') NFcoil, IsVaryCurrent, IsVaryGeometry
     endif
     ! initializations
     do icoil = 1, Ncoils
        ! general coil parameters;
        coil(icoil)%type = 1
        coil(icoil)%symm = IsSymmetric ! follow the general setting
        coil(icoil)%NS =  Nseg  
        coil(icoil)%I  =  init_current
        coil(icoil)%Ic =  IsVaryCurrent
        coil(icoil)%L  =  pi2*init_radius
        coil(icoil)%Lc =  IsVaryGeometry
        coil(icoil)%Lo =  target_length
        write(coil(icoil)%name,'("Mod_"I3.3)') icoil
        FATAL( rdcoils, coil(icoil)%Ic < 0 .or. coil(icoil)%Ic > 1, illegal )
        FATAL( rdcoils, coil(icoil)%Lc < 0 .or. coil(icoil)%Lc > 1, illegal )
        FATAL( rdcoils, coil(icoil)%Lo < zero                     , illegal )
        if(coil(icoil)%Ic == 0) Nfixcur = Nfixcur + 1
        if(coil(icoil)%Lc == 0) Nfixgeo = Nfixgeo + 1
        ! Fourier representation related;
        FouCoil(icoil)%NF = NFcoil
        SALLOCATE( FouCoil(icoil)%xc, (0:NFcoil), zero )
        SALLOCATE( FouCoil(icoil)%xs, (0:NFcoil), zero )
        SALLOCATE( FouCoil(icoil)%yc, (0:NFcoil), zero )
        SALLOCATE( FouCoil(icoil)%ys, (0:NFcoil), zero )
        SALLOCATE( FouCoil(icoil)%zc, (0:NFcoil), zero )
        SALLOCATE( FouCoil(icoil)%zs, (0:NFcoil), zero )
        ! get the geometry center
        zeta = (icoil-1+half) * pi2 / (Ncoils*surf_Nfp*2**symmetry)  ! put a half for a shift;       
        call surfcoord( plasma, zero, zeta, r1, z1)
        call surfcoord( plasma,   pi, zeta, r2, z2)
        Rmaj = half * (r1 + r2)
        z0   = half * (z1 + z2)              
        ! initilize with circular coils;
        FouCoil(icoil)%xc(0:1) = (/ Rmaj * cos(zeta), init_radius * cos(zeta) /)
        FouCoil(icoil)%xs(0:1) = (/ zero            , zero                    /)
        FouCoil(icoil)%yc(0:1) = (/ Rmaj * sin(zeta), init_radius * sin(zeta) /)
        FouCoil(icoil)%ys(0:1) = (/ zero            , zero                    /)
        FouCoil(icoil)%zc(0:1) = (/ z0              , zero                    /)
        Foucoil(icoil)%zs(0:1) = (/ zero            , init_radius             /)
     enddo ! end of do icoil;
  !------------- permanent dipoles and background magnetic field ----------------------------------------
  case( 2 ) ! averagely positioned permanent dipoles ; will be removed;  2020/01/17
     allocate( coil(1:Ncoils) )
     allocate(  DoF(1:Ncoils) )
     num_per_array = 16  ! number of dipoles at each toroidal cross-section
     num_tor = (Ncoils-1)/num_per_array ! number of toroidal arrangements
     if (myid == 0)  then
        write(ounit,'("rdcoils : initializing "i3" uniformly positioned magnetic dipoles with toroidal magnetif filed")') Ncoils-1
        if (IsQuiet < 1) write(ounit, '(8X,": Initialize "I4" X "I4" dipoles on r="ES12.5"m  with m="&
             ES12.5" A")') num_tor, num_per_array, init_radius, init_current
        if (IsQuiet < 0) write(ounit, '(8X,": IsVaryCurrent = "I1 " ; IsVaryGeometry = "I1)') &
             IsVaryCurrent, IsVaryGeometry
        FATAL( rdcoils, modulo(Ncoils-1, num_per_array) /= 0, Please provide a valid number )
     endif
     ! background magnetic field Bt Bz
     icoil = 1
     coil(icoil)%I  =  init_current
     coil(icoil)%Ic =  IsVaryCurrent
     coil(icoil)%L  =  pi2*init_radius
     coil(icoil)%Lc =  0               ! IsVaryGeometry ! ignore Bz first; 20190102
     coil(icoil)%Lo =  target_length
     !coil(icoil)%curv_k0 =  curv_k0
     coil(icoil)%Bz =  zero
     coil(icoil)%name = 'bg_BtBz_01'
     coil(icoil)%type = 3

     do itor = 1, num_tor
        zeta = (itor-1) * pi2 / num_tor  ! put a half for a shift;
        call surfcoord( plasma, zero, zeta, r1, z1)
        call surfcoord( plasma,   pi, zeta, r2, z2)
        Rmaj = half * (r1 + r2)
        z0   = half * (z1 + z2)     
        do ipol = 1, num_per_array
           icoil = icoil + 1
           !general coil parameters;
           coil(icoil)%type = 2
           coil(icoil)%Ic =  IsVaryCurrent
           coil(icoil)%I  =  init_current
           coil(icoil)%L  =  pi2*init_radius
           coil(icoil)%Lc =  IsVaryGeometry
           coil(icoil)%Lo =  target_length
           !coil(icoil)%curv_k0 =  curv_k0
           write(coil(icoil)%name,'("pm_"I6)') icoil
           FATAL( rdcoils, coil(icoil)%Ic < 0 .or. coil(icoil)%Ic > 1, illegal )
           FATAL( rdcoils, coil(icoil)%Lc < 0 .or. coil(icoil)%Lc > 1, illegal )
           FATAL( rdcoils, coil(icoil)%Lo < zero                     , illegal )
           if(coil(icoil)%Ic == 0) Nfixcur = Nfixcur + 1
           if(coil(icoil)%Lc == 0) Nfixgeo = Nfixgeo + 1

           teta = (ipol-1) * pi2 / num_per_array
           rtmp = Rmaj + init_radius * cos(teta)
           coil(icoil)%ox = rtmp * cos(zeta)
           coil(icoil)%oy = rtmp * sin(zeta)
           coil(icoil)%oz = z0 + init_radius * sin(teta)

!!$           ! toroidal direction
!!$           coil(icoil)%mx = - init_current * sin(zeta)
!!$           coil(icoil)%my =   init_current * cos(zeta)
!!$           coil(icoil)%mz = zero
!!$
!!$           ! poloidal direction
!!$           coil(icoil)%mx = - init_current * sin(teta) * cos(zeta)
!!$           coil(icoil)%my = - init_current * sin(teta) * sin(zeta)
!!$           coil(icoil)%mz =   init_current * cos(teta)

           ! poloidal and toroidal angle; in poloidal direction
           coil(icoil)%mt = -teta
           coil(icoil)%mp =  zeta
!!$
!!$           ! inward direction
!!$           coil(icoil)%mt =  teta + half * pi
!!$           coil(icoil)%mp =  zeta + pi
!!$
!!$           ! toroidal direction
!!$           coil(icoil)%mt = half * pi
!!$           coil(icoil)%mp = zeta + half * pi

        enddo ! enddo ipol
     enddo ! enddo itor
     FATAL( rdcoils, icoil .ne. Ncoils, counting coils wrong when initializing )
  end select
  FATAL( rdcoils, Nfixcur > Ncoils, error with fixed currents )
  FATAL( rdcoils, Nfixgeo > Ncoils, error with fixed geometry )

  !-----------------------normalize currents and geometries-------------------------------------
  ! sum the total currents;
  totalcurrent = zero
  do icoil = 1, Ncoils
     totalcurrent = totalcurrent + coil(icoil)%I

     ! Need to handle currents for multi-filaments

  enddo
  if(myid == 0 .and. IsQuiet <= 0) then
     write(ounit,'("        : "i3" fixed currents ; "i3" fixed geometries.")') &
          & Nfixcur, Nfixgeo
     !write( ounit,'("        : total current G ="es23.15" ; 2 . pi2 . G = "es23.15" ;")') &
     !     & totalcurrent, totalcurrent * pi2 * two
  endif

  !-----------------------allocate DoF arrays --------------------------------------------------  
  itmp = -1
  call AllocData(itmp)

  !-----------------------discretize coil data--------------------------------------------------
  if (myid == 0) then
     if (IsQuiet < 0) write(ounit, '(8X,": coils will be discretized in "I6" segments")') Nseg
  endif

  ifirst = 1
  call discoil(ifirst)
  ifirst = 0

  return

end subroutine rdcoils

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine discoil(ifirst)
!---------------------------------------------------------------------------------------------
! dicretize coils data;
! if ifirst = 1, it will update all the coils; otherwise, only update free coils;
! date: 20170314
!---------------------------------------------------------------------------------------------
  !use globals, only: dp, zero, pi2, myid, ounit, coil, FouCoil, Ncoils, DoF, coil_type_multi, &
  !                   Splines, coil_type_spline, MPI_COMM_FOCUS
  use globals
  use mpi
  implicit none

  INTEGER, intent(in) :: ifirst
  
  INTEGER             :: icoil, iseg, mm, NS, NCP, NF, ND, ip, i, j, k
  
  REAL                :: tt, C, leng, htrans, freq, xc, yc, zc
  REAL, allocatable   :: zeta(:), sinterms(:), costerms(:), rc(:,:), rcp(:,:), rcpp(:,:), rcppp(:,:), rcdof(:,:,:), rcpdof(:,:,:), rcppdof(:,:,:), &
                         xcdof(:), ycdof(:), zcdof(:), delta(:,:), deltadof(:,:,:), &
                         n(:,:), ndoff(:,:,:), np(:,:), npp(:,:), npdof(:,:,:), nhatp(:,:), nhatpp(:,:), normt(:), normtp(:), normtpp(:), normn(:), normtdof(:,:), normtpdof(:,:), normndof(:,:), &
                         that(:,:), nhat(:,:), bhat(:,:), thatdof(:,:,:), nhatdof(:,:,:), nhatpdof(:,:,:), bhatdof(:,:,:), bhatp(:,:), bhatpp(:,:), &
                         thatp(:,:), thatpp(:,:), nhata(:,:), bhata(:,:), nhatap(:,:), bhatap(:,:), nhatapp(:,:), bhatapp(:,:), nhatadof(:,:,:), bhatadof(:,:,:), &
                         alphac(:), alphas(:), alphacdof(:,:), alphasdof(:,:), Cdoff(:), g(:), w(:), h(:), gp(:), wp(:), hp(:), gpp(:), wpp(:), hpp(:), &
                         gdof(:,:), wdof(:,:), hdof(:,:), lengdof(:), &
                         nhatpfd(:,:), nhatpdoffd(:,:,:)
  !-------------------------------------------------------------------------------------------
  do icoil = 1, Ncoils
     ! first time or if Lc/=0, then need discretize;
     if( (coil(icoil)%Lc + ifirst) /= 0) then  
        !if( myid.ne.modulo(icoil-1,ncpu) ) cycle ! parallelization loop;
        select case (coil(icoil)%type)
        case( 1 )
           ! reset to zero for all the coils;
           coil(icoil)%xx = zero
           coil(icoil)%yy = zero
           coil(icoil)%zz = zero
           coil(icoil)%xt = zero
           coil(icoil)%yt = zero
           coil(icoil)%zt = zero
           coil(icoil)%xa = zero
           coil(icoil)%ya = zero
           coil(icoil)%za = zero
           coil(icoil)%xb = zero
           coil(icoil)%yb = zero
           coil(icoil)%zb = zero
           NS = coil(icoil)%NS
           NF = FouCoil(icoil)%NF  ! allias variable for simplicity;
           !-------------------------calculate coil data-------------------------------------------------  
           mm = 0
           coil(icoil)%xx(0:NS) = FouCoil(icoil)%cmt(0:NS,mm) * Foucoil(icoil)%xc(mm)
           coil(icoil)%yy(0:NS) = FouCoil(icoil)%cmt(0:NS,mm) * Foucoil(icoil)%yc(mm)
           coil(icoil)%zz(0:NS) = FouCoil(icoil)%cmt(0:NS,mm) * Foucoil(icoil)%zc(mm)    
           do mm = 1, NF   
              coil(icoil)%xx(0:NS) = coil(icoil)%xx(0:NS) + (   FouCoil(icoil)%cmt(0:NS,mm) * Foucoil(icoil)%xc(mm) &
                                                              + FouCoil(icoil)%smt(0:NS,mm) * Foucoil(icoil)%xs(mm) )
              coil(icoil)%yy(0:NS) = coil(icoil)%yy(0:NS) + (   FouCoil(icoil)%cmt(0:NS,mm) * Foucoil(icoil)%yc(mm) &
                                                              + FouCoil(icoil)%smt(0:NS,mm) * Foucoil(icoil)%ys(mm) )
              coil(icoil)%zz(0:NS) = coil(icoil)%zz(0:NS) + (   FouCoil(icoil)%cmt(0:NS,mm) * Foucoil(icoil)%zc(mm) &
                                                              + FouCoil(icoil)%smt(0:NS,mm) * Foucoil(icoil)%zs(mm) )
              coil(icoil)%xt(0:NS) = coil(icoil)%xt(0:NS) + ( - FouCoil(icoil)%smt(0:NS,mm) * Foucoil(icoil)%xc(mm) &
                                                              + FouCoil(icoil)%cmt(0:NS,mm) * Foucoil(icoil)%xs(mm) ) * mm
              coil(icoil)%yt(0:NS) = coil(icoil)%yt(0:NS) + ( - FouCoil(icoil)%smt(0:NS,mm) * Foucoil(icoil)%yc(mm) &
                                                              + FouCoil(icoil)%cmt(0:NS,mm) * Foucoil(icoil)%ys(mm) ) * mm
              coil(icoil)%zt(0:NS) = coil(icoil)%zt(0:NS) + ( - FouCoil(icoil)%smt(0:NS,mm) * Foucoil(icoil)%zc(mm) &
                                                              + FouCoil(icoil)%cmt(0:NS,mm) * Foucoil(icoil)%zs(mm) ) * mm
              coil(icoil)%xa(0:NS) = coil(icoil)%xa(0:NS) + ( - FouCoil(icoil)%cmt(0:NS,mm) * Foucoil(icoil)%xc(mm) &
                                                              - FouCoil(icoil)%smt(0:NS,mm) * Foucoil(icoil)%xs(mm) ) * mm*mm
              coil(icoil)%ya(0:NS) = coil(icoil)%ya(0:NS) + ( - FouCoil(icoil)%cmt(0:NS,mm) * Foucoil(icoil)%yc(mm) &
                                                              - FouCoil(icoil)%smt(0:NS,mm) * Foucoil(icoil)%ys(mm) ) * mm*mm
              coil(icoil)%za(0:NS) = coil(icoil)%za(0:NS) + ( - FouCoil(icoil)%cmt(0:NS,mm) * Foucoil(icoil)%zc(mm) &
                                                              - FouCoil(icoil)%smt(0:NS,mm) * Foucoil(icoil)%zs(mm) ) * mm*mm
              coil(icoil)%xb(0:NS) = coil(icoil)%xb(0:NS) + (   FouCoil(icoil)%smt(0:NS,mm) * Foucoil(icoil)%xc(mm) &
                                                              - FouCoil(icoil)%cmt(0:NS,mm) * Foucoil(icoil)%xs(mm) ) * mm*mm*mm
              coil(icoil)%yb(0:NS) = coil(icoil)%yb(0:NS) + (   FouCoil(icoil)%smt(0:NS,mm) * Foucoil(icoil)%yc(mm) &
                                                              - FouCoil(icoil)%cmt(0:NS,mm) * Foucoil(icoil)%ys(mm) ) * mm*mm*mm
              coil(icoil)%zb(0:NS) = coil(icoil)%zb(0:NS) + (   FouCoil(icoil)%smt(0:NS,mm) * Foucoil(icoil)%zc(mm) &
                                                              - FouCoil(icoil)%cmt(0:NS,mm) * Foucoil(icoil)%zs(mm) ) * mm*mm*mm
           enddo ! end of do mm; 

        case(2)

        case(3)

        case( coil_type_multi )
           NS = coil(icoil)%NS
           NF = FouCoil(icoil)%NF
           ND = 3 + 6*NF
           SALLOCATE( zeta, (0:NS), 0.0 )
           SALLOCATE( sinterms, (0:NS), 0.0 )
           SALLOCATE( costerms, (0:NS), 0.0 )
           SALLOCATE( rc, (0:NS,1:3), 0.0 )
           SALLOCATE( rcp, (0:NS,1:3), 0.0 )
           SALLOCATE( rcpp, (0:NS,1:3), 0.0 )
           SALLOCATE( rcppp, (0:NS,1:3), 0.0 )
           SALLOCATE( rcdof, (0:NS,1:3,1:ND), 0.0 )
           SALLOCATE( rcpdof, (0:NS,1:3,1:ND), 0.0 )
           SALLOCATE( rcppdof, (0:NS,1:3,1:ND), 0.0 )
           SALLOCATE( xcdof, (1:ND), 0.0 )
           SALLOCATE( ycdof, (1:ND), 0.0 )
           SALLOCATE( zcdof, (1:ND), 0.0 )
           SALLOCATE( delta, (0:NS,1:3), 0.0 )
           SALLOCATE( deltadof, (0:NS,1:3,1:ND), 0.0 )
           SALLOCATE( n, (0:NS,1:3), 0.0 )
           SALLOCATE( np, (0:NS,1:3), 0.0 )
           SALLOCATE( npp, (0:NS,1:3), 0.0 )
           SALLOCATE( ndoff, (0:NS,1:3,1:ND), 0.0 )
           SALLOCATE( npdof, (0:NS,1:3,1:ND), 0.0 )
           SALLOCATE( nhatp, (0:NS,1:3), 0.0 )
           SALLOCATE( nhatpp, (0:NS,1:3), 0.0 )
           SALLOCATE( nhatpdof, (0:NS,1:3,1:ND), 0.0 )
           SALLOCATE( normt, (0:NS), 0.0 )
           SALLOCATE( normtp, (0:NS), 0.0 )
           SALLOCATE( normtpp, (0:NS), 0.0 )
           SALLOCATE( normn, (0:NS), 0.0 )
           SALLOCATE( normtdof, (0:NS,1:ND), 0.0 )
           SALLOCATE( normtpdof, (0:NS,1:ND), 0.0 )
           SALLOCATE( normndof, (0:NS,1:ND), 0.0 )
           SALLOCATE( that, (0:NS,1:3), 0.0 )
           SALLOCATE( nhat, (0:NS,1:3), 0.0 )
           SALLOCATE( bhat, (0:NS,1:3), 0.0 )
           SALLOCATE( bhatp, (0:NS,1:3), 0.0 )
           SALLOCATE( bhatpp, (0:NS,1:3), 0.0 )
           SALLOCATE( thatp, (0:NS,1:3), 0.0 )
           SALLOCATE( thatpp, (0:NS,1:3), 0.0 )
           SALLOCATE( thatdof, (0:NS,1:3,1:ND), 0.0 )
           SALLOCATE( nhatdof, (0:NS,1:3,1:ND), 0.0 )
           SALLOCATE( bhatdof, (0:NS,1:3,1:ND), 0.0 )
           SALLOCATE( nhata, (0:NS,1:3), 0.0 )
           SALLOCATE( bhata, (0:NS,1:3), 0.0 )
           SALLOCATE( nhatadof, (0:NS,1:3,1:ND), 0.0 )
           SALLOCATE( bhatadof, (0:NS,1:3,1:ND), 0.0 )
           SALLOCATE( nhatap, (0:NS,1:3), 0.0 )
           SALLOCATE( bhatap, (0:NS,1:3), 0.0 )
           SALLOCATE( nhatapp, (0:NS,1:3), 0.0 )
           SALLOCATE( bhatapp, (0:NS,1:3), 0.0 )
           SALLOCATE( alphac, (0:Nalpha), 0.0 )
           SALLOCATE( alphas, (0:Nalpha), 0.0 )
           SALLOCATE( alphacdof, (0:Nalpha,1:ND), 0.0 )
           SALLOCATE( alphasdof, (0:Nalpha,1:ND), 0.0 )
           SALLOCATE( Cdoff, (1:ND), 0.0 )
           SALLOCATE( g, (0:NS), 0.0 )
           SALLOCATE( w, (0:NS), 0.0 )
           SALLOCATE( h, (0:NS), 0.0 )
           SALLOCATE( gp, (0:NS), 0.0 )
           SALLOCATE( wp, (0:NS), 0.0 )
           SALLOCATE( hp, (0:NS), 0.0 )
           SALLOCATE( gpp, (0:NS), 0.0 )
           SALLOCATE( wpp, (0:NS), 0.0 )
           SALLOCATE( hpp, (0:NS), 0.0 )
           SALLOCATE( gdof, (0:NS,1:ND), 0.0 )
           SALLOCATE( wdof, (0:NS,1:ND), 0.0 )
           SALLOCATE( hdof, (0:NS,1:ND), 0.0 )
           SALLOCATE( lengdof, (1:ND), 0.0 )
           SALLOCATE( nhatpfd, (0:NS,1:3), 0.0 )
           SALLOCATE( nhatpdoffd, (0:NS,1:3,1:ND), 0.0 )

           ! Add in zeta_0 phase shift later 
           do i = 0, NS
              zeta(i) = real(i)*pi2/real(NS)
           enddo
           mm = 0
           rc(0:NS,1) = Foucoil(icoil)%xc(mm)
           rc(0:NS,2) = Foucoil(icoil)%yc(mm)
           rc(0:NS,3) = Foucoil(icoil)%zc(mm)
           rcdof(0:NS,1,1) = 1.0
           rcdof(0:NS,2,2+2*NF) = 1.0
           rcdof(0:NS,3,3+4*NF) = 1.0
           do mm = 1, NF
              freq = real(mm*Nturns*Npancakes)
              costerms(0:NS) = cos(mm*Nturns*Npancakes*zeta(0:NS))
              sinterms(0:NS) = sin(mm*Nturns*Npancakes*zeta(0:NS))
              rc(0:NS,1)    = rc(0:NS,1)    +                Foucoil(icoil)%xc(mm)*costerms(0:NS) +                Foucoil(icoil)%xs(mm)*sinterms(0:NS)
              rc(0:NS,2)    = rc(0:NS,2)    +                Foucoil(icoil)%yc(mm)*costerms(0:NS) +                Foucoil(icoil)%ys(mm)*sinterms(0:NS)
              rc(0:NS,3)    = rc(0:NS,3)    +                Foucoil(icoil)%zc(mm)*costerms(0:NS) +                Foucoil(icoil)%zs(mm)*sinterms(0:NS)
              rcp(0:NS,1)   = rcp(0:NS,1)   -           freq*Foucoil(icoil)%xc(mm)*sinterms(0:NS) +           freq*Foucoil(icoil)%xs(mm)*costerms(0:NS)
              rcp(0:NS,2)   = rcp(0:NS,2)   -           freq*Foucoil(icoil)%yc(mm)*sinterms(0:NS) +           freq*Foucoil(icoil)%ys(mm)*costerms(0:NS)
              rcp(0:NS,3)   = rcp(0:NS,3)   -           freq*Foucoil(icoil)%zc(mm)*sinterms(0:NS) +           freq*Foucoil(icoil)%zs(mm)*costerms(0:NS)
              rcpp(0:NS,1)  = rcpp(0:NS,1)  -      freq*freq*Foucoil(icoil)%xc(mm)*costerms(0:NS) -      freq*freq*Foucoil(icoil)%xs(mm)*sinterms(0:NS)
              rcpp(0:NS,2)  = rcpp(0:NS,2)  -      freq*freq*Foucoil(icoil)%yc(mm)*costerms(0:NS) -      freq*freq*Foucoil(icoil)%ys(mm)*sinterms(0:NS)
              rcpp(0:NS,3)  = rcpp(0:NS,3)  -      freq*freq*Foucoil(icoil)%zc(mm)*costerms(0:NS) -      freq*freq*Foucoil(icoil)%zs(mm)*sinterms(0:NS)
              rcppp(0:NS,1) = rcppp(0:NS,1) + freq*freq*freq*Foucoil(icoil)%xc(mm)*sinterms(0:NS) - freq*freq*freq*Foucoil(icoil)%xs(mm)*costerms(0:NS)
              rcppp(0:NS,2) = rcppp(0:NS,2) + freq*freq*freq*Foucoil(icoil)%yc(mm)*sinterms(0:NS) - freq*freq*freq*Foucoil(icoil)%ys(mm)*costerms(0:NS)
              rcppp(0:NS,3) = rcppp(0:NS,3) + freq*freq*freq*Foucoil(icoil)%zc(mm)*sinterms(0:NS) - freq*freq*freq*Foucoil(icoil)%zs(mm)*costerms(0:NS)
              
              rcdof(0:NS,1,1+     mm) = costerms(0:NS)
              rcdof(0:NS,1,1+  NF+mm) = sinterms(0:NS)
              rcdof(0:NS,2,2+2*NF+mm) = costerms(0:NS)
              rcdof(0:NS,2,2+3*NF+mm) = sinterms(0:NS)
              rcdof(0:NS,3,3+4*NF+mm) = costerms(0:NS)
              rcdof(0:NS,3,3+5*NF+mm) = sinterms(0:NS)
              rcpdof(0:NS,1,1+     mm) = -1.0*sinterms(0:NS)*freq
              rcpdof(0:NS,1,1+  NF+mm) =      costerms(0:NS)*freq
              rcpdof(0:NS,2,2+2*NF+mm) = -1.0*sinterms(0:NS)*freq
              rcpdof(0:NS,2,2+3*NF+mm) =      costerms(0:NS)*freq
              rcpdof(0:NS,3,3+4*NF+mm) = -1.0*sinterms(0:NS)*freq
              rcpdof(0:NS,3,3+5*NF+mm) =      costerms(0:NS)*freq
              rcppdof(0:NS,1,1+     mm) = -1.0*costerms(0:NS)*freq*freq
              rcppdof(0:NS,1,1+  NF+mm) = -1.0*sinterms(0:NS)*freq*freq
              rcppdof(0:NS,2,2+2*NF+mm) = -1.0*costerms(0:NS)*freq*freq
              rcppdof(0:NS,2,2+3*NF+mm) = -1.0*sinterms(0:NS)*freq*freq
              rcppdof(0:NS,3,3+4*NF+mm) = -1.0*costerms(0:NS)*freq*freq
              rcppdof(0:NS,3,3+5*NF+mm) = -1.0*sinterms(0:NS)*freq*freq
           enddo
           normt(0:NS) = sqrt( rcp(0:NS,1)**2.0 + rcp(0:NS,2)**2.0 + rcp(0:NS,3)**2.0 )
           normtp(0:NS) = normt(0:NS)**-1.0*(rcp(0:NS,1)*rcpp(0:NS,1) + rcp(0:NS,2)*rcpp(0:NS,2) + rcp(0:NS,3)*rcpp(0:NS,3))
           normtpp(0:NS) = normt(0:NS)**-1.0*(rcpp(0:NS,1)*rcpp(0:NS,1)+rcpp(0:NS,2)*rcpp(0:NS,2)+rcpp(0:NS,3)*rcpp(0:NS,3)+rcp(0:NS,1)*rcppp(0:NS,1)+rcp(0:NS,2)*rcppp(0:NS,2)+rcp(0:NS,3)*rcppp(0:NS,3)) - &
                normt(0:NS)**-3.0*(rcp(0:NS,1)*rcpp(0:NS,1) + rcp(0:NS,2)*rcpp(0:NS,2) + rcp(0:NS,3)*rcpp(0:NS,3))**2.0
           do k = 1, ND
              if( myid.ne.modulo(k-1,ncpu) ) cycle
              normtdof(0:NS,k) = normt(0:NS)**-1.0*( rcp(0:NS,1)*rcpdof(0:NS,1,k) + rcp(0:NS,2)*rcpdof(0:NS,2,k) + rcp(0:NS,3)*rcpdof(0:NS,3,k) )
              normtpdof(0:NS,k) = -1.0*normt(0:NS)**-2.0*normtdof(0:NS,k)*(rcp(0:NS,1)*rcpp(0:NS,1) + rcp(0:NS,2)*rcpp(0:NS,2) + rcp(0:NS,3)*rcpp(0:NS,3)) + &
                   normt(0:NS)**-1.0*(rcpdof(0:NS,1,k)*rcpp(0:NS,1)+rcpdof(0:NS,2,k)*rcpp(0:NS,2)+rcpdof(0:NS,3,k)*rcpp(0:NS,3)+&
                   rcp(0:NS,1)*rcppdof(0:NS,1,k)+rcp(0:NS,2)*rcppdof(0:NS,2,k)+rcp(0:NS,3)*rcppdof(0:NS,3,k))
           enddo
           do j = 1, 3
              that(0:NS,j) = rcp(0:NS,j) / normt(0:NS)
              do k = 1, ND
                 if( myid.ne.modulo(k-1,ncpu) ) cycle
                 thatdof(0:NS,j,k) = rcpdof(0:NS,j,k) / normt(0:NS) - normtdof(0:NS,k)*rcp(0:NS,j) / normt(0:NS)**2.0
              enddo
           enddo
           xc = sum( rc(0:NS-1,1)*normt(0:NS-1) ) / sum( normt(0:NS-1) )
           yc = sum( rc(0:NS-1,2)*normt(0:NS-1) ) / sum( normt(0:NS-1) )
           zc = sum( rc(0:NS-1,3)*normt(0:NS-1) ) / sum( normt(0:NS-1) )
           do k = 1, ND
              if( myid.ne.modulo(k-1,ncpu) ) cycle
              xcdof(k) = sum( rcdof(0:NS-1,1,k)*normt(0:NS-1) + rc(0:NS-1,1)*normtdof(0:NS-1,k) ) / sum( normt(0:NS-1) ) - &
                         sum( normtdof(0:NS-1,k) )*sum( rc(0:NS-1,1)*normt(0:NS-1) ) / sum( normt(0:NS-1) )**2.0
              ycdof(k) = sum( rcdof(0:NS-1,2,k)*normt(0:NS-1) + rc(0:NS-1,2)*normtdof(0:NS-1,k) ) / sum( normt(0:NS-1) ) - &
                         sum( normtdof(0:NS-1,k) )*sum( rc(0:NS-1,2)*normt(0:NS-1) ) / sum( normt(0:NS-1) )**2.0
              zcdof(k) = sum( rcdof(0:NS-1,3,k)*normt(0:NS-1) + rc(0:NS-1,3)*normtdof(0:NS-1,k) ) / sum( normt(0:NS-1) ) - &
                         sum( normtdof(0:NS-1,k) )*sum( rc(0:NS-1,3)*normt(0:NS-1) ) / sum( normt(0:NS-1) )**2.0
           enddo
           delta(0:NS,1) = rc(0:NS,1) - xc
           delta(0:NS,2) = rc(0:NS,2) - yc
           delta(0:NS,3) = rc(0:NS,3) - zc
           do k = 1, ND
              if( myid.ne.modulo(k-1,ncpu) ) cycle
              deltadof(0:NS,1,k) = rcdof(0:NS,1,k) - xcdof(k)
              deltadof(0:NS,2,k) = rcdof(0:NS,2,k) - ycdof(k)
              deltadof(0:NS,3,k) = rcdof(0:NS,3,k) - zcdof(k)
           enddo
           do j = 1, 3
              n(0:NS,j) = delta(0:NS,j) - (delta(0:NS,1)*rcp(0:NS,1) + delta(0:NS,2)*rcp(0:NS,2) + delta(0:NS,3)*rcp(0:NS,3)) * rcp(0:NS,j) * normt(0:NS)**-2.0
              do k = 1, ND
                 if( myid.ne.modulo(k-1,ncpu) ) cycle
                 ndoff(0:NS,j,k) = deltadof(0:NS,j,k) - (deltadof(0:NS,1,k)*rcp(0:NS,1)+deltadof(0:NS,2,k)*rcp(0:NS,2)+deltadof(0:NS,3,k)*rcp(0:NS,3)+&
                      delta(0:NS,1)*rcpdof(0:NS,1,k)+delta(0:NS,2)*rcpdof(0:NS,2,k)+delta(0:NS,3)*rcpdof(0:NS,3,k))*rcp(0:NS,j)/(rcp(0:NS,1)**2.0+rcp(0:NS,2)**2.0+rcp(0:NS,3)**2.0) - &
                      (delta(0:NS,1)*rcp(0:NS,1)+delta(0:NS,2)*rcp(0:NS,2)+delta(0:NS,3)*rcp(0:NS,3))*rcpdof(0:NS,j,k)/(rcp(0:NS,1)**2.0+rcp(0:NS,2)**2.0+rcp(0:NS,3)**2.0) + &
                      2.0*(rcp(0:NS,1)*rcpdof(0:NS,1,k)+rcp(0:NS,2)*rcpdof(0:NS,2,k)+rcp(0:NS,3)*rcpdof(0:NS,3,k))*&
                      (delta(0:NS,1)*rcp(0:NS,1)+delta(0:NS,2)*rcp(0:NS,2)+delta(0:NS,3)*rcp(0:NS,3))*rcp(0:NS,j)/(rcp(0:NS,1)**2.0+rcp(0:NS,2)**2.0+rcp(0:NS,3)**2.0)**2.0
              enddo
           enddo
           normn(0:NS) = sqrt( n(0:NS,1)**2.0 + n(0:NS,2)**2.0 + n(0:NS,3)**2.0 )
           do k = 1, ND
              if( myid.ne.modulo(k-1,ncpu) ) cycle
              normndof(0:NS,k) = ( n(0:NS,1)**2.0 + n(0:NS,2)**2.0 + n(0:NS,3)**2.0 )**-0.5*( n(0:NS,1)*ndoff(0:NS,1,k) + n(0:NS,2)*ndoff(0:NS,2,k) + n(0:NS,3)*ndoff(0:NS,3,k) )
           enddo
           do j = 1, 3
              nhat(0:NS,j) = n(0:NS,j) / normn(0:NS)
              do k = 1, ND
                 if( myid.ne.modulo(k-1,ncpu) ) cycle
                 nhatdof(0:NS,j,k) = ndoff(0:NS,j,k) / normn(0:NS) - n(0:NS,j)*normndof(0:NS,k) / normn(0:NS)**2.0
              enddo
           enddo
           bhat(0:NS,1) = that(0:NS,2)*nhat(0:NS,3) - that(0:NS,3)*nhat(0:NS,2)
           bhat(0:NS,2) = that(0:NS,3)*nhat(0:NS,1) - that(0:NS,1)*nhat(0:NS,3)
           bhat(0:NS,3) = that(0:NS,1)*nhat(0:NS,2) - that(0:NS,2)*nhat(0:NS,1)
           do k = 1, ND
              if( myid.ne.modulo(k-1,ncpu) ) cycle
              bhatdof(0:NS,1,k) = thatdof(0:NS,2,k)*nhat(0:NS,3) - thatdof(0:NS,3,k)*nhat(0:NS,2) + that(0:NS,2)*nhatdof(0:NS,3,k) - that(0:NS,3)*nhatdof(0:NS,2,k)
              bhatdof(0:NS,2,k) = thatdof(0:NS,3,k)*nhat(0:NS,1) - thatdof(0:NS,1,k)*nhat(0:NS,3) + that(0:NS,3)*nhatdof(0:NS,1,k) - that(0:NS,1)*nhatdof(0:NS,3,k)
              bhatdof(0:NS,3,k) = thatdof(0:NS,1,k)*nhat(0:NS,2) - thatdof(0:NS,2,k)*nhat(0:NS,1) + that(0:NS,1)*nhatdof(0:NS,2,k) - that(0:NS,2)*nhatdof(0:NS,1,k)
           enddo
           do j = 1, 3
              thatp(0:NS,j) = rcpp(0:NS,j)/normt(0:NS) - ( rcp(0:NS,1)*rcpp(0:NS,1) + rcp(0:NS,2)*rcpp(0:NS,2) + rcp(0:NS,3)*rcpp(0:NS,3) )*rcp(0:NS,j) / normt(0:NS)**3.0
              thatpp(0:NS,j) = rcppp(0:NS,j)/normt(0:NS) - 2.0*( rcp(0:NS,1)*rcpp(0:NS,1) + rcp(0:NS,2)*rcpp(0:NS,2) + rcp(0:NS,3)*rcpp(0:NS,3) )*rcpp(0:NS,j) / normt(0:NS)**3.0 - &
                   (rcpp(0:NS,1)*rcpp(0:NS,1)+rcpp(0:NS,2)*rcpp(0:NS,2)+rcpp(0:NS,3)*rcpp(0:NS,3)+rcp(0:NS,1)*rcppp(0:NS,1)+rcp(0:NS,2)*rcppp(0:NS,2)+rcp(0:NS,3)*rcppp(0:NS,3))*rcp(0:NS,j) / normt(0:NS)**3.0 + &
                   3.0*( rcp(0:NS,1)*rcpp(0:NS,1) + rcp(0:NS,2)*rcpp(0:NS,2) + rcp(0:NS,3)*rcpp(0:NS,3) )*rcp(0:NS,j) / normt(0:NS)**4.0 * normtp(0:NS)
           enddo
           do j = 1, 3
              np(0:NS,j) = rcp(0:NS,j) - ( normt(0:NS)**2.0+delta(0:NS,1)*rcpp(0:NS,1)+delta(0:NS,2)*rcpp(0:NS,2)+delta(0:NS,3)*rcpp(0:NS,3) )*rcp(0:NS,j)*normt(0:NS)**-2.0 + &
                   ( delta(0:NS,1)*rcp(0:NS,1)+delta(0:NS,2)*rcp(0:NS,2)+delta(0:NS,3)*rcp(0:NS,3) )*( rcp(0:NS,j)*normt(0:NS)**-3.0*2.0*normtp(0:NS)-rcpp(0:NS,j)*normt(0:NS)**-2.0 )
              npp(0:NS,j) = rcpp(0:NS,j) - ( 2.0*normt(0:NS)*normtp(0:NS) + rcp(0:NS,1)*rcpp(0:NS,1)+rcp(0:NS,2)*rcpp(0:NS,2)+rcp(0:NS,3)*rcpp(0:NS,3) + &
                   delta(0:NS,1)*rcppp(0:NS,1)+delta(0:NS,2)*rcppp(0:NS,2)+delta(0:NS,3)*rcppp(0:NS,3) )*rcp(0:NS,j)*normt(0:NS)**-2.0 + &
                   2.0*( normt(0:NS)**2.0+delta(0:NS,1)*rcpp(0:NS,1)+delta(0:NS,2)*rcpp(0:NS,2)+delta(0:NS,3)*rcpp(0:NS,3) )*( rcp(0:NS,j)*normt(0:NS)**-3.0*2.0*normtp(0:NS)-rcpp(0:NS,j)*normt(0:NS)**-2.0 ) + &
                   ( delta(0:NS,1)*rcp(0:NS,1)+delta(0:NS,2)*rcp(0:NS,2)+delta(0:NS,3)*rcp(0:NS,3) )*( rcpp(0:NS,j)*normt(0:NS)**-3.0*2.0*normtp(0:NS) - &
                   3.0*rcp(0:NS,j)*normt(0:NS)**-4.0*2.0*normtp(0:NS)**2.0 + rcp(0:NS,j)*normt(0:NS)**-3.0*2.0*normtpp(0:NS) - &
                   rcppp(0:NS,j)*normt(0:NS)**-2.0 + 2.0*rcpp(0:NS,j)*normt(0:NS)**-3.0*normtp(0:NS) )
              do k = 1, ND
                 if( myid.ne.modulo(k-1,ncpu) ) cycle
                 npdof(0:NS,j,k) = rcpdof(0:NS,j,k) - ( 2.0*normt(0:NS)*normtdof(0:NS,k)+deltadof(0:NS,1,k)*rcpp(0:NS,1)+deltadof(0:NS,2,k)*rcpp(0:NS,2)+deltadof(0:NS,3,k)*rcpp(0:NS,3)+&
                      delta(0:NS,1)*rcppdof(0:NS,1,k)+delta(0:NS,2)*rcppdof(0:NS,2,k)+delta(0:NS,3)*rcppdof(0:NS,3,k) )*rcp(0:NS,j)*normt(0:NS)**-2.0 + &
                      ( normt(0:NS)**2.0+delta(0:NS,1)*rcpp(0:NS,1)+delta(0:NS,2)*rcpp(0:NS,2)+delta(0:NS,3)*rcpp(0:NS,3) )*&
                      ( 2.0*rcp(0:NS,j)*normt(0:NS)**-3.0*normtdof(0:NS,k)-rcpdof(0:NS,j,k)*normt(0:NS)**-2.0 ) + &
                      ( deltadof(0:NS,1,k)*rcp(0:NS,1)+deltadof(0:NS,2,k)*rcp(0:NS,2)+deltadof(0:NS,3,k)*rcp(0:NS,3)+delta(0:NS,1)*rcpdof(0:NS,1,k)+&
                      delta(0:NS,2)*rcpdof(0:NS,2,k)+delta(0:NS,3)*rcpdof(0:NS,3,k) )*( rcp(0:NS,j)*normt(0:NS)**-3.0*2.0*normtp(0:NS)-rcpp(0:NS,j)*normt(0:NS)**-2.0 ) + &
                      ( delta(0:NS,1)*rcp(0:NS,1)+delta(0:NS,2)*rcp(0:NS,2)+delta(0:NS,3)*rcp(0:NS,3) )*( rcpdof(0:NS,j,k)*normt(0:NS)**-3.0*2.0*normtp(0:NS)-&
                      3.0*rcp(0:NS,j)*normt(0:NS)**-4.0*normtdof(0:NS,k)*2.0*normtp(0:NS)+rcp(0:NS,j)*normt(0:NS)**-3.0*2.0*normtpdof(0:NS,k)-&
                      rcppdof(0:NS,j,k)*normt(0:NS)**-2.0+2.0*rcpp(0:NS,j)*normt(0:NS)**-3.0*normtdof(0:NS,k) )
              enddo
           enddo
           do j = 1, 3
              nhatp(0:NS,j) = np(0:NS,j)/normn(0:NS) - ( n(0:NS,1)*np(0:NS,1) + n(0:NS,2)*np(0:NS,2) + n(0:NS,3)*np(0:NS,3) )*n(0:NS,j) / normn(0:NS)**3.0
              nhatpp(0:NS,j) = npp(0:NS,j)/normn(0:NS) - 2.0*( n(0:NS,1)*np(0:NS,1) + n(0:NS,2)*np(0:NS,2) + n(0:NS,3)*np(0:NS,3) )*np(0:NS,j) / normn(0:NS)**3.0 - &
                   ( np(0:NS,1)*np(0:NS,1) + np(0:NS,2)*np(0:NS,2) + np(0:NS,3)*np(0:NS,3) + n(0:NS,1)*npp(0:NS,1) + n(0:NS,2)*npp(0:NS,2) + n(0:NS,3)*npp(0:NS,3) )*n(0:NS,j) / normn(0:NS)**3.0 + &
                   3.0*( n(0:NS,1)*np(0:NS,1) + n(0:NS,2)*np(0:NS,2) + n(0:NS,3)*np(0:NS,3) )*n(0:NS,j) / normn(0:NS)**5.0 * ( n(0:NS,1)*np(0:NS,1) + n(0:NS,2)*np(0:NS,2) + n(0:NS,3)*np(0:NS,3) )
              do k = 1, ND
                 if( myid.ne.modulo(k-1,ncpu) ) cycle
                 nhatpdof(0:NS,j,k) = npdof(0:NS,j,k)/normn(0:NS) - np(0:NS,j)*normn(0:NS)**-2.0*normndof(0:NS,k) - &
                      (ndoff(0:NS,1,k)*np(0:NS,1)+ndoff(0:NS,2,k)*np(0:NS,2)+ndoff(0:NS,3,k)*np(0:NS,3)+n(0:NS,1)*npdof(0:NS,1,k)+n(0:NS,2)*npdof(0:NS,2,k)+n(0:NS,3)*npdof(0:NS,3,k))*n(0:NS,j)/normn(0:NS)**3.0 - &
                      ( n(0:NS,1)*np(0:NS,1) + n(0:NS,2)*np(0:NS,2) + n(0:NS,3)*np(0:NS,3) )*( ndoff(0:NS,j,k) / normn(0:NS)**3.0 - 3.0*n(0:NS,j)*normn(0:NS)**-4.0*normndof(0:NS,k) )
              enddo
           enddo
           bhatp(0:NS,1) = thatp(0:NS,2)*nhat(0:NS,3) - thatp(0:NS,3)*nhat(0:NS,2) + that(0:NS,2)*nhatp(0:NS,3) - that(0:NS,3)*nhatp(0:NS,2)
           bhatp(0:NS,2) = thatp(0:NS,3)*nhat(0:NS,1) - thatp(0:NS,1)*nhat(0:NS,3) + that(0:NS,3)*nhatp(0:NS,1) - that(0:NS,1)*nhatp(0:NS,3)
           bhatp(0:NS,3) = thatp(0:NS,1)*nhat(0:NS,2) - thatp(0:NS,2)*nhat(0:NS,1) + that(0:NS,1)*nhatp(0:NS,2) - that(0:NS,2)*nhatp(0:NS,1)
           bhatpp(0:NS,1) = thatpp(0:NS,2)*nhat(0:NS,3) - thatpp(0:NS,3)*nhat(0:NS,2) + thatp(0:NS,2)*nhatp(0:NS,3) - thatp(0:NS,3)*nhatp(0:NS,2) + &
                thatp(0:NS,2)*nhatp(0:NS,3) - thatp(0:NS,3)*nhatp(0:NS,2) + that(0:NS,2)*nhatpp(0:NS,3) - that(0:NS,3)*nhatpp(0:NS,2)
           bhatpp(0:NS,2) = thatpp(0:NS,3)*nhat(0:NS,1) - thatpp(0:NS,1)*nhat(0:NS,3) + thatp(0:NS,3)*nhatp(0:NS,1) - thatp(0:NS,1)*nhatp(0:NS,3) + &
                thatp(0:NS,3)*nhatp(0:NS,1) - thatp(0:NS,1)*nhatp(0:NS,3) + that(0:NS,3)*nhatpp(0:NS,1) - that(0:NS,1)*nhatpp(0:NS,3)
           bhatpp(0:NS,3) = thatpp(0:NS,1)*nhat(0:NS,2) - thatpp(0:NS,2)*nhat(0:NS,1) + thatp(0:NS,1)*nhatp(0:NS,2) - thatp(0:NS,2)*nhatp(0:NS,1) + &
                thatp(0:NS,1)*nhatp(0:NS,2) - thatp(0:NS,2)*nhatp(0:NS,1) + that(0:NS,1)*nhatpp(0:NS,2) - that(0:NS,2)*nhatpp(0:NS,1)
           if( ifirst .eq. 1 .or. AlphaPert .eq. 1 ) then
              
              ! Improve robustness
              ! Add in increased NS
              !do i = 1, NS-1
              !   nhatpfd(i,1:3) = real(NS) * ( nhat(i+1,1:3) - nhat(i-1,1:3) ) / ( 2.0*pi2 )
              !enddo
              !nhatpfd(0,1:3) = real(NS) * ( nhat(1,1:3) - nhat(NS-1,1:3) ) / ( 2.0*pi2 )
              !nhatpfd(NS,1:3) = nhatpfd(0,1:3)
              !do k = 1, ND
              !   if( myid.ne.modulo(k-1,ncpu) ) cycle
              !   do i = 1, NS-1
              !      nhatpdoffd(i,1:3,k) = real(NS) * ( nhatdof(i+1,1:3,k) - nhatdof(i-1,1:3,k) ) / ( 2.0*pi2 )
              !   enddo
              !   nhatpdoffd(0,1:3,k) = real(NS) * ( nhatdof(1,1:3,k) - nhatdof(NS-1,1:3,k) ) / ( 2.0*pi2 )
              !   nhatpdoffd(NS,1:3,k) = nhatpdoffd(0,1:3,k)
              !enddo
              !if ( Alphafd .eq. 1 ) then
              !   nhatp(0:NS,1:3) = nhatpfd(0:NS,1:3)
              !   do k = 1, ND
              !      if( myid.ne.modulo(k-1,ncpu) ) cycle
              !      nhatpdof(0:NS,1:3,k) = nhatpdoffd(0:NS,1:3,k)
              !   enddo
              !endif

              C = sum(nhatp(0:NS-1,1)*bhat(0:NS-1,1)+nhatp(0:NS-1,2)*bhat(0:NS-1,2)+nhatp(0:NS-1,3)*bhat(0:NS-1,3))/(sum(normt(0:NS-1)))
              do k = 1, ND
                 if( myid.ne.modulo(k-1,ncpu) ) cycle
                 Cdoff(k) = sum(nhatpdof(0:NS-1,1,k)*bhat(0:NS-1,1)+nhatpdof(0:NS-1,2,k)*bhat(0:NS-1,2)+nhatpdof(0:NS-1,3,k)*bhat(0:NS-1,3)+&
                      nhatp(0:NS-1,1)*bhatdof(0:NS-1,1,k)+nhatp(0:NS-1,2)*bhatdof(0:NS-1,2,k)+nhatp(0:NS-1,3)*bhatdof(0:NS-1,3,k))/(sum(normt(0:NS-1))) - &
                      sum(nhatp(0:NS-1,1)*bhat(0:NS-1,1)+nhatp(0:NS-1,2)*bhat(0:NS-1,2)+nhatp(0:NS-1,3)*bhat(0:NS-1,3))*(sum(normt(0:NS-1)))**-2.0*(sum(normtdof(0:NS-1,k)))
              enddo

              do i = 1, Nalpha
                 alphac(i) = -2.0*sum(sin(Nturns*Npancakes*real(i)*zeta(0:NS-1))*(nhatp(0:NS-1,1)*bhat(0:NS-1,1)+&
                      nhatp(0:NS-1,2)*bhat(0:NS-1,2)+nhatp(0:NS-1,3)*bhat(0:NS-1,3)-C*normt(0:NS-1)))/(Nturns*Npancakes*real(NS)*real(i))
                 alphas(i) =  2.0*sum(cos(Nturns*Npancakes*real(i)*zeta(0:NS-1))*(nhatp(0:NS-1,1)*bhat(0:NS-1,1)+&
                      nhatp(0:NS-1,2)*bhat(0:NS-1,2)+nhatp(0:NS-1,3)*bhat(0:NS-1,3)-C*normt(0:NS-1)))/(Nturns*Npancakes*real(NS)*real(i))
                 do k = 1, ND
                    if( myid.ne.modulo(k-1,ncpu) ) cycle
                    alphacdof(i,k) = -2.0*sum(sin(Nturns*Npancakes*real(i)*zeta(0:NS-1))*(nhatpdof(0:NS-1,1,k)*bhat(0:NS-1,1)+nhatpdof(0:NS-1,2,k)*bhat(0:NS-1,2)+nhatpdof(0:NS-1,3,k)*bhat(0:NS-1,3)+&
                         nhatp(0:NS-1,1)*bhatdof(0:NS-1,1,k)+nhatp(0:NS-1,2)*bhatdof(0:NS-1,2,k)+nhatp(0:NS-1,3)*bhatdof(0:NS-1,3,k)-&
                         Cdoff(k)*normt(0:NS-1)-C*normtdof(0:NS-1,k)))/(Nturns*Npancakes*real(NS)*real(i))
                    alphasdof(i,k) =  2.0*sum(cos(Nturns*Npancakes*real(i)*zeta(0:NS-1))*(nhatpdof(0:NS-1,1,k)*bhat(0:NS-1,1)+nhatpdof(0:NS-1,2,k)*bhat(0:NS-1,2)+nhatpdof(0:NS-1,3,k)*bhat(0:NS-1,3)+&
                         nhatp(0:NS-1,1)*bhatdof(0:NS-1,1,k)+nhatp(0:NS-1,2)*bhatdof(0:NS-1,2,k)+nhatp(0:NS-1,3)*bhatdof(0:NS-1,3,k)-&
                         Cdoff(k)*normt(0:NS-1)-C*normtdof(0:NS-1,k)))/(Nturns*Npancakes*real(NS)*real(i))
                 enddo
              enddo
              coil(icoil)%alpha(0:NS) = 0.0
              coil(icoil)%alphap(0:NS) = 0.0
              coil(icoil)%alphapp(0:NS) = 0.0
              coil(icoil)%alphadof(0:NS,1:ND) = 0.0
              do i = 1, Nalpha
                 coil(icoil)%alpha(0:NS) = coil(icoil)%alpha(0:NS) + alphac(i)*cos(Nturns*Npancakes*real(i)*zeta(0:NS)) + alphas(i)*sin(Nturns*Npancakes*real(i)*zeta(0:NS))
                 coil(icoil)%alphap(0:NS) = coil(icoil)%alphap(0:NS) - alphac(i)*Nturns*Npancakes*real(i)*sin(Nturns*Npancakes*real(i)*zeta(0:NS)) + &
                      alphas(i)*Nturns*Npancakes*real(i)*cos(Nturns*Npancakes*real(i)*zeta(0:NS))
                 coil(icoil)%alphapp(0:NS) = coil(icoil)%alphapp(0:NS) - alphac(i)*(Nturns*Npancakes*real(i))**2.0*cos(Nturns*Npancakes*real(i)*zeta(0:NS)) - &
                      alphas(i)*(Nturns*Npancakes*real(i))**2.0*sin(Nturns*Npancakes*real(i)*zeta(0:NS))
                 do k = 1, ND
                    if( myid.ne.modulo(k-1,ncpu) ) cycle
                    coil(icoil)%alphadof(0:NS,k) = coil(icoil)%alphadof(0:NS,k) + alphacdof(i,k)*cos(Nturns*Npancakes*real(i)*zeta(0:NS)) + alphasdof(i,k)*sin(Nturns*Npancakes*real(i)*zeta(0:NS))
                 enddo
              enddo
           endif
           do j = 1, 3
              nhata(0:NS,j) = -1.0*sin(coil(icoil)%alpha(0:NS))*bhat(0:NS,j) + cos(coil(icoil)%alpha(0:NS))*nhat(0:NS,j)
              bhata(0:NS,j) =      sin(coil(icoil)%alpha(0:NS))*nhat(0:NS,j) + cos(coil(icoil)%alpha(0:NS))*bhat(0:NS,j)
              nhatap(0:NS,j) = -1.0*cos(coil(icoil)%alpha(0:NS))*coil(icoil)%alphap(0:NS)*bhat(0:NS,j) - sin(coil(icoil)%alpha(0:NS))*coil(icoil)%alphap(0:NS)*nhat(0:NS,j) - &
                              sin(coil(icoil)%alpha(0:NS))*bhatp(0:NS,j) + cos(coil(icoil)%alpha(0:NS))*nhatp(0:NS,j)
              bhatap(0:NS,j) =      cos(coil(icoil)%alpha(0:NS))*coil(icoil)%alphap(0:NS)*nhat(0:NS,j) - sin(coil(icoil)%alpha(0:NS))*coil(icoil)%alphap(0:NS)*bhat(0:NS,j) + &
                              sin(coil(icoil)%alpha(0:NS))*nhatp(0:NS,j) + cos(coil(icoil)%alpha(0:NS))*bhatp(0:NS,j)
              nhatapp(0:NS,j) = sin(coil(icoil)%alpha(0:NS))*coil(icoil)%alphap(0:NS)**2.0*bhat(0:NS,j) - cos(coil(icoil)%alpha(0:NS))*coil(icoil)%alphap(0:NS)**2.0*nhat(0:NS,j) - &
                   cos(coil(icoil)%alpha(0:NS))*coil(icoil)%alphapp(0:NS)*bhat(0:NS,j) - sin(coil(icoil)%alpha(0:NS))*coil(icoil)%alphapp(0:NS)*nhat(0:NS,j) - &
                   cos(coil(icoil)%alpha(0:NS))*coil(icoil)%alphap(0:NS)*bhatp(0:NS,j) - sin(coil(icoil)%alpha(0:NS))*coil(icoil)%alphap(0:NS)*nhatp(0:NS,j) - &
                   cos(coil(icoil)%alpha(0:NS))*coil(icoil)%alphap(0:NS)*bhatp(0:NS,j) - sin(coil(icoil)%alpha(0:NS))*bhatpp(0:NS,j) - &
                   sin(coil(icoil)%alpha(0:NS))*coil(icoil)%alphap(0:NS)*nhatp(0:NS,j) + cos(coil(icoil)%alpha(0:NS))*nhatpp(0:NS,j)
              bhatapp(0:NS,j) = -1.0*sin(coil(icoil)%alpha(0:NS))*coil(icoil)%alphap(0:NS)**2.0*nhat(0:NS,j) - cos(coil(icoil)%alpha(0:NS))*coil(icoil)%alphap(0:NS)**2.0*bhat(0:NS,j) + &
                   cos(coil(icoil)%alpha(0:NS))*coil(icoil)%alphapp(0:NS)*nhat(0:NS,j) - sin(coil(icoil)%alpha(0:NS))*coil(icoil)%alphapp(0:NS)*bhat(0:NS,j) + &
                   cos(coil(icoil)%alpha(0:NS))*coil(icoil)%alphap(0:NS)*nhatp(0:NS,j) - sin(coil(icoil)%alpha(0:NS))*coil(icoil)%alphap(0:NS)*bhatp(0:NS,j) + &
                   cos(coil(icoil)%alpha(0:NS))*coil(icoil)%alphap(0:NS)*nhatp(0:NS,j) + sin(coil(icoil)%alpha(0:NS))*nhatpp(0:NS,j) - &
                   sin(coil(icoil)%alpha(0:NS))*coil(icoil)%alphap(0:NS)*bhatp(0:NS,j) + cos(coil(icoil)%alpha(0:NS))*bhatpp(0:NS,j)
              do k = 1, ND
                 if( myid.ne.modulo(k-1,ncpu) ) cycle
                 if ( AlphaPert .eq. 0 ) then
                    nhatadof(0:NS,j,k) = -1.0*sin(coil(icoil)%alpha(0:NS))*bhatdof(0:NS,j,k) + cos(coil(icoil)%alpha(0:NS))*nhatdof(0:NS,j,k)
                    bhatadof(0:NS,j,k) =      sin(coil(icoil)%alpha(0:NS))*nhatdof(0:NS,j,k) + cos(coil(icoil)%alpha(0:NS))*bhatdof(0:NS,j,k)
                 else
                    nhatadof(0:NS,j,k) = -1.0*sin(coil(icoil)%alpha(0:NS))*bhatdof(0:NS,j,k) + cos(coil(icoil)%alpha(0:NS))*nhatdof(0:NS,j,k) - &
                         cos(coil(icoil)%alpha(0:NS))*coil(icoil)%alphadof(0:NS,k)*bhat(0:NS,j) - sin(coil(icoil)%alpha(0:NS))*coil(icoil)%alphadof(0:NS,k)*nhat(0:NS,j)
                    bhatadof(0:NS,j,k) =      sin(coil(icoil)%alpha(0:NS))*nhatdof(0:NS,j,k) + cos(coil(icoil)%alpha(0:NS))*bhatdof(0:NS,j,k) + &
                         cos(coil(icoil)%alpha(0:NS))*coil(icoil)%alphadof(0:NS,k)*nhat(0:NS,j) - sin(coil(icoil)%alpha(0:NS))*coil(icoil)%alphadof(0:NS,k)*bhat(0:NS,j)
                 endif
              enddo
           enddo
           
           FATAL( discoil , Npancakes .ne. 2, Npancakes for double pancake for now )
           leng = pi2 * sum( normt(0:NS-1) ) / real(NS)
           do j = 1, ND
              lengdof(j) = pi2 * sum( normtdof(0:NS-1,j) ) / real(NS)
           enddo
           htrans = Nturn_space*real(Nturns)/leng
           FATAL( discoil , htrans .ge. leng/2.0, Turn transition too long )
           do i = 1, NS
              g(i) = g(i-1) + pi2*0.5*(normt(i-1)+normt(i))/real(NS)
              do k = 1, ND
                 if( myid.ne.modulo(k-1,ncpu) ) cycle
                 gdof(i,k) = gdof(i-1,k) + pi2*0.5*(normtdof(i-1,k)+normtdof(i,k))/real(NS)
              enddo
           enddo
           gp(0:NS) = normt(0:NS)
           gpp(0:NS) = normtp(0:NS)
           trans_length = 0.5*trans_length
           do i = 0, NS
              if( g(i) .le. trans_length ) then
                 w(i) = -0.5*Npancake_space*sin(pi2*g(i)/(4.0*trans_length))
                 wp(i) = -0.5*Npancake_space*cos(pi2*g(i)/(4.0*trans_length))*pi2*gp(i)/(4.0*trans_length)
                 wpp(i) = 0.5*Npancake_space*sin(pi2*g(i)/(4.0*trans_length))*pi2*gp(i)/(4.0*trans_length)*pi2*gp(i)/(4.0*trans_length) - &
                      0.5*Npancake_space*cos(pi2*g(i)/(4.0*trans_length))*pi2*gpp(i)/(4.0*trans_length)
                 do k = 1, ND
                    if( myid.ne.modulo(k-1,ncpu) ) cycle
                    wdof(i,k) = -0.5*Npancake_space*cos(pi2*g(i)/(4.0*trans_length))*pi2*gdof(i,k)/(4.0*trans_length)
                 enddo
              endif
              if( g(i) .gt. trans_length .and. g(i) .lt. leng/2.0 - trans_length ) then
                 w(i) = -0.5*Npancake_space
                 wp(i) = 0.0
                 wpp(i) = 0.0
                 do k = 1, ND
                    if( myid.ne.modulo(k-1,ncpu) ) cycle
                    wdof(i,k) = 0.0
                 enddo
              endif
              if( g(i) .ge. leng/2.0 - trans_length .and. g(i) .le. leng/2.0 + trans_length ) then
                 w(i) = -0.5*Npancake_space*sin(pi2*(leng/2.0-g(i))/(4.0*trans_length))
                 wp(i) = 0.5*Npancake_space*cos(pi2*(leng/2.0-g(i))/(4.0*trans_length))*pi2*gp(i)/(4.0*trans_length)
                 wpp(i) = 0.5*Npancake_space*sin(pi2*(leng/2.0-g(i))/(4.0*trans_length))*pi2*gp(i)/(4.0*trans_length)*pi2*gp(i)/(4.0*trans_length) + &
                      0.5*Npancake_space*cos(pi2*(leng/2.0-g(i))/(4.0*trans_length))*pi2*gpp(i)/(4.0*trans_length)
                 do k = 1, ND
                    if( myid.ne.modulo(k-1,ncpu) ) cycle
                    wdof(i,k) = -0.5*Npancake_space*cos(pi2*(leng/2.0-g(i))/(4.0*trans_length))*pi2*(lengdof(k)/2.0-gdof(i,k))/(4.0*trans_length)
                 enddo
              endif
              if( g(i) .gt. leng/2.0 + trans_length .and. g(i) .lt. leng - trans_length ) then
                 w(i) = 0.5*Npancake_space
                 wp(i) = 0.0
                 wpp(i) = 0.0
                 do k = 1, ND
                    if( myid.ne.modulo(k-1,ncpu) ) cycle
                    wdof(i,k) = 0.0
                 enddo
              endif
              if( g(i) .ge. leng - trans_length ) then
                 w(i) = 0.5*Npancake_space*sin(pi2*(leng-g(i))/(4.0*trans_length))
                 wp(i) = -0.5*Npancake_space*cos(pi2*(leng-g(i))/(4.0*trans_length))*pi2*gp(i)/(4.0*trans_length)
                 wpp(i) = -0.5*Npancake_space*sin(pi2*(leng-g(i))/(4.0*trans_length))*pi2*gp(i)/(4.0*trans_length)*pi2*gp(i)/(4.0*trans_length) - &
                      0.5*Npancake_space*cos(pi2*(leng-g(i))/(4.0*trans_length))*pi2*gpp(i)/(4.0*trans_length)
                 do k = 1, ND
                    if( myid.ne.modulo(k-1,ncpu) ) cycle
                    wdof(i,k) = 0.5*Npancake_space*cos(pi2*(leng-g(i))/(4.0*trans_length))*pi2*(lengdof(k)-gdof(i,k))/(4.0*trans_length)
                 enddo
              endif
              if( g(i) .le. htrans ) then
                 h(i) = 0.5*Nturn_space*real(Nturns) - g(i)**2.0 - (Nturn_space*real(Nturns)/leng)**2.0
                 hp(i) = -2.0*g(i)*gp(i)
                 hpp(i) = -2.0*gp(i)**2.0 - 2.0*g(i)*gpp(i)
                 do k = 1, ND
                    if( myid.ne.modulo(k-1,ncpu) ) cycle
                    hdof(i,k) =  -2.0*g(i)*gdof(i,k) + 2.0*(Nturn_space*real(Nturns))**2.0*leng**-3.0*lengdof(k)
                 enddo
              endif
              if( g(i) .gt. htrans .and. g(i) .lt. leng/2.0 - htrans ) then
                 h(i) = 0.5*Nturn_space*real(Nturns) - 2.0*Nturn_space*real(Nturns)*g(i)/leng
                 hp(i) = -2.0*Nturn_space*real(Nturns)*gp(i)/leng
                 hpp(i) = -2.0*Nturn_space*real(Nturns)*gpp(i)/leng
                 do k = 1, ND
                    if( myid.ne.modulo(k-1,ncpu) ) cycle
                    hdof(i,k) = -2.0*Nturn_space*real(Nturns)*(gdof(i,k)/leng - g(i)*lengdof(k)/leng**2.0)
                 enddo
              endif
              if( g(i) .ge. leng/2.0 - htrans .and. g(i) .le. leng/2.0 + htrans ) then
                 h(i) = -0.5*Nturn_space*real(Nturns) + (leng/2.0-g(i))**2.0 + (Nturn_space*real(Nturns)/leng)**2.0
                 hp(i) = -2.0*(leng/2.0-g(i))*gp(i)
                 hpp(i) = 2.0*gp(i)**2.0 - 2.0*(leng/2.0-g(i))*gpp(i)
                 do k = 1, ND
                    if( myid.ne.modulo(k-1,ncpu) ) cycle
                    hdof(i,k) = 2.0*(leng/2.0-g(i))*(lengdof(k)/2.0-gdof(i,k)) - 2.0*(Nturn_space*real(Nturns))**2.0*(leng)**-3.0*lengdof(k)
                 enddo
              endif
              if( g(i) .gt. leng/2.0 + htrans .and. g(i) .lt. leng - htrans ) then
                 h(i) = 0.5*Nturn_space*real(Nturns) - 2.0*Nturn_space*real(Nturns)*(leng-g(i))/leng
                 hp(i) = 2.0*Nturn_space*real(Nturns)*gp(i)/leng
                 hpp(i) = 2.0*Nturn_space*real(Nturns)*gpp(i)/leng
                 do k = 1, ND
                    if( myid.ne.modulo(k-1,ncpu) ) cycle
                    hdof(i,k) = -2.0*Nturn_space*real(Nturns)*(-1.0*gdof(i,k)/leng+g(i)*lengdof(k)*leng**-2.0)
                 enddo
              endif
              if( g(i) .ge. leng - htrans ) then
                 h(i) = 0.5*Nturn_space*real(Nturns) - (leng-g(i))**2.0 - (Nturn_space*real(Nturns)/leng)**2.0
                 hp(i) = 2.0*(leng-g(i))*gp(i)
                 hpp(i) = -2.0*gp(i)**2.0 + 2.0*(leng-g(i))*gpp(i)
                 do k = 1, ND
                    if( myid.ne.modulo(k-1,ncpu) ) cycle
                    hdof(i,k) = -2.0*(leng-g(i))*(lengdof(k)-gdof(i,k)) + 2.0*(Nturn_space*real(Nturns))**2.0*(leng)**-3.0*lengdof(k)
                 enddo
              endif
           enddo
           trans_length = 2.0*trans_length

           coil(icoil)%xx(0:NS) = rc(0:NS,1) + h(0:NS)*nhata(0:NS,1) + w(0:NS)*bhata(0:NS,1)
           coil(icoil)%yy(0:NS) = rc(0:NS,2) + h(0:NS)*nhata(0:NS,2) + w(0:NS)*bhata(0:NS,2)
           coil(icoil)%zz(0:NS) = rc(0:NS,3) + h(0:NS)*nhata(0:NS,3) + w(0:NS)*bhata(0:NS,3)
           coil(icoil)%xt(0:NS) = rcp(0:NS,1) + hp(0:NS)*nhata(0:NS,1) + wp(0:NS)*bhata(0:NS,1) + h(0:NS)*nhatap(0:NS,1) + w(0:NS)*bhatap(0:NS,1)
           coil(icoil)%yt(0:NS) = rcp(0:NS,2) + hp(0:NS)*nhata(0:NS,2) + wp(0:NS)*bhata(0:NS,2) + h(0:NS)*nhatap(0:NS,2) + w(0:NS)*bhatap(0:NS,2)
           coil(icoil)%zt(0:NS) = rcp(0:NS,3) + hp(0:NS)*nhata(0:NS,3) + wp(0:NS)*bhata(0:NS,3) + h(0:NS)*nhatap(0:NS,3) + w(0:NS)*bhatap(0:NS,3)

           coil(icoil)%xa(0:NS) = rcpp(0:NS,1) + hpp(0:NS)*nhata(0:NS,1) + wpp(0:NS)*bhata(0:NS,1) + 2.0*hp(0:NS)*nhatap(0:NS,1) + 2.0*wp(0:NS)*bhatap(0:NS,1) + &
                h(0:NS)*nhatapp(0:NS,1) + w(0:NS)*bhatapp(0:NS,1)
           coil(icoil)%ya(0:NS) = rcpp(0:NS,2) + hpp(0:NS)*nhata(0:NS,2) + wpp(0:NS)*bhata(0:NS,2) + 2.0*hp(0:NS)*nhatap(0:NS,2) + 2.0*wp(0:NS)*bhatap(0:NS,2) + &
                h(0:NS)*nhatapp(0:NS,2) + w(0:NS)*bhatapp(0:NS,2)
           coil(icoil)%za(0:NS) = rcpp(0:NS,3) + hpp(0:NS)*nhata(0:NS,3) + wpp(0:NS)*bhata(0:NS,3) + 2.0*hp(0:NS)*nhatap(0:NS,3) + 2.0*wp(0:NS)*bhatap(0:NS,3) + &
                h(0:NS)*nhatapp(0:NS,3) + w(0:NS)*bhatapp(0:NS,3)
           
           DoF(icoil)%xof(0:NS-1,1:ND) = 0.0
           DoF(icoil)%yof(0:NS-1,1:ND) = 0.0
           DoF(icoil)%zof(0:NS-1,1:ND) = 0.0
           do k = 1, ND
              if( myid.ne.modulo(k-1,ncpu) ) cycle
              DoF(icoil)%xof(0:NS-1,k) = rcdof(0:NS-1,1,k) + hdof(0:NS-1,k)*nhata(0:NS-1,1) + wdof(0:NS-1,k)*bhata(0:NS-1,1) + &
                                            h(0:NS-1)*nhatadof(0:NS-1,1,k) + w(0:NS-1)*bhatadof(0:NS-1,1,k)
              DoF(icoil)%yof(0:NS-1,k) = rcdof(0:NS-1,2,k) + hdof(0:NS-1,k)*nhata(0:NS-1,2) + wdof(0:NS-1,k)*bhata(0:NS-1,2) + &
                                            h(0:NS-1)*nhatadof(0:NS-1,2,k) + w(0:NS-1)*bhatadof(0:NS-1,2,k)
              DoF(icoil)%zof(0:NS-1,k) = rcdof(0:NS-1,3,k) + hdof(0:NS-1,k)*nhata(0:NS-1,3) + wdof(0:NS-1,k)*bhata(0:NS-1,3) + &
                                            h(0:NS-1)*nhatadof(0:NS-1,3,k) + w(0:NS-1)*bhatadof(0:NS-1,3,k)
           enddo
           call MPI_BARRIER( MPI_COMM_FOCUS, ierr )
           call MPI_ALLREDUCE( MPI_IN_PLACE, DoF(icoil)%xof, NS*ND, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FOCUS, ierr )
           call MPI_ALLREDUCE( MPI_IN_PLACE, DoF(icoil)%yof, NS*ND, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FOCUS, ierr )
           call MPI_ALLREDUCE( MPI_IN_PLACE, DoF(icoil)%zof, NS*ND, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FOCUS, ierr )
           call MPI_BARRIER( MPI_COMM_FOCUS, ierr )
           
           DALLOCATE( zeta )
           DALLOCATE( costerms )
           DALLOCATE( sinterms )
           DALLOCATE( rc )
           DALLOCATE( rcp )
           DALLOCATE( rcpp )
           DALLOCATE( rcppp )
           DALLOCATE( rcdof )
           DALLOCATE( rcpdof )
           DALLOCATE( rcppdof )
           DALLOCATE( xcdof )
           DALLOCATE( ycdof )
           DALLOCATE( zcdof )
           DALLOCATE( delta )
           DALLOCATE( deltadof )
           DALLOCATE( n )
           DALLOCATE( np )
           DALLOCATE( npp )
           DALLOCATE( ndoff )
           DALLOCATE( npdof )
           DALLOCATE( nhatp )
           DALLOCATE( nhatpp )
           DALLOCATE( normt )
           DALLOCATE( normtp )
           DALLOCATE( normtpp )
           DALLOCATE( normn )
           DALLOCATE( normtdof )
           DALLOCATE( normtpdof )
           DALLOCATE( normndof )
           DALLOCATE( that )
           DALLOCATE( nhat )
           DALLOCATE( bhat )
           DALLOCATE( bhatp )
           DALLOCATE( bhatpp )
           DALLOCATE( thatp )
           DALLOCATE( thatpp )
           DALLOCATE( thatdof )
           DALLOCATE( nhatdof )
           DALLOCATE( nhatpdof )
           DALLOCATE( bhatdof )
           DALLOCATE( nhata )
           DALLOCATE( bhata )
           DALLOCATE( nhatap )
           DALLOCATE( bhatap )
           DALLOCATE( nhatapp )
           DALLOCATE( bhatapp )
           DALLOCATE( nhatadof )
           DALLOCATE( bhatadof )
           DALLOCATE( alphac )
           DALLOCATE( alphas )
           DALLOCATE( alphacdof )
           DALLOCATE( alphasdof )
           DALLOCATE( Cdoff )
           DALLOCATE( g )
           DALLOCATE( w )
           DALLOCATE( h )
           DALLOCATE( gp )
           DALLOCATE( wp )
           DALLOCATE( hp )
           DALLOCATE( gpp )
           DALLOCATE( wpp )
           DALLOCATE( hpp )
           DALLOCATE( gdof )
           DALLOCATE( wdof )
           DALLOCATE( hdof )
           DALLOCATE( lengdof )
           DALLOCATE( nhatpfd )
           DALLOCATE( nhatpdoffd )

        case( coil_type_spline )
           ! reset to zero for all the coils;
           coil(icoil)%xx = zero
           coil(icoil)%yy = zero
           coil(icoil)%zz = zero
           coil(icoil)%xt = zero
           coil(icoil)%yt = zero
           coil(icoil)%zt = zero
           coil(icoil)%xa = zero
           coil(icoil)%ya = zero
           coil(icoil)%za = zero
           NS = coil(icoil)%NS
           NCP = Splines(icoil)%NCP  ! allias variable for simplicity;

           !-------------------------enforce periodicity----------------------------------------------
           Splines(icoil)%Cpoints(NCP-3) = Splines(icoil)%Cpoints(0)
           Splines(icoil)%Cpoints(2*NCP-3) = Splines(icoil)%Cpoints(NCP) 
           Splines(icoil)%Cpoints(3*NCP-3) = Splines(icoil)%Cpoints(2*NCP)  

           Splines(icoil)%Cpoints(NCP-2) = Splines(icoil)%Cpoints(1)
           Splines(icoil)%Cpoints(2*NCP-2) = Splines(icoil)%Cpoints(NCP+1) 
           Splines(icoil)%Cpoints(3*NCP-2) = Splines(icoil)%Cpoints(2*NCP+1)  

           Splines(icoil)%Cpoints(NCP-1) = Splines(icoil)%Cpoints(2)
           Splines(icoil)%Cpoints(2*NCP-1) = Splines(icoil)%Cpoints(NCP+2) 
           Splines(icoil)%Cpoints(3*NCP-1) = Splines(icoil)%Cpoints(2*NCP+2)  
           !-------------------------calculate coil data-------------------------------------------------  
           do iseg=0,NS-1
                    coil(icoil)%xx(iseg) = sum(Splines(icoil)%Cpoints(0:NCP-1)*Splines(icoil)%basis_3(iseg,0:NCP-1))
                    coil(icoil)%yy(iseg) = sum(Splines(icoil)%Cpoints(NCP:2*NCP-1)*Splines(icoil)%basis_3(iseg,0:NCP-1))
                    coil(icoil)%zz(iseg) = sum(Splines(icoil)%Cpoints(2*NCP:3*NCP-1)*Splines(icoil)%basis_3(iseg,0:NCP-1))
                    coil(icoil)%xt(iseg) = sum(Splines(icoil)%Cpoints(0:NCP-1)*Splines(icoil)%db_dt(iseg,0:NCP-1))
                    coil(icoil)%yt(iseg) = sum(Splines(icoil)%Cpoints(NCP:2*NCP-1)*Splines(icoil)%db_dt(iseg,0:NCP-1))
                    coil(icoil)%zt(iseg) = sum(Splines(icoil)%Cpoints(2*NCP:3*NCP-1)*Splines(icoil)%db_dt(iseg,0:NCP-1))
                    coil(icoil)%xa(iseg) = sum(Splines(icoil)%Cpoints(0:NCP-1)*Splines(icoil)%db_dt_2(iseg,0:NCP-1))
                    coil(icoil)%ya(iseg) = sum(Splines(icoil)%Cpoints(NCP:2*NCP-1)*Splines(icoil)%db_dt_2(iseg,0:NCP-1))
                    coil(icoil)%za(iseg) = sum(Splines(icoil)%Cpoints(2*NCP:3*NCP-1)*Splines(icoil)%db_dt_2(iseg,0:NCP-1))
          enddo

          coil(icoil)%xx(NS) = coil(icoil)%xx(0)
          coil(icoil)%yy(NS) = coil(icoil)%yy(0)
          coil(icoil)%zz(NS) = coil(icoil)%zz(0)
          coil(icoil)%xt(NS) = coil(icoil)%xt(0)
          coil(icoil)%yt(NS) = coil(icoil)%yt(0)
          coil(icoil)%zt(NS) = coil(icoil)%zt(0)
          coil(icoil)%xa(NS) = coil(icoil)%xa(0)
          coil(icoil)%ya(NS) = coil(icoil)%ya(0)
          coil(icoil)%za(NS) = coil(icoil)%za(0)
        case default
           FATAL(discoil, .true., not supported coil types)
        end select

     endif

  enddo ! end of do icoil

  return
end subroutine discoil

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! I don't think this subroutine is used -Thomas

SUBROUTINE discfou2
  !---------------------------------------------------------------------------------------------
  ! Discretize coil data from Fourier harmonics to xx, yy, zz
  ! calling fouriermatrix for single set
  ! DATE: 2017/03/18
  !---------------------------------------------------------------------------------------------  
  use globals, only: dp, zero, pi2, myid, ncpu, ounit, coil, FouCoil, Ncoils, MPI_COMM_FOCUS
  use mpi
  implicit none

  INTEGER :: icoil, iorder, llmodnp, ierr, NS, NF
  !-------------------------call fouriermatr----------------------------------------------------  
  do icoil = 1, Ncoils

     NS = coil(icoil)%NS; NF = FouCoil(icoil)%NF  ! allias variable for simplicity;
     !reset to zero for all the coils;
     coil(icoil)%xx = zero
     coil(icoil)%yy = zero
     coil(icoil)%zz = zero
     coil(icoil)%xt = zero
     coil(icoil)%yt = zero
     coil(icoil)%zt = zero
     coil(icoil)%xa = zero
     coil(icoil)%ya = zero
     coil(icoil)%za = zero

     coil(icoil)%dd = pi2 / NS
     
     if( myid.ne.modulo(icoil-1,ncpu) ) cycle ! parallelization loop;

    iorder = 0
    call fouriermatrix( Foucoil(icoil)%xc, Foucoil(icoil)%xs, coil(icoil)%xx, NF, NS, iorder)
    call fouriermatrix( Foucoil(icoil)%yc, Foucoil(icoil)%ys, coil(icoil)%yy, NF, NS, iorder)
    call fouriermatrix( Foucoil(icoil)%zc, Foucoil(icoil)%zs, coil(icoil)%zz, NF, NS, iorder)

    iorder = 1  
    call fouriermatrix( Foucoil(icoil)%xc, Foucoil(icoil)%xs, coil(icoil)%xt, NF, NS, iorder)
    call fouriermatrix( Foucoil(icoil)%yc, Foucoil(icoil)%ys, coil(icoil)%yt, NF, NS, iorder)
    call fouriermatrix( Foucoil(icoil)%zc, Foucoil(icoil)%zs, coil(icoil)%zt, NF, NS, iorder)

    iorder = 2 
    call fouriermatrix( Foucoil(icoil)%xc, Foucoil(icoil)%xs, coil(icoil)%xa, NF, NS, iorder)
    call fouriermatrix( Foucoil(icoil)%yc, Foucoil(icoil)%ys, coil(icoil)%ya, NF, NS, iorder)
    call fouriermatrix( Foucoil(icoil)%zc, Foucoil(icoil)%zs, coil(icoil)%za, NF, NS, iorder)
   
  enddo
  !-------------------------broadcast coil data-------------------------------------------------  
  do icoil = 1, Ncoils ; llmodnp = modulo(icoil-1,ncpu)
     RlBCAST( coil(icoil)%xx(0:coil(icoil)%NS), coil(icoil)%NS+1, llmodnp )
     RlBCAST( coil(icoil)%yy(0:coil(icoil)%NS), coil(icoil)%NS+1, llmodnp )
     RlBCAST( coil(icoil)%zz(0:coil(icoil)%NS), coil(icoil)%NS+1, llmodnp )
     RlBCAST( coil(icoil)%xt(0:coil(icoil)%NS), coil(icoil)%NS+1, llmodnp )
     RlBCAST( coil(icoil)%yt(0:coil(icoil)%NS), coil(icoil)%NS+1, llmodnp )
     RlBCAST( coil(icoil)%zt(0:coil(icoil)%NS), coil(icoil)%NS+1, llmodnp )
     RlBCAST( coil(icoil)%xa(0:coil(icoil)%NS), coil(icoil)%NS+1, llmodnp )
     RlBCAST( coil(icoil)%ya(0:coil(icoil)%NS), coil(icoil)%NS+1, llmodnp )
     RlBCAST( coil(icoil)%za(0:coil(icoil)%NS), coil(icoil)%NS+1, llmodnp )
  enddo
  
  return
END SUBROUTINE discfou2

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! I don't think this subroutine is used -Thomas

subroutine fouriermatrix( xc, xs, xx, NF, ND, order )
  !---------------------------------------------------------------------------------------------
  ! This subroutine uses matrix operations to discretize data from Fourier harmonics.
  ! It's supposed to be the fastest method.
  ! DATE: 2017/03/18
  !---------------------------------------------------------------------------------------------
  use globals, only: dp, zero, pi2
  implicit none
  
  INTEGER, intent(in ) :: NF, ND, order
  REAL   , intent(in ) :: xc(0:NF), xs(0:NF)
  REAL   , intent(out) :: xx(0:ND)

  INTEGER              :: i
  REAL                 :: nn(0:NF, 1:1), tt(1:1, 0:ND), nt(0:NF, 0:ND), &
                       &  tc(1:1, 0:NF), ts(1:1, 0:NF), tx(1:1 , 0:ND), &
                       &  cnt(0:NF, 0:ND), snt(0:NF, 0:ND)

  !----------------------------data copy to matrix----------------------------------------------
  if ( size(xc) /= NF+1 ) STOP "Wrong input size for xc in subroutine fouriermatrix!"
  if ( size(xs) /= NF+1 ) STOP "Wrong input size for xs in subroutine fouriermatrix!"
  if ( size(xx) /= ND+1 ) STOP "Wrong input size for xx in subroutine fouriermatrix!"

  tc(1, 0:NF) = xc(0:NF) ! cos harmonics;
  ts(1, 0:NF) = xs(0:NF) ! sin harmonics;
  tx(1, 0:ND) = xx(0:ND) ! data coordinates;
  !----------------------------matrix assignmengt-----------------------------------------------
  nn(0:NF, 1) = (/ (i, i=0,NF) /) ! n;
  tt(1, 0:ND) = (/ (i*pi2/ND, i=0, ND) /) ! angle, t;
  nt = matmul(nn, tt)
  cnt = cos(nt)
  snt = sin(nt)
  !----------------------------select oder------------------------------------------------------
  select case (order)
  case (0)  ! 0-order

  case (1)  ! 1st-order
     do i = 0, ND
        cnt(0:NF, i) = - nn(0:NF, 1)     * snt(0:NF, i)
        snt(0:NF, i) =   nn(0:NF, 1)     * cnt(0:NF, i)
     enddo
  case (2)  ! 2nd-order
     do i = 0, ND
        cnt(0:NF, i) = -(nn(0:NF, 1)**2) * cnt(0:NF, i)
        snt(0:NF, i) = -(nn(0:NF, 1)**2) * snt(0:NF, i)
     enddo
  case default
     STOP "Invalid order in subroutine fouriermatrix"
  end select
  !----------------------------final multiplication---------------------------------------------
  tx = matmul(tc, cnt) + matmul(ts, snt)
  xx(0:ND) = tx(1, 0:ND)

  return
END subroutine fouriermatrix

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

SUBROUTINE readcoils(filename, maxnseg)
  use globals, only: dp, zero, coilsX, coilsY, coilsZ, coilsI, coilseg, coilname, Ncoils, ounit, myid, MPI_COMM_FOCUS
  use mpi
  implicit none

  INTEGER                    :: icoil, cunit, istat, astat, lstat, ierr, maxnseg, iseg
  CHARACTER*100              :: filename
  CHARACTER*200              :: line
  REAL                       :: tmp
  CHARACTER (LEN=200)        :: name

  cunit = 99
  
  ! check if file exists
  open(cunit,FILE=trim(filename),STATUS='old',IOSTAT=istat)
  if ( istat .ne. 0 ) then
     write(ounit,'("rdcoils : Error happens in reading "A, " with IOSTAT=", I6)') trim(filename), istat
     call MPI_ABORT( MPI_COMM_FOCUS, 1, ierr )
  endif
     
  ! read coils and segments data
  read(cunit,*)
  read(cunit,*)
  read(cunit,*)
  
  ! get the maximum number of segments first
  maxnseg = 0
  icoil = 0
  iseg = 0

  do
     read(cunit,'(a)', IOSTAT = istat) line
     if(istat .ne. 0 .or. line(1:3) == 'end') exit !detect EOF or end

     read(line, *, IOSTAT=lstat) tmp, tmp, tmp, tmp, tmp, name
     if (lstat .ne. 0) then
        iseg = iseg + 1
     else
        icoil = icoil + 1
        if (iseg .ge. maxnseg) maxnseg = iseg
        iseg = 0
     end if
  enddo
  close(cunit)

  ! ALLOCATE data
  Ncoils = icoil
  FATAL( readcoils, Ncoils .le. 0 , Errors in reading coils )
  FATAL( readcoils, maxnseg .le. 0 , Errors in reading coils )
  SALLOCATE( coilsX  , (1:maxnseg, 1:Ncoils), zero )
  SALLOCATE( coilsY  , (1:maxnseg, 1:Ncoils), zero )
  SALLOCATE( coilsZ  , (1:maxnseg, 1:Ncoils), zero )
  SALLOCATE( coilsI  , (1:Ncoils)           , zero )
  SALLOCATE( coilseg , (1:Ncoils)           ,    0 )
  SALLOCATE( coilname, (1:Ncoils)           , ''   )

  ! read and assign data
  open(cunit,FILE=trim(filename),STATUS='old',IOSTAT=istat)
  read(cunit,*)
  read(cunit,*)
  read(cunit,*)
  icoil = 1
  iseg = 0
  do
   read(cunit,'(a)', IOSTAT = istat) line
   if(istat .ne. 0 .or. line(1:3) == 'end') exit !detect EOF or end

   read(line, *, IOSTAT=lstat) tmp, tmp, tmp, tmp, tmp, coilname(icoil)
   if (lstat .ne. 0) then
      iseg = iseg + 1
      read(line, *) coilsX(iseg, icoil), coilsY(iseg, icoil), coilsZ(iseg, icoil), coilsI(icoil)
   else
      coilseg(icoil) = iseg
      icoil = icoil + 1
      iseg = 0
   end if
  enddo
  close(cunit)
  icoil = icoil - 1
  FATAL( readcoils, Ncoils .ne. icoil, These two should be equal)

#ifdef DEBUG
  write(ounit,'("rdcoils : Finding " I4 " coils and maximum segments number is " I6)') &
       Ncoils, maxnseg
  do icoil = 1, Ncoils
     write(ounit, '("rdcoils : Number of segments in coil " I4 " is " I6)') icoil, coilseg(icoil)
  enddo
#endif

  return

end SUBROUTINE readcoils

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

SUBROUTINE Fourier( X, XFC, XFS, Nsegs, NFcoil)
  use globals, only: dp, ounit, zero, pi2, half, myid, MPI_COMM_FOCUS
  use mpi
  implicit none

  REAL    :: X(1:Nsegs), XFC(0:NFcoil), XFS(0:NFcoil)
  INTEGER :: Nsegs, NFcoil, ifou, iseg, funit, ierr
  REAL, allocatable:: A(:), B(:)

  allocate(A(0:Nsegs-1))
  allocate(B(0:Nsegs-1))

  FATAL(Fourier, Nsegs < 2*NFcoil, Nsegs too small)
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

end SUBROUTINE Fourier

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
