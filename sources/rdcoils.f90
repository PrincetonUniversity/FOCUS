
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
               num_per_array, num_tor, ipol, itor, ns
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
        ns = coilseg(icoil)  
        coil(icoil)%symm = 0 ! no symmetry or periodicity
        coil(icoil)%NS =  ns
        coil(icoil)%I  =  coilsI(icoil)
        coil(icoil)%Ic =  IsVaryCurrent
        coil(icoil)%L  =  target_length ! irrelevant until re-computed;
        !coil(icoil)%k0 =  k0
        coil(icoil)%Lc =  IsVaryGeometry
        coil(icoil)%Lo =  target_length
        coil(icoil)%name = trim(coilname(icoil))
        SALLOCATE( coil(icoil)%xx, (0:coil(icoil)%NS), zero )
        SALLOCATE( coil(icoil)%yy, (0:coil(icoil)%NS), zero )
        SALLOCATE( coil(icoil)%zz, (0:coil(icoil)%NS), zero )     
        ! assignment
        coil(icoil)%xx(0:ns-1) = coilsX(1:ns, icoil); coil(icoil)%xx(ns) = coil(icoil)%xx(0)
        coil(icoil)%yy(0:ns-1) = coilsY(1:ns, icoil); coil(icoil)%yy(ns) = coil(icoil)%yy(0)
        coil(icoil)%zz(0:ns-1) = coilsZ(1:ns, icoil); coil(icoil)%zz(ns) = coil(icoil)%zz(0)
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
   end select  

  return

end subroutine rdcoils

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

SUBROUTINE readcoils(filename, maxnseg)
  use globals, only: dp, zero, coilsX, coilsY, coilsZ, coilsI, coilseg, coilname, Ncoils, ounit, myid, MPI_COMM_FOCUS
  implicit none
  include "mpif.h"

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

end SUBROUTINE READCOILS

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
