
!title (coils) ! Initialize coils data using Fourier series.

!latex \briefly{Initialize the coils data with Fourier series.}

!latex \calledby{\link{FAMUS}}
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
!latex  \item[2.] \inputvar{case\_init = 0} : Read coils data from {\bf ext.FAMUS} file. The format is as following. \red{This is the most flexible way, and
!latex             each coil can be different.}            
!latex  \begin{raw}
!latex   # Total number of coils
!latex              16
!latex   #------------1--------------------------------
!latex   #coil_type     coil_name
!latex       1    Mod_001   
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
!latex   #coil_type     coil_name
!latex       2    dipole_01  
!latex   #  Lc  ox   oy   oz  Ic  I  mt  mp
!latex      1   0.0  0.0  0.0  1 1.0E6  0.0  0.0
!latex   #-----------3--backgound Bt Bz----------------
!latex   #coil_type     coil_name
!latex       3    bg_BtBz_01
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

  use famus_globals
  use ncsx_ports_eval, only: in_ncsx_port

  implicit none

  include "mpif.h"

  LOGICAL   :: exist, in_port
  INTEGER   :: icoil, maxnseg, ifirst, NF, itmp, ip, icoef, total_coef, num_pm, num_bg, & 
               num_per_array, num_tor, ipol, itor, offset, icpu, iskip
  REAL      :: Rmaj, zeta, totalcurrent, z0, r1, r2, z1, z2, rtmp, teta
  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  Nfixcur = 0 ! fixed coil current number
  Nfixgeo = 0 ! fixed coil geometry number
  num_pm  = 0 ! number of permanent magnets
  num_bg  = 0 ! number of background field

  print *,' <----famus\rdcoils entered, myid=', myid, ' ncpu=', ncpu
  if(myid == 0) write(ounit, *) "-----------INITIALIZE COILS----------------------------------"

  ! Npc = 1

  select case( case_init )

     !-------------individual coil file---------------------------------------------------------------------
  case( 0 )
     
     ! myid = 0 read all the fixed coils
     if( myid==0 ) then
        print *, '<-----famus\rdcoils myid=0 is handling fixed coils'
        if (trim(fixed_coils) /= 'none') then  !get file number;
           open( runit, file=trim(fixed_coils), status="old", action='read')
           read( runit,*)
           read( runit,*) Ncoils
           write(ounit,'("rdcoils : identified "I10" normal coils in "A" ;")') Ncoils, trim(fixed_coils)
           if (allocated(FouCoil)) deallocate(FouCoil)
           if (allocated(coil)) deallocate(coil)
           if (allocated(dof)) deallocate(dof)
           allocate( FouCoil(1:Ncoils) )
           allocate(    coil(1:Ncoils) )
           allocate(     dof(1:Ncoils) )
           !if (.not. allocated(FouCoil)) allocate( FouCoil(1:Ncoils) )
           !if (.not. allocated(coil)) allocate(    coil(1:Ncoils) )
           !if (.not. allocated(dof)) allocate(     dof(1:Ncoils) )
           do icoil = 1, Ncoils
              read( runit,*)
              read( runit,*)
              read( runit,*) coil(icoil)%itype, coil(icoil)%symmetry, coil(icoil)%name
              if(coil(icoil)%itype == 1) then  ! Fourier representation
                 read( runit,*)
                 read( runit,*) coil(icoil)%NS, coil(icoil)%I, coil(icoil)%Ic, &
                      & coil(icoil)%L, coil(icoil)%Lc, coil(icoil)%Lo
                 FATAL( rdcoils, coil(icoil)%NS < 0                        , illegal )
                 FATAL( rdcoils, coil(icoil)%Ic < 0 .or. coil(icoil)%Ic > 1, illegal )
                 FATAL( rdcoils, coil(icoil)%Lc < 0 .or. coil(icoil)%Lc > 2, illegal )
                 FATAL( rdcoils, coil(icoil)%L  < zero                     , illegal )
                 FATAL( rdcoils, coil(icoil)%Lc < zero                     , illegal )
                 FATAL( rdcoils, coil(icoil)%Lo < zero                     , illegal )
                 read( runit,*)
                 read( runit,*) FouCoil(icoil)%NF
                 FATAL( rdcoils, Foucoil(icoil)%NF  < 0                    , illegal )
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
              else if (coil(icoil)%itype == 2) then  ! permanent magnets
                 read( runit,*)
                 read( runit,*) coil(icoil)%Lc, coil(icoil)%ox, coil(icoil)%oy, coil(icoil)%oz, &
                      coil(icoil)%Ic, coil(icoil)%I , coil(icoil)%mt, coil(icoil)%mp      
              else if (coil(icoil)%itype == 3) then  ! backgroud toroidal/vertical field
                 read( runit,*)
                 read( runit,*) coil(icoil)%Ic, coil(icoil)%I, coil(icoil)%Lc, coil(icoil)%Bz
              else
                 STOP " wrong coil type in rdcoils"
                 call MPI_ABORT(MPI_COMM_FAMUS, 1, ierr)
              endif

           enddo !end do icoil;
           close( runit )
        else ! initialize one fake coil;
           print *, '<-----rdcoils myid=0 is initializing one fake coil'
           Ncoils = 1
           if (allocated(coil)) deallocate(coil)
           if (allocated(dof)) deallocate(dof)
           allocate( coil(1:Ncoils) )
           allocate( dof(1:Ncoils) )
           !if (.not. allocated(coil)) allocate(    coil(1:Ncoils) )
           !if (.not. allocated(dof)) allocate(     dof(1:Ncoils) )
           coil(1)%itype = 3
           coil(1)%Ic = 0
           coil(1)%Lc = 0
           coil(1)%I  = zero
           coil(1)%Bz = zero
        endif
     endif

     print *, '<----famus\rdcoils is moving on to handle dipoles.'

     if( myid==0 ) then  !get file number;
        open( runit, file=trim(input_coils), status="old", action='read')
        read( runit,*)
        read( runit,*) Ncoils_total, momentq
        write(ounit,'("rdcoils : identified "I10" unique dipoles in "A" ;")') Ncoils_total, trim(input_coils)
        write(ounit,'("        : penalization coefficient q = "I2" ;")') momentq
        close( runit )
     endif

     print *, '<----famus\rdcoils is broadcasting Ncoils_total and momenq, myid=',myid
     IlBCAST( Ncoils_total, 1,  0 )
     IlBCAST( momentq, 1,  0 )
     if (Ncpu == 1) STOP ' At least use two cpus'
     if (myid/=0) Ncoils = Ncoils_total / (Ncpu-1)
     if (myid==Ncpu-1) Ncoils = Ncoils + mod(Ncoils_total, Ncpu-1) ! The last one read more
     !MPIOUT( Ncoils )
#ifdef TOPO
     offset = 3+(myid-1)*(Ncoils_total/(Ncpu-1))*1 ! only for dipoles
#else
     offset = 2+(myid-1)*(Ncoils_total/(Ncpu-1))*5 ! only for dipoles
#endif
     !MPIOUT( offset )
     print *, '<----famus myid=',myid,' calculated offset = ', offset
                               
     print *, '<----famus\rdcoils is allocaing coil & dof, myid=',myid
     if (myid /= 0) then
        if (allocated(coil)) deallocate(coil) ! necessary???
        allocate( coil(1:Ncoils) )
        if (allocated(dof)) deallocate(dof) ! necessary???
        allocate(  dof(1:Ncoils) )
        !if (.not. allocated(coil)) allocate(coil(1:Ncoils)) ! necessary???
        !if (.not. allocated(dof)) allocate(dof(1:Ncoils)) ! necessary???
        ! coil(1:Ncoils)%symmetry = 1
     endif

     print *, '<----famus\rdcoils is taking turns reading input_coils, myid=',myid,' offset=',offset, ' ncoils=', ncoils
     do icpu = 1, ncpu-1
        if (myid == icpu) then                              ! each cpu read the coils at the same time.
           open( runit, file=trim(input_coils), status="old", action='read')

           do iskip = 1, offset
              read( runit,*)
           enddo

           do icoil = 1, Ncoils
#ifdef TOPO
              read( runit,*) coil(icoil)%itype, coil(icoil)%symmetry, coil(icoil)%name, &
                   coil(icoil)%ox, coil(icoil)%oy, coil(icoil)%oz, &
                   coil(icoil)%Ic, coil(icoil)%moment, coil(icoil)%pho, &
                   coil(icoil)%Lc, coil(icoil)%mp, coil(icoil)%mt
              call in_ncsx_port(coil(icoil)%ox, coil(icoil)%oy, coil(icoil)%oz,&
                                in_port)
              if (in_port) then
                  coil(icoil)%Ic = 0
                  coil(icoil)%Lc = 0
                  coil(icoil)%pho = 0
              end if
#else
              read( runit,*)
              read( runit,*) coil(icoil)%itype, coil(icoil)%symmetry, coil(icoil)%name
              read( runit,*)
              read( runit,*) coil(icoil)%Lc, coil(icoil)%ox, coil(icoil)%oy, coil(icoil)%oz, &
                   coil(icoil)%Ic, coil(icoil)%I , coil(icoil)%mt, coil(icoil)%mp      
#endif
              !coil(icoil)%I = coil(icoil)%moment*sin(coil(icoil)%pho)**momentq
              coil(icoil)%I = coil(icoil)%moment*(coil(icoil)%pho)**momentq
           enddo !end do icoil;

           close(runit)
        endif ! end of if( myid == 0 )
     enddo

  end select

  !-----------------------allocate   coil data--------------------------------------------------
  print *, '<----famus\rdcoils allocating coil data, myid=',myid
!  if (myid == 0) then
     do icoil = 1, Ncoils
        SALLOCATE( coil(icoil)%xx, (0:coil(icoil)%NS), zero )
        SALLOCATE( coil(icoil)%yy, (0:coil(icoil)%NS), zero )
        SALLOCATE( coil(icoil)%zz, (0:coil(icoil)%NS), zero )
        SALLOCATE( coil(icoil)%xt, (0:coil(icoil)%NS), zero )
        SALLOCATE( coil(icoil)%yt, (0:coil(icoil)%NS), zero )
        SALLOCATE( coil(icoil)%zt, (0:coil(icoil)%NS), zero )
        SALLOCATE( coil(icoil)%xa, (0:coil(icoil)%NS), zero )
        SALLOCATE( coil(icoil)%ya, (0:coil(icoil)%NS), zero )
        SALLOCATE( coil(icoil)%za, (0:coil(icoil)%NS), zero )
        SALLOCATE( coil(icoil)%dl, (0:coil(icoil)%NS), zero )
        SALLOCATE( coil(icoil)%dd, (0:coil(icoil)%NS), zero )
     enddo
!  end if

  !-----------------------allocate DoF arrays --------------------------------------------------  
  print *, '<----famus\rdcoils allocating DoF arrays, myid=',myid

  itmp = -1
  call AllocData(itmp)

  !-----------------------discretize coil data--------------------------------------------------

  if (myid == 0) then
     if (IsQuiet < 0) write(ounit, '(8X,": coils will be discretized in "I6" segments")') Nseg
  endif

  ifirst = 1
  
  print *, '<----famus\rdcoils calling discoil, myid=',myid
  call discoil(ifirst)
  ifirst = 0

  print *, '<----famus\rdcoils calling mpi_barrier, myid=',myid
  call MPI_BARRIER( MPI_COMM_FAMUS, ierr )
  print *, '<----famus\rdcoils is exiting, myid=',myid

  return

end subroutine rdcoils

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine famus_rdcoils_mpi(coil_filename, Ncoils_total, momentq)

  use famus_globals, ONLY: MPI_COMM_FAMUS, arbitrarycoil, &
            coil_out => export_coil
  use ncsx_ports_eval, only: in_ncsx_port

  implicit none

  include "mpif.h"

  CHARACTER(100), INTENT(IN) :: coil_filename
  INTEGER, INTENT(OUT)       :: Ncoils_total
  INTEGER, INTENT(OUT)       :: momentq
!  type(arbitrarycoil) dimension(:), allocatable, INTENT(OUT) :: coil_out
  INTEGER   :: icoil, Ncoils,  & 
               myid, ncpu, runit, &
               ipol, offset, icpu, iskip, ierr
  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

   call MPI_COMM_RANK( MPI_COMM_FAMUS, myid, ierr )
   call MPI_COMM_SIZE( MPI_COMM_FAMUS, ncpu, ierr )

print *,'<----famus_rdcoils_mpi myid/rank=',myid,' ncpu/size=',ncpu

!  if(myid == 0) write(6, *) "-----------READING FAMUS COILS (MPI)  -----------------------------"
!     print *, '<----myid=', myid, ' has entered famus_rdcoils_mpi'
     !-------------individual coil file---------------------------------------------------------------------
     
     if( myid==0 ) then  !get file number;
        open( runit, file=trim(coil_filename), status="old", action='read')
        read( runit,*)
        read( runit,*) Ncoils_total, momentq
        write(6,'("famus_rdcoils_mpi : identified "I10" unique dipoles in "A" ;")') Ncoils_total, trim(coil_filename)
        write(6,'("        : penalization coefficient q = "I2" ;")') momentq
        close( runit )
     endif

     IlBCAST( Ncoils_total, 1,  0 )
     IlBCAST( momentq, 1,  0 )
     if (ncpu == 1) STOP ' At least use two cpus'
     Ncoils = Ncoils_total

!print *,'<----famus_rdcoils_mpi post broadcast myid/rank=',myid,' ncpu/size=',ncpu

     offset = 3                       
     if ( allocated(coil_out) ) then
       deallocate(coil_out)
     endif
     allocate( coil_out(Ncoils) )
     !if ( .not. allocated(coil_out) ) then
     !  allocate( coil_out(Ncoils) )
     !endif
!print *,'<----famus_rdcoils_mpi post allocate myid/rank=',myid,' ncpu/size=',ncpu

     do icpu = 0, ncpu-1
        if (myid == icpu) then                              ! each cpu read the coils at the same time.
           open( runit, file=trim(coil_filename), status="old", action='read')

           do iskip = 1, offset
              read( runit,*)
           enddo

           do icoil = 1, Ncoils
              read( runit,*) coil_out(icoil)%itype, coil_out(icoil)%symmetry, coil_out(icoil)%name, &
              coil_out(icoil)%ox, coil_out(icoil)%oy, coil_out(icoil)%oz, &
              coil_out(icoil)%Ic, coil_out(icoil)%moment , &
              coil_out(icoil)%pho, coil_out(icoil)%Lc, coil_out(icoil)%mp, coil_out(icoil)%mt      
              coil_out(icoil)%I = coil_out(icoil)%moment*(coil_out(icoil)%pho)**momentq
!              print *, '<----coilname: ',trim(coil_out(icoil)%name)
           enddo !end do icoil;

           close(runit)
        endif ! end of if( myid == 0 )
     enddo


  !-----------------------allocate   coil data--------------------------------------------------
!  print *, '<----checkup: coil_out(1)%ox=',coil_out(1)%ox
  print *, '<----famus_rdcoils_mpi is calling mpi_barrier prior to exit'
  call MPI_BARRIER( MPI_COMM_FAMUS, ierr )
  print *, '<----famus_rdcoils_mpi is exiting'

  return

end subroutine famus_rdcoils_mpi

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine discoil(ifirst)
!---------------------------------------------------------------------------------------------
! dicretize coils data;
! if ifirst = 1, it will update all the coils; otherwise, only update free coils;
! date: 20170314
!---------------------------------------------------------------------------------------------
  use famus_globals, only: dp, zero, pi2, myid, ounit, coil, FouCoil, Ncoils, DoF, MPI_COMM_FAMUS
  implicit none
  include "mpif.h"

  INTEGER, intent(in) :: ifirst

  INTEGER          :: icoil, iseg, mm, NS, NF, ierr, astat, ip
  REAL             :: tt
  REAL,allocatable :: cmt(:,:), smt(:,:)
  !-------------------------------------------------------------------------------------------

  do icoil = 1, Ncoils

     if( (coil(icoil)%Lc + coil(icoil)%Ic + ifirst) /= 0) then  !first time or if Lc/=0, then need discretize;

        !if( myid.ne.modulo(icoil-1,ncpu) ) cycle ! parallelization loop;

        select case (coil(icoil)%itype)
        case( 1 )

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

           NS = coil(icoil)%NS; NF = FouCoil(icoil)%NF  ! allias variable for simplicity;
           SALLOCATE( cmt, (0:NS, 0:NF), zero )
           SALLOCATE( smt, (0:NS, 0:NF), zero )

           do iseg = 0, NS ; tt = iseg * pi2 / NS
              do mm = 0, NF
                 cmt(iseg,mm) = cos( mm * tt )
                 smt(iseg,mm) = sin( mm * tt )
              enddo
           enddo

           !-------------------------calculate coil data-------------------------------------------------  
           mm = 0
           coil(icoil)%xx(0:NS) = cmt(0:NS,mm) * Foucoil(icoil)%xc(mm)
           coil(icoil)%yy(0:NS) = cmt(0:NS,mm) * Foucoil(icoil)%yc(mm)
           coil(icoil)%zz(0:NS) = cmt(0:NS,mm) * Foucoil(icoil)%zc(mm)    
           do mm = 1, NF   
              coil(icoil)%xx(0:NS) = coil(icoil)%xx(0:NS) + (   cmt(0:NS,mm) * Foucoil(icoil)%xc(mm) &
                                                              + smt(0:NS,mm) * Foucoil(icoil)%xs(mm) )
              coil(icoil)%yy(0:NS) = coil(icoil)%yy(0:NS) + (   cmt(0:NS,mm) * Foucoil(icoil)%yc(mm) &
                                                              + smt(0:NS,mm) * Foucoil(icoil)%ys(mm) )
              coil(icoil)%zz(0:NS) = coil(icoil)%zz(0:NS) + (   cmt(0:NS,mm) * Foucoil(icoil)%zc(mm) &
                                                              + smt(0:NS,mm) * Foucoil(icoil)%zs(mm) )

              coil(icoil)%xt(0:NS) = coil(icoil)%xt(0:NS) + ( - smt(0:NS,mm) * Foucoil(icoil)%xc(mm) &
                                                              + cmt(0:NS,mm) * Foucoil(icoil)%xs(mm) ) * mm
              coil(icoil)%yt(0:NS) = coil(icoil)%yt(0:NS) + ( - smt(0:NS,mm) * Foucoil(icoil)%yc(mm) &
                                                              + cmt(0:NS,mm) * Foucoil(icoil)%ys(mm) ) * mm
              coil(icoil)%zt(0:NS) = coil(icoil)%zt(0:NS) + ( - smt(0:NS,mm) * Foucoil(icoil)%zc(mm) &
                                                              + cmt(0:NS,mm) * Foucoil(icoil)%zs(mm) ) * mm

              coil(icoil)%xa(0:NS) = coil(icoil)%xa(0:NS) + ( - cmt(0:NS,mm) * Foucoil(icoil)%xc(mm) &
                                                              - smt(0:NS,mm) * Foucoil(icoil)%xs(mm) ) * mm*mm
              coil(icoil)%ya(0:NS) = coil(icoil)%ya(0:NS) + ( - cmt(0:NS,mm) * Foucoil(icoil)%yc(mm) &
                                                              - smt(0:NS,mm) * Foucoil(icoil)%ys(mm) ) * mm*mm
              coil(icoil)%za(0:NS) = coil(icoil)%za(0:NS) + ( - cmt(0:NS,mm) * Foucoil(icoil)%zc(mm) &
                                                              - smt(0:NS,mm) * Foucoil(icoil)%zs(mm) ) * mm*mm
           enddo ! end of do mm; 
           if(ifirst /= 0) then
              DoF(icoil)%xof(0:NS-1,      1:  NF+1) = cmt(0:NS-1, 0:NF)  !xc
              DoF(icoil)%xof(0:NS-1,   NF+2:2*NF+1) = smt(0:NS-1, 1:NF)  !xs
              DoF(icoil)%yof(0:NS-1, 2*NF+2:3*NF+2) = cmt(0:NS-1, 0:NF)  !yc
              DoF(icoil)%yof(0:NS-1, 3*NF+3:4*NF+2) = smt(0:NS-1, 1:NF)  !ys
              DoF(icoil)%zof(0:NS-1, 4*NF+3:5*NF+3) = cmt(0:NS-1, 0:NF)  !zc
              DoF(icoil)%zof(0:NS-1, 5*NF+4:6*NF+3) = smt(0:NS-1, 1:NF)  !zs
           endif
           coil(icoil)%dd = pi2 / NS  ! discretizing factor;

           DALLOCATE(cmt)
           DALLOCATE(smt)

        case(2)
           
           coil(icoil)%mx = sin(coil(icoil)%mt) * cos(coil(icoil)%mp) * coil(icoil)%I
           coil(icoil)%my = sin(coil(icoil)%mt) * sin(coil(icoil)%mp) * coil(icoil)%I
           coil(icoil)%mz = cos(coil(icoil)%mt) * coil(icoil)%I

        case(3)

        case default
           FATAL(discoil, .true., not supported coil types)
        end select

     endif

  enddo ! end of do icoil

  return
end subroutine discoil

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine fouriermatrix( xc, xs, xx, NF, ND, order )
  !---------------------------------------------------------------------------------------------
  ! This subroutine uses matrix operations to discretize data from Fourier harmonics.
  ! It's supposed to be the fastest method.
  ! DATE: 2017/03/18
  !---------------------------------------------------------------------------------------------
  use famus_globals, only: dp, zero, pi2
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
  use famus_globals, only: dp, zero, coilsX, coilsY, coilsZ, coilsI, coilseg, coilname, Ncoils, ounit, myid, &
                     MPI_COMM_FAMUS
  implicit none
  include "mpif.h"

  INTEGER, parameter         :: mcoil = 256, mseg = 1024 ! Largest coils and segments number
  INTEGER                    :: icoil, cunit, istat, astat, lstat, ierr, maxnseg, seg(1:mseg)
  REAL, dimension(mseg,mcoil):: x, y, z, I
  CHARACTER*100              :: filename
  CHARACTER*200              :: line
  REAL                       :: tmp
  CHARACTER (LEN=20), dimension(mcoil) :: name

  cunit = 99; I = 1.0; Ncoils= 1; maxnseg = 0; seg = 0;

  open(cunit,FILE=trim(filename),STATUS='old',IOSTAT=istat)
  if ( istat .ne. 0 ) then
     write(ounit,'("rdcoils : Reading coils data error in "A)') trim(filename)
     call MPI_ABORT( MPI_COMM_FAMUS, 1, ierr )
  endif
     
  ! read coils and segments data
  read(cunit,*)
  read(cunit,*)
  read(cunit,*)

  do
     read(cunit,'(a)', IOSTAT = istat) line
     if(istat .ne. 0 .or. line(1:3) == 'end') exit !detect EOF or end

     seg(Ncoils) = seg(Ncoils) + 1
     read(line,*, IOSTAT = lstat) x(seg(Ncoils), Ncoils), y(seg(Ncoils), Ncoils), z(seg(Ncoils), Ncoils), &
          I(seg(Ncoils), Ncoils)
     if ( I(seg(Ncoils), Ncoils) == zero) then
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

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

SUBROUTINE Fourier( X, XFC, XFS, Nsegs, NFcoil)
  use famus_globals, only: dp, ounit, zero, pi2, half, myid, MPI_COMM_FAMUS
  implicit none
  include "mpif.h"

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
