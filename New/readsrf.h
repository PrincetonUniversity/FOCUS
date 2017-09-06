!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

!title (main) ! Main program.

!latex \briefly{Main program.}

!latex \calledby{\link{}}
!latex \calls{\link{initial}, \link{rdsurf}, \link{rdcoils}, 
!latex \link{saving}, \link{solvers}}.

!latex \tableofcontents

!latex  \subsection{Brief introduction}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

subroutine readsrf

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  use globals, only : zero, half, pi2, sqrtmachprec, myid, ounit, &
                      Nteta, Nzeta, surf, case_surface, IsSymmetric, discretefactor, knotsurf
  
  implicit none
  
  include "mpif.h"
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  INTEGER :: astat, ierr, ii, jj, iteta, jzeta
  REAL    :: teta, zeta, xx(1:3), xt(1:3), xz(1:3), bn, ax(1:3), at(1:3), az(1:3), ds(1:3), dd

! INTEGER :: secs, mins, hrs, ii, jj
! REAL    :: teta, zeta, xx(1:3), xt(1:3), xz(1:3), bn, ax(1:3), at(1:3), az(1:3), ds(1:3), dd
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  select case( case_surface )
   
  case( 0 )
   
   write(ounit,'("readsrf : " 10x " : case_surface ="i2" ; Nteta ="i4", Nzeta ="i4" ;")') case_surface, Nteta, Nzeta
   
   call fousurf  ! VMEC  -style plasma boundary;
   
  case( 1 )

   write(ounit,'("readsrf : " 10x " : case_surface ="i2" ; Nteta ="i4", Nzeta ="i4" ; knotsurf ="f6.3" ;")') case_surface, Nteta, Nzeta, knotsurf

   call rdknot   ! axis  -style plasma boundary;

 !case( 2 )

  !call readwout ! Boozer-style plasma boundary (under construction);

  case default

   FATAL( focus , .true., selected value of case_surface not supported )

  end select

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! select case( IsSymmetric )
! case ( 0 ) ; Nfp = 1
! case ( 1 ) ;         ; Npc = 1
! case ( 2 ) ;         ; Npc = Nfp
! case default
!  FATAL( focus, .true., selected value of IsSymmetric not supported )
! end select
  
! discretefactor = discretefactor / Nfp ! discretefactor was initialized in initial; SRH; 28 Sep 17;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  allocate( surf(1:1) )
  
  surf(1)%Nteta = Nteta
  surf(1)%Nzeta = Nzeta
  
  SALLOCATE( surf(1)%xx, (0:Nteta  ,0:Nzeta  ), zero )
  SALLOCATE( surf(1)%yy, (0:Nteta  ,0:Nzeta  ), zero )
  SALLOCATE( surf(1)%zz, (0:Nteta  ,0:Nzeta  ), zero )
  
  SALLOCATE( surf(1)%nx, (0:Nteta  ,0:Nzeta  ), zero )
  SALLOCATE( surf(1)%ny, (0:Nteta  ,0:Nzeta  ), zero )
  SALLOCATE( surf(1)%nz, (0:Nteta  ,0:Nzeta  ), zero )
  
  SALLOCATE( surf(1)%ds, (0:Nteta  ,0:Nzeta  ), zero )
  
  SALLOCATE( surf(1)%xt, (0:Nteta  ,0:Nzeta  ), zero )
  SALLOCATE( surf(1)%yt, (0:Nteta  ,0:Nzeta  ), zero )
  SALLOCATE( surf(1)%zt, (0:Nteta  ,0:Nzeta  ), zero )
  
  SALLOCATE( surf(1)%Bx, (0:Nteta-1,0:Nzeta-1), zero )
  SALLOCATE( surf(1)%By, (0:Nteta-1,0:Nzeta-1), zero )
  SALLOCATE( surf(1)%Bz, (0:Nteta-1,0:Nzeta-1), zero )
  
  SALLOCATE( surf(1)%Bn, (0:Nteta-1,0:Nzeta-1), zero )

  SALLOCATE( surf(1)%pb, (0:Nteta  ,0:Nzeta  ), zero ) ! normal field; SRH; 28 Sep 17;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  CHECK( readsrf, Nteta.eq.0, error )
  CHECK( readsrf, Nzeta.eq.0, error )

  do iteta = 0, Nteta ; teta = ( iteta + half ) * pi2 / Nteta ! grid center;
   
   do jzeta = 0, Nzeta ; zeta = ( jzeta + half ) * pi2 / Nzeta
    
    select case( case_surface )
    case( 0 ) ; call surfxx(           teta, zeta,                            xx(1:3), xt(1:3), xz(1:3), bn )
    case( 1 ) ; call knotxx( knotsurf, teta, zeta, ax(1:3), at(1:3), az(1:3), xx(1:3), xt(1:3), xz(1:3)     ) ; bn = zero
    case default
     FATAL( focus , .true., selected case_surface not supported )
    end select
    
    ds(1:3) = - (/ xt(2) * xz(3) - xt(3) * xz(2), xt(3) * xz(1) - xt(1) * xz(3), xt(1) * xz(2) - xt(2) * xz(1) /) ! negative sign = counterclockwise;
    
    dd = sqrt( sum( ds(1:3)*ds(1:3) ) )

    surf(1)%xx(iteta,jzeta) = xx(1)
    surf(1)%yy(iteta,jzeta) = xx(2)
    surf(1)%zz(iteta,jzeta) = xx(3)
    
    surf(1)%xt(iteta,jzeta) = xt(1)
    surf(1)%yt(iteta,jzeta) = xt(2)
    surf(1)%zt(iteta,jzeta) = xt(3)

    CHECK( readsrf, dd.lt.sqrtmachprec, divide by zero )
    
    surf(1)%nx(iteta,jzeta) = ds(1) / dd
    surf(1)%ny(iteta,jzeta) = ds(2) / dd
    surf(1)%nz(iteta,jzeta) = ds(3) / dd
    
    surf(1)%ds(iteta,jzeta) =         dd
    
    surf(1)%pb(iteta,jzeta) = bn ! target normal field; SRH; 28 Sep 17;
    
   enddo ! end of do jzeta;
   
  enddo ! end of do iteta;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  return
  
end subroutine readsrf

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!



!title (boundary) ! The plasma boundary is read from file 

!latex \briefly{A Fourier representation for the plasma boundary is read from file }

!latex \calledby{\link{xfocus}}
!latex \calls{\link{}}

!latex \section{General representation (stellarator)}
!latex \subsection{overview}
!latex The general representation for plasma boundary is in \subroutine{generic}. The basic fomulation
!latex is
!latex \be
!latex \ds R &= \sum R_{mn}^c \, \cos(m\t - n\z) + R_{mn}^s \, \sin(m\t - n\z) \nonumber \\
!latex \ds Z &= \sum Z_{mn}^c \, \cos(m\t - n\z) + Z_{mn}^s \, \sin(m\t - n\z) \nonumber
!latex \ee
!latex Usually, if we imply stellarator symmetry, then $R_{mn}^s$ and $Z_{mn}^c$ would be zero.
!latex 
!latex The positive driection for poloidal angle $\t$ is \red{counterclockwise} and for toroidal angle is also
!latex \red{counterclockwise} from the top view. The positive surface normal should be pointed outwards.
!latex \subsection{Variables}
!latex The Fourier harmonics of the plasma boundary are  reqired in \verb+plasma.boundary+, 
!latex and the format of this file is as follows: 
!latex \begin{raw}
!latex Nfou       ! integer: number of Fourier harmonics for the plasma boundary;
!latex Nfp        ! integer: number of field periodicity;
!latex NBnf       ! integer: number of Fourier harmonics for Bn;
!latex ---------------------------------------------------------
!latex bim(1:bmn) ! integer: poloidal mode identification;
!latex bin(1:bmn) ! integer: toroidal mode identification;
!latex Bnim(1:bmn)! integer: poloidal mode identification, for Bn;
!latex Bnin(1:bmn)! integer: toroidal mode identification, for Bn;
!latex ---------------------------------------------------------
!latex Rbc(1:bmn) ! real   : cylindrical R cosine harmonics;
!latex Rbs(1:bmn) ! real   : cylindrical R   sine harmonics;
!latex Zbc(1:bmn) ! real   : cylindrical Z cosine harmonics;
!latex Zbs(1:bmn) ! real   : cylindrical Z   sine harmonics;
!latex Bns(1:nbf) ! real   : B normal sin harmonics;
!latex Bnc(1:nbf) ! real   : B normal cos harmonics; \end{raw}
!latex Note that immediately after reading (and broadcasting) 
!latex \verb+bin+, the field periodicity factor is included, i.e. \verb+bin = bin * Nfp+.
!latex \subsection{Sample file}
!latex Example of the plasma.boundary file:
!latex  { \begin{raw}
!latex #Nfou Nfp NBnf
!latex 4 2 1
!latex #plasma boundary
!latex # n m Rbc Rbs Zbc Zbs
!latex 0 0  3.00 0.0 0.0  0.00
!latex 0 1  0.30 0.0 0.0 -0.30
!latex 1 0  0.00 0.0 0.0 -0.06
!latex 1 1 -0.06 0.0 0.0 -0.06
!latex #Bn harmonics
!latex # n m bnc bns
!latex 0 0 0.0 0.0
!latex \end{raw}
!latex }
!latex \section{Knotran}
!latex The input surface file for knotrans is descriped in \code{knotxx}.
!latex \section{Tokamak}
!latex This part is reserved for later development of the interface for tokamaks.

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine fousurf
  
  use globals, only : zero, half, pi2, myid, ounit, runit, surffile, IsQuiet, IsSymmetric, &
                      Nfp, &
                      Nfou, NBnf, bim, bin, Bnim, Bnin, Rbc, Rbs, Zbc, Zbs, Bnc, Bns,  &
                      Nteta, Nzeta, surf, discretefactor
  
  implicit none
  
  include "mpif.h"
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOGICAL :: exist
  INTEGER :: iosta, astat, ierr, ii, jj, imn
  REAL    :: RR(0:2), ZZ(0:2), szeta, czeta, xx(1:3), xt(1:3), xz(1:3), ds(1:3), &
             teta, zeta, arg, dd
  
  !-------------read plasma.boundary---------------------------------------------------------------------  
  inquire( file=trim(surffile), exist=exist)  
  FATAL( surface, .not.exist, plasma.boundary does not exist ) 
  if( myid == 0 ) then
     open(runit, file=trim(surffile), status='old', action='read')
     read(runit,*) !empty line
     read(runit,*) Nfou, Nfp, NBnf !read dimensions
  endif
  
  !Broadcast the values
  IlBCAST( Nfou , 1, 0 )
  IlBCAST( Nfp  , 1, 0 )
  IlBCAST( NBnf , 1, 0 )  
  FATAL( surface, Nfou <= 0, invalid )
  FATAL( surface, Nfp  <= 0, invalid )
  FATAL( surface, NBnf <  0, invalid )

  !allocate arrays
  SALLOCATE( bim, (1:Nfou), 0 )
  SALLOCATE( bin, (1:Nfou), 0 )  
  SALLOCATE( Rbc, (1:Nfou), zero )
  SALLOCATE( Rbs, (1:Nfou), zero )
  SALLOCATE( Zbc, (1:Nfou), zero )
  SALLOCATE( Zbs, (1:Nfou), zero )

  if( myid == 0 ) then
   read(runit,*) !empty line
   read(runit,*) !empty line
   do imn = 1, Nfou
      read(runit,*) bin(imn), bim(imn), Rbc(imn), Rbs(imn), Zbc(imn), Zbs(imn)
   enddo
  endif  

  IlBCAST( bim(1:Nfou), Nfou, 0 )
  IlBCAST( bin(1:Nfou), Nfou, 0 )
 
  bin(1:Nfou) = bin(1:Nfou) * Nfp  !The full plasma;
     
  RlBCAST( Rbc(1:Nfou), Nfou, 0 )
  RlBCAST( Rbs(1:Nfou), Nfou, 0 )
  RlBCAST( Zbc(1:Nfou), Nfou, 0 )
  RlBCAST( Zbs(1:Nfou), Nfou, 0 )

  !read Bnormal ditributions
  if( NBnf  >   0) then
     SALLOCATE( Bnim, (1:NBnf), 0    )
     SALLOCATE( Bnin, (1:NBnf), 0    )
     SALLOCATE( Bnc , (1:NBnf), zero )
     SALLOCATE( Bns , (1:NBnf), zero )

     if( myid == 0 ) then
        read(runit,*) !empty line
        read(runit,*) !empty line
        do imn = 1, NBnf 
           read(runit,*) Bnin(imn), Bnim(imn), Bnc(imn), Bns(imn)
        enddo
     endif

     IlBCAST( Bnim(1:NBnf), NBnf, 0 )
     IlBCAST( Bnin(1:NBnf), NBnf, 0 )

     !if (IsSymmetric  ==  0)
     Bnin(1:NBnf) = Bnin(1:NBnf) * Nfp ! Disarde periodicity;
     ! This should be consistent with bnftran; Before fully constructed the stellarator symmetry,
     ! it's turned off;
     
     RlBCAST( Bnc(1:NBnf) , NBnf, 0 )
     RlBCAST( Bns(1:NBnf) , NBnf, 0 )
  endif

  if( myid == 0 ) close(runit,iostat=iosta)
  
  IlBCAST( iosta, 1, 0 )
  
  FATAL( surface, iosta.ne.0, error closing plasma.boundary )
  
  !-------------output for check-------------------------------------------------------------------------
  if( myid == 0 .and. IsQuiet <= 0) then
     write(ounit, *) "-----------Reading surface-----------------------------------"
     write(ounit, '("surface : Nfou = " I06 " ; Nfp = " I06 " ; NBnf = " I06 " ;" )') Nfou, Nfp, NBnf
  endif

  if( myid == 0 .and. IsQuiet <= -2) then !very detailed output;
     write(ounit,'("        : " 10x " : bim ="10i13   )') bim(1:Nfou)
     write(ounit,'("        : " 10x " : bin ="10i13   )') bin(1:Nfou)
     write(ounit,'("        : " 10x " : Rbc ="10es13.5)') Rbc(1:Nfou)
     write(ounit,'("        : " 10x " : Rbs ="10es13.5)') Rbs(1:Nfou)
     write(ounit,'("        : " 10x " : Zbc ="10es13.5)') Zbc(1:Nfou)
     write(ounit,'("        : " 10x " : Zbs ="10es13.5)') Zbs(1:Nfou)
     if(Nbnf > 0) then
        write(ounit,'("        : " 10x " : Bnim ="10i13  )') Bnim(1:NBnf)
        write(ounit,'("        : " 10x " : Bnin ="10i13  )') Bnin(1:NBnf)
        write(ounit,'("        : " 10x " : Bnc ="10es13.5)') Bnc (1:NBnf)
        write(ounit,'("        : " 10x " : Bns ="10es13.5)') Bns (1:NBnf)
     endif
  endif

  !-------------discretize surface data------------------------------------------------------------------  

! select case (IsSymmetric)
! case ( 0 )
!    Nfp = 1                          !reset Nfp to 1;
!    Npc = 1                          !number of coils periodicity 
! case ( 1 )                          !plasma periodicity enabled;
!    Npc = 1
! case ( 2 )                          !plasma and coil periodicity enabled;
!    Npc = Nfp
! end select
! discretefactor = discretefactor/Nfp

! allocate( surf(1:1) ) ! can allow for myltiple plasma boundaries 
                        ! if multiple currents are allowed; 14 Apr 16;
  
! surf(1)%Nteta = Nteta ! not used yet; used for multiple surfaces; 20170307;
! surf(1)%Nzeta = Nzeta ! not used yet; used for multiple surfaces; 20170307;
  
! SALLOCATE( surf(1)%xx, (0:Nteta-1,0:Nzeta-1), zero ) !x coordinates;
! SALLOCATE( surf(1)%yy, (0:Nteta-1,0:Nzeta-1), zero ) !y coordinates
! SALLOCATE( surf(1)%zz, (0:Nteta-1,0:Nzeta-1), zero ) !z coordinates
! SALLOCATE( surf(1)%nx, (0:Nteta-1,0:Nzeta-1), zero ) !unit nx;
! SALLOCATE( surf(1)%ny, (0:Nteta-1,0:Nzeta-1), zero ) !unit ny;
! SALLOCATE( surf(1)%nz, (0:Nteta-1,0:Nzeta-1), zero ) !unit nz;
! SALLOCATE( surf(1)%ds, (0:Nteta-1,0:Nzeta-1), zero ) !jacobian;
! SALLOCATE( surf(1)%xt, (0:Nteta-1,0:Nzeta-1), zero ) !dx/dtheta;
! SALLOCATE( surf(1)%yt, (0:Nteta-1,0:Nzeta-1), zero ) !dy/dtheta;
! SALLOCATE( surf(1)%zt, (0:Nteta-1,0:Nzeta-1), zero ) !dz/dtheta;
! SALLOCATE( surf(1)%pb, (0:Nteta-1,0:Nzeta-1), zero ) !target Bn;
 
! The center point value was used to discretize grid;
! do ii = 0, Nteta-1; teta = ( ii + half ) * pi2 / Nteta
!  do jj = 0, Nzeta-1; zeta = ( jj + half ) * pi2 / ( Nzeta*Nfp )
    
!   RR(0:2) = zero ; ZZ(0:2) = zero
    
!   do imn = 1, Nfou ; arg = bim(imn) * teta - bin(imn) * zeta
     
!    RR(0) =  RR(0) +     Rbc(imn) * cos(arg) + Rbs(imn) * sin(arg)
!    ZZ(0) =  ZZ(0) +     Zbc(imn) * cos(arg) + Zbs(imn) * sin(arg)
     
!    RR(1) =  RR(1) + ( - Rbc(imn) * sin(arg) + Rbs(imn) * cos(arg) ) * bim(imn)
!    ZZ(1) =  ZZ(1) + ( - Zbc(imn) * sin(arg) + Zbs(imn) * cos(arg) ) * bim(imn)
     
!    RR(2) =  RR(2) - ( - Rbc(imn) * sin(arg) + Rbs(imn) * cos(arg) ) * bin(imn)
!    ZZ(2) =  ZZ(2) - ( - Zbc(imn) * sin(arg) + Zbs(imn) * cos(arg) ) * bin(imn)
     
!   enddo ! end of do imn; 30 Oct 15;
    
!   szeta = sin(zeta)
!   czeta = cos(zeta)
    
!   xx(1:3) = (/   RR(0) * czeta,   RR(0) * szeta, ZZ(0) /)
!   xt(1:3) = (/   RR(1) * czeta,   RR(1) * szeta, ZZ(1) /)
!   xz(1:3) = (/   RR(2) * czeta,   RR(2) * szeta, ZZ(2) /) + (/ - RR(0) * szeta,   RR(0) * czeta, zero  /)

!   ds(1:3) = -(/ xt(2) * xz(3) - xt(3) * xz(2), & ! minus sign for theta counterclockwise direction;
!                 xt(3) * xz(1) - xt(1) * xz(3), &
!                 xt(1) * xz(2) - xt(2) * xz(1) /)

!   dd = sqrt( sum( ds(1:3)*ds(1:3) ) )

    ! x, y, z coordinates for the surface;
!   surf(1)%xx(ii,jj) = xx(1)
!   surf(1)%yy(ii,jj) = xx(2)
!   surf(1)%zz(ii,jj) = xx(3)

    ! dx/dt, dy/dt, dz/dt (dt for d theta)
!   surf(1)%xt(ii,jj) = xt(1)
!   surf(1)%yt(ii,jj) = xt(2)
!   surf(1)%zt(ii,jj) = xt(3)

    ! surface normal vectors and ds for the jacobian;
!   surf(1)%nx(ii,jj) = ds(1) / dd
!   surf(1)%ny(ii,jj) = ds(2) / dd
!   surf(1)%nz(ii,jj) = ds(3) / dd
!   surf(1)%ds(ii,jj) =         dd

!  enddo ! end of do jj; 14 Apr 16;
! enddo ! end of do ii; 14 Apr 16;

  !calculate target Bn with input harmonics; 05 Jan 17;
! if(NBnf >  0) then

!    do jj = 0, Nzeta-1 ; zeta = ( jj + half ) * pi2 / (Nzeta*Nfp)
!       do ii = 0, Nteta-1 ; teta = ( ii + half ) * pi2 / Nteta
!          do imn = 1, NBnf
!             arg = Bnim(imn) * teta - Bnin(imn) * zeta
!             surf(1)%pb(ii,jj) = surf(1)%pb(ii,jj) + Bnc(imn)*cos(arg) + Bns(imn)*sin(arg)
!          enddo
!       enddo
!    enddo

! endif
  
  
  return
  
end subroutine fousurf


subroutine surfxx( teta, zeta, xx, xt, xz, bn )
  
  use globals, only : zero, half, pi2, myid, ounit, runit, surffile, IsQuiet, IsSymmetric, &
                      Nfou, NBnf, bim, bin, Bnim, Bnin, Rbc, Rbs, Zbc, Zbs, Bnc, Bns,  &
                      Nteta, Nzeta, surf, discretefactor
  
  implicit none
  
  include "mpif.h"
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOGICAL :: exist
  INTEGER :: iosta, astat, ierr, ii, jj, imn
  REAL    :: RR(0:2), ZZ(0:2), szeta, czeta, xx(1:3), xt(1:3), xz(1:3), ds(1:3), bn(1:3), &
             teta, zeta, arg, dd
  
  RR(0:2) = zero ; ZZ(0:2) = zero
  
  do imn = 1, Nfou ; arg = bim(imn) * teta - bin(imn) * zeta
   
   RR(0) =  RR(0) +     Rbc(imn) * cos(arg) + Rbs(imn) * sin(arg)
   ZZ(0) =  ZZ(0) +     Zbc(imn) * cos(arg) + Zbs(imn) * sin(arg)
   
   RR(1) =  RR(1) + ( - Rbc(imn) * sin(arg) + Rbs(imn) * cos(arg) ) * bim(imn)
   ZZ(1) =  ZZ(1) + ( - Zbc(imn) * sin(arg) + Zbs(imn) * cos(arg) ) * bim(imn)
   
   RR(2) =  RR(2) - ( - Rbc(imn) * sin(arg) + Rbs(imn) * cos(arg) ) * bin(imn)
   ZZ(2) =  ZZ(2) - ( - Zbc(imn) * sin(arg) + Zbs(imn) * cos(arg) ) * bin(imn)
   
  enddo ! end of do imn; 30 Oct 15;
    
  szeta = sin(zeta)
  czeta = cos(zeta)
  
  xx(1:3) = (/   RR(0) * czeta,   RR(0) * szeta, ZZ(0) /)
  xt(1:3) = (/   RR(1) * czeta,   RR(1) * szeta, ZZ(1) /)
  xz(1:3) = (/   RR(2) * czeta,   RR(2) * szeta, ZZ(2) /) + (/ - RR(0) * szeta,   RR(0) * czeta, zero  /)
  
  ds(1:3) = -(/ xt(2) * xz(3) - xt(3) * xz(2), & ! minus sign for theta counterclockwise direction;
                xt(3) * xz(1) - xt(1) * xz(3), &
                xt(1) * xz(2) - xt(2) * xz(1) /)

  dd = sqrt( sum( ds(1:3)*ds(1:3) ) )


   bn = zero
  
! if(NBnf >  0) then


   
!    do jj = 0, Nzeta-1 ; zeta = ( jj + half ) * pi2 / (Nzeta*Nfp)
!       do ii = 0, Nteta-1 ; teta = ( ii + half ) * pi2 / Nteta
!          do imn = 1, NBnf
!             arg = Bnim(imn) * teta - Bnin(imn) * zeta
!             bn = bn + Bnc(imn)*cos(arg) + Bns(imn)*sin(arg)
!          enddo
!       enddo
!    enddo

! endif
 
  return
  
end subroutine surfxx

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!



!title (axis) ! Reads axis parameters from file.

!latex \briefly{The parameters describing a closed curve are read from file.}

!latex \calledby{\link{focus}}
!latex \calls{\link{knotxx}}

!latex \tableofcontents

!latex \subsection{overview}
!latex \bi
!latex \item[1.] A user-prescribed closed curve in three-dimensional space will serve as the proxy magnetic axis.
!latex \item[2.] The parameters describing the knot, namely \internal{axisxc}, \internal{axisyss} and \internal{axiszs}, are read from \verb+ext.axis+.
!latex \item[3.] A circular or elliptical cross section surface is constructed.
!latex \ei

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!!

subroutine rdknot
  
  use globals, only : zero, one, half, ten, pi2, myid, ncpu, ounit, runit, &
                      ext, &
                      axisNF, axisxc, axisxs, axisyc, axisys, axiszc, axiszs, Nfp, &
                      IsSymmetric, discretefactor, &
                      knotsurf, Nteta, Nzeta, surf
  
  implicit none
  
  include "mpif.h"
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!!
  
  LOGICAL              :: exist
  INTEGER              :: iostat, astat, ierr, ii, jj
  REAL                 :: teta, zeta, ax(1:3), at(1:3), az(1:3), xx(1:3), xt(1:3), xz(1:3), ds(1:3), dd
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!!
  
  inquire( file=trim(ext)//".axis", exist=exist )
  if( exist ) then
   if( myid.eq.0 ) then
    open(runit, file=trim(ext)//".axis", status="old", action='read', iostat=iostat)  
    write(ounit,'("readsrf : " 10x " : reading ext.axis ;")') 
   endif
  else
   inquire( file=".axis", exist=exist )
   if( exist ) then
    if( myid.eq.0 ) then
     open(runit, file=".axis", status="old", action='read', iostat=iostat)  
     write(ounit,'("rdcoil  : " 10x " : reading .axis ;")') 
    endif
   endif
  endif
  
  FATAL( rdknot, .not.exist, neither ext.axis nor .axis found )
  
  if( myid.eq.0 ) then
   read( runit, *, iostat=iostat )
   read( runit, *, iostat=iostat ) axisNF
  endif
  
  IlBCAST( axisNF, 1, 0 )
  
  SALLOCATE( axisxc, (0:axisNF), zero )
  SALLOCATE( axisxs, (0:axisNF), zero )
  SALLOCATE( axisyc, (0:axisNF), zero )
  SALLOCATE( axisys, (0:axisNF), zero )
  SALLOCATE( axiszc, (0:axisNF), zero )
  SALLOCATE( axiszs, (0:axisNF), zero )
  
  if( myid.eq.0 ) then
   
   read( runit, *, iostat=iostat ) 
   read( runit, *, iostat=iostat ) axisxc(0:axisNF)
   read( runit, *, iostat=iostat ) 
   read( runit, *, iostat=iostat ) axisxs(0:axisNF)
   read( runit, *, iostat=iostat ) 
   read( runit, *, iostat=iostat ) axisyc(0:axisNF)
   read( runit, *, iostat=iostat ) 
   read( runit, *, iostat=iostat ) axisys(0:axisNF)
   read( runit, *, iostat=iostat ) 
   read( runit, *, iostat=iostat ) axiszc(0:axisNF)
   read( runit, *, iostat=iostat ) 
   read( runit, *, iostat=iostat ) axiszs(0:axisNF)
   
   close(runit,iostat=iostat)
   
  endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!!
  
  RlBCAST( axisxc(0:axisNF), axisNF+1, 0 )
  RlBCAST( axisxs(0:axisNF), axisNF+1, 0 )
  RlBCAST( axisyc(0:axisNF), axisNF+1, 0 )
  RlBCAST( axisys(0:axisNF), axisNF+1, 0 )
  RlBCAST( axiszc(0:axisNF), axisNF+1, 0 )
  RlBCAST( axiszs(0:axisNF), axisNF+1, 0 )

  Nfp = 1
  
!  if( myid.eq.0 ) then
!   write(ounit,'("readsrf : " 10x " : axisxc =",    999es11.03)') axisxc(0:axisNF)
!   write(ounit,'("readsrf : " 10x " : axisxs =",11x,998es11.03)') axisxs(1:axisNF)
!   write(ounit,'("readsrf : " 10x " : axisyc =",    999es11.03)') axisyc(0:axisNF)
!   write(ounit,'("readsrf : " 10x " : axisys =",11x,999es11.03)') axisys(1:axisNF)
!   write(ounit,'("readsrf : " 10x " : axiszc =",    999es11.03)') axiszc(0:axisNF)
!   write(ounit,'("readsrf : " 10x " : axiszs =",11x,999es11.03)') axiszs(1:axisNF)
!  endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!!
  
! discretize the surface data; 2017/03/28; czhu;
  
! select case (IsSymmetric)
! case ( 0 )
!    Nfp = 1                          !reset Nfp to 1;
!    Npc = 1                          !number of coils periodicity 
! case ( 1 )                          !plasma periodicity enabled;
!    Npc = 1
! case ( 2 )                          !plasma and coil periodicity enabled;
!    Npc = Nfp
! end select
! discretefactor = discretefactor/Nfp

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!!
  return
  
end subroutine rdknot

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!!




subroutine knotxx( aa, teta, zeta, ax, at, az, xx, xt, xz )
  
  use globals, only : zero, one, pi2, small, myid, ounit, &
                      axisNF, axisxc, axisxs, axisyc, axisys, axiszc, axiszs
  
  implicit none
  
  include "mpif.h"
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!!
  
  REAL                 :: aa, teta, zeta, ax(1:3), at(1:3), az(1:3), xx(1:3), xt(1:3), xz(1:3)
  
  INTEGER              :: ierr, p(1:3), q(1:3), mm
  REAL                 :: cqz, sqz, cpz, spz, RR(0:3), ZZ(0:3), x0(1:3), x1(1:3), x2(1:3), x3(1:3)
  REAL                 :: a0, a1, a2, b0, b1, carg, sarg
  REAL                 :: tt(1:3), td(1:3), dd(1:3), xa, ya, za, ff, nn(1:3), nd(1:3), bb(1:3), bd(1:3)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!!
  
  x0(1:3) = (/ axisxc(0), axisyc(0), axiszc(0) /)
  x1(1:3) = (/ zero  , zero  , zero   /)
  x2(1:3) = (/ zero  , zero  , zero   /)
  x3(1:3) = (/ zero  , zero  , zero   /)
  
  if( aa.gt.zero ) then ! will need additional derivatives to construct normal and binormal; 14 Apr 16;
   
   do mm = 1, axisNF ; carg = cos(mm*zeta) ; sarg = sin(mm*zeta)
    x0(1:3) = x0(1:3) + ( + (/ axisxc(mm), axisyc(mm), axiszc(mm) /) * carg + &
                            (/ axisxs(mm), axisys(mm), axiszs(mm) /) * sarg ) 
    x1(1:3) = x1(1:3) + ( - (/ axisxc(mm), axisyc(mm), axiszc(mm) /) * sarg + &
                            (/ axisxs(mm), axisys(mm), axiszs(mm) /) * carg ) * mm   
    x2(1:3) = x2(1:3) + ( - (/ axisxc(mm), axisyc(mm), axiszc(mm) /) * carg - &
                            (/ axisxs(mm), axisys(mm), axiszs(mm) /) * sarg ) * mm**2
    x3(1:3) = x3(1:3) + ( + (/ axisxc(mm), axisyc(mm), axiszc(mm) /) * sarg - &
                            (/ axisxs(mm), axisys(mm), axiszs(mm) /) * carg ) * mm**3
   enddo
   
  else
   
   do mm = 1, axisNF ; carg = cos(mm*zeta) ; sarg = sin(mm*zeta)
    x0(1:3) = x0(1:3) + ( + (/ axisxc(mm), axisyc(mm), axiszc(mm) /) * carg + &
                            (/ axisxs(mm), axisys(mm), axiszs(mm) /) * sarg )       
    x1(1:3) = x1(1:3) + ( - (/ axisxc(mm), axisyc(mm), axiszc(mm) /) * sarg + &
                            (/ axisxs(mm), axisys(mm), axiszs(mm) /) * carg ) * mm   
   enddo
    
  endif ! end of if( aa.gt.zero ) ; 14 Apr 16;
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!!
  
  ax(1:3) = x0(1:3)                                                        ! Nov 12 15;
  at(1:3) = zero                                                           ! Nov 12 15;
  az(1:3) = x1(1:3)                                                        ! Nov 12 15;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!!
  
  if( aa.gt.zero ) then
   
   a0      = sqrt( x1(1)*x1(1) + x1(2)*x1(2) + x1(3)*x1(3) )                ! Nov 12 15;
   a1      =     ( x1(1)*x2(1) + x1(2)*x2(2) + x1(3)*x2(3) ) / a0           ! Nov 12 15;
   a2      =     ( x2(1)*x2(1) + x2(2)*x2(2) + x2(3)*x2(3)   &
               +   x1(1)*x3(1) + x1(2)*x3(2) + x1(3)*x3(3) - a1 * a1 ) / a0 ! Nov 12 15;

   tt(1:3) =   x1(1:3)                                                / a0  ! Nov 12 15;
   td(1:3) = ( x2(1:3) - tt(1:3) * a1                               ) / a0  ! Nov 12 15;
   dd(1:3) = ( x3(1:3) - td(1:3) * a1 - tt(1:3) * a2 - td(1:3) * a1 ) / a0  ! Nov 12 15;
   
   xa = ( x2(1) - tt(1) * a1 ) ! Nov 12 15;
   ya = ( x2(2) - tt(2) * a1 ) ! Nov 12 15;
   za = ( x2(3) - tt(3) * a1 ) ! Nov 12 15;
   
   ff = sqrt( xa**2 + ya**2 + za**2 ) ! Nov 12 15;
   
   b0 = ff / a0 ! Nov 12 15;
   
   b1 = ( ( xa * ( x3(1) - td(1) * a1 - tt(1) * a2 ) &
          + ya * ( x3(2) - td(2) * a1 - tt(2) * a2 ) &
          + za * ( x3(3) - td(3) * a1 - tt(3) * a2 ) ) / ff - b0 * a1 ) / a0 ! Nov 12 15;

   nn(1:3) =   td(1:3)                  / b0                                 ! Nov 12 15;
   nd(1:3) = ( dd(1:3) - nn(1:3) * b1 ) / b0                                 ! Nov 12 15;
   
   bb(1:3) = (/ tt(2)*nn(3)-tt(3)*nn(2), tt(3)*nn(1)-tt(1)*nn(3), tt(1)*nn(2)-tt(2)*nn(1) /)
   bd(1:3) = (/ td(2)*nn(3)-td(3)*nn(2), td(3)*nn(1)-td(1)*nn(3), td(1)*nn(2)-td(2)*nn(1) /) &
           + (/ tt(2)*nd(3)-tt(3)*nd(2), tt(3)*nd(1)-tt(1)*nd(3), tt(1)*nd(2)-tt(2)*nd(1) /)
   
   xx(1:3) = x0(1:3) + aa * (   cos(teta) * nn(1:3) - sin(teta) * bb(1:3) ) ! aa is minor radius;
   xt(1:3) = zero    + aa * ( - sin(teta) * nn(1:3) - cos(teta) * bb(1:3) )
   xz(1:3) = x1(1:3) + aa * (   cos(teta) * nd(1:3) - sin(teta) * bd(1:3) )
   
  else
   
   xx(1:3) = zero
   xt(1:3) = zero
   xz(1:3) = zero
   
  endif ! end of if( aa.gt.zero ) ; SRH; 28 Aug 17;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!!

  return
  
end subroutine knotxx

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!!
