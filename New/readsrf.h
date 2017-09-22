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
  
  use globals, only : zero, half, three, pi2, sqrtmachprec, myid, ounit, tstart, &
                      axisfile, &
                      Nt, Nz, surf, Isurface, minorrad
  
  implicit none
  
  include "mpif.h"
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  INTEGER :: astat, ierr, ii, jj
  REAL    :: teta, zeta, xx(1:3), xt(1:3), xz(1:3), bn, ax(1:3), at(1:3), az(1:3), ds(1:3), dd, tnow
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  tnow = MPI_WTIME()
  
  select case( Isurface )
   
  case( 0 )
   
   if( myid.eq.0 ) write(ounit,'("readsrf : ",f10.1," : Isurface ="i2" ; Nt ="i4", Nz ="i4" ;")') tnow-tstart, Isurface, Nt, Nz
   
   call fousurf  ! VMEC  -style plasma boundary;
   
  case( 1 )
   
   if( myid.eq.0 ) write(ounit,1000) tnow-tstart, Isurface, Nt, Nz, minorrad, trim(axisfile)
   
   call rdaxis   ! axis  -style plasma boundary;
   
  case default
   
   FATAL( focus , .true., selected value of Isurface is not supported )
   
  end select
  
1000 format("readsrf : ",f10.1," : Isurface ="i2" ; Nt ="i4", Nz ="i4" ; minorrad ="f6.3" ; reading ",a," ;")
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  surf%Nteta = Nt
  surf%Nzeta = Nz
  
  CHECK( readsrf, Nt.lt.0, error )
  CHECK( readsrf, Nz.lt.0, error )
  
  SALLOCATE( surf%xx, (0:Nt,0:Nz ), zero )
  SALLOCATE( surf%yy, (0:Nt,0:Nz ), zero )
  SALLOCATE( surf%zz, (0:Nt,0:Nz ), zero )
  
  SALLOCATE( surf%nx, (0:Nt,0:Nz ), zero )
  SALLOCATE( surf%ny, (0:Nt,0:Nz ), zero )
  SALLOCATE( surf%nz, (0:Nt,0:Nz ), zero )
  
  SALLOCATE( surf%ds, (0:Nt,0:Nz ), zero )
  
  SALLOCATE( surf%xt, (0:Nt,0:Nz ), zero )
  SALLOCATE( surf%yt, (0:Nt,0:Nz ), zero )
  SALLOCATE( surf%zt, (0:Nt,0:Nz ), zero )
  
! SALLOCATE( surf%Bn, (0:Nt,0:Nz ), zero ) ! normal field;
! SALLOCATE( surf%Tf, (     0:Nz ), zero ) ! toroidal flux;
  SALLOCATE( surf%Bp, (0:Nt,0:Nz ), zero ) ! normal field;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  CHECK( readsrf, Nt.le.0, error )
  CHECK( readsrf, Nz.le.0, error )
  
  do ii = 0, Nt ; teta = ( ii + half ) * pi2 / Nt ! grid center;
   do jj = 0, Nz ; zeta = ( jj + half ) * pi2 / Nz
    
    select case( Isurface )
    case( 0 ) ; call surfxx(           teta, zeta,                            xx(1:3), xt(1:3), xz(1:3), bn )
    case( 1 ) ; call knotxx( minorrad, teta, zeta, ax(1:3), at(1:3), az(1:3), xx(1:3), xt(1:3), xz(1:3)     ) ; bn = zero
    case default
     FATAL( focus , .true., selected Isurface is not supported )
    end select
    
    ds(1:3) = - (/ xt(2) * xz(3) - xt(3) * xz(2), &
                   xt(3) * xz(1) - xt(1) * xz(3), &
                   xt(1) * xz(2) - xt(2) * xz(1) /) ! negative sign = counterclockwise;
    
    dd = sqrt( sum( ds(1:3)*ds(1:3) ) )
    
    surf%xx(ii,jj) = xx(1)
    surf%yy(ii,jj) = xx(2)
    surf%zz(ii,jj) = xx(3)
    
    surf%xt(ii,jj) = xt(1)
    surf%yt(ii,jj) = xt(2)
    surf%zt(ii,jj) = xt(3)
    
    CHECK( readsrf, dd.lt.sqrtmachprec, divide by zero )
    
    surf%nx(ii,jj) = ds(1) / dd
    surf%ny(ii,jj) = ds(2) / dd
    surf%nz(ii,jj) = ds(3) / dd
    
    surf%ds(ii,jj) =         dd
    
    surf%Bp(ii,jj) = bn ! target/plasma normal field;
    
    surf%area = surf%area + sum(xx(1:3)*ds(1:3))

   enddo ! end of do jj;
  enddo ! end of do ii;

  surf%area = surf%area * (pi2/Nt) * (pi2/Nz) / three

 !if( myid.eq.0 ) write(ounit,'("readsrf : ", 10x ," : surface area =",e13.5," ;")') surf%area
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
! SALLOCATE( surf%area, (0:Nt-1), zero )
  
! do jj = 0, Nz-1
   
!  surf%area(jj) = sum( surf%xx(0:Nt-1,jj)*surf%yt(0:Nt-1,jj)-surf%yy(0:Nt-1,jj)*surf%xt(0:Nt-1,jj) ) * (pi2/Nt)
   
!  if( myid.eq.0 .and. jj.eq.0 ) write(ounit,'("readsrf : ", 10x ," : ",i4," cross-section area =",e13.5" ;")') jj, surf%area(jj)
   
! enddo
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  return
  
end subroutine readsrf

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

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

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

subroutine fousurf
  
  use globals, only : zero, half, pi2, myid, ounit, runit, surffile, &
                      Nfp, &
                      Nfou, NBnf, bim, bin, Bnim, Bnin, Rbc, Rbs, Zbc, Zbs, Bnc, Bns,  &
                      Nt, Nz, surf
  
  implicit none
  
  include "mpif.h"
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  LOGICAL :: exist
  INTEGER :: iosta, astat, ierr, imn
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
! if( myid == 0 .and. IsQuiet <= 0) then
!    write(ounit, *) "-----------Reading surface-----------------------------------"
!    write(ounit, '("surface : Nfou = " I06 " ; Nfp = " I06 " ; NBnf = " I06 " ;" )') Nfou, Nfp, NBnf
! endif

! if( myid == 0 .and. IsQuiet <= -2) then !very detailed output;
!    write(ounit,'("        : " 10x " : bim ="10i13   )') bim(1:Nfou)
!    write(ounit,'("        : " 10x " : bin ="10i13   )') bin(1:Nfou)
!    write(ounit,'("        : " 10x " : Rbc ="10es13.5)') Rbc(1:Nfou)
!    write(ounit,'("        : " 10x " : Rbs ="10es13.5)') Rbs(1:Nfou)
!    write(ounit,'("        : " 10x " : Zbc ="10es13.5)') Zbc(1:Nfou)
!    write(ounit,'("        : " 10x " : Zbs ="10es13.5)') Zbs(1:Nfou)
!    if(Nbnf > 0) then
!       write(ounit,'("        : " 10x " : Bnim ="10i13  )') Bnim(1:NBnf)
!       write(ounit,'("        : " 10x " : Bnin ="10i13  )') Bnin(1:NBnf)
!       write(ounit,'("        : " 10x " : Bnc ="10es13.5)') Bnc (1:NBnf)
!       write(ounit,'("        : " 10x " : Bns ="10es13.5)') Bns (1:NBnf)
!    endif
! endif

! allocate( surf(1:1) ) ! can allow for myltiple plasma boundaries 
                        ! if multiple currents are allowed; 14 Apr 16;
  
! surf%Nteta = Nt ! not used yet; used for multiple surfaces; 20170307;
! surf%Nzeta = Nz ! not used yet; used for multiple surfaces; 20170307;
  
! SALLOCATE( surf%xx, (0:Nt-1,0:Nz-1), zero ) !x coordinates;
! SALLOCATE( surf%yy, (0:Nt-1,0:Nz-1), zero ) !y coordinates
! SALLOCATE( surf%zz, (0:Nt-1,0:Nz-1), zero ) !z coordinates
! SALLOCATE( surf%nx, (0:Nt-1,0:Nz-1), zero ) !unit nx;
! SALLOCATE( surf%ny, (0:Nt-1,0:Nz-1), zero ) !unit ny;
! SALLOCATE( surf%nz, (0:Nt-1,0:Nz-1), zero ) !unit nz;
! SALLOCATE( surf%ds, (0:Nt-1,0:Nz-1), zero ) !jacobian;
! SALLOCATE( surf%xt, (0:Nt-1,0:Nz-1), zero ) !dx/dtheta;
! SALLOCATE( surf%yt, (0:Nt-1,0:Nz-1), zero ) !dy/dtheta;
! SALLOCATE( surf%zt, (0:Nt-1,0:Nz-1), zero ) !dz/dtheta;
! SALLOCATE( surf%pb, (0:Nt-1,0:Nz-1), zero ) !target Bn;
 
! The center point value was used to discretize grid;
! do ii = 0, Nt-1; teta = ( ii + half ) * pi2 / Nt
!  do jj = 0, Nz-1; zeta = ( jj + half ) * pi2 / ( Nz*Nfp )
    
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
!   surf%xx(ii,jj) = xx(1)
!   surf%yy(ii,jj) = xx(2)
!   surf%zz(ii,jj) = xx(3)

    ! dx/dt, dy/dt, dz/dt (dt for d theta)
!   surf%xt(ii,jj) = xt(1)
!   surf%yt(ii,jj) = xt(2)
!   surf%zt(ii,jj) = xt(3)

    ! surface normal vectors and ds for the jacobian;
!   surf%nx(ii,jj) = ds(1) / dd
!   surf%ny(ii,jj) = ds(2) / dd
!   surf%nz(ii,jj) = ds(3) / dd
!   surf%ds(ii,jj) =         dd

!  enddo ! end of do jj; 14 Apr 16;
! enddo ! end of do ii; 14 Apr 16;

  !calculate target Bn with input harmonics; 05 Jan 17;
! if(NBnf >  0) then

!    do jj = 0, Nz-1 ; zeta = ( jj + half ) * pi2 / (Nz*Nfp)
!       do ii = 0, Nt-1 ; teta = ( ii + half ) * pi2 / Nt
!          do imn = 1, NBnf
!             arg = Bnim(imn) * teta - Bnin(imn) * zeta
!             surf%pb(ii,jj) = surf%pb(ii,jj) + Bnc(imn)*cos(arg) + Bns(imn)*sin(arg)
!          enddo
!       enddo
!    enddo

! endif
  
  
  return
  
end subroutine fousurf

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

subroutine surfxx( teta, zeta, xx, xt, xz, bn )
  
  use globals, only : zero, half, pi2, myid, ounit, runit, surffile, &
                      Nfou, NBnf, bim, bin, Bnim, Bnin, Rbc, Rbs, Zbc, Zbs, Bnc, Bns,  &
                      Nt, Nz, surf
  
  implicit none
  
  include "mpif.h"
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

  LOGICAL :: exist
  INTEGER :: iosta, astat, ierr, imn
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


   
!    do jj = 0, Nz-1 ; zeta = ( jj + half ) * pi2 / (Nz*Nfp)
!       do ii = 0, Nt-1 ; teta = ( ii + half ) * pi2 / Nt
!          do imn = 1, NBnf
!             arg = Bnim(imn) * teta - Bnin(imn) * zeta
!             bn = bn + Bnc(imn)*cos(arg) + Bns(imn)*sin(arg)
!          enddo
!       enddo
!    enddo

! endif
 
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

  return
  
end subroutine surfxx

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

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

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

subroutine rdaxis
  
  use globals, only : zero, one, half, ten, pi2, myid, ncpu, ounit, runit, &
                      ext, axisfile, &
                      axisNF, axisxc, axisxs, axisyc, axisys, axiszc, axiszs, Nfp, &
                      minorrad, Nt, Nz, surf
  
  implicit none
  
  include "mpif.h"
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

  LOGICAL              :: exist
  INTEGER              :: iostat, astat, ierr
  REAL                 :: teta, zeta, ax(1:3), at(1:3), az(1:3), xx(1:3), xt(1:3), xz(1:3), ds(1:3), dd
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  inquire( file=trim(axisfile), exist=exist )
  
  FATAL( rdaxis , .not.exist, axisfile does not exist )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  if( myid.eq.0 ) then
   
   open( runit, file=trim(axisfile), status="old", action='read' )
   read( runit, * )
   read( runit, * ) axisNF ! write(ounit,'("readsrf : " 10x " : axisNF = ",i3," ;")') axisNF
   
  endif ! end of if( myid.eq.0 ) x;
  
  IlBCAST( axisNF, 1, 0 )
  
  SALLOCATE( axisxc, (0:axisNF), zero )
  SALLOCATE( axisxs, (0:axisNF), zero )
  SALLOCATE( axisyc, (0:axisNF), zero )
  SALLOCATE( axisys, (0:axisNF), zero )
  SALLOCATE( axiszc, (0:axisNF), zero )
  SALLOCATE( axiszs, (0:axisNF), zero )
  
  if( myid.eq.0 ) then
   
   read( runit, * ) 
   read( runit, * ) axisxc(0:axisNF) ! write(ounit,'("readsrf : " 10x " : axisxc = ",99es13.05," ;")') axisxc
   read( runit, * ) 
   read( runit, * ) axisxs(0:axisNF)
   read( runit, * ) 
   read( runit, * ) axisyc(0:axisNF)
   read( runit, * ) 
   read( runit, * ) axisys(0:axisNF)
   read( runit, * ) 
   read( runit, * ) axiszc(0:axisNF)
   read( runit, * ) 
   read( runit, * ) axiszs(0:axisNF)
   
   close(runit)
   
  endif ! end of if( myid.eq.0 ) ;
  
  
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
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

  return
  
end subroutine rdaxis

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

subroutine knotxx( aa, teta, zeta, ax, at, az, xx, xt, xz )
  
  use globals, only : zero, one, pi2, small, myid, ounit, &
                      axisNF, axisxc, axisxs, axisyc, axisys, axiszc, axiszs
  
  implicit none
  
  include "mpif.h"
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  REAL                 :: aa, teta, zeta, ax(1:3), at(1:3), az(1:3), xx(1:3), xt(1:3), xz(1:3)
  
  INTEGER              :: ierr, p(1:3), q(1:3), mm
  REAL                 :: cqz, sqz, cpz, spz, RR(0:3), ZZ(0:3), x0(1:3), x1(1:3), x2(1:3), x3(1:3)
  REAL                 :: a0, a1, a2, b0, b1, carg, sarg
  REAL                 :: tt(1:3), td(1:3), dd(1:3), xa, ya, za, ff, nn(1:3), nd(1:3), bb(1:3), bd(1:3)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  CHECK( readsrf, .not.allocated(axisxc), illegal )
  CHECK( readsrf, .not.allocated(axisxs), illegal )
  CHECK( readsrf, .not.allocated(axisyc), illegal )
  CHECK( readsrf, .not.allocated(axisys), illegal )
  CHECK( readsrf, .not.allocated(axiszc), illegal )
  CHECK( readsrf, .not.allocated(axiszs), illegal )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

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
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

  ax(1:3) = x0(1:3)                                                        ! Nov 12 15;
  at(1:3) = zero                                                           ! Nov 12 15;
  az(1:3) = x1(1:3)                                                        ! Nov 12 15;
    
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

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
   
  endif ! end of if( aa.gt.zero );
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

  return
  
end subroutine knotxx

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
