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
  
  use globals, only : zero, one, two, half, three, pi, pi2, sqrtmachprec, myid, ounit, tstart, &
                      axisfile, &
                      Nt, Nz, deltateta, deltazeta, surf, Isurface, minorrad, ellipticity, nrotate, zetaoff
  
  implicit none
  
  include "mpif.h"
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  INTEGER :: astat, ierr, ii, jj, isum, isurf
  REAL    :: teta, zeta, xx(1:3), xs(1:3), xu(1:3), xv(1:3), bn, ax(1:3), at(1:3), az(1:3), ds(1:3), dd, tnow, srho, v1(1:3), v2(1:3), w1(1:3), w2(1:3)
  REAL    :: xtt(1:3), xtz(1:3), xzz(1:3)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  tnow = MPI_WTIME()
  
  if( myid.eq.0 ) write(ounit,'("readsrf : ",f10.1," : Isurface ="i2" ; Nt ="i5", Nz ="i5" ;")') tnow-tstart, Isurface, Nt, Nz
  
  select case( Isurface )
   
  case( 0 )
   
   call fousurf  ! VMEC  -style plasma boundary;
   
  case( 1 )
   
   if( myid.eq.0 ) write(ounit,1000) minorrad, ellipticity, nrotate, zetaoff, trim(axisfile(1))//" & "//trim(axisfile(2))
   
1000 format("readsrf : ", 10x ," : minorrad ="f12.6" ; ellipticity ="f12.6" ; nrotate ="i3" ; zetaoff ="f12.6" ; reading ",a," ;")
   
   call rdaxis   ! axis  -style plasma boundary;
   
  case default
   
   FATAL( focus , .true., selected value of Isurface is not supported )
   
  end select
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  CHECK( readsrf, Nt.le.0, error )
  CHECK( readsrf, Nz.le.0, error )
  
  allocate( surf(1:2) )
  
  do isurf = 1, 2
   
   surf(isurf)%Nteta = Nt ; surf(isurf)%Nzeta = Nz
   
   SALLOCATE( surf(isurf)%xx , (1:3,0:Nt-1,0:Nz-1 ), zero )
   SALLOCATE( surf(isurf)%xs , (1:3,0:Nt-1,0:Nz-1 ), zero )
   SALLOCATE( surf(isurf)%xu , (1:3,0:Nt-1,0:Nz-1 ), zero )
   SALLOCATE( surf(isurf)%xv , (1:3,0:Nt-1,0:Nz-1 ), zero )
   SALLOCATE( surf(isurf)%sg , (    0:Nt-1,0:Nz-1 ), zero )
   SALLOCATE( surf(isurf)%gs , (1:3,0:Nt-1,0:Nz-1 ), zero )
   SALLOCATE( surf(isurf)%gu , (1:3,0:Nt-1,0:Nz-1 ), zero )
   SALLOCATE( surf(isurf)%gv , (1:3,0:Nt-1,0:Nz-1 ), zero )
   SALLOCATE( surf(isurf)%nn , (1:3,0:Nt-1,0:Nz-1 ), zero )
   SALLOCATE( surf(isurf)%ds , (    0:Nt-1,0:Nz-1 ), zero )
   SALLOCATE( surf(isurf)%guv, (1:6,0:Nt-1,0:Nz-1 ), zero )
   
! SALLOCATE( surf(isurf)%Bn , (    0:Nt-1,0:Nz-1 ), zero )
! SALLOCATE( surf(isurf)%Tf , (           0:Nz   ), zero )
   
   SALLOCATE( surf(isurf)%BB , (1:3,0:Nt-1,0:Nz-1 ), zero )
   SALLOCATE( surf(isurf)%Bn , (    0:Nt-1,0:Nz-1 ), zero )
   
   SALLOCATE( surf(isurf)%Bp , (    0:Nt-1,0:Nz-1 ), zero ) ! plasma normal field;
   
   SALLOCATE( surf(isurf)%EE , (    0:Nt-1,0:Nz-1 ), zero )
   SALLOCATE( surf(isurf)%FF , (    0:Nt-1,0:Nz-1 ), zero )
   SALLOCATE( surf(isurf)%GG , (    0:Nt-1,0:Nz-1 ), zero )
   SALLOCATE( surf(isurf)%LL , (    0:Nt-1,0:Nz-1 ), zero )
   SALLOCATE( surf(isurf)%MM , (    0:Nt-1,0:Nz-1 ), zero )
   SALLOCATE( surf(isurf)%PP , (    0:Nt-1,0:Nz-1 ), zero )
   SALLOCATE( surf(isurf)%HH , (    0:Nt-1,0:Nz-1 ), zero )
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
   
   srho = one
   
   surf(isurf)%area = zero ; surf(isurf)%vol = zero
   
   do ii = 0, Nt-1 ; teta = ( ii + half ) * deltateta ! please remind me why half-grid is required;
    do jj = 0, Nz-1 ; zeta = ( jj + half ) * deltazeta ! please remind me why half-grid is required;
     
     select case( Isurface )
     case( 0 ) ; call surfxx(              teta, zeta,             xx,     xu, xv, xtt, xtz, xzz,                    bn ) ; xs = zero
     case( 1 ) ; call knotxx( isurf, srho, teta, zeta, ax, at, az, xx, xs, xu, xv, xtt, xtz, xzz, v1, v2, w1, w2        ) ; bn = zero
     case default
      FATAL( focus , .true., selected Isurface is not supported )
     end select
     
     surf(isurf)%xx(1:3,ii,jj) = xx(1:3)
     
     surf(isurf)%xs(1:3,ii,jj) = xs(1:3)
     surf(isurf)%xu(1:3,ii,jj) = xu(1:3)
     surf(isurf)%xv(1:3,ii,jj) = xv(1:3)
     
     surf(isurf)%guv(1,ii,jj) = sum( xs(1:3)*xs(1:3) )
     surf(isurf)%guv(2,ii,jj) = sum( xs(1:3)*xu(1:3) )
     surf(isurf)%guv(3,ii,jj) = sum( xs(1:3)*xv(1:3) )
     surf(isurf)%guv(4,ii,jj) = sum( xu(1:3)*xu(1:3) )
     surf(isurf)%guv(5,ii,jj) = sum( xu(1:3)*xv(1:3) )
     surf(isurf)%guv(6,ii,jj) = sum( xv(1:3)*xv(1:3) )
     
     isum = 0 ; call cross( isum, xu(1:3), xv(1:3), ds(1:3) )
     
     surf(isurf)%sg(ii,jj) = sum( xs(1:3)*ds(1:3) )
     
     isum = 0 ; call cross( isum, xu(1:3), xv(1:3), surf(isurf)%gs(1:3,ii,jj) ) ; surf(isurf)%gs(1:3,ii,jj) = surf(isurf)%gs(1:3,ii,jj) / surf(isurf)%sg(ii,jj)
     isum = 0 ; call cross( isum, xv(1:3), xs(1:3), surf(isurf)%gu(1:3,ii,jj) ) ; surf(isurf)%gu(1:3,ii,jj) = surf(isurf)%gu(1:3,ii,jj) / surf(isurf)%sg(ii,jj)
     isum = 0 ; call cross( isum, xs(1:3), xu(1:3), surf(isurf)%gv(1:3,ii,jj) ) ; surf(isurf)%gv(1:3,ii,jj) = surf(isurf)%gv(1:3,ii,jj) / surf(isurf)%sg(ii,jj)
     
     dd = sqrt( sum( ds(1:3)*ds(1:3) ) )
     
     CHECK( readsrf, dd.lt.sqrtmachprec, divide by zero )
     
     surf(isurf)%nn(1:3,ii,jj) = ds(1:3) / dd
     surf(isurf)%ds(    ii,jj) =           dd
     surf(isurf)%Bp(    ii,jj) = bn           ! target/plasma normal field;
     
     surf(isurf)%area = surf(isurf)%area +              dd        ! area  ; 12 Nov 17;
     surf(isurf)%vol  = surf(isurf)%vol  + sum( xx(1:3)*ds(1:3) ) ! volume; 12 Nov 17;
     
     surf(isurf)%EE(ii,jj) = sum( xu(1:3) *  xu(1:3) )
     surf(isurf)%FF(ii,jj) = sum( xu(1:3) *  xv(1:3) )
     surf(isurf)%GG(ii,jj) = sum( xv(1:3) *  xv(1:3) )
     
     surf(isurf)%LL(ii,jj) = sum( ds(1:3) * xtt(1:3) ) / dd
     surf(isurf)%MM(ii,jj) = sum( ds(1:3) * xtz(1:3) ) / dd
     surf(isurf)%PP(ii,jj) = sum( ds(1:3) * xzz(1:3) ) / dd
     
    enddo ! end of do jj;
    
   enddo ! end of do ii;
   
   surf(isurf)%area = surf(isurf)%area * (pi2/Nt) * (pi2/Nz)
   surf(isurf)%vol  = surf(isurf)%vol  * (pi2/Nt) * (pi2/Nz) / three
   
   surf(isurf)%HH = ( surf(isurf)%LL*surf(isurf)%GG + surf(isurf)%PP*surf(isurf)%EE - two*surf(isurf)%MM*surf(isurf)%FF ) &
                  / ( surf(isurf)%EE*surf(isurf)%GG - surf(isurf)%FF*surf(isurf)%FF ) ! mean curvature;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
   
  enddo ! end of do isurf; 16 Nov 17;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

! do ii = 0, Nt-1
!  do jj = 0, Nz-1
!   surf(2)%xx(1:3,ii,jj) = surf(2)%xx(1:3,ii,jj) - surf(1)%xx(1:3,ii,jj) ! difference in position; 22 Nov 17;
!  enddo
! enddo
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  select case( Isurface )
  case( 0 ) ; write(ounit,'("readsrf : " 10x " : mean curvature not calculated for Isurface = 0 ; please revise surfxx ;")')
  case( 1 ) ; 
  case default
   FATAL( focus , .true., selected Isurface is not supported )
  end select
  
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

subroutine surfxx( teta, zeta, xx, xt, xz, xtt, xtz, xzz, bn )
  
  use globals, only : zero, half, pi2, myid, ounit, runit, surffile, &
                      Nfou, NBnf, bim, bin, Bnim, Bnin, Rbc, Rbs, Zbc, Zbs, Bnc, Bns,  &
                      Nt, Nz, surf
  
  implicit none
  
  include "mpif.h"
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

  LOGICAL :: exist
  INTEGER :: iosta, astat, ierr, imn
  REAL    :: RR(0:2), ZZ(0:2), szeta, czeta, xx(1:3), xt(1:3), xz(1:3), ds(1:3), bn(1:3), xtt(1:3), xtz(1:3), xzz(1:3), teta, zeta, arg, dd
  
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
  
  xtt(1:3) = zero ! second derivatives;
  xtz(1:3) = zero ! second derivatives;
  xzz(1:3) = zero ! second derivatives;
  
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
                      axisfile, &
                      axis, Nfp, &
                      minorrad, Nt, Nz, surf
  
  implicit none
  
  include "mpif.h"
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

  LOGICAL              :: exist
  INTEGER              :: iostat, astat, ierr, mm, NK , iaxis
  INTEGER, allocatable :: il(:), im(:), in(:)
  REAL                 :: teta, zeta, ax(1:3), at(1:3), az(1:3), xx(1:3), xt(1:3), xz(1:3), ds(1:3), dd, integratedtorsion
  REAL   , allocatable :: dx(:,:), dy(:,:), dz(:,:), a(:), b(:), c(:), d(:), tau(:)
    
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  allocate( axis(1:2) )
  
  do iaxis = 1, 2
   
   inquire( file=trim(axisfile(iaxis)), exist=exist )
   
   FATAL( rdaxis , iaxis.eq.1 .and. .not.exist, axisfile(1) does not exist )
   
   if( iaxis.eq.2 .and. .not.exist ) then
    write(ounit,'("readsrf : ", 10x ," : ",a," does not exist ; using axisfile(2) = ",a," ;")') trim(axisfile(iaxis)), trim(axisfile(1))
    axisfile(iaxis) = axisfile(1)
   endif
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
   
   if( myid.eq.0 ) then
    
    open( runit, file=trim(axisfile(iaxis)), status="old", action='read' )
    read( runit, * )
    read( runit, * ) axis(iaxis)%NF ! write(ounit,'("readsrf : " 10x " : axisNF = ",i3," ;")') axisNF
    
   endif ! end of if( myid.eq.0 ) ;
   
   IlBCAST( axis(iaxis)%NF, 1, 0 )
   
   SALLOCATE( axis(iaxis)%xc, (0:axis(iaxis)%NF), zero )
   SALLOCATE( axis(iaxis)%xs, (0:axis(iaxis)%NF), zero )
   SALLOCATE( axis(iaxis)%yc, (0:axis(iaxis)%NF), zero )
   SALLOCATE( axis(iaxis)%ys, (0:axis(iaxis)%NF), zero )
   SALLOCATE( axis(iaxis)%zc, (0:axis(iaxis)%NF), zero )
   SALLOCATE( axis(iaxis)%zs, (0:axis(iaxis)%NF), zero )
   
   SALLOCATE( axis(iaxis)%tc, (0:axis(iaxis)%NF), zero )
   SALLOCATE( axis(iaxis)%ts, (0:axis(iaxis)%NF), zero )
   
   if( myid.eq.0 ) then
    
    read( runit, * ) 
    read( runit, * ) axis(iaxis)%xc(0:axis(iaxis)%NF) ! write(ounit,'("readsrf : " 10x " : axisxc = ",99es13.05," ;")') axisxc
    read( runit, * ) 
    read( runit, * ) axis(iaxis)%xs(0:axis(iaxis)%NF)
    read( runit, * ) 
    read( runit, * ) axis(iaxis)%yc(0:axis(iaxis)%NF)
    read( runit, * ) 
    read( runit, * ) axis(iaxis)%ys(0:axis(iaxis)%NF)
    read( runit, * ) 
    read( runit, * ) axis(iaxis)%zc(0:axis(iaxis)%NF)
    read( runit, * ) 
    read( runit, * ) axis(iaxis)%zs(0:axis(iaxis)%NF)
    
    close(runit)
    
   endif ! end of if( myid.eq.0 ) ;
   
   RlBCAST( axis(iaxis)%xc(0:axis(iaxis)%NF), axis(iaxis)%NF+1, 0 )
   RlBCAST( axis(iaxis)%xs(0:axis(iaxis)%NF), axis(iaxis)%NF+1, 0 )
   RlBCAST( axis(iaxis)%yc(0:axis(iaxis)%NF), axis(iaxis)%NF+1, 0 )
   RlBCAST( axis(iaxis)%ys(0:axis(iaxis)%NF), axis(iaxis)%NF+1, 0 )
   RlBCAST( axis(iaxis)%zc(0:axis(iaxis)%NF), axis(iaxis)%NF+1, 0 )
   RlBCAST( axis(iaxis)%zs(0:axis(iaxis)%NF), axis(iaxis)%NF+1, 0 )
   
   Nfp = 1
   
!  if( myid.eq.0 ) then
!   write(ounit,'("rdaxis  : " 10x " : axisxc =",    999f11.07)') axis(iaxis)%xc(0:axis(iaxis)%NF)
!   write(ounit,'("rdaxis  : " 10x " : axisxs =",11x,998f11.07)') axis(iaxis)%xs(1:axis(iaxis)%NF)
!   write(ounit,'("rdaxis  : " 10x " : axisyc =",    999f11.07)') axis(iaxis)%yc(0:axis(iaxis)%NF)
!   write(ounit,'("rdaxis  : " 10x " : axisys =",11x,999f11.07)') axis(iaxis)%ys(1:axis(iaxis)%NF)
!   write(ounit,'("rdaxis  : " 10x " : axiszc =",    999f11.07)') axis(iaxis)%zc(0:axis(iaxis)%NF)
!   write(ounit,'("rdaxis  : " 10x " : axiszs =",11x,999f11.07)') axis(iaxis)%zs(1:axis(iaxis)%NF)
!  endif
   
  enddo ! end of do iaxis = 1, 2 ; 16 Nov 17;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  call inttorsion( axis(1)%NF, axis(1)%xc, axis(1)%xs, axis(1)%yc, axis(1)%ys, axis(1)%zc, axis(1)%zs, integratedtorsion )
  
!  SALLOCATE( il, (0:axis(1)%NF), (/ ( mm   , mm = 0, axis(1)%NF ) /) )
!  SALLOCATE( im, (0:axis(1)%NF), (/ ( mm**2, mm = 0, axis(1)%NF ) /) )
!  SALLOCATE( in, (0:axis(1)%NF), (/ ( mm**3, mm = 0, axis(1)%NF ) /) )
!  
!  NK = 3 * 4 * axis(1)%NF ! discrete resolution for resolving torsion;
!  
!  SALLOCATE( dx, (0:NK-1,0:3), zero )
!  SALLOCATE( dy, (0:NK-1,0:3), zero )
!  SALLOCATE( dz, (0:NK-1,0:3), zero )
!  
!  call fft( NK, dx(0:NK-1,0), axis(1)%NF,                      axis(1)%xc(0:axis(1)%NF),                      axis(1)%xs(0:axis(1)%NF), 1 )
!  call fft( NK, dy(0:NK-1,0), axis(1)%NF,                      axis(1)%yc(0:axis(1)%NF),                      axis(1)%ys(0:axis(1)%NF), 1 )
!  call fft( NK, dz(0:NK-1,0), axis(1)%NF,                      axis(1)%zc(0:axis(1)%NF),                      axis(1)%zs(0:axis(1)%NF), 1 )
!  
!  call fft( NK, dx(0:NK-1,1), axis(1)%NF, + il(0:axis(1)%NF) * axis(1)%xs(0:axis(1)%NF), - il(0:axis(1)%NF) * axis(1)%xc(0:axis(1)%NF), 1 )
!  call fft( NK, dy(0:NK-1,1), axis(1)%NF, + il(0:axis(1)%NF) * axis(1)%ys(0:axis(1)%NF), - il(0:axis(1)%NF) * axis(1)%yc(0:axis(1)%NF), 1 )
!  call fft( NK, dz(0:NK-1,1), axis(1)%NF, + il(0:axis(1)%NF) * axis(1)%zs(0:axis(1)%NF), - il(0:axis(1)%NF) * axis(1)%zc(0:axis(1)%NF), 1 )
!  
!  call fft( NK, dx(0:NK-1,2), axis(1)%NF, - im(0:axis(1)%NF) * axis(1)%xc(0:axis(1)%NF), - im(0:axis(1)%NF) * axis(1)%xs(0:axis(1)%NF), 1 )
!  call fft( NK, dy(0:NK-1,2), axis(1)%NF, - im(0:axis(1)%NF) * axis(1)%yc(0:axis(1)%NF), - im(0:axis(1)%NF) * axis(1)%ys(0:axis(1)%NF), 1 )
!  call fft( NK, dz(0:NK-1,2), axis(1)%NF, - im(0:axis(1)%NF) * axis(1)%zc(0:axis(1)%NF), - im(0:axis(1)%NF) * axis(1)%zs(0:axis(1)%NF), 1 )
!  
!  call fft( NK, dx(0:NK-1,3), axis(1)%NF, - in(0:axis(1)%NF) * axis(1)%xs(0:axis(1)%NF), + in(0:axis(1)%NF) * axis(1)%xc(0:axis(1)%NF), 1 )
!  call fft( NK, dy(0:NK-1,3), axis(1)%NF, - in(0:axis(1)%NF) * axis(1)%ys(0:axis(1)%NF), + in(0:axis(1)%NF) * axis(1)%yc(0:axis(1)%NF), 1 )
!  call fft( NK, dz(0:NK-1,3), axis(1)%NF, - in(0:axis(1)%NF) * axis(1)%zs(0:axis(1)%NF), + in(0:axis(1)%NF) * axis(1)%zc(0:axis(1)%NF), 1 )
!  
!  SALLOCATE( a, (0:NK-1), zero )
!  SALLOCATE( b, (0:NK-1), zero )
!  SALLOCATE( c, (0:NK-1), zero )
!  SALLOCATE( d, (0:NK-1), zero )
!  
!  a(0:NK-1) = dy(0:NK-1,1) * dz(0:NK-1,2) - dy(0:NK-1,2) * dz(0:NK-1,1)
!  b(0:NK-1) = dz(0:NK-1,1) * dx(0:NK-1,2) - dz(0:NK-1,2) * dx(0:NK-1,1)
!  c(0:NK-1) = dx(0:NK-1,1) * dy(0:NK-1,2) - dx(0:NK-1,2) * dy(0:NK-1,1)
!  
!  d(0:NK-1) = a(0:NK-1)**2 + b(0:NK-1)**2 + c(0:NK-1)**2
!  
!  SALLOCATE( tau, (0:NK-1), zero )
!  
!  tau(0:NK-1) = ( dx(0:NK-1,3) * a(0:NK-1) + dy(0:NK-1,3) * b(0:NK-1) + dz(0:NK-1,3) * c(0:NK-1) ) / d(0:NK-1)
!  
!  call fft( NK, tau(0:NK-1), axis(1)%NF, axis(1)%tc(0:axis(1)%NF), axis(1)%ts(0:axis(1)%NF), 0 )
!  
!  write(ounit,'("rdaxis  : " 10x " : integrated torsion =",es13.5," ;")') pi2 * axis(1)%tc(0)
!  
!  DALLOCATE( il )
!  DALLOCATE( im )
!  DALLOCATE( in )
!
!  DALLOCATE( dx )
!  DALLOCATE( dy )
!  DALLOCATE( dz )
!
!  DALLOCATE( a  )
!  DALLOCATE( b  )
!  DALLOCATE( c  )
!  DALLOCATE( d  )
!  
!  DALLOCATE( tau )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

  return
  
end subroutine rdaxis

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

subroutine knotxx( isurf, srho, teta, zeta, ax, at, az, xx, xs, xt, xz, xtt, xtz, xzz, v1, v2, w1, w2 )
  
  use globals, only : zero, one, two, half, pi2, small, myid, ounit, sqrtmachprec, &
                      axis, &
                      minorrad, ellipticity, nrotate, zetaoff
  
  implicit none
  
  include "mpif.h"
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  INTEGER              :: isurf
  REAL                 :: srho, teta, zeta, ax(1:3), at(1:3), az(1:3), xx(1:3), xs(1:3), xt(1:3), xz(1:3), v1(1:3), v2(1:3), w1(1:3), w2(1:3)
  REAL                 :: xtt(1:3), xtz(1:3), xzz(1:3), z1(1:3), z2(1:3)
  
  INTEGER              :: ierr, p(1:3), q(1:3), mm, isum
  REAL                 :: cqz, sqz, cpz, spz, RR(0:3), ZZ(0:3), x0(1:3), x1(1:3), x2(1:3), x3(1:3), x4(1:3), arg, darg
  REAL                 :: a0, a1, a2, a3, b0, b1, b2, carg, sarg, cost, sint
  REAL                 :: tt(1:3), td(1:3), dd(1:3), dt(1:3), xa, ya, za, xb, yb, zb, xc, yc, zc, nn(1:3), nz(1:3), bb(1:3), bz(1:3), nzz(1:3), bzz(1:3)
  REAL                 :: fa, fb, fc
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

  CHECK( readsrf, .not.allocated(axis(isurf)%xc), illegal )
  CHECK( readsrf, .not.allocated(axis(isurf)%xs), illegal )
  CHECK( readsrf, .not.allocated(axis(isurf)%yc), illegal )
  CHECK( readsrf, .not.allocated(axis(isurf)%ys), illegal )
  CHECK( readsrf, .not.allocated(axis(isurf)%zc), illegal )
  CHECK( readsrf, .not.allocated(axis(isurf)%zs), illegal )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  cost = cos(teta) ; sint = sin(teta)
  
  x0(1:3) = (/ axis(isurf)%xc(0), axis(isurf)%yc(0), axis(isurf)%zc(0) /)

  x1(1:3) = (/ zero             , zero             , zero              /)
  x2(1:3) = (/ zero             , zero             , zero              /)
  x3(1:3) = (/ zero             , zero             , zero              /)
  x4(1:3) = (/ zero             , zero             , zero              /)
  
  do mm = 1, axis(isurf)%NF ; carg = cos(mm*zeta) ; sarg = sin(mm*zeta)
   x0(1:3) = x0(1:3) + ( + (/ axis(isurf)%xc(mm), axis(isurf)%yc(mm), axis(isurf)%zc(mm) /) * carg &
                         + (/ axis(isurf)%xs(mm), axis(isurf)%ys(mm), axis(isurf)%zs(mm) /) * sarg ) 
   x1(1:3) = x1(1:3) + ( - (/ axis(isurf)%xc(mm), axis(isurf)%yc(mm), axis(isurf)%zc(mm) /) * sarg &
                         + (/ axis(isurf)%xs(mm), axis(isurf)%ys(mm), axis(isurf)%zs(mm) /) * carg ) * mm   
   x2(1:3) = x2(1:3) + ( - (/ axis(isurf)%xc(mm), axis(isurf)%yc(mm), axis(isurf)%zc(mm) /) * carg &
                         - (/ axis(isurf)%xs(mm), axis(isurf)%ys(mm), axis(isurf)%zs(mm) /) * sarg ) * mm**2
   x3(1:3) = x3(1:3) + ( + (/ axis(isurf)%xc(mm), axis(isurf)%yc(mm), axis(isurf)%zc(mm) /) * sarg &
                         - (/ axis(isurf)%xs(mm), axis(isurf)%ys(mm), axis(isurf)%zs(mm) /) * carg ) * mm**3
   x4(1:3) = x4(1:3) + ( + (/ axis(isurf)%xc(mm), axis(isurf)%yc(mm), axis(isurf)%zc(mm) /) * carg &
                         + (/ axis(isurf)%xs(mm), axis(isurf)%ys(mm), axis(isurf)%zs(mm) /) * sarg ) * mm**4
  enddo
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

  ax(1:3) = x0(1:3)
  at(1:3) = zero
  az(1:3) = x1(1:3)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  a0      = sqrt( x1(1)*x1(1) + x1(2)*x1(2) + x1(3)*x1(3) )
  
  CHECK( readsrf, abs(a0).lt.sqrtmachprec, divide by zero )
  
  a1      =     ( x1(1)*x2(1) + x1(2)*x2(2) + x1(3)*x2(3) ) / a0          
  a2      =     ( x2(1)*x2(1) + x2(2)*x2(2) + x2(3)*x2(3) &
              +   x1(1)*x3(1) + x1(2)*x3(2) + x1(3)*x3(3) - a1 * a1 ) / a0
  a3      = ( two * sum( x3(1:3)*x2(1:3) ) + sum( x2(1:3)*x3(1:3) ) + sum( x1(1:3)*x4(1:3) ) - 2 * a1 * a2 - a2 * a1 ) / a0

  tt(1:3) = ( x1(1:3)                                                                                 ) / a0 
  td(1:3) = ( x2(1:3)                                                                  - tt(1:3) * a1 ) / a0 ! curvature;
  dd(1:3) = ( x3(1:3) - td(1:3) * a1                          - tt(1:3) * a2           - td(1:3) * a1 ) / a0 
  dt(1:3) = ( x4(1:3) - two * ( dd(1:3) * a1 - td(1:3) * a2 ) - td(1:3) * a2 - tt * a3 - dd(1:3) * a1 ) / a0 ! please check; 12 Nov 17;
  
  xa = x2(1) - tt(1) * a1 ; xb = x3(1) - td(1) * a1 - tt(1) * a2 ; xc = x4(1) - dd(1) * a1 - two * td(1) * a2 - tt(1) * a3
  ya = x2(2) - tt(2) * a1 ; yb = x3(2) - td(2) * a1 - tt(2) * a2 ; yc = x4(2) - dd(2) * a1 - two * td(2) * a2 - tt(2) * a3
  za = x2(3) - tt(3) * a1 ; zb = x3(3) - td(3) * a1 - tt(3) * a2 ; zc = x4(3) - dd(3) * a1 - two * td(3) * a2 - tt(3) * a3
  
  fa = sqrt( xa**2+ya**2+za**2 ) ; fb = ( xa*xb+ya*yb+za*zb ) / fa ; fc = ( xb*xb+yb*yb+zb*zb + xa*xc+ya*yc+za*zc - fb*fb ) / fa
  
  CHECK( readsrf, abs(fa).lt.sqrtmachprec, divide by zero )
  
  b0 = ( fa                               ) / a0 ! magnitude of tangent;
  b1 = ( fb - b0 * a1                     ) / a0
  b2 = ( fc - b1 * a1 - b0 * a2 - b1 * a2 ) / a0

  CHECK( readsrf, abs(b0).lt.sqrtmachprec, divide by zero )
  
   nn(1:3) = ( td(1:3)                                              ) / b0
   nz(1:3) = ( dd(1:3) - nn(1:3) * b1                               ) / b0
  nzz(1:3) = ( dt(1:3) - nz(1:3) * b1 - nn(1:3) * b2 - nz(1:3) * b1 ) / b0
  
  !bb(1:3) = (/ tt(2)*nn(3)-tt(3)*nn(2), tt(3)*nn(1)-tt(1)*nn(3), tt(1)*nn(2)-tt(2)*nn(1) /)

   isum = 0 ; call cross( isum, tt(1:3), nn(1:3), bb(1:3) ) ! 12 Nov 17;

  !bz(1:3) = (/ td(2)*nn(3)-td(3)*nn(2), td(3)*nn(1)-td(1)*nn(3), td(1)*nn(2)-td(2)*nn(1) /) &
  !        + (/ tt(2)*nz(3)-tt(3)*nz(2), tt(3)*nz(1)-tt(1)*nz(3), tt(1)*nz(2)-tt(2)*nz(1) /)

   isum = 0 ; call cross( isum, td(1:3),  nn(1:3),  bz(1:3) )
   isum = 1 ; call cross( isum, tt(1:3),  nz(1:3),  bz(1:3) ) ! take care that cross accumulates; 12 Nov 17;
   
   isum = 0 ; call cross( isum, dd(1:3),  nn(1:3), bzz(1:3) )
   isum = 1 ; call cross( isum, td(1:3),  nz(1:3), bzz(1:3) ) ! take care that cross accumulates; 12 Nov 17;
   isum = 1 ; call cross( isum, td(1:3),  nz(1:3), bzz(1:3) ) ! take care that cross accumulates; 12 Nov 17;
   isum = 1 ; call cross( isum, tt(1:3), nzz(1:3), bzz(1:3) ) ! take care that cross accumulates; 12 Nov 17;
  
!#ifdef DEBUG
!  
!  ;                 ;  arg = (nrotate*half) * zeta + zetaoff - axistc(0) * zeta
!  ;                 ; darg = (nrotate*half)                  - axistc(0)
!  do mm = 1, axisNF ;  arg =  arg - ( axistc(mm) * sin(mm*zeta) - axists(mm) * cos(mm*zeta) ) / mm ! integrated torsion;
!   ;                ; darg = darg - ( axistc(mm) * cos(mm*zeta) + axists(mm) * sin(mm*zeta) )      ! 
!  enddo
!  
!#else
  
  ;                 ;  arg = (nrotate*half) * zeta + zetaoff
  ;                 ; darg = (nrotate*half)                 
  
!#endif

  v1(1:3)  =     cos(arg) *  nn(1:3) + sin(arg) *  bb(1:3)
  
  w1(1:3)  = ( - sin(arg) *  nn(1:3) + cos(arg) *  bb(1:3) ) * darg    &
           + (   cos(arg) *  nz(1:3) + sin(arg) *  bz(1:3) )

  z1(1:3)  = ( - cos(arg) *  nn(1:3) - sin(arg) *  bb(1:3) ) * darg**2 &
           + ( - sin(arg) *  nz(1:3) + cos(arg) *  bz(1:3) ) * darg    &
           + ( - sin(arg) *  nz(1:3) + cos(arg) *  bz(1:3) ) * darg    &
           + (   cos(arg) * nzz(1:3) + sin(arg) * bzz(1:3) )

! v1 = - v1 ; w1 = - w1 ; z1 = - z1 ! DEBUGGING TEST; 12 Nov 17;

  v2(1:3)  =   - sin(arg) *  nn(1:3) + cos(arg) *  bb(1:3)

  w2(1:3)  = ( - cos(arg) *  nn(1:3) - sin(arg) *  bb(1:3) ) * darg &
           + ( - sin(arg) *  nz(1:3) + cos(arg) *  bz(1:3) )

  z2(1:3)  = ( + sin(arg) *  nn(1:3) - cos(arg) *  bb(1:3) ) * darg**2 &
           + ( - cos(arg) *  nz(1:3) - sin(arg) *  bz(1:3) ) * darg    &
           + ( - cos(arg) *  nz(1:3) - sin(arg) *  bz(1:3) ) * darg    &
           + ( - sin(arg) * nzz(1:3) + cos(arg) * bzz(1:3) )

  xx(1:3)  = ax(1:3) + srho * minorrad * (   ellipticity * cost * v1(1:3) + sint * v2(1:3) )

  xs(1:3)  =                  minorrad * (   ellipticity * cost * v1(1:3) + sint * v2(1:3) )
  xt(1:3)  =           srho * minorrad * ( - ellipticity * sint * v1(1:3) + cost * v2(1:3) )
  xz(1:3)  = az(1:3) + srho * minorrad * (   ellipticity * cost * w1(1:3) + sint * w2(1:3) )

  xtt(1:3) =           srho * minorrad * ( - ellipticity * cost * v1(1:3) - sint * v2(1:3) )
  xtz(1:3) =           srho * minorrad * ( - ellipticity * sint * w1(1:3) + cost * w2(1:3) )
  xzz(1:3) = x2(1:3) + srho * minorrad * (   ellipticity * cost * z1(1:3) + sint * z2(1:3) )
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

!#ifdef DEBUG
!  write(ounit,'("knotxx  : " 10x " : xx(1:3), nn(1:3), bb(1:3) =",9f9.4)') xx(1:3), nn(1:3), bb(1:3)
!  FATAL(knotxx  , .true., just checking)
!#endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

  return
  
end subroutine knotxx

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine fft( NN, xi, NF, xc, xs, iflag )
  
  use globals, only : zero, one, two, half, myid
  
  implicit none
  
  include "mpif.h"
  
  INTEGER, intent(in) :: NN, NF, iflag
  REAL                :: xi(0:NN-1)
  REAL                :: xc(0:NF), xs(0:NF)
  
  INTEGER             :: ierr, ifail
  REAL                :: xm(0:NN-1), rwk(0:NN-1), sqrtN

  select case( iflag )
   
  case( 0 )
   
   sqrtN = one / sqrt(NN*one)

   xm(0:NN-1) = xi(0:NN-1)
  
   ifail = 0 ; call C06FAF( xm(0:NN-1) , NN , rwk(0:NN-1) , ifail )
   
   xc(0) = xm(0) * sqrtN ; xc(1:NF) = xm(1:NF) * two * sqrtN ; xs(1:NF) = - xm(NN-1:NN-NF:-1) * two * sqrtN
   
  case( 1 )

   sqrtN = sqrt(NN*one)
   
   xi(0:NN-1) = zero
   
   xi(0) = xc(0) * sqrtN ; xi(1:NF) = xc(1:NF) * sqrtN * half ; xi(NN-1:NN-NF:-1) = - xs(1:NF) * sqrtN * half

   ifail = 0 ; call C06GBF( xi(0:NN-1) , NN , ifail )
   ifail = 0 ; call C06EBF( xi(0:NN-1) , NN , ifail )
   
  case default
   
   FATAL( fft, .true., illegal iflag )
   
  end select
  
  return
  
end subroutine fft

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine cross( isum, aa, bb, cc )
  
  use globals, only : 
  
  implicit none
  
  include "mpif.h"
  
  INTEGER :: isum
  REAL    :: aa(1:3), bb(1:3), cc(1:3)
  
  select case( isum )
  case( 0 ) ; cc(1:3) =           (/ aa(2)*bb(3)-aa(3)*bb(2), aa(3)*bb(1)-aa(1)*bb(3), aa(1)*bb(2)-aa(2)*bb(1) /)   
  case( 1 ) ; cc(1:3) = cc(1:3) + (/ aa(2)*bb(3)-aa(3)*bb(2), aa(3)*bb(1)-aa(1)*bb(3), aa(1)*bb(2)-aa(2)*bb(1) /)   
  end select
  
  return
  
end subroutine cross

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine inttorsion( NF, xc, xs, yc, ys, zc, zs, integratedtorsion )
  
  use globals, only : zero, pi2, ounit
  
  implicit none
  
  include "mpif.h"

  INTEGER              :: NF
  REAL                 :: xc(0:NF), xs(0:NF), yc(0:NF), ys(0:NF), zc(0:NF), zs(0:NF), integratedtorsion
  
  INTEGER              :: mm, NK, astat
  INTEGER, allocatable :: il(:), im(:), in(:)
  REAL                 :: tc(0:NF), ts(0:NF)
  REAL   , allocatable :: dx(:,:), dy(:,:), dz(:,:), a(:), b(:), c(:), d(:), tau(:)
  
  SALLOCATE( il, (0:NF), (/ ( mm   , mm = 0, NF ) /) )
  SALLOCATE( im, (0:NF), (/ ( mm**2, mm = 0, NF ) /) )
  SALLOCATE( in, (0:NF), (/ ( mm**3, mm = 0, NF ) /) )
  
  NK = 3 * 4 * NF ! discrete resolution for resolving torsion;
  
  SALLOCATE( dx, (0:NK-1,0:3), zero )
  SALLOCATE( dy, (0:NK-1,0:3), zero )
  SALLOCATE( dz, (0:NK-1,0:3), zero )
  
  call fft( NK, dx(0:NK-1,0), NF,                      xc(0:NF),                      xs(0:NF), 1 )
  call fft( NK, dy(0:NK-1,0), NF,                      yc(0:NF),                      ys(0:NF), 1 )
  call fft( NK, dz(0:NK-1,0), NF,                      zc(0:NF),                      zs(0:NF), 1 )
  
  call fft( NK, dx(0:NK-1,1), NF, + il(0:NF) * xs(0:NF), - il(0:NF) * xc(0:NF), 1 )
  call fft( NK, dy(0:NK-1,1), NF, + il(0:NF) * ys(0:NF), - il(0:NF) * yc(0:NF), 1 )
  call fft( NK, dz(0:NK-1,1), NF, + il(0:NF) * zs(0:NF), - il(0:NF) * zc(0:NF), 1 )
  
  call fft( NK, dx(0:NK-1,2), NF, - im(0:NF) * xc(0:NF), - im(0:NF) * xs(0:NF), 1 )
  call fft( NK, dy(0:NK-1,2), NF, - im(0:NF) * yc(0:NF), - im(0:NF) * ys(0:NF), 1 )
  call fft( NK, dz(0:NK-1,2), NF, - im(0:NF) * zc(0:NF), - im(0:NF) * zs(0:NF), 1 )
  
  call fft( NK, dx(0:NK-1,3), NF, - in(0:NF) * xs(0:NF), + in(0:NF) * xc(0:NF), 1 )
  call fft( NK, dy(0:NK-1,3), NF, - in(0:NF) * ys(0:NF), + in(0:NF) * yc(0:NF), 1 )
  call fft( NK, dz(0:NK-1,3), NF, - in(0:NF) * zs(0:NF), + in(0:NF) * zc(0:NF), 1 )
  
  SALLOCATE( a, (0:NK-1), zero )
  SALLOCATE( b, (0:NK-1), zero )
  SALLOCATE( c, (0:NK-1), zero )
  SALLOCATE( d, (0:NK-1), zero )
  
  a(0:NK-1) = dy(0:NK-1,1) * dz(0:NK-1,2) - dy(0:NK-1,2) * dz(0:NK-1,1)
  b(0:NK-1) = dz(0:NK-1,1) * dx(0:NK-1,2) - dz(0:NK-1,2) * dx(0:NK-1,1)
  c(0:NK-1) = dx(0:NK-1,1) * dy(0:NK-1,2) - dx(0:NK-1,2) * dy(0:NK-1,1)
  
  d(0:NK-1) = a(0:NK-1)**2 + b(0:NK-1)**2 + c(0:NK-1)**2
  
  SALLOCATE( tau, (0:NK-1), zero )
  
  tau(0:NK-1) = ( dx(0:NK-1,3) * a(0:NK-1) + dy(0:NK-1,3) * b(0:NK-1) + dz(0:NK-1,3) * c(0:NK-1) ) / d(0:NK-1)
  
  call fft( NK, tau(0:NK-1), NF, tc(0:NF), ts(0:NF), 0 )
  
  integratedtorsion = pi2 * tc(0)

  write(ounit,'("rdaxis  : " 10x " : integrated torsion =",es13.5," ;")') integratedtorsion
  
  DALLOCATE( il )
  DALLOCATE( im )
  DALLOCATE( in )
  
  DALLOCATE( dx )
  DALLOCATE( dy )
  DALLOCATE( dz )
  
  DALLOCATE( a  )
  DALLOCATE( b  )
  DALLOCATE( c  )
  DALLOCATE( d  )
  
  DALLOCATE( tau )
  
  return
  
end subroutine inttorsion

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
