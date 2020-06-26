
!title (knot) ! Reads knot parameters from file. (Redundant?)

!latex \briefly{The parameters describing a knotted loop are read from file.}

!latex \calledby{\link{notopt}}
!latex \calls{\link{knotxx}}

!latex \tableofcontents

!latex \subsection{overview}
!latex \bi
!latex \item[1.] \red{This routine need to be updated}.
!latex \item[2.] A knot is a closed curve in three-dimensional space.
!latex \item[3.] The parameters describing the knot, namely \internal{xkc}, \internal{yks} and 
!latex \internal{zks}, are read from \verb+knot+.
!latex \item[4.] If \inputvar{kspring} $> 0$, the knot is ``relaxed'' to minimize the ``spring energy''
!latex  using 
!latex           \nag{http://www.nag.co.uk/numeric/FL/manual19/pdf/D02/d02bhf_fl19.pdf}{D02BHF}.
!latex \ei

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!!

subroutine rdknot
  
  use globals, only : dp, zero, one, half, ten, pi2, sqrtmachprec, myid, ncpu, ounit, runit, &
                      ext, input_surf, MPI_COMM_FOCUS, &
                      NFcoil, knotsurf, knotphase, &
                      xkc, xks, ykc, yks, zkc, zks!, kspring, tauend, itau， 
  
  implicit none
  
  include "mpif.h"
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!!
  
  LOGICAL              :: exist
  INTEGER              :: iostat, astat, ierr, lNFcoil

!  INTEGER, parameter   :: Ndof = 2, Ldf = Ndof, LRR = Ndof*(Ndof+1)/2
!  INTEGER              :: irevcm, mode, ic05pdf
!  REAL                 :: phi, teta, zeta
!  REAL                 :: xtol, factor, xdof(1:Ndof), ff(1:Ndof), df(1:Ldf,1:Ndof), diag(1:Ndof), RR(1:LRR), QTf(1:Ndof), work(1:Ndof,1:4)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!!
  
  inquire( file=trim(ext)//".op.knot", exist=exist )
  if( exist ) then
   if( myid.eq.0 ) then
    open(runit, file=trim(ext)//".op.knot", status="old", action='read', iostat=iostat)  
    write(ounit,'("rdknot : " 10x " : reading ext.op.knot ;")') 
   endif
  else
   inquire( file=trim(ext)//".knot", exist=exist )
   if( exist ) then
    if( myid.eq.0 ) then
     open(runit, file=trim(ext)//".knot", status="old", action='read', iostat=iostat)  
     write(ounit,'("rdcoil : " 10x " : reading ext.knot ; redundant naming convention ;")') 
    endif
   else
    inquire( file="knot", exist=exist )
    if( exist ) then
     if( myid.eq.0 ) then
      open(runit, file="knot", status="old", action='read', iostat=iostat)
      write(ounit,'("rdcoil : " 10x " : reading knot ;")') 
     endif
    endif
   endif
  endif
  
  FATAL( rdknot, .not.exist, input knot does not exist )
  
  SALLOCATE( xkc, (0:NFcoil), zero )
  SALLOCATE( xks, (0:NFcoil), zero )
  SALLOCATE( ykc, (0:NFcoil), zero )
  SALLOCATE( yks, (0:NFcoil), zero )
  SALLOCATE( zkc, (0:NFcoil), zero )
  SALLOCATE( zks, (0:NFcoil), zero )
  
  if( myid.eq.0 ) read( runit, *, iostat=iostat ) lNFcoil, knotphase
  
  IlBCAST( lNFcoil, 1, 0 )
  RlBCAST( knotphase, 1, 0 )
  
! NFcoil = qtorusknot + ptorusknot
!   
! FATAL(rdknot, ptorusknot.gt.qtorusknot, needs attention )
!   
! xkc(ptorusknot) =   Rtorusaxis ; xkc(qtorusknot-ptorusknot) = + half * atorusaxis ; xkc(qtorusknot+ptorusknot) = half * atorusaxis
! yks(ptorusknot) =   Rtorusaxis ; yks(qtorusknot-ptorusknot) = - half * atorusaxis ; yks(qtorusknot+ptorusknot) = half * atorusaxis
! zks(qtorusknot) = - atorusaxis
  
! suffix = 'square' & phase = pi2 / 4.0D
  
! xnc=[ 0.0,  1.767311692237854e+00,  6.828825699356500e-14, -9.444974660873413e-01,  3.381930152590584e-14, -2.645063698291779e-01 ]
! yns=[ 0.0,  6.908266991376877e-02, -2.093433065786243e-12, -9.791350364685059e-01,  7.350465722200106e-13, -2.231226861476898e-01 ]
! zns=[ 0.0, -9.695488214492798e-02, -1.645109395931321e-12,  5.140577629208565e-02, -3.898339984154120e-13, -4.254519343376160e-01 ]
  
! suffix = 'granny' & phase = 0.0D
  
! xkc(0:5) = (/ 0.0,  6.653696894645691e-01, -1.656940383509831e-11, -1.043724060058594e+00,  1.069502906375641e-11, -8.075384050607681e-02 /)
! yks(0:5) = (/ 0.0, -5.815748572349548e-01, -2.619011431337359e-11, -1.002896547317505e+00, -3.263673359343855e-11,  1.082220524549484e-01 /)
! zks(0:5) = (/ 0.0,  4.834449507384875e-11, -1.123250573873520e-01,  5.830259885986067e-11, -6.876130104064941e-01, -3.233802114976925e-12 /)
  
 !write(ounit,'("rdknot : " 10x " : lNFcoil, NFcoil=",2i6," ;")') lNFcoil, NFcoil

  if( myid.eq.0 ) read( runit, *, iostat=iostat ) xkc(0:min(lNFcoil,NFcoil))
  if( myid.eq.0 ) read( runit, *, iostat=iostat ) xks(0:min(lNFcoil,NFcoil))
  if( myid.eq.0 ) read( runit, *, iostat=iostat ) ykc(0:min(lNFcoil,NFcoil))
  if( myid.eq.0 ) read( runit, *, iostat=iostat ) yks(0:min(lNFcoil,NFcoil))
  if( myid.eq.0 ) read( runit, *, iostat=iostat ) zkc(0:min(lNFcoil,NFcoil))
  if( myid.eq.0 ) read( runit, *, iostat=iostat ) zks(0:min(lNFcoil,NFcoil))
  
  if( myid.eq.0 ) close(runit,iostat=iostat)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!!
  
  RlBCAST( xkc(0:NFcoil), NFcoil+1, 0 )
  RlBCAST( xks(0:NFcoil), NFcoil+1, 0 )
  RlBCAST( ykc(0:NFcoil), NFcoil+1, 0 )
  RlBCAST( yks(0:NFcoil), NFcoil+1, 0 )
  RlBCAST( zkc(0:NFcoil), NFcoil+1, 0 )
  RlBCAST( zks(0:NFcoil), NFcoil+1, 0 )
   
  if( myid.eq.0 ) then
   write(ounit,'("rdknot : " 10x " : xkc=",    999es11.03)') xkc(0:NFcoil)
   write(ounit,'("rdknot : " 10x " : xks=",11x,998es11.03)') xks(1:NFcoil)
   write(ounit,'("rdknot : " 10x " : ykc=",    999es11.03)') ykc(0:NFcoil)
   write(ounit,'("rdknot : " 10x " : yks=",11x,999es11.03)') yks(1:NFcoil)
   write(ounit,'("rdknot : " 10x " : zkc=",    999es11.03)') zkc(0:NFcoil)
   write(ounit,'("rdknot : " 10x " : zks=",11x,999es11.03)') zks(1:NFcoil)
  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!!
  
!  if( myid.eq.0 ) then
!
!   SALLOCATE( xsurf, (0:NDcoils,0:NDcoils), zero )
!   SALLOCATE( ysurf, (0:NDcoils,0:NDcoils), zero )
!   SALLOCATE( zsurf, (0:NDcoils,0:NDcoils), zero )
!   
!   xdof(1:Ndof) = (/ zero, zero + ten /)
!   
!   do jj = 0, NDcoils ; phi  = jj * (pi2/Nfp) / NDcoils
!    
!    do ii = 0, NDcoils ; teta = ii * pi2       / NDcoils
!     
!     ic05pdf = 1 ; irevcm = 0 ; xtol = sqrtmachprec ; mode = 1 ; factor = one
!     
!     do
!      
!      call C05PDF( irevcm, &
!                   Ndof, xdof(1:Ndof), ff(1:Ndof), df(1:Ldf,1:Ndof), &
!                   Ldf, xtol, diag(1:Ndof), mode, factor, RR(1:LRR), LRR, QTf(1:Ndof), work(1:Ndof,1:4), ic05pdf )
!
!      zeta = xdof(2) - ten ; call knotxx( knotsurf, xdof(1), zeta, ax(1:3), at(1:3), az(1:3), xx(1:3), xt(1:3), xz(1:3) )
!      
!      rx = sqrt( xx(1)*xx(1) + xx(2)*xx(2) ) ; rt = ( xx(1)*xt(1) + xx(2)*xt(2) ) / rx ; rz = ( xx(1)*xz(1) + xx(2)*xz(2) ) / rx
!      ox = sqrt( ax(1)*ax(1) + ax(2)*ax(2) ) ; ot = ( ax(1)*at(1) + ax(2)*at(2) ) / ox ; oz = ( ax(1)*az(1) + ax(2)*az(2) ) / ox
!      
!      select case( irevcm )
!      case( 0 ) ! if( ic05pdf.ne.0 ) write(ounit,2000) phi , teta, irevcm, ic05pdf, xx(1:3)
!       ;        ; exit
!      case( 1 ) ;
!      case( 2 ) ; ff(1  ) =   xx(1)        * sin(phi ) -   xx(2)           * cos(phi )
!       ;        ; ff(2  ) = ( rx    - ox ) * sin(teta) - ( xx(3) - ax(3) ) * cos(teta)
!      case( 3 ) ; df(1,1) =   xt(1)        * sin(phi ) -   xt(2)           * cos(phi )
!       ;        ; df(1,2) =   xz(1)        * sin(phi ) -   xz(2)           * cos(phi )
!       ;        ; df(2,1) = ( rt    - ot ) * sin(teta) - ( xt(3) - at(3) ) * cos(teta)
!       ;        ; df(2,2) = ( rz    - oz ) * sin(teta) - ( xz(3) - az(3) ) * cos(teta)
!      case default
!       FATAL(rdknot, .true., invalid irevcm )
!      end select
!      
!     enddo ! end of do; Nov 12 15;
!     
!2000 format("rdknot : " 10x " : phi ="es23.15" ; teta="es23.15" ; irevcm="i3" ; ic05pdf="i3" ; xx="3es23.15" ;")
!     
!     xsurf(ii,jj) = xx(1)
!     ysurf(ii,jj) = xx(2)
!     zsurf(ii,jj) = xx(3)
!     
!    enddo ! end of do ii; Nov 12 15;
!   enddo ! end of do jj; Nov 12 15;
!   
!  endif ! end of if( myid.eq.0 ) ; Nov 12 15;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!!
  
  return
  
end subroutine rdknot

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!!

!new knotxx subroutine
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!!

subroutine knotxx( aa, teta, zeta, ax, at, az, xx, xt, xz )
  
  use globals, only : dp, zero, one, pi2, small, myid, ounit, &
                      case_surface, NFcoil, knotphase, &
                      xkc, xks, ykc, yks, zkc, zks, MPI_COMM_FOCUS
  
  implicit none
  
  include "mpif.h"
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!!
  
  REAL                 :: aa, teta, zeta, ax(1:3), at(1:3), az(1:3), xx(1:3), xt(1:3), xz(1:3)
  
  INTEGER              :: ierr, p(1:3), q(1:3), mm
  REAL                 :: cqz, sqz, cpz, spz, RR(0:3), ZZ(0:3), x0(1:3), x1(1:3), x2(1:3), x3(1:3), &
                          a0, a1, a2, b0, b1, carg, sarg
  REAL                 :: tt(1:3), td(1:3), dd(1:3), xa, ya, za, ff, nn(1:3), nd(1:3), bb(1:3), bd(1:3)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!!
  
  select case( case_surface )
   
!  case( 1 ) ! torus knot; Nov 12 15;
!   
!   p(1:3) = (/ ptorusknot, ptorusknot**2, ptorusknot**3 /) ! Nov 12 15;
!   q(1:3) = (/ qtorusknot, qtorusknot**2, qtorusknot**3 /) ! Nov 12 15;
!   
!   cqz = cos( q(1)*zeta ) ; cpz = cos( p(1)*zeta ) ! Nov 12 15;
!   sqz = sin( q(1)*zeta ) ; spz = sin( p(1)*zeta ) ! Nov 12 15;
!   
!   RR(0) = Rtorusaxis + atorusaxis * cqz        ! Nov 12 15;
!   RR(1) =            - atorusaxis * sqz * q(1) ! Nov 12 15;
!   RR(2) =            - atorusaxis * cqz * q(2) ! Nov 12 15;
!   RR(3) =            + atorusaxis * sqz * q(3) ! Nov 12 15;
!   
!   ZZ(0) =            - atorusaxis * sqz        ! Nov 12 15;
!   ZZ(1) =            - atorusaxis * cqz * q(1) ! Nov 12 15;
!   ZZ(2) =            + atorusaxis * sqz * q(2) ! Nov 12 15;
!   ZZ(3) =            + atorusaxis * cqz * q(3) ! Nov 12 15;
!   
!   x0(1:3) = (/ + RR(0) * cpz, + RR(0) * spz, + ZZ(0) /) ! Nov 12 15;
!   
!   x1(1:3) = (/ + RR(1) * cpz, + RR(1) * spz, + ZZ(1) /) & 
!           + (/ - RR(0) * spz, + RR(0) * cpz,   zero  /) * p(1) ! Nov 12 15;
!
!   x2(1:3) = (/ + RR(2) * cpz, + RR(2) * spz, + ZZ(2) /)            & 
!           + (/ - RR(1) * spz, + RR(1) * cpz,   zero  /) * p(1) * 2 &
!           + (/ - RR(0) * cpz, - RR(0) * spz,   zero  /) * p(2)
!
!   x3(1:3) = (/ + RR(3) * cpz, + RR(3) * spz, + ZZ(3) /)            & 
!           + (/ - RR(2) * spz, + RR(2) * cpz,   zero  /) * p(1) * 3 &
!           + (/ - RR(1) * cpz, - RR(1) * spz,   zero  /) * p(2) * 3 &
!           + (/ + RR(0) * spz, - RR(0) * cpz,   zero  /) * p(3)       ! Nov 12 15;
   
  case( 1:2 ) ! arbitrary knot represented using Fourier series; Nov 12 15;
   
   FATAL( knotxx, abs(knotphase).gt.small, need to revise phase )
   
   x0(1:3) = (/ xkc(0), ykc(0), zkc(0) /)
   x1(1:3) = (/ zero  , zero  , zero   /)
   x2(1:3) = (/ zero  , zero  , zero   /)
   x3(1:3) = (/ zero  , zero  , zero   /)
   
   if( aa.gt.zero ) then ! will need additional derivatives to construct normal and binormal; 14 Apr 16;
    
    do mm = 1, NFcoil ; carg = cos(mm*zeta) ; sarg = sin(mm*zeta)
     x0(1:3) = x0(1:3) + ( + (/ xkc(mm), ykc(mm), zkc(mm) /) * carg + &
                             (/ xks(mm), yks(mm), zks(mm) /) * sarg ) 
     x1(1:3) = x1(1:3) + ( - (/ xkc(mm), ykc(mm), zkc(mm) /) * sarg + &
                             (/ xks(mm), yks(mm), zks(mm) /) * carg ) * mm   
     x2(1:3) = x2(1:3) + ( - (/ xkc(mm), ykc(mm), zkc(mm) /) * carg - &
                             (/ xks(mm), yks(mm), zks(mm) /) * sarg ) * mm**2
     x3(1:3) = x3(1:3) + ( + (/ xkc(mm), ykc(mm), zkc(mm) /) * sarg - &
                             (/ xks(mm), yks(mm), zks(mm) /) * carg ) * mm**3
    enddo
    
   else
    
    do mm = 1, NFcoil ; carg = cos(mm*zeta) ; sarg = sin(mm*zeta)
     x0(1:3) = x0(1:3) + ( + (/ xkc(mm), ykc(mm), zkc(mm) /) * carg + &
                             (/ xks(mm), yks(mm), zks(mm) /) * sarg )       
     x1(1:3) = x1(1:3) + ( - (/ xkc(mm), ykc(mm), zkc(mm) /) * sarg + &
                             (/ xks(mm), yks(mm), zks(mm) /) * carg ) * mm   
    !x2(1:3) = x2(1:3) + ( - (/ xkc(mm), ykc(mm), zkc(mm) /) * carg - &
    !                        (/ xks(mm), yks(mm), zks(mm) /) * sarg ) * mm**2
    !x3(1:3) = x3(1:3) + ( + (/ xkc(mm), ykc(mm), zkc(mm) /) * sarg - &
    !                        (/ xks(mm), yks(mm), zks(mm) /) * carg ) * mm**3
    enddo
    
   endif ! end of if( aa.gt.zero ) ; 14 Apr 16;
   
  case default
   
   FATAL( knotxx, .true., selected case_surface is not supported )
   
  end select
  
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
   
  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!!

  return
  
end subroutine knotxx

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!!
