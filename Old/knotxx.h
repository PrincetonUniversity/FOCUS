
!title (brief) ! Description.

!latex \briefly{Extended description.}

!latex \calledby{\link{notopt}, \link{plassf}, \link{rdknot}, \link{windsf}}
!latex \calls{\link{}}

!latex \tableofcontents

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine knotxx( aa, teta, zeta, ax, at, az, xx, xt, xz )
  
  use kmodule, only : zero, one, pi2, small, myid, ounit, &
                      Itopology, knotNF, knotsurf, ellipticity, nrotate, zetaoff, &
                      xkc, xks, ykc, yks, zkc, zks
  
  implicit none
  
  include "mpif.h"
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  REAL                 :: aa, teta, zeta, ax(1:3), at(1:3), az(1:3), xx(1:3), xt(1:3), xz(1:3)
  
  INTEGER              :: ierr, p(1:3), q(1:3), mm
  REAL                 :: rr, rt, rz
  REAL                 :: cqz, sqz, cpz, spz, x0(1:3), x1(1:3), x2(1:3), x3(1:3), a0, a1, a2, b0, b1, carg, sarg, ctz
  REAL                 :: tt(1:3), td(1:3), dd(1:3), xa, ya, za, ff, nn(1:3), nz(1:3), bb(1:3), bz(1:3), v1(1:3), v2(1:3), w1(1:3), w2(1:3), arg
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  select case( Itopology )
   
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
   
!  FATAL( knotxx, abs(knotphase).gt.small, need to revise phase )
   
   x0(1:3) = (/ xkc(0), ykc(0), zkc(0) /)
   x1(1:3) = (/ zero  , zero  , zero   /)
   x2(1:3) = (/ zero  , zero  , zero   /)
   x3(1:3) = (/ zero  , zero  , zero   /)
   
   if( aa.gt.zero ) then ! will need additional derivatives to construct normal and binormal; 14 Apr 16;
    
    do mm = 1, knotNF ; carg = cos(mm*zeta) ; sarg = sin(mm*zeta)
     x0(1:3) = x0(1:3) + ( + (/ xkc(mm), ykc(mm), zkc(mm) /) * carg + (/ xks(mm), yks(mm), zks(mm) /) * sarg ) 
     x1(1:3) = x1(1:3) + ( - (/ xkc(mm), ykc(mm), zkc(mm) /) * sarg + (/ xks(mm), yks(mm), zks(mm) /) * carg ) * mm   
     x2(1:3) = x2(1:3) + ( - (/ xkc(mm), ykc(mm), zkc(mm) /) * carg - (/ xks(mm), yks(mm), zks(mm) /) * sarg ) * mm**2
     x3(1:3) = x3(1:3) + ( + (/ xkc(mm), ykc(mm), zkc(mm) /) * sarg - (/ xks(mm), yks(mm), zks(mm) /) * carg ) * mm**3
    enddo
    
   else
    
    do mm = 1, knotNF ; carg = cos(mm*zeta) ; sarg = sin(mm*zeta)
     x0(1:3) = x0(1:3) + ( + (/ xkc(mm), ykc(mm), zkc(mm) /) * carg + (/ xks(mm), yks(mm), zks(mm) /) * sarg )       
     x1(1:3) = x1(1:3) + ( - (/ xkc(mm), ykc(mm), zkc(mm) /) * sarg + (/ xks(mm), yks(mm), zks(mm) /) * carg ) * mm   
    !x2(1:3) = x2(1:3) + ( - (/ xkc(mm), ykc(mm), zkc(mm) /) * carg - (/ xks(mm), yks(mm), zks(mm) /) * sarg ) * mm**2
    !x3(1:3) = x3(1:3) + ( + (/ xkc(mm), ykc(mm), zkc(mm) /) * sarg - (/ xks(mm), yks(mm), zks(mm) /) * carg ) * mm**3
    enddo
    
   endif ! end of if( aa.gt.zero ) ; 14 Apr 16;
   
  case default
   
   FATAL( knotxx, .true., selected Itopology is not supported )
   
  end select
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  ax(1:3) = x0(1:3)                                                        ! Nov 12 15;
  at(1:3) = zero                                                           ! Nov 12 15;
  az(1:3) = x1(1:3)                                                        ! Nov 12 15;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( aa.gt.zero ) then
   
   a0      = sqrt( x1(1)*x1(1) + x1(2)*x1(2) + x1(3)*x1(3) )                                   ! Nov 12 15;
   a1      =     ( x1(1)*x2(1) + x1(2)*x2(2) + x1(3)*x2(3) ) / a0                              ! Nov 12 15;
   a2      =     ( x2(1)*x2(1) + x2(2)*x2(2) + x2(3)*x2(3)   &
               +   x1(1)*x3(1) + x1(2)*x3(2) + x1(3)*x3(3) - a1 * a1 ) / a0                    ! Nov 12 15;

   tt(1:3) =   x1(1:3)                                                / a0                     ! Nov 12 15;
   td(1:3) = ( x2(1:3) - tt(1:3) * a1                               ) / a0                     ! Nov 12 15;
   dd(1:3) = ( x3(1:3) - td(1:3) * a1 - tt(1:3) * a2 - td(1:3) * a1 ) / a0                     ! Nov 12 15;
   
   xa = ( x2(1) - tt(1) * a1 ) ! Nov 12 15;
   ya = ( x2(2) - tt(2) * a1 ) ! Nov 12 15;
   za = ( x2(3) - tt(3) * a1 ) ! Nov 12 15;
   
   ff = sqrt( xa**2 + ya**2 + za**2 ) ! Nov 12 15;
   
   b0 = ff / a0 ! Nov 12 15;
   
   b1 = ( ( xa * ( x3(1) - td(1) * a1 - tt(1) * a2 ) &
          + ya * ( x3(2) - td(2) * a1 - tt(2) * a2 ) &
          + za * ( x3(3) - td(3) * a1 - tt(3) * a2 ) ) / ff - b0 * a1 ) / a0 ! Nov 12 15;

   nn(1:3) =   td(1:3)                  / b0                                                   ! Nov 12 15;
   nz(1:3) = ( dd(1:3) - nn(1:3) * b1 ) / b0                                                   ! Nov 12 15;
   
   bb(1:3) = (/ tt(2)*nn(3)-tt(3)*nn(2), tt(3)*nn(1)-tt(1)*nn(3), tt(1)*nn(2)-tt(2)*nn(1) /)   ! Nov 12 15;
   bz(1:3) = (/ td(2)*nn(3)-td(3)*nn(2), td(3)*nn(1)-td(1)*nn(3), td(1)*nn(2)-td(2)*nn(1) /) &
           + (/ tt(2)*nz(3)-tt(3)*nz(2), tt(3)*nz(1)-tt(1)*nz(3), tt(1)*nz(2)-tt(2)*nz(1) /)   ! Nov 12 15;

   arg = nrotate * zeta + zetaoff

   v1(1:3) =     cos(arg) * nn(1:3) + sin(arg) * bb(1:3)

   w1(1:3) = ( - sin(arg) * nn(1:3) + cos(arg) * bb(1:3) ) * nrotate &
           + (   cos(arg) * nz(1:3) + sin(arg) * bz(1:3) )

   v2(1:3) =   - sin(arg) * nn(1:3) + cos(arg) * bb(1:3)

   w2(1:3) = ( - cos(arg) * nn(1:3) - sin(arg) * bb(1:3) ) * nrotate &
           + ( - sin(arg) * nz(1:3) + cos(arg) * bz(1:3) )

   xx(1:3) = ax(1:3) + aa * knotsurf * (   ellipticity * cos(teta) * v1(1:3) + sin(teta) * v2(1:3) )

   xt(1:3) = at(1:3) + aa * knotsurf * ( - ellipticity * sin(teta) * v1(1:3) + cos(teta) * v2(1:3) )

   xz(1:3) = az(1:3) + aa * knotsurf * (   ellipticity * cos(teta) * w1(1:3) + sin(teta) * w2(1:3) )

!  xx(1:3) = ax(1:3) + aa * knotsurf * ( (   ellipticity * cos(teta) * cos(zeta) + sin(teta) * sin(zeta) ) * nn(1:3)   &
!                                      + (   ellipticity * cos(teta) * sin(zeta) - sin(teta) * cos(zeta) ) * bb(1:3) )   ! 31 May 17;

!  xt(1:3) = at(1:3) + aa * knotsurf * ( ( - ellipticity * sin(teta) * cos(zeta) + cos(teta) * sin(zeta) ) * nn(1:3)   &
!                                      + ( - ellipticity * sin(teta) * sin(zeta) - cos(teta) * cos(zeta) ) * bb(1:3) )   ! 31 May 17;

!  xz(1:3) = az(1:3) + aa * knotsurf * ( ( - ellipticity * cos(teta) * sin(zeta) + sin(teta) * cos(zeta) ) * nn(1:3)   &
!                                      + ( + ellipticity * cos(teta) * cos(zeta) + sin(teta) * sin(zeta) ) * bb(1:3) ) & 
!                    + aa * knotsurf * ( (   ellipticity * cos(teta) * cos(zeta) + sin(teta) * sin(zeta) ) * nz(1:3)   &
!                                      + (   ellipticity * cos(teta) * sin(zeta) - sin(teta) * cos(zeta) ) * bz(1:3) )


!  xx(1:3)  = ax(1:3) + aa * knotsurf * (   ellipticity * cos(teta) * nn(1:3) + sin(teta) * bb(1:3) )   ! 31 May 17;

!  xt(1:3)  = at(1:3) + aa * knotsurf * ( - ellipticity * sin(teta) * nn(1:3) + cos(teta) * bb(1:3) )   ! 31 May 17;

!  xz(1:3)  = az(1:3) + aa * knotsurf * ( + ellipticity * sin(teta) * nn(1:3) - cos(teta) * bb(1:3) ) & ! 31 May 17;
!                     + aa * knotsurf * (   ellipticity * cos(teta) * nz(1:3) + sin(teta) * bz(1:3) )   ! 31 May 17;

!  rr = aa * ( knotsurf + ellipticity * cos( teta - zeta )        )
!  rt = aa * (     zero - ellipticity * sin( teta - zeta ) * (+1) )
!  rz = aa * (     zero - ellipticity * sin( teta - zeta ) * (-1) )
   
!  xx(1:3) = ax(1:3) + rr * (   cos(teta) * nn(1:3) - sin(teta) * bb(1:3) ) ! aa is minor radius;


!  xt(1:3) = at(1:3) + rt * (   cos(teta) * nn(1:3) - sin(teta) * bb(1:3) ) &
!                    + rr * ( - sin(teta) * nn(1:3) - cos(teta) * bb(1:3) )

!  xz(1:3) = az(1:3) + rz * (   cos(teta) * nn(1:3) - sin(teta) * bb(1:3) ) &
!                    + rr * (   cos(teta) * nz(1:3) - sin(teta) * bz(1:3) )
   
!  xx(1:3) = x0(1:3) + aa * (   cos(teta) * nn(1:3) - sin(teta) * bb(1:3) ) ! aa is minor radius;
!  xt(1:3) = zero    + aa * ( - sin(teta) * nn(1:3) - cos(teta) * bb(1:3) )
!  xz(1:3) = x1(1:3) + aa * (   cos(teta) * nz(1:3) - sin(teta) * bz(1:3) )
   
  else
   
   xx(1:3) = zero
   xt(1:3) = zero
   xz(1:3) = zero
   
  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  return
  
end subroutine knotxx

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
