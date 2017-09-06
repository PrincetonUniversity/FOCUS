
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
