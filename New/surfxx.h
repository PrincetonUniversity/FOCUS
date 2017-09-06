subroutine surfxx( teta, zeta, xx, xt, xz, bn )
  
  use globals, only : zero, half, pi2, myid, ounit, runit, surffile, IsQuiet, IsSymmetric, &
                      Nfou, Nfp, NBnf, bim, bin, Bnim, Bnin, Rbc, Rbs, Zbc, Zbs, Bnc, Bns,  &
                      Nteta, Nzeta, surf, Npc, discretefactor
  
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
  
  if(NBnf >  0) then


   
     do jj = 0, Nzeta-1 ; zeta = ( jj + half ) * pi2 / (Nzeta*Nfp)
        do ii = 0, Nteta-1 ; teta = ( ii + half ) * pi2 / Nteta
           do imn = 1, NBnf
              arg = Bnim(imn) * teta - Bnin(imn) * zeta
              bn = bn + Bnc(imn)*cos(arg) + Bns(imn)*sin(arg)
           enddo
        enddo
     enddo

  endif
 
  return
  
end subroutine surfxx

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

