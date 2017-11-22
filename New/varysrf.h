!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

!title (force) ! Force.

!latex \briefly{Force.}

!latex \calledby{\link{}}
!latex \calls{\link{initial}, \link{rdsurf}, \link{rdcoils}, 
!latex \link{saving}, \link{solvers}}.

!latex \tableofcontents

!latex \subsection{Brief introduction}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

subroutine varysrf( Ndof, xdof )
  
  use globals, only : zero, two, three, half, pi2, sqrtmachprec, myid, ncpu, ounit, &
                      Ncoils, coil, cmt, smt, &
                      Nt, Nz, Ns, surf, discretecurve, deltateta, deltazeta, discretesurface, &
                      totlengt, Tfluxave, Bdotnsqd, &
                      wspectral, pspectral, Mdof, tdof, &
                      target_tflux, weight_length, weight_tflux, weight_bnorm
  
  implicit none
  
  include "mpif.h"
  
  INTEGER               :: Ndof
  REAL                  :: xdof(1:Ndof)
  
  INTEGER               :: ierr, icoil, mm, ii, jj, kk, idof, NF, icmodNc, isurf, isum
  REAL                  :: BB(0:Ns-1,1:3), CC, lBdotnsqd, ndotdx, BnRnH(1:3), dFds(1:3), Bs, Bu, Bv, gss, gst, gsz, BsgradRn(1:3), RsgradBn(1:3)
  REAL                  :: rr(1:3), rs, ru, rv, RRs(1:3), RRu(1:3), RRv(1:3), dist(1:5), dxdtau(1:3)
! REAL                  :: tx(0:Ns-1), ty(0:Ns-1), tz(0:Ns-1), ax(0:Ns-1), ay(0:Ns-1), az(0:Ns-1), txtx(0:Ns-1), txax(0:Ns-1), denom(0:Ns-1)
  CHARACTER             :: packorunpack
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  write(ounit,'("varysrf : " 10x " : return ;")') ; return
  
  packorunpack = 'U' ; call packdof( Ndof, xdof(1:Ndof), packorunpack )
  
  isurf = 1
  
  CHECK( varysrf, .not.allocated( surf(isurf)%BB ), fatal )
  
  surf(isurf)%BB(1:3,0:Nt-1,0:Nz-1) = zero ! initialize summation for components of total magnetic field; 16 Nov 17;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  do icoil = 1, Ncoils
   
   NF = coil(icoil)%NF ; CC = coil(icoil)%I
   
   do jj = 0, Nz-1
    do ii = 0, Nt-1
     
     call brfield( icoil, ii, jj, Ns, BB(0:Ns-1,1:3) )
     
     surf(isurf)%BB(1,ii,jj) = surf(isurf)%BB(1,ii,jj) + sum( BB(0:Ns-1,1) ) * CC
     surf(isurf)%BB(2,ii,jj) = surf(isurf)%BB(2,ii,jj) + sum( BB(0:Ns-1,2) ) * CC
     surf(isurf)%BB(3,ii,jj) = surf(isurf)%BB(3,ii,jj) + sum( BB(0:Ns-1,3) ) * CC
     
    enddo ! end do ii;
   enddo ! end do jj;
   
  enddo ! end of do icoil = 1, Ncoils;
  
  do jj = 0, Nz-1
   do ii = 0, Nt-1
    
    surf(isurf)%Bn(ii,jj) = surf(isurf)%BB(1,ii,jj) * surf(isurf)%nn(1,ii,jj) &
                          + surf(isurf)%BB(2,ii,jj) * surf(isurf)%nn(2,ii,jj) &
                          + surf(isurf)%BB(3,ii,jj) * surf(isurf)%nn(3,ii,jj)
    
   enddo
  enddo
  
! lBdotnsqd = sum( surf(isurf)%Bn(0:Nt-1,0:Nz-1) * surf(isurf)%Bn(0:Nt-1,0:Nz-1) * surf(isurf)%ds(0:Nt-1,0:Nz-1) ) * discretesurface * half
! write(ounit,'("varysrf : " 10x " : Bdotnsqd=",2es23.16" ; ratio =",f9.5," ;")') lBdotnsqd, Bdotnsqd(0), lBdotnsqd/Bdotnsqd(0) 
! return
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  do icoil = 1, Ncoils
   
   do kk = 0, Ns-1 ! line segments describing curve; 22 Nov 17;
    
    dFds(1:3) = zero ! initialize summation; 22 Nov 17;
    
    do ii = 0, Nt-1
     do jj = 0, Nz-1
      
      ndotdx = sum( surf(isurf)%nn(1:3,ii,jj) * surf(2)%xx(1:3,ii,jj) )

      
      Bs = sum( surf(isurf)%BB(1:3,ii,jj) * surf(isurf)%gs(1:3,ii,jj) )
      Bu = sum( surf(isurf)%BB(1:3,ii,jj) * surf(isurf)%gu(1:3,ii,jj) )
      Bv = sum( surf(isurf)%BB(1:3,ii,jj) * surf(isurf)%gv(1:3,ii,jj) )

      gss = surf(isurf)%guv(4,ii,jj) * surf(isurf)%guv(6,ii,jj) - surf(isurf)%guv(5,ii,jj) * surf(isurf)%guv(5,ii,jj) ! Jacobian factors cancel; 16 Nov 17;
      gst = surf(isurf)%guv(5,ii,jj) * surf(isurf)%guv(3,ii,jj) - surf(isurf)%guv(2,ii,jj) * surf(isurf)%guv(6,ii,jj) ! Jacobian factors cancel; 16 Nov 17;
      gsz = surf(isurf)%guv(2,ii,jj) * surf(isurf)%guv(5,ii,jj) - surf(isurf)%guv(4,ii,jj) * surf(isurf)%guv(3,ii,jj) ! Jacobian factors cancel; 16 Nov 17;
      
      BsgradRn(1:3) = ( Bu - Bs * gst/gss ) * ( coil(icoil)%Rn(1:3,ii+1,jj  ,kk) - coil(icoil)%Rn(1:3,ii-1,jj  ,kk) ) / ( two * deltateta ) & 
                    + ( Bv - Bs * gsz/gss ) * ( coil(icoil)%Rn(1:3,ii  ,jj+1,kk) - coil(icoil)%Rn(1:3,ii-1,jj-1,kk) ) / ( two * deltazeta )
      
      rr(1:3) = surf(isurf)%xx(1:3,ii,jj)-coil(icoil)%xx(1:3,kk)

      rs = sum( ( surf(isurf)%xx(1:3,ii,jj)-coil(icoil)%xx(1:3,kk) ) * surf(isurf)%gs(1:3,ii,jj) )
      ru = sum( ( surf(isurf)%xx(1:3,ii,jj)-coil(icoil)%xx(1:3,kk) ) * surf(isurf)%gu(1:3,ii,jj) )
      rv = sum( ( surf(isurf)%xx(1:3,ii,jj)-coil(icoil)%xx(1:3,kk) ) * surf(isurf)%gv(1:3,ii,jj) )
      
      dist(2) = sum((rr(1:3))**2) ; dist(1) = sqrt(dist(2)) ; dist(3) = dist(2)*dist(1) ; dist(5) = dist(3)*dist(2)

      RRs(1:3) = three * rr(1:3) * rs / dist(5) - surf(isurf)%gs(1:3,ii,jj) / dist(3)
      RRu(1:3) = three * rr(1:3) * ru / dist(5) - surf(isurf)%gu(1:3,ii,jj) / dist(3)
      RRv(1:3) = three * rr(1:3) * rv / dist(5) - surf(isurf)%gv(1:3,ii,jj) / dist(3)
      
      RsgradBn(1:3) = ( RRu(1:3) - RRs(1:3) * gst/gss ) * ( surf(isurf)%Bn(ii+1,jj  ) - surf(isurf)%Bn(ii-1,jj  ) ) / ( two * deltateta ) & 
                    + ( RRv(1:3) - RRs(1:3) * gsz/gss ) * ( surf(isurf)%Bn(ii  ,jj+1) - surf(isurf)%Bn(ii  ,jj-1) ) / ( two * deltateta )


      BnRnH(1:3) = surf(isurf)%Bn(ii,jj) * coil(icoil)%Rn(1:3,ii,jj,kk) * surf(isurf)%HH(ii,jj)
      

      dFds(1:3) = dFds(1:3) + ndotdx * ( RsgradBn(1:3) + BsgradRn(1:3) + BnRnH(1:3) )
      

     enddo ! end of do jj; 16 Nov 17;
    enddo ! end of do jj; 16 Nov 17;
    
    dFds(1:3) = dFds(1:3) * discretesurface
    
    isum = 0 ; call cross( isum, coil(icoil)%xt(1:3,kk), dFdS(1:3), dxdtau(1:3) ) !   now take cross product with tangent to coil; 22 Nov 17;
    
   enddo ! end of do kk; 22 Nov 17;
   
  enddo ! end of do icoil; 16 Nov 17;
  
  return

  do icoil = 1, Ncoils
   
   do kk = 0, Ns-1
    
    do ii = 0, Nt-1
     do jj = 0, Nz-1
      
!     B^s = sum( surf(isurf)%BB(1:3,ii,jj) * surf(isurf)%gs(1:3,ii,jj) )
!     B^u = sum( surf(isurf)%BB(1:3,ii,jj) * surf(isurf)%gu(1:3,ii,jj) )
!     B^v = sum( surf(isurf)%BB(1:3,ii,jj) * surf(isurf)%gv(1:3,ii,jj) )
      
!     r^s = sum( ( surf(isurf)%xx(1:3,ii,jj)-coil(icoil)%xx(1:3,kk) ) * surf(isurf)%gs(1:3,ii,jj) )
!     r^u = sum( ( surf(isurf)%xx(1:3,ii,jj)-coil(icoil)%xx(1:3,kk) ) * surf(isurf)%gu(1:3,ii,jj) )
!     r^v = sum( ( surf(isurf)%xx(1:3,ii,jj)-coil(icoil)%xx(1:3,kk) ) * surf(isurf)%gv(1:3,ii,jj) )
      
!     gss = surf(isurf)%guv(4,ii,jj) * surf(isurf)%guv(6,ii,jj) - surf(isurf)%guv(5,ii,jj) * surf(isurf)%guv(5,ii,jj) ! Jacobian factors cancel; 16 Nov 17;
!     gst = surf(isurf)%guv(5,ii,jj) * surf(isurf)%guv(3,ii,jj) - surf(isurf)%guv(2,ii,jj) * surf(isurf)%guv(6,ii,jj) ! Jacobian factors cancel; 16 Nov 17;
!     gsz = surf(isurf)%guv(2,ii,jj) * surf(isurf)%guv(5,ii,jj) - surf(isurf)%guv(4,ii,jj) * surf(isurf)%guv(3,ii,jj) ! Jacobian factors cancel; 16 Nov 17;
      
!      = ( B^u - B^s * gst / gss ) * ( coil(icoil)%Rn(1:3,ii+1,jj  ) - coil(icoil)%Rn(1:3,ii-1,jj  ) ) / ( two * deltateta ) & 
!      + ( B^v - B^s * gsz / gss ) * ( coil(icoil)%Rn(1:3,ii  ,jj+1) - coil(icoil)%Rn(1:3,ii-1,jj-1) ) / ( two * deltazeta )
      
!      R^s(1:3) = 3 * rr(1:3) * r^s / dist(5) - surf(isurf)%gs(1:3,ii,jj) / dist(3)
!      R^u(1:3) = 3 * rr(1:3) * r^u / dist(5) - surf(isurf)%gu(1:3,ii,jj) / dist(3)
!      R^v(1:3) = 3 * rr(1:3) * r^v / dist(5) - surf(isurf)%gv(1:3,ii,jj) / dist(3)

!     = ( R^u(1:3) - R^s(1:3) * gst / gss ) * ( surf(isurf)%Bn(ii+1,jj  ) - surf(isurf)%Bn(ii-1,jj  ) ) / ( two * deltateta ) & 
!     + ( R^z(1:3) - R^s(1:3) * gsz / gss ) * ( surf(isurf)%Bn(ii  ,jj+1) - surf(isurf)%Bn(ii  ,jj-1) ) / ( two * deltateta )
      
     enddo ! end of do jj; 16 Nov 17;
    enddo ! end of do ii; 16 Nov 17;
    
   enddo
   
  enddo ! end of do icoil; 16 Nov 17;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  return
  
end subroutine varysrf

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
