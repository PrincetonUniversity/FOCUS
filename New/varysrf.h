!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

!title (force) ! Force.

!latex \briefly{Force.}

!latex \calledby{\link{}}
!latex \calls{\link{initial}, \link{rdsurf}, \link{rdcoils}, 
!latex \link{saving}, \link{solvers}}.

!latex \tableofcontents

!latex \subsection{Brief introduction}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

subroutine varysrf( Ndof, xdof, fdof )
  
  use globals, only : zero, two, three, half, pi2, sqrtmachprec, myid, ncpu, ounit, tstart, &
                      Ncoils, coil, cmt, smt, &
                      Nt, Nz, Ns, surf, discretecurve, deltateta, deltazeta, discretesurface, &
                      totlengt, Tfluxave, Bdotnsqd, weight_bnorm
                     !wspectral, pspectral, Mdof, tdof, &
                     !target_tflux, weight_length, weight_tflux, weight_bnorm
  
  implicit none
  
  include "mpif.h"
  
  INTEGER               :: Ndof
  REAL                  :: xdof(1:Ndof), fdof(1:Ndof)
  
  INTEGER               :: ierr, icoil, mm, ii, jj, kk, idof, NF, icmodNc, isurf, isum, astat, il, iu, jl, ju
  REAL                  :: BB(0:Ns-1,1:3), CC, lBdotnsqd, ndotdx, BnRnH(1:3), dFds(1:3), Bs, Bu, Bv, gss, gst, gsz, gtt, gtz, gzz, BsgradRn(1:3), RsgradBn(1:3)
  REAL                  :: rr(1:3), rs, ru, rv, RRs(1:3), RRu(1:3), RRv(1:3), dist(1:5), tnow, told
  REAL   , allocatable  :: dxdtau(:,:)
  CHARACTER             :: packorunpack
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  told = MPI_WTIME()
  
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
     
     surf(isurf)%BB(1,ii,jj) = surf(isurf)%BB(1,ii,jj) + sum( BB(0:Ns-1,1) ) * CC ! * discretecurve ; 22 Nov 17;
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
  
  lBdotnsqd = sum( surf(isurf)%Bn(0:Nt-1,0:Nz-1) * surf(isurf)%Bn(0:Nt-1,0:Nz-1) * surf(isurf)%ds(0:Nt-1,0:Nz-1) ) * discretesurface * half
 
! write(ounit,'("varysrf : " 10x " : return ; 3 ;")') ; return
 
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

  SALLOCATE( dxdtau, (0:Ns-1,1:3), zero )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  do icoil = 1, Ncoils ! 28 Nov 17;
   
   do kk = 0, Ns-1 ! line segments describing curve; 22 Nov 17;
    
    dFds(1:3) = zero ! initialize summation; 22 Nov 17;
    
    do ii = 0, Nt-1 ; il = ii-1 ; iu = ii+1

     if( il.lt. 0 ) il = Nt-1
     if( iu.eq.Nt ) iu = 0

     do jj = 0, Nz-1 ; jl = jj-1 ; ju = jj+1
      
      if( jl.lt. 0 ) jl = Nz-1
      if( ju.eq.Nz ) ju = 0


     !ndotdx = sum( surf(isurf)%nn(1:3,ii,jj) * surf(2)%xx(1:3,ii,jj) )

      ndotdx = sum( surf(isurf)%nn(1:3,ii,jj) * ( surf(2)%xx(1:3,ii,jj) - surf(1)%xx(1:3,ii,jj) ) )

      
!     Bs = sum( surf(isurf)%BB(1:3,ii,jj) * surf(isurf)%gs(1:3,ii,jj) )
!     Bu = sum( surf(isurf)%BB(1:3,ii,jj) * surf(isurf)%gu(1:3,ii,jj) )
!     Bv = sum( surf(isurf)%BB(1:3,ii,jj) * surf(isurf)%gv(1:3,ii,jj) )

!     gss = surf(isurf)%guv(4,ii,jj) * surf(isurf)%guv(6,ii,jj) - surf(isurf)%guv(5,ii,jj) * surf(isurf)%guv(5,ii,jj) ! Jacobian factors cancel; 16 Nov 17;
!     gst = surf(isurf)%guv(5,ii,jj) * surf(isurf)%guv(3,ii,jj) - surf(isurf)%guv(2,ii,jj) * surf(isurf)%guv(6,ii,jj)
!     gsz = surf(isurf)%guv(2,ii,jj) * surf(isurf)%guv(5,ii,jj) - surf(isurf)%guv(4,ii,jj) * surf(isurf)%guv(3,ii,jj)

!     BsgradRn(1:3) = ( Bu - Bs * gst/gss ) * ( coil(icoil)%Rn(1:3,iu,jj,kk) - coil(icoil)%Rn(1:3,il,jj,kk) ) / ( two * deltateta ) & 
!                   + ( Bv - Bs * gsz/gss ) * ( coil(icoil)%Rn(1:3,ii,ju,kk) - coil(icoil)%Rn(1:3,ii,jl,kk) ) / ( two * deltazeta )
            
     !Bs = sum( surf(isurf)%BB(1:3,ii,jj) * surf(isurf)%gs(1:3,ii,jj) )
      Bu = sum( surf(isurf)%BB(1:3,ii,jj) * surf(isurf)%xu(1:3,ii,jj) )
      Bv = sum( surf(isurf)%BB(1:3,ii,jj) * surf(isurf)%xv(1:3,ii,jj) )

      gtt = surf(isurf)%guv(4,ii,jj) ; gtz = surf(isurf)%guv(5,ii,jj) ; gzz = surf(isurf)%guv(6,ii,jj) ; gss = gtt * gzz - gtz * gtz

     !gss = surf(isurf)%guv(4,ii,jj) * surf(isurf)%guv(6,ii,jj) - surf(isurf)%guv(5,ii,jj) * surf(isurf)%guv(5,ii,jj) ! Jacobian factors cancel; 16 Nov 17;
     !gst = surf(isurf)%guv(5,ii,jj) * surf(isurf)%guv(3,ii,jj) - surf(isurf)%guv(2,ii,jj) * surf(isurf)%guv(6,ii,jj)
     !gsz = surf(isurf)%guv(2,ii,jj) * surf(isurf)%guv(5,ii,jj) - surf(isurf)%guv(4,ii,jj) * surf(isurf)%guv(3,ii,jj)

      BsgradRn(1:3) = ( ( Bu * gzz - Bv * gtz ) * ( coil(icoil)%Rn(1:3,iu,jj,kk) - coil(icoil)%Rn(1:3,il,jj,kk) ) / ( two * deltateta )   & 
                      + ( Bv * gtt - Bu * gtz ) * ( coil(icoil)%Rn(1:3,ii,ju,kk) - coil(icoil)%Rn(1:3,ii,jl,kk) ) / ( two * deltazeta ) ) / gss

      rr(1:3) = surf(isurf)%xx(1:3,ii,jj) - coil(icoil)%xx(1:3,kk)

      rs = sum( rr(1:3) * surf(isurf)%gs(1:3,ii,jj) )
      ru = sum( rr(1:3) * surf(isurf)%gu(1:3,ii,jj) )
      rv = sum( rr(1:3) * surf(isurf)%gv(1:3,ii,jj) )
      
      dist(2) = sum((rr(1:3))**2) ; dist(1) = sqrt(dist(2)) ; dist(3) = dist(2)*dist(1) ; dist(5) = dist(3)*dist(2)

      RRs(1:3) = three * rr(1:3) * rs / dist(5) - surf(isurf)%gs(1:3,ii,jj) / dist(3)
      RRu(1:3) = three * rr(1:3) * ru / dist(5) - surf(isurf)%gu(1:3,ii,jj) / dist(3)
      RRv(1:3) = three * rr(1:3) * rv / dist(5) - surf(isurf)%gv(1:3,ii,jj) / dist(3)
      
      RsgradBn(1:3) = ( RRu(1:3) - RRs(1:3) * gst/gss ) * ( surf(isurf)%Bn(iu,jj) - surf(isurf)%Bn(il,jj) ) / ( two * deltateta ) & 
                    + ( RRv(1:3) - RRs(1:3) * gsz/gss ) * ( surf(isurf)%Bn(ii,ju) - surf(isurf)%Bn(ii,jl) ) / ( two * deltazeta )


      BnRnH(1:3) = surf(isurf)%Bn(ii,jj) * coil(icoil)%Rn(1:3,ii,jj,kk) * surf(isurf)%HH(ii,jj)
      

      dFds(1:3) = dFds(1:3) + ndotdx * ( RsgradBn(1:3) + BsgradRn(1:3) + BnRnH(1:3) ) * surf(isurf)%ds(ii,jj)
      

     enddo ! end of do jj; 16 Nov 17;
    enddo ! end of do jj; 16 Nov 17;
    
    dFds(1:3) = dFds(1:3) * discretesurface
    
    isum = 0 ; call cross( isum, coil(icoil)%xt(1:3,kk), dFdS(1:3), dxdtau(kk,1:3) )
    
   enddo ! end of do kk; 22 Nov 17;
   
   ;             ; idof = coil(icoil)%gdof
   ;  mm = 0     ; idof = idof + 1 ; fdof(idof) = sum( dxdtau(0:Ns-1,1) * cmt(0:Ns-1,mm) ) ! * discretecurve
   ;             ; idof = idof + 1 ; fdof(idof) = sum( dxdtau(0:Ns-1,2) * cmt(0:Ns-1,mm) )
   ;             ; idof = idof + 1 ; fdof(idof) = sum( dxdtau(0:Ns-1,3) * cmt(0:Ns-1,mm) )
   do mm = 1, NF ; idof = idof + 1 ; fdof(idof) = sum( dxdtau(0:Ns-1,1) * cmt(0:Ns-1,mm) )
    ;            ; idof = idof + 1 ; fdof(idof) = sum( dxdtau(0:Ns-1,1) * smt(0:Ns-1,mm) )
    ;            ; idof = idof + 1 ; fdof(idof) = sum( dxdtau(0:Ns-1,2) * cmt(0:Ns-1,mm) )
    ;            ; idof = idof + 1 ; fdof(idof) = sum( dxdtau(0:Ns-1,2) * smt(0:Ns-1,mm) )
    ;            ; idof = idof + 1 ; fdof(idof) = sum( dxdtau(0:Ns-1,3) * cmt(0:Ns-1,mm) )
    ;            ; idof = idof + 1 ; fdof(idof) = sum( dxdtau(0:Ns-1,3) * smt(0:Ns-1,mm) )
   enddo

  enddo ! end of do icoil; 16 Nov 17;
  
  fdof(1:Ndof) = fdof(1:Ndof) * weight_bnorm

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

  DALLOCATE( dxdtau )

  tnow = MPI_WTIME()

  write(ounit,'("varysrf : ",f10.1," : |B.n|^2 =",es23.15," = ",es23.15," ; time =",f8.2,"s ;" )') tnow-tstart, lBdotnsqd, Bdotnsqd(0), tnow-told

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  return
  
end subroutine varysrf

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
