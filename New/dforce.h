!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

!title (force) ! Force.

!latex \briefly{Force.}

!latex \calledby{\link{}}
!latex \calls{\link{initial}, \link{rdsurf}, \link{rdcoils}, 
!latex \link{saving}, \link{solvers}}.

!latex \tableofcontents

!latex \subsection{Brief introduction}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

subroutine dforce( Ndof, xdof, ff, fdof )
  
  use globals, only : zero, half, pi2, sqrtmachprec, myid, ncpu, ounit, &
                      Ncoils, coil, cmt, smt, &
                      Nt, Nz, Ns, surf, discretecurve, deltatheta, discretesurface, &
                      totlengt, Tfluxave, Bdotnsqd, &
                      wspectral, pspectral, Mdof, tdof, &
                      target_length, target_tflux, weight_length, weight_tflux, weight_bnorm
  
  implicit none
  
  include "mpif.h"
  
  INTEGER               :: Ndof
  REAL                  :: xdof(1:Ndof), ff, fdof(1:Ndof)
  
  INTEGER               :: ierr, icoil, mm, ii, jj, kk, idof, NF, icmodNc
  REAL                  :: dFdx(0:Ns,0:3), dBdx(0:Ns,0:3), CC
  REAL                  :: length(0:Ns), tx(0:Ns), ty(0:Ns), tz(0:Ns), ax(0:Ns), ay(0:Ns), az(0:Ns), txtx(0:Ns), txax(0:Ns), denom(0:Ns)
  CHARACTER             :: packorunpack
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  packorunpack = 'U' ; call packdof( Ndof, xdof(1:Ndof), packorunpack )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  do icoil = 1, Ncoils
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
   
   if( myid.ne.modulo(icoil-1,ncpu) ) cycle ! parallelization loop;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

   if( coil(icoil)%Lfree.eq.0 ) cycle ! geometry of coil is not allowed to vary;
   
   CHECK( dforce , coil(icoil)%gdof.lt.0, fatal )
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

   CHECK( dforce , .not.allocated( coil(icoil)%dL ), fatal )
   CHECK( dforce , .not.allocated( coil(icoil)%dT ), fatal )
   CHECK( dforce , .not.allocated( coil(icoil)%dB ), fatal )
   
   coil(icoil)%dL(              0:Ndof) = zero ! it is probably not required to set these to zero;
   coil(icoil)%dT(       0:Nz-1,0:Ndof) = zero
   coil(icoil)%dB(0:Nt-1,0:Nz-1,0:Ndof) = zero
   
   NF = coil(icoil)%NF ; CC = coil(icoil)%I ! shorthand; Fourier resolution; coil current;
   
   do jj = 0, Nz-1
    do ii = 0, Nt-1
     
     call abfield( icoil, ii, jj, Ns, dFdx(0:Ns,0:3), dBdx(0:Ns,0:3) ) ! compute toroidal flux and B.n produced by i-th coil on surface;
     
     ;             ; idof =        0         ; coil(icoil)%dT(   jj,idof) = coil(icoil)%dT(jj,idof) + sum( dFdx(1:Ns,0)                ) * CC
     ;             ;                         ; coil(icoil)%dB(ii,jj,idof) =                           sum( dBdx(1:Ns,0)                ) * CC
     ;             ; idof = coil(icoil)%gdof
     ;  mm = 0     ; idof = idof + 1         ; coil(icoil)%dT(   jj,idof) = coil(icoil)%dT(jj,idof) + sum( dFdx(1:Ns,1) * cmt(1:Ns,mm) ) * CC
     ;             ;                         ; coil(icoil)%dB(ii,jj,idof) =                           sum( dBdx(1:Ns,1) * cmt(1:Ns,mm) ) * CC
     ;             ; idof = idof + 1         ; coil(icoil)%dT(   jj,idof) = coil(icoil)%dT(jj,idof) + sum( dFdx(1:Ns,2) * cmt(1:Ns,mm) ) * CC
     ;             ;                         ; coil(icoil)%dB(ii,jj,idof) =                           sum( dBdx(1:Ns,2) * cmt(1:Ns,mm) ) * CC
     ;             ; idof = idof + 1         ; coil(icoil)%dT(   jj,idof) = coil(icoil)%dT(jj,idof) + sum( dFdx(1:Ns,3) * cmt(1:Ns,mm) ) * CC
     ;             ;                         ; coil(icoil)%dB(ii,jj,idof) =                           sum( dBdx(1:Ns,3) * cmt(1:Ns,mm) ) * CC
     do mm = 1, NF ; idof = idof + 1         ; coil(icoil)%dT(   jj,idof) = coil(icoil)%dT(jj,idof) + sum( dFdx(1:Ns,1) * cmt(1:Ns,mm) ) * CC
      ;            ;                         ; coil(icoil)%dB(ii,jj,idof) =                           sum( dBdx(1:Ns,1) * cmt(1:Ns,mm) ) * CC
      ;            ; idof = idof + 1         ; coil(icoil)%dT(   jj,idof) = coil(icoil)%dT(jj,idof) + sum( dFdx(1:Ns,1) * smt(1:Ns,mm) ) * CC
      ;            ;                         ; coil(icoil)%dB(ii,jj,idof) =                           sum( dBdx(1:Ns,1) * smt(1:Ns,mm) ) * CC
      ;            ; idof = idof + 1         ; coil(icoil)%dT(   jj,idof) = coil(icoil)%dT(jj,idof) + sum( dFdx(1:Ns,2) * cmt(1:Ns,mm) ) * CC
      ;            ;                         ; coil(icoil)%dB(ii,jj,idof) =                           sum( dBdx(1:Ns,2) * cmt(1:Ns,mm) ) * CC
      ;            ; idof = idof + 1         ; coil(icoil)%dT(   jj,idof) = coil(icoil)%dT(jj,idof) + sum( dFdx(1:Ns,2) * smt(1:Ns,mm) ) * CC
      ;            ;                         ; coil(icoil)%dB(ii,jj,idof) =                           sum( dBdx(1:Ns,2) * smt(1:Ns,mm) ) * CC
      ;            ; idof = idof + 1         ; coil(icoil)%dT(   jj,idof) = coil(icoil)%dT(jj,idof) + sum( dFdx(1:Ns,3) * cmt(1:Ns,mm) ) * CC
      ;            ;                         ; coil(icoil)%dB(ii,jj,idof) =                           sum( dBdx(1:Ns,3) * cmt(1:Ns,mm) ) * CC
      ;            ; idof = idof + 1         ; coil(icoil)%dT(   jj,idof) = coil(icoil)%dT(jj,idof) + sum( dFdx(1:Ns,3) * smt(1:Ns,mm) ) * CC
      ;            ;                         ; coil(icoil)%dB(ii,jj,idof) =                           sum( dBdx(1:Ns,3) * smt(1:Ns,mm) ) * CC
     enddo
     ;             ; idof = idof + 1         ; coil(icoil)%dT(   jj,idof) = coil(icoil)%dT(jj,idof) + sum( dFdx(1:Ns,0)                )
     ;             ;                         ; coil(icoil)%dB(ii,jj,idof) =                           sum( dBdx(1:Ns,0)                )
     
    enddo ! end do ii;
    
    coil(icoil)%dT(jj,0:Ndof) = coil(icoil)%dT(jj,0:Ndof) * discretecurve * deltatheta

   enddo ! end do jj;
   
   tx(1:Ns) = coil(icoil)%xt(1:Ns) ; ax(1:Ns) = coil(icoil)%xa(1:Ns) ! tangent to coil (shorthand);
   ty(1:Ns) = coil(icoil)%yt(1:Ns) ; ay(1:Ns) = coil(icoil)%ya(1:Ns) ! tangent to coil (shorthand);
   tz(1:Ns) = coil(icoil)%zt(1:Ns) ; az(1:Ns) = coil(icoil)%za(1:Ns) ! tangent to coil (shorthand);
   
   txtx(1:Ns) = tx(1:Ns)*tx(1:Ns) + ty(1:Ns)*ty(1:Ns) + tz(1:Ns)*tz(1:Ns) ; length(1:Ns) = sqrt(txtx(1:Ns)) ; denom(1:Ns) = length(1:Ns) * txtx(1:Ns)
   txax(1:Ns) = tx(1:Ns)*ax(1:Ns) + ty(1:Ns)*ay(1:Ns) + tz(1:Ns)*az(1:Ns)
   
#ifdef DEBUG
   do kk = 1, Ns
    FATAL( dforce , abs(denom(kk)).lt.sqrtmachprec, fatal )
   enddo
#endif

   ;             ; idof =        0         ; coil(icoil)%dL(idof) = sum(   length(1:Ns)                                                                 )
   ;             ; idof = coil(icoil)%gdof
   ;  mm = 0     ; idof = idof + 1         ; coil(icoil)%dL(idof) = sum( ( txax(1:Ns) * tx(1:Ns) - txtx(1:Ns) * ax(1:Ns) ) * cmt(1:Ns,mm) / denom(1:Ns) )
   ;             ; idof = idof + 1         ; coil(icoil)%dL(idof) = sum( ( txax(1:Ns) * ty(1:Ns) - txtx(1:Ns) * ay(1:Ns) ) * cmt(1:Ns,mm) / denom(1:Ns) )
   ;             ; idof = idof + 1         ; coil(icoil)%dL(idof) = sum( ( txax(1:Ns) * tz(1:Ns) - txtx(1:Ns) * az(1:Ns) ) * cmt(1:Ns,mm) / denom(1:Ns) )
   do mm = 1, NF ; idof = idof + 1         ; coil(icoil)%dL(idof) = sum( ( txax(1:Ns) * tx(1:Ns) - txtx(1:Ns) * ax(1:Ns) ) * cmt(1:Ns,mm) / denom(1:Ns) )
    ;            ; idof = idof + 1         ; coil(icoil)%dL(idof) = sum( ( txax(1:Ns) * tx(1:Ns) - txtx(1:Ns) * ax(1:Ns) ) * smt(1:Ns,mm) / denom(1:Ns) )
    ;            ; idof = idof + 1         ; coil(icoil)%dL(idof) = sum( ( txax(1:Ns) * ty(1:Ns) - txtx(1:Ns) * ay(1:Ns) ) * cmt(1:Ns,mm) / denom(1:Ns) )
    ;            ; idof = idof + 1         ; coil(icoil)%dL(idof) = sum( ( txax(1:Ns) * ty(1:Ns) - txtx(1:Ns) * ay(1:Ns) ) * smt(1:Ns,mm) / denom(1:Ns) )
    ;            ; idof = idof + 1         ; coil(icoil)%dL(idof) = sum( ( txax(1:Ns) * tz(1:Ns) - txtx(1:Ns) * az(1:Ns) ) * cmt(1:Ns,mm) / denom(1:Ns) )
    ;            ; idof = idof + 1         ; coil(icoil)%dL(idof) = sum( ( txax(1:Ns) * tz(1:Ns) - txtx(1:Ns) * az(1:Ns) ) * smt(1:Ns,mm) / denom(1:Ns) )
   enddo
   ;             ; idof = idof + 1         ; coil(icoil)%dL(idof) = zero ! length does not depend on current;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
   
  enddo ! end of do icoil = 1, Ncoils;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  do icoil = 1, Ncoils ; icmodNc = modulo(icoil-1,ncpu)
   RlBCAST( coil(icoil)%dL(              0:Ndof),       (Ndof+1), icmodNc )
   RlBCAST( coil(icoil)%dT(       0:Nz-1,0:Ndof),    Nz*(Ndof+1), icmodNc )
   RlBCAST( coil(icoil)%dB(0:Nt-1,0:Nz-1,0:Ndof), Nt*Nz*(Ndof+1), icmodNc )
  enddo
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  CHECK( dforce , .not.allocated(surf%dL), fatal )
  CHECK( dforce , .not.allocated(surf%dT), fatal )
  CHECK( dforce , .not.allocated(surf%dB), fatal )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  surf%dL(              0:Ndof) = zero
  surf%dT(       0:Nz-1,0:Ndof) = zero
  surf%dB(0:Nt-1,0:Nz-1,0:Ndof) = zero
  do icoil = 1, Ncoils
   surf%dL(              0:Ndof) = surf%dL(              0:Ndof) + coil(icoil)%dL(              0:Ndof)
   surf%dT(       0:Nz-1,0:Ndof) = surf%dT(       0:Nz-1,0:Ndof) + coil(icoil)%dT(       0:Nz-1,0:Ndof)
   surf%dB(0:Nt-1,0:Nz-1,0:Ndof) = surf%dB(0:Nt-1,0:Nz-1,0:Ndof) + coil(icoil)%dB(0:Nt-1,0:Nz-1,0:Ndof)
  enddo
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  CHECK( dforce , Ncoils.eq.0, divide by zero )
  CHECK( dforce , Nz    .eq.0, divide by zero )
  
  CHECK( dforce , .not.allocated(totlengt), fatal )
  CHECK( dforce , .not.allocated(Tfluxave), fatal )
  CHECK( dforce , .not.allocated(Bdotnsqd), fatal )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  do idof = 0, Ndof
   totlengt(idof) =      surf%dL(              idof)                                                       * discretecurve / Ncoils
   Tfluxave(idof) = sum( surf%dT(       0:Nz-1,idof)                                                     ) / Nz
   Bdotnsqd(idof) = sum( surf%dB(0:Nt-1,0:Nz-1,idof) * surf%dB(0:Nt-1,0:Nz-1,0) * surf%ds(0:Nt-1,0:Nz-1) ) * discretesurface
  enddo
  
  Bdotnsqd(0) = Bdotnsqd(0) * half
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  CHECK( dforce, target_length.lt.sqrtmachprec, illegal )

  ff = zero ; fdof(1:Ndof) = zero
  
! ff           =                    exp(totlengt(0)) * weight_length / exp(target_length)
! fdof(1:Ndof) = totlengt(1:Ndof) * exp(totlengt(0)) * weight_length / exp(target_length)

  ff           =                                    exp(totlengt(0)/target_length) * weight_length
  fdof(1:Ndof) = (totlengt(1:Ndof)/target_length) * exp(totlengt(0)/target_length) * weight_length
  
  ;ff         = ff         + weight_tflux * sum( (surf%dT(0:Nz-1,   0)-target_tflux) * (surf%dT(0:Nz-1,0)-target_tflux) ) * half
  do idof = 1, Ndof
   fdof(idof) = fdof(idof) + weight_tflux * sum( (surf%dT(0:Nz-1,idof)             ) * (surf%dT(0:Nz-1,0)-target_tflux) )
  enddo
  
  ff           = ff           + weight_bnorm * Bdotnsqd(0     )
  fdof(1:Ndof) = fdof(1:Ndof) + weight_bnorm * Bdotnsqd(1:Ndof)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

  fdof(1:Ndof) = fdof(1:Ndof) + wspectral * sum(tdof(1:Ndof)*Mdof(1:Ndof)) * tdof(1:Ndof)
  
! ff = totlengt(0) ; fdof(1:Ndof) = totlengt(1:Ndof) ! debugging;
! ff = Tfluxave(0) ; fdof(1:Ndof) = Tfluxave(1:Ndof) ! debugging;
! ff = Bdotnsqd(0) ; fdof(1:Ndof) = Bdotnsqd(1:Ndof) ! debugging;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

  return
  
end subroutine dforce

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
