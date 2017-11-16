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
                      target_tflux, weight_length, weight_tflux, weight_bnorm
  
  implicit none
  
  include "mpif.h"
  
  INTEGER               :: Ndof
  REAL                  :: xdof(1:Ndof), ff, fdof(1:Ndof)
  
  INTEGER               :: ierr, icoil, mm, ii, jj, kk, idof, NF, icmodNc
  REAL                  :: dFdx(0:Ns-1,0:3), dBdx(0:Ns-1,0:3), CC
  REAL                  :: length(0:Ns-1), tx(0:Ns-1), ty(0:Ns-1), tz(0:Ns-1), ax(0:Ns-1), ay(0:Ns-1), az(0:Ns-1), txtx(0:Ns-1), txax(0:Ns-1), dd(0:Ns-1)
  CHARACTER             :: packorunpack
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  packorunpack = 'U' ; call packdof( Ndof, xdof(1:Ndof), packorunpack )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  do icoil = 1, Ncoils
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
   
   if( myid.ne.modulo(icoil-1,ncpu) ) cycle ! parallelization loop;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

  !if( coil(icoil)%Lfree.eq.0 ) cycle ! geometry of coil is not allowed to vary; ! THIS IS AN ERROR!!! STILL NEED TO KNOW B_total \cdot \bn; 12 Nov 17;
   
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
     
     call abfield( icoil, ii, jj, Ns, dFdx(0:Ns-1,0:3), dBdx(0:Ns-1,0:3) ) ! compute toroidal flux and B.n produced by i-th coil on surface;
     
     ;             ; idof =        0         ; coil(icoil)%dT(   jj,idof) = coil(icoil)%dT(jj,idof) + sum( dFdx(0:Ns-1,0)                  ) * CC
     ;             ;                         ; coil(icoil)%dB(ii,jj,idof) =                           sum( dBdx(0:Ns-1,0)                  ) * CC
     ;             ; idof = coil(icoil)%gdof
     ;  mm = 0     ; idof = idof + 1         ; coil(icoil)%dT(   jj,idof) = coil(icoil)%dT(jj,idof) + sum( dFdx(0:Ns-1,1) * cmt(0:Ns-1,mm) ) * CC
     ;             ;                         ; coil(icoil)%dB(ii,jj,idof) =                           sum( dBdx(0:Ns-1,1) * cmt(0:Ns-1,mm) ) * CC
     ;             ; idof = idof + 1         ; coil(icoil)%dT(   jj,idof) = coil(icoil)%dT(jj,idof) + sum( dFdx(0:Ns-1,2) * cmt(0:Ns-1,mm) ) * CC
     ;             ;                         ; coil(icoil)%dB(ii,jj,idof) =                           sum( dBdx(0:Ns-1,2) * cmt(0:Ns-1,mm) ) * CC
     ;             ; idof = idof + 1         ; coil(icoil)%dT(   jj,idof) = coil(icoil)%dT(jj,idof) + sum( dFdx(0:Ns-1,3) * cmt(0:Ns-1,mm) ) * CC
     ;             ;                         ; coil(icoil)%dB(ii,jj,idof) =                           sum( dBdx(0:Ns-1,3) * cmt(0:Ns-1,mm) ) * CC
     do mm = 1, NF ; idof = idof + 1         ; coil(icoil)%dT(   jj,idof) = coil(icoil)%dT(jj,idof) + sum( dFdx(0:Ns-1,1) * cmt(0:Ns-1,mm) ) * CC
      ;            ;                         ; coil(icoil)%dB(ii,jj,idof) =                           sum( dBdx(0:Ns-1,1) * cmt(0:Ns-1,mm) ) * CC
      ;            ; idof = idof + 1         ; coil(icoil)%dT(   jj,idof) = coil(icoil)%dT(jj,idof) + sum( dFdx(0:Ns-1,1) * smt(0:Ns-1,mm) ) * CC
      ;            ;                         ; coil(icoil)%dB(ii,jj,idof) =                           sum( dBdx(0:Ns-1,1) * smt(0:Ns-1,mm) ) * CC
      ;            ; idof = idof + 1         ; coil(icoil)%dT(   jj,idof) = coil(icoil)%dT(jj,idof) + sum( dFdx(0:Ns-1,2) * cmt(0:Ns-1,mm) ) * CC
      ;            ;                         ; coil(icoil)%dB(ii,jj,idof) =                           sum( dBdx(0:Ns-1,2) * cmt(0:Ns-1,mm) ) * CC
      ;            ; idof = idof + 1         ; coil(icoil)%dT(   jj,idof) = coil(icoil)%dT(jj,idof) + sum( dFdx(0:Ns-1,2) * smt(0:Ns-1,mm) ) * CC
      ;            ;                         ; coil(icoil)%dB(ii,jj,idof) =                           sum( dBdx(0:Ns-1,2) * smt(0:Ns-1,mm) ) * CC
      ;            ; idof = idof + 1         ; coil(icoil)%dT(   jj,idof) = coil(icoil)%dT(jj,idof) + sum( dFdx(0:Ns-1,3) * cmt(0:Ns-1,mm) ) * CC
      ;            ;                         ; coil(icoil)%dB(ii,jj,idof) =                           sum( dBdx(0:Ns-1,3) * cmt(0:Ns-1,mm) ) * CC
      ;            ; idof = idof + 1         ; coil(icoil)%dT(   jj,idof) = coil(icoil)%dT(jj,idof) + sum( dFdx(0:Ns-1,3) * smt(0:Ns-1,mm) ) * CC
      ;            ;                         ; coil(icoil)%dB(ii,jj,idof) =                           sum( dBdx(0:Ns-1,3) * smt(0:Ns-1,mm) ) * CC
     enddo
     if( coil(icoil)%Ifree.ne.0 ) then
     ;             ; idof = idof + 1         ; coil(icoil)%dT(   jj,idof) = coil(icoil)%dT(jj,idof) + sum( dFdx(0:Ns-1,0)                )
     ;             ;                         ; coil(icoil)%dB(ii,jj,idof) =                           sum( dBdx(0:Ns-1,0)                )
     endif
     
    enddo ! end do ii;
    
    coil(icoil)%dT(jj,0:Ndof) = coil(icoil)%dT(jj,0:Ndof) * discretecurve * deltatheta

   enddo ! end do jj;
   
   tx(0:Ns-1) = coil(icoil)%xt(0:Ns-1) ; ax(0:Ns-1) = coil(icoil)%xa(0:Ns-1) ! tangent to coil (shorthand);
   ty(0:Ns-1) = coil(icoil)%yt(0:Ns-1) ; ay(0:Ns-1) = coil(icoil)%ya(0:Ns-1) ! tangent to coil (shorthand);
   tz(0:Ns-1) = coil(icoil)%zt(0:Ns-1) ; az(0:Ns-1) = coil(icoil)%za(0:Ns-1) ! tangent to coil (shorthand);
   
   txtx(0:Ns-1) = tx(0:Ns-1)*tx(0:Ns-1) + ty(0:Ns-1)*ty(0:Ns-1) + tz(0:Ns-1)*tz(0:Ns-1)
   txax(0:Ns-1) = tx(0:Ns-1)*ax(0:Ns-1) + ty(0:Ns-1)*ay(0:Ns-1) + tz(0:Ns-1)*az(0:Ns-1)

   length(0:Ns-1) = sqrt(txtx(0:Ns-1)) ; dd(0:Ns-1) = length(0:Ns-1) * txtx(0:Ns-1)
   
#ifdef DEBUG
   do kk = 0, Ns-1 ! 12 Nov 17;
    FATAL( dforce , abs(dd(kk)).lt.sqrtmachprec, fatal )
   enddo
#endif

   ;             ; idof =        0     ; coil(icoil)%dL(idof) = sum(   length(0:Ns-1)                                                                 )
   ;             ; idof = coil(icoil)%gdof
   ;  mm = 0     ; idof = idof + 1     ; coil(icoil)%dL(idof) = sum( ( txax(0:Ns-1) * tx(0:Ns-1) - txtx(0:Ns-1) * ax(0:Ns-1) ) * cmt(0:Ns-1,mm) / dd(0:Ns-1) )
   ;             ; idof = idof + 1     ; coil(icoil)%dL(idof) = sum( ( txax(0:Ns-1) * ty(0:Ns-1) - txtx(0:Ns-1) * ay(0:Ns-1) ) * cmt(0:Ns-1,mm) / dd(0:Ns-1) )
   ;             ; idof = idof + 1     ; coil(icoil)%dL(idof) = sum( ( txax(0:Ns-1) * tz(0:Ns-1) - txtx(0:Ns-1) * az(0:Ns-1) ) * cmt(0:Ns-1,mm) / dd(0:Ns-1) )
   do mm = 1, NF ; idof = idof + 1     ; coil(icoil)%dL(idof) = sum( ( txax(0:Ns-1) * tx(0:Ns-1) - txtx(0:Ns-1) * ax(0:Ns-1) ) * cmt(0:Ns-1,mm) / dd(0:Ns-1) )
    ;            ; idof = idof + 1     ; coil(icoil)%dL(idof) = sum( ( txax(0:Ns-1) * tx(0:Ns-1) - txtx(0:Ns-1) * ax(0:Ns-1) ) * smt(0:Ns-1,mm) / dd(0:Ns-1) )
    ;            ; idof = idof + 1     ; coil(icoil)%dL(idof) = sum( ( txax(0:Ns-1) * ty(0:Ns-1) - txtx(0:Ns-1) * ay(0:Ns-1) ) * cmt(0:Ns-1,mm) / dd(0:Ns-1) )
    ;            ; idof = idof + 1     ; coil(icoil)%dL(idof) = sum( ( txax(0:Ns-1) * ty(0:Ns-1) - txtx(0:Ns-1) * ay(0:Ns-1) ) * smt(0:Ns-1,mm) / dd(0:Ns-1) )
    ;            ; idof = idof + 1     ; coil(icoil)%dL(idof) = sum( ( txax(0:Ns-1) * tz(0:Ns-1) - txtx(0:Ns-1) * az(0:Ns-1) ) * cmt(0:Ns-1,mm) / dd(0:Ns-1) )
    ;            ; idof = idof + 1     ; coil(icoil)%dL(idof) = sum( ( txax(0:Ns-1) * tz(0:Ns-1) - txtx(0:Ns-1) * az(0:Ns-1) ) * smt(0:Ns-1,mm) / dd(0:Ns-1) )
   enddo
   if( coil(icoil)%Ifree.ne.0 ) then
   ;             ; idof = idof + 1     ; coil(icoil)%dL(idof) = zero ! length does not depend on current;
   endif

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
  
  Bdotnsqd(0) = Bdotnsqd(0) * half ! / surf%area
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  ff           = zero
  fdof(1:Ndof) = zero
  
! ff           =                                    exp(totlengt(0)/target_length) * weight_length
! fdof(1:Ndof) = (totlengt(1:Ndof)/target_length) * exp(totlengt(0)/target_length) * weight_length
  
  ff           =  totlengt(0     ) * weight_length ! 12 Nov 17;
  fdof(1:Ndof) =  totlengt(1:Ndof) * weight_length ! 12 Nov 17;
  
  ff           = ff         + weight_tflux * sum( (surf%dT(0:Nz-1,   0)-target_tflux) * (surf%dT(0:Nz-1,0)-target_tflux) ) * half
  do idof = 1, Ndof
  fdof(idof)   = fdof(idof) + weight_tflux * sum( (surf%dT(0:Nz-1,idof)             ) * (surf%dT(0:Nz-1,0)-target_tflux) )
  enddo
  
  ff           = ff           + weight_bnorm * Bdotnsqd(0     )
  fdof(1:Ndof) = fdof(1:Ndof) + weight_bnorm * Bdotnsqd(1:Ndof) ! / surf%area
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

  fdof(1:Ndof) = fdof(1:Ndof) + wspectral * sum(tdof(1:Ndof)*Mdof(1:Ndof)) * tdof(1:Ndof) ! spectral constraints; 12 Nov 17;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

  return
  
end subroutine dforce

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
