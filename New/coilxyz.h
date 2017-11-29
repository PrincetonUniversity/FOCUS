!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

!title (coilxyz) ! map coils to real space.

!latex \briefly{Initialize the coils data with Fourier series.}

!latex \calledby{\link{focus}}
!latex \calls{\link{}}

!latex \subsection{overview}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

subroutine coilxyz( icoil )
  
  use globals, only : zero, one, three, myid, coil, cmt, smt, Ns, surf, Nt, Nz
  
  implicit none
  
  include "mpif.h"
  
  INTEGER, intent(in) :: icoil
  
  INTEGER :: mm, ierr, idof, ii, jj, kk, isurf
  REAL    :: rx, ry, rz, dist, lx, ly, lz, rdotn
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  CHECK( coilxyz, .not.allocated(cmt), trigonometric factors have not been allocated )
  CHECK( coilxyz, .not.allocated(smt), trigonometric factors have not been allocated ) 

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  coil(icoil)%xx(1:3,0:Ns-1) = zero
 !coil(icoil)%yy = zero
 !coil(icoil)%zz = zero
  coil(icoil)%xt(1:3,0:Ns-1) = zero
 !coil(icoil)%yt = zero
 !coil(icoil)%zt = zero
  coil(icoil)%xa = zero
  coil(icoil)%ya = zero
  coil(icoil)%za = zero
  
  ;  mm = 0
  
  coil(icoil)%xx(1,0:Ns-1) = cmt(0:Ns-1,mm) * coil(icoil)%xc(mm)
  coil(icoil)%xx(2,0:Ns-1) = cmt(0:Ns-1,mm) * coil(icoil)%yc(mm)
  coil(icoil)%xx(3,0:Ns-1) = cmt(0:Ns-1,mm) * coil(icoil)%zc(mm)    
  
  do mm = 1, coil(icoil)%NF   
   
   coil(icoil)%xx(1,0:Ns-1) = coil(icoil)%xx(1,0:Ns-1) + (   cmt(0:Ns-1,mm) * coil(icoil)%xc(mm) + smt(0:Ns-1,mm) * coil(icoil)%xs(mm) )
   coil(icoil)%xx(2,0:Ns-1) = coil(icoil)%xx(2,0:Ns-1) + (   cmt(0:Ns-1,mm) * coil(icoil)%yc(mm) + smt(0:Ns-1,mm) * coil(icoil)%ys(mm) )
   coil(icoil)%xx(3,0:Ns-1) = coil(icoil)%xx(3,0:Ns-1) + (   cmt(0:Ns-1,mm) * coil(icoil)%zc(mm) + smt(0:Ns-1,mm) * coil(icoil)%zs(mm) )
   
   coil(icoil)%xt(1,0:Ns-1) = coil(icoil)%xt(1,0:Ns-1) + ( - smt(0:Ns-1,mm) * coil(icoil)%xc(mm) + cmt(0:Ns-1,mm) * coil(icoil)%xs(mm) ) * mm
   coil(icoil)%xt(2,0:Ns-1) = coil(icoil)%xt(2,0:Ns-1) + ( - smt(0:Ns-1,mm) * coil(icoil)%yc(mm) + cmt(0:Ns-1,mm) * coil(icoil)%ys(mm) ) * mm
   coil(icoil)%xt(3,0:Ns-1) = coil(icoil)%xt(3,0:Ns-1) + ( - smt(0:Ns-1,mm) * coil(icoil)%zc(mm) + cmt(0:Ns-1,mm) * coil(icoil)%zs(mm) ) * mm
   
   coil(icoil)%xa(0:Ns-1) = coil(icoil)%xa(0:Ns-1) + ( - cmt(0:Ns-1,mm) * coil(icoil)%xc(mm) - smt(0:Ns-1,mm) * coil(icoil)%xs(mm) ) * mm*mm
   coil(icoil)%ya(0:Ns-1) = coil(icoil)%ya(0:Ns-1) + ( - cmt(0:Ns-1,mm) * coil(icoil)%yc(mm) - smt(0:Ns-1,mm) * coil(icoil)%ys(mm) ) * mm*mm
   coil(icoil)%za(0:Ns-1) = coil(icoil)%za(0:Ns-1) + ( - cmt(0:Ns-1,mm) * coil(icoil)%zc(mm) - smt(0:Ns-1,mm) * coil(icoil)%zs(mm) ) * mm*mm
   
  enddo ! end of do mm; 
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  isurf = 1

  do kk = 0, Ns-1
   do ii = 0, Nt-1
    do jj = 0, Nz-1
     
     rx = surf(isurf)%xx(1,ii,jj) - coil(icoil)%xx(1,kk) ; lx = surf(isurf)%nn(1,ii,jj)
     ry = surf(isurf)%xx(2,ii,jj) - coil(icoil)%xx(2,kk) ; ly = surf(isurf)%nn(2,ii,jj)
     rz = surf(isurf)%xx(3,ii,jj) - coil(icoil)%xx(3,kk) ; lz = surf(isurf)%nn(3,ii,jj)
     
     dist = sqrt( rx*rx + ry*ry + rz*rz ) ; rdotn = rx*lx + ry*ly + rz*lz
     
     coil(icoil)%RR(1:9,ii,jj,kk) = three * (/ rx * rx, rx * ry, rx * rz, ry * rx, ry * ry, ry * rz, rz * rx, rz * ry, rz * rz /) / dist**5 &
                                  +         (/   one  ,  zero  ,  zero  ,  zero  ,   one  ,  zero  ,  zero  ,  zero  ,   one   /) / dist**3

     coil(icoil)%Rn(1:3,ii,jj,kk) = three * (/ rx, ry, rz /) * rdotn / dist**5 - (/ lx, ly, lz /) / dist**3 ! 28 Nov 17;
     
    enddo ! end of do jj; 16 Nov 17;
   enddo ! end of do ii; 16 Nov 17;
  enddo ! end of do kk; 16 Nov 17;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

  return

end subroutine coilxyz

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

!SUBROUTINE discfou2
!
!  use globals, only: zero, pi2, myid, ncpu, ounit, coil, FouCoil, Ncoils
!  implicit none
!  include "mpif.h"
!
!  INTEGER :: icoil, iorder, llmodnp, ierr, NS, NF
!  !-------------------------call fouriermatr----------------------------------------------------  
!  do icoil = 1, Ncoils
!
!     NS = coil(icoil)%NS; NF = FouCoil(icoil)%NF  ! allias variable for simplicity;
!     !reset to zero for all the coils;
!     coil(icoil)%xx = zero
!     coil(icoil)%yy = zero
!     coil(icoil)%zz = zero
!     coil(icoil)%xt = zero
!     coil(icoil)%yt = zero
!     coil(icoil)%zt = zero
!     coil(icoil)%xa = zero
!     coil(icoil)%ya = zero
!     coil(icoil)%za = zero
!
!     coil(icoil)%dd = pi2 / NS
!     
!     if( myid.ne.modulo(icoil-1,ncpu) ) cycle ! parallelization loop;
!
!    iorder = 0
!    call fouriermatrix( Foucoil(icoil)%xc, Foucoil(icoil)%xs, coil(icoil)%xx, NF, NS, iorder)
!    call fouriermatrix( Foucoil(icoil)%yc, Foucoil(icoil)%ys, coil(icoil)%yy, NF, NS, iorder)
!    call fouriermatrix( Foucoil(icoil)%zc, Foucoil(icoil)%zs, coil(icoil)%zz, NF, NS, iorder)
!
!    iorder = 1  
!    call fouriermatrix( Foucoil(icoil)%xc, Foucoil(icoil)%xs, coil(icoil)%xt, NF, NS, iorder)
!    call fouriermatrix( Foucoil(icoil)%yc, Foucoil(icoil)%ys, coil(icoil)%yt, NF, NS, iorder)
!    call fouriermatrix( Foucoil(icoil)%zc, Foucoil(icoil)%zs, coil(icoil)%zt, NF, NS, iorder)
!
!    iorder = 2 
!    call fouriermatrix( Foucoil(icoil)%xc, Foucoil(icoil)%xs, coil(icoil)%xa, NF, NS, iorder)
!    call fouriermatrix( Foucoil(icoil)%yc, Foucoil(icoil)%ys, coil(icoil)%ya, NF, NS, iorder)
!    call fouriermatrix( Foucoil(icoil)%zc, Foucoil(icoil)%zs, coil(icoil)%za, NF, NS, iorder)
!   
!  enddo
!  !-------------------------broadcast coil data-------------------------------------------------  
!  do icoil = 1, Ncoils ; llmodnp = modulo(icoil-1,ncpu)
!     RlBCAST( coil(icoil)%xx(0:coil(icoil)%NS), coil(icoil)%NS+1, llmodnp )
!     RlBCAST( coil(icoil)%yy(0:coil(icoil)%NS), coil(icoil)%NS+1, llmodnp )
!     RlBCAST( coil(icoil)%zz(0:coil(icoil)%NS), coil(icoil)%NS+1, llmodnp )
!     RlBCAST( coil(icoil)%xt(0:coil(icoil)%NS), coil(icoil)%NS+1, llmodnp )
!     RlBCAST( coil(icoil)%yt(0:coil(icoil)%NS), coil(icoil)%NS+1, llmodnp )
!     RlBCAST( coil(icoil)%zt(0:coil(icoil)%NS), coil(icoil)%NS+1, llmodnp )
!     RlBCAST( coil(icoil)%xa(0:coil(icoil)%NS), coil(icoil)%NS+1, llmodnp )
!     RlBCAST( coil(icoil)%ya(0:coil(icoil)%NS), coil(icoil)%NS+1, llmodnp )
!     RlBCAST( coil(icoil)%za(0:coil(icoil)%NS), coil(icoil)%NS+1, llmodnp )
!  enddo
!  
!  return
!END SUBROUTINE discfou2

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!subroutine fouriermatrix( xc, xs, xx, NF, ND, order )
!  !---------------------------------------------------------------------------------------------
!  ! This subroutine uses matrix operations to discretize data from Fourier harmonics.
!  ! It's supposed to be the fastest method.
!  ! DATE: 2017/03/18
!  !---------------------------------------------------------------------------------------------
!  use globals, only : zero, pi2
!  implicit none
!  
!  INTEGER, intent(in ) :: NF, ND, order
!  REAL   , intent(in ) :: xc(0:NF), xs(0:NF)
!  REAL   , intent(out) :: xx(0:ND)
!
!  INTEGER              :: i
!  REAL                 :: nn(0:NF, 1:1), tt(1:1, 0:ND), nt(0:NF, 0:ND), &
!                       &  tc(1:1, 0:NF), ts(1:1, 0:NF), tx(1:1 , 0:ND), &
!                       &  cnt(0:NF, 0:ND), snt(0:NF, 0:ND)
!
!  !----------------------------data copy to matrix----------------------------------------------
!  if ( size(xc) /= NF+1 ) STOP "Wrong input size for xc in subroutine fouriermatrix!"
!  if ( size(xs) /= NF+1 ) STOP "Wrong input size for xs in subroutine fouriermatrix!"
!  if ( size(xx) /= ND+1 ) STOP "Wrong input size for xx in subroutine fouriermatrix!"
!
!  tc(1, 0:NF) = xc(0:NF) ! cos harmonics;
!  ts(1, 0:NF) = xs(0:NF) ! sin harmonics;
!  tx(1, 0:ND) = xx(0:ND) ! data coordinates;
!  !----------------------------matrix assignmengt-----------------------------------------------
!  nn(0:NF, 1) = (/ (i, i=0,NF) /) ! n;
!  tt(1, 0:ND) = (/ (i*pi2/ND, i=0, ND) /) ! angle, t;
!  nt = matmul(nn, tt)
!  cnt = cos(nt)
!  snt = sin(nt)
!  !----------------------------select oder------------------------------------------------------
!
!
!
!  select case (order)
!  case (0)  ! 0-order
!  case (1)  ! 1st-order
!   do i = 0, ND
!    cnt(0:NF, i) = - nn(0:NF, 1)     * snt(0:NF, i)
!    snt(0:NF, i) =   nn(0:NF, 1)     * cnt(0:NF, i)
!   enddo
!  case (2)  ! 2nd-order
!   do i = 0, ND
!    cnt(0:NF, i) = -(nn(0:NF, 1)**2) * cnt(0:NF, i)
!    snt(0:NF, i) = -(nn(0:NF, 1)**2) * snt(0:NF, i)
!   enddo
!  case default
!   STOP "Invalid order in subroutine fouriermatrix"
!  end select
!
!  !----------------------------final multiplication---------------------------------------------
!  tx = matmul(tc, cnt) + matmul(ts, snt)
!  xx(0:ND) = tx(1, 0:ND)
!
!  return
!END subroutine fouriermatrix

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

