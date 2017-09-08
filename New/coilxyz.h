
!title (coilxyz) ! map coils to real space.

!latex \briefly{Initialize the coils data with Fourier series.}

!latex \calledby{\link{focus}}
!latex \calls{\link{}}

!latex \subsection{overview}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine coilxyz( icoil )
  
  use globals, only : zero, Nseg, coil, cmt, smt
  
  implicit none
  
  include "mpif.h"
  
  INTEGER, intent(in) :: icoil

! LOGICAL :: exist
  INTEGER :: mm !maxnseg, ifirst, NF, itmp, ip, icoef, jj
! REAL    :: zeta, totalcurrent, Ro, Zo, r1, r2, z1, z2, tt, ax(1:3), at(1:3), az(1:3), xx(1:3), xt(1:3), xz(1:3)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! INTEGER          :: iseg, mm, NS, ierr, astat, ip
! REAL             :: tt
! REAL,allocatable :: cmt(:,:), smt(:,:)

! FATAL( coilxyz, .true., under reconstruction )

!  do icoil = 1, Ncoils*Npc
  
!if( (coil(icoil)%Lc + ifirst) /= 0) then  !first time or if Lc/=0, then need discretize;
  
!reset to zero for all the coils;
  
  coil(icoil)%xx = zero
  coil(icoil)%yy = zero
  coil(icoil)%zz = zero
  coil(icoil)%xt = zero
  coil(icoil)%yt = zero
  coil(icoil)%zt = zero
  coil(icoil)%xa = zero
  coil(icoil)%ya = zero
  coil(icoil)%za = zero
!
!        !if( myid.ne.modulo(icoil-1,ncpu) ) cycle ! parallelization loop;
!
!        select case (coil(icoil)%itype)
!        case( 1 )
!
!           NS = coil(icoil)%NS; NF = FouCoil(icoil)%NF  ! allias variable for simplicity;
!           SALLOCATE( cmt, (0:NS, 0:NF), zero )
!           SALLOCATE( smt, (0:NS, 0:NF), zero )
!
!           do iseg = 0, NS ; tt = iseg * pi2 / NS
!              do mm = 0, NF
!                 cmt(iseg,mm) = cos( mm * tt )
!                 smt(iseg,mm) = sin( mm * tt )
!              enddo
!           enddo
!
!-------------------------calculate coil data-------------------------------------------------  
  ;  mm = 0

  coil(icoil)%xx(0:Nseg) = cmt(0:Nseg,mm) * coil(icoil)%xc(mm)
  coil(icoil)%yy(0:Nseg) = cmt(0:Nseg,mm) * coil(icoil)%yc(mm)
  coil(icoil)%zz(0:Nseg) = cmt(0:Nseg,mm) * coil(icoil)%zc(mm)    

  do mm = 1, coil(icoil)%NF   

   coil(icoil)%xx(0:Nseg) = coil(icoil)%xx(0:Nseg) + (   cmt(0:Nseg,mm) * coil(icoil)%xc(mm) + smt(0:Nseg,mm) * coil(icoil)%xs(mm) )
   coil(icoil)%yy(0:Nseg) = coil(icoil)%yy(0:Nseg) + (   cmt(0:Nseg,mm) * coil(icoil)%yc(mm) + smt(0:Nseg,mm) * coil(icoil)%ys(mm) )
   coil(icoil)%zz(0:Nseg) = coil(icoil)%zz(0:Nseg) + (   cmt(0:Nseg,mm) * coil(icoil)%zc(mm) + smt(0:Nseg,mm) * coil(icoil)%zs(mm) )
   
   coil(icoil)%xt(0:Nseg) = coil(icoil)%xt(0:Nseg) + ( - smt(0:Nseg,mm) * coil(icoil)%xc(mm) + cmt(0:Nseg,mm) * coil(icoil)%xs(mm) ) * mm
   coil(icoil)%yt(0:Nseg) = coil(icoil)%yt(0:Nseg) + ( - smt(0:Nseg,mm) * coil(icoil)%yc(mm) + cmt(0:Nseg,mm) * coil(icoil)%ys(mm) ) * mm
   coil(icoil)%zt(0:Nseg) = coil(icoil)%zt(0:Nseg) + ( - smt(0:Nseg,mm) * coil(icoil)%zc(mm) + cmt(0:Nseg,mm) * coil(icoil)%zs(mm) ) * mm
   
   coil(icoil)%xa(0:Nseg) = coil(icoil)%xa(0:Nseg) + ( - cmt(0:Nseg,mm) * coil(icoil)%xc(mm) - smt(0:Nseg,mm) * coil(icoil)%xs(mm) ) * mm*mm
   coil(icoil)%ya(0:Nseg) = coil(icoil)%ya(0:Nseg) + ( - cmt(0:Nseg,mm) * coil(icoil)%yc(mm) - smt(0:Nseg,mm) * coil(icoil)%ys(mm) ) * mm*mm
   coil(icoil)%za(0:Nseg) = coil(icoil)%za(0:Nseg) + ( - cmt(0:Nseg,mm) * coil(icoil)%zc(mm) - smt(0:Nseg,mm) * coil(icoil)%zs(mm) ) * mm*mm

  enddo ! end of do mm; 

!
!           if(ifirst /= 0) then
!              ip = (icoil-1)/Ncoils  ! the integer is the period number;
!              DoF(icoil)%xof(1:Nseg,      1:  NF+1) =  cosip(ip) * cmt(1:Nseg, 0:NF)  !x/xc
!              DoF(icoil)%xof(1:Nseg,   NF+2:2*NF+1) =  cosip(ip) * smt(1:Nseg, 1:NF)  !x/xs
!              DoF(icoil)%xof(1:NS, 2*NF+2:3*NF+2) = -sinip(ip) * cmt(1:NS, 0:NF)  !x/yc ; valid for ip>0 ;
!              DoF(icoil)%xof(1:NS, 3*NF+3:4*NF+2) = -sinip(ip) * smt(1:NS, 1:NF)  !x/ys ; valid for ip>0 ;
!              DoF(icoil)%yof(1:NS,      1:  NF+1) =  sinip(ip) * cmt(1:NS, 0:NF)  !y/xc ; valid for ip>0 ;
!              DoF(icoil)%yof(1:NS,   NF+2:2*NF+1) =  sinip(ip) * smt(1:NS, 1:NF)  !y/xs ; valid for ip>0 ;
!              DoF(icoil)%yof(1:NS, 2*NF+2:3*NF+2) =  cosip(ip) * cmt(1:NS, 0:NF)  !y/yc
!              DoF(icoil)%yof(1:NS, 3*NF+3:4*NF+2) =  cosip(ip) * smt(1:NS, 1:NF)  !y/ys
!              DoF(icoil)%zof(1:NS, 4*NF+3:5*NF+3) =              cmt(1:NS, 0:NF)  !z/zc
!              DoF(icoil)%zof(1:NS, 5*NF+4:6*NF+3) =              smt(1:NS, 1:NF)  !z/zs
!           endif
!
!           coil(icoil)%dd = pi2 / NS  ! discretizing factor;
!
!           DALLOCATE(cmt)
!           DALLOCATE(smt)
!
!        case default
!           FATAL(discoil, .true., not supported coil types)
!        end select
!
!     endif
!
!  enddo ! end of do icoil

  return
end subroutine coilxyz

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

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

