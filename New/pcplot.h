!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

!title (pcplot) ! Construct Poincare plot.

!latex \briefly{Construct Poincare plot.}

!latex \calledby{\link{solvers}}
!latex \calls{\link{solvers}}

!latex \section{overview}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

subroutine pcplot
  
  use globals, only : zero, one, ten, half, pi2, pi, ounit, myid, tstart, &
                      minorrad, ellipticity, nrotate, &
                      outplots, punit, &
                      Ntrj, Npts, odetol, iota, xyaxis
  
  implicit none  
  
  include "mpif.h"
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  INTEGER, parameter :: Node = 2, Nxyaxis = 2, Lrwork = Nxyaxis * ( 3 * Nxyaxis + 13 ) / 2
  INTEGER            :: ierr, NN, itrj, ipt, Lwa, ic05nbf, izeta
  REAL               :: Fxy(1:Nxyaxis), tnow, told, rwork(1:Lrwork), tol, xy(1:Node)

  REAL               :: srho, teta, zeta, ax(1:3), at(1:3), az(1:3), xx(1:3), xs(1:3), xt(1:3), xz(1:3), v1(1:3), v2(1:3), w1(1:3), w2(1:3)
  REAL               :: xtt(1:3), xtz(1:3), xzz(1:3)

  external           :: fxyaxis
  
  told = MPI_WTIME()
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
!#ifdef DEBUG
!  
!  srho = one ; teta = zero
!  
!  do izeta = 0, 1
!   zeta = izeta * pi2
!   call knotxx( srho, teta, zeta, ax(1:3), at(1:3), az(1:3), xx(1:3), xs(1:3), xt(1:3), xz(1:3), xtt, xtz, xzz, v1(1:3), v2(1:3), w1(1:3), w2(1:3) )
!   write(ounit,'("pcplot : " 10x " : "3es23.15)') ax(1:3) + minorrad * ( xy(1) * ellipticity * v1(1:3) + xy(2) * v2(1:3) )
!  enddo ! ;
!  
!#endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
!#ifdef DEBUG
!
!  xyaxis(1:Node) = zero ! this was already initialized in initial;
!
!  tnow = MPI_WTIME()
!
!#else

  xyaxis(1:Node) = zero ! this was already initialized in initial;
  
  NN = Nxyaxis ; tol = odetol * ten ; Lwa = Lrwork ; ic05nbf = 1
  
  ;            ; write(ounit,'("pcplot  : ",f10.1," : calling C05NBF ;")') told-tstart

 !FATAL( pcplot , .true., need to replace C05NBF )
  call C05NBF( fxyaxis, NN, xyaxis(1:NN), Fxy(1:NN), tol, rwork(1:Lwa), Lwa, ic05nbf )
  
  tnow = MPI_WTIME()

  select case( ic05nbf )
  case( 0 )    ; write(ounit,'("pcplot  : ",f10.1," : called  C05NBF ; success ; time ="f10.2"s ;")') tnow-tstart, tnow-told
  case default ; write(ounit,'("pcplot  : ",f10.1," : called  C05NBF ; error   ; time ="f10.2"s ;")') tnow-tstart, tnow-told ; xyaxis(1:Node) = zero
  end select
  
!#endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  if( myid.eq.0 ) write(ounit,'("pcplot  : ",f10.1," : Ntrj =",i4," ; Npts =",i6," ; odetol =",es12.5," ;")') tnow-tstart, Ntrj, Npts, odetol
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

  open( punit, file=trim(outplots), status="unknown", form="formatted" )
  
  write( punit,'(2i9)' ) Ntrj, Npts
  
  told = tnow
  
  do itrj = 0, Ntrj
   
   iota(itrj,1) = itrj * one / Ntrj
   
   xy(1:Node) = xyaxis(1:Node) + itrj * ( (/ one, zero /) - xyaxis(1:2) ) / Ntrj ! set starting point;
   if( myid.eq.0 ) write(punit,'(2es23.15)') xy(1:2) * minorrad * (/ ellipticity             , one        /)
!  xy(1:Node) = (/ iota(itrj,1), zero /)
!  if( myid.eq.0 ) write(punit,'(2es23.15)') xy(1) * minorrad * (/ ellipticity * cos(xy(2)), sin(xy(2)) /)
   
   do ipt = 1, Npts
    
    call pmap( Node, xy(1:Node) )
    
    if( myid.eq.0 ) write(punit,'(2es23.15)') xy(1:2) * minorrad * (/ ellipticity             , one        /)
!   if( myid.eq.0 ) write(punit,'(2es23.15)') xy(1) * minorrad * (/ ellipticity * cos(xy(2)), sin(xy(2)) /)
    
   enddo ! end of do ipt;
   
   iota(itrj,2) = zero
!  iota(itrj,2) = xy(2) / (Npts*pi2)
   
   tnow = MPI_WTIME()
   
   if( myid.eq.0 ) write(ounit,'("pcplot  : ",f10.1," : ",i4," : iota =",es23.15," ; time =",f10.2," ;")') tnow-tstart, itrj, iota(itrj,2), tnow-told
   
   told = tnow
   
  enddo ! end of do itrj;
  
  close(punit)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  return
  
end subroutine pcplot

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
!
!subroutine bstfield( zeta, st, Bst )
!  
!  use globals, only : zero, one, pi2, sqrtmachprec, ounit, myid, ncpu, Ncoils, Ns, coil, minorrad, ellipticity
!  
!  implicit none  
!  
!  include "mpif.h"
!  
!  REAL    :: zeta, st(*), Bst(*)
!  
!  INTEGER :: icoil, kk, ierr
!  REAL    :: ax(1:3), at(1:3), az(1:3), xx(1:3), xs(1:3), xt(1:3), xz(1:3), v1(1:3), v2(1:3), w1(1:3), w2(1:3)
!  REAL    :: a, b, c, d, e, f, g, h, i, idet, inversematrix(1:3,1:3)
!  REAL    :: dx, dy, dz, rr(1:3), dd(1:3), tx, ty, tz, Bt(1:3), Bi(1:3), Bs(1:3)
!  REAL    :: ssmax
!
!  ssmax = 1.1
!
!  if( st(1).lt.sqrtmachprec .or. st(1).gt.ssmax ) then ; Bst(1:2) = zero ; return
!  endif
!  
!  call knotxx( st(1), st(2), zeta, ax(1:3), at(1:3), az(1:3), xx(1:3), xs(1:3), xt(1:3), xz(1:3), v1(1:3), v2(1:3), w1(1:3), w2(1:3) )
!  
!  a = xs(1) ; b = xt(1) ; c = xz(1)
!  d = xs(2) ; e = xt(2) ; f = xz(2)
!  g = xs(3) ; h = xt(3) ; i = xz(3)
!
! !idet = one / ( - c * e * g + b * f * g + c * d * h - a * f * h - b * d * i + a * e * i ) ! this factor will cancel; not required;
!  
!  inversematrix(1,1) = ( -f * h + e * i ) ; inversematrix(1,2) = (  c * h - b * i ) ; inversematrix(1,3) = ( -c * e + b * f )
!  inversematrix(2,1) = (  f * g - d * i ) ; inversematrix(2,2) = ( -c * g + a * i ) ; inversematrix(2,3) = (  c * d - a * f )  
!  inversematrix(3,1) = ( -e * g + d * h ) ; inversematrix(3,2) = (  b * g - a * h ) ; inversematrix(3,3) = ( -b * d + a * e )
!  
!  Bt(1:3) = zero
!
!  do icoil = 1, Ncoils
!   
!  !if( myid.ne.modulo(icoil-1,ncpu) ) cycle ! parallelization loop;
!   
!   Bi(1:3) = zero
!   
!   do kk = 1, Ns
!    
!    dx = xx(1) - coil(icoil)%xx(kk) ! distance from evaluation point to curve;
!    dy = xx(2) - coil(icoil)%yy(kk)
!    dz = xx(3) - coil(icoil)%zz(kk)
!
!    rr(2) = dx**2 + dy**2 + dz**2 ; rr(1) = sqrt(rr(2)) ; rr(3) = rr(1) * rr(2)
!
!    FATAL( pcplot , abs(rr(3)).lt.sqrtmachprec, divide by zero )
!
!    ;                             ;                     ; dd(3) = one / rr(3)
!    
!    tx = coil(icoil)%xt(kk) ! tangent to coil (shorthand);
!    ty = coil(icoil)%yt(kk)
!    tz = coil(icoil)%zt(kk)
!    
!    Bi(1:3) = Bi(1:3) + (/ ty * dz - tz * dy, tz * dx - tx * dz, tx * dy - ty * dx /) * dd(3)
!    
!   enddo ! end of do kk;
!   
!  !Bi(1:3,icoil) = Bi(1:3,icoil) * coil(icoil)%I
!
!   Bt(1:3) = Bt(1:3) + Bi(1:3) * coil(icoil)%I
!   
!  enddo ! end of do icoil;
!  
! !call MPI_REDUCE( Bi(1:3), Bt(1:3), 3, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
!
! !RlBCAST( Bt(1:3), 3, 0 )
!
!  Bs(1:3) = matmul( inversematrix(1:3,1:3), Bt(1:3) )
!  
!  FATAL( pcplot , abs(Bs(3)).lt.sqrtmachprec, divide by zero )
!
!  Bst(1:2) = Bs(1:2) / Bs(3)
!
!  return
!  
!end subroutine bstfield
!
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

subroutine bxyfield( zeta, xy, Bxy )
  
  use globals, only : zero, one, pi2, sqrtmachprec, ounit, myid, ncpu, Ncoils, Ns, coil, minorrad, ellipticity
  
  implicit none  
  
  include "mpif.h"
  
  REAL    :: zeta, xy(*), Bxy(*)
  
  INTEGER :: icoil, kk, ierr, isurf
  REAL    :: st(1:2)
  REAL    :: ax(1:3), at(1:3), az(1:3), xx(1:3), xs(1:3), xt(1:3), xz(1:3), v1(1:3), v2(1:3), w1(1:3), w2(1:3)
  REAL    :: xtt(1:3), xtz(1:3), xzz(1:3)
  REAL    :: a, b, c, d, e, f, g, h, i, idet, inversematrix(1:3,1:3)
  REAL    :: dx, dy, dz, rr(1:3), dd(1:3), tx, ty, tz, Bt(1:3), Bi(1:3), Bs(1:3)
  REAL    :: ssmax, srho, teta

  st(1) = sqrt( xy(1)*xy(1) + xy(2)*xy(2) ) ; st(2) = atan2( xy(2), xy(1) )
  
#ifdef DEBUG

#else

  ssmax = 1.1
  if(                            st(1).gt.ssmax ) then ; Bxy(1:2) = zero ; return
  endif

#endif

  isurf = 1 ; srho = one ; teta = zero

  call knotxx( isurf, st(1), st(2), zeta, ax(1:3), at(1:3), az(1:3), xx(1:3), xs(1:3), xt(1:3), xz(1:3), xtt, xtz, xzz, v1(1:3), v2(1:3), w1(1:3), w2(1:3) )
  
  xx(1:3) = ax(1:3) + minorrad * ( xy(1) * ellipticity * v1(1:3) + xy(2) * v2(1:3) )
  
  a = minorrad * ellipticity * v1(1) ; b = minorrad * v2(1) ; c = az(1) + minorrad * ( xy(1) * ellipticity * w1(1) + xy(2) * w2(1) )
  d = minorrad * ellipticity * v1(2) ; e = minorrad * v2(2) ; f = az(2) + minorrad * ( xy(1) * ellipticity * w1(2) + xy(2) * w2(2) )
  g = minorrad * ellipticity * v1(3) ; h = minorrad * v2(3) ; i = az(3) + minorrad * ( xy(1) * ellipticity * w1(3) + xy(2) * w2(3) )
  
! idet = one / ( - c * e * g + b * f * g + c * d * h - a * f * h - b * d * i + a * e * i ) ! this factor will cancel; not required;
  
  inversematrix(1,1) = ( -f * h + e * i ) ; inversematrix(1,2) = (  c * h - b * i ) ; inversematrix(1,3) = ( -c * e + b * f )
  inversematrix(2,1) = (  f * g - d * i ) ; inversematrix(2,2) = ( -c * g + a * i ) ; inversematrix(2,3) = (  c * d - a * f )  
  inversematrix(3,1) = ( -e * g + d * h ) ; inversematrix(3,2) = (  b * g - a * h ) ; inversematrix(3,3) = ( -b * d + a * e )
  
  Bt(1:3) = zero

  do icoil = 1, Ncoils
   
  !if( myid.ne.modulo(icoil-1,ncpu) ) cycle ! parallelization loop;
   
   Bi(1:3) = zero
   
   do kk = 0, Ns-1 ! 12 Nov 17;
    
    dx = xx(1) - coil(icoil)%xx(1,kk) ! distance from evaluation point to curve;
    dy = xx(2) - coil(icoil)%xx(2,kk)
    dz = xx(3) - coil(icoil)%xx(3,kk)

    rr(2) = dx**2 + dy**2 + dz**2 ; rr(1) = sqrt(rr(2)) ; rr(3) = rr(1) * rr(2)

    FATAL( pcplot , abs(rr(3)).lt.sqrtmachprec, divide by zero )

    ;                             ;                     ; dd(3) = one / rr(3)
    
    tx = coil(icoil)%xt(1,kk) ! tangent to coil (shorthand);
    ty = coil(icoil)%xt(2,kk)
    tz = coil(icoil)%xt(3,kk)
    
    Bi(1:3) = Bi(1:3) + (/ ty * dz - tz * dy, tz * dx - tx * dz, tx * dy - ty * dx /) * dd(3)
    
   enddo ! end of do kk;
   
  !Bi(1:3,icoil) = Bi(1:3,icoil) * coil(icoil)%I

   Bt(1:3) = Bt(1:3) + Bi(1:3) * coil(icoil)%I
   
  enddo ! end of do icoil;
  
 !call MPI_REDUCE( Bi(1:3), Bt(1:3), 3, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )

 !RlBCAST( Bt(1:3), 3, 0 )

  Bs(1:3) = matmul( inversematrix(1:3,1:3), Bt(1:3) )
  
  FATAL( pcplot , abs(Bs(3)).lt.sqrtmachprec, divide by zero )

  Bxy(1:2) = Bs(1:2) / Bs(3)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  return
  
end subroutine bxyfield

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

subroutine pmap( Node, xy )
  
  use globals, only : zero, pi2, odetol, nrotate
  
  implicit none  
  
  include "mpif.h"
  
  INTEGER   :: Node, id02bjf
  REAL      :: xy(1:Node), zeta, zetaend, tol, rwork(1:20*Node)
  CHARACTER :: relabs
  
  external  :: bstfield, bxyfield, D02BJX, D02BJW
  
  zeta = zero ; zetaend = pi2 ; relabs = 'D' ; tol = odetol ; id02bjf = 1
  
  call D02BJF( zeta, zetaend, Node, xy(1:Node), bxyfield, odetol, relabs, D02BJX, D02BJW, rwork(1:20*Node), id02bjf )
  xy(1:2) = xy(1:2) * (-1)**nrotate

! call D02BJF( zeta, zetaend, Node, xy(1:Node), bstfield, odetol, relabs, D02BJX, D02BJW, rwork(1:20*Node), id02bjf )
! xy(2) = xy(2) + nrotate * pi
  
  return
  
end subroutine pmap

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

subroutine fxyaxis( Nxyaxis, xyi, Fxy, iflag )
  
  use globals, only : zero, pi2, ounit, odetol
  
  implicit none  
  
  include "mpif.h"
  
  INTEGER, parameter :: Node = 2
  INTEGER            :: Nxyaxis, iflag
  REAL               :: xyi(1:Nxyaxis), Fxy(1:Nxyaxis), xye(1:Nxyaxis)
  
  xye(1:Nxyaxis) = xyi(1:Nxyaxis)

  call pmap( Node, xye(1:Node) )

  Fxy(1:Nxyaxis) = xye(1:Node) - xyi(1:Node) ! fixed point = magnetic axis;

  write(ounit,'("fxyaxis : " 10x " : ("f13.9","f13.9" ) --> ("f13.9","f13.9" ) : F = ("es10.2","es10.2" ) ;")') xyi(1:Nxyaxis), xye(1:Nxyaxis), Fxy(1:Nxyaxis)

  return

end subroutine fxyaxis

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
