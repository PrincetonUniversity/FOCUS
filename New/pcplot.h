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
                      Ntrj, Etrj, Npts, odetol, iota, xyaxis, &
                      ga00aainput, ga01aaX, ga01aaO, tr00aainput, &
                      surf, pp, qq

 !use oculus, only : oculustr00aa

  implicit none  
  
  include "mpif.h"
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

  INTEGER, parameter :: Node = 2, Nxyaxis = 2, Lrwork = Nxyaxis * ( 3 * Nxyaxis + 13 ) / 2
  INTEGER            :: ierr, NN, itrj, ipt, Lwa, ic05nbf, izeta
  REAL               :: Fxy(1:Nxyaxis), tnow, told, rwork(1:Lrwork), tol, xy(1:Node)

  INTEGER            :: isurf, ifail , iguess, astat
  REAL               :: srho, teta, zeta, ax(1:3), at(1:3), az(1:3), xx(1:3), xs(1:3), xt(1:3), xz(1:3), xtt(1:3), xtz(1:3), xzz(1:3)
  REAL               :: v1(1:3), v2(1:3), w1(1:3), w2(1:3)
  
  external           :: fxyaxis
  
  told = MPI_WTIME()
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  isurf = 1 ; srho = one ; teta = zero ; zeta = zero
  
  call knotxx( isurf, srho , teta , zeta, ax, at, az, xx, xs, xt, xz, xtt, xtz, xzz, v1, v2, w1, w2 )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  ga00aainput%R = ax(1) ; ga00aainput%Z = ax(3) ; ifail = 1 ; call ga00aa( ga00aainput, ifail )
 
  write(ounit,1000) ga00aainput%R, ga00aainput%Z, ga00aainput%iota, ga00aainput%error, ifail
  
1000 format("pcplot  : " 10x " : oculus ga00aa axis  : (",f8.4,",",f8.4," ) ; iota=",f9.6," ; err=",es10.3," ; ifail=",i3," ;")
  
  itrj = 0 
  iota(itrj,1) = zero
  iota(itrj,2) = ga00aainput%iota

  ifail =  9 ; call ga00aa( ga00aainput, ifail )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  if( qq.gt.0 ) then
  
   ga01aaX%p = pp ; ga01aaX%q = qq ; ga01aaX%Ntor = 8 ; ga01aaX%R = 1.328 ; ga01aaX%Z = 0.000 ; ifail =  1 ; call ga01aa( ga01aaX, ifail )
   ga01aaO%p = pp ; ga01aaO%q = qq ; ga01aaO%Ntor = 8 ; ga01aaO%R = 1.293 ; ga01aaO%Z = 0.122 ; ifail =  1 ; call ga01aa( ga01aaO, ifail )
   
   write(ounit,1001) ga01aaX%p, ga01aaX%q, ga01aaX%R, ga01aaX%Z, ga01aaX%error, ifail
   write(ounit,1001) ga01aaO%p, ga01aaO%q, ga01aaO%R, ga01aaO%Z, ga01aaO%error, ifail
   
1001 format("pcplot  : ", 10x ," : oculus ga01aa (",i3,",",i3," ) : (",f8.4,",",f8.4," ) ; err=",es10.3," ; ifail=",i3," ;")
   
!  ifail =  9 ; call ga01aa( ga01aaX, ifail )
!  ifail =  9 ; call ga01aa( ga01aaO, ifail )
   
  endif ! end of if( qq.gt.0 ) ; 10 Dec 17;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  open( punit, file=trim(outplots), status="unknown", form="formatted" )
  
  write( punit,'(2i9)') Ntrj + Etrj - 1, Npts
  
  tr00aainput%Ppts = Npts
  
  do itrj = 1, Ntrj + Etrj
   
   tr00aainput%Ra     = ga00aainput%R
   tr00aainput%Za     = ga00aainput%Z
   tr00aainput%R      = ga00aainput%R + itrj * ( xx(1) - ga00aainput%R ) / ( Ntrj + Etrj )
   tr00aainput%Z      = ga00aainput%Z + itrj * ( xx(3) - ga00aainput%Z ) / ( Ntrj + Etrj )

   if( qq.gt.0 ) then
   tr00aainput%R      = ga01aaX%R + (itrj-1) * ( ga01aaO%R - ga01aaX%R ) / ( Ntrj + Etrj - 1 )
   tr00aainput%Z      = ga01aaX%Z + (itrj-1) * ( ga01aaO%Z - ga01aaX%Z ) / ( Ntrj + Etrj - 1 )
   endif

   iota(itrj,1) = sqrt((tr00aainput%R-tr00aainput%Ra)**2+(tr00aainput%Z-tr00aainput%Za)**2)
   
   ifail = 1 ; call tr00aa( tr00aainput, ifail )
   
   write(ounit,1002) itrj, tr00aainput%iota
   
1002 format("pcplot  : " 10x " : oculus tr00aa ",i4," : iota =",f9.6," ;")
   
   iota(itrj,2) = tr00aainput%iota

   write(punit,'(2es23.15)') tr00aainput%RZ(1:2,0:Npts)
   
  enddo
  
  close(punit)
  
  ifail = 9 ; call tr00aa( tr00aainput, ifail )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  return

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

  if( Npts.le.0 ) return

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  xyaxis(1:Node) = zero ! this was already initialized in initial;
  NN = Nxyaxis ; tol = odetol * ten ; Lwa = Lrwork ; ic05nbf = 1
  write(ounit,'("pcplot  : ",f10.1," : calling C05NBF ;")') told-tstart 
  call C05NBF( fxyaxis, NN, xyaxis(1:NN), Fxy(1:NN), tol, rwork(1:Lwa), Lwa, ic05nbf ) 
  
  tnow = MPI_WTIME()
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  if( myid.eq.0 ) write(ounit,'("pcplot  : ",f10.1," : Ntrj =",i4," ; Etrj ="i3" ; Npts =",i6," ; odetol =",es12.5," ;")') tnow-tstart, Ntrj, Etrj, Npts, odetol
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  open( punit, file=trim(outplots), status="unknown", form="formatted" )
  
  write( punit,'(2i9)' ) Ntrj + Etrj, Npts
  
  told = tnow
  
  do itrj = 0, Ntrj + Etrj
   
   iota(itrj,1) = itrj * one / Ntrj
   
!#ifdef DEBUG
!  xy(1:Node) = (/ iota(itrj,1), zero /) ! 28 Nov 17;
!  if( myid.eq.0 ) write(punit,'(2es23.15)') xy(1) * minorrad * (/ ellipticity * cos(xy(2)), sin(xy(2)) /) ! 28 Nov 17;
!#else
   xy(1:Node) = xyaxis(1:Node) + itrj * ( (/ one, zero /) - xyaxis(1:2) ) / Ntrj ! set starting point;
   if( myid.eq.0 ) write(punit,'(2es23.15)') xy(1:2) * minorrad * (/ ellipticity             , one        /)
!#endif   
   do ipt = 1, Npts
    
    call pmap( Node, xy(1:Node) )
    
!#ifdef DEBUG
!   if( myid.eq.0 ) write(punit,'(2es23.15)') xy(1) * minorrad * (/ ellipticity * cos(xy(2)), sin(xy(2)) /) ! 28 Nov 17;
!#else
    if( myid.eq.0 ) write(punit,'(2es23.15)') xy(1:2) * minorrad * (/ ellipticity             , one        /)
!   if( myid.eq.0 ) write(punit,'(2es23.15)') ax(1) + xy(1) * minorrad * ellipticity, &
!                                             ax(3) + xy(2) * minorrad * one
!#endif   
    
   enddo ! end of do ipt;
   
!#ifdef DEBUG
!  iota(itrj,2) = xy(2) / (Npts*pi2)
!#else
   iota(itrj,2) = zero
!#endif      
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

subroutine bstfield( zeta, st, Bst )
  
  use globals, only : zero, one, pi2, sqrtmachprec, ounit, myid, ncpu, Ncoils, Ns, coil, minorrad, ellipticity
  
  implicit none  
  
  include "mpif.h"
  
  REAL    :: zeta, st(*), Bst(*)
  
  INTEGER :: icoil, kk, ierr, isurf
  REAL    :: ax(1:3), at(1:3), az(1:3), xx(1:3), xs(1:3), xt(1:3), xz(1:3), v1(1:3), v2(1:3), w1(1:3), w2(1:3)
  REAL    :: xtt(1:3), xtz(1:3), xzz(1:3)
  REAL    :: a, b, c, d, e, f, g, h, i, idet, inversematrix(1:3,1:3)
  REAL    :: dx, dy, dz, rr(1:3), dd(1:3), tx, ty, tz, Bt(1:3), Bi(1:3), Bs(1:3)
  REAL    :: ssmax
  
  ssmax = 1.1 ; isurf = 1
  
  if( st(1).lt.sqrtmachprec .or. st(1).gt.ssmax ) then ; Bst(1:2) = zero ; return
  endif
  
  call knotxx( isurf, st(1), st(2), zeta, ax(1:3), at(1:3), az(1:3), xx(1:3), xs(1:3), xt(1:3), xz(1:3), xtt(1:3), xtz(1:3), xzz(1:3), &
v1(1:3), v2(1:3), w1(1:3), w2(1:3) )

  a = xs(1) ; b = xt(1) ; c = xz(1)
  d = xs(2) ; e = xt(2) ; f = xz(2)
  g = xs(3) ; h = xt(3) ; i = xz(3)
  
! idet = one / ( - c * e * g + b * f * g + c * d * h - a * f * h - b * d * i + a * e * i ) ! this factor will cancel; not required;
  
  inversematrix(1,1) = ( -f * h + e * i ) ; inversematrix(1,2) = (  c * h - b * i ) ; inversematrix(1,3) = ( -c * e + b * f )
  inversematrix(2,1) = (  f * g - d * i ) ; inversematrix(2,2) = ( -c * g + a * i ) ; inversematrix(2,3) = (  c * d - a * f )  
  inversematrix(3,1) = ( -e * g + d * h ) ; inversematrix(3,2) = (  b * g - a * h ) ; inversematrix(3,3) = ( -b * d + a * e )
  
  Bt(1:3) = zero
  
  do icoil = 1, Ncoils
   
!  if( myid.ne.modulo(icoil-1,ncpu) ) cycle ! parallelization loop;
   
   Bi(1:3) = zero
   
   do kk = 0, Ns-1 ! 28 Nov 17;
    
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
   
!  Bi(1:3,icoil) = Bi(1:3,icoil) * coil(icoil)%I
   
   Bt(1:3) = Bt(1:3) + Bi(1:3) * coil(icoil)%I
   
  enddo ! end of do icoil;
  
! call MPI_REDUCE( Bi(1:3), Bt(1:3), 3, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
  
! RlBCAST( Bt(1:3), 3, 0 )
  
  Bs(1:3) = matmul( inversematrix(1:3,1:3), Bt(1:3) )
  
  FATAL( pcplot , abs(Bs(3)).lt.sqrtmachprec, divide by zero )
  
  Bst(1:2) = Bs(1:2) / Bs(3)
  
  return
  
end subroutine bstfield

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

 !call knotxx( isurf, st(1), st(2), zeta, ax(1:3), at(1:3), az(1:3), xx(1:3), xs(1:3), xt(1:3), xz(1:3), xtt, xtz, xzz, v1(1:3), v2(1:3), w1(1:3), w2(1:3) )
  call knotxx( isurf, srho , teta , zeta, ax(1:3), at(1:3), az(1:3), xx(1:3), xs(1:3), xt(1:3), xz(1:3), xtt(1:3), xtz(1:3), xzz(1:3), &
v1(1:3), v2(1:3), w1(1:3), w2(1:3) )
  
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
  
  use globals, only : zero, pi, pi2, odetol, nrotate
  
  implicit none  
  
  include "mpif.h"
  
  INTEGER   :: Node, id02bjf
  REAL      :: xy(1:Node), zeta, zetaend, tol, rwork(1:20*Node)
  CHARACTER :: relabs
  
  external  :: bstfield, bxyfield, D02BJX, D02BJW
  
  zeta = zero ; zetaend = pi2 ; relabs = 'D' ; tol = odetol ; id02bjf = 1
  
!#ifdef DEBUG  
! call D02BJF( zeta, zetaend, Node, xy(1:Node), bstfield, odetol, relabs, D02BJX, D02BJW, rwork(1:20*Node), id02bjf ) ! 28 Nov 17;
! xy(2) = xy(2) + nrotate * pi ! 28 Nov 17;
!#else
  call D02BJF( zeta, zetaend, Node, xy(1:Node), bxyfield, odetol, relabs, D02BJX, D02BJW, rwork(1:20*Node), id02bjf )
  xy(1:2) = xy(1:2) * (-1)**nrotate
!#endif
  
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

subroutine bfield( RpZ, itangent, BRpZ, ifail )
  
  use globals, only : zero, one, pi2, three, sqrtmachprec, ounit, myid, ncpu, Ncoils, Ns, coil, minorrad, ellipticity
  
  implicit none  
  
  include "mpif.h"
  
  REAL    :: RpZ(1:3), BRpZ(1:3,0:3)
  INTEGER :: itangent, ifail
  
  INTEGER :: icoil, kk, ierr
  REAL    :: xx(1:3), R, zeta, czeta, szeta
  REAL    :: inversematrix(1:3,1:3), derivatmatrix(1:3,1:3)
  REAL    :: dx, dy, dz, rr(1:3), dd(1:3), tx, ty, tz, lBi(1:3), Bi(1:3,0:3), Bxyz(1:3,0:3)
  
  R = RpZ(1) ; zeta = RpZ(2) ; czeta = cos(zeta) ; szeta = sin(zeta)

  xx(1) = R * czeta ; xx(2) = R * szeta ; xx(3) = RpZ(3)

  Bxyz(1:3,0:3) = zero ! total magnetic field; 02 Dec 12;

  do icoil = 1, Ncoils
   
  !if( myid.ne.modulo(icoil-1,ncpu) ) cycle ! parallelization loop;
   
   Bi(1:3,0:3) = zero ! i-th magnetic field; 02 Dec 12;
   
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
    
    lBi(1:3) =                (/ ty * dz - tz * dy, tz * dx - tx * dz, tx * dy - ty * dx /) * dd(3)

    Bi(1:3,0) = Bi(1:3,0)                                                                   + lBi(1:3)
    Bi(1:3,1) = Bi(1:3,1) + ( (/       zero       , tz               ,         - ty      /) - lBi(1:3) * three * rr(1) * dx ) * dd(3)
    Bi(1:3,2) = Bi(1:3,2) + ( (/         - tz     ,       zero       , tx                /) - lBi(1:3) * three * rr(1) * dy ) * dd(3)
    Bi(1:3,3) = Bi(1:3,3) + ( (/ ty               ,         - tx     ,       zero        /) - lBi(1:3) * three * rr(1) * dz ) * dd(3)
    
   enddo ! end of do kk;
   
   Bxyz(1:3,0:3) = Bxyz(1:3,0:3) + Bi(1:3,0:3) * coil(icoil)%I
   
  enddo ! end of do icoil;
  
 !call MPI_REDUCE( Bi(1:3), Bxyz(1:3), 3, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )

 !RlBCAST( Bxyz(1:3), 3, 0 )
  
  BRpZ(1:3,0:3) = zero
  
  inversematrix(1,1) =   czeta         ; inversematrix(1,2) =   szeta         ; inversematrix(1,3) = zero
  inversematrix(2,1) = - szeta / R     ; inversematrix(2,2) =   czeta / R     ; inversematrix(2,3) = zero
  inversematrix(3,1) =    zero         ; inversematrix(3,2) =    zero         ; inversematrix(3,3) =  one
  
  BRpZ(1:3,0) = matmul( inversematrix(1:3,1:3), Bxyz(1:3,0) ) ! transform back to cylindrical; 02 Dec 12;
  
  derivatmatrix(1,1) =    zero         ; derivatmatrix(1,2) =    zero         ; derivatmatrix(1,3) = zero
  derivatmatrix(2,1) = + szeta / R**2  ; derivatmatrix(2,2) = - czeta / R**2  ; derivatmatrix(2,3) = zero
  derivatmatrix(3,1) =    zero         ; derivatmatrix(3,2) =    zero         ; derivatmatrix(3,3) = zero
  
  lBi(1:3) =     Bxyz(1:3,1) * czeta + Bxyz(1:3,2) * szeta
  
  BRpZ(1:3,1) = matmul( derivatmatrix(1:3,1:3), Bxyz(1:3,0) ) + matmul( inversematrix(1:3,1:3),   lBi(1:3) )
  
  derivatmatrix(1,1) = - szeta         ; derivatmatrix(1,2) =   czeta         ; derivatmatrix(1,3) = zero
  derivatmatrix(2,1) = - czeta / R     ; derivatmatrix(2,2) = - szeta / R     ; derivatmatrix(2,3) = zero
  derivatmatrix(3,1) =    zero         ; derivatmatrix(3,2) =    zero         ; derivatmatrix(3,3) = zero
  
  lBi(1:3) = ( - Bxyz(1:3,1) * szeta + Bxyz(1:3,2) * czeta                       ) * R
  
  BRpZ(1:3,2) = matmul( derivatmatrix(1:3,1:3), Bxyz(1:3,0) ) + matmul( inversematrix(1:3,1:3), lBi(1:3) )
  
  lBi(1:3) =                                                 Bxyz(1:3,3)
  
  BRpZ(1:3,3) =                                                 matmul( inversematrix(1:3,1:3), lBi(1:3) )
  
  ifail = 0 ! calculation was ok; 02 Dec 12;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  return
  
end subroutine bfield

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

!subroutine iccoil( t, x, y, z, ifail )
!
!  use globals, only : myid
!  
!  implicit none
!
!  include "mpif.h"
!  
!  REAL                 :: t, x(0:1), y(0:1), z(0:1)
!  INTEGER              :: ifail, ierr
!
!  REAL                 :: Rmaj, rmin
!
!  FATAL( pcplot, .true., what is iccoil )
!
!  return
!  
!end subroutine iccoil

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
