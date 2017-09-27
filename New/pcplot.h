!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

!title (pcplot) ! Minimize the target function via a differential flow.

!latex \briefly{The minimization problem is solved by integrating a system of ODEs.}

!latex \calledby{\link{solvers}}
!latex \calls{\link{solvers}}

!latex \section{overview}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

subroutine pcplot
  
  use globals, only : zero, one, half, pi2, ounit, myid, tstart, &
                      minorrad, ellipticity, nrotate, &
                      outplots, punit, &
                      Ntrj, Npts, odetol
  
  implicit none  
  
  include "mpif.h"
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  INTEGER, parameter :: Node = 2
  INTEGER            :: id02bjf, itrj, ipt
  REAL               :: zeta, zetaend, st(1:Node), tol, wk(1:20*Node), tnow, told
  CHARACTER          :: relabs
  
  external           :: bfield, D02BJX, D02BJW
  
  tnow = MPI_WTIME()

  if( myid.eq.0 ) write(ounit,'("pcplot  : ",f10.1," : Ntrj =",i4," ; Npts =",i6," ; odetol =",es12.5," ;")') tnow-tstart, Ntrj, Npts, odetol
  
  open( punit, file=trim(outplots), status="unknown", form="formatted" )
  
  write( punit,'(2i9)' ) Ntrj, Npts
  
  told = tnow

  do itrj = 1, Ntrj-1
   
   st(1:Node) = (/ itrj * one / Ntrj, zero /)
   
   if( myid.eq.0 ) write(punit,'(2es23.15)') st(1)*cos(st(2)), st(1)*sin(st(2))
   
   do ipt = 2, Npts
    
    zeta = zero ; zetaend = pi2 ; relabs = 'D' ; tol = odetol ; id02bjf = 1
    
    call D02BJF( zeta, zetaend, Node, st(1:Node), bfield, odetol, relabs, D02BJX, D02BJW, wk(1:20*Node), id02bjf )

    st(2) = st(2) - nrotate * half * pi2
    
!   if( myid.eq.0 ) write(ounit,'("pcplot  : ", 10x ," : "2i6" : ( s, t ) = (",es13.5,",",es13.5," ) ;")') itrj, ipt, st(1:Node)
    if( myid.eq.0 ) write(punit,'(2es23.15)') st(1) * minorrad * (/ ellipticity * cos(st(2)), sin(st(2)) /)
    
   enddo ! end of do ipt;
   
   tnow = MPI_WTIME()

   if( myid.eq.0 ) write(ounit,'("pcplot  : ",f10.1," : ",i4," : iota =",es23.15," ; time =",f10.2," ;")') tnow-tstart, itrj, st(2) / (Npts*pi2), tnow-told
   
   told = tnow
   
  enddo ! end of do itrj;
  
  close(punit)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  return
  
end subroutine pcplot

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

subroutine bfield( zeta, st, Bst )
  
  use globals, only : zero, one, pi2, sqrtmachprec, ounit, myid, ncpu, Ncoils, Ns, coil, minorrad
  
  implicit none  
  
  include "mpif.h"
  
  REAL    :: zeta, st(*), Bst(*)
  
  INTEGER :: icoil, kk, ierr
  REAL    :: ax(1:3), at(1:3), az(1:3), xx(1:3), xs(1:3), xt(1:3), xz(1:3)
  REAL    :: a, b, c, d, e, f, g, h, i, idet, inversematrix(1:3,1:3)
  REAL    :: dx, dy, dz, rr(1:3), dd(1:3), tx, ty, tz, Bt(1:3), Bi(1:3), Bs(1:3)

  if( st(1).lt.sqrtmachprec .or. st(1).gt.one ) then ; Bst(1:2) = zero ; return
  endif
  
  call knotxx( st(1), st(2), zeta, ax(1:3), at(1:3), az(1:3), xx(1:3), xs(1:3), xt(1:3), xz(1:3) )
  
  a = xs(1) ; b = xt(1) ; c = xz(1)
  d = xs(2) ; e = xt(2) ; f = xz(2)
  g = xs(3) ; h = xt(3) ; i = xz(3)
  
! idet = one / ( - c * e * g + b * f * g + c * d * h - a * f * h - b * d * i + a * e * i ) ! this factor will cancel; not required;
  
  inversematrix(1,1) = ( -f * h + e * i ) ; inversematrix(1,2) = (  c * h - b * i ) ; inversematrix(1,3) = ( -c * e + b * f )
  inversematrix(2,1) = (  f * g - d * i ) ; inversematrix(2,2) = ( -c * g + a * i ) ; inversematrix(2,3) = (  c * d - a * f )  
  inversematrix(3,1) = ( -e * g + d * h ) ; inversematrix(3,2) = (  b * g - a * h ) ; inversematrix(3,3) = ( -b * d + a * e )
  
  Bt(1:3) = zero

  do icoil = 1, Ncoils
   
  !if( myid.ne.modulo(icoil-1,ncpu) ) cycle ! parallelization loop;
   
   Bi(1:3) = zero
   
   do kk = 1, Ns
    
    dx = xx(1) - coil(icoil)%xx(kk) ! distance from evaluation point to curve;
    dy = xx(2) - coil(icoil)%yy(kk)
    dz = xx(3) - coil(icoil)%zz(kk)

    rr(2) = dx**2 + dy**2 + dz**2 ; rr(1) = sqrt(rr(2)) ; rr(3) = rr(1) * rr(2)

    FATAL( pcplot  , abs(rr(3)).lt.sqrtmachprec, divide by zero )

    ;                             ;                     ; dd(3) = one / rr(3)
    
    tx = coil(icoil)%xt(kk) ! tangent to coil (shorthand);
    ty = coil(icoil)%yt(kk)
    tz = coil(icoil)%zt(kk)
    
    Bi(1:3) = Bi(1:3) + (/ ty * dz - tz * dy, tz * dx - tx * dz, tx * dy - ty * dx /) * dd(3)
    
   enddo ! end of do kk;
   
  !Bi(1:3,icoil) = Bi(1:3,icoil) * coil(icoil)%I

   Bt(1:3) = Bt(1:3) + Bi(1:3) * coil(icoil)%I
   
  enddo ! end of do icoil;
  
 !call MPI_REDUCE( Bi(1:3), Bt(1:3), 3, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )

 !RlBCAST( Bt(1:3), 3, 0 )

  Bs(1:3) = matmul( inversematrix(1:3,1:3), Bt(1:3) )
  
  FATAL( pcplot  , abs(Bs(3)).lt.sqrtmachprec, divide by zero )

  Bst(1:2) = Bs(1:2) / Bs(3)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  return
  
end subroutine bfield

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
