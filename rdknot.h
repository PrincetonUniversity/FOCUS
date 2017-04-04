
!title (knot) ! Reads knot parameters from file. (Redundant?)

!latex \briefly{The parameters describing a knotted loop are read from file.}

!latex \calledby{\link{notopt}}
!latex \calls{\link{knotxx}}

!latex \tableofcontents

!latex \subsection{overview}
!latex \bi
!latex \item[1.] This routine is called by \link{notopt} only if \inputvar{Itopology} $= 2$.
!latex \item[2.] A knot is a closed curve in three-dimensional space.
!latex \item[3.] The parameters describing the knot, namely \internal{xkc}, \internal{yks} and \internal{zks}, are read from \verb+knot+.
!latex \item[4.] If \inputvar{kspring} $> 0$, the knot is ``relaxed'' to minimize the ``spring energy'' using 
!latex           \nag{http://www.nag.co.uk/numeric/FL/manual19/pdf/D02/d02bhf_fl19.pdf}{D02BHF}.
!latex \ei

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine rdknot
  
  use kmodule, only : zero, one, half, ten, pi2, sqrtmachprec, myid, ncpu, ounit, lunit, &
                      ext, &
                      knotNF, knotsurf, knotphase, surf, Nteta, Nzeta, bNfp,  &
                      xkc, xks, ykc, yks, zkc, zks!, kspring, tauend, itau
  
  implicit none
  
  include "mpif.h"
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOGICAL              :: exist
  INTEGER              :: iostat, astat, ierr
  INTEGER              :: ii, jj
  REAL                 :: teta, zeta, dd, ds(1:3), ax(1:3), at(1:3), az(1:3), xx(1:3), xt(1:3), xz(1:3) !for knotatron;

!  INTEGER, parameter   :: Ndof = 2, Ldf = Ndof, LRR = Ndof*(Ndof+1)/2
!  INTEGER              :: irevcm, mode, ic05pdf
!  REAL                 :: phi, teta, zeta
!  REAL                 :: xtol, factor, xdof(1:Ndof), ff(1:Ndof), df(1:Ldf,1:Ndof), diag(1:Ndof), RR(1:LRR), QTf(1:Ndof), work(1:Ndof,1:4)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  inquire( file=trim(ext)//".op.knot", exist=exist )
  if( exist ) then
   if( myid.eq.0 ) then
    open(lunit, file=trim(ext)//".op.knot", status="old", action='read', iostat=iostat)  
    write(ounit,'("rdknot : " 10x " : reading ext.op.knot ;")') 
   endif
  else
   inquire( file=trim(ext)//".knot", exist=exist )
   if( exist ) then
    if( myid.eq.0 ) then
     open(lunit, file=trim(ext)//".knot", status="old", action='read', iostat=iostat)  
     write(ounit,'("rdcoil : " 10x " : reading ext.knot ; redundant naming convention ;")') 
    endif
   else
    inquire( file="knot", exist=exist )
    if( exist ) then
     if( myid.eq.0 ) then
      open(lunit, file="knot", status="old", action='read', iostat=iostat)
      write(ounit,'("rdcoil : " 10x " : reading knot ;")') 
     endif
    endif
   endif
  endif
  
  FATAL( rdknot, .not.exist, input knot does not exist )

  if( myid.eq.0 ) read( lunit, *, iostat=iostat ) knotNF, knotphase
  
  IlBCAST( knotNF, 1, 0 )
  RlBCAST( knotphase, 1, 0 )
  
  SALLOCATE( xkc, (0:knotNF), zero )
  SALLOCATE( xks, (0:knotNF), zero )
  SALLOCATE( ykc, (0:knotNF), zero )
  SALLOCATE( yks, (0:knotNF), zero )
  SALLOCATE( zkc, (0:knotNF), zero )
  SALLOCATE( zks, (0:knotNF), zero )

  bNfp = 1
  
! NFcoil = qtorusknot + ptorusknot
!   
! FATAL(rdknot, ptorusknot.gt.qtorusknot, needs attention )
!   
! xkc(ptorusknot) =   Rtorusaxis ; xkc(qtorusknot-ptorusknot) = + half * atorusaxis ; xkc(qtorusknot+ptorusknot) = half * atorusaxis
! yks(ptorusknot) =   Rtorusaxis ; yks(qtorusknot-ptorusknot) = - half * atorusaxis ; yks(qtorusknot+ptorusknot) = half * atorusaxis
! zks(qtorusknot) = - atorusaxis
  
! suffix = 'square' & phase = pi2 / 4.0D
  
! xnc=[ 0.0,  1.767311692237854e+00,  6.828825699356500e-14, -9.444974660873413e-01,  3.381930152590584e-14, -2.645063698291779e-01 ]
! yns=[ 0.0,  6.908266991376877e-02, -2.093433065786243e-12, -9.791350364685059e-01,  7.350465722200106e-13, -2.231226861476898e-01 ]
! zns=[ 0.0, -9.695488214492798e-02, -1.645109395931321e-12,  5.140577629208565e-02, -3.898339984154120e-13, -4.254519343376160e-01 ]
  
! suffix = 'granny' & phase = 0.0D
  
! xkc(0:5) = (/ 0.0,  6.653696894645691e-01, -1.656940383509831e-11, -1.043724060058594e+00,  1.069502906375641e-11, -8.075384050607681e-02 /)
! yks(0:5) = (/ 0.0, -5.815748572349548e-01, -2.619011431337359e-11, -1.002896547317505e+00, -3.263673359343855e-11,  1.082220524549484e-01 /)
! zks(0:5) = (/ 0.0,  4.834449507384875e-11, -1.123250573873520e-01,  5.830259885986067e-11, -6.876130104064941e-01, -3.233802114976925e-12 /)
  
 !write(ounit,'("rdknot : " 10x " : lNFcoil, NFcoil=",2i6," ;")') lNFcoil, NFcoil

  if( myid.eq.0 ) read( lunit, *, iostat=iostat ) xkc(0:knotNF)
  if( myid.eq.0 ) read( lunit, *, iostat=iostat ) xks(0:knotNF)
  if( myid.eq.0 ) read( lunit, *, iostat=iostat ) ykc(0:knotNF)
  if( myid.eq.0 ) read( lunit, *, iostat=iostat ) yks(0:knotNF)
  if( myid.eq.0 ) read( lunit, *, iostat=iostat ) zkc(0:knotNF)
  if( myid.eq.0 ) read( lunit, *, iostat=iostat ) zks(0:knotNF)
  
  if( myid.eq.0 ) close(lunit,iostat=iostat)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  RlBCAST( xkc(0:knotNF), knotNF+1, 0 )
  RlBCAST( xks(0:knotNF), knotNF+1, 0 )
  RlBCAST( ykc(0:knotNF), knotNF+1, 0 )
  RlBCAST( yks(0:knotNF), knotNF+1, 0 )
  RlBCAST( zkc(0:knotNF), knotNF+1, 0 )
  RlBCAST( zks(0:knotNF), knotNF+1, 0 )
   
  if( myid.eq.0 ) then
   write(ounit,'("rdknot  : " 10x " : xkc=",    999es11.03)') xkc(0:knotNF)
   write(ounit,'("rdknot  : " 10x " : xks=",11x,998es11.03)') xks(1:knotNF)
   write(ounit,'("rdknot  : " 10x " : ykc=",    999es11.03)') ykc(0:knotNF)
   write(ounit,'("rdknot  : " 10x " : yks=",11x,999es11.03)') yks(1:knotNF)
   write(ounit,'("rdknot  : " 10x " : zkc=",    999es11.03)') zkc(0:knotNF)
   write(ounit,'("rdknot  : " 10x " : zks=",11x,999es11.03)') zks(1:knotNF)
  endif
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!discretize the knot surface data; 2017/03/28; czhu;
  allocate(surf(1:1))
  surf(1)%Nteta = Nteta
  surf(1)%Nzeta = Nzeta
  
  SALLOCATE( surf(1)%xx, (0:Nteta,0:Nzeta), zero )
  SALLOCATE( surf(1)%yy, (0:Nteta,0:Nzeta), zero )
  SALLOCATE( surf(1)%zz, (0:Nteta,0:Nzeta), zero )
  SALLOCATE( surf(1)%nx, (0:Nteta,0:Nzeta), zero )
  SALLOCATE( surf(1)%ny, (0:Nteta,0:Nzeta), zero )
  SALLOCATE( surf(1)%nz, (0:Nteta,0:Nzeta), zero )
  SALLOCATE( surf(1)%ds, (0:Nteta,0:Nzeta), zero )
  SALLOCATE( surf(1)%xt, (0:Nteta,0:Nzeta), zero )
  SALLOCATE( surf(1)%yt, (0:Nteta,0:Nzeta), zero )
  SALLOCATE( surf(1)%zt, (0:Nteta,0:Nzeta), zero )
  SALLOCATE( surf(1)%bnt,(0:Nteta,0:Nzeta), zero )
 
! The center point value was used to discretize grid 
  do ii = 0, Nteta ; teta = ( ii + half ) * pi2 / Nteta
   do jj = 0, Nzeta ; zeta = ( jj + half ) * pi2 / Nzeta

    call knotxx( knotsurf, teta, zeta, ax, at, az, xx, xt, xz )

    ds(1:3) = -(/ xt(2) * xz(3) - xt(3) * xz(2), xt(3) * xz(1) - xt(1) * xz(3), xt(1) * xz(2) - xt(2) * xz(1) /) !careful with the negative sign; means counterclockwise;

    dd = sqrt( sum( ds(1:3)*ds(1:3) ) )
   
    surf(1)%xx(ii,jj) = xx(1)
    surf(1)%yy(ii,jj) = xx(2)
    surf(1)%zz(ii,jj) = xx(3)

    surf(1)%xt(ii,jj) = xt(1)
    surf(1)%yt(ii,jj) = xt(2)
    surf(1)%zt(ii,jj) = xt(3)

    surf(1)%nx(ii,jj) = ds(1) / dd
    surf(1)%ny(ii,jj) = ds(2) / dd
    surf(1)%nz(ii,jj) = ds(3) / dd
    surf(1)%ds(ii,jj) =         dd

   enddo ! end of do jj;
  enddo ! end of do ii;


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!  if( myid.eq.0 ) then
!
!   SALLOCATE( xsurf, (0:NDcoils,0:NDcoils), zero )
!   SALLOCATE( ysurf, (0:NDcoils,0:NDcoils), zero )
!   SALLOCATE( zsurf, (0:NDcoils,0:NDcoils), zero )
!   
!   xdof(1:Ndof) = (/ zero, zero + ten /)
!   
!   do jj = 0, NDcoils ; phi  = jj * (pi2/Nfp) / NDcoils
!    
!    do ii = 0, NDcoils ; teta = ii * pi2       / NDcoils
!     
!     ic05pdf = 1 ; irevcm = 0 ; xtol = sqrtmachprec ; mode = 1 ; factor = one
!     
!     do
!      
!      call C05PDF( irevcm, &
!                   Ndof, xdof(1:Ndof), ff(1:Ndof), df(1:Ldf,1:Ndof), &
!                   Ldf, xtol, diag(1:Ndof), mode, factor, RR(1:LRR), LRR, QTf(1:Ndof), work(1:Ndof,1:4), ic05pdf )
!
!      zeta = xdof(2) - ten ; call knotxx( knotsurf, xdof(1), zeta, ax(1:3), at(1:3), az(1:3), xx(1:3), xt(1:3), xz(1:3) )
!      
!      rx = sqrt( xx(1)*xx(1) + xx(2)*xx(2) ) ; rt = ( xx(1)*xt(1) + xx(2)*xt(2) ) / rx ; rz = ( xx(1)*xz(1) + xx(2)*xz(2) ) / rx
!      ox = sqrt( ax(1)*ax(1) + ax(2)*ax(2) ) ; ot = ( ax(1)*at(1) + ax(2)*at(2) ) / ox ; oz = ( ax(1)*az(1) + ax(2)*az(2) ) / ox
!      
!      select case( irevcm )
!      case( 0 ) ! if( ic05pdf.ne.0 ) write(ounit,2000) phi , teta, irevcm, ic05pdf, xx(1:3)
!       ;        ; exit
!      case( 1 ) ;
!      case( 2 ) ; ff(1  ) =   xx(1)        * sin(phi ) -   xx(2)           * cos(phi )
!       ;        ; ff(2  ) = ( rx    - ox ) * sin(teta) - ( xx(3) - ax(3) ) * cos(teta)
!      case( 3 ) ; df(1,1) =   xt(1)        * sin(phi ) -   xt(2)           * cos(phi )
!       ;        ; df(1,2) =   xz(1)        * sin(phi ) -   xz(2)           * cos(phi )
!       ;        ; df(2,1) = ( rt    - ot ) * sin(teta) - ( xt(3) - at(3) ) * cos(teta)
!       ;        ; df(2,2) = ( rz    - oz ) * sin(teta) - ( xz(3) - az(3) ) * cos(teta)
!      case default
!       FATAL(rdknot, .true., invalid irevcm )
!      end select
!      
!     enddo ! end of do; Nov 12 15;
!     
!2000 format("rdknot : " 10x " : phi ="es23.15" ; teta="es23.15" ; irevcm="i3" ; ic05pdf="i3" ; xx="3es23.15" ;")
!     
!     xsurf(ii,jj) = xx(1)
!     ysurf(ii,jj) = xx(2)
!     zsurf(ii,jj) = xx(3)
!     
!    enddo ! end of do ii; Nov 12 15;
!   enddo ! end of do jj; Nov 12 15;
!   
!  endif ! end of if( myid.eq.0 ) ; Nov 12 15;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  return
  
end subroutine rdknot

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
