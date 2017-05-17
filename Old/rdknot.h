
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
  
  use kmodule, only : zero, one, half, ten, pi, pi2, sqrtmachprec, myid, ncpu, ounit, lunit, &
                      ext, &
                      knotNF, knotsurf, surf, Nteta, Nzeta, bNfp,  &
                      xkc, xks, ykc, yks, zkc, zks, &
                      Mpol, Ntor, &
                      Ntz, isr, trigm, trign, trigwk, ijreal, ijimag, mn, im, in, cfmn, sfmn, efmn, ofmn
  
  implicit none
  
  include "mpif.h"
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOGICAL              :: exist
  INTEGER              :: iostat, astat, ierr
  INTEGER              :: ii, jj, kk, jk, ifail
  REAL                 :: teta, zeta, alfa, dd, ds(1:3), ax(1:3), at(1:3), az(1:3), xx(1:3), xt(1:3), xz(1:3), xo(1:3), lo, lknotsurf !for knotatron;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

   INTEGER, parameter  :: Ndof = 2, Ldf = Ndof, LRR = Ndof * ( Ndof + 1 ) / 2
   INTEGER             :: irevcm, mode, ic05pdf
   REAL                :: gphi, rx, rt, rz, ox, ot, oz, xtol, factor
   REAL                :: alfazeta(1:Ndof), ff(1:Ndof), df(1:Ldf,1:Ndof), diag(1:Ndof), RR(1:LRR), QTf(1:Ndof), rwk(1:Ndof,1:4)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  inquire( file=trim(ext)//".fo.knot", exist=exist )
  if( exist ) then
   if( myid.eq.0 ) then
    open(lunit, file=trim(ext)//".fo.knot", status="old", action='read', iostat=iostat)  
    write(ounit,'("rdknot  : " 10x " : reading ext.fo.knot ;")') 
   endif
  else
   inquire( file=".fo.knot", exist=exist )
   if( exist ) then
    if( myid.eq.0 ) then
     open(lunit, file=trim(ext)//".knot", status="old", action='read', iostat=iostat)  
     write(ounit,'("rdknot  : " 10x " : reading .fo.knot ;")') 
    endif
  !else
  ! inquire( file="knot", exist=exist )
  ! if( exist ) then
  !  if( myid.eq.0 ) then
  !   open(lunit, file="knot", status="old", action='read', iostat=iostat)
  !   write(ounit,'("rdknot  : " 10x " : reading knot ;")') 
  !  endif
  ! endif
   endif
  endif
  
  FATAL( rdknot , .not.exist, input knot does not exist )
  
  if( myid.eq.0 ) read( lunit, *, iostat=iostat ) knotNF!, knotphase ! knotphase is redundant; 18 Apr 17;
  
  IlBCAST( knotNF, 1, 0 )
! RlBCAST( knotphase, 1, 0 ) ! knotphase is redundant; 18 Apr 17;
  
  SALLOCATE( xkc, (0:knotNF), zero )
  SALLOCATE( xks, (0:knotNF), zero )
  SALLOCATE( ykc, (0:knotNF), zero )
  SALLOCATE( yks, (0:knotNF), zero )
  SALLOCATE( zkc, (0:knotNF), zero )
  SALLOCATE( zks, (0:knotNF), zero )

  bNfp = 1
  
! NFcoil = qtorusknot + ptorusknot
!   
! FATAL(rdknot , ptorusknot.gt.qtorusknot, needs attention )
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
  
 !write(ounit,'("rdknot  : " 10x " : lNFcoil, NFcoil=",2i6," ;")') lNFcoil, NFcoil

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
   
 !if( myid.eq.0 ) then
 ! write(ounit,'("rdknot  : " 10x " : xkc=",    999es11.03)') xkc(0:knotNF)
 ! write(ounit,'("rdknot  : " 10x " : xks=",11x,998es11.03)') xks(1:knotNF)
 ! write(ounit,'("rdknot  : " 10x " : ykc=",    999es11.03)') ykc(0:knotNF)
 ! write(ounit,'("rdknot  : " 10x " : yks=",11x,999es11.03)') yks(1:knotNF)
 ! write(ounit,'("rdknot  : " 10x " : zkc=",    999es11.03)') zkc(0:knotNF)
 ! write(ounit,'("rdknot  : " 10x " : zks=",11x,999es11.03)') zks(1:knotNF)
 !endif

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

  SALLOCATE( surf(1)%rr, (1:Nteta*Nzeta), zero ) ! points on surface equally spaced in \theta and \phi; 18 Apr 17;
  SALLOCATE( surf(1)%rz, (1:Nteta*Nzeta), zero ) ! points on surface equally spaced in \theta and \phi; 18 Apr 17;
 
! The center point value was used to discretize grid;

  do ii = 0, Nteta ; teta = ( ii + half ) * pi2 / Nteta
   do jj = 0, Nzeta ; zeta = ( jj + half ) * pi2 / Nzeta
    
    lknotsurf = one ! knotsurf + ellipticity * cos( teta ) ! 11 May 17;

    call knotxx( lknotsurf, teta, zeta, ax, at, az, xx, xt, xz )
    
   !call knotxx( knotsurf, teta, zeta, ax, at, az, xx, xt, xz )
    
    ds(1:3) = -(/ xt(2) * xz(3) - xt(3) * xz(2), xt(3) * xz(1) - xt(1) * xz(3), xt(1) * xz(2) - xt(2) * xz(1) /) ! negative sign = counterclockwise;
    
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
  
!latex The transformation from $(\alpha,\zeta)$ to $(\theta,\phi)$ is performed using
!latex \nag{http://www.nag.co.uk/numeric/FL/manual19/pdf/C05/C05PDF_fl19.pdf}{C05PDF}.

  if( Mpol.gt.0 ) then
   
   if( myid.eq.0 ) then
    
    xo(1:3) = zero ; lo = zero ; alfa = zero
    
    do kk = 0, Nzeta-1 ; zeta = kk * (pi2/bNfp) / Nzeta
     lknotsurf = zero ; call knotxx( lknotsurf, alfa, zeta, ax(1:3), at(1:3), az(1:3), xx(1:3), xt(1:3), xz(1:3) )
     xo(1:3) = xo(1:3) + ax(1:3) * sqrt( sum(az(1:3)*az(1:3)) ) ! this is a very crude integration to find the curve center; 18 Apr 17;
     lo      = lo      +           sqrt( sum(az(1:3)*az(1:3)) ) ! this is a very crude integration to find the curve length; 18 Apr 17;
    enddo
    xo(1:3) = xo(1:3) / lo ; lo = lo * pi2 / Nzeta ! note that this is not broadcast; 18 Apr 17;
    
    write(ounit,'("rdknot  : " 10x " : (xo,yo,zo) = ("f12.7","f12.7","f12.7" ) ; lo ="f12.7" ;")') xo(1:3), lo
    

    alfazeta(1:2) = (/ zero, zero /)  + ten ! note offset in \teta & \zeta; 18 Apr 17;
    
    do kk = 0, Nzeta-1 ; gphi = kk * (pi2/bNfp) / Nzeta ! regularly spaced in \phi ; 27 Apr 17;
     
!    alfazeta(2) = gphi + ten ! reinitialize guess; 26 Apr 17;
     
     do jj = 0, Nteta-1 ; teta = jj *  pi2       / Nteta ; jk = 1 + jj + kk*Nteta ! regularly spaced in \teta; 27 Apr 17;
      
!     alfazeta(1) = teta + ten ! reinitialize guess; 26 Apr 17;
      
      ic05pdf = 1 ; irevcm = 0 ; xtol = sqrtmachprec ; mode = 1 ; factor = one
      
      do
       
       call C05PDF( irevcm, &
                    Ndof, alfazeta(1:Ndof), ff(1:Ndof), df(1:Ldf,1:Ndof), &
                    Ldf, xtol, diag(1:Ndof), mode, factor, RR(1:LRR), LRR, QTf(1:Ndof), rwk(1:Ndof,1:4), ic05pdf )
       
       alfa = mod( alfazeta(1) - ten + pi, pi2 ) ; zeta = alfazeta(2) - ten ! offsets must be consistent with above; 26 Apr 17;

       lknotsurf = one ! knotsurf + ellipticity * cos( alfa ) ! 11 May 17;

       call knotxx( lknotsurf, alfa, zeta, ax(1:3), at(1:3), az(1:3), xx(1:3), xt(1:3), xz(1:3) )
       
      !call knotxx( knotsurf, alfa, zeta, ax(1:3), at(1:3), az(1:3), xx(1:3), xt(1:3), xz(1:3) )
       
       rx = sqrt( xx(1)*xx(1) + xx(2)*xx(2) ) ; rt = ( xx(1)*xt(1) + xx(2)*xt(2) ) / rx ; rz = ( xx(1)*xz(1) + xx(2)*xz(2) ) / rx
       ox = sqrt( ax(1)*ax(1) + ax(2)*ax(2) ) ; ot = ( ax(1)*at(1) + ax(2)*at(2) ) / ox ; oz = ( ax(1)*az(1) + ax(2)*az(2) ) / ox
       
       select case( irevcm )
       case( 0 ) ;!if( ic05pdf.ne.0 ) write(ounit,2000) alfa, zeta, teta, gphi, ff(1:Ndof), ic05pdf
        ;        ; exit
       case( 1 ) ;                   !write(ounit,2000) alfa, zeta, teta, gphi, ff(1:Ndof)
      !case( 2 ) ; ff(1  ) = ( xx(1) - xo(1) ) * sin(gphi) - ( xx(2) - xo(2) ) * cos(gphi)
      ! ;        ; ff(2  ) = ( rx    - ox    ) * sin(teta) - ( xx(3) - ax(3) ) * cos(teta)
      !case( 3 ) ; df(1,1) = ( xt(1)         ) * sin(gphi) - ( xt(2)         ) * cos(gphi)
      ! ;        ; df(1,2) = ( xz(1)         ) * sin(gphi) - ( xz(2)         ) * cos(gphi)
      ! ;        ; df(2,1) = ( rt    - ot    ) * sin(teta) - ( xt(3) - at(3) ) * cos(teta)
      ! ;        ; df(2,2) = ( rz    - oz    ) * sin(teta) - ( xz(3) - az(3) ) * cos(teta)
       case( 2 ) ; ff(1  ) = ( xx(1)         ) * sin(gphi) - ( xx(2)         ) * cos(gphi) ! 27 Apr 17;
        ;        ; ff(2  ) = ( rx    - ox    ) * sin(teta) - ( xx(3) - ax(3) ) * cos(teta) ! 27 Apr 17;
       case( 3 ) ; df(1,1) = ( xt(1)         ) * sin(gphi) - ( xt(2)         ) * cos(gphi) ! 27 Apr 17;
        ;        ; df(1,2) = ( xz(1)         ) * sin(gphi) - ( xz(2)         ) * cos(gphi) ! 27 Apr 17;
        ;        ; df(2,1) = ( rt    - ot    ) * sin(teta) - ( xt(3) - at(3) ) * cos(teta) ! 27 Apr 17;
        ;        ; df(2,2) = ( rz    - oz    ) * sin(teta) - ( xz(3) - az(3) ) * cos(teta) ! 27 Apr 17;
       case default
        FATAL( rdknot , .true., invalid irevcm )
       end select
       
      enddo ! end of do; Nov 12 15;
      
2000  format("rdknot  : " 10x " : (\a,\z)=("f12.08" ,"f12.08" ) ; (\t,\p)=("f12.08" ,"f12.08" ) ; F="2es10.02" ; ":"ic05pdf="i3" ;")
      
      surf(1)%rr(jk) = sqrt( (xx(1)      )**2 + (xx(2)      )**2 )
      surf(1)%rz(jk) =        xx(3)!xo(3) ! 27 Apr 17;
      
     enddo ! end of do jj; Nov 12 15;
    enddo ! end of do kk; Nov 12 15;
    
    ijreal(1:Ntz) = surf(1)%rr(1:Ntz) ! 27 Apr 17;
    ijimag(1:Ntz) = surf(1)%rz(1:Ntz) ! 27 Apr 17;
    
    ifail = 0
    call tfft( Nteta, Nzeta, ijreal(1:Ntz), ijimag(1:Ntz), isr, trigm(1:2*Nteta), trign(1:2*Nzeta), trigwk(1:2*Ntz), &
               mn, im(1:mn), in(1:mn), cfmn(1:mn), sfmn(1:mn), efmn(1:mn), ofmn(1:mn), ifail )

    open(lunit, file=trim(ext)//".fo.RZ", status="unknown")
    do ii = 1, mn
     write(lunit,'(2i6,4es23.15)') in(ii), im(ii), cfmn(ii), sfmn(ii), efmn(ii), ofmn(ii)
    enddo
    close(lunit)
    
    open(lunit, file="input."//trim(ext), status="unknown")
    write(lunit,'("&INDATA")')
    write(lunit,'(" DELT = 0.5, NSTEP = 10000, TCON0 = 1.0,")')
    write(lunit,'(" NS_ARRAY    =      16,      25,      36,      49,")')
    write(lunit,'(" FTOL_ARRAY  = 1.0E-05, 1.0E-07, 1.0E-09, 1.0E-12,")')
    write(lunit,'(" NITER_ARRAY =    1000,   10000,  100000, 1000000,")')
! write(lunit,'(" PRECON_TYPE = 'none'")')
!   write(lunit,'(" PREC2D_THRESHOLD = 1.000000E-19")')
!   write(lunit,'(" LASYM = T, NFP = 1, MPOL ="i4", NTOR ="i3",")')  abs(Mpol), Ntor
    write(lunit,'(" LASYM = F, NFP = 1, MPOL ="i4", NTOR ="i3",")')          8,    4
    write(lunit,'(" PHIEDGE = 1.0, LFREEB = F, GAMMA = 0.0, SPRES_PED = 1.0,")')
!   write(lunit,'(" PRES_SCALE = 0.0, PMASS_TYPE = 'power_series' ")')
    write(lunit,'(" PRES_SCALE = 0.0, AM = 0.0, ")')
    write(lunit,'(" NCURR = 1, CURTOR = 0.0, AC = 0.0,")')
    write(lunit,'(" Raxis ="99(es13.05","))') cfmn(1:Ntor+1)
    write(lunit,'(" Zaxis ="99(es13.05","))') ofmn(1:Ntor+1)
    do ii = 1, mn
     write(lunit,1000) in(ii), im(ii), cfmn(ii), in(ii), im(ii), sfmn(ii), in(ii), im(ii), efmn(ii), in(ii), im(ii), ofmn(ii)
    enddo
    write(lunit,'("/")')
    close(lunit)
    
1000 format(" Rbc("i3","i2") ="es23.15", Rbs("i3","i2") ="es23.15", Zbc("i3","i2") ="es23.15", Zbs("i3","i2") ="es23.15",")

   endif ! end of if( myid.eq.0 ) ; Nov 12 15;
   
   RlBCAST( cfmn(1:mn), mn, 0 ) ! Rbc; 18 Apr 17;
   RlBCAST( sfmn(1:mn), mn, 0 ) ! Rbs; 18 Apr 17;
   RlBCAST( efmn(1:mn), mn, 0 ) ! Zbc; 18 Apr 17;
   RlBCAST( ofmn(1:mn), mn, 0 ) ! Zbs; 18 Apr 17;

   RlBCAST( surf(1)%rr(1:Ntz), Ntz, 0 ) ! 27 Apr 17;
   RlBCAST( surf(1)%rz(1:Ntz), Ntz, 0 ) ! 27 Apr 17;

  endif ! end of if( Mpol.gt.0 ) ; 27 Apr 17;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  return
  
end subroutine rdknot

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
