
!title (boundary) ! The plasma boundary is read from file.

!latex \briefly{A Fourier representation for the plasma boundary is read from file.}

!latex \calledby{\link{xfocus}}
!l tex \calls{\link{}}

!latex \tableofcontents

!latex \subsection{overview}
!latex \bi
!latex \item[1.] The Fourier harmonics of the plasma boundary are required in \verb+plasma.boundary+, 
!latex and the format of this file is as follows: \begin{verbatim}
!latex bmn        ! integers: number of Fourier harmonics for the plasma boundary;
!latex bNfp       ! integers: number of field periodicity;
!latex nbf        ! integers: number of Fourier harmonics for Bn;
!latex ---------------------------------------------------------
!latex bim(1:bmn) ! integers: poloidal mode identification;
!latex bin(1:bmn) ! integers: toroidal mode identification;
!latex bnim(1:bmn)! integers: poloidal mode identification, for Bn;
!latex bnin(1:bmn)! integers: toroidal mode identification, for Bn;
!latex ---------------------------------------------------------
!latex Rbc(1:bmn) ! real    : cylindrical R cosine harmonics;
!latex Rbs(1:bmn) ! real    : cylindrical R   sine harmonics;
!latex Zbc(1:bmn) ! real    : cylindrical Z cosine harmonics;
!latex Zbs(1:bmn) ! real    : cylindrical Z   sine harmonics;
!latex bns(1:nbf) ! real    : B normal sin harmonics;
!latex bnc(1:nbf) ! real    : B normal cos harmonics; \end{verbatim}
!latex \item[2.] Note that immediately after reading (and broadcasting) \verb+bin+, the field periodicity factor is included, i.e. \verb+bin = bin * bNfp+.
!latex \item[3.] Example of the plasma.boundary file:
!latex  { \begin{verbatim}
!latex #bmn bNfp nbf
!latex 4 2 1
!latex #plasma boundary
!latex # n m Rbc Rbs Zbc Zbs
!latex 0 0  3.00 0.0 0.0  0.00
!latex 0 1  0.30 0.0 0.0 -0.30
!latex 1 0  0.00 0.0 0.0 -0.06
!latex 1 1 -0.06 0.0 0.0 -0.06
!latex #Bn harmonics
!latex # n m bnc bns
!latex 0 0 0.0 0.0
!latex \end{verbatim}
!latex }
!latex \ei

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine surface
  
  use kmodule, only : zero, half, pi2, myid, ncpu, ounit, lunit, ext, Isymmetric, &
                      bmn, bNfp, nbf, bim, bin, bnim, bnin, Rbc, Rbs, Zbc, Zbs, bnc, bns, &
                      Nteta, Nzeta, surf
  
  implicit none
  
  include "mpif.h"
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOGICAL :: exist
  INTEGER :: iosta, iostb, iostc, iostd, ioste, iostf, astat, ierr, ii, jj, imn
  REAL    :: RR(0:2), ZZ(0:2), szeta, czeta, xx(1:3), xt(1:3), xz(1:3), ds(1:3), teta, zeta, arg, carg, sarg, dd, tmp
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  nbf = 0 !initialize
  
  inquire( file="plasma.boundary", exist=exist)
  
  FATAL( surface, .not.exist, plasma.boundary does not exist )
  
  if( myid.eq.0 ) open(lunit, file="plasma.boundary", status='old', action='read')
  
  if( myid.eq.0 ) read(lunit,*) !empty line
  if( myid.eq.0 ) read(lunit,*) bmn, bNfp, nbf !read dimensions
  
  !Broadcast the values
  IlBCAST( bmn , 1, 0 )
  IlBCAST( bNfp, 1, 0 )
  IlBCAST( nbf , 1, 0 )
  
  FATAL( surface, bmn .le.0, invalid )
  FATAL( surface, bNfp.le.0, invalid )
  
  SALLOCATE( bim, (1:bmn), 0 )
  SALLOCATE( bin, (1:bmn), 0 )
  
  SALLOCATE( Rbc, (1:bmn), zero )
  SALLOCATE( Rbs, (1:bmn), zero )
  SALLOCATE( Zbc, (1:bmn), zero )
  SALLOCATE( Zbs, (1:bmn), zero )

  if( myid.eq.0 ) then
   read(lunit,*) !empty line
   read(lunit,*) !empty line
   do imn = 1, bmn ; read(lunit,*) bin(imn), bim(imn), Rbc(imn), Rbs(imn), Zbc(imn), Zbs(imn)
   enddo
  endif  

  IlBCAST( bim(1:bmn), bmn, 0 )
  IlBCAST( bin(1:bmn), bmn, 0 )
  
  if (Isymmetric .eq. 0) bin(1:bmn) = bin(1:bmn) * bNfp  ! Disable periodicity
  
  RlBCAST( Rbc(1:bmn), bmn, 0 )
  RlBCAST( Rbs(1:bmn), bmn, 0 )
  RlBCAST( Zbc(1:bmn), bmn, 0 )
  RlBCAST( Zbs(1:bmn), bmn, 0 )

  if( nbf .gt. 0) then  !read Bn terms
     SALLOCATE( bnim, (1:nbf), 0    )
     SALLOCATE( bnin, (1:nbf), 0    )
     SALLOCATE( bnc , (1:nbf), zero )
     SALLOCATE( bns , (1:nbf), zero )

     if( myid.eq.0 ) then
        read(lunit,*) !empty line
        read(lunit,*) !empty line
        do imn = 1, nbf ; read(lunit,*) bnin(imn), bnim(imn), bnc(imn), bns(imn)
        enddo
     endif

     IlBCAST( bnim(1:nbf), nbf, 0 )
     IlBCAST( bnin(1:nbf), nbf, 0 )

     if (Isymmetric .eq. 0) bnin(1:nbf) = bnin(1:nbf) * bNfp ! Disable periodicity
     
     RlBCAST( bnc(1:nbf) , nbf, 0 )
     RlBCAST( bns(1:nbf) , nbf, 0 )
  endif

  if( myid.eq.0 ) close(lunit,iostat=iosta)
  
  IlBCAST( iosta, 1, 0 )
  
  FATAL( surface, iosta.ne.0, error closing plasma.boundary )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( myid.eq.0 ) then
   write(ounit,'("surface : " 10x " : bmn ="  i06   " ; nbf ="  i06   )') bmn, nbf
#ifdef DEBUG
   write(ounit,'("surface : " 10x " : bim ="99i13   )') bim(1:bmn)
   write(ounit,'("surface : " 10x " : bin ="99i13   )') bin(1:bmn)
   write(ounit,'("surface : " 10x " : Rbc ="99es13.5)') Rbc(1:bmn)
   write(ounit,'("surface : " 10x " : Rbs ="99es13.5)') Rbs(1:bmn)
   write(ounit,'("surface : " 10x " : Zbc ="99es13.5)') Zbc(1:bmn)
   write(ounit,'("surface : " 10x " : Zbs ="99es13.5)') Zbs(1:bmn)
  if (nbf .gt. 0) then
   write(ounit,'("surface : " 10x " : bnim ="99i13  )') bnim(1:nbf)
   write(ounit,'("surface : " 10x " : bnin ="99i13  )') bnin(1:nbf)
   write(ounit,'("surface : " 10x " : bnc ="99es13.5)') bnc(1:nbf)
   write(ounit,'("surface : " 10x " : bns ="99es13.5)') bns(1:nbf)
   endif
#endif
  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  allocate( surf(1:1) ) ! can allow for multiple target boundaries if multiple currents are allowed; 14 Apr 16;
  
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
    
    RR(0:2) = zero ; ZZ(0:2) = zero
    
    do imn = 1, bmn ; arg = bim(imn) * teta - bin(imn) * zeta ; carg = cos(arg) ; sarg = sin(arg)
     
     RR(0) =  RR(0) +     Rbc(imn) * carg + Rbs(imn) * sarg
     ZZ(0) =  ZZ(0) +     Zbc(imn) * carg + Zbs(imn) * sarg
     
     RR(1) =  RR(1) + ( - Rbc(imn) * sarg + Rbs(imn) * carg ) * bim(imn)
     ZZ(1) =  ZZ(1) + ( - Zbc(imn) * sarg + Zbs(imn) * carg ) * bim(imn)
     
     RR(2) =  RR(2) - ( - Rbc(imn) * sarg + Rbs(imn) * carg ) * bin(imn)
     ZZ(2) =  ZZ(2) - ( - Zbc(imn) * sarg + Zbs(imn) * carg ) * bin(imn)
     
    enddo ! end of do imn; 30 Oct 15;
    
    szeta = sin(zeta)
    czeta = cos(zeta)
    
    xx(1:3) = (/   RR(0) * czeta,   RR(0) * szeta, ZZ(0) /)
    xt(1:3) = (/   RR(1) * czeta,   RR(1) * szeta, ZZ(1) /)
    xz(1:3) = (/   RR(2) * czeta,   RR(2) * szeta, ZZ(2) /) + (/ - RR(0) * szeta,   RR(0) * czeta, zero  /)

    ds(1:3) = -(/ xt(2) * xz(3) - xt(3) * xz(2), xt(3) * xz(1) - xt(1) * xz(3), xt(1) * xz(2) - xt(2) * xz(1) /)!careful with the negative sign; means counterclockwise;

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

   enddo ! end of do jj; 14 Apr 16;
  enddo ! end of do ii; 14 Apr 16;

  !calculate target Bn with input harmonics; 05 Jan 17;
  if(nbf .gt. 0) then
     


     do ii = 0, Nteta ; teta = ( ii + half ) * pi2 / Nteta
        do jj = 0, Nzeta ; zeta = ( jj + half ) * pi2 / Nzeta

           tmp = 0.0
           do imn = 1, nbf

           !if( ii.eq.0 .and. jj.eq.0 ) write(ounit,'(2i6,2es23.15)') bnim(imn), bnin(imn), bnc(imn), bns(imn)
           !if( ii.eq.8 .and. jj.eq.8 ) write(ounit,'(2i6,2es23.15)') bnim(imn), bnin(imn), bnc(imn), bns(imn)

              arg = bnim(imn) * teta - bnin(imn) * zeta ; carg = cos(arg) ; sarg = sin(arg)
              tmp = tmp + bnc(imn) * carg + bns(imn) * sarg
              !if (ii.eq.0 .and. jj .eq. 0) write(*,*) tmp
              !surf(1)%bnt(ii, jj) = surf(1)%bnt(ii, jj) + bnc(imn) * carg + bns(imn) * sarg !there is bug for this; 05 Jan 17;
           enddo
           surf(1)%bnt(ii, jj) = tmp

          !if( ii.eq.0 .and. jj.eq.0 ) write(ounit,*) surf(1)%bnt(ii, jj)
          !if( ii.eq.8 .and. jj.eq.8 ) write(ounit,*) teta, zeta, surf(1)%bnt(ii, jj)

        enddo
     enddo

  endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  return
  
end subroutine surface

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine surfcoord( theta, zeta, r, z)
  use kmodule, only: zero, bmn, bNfp, bim, bin, Rbc, Rbs, Zbc, Zbs
  implicit none
  include "mpif.h"

  REAL, INTENT(in ) :: theta, zeta
  REAL, INTENT(out) :: r, z

  INTEGER           :: imn
  REAL              :: arg, carg, sarg
  
  if( .not. allocated(bim) ) STOP  "Please allocate surface data first!"

  r = zero; z = zero
  
  do imn = 1, bmn

     arg = bim(imn) * theta - bin(imn) * zeta
     R =  R + Rbc(imn) * cos(arg) + Rbs(imn) * sin(arg)
     Z =  Z + Zbc(imn) * cos(arg) + Zbs(imn) * sin(arg)

  enddo

  return
end subroutine surfcoord

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!t tle (brief) ! Description.

!latex \briefly{Extended description.}

!latex \calledby{\link{notopt}, \link{plassf}, \link{rdknot}, \link{windsf}}
!latex \calls{\link{}}

!latex \tableofcontents

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

#ifdef OLD

subroutine knotxx( aa, teta, zeta, ax, at, az, xx, xt, xz )
  
  use kmodule, only : zero, one, pi2, small, myid, ounit, &
                      Itopology, NFcoil, knotphase, &
                      xkc, xks, ykc, yks, zkc, zks
  
  implicit none
  
  include "mpif.h"
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  REAL                 :: aa, teta, zeta, ax(1:3), at(1:3), az(1:3), xx(1:3), xt(1:3), xz(1:3)
  
  INTEGER              :: ierr, p(1:3), q(1:3), mm
  REAL                 :: cqz, sqz, cpz, spz, RR(0:3), ZZ(0:3), x0(1:3), x1(1:3), x2(1:3), x3(1:3), a0, a1, a2, b0, b1, carg, sarg
  REAL                 :: tt(1:3), td(1:3), dd(1:3), xa, ya, za, ff, nn(1:3), nd(1:3), bb(1:3), bd(1:3)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  select case( Itopology )
   
!  case( 1 ) ! torus knot; Nov 12 15;
!   
!   p(1:3) = (/ ptorusknot, ptorusknot**2, ptorusknot**3 /) ! Nov 12 15;
!   q(1:3) = (/ qtorusknot, qtorusknot**2, qtorusknot**3 /) ! Nov 12 15;
!   
!   cqz = cos( q(1)*zeta ) ; cpz = cos( p(1)*zeta ) ! Nov 12 15;
!   sqz = sin( q(1)*zeta ) ; spz = sin( p(1)*zeta ) ! Nov 12 15;
!   
!   RR(0) = Rtorusaxis + atorusaxis * cqz        ! Nov 12 15;
!   RR(1) =            - atorusaxis * sqz * q(1) ! Nov 12 15;
!   RR(2) =            - atorusaxis * cqz * q(2) ! Nov 12 15;
!   RR(3) =            + atorusaxis * sqz * q(3) ! Nov 12 15;
!   
!   ZZ(0) =            - atorusaxis * sqz        ! Nov 12 15;
!   ZZ(1) =            - atorusaxis * cqz * q(1) ! Nov 12 15;
!   ZZ(2) =            + atorusaxis * sqz * q(2) ! Nov 12 15;
!   ZZ(3) =            + atorusaxis * cqz * q(3) ! Nov 12 15;
!   
!   x0(1:3) = (/ + RR(0) * cpz, + RR(0) * spz, + ZZ(0) /) ! Nov 12 15;
!   
!   x1(1:3) = (/ + RR(1) * cpz, + RR(1) * spz, + ZZ(1) /) & 
!           + (/ - RR(0) * spz, + RR(0) * cpz,   zero  /) * p(1) ! Nov 12 15;
!
!   x2(1:3) = (/ + RR(2) * cpz, + RR(2) * spz, + ZZ(2) /)            & 
!           + (/ - RR(1) * spz, + RR(1) * cpz,   zero  /) * p(1) * 2 &
!           + (/ - RR(0) * cpz, - RR(0) * spz,   zero  /) * p(2)
!
!   x3(1:3) = (/ + RR(3) * cpz, + RR(3) * spz, + ZZ(3) /)            & 
!           + (/ - RR(2) * spz, + RR(2) * cpz,   zero  /) * p(1) * 3 &
!           + (/ - RR(1) * cpz, - RR(1) * spz,   zero  /) * p(2) * 3 &
!           + (/ + RR(0) * spz, - RR(0) * cpz,   zero  /) * p(3)       ! Nov 12 15;
   
  case( 1:2 ) ! arbitrary knot represented using Fourier series; Nov 12 15;
   
   FATAL( knotxx, abs(knotphase).gt.small, need to revise phase )
   
   x0(1:3) = (/ xkc(0), ykc(0), zkc(0) /)
   x1(1:3) = (/ zero  , zero  , zero   /)
   x2(1:3) = (/ zero  , zero  , zero   /)
   x3(1:3) = (/ zero  , zero  , zero   /)
   
   if( aa.gt.zero ) then ! will need additional derivatives to construct normal and binormal; 14 Apr 16;
    
    do mm = 1, NFcoil ; carg = cos(mm*zeta) ; sarg = sin(mm*zeta)
     x0(1:3) = x0(1:3) + ( + (/ xkc(mm), ykc(mm), zkc(mm) /) * carg + (/ xks(mm), yks(mm), zks(mm) /) * sarg ) 
     x1(1:3) = x1(1:3) + ( - (/ xkc(mm), ykc(mm), zkc(mm) /) * sarg + (/ xks(mm), yks(mm), zks(mm) /) * carg ) * mm   
     x2(1:3) = x2(1:3) + ( - (/ xkc(mm), ykc(mm), zkc(mm) /) * carg - (/ xks(mm), yks(mm), zks(mm) /) * sarg ) * mm**2
     x3(1:3) = x3(1:3) + ( + (/ xkc(mm), ykc(mm), zkc(mm) /) * sarg - (/ xks(mm), yks(mm), zks(mm) /) * carg ) * mm**3
    enddo
    
   else
    
    do mm = 1, NFcoil ; carg = cos(mm*zeta) ; sarg = sin(mm*zeta)
     x0(1:3) = x0(1:3) + ( + (/ xkc(mm), ykc(mm), zkc(mm) /) * carg + (/ xks(mm), yks(mm), zks(mm) /) * sarg )       
     x1(1:3) = x1(1:3) + ( - (/ xkc(mm), ykc(mm), zkc(mm) /) * sarg + (/ xks(mm), yks(mm), zks(mm) /) * carg ) * mm   
    !x2(1:3) = x2(1:3) + ( - (/ xkc(mm), ykc(mm), zkc(mm) /) * carg - (/ xks(mm), yks(mm), zks(mm) /) * sarg ) * mm**2
    !x3(1:3) = x3(1:3) + ( + (/ xkc(mm), ykc(mm), zkc(mm) /) * sarg - (/ xks(mm), yks(mm), zks(mm) /) * carg ) * mm**3
    enddo
    
   endif ! end of if( aa.gt.zero ) ; 14 Apr 16;
   
  case default
   
   FATAL( knotxx, .true., selected Itopology is not supported )
   
  end select
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  ax(1:3) = x0(1:3)                                                        ! Nov 12 15;
  at(1:3) = zero                                                           ! Nov 12 15;
  az(1:3) = x1(1:3)                                                        ! Nov 12 15;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( aa.gt.zero ) then
   
   a0      = sqrt( x1(1)*x1(1) + x1(2)*x1(2) + x1(3)*x1(3) )                                   ! Nov 12 15;
   a1      =     ( x1(1)*x2(1) + x1(2)*x2(2) + x1(3)*x2(3) ) / a0                              ! Nov 12 15;
   a2      =     ( x2(1)*x2(1) + x2(2)*x2(2) + x2(3)*x2(3)   &
               +   x1(1)*x3(1) + x1(2)*x3(2) + x1(3)*x3(3) - a1 * a1 ) / a0                    ! Nov 12 15;

   tt(1:3) =   x1(1:3)                                                / a0                     ! Nov 12 15;
   td(1:3) = ( x2(1:3) - tt(1:3) * a1                               ) / a0                     ! Nov 12 15;
   dd(1:3) = ( x3(1:3) - td(1:3) * a1 - tt(1:3) * a2 - td(1:3) * a1 ) / a0                     ! Nov 12 15;
   
   xa = ( x2(1) - tt(1) * a1 ) ! Nov 12 15;
   ya = ( x2(2) - tt(2) * a1 ) ! Nov 12 15;
   za = ( x2(3) - tt(3) * a1 ) ! Nov 12 15;
   
   ff = sqrt( xa**2 + ya**2 + za**2 ) ! Nov 12 15;
   
   b0 = ff / a0 ! Nov 12 15;
   
   b1 = ( ( xa * ( x3(1) - td(1) * a1 - tt(1) * a2 ) &
          + ya * ( x3(2) - td(2) * a1 - tt(2) * a2 ) &
          + za * ( x3(3) - td(3) * a1 - tt(3) * a2 ) ) / ff - b0 * a1 ) / a0 ! Nov 12 15;

   nn(1:3) =   td(1:3)                  / b0                                                   ! Nov 12 15;
   nd(1:3) = ( dd(1:3) - nn(1:3) * b1 ) / b0                                                   ! Nov 12 15;
   
   bb(1:3) = (/ tt(2)*nn(3)-tt(3)*nn(2), tt(3)*nn(1)-tt(1)*nn(3), tt(1)*nn(2)-tt(2)*nn(1) /)   ! Nov 12 15;
   bd(1:3) = (/ td(2)*nn(3)-td(3)*nn(2), td(3)*nn(1)-td(1)*nn(3), td(1)*nn(2)-td(2)*nn(1) /) &
           + (/ tt(2)*nd(3)-tt(3)*nd(2), tt(3)*nd(1)-tt(1)*nd(3), tt(1)*nd(2)-tt(2)*nd(1) /)   ! Nov 12 15;
   
   xx(1:3) = x0(1:3) + aa * (   cos(teta) * nn(1:3) - sin(teta) * bb(1:3) ) ! aa is minor radius;
   xt(1:3) = zero    + aa * ( - sin(teta) * nn(1:3) - cos(teta) * bb(1:3) )
   xz(1:3) = x1(1:3) + aa * (   cos(teta) * nd(1:3) - sin(teta) * bd(1:3) )
   
  else
   
   xx(1:3) = zero
   xt(1:3) = zero
   xz(1:3) = zero
   
  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  return
  
end subroutine knotxx

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

#endif

