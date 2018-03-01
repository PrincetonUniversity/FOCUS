
!title (coils) ! Fourier description of coils.

!latex \briefly{The Fourier harmonics of the coils are read from file.}

!latex \calledby{\link{notopt}}
!l tex \calls{\link{}}

!latex \tableofcontents

!latex \subsection{overview}
!latex \bi
!latex \item[1.] The initial currents and geometries of the coil currents are read from \verb+.fo.coil.xxx+.
!latex \item[2.] The format of these files will be described in more detail upon request.
!latex \ei

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine rdcoils
  use kmodule, only : zero, half, one, two, pi, pi2, myid, ounit, lunit, ncpu, sqrtmachprec, &
       surf, Nteta, NFcoil, NDcoil, Ncoils, Ndof, coil, cmt, smt, itime, Ntauout, Tdof, &
       Linitialize, Itopology, Lnormalize, Rmaj, rmin, Ic, Io, Iw, Lc, Lo, Lw, Nfixcur, Nfixgeo,&
       coilspace, ext, coilsX, coilsY, coilsZ, coilsI, Nseg, bsconstant,antibscont, &
       Loptimize, weight_eqarc, deriv, norm, Inorm, Gnorm
  implicit none

  include "mpif.h"

  LOGICAL   :: exist
  INTEGER   :: ierr, astat, ii, icoil, mm, jj, ifail, maxnseg, idof, tmp
  REAL      :: zeta, tt, totalcurrent, r1, r2, z1, z2, z0
  REAL      :: ax(1:3), at(1:3), az(1:3), xx(1:3), xt(1:3), xz(1:3) !for knotatron;
  CHARACTER :: suffix*3, coilsfile*40

  Nfixcur = 0 ! fixed coil current number
  Nfixgeo = 0 ! fixed coil geometry number
  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  select case( Linitialize )

     !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!!$  case(: -2 )
!!$
!!$     FATAL( rdcoils, Itopology .ne. 1, this is for knotatron only.)
!!$
!!$     Ncoils = abs(Linitialize)
!!$     if (myid .eq. 0) then
!!$        write(ounit,'("rdcoils : " 10x " : Initialize" I3" coils for knotatrans.")') Ncoils
!!$     endif
!!$
!!$     if( .not. allocated(coilsX) ) then
!!$        SALLOCATE( coilsX, (1:NDcoil,1:Ncoils), zero )    !allocate on other nodes
!!$        SALLOCATE( coilsY, (1:NDcoil,1:Ncoils), zero )
!!$        SALLOCATE( coilsZ, (1:NDcoil,1:Ncoils), zero )
!!$     endif
!!$
!!$     !preparing for coilsX, coilsY and coilsZ;
!!$     do ii = 1, Ncoils
!!$        zeta = (ii-1)*pi2/Ncoils ! toroidal angle;
!!$        do jj = 1, NDcoil
!!$           tt = (jj-1)*pi2/NDcoil !poloidal angle;
!!$
!!$           FATAL( rdcoils, rmin.lt.one, definition of rmin has changed ask SRH )  ! 11 May 17;
!!$
!!$           call knotxx( rmin, tt, zeta, ax, at, az, xx, xt, xz )
!!$
!!$           coilsX(jj, ii) = xx(1)
!!$           coilsY(jj, ii) = xx(2)
!!$           coilsZ(jj, ii) = xx(3)
!!$
!!$        enddo
!!$     enddo
!!$
!!$     allocate( coil(1:Ncoils) )
!!$
!!$     icoil = 0
!!$     do icoil = 1, Ncoils
!!$
!!$        coil(icoil)%N  =  NFcoil
!!$        coil(icoil)%D  =  NDcoil
!!$        coil(icoil)%I  =  Io 
!!$        coil(icoil)%Ic =  Ic
!!$        coil(icoil)%Io =  Io
!!$        coil(icoil)%Iw =  Iw
!!$
!!$        coil(icoil)%L  =  Lo ! irrelevant until re-computed; 14 Apr 16;
!!$        coil(icoil)%Lc =  Lc
!!$        coil(icoil)%Lo =  Lo
!!$        coil(icoil)%Lw =  Lw
!!$
!!$        FATAL( rdcoils, coil(icoil)%Ic.lt.0 .or. coil(icoil)%Ic.gt.1, illegal )
!!$        FATAL( rdcoils, coil(icoil)%Iw.lt.zero                      , illegal )
!!$        FATAL( rdcoils, coil(icoil)%Lc.lt.0                         , illegal )
!!$        FATAL( rdcoils, coil(icoil)%Lo.lt.zero                      , illegal )
!!$        FATAL( rdcoils, coil(icoil)%Lw.lt.zero                      , illegal )
!!$
!!$        SALLOCATE( coil(icoil)%xc, (0:NFcoil), zero )
!!$        SALLOCATE( coil(icoil)%xs, (0:NFcoil), zero )
!!$        SALLOCATE( coil(icoil)%yc, (0:NFcoil), zero )
!!$        SALLOCATE( coil(icoil)%ys, (0:NFcoil), zero )
!!$        SALLOCATE( coil(icoil)%zc, (0:NFcoil), zero )
!!$        SALLOCATE( coil(icoil)%zs, (0:NFcoil), zero )
!!$
!!$        SALLOCATE( coil(icoil)%xx, (0:coil(icoil)%D), zero )
!!$        SALLOCATE( coil(icoil)%yy, (0:coil(icoil)%D), zero )
!!$        SALLOCATE( coil(icoil)%zz, (0:coil(icoil)%D), zero )
!!$        SALLOCATE( coil(icoil)%xt, (0:coil(icoil)%D), zero )
!!$        SALLOCATE( coil(icoil)%yt, (0:coil(icoil)%D), zero )
!!$        SALLOCATE( coil(icoil)%zt, (0:coil(icoil)%D), zero )
!!$        SALLOCATE( coil(icoil)%xa, (0:coil(icoil)%D), zero )
!!$        SALLOCATE( coil(icoil)%ya, (0:coil(icoil)%D), zero )
!!$        SALLOCATE( coil(icoil)%za, (0:coil(icoil)%D), zero )
!!$
!!$        write(coil(icoil)%name,'(i3.3"th-coil")') icoil
!!$
!!$        if(myid .ne. modulo(icoil-1, ncpu)) cycle
!!$
!!$        call Fourier( coilsX(1:NDcoil,icoil), coil(icoil)%xc, coil(icoil)%xs, NDcoil, NFcoil )
!!$        call Fourier( coilsY(1:NDcoil,icoil), coil(icoil)%yc, coil(icoil)%ys, NDcoil, NFcoil )
!!$        call Fourier( coilsZ(1:NDcoil,icoil), coil(icoil)%zc, coil(icoil)%zs, NDcoil, NFcoil )
!!$
!!$        if(coil(icoil)%Ic .eq. 0) Nfixcur = Nfixcur + 1
!!$        if(coil(icoil)%Lc .eq. 0) Nfixgeo = Nfixgeo + 1
!!$
!!$     enddo
!!$
!!$     do icoil = 1, NCoils
!!$        RlBCAST( coil(icoil)%xc(0:NFcoil) , 1+NFcoil ,  modulo(icoil-1, ncpu) )
!!$        RlBCAST( coil(icoil)%xs(0:NFcoil) , 1+NFcoil ,  modulo(icoil-1, ncpu) )
!!$        RlBCAST( coil(icoil)%yc(0:NFcoil) , 1+NFcoil ,  modulo(icoil-1, ncpu) )
!!$        RlBCAST( coil(icoil)%ys(0:NFcoil) , 1+NFcoil ,  modulo(icoil-1, ncpu) )
!!$        RlBCAST( coil(icoil)%zc(0:NFcoil) , 1+NFcoil ,  modulo(icoil-1, ncpu) )
!!$        RlBCAST( coil(icoil)%zs(0:NFcoil) , 1+NFcoil ,  modulo(icoil-1, ncpu) )
!!$     enddo
!!$
!!$     FATAL(rdcoils, .not. allocated(coilsX), coils file data not allocated)
!!$     deallocate(coilsX, coilsY, coilsZ)
! comment out on 20180228 for NAG incompative
     !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  case(-1 )

     if (myid .eq. 0) then
        write(ounit,'("rdcoils : " 10x " : Reading coils data from coils."A)') trim(ext)
        coilsfile = 'coils.'//trim(ext)
        call readcoils(coilsfile, maxnseg)
        write(ounit,'("rdcoils : " 10x " : Read ",i6," coils")') Ncoils
     endif

     IlBCAST( Ncoils   ,      1, 0 )
     IlBCAST( maxnseg  ,      1, 0 )

     if( .not. allocated(coilsX) ) then
        SALLOCATE( coilsX, (1:maxnseg,1:Ncoils), zero )    !allocate on other nodes
        SALLOCATE( coilsY, (1:maxnseg,1:Ncoils), zero )
        SALLOCATE( coilsZ, (1:maxnseg,1:Ncoils), zero )
        SALLOCATE( coilsI, (          1:Ncoils), zero )
        SALLOCATE(   Nseg, (          1:Ncoils),    0 )
     endif

     RlBCAST( coilsX, maxnseg*Ncoils, 0 )
     RlBCAST( coilsY, maxnseg*Ncoils, 0 )
     RlBCAST( coilsZ, maxnseg*Ncoils, 0 )
     RlBCAST( coilsI,         Ncoils, 0 )   
     IlBCAST( Nseg  ,         Ncoils, 0 )

     !call identfy  !in identfy.h

     allocate( coil(1:Ncoils) )

     icoil = 0
     do icoil = 1, Ncoils

        coil(icoil)%N  =  NFcoil
        coil(icoil)%D  =  NDcoil
        coil(icoil)%I  =  coilsI(icoil) * antibscont       ! a scale constant on current to eliminate the large order of currents; only for coils file yet; 07/14/2016
        coil(icoil)%Ic =  Ic
        coil(icoil)%Io =  Io
        coil(icoil)%Iw =  Iw

        coil(icoil)%L  =  Lo ! irrelevant until re-computed; 14 Apr 16;
        coil(icoil)%Lc =  Lc
        coil(icoil)%Lo =  Lo
        coil(icoil)%Lw =  Lw

        FATAL( rdcoils, coil(icoil)%Ic.lt.0 .or. coil(icoil)%Ic.gt.1, illegal )
        FATAL( rdcoils, coil(icoil)%Iw.lt.zero                      , illegal )
        FATAL( rdcoils, coil(icoil)%Lc.lt.0                         , illegal )
        FATAL( rdcoils, coil(icoil)%Lo.lt.zero                      , illegal )
        FATAL( rdcoils, coil(icoil)%Lw.lt.zero                      , illegal )

        SALLOCATE( coil(icoil)%xc, (0:NFcoil), zero )
        SALLOCATE( coil(icoil)%xs, (0:NFcoil), zero )
        SALLOCATE( coil(icoil)%yc, (0:NFcoil), zero )
        SALLOCATE( coil(icoil)%ys, (0:NFcoil), zero )
        SALLOCATE( coil(icoil)%zc, (0:NFcoil), zero )
        SALLOCATE( coil(icoil)%zs, (0:NFcoil), zero )

        SALLOCATE( coil(icoil)%xx, (0:coil(icoil)%D), zero )
        SALLOCATE( coil(icoil)%yy, (0:coil(icoil)%D), zero )
        SALLOCATE( coil(icoil)%zz, (0:coil(icoil)%D), zero )
        SALLOCATE( coil(icoil)%xt, (0:coil(icoil)%D), zero )
        SALLOCATE( coil(icoil)%yt, (0:coil(icoil)%D), zero )
        SALLOCATE( coil(icoil)%zt, (0:coil(icoil)%D), zero )
        SALLOCATE( coil(icoil)%xa, (0:coil(icoil)%D), zero )
        SALLOCATE( coil(icoil)%ya, (0:coil(icoil)%D), zero )
        SALLOCATE( coil(icoil)%za, (0:coil(icoil)%D), zero )

        write(coil(icoil)%name,'(i3.3"th-coil")') icoil

        if(myid .ne. modulo(icoil-1, ncpu)) cycle

        call Fourier( coilsX(1:nseg(icoil),icoil), coil(icoil)%xc, coil(icoil)%xs, nseg(icoil), NFcoil )
        call Fourier( coilsY(1:nseg(icoil),icoil), coil(icoil)%yc, coil(icoil)%ys, nseg(icoil), NFcoil )
        call Fourier( coilsZ(1:nseg(icoil),icoil), coil(icoil)%zc, coil(icoil)%zs, nseg(icoil), NFcoil )

        !change theta direction if needed; 2016/08/04
        !comment out on 01/23/2017
!!$    if ( coil(icoil)%zs(1) .gt. zero ) then
!!$       coil(icoil)%xs = - coil(icoil)%xs
!!$       coil(icoil)%ys = - coil(icoil)%ys
!!$       coil(icoil)%zs = - coil(icoil)%zs
!!$#ifdef DEBUG
!!$       TMPOUT("reset theta direction after reading coils")
!!$#endif
!!$    endif

        if(coil(icoil)%Ic .eq. 0) Nfixcur = Nfixcur + 1
        if(coil(icoil)%Lc .eq. 0) Nfixgeo = Nfixgeo + 1

     enddo

     do icoil = 1, NCoils
        RlBCAST( coil(icoil)%xc(0:NFcoil) , 1+NFcoil ,  modulo(icoil-1, ncpu) )
        RlBCAST( coil(icoil)%xs(0:NFcoil) , 1+NFcoil ,  modulo(icoil-1, ncpu) )
        RlBCAST( coil(icoil)%yc(0:NFcoil) , 1+NFcoil ,  modulo(icoil-1, ncpu) )
        RlBCAST( coil(icoil)%ys(0:NFcoil) , 1+NFcoil ,  modulo(icoil-1, ncpu) )
        RlBCAST( coil(icoil)%zc(0:NFcoil) , 1+NFcoil ,  modulo(icoil-1, ncpu) )
        RlBCAST( coil(icoil)%zs(0:NFcoil) , 1+NFcoil ,  modulo(icoil-1, ncpu) )
     enddo

     FATAL(rdcoils, .not. allocated(coilsX), coils file data not allocated)
     deallocate(coilsX, coilsY, coilsZ, coilsI)
     !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  case( 0 )

     if( myid==0 ) then  !get file number;
        inquire( file=trim(ext)//".focus", exist=exist )
        if ( .not. exist ) then
           STOP "ext.focus NOT existed"
           call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
        endif
        open( lunit, file=trim(ext)//".focus", status="old" )
        read( lunit,*)
        read( lunit,*) Ncoils
        write(ounit,'("rdcoils : " 10X " : identified "i3" coils in ext.focus ;")') Ncoils
     endif

     IlBCAST( Ncoils        ,    1,  0 )
     allocate(    coil(1:Ncoils) )

     if( myid==0 ) then
        do icoil = 1, Ncoils
           read( lunit,*)
           read( lunit,*)
           read( lunit,*) tmp, coil(icoil)%name  !tmp space is reserved for new version;
           read( lunit,*)
           read( lunit,*) coil(icoil)%D, coil(icoil)%I, coil(icoil)%Ic, &
                & coil(icoil)%L, coil(icoil)%Lc, coil(icoil)%Lo
           FATAL( coilfou, coil(icoil)%D < 0                         , illegal )
           FATAL( coilfou, coil(icoil)%Ic < 0 .or. coil(icoil)%Ic > 1, illegal )
           FATAL( coilfou, coil(icoil)%Lc < 0 .or. coil(icoil)%Lc > 2, illegal )
           FATAL( coilfou, coil(icoil)%L  < zero                     , illegal )
           FATAL( coilfou, coil(icoil)%Lc < zero                     , illegal )
           FATAL( coilfou, coil(icoil)%Lo < zero                     , illegal )
           read( lunit,*)
           read( lunit,*) coil(icoil)%N
           FATAL( coilfou, coil(icoil)%N  < 0                    , illegal )
           SALLOCATE( coil(icoil)%xc, (0:NFcoil), zero )
           SALLOCATE( coil(icoil)%xs, (0:NFcoil), zero )
           SALLOCATE( coil(icoil)%yc, (0:NFcoil), zero )
           SALLOCATE( coil(icoil)%ys, (0:NFcoil), zero )
           SALLOCATE( coil(icoil)%zc, (0:NFcoil), zero )
           SALLOCATE( coil(icoil)%zs, (0:NFcoil), zero )
           read( lunit,*)
           read( lunit,*) coil(icoil)%xc(0:min(coil(icoil)%N, NFcoil))
           read( lunit,*) coil(icoil)%xs(0:min(coil(icoil)%N, NFcoil))
           read( lunit,*) coil(icoil)%yc(0:min(coil(icoil)%N, NFcoil))
           read( lunit,*) coil(icoil)%ys(0:min(coil(icoil)%N, NFcoil))
           read( lunit,*) coil(icoil)%zc(0:min(coil(icoil)%N, NFcoil))
           read( lunit,*) coil(icoil)%zs(0:min(coil(icoil)%N, NFcoil))

        enddo !end do icoil;

        close( lunit )
     endif ! end of if( myid==0 );

     do icoil = 1, Ncoils

        ClBCAST( coil(icoil)%name         , 10       ,  0 )
        IlBCAST( coil(icoil)%D            , 1        ,  0 )
        RlBCAST( coil(icoil)%I            , 1        ,  0 )
        IlBCAST( coil(icoil)%Ic           , 1        ,  0 )
        RlBCAST( coil(icoil)%L            , 1        ,  0 )
        IlBCAST( coil(icoil)%Lc           , 1        ,  0 )
        RlBCAST( coil(icoil)%Lo           , 1        ,  0 )
        IlBCAST( coil(icoil)%N            , 1        ,  0 )

        coil(icoil)%N = NFcoil
        coil(icoil)%D = NDcoil

        coil(icoil)%Io =  Io
        coil(icoil)%Iw =  Iw
        coil(icoil)%Lw =  Lw

        if (.not. allocated(coil(icoil)%xc) ) then
           SALLOCATE( coil(icoil)%xc, (0:coil(icoil)%N), zero )
           SALLOCATE( coil(icoil)%xs, (0:coil(icoil)%N), zero )
           SALLOCATE( coil(icoil)%yc, (0:coil(icoil)%N), zero )
           SALLOCATE( coil(icoil)%ys, (0:coil(icoil)%N), zero )
           SALLOCATE( coil(icoil)%zc, (0:coil(icoil)%N), zero )
           SALLOCATE( coil(icoil)%zs, (0:coil(icoil)%N), zero ) 
        endif
        RlBCAST( coil(icoil)%xc(0:coil(icoil)%N) , 1+coil(icoil)%N ,  0 )
        RlBCAST( coil(icoil)%xs(0:coil(icoil)%N) , 1+coil(icoil)%N ,  0 )
        RlBCAST( coil(icoil)%yc(0:coil(icoil)%N) , 1+coil(icoil)%N ,  0 )
        RlBCAST( coil(icoil)%ys(0:coil(icoil)%N) , 1+coil(icoil)%N ,  0 )
        RlBCAST( coil(icoil)%zc(0:coil(icoil)%N) , 1+coil(icoil)%N ,  0 )
        RlBCAST( coil(icoil)%zs(0:coil(icoil)%N) , 1+coil(icoil)%N ,  0 )

        if(coil(icoil)%Ic == 0) Nfixcur = Nfixcur + 1
        if(coil(icoil)%Lc == 0) Nfixgeo = Nfixgeo + 1

        SALLOCATE( coil(icoil)%xx, (0:coil(icoil)%D), zero )
        SALLOCATE( coil(icoil)%yy, (0:coil(icoil)%D), zero )
        SALLOCATE( coil(icoil)%zz, (0:coil(icoil)%D), zero )
        SALLOCATE( coil(icoil)%xt, (0:coil(icoil)%D), zero )
        SALLOCATE( coil(icoil)%yt, (0:coil(icoil)%D), zero )
        SALLOCATE( coil(icoil)%zt, (0:coil(icoil)%D), zero )
        SALLOCATE( coil(icoil)%xa, (0:coil(icoil)%D), zero )
        SALLOCATE( coil(icoil)%ya, (0:coil(icoil)%D), zero )
        SALLOCATE( coil(icoil)%za, (0:coil(icoil)%D), zero )

     enddo

!!$
!!$   if( myid.eq.0 ) then
!!$
!!$    Ncoils = 0
!!$    do ii = 1, 999
!!$     write(suffix,'(i3.3)') ii
!!$     inquire( file=".fo.coil."//suffix, exist=exist )
!!$     if( exist ) Ncoils = Ncoils + 1
!!$    enddo
!!$
!!$    write(ounit,'("rdcoils : " 10x " : identified "i3" .fo.coil.xxx files ;")') Ncoils
!!$
!!$   endif
!!$
!!$   IlBCAST( Ncoils        ,    1,  0 )
!!$
!!$   allocate( coil(1:Ncoils) )
!!$
!!$   icoil = 0 
!!$
!!$   do ii = 1, 999
!!$
!!$    write(suffix,'(i3.3)') ii
!!$    inquire( file=".fo.coil."//suffix, exist=exist )
!!$
!!$    if( exist ) then
!!$
!!$     icoil = icoil + 1
!!$
!!$     SALLOCATE( coil(icoil)%xc, (0:NFcoil), zero )
!!$     SALLOCATE( coil(icoil)%xs, (0:NFcoil), zero )
!!$     SALLOCATE( coil(icoil)%yc, (0:NFcoil), zero )
!!$     SALLOCATE( coil(icoil)%ys, (0:NFcoil), zero )
!!$     SALLOCATE( coil(icoil)%zc, (0:NFcoil), zero )
!!$     SALLOCATE( coil(icoil)%zs, (0:NFcoil), zero )
!!$
!!$     if( myid.eq.0 ) then
!!$      open( lunit, file=".fo.coil."//suffix, status="old" )
!!$      read( lunit,*) coil(icoil)%N, coil(icoil)%D
!!$      read( lunit,*) coil(icoil)%I, coil(icoil)%Ic, coil(icoil)%Io, coil(icoil)%Iw
!!$      read( lunit,*) coil(icoil)%L, coil(icoil)%Lc, coil(icoil)%Lo, coil(icoil)%Lw
!!$      read( lunit,*) coil(icoil)%xc(0:min(coil(icoil)%N,NFcoil))
!!$      read( lunit,*) coil(icoil)%xs(0:min(coil(icoil)%N,NFcoil))
!!$      read( lunit,*) coil(icoil)%yc(0:min(coil(icoil)%N,NFcoil))
!!$      read( lunit,*) coil(icoil)%ys(0:min(coil(icoil)%N,NFcoil))
!!$      read( lunit,*) coil(icoil)%zc(0:min(coil(icoil)%N,NFcoil))
!!$      read( lunit,*) coil(icoil)%zs(0:min(coil(icoil)%N,NFcoil))
!!$      close( lunit )
!!$     endif ! end of if( myid.eq.0 ) ; 14 Apr 16;
!!$
!!$     IlBCAST( coil(icoil)%N            , 1        ,  0 )
!!$
!!$     IlBCAST( coil(icoil)%D            , 1        ,  0 )
!!$
!!$     RlBCAST( coil(icoil)%I            , 1        ,  0 )
!!$     IlBCAST( coil(icoil)%Ic           , 1        ,  0 )
!!$     RlBCAST( coil(icoil)%Io           , 1        ,  0 )
!!$     RlBCAST( coil(icoil)%Iw           , 1        ,  0 )
!!$
!!$     FATAL( rdcoils, coil(icoil)%Ic.lt.0 .or. coil(icoil)%Ic.gt.1, illegal )
!!$     FATAL( rdcoils, coil(icoil)%Iw.lt.zero                      , illegal )
!!$
!!$     RlBCAST( coil(icoil)%L            , 1        ,  0 )
!!$     IlBCAST( coil(icoil)%Lc           , 1        ,  0 )
!!$     RlBCAST( coil(icoil)%Lo           , 1        ,  0 )
!!$     RlBCAST( coil(icoil)%Lw           , 1        ,  0 )
!!$
!!$     FATAL( rdcoils, coil(icoil)%Lc.lt.0 .or. coil(icoil)%Lc.gt.2, illegal )
!!$     FATAL( rdcoils, coil(icoil)%Lo.lt.zero                      , illegal )
!!$     FATAL( rdcoils, coil(icoil)%Lw.lt.zero                      , illegal )
!!$
!!$     RlBCAST( coil(icoil)%xc(0:NFcoil) , 1+NFcoil ,  0 )
!!$     RlBCAST( coil(icoil)%xs(0:NFcoil) , 1+NFcoil ,  0 )
!!$     RlBCAST( coil(icoil)%yc(0:NFcoil) , 1+NFcoil ,  0 )
!!$     RlBCAST( coil(icoil)%ys(0:NFcoil) , 1+NFcoil ,  0 )
!!$     RlBCAST( coil(icoil)%zc(0:NFcoil) , 1+NFcoil ,  0 )
!!$     RlBCAST( coil(icoil)%zs(0:NFcoil) , 1+NFcoil ,  0 )
!!$
!!$     coil(icoil)%N = NFcoil ! all coils have same Fourier  resolution; should allow each coil to have different Fourier resolution; 14 Apr 16;
!!$     coil(icoil)%D = NDcoil ! all coils have same discrete resolution; should allow each coil to have different Fourier resolution; 14 Apr 16;
!!$
!!$     !   if( myid.eq.0 ) write(ounit,1000) icoil, coil(icoil)%I, coil(icoil)%Ic, coil(icoil)%Io, coil(icoil)%Iw
!!$
!!$     SALLOCATE( coil(icoil)%xx, (0:coil(icoil)%D), zero )
!!$     SALLOCATE( coil(icoil)%yy, (0:coil(icoil)%D), zero )
!!$     SALLOCATE( coil(icoil)%zz, (0:coil(icoil)%D), zero )
!!$     SALLOCATE( coil(icoil)%xt, (0:coil(icoil)%D), zero )
!!$     SALLOCATE( coil(icoil)%yt, (0:coil(icoil)%D), zero )
!!$     SALLOCATE( coil(icoil)%zt, (0:coil(icoil)%D), zero )
!!$     SALLOCATE( coil(icoil)%xa, (0:coil(icoil)%D), zero )
!!$     SALLOCATE( coil(icoil)%ya, (0:coil(icoil)%D), zero )
!!$     SALLOCATE( coil(icoil)%za, (0:coil(icoil)%D), zero )
!!$
!!$     write(coil(icoil)%name,'(i3.3"th-coil")') icoil
!!$
!!$     if(coil(icoil)%Ic .eq. 0) Nfixcur = Nfixcur + 1
!!$     if(coil(icoil)%Lc .eq. 0) Nfixgeo = Nfixgeo + 1
!!$
!!$    endif ! end of if( exist ) ; 14 Apr 16;
!!$
!!$   enddo ! matches do ii; 14 Apr 16;

     if(myid .eq. 0) write(ounit,'("rdcoils : " 10x " : "i3" fixed currents ; "i3" fixed geometries.")') Nfixcur, Nfixgeo

     !1000 format("rdcoils : " 10x " : "i3") I ="es13.5" ; Ic ="i2" ; Io ="es13.5" ; Iw ="es12.5" ;")

     !FATAL( rdcoils, icoil.ne.Ncoils, counting error )

     !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  case( 1: ) ! Linitialize; 14 Apr 16;

     !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

     !latex \subsection{initializing coils}
     !latex If suitable coils are not pre-determined, then various options for automatically constructing a suitable coil arrangement are available:
     !latex \bi
     !latex \item[1.] \inputvar{Linitialize} $= 1$:
     !latex \ei

     !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

     !if( myid.eq.0 ) write(ounit,1010) Linitialize, Ic, Io, Iw, Lc, Lo, Lw

1010 format("rdcoils : " 10x " : Linitialize ="i4" ; ":"Ic ="i2" ; Io ="es13.5" ; Iw ="es12.5" ; Lc ="i2" ; Lo ="es12.5" ; Lw ="es12.5" ;")

     Ncoils = Linitialize

     allocate( coil(1:Ncoils) )

     do icoil = 1, Ncoils

        coil(icoil)%N  =  NFcoil
        coil(icoil)%D  =  NDcoil

        !    Io = 8.1E7/Ncoils
        coil(icoil)%I  =  Io
        coil(icoil)%Ic =  Ic
        coil(icoil)%Io =  Io
        coil(icoil)%Iw =  Iw

        coil(icoil)%L  =  Lo ! irrelevant until re-computed; 14 Apr 16;
        coil(icoil)%Lc =  Lc
        coil(icoil)%Lo =  Lo
        coil(icoil)%Lw =  Lw

        SALLOCATE( coil(icoil)%xc, (0:NFcoil), zero )
        SALLOCATE( coil(icoil)%xs, (0:NFcoil), zero )
        SALLOCATE( coil(icoil)%yc, (0:NFcoil), zero )
        SALLOCATE( coil(icoil)%ys, (0:NFcoil), zero )
        SALLOCATE( coil(icoil)%zc, (0:NFcoil), zero )
        SALLOCATE( coil(icoil)%zs, (0:NFcoil), zero )

        zeta = (icoil-1+half) * pi2 / Ncoils ! +half for different initialize option; 2017/04/17 

        call surfcoord( zero, zeta, r1, z1)
        call surfcoord(   pi, zeta, r2, z2)

        Rmaj = half * (r1 + r2)
        z0   = half * (z1 + z2)

        ! shudson's representation;

!!$    coil(icoil)%xc(0:1) = (/ Rmaj * cos(zeta), rmin * cos(zeta) /)
!!$    coil(icoil)%xs(0:1) = (/ 0.0             , 0.0              /)
!!$    coil(icoil)%yc(0:1) = (/ Rmaj * sin(zeta), rmin * sin(zeta) /)
!!$    coil(icoil)%ys(0:1) = (/ 0.0             , 0.0              /)
!!$    coil(icoil)%zc(0:1) = (/ 0.0             , 0.0              /)
!!$    coil(icoil)%zs(0:1) = (/ 0.0             ,-1.0              /)

        ! czhu's representation;  07/09/2016

        coil(icoil)%xc(0:1) = (/ Rmaj * cos(zeta), rmin * cos(zeta) /)
        coil(icoil)%xs(0:1) = (/ 0.0             , 0.0              /)
        coil(icoil)%yc(0:1) = (/ Rmaj * sin(zeta), rmin * sin(zeta) /)
        coil(icoil)%ys(0:1) = (/ 0.0             , 0.0              /)
        coil(icoil)%zc(0:1) = (/ z0              , 0.0              /)
        coil(icoil)%zs(0:1) = (/ 0.0             , rmin             /)
!!$    ! angle difference in coil representation;
!!$    coil(icoil)%xc(0:1) = (/ Rmaj * cos(zeta), sqrt(2.0)/2 * rmin * cos(zeta) /)
!!$    coil(icoil)%xs(0:1) = (/ 0.0             ,-sqrt(2.0)/2 * rmin * cos(zeta) /)
!!$    coil(icoil)%yc(0:1) = (/ Rmaj * sin(zeta), sqrt(2.0)/2 * rmin * sin(zeta) /)
!!$    coil(icoil)%ys(0:1) = (/ 0.0             ,-sqrt(2.0)/2 * rmin * sin(zeta) /)
!!$    coil(icoil)%zc(0:1) = (/ z0              , sqrt(2.0)/2 * rmin             /)
!!$    coil(icoil)%zs(0:1) = (/ 0.0             , sqrt(2.0)/2 * rmin             /)

        SALLOCATE( coil(icoil)%xx, (0:coil(icoil)%D), zero )
        SALLOCATE( coil(icoil)%yy, (0:coil(icoil)%D), zero )
        SALLOCATE( coil(icoil)%zz, (0:coil(icoil)%D), zero )
        SALLOCATE( coil(icoil)%xt, (0:coil(icoil)%D), zero )
        SALLOCATE( coil(icoil)%yt, (0:coil(icoil)%D), zero )
        SALLOCATE( coil(icoil)%zt, (0:coil(icoil)%D), zero )
        SALLOCATE( coil(icoil)%xa, (0:coil(icoil)%D), zero )
        SALLOCATE( coil(icoil)%ya, (0:coil(icoil)%D), zero )
        SALLOCATE( coil(icoil)%za, (0:coil(icoil)%D), zero )

        write(coil(icoil)%name,'(i3.3"th-coil")') icoil

        if(coil(icoil)%Ic .eq. 0) Nfixcur = Nfixcur + 1
        if(coil(icoil)%Lc .eq. 0) Nfixgeo = Nfixgeo + 1

     enddo ! end of do icoil; 14 Apr 16;

  end select

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!  
  !if( myid.eq.0 ) then
  ! write(ounit,'("rdcoils : " 10x " : I ="999(es10.2","))') ( coil(icoil)%I, icoil = 1, Ncoils )
  ! write(ounit,'("rdcoils : " 10x " : L ="999(es10.2","))') ( coil(icoil)%L, icoil = 1, Ncoils )
  !endif
  totalcurrent = zero
  do icoil = 1, Ncoils
     totalcurrent = totalcurrent + coil(icoil)%I
     !write( ounit,'("rdcoils : ",i5," : I =",es23.15," ; \sum I ="es23.15" ;")') icoil, coil(icoil)%I, totalcurrent
  enddo

  if( myid.eq.0 ) write( ounit,'("rdcoils : "10x" : total current G ="es23.15" ; 2 . pi2 . G = "es23.15" ;")') totalcurrent, totalcurrent * pi2 * two

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  SALLOCATE( cmt, (0:NDcoil,0:NFcoil), zero )
  SALLOCATE( smt, (0:NDcoil,0:NFcoil), zero )

  do itime = 0, NDcoil ; tt = itime * pi2 / NDcoil
     do mm = 0, NFcoil
        cmt(itime,mm) = cos( mm * tt )
        smt(itime,mm) = sin( mm * tt )
     enddo
  enddo

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  Ndof = ( 1 + 3 + 6*NFcoil ) * Ncoils - Nfixcur - (3+6*NFcoil)*Nfixgeo

!  if ( Loptimize .eq. 3 ) Ndof = ( 1 + 3 + 8*NFcoil ) * Ncoils - Nfixcur - (3+8*NFcoil)*Nfixgeo ! for newton method only;

  call discretecoil

  Tdof = ( 1 + 3 + 6*NFcoil ) * Ncoils !total dimension
  SALLOCATE( coilspace, (0:Ntauout,1:Tdof), zero ) ! coils' current and fourier harmonics data;
!  if (Loptimize .ne. 3) then
     SALLOCATE( deriv    , (1:Ndof,0:5), zero )
!  endif

  do icoil = 1, Ncoils
     SALLOCATE( coil(icoil)%lmdc, (0:NFcoil), one )
     SALLOCATE( coil(icoil)%lmds, (0:NFcoil), one )   ! initialized as one; 07/26/2016
  enddo

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  SALLOCATE(   norm, (1:Ndof)        , zero )  !allocate normalization vector; 20170510;
  if (Lnormalize .ge. 2) then

     totalcurrent = zero
     do icoil = 1, Ncoils
        totalcurrent = totalcurrent + abs(coil(icoil)%I)
     enddo
     Inorm = totalcurrent/Ncoils !mean of abs values;

     r1 = sqrt( surf(1)%xx(      0,0)**2 + surf(1)%yy(      0,0)**2 ) ! R at (0 ,0)
     r2 = sqrt( surf(1)%xx(Nteta/2,0)**2 + surf(1)%yy(Nteta/2,0)**2 ) ! R at (pi,0)
     Gnorm = half * (r1 + r2) ! something like the major radius; need to be changded;

     idof = 0
     do icoil = 1, Ncoils

        if(coil(icoil)%Ic.ne. 0) then 
           idof = idof + 1 ; norm(idof) = Inorm
        endif

        if(coil(icoil)%Lc.ne. 0) then  
           idof = idof + 1 ; norm(idof) = Gnorm
           idof = idof + 1 ; norm(idof) = Gnorm
           idof = idof + 1 ; norm(idof) = Gnorm
           do mm = 1, NFcoil 
              idof = idof + 1 ; norm(idof) = Gnorm
              idof = idof + 1 ; norm(idof) = Gnorm
              idof = idof + 1 ; norm(idof) = Gnorm
              idof = idof + 1 ; norm(idof) = Gnorm
              idof = idof + 1 ; norm(idof) = Gnorm
              idof = idof + 1 ; norm(idof) = Gnorm
           enddo
        endif

     enddo
     FATAL( rdcoils , idof .ne. Ndof, counting error in packing )

     if( myid.eq.0 ) write( ounit,'("rdcoils : "10x" : currents are normalized by " ES23.15 " ; and geometry by "ES23.15 " .")') Inorm, Gnorm

  else

     Inorm = one; Gnorm = one
     norm = one

  endif

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  return

end subroutine rdcoils

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!


SUBROUTINE readcoils(filename, maxnseg)
  use kmodule, only : zero, coilsX, coilsY, coilsZ, coilsI, Nseg, Ncoils, ounit, myid
  implicit none
  include "mpif.h"

  INTEGER, parameter         :: mcoil = 256, mseg = 1024 ! Largest coils and segments number
  INTEGER                    :: cunit, istat, astat, lstat, ierr, maxnseg, seg(1:mseg), icoil
  REAL, dimension(mseg,mcoil):: x, y, z, I
  CHARACTER*40               :: filename
  CHARACTER*200              :: line

  cunit = 99; I = 1.0; Ncoils= 1; maxnseg = 0; seg = 0;

  open(cunit,FILE=filename,STATUS='old',IOSTAT=istat)
  if ( istat .ne. 0 ) stop "Reading coils error!"

! read coils and segments data
  read(cunit,*)
  read(cunit,*)
  read(cunit,*)

  do
     read(cunit,'(a)', IOSTAT = istat) line
     if(istat .ne. 0 .or. line(1:3) .eq. 'end') exit !detect EOF or end

     seg(Ncoils) = seg(Ncoils) + 1
     read(line,*, IOSTAT = lstat) x(seg(Ncoils), Ncoils), y(seg(Ncoils), Ncoils), z(seg(Ncoils), Ncoils), &
                                  I(seg(Ncoils), Ncoils)
     if ( I(seg(Ncoils), Ncoils) .eq. 0.0 ) then
        seg(Ncoils) = seg(Ncoils) - 1  !remove the duplicated last point
        !if ( seg(Ncoils) .ge. maxnseg ) maxnseg = seg(Ncoils)
        Ncoils = Ncoils + 1
     endif    
  enddo

  close(cunit)

  Ncoils = Ncoils - 1
  maxnseg = maxval(seg)
  FATAL( readcoils, Ncoils  .le. 0 .or. Ncoils  .gt. mcoil, illegal )
  FATAL( readcoils, maxnseg .le. 0 .or. maxnseg .gt. mseg , illegal )

#ifdef DEBUG
  write(ounit,'("rdcoils : " 10X " : Finding " I4 " coils and maximum segments number is " I6)') Ncoils, maxnseg
  do icoil = 1, Ncoils
     write(ounit, *) icoil, seg(icoil)
  enddo
#endif

  SALLOCATE( coilsX, (1:maxnseg, 1:Ncoils), zero )
  SALLOCATE( coilsY, (1:maxnseg, 1:Ncoils), zero )
  SALLOCATE( coilsZ, (1:maxnseg, 1:Ncoils), zero )
  SALLOCATE( coilsI, (1:Ncoils)           , zero )
  SALLOCATE(   Nseg, (1:Ncoils)           ,    0 )

  coilsX( 1:maxnseg, 1:Ncoils) = x( 1:maxnseg, 1:Ncoils)
  coilsY( 1:maxnseg, 1:Ncoils) = y( 1:maxnseg, 1:Ncoils)
  coilsZ( 1:maxnseg, 1:Ncoils) = z( 1:maxnseg, 1:Ncoils)
  coilsI(            1:Ncoils) = I( 1        , 1:Ncoils)
    Nseg(            1:Ncoils) = seg(          1:Ncoils)
  
  return

end SUBROUTINE READCOILS

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!!$
!!$SUBROUTINE FFT( X, XFC, XFS, Nsegs, NFcoil)
!!$  use kmodule, only: ounit
!!$  implicit none
!!$  include "mpif.h"
!!$
!!$  REAL    :: X(1:Nsegs), XFC(0:NFcoil), XFS(0:NFcoil), Xtmp(0:Nsegs-1)
!!$  INTEGER :: Nsegs, NFcoil, ifail, ii, m, n
!!$  REAL, allocatable:: A(:), B(:)
!!$  LOGICAL :: IsOver
!!$
!!$  allocate(A(0:Nsegs-1))
!!$  allocate(B(0:Nsegs-1))
!!$  
!!$
!!$  if( IsOver(Nsegs) ) then
!!$#ifdef DEBUG
!!$   write(ounit,'(A, i3, A)') "The segments number ", nsegs, " has a prime number larger than 19. The last point was removed."
!!$#endif
!!$   Nsegs = Nsegs -1
!!$  endif
!!$ 
!!$  do n = 0, Nsegs-1
!!$   Xtmp(n) = X(n+1)
!!$  enddo
!!$
!!$  ifail = 0
!!$  call C06EAF( Xtmp(0:Nsegs-1), Nsegs, ifail ) 
!!$
!!$  ifail = 0
!!$  call C06GBF( Xtmp(0:Nsegs-1), Nsegs, ifail )
!!$
!!$  m = 1
!!$  ifail = 0
!!$  call c06gsf(m,Nsegs,Xtmp,A,B,ifail)
!!$
!!$  XFC(0) = A(0)
!!$  XFS(0) = B(0)
!!$
!!$  do ii = 1, NFcoil
!!$   XFC(ii) = 2.0 * A(ii)
!!$   XFS(ii) = 2.0 * B(ii)
!!$  enddo
!!$
!!$  XFC = XFC / sqrt(real(Nsegs))
!!$  XFS = XFS / sqrt(real(Nsegs))
!!$
!!$  deallocate(A, B)
!!$
!!$  return
!!$
!!$end SUBROUTINE FFT


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!find if a number's largest prime number is bigger than 19
LOGICAL function  IsOver( num )
  IMPLICIT  NONE
  include "mpif.h"

  INTEGER, intent(in):: num
  INTEGER :: Input
  INTEGER :: Divisor, Count, maxprime

  if(num .eq. 0) stop "Number of segments can't be zero!"
  input = num
  maxprime = 19
  IsOver = .false.
  Count = 0
  DO                         ! here, we try to remove all factors of 2
   IF (MOD(Input,2) /= 0 .OR. Input == 1)  EXIT
   Count = Count + 1       ! increase count
   Input = Input / 2       ! remove this factor from Input
  END DO

  Divisor = 3                ! now we only worry about odd factors
  DO                         ! 3, 5, 7, .... will be tried
   IF (Divisor > Input) EXIT    ! if a factor is too large, exit and done
   DO                      ! try this factor repeatedly
    IF (MOD(Input,Divisor) /= 0 .OR. Input == 1)  EXIT
    Count = Count + 1
    Input = Input / Divisor   ! remove this factor from Input
   END DO
   Divisor = Divisor + 2   ! move to next odd number
   if ( Divisor .gt. 19 ) IsOver = .true.
  END DO

  return
END function IsOver

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

SUBROUTINE Fourier( X, XFC, XFS, Nsegs, NFcoil)
  use kmodule, only: ounit, zero, pi2, half, myid
  implicit none
  include "mpif.h"

  REAL    :: X(1:Nsegs), XFC(0:NFcoil), XFS(0:NFcoil)
  INTEGER :: Nsegs, NFcoil, ifou, iseg, funit
  REAL, allocatable:: A(:), B(:)
  LOGICAL :: IsOver

  allocate(A(0:Nsegs-1))
  allocate(B(0:Nsegs-1))

  A = zero; B = zero

  do ifou = 0, Nsegs-1

     do iseg = 1, Nsegs
        A(ifou) = A(ifou) + X(iseg)*cos(ifou*pi2*(iseg-1)/Nsegs)
        B(ifou) = B(ifou) + X(iseg)*sin(ifou*pi2*(iseg-1)/Nsegs)
     enddo

  enddo

  A = 2.0/Nsegs * A; A(0) = half*A(0)
  B = 2.0/Nsegs * B

  XFC(0:NFcoil) = A(0:NFcoil)
  XFS(0:NFcoil) = B(0:NFcoil)
  
  deallocate(A, B)

  return

end SUBROUTINE Fourier

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine discretecoil
  use kmodule, only: zero, &
                     coil, cmt, smt, NFcoil, NDcoil, Ncoils, &
                     ncpu, myid
  implicit none
  include "mpif.h"


  INTEGER           :: astat, ierr
  INTEGER           :: icoil,llmodnp, mm

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

  !coil discresization; xx, xt, xa are 0, 1st and 2nd derivatives;
  do icoil = 1, Ncoils

  if( myid.ne.modulo(icoil-1,ncpu) ) cycle ! parallelization loop;

  mm = 0
  ;coil(icoil)%xx(0:NDcoil) =                            cmt(0:NDcoil,mm) * coil(icoil)%xc(mm)
  ;coil(icoil)%yy(0:NDcoil) =                            cmt(0:NDcoil,mm) * coil(icoil)%yc(mm)
  ;coil(icoil)%zz(0:NDcoil) =                            cmt(0:NDcoil,mm) * coil(icoil)%zc(mm)
  ;coil(icoil)%xt(0:NDcoil) =                                               zero
  ;coil(icoil)%yt(0:NDcoil) =                                               zero
  ;coil(icoil)%zt(0:NDcoil) =                                               zero
  ;coil(icoil)%xa(0:NDcoil) =                                               zero
  ;coil(icoil)%ya(0:NDcoil) =                                               zero
  ;coil(icoil)%za(0:NDcoil) =                                               zero
  do mm = 1, NFcoil

     coil(icoil)%xx(0:NDcoil) = coil(icoil)%xx(0:NDcoil) + (   cmt(0:NDcoil,mm) * coil(icoil)%xc(mm) + smt(0:NDcoil,mm) * coil(icoil)%xs(mm) )
     coil(icoil)%yy(0:NDcoil) = coil(icoil)%yy(0:NDcoil) + (   cmt(0:NDcoil,mm) * coil(icoil)%yc(mm) + smt(0:NDcoil,mm) * coil(icoil)%ys(mm) )
     coil(icoil)%zz(0:NDcoil) = coil(icoil)%zz(0:NDcoil) + (   cmt(0:NDcoil,mm) * coil(icoil)%zc(mm) + smt(0:NDcoil,mm) * coil(icoil)%zs(mm) )
     coil(icoil)%xt(0:NDcoil) = coil(icoil)%xt(0:NDcoil) + ( - smt(0:NDcoil,mm) * coil(icoil)%xc(mm) + cmt(0:NDcoil,mm) * coil(icoil)%xs(mm) ) * mm
     coil(icoil)%yt(0:NDcoil) = coil(icoil)%yt(0:NDcoil) + ( - smt(0:NDcoil,mm) * coil(icoil)%yc(mm) + cmt(0:NDcoil,mm) * coil(icoil)%ys(mm) ) * mm
     coil(icoil)%zt(0:NDcoil) = coil(icoil)%zt(0:NDcoil) + ( - smt(0:NDcoil,mm) * coil(icoil)%zc(mm) + cmt(0:NDcoil,mm) * coil(icoil)%zs(mm) ) * mm
     coil(icoil)%xa(0:NDcoil) = coil(icoil)%xa(0:NDcoil) + ( - cmt(0:NDcoil,mm) * coil(icoil)%xc(mm) - smt(0:NDcoil,mm) * coil(icoil)%xs(mm) ) * mm * mm
     coil(icoil)%ya(0:NDcoil) = coil(icoil)%ya(0:NDcoil) + ( - cmt(0:NDcoil,mm) * coil(icoil)%yc(mm) - smt(0:NDcoil,mm) * coil(icoil)%ys(mm) ) * mm * mm
     coil(icoil)%za(0:NDcoil) = coil(icoil)%za(0:NDcoil) + ( - cmt(0:NDcoil,mm) * coil(icoil)%zc(mm) - smt(0:NDcoil,mm) * coil(icoil)%zs(mm) ) * mm * mm

  enddo ! end of do ii; 03 May 16;
  enddo ! end of do icoil

  do icoil = 1, Ncoils ; llmodnp = modulo(icoil-1,ncpu)
   RlBCAST( coil(icoil)%xx(0:NDcoil), NDcoil+1, llmodnp )
   RlBCAST( coil(icoil)%yy(0:NDcoil), NDcoil+1, llmodnp )
   RlBCAST( coil(icoil)%zz(0:NDcoil), NDcoil+1, llmodnp )
   RlBCAST( coil(icoil)%xt(0:NDcoil), NDcoil+1, llmodnp )
   RlBCAST( coil(icoil)%yt(0:NDcoil), NDcoil+1, llmodnp )
   RlBCAST( coil(icoil)%zt(0:NDcoil), NDcoil+1, llmodnp )
   RlBCAST( coil(icoil)%xa(0:NDcoil), NDcoil+1, llmodnp )
   RlBCAST( coil(icoil)%ya(0:NDcoil), NDcoil+1, llmodnp )
   RlBCAST( coil(icoil)%za(0:NDcoil), NDcoil+1, llmodnp )
  enddo

return
end subroutine discretecoil

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine pack( xdof )

  use kmodule, only : NFcoil, Ncoils, coil, Ndof, myid, ounit, norm
  implicit none
  include "mpif.h"

  REAL    :: xdof(Ndof)

  INTEGER :: idof, icoil, mm, ierr

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  idof = 0
  do icoil = 1, Ncoils
     
     if(coil(icoil)%Ic.ne. 0) then 
        idof = idof + 1 ; xdof(idof) = coil(icoil)%I      / norm(idof)
     endif

     if(coil(icoil)%Lc.ne. 0) then  
        idof = idof + 1 ; xdof(idof) = coil(icoil)%xc( 0) / norm(idof)
        idof = idof + 1 ; xdof(idof) = coil(icoil)%yc( 0) / norm(idof)
        idof = idof + 1 ; xdof(idof) = coil(icoil)%zc( 0) / norm(idof)
      do mm = 1, NFcoil 
        idof = idof + 1 ; xdof(idof) = coil(icoil)%xc(mm) / norm(idof)
        idof = idof + 1 ; xdof(idof) = coil(icoil)%yc(mm) / norm(idof)
        idof = idof + 1 ; xdof(idof) = coil(icoil)%zc(mm) / norm(idof)
        idof = idof + 1 ; xdof(idof) = coil(icoil)%xs(mm) / norm(idof)
        idof = idof + 1 ; xdof(idof) = coil(icoil)%ys(mm) / norm(idof)
        idof = idof + 1 ; xdof(idof) = coil(icoil)%zs(mm) / norm(idof)
      enddo
     endif

  enddo

  FATAL( pack , idof .ne. Ndof, counting error in packing )

  return

end subroutine pack

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!  

subroutine unpack( xdof )

  use kmodule, only : NFcoil, Ncoils, coil, Ndof, myid, zero, norm
  implicit none
  include "mpif.h"

  REAL    :: xdof(Ndof)

  INTEGER :: idof, icoil, mm, ierr

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!  
  idof = 0
  
  do icoil = 1, Ncoils

     if(coil(icoil)%Ic.ne. 0) then
        idof = idof + 1 ; coil(icoil)%I      = xdof(idof) * norm(idof)
     endif

     if(coil(icoil)%Lc.ne. 0) then
        idof = idof + 1 ; coil(icoil)%xc( 0) = xdof(idof) * norm(idof)
                        ; coil(icoil)%xs( 0) = zero
        idof = idof + 1 ; coil(icoil)%yc( 0) = xdof(idof) * norm(idof)
                        ; coil(icoil)%ys( 0) = zero
        idof = idof + 1 ; coil(icoil)%zc( 0) = xdof(idof) * norm(idof)
                        ; coil(icoil)%zs( 0) = zero        
      do mm = 1, NFcoil
        idof = idof + 1 ; coil(icoil)%xc(mm) = xdof(idof) * norm(idof)
        idof = idof + 1 ; coil(icoil)%yc(mm) = xdof(idof) * norm(idof)
        idof = idof + 1 ; coil(icoil)%zc(mm) = xdof(idof) * norm(idof)
        idof = idof + 1 ; coil(icoil)%xs(mm) = xdof(idof) * norm(idof)
        idof = idof + 1 ; coil(icoil)%ys(mm) = xdof(idof) * norm(idof)
        idof = idof + 1 ; coil(icoil)%zs(mm) = xdof(idof) * norm(idof)
      enddo
     endif

  enddo

  FATAL( unpack , idof .ne. Ndof, counting error in unpacking )
  return

end subroutine unpack

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!  

SUBROUTINE DoFconvert(idof, icoil, inf)  !not applied to fixed DoFs; 01/09/2017
  use kmodule, only: shudson, Ncoils, NFcoil, ounit, myid
  implicit none
  include "mpif.h"

  INTEGER, intent(in ) :: idof
  INTEGER, intent(out) :: icoil, inf

  INTEGER              :: imod, ncdof, ierr

  FATAL( DoFconvert, .not. allocated(shudson), need allocate shudson array first)

  ncdof = 1 + 3 + 6 * NFcoil
  if( (idof .le. 0) .or. (idof .gt. Ncoils*ncdof) ) then
   call MPI_ABORT( MPI_COMM_WORLD, 1, ierr )
   stop "index out of range in DoFconvert!"
  endif
  
  icoil = (idof-1) / ncdof + 1
  imod  = modulo(idof,ncdof)
  if ( imod .eq. 0 ) then
   inf = shudson(ncdof)
  else
   inf = shudson(imod)
  endif

  return
end SUBROUTINE DoFconvert
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine ntpack( xdof )

  use kmodule, only : NFcoil, Ncoils, coil, Ndof, myid
  implicit none
  include "mpif.h"

  REAL    :: xdof(Ndof)

  INTEGER :: idof, icoil, mm, ierr

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  FATAL( ntpack , Ndof .ne. (4+8*NFcoil)*Ncoils, change Ndof first )
  
  idof = 0
  do icoil = 1, Ncoils
   ;                 ; idof = idof + 1 ; xdof(idof) = coil(icoil)%I
   ;  mm = 0         ; idof = idof + 1 ; xdof(idof) = coil(icoil)%xc(  mm)
   ;                 ; idof = idof + 1 ; xdof(idof) = coil(icoil)%yc(  mm)
   ;                !; idof = idof + 1 ; xdof(idof) = coil(icoil)%zc(  mm)   !comment out for fixed poloidal angel origin
   ;                 ; idof = idof + 1 ; xdof(idof) = coil(icoil)%lmdc(mm)
   do mm = 1, NFcoil ; idof = idof + 1 ; xdof(idof) = coil(icoil)%xc(  mm)
    ;                ; idof = idof + 1 ; xdof(idof) = coil(icoil)%yc(  mm)
    ;                ; idof = idof + 1 ; xdof(idof) = coil(icoil)%zc(  mm)
    ;                ; idof = idof + 1 ; xdof(idof) = coil(icoil)%lmdc(mm)
    ;                ; idof = idof + 1 ; xdof(idof) = coil(icoil)%xs(  mm)
    ;                ; idof = idof + 1 ; xdof(idof) = coil(icoil)%ys(  mm)
    ;                ; idof = idof + 1 ; xdof(idof) = coil(icoil)%zs(  mm)
    ;                ; idof = idof + 1 ; xdof(idof) = coil(icoil)%lmds(mm)
   enddo
  enddo

  FATAL( ntpack , idof .ne. Ndof, counting error in newton method packing )
  return

end subroutine ntpack

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!  

subroutine ntunpack( xdof )

  use kmodule, only : NFcoil, Ncoils, coil, Ndof, myid, zero
  implicit none
  include "mpif.h"

  REAL    :: xdof(Ndof)

  INTEGER :: idof, icoil, mm, ierr

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-! 

  FATAL( ntunpack , Ndof .ne. (4+8*NFcoil)*Ncoils, change Ndof first )

  idof = 0
  

  do icoil = 1, Ncoils
   ;                 ; idof = idof + 1 ; coil(icoil)%I        = xdof(idof)
   ;  mm = 0         ; idof = idof + 1 ; coil(icoil)%xc(  mm) = xdof(idof)
                                       ; coil(icoil)%xs(  mm) = zero      
   ;                 ; idof = idof + 1 ; coil(icoil)%yc(  mm) = xdof(idof)
                                       ; coil(icoil)%ys(  mm) = zero
   ;                !; idof = idof + 1 ; coil(icoil)%zc(  mm) = xdof(idof)
                                       ; coil(icoil)%zs(  mm) = zero
   ;                 ; idof = idof + 1 ; coil(icoil)%lmdc(mm) = xdof(idof)
                                       ; coil(icoil)%lmds(mm) = zero
   do mm = 1, NFcoil ; idof = idof + 1 ; coil(icoil)%xc(  mm) = xdof(idof)
    ;                ; idof = idof + 1 ; coil(icoil)%yc(  mm) = xdof(idof)
    ;                ; idof = idof + 1 ; coil(icoil)%zc(  mm) = xdof(idof)
    ;                ; idof = idof + 1 ; coil(icoil)%lmdc(mm) = xdof(idof)
    ;                ; idof = idof + 1 ; coil(icoil)%xs(  mm) = xdof(idof)
    ;                ; idof = idof + 1 ; coil(icoil)%ys(  mm) = xdof(idof)
    ;                ; idof = idof + 1 ; coil(icoil)%zs(  mm) = xdof(idof)
    ;                ; idof = idof + 1 ; coil(icoil)%lmds(mm) = xdof(idof)
   enddo
  enddo

  FATAL( ntunpack , idof .ne. Ndof, counting error in newton unpacking )
  return

end subroutine ntunpack

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!  

SUBROUTINE ntconvert(idof, icoil, inf)
  use kmodule, only: newton, Ncoils, NFcoil, ounit, myid
  implicit none
  include "mpif.h"

  INTEGER, intent(in ) :: idof
  INTEGER, intent(out) :: icoil, inf

  INTEGER              :: imod, ncdof, ierr

  FATAL( ntconvert, .not. allocated(newton), need allocate newton array first)

  ncdof = 1 + 3 + 8 * NFcoil
  if( (idof .le. 0) .or. (idof .gt. Ncoils*ncdof) ) then
   call MPI_ABORT( MPI_COMM_WORLD, 1, ierr )
   stop "index out of range in DoFconvert!"
  endif
  
  icoil = (idof-1) / ncdof + 1
  imod  = modulo(idof,ncdof)
  if ( imod .eq. 0 ) then
   inf = newton(ncdof)
  else
   inf = newton(imod)
  endif

  return

END SUBROUTINE ntconvert

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!  

! test shudson and newton array
!!$SUBROUTINE test
!!$  use kmodule, only: Ndof, NFcoil, Ncoils, Ndof, myid, ounit
!!$  implicit none
!!$  include "mpif.h"
!!$
!!$  INTEGER :: idof, icoil, inf, imod, ncdof, ierr
!!$  
!!$  ncdof = 1 + 3 + 8 * NFcoil
!!$
!!$  Ndof = ncdof * Ncoils
!!$  
!!$  do idof = 1, Ndof
!!$
!!$     call ntconvert(idof, icoil, inf)
!!$
!!$     if(myid .eq. 0) write(ounit,'("idof = "I6" : icoil = " I3 " ; inf = "I3)') idof, icoil, inf
!!$
!!$  enddo
!!$
!!$  return
!!$
!!$END SUBROUTINE test

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
