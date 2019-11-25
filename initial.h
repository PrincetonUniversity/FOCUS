
!title (initialize) ! Read input file, initialize global variables.

!latex \briefly{Reads input file, broadcasts, and allocates/initializes global variables.}

!latex \calledby{\link{focus}}
!l tex \calls{\link{}}

!latex \tableofcontents

!latex \subsection{overview}
!latex \bi
!latex \item[1.] 
!latex \ei

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine initial

  use kmodule

  implicit none

  include "mpif.h"

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  LOGICAL :: exist
  INTEGER :: ierr, astat, idof, mm, ncdof, imn, nn, Nfp, ifail
  REAL    :: X02AJF

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  machprec = X02AJF() ; sqrtmachprec = sqrt(machprec) ; vsmall = ten * machprec ; small = thousand * machprec

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( myid.eq.0 ) then ! only the master node reads the input; 25 Mar 15;
     call getarg(1,ext)
     write(ounit,'("initial : " 10x " : machprec ="es12.5" ; sqrtmachprec ="es12.5" ; ext = "a)') machprec, sqrtmachprec, trim(ext)
     inquire( file=trim(ext)//".fo", exist=exist )
  endif

  LlBCAST( exist, 1, 0 )
  FATAL( initial, .not.exist, input file ext.fo not provided )

  if( myid.eq.0 ) then
     open(lunit, file=trim(ext)//".fo", status="unknown")
     read(lunit, focusin)
     close(lunit)
  endif ! end of if( myid.eq.0 ) ; 25 Mar 15;

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  ClBCAST( ext           ,  100,  0 )
  IlBCAST( Idisplay      ,    1,  0 )
  IlBCAST( Isymmetric    ,    1,  0 )
  IlBCAST( Itopology     ,    1,  0 )
  RlBCAST( knotsurf      ,    1,  0 )
  RlBCAST( ellipticity   ,    1,  0 )
  IlBCAST( nrotate       ,    1,  0 )
  RlBCAST( zetaoff       ,    1,  0 )
  IlBCAST( NFcoil        ,    1,  0 )
  IlBCAST( NDcoil        ,    1,  0 )
  IlBCAST( Linitialize   ,    1,  0 )
  RlBCAST( Rmaj          ,    1,  0 )
  RlBCAST( rmin          ,    1,  0 )
  IlBCAST( Ic            ,    1,  0 )
  RlBCAST( Io            ,    1,  0 )
  RlBCAST( Iw            ,    1,  0 )
  IlBCAST( Lc            ,    1,  0 )
  RlBCAST( Lo            ,    1,  0 )
  RlBCAST( Lw            ,    1,  0 )
  IlBCAST( Loptimize     ,    1,  0 )
  IlBCAST( Lnormalize    ,    1,  0 )
  RlBCAST( weight_bnorm  ,    1,  0 )
  RlBCAST( weight_tflux  ,    1,  0 )
  RlBCAST( target_tflux  ,    1,  0 )
  RlBCAST( weight_ttlen  ,    1,  0 )
  RlBCAST( weight_eqarc  ,    1,  0 )
  RlBCAST( weight_ccsep  ,    1,  0 )
  RlBCAST( weight_qasym  ,    1,  0 )
  RlBCAST( weight_resbn  ,    1,  0 )
  RlBCAST( tauend        ,    1,  0 )
  RlBCAST( tautol        ,    1,  0 )
  IlBCAST( Ntauout       ,    1,  0 )
  IlBCAST( Savfreq       ,    1,  0 )
  IlBCAST( Nteta         ,    1,  0 )
  IlBCAST( Nzeta         ,    1,  0 )
  RlBCAST( absacc        ,    1,  0 )
  RlBCAST( absreq        ,    1,  0 )
  RlBCAST( relreq        ,    1,  0 )
  RlBCAST( xtol          ,    1,  0 )
  RlBCAST( eta           ,    1,  0 )
  RlBCAST( stepmx        ,    1,  0 )
  IlBCAST( Mpol          ,    1,  0 )
  IlBCAST( Ntor          ,    1,  0 )
  IlBCAST( Lpoincare     ,    1,  0 )
  RlBCAST( odetol        ,    1,  0 )
  IlBCAST( Ppts          ,    1,  0 )
  IlBCAST( Ptrj          ,    1,  0 )
  IlBCAST( iphi          ,    1,  0 )
  RlBCAST( bstol         ,    1,  0 )
  IlBCAST( bsnlimit      ,    1,  0 )
  
  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  !bsfield%tol = bstol ; bsfield%N = bsnlimit ! bsfield is global; passed through to Oculus:bs00aa; 30 Oct 15;

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( myid.eq.0 ) then

     if( Itopology   .eq. 0 ) write(ounit,0100) Idisplay, Isymmetric, Itopology, NFcoil, NDcoil
     if( Itopology   .eq. 1 ) write(ounit,0101) Idisplay, Isymmetric, Itopology, knotsurf, ellipticity, NFcoil, NDcoil
     if( Itopology   .eq. 2 ) write(ounit,0100) Idisplay, Isymmetric, Itopology, NFcoil, NDcoil

     if( Linitialize .le. 0 ) write(ounit,0102) Linitialize
     if( Linitialize .gt. 0 ) write(ounit,0103) Linitialize, Rmaj, rmin, Ic, Io, Iw, Lc, Lo, Lw

     if( Loptimize   .eq. 0 ) write(ounit,0104) Loptimize
     if( Loptimize   .ne. 0 ) then 
         write(ounit,0105) Loptimize, tauend, tautol, Ntauout, Nteta
         write(ounit,0106) Nzeta, xtol, eta, stepmx
         write(ounit,0111) Lnormalize
         write(ounit,0107) weight_bnorm, weight_tflux, target_tflux
#ifdef NORM
         write(ounit, '(36X " ; Bn normalized to Bmod.")')
#endif
         write(ounit,0108) weight_ttlen, weight_eqarc, weight_ccsep
      endif
      
     if( Lpoincare   .eq. 0 ) write(ounit,0109) Lpoincare
     if( Lpoincare   .ne. 0 ) write(ounit,0110) Lpoincare, odetol, Ppts, Ptrj, iphi, bstol, bsnlimit

  endif

0100 format("initial : " 10x " : Idisplay ="i2 " ; Isymmetric ="i2" ; Itopology ="i2" ; NFcoil="i3" ; NDcoil ="i4" ;")
0101 format("initial : " 10x " : Idisplay ="i2 " ; Isymmetric ="i2" ; Itopology ="i2" ; knotsurf ="f7.3" ; ellipticity="f7.3" ; NFcoil="i3" ; NDcoil ="i4" ;")

0102 format("initial : " 10x " : Linitialize ="i4" ;")
0103 format("initial : " 10x " : Linitialize ="i4" ; Rmaj ="f5.2" ; rmin ="f5.2" ; Ic ="i2" ; Io ="es13.5" ; Iw ="es12.5" ; Lc ="i2" ; Lo ="es12.5" ; Lw ="es12.5" ;")

0104 format("initial : " 10x " : Loptimize ="i2" ;")
0105 format("initial : " 10x " : Loptimize ="i2" ; tauend ="es12.5" ; tautol ="es9.2" ; Ntauout ="i6" ; Nteta ="i5)
0106 format( 36X                               " ; Nzeta ="i5" ; xtol ="es12.5" ; eta =" es12.5" ; stepmx ="es12.5)
0111 format( 36X                               " ; Lnormalize = "I2)
0107 format( 36X                               " ; weight_bnorm = "es12.5" ; weight_tflux ="es12.5" ; target_tflux ="es12.5)
0108 format( 36X                               " ; weight_ttlen = "es12.5" ; weight_eqarc ="es12.5" ; weight_ccsep ="es12.5)

0109 format("initial : " 10x " : Lpoincare ="i2" ;")
0110 format("initial : " 10x " : Lpoincare ="i2" ; odetol ="es8.1" ; Ppts ="i6" ; Ptrj ="i4" ; iphi ="i4" ; bstol ="es8.1" ; bsnlimit ="i9" ;")

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  select case( Itopology )
  case( 0 )
  case( 1 )
     FATAL( initial, knotsurf.lt.  zero, illegal)
  case( 2 )
  case default
     FATAL( initial, .true., selected Itopology is not supported )
  end select

  FATAL( initial, NFcoil .le. 0, illegal )
  FATAL( initial, NDcoil .le. 0, illegal )

  select case( Linitialize )
  case(:-2 )
  case( -1 )
     inquire( file="coils."//trim(ext), exist=exist )
     FATAL( initial, .not.exist, coils file coils.ext not provided )
  case( 0  )
  case( 1: )
     FATAL( initial, Ic.lt.0    .or. Ic.gt.1, illegal )
     FATAL( initial, Iw.lt.zero             , illegal )
     FATAL( initial, Lc.lt.0                , illegal )
     FATAL( initial, Lo.lt.zero             , illegal )
     FATAL( initial, Lw.lt.zero             , illegal )
  case default
     FATAL( initial, .true., selected Linitialize is not supported )
  end select

  select case( Loptimize )
  case( -2 )
     FATAL( initial, weight_bnorm  .lt.zero, illegal )
     FATAL( initial, weight_tflux  .lt.zero, illegal )
!     FATAL( initial, target_tflux  .lt.zero, illegal )
     FATAL( initial, weight_ttlen  .lt.zero, illegal )
     FATAL( initial, Nteta   .le.   0, illegal )
     FATAL( initial, Nzeta   .le.   0, illegal )
     !FATAL( initial, absacc  .le.zero, illegal )
     FATAL( initial, absreq  .le.zero, illegal )
     FATAL( initial, relreq  .le.zero, illegal )
  case( -1 )
     FATAL( initial, weight_bnorm  .lt.zero, illegal )
     FATAL( initial, weight_tflux  .lt.zero, illegal )
!     FATAL( initial, target_tflux  .lt.zero, illegal )
     FATAL( initial, weight_ttlen  .lt.zero, illegal )
     FATAL( initial, Nteta   .le.   0, illegal )
     FATAL( initial, Nzeta   .le.   0, illegal )
     !FATAL( initial, absacc  .le.zero, illegal )
     FATAL( initial, absreq  .le.zero, illegal )
     FATAL( initial, relreq  .le.zero, illegal )
  case( 0 )
  case( 1 )
     FATAL( initial, tauend  .le.zero, illegal )
     FATAL( initial, tautol  .le.zero, illegal )
     FATAL( initial, Ntauout .le.   0, illegal )
     FATAL( initial, Nteta   .le.   0, illegal )
     FATAL( initial, Nzeta   .le.   0, illegal )
     FATAL( initial, absacc  .le.zero, illegal )
     FATAL( initial, absreq  .le.zero, illegal )
     FATAL( initial, relreq  .le.zero, illegal )
  case( 2 )
     FATAL( initial, weight_bnorm  .lt.zero, illegal )
     FATAL( initial, weight_tflux  .lt.zero, illegal )
!     FATAL( initial, target_tflux  .lt.zero, illegal )
     FATAL( initial, weight_ttlen  .lt.zero, illegal )
     FATAL( initial, weight_eqarc  .lt.zero, illegal )
     FATAL( initial, weight_ccsep  .lt.zero, illegal )
     FATAL( initial, tauend  .le.zero, illegal )
     FATAL( initial, tautol  .le.zero, illegal )
     FATAL( initial, Ntauout .le.   0, illegal )
     FATAL( initial, Nteta   .le.   0, illegal )
     FATAL( initial, Nzeta   .le.   0, illegal )
     FATAL( initial, absacc  .le.zero, illegal )
     FATAL( initial, absreq  .le.zero, illegal )
     FATAL( initial, relreq  .le.zero, illegal )
  case ( 3 )
     FATAL( initial, weight_bnorm  .lt.zero, illegal )
     FATAL( initial, weight_tflux  .lt.zero, illegal )
!     FATAL( initial, target_tflux  .lt.zero, illegal )
     FATAL( initial, weight_ttlen  .lt.zero, illegal )
     FATAL( initial, weight_eqarc  .lt.zero, illegal )
     FATAL( initial, weight_ccsep  .lt.zero, illegal )
     FATAL( initial, Ntauout .le.   0, illegal )
     FATAL( initial, Nteta   .le.   0, illegal )
     FATAL( initial, Nzeta   .le.   0, illegal )
  case ( 4 )
     FATAL( initial, xtol .lt. zero  , illegal )
     FATAL( initial, eta .lt. zero .or. eta .ge. one, illegal )
     FATAL( initial, stepmx .lt. xtol, illegal )
     NT_Niter = Ntauout
  case ( 5 )
     FATAL( initial, weight_bnorm  .lt.zero, illegal )
     FATAL( initial, weight_tflux  .lt.zero, illegal )
     FATAL( initial, weight_ttlen  .lt.zero, illegal )
     FATAL( initial, weight_eqarc  .lt.zero, illegal )
     FATAL( initial, weight_ccsep  .lt.zero, illegal )
     FATAL( initial, Ntauout .le.   0, illegal )
     FATAL( initial, Nteta   .le.   0, illegal )
     FATAL( initial, Nzeta   .le.   0, illegal )
     CG_Niter = Ntauout
  case ( 6 )
     FATAL( initial, weight_bnorm  .lt.zero, illegal )
     FATAL( initial, weight_tflux  .lt.zero, illegal )
     FATAL( initial, weight_ttlen  .lt.zero, illegal )
     FATAL( initial, weight_eqarc  .lt.zero, illegal )
     FATAL( initial, weight_ccsep  .lt.zero, illegal )
     FATAL( initial, Ntauout .le.   0, illegal )
     FATAL( initial, Nteta   .le.   0, illegal )
     FATAL( initial, Nzeta   .le.   0, illegal )
     CG_Niter = Ntauout
     NT_Niter = Ntauout
     Ntauout = CG_Niter + NT_Niter
  case ( 9 )
     FATAL( initial, weight_bnorm  .lt.zero, illegal )
     FATAL( initial, weight_tflux  .lt.zero, illegal )
     FATAL( initial, weight_ttlen  .lt.zero, illegal )
     FATAL( initial, weight_eqarc  .lt.zero, illegal )
     FATAL( initial, weight_ccsep  .lt.zero, illegal )
     FATAL( initial, Ntauout .le.   0, illegal )
     FATAL( initial, Nteta   .le.   0, illegal )
     FATAL( initial, Nzeta   .le.   0, illegal )
     CG_Niter = Ntauout
     NT_Niter = Ntauout
     Ntauout = CG_Niter + NT_Niter

  case ( 10 )
     FATAL( initial, weight_bnorm  .lt.zero, illegal )
     FATAL( initial, weight_tflux  .lt.zero, illegal )
     FATAL( initial, weight_ttlen  .lt.zero, illegal )
     FATAL( initial, weight_eqarc  .lt.zero, illegal )
     FATAL( initial, weight_ccsep  .lt.zero, illegal )
     FATAL( initial, Ntauout .le.   0, illegal )
     FATAL( initial, Nteta   .le.   0, illegal )
     FATAL( initial, Nzeta   .le.   0, illegal )

  case default
     FATAL( initial, .true., selected Loptimize is not supported )
  end select

  FATAL( initial, lc .eq. 0 .and. weight_ttlen .ne. zero, conflicts between lc and weight_ttlen)

  if( Lpoincare.ne.0 ) then
   FATAL( initial, odetol   .le.zero  , illegal )
   FATAL( initial, Ppts     .lt.0     , illegal )
   FATAL( initial, Ptrj     .lt.0     , illegal )
   FATAL( initial, bstol    .le.zero  , illegal )
   FATAL( initial, bsnlimit .le.   0  , illegal )
   FATAL( initial, iphi     .gt. Nzeta, illegal )
  endif

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  write(nodelabel,'(i3.3)') myid ! nodelabel is global; 30 Oct 15;

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  !SALLOCATE( xsurf, (0:NDcoil -1,0:NDcoil -1), zero )
  !SALLOCATE( ysurf, (0:NDcoil -1,0:NDcoil -1), zero )
  !SALLOCATE( zsurf, (0:NDcoil -1,0:NDcoil -1), zero )

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  itau = 0

  nrestart = -1 ; bmn = 0 ; discretefactor = (pi2/Nteta) * (pi2/Nzeta); Cdof = 6*NFcoil + 6 ! DoF of per coil

  SALLOCATE( evolution, (0:Ntauout,0:9), zero )

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  ! construct a array shudson/newton to convert different representation of array index in each coil
!!$  if ( Loptimize .eq. 3 ) then
!!$     ncdof = 1 + 3 + 8 * NFcoil
!!$     SALLOCATE( newton, (1:ncdof), zero)
!!$     idof = 0
!!$     ;            ; idof = idof + 1 ; newton(idof) =             0 ! coil(icoil)%I
!!$     ; mm = 0     ; idof = idof + 1 ; newton(idof) =             1 ! coil(icoil)%xc( 0)
!!$     ;            ; idof = idof + 1 ; newton(idof) =    2*NFcoil+3 ! coil(icoil)%yc( 0)
!!$!    ;            ; idof = idof + 1 ; newton(idof) =    4*NFcoil+5 ! coil(icoil)%zc( 0)
!!$     ;            ; idof = idof + 1 ; newton(idof) =    6*NFcoil+7 ! coil(icoil)%lmdc( 0)
!!$     do mm = 1, NFcoil ; idof = idof + 1 ; newton(idof) = mm+         1 ! coil(icoil)%xc(mm)
!!$        ;              ; idof = idof + 1 ; newton(idof) = mm+2*NFcoil+3 ! coil(icoil)%yc(mm)
!!$        ;              ; idof = idof + 1 ; newton(idof) = mm+4*NFcoil+5 ! coil(icoil)%zc(mm)
!!$        ;              ; idof = idof + 1 ; newton(idof) = mm+6*NFcoil+7 ! coil(icoil)%lmdc(mm)
!!$        ;              ; idof = idof + 1 ; newton(idof) = mm+  NFcoil+2 ! coil(icoil)%xs(mm)
!!$        ;              ; idof = idof + 1 ; newton(idof) = mm+3*NFcoil+4 ! coil(icoil)%ys(mm)
!!$        ;              ; idof = idof + 1 ; newton(idof) = mm+5*NFcoil+6 ! coil(icoil)%zs(mm)
!!$        ;              ; idof = idof + 1 ; newton(idof) = mm+7*NFcoil+8 ! coil(icoil)%lmds(mm)
!!$     enddo ! end of do mm
!!$     FATAL( initial, idof .ne. ncdof, counting error in constructing newton array )
!!$  else

     ncdof = 1 + 3 + 6 * NFcoil;
     SALLOCATE( shudson, (1:ncdof), zero)
     idof = 0
     ;            ; idof = idof + 1 ; shudson(idof) =             0 ! coil(icoil)%I
     ; mm = 0     ; idof = idof + 1 ; shudson(idof) =             1 ! coil(icoil)%xc( 0)
     ;            ; idof = idof + 1 ; shudson(idof) =    2*NFcoil+3 ! coil(icoil)%yc( 0)
     ;            ; idof = idof + 1 ; shudson(idof) =    4*NFcoil+5 ! coil(icoil)%zc( 0)
     do mm = 1, NFcoil ; idof = idof + 1 ; shudson(idof) = mm+         1 ! coil(icoil)%xc(mm)
        ;            ; idof = idof + 1 ; shudson(idof) = mm+2*NFcoil+3 ! coil(icoil)%yc(mm)
        ;            ; idof = idof + 1 ; shudson(idof) = mm+4*NFcoil+5 ! coil(icoil)%zc(mm)
        ;            ; idof = idof + 1 ; shudson(idof) = mm+  NFcoil+2 ! coil(icoil)%xs(mm)
        ;            ; idof = idof + 1 ; shudson(idof) = mm+3*NFcoil+4 ! coil(icoil)%ys(mm)
        ;            ; idof = idof + 1 ; shudson(idof) = mm+5*NFcoil+6 ! coil(icoil)%zs(mm)
     enddo ! end of do mm
     FATAL( initial, idof .ne. ncdof, counting error in constructing shudson array )
!  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( Loptimize .eq. 3 .and. weight_eqarc .le. zero ) then

     weight_eqarc = one
     if(myid .eq. 0) write(ounit,'("initial : " 10x " : When using Newton method, the eqarc constraint is automaticallly turned on.")')

  endif

  !save weights before normalized
  tmpw_bnorm = weight_bnorm
  tmpw_tflux = weight_tflux
  tmpt_tflux = target_tflux
  tmpw_ttlen = weight_ttlen
  tmpw_eqarc = weight_eqarc
  tmpw_ccsep = weight_ccsep
  tmpw_qasym = weight_qasym

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
 !Nt = max( Ndiscrete*4*Mpol, 1 ) ; Nz = max( Ndiscrete*4*Ntor, 1 ) ; Ntz = Nt*Nz ; soNtz = one / sqrt( one*Ntz ) ! exaggerated discrete resolution;

  Ntz = Nteta * Nzeta

  mn = 1 + Ntor + abs(Mpol) * ( 2 * Ntor + 1 ) ! Fourier resolution of interface geometry & vector potential;
  
  SALLOCATE( im, (1:mn), 0 )
  SALLOCATE( in, (1:mn), 0 )
  
  imn = 0 ; Nfp = 1 ! TAKE CARE WITH NFP; 18 Apr 17;
  
  ;  mm = 0
  ;do nn = 0, Ntor
  ; imn = imn+1 ; im(imn) = mm ; in(imn) = nn * Nfp
  ;enddo
  ;
  
  do mm = 1, abs(Mpol)
   do nn = -Ntor, Ntor
    imn = imn+1 ; im(imn) = mm ; in(imn) = nn * Nfp
   enddo
  enddo

  SALLOCATE( efmn, (1:mn), zero ) ! workspace for Fourier harmonics; 18 Apr 17;
  SALLOCATE( ofmn, (1:mn), zero )
  SALLOCATE( cfmn, (1:mn), zero )
  SALLOCATE( sfmn, (1:mn), zero )
  
  SALLOCATE( ijreal, (1:Ntz), zero ) ! workspace for real space grid; 18 Apr 17;
  SALLOCATE( ijimag, (1:Ntz), zero )
  SALLOCATE( jireal, (1:Ntz), zero )
  SALLOCATE( jiimag, (1:Ntz), zero )

  SALLOCATE( trigm , (1:2*Nteta) , zero ) ! trignometric factors required for fast Fourier transform;
  SALLOCATE( trign , (1:2*Nzeta) , zero )
  SALLOCATE( trigwk, (1:2*Ntz), zero )

  isr = 'I' ; ifail = 0

  ! call C06FUF( Nteta, Nzeta, ijreal(1:Ntz), ijimag(1:Ntz), isr, trigm(1:2*Nteta), trign(1:2*Nzeta), trigwk(1:2*Ntz), ifail )
  ! comment out on 20180228 for NAG incompative

  isr = 'S'

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  return

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

end subroutine initial
