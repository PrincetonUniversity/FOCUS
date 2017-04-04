
!title (optimize) ! Nonlinear search.

!latex \briefly{Nonlinear search.}

!latex \tableofcontents

!latex \subsection{overview}
!latex \bi
!latex \item[1.] First, 
!latex \ei

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine search
  
  use knotopt
  
 !use oculus, only : pp00aa
  
  implicit none
  
  include "mpif.h"

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOGICAL              :: exist 
  INTEGER              :: ierr, astat, iostat, itrj, ipp00aa
  REAL                 :: aa, teta, zeta

  REAL                 :: ax(1:3), at(1:3), az(1:3), xx(1:3), xt(1:3), xz(1:3)

  INTEGER              :: Ibound, LIwork, LRwork, Iuser(1:1), ie04jyf ! optimization routine; 30 Oct 15;
  INTEGER, allocatable :: Iwork(:)
  REAL                 :: fBnsqd, Ruser(1:1)
  REAL   , allocatable :: xcdof(:), lowerbound(:), upperbound(:), Rwork(:)

  REAL                 :: ofunct
  external             :: ofunct

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! THE FOLLOWING MUST BE CONSISTENT WITH SETDOF; 04 Nov 15;
   
!                         xcu[i]    + xtc[0 ,i] + xzc[0 ,i] + xts[1:,i] + xtc[1:,i] + xzs[1:,i] + xzc[1:,i]
   
!  Ncdof = Osymmetric * ( Vcurrents + 0         + 0         + 0         + 0         + NFourier  + NFourier  )     &
!        + Hsymmetric * ( Vcurrents + 0         + 0         + 0         + 0         + NFourier  + NFourier  )     &
!        + Nsympairs  * ( Vcurrents + 0         + Vzetaoff  + 0         + 0         + NFourier  + NFourier  )

! THE ABOVE     MUST BE CONSISTENT WITH SETDOF; 04 Nov 15;
   
!   FATAL( notopt, Ncdof.le.0, invalid )
!   
!   SALLOCATE( xcdof, (1:Ncdof), zero )
!   
!   SALLOCATE( lowerbound, (1:Ncdof), zero ) ! lower bounds;
!   SALLOCATE( upperbound, (1:Ncdof), zero ) ! upper bounds;
!   
!   call setdof( Ncdof, xcdof(1:Ncdof), 'C', lowerbound(1:Ncdof), upperbound(1:Ncdof) ) ! label degrees-of-freedom; provide default values; 30 Oct 15;
!   
!   call setdof( Ncdof, xcdof(1:Ncdof), 'P', lowerbound(1:Ncdof), upperbound(1:Ncdof) ) ! "pack" degrees-of-freedom into single vector; 30 Oct 15;
!   
!   select case( Vcurrents )
!   case(  0 ) ; Ibound = 1 ; lowerbound(1:Ncdof) = zero ; upperbound(1:Ncdof) = zero ! no bounds;
!   case(  1 ) ; Ibound = 1 ; lowerbound(1      ) = zero ; upperbound(1      ) = zero ! user-supplied bounds; other bounds set elsewhere; 
!   end select
!      
!   if( myid.eq.0 ) open(funit,file="."//trim(ext)//".op.F.dat", status="unknown" ) ! this will save the |F| as a function of iteration; 30 Oct 15;
!   
!   
!   LIwork = Ncdof + 2
!   SALLOCATE(Iwork,(1:LIwork),   0) ! integer workspace; 30 Oct 15;
!   
!   LRwork = Ncdof * ( Ncdof - 1 ) / 2 + 12 * Ncdof
!   SALLOCATE(Rwork,(1:LRwork),zero) ! real    workspace; 30 Oct 15;
!   
!   Iuser(1:1)  =  (/ 0 /) ; Ruser(1:1) = (/ zero /) ! these are not used; 30 Oct 15;
!   
!   if( myid.eq.0 ) write(ounit,1000) nrestart, Ncdof
!   if( myid.eq.0 ) write(ounit,1001) nrestart, Ncdof, Wbnormal, Wclength, Wepotent, Wtorflux
!   
!1000 format("notopt : " 10x " : "i10" /"i6" :    "  22x  " ;"  23x  " ;          "3x" ; calling E04JYF ;")
!1001 format("notopt : " 10x " : "i10" /"i6" :    "  22x  " ;"  23x  " ; weights  "3x" ;    "es22.15",    "es22.15",    "es22.15",    "es22.15" ;")
!   
!   ie04jyf = 1
!   
!   call E04JYF( Ncdof,                                    &                                 ! NAG minimization algorithm ; 30 Oct 15;
!                Ibound,                                   &
!                ofunct,                                   &                                 ! |F| = \int_ds (B.n)^2 ds   ; 30 Oct 15;
!                lowerbound(1:Ncdof), upperbound(1:Ncdof), &
!                xcdof(1:Ncdof),                           &                                 ! "packed" degrees-of-freedom; 30 Oct 15;
!                fBnsqd,                                   &                                 !                            ; 04 Nov 15;
!                Iwork(1:LIwork), LIwork, Rwork(1:LRwork), LRwork, Iuser(1:1), Ruser(1:1), & ! workspace                  ; 30 Oct 15;
!                ie04jyf )
!   
!   DALLOCATE(Rwork)
!   DALLOCATE(Iwork)
!   
!   if( myid.eq.0 ) write(ounit,1002) nrestart, Ncdof, fBnsqd, fBnsqd-oBnsqd, ie04jyf
!   
!1002 format("notopt : " 10x " : "i10" /"i4" : F ="es22.15" ;"es23.15" ; ie04jyf ="i3" ;")
!   
!   if( myid.eq.0 ) close(funit)
!   
!   call setdof( Ncdof, xcdof(1:Ncdof), 'U', lowerbound(1:Ncdof), upperbound(1:Ncdof) ) ! "unpack" degrees-of-freedom; 30 Oct 15;
!   
!9000 continue
!   
!   DALLOCATE(lowerbound)
!   DALLOCATE(upperbound)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
   return

 end subroutine search

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
