
!title (coil freedom) ! Determines degrees-of-freedom in coil representation.

!latex \briefly{Determines degrees-of-freedom in coil representation.}

!latex \calledby{\link{xoptim}}
!l tex \calls{\link{}}

!latex \tableofcontents

!latex \subsection{overview}


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine setdof( Ndof, Xdof, packorunpack, lowerbound, upperbound )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  use knotopt , only : zero, half, one, two, myid, ncpu, ounit, &
                       Nfp, &
                       AddCurrent, Osymmetric, Hsymmetric, Nsympairs, NFourier, &
                       Ncoils, &
                       Vcurrents, Vzetaoff, &
                       itc, its, izc, izs, icu, &
                       xtc, xts, xzc, xzs, xcu
  
  implicit none
  
  include "mpif.h"
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  INTEGER               :: Ndof
  REAL                  :: Xdof(1:Ndof), lowerbound(1:Ndof), upperbound(1:Ndof)
  CHARACTER             :: packorunpack
  
  INTEGER               :: astat, ierr, idof, nn, icoil, ipair, sdof, ifp
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  select case( packorunpack )
   
  case( 'C' ) ! counting degrees-of-freedom; 30 Oct 15;
   
   SALLOCATE( icu, (           1:Ncoils), 0 )
   SALLOCATE( itc, (0:NFourier,1:Ncoils), 0 )
   SALLOCATE( its, (0:NFourier,1:Ncoils), 0 )
   SALLOCATE( izc, (0:NFourier,1:Ncoils), 0 )
   SALLOCATE( izs, (0:NFourier,1:Ncoils), 0 )
   
   idof = 0 ; icoil = AddCurrent
   
   if( Osymmetric.eq.1 ) then ! first  stellarator-symmetric coil; 30 Oct 15;
    
    if( Vcurrents.eq.1 ) then
     idof = idof + 1 ; icu(   icoil+1:icoil+Nfp) = + idof ;                                            ; lowerbound(idof) = zero ; upperbound(idof) = zero
    endif

    do nn = 1, NFourier
     idof = idof + 1 ; izs(nn,icoil+1:icoil+Nfp) = + idof
     idof = idof + 1 ; izc(nn,icoil+1:icoil+Nfp) = + idof
    enddo
    
    icoil = icoil + Nfp
    
   endif! end of if( Osymmetric.eq.1 ) ; 04 Nov 15;
   
   if( Hsymmetric.eq.1 ) then ! second stellarator-symmetric coil; 30 Oct 15;
    
    if( Vcurrents.eq.1 ) then
     idof = idof + 1 ; icu(   icoil+1:icoil+Nfp) = + idof ;                                            ; lowerbound(idof) = zero ; upperbound(idof) = zero
    endif

    do nn = 1, NFourier
     idof = idof + 1 ; izs(nn,icoil+1:icoil+Nfp) = + idof
     idof = idof + 1 ; izc(nn,icoil+1:icoil+Nfp) = + idof
    enddo
    
    icoil = icoil + Nfp
    
   endif ! end of if( Hsymmetric.eq.1 ) ; 04 Nov 15;
   
   do ipair = 1, Nsympairs
    
    if( Vcurrents.eq.1 ) then
     idof = idof + 1 ; icu(   icoil+1:icoil+Nfp) = + idof ; icu(   icoil+1+Nfp:icoil+Nfp+Nfp) = + idof ; lowerbound(idof) = zero ; upperbound(idof) = zero
    endif

    if( Vzetaoff.eq.1 ) then
     ; nn = 0
     idof = idof + 1 ; izc(nn,icoil+1:icoil+Nfp) = + idof ; izc(nn,icoil+1+Nfp:icoil+Nfp+Nfp) = - idof
    endif
    
    do nn = 1, NFourier
     idof = idof + 1 ; izs(nn,icoil+1:icoil+Nfp) = + idof ; izs(nn,icoil+1+Nfp:icoil+Nfp+Nfp) = + idof
     idof = idof + 1 ; izc(nn,icoil+1:icoil+Nfp) = + idof ; izc(nn,icoil+1+Nfp:icoil+Nfp+Nfp) = - idof
    enddo
    
    icoil = icoil + Nfp + Nfp
    
   enddo ! end of do ipair; 30 Oct 15;

   if( myid.eq.0 ) write(ounit,'("setdof : " 10x " : labelled  degrees-of-freedom ; idof ="i6" ; Ndof ="i6" ;")') idof, Ndof

   FATAL( setdof, idof .ne.Ndof  , counting error )
   FATAL( setdof, icoil.ne.Ncoils, counting error )
   
  case( 'P' ) ! packing; 30 Oct 15;
   
   do icoil = 1, Ncoils
    
    ;idof = icu(   icoil) ; sdof = sign(1,idof) ; if( idof.ne.0 ) Xdof(abs(idof)) =   xcu(   icoil) * sdof
    do nn = 0, NFourier
     idof = its(nn,icoil) ; sdof = sign(1,idof) ; if( idof.ne.0 ) Xdof(abs(idof)) =   xts(nn,icoil) * sdof
     idof = itc(nn,icoil) ; sdof = sign(1,idof) ; if( idof.ne.0 ) Xdof(abs(idof)) =   xtc(nn,icoil) * sdof
     idof = izs(nn,icoil) ; sdof = sign(1,idof) ; if( idof.ne.0 ) Xdof(abs(idof)) =   xzs(nn,icoil) * sdof
     idof = izc(nn,icoil) ; sdof = sign(1,idof) ; if( idof.ne.0 ) Xdof(abs(idof)) =   xzc(nn,icoil) * sdof
    enddo ! end of do nn; 30 Oct 15;

   enddo ! end of do icoil; 30 Oct 15;

   if( myid.eq.0 ) write(ounit,'("setdof : " 10x " : packed    degrees-of-freedom ;")') 
   
  case( 'U' ) ! un-packing; 30 Oct 15;
   
   do icoil = 1, Ncoils
    
    ;idof = icu(   icoil) ; sdof = sign(1,idof) ; if( idof.ne.0 ) xcu(   icoil) = Xdof(abs(idof)) / sdof
    do nn = 0, NFourier
     idof = its(nn,icoil) ; sdof = sign(1,idof) ; if( idof.ne.0 ) xts(nn,icoil) = Xdof(abs(idof)) / sdof
     idof = itc(nn,icoil) ; sdof = sign(1,idof) ; if( idof.ne.0 ) xtc(nn,icoil) = Xdof(abs(idof)) / sdof
     idof = izs(nn,icoil) ; sdof = sign(1,idof) ; if( idof.ne.0 ) xzs(nn,icoil) = Xdof(abs(idof)) / sdof
     idof = izc(nn,icoil) ; sdof = sign(1,idof) ; if( idof.ne.0 ) xzc(nn,icoil) = Xdof(abs(idof)) / sdof
    enddo ! end of do nn; 30 Oct 15;

   enddo ! end of do icoil; 30 Oct 15;
   
  case default
   
   FATAL( setdof, .true., illegal packorunpack )
   
  end select
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
end subroutine setdof

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
