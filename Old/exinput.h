
!title (diagnostic) ! Prepare input information for other codes.

!latex \briefly{Prepare input for M3D-$C^1$, SPEC, mgrid, etc.}

!latex \calledby{\link{notopt}}

!latex \tableofcontents

!latex \subsection{overview}
!latex \bi
!latex \item[1.] SPEC input information is constructed using \oculus{bn00aa}.
!latex \ei

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine exinput
  
  use kmodule, only : zero, half, pi2, ext, ncpu, myid, ounit, lunit, funit, Lpoincare                     
  
! use oculus , only : bnormal, bn00aa
  
  implicit none
  
  include "mpif.h"
  
!!$!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!!$  
!!$  LOGICAL              :: exist
!!$  INTEGER              :: ierr, astat, iostat, imn, ibn00aa, Nt, Nz, Ntz, Mvol, Igeometry , vvol, Lrad, ii, jj, kk, jk, itangent, ibfield, nijk
!!$  REAL                 :: pi2nfp, teta, zeta, RpZ(1:3), dBRpZ(1:3,0:3), gdB
!!$  REAL, allocatable    :: RR(:), ZZ(:), sg(:), BR(:), Bp(:), BZ(:), le(:)
!!$  type(bnormal)        :: bn
!!$  
!!$!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!!$  
!!$  if( myid.ne.0 ) return
!!$  
!!$!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!!$  
!!$  inquire( file="."//trim(ext)//".sp.grid", exist=exist)
!!$  
!!$  if( exist ) then
!!$   
!!$   write(ounit,'("exinput : " 10x " : reading .ext.sp.grid ; computing Biot-Savart field ; ")') 
!!$   
!!$   open(lunit, file="."//trim(ext)//".sp.grid", status="old"    , form="unformatted" ) ! coordinate grid;
!!$
!!$   open(funit, file="."//trim(ext)//".fo.grid", status="unknown", form="unformatted" ) ! coordinate grid;
!!$
!!$   read(lunit) Nt, Nz, Ntz, Mvol, Igeometry, pi2nfp
!!$   
!!$   SALLOCATE( RR, (1:Ntz), zero )
!!$   SALLOCATE( ZZ, (1:Ntz), zero )
!!$   SALLOCATE( sg, (1:Ntz), zero )
!!$   SALLOCATE( BR, (1:Ntz), zero )
!!$   SALLOCATE( Bp, (1:Ntz), zero )
!!$   SALLOCATE( BZ, (1:Ntz), zero )
!!$
!!$   SALLOCATE( le, (1:Ntz), zero ) ! local error; 05 Jan 17;
!!$   
!!$   itangent = 0 
!!$   
!!$   nijk = 0 ; gdB = zero
!!$   
!!$   do vvol = 1, Mvol
!!$
!!$    read(lunit) Lrad
!!$
!!$    do ii = 0, Lrad
!!$
!!$     write(ounit,'("exinput : ", 10x ," : vvol =",i3," ; ii =",i3," /",i3," ;")') vvol, ii, Lrad
!!$
!!$     read(lunit) RR(1:Ntz)
!!$     read(lunit) ZZ(1:Ntz)
!!$     read(lunit) sg(1:Ntz)
!!$     read(lunit) BR(1:Ntz)
!!$     read(lunit) Bp(1:Ntz) ; Bp(1:Ntz) = Bp(1:Ntz) * RR(1:Ntz)
!!$     read(lunit) BZ(1:Ntz)
!!$     
!!$     if( ii.eq.0 ) cycle
!!$     
!!$     do kk = 0, Nz-1 ; zeta = kk * pi2nfp / Nz
!!$      do jj = 0, Nt-1 ; teta = jj * pi2    / Nt ; jk = 1 + jj + kk*Nt
!!$       
!!$       RpZ(1:3) = (/ RR(jk), zeta, ZZ(jk) /)
!!$       
!!$       call bfield( RpZ(1:3), itangent, dBRpZ(1:3,0:3), ibfield ) ; dBRpZ(2,0) = dBRpZ(2,0) * RR(jk)
!!$       
!!$      !write(ounit,'("exinput : " 10x " : BR="2es13.5" ; |dBR|="es13.5" ; Bp="2es13.5" ; |dBp|="es13.5" ; BZ="2es13.5" ; |dBZ|="es13.5" ;")') &
!!$   !BR(jk), dBRpZ(1,0), BR(jk)-dBRpZ(1,0), Bp(jk), dBRpZ(2,0), Bp(jk)-dBRpZ(2,0), BZ(jk), dBRpZ(3,0), BZ(jk)-dBRpZ(3,0)
!!$       
!!$       nijk = nijk + 1
!!$
!!$       le(jk) = sqrt( (BR(jk)-dBRpZ(1,0))**2 + (Bp(jk)-dBRpZ(2,0))**2 + (BZ(jk)-dBRpZ(3,0))**2 ) / sqrt( BR(jk)**2 + Bp(jk)**2 + BZ(jk)**2 )
!!$       
!!$       gdB = gdB + sg(jk) * le(jk)
!!$
!!$      enddo ! end of do jj; 05 Jan 17;
!!$     enddo ! end of do kk; end of jk; 05 Jan 17;
!!$     
!!$     write(funit) le(1:Ntz)     
!!$
!!$    enddo ! end of do ii; 29 Jan 13;
!!$
!!$   enddo ! end of do vvol; 29 Jan 13;
!!$   
!!$   close(lunit)
!!$   close(funit)
!!$   
!!$   deallocate( RR )
!!$   deallocate( ZZ )
!!$   deallocate( sg )
!!$   deallocate( BR )
!!$   deallocate( Bp )
!!$   deallocate( BZ )
!!$   
!!$   write(ounit,'("exinput : " 10x " : |B_s-B_v| =",es22.15," ; wrote .ext.fo.grid ;")') gdB / nijk
!!$
!!$  else
!!$
!!$   write(ounit,'("exinput : " 10x " : .ext.sp.grid not found ;")')
!!$
!!$  endif
!!$  
!!$!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!!$  
!!$  inquire( file=trim(ext)//".control.surface", exist=exist)
!!$  
!!$  FATAL( exinput, .not.exist, ext.control.surface does not exist )
!!$  
!!$  write( ounit,'("exinput : " 10x " : reading ext.control.surface ;")') 
!!$  
!!$  open( lunit, file=trim(ext)//".control.surface", status='old', action='read' )
!!$  
!!$  read( lunit, * ) bn%cmn, bn%Nfp
!!$  
!!$  FATAL( exinput, bn%cmn.le.0, invalid )
!!$  FATAL( exinput, bn%Nfp.le.0, invalid )
!!$  
!!$  SALLOCATE( bn%cim, (1:bn%cmn), 0 )
!!$  SALLOCATE( bn%cin, (1:bn%cmn), 0 )
!!$  
!!$  SALLOCATE( bn%Rcc, (1:bn%cmn), zero )
!!$  SALLOCATE( bn%Rcs, (1:bn%cmn), zero )
!!$  SALLOCATE( bn%Zcc, (1:bn%cmn), zero )
!!$  SALLOCATE( bn%Zcs, (1:bn%cmn), zero )
!!$  
!!$  do imn = 1, bn%cmn ; read( lunit, * ) bn%cin(imn), bn%cim(imn), bn%Rcc(imn), bn%Rcs(imn), bn%Zcc(imn), bn%Zcs(imn)
!!$  enddo
!!$  
!!$  bn%cin(1:bn%cmn) = bn%cin(1:bn%cmn) * bn%Nfp
!!$  
!!$  close(lunit)
!!$  
!!$  write( ounit,'("exinput : " 10x " : calling Oculus:bn00aa ;")')
!!$  
!!$  bn%Mpol = 20 ; bn%Ntor = 10 ; bn%tol = 1.0e-14 ! this should be user input; 04 Aug 16;
!!$  
!!$  ibn00aa =  0 ; call bn00aa( bn, ibn00aa )
!!$  
!!$  write( ounit,'("exinput : " 10x " : writing ext.Vn ; input for SPEC (free-boundary) ;")')
!!$  
!!$  open( lunit, file=trim(ext)//".Vn", status="unknown", form="formatted", iostat=iostat )
!!$  FATAL( exinput, iostat.ne.0, error opening ext.Vn ) 
!!$  write( lunit,'(2i6         )') bn%cmn, bn%Nfp
!!$  write( lunit,'(2i6,4es23.15)') ( bn%cim(imn), bn%cin(imn)/bn%Nfp, bn%Rcc(imn), bn%Rcs(imn), bn%Zcc(imn), bn%Zcs(imn), imn = 1, bn%cmn )
!!$  write( lunit,'(2i6,2es23.15)') bn%mn , bn%Nfp, bn%Itor, bn%Gpol
!!$  write( lunit,'(2i6,2es23.15)') ( bn%im(imn) , bn%in(imn)/bn%Nfp , bn%gBc(imn), bn%gBs(imn),                           imn = 1, bn%mn  )
!!$  close( lunit )
!!$  
!!$ !open( lunit, file="bnorm."//trim(ext), status="unknown", form="formatted", iostat=iostat )
!!$ !FATAL( exinput, iostat.ne.0, error opening bnorm.ext ) 
!!$ !write( lunit,'(2i6,2es23.15)') ( bn%im(imn) , bn%in(imn)/bn%Nfp , bn%gBc(imn), bn%gBs(imn),                           imn = 1, bn%mn  )
!!$ !close( lunit )
!!$  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  Lpoincare = 0

  return
  
end subroutine exinput

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
