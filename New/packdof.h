!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

!title (packdof) ! paking degree of freedom (dof) into one vector

!latex  \briefly{Packing all the free coil parameters into a one rank vector.}
!latex 
!latex  \subsection{Overview}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

subroutine packdof( Ndof, xdof, packorunpack )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  use globals, only : zero, myid, ounit, &
                      Ncoils, coil, tdof, pspectral, Mdof
  
  implicit none
  
  include "mpif.h"
  
  INTEGER  , intent(in) :: Ndof
  REAL                  :: xdof(1:Ndof)
  CHARACTER, intent(in) :: packorunpack
  
  INTEGER :: icoil, idof, mm, NF, ierr, ideriv
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  select case( packorunpack )
   
  case('P') ! packing curve parameters that are degrees-of-freedom into single vector;
   
   idof = 0 ! labels degree-of-freedom;
   
   do icoil = 1, Ncoils ; NF = coil(icoil)%NF
    
    if( coil(icoil)%Lfree.eq.1 ) then ; coil(icoil)%gdof = idof ! geometry of curve is allowed to vary;
     
     ;  mm = 0     ; idof = idof + 1 ; xdof(idof) = coil(icoil)%xc(mm)
     ;             ; idof = idof + 1 ; xdof(idof) = coil(icoil)%yc(mm)
     ;             ; idof = idof + 1 ; xdof(idof) = coil(icoil)%zc(mm)
     do mm = 1, NF ; idof = idof + 1 ; xdof(idof) = coil(icoil)%xc(mm)
      ;            ; idof = idof + 1 ; xdof(idof) = coil(icoil)%xs(mm)
      ;            ; idof = idof + 1 ; xdof(idof) = coil(icoil)%yc(mm)
      ;            ; idof = idof + 1 ; xdof(idof) = coil(icoil)%ys(mm)
      ;            ; idof = idof + 1 ; xdof(idof) = coil(icoil)%zc(mm)
      ;            ; idof = idof + 1 ; xdof(idof) = coil(icoil)%zs(mm)
     enddo
     
    !xdof(idof+     1:idof+  NF+1) = coil(icoil)%xc(0:NF)
    !xdof(idof+  NF+2:idof+2*NF+1) = coil(icoil)%xs(1:NF)
    !xdof(idof+2*NF+2:idof+3*NF+2) = coil(icoil)%yc(0:NF)
    !xdof(idof+3*NF+3:idof+3*NF+2) = coil(icoil)%ys(1:NF)
    !xdof(idof+3*NF+3:idof+4*NF+3) = coil(icoil)%zc(0:NF)
    !xdof(idof+4*NF+4:idof+5*NF+2) = coil(icoil)%zs(1:NF)
    !idof = idof + 6*NF+3
     
    endif
    
    if( coil(icoil)%Ifree.eq.1 ) then ; coil(icoil)%idof = idof ! current of curve is allowed to vary;
     
     ;             ; idof = idof + 1 ; xdof(idof) = coil(icoil)%I

    !xdof(idof+     1:idof+     1) = coil(icoil)%I
    !idof = idof + 1
     
    endif
    
   enddo ! end of do icoil;
   
  case('U') ! unpacking degrees-of-freedom to curve parameters;
   
   idof = 0
   
   do icoil = 1, Ncoils ; NF = coil(icoil)%NF
    
    if( coil(icoil)%Lfree.eq.1 ) then ! geometry of curve is allowed to vary;
     
     ;  mm = 0     ; idof = idof + 1 ; coil(icoil)%xc(mm) = xdof(idof)
     ;             ; idof = idof + 1 ; coil(icoil)%yc(mm) = xdof(idof)
     ;             ; idof = idof + 1 ; coil(icoil)%zc(mm) = xdof(idof)
     do mm = 1, NF ; idof = idof + 1 ; coil(icoil)%xc(mm) = xdof(idof)
      ;            ; idof = idof + 1 ; coil(icoil)%xs(mm) = xdof(idof)
      ;            ; idof = idof + 1 ; coil(icoil)%yc(mm) = xdof(idof)
      ;            ; idof = idof + 1 ; coil(icoil)%ys(mm) = xdof(idof)
      ;            ; idof = idof + 1 ; coil(icoil)%zc(mm) = xdof(idof)
      ;            ; idof = idof + 1 ; coil(icoil)%zs(mm) = xdof(idof)
     enddo

    !coil(icoil)%xc(0:NF) = xdof(idof+     1:idof+  NF+1)
    !coil(icoil)%xs(1:NF) = xdof(idof+  NF+2:idof+2*NF+1)
    !coil(icoil)%yc(0:NF) = xdof(idof+2*NF+2:idof+3*NF+2)
    !coil(icoil)%ys(1:NF) = xdof(idof+3*NF+3:idof+3*NF+2)
    !coil(icoil)%zc(0:NF) = xdof(idof+3*NF+3:idof+4*NF+3)
    !coil(icoil)%zs(1:NF) = xdof(idof+4*NF+4:idof+5*NF+2)
    !idof = idof + 6*NF+3

     call coilxyz( icoil ) ! map curve to real space;
     
    endif
    
    if( coil(icoil)%Ifree.eq.1 ) then ! current of curve is allowed to vary;
     
     ;             ; idof = idof + 1 ; coil(icoil)%I      = xdof(idof)

    !coil(icoil)%I        = xdof(idof+     1            )
    !idof = idof + 1
     
    endif
    
   enddo ! end of do icoil;

   CHECK( packdof, .not.allocated(tdof), error )
   
   idof = 0
   
   do icoil = 1, Ncoils ; NF = coil(icoil)%NF
    
    if( coil(icoil)%Lfree.eq.1 ) then ! geometry of curve is allowed to vary;

     ;  mm = 0     ; idof = idof + 1 ; tdof(idof) = zero                      ; Mdof(idof) = zero
     ;             ; idof = idof + 1 ; tdof(idof) = zero                      ; Mdof(idof) = zero
     ;             ; idof = idof + 1 ; tdof(idof) = zero                      ; Mdof(idof) = zero
     do mm = 1, NF ; idof = idof + 1 ; tdof(idof) = + mm * coil(icoil)%xs(mm) ; Mdof(idof) = mm**pspectral * coil(icoil)%xc(mm) 
      ;            ; idof = idof + 1 ; tdof(idof) = - mm * coil(icoil)%xc(mm) ; Mdof(idof) = mm**pspectral * coil(icoil)%xs(mm) 
      ;            ; idof = idof + 1 ; tdof(idof) = + mm * coil(icoil)%ys(mm) ; Mdof(idof) = mm**pspectral * coil(icoil)%yc(mm) 
      ;            ; idof = idof + 1 ; tdof(idof) = - mm * coil(icoil)%yc(mm) ; Mdof(idof) = mm**pspectral * coil(icoil)%ys(mm) 
      ;            ; idof = idof + 1 ; tdof(idof) = + mm * coil(icoil)%zs(mm) ; Mdof(idof) = mm**pspectral * coil(icoil)%zc(mm) 
      ;            ; idof = idof + 1 ; tdof(idof) = - mm * coil(icoil)%zc(mm) ; Mdof(idof) = mm**pspectral * coil(icoil)%zs(mm) 
     enddo

    endif

    if( coil(icoil)%Ifree.eq.1 ) then ! current of curve is allowed to vary;
     
     ;             ; idof = idof + 1 ; tdof(idof) = zero                      ; Mdof(idof) = zero

    !coil(icoil)%I        = xdof(idof+     1            )
    !idof = idof + 1
     
    endif
   
   enddo ! end of do icoil;
 
  case default
   
   FATAL( packdof, .true., selected packorunpack not supported )

  end select
  
  return
  
end subroutine packdof

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
