!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

!title (packdof) ! paking degree of freedom (dof) into one vector

!latex  \briefly{Packing all the free coil parameters into a one rank vector.}
!latex 
!latex  \subsection{Overview}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

SUBROUTINE packdof( packorunpack )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  use globals, only : zero, myid, ounit, Ncoils, coil, DoF, Ndof, xdof, Inorm, Gnorm
  
  implicit none
  
  include "mpif.h"

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

  CHARACTER, intent(in) :: packorunpack

  INTEGER :: icoil, NF, idof, ierr
  
! REAL    :: lxdof(1:Ndof)
! INTEGER :: idof, icoil, ND, astat, ierr

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  select case( packorunpack )
   
  case('P')
   
   idof = 0 
   
   do icoil = 1, Ncoils ; NF = coil(icoil)%NF
    
    if( coil(icoil)%Lc == 0 ) then
     
     xdof(idof+     1:idof+  NF+1) = coil(icoil)%xc(0:NF)
     xdof(idof+  NF+2:idof+2*NF+1) = coil(icoil)%xs(1:NF)
     xdof(idof+2*NF+2:idof+3*NF+2) = coil(icoil)%yc(0:NF)
     xdof(idof+3*NF+3:idof+3*NF+2) = coil(icoil)%ys(1:NF)
     xdof(idof+3*NF+3:idof+4*NF+3) = coil(icoil)%zc(0:NF)
     xdof(idof+4*NF+4:idof+5*NF+2) = coil(icoil)%zs(1:NF)
     
     idof = idof + 6*NF+3
     
    endif
    
    if( coil(icoil)%Ic == 0 ) then
     
     xdof(idof+     1:idof+     1) = coil(icoil)%I
     
     idof = idof + 1
     
    endif
    
   enddo ! end of do icoil; 04 Sep 17;
   
  case('U')
   
   idof = 0 
   
   do icoil = 1, Ncoils ; NF = coil(icoil)%NF
    
    if( coil(icoil)%Lc == 0 ) then
     
     coil(icoil)%xc(0:NF) = xdof(idof+     1:idof+  NF+1)
     coil(icoil)%xs(1:NF) = xdof(idof+  NF+2:idof+2*NF+1)
     coil(icoil)%yc(0:NF) = xdof(idof+2*NF+2:idof+3*NF+2)
     coil(icoil)%ys(1:NF) = xdof(idof+3*NF+3:idof+3*NF+2)
     coil(icoil)%zc(0:NF) = xdof(idof+3*NF+3:idof+4*NF+3)
     coil(icoil)%zs(1:NF) = xdof(idof+4*NF+4:idof+5*NF+2)
     
     idof = idof + 6*NF+3

     call coilxyz( icoil )
     
    endif
    
    if( coil(icoil)%Ic == 0 ) then
     
     coil(icoil)%I        = xdof(idof+     1            )
     
     idof = idof + 1
     
    endif
    


   enddo ! end of do icoil; 04 Sep 17;

  case default
   
   FATAL( packdof, .true., selected packorunpack not supported )

  end select
  
!  lxdof(1:Ndof) = zero ! reset xdof;
!  
!  call packcoil !pack coil parameters into DoF;
!  
!  idof = 0
!  
!  do icoil = 1, Ncoils
!   
!   if( coil(icoil)%Ic /= 0 ) then 
!    lxdof(idof+1) = coil(icoil)%I / Inorm
!    idof = idof + 1
!   endif
!   
!   ND = DoF(icoil)%ND
!   
!   if( coil(icoil)%Lc /= 0 ) then
!    lxdof(idof+1:idof+ND) = DoF(icoil)%xdof(1:ND) / Gnorm
!    idof = idof + ND
!   endif
!   
!  enddo !end do icoil;
!  
!  FATAL( packdof , idof .ne. Ndof, counting error in packing )
!  
!  call mpi_barrier(MPI_COMM_WORLD, ierr)
  
  return
  
END SUBROUTINE packdof

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

!SUBROUTINE unpacking(lxdof)
!  
!  use globals, only : zero, myid, ounit, Ncoils, coil, DoF, Ndof, Inorm, Gnorm
!  
!  implicit none
!  
!  include "mpif.h"
!  
!  REAL    :: lxdof(1:Ndof)
!  INTEGER :: idof, icoil, ND, astat, ierr, ifirst
!  
!  idof = 0 ; ifirst = 0
!  
!  do icoil = 1, Ncoils
!   
!   if( coil(icoil)%Ic /= 0 ) then
!    coil(icoil)%I = lxdof(idof+1) * Inorm
!    idof = idof + 1
!   endif
!   
!   ND = DoF(icoil)%ND
!   
!   if( coil(icoil)%Lc /= 0 ) then
!    DoF(icoil)%xdof(1:ND) = lxdof(idof+1:idof+ND) * Gnorm
!    idof = idof + ND
!   endif
!   
!  enddo !end do icoil;
!  
!  FATAL( unpacking , idof .ne. Ndof, counting error in unpacking )
!  
!  call unpackcoil !unpack DoF to coil parameters;
!  
! !call discoil(ifirst)
!  
!  call mpi_barrier(MPI_COMM_WORLD, ierr)
!  
!  return
!  
!END SUBROUTINE unpacking

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!SUBROUTINE packcoil
!  !---------------------------------------------------------------------------------------------
!  ! pack coil representation variables into DoF (only geometries without currents);
!  ! DATE: 2017/03/25
!  !--------------------------------------------------------------------------------------------- 
!  use globals, only : zero, myid, ounit, Ncoils, coil, FouCoil, DoF
!  implicit none
!  include "mpif.h"
!
!  INTEGER  :: icoil, idof, NF, ierr, astat
!
!  FATAL( packcoil, .not. allocated(coil)   , illegal )
!  FATAL( packcoil, .not. allocated(FouCoil), illegal )
!  FATAL( packcoil, .not. allocated(DoF)    , illegal )
!
!  do icoil = 1, Ncoils
!
!     select case (coil(icoil)%itype)
!     !--------------------------------------------------------------------------------------------- 
!     case(1)
!        ! get number of DoF for each coil and allocate arrays;
!        NF = FouCoil(icoil)%NF
!        DoF(icoil)%xdof = zero
!
!        !pack Fourier series;
!        idof = 0
!        if(coil(icoil)%Lc /= 0) then
!           DoF(icoil)%xdof(idof+1:idof+NF+1) = FouCoil(icoil)%xc(0:NF) ; idof = idof + NF+1
!           DoF(icoil)%xdof(idof+1:idof+NF  ) = FouCoil(icoil)%xs(1:NF) ; idof = idof + NF    
!           DoF(icoil)%xdof(idof+1:idof+NF+1) = FouCoil(icoil)%yc(0:NF) ; idof = idof + NF+1
!           DoF(icoil)%xdof(idof+1:idof+NF  ) = FouCoil(icoil)%ys(1:NF) ; idof = idof + NF    
!           DoF(icoil)%xdof(idof+1:idof+NF+1) = FouCoil(icoil)%zc(0:NF) ; idof = idof + NF+1
!           DoF(icoil)%xdof(idof+1:idof+NF  ) = FouCoil(icoil)%zs(1:NF) ; idof = idof + NF    
!        endif
!        FATAL( packcoil , idof .ne. DoF(icoil)%ND, counting error in packing )
!        
!     !---------------------------------------------------------------------------------------------
!     case default
!        FATAL(packcoil, .true., not supported coil types)
!     end select     
!
!  enddo ! end do icoil;
!
!END SUBROUTINE packcoil

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!SUBROUTINE unpackcoil
!  
!  use globals, only : zero, myid, ounit, Ncoils, coil, FouCoil, DoF
!  
!  implicit none
!  
!  include "mpif.h"
!  
!  INTEGER  :: icoil, idof, NF, ierr, astat
!  
!  FATAL( unpackcoil, .not. allocated(coil)   , illegal )
!  FATAL( unpackcoil, .not. allocated(FouCoil), illegal )
!  FATAL( unpackcoil, .not. allocated(DoF)    , illegal )
!  
!  do icoil = 1, Ncoils
!   
!   select case (coil(icoil)%itype)
!    
!   case(1)
!    
!    NF = FouCoil(icoil)%NF
!    
!    idof = 0
!    
!    if( coil(icoil)%Lc /= 0 ) then
!     
!     FouCoil(icoil)%xc(0:NF) = DoF(icoil)%xdof(idof+1 : idof+NF+1) ; idof = idof + NF + 1
!     FouCoil(icoil)%xs(1:NF) = DoF(icoil)%xdof(idof+1 : idof+NF  ) ; idof = idof + NF    
!     FouCoil(icoil)%yc(0:NF) = DoF(icoil)%xdof(idof+1 : idof+NF+1) ; idof = idof + NF + 1
!     FouCoil(icoil)%ys(1:NF) = DoF(icoil)%xdof(idof+1 : idof+NF  ) ; idof = idof + NF    
!     FouCoil(icoil)%zc(0:NF) = DoF(icoil)%xdof(idof+1 : idof+NF+1) ; idof = idof + NF + 1
!     FouCoil(icoil)%zs(1:NF) = DoF(icoil)%xdof(idof+1 : idof+NF  ) ; idof = idof + NF   
!    endif
!
!    FATAL( packcoil , idof .ne. DoF(icoil)%ND, counting error in packing )
!    
!   case default
!   
!    FATAL(packcoil, .true., not supported coil types)
!    
!   end select
!   
!  enddo ! end do icoil;
!
!return
!
!END SUBROUTINE unpackcoil

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
