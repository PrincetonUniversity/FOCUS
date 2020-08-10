!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

SUBROUTINE diagnos
!------------------------------------------------------------------------------------------------------ 
! DATE: 07/13/2017
! diagonose the coil performance
!------------------------------------------------------------------------------------------------------   
  use globals, only: dp, zero, one, myid, ounit, sqrtmachprec, IsQuiet, case_optimize, coil, surf, Ncoils, &
       Nteta, Nzeta, bnorm, bharm, tflux, ttlen, specw, ccsep, coilspace, FouCoil, iout, Tdof, case_length, &
       cssep, Bmnc, Bmns, tBmnc, tBmns, weight_bharm, coil_importance, weight_bnorm, overlap, &
       pmsum, total_moment, magtorque, ext, MPI_COMM_FAMUS
                     
  implicit none
  include "mpif.h"

  INTEGER           :: icoil, itmp=0, astat, ierr, NF, idof, i, j,iteta, jzeta
  LOGICAL           :: lwbnorm = .True. , l_raw = .False.!if use raw coils data
  REAL              :: MaxCurv, AvgLength, MinCCdist, MinCPdist, tmp_dist, ReDot, ImDot, B(3), x, y, z
  REAL, parameter   :: infmax = 1.0E6
  REAL, allocatable :: Atmp(:,:), Btmp(:,:)

  if (myid == 0 .and. IsQuiet < 0) write(ounit, *) "-----------COIL DIAGNOSTICS----------------------------------"

  !--------------------------------cost functions-------------------------------------------------------  
  if (case_optimize == 0) call AllocData(0) ! if not allocate data;
  call costfun(case_optimize)

  if (myid == 0) write(ounit, '("diagnos : "5(A12," ; "))') , &
       "Bnormal", "Bmn harmonics", "tor. flux", "coil length", "c-s sep." 
  if (myid == 0) write(ounit, '("        : "6(ES12.5," ; "))') bnorm, bharm, tflux, ttlen, cssep

  ! compute Bx, By, Bz for later calculations
  do jzeta = 0, Nzeta - 1
     do iteta = 0, Nteta - 1
        ! x = 1.49; y = 0 ; z=0; B=0
        x = surf(1)%xx(iteta, jzeta)
        y = surf(1)%yy(iteta, jzeta)
        z = surf(1)%zz(iteta, jzeta)
        B = 0
        call coils_bfield(B,x,y,z)
        surf(1)%Bx(iteta, jzeta) = surf(1)%Bx(iteta, jzeta) + B(1)
        surf(1)%By(iteta, jzeta) = surf(1)%By(iteta, jzeta) + B(2)
        surf(1)%Bz(iteta, jzeta) = surf(1)%Bz(iteta, jzeta) + B(2)
     enddo
  enddo
  

  !--------------------------------calculate the average Bn error-------------------------------
  if (allocated(surf(1)%bn)) then
     ! \sum{ |Bn| / |B| }/ (Nt*Nz)
     if(myid .eq. 0) then 
        write(ounit, '(8X": Ave. relative absolute Bn error |Bn|/B : " ES12.5"; max(|Bn|)="ES12.5)') &
          sum(abs(surf(1)%bn/sqrt(surf(1)%Bx**2 + surf(1)%By**2 + surf(1)%Bz**2))) / (Nteta*Nzeta), &
          maxval(abs(surf(1)%bn))
     endif
  endif

  !--------------------------------calculate dipole effective volume------------------------------------  
  call minvol(0)
  if(myid .eq. 0)  write(ounit, '(8X": Total free magnetic moment M (L-2 norm):", ES12.5, &
       "; Effective ratio=", ES12.5)') sqrt(total_moment), sqrt(pmsum)
  !--------------------------------------------------------------------------------------------- 

  if (magtorque) then
     if (myid==0) write(ounit, '(8X": Mag. torque removed!")')
  endif
  
  return

  !--------------------------------calculate coil importance------------------------------------  
  if (.not. allocated(coil_importance)) then
     SALLOCATE( coil_importance, (1:Ncoils), zero )
  endif

  if (weight_bnorm > sqrtmachprec .or. weight_bharm > sqrtmachprec) then  ! make sure data_allocated
     do icoil = 1, Ncoils
        call importance(icoil)
     enddo
  
     if(myid .eq. 0) write(ounit, '(8X": The most and least important coils are :  " & 
          ES12.5" at coil" I4 " ; " ES12.5" at coil "I4)')      &
      maxval(coil_importance), maxloc(coil_importance), &
      minval(coil_importance), minloc(coil_importance)

  endif

  if (myid == 0 .and. IsQuiet < 0) write(ounit, *) "--------------------------------------------"
 
  return

END SUBROUTINE diagnos

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine curvature(icoil)

  use globals, only: dp, zero, pi2, ncpu, astat, ierr, myid, ounit, coil, NFcoil, Nseg, Ncoils
  implicit none
  include "mpif.h"

  INTEGER, INTENT(in) :: icoil

  REAL,allocatable    :: curv(:)

  SALLOCATE(curv, (0:coil(icoil)%NS), zero)

  curv = sqrt( (coil(icoil)%za*coil(icoil)%yt-coil(icoil)%zt*coil(icoil)%ya)**2  &
             + (coil(icoil)%xa*coil(icoil)%zt-coil(icoil)%xt*coil(icoil)%za)**2  & 
             + (coil(icoil)%ya*coil(icoil)%xt-coil(icoil)%yt*coil(icoil)%xa)**2 )& 
             / ((coil(icoil)%xt)**2+(coil(icoil)%yt)**2+(coil(icoil)%zt)**2)**(1.5)
  coil(icoil)%maxcurv = maxval(curv)

  return
end subroutine curvature

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine mindist(array_A, dim_A, array_B, dim_B, minimum)
 
  use globals, only: dp
  implicit none

  INTEGER, INTENT(IN ) :: dim_A, dim_B
  REAL   , INTENT(IN ) :: array_A(1:3,1:dim_A), array_B(1:3,1:dim_B)
  REAL   , INTENT(OUT) :: minimum

  INTEGER :: ipoint, jpoint, itmp, jtmp
  REAL, parameter :: infmax = 1.0E6
  REAL    :: distance

  minimum = infmax
  do ipoint = 1, dim_A
     do jpoint = 1, dim_B

        distance = ( array_A(1, ipoint) - array_B(1, jpoint) )**2 &
                  +( array_A(2, ipoint) - array_B(2, jpoint) )**2 &
                  +( array_A(3, ipoint) - array_B(3, jpoint) )**2
        if (distance .le. minimum) minimum = distance
        
     enddo
  enddo

  minimum = sqrt(minimum)
  
  return

end subroutine mindist

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine importance(icoil)
  use globals, only: dp,  zero, pi2, ncpu, astat, ierr, myid, ounit, coil, NFcoil, Nseg, Ncoils, &
                     surf, Nteta, Nzeta, bsconstant, coil_importance, MPI_COMM_FAMUS
  implicit none
  include "mpif.h"

  INTEGER, INTENT(in) :: icoil  

  INTEGER                               :: iteta, jzeta, NumGrid
  REAL                                  :: dBx, dBy, dBz
  REAL, dimension(0:Nteta-1, 0:Nzeta-1) :: tbx, tby, tbz        ! summed Bx, By and Bz

  !--------------------------initialize and allocate arrays------------------------------------- 

  NumGrid = Nteta*Nzeta
  tbx = zero; tby = zero; tbz = zero        !already allocted; reset to zero;

  do jzeta = 0, Nzeta - 1
     do iteta = 0, Nteta - 1

        if( myid.ne.modulo(jzeta*Nteta+iteta,ncpu) ) cycle ! parallelization loop;
        call bfield0(icoil, surf(1)%xx(iteta, jzeta), surf(1)%yy(iteta, jzeta), &
                   & surf(1)%zz(iteta, jzeta), tbx(iteta, jzeta), tby(iteta, jzeta), tbz(iteta, jzeta))

     enddo ! end do iteta
  enddo ! end do jzeta

  call MPI_BARRIER( MPI_COMM_FAMUS, ierr )     
  call MPI_ALLREDUCE( MPI_IN_PLACE, tbx, NumGrid, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_FAMUS, ierr )
  call MPI_ALLREDUCE( MPI_IN_PLACE, tby, NumGrid, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_FAMUS, ierr )
  call MPI_ALLREDUCE( MPI_IN_PLACE, tbz, NumGrid, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_FAMUS, ierr )

  coil_importance(icoil) = sum( (tbx*surf(1)%Bx + tby*surf(1)%By + tbz*surf(1)%Bz) / &
                                (surf(1)%Bx**2 + surf(1)%By**2 + surf(1)%Bz**2) ) / NumGrid

  return

end subroutine importance
