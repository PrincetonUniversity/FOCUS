!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

SUBROUTINE diagnos
!------------------------------------------------------------------------------------------------------ 
! DATE: 07/13/2017
! diagonose the coil performance
!------------------------------------------------------------------------------------------------------   
  use focus_globals, only: dp, zero, one, myid, ounit, sqrtmachprec, IsQuiet, case_optimize, coil, surf, Ncoils, &
       Nteta, Nzeta, bnorm, bharm, tflux, ttlen, specw, ccsep, coilspace, FouCoil, iout, Tdof, case_length, &
       cssep, Bmnc, Bmns, tBmnc, tBmns, weight_bharm, coil_importance, Npc, weight_bnorm, overlap, &
       pmsum, total_moment, magtorque, ext
                     
  implicit none
  include "mpif.h"

  INTEGER           :: icoil, itmp=0, astat, ierr, NF, idof, i, j
  LOGICAL           :: lwbnorm = .True. , l_raw = .False.!if use raw coils data
  REAL              :: MaxCurv, AvgLength, MinCCdist, MinCPdist, tmp_dist, ReDot, ImDot, B(3), x, y, z
  REAL, parameter   :: infmax = 1.0E6
  REAL, allocatable :: Atmp(:,:), Btmp(:,:)

  if (myid == 0 .and. IsQuiet < 0) write(ounit, *) "-----------COIL DIAGNOSTICS----------------------------------"

  !--------------------------------cost functions-------------------------------------------------------  
  if (case_optimize == 0) call AllocData(0) ! if not allocate data;
  call costfun(0)

  if (myid == 0) write(ounit, '("diagnos : "5(A12," ; "))') , &
       "Bnormal", "Bmn harmonics", "tor. flux", "coil length", "c-s sep." 
  if (myid == 0) write(ounit, '("        : "6(ES12.5," ; "))') bnorm, bharm, tflux, ttlen, cssep

!!$  ! output origin magnetic field
!!$  x = 1.5_dp ; y = 0 ; z = 0 ; B = 0
!!$  call coils_bfield(B,x,y,z)
!!$  if (myid == 0) write(ounit, '("        : B("3(F4.2,", ")") = ("3(ES12.5", ")")")') x,y,z,B(1),B(2),B(3) 

  !save all the coil parameters;
  if (allocated(coilspace)) then
     idof = 0
     do icoil = 1, Ncoils
        coilspace(iout, idof+1 ) = coil(icoil)%I ;  idof = idof + 1

        select case (coil(icoil)%itype)
        case (1)
           NF = FouCoil(icoil)%NF
           coilspace(iout, idof+1:idof+NF+1) = FouCoil(icoil)%xc(0:NF) ; idof = idof + NF +1
           coilspace(iout, idof+1:idof+NF  ) = FouCoil(icoil)%xs(1:NF) ; idof = idof + NF
           coilspace(iout, idof+1:idof+NF+1) = FouCoil(icoil)%yc(0:NF) ; idof = idof + NF +1
           coilspace(iout, idof+1:idof+NF  ) = FouCoil(icoil)%ys(1:NF) ; idof = idof + NF
           coilspace(iout, idof+1:idof+NF+1) = FouCoil(icoil)%zc(0:NF) ; idof = idof + NF +1
           coilspace(iout, idof+1:idof+NF  ) = FouCoil(icoil)%zs(1:NF) ; idof = idof + NF
!!$        case default
!!$           FATAL(descent, .true., not supported coil types)
        end select
     enddo
!!$     FATAL( output , idof .ne. Tdof, counting error in restart )
  endif

  !-------------------------------coil maximum curvature----------------------------------------------------  
  MaxCurv = zero
  do icoil = 1, Ncoils
     if(coil(icoil)%itype .ne. 1) exit ! only for Fourier
     call curvature(icoil)
     if (coil(icoil)%maxcurv .ge. MaxCurv) then
        MaxCurv = coil(icoil)%maxcurv
        itmp = icoil !record the number
     endif
#ifdef DEBUG
     if(myid .eq. 0) write(ounit, '(8X": Maximum curvature of "I3 "-th coil is : " ES23.15)') &
          icoil, coil(icoil)%maxcurv
#endif
  enddo
  if(myid .eq. 0) write(ounit, '(8X": Maximum curvature of all the coils is  :" ES23.15 &
       " ; at coil " I3)') MaxCurv, itmp

  !-------------------------------average coil length-------------------------------------------------------  
  AvgLength = zero
  if ( (case_length == 1) .and. (sum(coil(1:Ncoils)%Lo) < sqrtmachprec) ) coil(1:Ncoils)%Lo = one
  call length(0)
  do icoil = 1, Ncoils
     if(coil(icoil)%itype .ne. 1) exit ! only for Fourier
     AvgLength = AvgLength + coil(icoil)%L
  enddo
  AvgLength = AvgLength / Ncoils
  if(myid .eq. 0) write(ounit, '(8X": Average length of the coils is"8X" :" ES23.15)') AvgLength

  !-----------------------------minimum coil coil separation------------------------------------  
  ! coils are supposed to be placed in order
  minCCdist = infmax
  do icoil = 1, Ncoils

     if(coil(icoil)%itype .ne. 1) exit ! only for Fourier

     if(Ncoils .eq. 1) exit !if only one coil
     itmp = icoil + 1
     if(icoil .eq. Ncoils) itmp = 1

     SALLOCATE(Atmp, (1:3,0:coil(icoil)%NS-1), zero)
     SALLOCATE(Btmp, (1:3,0:coil(itmp )%NS-1), zero)

     Atmp(1, 0:coil(icoil)%NS-1) = coil(icoil)%xx(0:coil(icoil)%NS-1)
     Atmp(2, 0:coil(icoil)%NS-1) = coil(icoil)%yy(0:coil(icoil)%NS-1)
     Atmp(3, 0:coil(icoil)%NS-1) = coil(icoil)%zz(0:coil(icoil)%NS-1)

     Btmp(1, 0:coil(itmp )%NS-1) = coil(itmp)%xx(0:coil(itmp )%NS-1)
     Btmp(2, 0:coil(itmp )%NS-1) = coil(itmp)%yy(0:coil(itmp )%NS-1)
     Btmp(3, 0:coil(itmp )%NS-1) = coil(itmp)%zz(0:coil(itmp )%NS-1)
     
     call mindist(Atmp, coil(icoil)%NS, Btmp, coil(itmp)%NS, tmp_dist)

     if (minCCdist .ge. tmp_dist) minCCdist=tmp_dist

     DALLOCATE(Atmp)
     DALLOCATE(Btmp)

  enddo

  if(myid .eq. 0) write(ounit, '(8X": The minimum coil-coil distance is "4X" :" ES23.15)') minCCdist

  !--------------------------------minimum coil plasma separation-------------------------------  
  minCPdist = infmax
  do icoil = 1, Ncoils

     if(coil(icoil)%itype .ne. 1) exit ! only for Fourier

     SALLOCATE(Atmp, (1:3,0:coil(icoil)%NS-1), zero)
     SALLOCATE(Btmp, (1:3,1:(Nteta*Nzeta)), zero)

     Atmp(1, 0:coil(icoil)%NS-1) = coil(icoil)%xx(0:coil(icoil)%NS-1)
     Atmp(2, 0:coil(icoil)%NS-1) = coil(icoil)%yy(0:coil(icoil)%NS-1)
     Atmp(3, 0:coil(icoil)%NS-1) = coil(icoil)%zz(0:coil(icoil)%NS-1)

     Btmp(1, 1:(Nteta*Nzeta)) = reshape(surf(1)%xx(0:Nteta-1, 0:Nzeta-1), (/Nteta*Nzeta/))
     Btmp(2, 1:(Nteta*Nzeta)) = reshape(surf(1)%yy(0:Nteta-1, 0:Nzeta-1), (/Nteta*Nzeta/))
     Btmp(3, 1:(Nteta*Nzeta)) = reshape(surf(1)%zz(0:Nteta-1, 0:Nzeta-1), (/Nteta*Nzeta/))

     call mindist(Atmp, coil(icoil)%NS, Btmp, Nteta*Nzeta, tmp_dist)

     if (minCPdist .ge. tmp_dist) then 
        minCPdist=tmp_dist
        itmp = icoil
     endif

     DALLOCATE(Atmp)
     DALLOCATE(Btmp)

  enddo
  if(myid .eq. 0) then
     write(ounit, '(8X": The minimum coil-plasma distance is    :" ES23.15 &
          " ; at coil " I3)') minCPdist, itmp
     if (minCPdist < 0.1) then
        write(ounit, '(8X":-----------WARNING!!!-------------------------")')
        write(ounit, '(8X": Initial coils might intersect with the plasma")')
        write(ounit, '(8X":----------------------------------------------")') 
     endif
  endif
  
  !--------------------------------overlap of Bn_mn harmonics-----------------------------------  
  ReDot = zero ; ImDot = zero
  if (weight_bharm > sqrtmachprec) then  ! do not care n,m; use cos - i sin
     ReDot = sum(Bmnc*tBmnc) + sum(Bmns*tBmns)
     ImDot = sum(Bmnc*tBmns) - sum(Bmns*tBmnc)
     overlap = sqrt( (ReDot*ReDot + ImDot*ImDot) &
                 / ( (sum(Bmnc*Bmnc)+sum(Bmns*Bmns)) * (sum(tBmnc*tBmnc)+sum(tBmns*tBmns)) ) )
     if(myid .eq. 0) write(ounit, '(8X": Overlap of the realized Bn harmonics is:" F8.3 " %")') 100*overlap
  endif

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
  if(myid .eq. 0)  write(ounit, '(8X": Total free magnetic moment M           :", ES12.5, &
       "; Effective ratio=", ES12.5)') total_moment, pmsum
  !--------------------------------------------------------------------------------------------- 

  if (magtorque) then
     call mag_torque()
     if (myid==0) write(ounit, '(8X": magnetic field and moment data are saved in "A)') trim(ext)//".mag"
  endif
  
  return

  !--------------------------------calculate coil importance------------------------------------  
  if (.not. allocated(coil_importance)) then
     SALLOCATE( coil_importance, (1:Ncoils*Npc), zero )
  endif

  if (weight_bnorm > sqrtmachprec .or. weight_bharm > sqrtmachprec) then  ! make sure data_allocated
     do icoil = 1, Ncoils*Npc
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

  use focus_globals, only: dp, zero, pi2, ncpu, astat, ierr, myid, ounit, coil, NFcoil, Nseg, Ncoils
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
 
  use focus_globals, only: dp
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
  use focus_globals, only: dp,  zero, pi2, ncpu, astat, ierr, myid, ounit, coil, NFcoil, Nseg, Ncoils, &
                     surf, Nteta, Nzeta, bsconstant, coil_importance
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

  call MPI_BARRIER( MPI_COMM_WORLD, ierr )     
  call MPI_ALLREDUCE( MPI_IN_PLACE, tbx, NumGrid, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
  call MPI_ALLREDUCE( MPI_IN_PLACE, tby, NumGrid, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
  call MPI_ALLREDUCE( MPI_IN_PLACE, tbz, NumGrid, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )

  coil_importance(icoil) = sum( (tbx*surf(1)%Bx + tby*surf(1)%By + tbz*surf(1)%Bz) / &
                                (surf(1)%Bx**2 + surf(1)%By**2 + surf(1)%Bz**2) ) / NumGrid

  return

end subroutine importance
