!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

SUBROUTINE diagnos
!------------------------------------------------------------------------------------------------------ 
! DATE: 07/13/2017
! diagonose the coil performance
!------------------------------------------------------------------------------------------------------   
  use globals, only: dp, zero, one, myid, ounit, sqrtmachprec, IsQuiet, case_optimize, coil, surf, Ncoils, &
       Nteta, Nzeta, bnorm, bharm, tflux, ttlen, specw, ccsep, coilspace, FouCoil, iout, Tdof, case_length, &
       cssep, Bmnc, Bmns, tBmnc, tBmns, weight_bharm, coil_importance, Nfp, weight_bnorm, overlap, plasma, &
       cosnfp, sinnfp, symmetry, discretefactor, MPI_COMM_FOCUS, surf_Nfp, curv, case_curv, tors, nis, &
       lambda_alpha, weight_nis
  use mpi
  implicit none

  INTEGER           :: icoil, itmp, astat, ierr, NF, idof, i, j, isurf, cs, ip, is, Npc, coilInd0, coilInd1
  LOGICAL           :: lwbnorm, l_raw
  REAL              :: MaxCurv, AvgLength, MinCCdist, MinCPdist, tmp_dist, ReDot, ImDot, dum, AvgCurv, AvgTors, &
                       MinLambda, MaxS, torsRet
  REAL, parameter   :: infmax = 1.0E6
  REAL, allocatable :: Atmp(:,:), Btmp(:,:)

  isurf = plasma
  itmp = 0
  lwbnorm = .True. 
  l_raw = .False. ! if use raw coils data
  if (myid == 0 .and. IsQuiet < 0) write(ounit, *) "-----------COIL DIAGNOSTICS----------------------------------"

  !--------------------------------cost functions-------------------------------------------------------  
  if (case_optimize == 0) call AllocData(0) ! if not allocate data;
  call costfun(0)

  !if (myid == 0) write(ounit, '("diagnos : "8(A12," ; "))') , &
  !     "Bnormal", "Bmn harmonics", "tor. flux", "coil length", "c-s sep." , "curvature", "c-c sep.", "torsion"
  !if (myid == 0) write(ounit, '("        : "8(ES12.5," ; "))') bnorm, bharm, tflux, ttlen, cssep, curv, ccsep, tors
  if (myid == 0) write(ounit, '("diagnos : "9(A12," ; "))') , &
       "Bnormal", "Bmn harmonics", "tor. flux", "coil length", "c-s sep." , "curvature", "c-c sep.", "torsion", "nissin"
  if (myid == 0) write(ounit, '("        : "9(ES12.5," ; "))') bnorm, bharm, tflux, ttlen, cssep, curv, ccsep, tors, nis

  !save all the coil parameters;
  if (allocated(coilspace)) then
     idof = 0
     do icoil = 1, Ncoils
        coilspace(iout, idof+1 ) = coil(icoil)%I ;  idof = idof + 1

        select case (coil(icoil)%type)
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
     if(coil(icoil)%type .ne. 1) exit ! only for Fourier
     call CurvDeriv0(icoil,dum) !dummy return
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

  !-------------------------------average coil curvature-------------------------------------------------------
  AvgCurv = zero
  do icoil = 1, Ncoils
     if(coil(icoil)%type .ne. 1) exit ! only for Fourier
     call avgcurvature(icoil)
     AvgCurv = AvgCurv + coil(icoil)%avgcurv
#ifdef DEBUG
     if(myid .eq. 0) write(ounit, '(8X": Average curvature of "I3 "-th coil is : " ES23.15)') &
          icoil, coil(icoil)%avgcurv
#endif
  enddo
  AvgCurv = AvgCurv / Ncoils
  if(myid .eq. 0) write(ounit, '(8X": Average curvature of the coils is"5X" :" ES23.15)') AvgCurv

  !-------------------------------average coil torsion-----------------------------------------------------
  AvgTors = zero
  do icoil = 1, Ncoils
     if(coil(icoil)%type .ne. 1) exit ! only for Fourier 
     call TorsDeriv0(icoil,torsRet,dum) !dummy return
     AvgTors = AvgTors + abs(torsRet) 
#ifdef DEBUG
     if(myid .eq. 0) write(ounit, '(8X": Average torsion of "I3 "-th coil is   : " ES23.15)') &
        icoil, abs(torsRet)
#endif
  enddo
  AvgTors = AvgTors / Ncoils
  if(myid .eq. 0) write(ounit, '(8X": Average torsion of the coils is"5X"   :" ES23.15)') AvgTors

  !-------------------------------minimum coil lambda-----------------------------------------------------
  if( lambda_alpha .ne. 0 ) then
  call TorsDeriv0(1,dum,dum)
  MinLambda = coil(1)%minlambda
  do icoil = 1, Ncoils
     if(coil(icoil)%type .ne. 1) exit ! only for Fourier
     call TorsDeriv0(icoil,dum,dum) !dummy return
     if( coil(icoil)%minlambda .le. MinLambda) then
        MinLambda = coil(icoil)%minlambda
     endif
#ifdef DEBUG
     if(myid .eq. 0) write(ounit, '(8X": Minimum lambda of "I3 "-th coil is    : " ES23.15)') &
        icoil, coil(icoil)%minlambda
#endif
  enddo
  if(myid .eq. 0) write(ounit, '(8X": Minimum lambda of the coils is"5X"    :" ES23.15)') MinLambda
  endif

  !-------------------------------maximum coil S------------------------------------------------------------
  if( weight_nis .ne. 0 ) then
  call NisDeriv0(1,dum)
  MaxS = coil(1)%maxs
  do icoil = 1, Ncoils
     if(coil(icoil)%type .ne. 1) exit ! only for Fourier
     call NisDeriv0(icoil,dum) !dummy return
     if( coil(icoil)%maxs .ge. MaxS) then
        MaxS = coil(icoil)%maxs
     endif
#ifdef DEBUG
     if(myid .eq. 0) write(ounit, '(8X": Maximum S of "I3 "-th coil is         : " ES23.15)') &
        icoil, coil(icoil)%maxs
#endif
  enddo
  if(myid .eq. 0) write(ounit, '(8X": Maximum S of the coils is"5X"         :" ES23.15)') MaxS
  endif

  !-------------------------------average coil length-------------------------------------------------------  
  AvgLength = zero
  if ( (case_length == 1) .and. (sum(coil(1:Ncoils)%Lo) < sqrtmachprec) ) coil(1:Ncoils)%Lo = one
  call length(0)
  do icoil = 1, Ncoils
     if(coil(icoil)%type .ne. 1) exit ! only for Fourier
     AvgLength = AvgLength + coil(icoil)%L
  enddo
  AvgLength = AvgLength / Ncoils
  if(myid .eq. 0) write(ounit, '(8X": Average length of the coils is"8X" :" ES23.15)') AvgLength

  !-----------------------------minimum coil coil separation------------------------------------  
  ! coils are supposed to be placed in order
  minCCdist = infmax
  coilInd0 = 0
  coilInd1 = 1
  do icoil = 1, Ncoils
     if(coil(icoil)%type .ne. 1) exit ! only for Fourier
     if(Ncoils .eq. 1) exit ! if only one coil
     ! Data for the first coil
     SALLOCATE(Atmp, (1:3,0:coil(icoil)%NS-1), zero)
     Atmp(1, 0:coil(icoil)%NS-1) = coil(icoil)%xx(0:coil(icoil)%NS-1)
     Atmp(2, 0:coil(icoil)%NS-1) = coil(icoil)%yy(0:coil(icoil)%NS-1)
     Atmp(3, 0:coil(icoil)%NS-1) = coil(icoil)%zz(0:coil(icoil)%NS-1)
     ! Check distances for coil and its symmetric coils 
     if ( coil(icoil)%symm /= 0 ) then
        SALLOCATE(Btmp, (1:3,0:coil(icoil)%NS-1), zero)
        if ( coil(icoil)%symm .eq. 1 ) then
           cs = 0
           Npc = Nfp
        elseif ( coil(icoil)%symm .eq. 2 ) then
           cs = 1
           Npc = Nfp
        endif
        do ip = 1, Npc
           do is = 0, cs
              if ( ip .eq. 1 .and. is .eq. 0 ) cycle
              Btmp(1, 0:coil(icoil)%NS-1) = (coil(icoil)%xx(0:coil(icoil)%NS-1)*cosnfp(ip) &
                                        & - coil(icoil)%yy(0:coil(icoil)%NS-1)*sinnfp(ip) ) 
              Btmp(2, 0:coil(icoil)%NS-1) = (-1)**is * (coil(icoil)%xx(0:coil(icoil)%NS-1)*sinnfp(ip) &
                                                   & + coil(icoil)%yy(0:coil(icoil)%NS-1)*cosnfp(ip) )
              Btmp(3, 0:coil(icoil)%NS-1) = (-1)**is * (coil(icoil)%zz(0:coil(icoil)%NS-1))
              call mindist(Atmp, coil(icoil)%NS, Btmp, coil(icoil)%NS, tmp_dist)
#ifdef DEBUG
              if(myid .eq. 0) write(ounit, '(8X": distance between  "I3.3"-th and "I3.3"-th coil (ip="I2.2 &
                   ", is="I1") is : " ES23.15)') icoil, icoil, ip, is, tmp_dist
#endif
              if (minCCdist .ge. tmp_dist) then
                 minCCdist = tmp_dist
                 coilInd0 = icoil
                 coilInd1 = icoil
              endif
           enddo
        enddo
        DALLOCATE(Btmp)
     endif
     do itmp = 1, Ncoils
        ! skip self and non-Fourier coils
        if (itmp == icoil .or. coil(icoil)%type /= 1) cycle
        SALLOCATE(Btmp, (1:3,0:coil(itmp)%NS-1), zero)
        ! check if the coil is stellarator symmetric
        !select case (coil(icoil)%symm)
        select case (coil(itmp)%symm) 
        case ( 0 )
           cs  = 0
           Npc = 1
        case ( 1 )
           cs  = 0
           Npc = Nfp
        case ( 2 ) 
           cs  = 1
           Npc = Nfp
        end select
        ! periodicity and stellarator symmetry
        do ip = 1, Npc
           do is = 0, cs
              Btmp(1, 0:coil(itmp)%NS-1) =            (coil(itmp)%xx(0:coil(itmp)%NS-1)*cosnfp(ip) &
                                                   & - coil(itmp)%yy(0:coil(itmp)%NS-1)*sinnfp(ip) )
              Btmp(2, 0:coil(itmp)%NS-1) = (-1)**is * (coil(itmp)%xx(0:coil(itmp)%NS-1)*sinnfp(ip) &
                                                   & + coil(itmp)%yy(0:coil(itmp)%NS-1)*cosnfp(ip) )
              Btmp(3, 0:coil(itmp)%NS-1) = (-1)**is * (coil(itmp)%zz(0:coil(itmp)%NS-1))
              call mindist(Atmp, coil(icoil)%NS, Btmp, coil(itmp)%NS, tmp_dist)
#ifdef DEBUG
              if(myid .eq. 0) write(ounit, '(8X": distance between  "I3.3"-th and "I3.3"-th coil (ip="I2.2 & 
                   ", is="I1") is : " ES23.15)') icoil, itmp, ip, is, tmp_dist
#endif
              if (minCCdist .ge. tmp_dist) then 
                 minCCdist = tmp_dist
                 coilInd0 = icoil
                 coilInd1 = itmp
              endif
           enddo
        enddo
        DALLOCATE(Btmp)
     enddo
     DALLOCATE(Atmp)
  enddo

  if(myid .eq. 0) write(ounit, '(8X": The minimum coil-coil distance is "4X" :" ES23.15" ; at coils"I3" &
          ,"I3)') minCCdist, coilInd0, coilInd1

  !--------------------------------minimum coil plasma separation-------------------------------  
  minCPdist = infmax
  do icoil = 1, Ncoils

     if(coil(icoil)%type .ne. 1) exit ! only for Fourier

     SALLOCATE(Atmp, (1:3,0:coil(icoil)%NS-1), zero)
     SALLOCATE(Btmp, (1:3,1:(Nteta*Nzeta)), zero)

     Atmp(1, 0:coil(icoil)%NS-1) = coil(icoil)%xx(0:coil(icoil)%NS-1)
     Atmp(2, 0:coil(icoil)%NS-1) = coil(icoil)%yy(0:coil(icoil)%NS-1)
     Atmp(3, 0:coil(icoil)%NS-1) = coil(icoil)%zz(0:coil(icoil)%NS-1)

     Btmp(1, 1:(Nteta*Nzeta)) = reshape(surf(isurf)%xx(0:Nteta-1, 0:Nzeta-1), (/Nteta*Nzeta/))
     Btmp(2, 1:(Nteta*Nzeta)) = reshape(surf(isurf)%yy(0:Nteta-1, 0:Nzeta-1), (/Nteta*Nzeta/))
     Btmp(3, 1:(Nteta*Nzeta)) = reshape(surf(isurf)%zz(0:Nteta-1, 0:Nzeta-1), (/Nteta*Nzeta/))

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
  if (allocated(surf(isurf)%bn)) then
     ! \sum{ |Bn| / |B| }/ (Nt*Nz)
     if(myid .eq. 0) then
        write(ounit, '(8X": Ave. relative absolute Bn error |Bn|/B : " ES12.5"; max(|Bn|)="ES12.5)') &
             sum(abs(surf(plasma)%bn/sqrt(surf(plasma)%Bx**2+surf(plasma)%By**2+surf(plasma)%Bz**2))) &
             / (Nteta*Nzeta), maxval(abs(surf(plasma)%bn))
        write(ounit, '(8X": Surface area normalized Bn error int(|Bn|/B*ds)/A : "ES23.15)') &
             sum(abs(surf(plasma)%bn)/sqrt(surf(plasma)%Bx**2+surf(plasma)%By**2+surf(plasma)%Bz**2) &
             *surf(plasma)%ds)*discretefactor/(surf(plasma)%area/(surf_Nfp*2**symmetry))
     endif
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
  !--------------------------------------------------------------------------------------------- 

  if (myid == 0 .and. IsQuiet < 0) write(ounit, *) "--------------------------------------------"
 
  return

END SUBROUTINE diagnos

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine avgcurvature(icoil)

  use globals, only: dp, zero, pi2, ncpu, astat, ierr, myid, ounit, coil, NFcoil, Nseg, Ncoils
  implicit none
  include "mpif.h"

  INTEGER, INTENT(in) :: icoil

  REAL,allocatable    :: davgcurv(:)

  SALLOCATE(davgcurv, (0:coil(icoil)%NS-1), zero)

  davgcurv = sqrt( (coil(icoil)%za*coil(icoil)%yt-coil(icoil)%zt*coil(icoil)%ya)**2  &
             + (coil(icoil)%xa*coil(icoil)%zt-coil(icoil)%xt*coil(icoil)%za)**2  & 
             + (coil(icoil)%ya*coil(icoil)%xt-coil(icoil)%yt*coil(icoil)%xa)**2 )& 
             / ((coil(icoil)%xt)**2+(coil(icoil)%yt)**2+(coil(icoil)%zt)**2)**(1.5)
  davgcurv = davgcurv*sqrt(coil(icoil)%xt**2+coil(icoil)%yt**2+coil(icoil)%zt**2)
  coil(icoil)%avgcurv = pi2*sum(davgcurv)/size(davgcurv)
  coil(icoil)%avgcurv = coil(icoil)%avgcurv/coil(icoil)%L

  return
end subroutine avgcurvature

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
                     surf, Nteta, Nzeta, bsconstant, coil_importance, plasma, MPI_COMM_FOCUS
  implicit none
  include "mpif.h"

  INTEGER, INTENT(in) :: icoil  

  INTEGER                               :: iteta, jzeta, NumGrid, isurf
  REAL                                  :: dBx, dBy, dBz
  REAL, dimension(0:Nteta-1, 0:Nzeta-1) :: tbx, tby, tbz        ! summed Bx, By and Bz

  !--------------------------initialize and allocate arrays------------------------------------- 

  isurf = plasma
  NumGrid = Nteta*Nzeta
  tbx = zero; tby = zero; tbz = zero        !already allocted; reset to zero;

  do jzeta = 0, Nzeta - 1
     do iteta = 0, Nteta - 1

        if( myid.ne.modulo(jzeta*Nteta+iteta,ncpu) ) cycle ! parallelization loop;
        call bfield0(icoil, surf(isurf)%xx(iteta, jzeta), surf(isurf)%yy(iteta, jzeta), &
                   & surf(isurf)%zz(iteta, jzeta), tbx(iteta, jzeta), tby(iteta, jzeta), tbz(iteta, jzeta))

     enddo ! end do iteta
  enddo ! end do jzeta

  call MPI_BARRIER( MPI_COMM_FOCUS, ierr )     
  call MPI_ALLREDUCE( MPI_IN_PLACE, tbx, NumGrid, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_FOCUS, ierr )
  call MPI_ALLREDUCE( MPI_IN_PLACE, tby, NumGrid, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_FOCUS, ierr )
  call MPI_ALLREDUCE( MPI_IN_PLACE, tbz, NumGrid, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_FOCUS, ierr )

  coil_importance(icoil) = sum( (tbx*surf(isurf)%Bx + tby*surf(isurf)%By + tbz*surf(isurf)%Bz) / &
                                (surf(isurf)%Bx**2 + surf(isurf)%By**2 + surf(isurf)%Bz**2) ) / NumGrid

  return

end subroutine importance
