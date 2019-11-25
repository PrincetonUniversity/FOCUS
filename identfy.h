
!title (identify) ! Poincar&eacute; plot.

!latex \briefly{Calculate actual Bnormal value from coils file; Analyze coils properties.}

!latex \calledby{\link{notopt}}

!latex \tableofcontents

!latex \subsection{overview}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine identfy
  use kmodule, only: ounit, myid, zero, coil, surf, Nteta, Nzeta, coilsI, antibscont, Ncoils, &
       bnorm, tflux, ttlen, tbn, target_tflux, isign, Lnormalize, qasym
  implicit none
  include "mpif.h"

  INTEGER           :: icoil, itmp, astat, ierr, mf=20, nf=20
  LOGICAL           :: lwbnorm = .True. , l_raw = .False.!if use raw coils data
  REAL              :: MaxCurv, AvgLength, MinCCdist, MinCPdist, tmp_dist
  REAL, parameter   :: infmax = 1.0E6
  REAL, allocatable :: Atmp(:,:), Btmp(:,:)

  if( .not. allocated(coil) ) l_raw = .True.

  if(l_raw) then
     allocate( coil(1:Ncoils) )

     call coilsprocess

     do icoil = 1, Ncoils
        coil(icoil)%I  =  coilsI(icoil) * antibscont
     enddo

  endif

  ! evaluate object functions
#ifdef NORM
    call bnormal(0)
#else
    call bnormal2(0)
#endif

   if ( target_tflux .eq. 0.0 ) then
    call torflux(0)
    target_tflux = isign*sqrt(2.0*tflux)             
    if(myid .eq. 0) write(ounit,'("Costfun :"11X" : Reset target toroidal flux to"ES23.15)') target_tflux
   endif

  call torflux(0)
  call tlengthExp(0)
  if(myid .eq. 0) write(ounit, '("evaluate: "10X" : bnorm = "ES23.15" ; tflux = "ES23.15" ; ttlen = "ES23.15" ; ")') bnorm, tflux, ttlen
  if (myid==0) write(ounit, '("evaluate: "10X" : Quasi-axisymmetric functional = "ES23.15)'), qasym

  !if(l_raw) deallocate( coil ) 
  !return        !for avoiding evaluating curvature for straight lines; 2017/02/10
  
  ! maximum curvature of all the coils
  MaxCurv = zero
  do icoil = 1, Ncoils
     call curvature(icoil)
     if (coil(icoil)%maxcurv .ge. MaxCurv) then
        MaxCurv = coil(icoil)%maxcurv
        itmp = icoil !record the number
     endif
#ifdef DEBUG
     if(myid .eq. 0) write(ounit, '("evaluate: "10X" : Maximum curvature of "I3 "-th coil is : " ES23.15)') icoil, coil(icoil)%maxcurv
#endif
  enddo
  if(myid .eq. 0) write(ounit, '("evaluate: "10X" : Maximum curvature of all the coils is  :" ES23.15 " ; at coil " I3)') MaxCurv, itmp
  
  ! calculate the average length
  AvgLength = zero
  do icoil = 1, Ncoils
     AvgLength = AvgLength + coil(icoil)%L
  enddo
  AvgLength = AvgLength / Ncoils
  if(myid .eq. 0) write(ounit, '("evaluate: "10X" : Average length of the coils is"8X" :" ES23.15)') AvgLength

  ! calculate the minimum distance of coil-coil separation
  ! coils are supporsed to be placed in order
  minCCdist = infmax
  do icoil = 1, Ncoils

     if(Ncoils .eq. 1) exit !if only one coil
     itmp = icoil + 1
     if(icoil .eq. Ncoils) itmp = 1

     SALLOCATE(Atmp, (1:3,1:coil(icoil)%D), zero)
     SALLOCATE(Btmp, (1:3,1:coil(itmp )%D), zero)

     Atmp(1, 1:coil(icoil)%D) = coil(icoil)%xx(1:coil(icoil)%D)
     Atmp(2, 1:coil(icoil)%D) = coil(icoil)%yy(1:coil(icoil)%D)
     Atmp(3, 1:coil(icoil)%D) = coil(icoil)%zz(1:coil(icoil)%D)

     Btmp(1, 1:coil(itmp)%D) = coil(itmp)%xx(1:coil(itmp)%D)
     Btmp(2, 1:coil(itmp)%D) = coil(itmp)%yy(1:coil(itmp)%D)
     Btmp(3, 1:coil(itmp)%D) = coil(itmp)%zz(1:coil(itmp)%D)
     
     call mindist(Atmp, coil(icoil)%D, Btmp, coil(itmp)%D, tmp_dist)

     if (minCCdist .ge. tmp_dist) minCCdist=tmp_dist

     DALLOCATE(Atmp)
     DALLOCATE(Btmp)

  enddo

  if(myid .eq. 0) write(ounit, '("evaluate: "10X" : The minimum coil-coil distance is"5X" :" ES23.15)') minCCdist

  ! calculate the minimum distance of coil-plasma separation
  minCPdist = infmax
  do icoil = 1, Ncoils

     SALLOCATE(Atmp, (1:3,1:coil(icoil)%D), zero)
     SALLOCATE(Btmp, (1:3,1:(Nteta*Nzeta)), zero)

     Atmp(1, 1:coil(icoil)%D) = coil(icoil)%xx(1:coil(icoil)%D)
     Atmp(2, 1:coil(icoil)%D) = coil(icoil)%yy(1:coil(icoil)%D)
     Atmp(3, 1:coil(icoil)%D) = coil(icoil)%zz(1:coil(icoil)%D)

     Btmp(1, 1:(Nteta*Nzeta)) = reshape(surf(1)%xx(0:Nteta-1, 0:Nzeta-1), (/Nteta*Nzeta/))
     Btmp(2, 1:(Nteta*Nzeta)) = reshape(surf(1)%yy(0:Nteta-1, 0:Nzeta-1), (/Nteta*Nzeta/))
     Btmp(3, 1:(Nteta*Nzeta)) = reshape(surf(1)%zz(0:Nteta-1, 0:Nzeta-1), (/Nteta*Nzeta/))

     call mindist(Atmp, coil(icoil)%D, Btmp, Nteta*Nzeta, tmp_dist)

     if (minCPdist .ge. tmp_dist) then 
        minCPdist=tmp_dist
        itmp = icoil
     endif

     DALLOCATE(Atmp)
     DALLOCATE(Btmp)

  enddo
  if(myid .eq. 0) write(ounit, '("evaluate: "10X" : The minimum coil-plasma distance is    :" ES23.15 " ; at coil " I3)') minCPdist, itmp

! if(myid .eq. 0) then
!  FATAL( identfy, 8.gt.Nteta .or. 8.gt.Nzeta, illegal )
!  write(*,*) surf(1)%bnt(0,0), surf(1)%bnt(8,8)
! endif

  !if(l_raw) surf(1)%bnt = -tbn ! for substracting the bnorm error from the existing coils

  if(l_raw) deallocate( coil )

  !TMPOUT(surf(1)%bnt(0,0))
  !TMPOUT(tbn(0,0))

  !TMPOUT(surf(1)%bnt(Nteta/2, Nzeta/2))

  if(lwbnorm) call bnftran(mf, nf) ! need debug ; 2017/01/11

end subroutine identfy

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
subroutine coilsprocess
!directly use the data in coil.XX file. not so accurate, better use Fourier series first.

  use kmodule, only: zero, pi2, &
       coil, cmt, smt, NFcoil, NDcoil, Ncoils, coilsX, coilsY, coilsZ, coilsI, Nseg, Lw, &
       ncpu, myid, ounit
  implicit none
  include "mpif.h"

  INTEGER           :: astat, ierr
  INTEGER           :: icoil, iseg
  REAL              :: dlength
  REAL,allocatable  :: tmpx(:), tmpy(:), tmpz(:)

  icoil = 0; iseg = 0
  do icoil = 1, Ncoils

     coil(icoil)%D = Nseg(icoil)
     coil(icoil)%Lw = Lw

     SALLOCATE( coil(icoil)%xx, (0:coil(icoil)%D), zero )
     SALLOCATE( coil(icoil)%yy, (0:coil(icoil)%D), zero )
     SALLOCATE( coil(icoil)%zz, (0:coil(icoil)%D), zero )
     SALLOCATE( coil(icoil)%xt, (0:coil(icoil)%D), zero )
     SALLOCATE( coil(icoil)%yt, (0:coil(icoil)%D), zero )
     SALLOCATE( coil(icoil)%zt, (0:coil(icoil)%D), zero )

     SALLOCATE(           tmpx, (0:coil(icoil)%D), zero ) 
     SALLOCATE(           tmpy, (0:coil(icoil)%D), zero )
     SALLOCATE(           tmpz, (0:coil(icoil)%D), zero )

     coil(icoil)%xx(0:coil(icoil)%D-1) = coilsX(1:Nseg(icoil), icoil); coil(icoil)%xx(coil(icoil)%D) = coil(icoil)%xx(0)
     coil(icoil)%yy(0:coil(icoil)%D-1) = coilsY(1:Nseg(icoil), icoil); coil(icoil)%yy(coil(icoil)%D) = coil(icoil)%yy(0)
     coil(icoil)%zz(0:coil(icoil)%D-1) = coilsZ(1:Nseg(icoil), icoil); coil(icoil)%zz(coil(icoil)%D) = coil(icoil)%zz(0)

     dlength = pi2/(Nseg(icoil)) ! step size
     tmpx = cshift(coil(icoil)%xx, 1) !shift 1 dimension
     tmpy = cshift(coil(icoil)%yy, 1) !shift 1 dimension
     tmpz = cshift(coil(icoil)%zz, 1) !shift 1 dimension

     coil(icoil)%xt = (tmpx-coil(icoil)%xx)/dlength  ![f(x+h) - f(x)] / h
     coil(icoil)%yt = (tmpy-coil(icoil)%yy)/dlength  ![f(x+h) - f(x)] / h
     coil(icoil)%zt = (tmpz-coil(icoil)%zz)/dlength  ![f(x+h) - f(x)] / h

     coil(icoil)%xt(Nseg(icoil)) = coil(icoil)%xt(0)
     coil(icoil)%yt(Nseg(icoil)) = coil(icoil)%yt(0)
     coil(icoil)%zt(Nseg(icoil)) = coil(icoil)%zt(0)

     DALLOCATE(tmpx)
     DALLOCATE(tmpy)
     DALLOCATE(tmpz)
  enddo

end subroutine coilsprocess
 
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine curvature(icoil)

  use kmodule, only : zero, pi2, &
       coil, NFcoil, NDcoil, Ncoils,  &
       ncpu, myid, ounit
  implicit none
  include "mpif.h"

  INTEGER           :: icoil

  
  INTEGER           :: astat, ierr
  REAL,allocatable  :: curv(:)

  SALLOCATE(curv, (0:coil(icoil)%D), zero)

  curv = sqrt( (coil(icoil)%za*coil(icoil)%yt-coil(icoil)%zt*coil(icoil)%ya)**2  &
             + (coil(icoil)%xa*coil(icoil)%zt-coil(icoil)%xt*coil(icoil)%za)**2  & 
             + (coil(icoil)%ya*coil(icoil)%xt-coil(icoil)%yt*coil(icoil)%xa)**2 )& 
             / ((coil(icoil)%xt)**2+(coil(icoil)%yt)**2+(coil(icoil)%zt)**2)**(1.5)
  coil(icoil)%maxcurv = maxval(curv)

  return
end subroutine curvature

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine mindist(array_A, dim_A, array_B, dim_B, minimum)

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

  
