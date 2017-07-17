
!title (diagnostic) ! Poincar&eacute; plot.

!latex \briefly{\Poincare plot.}

!latex \calledby{\link{notopt}}

!latex \tableofcontents

!latex \subsection{overview}
!latex \bi
!latex \item[1.] First,
!latex \item[5.] If \inputvar{Lpoincare} $\ne 0$, then 
!latex \item[i.] a \Poincare plot is constructed using \oculus{pp00aa}.
!latex \ei

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine diagnos
  
  use kmodule, only : zero, half, pi, pi2, small, ext, ncpu, myid, ounit, nunit, lunit, nodelabel, &
                      Itopology, knotsurf, Nteta, Nzeta, surf, &
                      bNfp, &
                      Lpoincare, phi, Ppts, Ptrj, odetol, poincare
  
  use oculus , only : pp00aa
  
  implicit none
  
  include "mpif.h"
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOGICAL              :: exist
  INTEGER              :: ierr, astat, iostat, itrj, ipp00aa, Nfp, ii(1:1)
  REAL                 :: lRl, lRu, lZl, lZu
  REAL, allocatable    :: RZ(:,:)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( myid.eq.0 ) write(ounit,'("diagnos : " 10x " : Ppts ="i6" ; Ptrj ="i4" ; odetol ="es8.1" ;")') Ppts, Ptrj, odetol
  
  if( Ppts.le.0 ) return
  
  poincare%Nfp = bNfp ; poincare%phi = phi !poincare%phi = zero * (pi2/Nfp) !try to plot PP at input phi; czhu 2016/04/22
  
  open(nunit+myid, file="."//trim(ext)//".fo.P."//nodelabel//".dat", status="unknown", form="unformatted", iostat=iostat )
  FATAL( diagnos, iostat.ne.0, error opening .ext.fo.P.xxx )

  poincare%idirection = +1 ; poincare%Ppts = Ppts ; poincare%odetol = odetol ; poincare%iLyapunov = 0 ; poincare%Ltol = zero

  if( LPoincare.ge. 1 ) poincare%flparameter = 0
  if( LPoincare.eq.-1 ) poincare%flparameter = 1
  
 !ii(1:1) = minloc(surf(1)%xx(0:Nteta,0)) ; lRl = surf(1)%xx(ii(1)-1,0) ; lZl = surf(1)%zz(ii(1)-1,0) ! THIS ONLY WORKS AT PHI = ZERO; 14 Apr 16;
 !ii(1:1) = maxloc(surf(1)%xx(0:Nteta,0)) ; lRu = surf(1)%xx(ii(1)-1,0) ; lZu = surf(1)%zz(ii(1)-1,0) ! THIS ONLY WORKS AT PHI = ZERO; 14 Apr 16;

  lRl = surf(1)%xx(Nteta/2,0) ; lZl = zero ; !lZl = surf(1)%zz(      0,0)
  lRu = surf(1)%xx(Nteta  ,0) ; lZu = zero ; !lZu = surf(1)%zz(Nteta/2,0)

  lRl = ( lRl + lRu ) * half ! re-center to mid-point; 04 Aug 16;
  lZl = ( lZl + lZu ) * half     

  inquire( file = trim(ext)//".fo.starting.points", exist=exist )

  if( exist ) then
   if( myid.eq.0 ) open( lunit, file = trim(ext)//".fo.starting.points", status="old" )
   if( myid.eq.0 ) read( lunit, * ) Ptrj
   if( myid.eq.0 ) write(ounit,'("diagnos : " 10x " : reading ",i3," starting points from ext.fo.starting.points ;")') Ptrj
   IlBCAST( Ptrj, 1, 0 )
   if( Ptrj.le.0 ) then
    exist = .false.
   else
    SALLOCATE( RZ, (1:2,1:Ptrj), zero )
    if( myid.eq.0 ) read( lunit, * ) RZ(1:2,1:Ptrj)
    RlBCAST( RZ(1:2,1:Ptrj), 2*Ptrj, 0 )
   endif ! end of if( Ptrj.le.0 ) ; 04 Aug 16;
   if( myid.eq.0 ) close( lunit )
  else
   if( myid.eq.0 ) write(ounit,'("diagnos : " 10x " : ext.fo.starting.points does not exist ;")')
  endif ! end of if( exist ) ; 04 Aug 16;

  if( Ptrj.le.0 ) return

  do itrj = 1, Ptrj

   if( myid.ne.modulo(itrj-1,ncpu) ) cycle ! construction of Poincare plot is in parallel; 30 Oct 15;
    
   poincare%R = lRl + (itrj-1) * ( lRu - lRl ) / max(1,Ptrj-1)
   poincare%Z = lZl + (itrj-1) * ( lZu - lZl ) / max(1,Ptrj-1)
    
   if( exist ) then ; poincare%R = RZ(1,itrj) ! starting points have been read from file; 04 Aug 16;
    ;               ; poincare%Z = RZ(2,itrj)
   endif

   FATAL( diagnos, poincare%R.lt.small, illegal R )

   ipp00aa = 0 ; call pp00aa( poincare, ipp00aa ) ! problem happens here; czhu 2016/04/22
   write(ounit,1200) ipp00aa, myid, "called ", itrj, poincare%phi, poincare%R, poincare%Z, poincare%Ppts, poincare%odetol, poincare%ipts
   
   if( poincare%Lallocated.eq.1 ) then
    write(nunit+myid)                   poincare%ipts
    write(nunit+myid) poincare%RZ(1:2,0:poincare%ipts)
   endif
  enddo ! end of do itrj = 0, Ptrj-1 ; 25 Mar 15;
   
  close(nunit+myid)
  
1200 format("diagnos : ipp00aa="i2" : myid="i3" ; "a7" pp00aa ; itrj="i5" ; phi="f9.5" ; (R,Z)=("f9.5" ,"f9.5" ) ; Ppts="i7" ; odetol="es8.1" ;"i6" ;")
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  return
  
end subroutine diagnos

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine wtmgrid
  use kmodule, only : zero, half, pi2, ext, ncpu, myid, ounit, lunit, Lpoincare
  implicit none

  include "mpif.h"

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  LOGICAL              :: exist
  INTEGER              :: ierr, astat, iostat, imn, ibn00aa, ip, iz, ir, np, nz, nr, itangent, ibfield, Mfp, nextcur
  REAL                 :: RpZ(1:3), R, P, Z, Pmin, Pmax, Zmin, Zmax, Rmin, Rmax, Btmp(1:3,0:3), pressure, gap
  REAL, allocatable    :: BRZp(:,:,:,:), dBRZp(:,:,:,:), BRpZ(:,:,:,:), dBRpZ(:,:,:,:)
  CHARACTER*13         :: suffix
  CHARACTER*30         :: curlabel(1:1)

  np = 36 ; nz = 128 ; nr = 128 ; Mfp = 1 ! SHOULD BE USER INPUT; 04 Aug 16;

  SALLOCATE( BRZp, (1:3,1:Nr,1:Nz,1:Np), zero )
  SALLOCATE(dBRZp, (1:3,1:Nr,1:Nz,1:Np), zero )
  SALLOCATE( BRpZ, (1:3,1:Nr,1:Np,1:Nz), zero )
  SALLOCATE(dBRpZ, (1:3,1:Nr,1:Np,1:Nz), zero )

  Pmin = zero ; Pmax = pi2 ! DO NOT CHANGE; 04 Aug 16;

  call plasdim(Rmin, Rmax, Zmin, Zmax)  !calculate plasma surface boundary  ;09/11/2016
 !call coildim(Rmin, Rmax, Zmin, Zmax)

  gap = 0.05
  Rmin = Rmin - gap; Rmax = Rmax + gap
  Zmin = Zmin - gap; Zmax = Zmax + gap

 !Rmin =  2.4 ; Rmax =  3.6
 !Zmin = -0.6 ; Zmax =  0.6

  if( myid.eq.0 ) write( ounit,'("wtmgrid : " 10x " : writing mgrid or m3dc1 file at grid of [ "4(ES23.15,2X)" ]",3i6)') Rmin, Rmax, Zmin, Zmax, np, nr, nz

  do ip = 1, np ; RpZ(2) = Pmin + ( Pmax - Pmin ) * ( ip - 1 ) / ( np - 0 ) / Mfp

     if ( myid.ne.modulo(ip,ncpu) ) cycle

     do iz = 1, nz ; RpZ(3) = Zmin + ( Zmax - Zmin ) * ( iz - 1 ) / ( nz - 1 )

        do ir = 1, nr ; RpZ(1) = Rmin + ( Rmax - Rmin ) * ( ir - 1 ) / ( nr - 1 )

           itangent = 0 ; call bfield( RpZ(1:3), itangent, Btmp(1:3,0:3), ibfield ) 

           Btmp(2,0) = Btmp(2,0) * RpZ(1)

           dBRpZ(1:3,ir,ip,iz) = (/ Btmp(1,0), Btmp(2,0), Btmp(3,0) /)
           dBRZp(1:3,ir,iz,ip) = (/ Btmp(1,0), Btmp(3,0), Btmp(2,0) /)

        enddo

     enddo

  enddo

  call MPI_Reduce(dBRpZ, BRpZ, 3*nr*nz*np, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_Reduce(dBRZp, BRZp, 3*nr*nz*np, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

  if( myid.eq.0 .and. Lpoincare.ne.8 ) then

     nextcur = 1 ; curlabel(1) = "focus-space-curves"

     write(suffix,'(i3.3,".",i4.4,".",i4.4)') Np, Nr, Nz

     write( ounit,'("wtmgrid : " 10x " : writing mgrid.ext.f."i3.3"."i4.4"."i4.4" ; Mfp="i3" ;")') np, nr, nz, Mfp

    !open( lunit, file=trim(ext)//".fo.mgrid", status="unknown", form="unformatted", iostat=iostat )
     open( lunit, file="mgrid."//trim(ext)//".f."//suffix, status="unknown", form="unformatted", iostat=iostat )
     FATAL( wtmgrid, iostat.ne.0, error opening ext.fo.mgrid )
     write(lunit) Nr, Nz, Np, Mfp, nextcur
     write(lunit) Rmin, Zmin, Rmax, Zmax
     write(lunit) curlabel(1:nextcur)
     write(lunit) BRZp(1:3,1:Nr,1:Nz,1:Np)
     close(lunit)

  endif

  ! write 4m3dc1 only when Lpoicare >= 8
  if(myid .eq. 0 .and. Lpoincare .ge. 8) then
     write( ounit,'("wtgrid : " 10x " : writing ext.fo.4m3dc1 ; input for M3D-C^1 ; nr="i4" ; np="i4" ; nz="i4" ;")') nr, np, nz

     open( lunit, file=trim(ext)//".fo.4m3dc1", status="unknown", form="formatted", iostat=iostat )
     FATAL( exinput, iostat.ne.0, error opening ext.4m3d )

     write( lunit,'("# NR     Nphi    NZ")')
     write( lunit,'("# ",3i6)') nr, np, nz
     write( lunit,'("# R[m]                   phi[rad]           Z[m]               BR[T]              Bphi[T]            BZ[T]              Pres[Pa]")')
     
     pressure = zero
     do ip = 1, np ; RpZ(2) = Pmin + ( Pmax - Pmin ) * ( ip - 1 ) / ( np - 0 )
        do iz = 1, nz ; RpZ(3) = Zmin + ( Zmax - Zmin ) * ( iz - 1 ) / ( nz - 1 )
           do ir = 1, nr ; RpZ(1) = Rmin + ( Rmax - Rmin ) * ( ir - 1 ) / ( nr - 1 )
              write(lunit,'(7es23.15)') RpZ(1:3), BRpZ(1:3,ir, iz, ip), pressure
           enddo
        enddo
     enddo
  
     close( lunit )
  endif

  DEALLOCATE(dBRZp)
  DEALLOCATE(dBRpZ)
  DEALLOCATE( BRZp)
  DEALLOCATE( BRpZ)


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  return

end subroutine wtmgrid

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine plasdim(Rmin, Rmax, Zmin, Zmax)

  use kmodule, only: zero, thousand, nteta, nzeta, surf
  implicit none
  include "mpif.h"

  REAL    :: Rmin, Rmax, Zmin, Zmax

  INTEGER :: ii, jj
  REAL    :: R, Z

  Rmin=thousand; Rmax=zero
  Zmin=thousand; Zmax=zero

  do jj = 1, nzeta
     do ii = 1, nteta

        R = sqrt(surf(1)%xx(ii,jj)**2 + surf(1)%yy(ii,jj)**2)
        Z = surf(1)%zz(ii,jj)
        if( R .ge. Rmax) Rmax = R; if( R .le. Rmin) Rmin = R
        if( Z .ge. Zmax) Zmax = Z; if( Z .le. Zmin) Zmin = Z

     enddo
  enddo

  return
end subroutine plasdim

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!


subroutine coildim(Rmin, Rmax, Zmin, Zmax)

  use kmodule, only: zero, thousand, Ncoils, NDcoil, coil
  implicit none
  include "mpif.h"

  REAL    :: Rmin, Rmax, Zmin, Zmax

  INTEGER :: ii, nn
  REAL    :: R, Z


  call discretecoil

  Rmin=thousand; Rmax=zero
  Zmin=thousand; Zmax=zero

  do ii = 1, Ncoils
     do nn = 1, NDcoil
     
        R = sqrt(coil(ii)%xx(nn)**2 + coil(ii)%yy(nn)**2)
        Z = coil(ii)%zz(nn)
        if( R .ge. Rmax) Rmax = R; if( R .le. Rmin) Rmin = R
        if( Z .ge. Zmax) Zmax = Z; if( Z .le. Zmin) Zmin = Z

     enddo
  enddo

  return
end subroutine coildim
