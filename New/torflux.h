!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

!title (toroidal flux) ! Compute toroidal flux.

!latex \briefly{The toroidal magnetic flux . . }

!latex \calledby{\link{solvers}}

!latex \section{Outline}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

subroutine torflux( ideriv )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  use globals, only : zero, half, one, pi2, sqrtmachprec, ncpu, myid, ounit, &
                      bsconstant, &
                      coil, &
                      surf, Ncoils, Nteta, Nzeta, &
                      toroidalfluxerror, toroidalfluxaverage, target_tflux, &
!                     t1F, t2F, &
                      deltatheta
  
  implicit none
  
  include "mpif.h"
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  INTEGER, INTENT(in)                :: ideriv
  
  INTEGER                            :: astat, ierr, icoil, iteta, jzeta
  REAL                               :: toroidalflux, summedtoroidalflux, summedtoroidalfluxerror
  REAL                               :: lax, lay, laz
  REAL                               :: dAx, dAy, dAz
! REAL, dimension(1:Ndof, 0:Nzeta-1) :: ldF, dF
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  summedtoroidalflux      = zero
  summedtoroidalfluxerror = zero
  
  do jzeta = 0, Nzeta-1
   
   if( myid.ne.modulo(jzeta,ncpu) ) cycle
   
   toroidalflux = zero
   
   do iteta = 0, Nteta-1
    
    lax = zero ; lay = zero ; laz = zero
    
    do icoil = 1, Ncoils
     
    !call abfield( icoil, iteta, jzeta, dAx, dAy, dAz )
     
     lax = lax + dAx * coil(icoil)%I
     lay = lay + dAy * coil(icoil)%I
     laz = laz + dAz * coil(icoil)%I
     
    enddo ! end of do icoil;
    
    toroidalflux = toroidalflux + lax*surf%xt(iteta,jzeta) + lay*surf%yt(iteta,jzeta) + laz*surf%zt(iteta,jzeta)
    
   enddo ! end of do iteta;
   
   toroidalflux = toroidalflux * deltatheta * bsconstant ! normalize; 04 Sep 17;
   
   summedtoroidalflux      = summedtoroidalflux      +   toroidalflux
   summedtoroidalfluxerror = summedtoroidalfluxerror + ( toroidalflux-target_tflux )**2
   
  enddo ! end of do jzeta;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  call MPI_BARRIER( MPI_COMM_WORLD, ierr )
  
  call MPI_REDUCE( summedtoroidalflux     , toroidalfluxaverage, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
  call MPI_REDUCE( summedtoroidalfluxerror, toroidalfluxerror  , 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
  
  RlBCAST( toroidalfluxerror   , 1    , 0)
  RlBCAST( toroidalfluxaverage , 1    , 0)
  
  CHECK( torflux, Nzeta.eq.0, divide by zero)

  toroidalfluxerror   = half * toroidalfluxerror   / Nzeta 
  toroidalfluxaverage =        toroidalfluxaverage / Nzeta
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  if ( ideriv >= 1 ) then
   
   FATAL( torflux, .true., under construction )
   
!     ldF = zero; dF = zero; t1F = zero
!
!     do jzeta = 0, Nzeta - 1
!        if( myid.ne.modulo(jzeta,ncpu) ) cycle ! parallelization loop;
!
!        do iteta = 0, Nteta - 1             
!           
!          !do ip = 1, Npc
!              idof = 0
!              do icoil = 1, Ncoils
!                 ND = DoF(icoil)%ND
!                 if ( coil(icoil)%Ic /= 0 ) then !if current is free;
!                    call bpotential0(icoil, iteta, jzeta, &
!                         & dAx(0,0), dAy(0,0), dAz(0,0))
!
!                    ldF(idof+1, jzeta) = ldF(idof+1, jzeta) &
!                         & + bsconstant * ( dAx(0,0)*surf%xt(iteta,jzeta)   &
!                         &                + dAy(0,0)*surf%yt(iteta,jzeta)   &
!                         &                + dAz(0,0)*surf%zt(iteta,jzeta) )
!                    idof = idof +1
!                 endif
!
!                 if ( coil(icoil)%Lc /= 0 ) then !if geometry is free;
!                    call bpotential1(icoil, iteta, jzeta, &
!                         &       dAx(1:ND,0), dAy(1:ND,0), dAz(1:ND,0), ND)
!
!                    ldF(idof+1:idof+ND, jzeta) = ldF(idof+1:idof+ND, jzeta) &
!                         & + bsconstant * coil(icoil)%I * ( dAx(1:ND,0)*surf%xt(iteta,jzeta)   &
!                         &                                + dAy(1:ND,0)*surf%yt(iteta,jzeta)   &
!                         &                                + dAz(1:ND,0)*surf%zt(iteta,jzeta) )
!
!                    idof = idof + ND
!                 endif
!
!              enddo !end icoil;
!              FATAL( torflux , idof .ne. Ndof, counting error in packing )
!          !enddo  ! end do ip;
!
!        enddo !end iteta;
!     enddo !end jzeta
!
!     ldF = ldF * pi2/Nteta
!
!     call MPI_BARRIER( MPI_COMM_WORLD, ierr )
!     call MPI_REDUCE(ldF, dF, Ndof*Nzeta, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
!     RlBCAST( dF, Ndof*Nzeta, 0 )
!
!     do idof = 1, Ndof
!        t1F(idof) = sum( psi_diff(0:Nzeta-1) * dF(idof, 0:Nzeta-1) ) / Nzeta
!     enddo
     
  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

  call MPI_barrier( MPI_COMM_WORLD, ierr )

  return

end subroutine torflux

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

!subroutine afield( icoil, iteta, jzeta, NS, dAdx )
!  
!  use globals, only : zero, one, sqrtmachprec, myid, myid, ounit, &
!                      surf, Nteta, Nzeta, &
!                      Ncoils, coil, Nseg, deltacurveparameter
!  
!  implicit none
!  
!  include "mpif.h"
!  
!  INTEGER, intent(in ) :: icoil, iteta, jzeta, NS
!  REAL   , intent(out) :: dAdx(0:NS,1:3,0:3)
!  
!!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
!  
!  INTEGER :: ierr, astat, ii
!  REAL    :: dx, dy, dz, rr(1:5), invr
!  
!!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
!  
!  CHECK( afield, icoil.lt.1 .or. icoil.gt.Ncoils, icoil not in range )
!  CHECK( afield, iteta.lt.0 .or. iteta.gt.Nteta , iteta not in range )
!  CHECK( afield, jzeta.lt.0 .or. jzeta.gt.Nzeta , jzeta not in range )
!  
!!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
!  
!  dAdx(1:NS,1:3,0:3) = zero ! perhaps this is not required; 04 Sep 17;
!
!  do ii = 1, Nseg
!   
!   dx = surf%xx(iteta,jzeta) - coil(icoil)%xx(ii)
!   dy = surf%yy(iteta,jzeta) - coil(icoil)%yy(ii)
!   dz = surf%zz(iteta,jzeta) - coil(icoil)%zz(ii)
!   
!   rr(2) = dx**2 + dy**2 + dz**2 ; rr(1) = sqrt(rr(2)) ! rr(3) = rr(1) * rr(2) ; rr(5) = rr(3) * rr(2)
!   
!   CHECK( afield , rr(1).lt.sqrtmachprec, divide by zero )
!   
!   invr = one / rr(1)
!   
!   dAdx(ii,1,0) = coil(icoil)%xt(ii) * invr ! * deltacurveparameter
!   dAdx(ii,2,0) = coil(icoil)%yt(ii) * invr
!   dAdx(ii,3,0) = coil(icoil)%zt(ii) * invr
!   
!  enddo ! end do ii;
!
!  dAdx(0,1,0) = sum( dAdx(1:NS,1,0) ) ! * deltacurveparameter
!  dAdx(0,2,0) = sum( dAdx(1:NS,2,0) ) ! * deltacurveparameter
!  dAdx(0,3,0) = sum( dAdx(1:NS,3,0) ) ! * deltacurveparameter
!
!!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
!  
!  return
!
!!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
!
!end subroutine afield

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!subroutine bpotential1(icoil, iteta, jzeta, Ax, Ay, Az, ND)
!
!  use globals, only: coil, surf, NFcoil, Ncoils, Nteta, Nzeta,      &
!                     zero, myid, ounit
!  implicit none
!  include "mpif.h"
!
!  INTEGER, intent(in ) :: icoil, iteta, jzeta, ND
!  REAL, dimension(1:1, 1:ND), intent(inout) :: Ax, Ay, Az
!

!
!  INTEGER              :: ierr, astat, kseg, NS
!  REAL                 :: dlx, dly, dlz, r, rm3, ltx, lty, ltz
!  REAL, dimension(1:1, 1:coil(icoil)%NS)   :: dAxx, dAxy, dAxz, dAyx, dAyy, dAyz, dAzx, dAzy, dAzz
!

!  
!  FATAL( bpotential1, icoil .lt. 1 .or. icoil .gt. Ncoils    , icoil not in right range )
!  FATAL( bpotential1, iteta .lt. 0 .or. iteta .gt. Nteta     , iteta not in right range )
!  FATAL( bpotential1, jzeta .lt. 0 .or. jzeta .gt. Nzeta     , jzeta not in right range )
!  FATAL( bpotential1, ND <= 0, wrong inout dimension of ND )
!  
!  NS = coil(icoil)%NS
!
!  dlx = zero; ltx = zero; Ax = zero
!  dly = zero; lty = zero; Ay = zero
!  dlz = zero; ltz = zero; Az = zero
!
!  do kseg = 1, NS
!     
!     dlx = surf%xx(iteta,jzeta) - coil(icoil)%xx(kseg)
!     dly = surf%yy(iteta,jzeta) - coil(icoil)%yy(kseg)
!     dlz = surf%zz(iteta,jzeta) - coil(icoil)%zz(kseg)
!
!     r = sqrt(dlx**2 + dly**2 + dlz**2); rm3 = r**(-3)
!
!     ltx = coil(icoil)%xt(kseg)
!     lty = coil(icoil)%yt(kseg)
!     ltz = coil(icoil)%zt(kseg)
!
!     dAxx(1,kseg) = - (dly*lty + dlz*ltz) * rm3 * coil(icoil)%dd(kseg) !Ax/x
!     dAxy(1,kseg) =    dly*ltx            * rm3 * coil(icoil)%dd(kseg) !Ax/y
!     dAxz(1,kseg) =    dlz*ltx            * rm3 * coil(icoil)%dd(kseg) !Ax/z
!
!     dAyx(1,kseg) =    dlx*lty            * rm3 * coil(icoil)%dd(kseg) !Ay/x
!     dAyy(1,kseg) = - (dlx*ltx + dlz*ltz) * rm3 * coil(icoil)%dd(kseg) !Ay/y
!     dAyz(1,kseg) =    dlz*lty            * rm3 * coil(icoil)%dd(kseg) !Ay/z
!
!     dAzx(1,kseg) =    dlx*ltz            * rm3 * coil(icoil)%dd(kseg) !Az/x
!     dAzy(1,kseg) =    dly*ltz            * rm3 * coil(icoil)%dd(kseg) !Az/y
!     dAzz(1,kseg) = - (dlx*ltx + dly*lty) * rm3 * coil(icoil)%dd(kseg) !Az/z
!
!  enddo    ! enddo kseg
!
!  Ax(1:1, 1:ND) = matmul(dAxx, DoF(icoil)%xof) + matmul(dAxy, DoF(icoil)%yof) + matmul(dAxz, DoF(icoil)%zof)
!  Ay(1:1, 1:ND) = matmul(dAyx, DoF(icoil)%xof) + matmul(dAyy, DoF(icoil)%yof) + matmul(dAyz, DoF(icoil)%zof)
!  Az(1:1, 1:ND) = matmul(dAzx, DoF(icoil)%xof) + matmul(dAzy, DoF(icoil)%yof) + matmul(dAzz, DoF(icoil)%zof)
!
!  return
!
!end subroutine bpotential1

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
