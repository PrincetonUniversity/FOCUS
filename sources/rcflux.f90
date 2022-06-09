
!title (bnormal) ! Calculate total bnormal and its derivatives.

!latex \briefly{Calculate the total bnormal of all coils on plasma surface and the derivatives with respect to coil geometry and currents, including the first and second dirivatives. 
!latex          Calling \emph{bnormal(0), bnormal(1), bnormal(2)} calculates the $0-order$, $1^{st}-order$ and $2^{nd}-order$ derivatives respectively.}

!latex \calledby{\link{costfun}}
!latex \calls{\link{bfield}}

!latex \tableofcontents

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

module bnorm_mod
  ! contains some common variables used in subroutine bnormal
  ! allocating once and re-using them will save allocation time
  use globals, only : dp
  implicit none

  ! 0-order
  REAL, allocatable :: dBx(:,:), dBy(:,:), dBz(:,:), Bm(:,:), dAx(:,:), dAy(:,:), dAz(:,:)
  ! 1st-order
  REAL, allocatable :: dBn(:), dBm(:), d1B(:,:,:)

end module bnorm_mod
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine rcflux( ideriv )
!------------------------------------------------------------------------------------------------------ 
! DATE:  04/02/2017;
! Calculate the Bn surface integral and its derivatives;
! ideriv = 0 -> only calculate the Bn surface integral;
! ideriv = 1 -> calculate the Bn surface integral and its first derivatives;
! ideriv = 2 -> calculate the Bn surface integral and its first & second derivatives;
!------------------------------------------------------------------------------------------------------   
  use globals, only: dp, zero, half, one, pi2, sqrtmachprec, bsconstant, ncpu, myid, ounit, &
       coil, DoF, surf, Ncoils, Nteta, Nzeta, discretefactor, plasma, &
       bnorm, t1B, t2B, bn, Ndof, Cdof, weight_bharm, case_bnormal, &
       weight_bnorm, ibnorm, mbnorm, ibharm, mbharm, LM_fvec, LM_fjac, &
       bharm, t1H, Bmnc, Bmns, wBmn, tBmnc, tBmns, Bmnim, Bmnin, NBmn, &
       weight_resbn, target_resbn, resbn, resbn_m, resbn_n, t1R, b1s, b1c, resbn_bnc, resbn_bns, &
       gsurf, ghost_use, ghost_call, ghost_once, rcflux_use, machprec, MPI_COMM_FOCUS
  use bnorm_mod
  use bharm_mod
  use mpi
  implicit none

  INTEGER, INTENT(in)                   :: ideriv
  !--------------------------------------------------------------------------------------------
  INTEGER                               :: astat, ierr, suc
  INTEGER                               :: icoil, iteta, jzeta, idof, ND, NumGrid, isurf
  REAL                                  :: arg, teta, zeta, bnc, bns, shift, rcflux, psmall, &
                                           small, negvalue, posvalue, tmp_xdof(1:Ndof)
  REAL, allocatable                     :: cosarg(:,:), sinarg(:,:)
  !--------------------------initialize and allocate arrays-------------------------------------
  isurf = plasma 
  NumGrid = Nteta*Nzeta
  ! reset to zero;
  bnorm = zero 
  surf(isurf)%Bx = zero; surf(isurf)%By = zero; surf(isurf)%Bz = zero; surf(isurf)%Bn = zero     
  dBx = zero; dBy = zero; dBz = zero; Bm = zero; dAx = zero; dAy = zero; dAz = zero
  bn = zero
  suc = 1
  if (weight_resbn .gt. sqrtmachprec) then
      FATAL( bnormal, resbn_m .le. 0, wrong poloidal mode number)
      FATAL( bnormal, resbn_n .le. 0, wrong toroidal mode number)
      resbn = zero ; bnc = zero ;  bns = zero
      shift = half
      if (ghost_use .eq. 1 .and. ghost_call .eq. 1) then
         call ghost(1, suc)
         if ( suc .eq. 0 ) then
            resbn = 10.0
            ! return
         endif
         if ( ghost_once .eq. 1 ) then
            ghost_call = 0
         endif
      endif
  endif
  !-------------------------------calculate Bn-------------------------------------------------- 
  if( ideriv >= 0 ) then
     if (weight_resbn .gt. sqrtmachprec) then ! resonant Bn perturbation
         if ( suc .eq. 1 ) then ! No indent
         rcflux = 0.0
         do jzeta = 1, gsurf(1)%Nseg_stable-1
            if( myid.ne.modulo(jzeta,ncpu) ) cycle ! parallelization loop;
            do icoil = 1, Ncoils
               call afield0(icoil, gsurf(1)%ox(jzeta), gsurf(1)%oy(jzeta), gsurf(1)%oz(jzeta), dAx(0,0), dAy(0,0), dAz(0,0) )
               rcflux = rcflux + dAx(0,0)*gsurf(1)%oxdot(jzeta)
               rcflux = rcflux + dAy(0,0)*gsurf(1)%oydot(jzeta)
               rcflux = rcflux + dAz(0,0)*gsurf(1)%ozdot(jzeta)
               call afield0(icoil, gsurf(1)%xx(jzeta), gsurf(1)%xy(jzeta), gsurf(1)%xz(jzeta), dAx(0,0), dAy(0,0), dAz(0,0) )
               rcflux = rcflux - dAx(0,0)*gsurf(1)%xxdot(jzeta)
               rcflux = rcflux - dAy(0,0)*gsurf(1)%xydot(jzeta)
               rcflux = rcflux - dAz(0,0)*gsurf(1)%xzdot(jzeta)
            enddo
         enddo
         call MPI_BARRIER( MPI_COMM_FOCUS, ierr )
         call MPI_ALLREDUCE( MPI_IN_PLACE, rcflux, 1  , MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FOCUS, ierr )
         rcflux = rcflux*pi2*resbn_m/(gsurf(1)%Nseg_stable-1)
         resbn = rcflux**2
         endif
         ! Include negative for maximizing islands
     endif
  endif

  !-------------------------------calculate Bn/x------------------------------------------------
  if ( ideriv >= 1 ) then
     ! reset data
     if (weight_resbn .gt. sqrtmachprec) then
         t1R = zero
     endif
     if (weight_resbn .gt. sqrtmachprec) then
        do jzeta = 1, gsurf(1)%Nseg_stable-1
           if( myid.ne.modulo(jzeta,ncpu) ) cycle ! parallelization loop;
           idof = 0
           do icoil = 1, Ncoils
              ND = DoF(icoil)%ND
              if ( coil(icoil)%Ic /= 0 ) then
                 call afield0(icoil, gsurf(1)%ox(jzeta), gsurf(1)%oy(jzeta), gsurf(1)%oz(jzeta), dAx(0,0), dAy(0,0), dAz(0,0) )
                 t1R(idof+1) = t1R(idof+1) + dAx(0,0)*gsurf(1)%oxdot(jzeta) / coil(icoil)%I
                 t1R(idof+1) = t1R(idof+1) + dAy(0,0)*gsurf(1)%oydot(jzeta) / coil(icoil)%I
                 t1R(idof+1) = t1R(idof+1) + dAz(0,0)*gsurf(1)%ozdot(jzeta) / coil(icoil)%I
                 call afield0(icoil, gsurf(1)%xx(jzeta), gsurf(1)%xy(jzeta), gsurf(1)%xz(jzeta), dAx(0,0), dAy(0,0), dAz(0,0) )
                 t1R(idof+1) = t1R(idof+1) - dAx(0,0)*gsurf(1)%xxdot(jzeta) / coil(icoil)%I
                 t1R(idof+1) = t1R(idof+1) - dAy(0,0)*gsurf(1)%xydot(jzeta) / coil(icoil)%I
                 t1R(idof+1) = t1R(idof+1) - dAz(0,0)*gsurf(1)%xzdot(jzeta) / coil(icoil)%I
                 idof = idof + 1
              endif
              if ( coil(icoil)%Lc /= 0 ) then
                 call afield1(icoil, gsurf(1)%ox(jzeta), gsurf(1)%oy(jzeta), gsurf(1)%oz(jzeta), dAx(1:ND,0), dAy(1:ND,0), dAz(1:ND,0), ND )
                 t1R(idof+1:idof+ND) = t1R(idof+1:idof+ND) + dAx(1:ND,0)*gsurf(1)%oxdot(jzeta)
                 t1R(idof+1:idof+ND) = t1R(idof+1:idof+ND) + dAy(1:ND,0)*gsurf(1)%oydot(jzeta)
                 t1R(idof+1:idof+ND) = t1R(idof+1:idof+ND) + dAz(1:ND,0)*gsurf(1)%ozdot(jzeta)
                 call afield1(icoil, gsurf(1)%xx(jzeta), gsurf(1)%xy(jzeta), gsurf(1)%xz(jzeta), dAx(1:ND,0), dAy(1:ND,0), dAz(1:ND,0), ND )
                 t1R(idof+1:idof+ND) = t1R(idof+1:idof+ND) - dAx(1:ND,0)*gsurf(1)%xxdot(jzeta)
                 t1R(idof+1:idof+ND) = t1R(idof+1:idof+ND) - dAy(1:ND,0)*gsurf(1)%xydot(jzeta)
                 t1R(idof+1:idof+ND) = t1R(idof+1:idof+ND) - dAz(1:ND,0)*gsurf(1)%xzdot(jzeta)
                 idof = idof + ND
              endif
           enddo
           FATAL( bnormal, idof .ne. Ndof, counting error in packing ) 
        enddo
        call MPI_BARRIER( MPI_COMM_FOCUS, ierr )
        call MPI_ALLREDUCE( MPI_IN_PLACE, t1R, Ndof  , MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FOCUS, ierr )
        t1R(1:Ndof) = t1R(1:Ndof)*pi2*resbn_m/(gsurf(1)%Nseg_stable-1)
        t1R(1:Ndof) = 2*rcflux*t1R(1:Ndof)
     endif
  endif
  !--------------------------------------------------------------------------------------------
  call MPI_barrier( MPI_COMM_FOCUS, ierr )
  return
end subroutine rcflux
