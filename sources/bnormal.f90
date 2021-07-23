
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
  REAL, allocatable :: dBx(:,:), dBy(:,:), dBz(:,:), Bm(:,:)
  REAL, allocatable :: gx(:,:,:), gy(:,:,:), gz(:,:,:) ! inductance matrix
  ! 1st-order
  REAL, allocatable :: dBn(:), dBm(:), d1B(:,:,:)
  REAL, allocatable :: b1c(:), b1s(:)

end module bnorm_mod
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine bnormal( ideriv )
!------------------------------------------------------------------------------------------------------ 
! DATE:  10/28/2019;
! Calculate the Bn surface integral and its derivatives;
! ideriv = 0 -> only calculate the Bn surface integral;
! ideriv = 1 -> calculate the Bn surface integral and its first derivatives;
!------------------------------------------------------------------------------------------------------   
  use globals, only: dp, zero, half, one, two, pi2, sqrtmachprec, bsconstant, ncpu, myid, ounit, &
       coil, DoF, surf, Ncoils, Nteta, Nzeta, discretefactor, cosnfp, sinnfp, &
       bnorm, t1B, t2B, bn, Ndof, Cdof, weight_bharm, case_bnormal, &
       weight_bnorm, ibnorm, mbnorm, ibharm, mbharm, LM_fvec, LM_fjac, &
       bharm, t1H, Bmnc, Bmns, wBmn, tBmnc, tBmns, Bmnim, Bmnin, NBmn, dof_offset, ldof, momentq, &
       MPI_COMM_FAMUS, &
       weight_resbn, target_resbn, resbn, resbn_m, resbn_n, T1R, shift, resbn_bnc, resbn_bns
  use bnorm_mod
  use bharm_mod
  use mpi
  implicit none

  INTEGER, INTENT(in)                   :: ideriv
  !--------------------------------------------------------------------------------------------
  INTEGER                               :: astat, ierr
  INTEGER                               :: icoil, iteta, jzeta, idof, ND, NumGrid, ip, is, index
  REAL                                  :: sinp, sint, cosp, cost, drho, dmx, dmy, dmz
  REAL                                  :: arg, teta, zeta, bnc, bns
  REAL   , allocatable                  :: cosarg(:,:), sinarg(:,:)

  !--------------------------initialize and allocate arrays------------------------------------- 

  NumGrid = Nteta*Nzeta
  ! reset to zero;
  bnorm = zero 
  surf(1)%Bn = zero     
  dBx = zero; dBy = zero; dBz = zero; Bm = zero
  bn = zero

  ! resonant Bn
   if (weight_resbn .gt. sqrtmachprec) then
      FATAL( bnormal, resbn_m .le. 0, wrong poloidal mode number)
      FATAL( bnormal, resbn_n .le. 0, wrong toroidal mode number)
      resbn = zero ; bnc = zero ;  bns = zero 
      if (.not. allocated(cosarg) ) then
         SALLOCATE( cosarg, (0:Nteta-1, 0:Nzeta-1), zero )
         SALLOCATE( sinarg, (0:Nteta-1, 0:Nzeta-1), zero )
         do jzeta = 0, Nzeta - 1
            zeta = ( jzeta + shift ) * pi2 / surf(1)%Nzeta
            do iteta = 0, Nteta - 1
               teta = ( iteta + shift ) * pi2 / surf(1)%Nteta
               arg = resbn_m*teta - resbn_n*zeta
               cosarg(iteta, jzeta) = cos(arg)
               sinarg(iteta, jzeta) = sin(arg)
            enddo
         enddo
      endif
   endif

  !-------------------------------calculate Bn-------------------------------------------------- 
  if( ideriv >= 0 ) then

     do jzeta = 0, Nzeta - 1
        do iteta = 0, Nteta - 1
           
           if (myid==0) then ! contribution from the master cpu and plasma currents
              surf(1)%Bn(iteta, jzeta) = surf(1)%Bx(iteta, jzeta)*surf(1)%nx(iteta, jzeta) &
                   &                   + surf(1)%By(iteta, jzeta)*surf(1)%ny(iteta, jzeta) &
                   &                   + surf(1)%Bz(iteta, jzeta)*surf(1)%nz(iteta, jzeta) &
                   &                   + surf(1)%pb(iteta, jzeta)
           else ! contribution from dipoles
              surf(1)%Bn(iteta, jzeta) = 1.0/surf(1)%ds(iteta, jzeta)*(sum(gx(iteta, jzeta, 1:Ncoils) * coil(1:Ncoils)%mx) &
                   & + sum(gy(iteta, jzeta, 1:Ncoils) * coil(1:Ncoils)%my) + sum(gz(iteta, jzeta, 1:Ncoils) * coil(1:Ncoils)%mz))
           endif

           ! gather all the data
           call MPI_ALLREDUCE( MPI_IN_PLACE, surf(1)%Bn(iteta, jzeta), 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FAMUS, ierr )
           
           select case (case_bnormal)
           case (0)     ! no normalization over |B|;
              bnorm = bnorm + surf(1)%Bn(iteta, jzeta) * surf(1)%Bn(iteta, jzeta) * surf(1)%ds(iteta, jzeta)
           case default
              FATAL( bnorm, .true., case_bnormal can only be 0 )
           end select

        enddo ! end do iteta
     enddo ! end do jzeta
     
     bnorm = bnorm * half * discretefactor
     bn = surf(1)%Bn -  surf(1)%pb  ! bn is B.n from coils

     ! resonant Bn
      if (weight_resbn .gt. sqrtmachprec) then ! resonant Bn perturbation
         bnc = sum(surf(1)%Bn * cosarg) * discretefactor
         bns = sum(surf(1)%Bn * sinarg) * discretefactor
         resbn_bnc = bnc
         resbn_bns = bns
         resbn = - abs(sqrt(bnc*bnc + bns*bns) - target_resbn)
         ! if(myid==0) write(ounit, '("Resonant Bmn spectrum = "ES23.15)'), sqrt(bnc*bnc + bns*bns)
      endif  

  endif
  
  !-------------------------------calculate Bn/x------------------------------------------------
  if ( ideriv >= 1 ) then

     t1B = zero 
     if (mbnorm > 0 .or. weight_bharm > sqrtmachprec)  d1B = zero
     dBn = zero ; dBm = zero

      if (weight_resbn .gt. sqrtmachprec) then
         t1R = zero ; b1c = zero ; b1s = zero
      endif

     do jzeta = 0, Nzeta - 1
        do iteta = 0, Nteta - 1           
           idof = dof_offset
           dBn = zero
           if (myid==0) then
              continue
           else
              do icoil = 1, Ncoils
                 cost = cos(coil(icoil)%mt) ; sint = sin(coil(icoil)%mt)
                 cosp = cos(coil(icoil)%mp) ; sinp = sin(coil(icoil)%mp) 
                 if ( coil(icoil)%Ic /= 0 ) then !if current is free;
                    ! dBn/drho
                    drho = coil(icoil)%moment*momentq*(coil(icoil)%pho)**(momentq-1)
                    dmx = drho*sint*cosp
                    dmy = drho*sint*sinp
                    dmz = drho*cost
                    dBn(idof+1) = ( gx(iteta, jzeta, icoil) * dmx &
                         &        + gy(iteta, jzeta, icoil) * dmy &
                         &        + gz(iteta, jzeta, icoil) * dmz ) !/ surf(1)%ds(iteta, jzeta)                   
                    idof = idof +1
                 endif
                 if ( coil(icoil)%Lc /= 0 ) then !if geometry is free;
                    drho = coil(icoil)%I
                    ! dBn/dmt
                    dmx =  drho*cost*cosp
                    dmy =  drho*cost*sinp
                    dmz = -drho*sint
                    dBn(idof+1) = ( gx(iteta, jzeta, icoil) * dmx &
                         &        + gy(iteta, jzeta, icoil) * dmy &
                         &        + gz(iteta, jzeta, icoil) * dmz ) !/ surf(1)%ds(iteta, jzeta)
                    idof = idof +1
                    ! dBn/dmp
                    dmx = -drho*sint*sinp
                    dmy =  drho*sint*cosp
                    dmz =  0
                    dBn(idof+1) = ( gx(iteta, jzeta, icoil) * dmx &
                         &        + gy(iteta, jzeta, icoil) * dmy &
                         &        + gz(iteta, jzeta, icoil) * dmz ) !/ surf(1)%ds(iteta, jzeta)
                    idof = idof +1
                 endif           
              enddo  !end icoil;
           endif
           FATAL( bnormal , idof-dof_offset .ne. ldof, counting error in packing )             
           select case (case_bnormal)
           case (0)     ! no normalization over |B|;
              t1B(1:Ndof) = t1B(1:Ndof) + surf(1)%bn(iteta,jzeta) * dBn(1:Ndof) !* surf(1)%ds(iteta,jzeta)
           case default
              FATAL( bnorm, .true., case_bnormal can only be 0 )
           end select

           if (weight_resbn .gt. sqrtmachprec) then  ! resonant Bn error
               b1c = b1c + dBn(1:Ndof) * cosarg(iteta, jzeta) / surf(1)%ds(iteta, jzeta)
               b1s = b1s + dBn(1:Ndof) * sinarg(iteta, jzeta) / surf(1)%ds(iteta, jzeta)
           endif 

        enddo  !end iteta;
     enddo  !end jzeta;

     call MPI_ALLREDUCE( MPI_IN_PLACE, t1B, Ndof, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FAMUS, ierr )

     t1B = t1B * discretefactor

      if (weight_resbn .gt. sqrtmachprec) then 
         b1c = b1c * discretefactor
         b1s = b1s * discretefactor
         t1R = -sign(1.0_dp, sqrt(bnc*bnc+bns*bns)-target_resbn) * (bnc*bnc + bns*bns)**(-0.5)*(bnc * b1c + bns * b1s)
         call MPI_ALLREDUCE( MPI_IN_PLACE, t1R, Ndof, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FAMUS, ierr )
      endif 

  endif

  !--------------------------------------------------------------------------------------------

  call MPI_barrier( MPI_COMM_FAMUS, ierr )

  return
  
end subroutine bnormal

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine prepare_inductance()
  use globals, only: dp, ierr, iout, myid, ounit, zero, Ncpu, Ncoils_total, Nfp_raw, &
       coil, surf, bsconstant, cosnfp, sinnfp, Ncoils, Nzeta, Nteta, one, three
  use bnorm_mod
  use mpi
  implicit none

  ! local variables
  INTEGER :: is, ip, iteta, jzeta, icoil, astat, symmetry, npc
  REAL :: ox, oy, oz, rx, ry, rz, rr, rm3, rm5, rdotN

  ! return if allocated
  if (allocated(gx)) return

  SALLOCATE( gx, (0:Nteta-1,0:Nzeta-1,1:Ncoils), zero ) ! inductance matrix for calculating B.n
  SALLOCATE( gy, (0:Nteta-1,0:Nzeta-1,1:Ncoils), zero ) ! inductance matrix for calculating B.n
  SALLOCATE( gz, (0:Nteta-1,0:Nzeta-1,1:Ncoils), zero ) ! inductance matrix for calculating B.n
  
  if (myid==0) then   ! calculate the field from coil on the master cpu
     surf(1)%Bx = zero; surf(1)%By = zero; surf(1)%Bz = zero
     do jzeta = 0, Nzeta-1
        do iteta = 0, Nteta-1
           do icoil = 1, Ncoils         
              call bfield0(icoil, surf(1)%xx(iteta, jzeta), surf(1)%yy(iteta, jzeta), &
                   & surf(1)%zz(iteta, jzeta), dBx(0,0), dBy(0,0), dBz(0,0))
              surf(1)%Bx(iteta, jzeta) = surf(1)%Bx(iteta, jzeta) + dBx( 0, 0) * coil(icoil)%I * bsconstant
              surf(1)%By(iteta, jzeta) = surf(1)%By(iteta, jzeta) + dBy( 0, 0) * coil(icoil)%I * bsconstant 
              surf(1)%Bz(iteta, jzeta) = surf(1)%Bz(iteta, jzeta) + dBz( 0, 0) * coil(icoil)%I * bsconstant 
           enddo ! end do icoil
        enddo
     enddo
  else ! each slave cpu calculates their contribution (0:Ntheta-1, 0:Nzeta-1, 1:Ncoils*Npc*Symmetry).     
     do icoil = 1, Ncoils
        ! check if stellarator symmetric
        npc = Nfp_raw
        symmetry = 0 
        if (coil(icoil)%symmetry == 0) Npc = 1
        if (coil(icoil)%symmetry == 2) symmetry = 1
        do jzeta = 0, Nzeta-1
           do iteta = 0, Nteta-1
              do ip = 1, Npc
                 do is = 0, symmetry
                    ! position vector
                    ox = coil(icoil)%ox*cosnfp(ip) - coil(icoil)%oy*(-1)**is*sinnfp(ip)
                    oy = coil(icoil)%ox*sinnfp(ip) + coil(icoil)%oy*(-1)**is*cosnfp(ip)
                    oz = coil(icoil)%oz*(-1)**is
                    rx = surf(1)%xx(iteta, jzeta) - ox
                    ry = surf(1)%yy(iteta, jzeta) - oy
                    rz = surf(1)%zz(iteta, jzeta) - oz
                    rr = sqrt(rx*rx+ry*ry+rz*rz)
                    rm3 = one/(rr**3)
                    rm5 = one/(rr**5)
                    rdotN = three*(rx*surf(1)%nx(iteta, jzeta)+ry*surf(1)%ny(iteta, jzeta)+rz*surf(1)%nz(iteta, jzeta))*surf(1)%ds(iteta, jzeta)
                    ! compute gx, gy, gz
                    gx(iteta, jzeta, icoil) = gx(iteta, jzeta, icoil) &
                         & + (rdotN*rm5*rx - rm3*surf(1)%nx(iteta, jzeta)*surf(1)%ds(iteta, jzeta))*(-1)**is*cosnfp(ip) &
                         & + (rdotN*rm5*ry - rm3*surf(1)%ny(iteta, jzeta)*surf(1)%ds(iteta, jzeta))*(-1)**is*sinnfp(ip)
                    gy(iteta, jzeta, icoil) = gy(iteta, jzeta, icoil) &
                         & + (rdotN*rm5*ry - rm3*surf(1)%ny(iteta, jzeta)*surf(1)%ds(iteta, jzeta))*cosnfp(ip) &
                         & - (rdotN*rm5*rx - rm3*surf(1)%nx(iteta, jzeta)*surf(1)%ds(iteta, jzeta))*sinnfp(ip)                    
                    gz(iteta, jzeta, icoil) = gz(iteta, jzeta, icoil) + rdotN*rm5*rz - rm3*surf(1)%nz(iteta, jzeta)*surf(1)%ds(iteta, jzeta)
                 enddo
              enddo
           enddo ! end iteta
        enddo ! end jzeta
     enddo ! end icoil
     
     ! add contant
     gx = gx*bsconstant
     gy = gy*bsconstant
     gz = gz*bsconstant
  endif

  return
end subroutine prepare_inductance
