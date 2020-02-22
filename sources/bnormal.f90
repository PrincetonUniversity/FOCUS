
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

end module bnorm_mod
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine bnormal_old( ideriv )
!------------------------------------------------------------------------------------------------------ 
! DATE:  04/02/2017;
! Calculate the Bn surface integral and its derivatives;
! ideriv = 0 -> only calculate the Bn surface integral;
! ideriv = 1 -> calculate the Bn surface integral and its first derivatives;
! ideriv = 2 -> calculate the Bn surface integral and its first & second derivatives;
!------------------------------------------------------------------------------------------------------   
  use globals, only: dp, zero, half, one, pi2, sqrtmachprec, bsconstant, ncpu, myid, ounit, &
       coil, DoF, surf, Ncoils, Nteta, Nzeta, discretefactor, &
       bnorm, t1B, t2B, bn, Ndof, Npc, Cdof, weight_bharm, case_bnormal, &
       weight_bnorm, ibnorm, mbnorm, ibharm, mbharm, LM_fvec, LM_fjac, &
       bharm, t1H, Bmnc, Bmns, wBmn, tBmnc, tBmns, Bmnim, Bmnin, NBmn, dof_offset, ldof, momentq
  use bnorm_mod
  use bharm_mod
  implicit none
  include "mpif.h"

  INTEGER, INTENT(in)                   :: ideriv
  !--------------------------------------------------------------------------------------------
  INTEGER                               :: astat, ierr
  INTEGER                               :: icoil, iteta, jzeta, idof, ND, NumGrid, ip

  !--------------------------initialize and allocate arrays------------------------------------- 

  NumGrid = Nteta*Nzeta
  ! reset to zero;
  bnorm = zero 
  surf(1)%Bx = zero; surf(1)%By = zero; surf(1)%Bz = zero; surf(1)%Bn = zero     
  dBx = zero; dBy = zero; dBz = zero; Bm = zero

  bn = zero

  !-------------------------------calculate Bn-------------------------------------------------- 
  if( ideriv >= 0 ) then

     do jzeta = 0, Nzeta - 1
        do iteta = 0, Nteta - 1
           !if( myid.ne.modulo(jzeta*Nteta+iteta,ncpu) ) cycle ! parallelization loop;

           do icoil = 1, Ncoils
              call bfield0(icoil, surf(1)%xx(iteta, jzeta), surf(1)%yy(iteta, jzeta), &
                   & surf(1)%zz(iteta, jzeta), dBx(0,0), dBy(0,0), dBz(0,0))
              surf(1)%Bx(iteta, jzeta) = surf(1)%Bx(iteta, jzeta) + dBx( 0, 0) * coil(icoil)%I * bsconstant
              surf(1)%By(iteta, jzeta) = surf(1)%By(iteta, jzeta) + dBy( 0, 0) * coil(icoil)%I * bsconstant 
              surf(1)%Bz(iteta, jzeta) = surf(1)%Bz(iteta, jzeta) + dBz( 0, 0) * coil(icoil)%I * bsconstant 
           enddo ! end do icoil

           call MPI_ALLREDUCE( MPI_IN_PLACE, surf(1)%Bx(iteta, jzeta), 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )
           call MPI_ALLREDUCE( MPI_IN_PLACE, surf(1)%By(iteta, jzeta), 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )
           call MPI_ALLREDUCE( MPI_IN_PLACE, surf(1)%Bz(iteta, jzeta), 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )

           surf(1)%Bn(iteta, jzeta) = surf(1)%Bx(iteta, jzeta)*surf(1)%nx(iteta, jzeta) &
                &            + surf(1)%By(iteta, jzeta)*surf(1)%ny(iteta, jzeta) &
                &            + surf(1)%Bz(iteta, jzeta)*surf(1)%nz(iteta, jzeta) &
                &            + surf(1)%pb(iteta, jzeta)

           select case (case_bnormal)
           case (0)     ! no normalization over |B|;
              bnorm = bnorm + surf(1)%Bn(iteta, jzeta) * surf(1)%Bn(iteta, jzeta) * surf(1)%ds(iteta, jzeta)
           case (1)    ! normalized over |B|;
              Bm(iteta, jzeta) = surf(1)%Bx(iteta, jzeta)*surf(1)%Bx(iteta, jzeta) &
                &              + surf(1)%By(iteta, jzeta)*surf(1)%By(iteta, jzeta) &
                &              + surf(1)%Bz(iteta, jzeta)*surf(1)%Bz(iteta, jzeta)
              bnorm = bnorm + surf(1)%Bn(iteta, jzeta) * surf(1)%Bn(iteta, jzeta) &
                &             / Bm(iteta, jzeta) * surf(1)%ds(iteta, jzeta)
           case default
              FATAL( bnorm, .true., case_bnormal can only be 0/1 )
           end select

        enddo ! end do iteta
     enddo ! end do jzeta
     

     bnorm = bnorm * half * discretefactor
     bn = surf(1)%Bn -  surf(1)%pb  ! bn is B.n from coils
     ! bn = surf(1)%Bx * surf(1)%nx + surf(1)%By * surf(1)%ny + surf(1)%Bz * surf(1)%nz
     !! if (case_bnormal == 0) bnorm = bnorm * bsconstant * bsconstant ! take bsconst back
     
     ! Another type of target functions
     if (mbnorm > 0) then
        select case (case_bnormal)
        case (0)     ! no normalization over |B|;
           LM_fvec(ibnorm+1:ibnorm+mbnorm) =  weight_bnorm &
                &  * reshape(surf(1)%bn(0:Nteta-1, 0:Nzeta-1)                               , (/Nteta*Nzeta/))
        case (1)    ! normalized over |B|;
           LM_fvec(ibnorm+1:ibnorm+mbnorm) =  weight_bnorm &
                &  * reshape(surf(1)%bn(0:Nteta-1, 0:Nzeta-1)/sqrt(bm(0:Nteta-1, 0:Nzeta-1)), (/Nteta*Nzeta/))
        case default
           FATAL( bnorm, .true., case_bnormal can only be 0/1 )
        end select
           
     endif

     ! Bn harmonics related
     if (weight_bharm > sqrtmachprec) then
        call twodft( bn, Bmns, Bmnc, Bmnim, Bmnin, NBmn ) ! Bn from coils
        bharm = half * sum( wBmn * ((Bmnc - tBmnc)**2 + (Bmns - tBmns)**2) )

        if (mbharm > 0) then
           LM_fvec(ibharm+1:ibharm+mbharm/2) = weight_bharm * wBmn * (Bmnc - tBmnc)
           LM_fvec(ibharm+mbharm/2+1:ibharm+mbharm) = weight_bharm * wBmn * (Bmns - tBmns)
        endif

     endif
  endif

  !-------------------------------calculate Bn/x------------------------------------------------
  if ( ideriv >= 1 ) then

     t1B = zero 
     if (mbnorm > 0 .or. weight_bharm > sqrtmachprec)  d1B = zero
     dBn = zero ; dBm = zero

     do jzeta = 0, Nzeta - 1
        do iteta = 0, Nteta - 1
           
           idof = dof_offset
           do icoil = 1, Ncoils
              ND = DoF(icoil)%ND
              if ( coil(icoil)%Ic /= 0 ) then !if current is free;
                 call bfield0(icoil, surf(1)%xx(iteta, jzeta), surf(1)%yy(iteta, jzeta), &
                      & surf(1)%zz(iteta, jzeta), dBx(0,0), dBy(0,0), dBz(0,0))
                 if (coil(icoil)%itype == 3) dBz(0,0) = zero  ! Bz doesn't change in itype=3
                 dBn(idof+1) = ( dBx(0,0)*surf(1)%nx(iteta,jzeta) &
                      &        + dBy(0,0)*surf(1)%ny(iteta,jzeta) &
                      &        + dBz(0,0)*surf(1)%nz(iteta,jzeta) ) * bsconstant
                ! if (coil(icoil)%itype == 2) dBn(idof+1) = dBn(idof+1)*(coil(icoil)%moment*momentq*sin(coil(icoil)%pho)**(momentq-1)*cos(coil(icoil)%pho)) 
                 if (coil(icoil)%itype == 2) dBn(idof+1) = dBn(idof+1)*(coil(icoil)%moment*momentq*(coil(icoil)%pho)**(momentq-1))
                 if (case_bnormal == 1) then  ! normalized over |B|;
                    dBm(idof+1) = ( dBx(0,0)*surf(1)%Bx(iteta,jzeta) &
                         &        + dBy(0,0)*surf(1)%By(iteta,jzeta) &
                         &        + dBz(0,0)*surf(1)%Bz(iteta,jzeta) ) * bsconstant ! two is canceled below
                   ! if (coil(icoil)%itype == 2) dBm(idof+1) = dBm(idof+1)*(coil(icoil)%moment*momentq*sin(coil(icoil)%pho)**(momentq-1)*cos(coil(icoil)%pho))
                    if (coil(icoil)%itype == 2) dBm(idof+1) = dBm(idof+1)*(coil(icoil)%moment*momentq*(coil(icoil)%pho)**(momentq-1))
                 endif

                 idof = idof +1
              endif

              if ( coil(icoil)%Lc /= 0 ) then !if geometry is free;
                 call bfield1(icoil, surf(1)%xx(iteta, jzeta), surf(1)%yy(iteta, jzeta), &
                      & surf(1)%zz(iteta, jzeta), dBx(1:ND,0), dBy(1:ND,0), dBz(1:ND,0), ND)
                 dBn(idof+1:idof+ND) = ( dBx(1:ND,0)*surf(1)%nx(iteta,jzeta) &
                      &                + dBy(1:ND,0)*surf(1)%ny(iteta,jzeta) &
                      &                + dBz(1:ND,0)*surf(1)%nz(iteta,jzeta) )
                 if (case_bnormal == 1) then  ! normalized over |B|;
                    dBm(idof+1:idof+ND) = ( dBx(1:ND,0)*surf(1)%Bx(iteta,jzeta) &
                         &                + dBy(1:ND,0)*surf(1)%By(iteta,jzeta) &
                         &                + dBz(1:ND,0)*surf(1)%Bz(iteta,jzeta) )
                 endif

                 idof = idof + ND

              endif
              
!!$              if (myid==0 .and. iteta==0 .and. jzeta==0) then
!!$                 if (icoil==1) then
!!$                    print *, dBx(1:3,0)
!!$                    print *, dBy(1:3,0)
!!$                    print *, dBz(1:3,0)
!!$                    print *, dBn(1:3)
!!$                 endif
!!$              endif
           enddo  !end icoil;
           FATAL( bnormal , idof-dof_offset .ne. ldof, counting error in packing )

!!$           call MPI_ALLREDUCE( MPI_IN_PLACE, dBn, Ndof, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )
!!$           if (case_bnormal == 1) then
!!$              call MPI_ALLREDUCE( MPI_IN_PLACE, dBm, Ndof, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )
!!$           endif
              
           select case (case_bnormal)
           case (0)     ! no normalization over |B|;
              t1B(1:Ndof) = t1B(1:Ndof) + surf(1)%bn(iteta,jzeta) * surf(1)%ds(iteta,jzeta) * dBn(1:Ndof)
!!$              if (myid==0 .and. iteta==0 .and. jzeta==0) then
!!$                 print *, surf(1)%bn(iteta,jzeta), surf(1)%ds(iteta,jzeta)
!!$                 print *, dBn(1:10)
!!$              endif
              if (mbnorm > 0) then ! L-M
                 d1B(1:Ndof, iteta, jzeta) =  d1B(1:Ndof, iteta, jzeta) + dBn(1:Ndof)
              endif 
           case (1)     ! normalized over |B|;
              t1B(1:Ndof) = t1B(1:Ndof) + ( surf(1)%Bn(iteta,jzeta) * dBn(1:Ndof) &
                   &                       / bm(iteta, jzeta)             &
                   &                       - surf(1)%Bn(iteta,jzeta) * surf(1)%Bn(iteta,jzeta) &
                   &                       / (bm(iteta, jzeta)*bm(iteta, jzeta)) &  
                   &                       * dBm(1:Ndof) ) * surf(1)%ds(iteta,jzeta)
              if (mbnorm > 0) then ! L-M
                 d1B(1:Ndof, iteta, jzeta) = d1B(1:Ndof, iteta, jzeta) + dBn(1:Ndof) & 
                      &                    / sqrt(bm(iteta, jzeta)) &
                      &                    - surf(1)%Bn(iteta,jzeta) * dBm(1:Ndof) &
                      &                    / (bm(iteta, jzeta) * sqrt(bm(iteta, jzeta)))
              endif
           case default
              FATAL( bnorm, .true., case_bnormal can only be 0/1 )
           end select

        enddo  !end iteta;
     enddo  !end jzeta;

     call MPI_ALLREDUCE( MPI_IN_PLACE, t1B, Ndof, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )

     t1B = t1B * discretefactor

     ! Another type of target functions
     if (mbnorm > 0) then
        call MPI_ALLREDUCE( MPI_IN_PLACE, d1B, Ndof*Nteta*Nzeta, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )
        do idof = 1, Ndof
           LM_fjac(ibnorm+1:ibnorm+mbnorm, idof) = weight_bnorm &
                &  * reshape(d1B(idof, 0:Nteta-1, 0:Nzeta-1), (/Nteta*Nzeta/))
        enddo
     endif

     if (weight_bharm > sqrtmachprec) then
        dBc = zero ; dBs = zero
        do idof = 1, Ndof
           call twodft( d1B(idof,  0:Nteta-1, 0:Nzeta-1), dBs, dBc, Bmnim, Bmnin, NBmn )
           t1H(idof) = sum( wBmn * ( (Bmnc - tBmnc)*dBc + (Bmns - tBmns)*dBs ) )
           if (mbharm > 0) then
              LM_fjac(ibharm+1         :ibharm+mbharm/2, idof) = weight_bharm * wBmn * dBc
              LM_fjac(ibharm+mbharm/2+1:ibharm+mbharm  , idof) = weight_bharm * wBmn * dBs
           endif
        enddo    
     endif

  endif

  !--------------------------------------------------------------------------------------------

  call MPI_barrier( MPI_COMM_WORLD, ierr )

  return
end subroutine bnormal_old

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine bnormal( ideriv )
!------------------------------------------------------------------------------------------------------ 
! DATE:  10/28/2019;
! Calculate the Bn surface integral and its derivatives;
! ideriv = 0 -> only calculate the Bn surface integral;
! ideriv = 1 -> calculate the Bn surface integral and its first derivatives;
!------------------------------------------------------------------------------------------------------   
  use globals, only: dp, zero, half, one, pi2, sqrtmachprec, bsconstant, ncpu, myid, ounit, &
       coil, DoF, surf, Ncoils, Nteta, Nzeta, discretefactor, npc, cosnfp, sinnfp, &
       bnorm, t1B, t2B, bn, Ndof, Npc, Cdof, weight_bharm, case_bnormal, &
       weight_bnorm, ibnorm, mbnorm, ibharm, mbharm, LM_fvec, LM_fjac, &
       bharm, t1H, Bmnc, Bmns, wBmn, tBmnc, tBmns, Bmnim, Bmnin, NBmn, dof_offset, ldof, momentq
  use bnorm_mod
  use bharm_mod
  use mpi
  implicit none

  INTEGER, INTENT(in)                   :: ideriv
  !--------------------------------------------------------------------------------------------
  INTEGER                               :: astat, ierr
  INTEGER                               :: icoil, iteta, jzeta, idof, ND, NumGrid, ip, is, index
  REAL                                  :: sinp, sint, cosp, cost, drho, dmx, dmy, dmz

  !--------------------------initialize and allocate arrays------------------------------------- 

  NumGrid = Nteta*Nzeta
  ! reset to zero;
  bnorm = zero 
  surf(1)%Bn = zero     
  dBx = zero; dBy = zero; dBz = zero; Bm = zero
  bn = zero

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
           call MPI_ALLREDUCE( MPI_IN_PLACE, surf(1)%Bn(iteta, jzeta), 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )
           
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

  endif
  
  !-------------------------------calculate Bn/x------------------------------------------------
  if ( ideriv >= 1 ) then

     t1B = zero 
     if (mbnorm > 0 .or. weight_bharm > sqrtmachprec)  d1B = zero
     dBn = zero ; dBm = zero

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

        enddo  !end iteta;
     enddo  !end jzeta;

     call MPI_ALLREDUCE( MPI_IN_PLACE, t1B, Ndof, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )

     t1B = t1B * discretefactor

  endif

  !--------------------------------------------------------------------------------------------

  call MPI_barrier( MPI_COMM_WORLD, ierr )

  return
  
end subroutine bnormal

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine prepare_inductance()
  use globals, only: dp, ierr, iout, myid, ounit, zero, Ncpu, Ncoils_total, Npc, &
       coil, surf, bsconstant, cosnfp, sinnfp, Ncoils, Nzeta, Nteta, one, three
  use bnorm_mod
  use mpi
  implicit none

  ! local variables
  INTEGER :: is, ip, iteta, jzeta, icoil, astat, symmetry
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
        if (coil(icoil)%symmetry == 2) then
           symmetry = 1
        else
           symmetry = 0
        endif
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
