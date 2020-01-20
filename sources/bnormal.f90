
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
  ! 1st-order
  REAL, allocatable :: dBn(:), dBm(:), d1B(:,:,:)

end module bnorm_mod
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine bnormal( ideriv )
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
       bharm, t1H, Bmnc, Bmns, wBmn, tBmnc, tBmns, Bmnim, Bmnin, NBmn
  use bnorm_mod
  use bharm_mod
  implicit none
  include "mpif.h"

  INTEGER, INTENT(in)                   :: ideriv
  !--------------------------------------------------------------------------------------------
  INTEGER                               :: astat, ierr
  INTEGER                               :: icoil, iteta, jzeta, idof, ND, NumGrid, ip, isurf

  !--------------------------initialize and allocate arrays-------------------------------------
  
  isurf = plasma 
  NumGrid = Nteta*Nzeta
  ! reset to zero;
  bnorm = zero 
  surf(isurf)%Bx = zero; surf(isurf)%By = zero; surf(isurf)%Bz = zero; surf(isurf)%Bn = zero     
  dBx = zero; dBy = zero; dBz = zero; Bm = zero

  bn = zero

  !-------------------------------calculate Bn-------------------------------------------------- 
  if( ideriv >= 0 ) then

     do jzeta = 0, Nzeta - 1
        do iteta = 0, Nteta - 1
           if( myid.ne.modulo(jzeta*Nteta+iteta,ncpu) ) cycle ! parallelization loop;

           do icoil = 1, Ncoils
              call bfield0(icoil, surf(isurf)%xx(iteta, jzeta), surf(isurf)%yy(iteta, jzeta), &
                   & surf(isurf)%zz(iteta, jzeta), dBx(0,0), dBy(0,0), dBz(0,0))
              surf(isurf)%Bx(iteta, jzeta) = surf(isurf)%Bx(iteta, jzeta) + dBx( 0, 0) 
              surf(isurf)%By(iteta, jzeta) = surf(isurf)%By(iteta, jzeta) + dBy( 0, 0) 
              surf(isurf)%Bz(iteta, jzeta) = surf(isurf)%Bz(iteta, jzeta) + dBz( 0, 0) 
           enddo ! end do icoil

           surf(isurf)%Bn(iteta, jzeta) = surf(isurf)%Bx(iteta, jzeta)*surf(isurf)%nx(iteta, jzeta) &
                &            + surf(isurf)%By(iteta, jzeta)*surf(isurf)%ny(iteta, jzeta) &
                &            + surf(isurf)%Bz(iteta, jzeta)*surf(isurf)%nz(iteta, jzeta) &
                &            - surf(isurf)%pb(iteta, jzeta)

           select case (case_bnormal)
           case (0)     ! no normalization over |B|;
              bnorm = bnorm + surf(isurf)%Bn(iteta, jzeta) * surf(isurf)%Bn(iteta, jzeta) * surf(isurf)%ds(iteta, jzeta)
           case (1)    ! normalized over |B|;
              Bm(iteta, jzeta) = surf(isurf)%Bx(iteta, jzeta)*surf(isurf)%Bx(iteta, jzeta) &
                &              + surf(isurf)%By(iteta, jzeta)*surf(isurf)%By(iteta, jzeta) &
                &              + surf(isurf)%Bz(iteta, jzeta)*surf(isurf)%Bz(iteta, jzeta)
              bnorm = bnorm + surf(isurf)%Bn(iteta, jzeta) * surf(isurf)%Bn(iteta, jzeta) &
                &             / Bm(iteta, jzeta) * surf(isurf)%ds(iteta, jzeta)
           case default
              FATAL( bnorm, .true., case_bnormal can only be 0/1 )
           end select

        enddo ! end do iteta
     enddo ! end do jzeta
     
     call MPI_BARRIER( MPI_COMM_WORLD, ierr )     
     call MPI_ALLREDUCE( MPI_IN_PLACE, surf(isurf)%Bx, NumGrid, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )
     call MPI_ALLREDUCE( MPI_IN_PLACE, surf(isurf)%By, NumGrid, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )
     call MPI_ALLREDUCE( MPI_IN_PLACE, surf(isurf)%Bz, NumGrid, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )
     call MPI_ALLREDUCE( MPI_IN_PLACE, surf(isurf)%Bn, NumGrid, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )
     call MPI_ALLREDUCE( MPI_IN_PLACE, bnorm, 1  , MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )

     bnorm = bnorm * half * discretefactor
     bn = surf(isurf)%Bn +  surf(isurf)%pb  ! bn is B.n from coils
     ! bn = surf(isurf)%Bx * surf(isurf)%nx + surf(isurf)%By * surf(isurf)%ny + surf(isurf)%Bz * surf(isurf)%nz
     !! if (case_bnormal == 0) bnorm = bnorm * bsconstant * bsconstant ! take bsconst back

     if (case_bnormal == 1) then    ! collect |B|
        call MPI_ALLREDUCE( MPI_IN_PLACE, Bm, NumGrid, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )
        !! bm = bm * bsconstant * bsconstant 
     endif
     
     ! Another type of target functions
     if (mbnorm > 0) then
        select case (case_bnormal)
        case (0)     ! no normalization over |B|;
           LM_fvec(ibnorm+1:ibnorm+mbnorm) =  weight_bnorm &
                &  * reshape(surf(isurf)%bn(0:Nteta-1, 0:Nzeta-1)                               , (/Nteta*Nzeta/))
        case (1)    ! normalized over |B|;
           LM_fvec(ibnorm+1:ibnorm+mbnorm) =  weight_bnorm &
                &  * reshape(surf(isurf)%bn(0:Nteta-1, 0:Nzeta-1)/sqrt(bm(0:Nteta-1, 0:Nzeta-1)), (/Nteta*Nzeta/))
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

     t1B = zero ; d1B = zero
     dBn = zero ; dBm = zero

     do jzeta = 0, Nzeta - 1
        do iteta = 0, Nteta - 1

           if( myid.ne.modulo(jzeta*Nteta+iteta,ncpu) ) cycle ! parallelization loop;

           do ip = 1, 1

              idof = 0
              do icoil = 1, Ncoils
                 ND = DoF(icoil)%ND
                 if ( coil(icoil)%Ic /= 0 ) then !if current is free;

                    call bfield0(icoil+(ip-1)*Ncoils, surf(isurf)%xx(iteta, jzeta), surf(isurf)%yy(iteta, jzeta), &
                   & surf(isurf)%zz(iteta, jzeta), dBx(0,0), dBy(0,0), dBz(0,0))
                    if (coil(icoil+(ip-1)*Ncoils)%type == 3) dBz(0,0) = zero  ! Bz doesn't change in type=3
                    dBn(idof+1) = ( dBx(0,0)*surf(isurf)%nx(iteta,jzeta) &
                         &        + dBy(0,0)*surf(isurf)%ny(iteta,jzeta) &
                         &        + dBz(0,0)*surf(isurf)%nz(iteta,jzeta) ) / coil(icoil+(ip-1)*Ncoils)%I

                    if (case_bnormal == 1) then  ! normalized over |B|;
                    dBm(idof+1) = ( dBx(0,0)*surf(isurf)%Bx(iteta,jzeta) &
                         &        + dBy(0,0)*surf(isurf)%By(iteta,jzeta) &
                         &        + dBz(0,0)*surf(isurf)%Bz(iteta,jzeta) ) / coil(icoil+(ip-1)*Ncoils)%I 
                    endif

                    idof = idof +1
                 endif

                 if ( coil(icoil)%Lc /= 0 ) then !if geometry is free;
                    call bfield1(icoil+(ip-1)*Ncoils, surf(isurf)%xx(iteta, jzeta), surf(isurf)%yy(iteta, jzeta), &
                   & surf(isurf)%zz(iteta, jzeta), dBx(1:ND,0), dBy(1:ND,0), dBz(1:ND,0), ND)
                    dBn(idof+1:idof+ND) = ( dBx(1:ND,0)*surf(isurf)%nx(iteta,jzeta) &
                         &                + dBy(1:ND,0)*surf(isurf)%ny(iteta,jzeta) &
                         &                + dBz(1:ND,0)*surf(isurf)%nz(iteta,jzeta) )
                    if (case_bnormal == 1) then  ! normalized over |B|;
                    dBm(idof+1:idof+ND) = ( dBx(1:ND,0)*surf(isurf)%Bx(iteta,jzeta) &
                         &                + dBy(1:ND,0)*surf(isurf)%By(iteta,jzeta) &
                         &                + dBz(1:ND,0)*surf(isurf)%Bz(iteta,jzeta) )
                    endif

                    idof = idof + ND

                 endif

              enddo  !end icoil;
              FATAL( bnormal , idof .ne. Ndof, counting error in packing )

              select case (case_bnormal)
              case (0)     ! no normalization over |B|;
                 t1B(1:Ndof) = t1B(1:Ndof) + surf(isurf)%bn(iteta,jzeta) * surf(isurf)%ds(iteta,jzeta) * dBn(1:Ndof)
                 d1B(1:Ndof, iteta, jzeta) =  d1B(1:Ndof, iteta, jzeta) + dBn(1:Ndof)
              case (1)     ! normalized over |B|;
                 t1B(1:Ndof) = t1B(1:Ndof) + ( surf(isurf)%Bn(iteta,jzeta) * dBn(1:Ndof) &
                      &                       / bm(iteta, jzeta)             &
                      &                       - surf(isurf)%Bn(iteta,jzeta) * surf(isurf)%Bn(iteta,jzeta) &
                      &                       / (bm(iteta, jzeta)*bm(iteta, jzeta)) &  
                      &                       * dBm(1:Ndof) ) * surf(isurf)%ds(iteta,jzeta)
                 d1B(1:Ndof, iteta, jzeta) = d1B(1:Ndof, iteta, jzeta) + dBn(1:Ndof) & 
                      &                    / sqrt(bm(iteta, jzeta)) &
                      &                    - surf(isurf)%Bn(iteta,jzeta) * dBm(1:Ndof) &
                      &                    / (bm(iteta, jzeta) * sqrt(bm(iteta, jzeta)))
              case default
                 FATAL( bnorm, .true., case_bnormal can only be 0/1 )
              end select
                 
           enddo  !end ip;

        enddo  !end iteta;
     enddo  !end jzeta;

     call MPI_BARRIER( MPI_COMM_WORLD, ierr )
     call MPI_ALLREDUCE( MPI_IN_PLACE, t1B, Ndof        , MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )
     call MPI_ALLREDUCE( MPI_IN_PLACE, d1B, Ndof*NumGrid, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )

     t1B = t1B * discretefactor

     ! Another type of target functions
     if (mbnorm > 0) then
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
end subroutine bnormal

