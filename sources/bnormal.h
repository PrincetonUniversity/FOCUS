
!title (bnormal) ! Calculate total bnormal and its derivatives.

!latex \briefly{Calculate the total bnormal of all coils on plasma surface and the derivatives with respect to coil geometry and currents, including the first and second dirivatives. 
!latex          Calling \emph{bnormal(0), bnormal(1), bnormal(2)} calculates the $0-order$, $1^{st}-order$ and $2^{nd}-order$ derivatives respectively.}

!latex \calledby{\link{costfun}}
!latex \calls{\link{bfield}}

!latex \tableofcontents
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
subroutine bnormal( ideriv )
!------------------------------------------------------------------------------------------------------ 
! DATE:  04/02/2017;
! Calculate the Bn surface integral and its derivatives;
! ideriv = 0 -> only calculate the Bn surface integral;
! ideriv = 1 -> calculate the Bn surface integral and its first derivatives;
! ideriv = 2 -> calculate the Bn surface integral and its first & second derivatives;
!------------------------------------------------------------------------------------------------------   
  use globals, only: zero, half, one, pi2, sqrtmachprec, bsconstant, ncpu, myid, ounit, &
       coil, DoF, surf, Ncoils, Nteta, Nzeta, discretefactor, &
       bnorm, t1B, t2B, bn, Ndof, Npc, dB, Cdof, weight_bharm

  implicit none
  include "mpif.h"

  INTEGER, INTENT(in)                   :: ideriv
  !--------------------------------------------------------------------------------------------
  INTEGER                               :: astat, ierr
  INTEGER                               :: icoil, iteta, jzeta, idof, ND, NumGrid, ip
  REAL, dimension(0:Nteta-1, 0:Nzeta-1) :: lbx, lby, lbz        ! local Bx, By and Bz
  REAL, dimension(0:Cdof, 0:Cdof)       :: dBx, dBy, dBz        ! dB of each coil;
  REAL, dimension(1:Ndof)               :: l1B
  REAL, allocatable                     :: ldB(:,:,:)

  !--------------------------initialize and allocate arrays------------------------------------- 

  NumGrid = Nteta*Nzeta
  bnorm = zero
  lbx = zero; lby = zero; lbz = zero        !already allocted; reset to zero;
  dBx = zero; dBy = zero; dBz = zero

  bn = zero
  surf(1)%bn = zero; surf(1)%Bx = zero; surf(1)%By = zero; surf(1)%Bz = zero

  !-------------------------------calculate Bn-------------------------------------------------- 
  if( ideriv >= 0 ) then

 
     do jzeta = 0, Nzeta - 1
        do iteta = 0, Nteta - 1
           if( myid.ne.modulo(jzeta*Nteta+iteta,ncpu) ) cycle ! parallelization loop;

           do icoil = 1, Ncoils*Npc
              call bfield0(icoil, iteta, jzeta, dBx(0,0), dBy(0,0), dBz(0,0))
              lbx(iteta, jzeta) = lbx(iteta, jzeta) + dBx( 0, 0) * coil(icoil)%I * bsconstant
              lby(iteta, jzeta) = lby(iteta, jzeta) + dBy( 0, 0) * coil(icoil)%I * bsconstant
              lbz(iteta, jzeta) = lbz(iteta, jzeta) + dBz( 0, 0) * coil(icoil)%I * bsconstant
           enddo ! end do icoil

        enddo ! end do iteta
     enddo ! end do jzeta
     
     call MPI_BARRIER( MPI_COMM_WORLD, ierr )     
     call MPI_REDUCE( lbx, surf(1)%Bx, NumGrid, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
     call MPI_REDUCE( lby, surf(1)%By, NumGrid, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
     call MPI_REDUCE( lbz, surf(1)%Bz, NumGrid, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
     
     RlBCAST( surf(1)%Bx, NumGrid, 0 )  ! total Bx from coils;
     RlBCAST( surf(1)%By, NumGrid, 0 )  ! total By from coils;
     RlBCAST( surf(1)%Bz, NumGrid, 0 )  ! total Bz from coils;

     bn =  surf(1)%Bx * surf(1)%nx + surf(1)%By * surf(1)%ny + surf(1)%Bz * surf(1)%nz

     surf(1)%bn = bn - surf(1)%pb       ! total Bn; change to minus on 05/28/2018

     bnorm = sum( surf(1)%bn * surf(1)%bn * surf(1)%ds ) * half * discretefactor

  endif

  !-------------------------------calculate Bn/x------------------------------------------------
  if ( ideriv >= 1 ) then

     t1B = zero ; l1B = zero
     if (weight_bharm > sqrtmachprec) then
        SALLOCATE( ldB, (1:Ndof, 0:Nteta-1, 0:Nzeta-1), zero)
        dB = zero
     endif

     do jzeta = 0, Nzeta - 1
        do iteta = 0, Nteta - 1

           if( myid.ne.modulo(jzeta*Nteta+iteta,ncpu) ) cycle ! parallelization loop;

           do ip = 1, Npc

              idof = 0
              do icoil = 1, Ncoils
                 ND = DoF(icoil)%ND
                 if ( coil(icoil)%Ic /= 0 ) then !if current is free;
                    call bfield0(icoil+(ip-1)*Ncoils, iteta, jzeta, dBx(0,0), dBy(0,0), dBz(0,0))
                    l1B(idof+1) = l1B(idof+1) + surf(1)%bn(iteta,jzeta) * surf(1)%ds(iteta,jzeta) &
                         & * bsconstant * ( dBx(0,0)*surf(1)%nx(iteta,jzeta)   &
                         &                + dBy(0,0)*surf(1)%ny(iteta,jzeta)   &
                         &                + dBz(0,0)*surf(1)%nz(iteta,jzeta) )

                    if (weight_bharm > sqrtmachprec) then
                    ldB(idof+1, iteta, jzeta) = ldB(idof+1, iteta, jzeta)                      & 
                         &                    + bsconstant * (dBx(0,0)*surf(1)%nx(iteta,jzeta) &
                         &                                   +dBy(0,0)*surf(1)%ny(iteta,jzeta) &
                         &                                   +dBz(0,0)*surf(1)%nz(iteta,jzeta) )
                    endif

                    idof = idof +1
                 endif

                 if ( coil(icoil)%Lc /= 0 ) then !if geometry is free;
                    call bfield1(icoil+(ip-1)*Ncoils, iteta, jzeta, dBx(1:ND,0), dBy(1:ND,0), dBz(1:ND,0), ND)
                    l1B(idof+1:idof+ND) = l1B(idof+1:idof+ND) + surf(1)%bn(iteta,jzeta)          &
                         &              * surf(1)%ds(iteta,jzeta) * bsconstant * coil(icoil)%I   &
                         &                            * ( dBx(1:ND,0)*surf(1)%nx(iteta,jzeta)    &
                         &                              + dBy(1:ND,0)*surf(1)%ny(iteta,jzeta)    &
                         &                              + dBz(1:ND,0)*surf(1)%nz(iteta,jzeta) )

                    if (weight_bharm > sqrtmachprec) then
                    ldB(idof+1:idof+ND, iteta, jzeta) = ldB(idof+1:idof+ND, iteta, jzeta)      &
                         &                            + bsconstant * coil(icoil)%I             &
                         &                            * (dBx(1:ND,0)*surf(1)%nx(iteta,jzeta)   &
                         &                              +dBy(1:ND,0)*surf(1)%ny(iteta,jzeta)   &
                         &                              +dBz(1:ND,0)*surf(1)%nz(iteta,jzeta) )
                    endif

                    idof = idof + ND

                 endif

              enddo  !end icoil;
              FATAL( bnormal , idof .ne. Ndof, counting error in packing )

           enddo  !end ip;

        enddo  !end iteta;
     enddo  !end jzeta;

     call MPI_BARRIER( MPI_COMM_WORLD, ierr )
     call MPI_REDUCE(l1B, t1B, Ndof, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )

     if (weight_bharm > sqrtmachprec) then
        call MPI_REDUCE(ldB, dB, Ndof*NumGrid, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
        RlBCAST( dB, Ndof*NumGrid, 0 )
        DALLOCATE( ldB )
     endif

     RlBCAST( t1B, Ndof, 0 )
     t1B = t1B * discretefactor

  endif

  !--------------------------------------------------------------------------------------------

  call MPI_barrier( MPI_COMM_WORLD, ierr )

  return
end subroutine bnormal