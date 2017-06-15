
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
       bnorm, t1B, t2B, bn, Ndof

  implicit none
  include "mpif.h"

  INTEGER, INTENT(in)                   :: ideriv
  !--------------------------------------------------------------------------------------------
  INTEGER                               :: astat, ierr
  INTEGER                               :: icoil, iteta, jzeta, idof, ND, NumGrid
  REAL, dimension(0:Nteta-1, 0:Nzeta-1) :: lbx, lby, lbz, lbn         ! local Bx, By and Bz
  REAL, dimension(1:Ndof, 0:Nteta-1, 0:Nzeta-1) :: ldB, dB

  !--------------------------initialize and allocate arrays------------------------------------- 

  NumGrid = Nteta*Nzeta
  bnorm = zero
  lbx = zero; lby = zero; lbz = zero; lbn =  zero        !already allocted; reset to zero;

  bn = zero
  surf(1)%bn = zero; surf(1)%Bx = zero; surf(1)%By = zero; surf(1)%Bz = zero
  do icoil = 1, Ncoils
     coil(icoil)%Bx = zero; coil(icoil)%By = zero; coil(icoil)%Bz = zero
  enddo

  !-------------------------------calculate Bn-------------------------------------------------- 
  if( ideriv >= 0 ) then

     do jzeta = 0, Nzeta - 1
        do iteta = 0, Nteta - 1
           if( myid.ne.modulo(jzeta*Nteta+iteta,ncpu) ) cycle ! parallelization loop;

           do icoil = 1, Ncoils
              call bfield0(icoil, iteta, jzeta, coil(icoil)%Bx(0,0), coil(icoil)%By(0,0), coil(icoil)%Bz(0,0))
              lbx(iteta, jzeta) = lbx(iteta, jzeta) + coil(icoil)%Bx( 0, 0) * coil(icoil)%I * bsconstant
              lby(iteta, jzeta) = lby(iteta, jzeta) + coil(icoil)%By( 0, 0) * coil(icoil)%I * bsconstant
              lbz(iteta, jzeta) = lbz(iteta, jzeta) + coil(icoil)%Bz( 0, 0) * coil(icoil)%I * bsconstant
           enddo ! end do icoil
           lbn(iteta, jzeta) = lbx(iteta, jzeta) * surf(1)%nx(iteta,jzeta)  &
                &            + lby(iteta, jzeta) * surf(1)%ny(iteta,jzeta)  &
                &            + lbz(iteta, jzeta) * surf(1)%nz(iteta,jzeta)
           !surf(1)%bn(iteta, jzeta) = lbn(iteta, jzeta) + surf(1)%pb(iteta, jzeta) !coilBn - targetBn;
        enddo ! end do iteta
     enddo ! end do jzeta
  
     call MPI_REDUCE( lbn, bn        , NumGrid, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
     call MPI_REDUCE( lbx, surf(1)%Bx, NumGrid, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
     call MPI_REDUCE( lby, surf(1)%By, NumGrid, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
     call MPI_REDUCE( lbz, surf(1)%Bz, NumGrid, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )

     RlBCAST( bn        , NumGrid, 0 )  ! coil Bn distribution;
     RlBCAST( surf(1)%Bx, NumGrid, 0 )  ! total Bx from coils;
     RlBCAST( surf(1)%By, NumGrid, 0 )  ! total By from coils;
     RlBCAST( surf(1)%Bz, NumGrid, 0 )  ! total Bz from coils;

     surf(1)%bn = bn + surf(1)%pb       ! total Bn;

     bnorm = sum( surf(1)%bn * surf(1)%bn * surf(1)%ds ) * half * discretefactor

     ! mapping the last  row & column; only used for plotting;
!!$             bn(Nteta, 0:Nzeta-1) =         bn(0, 0:Nzeta-1)
!!$     surf(1)%bn(Nteta, 0:Nzeta-1) = surf(1)%bn(0, 0:Nzeta-1)
!!$     surf(1)%Bx(Nteta, 0:Nzeta-1) = surf(1)%Bx(0, 0:Nzeta-1)
!!$     surf(1)%By(Nteta, 0:Nzeta-1) = surf(1)%By(0, 0:Nzeta-1)
!!$     surf(1)%Bz(Nteta, 0:Nzeta-1) = surf(1)%Bz(0, 0:Nzeta-1)
!!$
!!$             bn(0:Nteta-1, Nzeta) =         bn(0:Nteta-1, 0)
!!$     surf(1)%bn(0:Nteta-1, Nzeta) = surf(1)%bn(0:Nteta-1, 0)
!!$     surf(1)%Bx(0:Nteta-1, Nzeta) = surf(1)%Bx(0:Nteta-1, 0)
!!$     surf(1)%By(0:Nteta-1, Nzeta) = surf(1)%By(0:Nteta-1, 0)
!!$     surf(1)%Bz(0:Nteta-1, Nzeta) = surf(1)%Bz(0:Nteta-1, 0)
!!$
!!$             bn(Nteta   ,  Nzeta) =         bn(0    ,     0)
!!$     surf(1)%bn(Nteta   ,  Nzeta) = surf(1)%bn(0    ,     0)
!!$     surf(1)%Bx(Nteta   ,  Nzeta) = surf(1)%Bx(0    ,     0)
!!$     surf(1)%By(Nteta   ,  Nzeta) = surf(1)%By(0    ,     0)
!!$     surf(1)%Bz(Nteta   ,  Nzeta) = surf(1)%Bz(0    ,     0) 

  endif

  !-------------------------------calculate Bn/x------------------------------------------------
  if ( ideriv >= 1 ) then

     t1B = zero ; ldB = zero ; dB = zero

     do jzeta = 0, Nzeta - 1
        do iteta = 0, Nteta - 1

           if( myid.ne.modulo(jzeta*Nteta+iteta,ncpu) ) cycle ! parallelization loop;

           idof = 0
           
           do icoil = 1, Ncoils
              ND = DoF(icoil)%ND
              if ( coil(icoil)%Ic /= 0 ) then !if current is free;
                 call bfield0(icoil, iteta, jzeta, &
                      & coil(icoil)%Bx(0,0), coil(icoil)%By(0,0), coil(icoil)%Bz(0,0))
                 ldB(idof+1, iteta, jzeta) = bsconstant * ( coil(icoil)%Bx(0,0)*surf(1)%nx(iteta,jzeta)   &
                      &                                   + coil(icoil)%By(0,0)*surf(1)%ny(iteta,jzeta)   &
                      &                                   + coil(icoil)%Bz(0,0)*surf(1)%nz(iteta,jzeta) )
                 idof = idof +1
              endif

              if ( coil(icoil)%Lc /= 0 ) then !if geometry is free;
                 call bfield1(icoil, iteta, jzeta, &
                      &       coil(icoil)%Bx(1:ND,0), coil(icoil)%By(1:ND,0), coil(icoil)%Bz(1:ND,0), ND)

                 ldB(idof+1:idof+ND, iteta, jzeta) = bsconstant * coil(icoil)%I          &
                      &                            * ( coil(icoil)%Bx(1:ND,0)*surf(1)%nx(iteta,jzeta)   &
                      &                              + coil(icoil)%By(1:ND,0)*surf(1)%ny(iteta,jzeta)   &
                      &                              + coil(icoil)%Bz(1:ND,0)*surf(1)%nz(iteta,jzeta) )

                 idof = idof + ND
              endif
              
           enddo
           FATAL( bnormal , idof .ne. Ndof, counting error in packing )

        enddo
     enddo

     call MPI_REDUCE(ldB, dB, Ndof*NumGrid, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
     RlBCAST( dB, Ndof*NumGrid, 0 )

     do idof = 1, Ndof
        t1B(idof) = sum( surf(1)%bn(0:Nteta-1, 0:Nzeta-1) * dB(idof, 0:Nteta-1, 0:Nzeta-1) &
                        * surf(1)%ds(0:Nteta-1, 0:Nzeta-1) ) * discretefactor
     enddo


  endif

  !--------------------------------------------------------------------------------------------

  call MPI_barrier( MPI_COMM_WORLD, ierr )
  
  return
end subroutine bnormal
