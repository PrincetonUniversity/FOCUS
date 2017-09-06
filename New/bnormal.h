!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

!title (bnormal) ! Calculate total bnormal and its derivatives.

!latex \briefly{Calculate the total bnormal of all coils on plasma surface and the derivatives}

!latex \calledby{\link{costfun}}
!latex \calls{\link{bfield}}

!latex \tableofcontents

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

subroutine bnormal( ideriv )
  
  use globals, only : zero, half, one, pi2, sqrtmachprec, ncpu, myid, ounit, &
                      bsconstant, &
                      surf, Nteta, Nzeta, Ntz, Ncoils, coil, DoF, discretefactor, &
                      Bdotnsquared, &
                      t1B, t2B, Ndof, dB, Cdof, weight_bharm
  
  implicit none

  include "mpif.h"

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  INTEGER, INTENT(in)                   :: ideriv

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

  INTEGER                              :: astat, ierr, icoil, iteta, jzeta! idof, ND, ip
  REAL, dimension(0:Nteta-1,0:Nzeta-1) :: lbx, lby, lbz        ! local Bx, By and Bz
  REAL                                 :: dBx, dBy, dBz        ! dB of each coil;
! REAL, dimension(0:Cdof, 0:Cdof)      :: dBx, dBy, dBz        ! dB of each coil;
! REAL, dimension(1:Ndof)              :: l1B
! REAL, allocatable                    :: ldB(:,:,:)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  Bdotnsquared = zero ! initialize intent out; 04 Sep 17;
  
  lbx(0:Nteta-1,0:Nzeta-1) = zero
  lby(0:Nteta-1,0:Nzeta-1) = zero
  lbz(0:Nteta-1,0:Nzeta-1) = zero
  
  do jzeta = 0, Nzeta-1
   do iteta = 0, Nteta-1
    
    if( myid.ne.modulo(jzeta*Nteta+iteta,ncpu) ) cycle ! parallelization loop;
    
    do icoil = 1, Ncoils
     
     call bfield( icoil, iteta, jzeta, dBx, dBy, dBz )

     lbx(iteta,jzeta) = lbx(iteta,jzeta) + dBx * coil(icoil)%I! * bsconstant
     lby(iteta,jzeta) = lby(iteta,jzeta) + dBy * coil(icoil)%I! * bsconstant
     lbz(iteta,jzeta) = lbz(iteta,jzeta) + dBz * coil(icoil)%I! * bsconstant
     
    enddo ! end do icoil
    
   enddo ! end do iteta
  enddo ! end do jzeta
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  call MPI_BARRIER( MPI_COMM_WORLD, ierr )
  
  call MPI_REDUCE( lbx(0:Nteta-1,0:Nzeta-1), surf(1)%Bx(0:Nteta-1,0:Nzeta-1), Ntz, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
  call MPI_REDUCE( lby(0:Nteta-1,0:Nzeta-1), surf(1)%By(0:Nteta-1,0:Nzeta-1), Ntz, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
  call MPI_REDUCE( lbz(0:Nteta-1,0:Nzeta-1), surf(1)%Bz(0:Nteta-1,0:Nzeta-1), Ntz, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
  
  RlBCAST( surf(1)%Bx(0:Nteta-1,0:Nzeta-1), Ntz, 0 )  ! total Bx from coils;
  RlBCAST( surf(1)%By(0:Nteta-1,0:Nzeta-1), Ntz, 0 )  ! total By from coils;
  RlBCAST( surf(1)%Bz(0:Nteta-1,0:Nzeta-1), Ntz, 0 )  ! total Bz from coils;
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
   
  surf(1)%bn(0:Nteta-1,0:Nzeta-1) = surf(1)%Bx*surf(1)%nx + surf(1)%By*surf(1)%ny + surf(1)%Bz*surf(1)%nz + surf(1)%pb ! total Bn, including plasma component;
  
  Bdotnsquared = sum( surf(1)%bn * surf(1)%bn * surf(1)%ds ) * half * discretefactor * bsconstant
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  if ( ideriv.ge.1 ) then
   
   FATAL( bnormal, .true., under reconstruction )

!     t1B = zero ; l1B = zero
!     if (weight_bharm > sqrtmachprec) then
!        SALLOCATE( ldB, (1:Ndof, 0:Nteta-1, 0:Nzeta-1), zero)
!        dB = zero
!     endif
 
      do jzeta = 0, Nzeta-1
         do iteta = 0, Nteta-1
 
            if( myid.ne.modulo(jzeta*Nteta+iteta,ncpu) ) cycle ! parallelization loop;
 
!           do ip = 1, Npc 
 
!              idof = 0
!              do icoil = 1, Ncoils
!                 ND = DoF(icoil)%ND
!                 if ( coil(icoil)%Ic /= 0 ) then !if current is free;
!                    call bfield( icoil, iteta, jzeta, dBx(0,0), dBy(0,0), dBz(0,0) )
!                    l1B(idof+1) = l1B(idof+1) + surf(1)%bn(iteta,jzeta) * surf(1)%ds(iteta,jzeta) &
!                         & * bsconstant * ( dBx(0,0)*surf(1)%nx(iteta,jzeta)   &
!                         &                + dBy(0,0)*surf(1)%ny(iteta,jzeta)   &
!                         &                + dBz(0,0)*surf(1)%nz(iteta,jzeta) )
!
!                    if (weight_bharm > sqrtmachprec) then
!                    ldB(idof+1, iteta, jzeta) = ldB(idof+1, iteta, jzeta)                      & 
!                         &                    + bsconstant * (dBx(0,0)*surf(1)%nx(iteta,jzeta) &
!                         &                                   +dBy(0,0)*surf(1)%ny(iteta,jzeta) &
!                         &                                   +dBz(0,0)*surf(1)%nz(iteta,jzeta) )
!                    endif
!
!                    idof = idof +1
!                 endif
!
!                 if ( coil(icoil)%Lc /= 0 ) then !if geometry is free;
!                    call bfield1(icoil+(ip-1)*Ncoils, iteta, jzeta, dBx(1:ND,0), dBy(1:ND,0), dBz(1:ND,0), ND)
!                    l1B(idof+1:idof+ND) = l1B(idof+1:idof+ND) + surf(1)%bn(iteta,jzeta)          &
!                         &              * surf(1)%ds(iteta,jzeta) * bsconstant * coil(icoil)%I   &
!                         &                            * ( dBx(1:ND,0)*surf(1)%nx(iteta,jzeta)    &
!                         &                              + dBy(1:ND,0)*surf(1)%ny(iteta,jzeta)    &
!                         &                              + dBz(1:ND,0)*surf(1)%nz(iteta,jzeta) )
!
!                    if (weight_bharm > sqrtmachprec) then
!                    ldB(idof+1:idof+ND, iteta, jzeta) = ldB(idof+1:idof+ND, iteta, jzeta)      &
!                         &                            + bsconstant * coil(icoil)%I             &
!                         &                            * (dBx(1:ND,0)*surf(1)%nx(iteta,jzeta)   &
!                         &                              +dBy(1:ND,0)*surf(1)%ny(iteta,jzeta)   &
!                         &                              +dBz(1:ND,0)*surf(1)%nz(iteta,jzeta) )
!                    endif
!
!                    idof = idof + ND
!
!                 endif
!
!              enddo  !end icoil;
!              FATAL( bnormal , idof .ne. Ndof, counting error in packing )
!
!           enddo  !end ip;
!
           enddo  !end iteta;
          enddo  !end jzeta;
!
!     call MPI_BARRIER( MPI_COMM_WORLD, ierr )
!     call MPI_REDUCE(l1B, t1B, Ndof, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
!
!     if (weight_bharm > sqrtmachprec) then
!        call MPI_REDUCE(ldB, dB, Ndof*Ntz, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
!        RlBCAST( dB, Ndof*Ntz, 0 )
!        DALLOCATE( ldB )
!     endif
!
!     RlBCAST( t1B, Ndof, 0 )
!     t1B = t1B * discretefactor

    endif ! end of if( ideriv.ge.1 ) ; 04 Sep 17;
    
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

  call MPI_barrier( MPI_COMM_WORLD, ierr )

  return

end subroutine bnormal

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
