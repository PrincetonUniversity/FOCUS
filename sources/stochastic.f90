!title (tlength) ! Calculate coil variations and field error from variations.. (tkruger)

!latex \briefly{Write something...}

!latex \calledby{\link{solvers}}

!latex  \section{General}
!latex  Write something...

subroutine stochastic(ideriv,bnormmax,bnormavg)
  use globals, only: dp, zero, half, pi2, machprec, ncpu, myid, ounit, &
       coil, DoF, Ncoils, Nfixgeo, Ndof, &
       ittlen, mttlen, LM_fvec, LM_fjac, bnorm, Npert, sdelta, MPI_COMM_FOCUS

  use mpi
  implicit none
  INTEGER, INTENT(in) :: ideriv
  REAL, INTENT(out)   :: bnormmax, bnormavg

  INTEGER             :: astat, ierr, icoil, j, idof, ND, ivec, nxx(1:Ncoils,1:Npert), nxxhold
  REAL                :: d1L(1:Ndof), bnormhold, pertxx(0:coil(1)%NS-1), pertxt(0:coil(1)%NS-1), psx(1:Ncoils,1:Npert), psxhold, torsRet
  REAL, allocatable   :: xxhold(:,:), xthold(:,:), g(:,:), gp(:,:)

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  ! Stochastic variables set here
  do icoil = 1, Ncoils
     do j = 1, Npert
        if (myid .eq. 0) then
           call random_number(torsRet)
           psx(icoil,j) = pi2*torsRet
           call random_number(torsRet)
           nxx(icoil,j) = FLOOR(3.0*torsRet) + 1 ! Should give int between 1 and 3
        endif
        call MPI_BARRIER( MPI_COMM_FOCUS, ierr )
        RlBCAST( psx(icoil,j), 1, 0 )
        IlBCAST( nxx(icoil,j), 1, 0 )
        call MPI_BARRIER( MPI_COMM_FOCUS, ierr )
     enddo
  enddo

  call bnormal( 0 )
  bnormhold = bnorm
  bnormmax = zero
  bnormavg = zero

  ! Assume all coils have same resolution, CHANGE LATER
  SALLOCATE(xxhold, ( 1:Ncoils , 0:coil(1)%NS-1) , zero)
  SALLOCATE(xthold, ( 1:Ncoils , 0:coil(1)%NS-1) , zero)
  do icoil = 1, Ncoils
     xxhold(icoil,0:coil(1)%NS-1) = coil(icoil)%xx(0:coil(1)%NS-1)
     xthold(icoil,0:coil(1)%NS-1) = coil(icoil)%xt(0:coil(1)%NS-1)
     ! Do the same for y and z
  enddo

  SALLOCATE(g, ( 1:Ncoils , 0:coil(1)%NS-1) , zero)
  SALLOCATE(gp, ( 1:Ncoils , 0:coil(1)%NS-1) , zero)
  do icoil = 1, Ncoils
     call lenDeriv0( icoil, coil(icoil)%L )
     do j = 1, coil(icoil)%NS-1
        g(icoil,j) = sqrt( coil(icoil)%xt(j-1)**2 + coil(icoil)%yt(j-1)**2 + coil(icoil)%zt(j-1)**2 ) + g(icoil,j-1)
     enddo
     g(icoil,0:coil(1)%NS-1) = pi2*g(icoil,0:coil(1)%NS-1)/(coil(icoil)%NS*coil(icoil)%L)
     gp(icoil,0:coil(1)%NS-1) = sqrt( coil(icoil)%xt(0:coil(1)%NS-1)**2 + coil(icoil)%yt(0:coil(1)%NS-1)**2 + coil(icoil)%zt(0:coil(1)%NS-1)**2 )/coil(icoil)%L
  enddo
  ! Current perturbations are not parameterization dependent
  ! which can influence a stochastic optimization, but probably won't

  ! Set sdelta to sdelta/sqrt(3.0)
  do j = 1, Npert
     do icoil = 1, Ncoils
        if( coil(icoil)%type .ne. 1 ) cycle
        pertxx(0:coil(1)%NS-1) = sdelta*cos(pi2*nxx(icoil,j)*g(icoil,0:coil(1)%NS-1) + psx(icoil,j))
        pertxt(0:coil(1)%NS-1) = -pi2*nxx(icoil,j)*sdelta*sin(pi2*nxx(icoil,j)*g(icoil,0:coil(1)%NS-1) + psx(icoil,j))*gp(icoil,0:coil(1)%NS-1)
        coil(icoil)%xx(0:coil(1)%NS-1) = xxhold(icoil,0:coil(1)%NS-1)
        coil(icoil)%xt(0:coil(1)%NS-1) = xthold(icoil,0:coil(1)%NS-1)
        coil(icoil)%xx(0:coil(1)%NS-1) = coil(icoil)%xx(0:coil(1)%NS-1) + pertxx(0:coil(1)%NS-1)
        coil(icoil)%xt(0:coil(1)%NS-1) = coil(icoil)%xt(0:coil(1)%NS-1) + pertxt(0:coil(1)%NS-1)
     enddo
     call MPI_BARRIER( MPI_COMM_FOCUS, ierr ) ! Might not need this
     call bnormal( 0 )
     if ( bnormmax .le. bnorm ) bnormmax = bnorm
     bnormavg = bnormavg + bnorm
  enddo
  bnormavg = bnormavg / real(Npert)

  bnorm = bnormhold
  do icoil = 1, Ncoils
     coil(icoil)%xx(0:coil(1)%NS-1) = xxhold(icoil,0:coil(1)%NS-1)
     coil(icoil)%xt(0:coil(1)%NS-1) = xthold(icoil,0:coil(1)%NS-1)
  enddo

  DALLOCATE(xxhold)
  DALLOCATE(xthold)
  DALLOCATE(g)
  DALLOCATE(gp)

  return
end subroutine stochastic
