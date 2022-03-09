!title (tlength) ! Calculate coil variations and field error from variations.. (tkruger)

!latex \briefly{Write something...}

!latex \calledby{\link{solvers}}

!latex  \section{General}
!latex  Write something...

! not parallelized; communications may take more time;
subroutine stochastic(ideriv,psx,nxx,bnormmax,bnormavg)
  use globals, only: dp, zero, half, pi2, machprec, ncpu, myid, ounit, &
       coil, DoF, Ncoils, Nfixgeo, Ndof, &
       ittlen, mttlen, LM_fvec, LM_fjac, bnorm, Npert, sdelta, MPI_COMM_FOCUS

  use mpi
  implicit none
  INTEGER, INTENT(in) :: ideriv, nxx(Ncoils,Npert)
  REAL, INTENT(in)    :: psx(Ncoils,Npert)
  REAL, INTENT(out)   :: bnormmax, bnormavg

  INTEGER             :: astat, ierr, icoil, j, idof, ND, ivec
  REAL                :: d1L(1:Ndof), bnormhold, pertxx(0:coil(1)%NS-1), pertxt(0:coil(1)%NS-1)
  REAL, allocatable   :: xxhold(:,:), xthold(:,:), g(:,:), gp(:,:)

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  ! Call this routine from and set random variables from diagnos
  ! Any parallelization will mess up random numbers

  call bnormal( 0 )
  bnormhold = bnorm
  bnormmax = zero
  bnormavg = zero

  ! Assume all coils have same resolution, CHANGE LATER
  SALLOCATE(xxhold, ( 1:Ncoils , 0:coil(1)%NS-1) , zero)
  SALLOCATE(xthold, ( 1:Ncoils , 0:coil(1)%NS-1) , zero)
  do icoil = 1, Ncoils
     xxhold(icoil,:) = coil(icoil)%xx(:)
     xthold(icoil,:) = coil(icoil)%xt(:)
     ! Do the same for y and z
  enddo

  SALLOCATE(g, ( 1:Ncoils , 0:coil(1)%NS-1) , zero)
  SALLOCATE(gp, ( 1:Ncoils , 0:coil(1)%NS-1) , zero)
  do icoil = 1, Ncoils
     call lenDeriv0( icoil, coil(icoil)%L )
     do j = 1, coil(icoil)%NS-1
        g(icoil,j) = sqrt( coil(icoil)%xt(j-1)**2 + coil(icoil)%yt(j-1)**2 + coil(icoil)%zt(j-1)**2 ) + g(icoil,j-1)
     enddo
     g(icoil,:) = pi2*g(icoil,:)/(coil(icoil)%NS*coil(icoil)%L)
     gp(icoil,:) = sqrt( coil(icoil)%xt(:)**2 + coil(icoil)%yt(:)**2 + coil(icoil)%zt(:)**2 )/coil(icoil)%L
  enddo
  ! Current perturbations are not neccessarily parameterization dependent
  ! which can influence a stochastic optimization 

  ! This should be parallelized
  ! Set sdelta to sdelta/sqrt(3.0)
  do j = 1, Npert
     do icoil = 1, Ncoils
        if( coil(icoil)%type .ne. 1 ) cycle
        pertxx(:) = sdelta*cos(pi2*nxx(icoil,j)*g(icoil,:) + psx(icoil,j))
        pertxt(:) = -pi2*nxx(icoil,j)*sdelta*sin(pi2*nxx(icoil,j)*g(icoil,:) + psx(icoil,j))*gp(icoil,:)
        coil(icoil)%xx(:) = xxhold(icoil,:)
        coil(icoil)%xt(:) = xthold(icoil,:)
        coil(icoil)%xx(:) = coil(icoil)%xx(:) + pertxx(:)
        coil(icoil)%xt(:) = coil(icoil)%xt(:) + pertxt(:)
     enddo
     call bnormal( 0 )
     if ( bnormmax .le. bnorm ) bnormmax = bnorm
     bnormavg = bnormavg + bnorm
     ! Stop until perturbations are made and field is called
  enddo
  bnormavg = bnormavg / real(Npert)

  bnorm = bnormhold
  do icoil = 1, Ncoils
     coil(icoil)%xx(:) = xxhold(icoil,:)
     coil(icoil)%xt(:) = xthold(icoil,:)
  enddo

  DALLOCATE(xxhold)
  DALLOCATE(xthold)
  DALLOCATE(g)
  DALLOCATE(gp)

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!  ttlen = zero
  !ivec = 1

!  if( ideriv >= 0 ) then
!     do icoil = 1, Ncoils     !only care about unique coils;
!        if(coil(icoil)%type == 1) then  ! only for Fourier
           !if( myid.ne.modulo(icoil-1,ncpu) ) cycle ! parallelization loop;
!           call LenDeriv0(icoil, coil(icoil)%L)
!           if ( coil(icoil)%Lc /= 0 ) then
!              ttlen = ttlen +  half * (coil(icoil)%L - coil(icoil)%Lo)**2 / coil(icoil)%Lo**2
              !if (mttlen > 0) then ! L-M format of targets
                 !LM_fvec(ittlen+ivec) = weight_ttlen * (coil(icoil)%L - coil(icoil)%Lo)
                 !ivec = ivec + 1
              !endif
!           endif
!        endif
!     enddo
     !if (mttlen > 0) then ! L-M format of targets
     !   FATAL( length, ivec == mttlen, Errors in counting ivec for L-M )
     !endif
!     ttlen = ttlen / (Ncoils - Nfixgeo + machprec)
!  endif

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!  if ( ideriv >= 1 ) then
!     t1A = zero ; d1L = zero
!     idof = 0 !; ivec = 1
!     do icoil = 1, Ncoils
!        ND = DoF(icoil)%ND
!        if ( coil(icoil)%Ic /= 0 ) then !if current is free;
!           idof = idof +1
!        endif
!        if ( coil(icoil)%Lc /= 0 ) then !if geometry is free;
!           if(coil(icoil)%type .eq. 1) then ! only for Fourier
!              call lenDeriv1( icoil, d1L(idof+1:idof+ND), ND )
!              t1A(idof+1:idof+ND) = d1L(idof+1:idof+ND)
              !if (mttlen > 0) then ! L-M format of targets
                 !LM_fjac(ittlen+ivec, idof+1:idof+ND) = weight_ttlen * d1L(idof+1:idof+ND)
                 !ivec = ivec + 1
              !endif
!           endif 
!           idof = idof + ND
!        endif
!     enddo !end icoil;
!     FATAL( length , idof .ne. Ndof, counting error in packing )
     !if (mttlen > 0) then ! L-M format of targets
     !   FATAL( length, ivec == mttlen, Errors in counting ivec for L-M )
     !endif
!     t1A = t1A / (Ncoils - Nfixgeo + machprec)
!  endif

  return
end subroutine stochastic
