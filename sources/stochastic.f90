!title (tlength) ! Calculate coil variations and field error from variations.. (tkruger)

!latex \briefly{Write something...}

!latex \calledby{\link{solvers}}

!latex  \section{General}
!latex  Write something...

subroutine sbnormal(ideriv)
  use globals, only: dp, zero, half, pi2, machprec, ncpu, myid, ounit, &
       coil, DoF, Ncoils, Nfixgeo, Ndof, bnorm, t1B, t1Bavg, Npert, sdelta, Nmax, &
       bnormavg, bnormmax, surf, plasma, DoF, MPI_COMM_FOCUS

  use mpi
  implicit none
  INTEGER, INTENT(in) :: ideriv
  INTEGER             :: astat, ierr, icoil, j, idof, ivec, NSmax, NDmax, jj, NS
  REAL, allocatable   :: xxhold(:,:), xthold(:,:), yyhold(:,:), ythold(:,:), &
                         zzhold(:,:), zthold(:,:), g(:,:), gp(:,:), dgddof(:,:,:), &
                         d1L(:,:), xofhold(:,:,:), yofhold(:,:,:), zofhold(:,:,:)

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
 
  NSmax = zero
  NDmax = zero 
  do icoil = 1, Ncoils
     if ( coil(icoil)%NS .gt. NSmax ) NSmax = coil(icoil)%NS
     if ( DoF(icoil)%ND .gt. NDmax ) NDmax = DoF(icoil)%ND
  enddo

  call bnormal( ideriv )

  bnormmax = zero
  bnormavg = zero
  if ( ideriv .eq. 1 ) t1Bavg = zero

  SALLOCATE(xxhold, ( 1:Ncoils , 0:NSmax) , zero)
  SALLOCATE(xthold, ( 1:Ncoils , 0:NSmax) , zero)
  SALLOCATE(yyhold, ( 1:Ncoils , 0:NSmax) , zero)
  SALLOCATE(ythold, ( 1:Ncoils , 0:NSmax) , zero)
  SALLOCATE(zzhold, ( 1:Ncoils , 0:NSmax) , zero)
  SALLOCATE(zthold, ( 1:Ncoils , 0:NSmax) , zero)
  SALLOCATE(xofhold, ( 1:Ncoils, 0:NSmax-1, 1:NDmax), zero)
  SALLOCATE(yofhold, ( 1:Ncoils, 0:NSmax-1, 1:NDmax), zero)
  SALLOCATE(zofhold, ( 1:Ncoils, 0:NSmax-1, 1:NDmax), zero)
  do icoil = 1, Ncoils
     xxhold(icoil,0:coil(icoil)%NS) = coil(icoil)%xx(0:coil(icoil)%NS)
     xthold(icoil,0:coil(icoil)%NS) = coil(icoil)%xt(0:coil(icoil)%NS)
     yyhold(icoil,0:coil(icoil)%NS) = coil(icoil)%yy(0:coil(icoil)%NS)
     ythold(icoil,0:coil(icoil)%NS) = coil(icoil)%yt(0:coil(icoil)%NS)
     zzhold(icoil,0:coil(icoil)%NS) = coil(icoil)%zz(0:coil(icoil)%NS)
     zthold(icoil,0:coil(icoil)%NS) = coil(icoil)%zt(0:coil(icoil)%NS)
     xofhold(icoil,0:coil(icoil)%NS-1,1:DoF(icoil)%ND) = DoF(icoil)%xof(0:coil(icoil)%NS-1,1:DoF(icoil)%ND)
     yofhold(icoil,0:coil(icoil)%NS-1,1:DoF(icoil)%ND) = DoF(icoil)%yof(0:coil(icoil)%NS-1,1:DoF(icoil)%ND)
     zofhold(icoil,0:coil(icoil)%NS-1,1:DoF(icoil)%ND) = DoF(icoil)%zof(0:coil(icoil)%NS-1,1:DoF(icoil)%ND)
  enddo

  SALLOCATE(g,  ( 1:Ncoils , 0:NSmax) , zero)
  SALLOCATE(gp, ( 1:Ncoils , 0:NSmax) , zero)
  SALLOCATE(dgddof, ( 1:Ncoils, 0:NSmax, 1:NDmax ) , zero)
  SALLOCATE(d1L, (1:1,1:NDmax), zero) 
  do icoil = 1, Ncoils
     call lenDeriv0( icoil, coil(icoil)%L )
     call lenDeriv1( icoil, d1L(1:1,1:DoF(icoil)%ND), DoF(icoil)%ND )
     do j = 1, coil(icoil)%NS
        g(icoil,j) = sqrt( coil(icoil)%xt(j-1)**2 + coil(icoil)%yt(j-1)**2 + coil(icoil)%zt(j-1)**2 )
        g(icoil,j) = g(icoil,j) + sqrt( coil(icoil)%xt(j)**2 + coil(icoil)%yt(j)**2 + coil(icoil)%zt(j)**2 )
        g(icoil,j) = 0.5*g(icoil,j) + g(icoil,j-1)
        
        dgddof(icoil, j, 1:DoF(icoil)%ND) = ( coil(icoil)%xt(j-1)*DoF(icoil)%xtof(j-1,1:DoF(icoil)%ND) + &
             coil(icoil)%yt(j-1)*DoF(icoil)%ytof(j-1,1:DoF(icoil)%ND) + coil(icoil)%zt(j-1)*DoF(icoil)%ztof(j-1,1:DoF(icoil)%ND) ) / &
             sqrt( coil(icoil)%xt(j-1)**2 + coil(icoil)%yt(j-1)**2 + coil(icoil)%zt(j-1)**2 )
        
        dgddof(icoil, j, 1:DoF(icoil)%ND) = dgddof(icoil, j, 1:DoF(icoil)%ND) - &
             d1L(1,1:DoF(icoil)%ND)*sqrt( coil(icoil)%xt(j-1)**2 + coil(icoil)%yt(j-1)**2 + coil(icoil)%zt(j-1)**2 )/coil(icoil)%L

        dgddof(icoil, j, 1:DoF(icoil)%ND) = dgddof(icoil, j, 1:DoF(icoil)%ND) + ( coil(icoil)%xt(j)*DoF(icoil)%xtof(j,1:DoF(icoil)%ND) + &
             coil(icoil)%yt(j)*DoF(icoil)%ytof(j,1:DoF(icoil)%ND) + coil(icoil)%zt(j)*DoF(icoil)%ztof(j,1:DoF(icoil)%ND) ) / sqrt( coil(icoil)%xt(j)**2 + coil(icoil)%yt(j)**2 + coil(icoil)%zt(j)**2 )
        
        dgddof(icoil, j, 1:DoF(icoil)%ND) = dgddof(icoil, j, 1:DoF(icoil)%ND) - &
             d1L(1,1:DoF(icoil)%ND)*sqrt( coil(icoil)%xt(j)**2 + coil(icoil)%yt(j)**2 + coil(icoil)%zt(j)**2 )/coil(icoil)%L
        
        dgddof(icoil, j, 1:DoF(icoil)%ND) = 0.5*dgddof(icoil, j, 1:DoF(icoil)%ND) + dgddof(icoil, j-1, 1:DoF(icoil)%ND)
     enddo
     g(icoil, 0:coil(icoil)%NS) = g(icoil, 0:coil(icoil)%NS)*pi2 / (coil(icoil)%L*real(coil(icoil)%NS))
     gp(icoil,0:coil(icoil)%NS) = sqrt( coil(icoil)%xt(0:coil(icoil)%NS)**2 + coil(icoil)%yt(0:coil(icoil)%NS)**2 + coil(icoil)%zt(0:coil(icoil)%NS)**2 )/coil(icoil)%L
     dgddof(icoil, 0:coil(icoil)%NS, 1:DoF(icoil)%ND) = dgddof(icoil, 0:coil(icoil)%NS, 1:DoF(icoil)%ND)*pi2 / (coil(icoil)%L*real(coil(icoil)%NS))
  enddo

  sdelta = sdelta/sqrt(3.0)

  do j = 1, Npert
     do icoil = 1, Ncoils
        if( coil(icoil)%type .ne. 1 ) cycle
        NS = coil(icoil)%NS
        coil(icoil)%xx(0:NS) = xxhold(icoil,0:NS) + sdelta*cos(pi2*real(coil(icoil)%nxx(j))*g(icoil,0:NS) + coil(icoil)%psx(j))
        coil(icoil)%xt(0:NS) = xthold(icoil,0:NS) - pi2*real(coil(icoil)%nxx(j))*sdelta*sin(pi2*real(coil(icoil)%nxx(j))*g(icoil,0:NS) + &
                                                    coil(icoil)%psx(j))*gp(icoil,0:NS)

        coil(icoil)%yy(0:NS) = yyhold(icoil,0:NS) + sdelta*cos(pi2*real(coil(icoil)%nyy(j))*g(icoil,0:NS) + coil(icoil)%psy(j))
        coil(icoil)%yt(0:NS) = ythold(icoil,0:NS) - pi2*real(coil(icoil)%nyy(j))*sdelta*sin(pi2*real(coil(icoil)%nyy(j))*g(icoil,0:NS) + &
                                                    coil(icoil)%psy(j))*gp(icoil,0:NS)
        
        coil(icoil)%zz(0:NS) = zzhold(icoil,0:NS) + sdelta*cos(pi2*real(coil(icoil)%nzz(j))*g(icoil,0:NS) + coil(icoil)%psz(j))
        coil(icoil)%zt(0:NS) = zthold(icoil,0:NS) - pi2*real(coil(icoil)%nzz(j))*sdelta*sin(pi2*real(coil(icoil)%nzz(j))*g(icoil,0:NS) + &
                                                    coil(icoil)%psz(j))*gp(icoil,0:NS)

        do jj = 1,DoF(icoil)%ND
           DoF(icoil)%xof(0:NS-1,jj) = xofhold(icoil,0:NS-1,jj) - & ! dx/ddof
                pi2*real(coil(icoil)%nxx(j))*sdelta*sin(pi2*real(coil(icoil)%nxx(j))*g(icoil,0:NS-1) + coil(icoil)%psx(j))*dgddof(icoil,0:NS-1,jj)
           DoF(icoil)%yof(0:NS-1,jj) = yofhold(icoil,0:NS-1,jj) - &
                pi2*real(coil(icoil)%nyy(j))*sdelta*sin(pi2*real(coil(icoil)%nyy(j))*g(icoil,0:NS-1) + coil(icoil)%psy(j))*dgddof(icoil,0:NS-1,jj)
           DoF(icoil)%zof(0:NS-1,jj) = zofhold(icoil,0:NS-1,jj) - &
                pi2*real(coil(icoil)%nzz(j))*sdelta*sin(pi2*real(coil(icoil)%nzz(j))*g(icoil,0:NS-1) + coil(icoil)%psz(j))*dgddof(icoil,0:NS-1,jj)
        enddo
     enddo
     call bnormal(ideriv)
     ! Could save every evalutation
     bnormavg = bnormavg + bnorm
     if ( ideriv .eq. 1 ) t1Bavg(1:Ndof) = t1Bavg(1:Ndof) + t1B(1:Ndof)
     if ( bnormmax .le. bnorm ) bnormmax = bnorm
  enddo

  sdelta = sdelta*sqrt(3.0)
  bnormavg = bnormavg / real(Npert)
  if ( ideriv .eq. 1 ) t1Bavg = t1Bavg / real(Npert)

  do icoil = 1, Ncoils
     coil(icoil)%xx(0:coil(icoil)%NS) = xxhold(icoil,0:coil(icoil)%NS)
     coil(icoil)%xt(0:coil(icoil)%NS) = xthold(icoil,0:coil(icoil)%NS)
     coil(icoil)%yy(0:coil(icoil)%NS) = yyhold(icoil,0:coil(icoil)%NS)
     coil(icoil)%yt(0:coil(icoil)%NS) = ythold(icoil,0:coil(icoil)%NS)
     coil(icoil)%zz(0:coil(icoil)%NS) = zzhold(icoil,0:coil(icoil)%NS)
     coil(icoil)%zt(0:coil(icoil)%NS) = zthold(icoil,0:coil(icoil)%NS)
     DoF(icoil)%xof(0:coil(icoil)%NS-1,1:DoF(icoil)%ND) = xofhold(icoil,0:coil(icoil)%NS-1,1:DoF(icoil)%ND)
     DoF(icoil)%yof(0:coil(icoil)%NS-1,1:DoF(icoil)%ND) = yofhold(icoil,0:coil(icoil)%NS-1,1:DoF(icoil)%ND)
     DoF(icoil)%zof(0:coil(icoil)%NS-1,1:DoF(icoil)%ND) = zofhold(icoil,0:coil(icoil)%NS-1,1:DoF(icoil)%ND)
  enddo
  call bnormal(ideriv)

  DALLOCATE(xxhold)
  DALLOCATE(xthold)
  DALLOCATE(yyhold)
  DALLOCATE(ythold)
  DALLOCATE(zzhold)
  DALLOCATE(zthold)
  DALLOCATE(g)
  DALLOCATE(gp)
  DALLOCATE(dgddof)
  DALLOCATE(d1L)
  DALLOCATE(xofhold)
  DALLOCATE(yofhold)
  DALLOCATE(zofhold)

  return
end subroutine sbnormal

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine perturbation(ideriv)
 
  use globals, only: dp, zero, half, pi2, machprec, ncpu, myid, ounit, &
  coil, DoF, Ncoils, Nfixgeo, Ndof, bnorm, Npert, Nmax, MPI_COMM_FOCUS
  
  use mpi
 
  implicit none
  INTEGER, INTENT(in) :: ideriv
  INTEGER             :: astat, ierr, icoil, j
  REAL                :: arb
 
  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
 
  ! Stochastic variables set here
  do icoil = 1, Ncoils
     SALLOCATE( coil(icoil)%psx, (1:Npert), zero )
     SALLOCATE( coil(icoil)%psy, (1:Npert), zero )
     SALLOCATE( coil(icoil)%psz, (1:Npert), zero )
     SALLOCATE( coil(icoil)%nxx, (1:Npert), 0 )
     SALLOCATE( coil(icoil)%nyy, (1:Npert), 0 )
     SALLOCATE( coil(icoil)%nzz, (1:Npert), 0 )
     if (myid .eq. 0) then
        do j = 1, Npert
           call random_number(arb)
           !arb = rand(0)                ! Used for deterministic debugging 
           coil(icoil)%psx(j) = pi2*arb
           call random_number(arb)
           coil(icoil)%nxx(j) = FLOOR(real(Nmax)*arb) + 1 ! Gives int between 1 and 3
           call random_number(arb)
           coil(icoil)%psy(j) = pi2*arb
           call random_number(arb)
           coil(icoil)%nyy(j) = FLOOR(real(Nmax)*arb) + 1
           call random_number(arb)
           coil(icoil)%psz(j) = pi2*arb
           call random_number(arb)
           coil(icoil)%nzz(j) = FLOOR(real(Nmax)*arb) + 1
        enddo
     endif
     RlBCAST( coil(icoil)%psx, Npert, 0 )
     IlBCAST( coil(icoil)%nxx, Npert, 0 )
     RlBCAST( coil(icoil)%psy, Npert, 0 )
     IlBCAST( coil(icoil)%nyy, Npert, 0 )
     RlBCAST( coil(icoil)%psz, Npert, 0 )
     IlBCAST( coil(icoil)%nzz, Npert, 0 )
  enddo
  
end subroutine perturbation
