!title (tlength) ! Calculate coil variations and field error from variations.. (tkruger)

!latex \briefly{Write something...}

!latex \calledby{\link{solvers}}

!latex  \section{General}
!latex  Write something...

subroutine sbnormal(ideriv)
  use globals, only: dp, zero, half, pi2, machprec, ncpu, myid, ounit, &
       coil, DoF, Ncoils, Nfixgeo, Ndof, bnorm, t1B, t1R, t1Bavg, t1Ravg, resbn, Npert, sdelta, Nmax, &
       bnormavg, bnormmax, resbnavg, abspsimax, surf, plasma, DoF, weight_sbnorm, weight_sresbn, &
       rcflux_target, psi, stochpsi, stochpsipred, case_optimize, CG_maxiter, stoch_done, MPI_COMM_FOCUS

  use mpi
  implicit none
  INTEGER, INTENT(in) :: ideriv
  INTEGER             :: astat, ierr, icoil, j, idof, ivec, NSmax, NDmax, jj, NS, failure
  REAL                :: psipred, psi0, dummy
  REAL, allocatable   :: xxhold(:,:), xthold(:,:), yyhold(:,:), ythold(:,:), &
                         zzhold(:,:), zthold(:,:), g(:,:), gp(:,:), dgddof(:,:,:), &
                         d1L(:,:), xofhold(:,:,:), yofhold(:,:,:), zofhold(:,:,:)

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( stoch_done .eq. 1) return
  if( case_optimize .eq. 0 .or. CG_maxiter .eq. 0 ) then
     stoch_done = 1
  endif
  
  NSmax = zero
  NDmax = zero 
  do icoil = 1, Ncoils
     if ( coil(icoil)%NS .gt. NSmax ) NSmax = coil(icoil)%NS
     if ( DoF(icoil)%ND .gt. NDmax ) NDmax = DoF(icoil)%ND
  enddo

  if ( weight_sbnorm .ne. 0.0 ) then
     call bnormal( ideriv )
     if ( ideriv .eq. 1 ) t1Bavg = zero
  endif
  if ( weight_sresbn .ne. 0.0 ) then
     call rcflux( ideriv )
     if ( ideriv .eq. 1 ) t1Ravg = zero
  endif
  psi0 = psi

  bnormmax = zero
  bnormavg = zero
  abspsimax = zero
  resbnavg = zero
  failure = zero

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
        ! Use global pertx ...
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
     
     if ( weight_sbnorm .ne. 0.0 ) then
        call bnormal(ideriv)
        bnormavg = bnormavg + bnorm
        if ( bnormmax .le. bnorm ) bnormmax = bnorm
        if ( ideriv .eq. 1 ) t1Bavg(1:Ndof) = t1Bavg(1:Ndof) + t1B(1:Ndof)
     endif
     
     if ( weight_sresbn .ne. 0.0 ) then
        call rcflux(ideriv)
        if ( resbn .ne. 1.0 ) then
           psipred = psi0
           do icoil = 1, Ncoils
              NS = coil(icoil)%NS
              call dpsi0(icoil,dummy)
              ! Use coil symmetry instead of 8
              psipred = psipred + 8.0*sum(coil(icoil)%dpsidx(1:NS)*sdelta*cos(pi2*real(coil(icoil)%nxx(j))*g(icoil,1:NS)+coil(icoil)%psx(j)) + &
                                          coil(icoil)%dpsidy(1:NS)*sdelta*cos(pi2*real(coil(icoil)%nyy(j))*g(icoil,1:NS)+coil(icoil)%psy(j)) + &
                                          coil(icoil)%dpsidz(1:NS)*sdelta*cos(pi2*real(coil(icoil)%nzz(j))*g(icoil,1:NS)+coil(icoil)%psz(j)))
           enddo
           resbnavg = resbnavg + (psi-rcflux_target)**2.0
           if ( ideriv .eq. 1 ) t1Ravg(1:Ndof) = t1Ravg(1:Ndof) + t1R(1:Ndof)
           stochpsi(j) = psi
           stochpsipred(j) = psipred
           if ( abspsimax .le. abs(psi) ) abspsimax = abs(psi)
        else
           failure = failure + 1
           stochpsi(j) = HUGE(psi)
           stochpsipred(j) = HUGE(psi)
        endif
     endif
  enddo

  sdelta = sdelta*sqrt(3.0)

  if ( weight_sbnorm .ne. 0.0 ) then
     bnormavg = bnormavg / real(Npert)
     if ( ideriv .eq. 1 ) t1Bavg = t1Bavg / real(Npert)
  endif  

  if ( weight_sresbn .ne. 0.0 ) then
     if ( Npert-failure .ne. 0 ) then
        resbnavg = resbnavg / real(Npert-failure)
        if ( ideriv .eq. 1 ) t1Ravg = t1Ravg / real(Npert-failure)
     endif
  endif
     
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
  if ( weight_sbnorm .ne. 0.0 ) then
     call bnormal(ideriv)
  endif
  if ( weight_sresbn .ne. 0.0 ) then
     call rcflux(ideriv)
  endif

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
  coil, DoF, Ncoils, Nfixgeo, Ndof, bnorm, Npert, Nmax, stochpsi, stochpsipred, &
  sdelta, Nturns, Npancakes, coil_type_multi, MPI_COMM_FOCUS
  
  use mpi
 
  implicit none
  INTEGER, INTENT(in) :: ideriv
  INTEGER             :: astat, ierr, icoil, i, j, NS
  REAL                :: arb, freq
  REAL, allocatable   :: g(:), gp(:)
 
  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
 
  ! Stochastic variables set here
  SALLOCATE( stochpsi, (1:Npert), zero )
  SALLOCATE( stochpsipred, (1:Npert), zero )
  do icoil = 1, Ncoils
     SALLOCATE( coil(icoil)%psx, (1:Npert), zero )
     SALLOCATE( coil(icoil)%psy, (1:Npert), zero )
     SALLOCATE( coil(icoil)%psz, (1:Npert), zero )
     SALLOCATE( coil(icoil)%nxx, (1:Npert), 0 )
     SALLOCATE( coil(icoil)%nyy, (1:Npert), 0 )
     SALLOCATE( coil(icoil)%nzz, (1:Npert), 0 )
     if (myid .eq. 0) then
        do j = 1, Npert
           !arb = rand(0) ! Deterministic
           call random_number(arb)
           coil(icoil)%psx(j) = pi2*arb
           call random_number(arb)
           coil(icoil)%nxx(j) = FLOOR(real(Nmax)*arb) + 1
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
     
     NS = coil(icoil)%NS
     SALLOCATE( coil(icoil)%pertx, (0:NS,1:Npert), 0.0 )
     SALLOCATE( coil(icoil)%perty, (0:NS,1:Npert), 0.0 )
     SALLOCATE( coil(icoil)%pertz, (0:NS,1:Npert), 0.0 )
     SALLOCATE( coil(icoil)%pertxp, (0:NS,1:Npert), 0.0 )
     SALLOCATE( coil(icoil)%pertyp, (0:NS,1:Npert), 0.0 )
     SALLOCATE( coil(icoil)%pertzp, (0:NS,1:Npert), 0.0 )
    
     if ( coil(icoil)%type .eq. 1 ) then
        freq = 1.0
     else if ( coil(icoil)%type .eq. coil_type_multi ) then
        freq = real(Nturns*Npancakes) ! Need to use derivative of single-filament
     else
        cycle
     endif
     sdelta = sdelta/sqrt(3.0)
     SALLOCATE( g, (0:NS), 0.0 )
     SALLOCATE( gp, (0:NS), 0.0 )
     do i = 1, NS
        g(i) = sqrt( coil(icoil)%xt(i-1)**2 + coil(icoil)%yt(i-1)**2 + coil(icoil)%zt(i-1)**2 )
        g(i) = g(i) + sqrt( coil(icoil)%xt(i)**2 + coil(icoil)%yt(i)**2 + coil(icoil)%zt(i)**2 )
        g(i) = 0.5*g(i) + g(i-1)
     enddo
     call lenDeriv0( icoil, coil(icoil)%L )
     g(0:NS) = g(0:NS)*pi2 / ( coil(icoil)%L*real(NS) )
     gp(0:NS) = sqrt( coil(icoil)%xt(0:NS)**2 + coil(icoil)%yt(0:NS)**2 + coil(icoil)%zt(0:NS)**2 ) / ( coil(icoil)%L )
     do j = 1, Npert
        coil(icoil)%pertx(0:NS,j) = sdelta*cos( freq*pi2*real(coil(icoil)%nxx(j))*g(0:NS) + coil(icoil)%psx(j) )
        coil(icoil)%perty(0:NS,j) = sdelta*cos( freq*pi2*real(coil(icoil)%nyy(j))*g(0:NS) + coil(icoil)%psy(j) )
        coil(icoil)%pertz(0:NS,j) = sdelta*cos( freq*pi2*real(coil(icoil)%nzz(j))*g(0:NS) + coil(icoil)%psz(j) )
        coil(icoil)%pertxp(0:NS,j) = -1.0*sdelta*sin(freq*pi2*real(coil(icoil)%nxx(j))*g(0:NS)+coil(icoil)%psx(j))*freq*pi2*real(coil(icoil)%nxx(j))*gp(0:NS)
        coil(icoil)%pertyp(0:NS,j) = -1.0*sdelta*sin(freq*pi2*real(coil(icoil)%nyy(j))*g(0:NS)+coil(icoil)%psy(j))*freq*pi2*real(coil(icoil)%nyy(j))*gp(0:NS)
        coil(icoil)%pertzp(0:NS,j) = -1.0*sdelta*sin(freq*pi2*real(coil(icoil)%nzz(j))*g(0:NS)+coil(icoil)%psz(j))*freq*pi2*real(coil(icoil)%nzz(j))*gp(0:NS)
     enddo
     DALLOCATE( g )
     DALLOCATE( gp )
     sdelta = sdelta*sqrt(3.0)

  enddo

end subroutine perturbation
