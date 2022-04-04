!title (tlength) ! Calculate coil variations and field error from variations.. (tkruger)

!latex \briefly{Write something...}

!latex \calledby{\link{solvers}}

!latex  \section{General}
!latex  Write something...

subroutine sbnormal(ideriv)
  use globals, only: dp, zero, half, pi2, machprec, ncpu, myid, ounit, &
       coil, DoF, Ncoils, Nfixgeo, Ndof, bnorm, t1B, t1Bavg, Npert, sdelta, Nmax, &
       bnormavg, bnormmax, MPI_COMM_FOCUS

  use mpi
  implicit none
  INTEGER, INTENT(in) :: ideriv

  INTEGER             :: astat, ierr, icoil, j, idof, ND, ivec, nxx(1:Ncoils,1:Npert), nxxhold, &
                         nyy(1:Ncoils,1:Npert), nyyhold, nzz(1:Ncoils,1:Npert), nzzhold
  REAL                :: d1L(1:Ndof), bnormhold, &
                         pertxx(0:coil(1)%NS-1), pertxt(0:coil(1)%NS-1), psx(1:Ncoils,1:Npert), &
                         pertyy(0:coil(1)%NS-1), pertyt(0:coil(1)%NS-1), psy(1:Ncoils,1:Npert), &
                         pertzz(0:coil(1)%NS-1), pertzt(0:coil(1)%NS-1), psz(1:Ncoils,1:Npert)
  REAL, allocatable   :: t1Bhold(:), xxhold(:,:), xthold(:,:), yyhold(:,:), ythold(:,:), &
                         zzhold(:,:), zthold(:,:), g(:,:), gp(:,:)

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  do icoil = 1, Ncoils
     psx(icoil,1:Npert) = coil(icoil)%psx(1:Npert)
     psy(icoil,1:Npert) = coil(icoil)%psy(1:Npert)
     psz(icoil,1:Npert) = coil(icoil)%psz(1:Npert)
     nxx(icoil,1:Npert) = coil(icoil)%nxx(1:Npert)
     nyy(icoil,1:Npert) = coil(icoil)%nyy(1:Npert)
     nzz(icoil,1:Npert) = coil(icoil)%nzz(1:Npert)
  enddo

  if ( ideriv .eq. 0 ) then
     call bnormal( 0 )
  else
     call bnormal( 1 )
  endif
  bnormhold = bnorm
  t1Bhold = t1B
  bnormmax = zero
  bnormavg = zero
  t1Bavg = zero

  SALLOCATE(xxhold, ( 1:Ncoils , 0:coil(1)%NS-1) , zero)
  SALLOCATE(xthold, ( 1:Ncoils , 0:coil(1)%NS-1) , zero)
  SALLOCATE(yyhold, ( 1:Ncoils , 0:coil(1)%NS-1) , zero)
  SALLOCATE(ythold, ( 1:Ncoils , 0:coil(1)%NS-1) , zero)
  SALLOCATE(zzhold, ( 1:Ncoils , 0:coil(1)%NS-1) , zero)
  SALLOCATE(zthold, ( 1:Ncoils , 0:coil(1)%NS-1) , zero)
  do icoil = 1, Ncoils
     if( coil(icoil)%NS .ne. coil(1)%NS ) then
        if (myid .eq. 0) write(ounit, '(8X"WARNING: Stochastic needs NS constant")')
     endif
     xxhold(icoil,0:coil(1)%NS-1) = coil(icoil)%xx(0:coil(1)%NS-1)
     xthold(icoil,0:coil(1)%NS-1) = coil(icoil)%xt(0:coil(1)%NS-1)
     yyhold(icoil,0:coil(1)%NS-1) = coil(icoil)%yy(0:coil(1)%NS-1)
     ythold(icoil,0:coil(1)%NS-1) = coil(icoil)%yt(0:coil(1)%NS-1)
     zzhold(icoil,0:coil(1)%NS-1) = coil(icoil)%zz(0:coil(1)%NS-1)
     zthold(icoil,0:coil(1)%NS-1) = coil(icoil)%zt(0:coil(1)%NS-1)
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

  sdelta = sdelta/sqrt(3.0)

  do j = 1, Npert
     do icoil = 1, Ncoils
        if( coil(icoil)%type .ne. 1 ) cycle
        !pertxx(0:coil(1)%NS-1) = sdelta*cos(pi2*nxx(icoil,j)*g(icoil,0:coil(1)%NS-1) + psx(icoil,j))
        !pertxt(0:coil(1)%NS-1) = -1.0*pi2*nxx(icoil,j)*sdelta*sin(pi2*nxx(icoil,j)*g(icoil,0:coil(1)%NS-1) + psx(icoil,j))*gp(icoil,0:coil(1)%NS-1)
        !coil(icoil)%xx(0:coil(1)%NS-1) = xxhold(icoil,0:coil(1)%NS-1)
        !coil(icoil)%xt(0:coil(1)%NS-1) = xthold(icoil,0:coil(1)%NS-1)
        !coil(icoil)%xx(0:coil(1)%NS-1) = coil(icoil)%xx(0:coil(1)%NS-1) + pertxx(0:coil(1)%NS-1)
        !coil(icoil)%xt(0:coil(1)%NS-1) = coil(icoil)%xt(0:coil(1)%NS-1) + pertxt(0:coil(1)%NS-1)
        coil(icoil)%xx(0:coil(1)%NS-1) = xxhold(icoil,0:coil(1)%NS-1) + sdelta*cos(pi2*nxx(icoil,j)*g(icoil,0:coil(1)%NS-1) + psx(icoil,j))
        coil(icoil)%xt(0:coil(1)%NS-1) = xthold(icoil,0:coil(1)%NS-1) - pi2*nxx(icoil,j)*sdelta*sin(pi2*nxx(icoil,j)*g(icoil,0:coil(1)%NS-1) + psx(icoil,j))*gp(icoil,0:coil(1)%NS-1)

        !pertyy(0:coil(1)%NS-1) = sdelta*cos(pi2*nyy(icoil,j)*g(icoil,0:coil(1)%NS-1) + psy(icoil,j))
        !pertyt(0:coil(1)%NS-1) = -1.0*pi2*nyy(icoil,j)*sdelta*sin(pi2*nyy(icoil,j)*g(icoil,0:coil(1)%NS-1) + psy(icoil,j))*gp(icoil,0:coil(1)%NS-1)
        !coil(icoil)%yy(0:coil(1)%NS-1) = yyhold(icoil,0:coil(1)%NS-1)
        !coil(icoil)%yt(0:coil(1)%NS-1) = ythold(icoil,0:coil(1)%NS-1)
        !coil(icoil)%yy(0:coil(1)%NS-1) = coil(icoil)%yy(0:coil(1)%NS-1) + pertyy(0:coil(1)%NS-1)
        !coil(icoil)%yt(0:coil(1)%NS-1) = coil(icoil)%yt(0:coil(1)%NS-1) + pertyt(0:coil(1)%NS-1)
        coil(icoil)%yy(0:coil(1)%NS-1) = yyhold(icoil,0:coil(1)%NS-1) + sdelta*cos(pi2*nyy(icoil,j)*g(icoil,0:coil(1)%NS-1) + psy(icoil,j))
        coil(icoil)%yt(0:coil(1)%NS-1) = ythold(icoil,0:coil(1)%NS-1) - pi2*nyy(icoil,j)*sdelta*sin(pi2*nyy(icoil,j)*g(icoil,0:coil(1)%NS-1) + psy(icoil,j))*gp(icoil,0:coil(1)%NS-1)
        
        !pertzz(0:coil(1)%NS-1) = sdelta*cos(pi2*nzz(icoil,j)*g(icoil,0:coil(1)%NS-1) + psz(icoil,j))
        !pertzt(0:coil(1)%NS-1) = -1.0*pi2*nzz(icoil,j)*sdelta*sin(pi2*nzz(icoil,j)*g(icoil,0:coil(1)%NS-1) + psz(icoil,j))*gp(icoil,0:coil(1)%NS-1)
        !coil(icoil)%zz(0:coil(1)%NS-1) = zzhold(icoil,0:coil(1)%NS-1)
        !coil(icoil)%zt(0:coil(1)%NS-1) = zthold(icoil,0:coil(1)%NS-1)
        !coil(icoil)%zz(0:coil(1)%NS-1) = coil(icoil)%zz(0:coil(1)%NS-1) + pertzz(0:coil(1)%NS-1)
        !coil(icoil)%zt(0:coil(1)%NS-1) = coil(icoil)%zt(0:coil(1)%NS-1) + pertzt(0:coil(1)%NS-1)
        coil(icoil)%zz(0:coil(1)%NS-1) = zzhold(icoil,0:coil(1)%NS-1) + sdelta*cos(pi2*nzz(icoil,j)*g(icoil,0:coil(1)%NS-1) + psz(icoil,j))
        coil(icoil)%zt(0:coil(1)%NS-1) = zthold(icoil,0:coil(1)%NS-1) - pi2*nzz(icoil,j)*sdelta*sin(pi2*nzz(icoil,j)*g(icoil,0:coil(1)%NS-1) + psz(icoil,j))*gp(icoil,0:coil(1)%NS-1)
     enddo
     !call MPI_BARRIER( MPI_COMM_FOCUS, ierr )
     if ( ideriv .eq. 0 ) then
        call bnormal( 0 )
     else
        call bnormal( 1 )
        t1Bavg = t1Bavg + t1B
     endif
     if ( bnormmax .le. bnorm ) bnormmax = bnorm
     bnormavg = bnormavg + bnorm
  enddo

  sdelta = sdelta*sqrt(3.0)
  bnormavg = bnormavg / real(Npert)
  if ( ideriv .eq. 1 ) t1Bavg = t1Bavg / real(Npert)
  bnorm = bnormhold
  if ( ideriv .eq. 1 ) t1B = t1Bhold

  do icoil = 1, Ncoils
     coil(icoil)%xx(0:coil(1)%NS-1) = xxhold(icoil,0:coil(1)%NS-1)
     coil(icoil)%xt(0:coil(1)%NS-1) = xthold(icoil,0:coil(1)%NS-1)
     coil(icoil)%yy(0:coil(1)%NS-1) = yyhold(icoil,0:coil(1)%NS-1)
     coil(icoil)%yt(0:coil(1)%NS-1) = ythold(icoil,0:coil(1)%NS-1)
     coil(icoil)%zz(0:coil(1)%NS-1) = zzhold(icoil,0:coil(1)%NS-1)
     coil(icoil)%zt(0:coil(1)%NS-1) = zthold(icoil,0:coil(1)%NS-1)
  enddo

  DALLOCATE(xxhold)
  DALLOCATE(xthold)
  DALLOCATE(yyhold)
  DALLOCATE(ythold)
  DALLOCATE(zzhold)
  DALLOCATE(zthold)
  DALLOCATE(g)
  DALLOCATE(gp)

  return
end subroutine sbnormal

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine perturbation(ideriv)
 
  use globals, only: dp, zero, half, pi2, machprec, ncpu, myid, ounit, &
  coil, DoF, Ncoils, Nfixgeo, Ndof, bnorm, Npert, sdelta, Nmax, MPI_COMM_FOCUS
  
  use mpi
 
  implicit none
  INTEGER, INTENT(in) :: ideriv
  INTEGER             :: astat, ierr, icoil, j
  REAL                :: arb
 
  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
 
  ! Stochastic variables set here, not neccessarily parameterization independent
  do icoil = 1, Ncoils
     SALLOCATE( coil(icoil)%psx, (1:Npert), zero )
     SALLOCATE( coil(icoil)%psy, (1:Npert), zero )
     SALLOCATE( coil(icoil)%psz, (1:Npert), zero )
     SALLOCATE( coil(icoil)%nxx, (1:Npert), zero )
     SALLOCATE( coil(icoil)%nyy, (1:Npert), zero )
     SALLOCATE( coil(icoil)%nzz, (1:Npert), zero )
     do j = 1, Npert
        if (myid .eq. 0) then
           call random_number(arb)
           coil(icoil)%psx(j) = pi2*arb
           call random_number(arb)
           coil(icoil)%nxx(j) = FLOOR(real(Nmax)*arb) + 1 ! Should give int between 1 and 3
           call random_number(arb)
           coil(icoil)%psy(j) = pi2*arb
           call random_number(arb)
           coil(icoil)%nyy(j) = FLOOR(real(Nmax)*arb) + 1
           call random_number(arb)
           coil(icoil)%psz(j) = pi2*arb
           call random_number(arb)
           coil(icoil)%nzz(j) = FLOOR(real(Nmax)*arb) + 1
        endif
        call MPI_BARRIER( MPI_COMM_FOCUS, ierr )
        RlBCAST( coil(icoil)%psx(j), 1, 0 )
        IlBCAST( coil(icoil)%nxx(j), 1, 0 )
        RlBCAST( coil(icoil)%psy(j), 1, 0 )
        IlBCAST( coil(icoil)%nyy(j), 1, 0 )
        RlBCAST( coil(icoil)%psz(j), 1, 0 )
        IlBCAST( coil(icoil)%nzz(j), 1, 0 )
        call MPI_BARRIER( MPI_COMM_FOCUS, ierr )
     enddo
  enddo

end subroutine perturbation
