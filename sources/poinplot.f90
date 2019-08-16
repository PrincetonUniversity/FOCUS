SUBROUTINE poinplot
  !------------------------------------------------------------------------------------------------------ 
  ! DATE:  12/12/2018
  ! Poincare plots of the vacuum flux surfaces and calculate the rotational transform
  !------------------------------------------------------------------------------------------------------ 
  USE globals, only : dp, myid, ncpu, zero, half, pi, pi2, ounit, pi, sqrtmachprec, pp_maxiter, &
                      pp_phi, pp_raxis, pp_zaxis, pp_xtol, pp_rmax, pp_zmax, ppr, ppz, pp_ns, iota, nfp_raw, &
                      XYZB, lboozmn, booz_mnc, booz_mns, booz_mn, total_num, &
                      master, nmaster, nworker, masterid, color, myworkid, MPI_COMM_MASTERS, MPI_COMM_MYWORLD, MPI_COMM_WORKERS
  USE mpi
  IMPLICIT NONE

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  INTEGER              :: ierr, astat, iflag
  INTEGER              :: ip, is, niter, icommand
  REAL                 :: theta, zeta, r, RZ(2), r1, z1, rzrzt(5), x, y, z
  REAL                 :: B(3), start, finish

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  FATAL( poinplot, pp_ns < 1    , not enough starting points )
  FATAL( poinplot, pp_maxiter<1 , not enough max. iterations )

  pp_phi = pp_phi * pi  ! pp_phi=0.5 -> pi/2
  x = zero ; y = zero ; z = zero ; B = zero

  ! if raxis, zaxis not provided
  if ( (abs(pp_raxis) + abs(pp_zaxis)) < sqrtmachprec) then
     zeta = pp_phi
     theta = zero ; call surfcoord( theta, zeta, r , z )
     theta = pi   ; call surfcoord( theta, zeta, r1, z1)
     
     pp_raxis = (r+r1)*half
     pp_zaxis = (z+z1)*half
  endif
  
  
  RZ(1) = pp_raxis ; RZ(2) = pp_zaxis
  start = MPI_Wtime()
  call find_axis(RZ, pp_maxiter, pp_xtol)
  pp_raxis = RZ(1) ; pp_zaxis = RZ(2)
  RlBCAST( pp_raxis, 1 , 0 )
  RlBCAST( pp_zaxis, 1 , 0 )

  ! poincare plot and calculate iota
  SALLOCATE( ppr , (1:pp_ns, 0:pp_maxiter), zero )
  SALLOCATE( ppz , (1:pp_ns, 0:pp_maxiter), zero )
  SALLOCATE( iota, (1:pp_ns)              , zero )

  ! if pp_rmax and pp_zmax not provied 
  if ( (abs(pp_rmax) + abs(pp_zmax)) < sqrtmachprec) then
     zeta = pp_phi
     theta = zero ; call surfcoord( theta, zeta, r , z )
     pp_rmax = r*1.0 ; pp_zmax = z*1.0
  endif

  if(myid==0) write(ounit, '("poinplot: following fieldlines between ("ES12.5 &
       ","ES12.5" ) and ("ES12.5","ES12.5" )")') pp_raxis, pp_zaxis, pp_rmax, pp_zmax

  do is = 1, pp_ns     ! pp_ns is the number of eavaluation surfaces
     niter = 0    ! number of successful iterations
     !if ( modulo(myid, pp_ns) /= modulo((is-1), nmaster))  cycle  ! MPI
     rzrzt(1:5) = (/ pp_raxis + is*(pp_rmax-pp_raxis)/pp_ns, &
          pp_zaxis + is*(pp_zmax-pp_zaxis)/pp_ns, &
          pp_raxis, pp_zaxis, zero                   /)
     ppr(is, 0) = rzrzt(1) ; ppz(is, 0) = rzrzt(2)

     do ip = 1, pp_maxiter
        iflag = 1
        call ppiota(rzrzt, iflag)
        if (iflag >= 0) niter = niter + 1   ! counting
        ppr(is, ip) = rzrzt(1)
        ppz(is, ip) = rzrzt(2)
        ! FATAL( poinplot, abs((rzrzt(3)-pp_raxis)/pp_raxis)>pp_xtol, magnetic axis is not coming back )
     enddo

     if (niter==0) then
        iota(is) = zero
     else 
        iota(is) = rzrzt(5) / (niter*pi2/Nfp_raw)
     endif

     if (myid == 0) write(ounit, '(8X": order="I6" ; masterid="I6" ; (R,Z)=("ES12.5","ES12.5 & 
          " ) ; iota="ES12.5" ; niter="I6" .")') is, myid, ppr(is,0), ppz(is,0), iota(is), niter

     if(lboozmn .and. abs(iota(is))>sqrtmachprec) then
        x = ppr(is, 0) * cos(pp_phi) ; y = ppr(is, 0) * sin(pp_phi) ; z = ppz(is, 0)
        call boozsurf( XYZB(1:total_num, 1:4, is), x, y, z, iota(is), is)
     endif
  enddo

return

END SUBROUTINE poinplot

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

SUBROUTINE find_axis(RZ, MAXFEV, XTOL)
  USE globals, only : dp, myid, ounit, zero, pp_phi, Nfp_raw
  USE mpi
  IMPLICIT NONE

  REAL, INTENT(INOUT)  :: RZ(2)
  REAL, INTENT(IN   )  :: XTOL
  INTEGER, INTENT(IN)  :: MAXFEV

  INTEGER, parameter   :: n=2
  INTEGER              :: ml,mu,mode,nprint,info,nfev,ldfjac,lr
  REAL     :: epsfcn,factor
  REAL     :: fvec(n),diag(n),qtf(n),wa1(n),wa2(n),wa3(n),wa4(n)
  REAL, allocatable :: fjac(:,:),r(:)
  external :: axis_fcn
 
  LR = N*(N+1)/2
  LDFJAC = N
  ml = n-1
  mu = n-1
  epsfcn = 1.0E-4
  mode = 1
  factor = 100.0
  nprint = -1
  
  allocate(fjac(ldfjac,n))
  allocate(r(lr))

  call hybrd(axis_fcn,n,RZ,fvec,xtol,maxfev,ml,mu,epsfcn,diag, &
       mode,factor,nprint,info,nfev,fjac,ldfjac,r,lr,qtf,wa1,wa2,wa3,wa4)

  if (myid == 0) then
     write(ounit,'("findaxis: Finding axis at phi = "ES12.5" with (R,Z) = ( "ES12.5,","ES12.5" ).")') &
          pp_phi, RZ(1), RZ(2)
     select case (info)
     case (0)
        write(ounit,'("findaxis: info=0, improper input parameters.")')
     case (1)
        write(ounit,'("findaxis: info=1, relative error between two consecutive iterates is at most xtol.")')
     case (2)
        write(ounit,'("findaxis: info=2, number of calls to fcn has reached or exceeded maxfev.")')
     case (3)
        write(ounit,'("findaxis: info=3, xtol is too small.")')    
     case (4)
        write(ounit,'("findaxis: info=4, iteration is not making good progress, jacobian.")')    
     case (5)
        write(ounit,'("findaxis: info=5, iteration is not making good progress, function.")')
     case default
        write(ounit,'("findaxis: info="I2", something wrong with the axis finding subroutine.")') info
     end select
  endif

  return

END SUBROUTINE find_axis

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

SUBROUTINE axis_fcn(n,x,fvec,iflag)
  USE globals, only : dp, myid, IsQuiet, ounit, zero, pi2, sqrtmachprec, pp_phi, Nfp_raw, pp_xtol
  USE mpi
  IMPLICIT NONE

  INTEGER  :: n, iflag
  REAL :: x(n), fvec(n)

  INTEGER  :: iwork(5), ierr, ifail
  REAL     :: rz_end(n), phi_init, phi_stop, relerr, abserr, work(100+21*N)
  EXTERNAL :: BRpZ
  
  ifail = 1
  relerr = pp_xtol
  abserr = sqrtmachprec
  phi_init = pp_phi
  phi_stop = pp_phi + pi2/Nfp_raw
  rz_end = x

  call ode( BRpZ, n, rz_end, phi_init, phi_stop, relerr, abserr, ifail, work, iwork )
  if ( ifail /= 2 ) then     
     if ( myid==0 .and. IsQuiet < 0 ) then
        write ( ounit, '(A,I3)' ) 'axis_fcn: ODE solver ERROR; returned IFAIL = ', ifail
        select case ( ifail )
        case ( 3 )
           write(ounit, '("axis_fcn: DF_xtol or abserr too small.")')
        case ( 4 )
           write(ounit, '("axis_fcn: tau not reached after 500 steps.")')
        case ( 5 )
           write(ounit, '("axis_fcn: tau not reached because equation to be stiff.")')
        case ( 6 )
           write(ounit, '("axis_fcn: INVALID input parameters.")')
        end select
     end if
     iflag = -1
     ! call MPI_ABORT( MPI_COMM_WORLD, 1, ierr )
  end if

  fvec = rz_end - x

  return
END SUBROUTINE axis_fcn

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

SUBROUTINE ppiota(rzrzt,iflag)
  USE globals, only : dp, myid, IsQuiet, ounit, zero, pi2, sqrtmachprec, pp_phi, Nfp_raw, pp_xtol
  USE mpi
  IMPLICIT NONE

  INTEGER, parameter  :: n = 5
  INTEGER  :: iflag
  REAL     :: rzrzt(n)

  INTEGER  :: iwork(5), ierr, ifail
  REAL     :: phi_init, phi_stop, relerr, abserr, work(100+21*N)
  EXTERNAL :: BRpZ_iota
  
  ifail = 1
  relerr = pp_xtol
  abserr = sqrtmachprec
  phi_init = pp_phi
  phi_stop = pp_phi + pi2/Nfp_raw

  call ode( BRpZ_iota, n, rzrzt, phi_init, phi_stop, relerr, abserr, ifail, work, iwork )
  if ( ifail /= 2 ) then     
     if ( IsQuiet < -1 ) then
        write ( ounit, '(A,I3)' ) 'ppiota  : ODE solver ERROR; returned IFAIL = ', ifail
        select case ( ifail )
        case ( 3 )
           write(ounit, '("ppiota  : DF_xtol or abserr too small.")')
        case ( 4 )
           write(ounit, '("ppiota  : tau not reached after 500 steps.")')
        case ( 5 )
           write(ounit, '("ppiota  : tau not reached because equation to be stiff.")')
        case ( 6 )
           write(ounit, '("ppiota  : INVALID input parameters.")')
        end select
     end if
     iflag = -1
     ! call MPI_ABORT( MPI_COMM_WORLD, 1, ierr )
  end if

  return
END SUBROUTINE ppiota

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

SUBROUTINE BRpZ( t, x, dx )
  !----------------------
  ! dR/dphi = BR / Bphi
  ! dZ/dphi = BZ / Bphi
  !----------------------
  use globals, only : dp, zero, ounit, myid, ierr, myworkid, nworker, MPI_COMM_MYWORLD
  USE MPI
  implicit none

  !---------------------------------------------------------------------------------------------   
  INTEGER, parameter   :: n=2
  REAL, INTENT( IN)    :: t, x(n)
  REAL, INTENT(OUT)    :: dx(n)

  REAL                 :: RR, ZZ, XX, YY, BR, BP, BZ, B(3)
  external             :: coils_bfield
  !---------------------------------------------------------------------------------------------
  
  RR = x(1); ZZ = x(2)           ! cylindrical coordinate
  XX = RR*cos(t); YY = RR*sin(t) ! cartesian   coordinate
  B = zero

  call coils_bfield(B, XX, YY, ZZ)
  
  BR =     B(1)*cos(t) + B(2)*sin(t)
  BP = ( - B(1)*sin(t) + B(2)*cos(t) ) / RR
  BZ =     B(3)

  dx(1) = BR/BP
  dx(2) = BZ/BP

  return
END SUBROUTINE BRpZ
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

SUBROUTINE BRpZ_iota( t, x, dx )
  !----------------------
  ! dR/dphi = BR / Bphi
  ! dZ/dphi = BZ / Bphi
  !----------------------
  use globals, only : dp, zero, ounit, myid, ierr
  USE MPI
  implicit none

  !---------------------------------------------------------------------------------------------   
  INTEGER, parameter   :: n=5
  REAL, INTENT( IN)    :: t, x(n)
  REAL, INTENT(OUT)    :: dx(n)

  REAL                 :: RR, ZZ, XX, YY, BR, BP, BZ, B(3), length
  external             :: coils_bfield
  !---------------------------------------------------------------------------------------------

  ! field line
  RR = x(1); ZZ = x(2)           ! cylindrical coordinate
  XX = RR*cos(t); YY = RR*sin(t) ! cartesian   coordinate
  B = zero

  call coils_bfield(B, XX, YY, ZZ)
  
  BR =     B(1)*cos(t) + B(2)*sin(t)
  BP = ( - B(1)*sin(t) + B(2)*cos(t) ) / RR
  BZ =     B(3)

  dx(1) = BR/BP
  dx(2) = BZ/BP

  ! magnetic axis
  RR = x(3); ZZ = x(4)           ! cylindrical coordinate
  XX = RR*cos(t); YY = RR*sin(t) ! cartesian   coordinate
  B = zero

  call coils_bfield(B, XX, YY, ZZ)
  
  BR =     B(1)*cos(t) + B(2)*sin(t)
  BP = ( - B(1)*sin(t) + B(2)*cos(t) ) / RR
  BZ =     B(3)

  dx(3) = BR/BP
  dx(4) = BZ/BP

  ! integrate theta
  length = (x(1) - x(3))**2 + (x(2)-x(4))**2  ! delta R^2 + delta Z^2
  dx(5) = ( (x(1) - x(3))*(dx(2)-dx(4)) - (x(2)-x(4))*(dx(1)-dx(3)) ) / length

  return
END SUBROUTINE BRpZ_iota
    