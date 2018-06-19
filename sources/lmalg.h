!title (lmalg) ! levenberg-marquardt algorithm

!latex \briefly{This subroutine will call the function lmder1 from MINPACK.
!latex \em{lmder} is a modified levenberg-marquardt algorithm with analytically calculated jacobian.}

!latex \calledby{\link{solvers}}
!latex \calls{\link{packdof}}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine lmalg
use globals, only: sqrtmachprec, zero, myid, ounit, Ncoils, Ndof, t1E, iout, xdof, &
     tstart, tfinish, NBmn, Nzeta, Nteta, tstart, tfinish, &
     LM_maxiter, LM_xtol, LM_ftol, LM_iter, LM_factor, LM_mfvec
  use mpi
  implicit none
  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  INTEGER                 :: m, n, ldfjac, info, lwa, idof, ierr, astat, iter, ipvt(1:Ndof)
  INTEGER                 :: maxfev, mode, nprint, nfev, njev
  REAL                    :: ftol, xtol, gtol, factor, x(1:Ndof)
  REAL   , allocatable    :: wa(:), fvec(:), fjac(:,:)
  EXTERNAL                :: focus_fcn
  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  tstart = MPI_Wtime()

  m = LM_mfvec
  n = Ndof
  FATAL( lmalg, m <= 0, wrong number of targets )
  FATAL( lmalg, n <= 0, wrong number of parameters )

  ldfjac = m
  lwa = 5*Ndof+m
  info = 0
  ipvt = 0
  nfev = 0
  njev = 0
  maxfev = max(100*(n + 1), LM_maxiter)
  factor = LM_factor
  ftol = LM_ftol
  xtol = LM_xtol
  gtol = zero
  mode = 1
  nprint = 1

  LM_iter = 0 

  SALLOCATE(fvec, (1:m), zero)
  SALLOCATE(fjac, (1:ldfjac, 1:Ndof), zero)
  SALLOCATE(wa, (1:lwa), zero)
  
  call packdof(x(1:Ndof)) ! initial xdof;
  if (myid == 0) write(ounit, '("output  : "A6" : "8(A12," ; "))') "iout", "time (s)", "chi", "dE_norm", &
       "Bnormal", "Bmn harmonics", "tor. flux", "coil length", "c-s sep." 
  !call lmder1(focus_fcn,m,Ndof,x,fvec,fjac,ldfjac,tol,info,ipvt,wa,lwa)
  call lmder(focus_fcn,m,n,x,fvec,fjac,ldfjac,ftol,xtol,gtol,maxfev, &
       wa(1),mode,factor,nprint,info,nfev,njev,ipvt,wa(n+1), &
       wa(2*n+1),wa(3*n+1),wa(4*n+1),wa(5*n+1))

  if (myid == 0) then
     select case (info)
     case (0)
        write(ounit,1000) info, "improper input parameters."
     case (1)
        write(ounit,1000) info, "error in the sum of squares is at most tol."
     case (2)
        write(ounit,1000) info, "error between x and the solution is at most tol."
     case (3)
        write(ounit,1000) info, "conditions for info = 1 and info = 2 both hold."
     case (4)
        write(ounit,1000) info, "fvec is orthogonal to the columns of the jacobian to machine precision."
     case (5)
        write(ounit,1000) info, "number of calls to fcn with iflag = 1 has reached 100*(n+1)."
     case (6)
        write(ounit,1000) info, "ftol is too small, no further reduction in the sum of squares."
     case (7)
        write(ounit,1000) info, "xtol is too small, no further improvement in the approximate solution x."
     case (8)
        write(ounit,1000) info, "gtol is too small, fvec is orthogonal to the columns of the jacobian."
     case default
        write(ounit,1000) info, "unsupported info, due to manually exit or others."
     end select
  endif
         
1000 format(8X, ": ", "info = ", I, ", ", A)

  DALLOCATE( fvec )
  DALLOCATE( fjac )
  DALLOCATE( wa   )
  
 
  return
end subroutine lmalg

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine focus_fcn(m,n,x,fvec,fjac,ldfjac,iflag)
  use globals, only: zero, myid, ounit, tstart, tfinish, LM_iter, LM_maxiter, exit_signal, LM_fvec, LM_fjac
  use mpi
  implicit none
  
  INTEGER, INTENT(in)             :: m,n,ldfjac
  INTEGER, INTENT(inout)          :: iflag
  DOUBLE PRECISION, INTENT(in)    :: x(n)
  DOUBLE PRECISION, INTENT(out)   :: fvec(m),fjac(ldfjac,n)
  INTEGER                         :: idof, ierr, astat


  call unpacking(x(1:n))
  select case (iflag)
  case ( 0 )
     LM_iter = LM_iter + 1

     if (LM_iter > LM_maxiter) then
        iflag = -1
        if(myid .eq. 0) write(ounit, '("lmalg   : EXITING--------maximum iterations reached.")')
        return
     endif

     if ( exit_signal ) then
        if(myid .eq. 0) write(ounit, '("lmalg   : EXITING-------No obvious change in last 5 outputs!")')
        iflag = -1
        return
     endif

     call costfun(0)
     tfinish = MPI_Wtime()
     call output(tfinish-tstart)     
  case ( 1 )
     call costfun(0)
     fvec(1:m) = LM_fvec(1:m)
  case ( 2 )
     call costfun(1)
     fjac(1:ldfjac, 1:n) = LM_fjac(1:ldfjac, 1:n)
  case default
     FATAL( lmalg, .true. , unsupported iflag value )
  end select
  
  return
end subroutine focus_fcn
     
