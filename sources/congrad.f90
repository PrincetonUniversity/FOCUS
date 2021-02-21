!title (congrad) ! Nonlinear conjugate gradient method

!latex \briefly{The nonlinear conjugate gradient method is an iterative method using the conjugacy property 
!latex          of the matrix to solve a nonlinear system of equations. Details of the method can be found 
!latex          in \emph{Numerical Optimization}(Nocedal, J., \& Wright, S. (2006). Numerical optimization. 
!latex          Springer Science \& Business Media.).}

!latex \calledby{\link{solvers}}
!latex \calls{\link{packdof}}

!latex \section{Basic algorithm}
!latex The basic nonlinear conjugate gradient algorithm used in FAMUS is a hybrid of {\bf Algorithm 5.4} 
!latex and {\bf Equation (5.49)} from \emph{Numerical Optimization} and 
!latex \doilink{10.1137/S1052623497318992}{Dai \& Yuan}. This version of conjugate gradient method converges 
!latex globally, provided the line search satisfies the standard Wolfe conditions. 
!latex 
!latex Our target function is $\chi^2(\vect{X})$, while $\vect{X}$ is the variables vector. As we mentioned 
!latex before, we can calculate the gradient $G(\vect{X}) = \pdv{\chi^2}{\vect{X}}$ accurately with analytical
!latex expressions. The structure of the algorithm is as below.
!latex \begin{tcolorbox}
!latex  $k=0$: for initial $\vect{X}_0$, evaluate $\chi^2(\vect{X}_0)$ and $G(\vect{X}_0)$; $p_0 = -G_0$; \\
!latex {\bf while} $G_k>\epsilon$ : \\ 
!latex               $\vect{X}_{k+1} = \vect{X}_k + \a_k p_k$ ($\a_k$ satisfies strong Wolfe condition); \\
!latex               $\b_{k+1} = \frac{|G_{k+1}|^2}{(G_{k+1} - G_k)^T p_k}$; \\
!latex               $p_{k+1} = -G_{k+1} + \b_{k+1} p_k$ ; \\
!latex               $k = k + 1$ ; \\
!latex {\bf end(while)}
!latex \end{tcolorbox}
!latex 
!latex The line search algorithm is applying the {\bf Algorithm 3.5 \& 3.6} in the book to satisfy 
!latex the strong Wolfe conditions.

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
#ifdef oldcg
SUBROUTINE congrad
  use globals, only: dp, sqrtmachprec, myid, ounit, Ncoils, Ndof, t1E, iout, CG_maxiter, CG_xtol, xdof, &
       exit_signal, tstart, tfinish
  use mpi
  implicit none

  INTEGER                 :: idof, icoil, c1, n1, ierr, astat, iter, iflag
  REAL                    :: alpha, beta, f
  REAL, dimension(1:Ndof) :: lxdof, p, gradk, gradf

  tfinish = MPI_Wtime()
  iter = 0
  call packdof(lxdof(1:Ndof)) ! initial xdof;
  call getdf(lxdof, f, gradk)
  p(1:Ndof) = -gradk ! initial step direction;
  alpha = 1.0 ! initial step size;

  if (myid == 0) write(ounit, '("output  : "A6" : "8(A12," ; "))') "iout", "time (s)", "chi", "dE_norm", &
       "Bnormal", "PM eff. vol.", "Dip. Binary", "Total PM Vol.", "Res. Bn_mn"

  do

     iter = iter + 1
     call wolfe(lxdof, p, alpha, iflag) ! find a step size matching the Wolfe condiction;
     if ( iflag == -1 ) then
        if(myid .eq. 0) write(ounit, '("congrad : EXITING--------Unconverged line search!")')
        exit  ! reach maximum iterations;
     endif
        
     lxdof = lxdof + alpha*p ! next xdof

     call getdf(lxdof, f, gradf)

     tstart = MPI_Wtime()
     call output(tstart-tfinish)

     if ( sqrt(sum(gradf**2)) < CG_xtol ) then
        if(myid .eq. 0) write(ounit, '("congrad : EXITING--------Local minimum reached!")')
        exit  ! reach minimum
     endif

     beta = sum(gradf**2) / sum( (gradf-gradk)*p )
     p = -gradf + beta*p ! direction for next step;
     gradk = gradf  !save for the current step;

     alpha = 1.0  ! reset alpha;

     if (iter .ge. CG_maxiter) then
        if(myid .eq. 0) write(ounit, '("congrad : EXITING--------Maximum iterations reached!")')
        exit  ! reach maximum iterations;
     endif

     if ( exit_signal ) then
        if(myid .eq. 0) write(ounit, '("congrad : EXITING-------No obvious change in last 5 outputs!")')
        exit         ! no obvious changes in past 5 iterations; 07/20/2017
     endif

  enddo

  if(myid .eq. 0) write(ounit, '("congrad : Computation using conjugate gradient finished.")')

  return
END SUBROUTINE congrad
#endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

SUBROUTINE wolfe( x0, p, alpha, iflag )

  use globals, only: dp, zero, sqrtmachprec, ounit, myid, Ndof, CG_wolfe_c1, CG_wolfe_c2

  implicit none
  include "mpif.h"

  REAL   , INTENT( in)     :: x0(1:Ndof), p(1:Ndof)
  INTEGER, INTENT(out)     :: iflag
  REAL   , INTENT(inout)   :: alpha

  REAL                    :: zoom
  INTEGER                 :: i, maxiter
  REAL                    :: a0, ap, ac, f0, fc, fp, rd, c1, c2
  REAL, dimension(1:Ndof) :: xc, g0, gc, gp

  iflag = 0 ! normal value;

  c1 = CG_wolfe_c1
  c2 = CG_wolfe_c2     ! c1 & c2
  i = 0 ; maxiter = 10
  a0 = 0.0

  call getdf(x0, f0, g0)
  
  ap = a0            ! previous alpha;
  fp = f0            ! previous fnction;
  gp = g0            ! previous gradient;

  ac = alpha

  do i = 1, maxiter

     xc = x0 + ac*p  ! current xdof

     call getdf(xc, fc, gc)

     if ( (fc > f0 + c1*ac*sum(p*g0)) .or. (fc >= fp .and. i > 1) ) then
#ifdef DEBUG
        if (myid == 0) write(ounit,'("congrad : option 1 ; ap = "F12.5" ; ac = "F12.5)') ap, ac
#endif
        alpha = zoom( x0, p, ap, ac )
        if (abs(alpha) < sqrtmachprec) exit    ! if zoom failed;
#ifdef DEBUG
        if (myid == 0) write(ounit,'("congrad : zoom failed to find next alpha.")')
#endif
        return
     endif

     if ( abs(sum(p*gc)) <= -c2*sum(p*g0) ) then
        alpha = ac
        return
     endif

     if ( sum(p*gc) >= zero ) then
#ifdef DEBUG
        if (myid == 0) write(ounit,'("congrad : option 3 ; ap = "F12.5" ; ac = "F12.5)') ap, ac
#endif
        alpha = zoom( x0, p, ac, ap )
        if (abs(alpha) < sqrtmachprec) exit    ! if zoom failed;
#ifdef DEBUG
        if (myid == 0) write(ounit,'("congrad : zoom failed to find next alpha.")')
#endif
        return
     endif          

     ap = ac
     fp = fc
     gp = gc

     ! if too many iterations then increase alpha_max;
     ac = ac*2.0
     
  end do

#ifdef DEBUG
  if (myid == 0) write(ounit, '("congrad : line search does not converge. Please change c1, c2, or weights.")')
#endif

  iflag = -1 ! not converge; terminate the optimization;

  return

END SUBROUTINE wolfe

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

REAL FUNCTION zoom( x0, p, alo, ahi )

  use globals, only : dp, zero, ounit, myid, Ndof, CG_wolfe_c1, CG_wolfe_c2

  implicit none
  include "mpif.h"

  REAL, INTENT(   in)     :: x0(1:Ndof), p(1:Ndof)
  REAL, INTENT(inout)     :: alo, ahi

  INTEGER                 :: iter, maxiter = 20 ! maximum iteration;
  REAL                    :: f0, fc, fl, fp, alpha, c1, c2
  REAL, dimension(1:Ndof) :: xc, xl, g0, gc, gp, gl

  iter = 0
  c1 = CG_wolfe_c1
  c2 = CG_wolfe_c2     ! c1 & c2

  call getdf(x0, f0, g0)

  do
     alpha = 0.5*(alo + ahi) ! linearly intepolation;

     xc = x0 + alpha*p
     call getdf(xc, fc, gc)

     xl = x0 + alo*p
     call getdf(xl, fl, gl)

     if ( (fc > f0 + c1*alpha*sum(p*g0)) .or. (fc >= fl) ) then 
        ahi = alpha
     else
        if ( abs(sum(p*gc)) <= -c2*sum(p*g0) ) then
           zoom = alpha
           return
        endif

        if ( (ahi-alo)*sum(p*gc) >= zero ) then
           ahi = alo
        endif

        alo = alpha

     endif

     iter = iter + 1
     if (iter >= maxiter) then 
#ifdef DEBUG
        if (myid == 0) write(ounit,'("congrad : zoom failed to find next alpha, reset to zero.")')
#endif
        zoom = zero
        return
     endif

  end do

  return
END FUNCTION zoom

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!


SUBROUTINE getdf(lxdof, f, g)
  use globals, only: dp, myid, ounit, ierr, Ndof, chi, t1E, MPI_COMM_FAMUS
  implicit none
  include "mpif.h"

  REAL, INTENT(in ) :: lxdof(1:Ndof)
  REAL, INTENT(out) :: f, g(1:Ndof)

  
  call MPI_BARRIER( MPI_COMM_FAMUS, ierr ) ! wait all cpus;

  call unpacking(lxdof)
  call costfun(1)
  f = chi
  g = t1E

  return
END SUBROUTINE getdf

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
#ifndef oldcg

SUBROUTINE congrad
  use globals, only: dp, sqrtmachprec, myid, ounit, Ncoils, Ndof, t1E, iout, CG_maxiter, CG_xtol, xdof, &
       exit_signal, tstart, tfinish
  use mpi
  implicit none

  INTEGER                 :: ierr, astat, iter, n, nfunc, ngrad, status
  REAL                    :: f, gnorm
  REAL, dimension(1:Ndof) :: x, g, d, xtemp, gtemp
  EXTERNAL                :: myvalue, mygrad

  tfinish = MPI_Wtime()
  iter = 0
  n = Ndof
  call packdof(x(1:n)) ! initial xdof;
  call cg_descent (CG_xtol, x, n, myvalue, mygrad, status, gnorm, f, iter, nfunc, ngrad, d, g, xtemp, gtemp)

  if (myid == 0) then
     select case (status)
     case (0) 
        write(ounit, '("congrad : status="I1": convergence tolerance satisfied.")')  status
     case (1) 
        write(ounit, '("congrad : status="I1": change in func <= feps*|f|.")')  status
     case (2) 
        write(ounit, '("congrad : status="I1": total iterations exceeded maxit.")')  status
     case (3) 
        write(ounit, '("congrad : status="I1": slope always negative in line search.")')  status
     case (4) 
        write(ounit, '("congrad : status="I1": number secant iterations exceed nsecant.")')  status
     case (5) 
        write(ounit, '("congrad : status="I1": search direction not a descent direction.")')  status
     case (6) 
        write(ounit, '("congrad : status="I1": line search fails in initial interval.")')  status
     case (7) 
        write(ounit, '("congrad : status="I1": line search fails during bisection.")')  status
     case (8) 
        write(ounit, '("congrad : status="I1": line search fails during interval update.")')  status
     case default 
        write(ounit, '("congrad : status="I1": unknow options!")')  status
     end select
  end if

  if(myid .eq. 0) write(ounit, '("congrad : Computation using conjugate gradient finished.")')
     
  return
END SUBROUTINE congrad
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
SUBROUTINE myvalue(f, x, n)
  use globals, only: dp, myid, ounit, ierr, chi, MPI_COMM_FAMUS
  implicit none
  include "mpif.h"

  INTEGER, INTENT(in) :: n
  REAL, INTENT(in)    :: x(n) 
  REAL, INTENT(out)   :: f

  call MPI_BARRIER( MPI_COMM_FAMUS, ierr ) ! wait all cpus;
  call unpacking(x)
  call costfun(0)
  f = chi
  return

END SUBROUTINE myvalue
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
SUBROUTINE mygrad(g, x, n)
  use globals, only: dp, myid, ounit, ierr, t1E, MPI_COMM_FAMUS
  implicit none
  include "mpif.h"

  INTEGER, INTENT(in) :: n
  REAL, INTENT(in)    :: x(n)
  REAL, INTENT(out)   :: g(n)

  call MPI_BARRIER( MPI_COMM_FAMUS, ierr ) ! wait all cpus;
  call unpacking(x)
  call costfun(1)
  g = t1E
END SUBROUTINE mygrad
#endif
