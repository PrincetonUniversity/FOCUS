!title (congrad) ! Nonlinear conjugate gradient method

!latex \briefly{The nonlinear conjugate gradient method is an iterative method using the conjugacy property 
!latex          of the matrix to solve a nonlinear system of equations. Details of the method can be found 
!latex          in \emph{Numerical Optimization}(Nocedal, J., \& Wright, S. (2006). Numerical optimization. 
!latex          Springer Science \& Business Media.).}

!latex \calledby{\link{solvers}}
!latex \calls{packdof}

!latex \section{Basic algorithm}
!latex The basic nonlinear conjugate gradient algorithm used in FOCUS is a hybrid of {\bf Algorithm 5.4} 
!latex and {\bf Equation (5.49)} from \emph{Numerical Optimization} and 
!latex \doilink{10.1137/S1052623497318992}{Dai & Yuan}. This version of conjugate gradient method converges 
!latex globally, provided the line search satisfies the standard Wolfe conditions. 
!latex 
!latex Our target function is $\chi^2(\vect{X})$, while $\vect{X}$ is the variables vector. As we mentioned 
!latex before, we can calculate the gradient $G(\vect{X}) = \pdv{\chi^2}{\vect{X}}$ accurately with analytical
!latex expressions. The structure of the algorithm is as below.
!latex \fbox{ \beigin{minipage}{\textwidth}
!latex \begin{itemize}
!latex \item[$k=0$:] for initial $\vect{X}_0$, evaluate $\chi^2(\vect{X}_0)$ and $G(\vect{X}_0)$;
!latex               $p_0 =  -G_0$
!latex \item[{\bf while} $G_k>\epsilon$] \\
!latex               $\vect{X}_{k+1} = \vect{X}_k + \a_k p_k$ ($\a_k$ satisfies strong Wolfe condition) \\
!latex               $\b_{k+1} = \frac{\abs{G_{k+1}}^2}{(G_{k+1} - G_k)^T p_k}$ \\
!latex               $p_{k+1} = -G_{k+1} + \b_{k+1} p_k$
!latex               $k = k + 1$
!latex \item[{\bf end(while)}]
!latex \end{itemize}
!latex \end{minipage}}
!latex 
!latex The line search algorithm is applying the {\bf Algorithm 3.5 & 3.6} in the book to satisfy 
!latex the strong Wolfe conditions.

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

SUBROUTINE congrad
  use globals, only: sqrtmachprec, myid, ounit, Ncoils, Ndof, t1E, iout, CG_maxiter, CG_xtol, &
       chi, t1E
  implicit none
  include "mpif.h"

  INTEGER                 :: idof, icoil, c1, n1, ierr, astat, iter
  REAL                    :: alpha, beta, f
  REAL, dimension(1:Ndof) :: lxdof, p, gradk, gradf

  iter = 0
  call packdof(lxdof(1:Ndof)) ! initial xdof;  
  call costfun(1)
  p(1:Ndof) = -t1E ! initial step direction;
  alpha = 1.0 ! initial step size;

  if (myid == 0) write(ounit, '("output  : "A6" : "9(A12," ; "))') "iout", "iter", "chi", "dE_norm", &
       "Bnormal", "Bmn harmonics", "tor. flux", "coil length", "spectral", "c-c sep." 
  call output(real(iter))

  do

     iter = iter + 1
     call wolfe(lxdof, p, alpha) ! find a step size matching the Wolfe condiction;
     lxdof = lxdof + alpha*p(1:Ndof) ! next xdof

     call unpacking(lxdof)
     call costfun(1)
     f = chi
     gradf = t1E

     call output(real(iter))

     if ( sum(gradf**2) .lt. sqrtmachprec) exit  ! reach minimum
 
     beta = sum(gradf**2) / sum( (gradf-gradk)*p )
     p = -gradf + beta*p ! direction for next step;
     gradk = gradf  !save for the current step;

     alpha = 1.0  ! reset alpha;

     if (iter .ge. CG_maxiter) exit  ! reach maximum iterations;

  enddo

  if(myid .eq. 0) write(ounit, '("congrad : Computation using conjugate gradient finished.")')

  return
END SUBROUTINE congrad

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

SUBROUTINE wolfe( x0, p, alpha )

  use globals, only : zero, sqrtmachprec, ounit, myid, Ndof, CG_wolfe_c1, CG_wolfe_c2, chi, t1E

  implicit none
  include "mpif.h"

  REAL, INTENT( in)       :: x0(1:Ndof), p(1:Ndof)
  REAL, INTENT(out)       :: alpha

  REAL                    :: zoom
  INTEGER                 :: i
  REAL                    :: a0, am, ap, ac, f0, fc, fp, rd, c1, c2
  REAL, dimension(1:Ndof) :: xc, g0, gc, gp


  c1 = CG_wolfe_c1
  c2 = CG_wolfe_c2     ! c1 & c2
  i = 0
  a0 = 0.0
  am = alpha
  if (am <= 0.0) STOP "alpha_max should be larger than alpha_0"

  call unpacking(x0)
  call costfun(1)
  f0 = chi
  g0 = t1E
  
  ap = a0            ! previous alpha;
  fp = f0            ! previous fnction;
  gp = g0            ! previous gradient;

  call RANDOM_NUMBER(rd)
  ac = am*rd   ! current alpha;

  do
     xc = x0 + ac*p  ! current xdof

     call unpacking(xc)
     call costfun(1)
     fc = chi
     gc = t1E

     if ( (fc > f0 + c1*ac*sum(p*g0)) .or. (fc >= fp .and. i > 1) ) then
        alpha = zoom( x0, p, ap, ac )
        return
     endif

     if ( abs(sum(p*gc)) <= -c2*sum(p*g0) ) then
        alpha = ac
        return
     endif

     if ( sum(p*gc) >= zero ) then
        alpha = zoom( x0, p, ac, ap )
        return
     endif        

     ap = ac
     fp = fc
     gp = gc

     ! if too many iterations then increase alpha_max;
     if ( i > 1000 .or. (am-ac) .lt. sqrtmachprec ) then 
        am = 2.0 * am
        i = 0
#ifdef DEBUG
        if (myid .eq. 0) print '("congrad : ", A, F10.5)', "change am to", am
#endif
     endif

     call RANDOM_NUMBER(rd)
     ac = ac + (am-ac)*rd

     i = i+1
     
  end do

  return

END SUBROUTINE wolfe

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

REAL FUNCTION zoom( x0, p, alo, ahi )

  use globals, only : zero, ounit, myid, Ndof, CG_wolfe_c1, CG_wolfe_c2, chi, t1E

  implicit none
  include "mpif.h"

  REAL, INTENT(   in)     :: x0(1:Ndof), p(1:Ndof)
  REAL, INTENT(inout)     :: alo, ahi

  REAL                    :: f0, fc, fl, fp, alpha, c1, c2
  REAL, dimension(1:Ndof) :: xc, xl, g0, gc, gp, gl


  c1 = CG_wolfe_c1
  c2 = CG_wolfe_c2     ! c1 & c2

  call unpacking(x0)
  call costfun(1)
  f0 = chi
  g0 = t1E

  do
     alpha = 0.5*(alo + ahi)
     xc = x0 + alpha*p
     call unpacking(xc)
     call costfun(1)
     fc = chi
     gc = t1E

     xl = x0 + alo*p
     call unpacking(xl)
     fl = chi
     gl = t1E

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
  end do

  return
END FUNCTION zoom

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
