!title (lmalg) ! levenberg-marquardt algorithm

!latex \briefly{ This subroutine has the interface for calling the function {\em lmder} from MINPACK.
!latex {\em{lmder}} is a modified levenberg-marquardt algorithm with analytically calculated jacobian.}

!latex \calledby{\link{solvers}}
!latex \calls{\link{packdof}}

!latex  \section{General}
!latex  Levenberg-Marquardt algorithm is one of the most famous minimization algorithms.
!latex  There is a brief introduction on \href{https://en.wikipedia.org/wiki/Levenberg%E2%80%93Marquardt_algorithm}{Wikipedia}.
!latex  In FOCUS, we use the subroutine {\em{lmder}} from \href{http://www.netlib.org/minpack/}{MINPACK}.

!latex  \section{Documentation}
!latex  The cost functions (targets) are stored in {\bf LM\_fvec(1:LM\_mfvec)}, 
!latex  and the jacobian is stored in {\bf LM\_fjac(1:LM\_mfvec, 1:Ndof)}. 
!latex  The number of targets ({\bf LM\_mfvec}) must be not be smaller than the number of parameters ({\bf Ndof}).
!latex  
!latex  The targets are consistant of different terms. The following table lists the details.
!latex \begin{table}[h!]
!latex \centering
!latex \caption{{\bf LM\_fvec} components}
!latex \label{my-label}
!latex \begin{tabular}{cccc}
!latex \hline
!latex cost functions &  physiccal meaning &  switch & length  \\
!latex \hline
!latex \link{bnormal} & $\vect{B} \cdot \vect{n}$ on each surface element & weight\_bnorm $> 0$ & \inputvar{Nteta}*\inputvar{Nzeta}  \\
!latex \link{bmnharm} & $wBmn_i (Bmn_i - Bmn^o_i)$ for each reasonant harmonics & weight\_bharm $> 0$ & 2*NBmn in {\em target.harmonics} \\
!latex \link{torflux} & $\Psi_i - \Psi_o$ at each toroidal cross-sections & weight\_tflux $> 0$ & \inputvar{Nzeta}  \\
!latex \link{length} & length penalty of each coil, $L_i - L_o$ or $\exp(L_i)/\exp{L_o}$ & weight\_ttlen $> 0$ & Ncoils - Nfixgeo \\
!latex \link{surfsep} & potential energy between each coil and a control surface & weight\_cssep $> 0$ & Ncoils - Nfixgeo 
!latex \hline
!latex \end{tabular}
!latex \end{table}

!latex Here are the comments from {\em lmder}. The four input parameters for L-M optimizer are \inputvar{LM\_maxiter}, \inputvar{LM\_xtol},  \inputvar{LM\_ftol} and  \inputvar{LM\_factor}. Please look at \link{initial} for more details.
!latex \lstset{language=Fortran}
!latex  \begin{lstlisting}[frame=single]
!latex  c     the subroutine statement is
!latex  c
!latex  c       subroutine lmder(fcn,m,n,x,fvec,fjac,ldfjac,ftol,xtol,gtol,
!latex  c                        maxfev,diag,mode,factor,nprint,info,nfev,
!latex  c                        njev,ipvt,qtf,wa1,wa2,wa3,wa4)
!latex  c
!latex  c     where
!latex  c
!latex  c       fcn is the name of the user-supplied subroutine which
!latex  c         calculates the functions and the jacobian. fcn must
!latex  c         be declared in an external statement in the user
!latex  c         calling program, and should be written as follows.
!latex  c
!latex  c         subroutine fcn(m,n,x,fvec,fjac,ldfjac,iflag)
!latex  c         integer m,n,ldfjac,iflag
!latex  c         double precision x(n),fvec(m),fjac(ldfjac,n)
!latex  c         ----------
!latex  c         if iflag = 1 calculate the functions at x and
!latex  c         return this vector in fvec. do not alter fjac.
!latex  c         if iflag = 2 calculate the jacobian at x and
!latex  c         return this matrix in fjac. do not alter fvec.
!latex  c         ----------
!latex  c         return
!latex  c         end
!latex  c
!latex  c         the value of iflag should not be changed by fcn unless
!latex  c         the user wants to terminate execution of lmder.
!latex  c         in this case set iflag to a negative integer.
!latex  c
!latex  c       m is a positive integer input variable set to the number
!latex  c         of functions.
!latex  c
!latex  c       n is a positive integer input variable set to the number
!latex  c         of variables. n must not exceed m.
!latex  c
!latex  c       x is an array of length n. on input x must contain
!latex  c         an initial estimate of the solution vector. on output x
!latex  c         contains the final estimate of the solution vector.
!latex  c
!latex  c       fvec is an output array of length m which contains
!latex  c         the functions evaluated at the output x.
!latex  c
!latex  c       fjac is an output m by n array. the upper n by n submatrix
!latex  c         of fjac contains an upper triangular matrix r with
!latex  c         diagonal elements of nonincreasing magnitude such that
!latex  c
!latex  c                t     t           t
!latex  c               p *(jac *jac)*p = r *r,
!latex  c
!latex  c         where p is a permutation matrix and jac is the final
!latex  c         calculated jacobian. column j of p is column ipvt(j)
!latex  c         (see below) of the identity matrix. the lower trapezoidal
!latex  c         part of fjac contains information generated during
!latex  c         the computation of r.
!latex  c
!latex  c       ldfjac is a positive integer input variable not less than m
!latex  c         which specifies the leading dimension of the array fjac.
!latex  c
!latex  c       ftol is a nonnegative input variable. termination
!latex  c         occurs when both the actual and predicted relative
!latex  c         reductions in the sum of squares are at most ftol.
!latex  c         therefore, ftol measures the relative error desired
!latex  c         in the sum of squares.
!latex  c
!latex  c       xtol is a nonnegative input variable. termination
!latex  c         occurs when the relative error between two consecutive
!latex  c         iterates is at most xtol. therefore, xtol measures the
!latex  c         relative error desired in the approximate solution.
!latex  c
!latex  c       gtol is a nonnegative input variable. termination
!latex  c         occurs when the cosine of the angle between fvec and
!latex  c         any column of the jacobian is at most gtol in absolute
!latex  c         value. therefore, gtol measures the orthogonality
!latex  c         desired between the function vector and the columns
!latex  c         of the jacobian.
!latex  c
!latex  c       maxfev is a positive integer input variable. termination
!latex  c         occurs when the number of calls to fcn with iflag = 1
!latex  c         has reached maxfev.
!latex  c
!latex  c       diag is an array of length n. if mode = 1 (see
!latex  c         below), diag is internally set. if mode = 2, diag
!latex  c         must contain positive entries that serve as
!latex  c         multiplicative scale factors for the variables.
!latex  c
!latex  c       mode is an integer input variable. if mode = 1, the
!latex  c         variables will be scaled internally. if mode = 2,
!latex  c         the scaling is specified by the input diag. other
!latex  c         values of mode are equivalent to mode = 1.
!latex  c
!latex  c       factor is a positive input variable used in determining the
!latex  c         initial step bound. this bound is set to the product of
!latex  c         factor and the euclidean norm of diag*x if nonzero, or else
!latex  c         to factor itself. in most cases factor should lie in the
!latex  c         interval (.1,100.).100. is a generally recommended value.
!latex  c
!latex  c       nprint is an integer input variable that enables controlled
!latex  c         printing of iterates if it is positive. in this case,
!latex  c         fcn is called with iflag = 0 at the beginning of the first
!latex  c         iteration and every nprint iterations thereafter and
!latex  c         immediately prior to return, with x, fvec, and fjac
!latex  c         available for printing. fvec and fjac should not be
!latex  c         altered. if nprint is not positive, no special calls
!latex  c         of fcn with iflag = 0 are made.
!latex  c
!latex  c       info is an integer output variable. if the user has
!latex  c         terminated execution, info is set to the (negative)
!latex  c         value of iflag. see description of fcn. otherwise,
!latex  c         info is set as follows.
!latex  c
!latex  c         info = 0  improper input parameters.
!latex  c
!latex  c         info = 1  both actual and predicted relative reductions
!latex  c                   in the sum of squares are at most ftol.
!latex  c
!latex  c         info = 2  relative error between two consecutive iterates
!latex  c                   is at most xtol.
!latex  c
!latex  c         info = 3  conditions for info = 1 and info = 2 both hold.
!latex  c
!latex  c         info = 4  the cosine of the angle between fvec and any
!latex  c                   column of the jacobian is at most gtol in
!latex  c                   absolute value.
!latex  c
!latex  c         info = 5  number of calls to fcn with iflag = 1 has
!latex  c                   reached maxfev.
!latex  c
!latex  c         info = 6  ftol is too small. no further reduction in
!latex  c                   the sum of squares is possible.
!latex  c
!latex  c         info = 7  xtol is too small. no further improvement in
!latex  c                   the approximate solution x is possible.
!latex  c
!latex  c         info = 8  gtol is too small. fvec is orthogonal to the
!latex  c                   columns of the jacobian to machine precision.
!latex  c
!latex  c       nfev is an integer output variable set to the number of
!latex  c         calls to fcn with iflag = 1.
!latex  c
!latex  c       njev is an integer output variable set to the number of
!latex  c         calls to fcn with iflag = 2.
!latex  c
!latex  c       ipvt is an integer output array of length n. ipvt
!latex  c         defines a permutation matrix p such that jac*p = q*r,
!latex  c         where jac is the final calculated jacobian, q is
!latex  c         orthogonal (not stored), and r is upper triangular
!latex  c         with diagonal elements of nonincreasing magnitude.
!latex  c         column j of p is column ipvt(j) of the identity matrix.
!latex  c
!latex  c       qtf is an output array of length n which contains
!latex  c         the first n elements of the vector (q transpose)*fvec.
!latex  c
!latex  c       wa1, wa2, and wa3 are work arrays of length n.
!latex  c
!latex  c       wa4 is a work array of length m.
!latex  c
!latex  c     subprograms called
!latex  c
!latex  c       user-supplied ...... fcn
!latex  c
!latex  c       minpack-supplied ... dpmpar,enorm,lmpar,qrfac
!latex  c
!latex  c       fortran-supplied ... dabs,dmax1,dmin1,dsqrt,mod
!latex  c
!latex  c     argonne national laboratory. minpack project. march 1980.
!latex  c     burton s. garbow, kenneth e. hillstrom, jorge j. more
!latex  \end{lstlisting}

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
     
