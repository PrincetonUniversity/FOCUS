c      ________________________________________________________________
c     |      A conjugate gradient method with guaranteed descent       |
c     |                                                                |
c     |             Version 1.1  (December 10, 2004)                   |
c     |             Version 1.2  (June 4, 2005)                        |
c     |             Version 1.3  (October 6, 2005)                     |
c     |             Version 1.4  (November 14, 2005)                   |
c     |                                                                |
c     |           William W. Hager    and   Hongchao Zhang             |
c     |          hager@math.ufl.edu       hzhang@math.ufl.edu          |
c     |                   Department of Mathematics                    |
c     |                     University of Florida                      |
c     |                 Gainesville, Florida 32611 USA                 |
c     |                      352-392-0281 x 244                        |
c     |                                                                |
c     |              Copyright 2004 by William W. Hager                |
c     |                                                                |
c     |This program is free software; you can redistribute it and/or   |
c     |modify it under the terms of the GNU General Public License as  |
c     |published by the Free Software Foundation; either version 2 of  |
c     |the License, or (at your option) any later version.             |
c     |This program is distributed in the hope that it will be useful, |
c     |but WITHOUT ANY WARRANTY; without even the implied warranty of  |
c     |MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the   |
c     |GNU General Public License for more details.                    |
c     |                                                                |
c     |You should have received a copy of the GNU General Public       |
c     |License along with this program; if not, write to the Free      |
c     |Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, |
c     |MA  02110-1301  USA                                             |
c     |                                                                |
c     |          http://www.math.ufl.edu/~hager/papers/CG              |
c     |                                                                |
c     |    INPUT:                                                      |
c     |                                                                |
c     |(double) grad_tol-- StopRule = T: |g|_infty <= max (grad_tol,   |
c     |                          StopFac*initial |g|_infty) [default]  |
c     |                    StopRule = F: |g|_infty <= grad_tol(1+|f|)  |
c     |                                                                |
c     |(double) x       --starting guess (length n)                    |
c     |                                                                |
c     |(int)    dim     --problem dimension (also denoted n)           |
c     |                                                                |
c     |         cg_value--name of cost evaluation subroutine           |
c     |                  (external in main program, cg_value(f, x, n)  |
c     |                   puts value of cost function at x in f        |
c     |                   f is double precision scalar and x is        |
c     |                   double precision array of length n)          |
c     |                                                                |
c     |         cg_grad --name gradient evaluation subroutine          |
c     |                  (external in main program, cg_grad (g, x, n)  |
c     |                   puts gradient at x in g, g and x are         |
c     |                   double precision arrays of length n)         |
c     |                                                                |
c     |(double) gnorm   --if the parameter Step in cg_descent.parm is  |
c     |                   .true., then gnorm contains the initial step |
c     |                   used at iteration 0 in the line search       |
c     |                                                                |
c     |(double) d       --direction (work array, length n)             |
c     |                                                                |
c     |(double) g       --gradient (work array, length n)              |
c     |                                                                |
c     |(double) xtemp   --work array (work array, length n)            |
c     |                                                                |
c     |(double) gtemp   --work array (work array, length n)            |
c     |                                                                |
c     |    OUTPUT:                                                     |
c     |                                                                |
c     |(int)    status  -- 0 (convergence tolerance satisfied)         |
c     |                    1 (change in func <= feps*|f|)              |
c     |                    2 (total iterations exceeded maxit)         |
c     |                    3 (slope always negative in line search)    |
c     |                    4 (number secant iterations exceed nsecant) |
c     |                    5 (search direction not a descent direction)|
c     |                    6 (line search fails in initial interval)   |
c     |                    7 (line search fails during bisection)      |
c     |                    8 (line search fails during interval update)|
c     |                                                                |
c     |(double) gnorm   --max abs component of gradient                |
c     |                                                                |
c     |(double) f       --function value at solution                   |
c     |                                                                |
c     |(double) x       --solution (length n)                          |
c     |                                                                |
c     |(int)    iter    --number of iterations                         |
c     |                                                                |
c     |(int)    nfunc   --number of function evaluations               |
c     |                                                                |
c     |(int)    ngrad   --number of gradient evaluations               |
c     |                                                                |
c     |Note: The file cg_descent.parm must be placed in the directory  |
c     |      where the code is run                                     |
c     |________________________________________________________________|
c
      subroutine cg_descent (grad_tol, x, dim, cg_value, cg_grad,
     &                       status, gnorm, f, iter, nfunc, ngrad,
     &                       d, g, xtemp, gtemp)

       use globals, only: dp, myid, ounit, IsQuiet, output_use, tstart, 
     &                    tfinish
       use mpi

      double precision x (*), d (*), g (*), xtemp (*), gtemp (*),
     &                 delta, sigma, eps,
     &                 gamma, rho, tol, eta, fpert, f0, Ck, Qdecay,
     &                 wolfe_hi, wolfe_lo, awolfe_hi,
     &                 QuadCutOff, StopFac, AWolfeFac,
     &                 zero, feps, psi0, psi1, psi2,
     &                 grad_tol, delta2, eta_sq, Qk,
     &                 f, ftemp, gnorm, xnorm, gnorm2, dnorm2, denom,
     &                 t, t1, t2, t3, t4, dphi, dphi0, alpha, talpha,
     &                 yk, yk2, ykgk, dkyk, beta

      integer          n, n5, n6, nf, ng, info, nrestart,
     &                 nexpand, nsecant, maxit,
     &                 iter, status, nfunc, ngrad,
     &                 i, j, i1, i2, i3, i4, dim

      logical          PertRule, QuadOK, QuadStep, PrintLevel,
     &                 PrintFinal, StopRule, AWolfe, Step, debug,
     &                 cg_tol

      external         cg_value, cg_grad

      common /cgparms/delta, sigma, eps,
     &                gamma, rho, tol, eta, fpert, f0, Ck, Qdecay,
     &                wolfe_hi, wolfe_lo, awolfe_hi,
     &                QuadCutOff, StopFac, AWolfeFac,
     &                zero, feps, psi0, psi1, psi2,
     &                n, n5, n6, nf, ng, info,
     &                nrestart, nexpand, nsecant, maxit,
     &                PertRule, QuadOK, QuadStep, PrintLevel,
     &                PrintFinal, StopRule, AWolfe, Step, debug

c initialize the parameters

      call cg_init (grad_tol, dim)

      if ( Step ) then
          alpha = gnorm
      endif
      delta2 = 2*delta - 1
      eta_sq = eta*eta
      iter = 0
      Ck = 0
      Qk = 0

c initial function and gradient evaluations, initial direction

      call cg_value (f, x, n)
      nf = nf + 1
      call cg_grad (g, x, n)
      ng = ng + 1
      f0 = f + f
      gnorm = zero
      xnorm = zero
      gnorm2 = zero
      do i = 1, n5
          xnorm = dmax1 (xnorm, dabs (x (i)))
          t = g (i)
          d (i) = -t
          gnorm = dmax1 (gnorm, dabs(t))
          gnorm2 = gnorm2 + t*t
      enddo
      do i = n6, n, 5
          xnorm = dmax1 (xnorm, dabs (x (i)))
          t = g (i)
          gnorm = dmax1 (gnorm, dabs (t))
          d (i)   = -t
          j = i + 1
          t1 = g (j)
          d (j) = -t1
          gnorm = dmax1 (gnorm, dabs (t1))
          xnorm = dmax1 (xnorm, dabs (x (j)))
          j = i + 2
          t2 = g (j)
          d (j) = -t2
          gnorm = dmax1 (gnorm, dabs (t2))
          xnorm = dmax1 (xnorm, dabs (x (j)))
          j = i + 3
          t3 = g (j)
          d (j) = -t3
          gnorm = dmax1 (gnorm, dabs (t3)) 
          xnorm = dmax1 (xnorm, dabs (x (j)))
          j = i + 4
          t4 = g (j)
          d (j) = -t4
          gnorm = dmax1 (gnorm, dabs (t4)) 
          xnorm = dmax1 (xnorm, dabs (x (j)))
          gnorm2 = gnorm2 + t*t + t1*t1 + t2*t2 + t3*t3 + t4*t4
      enddo

      if ( StopRule ) then
          tol = dmax1 (gnorm*StopFac, tol)
      endif

      if ( PrintLevel ) then        
          write (*, 10) iter, f, gnorm, AWolfe
10        format ('iter: ', i5, ' f= ', e14.6,
     &            ' gnorm= ', e14.6, ' AWolfe= ', l2)
      endif

      if ( cg_tol (f, gnorm) ) goto 100

      dphi0 = -gnorm2
      if ( .not.Step ) then
          alpha = psi0*xnorm/gnorm
          if ( xnorm .eq. zero ) then
              if ( f .ne. zero ) then
                  alpha = psi0*dabs (f)/gnorm2
              else
                  alpha = 1.d0
              endif
          endif
      endif
 
c start the conjugate gradient iteration

c
c   alpha starts as old step, ends as initial step for next iteration
c   f is function value for alpha = 0
c   QuadOK = .true. means that a quadratic step was taken
c

      do iter = 1, maxit
          QuadOK = .false.
          alpha = psi2*alpha
          if ( QuadStep ) then               
              if ( f .ne. zero ) then
                  t = dabs ((f-f0)/f)
              else 
                  t = 1.d0
              endif
              if ( t .gt. QuadCutOff ) then
                  talpha = psi1*alpha  
                  call cg_step (xtemp, x, d, talpha)
                  call cg_value (ftemp, xtemp, n)
                  nf = nf + 1
                  if ( ftemp .lt. f ) then
                     denom = 2.0d0*(((ftemp-f)/talpha)-dphi0)
                     if ( denom .gt. zero ) then
                         QuadOK = .true.
                         alpha = -dphi0*talpha/denom
                     endif
                  endif
              endif
          endif
          f0 = f

          if ( PrintLevel .and. IsQuiet<0 ) then
              if (myid .eq. 0) write (*, 20) QuadOK, alpha, f0, dphi0
20            format ('QuadOK:', l2, ' initial a:',
     &                 e14.6,' f0:', e14.6, ' dphi', e14.6)
          endif

c parameters in Wolfe and approximiate Wolfe conditions, and in update

          Qk = Qdecay*Qk + 1.
          Ck = Ck + (dabs (f) - Ck)/Qk

          if ( PertRule ) then
              fpert = f + eps*Ck
          else
              fpert = f + eps
          endif

          wolfe_hi = delta*dphi0
          wolfe_lo = sigma*dphi0
          awolfe_hi = delta2*dphi0
          if ( AWolfe ) then
              call cg_line  (alpha, f, dphi, dphi0, x, xtemp, d, gtemp,
     &                     cg_value, cg_grad)
          else
              call cg_lineW (alpha, f, dphi, dphi0, x, xtemp, d, gtemp,
     &                     cg_value, cg_grad)
          endif

          if ( info .gt. 0 ) goto 100
c
c Test for convergence to within machine epsilon
c (set feps to zero to remove this test)
c
          if ( -alpha*dphi0 .le. feps*dabs (f) ) then
              info = 1
              goto 100
          endif

c compute beta, yk2, gnorm, gnorm2, dnorm2, update x and g, 

          if ( mod (iter, nrestart) .ne. 0 ) then
              gnorm = zero
              dnorm2 = zero
              yk2 = zero
              ykgk = zero
              do i = 1, n5
                  x (i) = xtemp (i)
                  t = gtemp (i)
                  yk = t - g (i)
                  yk2 = yk2 + yk**2
                  ykgk = ykgk + yk*t
                  g (i) = t
                  gnorm = dmax1 (gnorm, dabs (t))
                  dnorm2 = dnorm2 + d (i)**2
              enddo
              do i = n6, n, 5
                  x (i) = xtemp (i)
                  t = gtemp (i)
                  yk = t - g (i)
                  yk2 = yk2 + yk**2
                  ykgk = ykgk + yk*t
                  i1 = i + 1
                  x (i1) = xtemp (i1)
                  t1 = gtemp (i1)
                  i2 = i + 2
                  x (i2) = xtemp (i2)
                  t2 = gtemp (i2)
                  i3 = i + 3
                  x (i3) = xtemp (i3)
                  t3 = gtemp (i3)
                  i4 = i + 4
                  x (i4) = xtemp (i4)
                  t4 = gtemp (i4)
                  yk2 = yk2 + (t1-g (i1))**2 + (t2-g (i2))**2
     &                      + (t3-g (i3))**2 + (t4-g (i4))**2
                  ykgk = ykgk + (t1-g (i1))*t1 + (t2-g (i2))*t2
     &                        + (t3-g (i3))*t3 + (t4-g (i4))*t4
                  g (i) = t
                  gnorm = dmax1 (gnorm, dabs (t))
                  g (i1) = t1
                  gnorm = dmax1 (gnorm, dabs (t1))
                  g (i2) = t2
                  gnorm = dmax1 (gnorm, dabs (t2))
                  g (i3) = t3
                  gnorm = dmax1 (gnorm, dabs (t3)) 
                  g (i4) = t4
                  gnorm = dmax1 (gnorm, dabs (t4)) 
                  dnorm2 = dnorm2 + d (i)**2 + d (i1)**2 + d (i2)**2
     &                                       + d (i3)**2 + d (i4)**2
              enddo
              if ( cg_tol (f, gnorm) ) goto 100
              dkyk = dphi - dphi0
              beta = (ykgk - 2.d0*dphi*yk2/dkyk)/dkyk

c   faster: initialize dnorm2 = gnorm2 at start, then
c             dnorm2 = gnorm2 + beta**2*dnorm2 - 2.d0*beta*dphi
c             gnorm2 = ||g_{k+1}||^2
c             dnorm2 = ||d_{k+1}||^2
c             dpi = g_{k+1}' d_k

              beta = dmax1 (beta,
     &               -1.d0/dsqrt (dmin1 (eta_sq, gnorm2)*dnorm2))

c     update search direction d = -g + beta*dold

              gnorm2 = zero
              do i = 1, n5
                  t = g (i)
                  d (i) = -t + beta*d (i)
                  gnorm2 = gnorm2 + t*t
              enddo
              do i = n6, n, 5
                  d (i) = -g (i) + beta*d (i)
                  i1 = i + 1
                  d (i1) = -g (i1) + beta*d (i1)
                  i2 = i + 2
                  d (i2) = -g (i2) + beta*d (i2)
                  i3 = i + 3
                  d (i3) = -g (i3) + beta*d (i3)
                  i4 = i + 4
                  d (i4) = -g (i4) + beta*d (i4)
                  gnorm2 = gnorm2 + g (i)**2 + g (i1)**2 + g (i2)**2
     &                                       + g (i3)**2 + g (i4)**2
              enddo
              dphi0 = -gnorm2 + beta*dphi

          else

c     search direction d = -g

              if ( PrintLevel .and. IsQuiet <0 ) then
                  if (myid .eq. 0) write (*, *) "RESTART CG"
              endif
              gnorm = zero
              gnorm2 = zero
              do i = 1, n5
                  x (i) = xtemp (i)
                  t = gtemp (i)
                  g (i) = t
                  d (i) = -t
                  gnorm = dmax1 (gnorm, dabs(t))
                  gnorm2 = gnorm2 + t*t
              enddo
              do i = n6, n, 5
                  x (i) = xtemp (i)
                  t = gtemp (i)
                  g (i) = t
                  d (i) = -t
                  gnorm = dmax1 (gnorm, dabs(t))
                  j = i + 1
                  x (j) = xtemp (j)
                  t1 = gtemp (j)
                  g (j) = t1
                  d (j) = -t1
                  gnorm = dmax1 (gnorm, dabs(t1))
                  j = i + 2
                  x (j) = xtemp (j)
                  t2 = gtemp (j)
                  g (j) = t2
                  d (j) = -t2
                  gnorm = dmax1 (gnorm, dabs(t2))
                  j = i + 3
                  x (j) = xtemp (j)
                  t3 = gtemp (j)
                  g (j) = t3
                  d (j) = -t3
                  gnorm = dmax1 (gnorm, dabs(t3)) 
                  j = i + 4
                  x (j) = xtemp (j)
                  t4 = gtemp (j)
                  g (j) = t4
                  d (j) = -t4
                  gnorm = dmax1 (gnorm, dabs(t4)) 
                  gnorm2 = gnorm2 + t*t + t1*t1 + t2*t2 + t3*t3 + t4*t4
              enddo
              if ( cg_tol (f, gnorm) ) goto 100
              dphi0 = -gnorm2
          endif
          if ( .not.AWolfe ) then
              if ( dabs (f-f0) .lt. AWolfeFac*Ck ) then
                  AWolfe = .true.
              endif
          endif
      
          if ( PrintLevel .or. PrintFinal .and. output_use .eq. 1) then
             tstart = MPI_Wtime()
             call output(tstart-tfinish)  
c              write (*, 10) iter, f, gnorm, AWolfe
          endif

          if ( debug ) then
              if ( f .gt. f0 + 1.e-10*Ck ) then
                  write (*, 270)
                  write (*, 50) f, f0
50                format (' new value:', e30.16, 'old value:', e30.16)
                  stop
              endif
          endif
                  
          if ( dphi0 .gt. zero ) then
             info = 5
             goto 100
          endif
      enddo
      info = 2
100   nfunc = nf
      ngrad = ng
      status = info
      if ( info .gt. 2 ) then
          gnorm = zero
          do i = 1, n
              x (i) = xtemp (i)
              g (i) = gtemp (i)
              gnorm = dmax1 (gnorm, dabs(g (i))) 
          enddo
      endif
      if ( PrintFinal .and. .false. ) then
          write (6, *) 'Termination status:', status
          if ( status .eq. 0 ) then
              write (6, 200)
          else if ( status .eq. 1 ) then
              write (6, 210)
          else if ( status .eq. 2 ) then
              write (6, 220) maxit
              write (6, 300)
              write (6, 400) grad_tol
          else if ( status .eq. 3 ) then
              write (6, 230)
              write (6, 300)
              write (6, 430)
              write (6, 410)
          else if ( status .eq. 4 ) then
              write (6, 240)
              write (6, 300)
              write (6, 400) grad_tol
          else if ( status .eq. 5 ) then
              write (6, 250)
          else if ( status .eq. 6 ) then
              write (6, 260)
              write (6, 300)
              write (6, 400) grad_tol
              write (6, 410)
              write (6, 420)
          else if ( status .eq. 7 ) then
              write (6, 260)
              write (6, 300)
              write (6, 400) grad_tol
          else if ( status .eq. 8 ) then
              write (6, 260)
              write (6, 300)
              write (6, 400) grad_tol
              write (6, 410)
              write (6, 420)
          endif
          write (6, 500) gnorm
          write (6, *) 'function value:', f
          write (6, *) 'cg iterations:', iter
          write (6, *) 'function evaluations:', nfunc
          write (6, *) 'gradient evaluations:', ngrad
      endif
      return
200   format (' Convergence tolerance for gradient satisfied')
210   format (' Terminating since change in function value <= feps*|f|')
220   format (' Total number of iterations exceed max allow:', i10)
230   format (' Slope always negative in line search')
240   format (' Line search fails, too many secant steps')
250   format (' Search direction not a descent direction')
260   format (' Line search fails')
270   format (' Debugger is on, function value does not improve')
300   format (' Possible causes of this error message:')
400   format ('   - your tolerance (grad_tol = ', d12.4,
     &        ') may be too strict')
410   format ('   - your gradient routine has an error')
420   format ('   - parameter epsilon in cg_descent.parm is too small')
430   format ('   - your cost function has an error')
500   format (' absolute largest component of gradient: ', d12.4)
      end

c     PARAMETERS:
c
c     delta - range (0, .5), used in the Wolfe conditions
c     sigma - range [delta, 1), used in the Wolfe conditions
c     eps - range [0, infty), used to compute line search perturbation
c     gamma - range (0,1), determines when to perform bisection step
c     rho   - range (1, infty), growth factor when finding initial interval
c     eta   - range (0, infty), used in lower bound for beta
c     psi0  - range (0, 1), factor used in very initial starting guess
c     psi1  - range (0, 1), factor previous step multiplied by in QuadStep
c     psi2  - range (1, infty), factor previous step is multipled by for startup
c     QuadCutOff - perform QuadStep if relative change in f > QuadCutOff
c     StopFac - used in StopRule
c     AWolfeFac - used to decide when to switch from Wolfe to AWolfe
c     restart_fac - range (0, infty) restart cg when iter = n*restart 
c     maxit_fac - range (0, infty) terminate in maxit = maxit_fac*n iterations
c     feps - stop when -alpha*dphi0 (est. change in value) <= feps*|f|
c            (feps = 0 removes this test, example: feps = eps*1.e-5
c             where eps is machine epsilon)
c     tol   - range (0, infty), convergence tolerance
c     nexpand - range [0, infty), number of grow/shrink allowed in bracket
c     nsecant - range [0, infty), maximum number of secant steps
c     PertRule - gives the rule used for the perturbation in f
c                 F => fpert = eps
c                 T => fpert = eps*Ck, Ck is an average of prior |f|
c                             Ck is an average of prior |f|
c     QuadStep- .true. (use quadratic step) .false. (no quadratic step)
c     PrintLevel- .false. (no printout) .true. (print intermediate results)
c     PrintFinal- .false. (no printout) .true. (print messages, final error)
c     StopRule - .true. (max abs grad <= max (tol, StopFac*initial abs grad))
c                .false. (... <= tol*(1+|f|))
c     AWolfe - .false. (use standard Wolfe initially)
c            - .true. (use approximate + standard Wolfe)
c     Step - .false. (program computing starting step at iteration 0)
c          - .true. (user provides starting step in gnorm argument of cg_descent
c     debug - .false. (no debugging)
c           - .true. (check that function values do not increase)
c     info  - same as status
c
c     DEFAULT PARAMETER VALUES:
c
c         delta : 0.1
c         sigma : 0.9
c         eps : 1.e-6
c         gamma : 0.66
c         rho   : 5.0
c         restart: 1.0
c         eta   : 0.01
c         psi0  : 0.01
c         psi1  : 0.1 
c         psi2  : 2.0 
c         QuadCutOff: 1.d-12
c         StopFac: 0.d0
c         AWolfeFac: 1.d-3
c         tol   : grad_tol
c         nrestart: n (restart_fac = 1)
c         maxit : 500*n (maxit_fac = 500)
c         feps : 0.0
c         Qdecay : 0.7
c         nexpand: 50
c         nsecant: 50
c         PertRule: .true.
c         QuadStep: .true.
c         PrintLevel: .false.
c         PrintFinal: .true.
c         StopRule: .true.
c         AWolfe: .false.
c         Step: .false.
c         debug: .false.
c         info  : 0
c         feps  : 0.0
c

c      (double) grad_tol-- used in stopping rule
c      (int)    dim     --problem dimension (also denoted n)

      subroutine cg_init (grad_tol, dim)
      use globals, only : cg_maxiter, CG_wolfe_c1, CG_wolfe_c2, cg_xtol,
     &     IsQuiet
      double precision delta, sigma, eps,
     &                 gamma, rho, tol, eta, fpert, f0, Ck, Qdecay,
     &                 wolfe_hi, wolfe_lo, awolfe_hi,
     &                 QuadCutOff, StopFac, AWolfeFac,
     &                 zero, feps, psi0, psi1, psi2,
     &                 grad_tol, restart_fac, maxit_fac

      integer          n, n5, n6, nf, ng, info, nrestart,
     &                 nexpand, nsecant, maxit, dim

      logical          PertRule, QuadOK, QuadStep, PrintLevel,
     &                 PrintFinal, StopRule, AWolfe, Step, debug

      common /cgparms/delta, sigma, eps,
     &                gamma, rho, tol, eta, fpert, f0, Ck, Qdecay,
     &                wolfe_hi, wolfe_lo, awolfe_hi,
     &                QuadCutOff, StopFac, AWolfeFac,
     &                zero, feps, psi0, psi1, psi2,
     &                n, n5, n6, nf, ng, info,
     &                nrestart, nexpand, nsecant, maxit,
     &                PertRule, QuadOK, QuadStep, PrintLevel,
     &                PrintFinal, StopRule, AWolfe, Step, debug

      n = dim
      tol = grad_tol
c$$$      open (10, file='cg_descent_f.parm')
c$$$      read (10, *) delta
c$$$      read (10, *) sigma
c$$$      read (10, *) eps
c$$$      read (10, *) gamma
c$$$      read (10, *) rho
c$$$      read (10, *) eta
c$$$      read (10, *) psi0
c$$$      read (10, *) psi1
c$$$      read (10, *) psi2
c$$$      read (10, *) QuadCutOff
c$$$      read (10, *) StopFac
c$$$      read (10, *) AWolfeFac
c$$$      read (10, *) restart_fac
c$$$      read (10, *) maxit_fac
c$$$      read (10, *) feps
c$$$      read (10, *) Qdecay
c$$$      read (10, *) nexpand
c$$$      read (10, *) nsecant
c$$$      read (10, *) PertRule
c$$$      read (10, *) QuadStep
c$$$      read (10, *) PrintLevel
c$$$      read (10, *) PrintFinal
c$$$      read (10, *) StopRule
c$$$      read (10, *) AWolfe
c$$$      read (10, *) Step
c$$$      read (10, *) debug
      delta      = cg_Wolfe_c1   
      sigma      = cg_Wolfe_c2
      eps        = 1.d-6  
      gamma      = .66d0  
      rho        = 5.0d0  
      eta        = .01d0  
      psi0       = .01d0  
      psi1       = .1d0   
      psi2       = 2.d0   
      QuadCutOff = 1.d-12 
      StopFact   = 0.d-12 
      AWolfeFac  = 1.d-3  
      restart_fac= 1.0d0  
      maxit_fac  = 500.d0 
      feps       = 0.d0   
      Qdecay     = .7d0   
      nexpand    = 50     
      nsecant    = 50     
      PertRule   = .true. 
      QuadStep   = .true.
      if (myid==0 .and. IsQuiet<-1) then
         PrintLevel = .true.
      else 
         PrintLevel = .false.
      end if 
      if (myid==0) then
         PrintFinal = .true. 
      else 
         PrintFinal = .false. 
      end if 
      StopRule   = .true. 
      AWolfe     = .false.
      Step       = .false.
      if (myid==0 .and. IsQuiet<-1) then
         debug      = .true.
      else 
         debug      = .false.
      end if 
      nrestart = n*restart_fac
c      maxit = n*maxit_fac
      maxit = cg_maxiter
      zero = 0.d0
      info = 0
      n5 = mod (n, 5)
      n6 = n5 + 1
      nf = 0
      ng = 0
c     close (10)
      return
      end

c check whether the Wolfe or the approximate Wolfe conditions
c     are satisfied

c      (double) alpha   -- stepsize
c      (double) f       -- function value associated with stepsize alpha
c      (double) dphi    -- derivative value associated with stepsize alpha

      logical function cg_Wolfe (alpha, f, dphi)

      double precision delta, sigma, eps,
     &                 gamma, rho, tol, eta, fpert, f0, Ck, Qdecay,
     &                 wolfe_hi, wolfe_lo, awolfe_hi,
     &                 QuadCutOff, StopFac, AWolfeFac,
     &                 zero, feps, psi0, psi1, psi2,
     &                 alpha, f, dphi

      integer          n, n5, n6, nf, ng, info, nrestart,
     &                 nexpand, nsecant, maxit

      logical          PertRule, QuadOK, QuadStep, PrintLevel,
     &                 PrintFinal, StopRule, AWolfe, Step, debug

      common /cgparms/delta, sigma, eps,
     &                gamma, rho, tol, eta, fpert, f0, Ck, Qdecay,
     &                wolfe_hi, wolfe_lo, awolfe_hi,
     &                QuadCutOff, StopFac, AWolfeFac,
     &                zero, feps, psi0, psi1, psi2,
     &                n, n5, n6, nf, ng, info,
     &                nrestart, nexpand, nsecant, maxit,
     &                PertRule, QuadOK, QuadStep, PrintLevel,
     &                PrintFinal, StopRule, AWolfe, Step, debug

      if ( dphi .ge. wolfe_lo ) then

c test original Wolfe conditions

          if ( f-f0 .le. alpha*wolfe_hi ) then
              cg_Wolfe = .true.
              if ( PrintLevel) then
                 write (*, 10) f, f0, alpha*wolfe_hi, dphi
10                format (' wolfe f:', e14.6, ' f0: ',
     &                    e14.6, e14.6, ' dphi:', e14.6)
              endif
              return

c test approximate Wolfe conditions

          elseif ( AWolfe ) then
              if ( (f .le. fpert).and.(dphi .le. awolfe_hi) ) then
                  cg_Wolfe = .true.
                  if ( PrintLevel ) then
                     write (*, 20) f, fpert, dphi, awolfe_hi
20                        format ('f:', e14.6, ' fpert:', e14.6,
     &                            ' dphi: ', e14.6, ' fappx:', e14.6)
                  endif
                  return
              endif
          endif
      endif
      cg_Wolfe = .false.
      return
      end

c check for convergence of the cg iterations
c      (double) f       -- function value associated with stepsize
c      (double) gnorm   -- gradient (infinity) norm

      logical function cg_tol (f, gnorm)

      double precision delta, sigma, eps,
     &                 gamma, rho, tol, eta, fpert, f0, Ck, Qdecay,
     &                 wolfe_hi, wolfe_lo, awolfe_hi,
     &                 QuadCutOff, StopFac, AWolfeFac,
     &                 zero, feps, psi0, psi1, psi2,
     &                 f, gnorm

      integer          n, n5, n6, nf, ng, info, nrestart,
     &                 nexpand, nsecant, maxit

      logical          PertRule, QuadOK, QuadStep, PrintLevel,
     &                 PrintFinal, StopRule, AWolfe, Step, debug

      common /cgparms/delta, sigma, eps,
     &                gamma, rho, tol, eta, fpert, f0, Ck, Qdecay,
     &                wolfe_hi, wolfe_lo, awolfe_hi,
     &                QuadCutOff, StopFac, AWolfeFac,
     &                zero, feps, psi0, psi1, psi2,
     &                n, n5, n6, nf, ng, info,
     &                nrestart, nexpand, nsecant, maxit,
     &                PertRule, QuadOK, QuadStep, PrintLevel,
     &                PrintFinal, StopRule, AWolfe, Step, debug

      if ( StopRule ) then
          if ( gnorm .le. tol ) then
              cg_tol = .true.
              return
          endif
      else
          if ( gnorm .le. tol*(1.0 + dabs (f)) ) then
              cg_tol = .true.
              return
          endif
      endif
      cg_tol = .false.
      return
      end

c compute dot product of x and y, vectors of length n
c      (double) x       -- first vector
c      (double) y       -- second vector

      double precision function cg_dot (x, y)

      double precision delta, sigma, eps,
     &                 gamma, rho, tol, eta, fpert, f0, Ck, Qdecay,
     &                 wolfe_hi, wolfe_lo, awolfe_hi,
     &                 QuadCutOff, StopFac, AWolfeFac,
     &                 zero, feps, psi0, psi1, psi2,
     &                 x (*), y(*), t

      integer          n, n5, n6, nf, ng, info, nrestart,
     &                 nexpand, nsecant, maxit, i

      logical          PertRule, QuadOK, QuadStep, PrintLevel,
     &                 PrintFinal, StopRule, AWolfe, Step, debug

      common /cgparms/delta, sigma, eps,
     &                gamma, rho, tol, eta, fpert, f0, Ck, Qdecay,
     &                wolfe_hi, wolfe_lo, awolfe_hi,
     &                QuadCutOff, StopFac, AWolfeFac,
     &                zero, feps, psi0, psi1, psi2,
     &                n, n5, n6, nf, ng, info,
     &                nrestart, nexpand, nsecant, maxit,
     &                PertRule, QuadOK, QuadStep, PrintLevel,
     &                PrintFinal, StopRule, AWolfe, Step, debug

      t = zero
      do i = 1, n5
          t = t + x (i)*y (i)
      enddo
      do i = n6, n, 5
          t = t + x (i)*y(i) + x (i+1)*y (i+1) + x (i+2)*y (i+2)
     &                       + x (i+3)*y (i+3) + x (i+4)*y (i+4)
      enddo
      cg_dot = t
      return
      end

c
c  compute xtemp = x + alpha d
c
c      (double) xtemp   -- output vector
c      (double) x       -- initial vector
c      (double) d       -- search direction vector
c      (double) alpha   -- stepsize along search direction vector

      subroutine cg_step (xtemp, x, d, alpha)

      double precision delta, sigma, eps,
     &                 gamma, rho, tol, eta, fpert, f0, Ck, Qdecay,
     &                 wolfe_hi, wolfe_lo, awolfe_hi,
     &                 QuadCutOff, StopFac, AWolfeFac,
     &                 zero, feps, psi0, psi1, psi2,
     &                 xtemp (*), x (*), d (*), alpha

      integer          n, n5, n6, nf, ng, info, nrestart,
     &                 nexpand, nsecant, maxit, i, j

      logical          PertRule, QuadOK, QuadStep, PrintLevel,
     &                 PrintFinal, StopRule, AWolfe, Step, debug

      common /cgparms/delta, sigma, eps,
     &                gamma, rho, tol, eta, fpert, f0, Ck, Qdecay,
     &                wolfe_hi, wolfe_lo, awolfe_hi,
     &                QuadCutOff, StopFac, AWolfeFac,
     &                zero, feps, psi0, psi1, psi2,
     &                n, n5, n6, nf, ng, info,
     &                nrestart, nexpand, nsecant, maxit,
     &                PertRule, QuadOK, QuadStep, PrintLevel,
     &                PrintFinal, StopRule, AWolfe, Step, debug

      do i = 1, n5
          xtemp (i) = x(i) + alpha*d(i)
      enddo
      do i = n6, n, 5
          xtemp (i) = x (i) + alpha*d (i)
          j = i + 1
          xtemp (j) = x (j) + alpha*d (j)
          j = i + 2
          xtemp (j) = x (j) + alpha*d (j)
          j = i + 3
          xtemp (j) = x (j) + alpha*d (j)
          j = i + 4
          xtemp (j) = x (j) + alpha*d (j)
      enddo
      end

c      (double) alpha   -- stepsize along search direction vector
c      (double) phi     -- function value for step alpha
c      (double) dphi    -- function derivative for step alpha
c      (double) dphi0   -- function derivative at starting point (alpha = 0)
c      (double) x       -- current iterate
c      (double) xtemp   -- x + alpha*d
c      (double) d       -- current search direction
c      (double) gtemp   -- gradient at x + alpha*d
c      (external) cg_value -- routine to evaluate function value
c      (external) cg_grad  -- routine to evaluate function gradient

      subroutine cg_line (alpha, phi, dphi, dphi0, x, xtemp, d, gtemp,
     &                    cg_value, cg_grad)

      double precision delta, sigma, eps,
     &                 gamma, rho, tol, eta, fpert, f0, Ck, Qdecay,
     &                 wolfe_hi, wolfe_lo, awolfe_hi,
     &                 QuadCutOff, StopFac, AWolfeFac,
     &                 zero, feps, psi0, psi1, psi2,
     &                 x (*), xtemp (*), d (*), gtemp (*),
     &                 a, dphia, b, dphib, alpha, phi, dphi, c,
     &                 a0, da0, b0, db0, width, fquad, dphi0,
     &                 cg_dot

      integer          n, n5, n6, nf, ng, info, nrestart,
     &                 nexpand, nsecant, maxit,
     &                 ngrow, nshrink, cg_update, iter, flag

      logical          PertRule, QuadOK, QuadStep, PrintLevel,
     &                 PrintFinal, StopRule, AWolfe, Step, debug,
     &                 cg_Wolfe

      external         cg_value, cg_grad

      common /cgparms/delta, sigma, eps,
     &                gamma, rho, tol, eta, fpert, f0, Ck, Qdecay,
     &                wolfe_hi, wolfe_lo, awolfe_hi,
     &                QuadCutOff, StopFac, AWolfeFac,
     &                zero, feps, psi0, psi1, psi2,
     &                n, n5, n6, nf, ng, info,
     &                nrestart, nexpand, nsecant, maxit,
     &                PertRule, QuadOK, QuadStep, PrintLevel,
     &                PrintFinal, StopRule, AWolfe, Step, debug

      call cg_step (xtemp, x, d, alpha)
      call cg_grad (gtemp, xtemp, n) 
      ng = ng + 1
      dphi = cg_dot (gtemp, d)
c
c Find initial interval [a,b] such that dphia < 0, dphib >= 0,
c        and phia <= phi0 + feps*dabs (phi0)
c
      a = zero
      dphia = dphi0 
      ngrow = 0
      nshrink = 0
      do while ( dphi .lt. zero )
          call cg_value (phi, xtemp, n)
          nf = nf + 1
c
c if quadstep in effect and quadratic conditions hold, check wolfe condition
c
          if ( QuadOK ) then
              if ( ngrow .eq. 0 ) fquad = dmin1 (phi, f0)
              if ( phi .le. fquad ) then
                  if ( PrintLevel ) then
                      write (*, 10) alpha, phi, fquad
10                    format ('alpha:', e14.6, ' phi:', e14.6,
     &                        ' fquad:', e14.6)
                  endif
                  if ( cg_Wolfe (alpha, phi, dphi) ) return
              endif
          endif
          if ( phi .le. fpert ) then
              a = alpha
              dphia = dphi
          else
c
c contraction phase
c
              b = alpha
              do while ( .true. )
                  alpha = .5d0*(a+b)
                  nshrink = nshrink + 1
                  if ( nshrink .gt. nexpand ) then
                      info = 6
                      return
                  endif
                  call cg_step (xtemp, x, d, alpha)
                  call cg_grad (gtemp, xtemp, n) 
                  ng = ng + 1
                  dphi = cg_dot (gtemp, d)
                  if ( dphi .ge. zero ) goto 100
                  call cg_value (phi, xtemp, n)
                  nf = nf + 1
                  if ( PrintLevel ) then
                      write (6, 20) a, b, alpha, phi, dphi
20                    format ('contract, a:', e14.6,
     &                        ' b:', e14.6, ' alpha:', e14.6,
     &                        ' phi:', e14.6, ' dphi:', e14.6)
                  endif
                  if ( QuadOK .and. (phi .le. fquad) ) then
                      if ( cg_Wolfe (alpha, phi, dphi) ) return
                  endif
                  if ( phi .le. fpert ) then
                      a = alpha
                      dphia = dphi
                  else
                      b = alpha
                  endif
              enddo
          endif
c
c expansion phase
c
          ngrow = ngrow + 1
          if ( ngrow .gt. nexpand ) then
              info = 3
              return
          endif
          alpha = rho*alpha
          call cg_step (xtemp, x, d, alpha)
          call cg_grad (gtemp, xtemp, n) 
          ng = ng + 1
          dphi = cg_dot (gtemp, d)
          if ( PrintLevel ) then
              write (*, 30) a, alpha, phi, dphi
30            format ('expand,   a:', e14.6, ' alpha:', e14.6,
     &                 ' phi:', e14.6, ' dphi:', e14.6)
          endif
      enddo
100   continue
      b = alpha
      dphib = dphi
      if ( QuadOK ) then
          call cg_value (phi, xtemp, n)
          nf = nf + 1
          if ( ngrow + nshrink .eq. 0 ) fquad = dmin1 (phi, f0)
          if ( phi .le. fquad ) then
              if ( cg_Wolfe (alpha, phi, dphi) ) return
          endif
      endif
      do iter = 1, nsecant
          if ( PrintLevel ) then
              write (*, 40) a, b, dphia, dphib
40            format ('secant, a:', e14.6, ' b:', e14.6,
     &                 ' da:', e14.6, ' db:', e14.6)
          endif
          width = gamma*(b - a)
          if ( -dphia .le. dphib ) then
              alpha = a - (a-b)*(dphia/(dphia-dphib))
          else
              alpha = b - (a-b)*(dphib/(dphia-dphib))
          endif
          c = alpha
          a0 = a
          b0 = b
          da0 = dphia
          db0 = dphib
          flag = cg_update (a, dphia, b, dphib, alpha, phi,
     &              dphi, x, xtemp, d, gtemp, cg_value, cg_grad)
          if ( flag .gt. 0 ) then
              return
          else if ( flag .eq. 0 ) then
              if ( c .eq. a ) then
                  if ( dphi .gt. da0 ) then
                      alpha = c - (c-a0)*(dphi/(dphi-da0))
                  else
                      alpha = a
                  endif
              else
                  if ( dphi .lt. db0 ) then
                      alpha = c - (c-b0)*(dphi/(dphi-db0))
                  else
                      alpha = b
                  endif
              endif
              if ( (alpha .gt. a) .and. (alpha .lt. b) ) then
                  if ( PrintLevel ) write (*, *) "2nd secant"
                  flag = cg_update (a, dphia, b, dphib, alpha, phi,
     &                      dphi, x, xtemp, d, gtemp, cg_value, cg_grad)
                  if ( flag .gt. 0 ) return
              endif
          endif
c
c    bisection iteration
c  
          if ( (b-a) .ge. width ) then
              alpha = .5d0*(b+a)
              if ( PrintLevel ) write (*, *) "bisection"
              flag = cg_update (a, dphia, b, dphib, alpha, phi,
     &                  dphi, x, xtemp, d, gtemp, cg_value, cg_grad)
              if ( flag .gt. 0 ) return
          else
              if ( b .le. a ) then
                  info = 7
                  return
              endif
          endif
      end do
      info = 4
      return
      end

c  This routine is identical to cg_line except that the function
c  psi (a) = phi (a) - phi (0) - a*delta*dphi (0) is miniminized instead of
c  the function phi

c      (double) alpha   -- stepsize along search direction vector
c      (double) phi     -- function value for step alpha
c      (double) dphi    -- function derivative for step alpha
c      (double) dphi0   -- function derivative at starting point (alpha = 0)
c      (double) x       -- current iterate
c      (double) xtemp   -- x + alpha*d
c      (double) d       -- current search direction
c      (double) gtemp   -- gradient at x + alpha*d
c      (external) cg_value -- routine to evaluate function value
c      (external) cg_grad  -- routine to evaluate function gradient

      subroutine cg_lineW (alpha, phi, dphi, dphi0, x, xtemp, d, gtemp,
     &                    cg_value, cg_grad)

      double precision delta, sigma, eps,
     &                 gamma, rho, tol, eta, fpert, f0, Ck, Qdecay,
     &                 wolfe_hi, wolfe_lo, awolfe_hi,
     &                 QuadCutOff, StopFac, AWolfeFac,
     &                 zero, feps, psi0, psi1, psi2,
     &                 x (*), xtemp (*), d (*), gtemp (*),
     &                 a, dpsia, b, dpsib, alpha, phi, dphi, c,
     &                 a0, da0, b0, db0, width, fquad, dphi0,
     &                 cg_dot, psi, dpsi

      integer          n, n5, n6, nf, ng, info, nrestart,
     &                 nexpand, nsecant, maxit,
     &                 ngrow, nshrink, cg_updateW, iter, flag

      logical          PertRule, QuadOK, QuadStep, PrintLevel,
     &                 PrintFinal, StopRule, AWolfe, Step, debug,
     &                 cg_Wolfe

      external         cg_value, cg_grad

      common /cgparms/delta, sigma, eps,
     &                gamma, rho, tol, eta, fpert, f0, Ck, Qdecay,
     &                wolfe_hi, wolfe_lo, awolfe_hi,
     &                QuadCutOff, StopFac, AWolfeFac,
     &                zero, feps, psi0, psi1, psi2,
     &                n, n5, n6, nf, ng, info,
     &                nrestart, nexpand, nsecant, maxit,
     &                PertRule, QuadOK, QuadStep, PrintLevel,
     &                PrintFinal, StopRule, AWolfe, Step, debug

      call cg_step (xtemp, x, d, alpha)
      call cg_grad (gtemp, xtemp, n) 
      ng = ng + 1
      dphi = cg_dot (gtemp, d)
      dpsi = dphi - wolfe_hi
c
c Find initial interval [a,b] such that dpsia < 0, dpsib >= 0,
c        and psia <= phi0 + feps*dabs (phi0)
c
      a = zero
      dpsia = dphi0 - wolfe_hi
      ngrow = 0
      nshrink = 0
      do while ( dpsi .lt. zero )
          call cg_value (phi, xtemp, n)
          psi = phi - alpha*wolfe_hi
          
          nf = nf + 1
c
c if quadstep in effect and quadratic conditions hold, check wolfe condition
c
          if ( QuadOK ) then
              if ( ngrow .eq. 0 ) fquad = dmin1 (phi, f0)
              if ( phi .le. fquad ) then
                  if ( PrintLevel ) then
                      write (*, 10) alpha, phi, fquad
10                    format ('alpha:', e14.6, ' phi:', e14.6,
     &                        ' fquad:', e14.6)
                  endif
                  if ( cg_Wolfe (alpha, phi, dphi) ) return
              endif
          endif
          if ( psi .le. fpert ) then
              a = alpha
              dpsia = dpsi
          else
c
c contraction phase
c
              b = alpha
              do while ( .true. )
                  alpha = .5d0*(a+b)
                  nshrink = nshrink + 1
                  if ( nshrink .gt. nexpand ) then
                      info = 6
                      return
                  endif
                  call cg_step (xtemp, x, d, alpha)
                  call cg_grad (gtemp, xtemp, n) 
                  ng = ng + 1
                  dphi = cg_dot (gtemp, d)
                  dpsi = dphi - wolfe_hi
                  if ( dpsi .ge. zero ) goto 100
                  call cg_value (phi, xtemp, n)
                  psi = phi - alpha*wolfe_hi
                  nf = nf + 1
                  if ( PrintLevel ) then
                      write (6, 20) a, b, alpha, phi, dphi
20                    format ('contract, a:', e14.6,
     &                        ' b:', e14.6, ' alpha:', e14.6,
     &                        ' phi:', e14.6, ' dphi:', e14.6)
                  endif
                  if ( QuadOK .and. (phi .le. fquad) ) then
                      if ( cg_Wolfe (alpha, phi, dphi) ) return
                  endif
                  if ( psi .le. fpert ) then
                      a = alpha
                      dpsia = dpsi
                  else
                      b = alpha
                  endif
              enddo
          endif
c
c expansion phase
c
          ngrow = ngrow + 1
          if ( ngrow .gt. nexpand ) then
              info = 3
              return
          endif
          alpha = rho*alpha
          call cg_step (xtemp, x, d, alpha)
          call cg_grad (gtemp, xtemp, n) 
          ng = ng + 1
          dphi = cg_dot (gtemp, d)
          dpsi = dphi - wolfe_hi
          if ( PrintLevel ) then
              write (*, 30) a, alpha, phi, dphi
30            format ('expand,   a:', e14.6, ' alpha:', e14.6,
     &                 ' phi:', e14.6, ' dphi:', e14.6)
              write (6, *) "expand, alpha:", alpha, "dphi:", dphi
          endif
      enddo
100   continue
      b = alpha
      dpsib = dpsi
      if ( QuadOK ) then
          call cg_value (phi, xtemp, n)
          nf = nf + 1
          if ( ngrow + nshrink .eq. 0 ) fquad = dmin1 (phi, f0)
          if ( phi .le. fquad ) then
              if ( cg_Wolfe (alpha, phi, dphi) ) return
          endif
      endif
      do iter = 1, nsecant
          if ( PrintLevel ) then
              write (*, 40) a, b, dpsia, dpsib
40            format ('secant, a:', e14.6, ' b:', e14.6,
     &                 ' da:', e14.6, ' db:', e14.6)
          endif
          width = gamma*(b - a)
          if ( -dpsia .le. dpsib ) then
              alpha = a - (a-b)*(dpsia/(dpsia-dpsib))
          else
              alpha = b - (a-b)*(dpsib/(dpsia-dpsib))
          endif
          c = alpha
          a0 = a
          b0 = b
          da0 = dpsia
          db0 = dpsib
          flag = cg_updateW (a, dpsia, b, dpsib, alpha,
     &               phi, dphi, dpsi, x, xtemp, d, gtemp,
     &               cg_value, cg_grad)
          if ( flag .gt. 0 ) then
              return
          else if ( flag .eq. 0 ) then
              if ( c .eq. a ) then
                  if ( dpsi .gt. da0 ) then
                      alpha = c - (c-a0)*(dpsi/(dpsi-da0))
                  else
                      alpha = a
                  endif
              else
                  if ( dpsi .lt. db0 ) then
                      alpha = c - (c-b0)*(dpsi/(dpsi-db0))
                  else
                      alpha = b
                  endif
              endif
              if ( (alpha .gt. a) .and. (alpha .lt. b) ) then
                  if ( PrintLevel ) write (*, *) "2nd secant"
                  flag = cg_updateW (a, dpsia, b, dpsib, alpha,
     &                       phi, dphi, dpsi, x, xtemp, d, gtemp,
     &                       cg_value, cg_grad)
                  if ( flag .gt. 0 ) return
              endif
          endif
c
c    bisection iteration
c  
          if ( (b-a) .ge. width ) then
              alpha = .5d0*(b+a)
              if ( PrintLevel ) write (*, *) "bisection"
              flag = cg_updateW (a, dpsia, b, dpsib, alpha,
     &                   phi, dphi, dpsi, x, xtemp, d, gtemp,
     &                   cg_value, cg_grad)
              if ( flag .gt. 0 ) return
          else
              if ( b .le. a ) then
                  info = 7
                  return
              endif
          endif
      end do
      info = 4
      return
      end
c
c update returns 1 if Wolfe condition is satisfied or too many iterations
c        returns  0 if the interval updated successfully
c        returns -1 if search done
c
c      (double) a       -- left side of bracketting interval
c      (double) dphia   -- derivative at a
c      (double) b       -- right side of bracketting interval
c      (double) dphib   -- derivative at b
c      (double) alpha   -- trial step (between a and b)
c      (double) phi     -- function value at alpha (returned)
c      (double) dphi    -- function derivative at alpha (returned)
c      (double) x       -- current iterate
c      (double) xtemp   -- x + alpha*d
c      (double) d       -- current search direction
c      (double) gtemp   -- gradient at x + alpha*d
c      (external) cg_value -- routine to evaluate function value
c      (external) cg_grad  -- routine to evaluate function gradient

      integer function cg_update (a, dphia, b, dphib, alpha, phi,
     &                    dphi, x, xtemp, d, gtemp, cg_value, cg_grad)

      double precision delta, sigma, eps,
     &                 gamma, rho, tol, eta, fpert, f0, Ck, Qdecay,
     &                 wolfe_hi, wolfe_lo, awolfe_hi,
     &                 QuadCutOff, StopFac, AWolfeFac,
     &                 zero, feps, psi0, psi1, psi2,
     &                 a, dphia, b, dphib, alpha, phi, dphi,
     &                 x (*), xtemp (*), d (*), gtemp (*),
     &                 cg_dot

      integer          n, n5, n6, nf, ng, info, nrestart,
     &                 nexpand, nsecant, maxit,
     &                 nshrink

      logical          PertRule, QuadOK, QuadStep, PrintLevel,
     &                 PrintFinal, StopRule, AWolfe, Step, debug,
     &                 cg_Wolfe

      external         cg_value, cg_grad

      common /cgparms/delta, sigma, eps,
     &                gamma, rho, tol, eta, fpert, f0, Ck, Qdecay,
     &                wolfe_hi, wolfe_lo, awolfe_hi,
     &                QuadCutOff, StopFac, AWolfeFac,
     &                zero, feps, psi0, psi1, psi2,
     &                n, n5, n6, nf, ng, info,
     &                nrestart, nexpand, nsecant, maxit,
     &                PertRule, QuadOK, QuadStep, PrintLevel,
     &                PrintFinal, StopRule, AWolfe, Step, debug

      call cg_step (xtemp, x, d, alpha)
      call cg_value (phi, xtemp, n)
      nf = nf + 1
      call cg_grad (gtemp, xtemp, n) 
      ng = ng + 1
      dphi = cg_dot (gtemp, d)
      if ( PrintLevel ) then
          write (*, 10) alpha, phi, dphi
10        format ('update alpha:', e14.6, ' phi:', e14.6,
     &            ' dphi:', e14.6)
      endif
      cg_update = 0
      if ( cg_Wolfe (alpha, phi, dphi) ) then
          cg_update = 1
          goto 110
      endif
      if ( dphi .ge. zero ) then
          b = alpha
          dphib = dphi
          goto 110
      else
          if ( phi .le. fpert ) then
              a = alpha
              dphia = dphi
              goto 110
          endif
      endif
      nshrink = 0
      b = alpha
      do while ( .true. )
          alpha = .5d0*(a+b)
          nshrink = nshrink + 1
          if ( nshrink .gt. nexpand ) then
              info = 8
              cg_update = 1
              goto 110
          endif
          call cg_step (xtemp, x, d, alpha)
          call cg_grad (gtemp, xtemp, n) 
          ng = ng + 1
          dphi = cg_dot (gtemp, d)
          call cg_value (phi, xtemp, n)
          nf = nf + 1
          if ( PrintLevel ) then
              write (6, 20) a, alpha, phi, dphi
20            format ('contract, a:', e14.6, ' alpha:', e14.6,
     &                 ' phi:', e14.6, ' dphi:', e14.6)
          endif
          if ( cg_Wolfe (alpha, phi, dphi) ) then
              cg_update = 1
              goto 110
          endif
          if ( dphi .ge. zero ) then
              b = alpha
              dphib = dphi
              goto 100
          endif
          if ( phi .le. fpert ) then
              if ( PrintLevel ) then
                  write (6, *) "update a:", alpha, "dphia:", dphi
              endif
              a = alpha
              dphia = dphi
          else
              b = alpha
          endif
      enddo
100   continue
      cg_update = -1
110   continue
      if ( PrintLevel ) then
          write (*, 200) a, b, dphia, dphib, cg_update
200       format ('UP a:', e14.6, ' b:', e14.6,
     &             ' da:', e14.6, ' db:', e14.6, ' up:', i2)
      endif
      return
      end

c  This routine is identical to cg_update except that the function
c  psi (a) = phi (a) - phi (0) - a*delta*dphi (0) is miniminized instead of
c  the function phi
c
c update returns 1 if Wolfe condition is satisfied or too many iterations
c        returns  0 if the interval updated successfully
c        returns -1 if search done
c
c      (double) a       -- left side of bracketting interval
c      (double) dpsia   -- derivative at a
c      (double) b       -- right side of bracketting interval
c      (double) dpsib   -- derivative at b
c      (double) alpha   -- trial step (between a and b)
c      (double) phi     -- function value at alpha (returned)
c      (double) dphi    -- derivative of phi at alpha (returned)
c      (double) dpsi    -- derivative of psi at alpha (returned)
c      (double) x       -- current iterate
c      (double) xtemp   -- x + alpha*d
c      (double) d       -- current search direction
c      (double) gtemp   -- gradient at x + alpha*d
c      (external) cg_value -- routine to evaluate function value
c      (external) cg_grad  -- routine to evaluate function gradient

      integer function cg_updateW (a, dpsia, b, dpsib, alpha, phi, dphi,
     &                      dpsi, x, xtemp, d, gtemp, cg_value, cg_grad)

      double precision delta, sigma, eps,
     &                 gamma, rho, tol, eta, fpert, f0, Ck, Qdecay,
     &                 wolfe_hi, wolfe_lo, awolfe_hi,
     &                 QuadCutOff, StopFac, AWolfeFac,
     &                 zero, feps, psi0, psi1, psi2,
     &                 a, dpsia, b, dpsib, alpha, phi, dphi,
     &                 x (*), xtemp (*), d (*), gtemp (*),
     &                 cg_dot, psi, dpsi

      integer          n, n5, n6, nf, ng, info, nrestart,
     &                 nexpand, nsecant, maxit, nshrink

      logical          PertRule, QuadOK, QuadStep, PrintLevel,
     &                 PrintFinal, StopRule, AWolfe, Step, debug,
     &                 cg_Wolfe

      external         cg_value, cg_grad

      common /cgparms/delta, sigma, eps,
     &                gamma, rho, tol, eta, fpert, f0, Ck, Qdecay,
     &                wolfe_hi, wolfe_lo, awolfe_hi,
     &                QuadCutOff, StopFac, AWolfeFac,
     &                zero, feps, psi0, psi1, psi2,
     &                n, n5, n6, nf, ng, info,
     &                nrestart, nexpand, nsecant, maxit,
     &                PertRule, QuadOK, QuadStep, PrintLevel,
     &                PrintFinal, StopRule, AWolfe, Step, debug

      call cg_step (xtemp, x, d, alpha)
      call cg_value (phi, xtemp, n)
      psi = phi - alpha*wolfe_hi
      nf = nf + 1
      call cg_grad (gtemp, xtemp, n) 
      ng = ng + 1
      dphi = cg_dot (gtemp, d)
      dpsi = dphi - wolfe_hi
      if ( PrintLevel ) then
          write (*, 10) alpha, psi, dpsi
10        format ('update alpha:', e14.6, ' psi:', e14.6,
     &            ' dpsi:', e14.6)
      endif
      cg_updateW = 0
      if ( cg_Wolfe (alpha, phi, dphi) ) then
          cg_updateW = 1
          goto 110
      endif
      if ( dpsi .ge. zero ) then
          b = alpha
          dpsib = dpsi
          goto 110
      else
          if ( psi .le. fpert ) then
              a = alpha
              dpsia = dpsi
              goto 110
          endif
      endif
      nshrink = 0
      b = alpha
      do while ( .true. )
          alpha = .5d0*(a+b)
          nshrink = nshrink + 1
          if ( nshrink .gt. nexpand ) then
              info = 8
              cg_updateW = 1
              goto 110
          endif
          call cg_step (xtemp, x, d, alpha)
          call cg_grad (gtemp, xtemp, n) 
          ng = ng + 1
          dphi = cg_dot (gtemp, d)
          dpsi = dphi - wolfe_hi
          call cg_value (phi, xtemp, n)
          psi = phi - alpha*wolfe_hi
          nf = nf + 1
          if ( PrintLevel ) then
              write (6, 20) a, alpha, phi, dphi
20            format ('contract, a:', e14.6, ' alpha:', e14.6,
     &                 ' phi:', e14.6, ' dphi:', e14.6)
          endif
          if ( cg_Wolfe (alpha, phi, dphi) ) then
              cg_updateW = 1
              goto 110
          endif
          if ( dpsi .ge. zero ) then
              b = alpha
              dpsib = dpsi
              goto 100
          endif
          if ( psi .le. fpert ) then
              if ( PrintLevel ) then
                  write (6, *) "update a:", alpha, "dpsia:", dpsi
              endif
              a = alpha
              dpsia = dpsi
          else
              b = alpha
          endif
      enddo
100   continue
      cg_updateW = -1
110   continue
      if ( PrintLevel ) then
          write (*, 200) a, b, dpsia, dpsib, cg_updateW
200       format ('UP a:', e14.6, ' b:', e14.6,
     &             ' da:', e14.6, ' db:', e14.6, ' up:', i2)
      endif
      return
      end
c Version 1.2 Changes:
c
c   1. Fix problem with user specified initial step (overwriting step)
c   2. Change dphi to dpsi at lines 1228 and 1234 in cg_lineW
c   3. Add comment about how to compute dnorm2 by an update of previous dnorm2
c   4. In comment statements for cg_lineW and cg_updateW, insert "delta"
c      in definition of psi (a)
c   5. In dimension statements, change "(1)" to "(*)"

c Version 1.3 Changes:
c   1. Remove extraneous write in line 985 (same thing written out twice)
c   2. Remove the parameter theta from cg_descent.parm and from the code
c      (we use theta = .5 in the cg_update)

c Version 1.4 Change:
c   1. The variable dpsi needs to be included in the argument list for
c      subroutine updateW (update of a Wolfe line search)
