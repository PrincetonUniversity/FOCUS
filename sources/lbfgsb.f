c                                                                                      
c  L-BFGS-B is released under the “New BSD License” (aka “Modified BSD License”        
c  or “3-clause license”)                                                              
c  Please read attached file License.txt                                               
c                                        
      subroutine dpofa(a,lda,n,info)
      integer lda,n,info
      double precision a(lda,*)
c
c     dpofa factors a double precision symmetric positive definite
c     matrix.
c
c     dpofa is usually called by dpoco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c     (time for dpoco) = (1 + 18/n)*(time for dpofa) .
c
c     on entry
c
c        a       double precision(lda, n)
c                the symmetric matrix to be factored.  only the
c                diagonal and upper triangle are used.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix  r  so that  a = trans(r)*r
c                where  trans(r)  is the transpose.
c                the strict lower triangle is unaltered.
c                if  info .ne. 0 , the factorization is not complete.
c
c        info    integer
c                = 0  for normal return.
c                = k  signals an error condition.  the leading minor
c                     of order  k  is not positive definite.
c
c     linpack.  this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas ddot
c     fortran sqrt
c
c     internal variables
c
      double precision ddot,t
      double precision s
      integer j,jm1,k
c     begin block with ...exits to 40
c
c
         do 30 j = 1, n
            info = j
            s = 0.0d0
            jm1 = j - 1
            if (jm1 .lt. 1) go to 20
            do 10 k = 1, jm1
               t = a(k,j) - ddot(k-1,a(1,k),1,a(1,j),1)
               t = t/a(k,k)
               a(k,j) = t
               s = s + t*t
   10       continue
   20       continue
            s = a(j,j) - s
c     ......exit
            if (s .le. 0.0d0) go to 40
            a(j,j) = sqrt(s)
   30    continue
         info = 0
   40 continue
      return
      end
      
c====================== The end of dpofa ===============================

      subroutine dtrsl(t,ldt,n,b,job,info)
      integer ldt,n,job,info
      double precision t(ldt,*),b(*)
c
c
c     dtrsl solves systems of the form
c
c                   t * x = b
c     or
c                   trans(t) * x = b
c
c     where t is a triangular matrix of order n. here trans(t)
c     denotes the transpose of the matrix t.
c
c     on entry
c
c         t         double precision(ldt,n)
c                   t contains the matrix of the system. the zero
c                   elements of the matrix are not referenced, and
c                   the corresponding elements of the array can be
c                   used to store other information.
c
c         ldt       integer
c                   ldt is the leading dimension of the array t.
c
c         n         integer
c                   n is the order of the system.
c
c         b         double precision(n).
c                   b contains the right hand side of the system.
c
c         job       integer
c                   job specifies what kind of system is to be solved.
c                   if job is
c
c                        00   solve t*x=b, t lower triangular,
c                        01   solve t*x=b, t upper triangular,
c                        10   solve trans(t)*x=b, t lower triangular,
c                        11   solve trans(t)*x=b, t upper triangular.
c
c     on return
c
c         b         b contains the solution, if info .eq. 0.
c                   otherwise b is unaltered.
c
c         info      integer
c                   info contains zero if the system is nonsingular.
c                   otherwise info contains the index of
c                   the first zero diagonal element of t.
c
c     linpack. this version dated 08/14/78 .
c     g. w. stewart, university of maryland, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,ddot
c     fortran mod
c
c     internal variables
c
      double precision ddot,temp
      integer case,j,jj
c
c     begin block permitting ...exits to 150
c
c        check for zero diagonal elements.
c
         do 10 info = 1, n
c     ......exit
            if (t(info,info) .eq. 0.0d0) go to 150
   10    continue
         info = 0
c
c        determine the task and go to it.
c
         case = 1
         if (mod(job,10) .ne. 0) case = 2
         if (mod(job,100)/10 .ne. 0) case = case + 2
         go to (20,50,80,110), case
c
c        solve t*x=b for t lower triangular
c
   20    continue
            b(1) = b(1)/t(1,1)
            if (n .lt. 2) go to 40
            do 30 j = 2, n
               temp = -b(j-1)
               call daxpy(n-j+1,temp,t(j,j-1),1,b(j),1)
               b(j) = b(j)/t(j,j)
   30       continue
   40       continue
         go to 140
c
c        solve t*x=b for t upper triangular.
c
   50    continue
            b(n) = b(n)/t(n,n)
            if (n .lt. 2) go to 70
            do 60 jj = 2, n
               j = n - jj + 1
               temp = -b(j+1)
               call daxpy(j,temp,t(1,j+1),1,b(1),1)
               b(j) = b(j)/t(j,j)
   60       continue
   70       continue
         go to 140
c
c        solve trans(t)*x=b for t lower triangular.
c
   80    continue
            b(n) = b(n)/t(n,n)
            if (n .lt. 2) go to 100
            do 90 jj = 2, n
               j = n - jj + 1
               b(j) = b(j) - ddot(jj-1,t(j+1,j),1,b(j+1),1)
               b(j) = b(j)/t(j,j)
   90       continue
  100       continue
         go to 140
c
c        solve trans(t)*x=b for t upper triangular.
c
  110    continue
            b(1) = b(1)/t(1,1)
            if (n .lt. 2) go to 130
            do 120 j = 2, n
               b(j) = b(j) - ddot(j-1,t(1,j),1,b(1),1)
               b(j) = b(j)/t(j,j)
  120       continue
  130       continue
  140    continue
  150 continue
      return
      end
      
c====================== The end of dtrsl ===============================


c                                                                                      
c  L-BFGS-B is released under the “New BSD License” (aka “Modified BSD License”        
c  or “3-clause license”)                                                              
c  Please read attached file License.txt                                               
c                                        

      double precision function dnrm2(n,x,incx)
      integer n,incx
      double precision x(n)
c     **********
c
c     Function dnrm2
c
c     Given a vector x of length n, this function calculates the
c     Euclidean norm of x with stride incx.
c
c     The function statement is
c
c       double precision function dnrm2(n,x,incx)
c
c     where
c
c       n is a positive integer input variable.
c
c       x is an input array of length n.
c
c       incx is a positive integer variable that specifies the 
c         stride of the vector.
c
c     Subprograms called
c
c       FORTRAN-supplied ... abs, max, sqrt
c
c     MINPACK-2 Project. February 1991.
c     Argonne National Laboratory.
c     Brett M. Averick.
c
c     **********
      integer i
      double precision scale

      dnrm2 = 0.0d0
      scale = 0.0d0

      do 10 i = 1, n, incx
         scale = max(scale, abs(x(i)))
   10 continue

      if (scale .eq. 0.0d0) return

      do 20 i = 1, n, incx
         dnrm2 = dnrm2 + (x(i)/scale)**2
   20 continue

      dnrm2 = scale*sqrt(dnrm2)

 
      return

      end
      
c====================== The end of dnrm2 ===============================

      subroutine daxpy(n,da,dx,incx,dy,incy)
c
c     constant times a vector plus a vector.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(*),dy(*),da
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if (da .eq. 0.0d0) return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dy(i) + da*dx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
        dy(i) = dy(i) + da*dx(i)
        dy(i + 1) = dy(i + 1) + da*dx(i + 1)
        dy(i + 2) = dy(i + 2) + da*dx(i + 2)
        dy(i + 3) = dy(i + 3) + da*dx(i + 3)
   50 continue
      return
      end
      
c====================== The end of daxpy ===============================

      subroutine dcopy(n,dx,incx,dy,incy)
c
c     copies a vector, x, to a vector, y.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(*),dy(*)
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,7)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dx(i)
   30 continue
      if( n .lt. 7 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,7
        dy(i) = dx(i)
        dy(i + 1) = dx(i + 1)
        dy(i + 2) = dx(i + 2)
        dy(i + 3) = dx(i + 3)
        dy(i + 4) = dx(i + 4)
        dy(i + 5) = dx(i + 5)
        dy(i + 6) = dx(i + 6)
   50 continue
      return
      end
      
c====================== The end of dcopy ===============================

      double precision function ddot(n,dx,incx,dy,incy)
c
c     forms the dot product of two vectors.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(*),dy(*),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
c
      ddot = 0.0d0
      dtemp = 0.0d0
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = dtemp + dx(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      ddot = dtemp
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dtemp + dx(i)*dy(i)
   30 continue
      if( n .lt. 5 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) +
     *   dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
   50 continue
   60 ddot = dtemp
      return
      end
      
c====================== The end of ddot ================================

      subroutine  dscal(n,da,dx,incx)
c
c     scales a vector by a constant.
c     uses unrolled loops for increment equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c
      double precision da,dx(*)
      integer i,incx,m,mp1,n,nincx
c
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        dx(i) = da*dx(i)
   10 continue
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dx(i) = da*dx(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dx(i) = da*dx(i)
        dx(i + 1) = da*dx(i + 1)
        dx(i + 2) = da*dx(i + 2)
        dx(i + 3) = da*dx(i + 3)
        dx(i + 4) = da*dx(i + 4)
   50 continue
      return
      end
      
c====================== The end of dscal ===============================

c                                                                                      
c  L-BFGS-B is released under the “New BSD License” (aka “Modified BSD License”        
c  or “3-clause license”)                                                              
c  Please read attached file License.txt                                               
c                                        
      subroutine timer(ttime)
      double precision ttime
c
      real temp
c
c     This routine computes cpu time in double precision; it makes use of 
c     the intrinsic f90 cpu_time therefore a conversion type is
c     needed.
c
c           J.L Morales  Departamento de Matematicas, 
c                        Instituto Tecnologico Autonomo de Mexico
c                        Mexico D.F.
c
c           J.L Nocedal  Department of Electrical Engineering and
c                        Computer Science.
c                        Northwestern University. Evanston, IL. USA
c                         
c                        January 21, 2011
c
      temp = sngl(ttime)
      call cpu_time(temp)
      ttime = dble(temp) 

      return

      end
      
c                                                                                      
c  L-BFGS-B is released under the “New BSD License” (aka “Modified BSD License”        
c  or “3-clause license”)                                                              
c  Please read attached file License.txt                                               
c                                        
c===========   L-BFGS-B (version 3.0.  April 25, 2011  ===================
c
c     This is a modified version of L-BFGS-B. Minor changes in the updated 
c     code appear preceded by a line comment as follows 
c  
c     c-jlm-jn 
c
c     Major changes are described in the accompanying paper:
c
c         Jorge Nocedal and Jose Luis Morales, Remark on "Algorithm 778: 
c         L-BFGS-B: Fortran Subroutines for Large-Scale Bound Constrained 
c         Optimization"  (2011). To appear in  ACM Transactions on 
c         Mathematical Software,
c
c     The paper describes an improvement and a correction to Algorithm 778. 
c     It is shown that the performance of the algorithm can be improved 
c     significantly by making a relatively simple modication to the subspace 
c     minimization phase. The correction concerns an error caused by the use 
c     of routine dpmeps to estimate machine precision. 
c
c     The total work space **wa** required by the new version is 
c 
c                  2*m*n + 11m*m + 5*n + 8*m 
c
c     the old version required 
c
c                  2*m*n + 12m*m + 4*n + 12*m 
c
c
c            J. Nocedal  Department of Electrical Engineering and
c                        Computer Science.
c                        Northwestern University. Evanston, IL. USA
c
c
c           J.L Morales  Departamento de Matematicas, 
c                        Instituto Tecnologico Autonomo de Mexico
c                        Mexico D.F. Mexico.
c
c                        March  2011    
c                                                 
c============================================================================= 
      subroutine setulb(n, m, x, l, u, nbd, f, g, factr, pgtol, wa, iwa,
     +                 task, iprint, csave, lsave, isave, dsave)
 
      character*60     task, csave
      logical          lsave(4)
      integer          n, m, iprint, 
     +                 nbd(n), iwa(3*n), isave(44)
      double precision f, factr, pgtol, x(n), l(n), u(n), g(n),
c
c-jlm-jn
     +                 wa(2*m*n + 5*n + 11*m*m + 8*m), dsave(29)
 
c     ************
c
c     Subroutine setulb
c
c     This subroutine partitions the working arrays wa and iwa, and 
c       then uses the limited memory BFGS method to solve the bound
c       constrained optimization problem by calling mainlb.
c       (The direct method will be used in the subspace minimization.)
c
c     n is an integer variable.
c       On entry n is the dimension of the problem.
c       On exit n is unchanged.
c
c     m is an integer variable.
c       On entry m is the maximum number of variable metric corrections
c         used to define the limited memory matrix.
c       On exit m is unchanged.
c
c     x is a double precision array of dimension n.
c       On entry x is an approximation to the solution.
c       On exit x is the current approximation.
c
c     l is a double precision array of dimension n.
c       On entry l is the lower bound on x.
c       On exit l is unchanged.
c
c     u is a double precision array of dimension n.
c       On entry u is the upper bound on x.
c       On exit u is unchanged.
c
c     nbd is an integer array of dimension n.
c       On entry nbd represents the type of bounds imposed on the
c         variables, and must be specified as follows:
c         nbd(i)=0 if x(i) is unbounded,
c                1 if x(i) has only a lower bound,
c                2 if x(i) has both lower and upper bounds, and
c                3 if x(i) has only an upper bound.
c       On exit nbd is unchanged.
c
c     f is a double precision variable.
c       On first entry f is unspecified.
c       On final exit f is the value of the function at x.
c
c     g is a double precision array of dimension n.
c       On first entry g is unspecified.
c       On final exit g is the value of the gradient at x.
c
c     factr is a double precision variable.
c       On entry factr >= 0 is specified by the user.  The iteration
c         will stop when
c
c         (f^k - f^{k+1})/max{|f^k|,|f^{k+1}|,1} <= factr*epsmch
c
c         where epsmch is the machine precision, which is automatically
c         generated by the code. Typical values for factr: 1.d+12 for
c         low accuracy; 1.d+7 for moderate accuracy; 1.d+1 for extremely
c         high accuracy.
c       On exit factr is unchanged.
c
c     pgtol is a double precision variable.
c       On entry pgtol >= 0 is specified by the user.  The iteration
c         will stop when
c
c                 max{|proj g_i | i = 1, ..., n} <= pgtol
c
c         where pg_i is the ith component of the projected gradient.   
c       On exit pgtol is unchanged.
c
c     wa is a double precision working array of length 
c       (2mmax + 5)nmax + 12mmax^2 + 12mmax.
c
c     iwa is an integer working array of length 3nmax.
c
c     task is a working string of characters of length 60 indicating
c       the current job when entering and quitting this subroutine.
c
c     iprint is an integer variable that must be set by the user.
c       It controls the frequency and type of output generated:
c        iprint<0    no output is generated;
c        iprint=0    print only one line at the last iteration;
c        0<iprint<99 print also f and |proj g| every iprint iterations;
c        iprint=99   print details of every iteration except n-vectors;
c        iprint=100  print also the changes of active set and final x;
c        iprint>100  print details of every iteration including x and g;
c       When iprint > 0, the file iterate.dat will be created to
c                        summarize the iteration.
c
c     csave is a working string of characters of length 60.
c
c     lsave is a logical working array of dimension 4.
c       On exit with 'task' = NEW_X, the following information is 
c                                                             available:
c         If lsave(1) = .true.  then  the initial X has been replaced by
c                                     its projection in the feasible set;
c         If lsave(2) = .true.  then  the problem is constrained;
c         If lsave(3) = .true.  then  each variable has upper and lower
c                                     bounds;
c
c     isave is an integer working array of dimension 44.
c       On exit with 'task' = NEW_X, the following information is 
c                                                             available:
c         isave(22) = the total number of intervals explored in the 
c                         search of Cauchy points;
c         isave(26) = the total number of skipped BFGS updates before 
c                         the current iteration;
c         isave(30) = the number of current iteration;
c         isave(31) = the total number of BFGS updates prior the current
c                         iteration;
c         isave(33) = the number of intervals explored in the search of
c                         Cauchy point in the current iteration;
c         isave(34) = the total number of function and gradient 
c                         evaluations;
c         isave(36) = the number of function value or gradient
c                                  evaluations in the current iteration;
c         if isave(37) = 0  then the subspace argmin is within the box;
c         if isave(37) = 1  then the subspace argmin is beyond the box;
c         isave(38) = the number of free variables in the current
c                         iteration;
c         isave(39) = the number of active constraints in the current
c                         iteration;
c         n + 1 - isave(40) = the number of variables leaving the set of
c                           active constraints in the current iteration;
c         isave(41) = the number of variables entering the set of active
c                         constraints in the current iteration.
c
c     dsave is a double precision working array of dimension 29.
c       On exit with 'task' = NEW_X, the following information is
c                                                             available:
c         dsave(1) = current 'theta' in the BFGS matrix;
c         dsave(2) = f(x) in the previous iteration;
c         dsave(3) = factr*epsmch;
c         dsave(4) = 2-norm of the line search direction vector;
c         dsave(5) = the machine precision epsmch generated by the code;
c         dsave(7) = the accumulated time spent on searching for
c                                                         Cauchy points;
c         dsave(8) = the accumulated time spent on
c                                                 subspace minimization;
c         dsave(9) = the accumulated time spent on line search;
c         dsave(11) = the slope of the line search function at
c                                  the current point of line search;
c         dsave(12) = the maximum relative step length imposed in
c                                                           line search;
c         dsave(13) = the infinity norm of the projected gradient;
c         dsave(14) = the relative step length in the line search;
c         dsave(15) = the slope of the line search function at
c                                 the starting point of the line search;
c         dsave(16) = the square of the 2-norm of the line search
c                                                      direction vector.
c
c     Subprograms called:
c
c       L-BFGS-B Library ... mainlb.    
c
c
c     References:
c
c       [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited
c       memory algorithm for bound constrained optimization'',
c       SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.
c
c       [2] C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, ``L-BFGS-B: a
c       limited memory FORTRAN code for solving bound constrained
c       optimization problems'', Tech. Report, NAM-11, EECS Department,
c       Northwestern University, 1994.
c
c       (Postscript files of these papers are available via anonymous
c        ftp to eecs.nwu.edu in the directory pub/lbfgs/lbfgs_bcm.)
c
c                           *  *  *
c
c     NEOS, November 1994. (Latest revision June 1996.)
c     Optimization Technology Center.
c     Argonne National Laboratory and Northwestern University.
c     Written by
c                        Ciyou Zhu
c     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
c
c
c     ************
c-jlm-jn 
      integer   lws,lr,lz,lt,ld,lxp,lwa,
     +          lwy,lsy,lss,lwt,lwn,lsnd

      if (task .eq. 'START') then
         isave(1)  = m*n
         isave(2)  = m**2
         isave(3)  = 4*m**2
         isave(4)  = 1                      ! ws      m*n
         isave(5)  = isave(4)  + isave(1)   ! wy      m*n
         isave(6)  = isave(5)  + isave(1)   ! wsy     m**2
         isave(7)  = isave(6)  + isave(2)   ! wss     m**2
         isave(8)  = isave(7)  + isave(2)   ! wt      m**2
         isave(9)  = isave(8)  + isave(2)   ! wn      4*m**2
         isave(10) = isave(9)  + isave(3)   ! wsnd    4*m**2
         isave(11) = isave(10) + isave(3)   ! wz      n
         isave(12) = isave(11) + n          ! wr      n
         isave(13) = isave(12) + n          ! wd      n
         isave(14) = isave(13) + n          ! wt      n
         isave(15) = isave(14) + n          ! wxp     n
         isave(16) = isave(15) + n          ! wa      8*m
      endif
      lws  = isave(4)
      lwy  = isave(5)
      lsy  = isave(6)
      lss  = isave(7)
      lwt  = isave(8)
      lwn  = isave(9)
      lsnd = isave(10)
      lz   = isave(11)
      lr   = isave(12)
      ld   = isave(13)
      lt   = isave(14)
      lxp  = isave(15)
      lwa  = isave(16)

      call mainlb(n,m,x,l,u,nbd,f,g,factr,pgtol,
     +  wa(lws),wa(lwy),wa(lsy),wa(lss), wa(lwt),
     +  wa(lwn),wa(lsnd),wa(lz),wa(lr),wa(ld),wa(lt),wa(lxp),
     +  wa(lwa),
     +  iwa(1),iwa(n+1),iwa(2*n+1),task,iprint, 
     +  csave,lsave,isave(22),dsave)

      return

      end

c======================= The end of setulb =============================
 
      subroutine mainlb(n, m, x, l, u, nbd, f, g, factr, pgtol, ws, wy,
     +                  sy, ss, wt, wn, snd, z, r, d, t, xp, wa, 
     +                  index, iwhere, indx2, task,
     +                  iprint, csave, lsave, isave, dsave)
      implicit none
      character*60     task, csave
      logical          lsave(4)
      integer          n, m, iprint, nbd(n), index(n),
     +                 iwhere(n), indx2(n), isave(23)
      double precision f, factr, pgtol,
     +                 x(n), l(n), u(n), g(n), z(n), r(n), d(n), t(n), 
c-jlm-jn
     +                 xp(n), 
     +                 wa(8*m), 
     +                 ws(n, m), wy(n, m), sy(m, m), ss(m, m), 
     +                 wt(m, m), wn(2*m, 2*m), snd(2*m, 2*m), dsave(29)

c     ************
c
c     Subroutine mainlb
c
c     This subroutine solves bound constrained optimization problems by
c       using the compact formula of the limited memory BFGS updates.
c       
c     n is an integer variable.
c       On entry n is the number of variables.
c       On exit n is unchanged.
c
c     m is an integer variable.
c       On entry m is the maximum number of variable metric
c          corrections allowed in the limited memory matrix.
c       On exit m is unchanged.
c
c     x is a double precision array of dimension n.
c       On entry x is an approximation to the solution.
c       On exit x is the current approximation.
c
c     l is a double precision array of dimension n.
c       On entry l is the lower bound of x.
c       On exit l is unchanged.
c
c     u is a double precision array of dimension n.
c       On entry u is the upper bound of x.
c       On exit u is unchanged.
c
c     nbd is an integer array of dimension n.
c       On entry nbd represents the type of bounds imposed on the
c         variables, and must be specified as follows:
c         nbd(i)=0 if x(i) is unbounded,
c                1 if x(i) has only a lower bound,
c                2 if x(i) has both lower and upper bounds,
c                3 if x(i) has only an upper bound.
c       On exit nbd is unchanged.
c
c     f is a double precision variable.
c       On first entry f is unspecified.
c       On final exit f is the value of the function at x.
c
c     g is a double precision array of dimension n.
c       On first entry g is unspecified.
c       On final exit g is the value of the gradient at x.
c
c     factr is a double precision variable.
c       On entry factr >= 0 is specified by the user.  The iteration
c         will stop when
c
c         (f^k - f^{k+1})/max{|f^k|,|f^{k+1}|,1} <= factr*epsmch
c
c         where epsmch is the machine precision, which is automatically
c         generated by the code.
c       On exit factr is unchanged.
c
c     pgtol is a double precision variable.
c       On entry pgtol >= 0 is specified by the user.  The iteration
c         will stop when
c
c                 max{|proj g_i | i = 1, ..., n} <= pgtol
c
c         where pg_i is the ith component of the projected gradient.
c       On exit pgtol is unchanged.
c
c     ws, wy, sy, and wt are double precision working arrays used to
c       store the following information defining the limited memory
c          BFGS matrix:
c          ws, of dimension n x m, stores S, the matrix of s-vectors;
c          wy, of dimension n x m, stores Y, the matrix of y-vectors;
c          sy, of dimension m x m, stores S'Y;
c          ss, of dimension m x m, stores S'S;
c          yy, of dimension m x m, stores Y'Y;
c          wt, of dimension m x m, stores the Cholesky factorization
c                                  of (theta*S'S+LD^(-1)L'); see eq.
c                                  (2.26) in [3].
c
c     wn is a double precision working array of dimension 2m x 2m
c       used to store the LEL^T factorization of the indefinite matrix
c                 K = [-D -Y'ZZ'Y/theta     L_a'-R_z'  ]
c                     [L_a -R_z           theta*S'AA'S ]
c
c       where     E = [-I  0]
c                     [ 0  I]
c
c     snd is a double precision working array of dimension 2m x 2m
c       used to store the lower triangular part of
c                 N = [Y' ZZ'Y   L_a'+R_z']
c                     [L_a +R_z  S'AA'S   ]
c            
c     z(n),r(n),d(n),t(n), xp(n),wa(8*m) are double precision working arrays.
c       z  is used at different times to store the Cauchy point and
c          the Newton point.
c       xp is used to safeguard the projected Newton direction
c
c     sg(m),sgo(m),yg(m),ygo(m) are double precision working arrays. 
c
c     index is an integer working array of dimension n.
c       In subroutine freev, index is used to store the free and fixed
c          variables at the Generalized Cauchy Point (GCP).
c
c     iwhere is an integer working array of dimension n used to record
c       the status of the vector x for GCP computation.
c       iwhere(i)=0 or -3 if x(i) is free and has bounds,
c                 1       if x(i) is fixed at l(i), and l(i) .ne. u(i)
c                 2       if x(i) is fixed at u(i), and u(i) .ne. l(i)
c                 3       if x(i) is always fixed, i.e.,  u(i)=x(i)=l(i)
c                -1       if x(i) is always free, i.e., no bounds on it.
c
c     indx2 is an integer working array of dimension n.
c       Within subroutine cauchy, indx2 corresponds to the array iorder.
c       In subroutine freev, a list of variables entering and leaving
c       the free set is stored in indx2, and it is passed on to
c       subroutine formk with this information.
c
c     task is a working string of characters of length 60 indicating
c       the current job when entering and leaving this subroutine.
c
c     iprint is an INTEGER variable that must be set by the user.
c       It controls the frequency and type of output generated:
c        iprint<0    no output is generated;
c        iprint=0    print only one line at the last iteration;
c        0<iprint<99 print also f and |proj g| every iprint iterations;
c        iprint=99   print details of every iteration except n-vectors;
c        iprint=100  print also the changes of active set and final x;
c        iprint>100  print details of every iteration including x and g;
c       When iprint > 0, the file iterate.dat will be created to
c                        summarize the iteration.
c
c     csave is a working string of characters of length 60.
c
c     lsave is a logical working array of dimension 4.
c
c     isave is an integer working array of dimension 23.
c
c     dsave is a double precision working array of dimension 29.
c
c
c     Subprograms called
c
c       L-BFGS-B Library ... cauchy, subsm, lnsrlb, formk, 
c
c        errclb, prn1lb, prn2lb, prn3lb, active, projgr,
c
c        freev, cmprlb, matupd, formt.
c
c       Minpack2 Library ... timer
c
c       Linpack Library ... dcopy, ddot.
c
c
c     References:
c
c       [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited
c       memory algorithm for bound constrained optimization'',
c       SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.
c
c       [2] C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, ``L-BFGS-B: FORTRAN
c       Subroutines for Large Scale Bound Constrained Optimization''
c       Tech. Report, NAM-11, EECS Department, Northwestern University,
c       1994.
c 
c       [3] R. Byrd, J. Nocedal and R. Schnabel "Representations of
c       Quasi-Newton Matrices and their use in Limited Memory Methods'',
c       Mathematical Programming 63 (1994), no. 4, pp. 129-156.
c
c       (Postscript files of these papers are available via anonymous
c        ftp to eecs.nwu.edu in the directory pub/lbfgs/lbfgs_bcm.)
c
c                           *  *  *
c
c     NEOS, November 1994. (Latest revision June 1996.)
c     Optimization Technology Center.
c     Argonne National Laboratory and Northwestern University.
c     Written by
c                        Ciyou Zhu
c     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
c
c
c     ************
 
      logical          prjctd,cnstnd,boxed,updatd,wrk
      character*3      word
      integer          i,k,nintol,itfile,iback,nskip,
     +                 head,col,iter,itail,iupdat,
     +                 nseg,nfgv,info,ifun,
     +                 iword,nfree,nact,ileave,nenter
      double precision theta,fold,ddot,dr,rr,tol,
     +                 xstep,sbgnrm,ddum,dnorm,dtd,epsmch,
     +                 cpu1,cpu2,cachyt,sbtime,lnscht,time1,time2,
     +                 gd,gdold,stp,stpmx,time
      double precision one,zero
      parameter        (one=1.0d0,zero=0.0d0)
      
      if (task .eq. 'START') then

         epsmch = epsilon(one)

         call timer(time1)

c        Initialize counters and scalars when task='START'.

c           for the limited memory BFGS matrices:
         col    = 0
         head   = 1
         theta  = one
         iupdat = 0
         updatd = .false.
         iback  = 0
         itail  = 0
         iword  = 0
         nact   = 0
         ileave = 0
         nenter = 0
         fold   = zero
         dnorm  = zero
         cpu1   = zero
         gd     = zero
         stpmx  = zero
         sbgnrm = zero
         stp    = zero
         gdold  = zero
         dtd    = zero

c           for operation counts:
         iter   = 0
         nfgv   = 0
         nseg   = 0
         nintol = 0
         nskip  = 0
         nfree  = n
         ifun   = 0
c           for stopping tolerance:
         tol = factr*epsmch

c           for measuring running time:
         cachyt = 0
         sbtime = 0
         lnscht = 0
 
c           'word' records the status of subspace solutions.
         word = '---'

c           'info' records the termination information.
         info = 0

         itfile = 8
         if (iprint .ge. 1) then
c                                open a summary file 'iterate.dat'
            open (8, file = 'iterate.dat', status = 'unknown')
         endif            

c        Check the input arguments for errors.

         call errclb(n,m,factr,l,u,nbd,task,info,k)
         if (task(1:5) .eq. 'ERROR') then
            call prn3lb(n,x,f,task,iprint,info,itfile,
     +                  iter,nfgv,nintol,nskip,nact,sbgnrm,
     +                  zero,nseg,word,iback,stp,xstep,k,
     +                  cachyt,sbtime,lnscht)
            return
         endif

         call prn1lb(n,m,l,u,x,iprint,itfile,epsmch)
 
c        Initialize iwhere & project x onto the feasible set.
 
         call active(n,l,u,nbd,x,iwhere,iprint,prjctd,cnstnd,boxed) 

c        The end of the initialization.

      else
c          restore local variables.

         prjctd = lsave(1)
         cnstnd = lsave(2)
         boxed  = lsave(3)
         updatd = lsave(4)

         nintol = isave(1)
         itfile = isave(3)
         iback  = isave(4)
         nskip  = isave(5)
         head   = isave(6)
         col    = isave(7)
         itail  = isave(8)
         iter   = isave(9)
         iupdat = isave(10)
         nseg   = isave(12)
         nfgv   = isave(13)
         info   = isave(14)
         ifun   = isave(15)
         iword  = isave(16)
         nfree  = isave(17)
         nact   = isave(18)
         ileave = isave(19)
         nenter = isave(20)

         theta  = dsave(1)
         fold   = dsave(2)
         tol    = dsave(3)
         dnorm  = dsave(4)
         epsmch = dsave(5)
         cpu1   = dsave(6)
         cachyt = dsave(7)
         sbtime = dsave(8)
         lnscht = dsave(9)
         time1  = dsave(10)
         gd     = dsave(11)
         stpmx  = dsave(12)
         sbgnrm = dsave(13)
         stp    = dsave(14)
         gdold  = dsave(15)
         dtd    = dsave(16)
   
c        After returning from the driver go to the point where execution
c        is to resume.

         if (task(1:5) .eq. 'FG_LN') goto 666
         if (task(1:5) .eq. 'NEW_X') goto 777
         if (task(1:5) .eq. 'FG_ST') goto 111
         if (task(1:4) .eq. 'STOP') then
            if (task(7:9) .eq. 'CPU') then
c                                          restore the previous iterate.
               call dcopy(n,t,1,x,1)
               call dcopy(n,r,1,g,1)
               f = fold
            endif
            goto 999
         endif
      endif 

c     Compute f0 and g0.

      task = 'FG_START' 
c          return to the driver to calculate f and g; reenter at 111.
      goto 1000
 111  continue
      nfgv = 1
 
c     Compute the infinity norm of the (-) projected gradient.
 
      call projgr(n,l,u,nbd,x,g,sbgnrm)
  
      if (iprint .ge. 1) then
         write (6,1002) iter,f,sbgnrm
         write (itfile,1003) iter,nfgv,sbgnrm,f
      endif
      if (sbgnrm .le. pgtol) then
c                                terminate the algorithm.
         task = 'CONVERGENCE: NORM_OF_PROJECTED_GRADIENT_<=_PGTOL'
         goto 999
      endif 
 
c ----------------- the beginning of the loop --------------------------
 
 222  continue
      if (iprint .ge. 99) write (6,1001) iter + 1
      iword = -1
c
      if (.not. cnstnd .and. col .gt. 0) then 
c                                            skip the search for GCP.
         call dcopy(n,x,1,z,1)
         wrk = updatd
         nseg = 0
         goto 333
      endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the Generalized Cauchy Point (GCP).
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call timer(cpu1) 
      call cauchy(n,x,l,u,nbd,g,indx2,iwhere,t,d,z,
     +            m,wy,ws,sy,wt,theta,col,head,
     +            wa(1),wa(2*m+1),wa(4*m+1),wa(6*m+1),nseg,
     +            iprint, sbgnrm, info, epsmch)
      if (info .ne. 0) then 
c         singular triangular system detected; refresh the lbfgs memory.
         if(iprint .ge. 1) write (6, 1005)
         info   = 0
         col    = 0
         head   = 1
         theta  = one
         iupdat = 0
         updatd = .false.
         call timer(cpu2) 
         cachyt = cachyt + cpu2 - cpu1
         goto 222
      endif
      call timer(cpu2) 
      cachyt = cachyt + cpu2 - cpu1
      nintol = nintol + nseg

c     Count the entering and leaving variables for iter > 0; 
c     find the index set of free and active variables at the GCP.

      call freev(n,nfree,index,nenter,ileave,indx2,
     +           iwhere,wrk,updatd,cnstnd,iprint,iter)
      nact = n - nfree

 333  continue
 
c     If there are no free variables or B=theta*I, then
c                                        skip the subspace minimization.
 
      if (nfree .eq. 0 .or. col .eq. 0) goto 555
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Subspace minimization.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call timer(cpu1) 

c     Form  the LEL^T factorization of the indefinite
c       matrix    K = [-D -Y'ZZ'Y/theta     L_a'-R_z'  ]
c                     [L_a -R_z           theta*S'AA'S ]
c       where     E = [-I  0]
c                     [ 0  I]

      if (wrk) call formk(n,nfree,index,nenter,ileave,indx2,iupdat,
     +                 updatd,wn,snd,m,ws,wy,sy,theta,col,head,info)
      if (info .ne. 0) then
c          nonpositive definiteness in Cholesky factorization;
c          refresh the lbfgs memory and restart the iteration.
         if(iprint .ge. 1) write (6, 1006)
         info   = 0
         col    = 0
         head   = 1
         theta  = one
         iupdat = 0
         updatd = .false.
         call timer(cpu2) 
         sbtime = sbtime + cpu2 - cpu1 
         goto 222
      endif 

c        compute r=-Z'B(xcp-xk)-Z'g (using wa(2m+1)=W'(xcp-x)
c                                                   from 'cauchy').
      call cmprlb(n,m,x,g,ws,wy,sy,wt,z,r,wa,index,
     +           theta,col,head,nfree,cnstnd,info)
      if (info .ne. 0) goto 444

c-jlm-jn   call the direct method. 

      call subsm( n, m, nfree, index, l, u, nbd, z, r, xp, ws, wy,
     +           theta, x, g, col, head, iword, wa, wn, iprint, info)
 444  continue
      if (info .ne. 0) then 
c          singular triangular system detected;
c          refresh the lbfgs memory and restart the iteration.
         if(iprint .ge. 1) write (6, 1005)
         info   = 0
         col    = 0
         head   = 1
         theta  = one
         iupdat = 0
         updatd = .false.
         call timer(cpu2) 
         sbtime = sbtime + cpu2 - cpu1 
         goto 222
      endif
 
      call timer(cpu2) 
      sbtime = sbtime + cpu2 - cpu1 
 555  continue
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Line search and optimality tests.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
c     Generate the search direction d:=z-x.

      do 40 i = 1, n
         d(i) = z(i) - x(i)
  40  continue
      call timer(cpu1) 
 666  continue
      call lnsrlb(n,l,u,nbd,x,f,fold,gd,gdold,g,d,r,t,z,stp,dnorm,
     +            dtd,xstep,stpmx,iter,ifun,iback,nfgv,info,task,
     +            boxed,cnstnd,csave,isave(22),dsave(17))
      if (info .ne. 0 .or. iback .ge. 20) then
c          restore the previous iterate.
         call dcopy(n,t,1,x,1)
         call dcopy(n,r,1,g,1)
         f = fold
         if (col .eq. 0) then
c             abnormal termination.
            if (info .eq. 0) then
               info = -9
c                restore the actual number of f and g evaluations etc.
               nfgv = nfgv - 1
               ifun = ifun - 1
               iback = iback - 1
            endif
            task = 'ABNORMAL_TERMINATION_IN_LNSRCH'
            iter = iter + 1
            goto 999
         else
c             refresh the lbfgs memory and restart the iteration.
            if(iprint .ge. 1) write (6, 1008)
            if (info .eq. 0) nfgv = nfgv - 1
            info   = 0
            col    = 0
            head   = 1
            theta  = one
            iupdat = 0
            updatd = .false.
            task   = 'RESTART_FROM_LNSRCH'
            call timer(cpu2)
            lnscht = lnscht + cpu2 - cpu1
            goto 222
         endif
      else if (task(1:5) .eq. 'FG_LN') then
c          return to the driver for calculating f and g; reenter at 666.
         goto 1000
      else 
c          calculate and print out the quantities related to the new X.
         call timer(cpu2) 
         lnscht = lnscht + cpu2 - cpu1
         iter = iter + 1
 
c        Compute the infinity norm of the projected (-)gradient.
 
         call projgr(n,l,u,nbd,x,g,sbgnrm)
 
c        Print iteration information.

         call prn2lb(n,x,f,g,iprint,itfile,iter,nfgv,nact,
     +               sbgnrm,nseg,word,iword,iback,stp,xstep)
         goto 1000
      endif
 777  continue

c     Test for termination.

      if (sbgnrm .le. pgtol) then
c                                terminate the algorithm.
         task = 'CONVERGENCE: NORM_OF_PROJECTED_GRADIENT_<=_PGTOL'
         goto 999
      endif 

      ddum = max(abs(fold), abs(f), one)
      if ((fold - f) .le. tol*ddum) then
c                                        terminate the algorithm.
         task = 'CONVERGENCE: REL_REDUCTION_OF_F_<=_FACTR*EPSMCH'
         if (iback .ge. 10) info = -5
c           i.e., to issue a warning if iback>10 in the line search.
         goto 999
      endif 

c     Compute d=newx-oldx, r=newg-oldg, rr=y'y and dr=y's.
 
      do 42 i = 1, n
         r(i) = g(i) - r(i)
  42  continue
      rr = ddot(n,r,1,r,1)
      if (stp .eq. one) then  
         dr = gd - gdold
         ddum = -gdold
      else
         dr = (gd - gdold)*stp
         call dscal(n,stp,d,1)
         ddum = -gdold*stp
      endif
 
      if (dr .le. epsmch*ddum) then
c                            skip the L-BFGS update.
         nskip = nskip + 1
         updatd = .false.
         if (iprint .ge. 1) write (6,1004) dr, ddum
         goto 888
      endif 
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Update the L-BFGS matrix.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
      updatd = .true.
      iupdat = iupdat + 1

c     Update matrices WS and WY and form the middle matrix in B.

      call matupd(n,m,ws,wy,sy,ss,d,r,itail,
     +            iupdat,col,head,theta,rr,dr,stp,dtd)

c     Form the upper half of the pds T = theta*SS + L*D^(-1)*L';
c        Store T in the upper triangular of the array wt;
c        Cholesky factorize T to J*J' with
c           J' stored in the upper triangular of wt.

      call formt(m,wt,sy,ss,col,theta,info)
 
      if (info .ne. 0) then 
c          nonpositive definiteness in Cholesky factorization;
c          refresh the lbfgs memory and restart the iteration.
         if(iprint .ge. 1) write (6, 1007)
         info = 0
         col = 0
         head = 1
         theta = one
         iupdat = 0
         updatd = .false.
         goto 222
      endif

c     Now the inverse of the middle matrix in B is

c       [  D^(1/2)      O ] [ -D^(1/2)  D^(-1/2)*L' ]
c       [ -L*D^(-1/2)   J ] [  0        J'          ]

 888  continue
 
c -------------------- the end of the loop -----------------------------
 
      goto 222
 999  continue
      call timer(time2)
      time = time2 - time1
      call prn3lb(n,x,f,task,iprint,info,itfile,
     +            iter,nfgv,nintol,nskip,nact,sbgnrm,
     +            time,nseg,word,iback,stp,xstep,k,
     +            cachyt,sbtime,lnscht)
 1000 continue

c     Save local variables.

      lsave(1)  = prjctd
      lsave(2)  = cnstnd
      lsave(3)  = boxed
      lsave(4)  = updatd

      isave(1)  = nintol 
      isave(3)  = itfile 
      isave(4)  = iback 
      isave(5)  = nskip 
      isave(6)  = head 
      isave(7)  = col 
      isave(8)  = itail 
      isave(9)  = iter 
      isave(10) = iupdat 
      isave(12) = nseg
      isave(13) = nfgv 
      isave(14) = info 
      isave(15) = ifun 
      isave(16) = iword 
      isave(17) = nfree 
      isave(18) = nact 
      isave(19) = ileave 
      isave(20) = nenter 

      dsave(1)  = theta 
      dsave(2)  = fold 
      dsave(3)  = tol 
      dsave(4)  = dnorm 
      dsave(5)  = epsmch 
      dsave(6)  = cpu1 
      dsave(7)  = cachyt 
      dsave(8)  = sbtime 
      dsave(9)  = lnscht 
      dsave(10) = time1 
      dsave(11) = gd 
      dsave(12) = stpmx 
      dsave(13) = sbgnrm
      dsave(14) = stp
      dsave(15) = gdold
      dsave(16) = dtd  

 1001 format (//,'ITERATION ',i5)
 1002 format
     +  (/,'At iterate',i5,4x,'f= ',1p,d12.5,4x,'|proj g|= ',1p,d12.5)
 1003 format (2(1x,i4),5x,'-',5x,'-',3x,'-',5x,'-',5x,'-',8x,'-',3x,
     +        1p,2(1x,d10.3))
 1004 format ('  ys=',1p,e10.3,'  -gs=',1p,e10.3,' BFGS update SKIPPED')
 1005 format (/, 
     +' Singular triangular system detected;',/,
     +'   refresh the lbfgs memory and restart the iteration.')
 1006 format (/, 
     +' Nonpositive definiteness in Cholesky factorization in formk;',/,
     +'   refresh the lbfgs memory and restart the iteration.')
 1007 format (/, 
     +' Nonpositive definiteness in Cholesky factorization in formt;',/,
     +'   refresh the lbfgs memory and restart the iteration.')
 1008 format (/, 
     +' Bad direction in the line search;',/,
     +'   refresh the lbfgs memory and restart the iteration.')

      return   

      end
 
c======================= The end of mainlb =============================

      subroutine active(n, l, u, nbd, x, iwhere, iprint,
     +                  prjctd, cnstnd, boxed)

      logical          prjctd, cnstnd, boxed
      integer          n, iprint, nbd(n), iwhere(n)
      double precision x(n), l(n), u(n)

c     ************
c
c     Subroutine active
c
c     This subroutine initializes iwhere and projects the initial x to
c       the feasible set if necessary.
c
c     iwhere is an integer array of dimension n.
c       On entry iwhere is unspecified.
c       On exit iwhere(i)=-1  if x(i) has no bounds
c                         3   if l(i)=u(i)
c                         0   otherwise.
c       In cauchy, iwhere is given finer gradations.
c
c
c                           *  *  *
c
c     NEOS, November 1994. (Latest revision June 1996.)
c     Optimization Technology Center.
c     Argonne National Laboratory and Northwestern University.
c     Written by
c                        Ciyou Zhu
c     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
c
c
c     ************

      integer          nbdd,i
      double precision zero
      parameter        (zero=0.0d0)

c     Initialize nbdd, prjctd, cnstnd and boxed.

      nbdd = 0
      prjctd = .false.
      cnstnd = .false.
      boxed = .true.

c     Project the initial x to the easible set if necessary.

      do 10 i = 1, n
         if (nbd(i) .gt. 0) then
            if (nbd(i) .le. 2 .and. x(i) .le. l(i)) then
               if (x(i) .lt. l(i)) then
                  prjctd = .true.
                  x(i) = l(i)
               endif
               nbdd = nbdd + 1
            else if (nbd(i) .ge. 2 .and. x(i) .ge. u(i)) then
               if (x(i) .gt. u(i)) then
                  prjctd = .true.
                  x(i) = u(i)
               endif
               nbdd = nbdd + 1
            endif
         endif
  10  continue

c     Initialize iwhere and assign values to cnstnd and boxed.

      do 20 i = 1, n
         if (nbd(i) .ne. 2) boxed = .false.
         if (nbd(i) .eq. 0) then
c                                this variable is always free
            iwhere(i) = -1

c           otherwise set x(i)=mid(x(i), u(i), l(i)).
         else
            cnstnd = .true.
            if (nbd(i) .eq. 2 .and. u(i) - l(i) .le. zero) then
c                   this variable is always fixed
               iwhere(i) = 3
            else 
               iwhere(i) = 0
            endif
         endif
  20  continue

      if (iprint .ge. 0) then
         if (prjctd) write (6,*)
     +   'The initial X is infeasible.  Restart with its projection.'
         if (.not. cnstnd)
     +      write (6,*) 'This problem is unconstrained.'
      endif

      if (iprint .gt. 0) write (6,1001) nbdd

 1001 format (/,'At X0 ',i9,' variables are exactly at the bounds') 

      return

      end

c======================= The end of active =============================
 
      subroutine bmv(m, sy, wt, col, v, p, info)

      integer m, col, info
      double precision sy(m, m), wt(m, m), v(2*col), p(2*col)

c     ************
c
c     Subroutine bmv
c
c     This subroutine computes the product of the 2m x 2m middle matrix 
c       in the compact L-BFGS formula of B and a 2m vector v;  
c       it returns the product in p.
c       
c     m is an integer variable.
c       On entry m is the maximum number of variable metric corrections
c         used to define the limited memory matrix.
c       On exit m is unchanged.
c
c     sy is a double precision array of dimension m x m.
c       On entry sy specifies the matrix S'Y.
c       On exit sy is unchanged.
c
c     wt is a double precision array of dimension m x m.
c       On entry wt specifies the upper triangular matrix J' which is 
c         the Cholesky factor of (thetaS'S+LD^(-1)L').
c       On exit wt is unchanged.
c
c     col is an integer variable.
c       On entry col specifies the number of s-vectors (or y-vectors)
c         stored in the compact L-BFGS formula.
c       On exit col is unchanged.
c
c     v is a double precision array of dimension 2col.
c       On entry v specifies vector v.
c       On exit v is unchanged.
c
c     p is a double precision array of dimension 2col.
c       On entry p is unspecified.
c       On exit p is the product Mv.
c
c     info is an integer variable.
c       On entry info is unspecified.
c       On exit info = 0       for normal return,
c                    = nonzero for abnormal return when the system
c                                to be solved by dtrsl is singular.
c
c     Subprograms called:
c
c       Linpack ... dtrsl.
c
c
c                           *  *  *
c
c     NEOS, November 1994. (Latest revision June 1996.)
c     Optimization Technology Center.
c     Argonne National Laboratory and Northwestern University.
c     Written by
c                        Ciyou Zhu
c     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
c
c
c     ************
 
      integer          i,k,i2
      double precision sum
 
      if (col .eq. 0) return
 
c     PART I: solve [  D^(1/2)      O ] [ p1 ] = [ v1 ]
c                   [ -L*D^(-1/2)   J ] [ p2 ]   [ v2 ].

c       solve Jp2=v2+LD^(-1)v1.
      p(col + 1) = v(col + 1)
      do 20 i = 2, col
         i2 = col + i
         sum = 0.0d0
         do 10 k = 1, i - 1
            sum = sum + sy(i,k)*v(k)/sy(k,k)
  10     continue
         p(i2) = v(i2) + sum
  20  continue  
c     Solve the triangular system
      call dtrsl(wt,m,col,p(col+1),11,info)
      if (info .ne. 0) return
 
c       solve D^(1/2)p1=v1.
      do 30 i = 1, col
         p(i) = v(i)/sqrt(sy(i,i))
  30  continue 
 
c     PART II: solve [ -D^(1/2)   D^(-1/2)*L'  ] [ p1 ] = [ p1 ]
c                    [  0         J'           ] [ p2 ]   [ p2 ]. 
 
c       solve J^Tp2=p2. 
      call dtrsl(wt,m,col,p(col+1),01,info)
      if (info .ne. 0) return
 
c       compute p1=-D^(-1/2)(p1-D^(-1/2)L'p2)
c                 =-D^(-1/2)p1+D^(-1)L'p2.  
      do 40 i = 1, col
         p(i) = -p(i)/sqrt(sy(i,i))
  40  continue
      do 60 i = 1, col
         sum = 0.d0
         do 50 k = i + 1, col
            sum = sum + sy(k,i)*p(col+k)/sy(i,i)
  50     continue
         p(i) = p(i) + sum
  60  continue

      return

      end

c======================== The end of bmv ===============================

      subroutine cauchy(n, x, l, u, nbd, g, iorder, iwhere, t, d, xcp, 
     +                  m, wy, ws, sy, wt, theta, col, head, p, c, wbp, 
     +                  v, nseg, iprint, sbgnrm, info, epsmch)
      implicit none
      integer          n, m, head, col, nseg, iprint, info, 
     +                 nbd(n), iorder(n), iwhere(n)
      double precision theta, epsmch,
     +                 x(n), l(n), u(n), g(n), t(n), d(n), xcp(n),
     +                 wy(n, col), ws(n, col), sy(m, m),
     +                 wt(m, m), p(2*m), c(2*m), wbp(2*m), v(2*m)

c     ************
c
c     Subroutine cauchy
c
c     For given x, l, u, g (with sbgnrm > 0), and a limited memory
c       BFGS matrix B defined in terms of matrices WY, WS, WT, and
c       scalars head, col, and theta, this subroutine computes the
c       generalized Cauchy point (GCP), defined as the first local
c       minimizer of the quadratic
c
c                  Q(x + s) = g's + 1/2 s'Bs
c
c       along the projected gradient direction P(x-tg,l,u).
c       The routine returns the GCP in xcp. 
c       
c     n is an integer variable.
c       On entry n is the dimension of the problem.
c       On exit n is unchanged.
c
c     x is a double precision array of dimension n.
c       On entry x is the starting point for the GCP computation.
c       On exit x is unchanged.
c
c     l is a double precision array of dimension n.
c       On entry l is the lower bound of x.
c       On exit l is unchanged.
c
c     u is a double precision array of dimension n.
c       On entry u is the upper bound of x.
c       On exit u is unchanged.
c
c     nbd is an integer array of dimension n.
c       On entry nbd represents the type of bounds imposed on the
c         variables, and must be specified as follows:
c         nbd(i)=0 if x(i) is unbounded,
c                1 if x(i) has only a lower bound,
c                2 if x(i) has both lower and upper bounds, and
c                3 if x(i) has only an upper bound. 
c       On exit nbd is unchanged.
c
c     g is a double precision array of dimension n.
c       On entry g is the gradient of f(x).  g must be a nonzero vector.
c       On exit g is unchanged.
c
c     iorder is an integer working array of dimension n.
c       iorder will be used to store the breakpoints in the piecewise
c       linear path and free variables encountered. On exit,
c         iorder(1),...,iorder(nleft) are indices of breakpoints
c                                which have not been encountered; 
c         iorder(nleft+1),...,iorder(nbreak) are indices of
c                                     encountered breakpoints; and
c         iorder(nfree),...,iorder(n) are indices of variables which
c                 have no bound constraits along the search direction.
c
c     iwhere is an integer array of dimension n.
c       On entry iwhere indicates only the permanently fixed (iwhere=3)
c       or free (iwhere= -1) components of x.
c       On exit iwhere records the status of the current x variables.
c       iwhere(i)=-3  if x(i) is free and has bounds, but is not moved
c                 0   if x(i) is free and has bounds, and is moved
c                 1   if x(i) is fixed at l(i), and l(i) .ne. u(i)
c                 2   if x(i) is fixed at u(i), and u(i) .ne. l(i)
c                 3   if x(i) is always fixed, i.e.,  u(i)=x(i)=l(i)
c                 -1  if x(i) is always free, i.e., it has no bounds.
c
c     t is a double precision working array of dimension n. 
c       t will be used to store the break points.
c
c     d is a double precision array of dimension n used to store
c       the Cauchy direction P(x-tg)-x.
c
c     xcp is a double precision array of dimension n used to return the
c       GCP on exit.
c
c     m is an integer variable.
c       On entry m is the maximum number of variable metric corrections 
c         used to define the limited memory matrix.
c       On exit m is unchanged.
c
c     ws, wy, sy, and wt are double precision arrays.
c       On entry they store information that defines the
c                             limited memory BFGS matrix:
c         ws(n,m) stores S, a set of s-vectors;
c         wy(n,m) stores Y, a set of y-vectors;
c         sy(m,m) stores S'Y;
c         wt(m,m) stores the
c                 Cholesky factorization of (theta*S'S+LD^(-1)L').
c       On exit these arrays are unchanged.
c
c     theta is a double precision variable.
c       On entry theta is the scaling factor specifying B_0 = theta I.
c       On exit theta is unchanged.
c
c     col is an integer variable.
c       On entry col is the actual number of variable metric
c         corrections stored so far.
c       On exit col is unchanged.
c
c     head is an integer variable.
c       On entry head is the location of the first s-vector (or y-vector)
c         in S (or Y).
c       On exit col is unchanged.
c
c     p is a double precision working array of dimension 2m.
c       p will be used to store the vector p = W^(T)d.
c
c     c is a double precision working array of dimension 2m.
c       c will be used to store the vector c = W^(T)(xcp-x).
c
c     wbp is a double precision working array of dimension 2m.
c       wbp will be used to store the row of W corresponding
c         to a breakpoint.
c
c     v is a double precision working array of dimension 2m.
c
c     nseg is an integer variable.
c       On exit nseg records the number of quadratic segments explored
c         in searching for the GCP.
c
c     sg and yg are double precision arrays of dimension m.
c       On entry sg  and yg store S'g and Y'g correspondingly.
c       On exit they are unchanged. 
c 
c     iprint is an INTEGER variable that must be set by the user.
c       It controls the frequency and type of output generated:
c        iprint<0    no output is generated;
c        iprint=0    print only one line at the last iteration;
c        0<iprint<99 print also f and |proj g| every iprint iterations;
c        iprint=99   print details of every iteration except n-vectors;
c        iprint=100  print also the changes of active set and final x;
c        iprint>100  print details of every iteration including x and g;
c       When iprint > 0, the file iterate.dat will be created to
c                        summarize the iteration.
c
c     sbgnrm is a double precision variable.
c       On entry sbgnrm is the norm of the projected gradient at x.
c       On exit sbgnrm is unchanged.
c
c     info is an integer variable.
c       On entry info is 0.
c       On exit info = 0       for normal return,
c                    = nonzero for abnormal return when the the system
c                              used in routine bmv is singular.
c
c     Subprograms called:
c 
c       L-BFGS-B Library ... hpsolb, bmv.
c
c       Linpack ... dscal dcopy, daxpy.
c
c
c     References:
c
c       [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited
c       memory algorithm for bound constrained optimization'',
c       SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.
c
c       [2] C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, ``L-BFGS-B: FORTRAN
c       Subroutines for Large Scale Bound Constrained Optimization''
c       Tech. Report, NAM-11, EECS Department, Northwestern University,
c       1994.
c
c       (Postscript files of these papers are available via anonymous
c        ftp to eecs.nwu.edu in the directory pub/lbfgs/lbfgs_bcm.)
c
c                           *  *  *
c
c     NEOS, November 1994. (Latest revision June 1996.)
c     Optimization Technology Center.
c     Argonne National Laboratory and Northwestern University.
c     Written by
c                        Ciyou Zhu
c     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
c
c
c     ************

      logical          xlower,xupper,bnded
      integer          i,j,col2,nfree,nbreak,pointr,
     +                 ibp,nleft,ibkmin,iter
      double precision f1,f2,dt,dtm,tsum,dibp,zibp,dibp2,bkmin,
     +                 tu,tl,wmc,wmp,wmw,ddot,tj,tj0,neggi,sbgnrm,
     +                 f2_org
      double precision one,zero
      parameter        (one=1.0d0,zero=0.0d0)
 
c     Check the status of the variables, reset iwhere(i) if necessary;
c       compute the Cauchy direction d and the breakpoints t; initialize
c       the derivative f1 and the vector p = W'd (for theta = 1).
 
      if (sbgnrm .le. zero) then
         if (iprint .ge. 0) write (6,*) 'Subgnorm = 0.  GCP = X.'
         call dcopy(n,x,1,xcp,1)
         return
      endif 
      bnded = .true.
      nfree = n + 1
      nbreak = 0
      ibkmin = 0
      bkmin = zero
      col2 = 2*col
      f1 = zero
      if (iprint .ge. 99) write (6,3010)

c     We set p to zero and build it up as we determine d.

      do 20 i = 1, col2
         p(i) = zero
  20  continue 

c     In the following loop we determine for each variable its bound
c        status and its breakpoint, and update p accordingly.
c        Smallest breakpoint is identified.

      do 50 i = 1, n 
         neggi = -g(i)      
         if (iwhere(i) .ne. 3 .and. iwhere(i) .ne. -1) then
c             if x(i) is not a constant and has bounds,
c             compute the difference between x(i) and its bounds.
            if (nbd(i) .le. 2) tl = x(i) - l(i)
            if (nbd(i) .ge. 2) tu = u(i) - x(i)

c           If a variable is close enough to a bound
c             we treat it as at bound.
            xlower = nbd(i) .le. 2 .and. tl .le. zero
            xupper = nbd(i) .ge. 2 .and. tu .le. zero

c              reset iwhere(i).
            iwhere(i) = 0
            if (xlower) then
               if (neggi .le. zero) iwhere(i) = 1
            else if (xupper) then
               if (neggi .ge. zero) iwhere(i) = 2
            else
               if (abs(neggi) .le. zero) iwhere(i) = -3
            endif
         endif 
         pointr = head
         if (iwhere(i) .ne. 0 .and. iwhere(i) .ne. -1) then
            d(i) = zero
         else
            d(i) = neggi
            f1 = f1 - neggi*neggi
c             calculate p := p - W'e_i* (g_i).
            do 40 j = 1, col
               p(j) = p(j) +  wy(i,pointr)* neggi
               p(col + j) = p(col + j) + ws(i,pointr)*neggi
               pointr = mod(pointr,m) + 1
  40        continue 
            if (nbd(i) .le. 2 .and. nbd(i) .ne. 0
     +                        .and. neggi .lt. zero) then
c                                 x(i) + d(i) is bounded; compute t(i).
               nbreak = nbreak + 1
               iorder(nbreak) = i
               t(nbreak) = tl/(-neggi)
               if (nbreak .eq. 1 .or. t(nbreak) .lt. bkmin) then
                  bkmin = t(nbreak)
                  ibkmin = nbreak
               endif
            else if (nbd(i) .ge. 2 .and. neggi .gt. zero) then
c                                 x(i) + d(i) is bounded; compute t(i).
               nbreak = nbreak + 1
               iorder(nbreak) = i
               t(nbreak) = tu/neggi
               if (nbreak .eq. 1 .or. t(nbreak) .lt. bkmin) then
                  bkmin = t(nbreak)
                  ibkmin = nbreak
               endif
            else
c                x(i) + d(i) is not bounded.
               nfree = nfree - 1
               iorder(nfree) = i
               if (abs(neggi) .gt. zero) bnded = .false.
            endif
         endif
  50  continue 
 
c     The indices of the nonzero components of d are now stored
c       in iorder(1),...,iorder(nbreak) and iorder(nfree),...,iorder(n).
c       The smallest of the nbreak breakpoints is in t(ibkmin)=bkmin.
 
      if (theta .ne. one) then
c                   complete the initialization of p for theta not= one.
         call dscal(col,theta,p(col+1),1)
      endif
 
c     Initialize GCP xcp = x.

      call dcopy(n,x,1,xcp,1)

      if (nbreak .eq. 0 .and. nfree .eq. n + 1) then
c                  is a zero vector, return with the initial xcp as GCP.
         if (iprint .gt. 100) write (6,1010) (xcp(i), i = 1, n)
         return
      endif    
 
c     Initialize c = W'(xcp - x) = 0.
  
      do 60 j = 1, col2
         c(j) = zero
  60  continue 
 
c     Initialize derivative f2.
 
      f2 =  -theta*f1 
      f2_org  =  f2
      if (col .gt. 0) then
         call bmv(m,sy,wt,col,p,v,info)
         if (info .ne. 0) return
         f2 = f2 - ddot(col2,v,1,p,1)
      endif
      dtm = -f1/f2
      tsum = zero
      nseg = 1
      if (iprint .ge. 99) 
     +   write (6,*) 'There are ',nbreak,'  breakpoints '
 
c     If there are no breakpoints, locate the GCP and return. 
 
      if (nbreak .eq. 0) goto 888
             
      nleft = nbreak
      iter = 1
 
 
      tj = zero
 
c------------------- the beginning of the loop -------------------------
 
 777  continue
 
c     Find the next smallest breakpoint;
c       compute dt = t(nleft) - t(nleft + 1).
 
      tj0 = tj
      if (iter .eq. 1) then
c         Since we already have the smallest breakpoint we need not do
c         heapsort yet. Often only one breakpoint is used and the
c         cost of heapsort is avoided.
         tj = bkmin
         ibp = iorder(ibkmin)
      else
         if (iter .eq. 2) then
c             Replace the already used smallest breakpoint with the
c             breakpoint numbered nbreak > nlast, before heapsort call.
            if (ibkmin .ne. nbreak) then
               t(ibkmin) = t(nbreak)
               iorder(ibkmin) = iorder(nbreak)
            endif 
c        Update heap structure of breakpoints
c           (if iter=2, initialize heap).
         endif
         call hpsolb(nleft,t,iorder,iter-2)
         tj = t(nleft)
         ibp = iorder(nleft)  
      endif 
         
      dt = tj - tj0
 
      if (dt .ne. zero .and. iprint .ge. 100) then
         write (6,4011) nseg,f1,f2
         write (6,5010) dt
         write (6,6010) dtm
      endif          
 
c     If a minimizer is within this interval, locate the GCP and return. 
 
      if (dtm .lt. dt) goto 888
 
c     Otherwise fix one variable and
c       reset the corresponding component of d to zero.
    
      tsum = tsum + dt
      nleft = nleft - 1
      iter = iter + 1
      dibp = d(ibp)
      d(ibp) = zero
      if (dibp .gt. zero) then
         zibp = u(ibp) - x(ibp)
         xcp(ibp) = u(ibp)
         iwhere(ibp) = 2
      else
         zibp = l(ibp) - x(ibp)
         xcp(ibp) = l(ibp)
         iwhere(ibp) = 1
      endif
      if (iprint .ge. 100) write (6,*) 'Variable  ',ibp,'  is fixed.'
      if (nleft .eq. 0 .and. nbreak .eq. n) then
c                                             all n variables are fixed,
c                                                return with xcp as GCP.
         dtm = dt
         goto 999
      endif
 
c     Update the derivative information.
 
      nseg = nseg + 1
      dibp2 = dibp**2
 
c     Update f1 and f2.
 
c        temporarily set f1 and f2 for col=0.
      f1 = f1 + dt*f2 + dibp2 - theta*dibp*zibp
      f2 = f2 - theta*dibp2

      if (col .gt. 0) then
c                          update c = c + dt*p.
         call daxpy(col2,dt,p,1,c,1)
 
c           choose wbp,
c           the row of W corresponding to the breakpoint encountered.
         pointr = head
         do 70 j = 1,col
            wbp(j) = wy(ibp,pointr)
            wbp(col + j) = theta*ws(ibp,pointr)
            pointr = mod(pointr,m) + 1
  70     continue 
 
c           compute (wbp)Mc, (wbp)Mp, and (wbp)M(wbp)'.
         call bmv(m,sy,wt,col,wbp,v,info)
         if (info .ne. 0) return
         wmc = ddot(col2,c,1,v,1)
         wmp = ddot(col2,p,1,v,1) 
         wmw = ddot(col2,wbp,1,v,1)
 
c           update p = p - dibp*wbp. 
         call daxpy(col2,-dibp,wbp,1,p,1)
 
c           complete updating f1 and f2 while col > 0.
         f1 = f1 + dibp*wmc
         f2 = f2 + 2.0d0*dibp*wmp - dibp2*wmw
      endif

      f2 = max(epsmch*f2_org,f2)
      if (nleft .gt. 0) then
         dtm = -f1/f2
         goto 777
c                 to repeat the loop for unsearched intervals. 
      else if(bnded) then
         f1 = zero
         f2 = zero
         dtm = zero
      else
         dtm = -f1/f2
      endif 

c------------------- the end of the loop -------------------------------
 
 888  continue
      if (iprint .ge. 99) then
         write (6,*)
         write (6,*) 'GCP found in this segment'
         write (6,4010) nseg,f1,f2
         write (6,6010) dtm
      endif 
      if (dtm .le. zero) dtm = zero
      tsum = tsum + dtm
 
c     Move free variables (i.e., the ones w/o breakpoints) and 
c       the variables whose breakpoints haven't been reached.
 
      call daxpy(n,tsum,d,1,xcp,1)
 
 999  continue
 
c     Update c = c + dtm*p = W'(x^c - x) 
c       which will be used in computing r = Z'(B(x^c - x) + g).
 
      if (col .gt. 0) call daxpy(col2,dtm,p,1,c,1)
      if (iprint .gt. 100) write (6,1010) (xcp(i),i = 1,n)
      if (iprint .ge. 99) write (6,2010)

 1010 format ('Cauchy X =  ',/,(4x,1p,6(1x,d11.4)))
 2010 format (/,'---------------- exit CAUCHY----------------------',/)
 3010 format (/,'---------------- CAUCHY entered-------------------')
 4010 format ('Piece    ',i3,' --f1, f2 at start point ',1p,2(1x,d11.4))
 4011 format (/,'Piece    ',i3,' --f1, f2 at start point ',
     +        1p,2(1x,d11.4))
 5010 format ('Distance to the next break point =  ',1p,d11.4)
 6010 format ('Distance to the stationary point =  ',1p,d11.4) 
 
      return
 
      end

c====================== The end of cauchy ==============================

      subroutine cmprlb(n, m, x, g, ws, wy, sy, wt, z, r, wa, index, 
     +                 theta, col, head, nfree, cnstnd, info)
 
      logical          cnstnd
      integer          n, m, col, head, nfree, info, index(n)
      double precision theta, 
     +                 x(n), g(n), z(n), r(n), wa(4*m), 
     +                 ws(n, m), wy(n, m), sy(m, m), wt(m, m)

c     ************
c
c     Subroutine cmprlb 
c
c       This subroutine computes r=-Z'B(xcp-xk)-Z'g by using 
c         wa(2m+1)=W'(xcp-x) from subroutine cauchy.
c
c     Subprograms called:
c
c       L-BFGS-B Library ... bmv.
c
c
c                           *  *  *
c
c     NEOS, November 1994. (Latest revision June 1996.)
c     Optimization Technology Center.
c     Argonne National Laboratory and Northwestern University.
c     Written by
c                        Ciyou Zhu
c     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
c
c
c     ************
 
      integer          i,j,k,pointr
      double precision a1,a2

      if (.not. cnstnd .and. col .gt. 0) then 
         do 26 i = 1, n
            r(i) = -g(i)
  26     continue
      else
         do 30 i = 1, nfree
            k = index(i)
            r(i) = -theta*(z(k) - x(k)) - g(k)
  30     continue
         call bmv(m,sy,wt,col,wa(2*m+1),wa(1),info)
         if (info .ne. 0) then
            info = -8
            return
         endif
         pointr = head 
         do 34 j = 1, col
            a1 = wa(j)
            a2 = theta*wa(col + j)
            do 32 i = 1, nfree
               k = index(i)
               r(i) = r(i) + wy(k,pointr)*a1 + ws(k,pointr)*a2
  32        continue
            pointr = mod(pointr,m) + 1
  34     continue
      endif

      return

      end

c======================= The end of cmprlb =============================

      subroutine errclb(n, m, factr, l, u, nbd, task, info, k)
 
      character*60     task
      integer          n, m, info, k, nbd(n)
      double precision factr, l(n), u(n)

c     ************
c
c     Subroutine errclb
c
c     This subroutine checks the validity of the input data.
c
c
c                           *  *  *
c
c     NEOS, November 1994. (Latest revision June 1996.)
c     Optimization Technology Center.
c     Argonne National Laboratory and Northwestern University.
c     Written by
c                        Ciyou Zhu
c     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
c
c
c     ************

      integer          i
      double precision one,zero
      parameter        (one=1.0d0,zero=0.0d0)

c     Check the input arguments for errors.

      if (n .le. 0) task = 'ERROR: N .LE. 0'
      if (m .le. 0) task = 'ERROR: M .LE. 0'
      if (factr .lt. zero) task = 'ERROR: FACTR .LT. 0'

c     Check the validity of the arrays nbd(i), u(i), and l(i).

      do 10 i = 1, n
         if (nbd(i) .lt. 0 .or. nbd(i) .gt. 3) then
c                                                   return
            task = 'ERROR: INVALID NBD'
            info = -6
            k = i
         endif
         if (nbd(i) .eq. 2) then
            if (l(i) .gt. u(i)) then
c                                    return
               task = 'ERROR: NO FEASIBLE SOLUTION'
               info = -7
               k = i
            endif
         endif
  10  continue

      return

      end

c======================= The end of errclb =============================
 
      subroutine formk(n, nsub, ind, nenter, ileave, indx2, iupdat, 
     +                 updatd, wn, wn1, m, ws, wy, sy, theta, col,
     +                 head, info)

      integer          n, nsub, m, col, head, nenter, ileave, iupdat,
     +                 info, ind(n), indx2(n)
      double precision theta, wn(2*m, 2*m), wn1(2*m, 2*m),
     +                 ws(n, m), wy(n, m), sy(m, m)
      logical          updatd

c     ************
c
c     Subroutine formk 
c
c     This subroutine forms  the LEL^T factorization of the indefinite
c
c       matrix    K = [-D -Y'ZZ'Y/theta     L_a'-R_z'  ]
c                     [L_a -R_z           theta*S'AA'S ]
c                                                    where E = [-I  0]
c                                                              [ 0  I]
c     The matrix K can be shown to be equal to the matrix M^[-1]N
c       occurring in section 5.1 of [1], as well as to the matrix
c       Mbar^[-1] Nbar in section 5.3.
c
c     n is an integer variable.
c       On entry n is the dimension of the problem.
c       On exit n is unchanged.
c
c     nsub is an integer variable
c       On entry nsub is the number of subspace variables in free set.
c       On exit nsub is not changed.
c
c     ind is an integer array of dimension nsub.
c       On entry ind specifies the indices of subspace variables.
c       On exit ind is unchanged. 
c
c     nenter is an integer variable.
c       On entry nenter is the number of variables entering the 
c         free set.
c       On exit nenter is unchanged. 
c
c     ileave is an integer variable.
c       On entry indx2(ileave),...,indx2(n) are the variables leaving
c         the free set.
c       On exit ileave is unchanged. 
c
c     indx2 is an integer array of dimension n.
c       On entry indx2(1),...,indx2(nenter) are the variables entering
c         the free set, while indx2(ileave),...,indx2(n) are the
c         variables leaving the free set.
c       On exit indx2 is unchanged. 
c
c     iupdat is an integer variable.
c       On entry iupdat is the total number of BFGS updates made so far.
c       On exit iupdat is unchanged. 
c
c     updatd is a logical variable.
c       On entry 'updatd' is true if the L-BFGS matrix is updatd.
c       On exit 'updatd' is unchanged. 
c
c     wn is a double precision array of dimension 2m x 2m.
c       On entry wn is unspecified.
c       On exit the upper triangle of wn stores the LEL^T factorization
c         of the 2*col x 2*col indefinite matrix
c                     [-D -Y'ZZ'Y/theta     L_a'-R_z'  ]
c                     [L_a -R_z           theta*S'AA'S ]
c
c     wn1 is a double precision array of dimension 2m x 2m.
c       On entry wn1 stores the lower triangular part of 
c                     [Y' ZZ'Y   L_a'+R_z']
c                     [L_a+R_z   S'AA'S   ]
c         in the previous iteration.
c       On exit wn1 stores the corresponding updated matrices.
c       The purpose of wn1 is just to store these inner products
c       so they can be easily updated and inserted into wn.
c
c     m is an integer variable.
c       On entry m is the maximum number of variable metric corrections
c         used to define the limited memory matrix.
c       On exit m is unchanged.
c
c     ws, wy, sy, and wtyy are double precision arrays;
c     theta is a double precision variable;
c     col is an integer variable;
c     head is an integer variable.
c       On entry they store the information defining the
c                                          limited memory BFGS matrix:
c         ws(n,m) stores S, a set of s-vectors;
c         wy(n,m) stores Y, a set of y-vectors;
c         sy(m,m) stores S'Y;
c         wtyy(m,m) stores the Cholesky factorization
c                                   of (theta*S'S+LD^(-1)L')
c         theta is the scaling factor specifying B_0 = theta I;
c         col is the number of variable metric corrections stored;
c         head is the location of the 1st s- (or y-) vector in S (or Y).
c       On exit they are unchanged.
c
c     info is an integer variable.
c       On entry info is unspecified.
c       On exit info =  0 for normal return;
c                    = -1 when the 1st Cholesky factorization failed;
c                    = -2 when the 2st Cholesky factorization failed.
c
c     Subprograms called:
c
c       Linpack ... dcopy, dpofa, dtrsl.
c
c
c     References:
c       [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited
c       memory algorithm for bound constrained optimization'',
c       SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.
c
c       [2] C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, ``L-BFGS-B: a
c       limited memory FORTRAN code for solving bound constrained
c       optimization problems'', Tech. Report, NAM-11, EECS Department,
c       Northwestern University, 1994.
c
c       (Postscript files of these papers are available via anonymous
c        ftp to eecs.nwu.edu in the directory pub/lbfgs/lbfgs_bcm.)
c
c                           *  *  *
c
c     NEOS, November 1994. (Latest revision June 1996.)
c     Optimization Technology Center.
c     Argonne National Laboratory and Northwestern University.
c     Written by
c                        Ciyou Zhu
c     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
c
c
c     ************

      integer          m2,ipntr,jpntr,iy,is,jy,js,is1,js1,k1,i,k,
     +                 col2,pbegin,pend,dbegin,dend,upcl
      double precision ddot,temp1,temp2,temp3,temp4
      double precision one,zero
      parameter        (one=1.0d0,zero=0.0d0)

c     Form the lower triangular part of
c               WN1 = [Y' ZZ'Y   L_a'+R_z'] 
c                     [L_a+R_z   S'AA'S   ]
c        where L_a is the strictly lower triangular part of S'AA'Y
c              R_z is the upper triangular part of S'ZZ'Y.
      
      if (updatd) then
         if (iupdat .gt. m) then 
c                                 shift old part of WN1.
            do 10 jy = 1, m - 1
               js = m + jy
               call dcopy(m-jy,wn1(jy+1,jy+1),1,wn1(jy,jy),1)
               call dcopy(m-jy,wn1(js+1,js+1),1,wn1(js,js),1)
               call dcopy(m-1,wn1(m+2,jy+1),1,wn1(m+1,jy),1)
  10        continue
         endif
 
c          put new rows in blocks (1,1), (2,1) and (2,2).
         pbegin = 1
         pend = nsub
         dbegin = nsub + 1
         dend = n
         iy = col
         is = m + col
         ipntr = head + col - 1
         if (ipntr .gt. m) ipntr = ipntr - m    
         jpntr = head
         do 20 jy = 1, col
            js = m + jy
            temp1 = zero
            temp2 = zero
            temp3 = zero
c             compute element jy of row 'col' of Y'ZZ'Y
            do 15 k = pbegin, pend
               k1 = ind(k)
               temp1 = temp1 + wy(k1,ipntr)*wy(k1,jpntr)
  15        continue
c             compute elements jy of row 'col' of L_a and S'AA'S
            do 16 k = dbegin, dend
               k1 = ind(k)
               temp2 = temp2 + ws(k1,ipntr)*ws(k1,jpntr)
               temp3 = temp3 + ws(k1,ipntr)*wy(k1,jpntr)
  16        continue
            wn1(iy,jy) = temp1
            wn1(is,js) = temp2
            wn1(is,jy) = temp3
            jpntr = mod(jpntr,m) + 1
  20     continue
 
c          put new column in block (2,1).
         jy = col       
         jpntr = head + col - 1
         if (jpntr .gt. m) jpntr = jpntr - m
         ipntr = head
         do 30 i = 1, col
            is = m + i
            temp3 = zero
c             compute element i of column 'col' of R_z
            do 25 k = pbegin, pend
               k1 = ind(k)
               temp3 = temp3 + ws(k1,ipntr)*wy(k1,jpntr)
  25        continue 
            ipntr = mod(ipntr,m) + 1
            wn1(is,jy) = temp3
  30     continue
         upcl = col - 1
      else
         upcl = col
      endif
 
c       modify the old parts in blocks (1,1) and (2,2) due to changes
c       in the set of free variables.
      ipntr = head      
      do 45 iy = 1, upcl
         is = m + iy
         jpntr = head
         do 40 jy = 1, iy
            js = m + jy
            temp1 = zero
            temp2 = zero
            temp3 = zero
            temp4 = zero
            do 35 k = 1, nenter
               k1 = indx2(k)
               temp1 = temp1 + wy(k1,ipntr)*wy(k1,jpntr)
               temp2 = temp2 + ws(k1,ipntr)*ws(k1,jpntr)
  35        continue
            do 36 k = ileave, n
               k1 = indx2(k)
               temp3 = temp3 + wy(k1,ipntr)*wy(k1,jpntr)
               temp4 = temp4 + ws(k1,ipntr)*ws(k1,jpntr)
  36        continue
            wn1(iy,jy) = wn1(iy,jy) + temp1 - temp3 
            wn1(is,js) = wn1(is,js) - temp2 + temp4 
            jpntr = mod(jpntr,m) + 1
  40     continue
         ipntr = mod(ipntr,m) + 1
  45  continue
 
c       modify the old parts in block (2,1).
      ipntr = head      
      do 60 is = m + 1, m + upcl
         jpntr = head 
         do 55 jy = 1, upcl
            temp1 = zero
            temp3 = zero
            do 50 k = 1, nenter
               k1 = indx2(k)
               temp1 = temp1 + ws(k1,ipntr)*wy(k1,jpntr)
  50        continue
            do 51 k = ileave, n
               k1 = indx2(k)
               temp3 = temp3 + ws(k1,ipntr)*wy(k1,jpntr)
  51        continue
         if (is .le. jy + m) then
               wn1(is,jy) = wn1(is,jy) + temp1 - temp3  
            else
               wn1(is,jy) = wn1(is,jy) - temp1 + temp3  
            endif
            jpntr = mod(jpntr,m) + 1
  55     continue
         ipntr = mod(ipntr,m) + 1
  60  continue
 
c     Form the upper triangle of WN = [D+Y' ZZ'Y/theta   -L_a'+R_z' ] 
c                                     [-L_a +R_z        S'AA'S*theta]

      m2 = 2*m
      do 70 iy = 1, col
         is = col + iy
         is1 = m + iy
         do 65 jy = 1, iy
            js = col + jy
            js1 = m + jy
            wn(jy,iy) = wn1(iy,jy)/theta
            wn(js,is) = wn1(is1,js1)*theta
  65     continue
         do 66 jy = 1, iy - 1
            wn(jy,is) = -wn1(is1,jy)
  66     continue
         do 67 jy = iy, col
            wn(jy,is) = wn1(is1,jy)
  67     continue
         wn(iy,iy) = wn(iy,iy) + sy(iy,iy)
  70  continue

c     Form the upper triangle of WN= [  LL'            L^-1(-L_a'+R_z')] 
c                                    [(-L_a +R_z)L'^-1   S'AA'S*theta  ]

c        first Cholesky factor (1,1) block of wn to get LL'
c                          with L' stored in the upper triangle of wn.
      call dpofa(wn,m2,col,info)
      if (info .ne. 0) then
         info = -1
         return
      endif
c        then form L^-1(-L_a'+R_z') in the (1,2) block.
      col2 = 2*col
      do 71 js = col+1 ,col2
         call dtrsl(wn,m2,col,wn(1,js),11,info)
  71  continue

c     Form S'AA'S*theta + (L^-1(-L_a'+R_z'))'L^-1(-L_a'+R_z') in the
c        upper triangle of (2,2) block of wn.
                      

      do 72 is = col+1, col2
         do 74 js = is, col2
               wn(is,js) = wn(is,js) + ddot(col,wn(1,is),1,wn(1,js),1)
  74        continue
  72     continue

c     Cholesky factorization of (2,2) block of wn.

      call dpofa(wn(col+1,col+1),m2,col,info)
      if (info .ne. 0) then
         info = -2
         return
      endif

      return

      end

c======================= The end of formk ==============================

      subroutine formt(m, wt, sy, ss, col, theta, info)
 
      integer          m, col, info
      double precision theta, wt(m, m), sy(m, m), ss(m, m)

c     ************
c
c     Subroutine formt
c
c       This subroutine forms the upper half of the pos. def. and symm.
c         T = theta*SS + L*D^(-1)*L', stores T in the upper triangle
c         of the array wt, and performs the Cholesky factorization of T
c         to produce J*J', with J' stored in the upper triangle of wt.
c
c     Subprograms called:
c
c       Linpack ... dpofa.
c
c
c                           *  *  *
c
c     NEOS, November 1994. (Latest revision June 1996.)
c     Optimization Technology Center.
c     Argonne National Laboratory and Northwestern University.
c     Written by
c                        Ciyou Zhu
c     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
c
c
c     ************

      integer          i,j,k,k1
      double precision ddum
      double precision zero
      parameter        (zero=0.0d0)


c     Form the upper half of  T = theta*SS + L*D^(-1)*L',
c        store T in the upper triangle of the array wt.
 
      do 52 j = 1, col
         wt(1,j) = theta*ss(1,j)
  52  continue
      do 55 i = 2, col
         do 54 j = i, col
            k1 = min(i,j) - 1
            ddum  = zero
            do 53 k = 1, k1
               ddum  = ddum + sy(i,k)*sy(j,k)/sy(k,k)
  53        continue
            wt(i,j) = ddum + theta*ss(i,j)
  54     continue
  55  continue
 
c     Cholesky factorize T to J*J' with 
c        J' stored in the upper triangle of wt.
 
      call dpofa(wt,m,col,info)
      if (info .ne. 0) then
         info = -3
      endif

      return

      end

c======================= The end of formt ==============================
 
      subroutine freev(n, nfree, index, nenter, ileave, indx2, 
     +                 iwhere, wrk, updatd, cnstnd, iprint, iter)

      integer n, nfree, nenter, ileave, iprint, iter, 
     +        index(n), indx2(n), iwhere(n)
      logical wrk, updatd, cnstnd

c     ************
c
c     Subroutine freev 
c
c     This subroutine counts the entering and leaving variables when
c       iter > 0, and finds the index set of free and active variables
c       at the GCP.
c
c     cnstnd is a logical variable indicating whether bounds are present
c
c     index is an integer array of dimension n
c       for i=1,...,nfree, index(i) are the indices of free variables
c       for i=nfree+1,...,n, index(i) are the indices of bound variables
c       On entry after the first iteration, index gives 
c         the free variables at the previous iteration.
c       On exit it gives the free variables based on the determination
c         in cauchy using the array iwhere.
c
c     indx2 is an integer array of dimension n
c       On entry indx2 is unspecified.
c       On exit with iter>0, indx2 indicates which variables
c          have changed status since the previous iteration.
c       For i= 1,...,nenter, indx2(i) have changed from bound to free.
c       For i= ileave+1,...,n, indx2(i) have changed from free to bound.
c 
c
c                           *  *  *
c
c     NEOS, November 1994. (Latest revision June 1996.)
c     Optimization Technology Center.
c     Argonne National Laboratory and Northwestern University.
c     Written by
c                        Ciyou Zhu
c     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
c
c
c     ************
 
      integer iact,i,k

      nenter = 0
      ileave = n + 1
      if (iter .gt. 0 .and. cnstnd) then
c                           count the entering and leaving variables.
         do 20 i = 1, nfree
            k = index(i)

c            write(6,*) ' k  = index(i) ', k
c            write(6,*) ' index = ', i

            if (iwhere(k) .gt. 0) then
               ileave = ileave - 1
               indx2(ileave) = k
               if (iprint .ge. 100) write (6,*)
     +             'Variable ',k,' leaves the set of free variables'
            endif
  20     continue
         do 22 i = 1 + nfree, n
            k = index(i)
            if (iwhere(k) .le. 0) then
               nenter = nenter + 1
               indx2(nenter) = k
               if (iprint .ge. 100) write (6,*)
     +             'Variable ',k,' enters the set of free variables'
            endif
  22     continue
         if (iprint .ge. 99) write (6,*)
     +       n+1-ileave,' variables leave; ',nenter,' variables enter'
      endif
      wrk = (ileave .lt. n+1) .or. (nenter .gt. 0) .or. updatd
 
c     Find the index set of free and active variables at the GCP.
 
      nfree = 0 
      iact = n + 1
      do 24 i = 1, n
         if (iwhere(i) .le. 0) then
            nfree = nfree + 1
            index(nfree) = i
         else
            iact = iact - 1
            index(iact) = i
         endif
  24  continue
      if (iprint .ge. 99) write (6,*)
     +      nfree,' variables are free at GCP ',iter + 1  

      return

      end

c======================= The end of freev ==============================

      subroutine hpsolb(n, t, iorder, iheap)
      integer          iheap, n, iorder(n)
      double precision t(n)
  
c     ************
c
c     Subroutine hpsolb 
c
c     This subroutine sorts out the least element of t, and puts the
c       remaining elements of t in a heap.
c 
c     n is an integer variable.
c       On entry n is the dimension of the arrays t and iorder.
c       On exit n is unchanged.
c
c     t is a double precision array of dimension n.
c       On entry t stores the elements to be sorted,
c       On exit t(n) stores the least elements of t, and t(1) to t(n-1)
c         stores the remaining elements in the form of a heap.
c
c     iorder is an integer array of dimension n.
c       On entry iorder(i) is the index of t(i).
c       On exit iorder(i) is still the index of t(i), but iorder may be
c         permuted in accordance with t.
c
c     iheap is an integer variable specifying the task.
c       On entry iheap should be set as follows:
c         iheap .eq. 0 if t(1) to t(n) is not in the form of a heap,
c         iheap .ne. 0 if otherwise.
c       On exit iheap is unchanged.
c
c
c     References:
c       Algorithm 232 of CACM (J. W. J. Williams): HEAPSORT.
c
c                           *  *  *
c
c     NEOS, November 1994. (Latest revision June 1996.)
c     Optimization Technology Center.
c     Argonne National Laboratory and Northwestern University.
c     Written by
c                        Ciyou Zhu
c     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
c
c     ************
  
      integer          i,j,k,indxin,indxou
      double precision ddum,out

      if (iheap .eq. 0) then

c        Rearrange the elements t(1) to t(n) to form a heap.

         do 20 k = 2, n
            ddum  = t(k)
            indxin = iorder(k)

c           Add ddum to the heap.
            i = k
   10       continue
            if (i.gt.1) then
               j = i/2
               if (ddum .lt. t(j)) then
                  t(i) = t(j)
                  iorder(i) = iorder(j)
                  i = j
                  goto 10 
               endif  
            endif  
            t(i) = ddum
            iorder(i) = indxin
   20    continue
      endif
 
c     Assign to 'out' the value of t(1), the least member of the heap,
c        and rearrange the remaining members to form a heap as
c        elements 1 to n-1 of t.
 
      if (n .gt. 1) then
         i = 1
         out = t(1)
         indxou = iorder(1)
         ddum  = t(n)
         indxin  = iorder(n)

c        Restore the heap 
   30    continue
         j = i+i
         if (j .le. n-1) then
            if (t(j+1) .lt. t(j)) j = j+1
            if (t(j) .lt. ddum ) then
               t(i) = t(j)
               iorder(i) = iorder(j)
               i = j
               goto 30
            endif 
         endif 
         t(i) = ddum
         iorder(i) = indxin
 
c     Put the least member in t(n). 

         t(n) = out
         iorder(n) = indxou
      endif 

      return

      end

c====================== The end of hpsolb ==============================

      subroutine lnsrlb(n, l, u, nbd, x, f, fold, gd, gdold, g, d, r, t,
     +                  z, stp, dnorm, dtd, xstep, stpmx, iter, ifun,
     +                  iback, nfgv, info, task, boxed, cnstnd, csave,
     +                  isave, dsave)

      character*60     task, csave
      logical          boxed, cnstnd
      integer          n, iter, ifun, iback, nfgv, info,
     +                 nbd(n), isave(2)
      double precision f, fold, gd, gdold, stp, dnorm, dtd, xstep,
     +                 stpmx, x(n), l(n), u(n), g(n), d(n), r(n), t(n),
     +                 z(n), dsave(13)
c     **********
c
c     Subroutine lnsrlb
c
c     This subroutine calls subroutine dcsrch from the Minpack2 library
c       to perform the line search.  Subroutine dscrch is safeguarded so
c       that all trial points lie within the feasible region.
c
c     Subprograms called:
c
c       Minpack2 Library ... dcsrch.
c
c       Linpack ... dtrsl, ddot.
c
c
c                           *  *  *
c
c     NEOS, November 1994. (Latest revision June 1996.)
c     Optimization Technology Center.
c     Argonne National Laboratory and Northwestern University.
c     Written by
c                        Ciyou Zhu
c     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
c
c
c     **********

      integer          i
      double           precision ddot,a1,a2
      double precision one,zero,big
      parameter        (one=1.0d0,zero=0.0d0,big=1.0d+10)
      double precision ftol,gtol,xtol
      parameter        (ftol=1.0d-3,gtol=0.9d0,xtol=0.1d0)

      if (task(1:5) .eq. 'FG_LN') goto 556

      dtd = ddot(n,d,1,d,1)
      dnorm = sqrt(dtd)

c     Determine the maximum step length.

      stpmx = big
      if (cnstnd) then
         if (iter .eq. 0) then
            stpmx = one
         else
            do 43 i = 1, n
               a1 = d(i)
               if (nbd(i) .ne. 0) then
                  if (a1 .lt. zero .and. nbd(i) .le. 2) then
                     a2 = l(i) - x(i)
                     if (a2 .ge. zero) then
                        stpmx = zero
                     else if (a1*stpmx .lt. a2) then
                        stpmx = a2/a1
                     endif
                  else if (a1 .gt. zero .and. nbd(i) .ge. 2) then
                     a2 = u(i) - x(i)
                     if (a2 .le. zero) then
                        stpmx = zero
                     else if (a1*stpmx .gt. a2) then
                        stpmx = a2/a1
                     endif
                  endif
               endif
  43        continue
         endif
      endif
 
      if (iter .eq. 0 .and. .not. boxed) then
         stp = min(one/dnorm, stpmx)
      else
         stp = one
      endif 

      call dcopy(n,x,1,t,1)
      call dcopy(n,g,1,r,1)
      fold = f
      ifun = 0
      iback = 0
      csave = 'START'
 556  continue
      gd = ddot(n,g,1,d,1)
      if (ifun .eq. 0) then
         gdold=gd
         if (gd .ge. zero) then
c                               the directional derivative >=0.
c                               Line search is impossible.
            write(6,*)' ascent direction in projection gd = ', gd
            info = -4
            return
         endif
      endif

      call dcsrch(f,gd,stp,ftol,gtol,xtol,zero,stpmx,csave,isave,dsave)

      xstep = stp*dnorm
      if (csave(1:4) .ne. 'CONV' .and. csave(1:4) .ne. 'WARN') then
         task = 'FG_LNSRCH'
         ifun = ifun + 1
         nfgv = nfgv + 1
         iback = ifun - 1 
         if (stp .eq. one) then
            call dcopy(n,z,1,x,1)
         else
            do 41 i = 1, n
               x(i) = stp*d(i) + t(i)
  41        continue
         endif
      else
         task = 'NEW_X'
      endif

      return

      end

c======================= The end of lnsrlb =============================

      subroutine matupd(n, m, ws, wy, sy, ss, d, r, itail, 
     +                  iupdat, col, head, theta, rr, dr, stp, dtd)
 
      integer          n, m, itail, iupdat, col, head
      double precision theta, rr, dr, stp, dtd, d(n), r(n), 
     +                 ws(n, m), wy(n, m), sy(m, m), ss(m, m)

c     ************
c
c     Subroutine matupd
c
c       This subroutine updates matrices WS and WY, and forms the
c         middle matrix in B.
c
c     Subprograms called:
c
c       Linpack ... dcopy, ddot.
c
c
c                           *  *  *
c
c     NEOS, November 1994. (Latest revision June 1996.)
c     Optimization Technology Center.
c     Argonne National Laboratory and Northwestern University.
c     Written by
c                        Ciyou Zhu
c     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
c
c
c     ************
 
      integer          j,pointr
      double precision ddot
      double precision one
      parameter        (one=1.0d0)

c     Set pointers for matrices WS and WY.
 
      if (iupdat .le. m) then
         col = iupdat
         itail = mod(head+iupdat-2,m) + 1
      else
         itail = mod(itail,m) + 1
         head = mod(head,m) + 1
      endif
 
c     Update matrices WS and WY.

      call dcopy(n,d,1,ws(1,itail),1)
      call dcopy(n,r,1,wy(1,itail),1)
 
c     Set theta=yy/ys.
 
      theta = rr/dr
 
c     Form the middle matrix in B.
 
c        update the upper triangle of SS,
c                                         and the lower triangle of SY:
      if (iupdat .gt. m) then
c                              move old information
         do 50 j = 1, col - 1
            call dcopy(j,ss(2,j+1),1,ss(1,j),1)
            call dcopy(col-j,sy(j+1,j+1),1,sy(j,j),1)
  50     continue
      endif
c        add new information: the last row of SY
c                                             and the last column of SS:
      pointr = head
      do 51 j = 1, col - 1
         sy(col,j) = ddot(n,d,1,wy(1,pointr),1)
         ss(j,col) = ddot(n,ws(1,pointr),1,d,1)
         pointr = mod(pointr,m) + 1
  51  continue
      if (stp .eq. one) then
         ss(col,col) = dtd
      else
         ss(col,col) = stp*stp*dtd
      endif
      sy(col,col) = dr
 
      return

      end

c======================= The end of matupd =============================

      subroutine prn1lb(n, m, l, u, x, iprint, itfile, epsmch)
 
      integer n, m, iprint, itfile
      double precision epsmch, x(n), l(n), u(n)

c     ************
c
c     Subroutine prn1lb
c
c     This subroutine prints the input data, initial point, upper and
c       lower bounds of each variable, machine precision, as well as 
c       the headings of the output.
c
c
c                           *  *  *
c
c     NEOS, November 1994. (Latest revision June 1996.)
c     Optimization Technology Center.
c     Argonne National Laboratory and Northwestern University.
c     Written by
c                        Ciyou Zhu
c     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
c
c
c     ************

      integer i

      if (iprint .ge. 0) then
         write (6,7001) epsmch
         write (6,*) 'N = ',n,'    M = ',m
         if (iprint .ge. 1) then
            write (itfile,2001) epsmch
            write (itfile,*)'N = ',n,'    M = ',m
            write (itfile,9001)
            if (iprint .gt. 100) then
               write (6,1004) 'L =',(l(i),i = 1,n)
               write (6,1004) 'X0 =',(x(i),i = 1,n)
               write (6,1004) 'U =',(u(i),i = 1,n)
            endif 
         endif
      endif 

 1004 format (/,a4, 1p, 6(1x,d11.4),/,(4x,1p,6(1x,d11.4)))
 2001 format ('RUNNING THE L-BFGS-B CODE',/,/,
     + 'it    = iteration number',/,
     + 'nf    = number of function evaluations',/,
     + 'nseg  = number of segments explored during the Cauchy search',/,
     + 'nact  = number of active bounds at the generalized Cauchy point'
     + ,/,
     + 'sub   = manner in which the subspace minimization terminated:'
     + ,/,'        con = converged, bnd = a bound was reached',/,
     + 'itls  = number of iterations performed in the line search',/,
     + 'stepl = step length used',/,
     + 'tstep = norm of the displacement (total step)',/,
     + 'projg = norm of the projected gradient',/,
     + 'f     = function value',/,/,
     + '           * * *',/,/,
     + 'Machine precision =',1p,d10.3)
 7001 format ('RUNNING THE L-BFGS-B CODE',/,/,
     + '           * * *',/,/,
     + 'Machine precision =',1p,d10.3)
 9001 format (/,3x,'it',3x,'nf',2x,'nseg',2x,'nact',2x,'sub',2x,'itls',
     +        2x,'stepl',4x,'tstep',5x,'projg',8x,'f')

      return

      end

c======================= The end of prn1lb =============================

      subroutine prn2lb(n, x, f, g, iprint, itfile, iter, nfgv, nact, 
     +                  sbgnrm, nseg, word, iword, iback, stp, xstep)
 
      character*3      word
      integer          n, iprint, itfile, iter, nfgv, nact, nseg,
     +                 iword, iback
      double precision f, sbgnrm, stp, xstep, x(n), g(n)

c     ************
c
c     Subroutine prn2lb
c
c     This subroutine prints out new information after a successful
c       line search. 
c
c
c                           *  *  *
c
c     NEOS, November 1994. (Latest revision June 1996.)
c     Optimization Technology Center.
c     Argonne National Laboratory and Northwestern University.
c     Written by
c                        Ciyou Zhu
c     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
c
c
c     ************

      integer i,imod

c           'word' records the status of subspace solutions.
      if (iword .eq. 0) then
c                            the subspace minimization converged.
         word = 'con'
      else if (iword .eq. 1) then
c                          the subspace minimization stopped at a bound.
         word = 'bnd'
      else if (iword .eq. 5) then
c                             the truncated Newton step has been used.
         word = 'TNT'
      else
         word = '---'
      endif
      if (iprint .ge. 99) then
         write (6,*) 'LINE SEARCH',iback,' times; norm of step = ',xstep
         write (6,2001) iter,f,sbgnrm
         if (iprint .gt. 100) then      
            write (6,1004) 'X =',(x(i), i = 1, n)
            write (6,1004) 'G =',(g(i), i = 1, n)
         endif
      else if (iprint .gt. 0) then 
         imod = mod(iter,iprint)
         if (imod .eq. 0) write (6,2001) iter,f,sbgnrm
      endif
      if (iprint .ge. 1) write (itfile,3001)
     +          iter,nfgv,nseg,nact,word,iback,stp,xstep,sbgnrm,f

 1004 format (/,a4, 1p, 6(1x,d11.4),/,(4x,1p,6(1x,d11.4)))
 2001 format
     +  (/,'At iterate',i5,4x,'f= ',1p,d12.5,4x,'|proj g|= ',1p,d12.5)
 3001 format(2(1x,i4),2(1x,i5),2x,a3,1x,i4,1p,2(2x,d8.1),1p,2(1x,d10.3))

      return

      end

c======================= The end of prn2lb =============================

      subroutine prn3lb(n, x, f, task, iprint, info, itfile, 
     +                  iter, nfgv, nintol, nskip, nact, sbgnrm, 
     +                  time, nseg, word, iback, stp, xstep, k, 
     +                  cachyt, sbtime, lnscht)
 
      character*60     task
      character*3      word
      integer          n, iprint, info, itfile, iter, nfgv, nintol,
     +                 nskip, nact, nseg, iback, k
      double precision f, sbgnrm, time, stp, xstep, cachyt, sbtime,
     +                 lnscht, x(n)

c     ************
c
c     Subroutine prn3lb
c
c     This subroutine prints out information when either a built-in
c       convergence test is satisfied or when an error message is
c       generated.
c       
c
c                           *  *  *
c
c     NEOS, November 1994. (Latest revision June 1996.)
c     Optimization Technology Center.
c     Argonne National Laboratory and Northwestern University.
c     Written by
c                        Ciyou Zhu
c     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
c
c
c     ************

      integer i

      if (task(1:5) .eq. 'ERROR') goto 999

      if (iprint .ge. 0) then
         write (6,3003)
         write (6,3004)
         write(6,3005) n,iter,nfgv,nintol,nskip,nact,sbgnrm,f
         if (iprint .ge. 100) then
            write (6,1004) 'X =',(x(i),i = 1,n)
         endif  
         if (iprint .ge. 1) write (6,*) ' F =',f
      endif 
 999  continue
      if (iprint .ge. 0) then
         write (6,3009) task
         if (info .ne. 0) then
            if (info .eq. -1) write (6,9011)
            if (info .eq. -2) write (6,9012)
            if (info .eq. -3) write (6,9013)
            if (info .eq. -4) write (6,9014)
            if (info .eq. -5) write (6,9015)
            if (info .eq. -6) write (6,*)' Input nbd(',k,') is invalid.'
            if (info .eq. -7) 
     +      write (6,*)' l(',k,') > u(',k,').  No feasible solution.'
            if (info .eq. -8) write (6,9018)
            if (info .eq. -9) write (6,9019)
         endif
         if (iprint .ge. 1) write (6,3007) cachyt,sbtime,lnscht
         write (6,3008) time
         if (iprint .ge. 1) then
            if (info .eq. -4 .or. info .eq. -9) then
               write (itfile,3002)
     +             iter,nfgv,nseg,nact,word,iback,stp,xstep
            endif
            write (itfile,3009) task
            if (info .ne. 0) then
               if (info .eq. -1) write (itfile,9011)
               if (info .eq. -2) write (itfile,9012)
               if (info .eq. -3) write (itfile,9013)
               if (info .eq. -4) write (itfile,9014)
               if (info .eq. -5) write (itfile,9015)
               if (info .eq. -8) write (itfile,9018)
               if (info .eq. -9) write (itfile,9019)
            endif
            write (itfile,3008) time
         endif
      endif

 1004 format (/,a4, 1p, 6(1x,d11.4),/,(4x,1p,6(1x,d11.4)))
 3002 format(2(1x,i4),2(1x,i5),2x,a3,1x,i4,1p,2(2x,d8.1),6x,'-',10x,'-')
 3003 format (/,
     + '           * * *',/,/,
     + 'Tit   = total number of iterations',/,
     + 'Tnf   = total number of function evaluations',/,
     + 'Tnint = total number of segments explored during',
     +           ' Cauchy searches',/,
     + 'Skip  = number of BFGS updates skipped',/,
     + 'Nact  = number of active bounds at final generalized',
     +          ' Cauchy point',/,
     + 'Projg = norm of the final projected gradient',/,
     + 'F     = final function value',/,/,
     + '           * * *')
 3004 format (/,3x,'N',4x,'Tit',5x,'Tnf',2x,'Tnint',2x,
     +       'Skip',2x,'Nact',5x,'Projg',8x,'F')
 3005 format (i5,2(1x,i6),(1x,i6),(2x,i4),(1x,i5),1p,2(2x,d10.3))
 3007 format (/,' Cauchy                time',1p,e10.3,' seconds.',/ 
     +        ' Subspace minimization time',1p,e10.3,' seconds.',/
     +        ' Line search           time',1p,e10.3,' seconds.')
 3008 format (/,' Total User time',1p,e10.3,' seconds.',/)
 3009 format (/,a60)
 9011 format (/,
     +' Matrix in 1st Cholesky factorization in formk is not Pos. Def.')
 9012 format (/,
     +' Matrix in 2st Cholesky factorization in formk is not Pos. Def.')
 9013 format (/,
     +' Matrix in the Cholesky factorization in formt is not Pos. Def.')
 9014 format (/,
     +' Derivative >= 0, backtracking line search impossible.',/,
     +'   Previous x, f and g restored.',/,
     +' Possible causes: 1 error in function or gradient evaluation;',/,
     +'                  2 rounding errors dominate computation.')
 9015 format (/,
     +' Warning:  more than 10 function and gradient',/,
     +'   evaluations in the last line search.  Termination',/,
     +'   may possibly be caused by a bad search direction.')
 9018 format (/,' The triangular system is singular.')
 9019 format (/,
     +' Line search cannot locate an adequate point after 20 function',/
     +,'  and gradient evaluations.  Previous x, f and g restored.',/,
     +' Possible causes: 1 error in function or gradient evaluation;',/,
     +'                  2 rounding error dominate computation.')

      return

      end

c======================= The end of prn3lb =============================

      subroutine projgr(n, l, u, nbd, x, g, sbgnrm)

      integer          n, nbd(n)
      double precision sbgnrm, x(n), l(n), u(n), g(n)

c     ************
c
c     Subroutine projgr
c
c     This subroutine computes the infinity norm of the projected
c       gradient.
c
c
c                           *  *  *
c
c     NEOS, November 1994. (Latest revision June 1996.)
c     Optimization Technology Center.
c     Argonne National Laboratory and Northwestern University.
c     Written by
c                        Ciyou Zhu
c     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
c
c
c     ************

      integer i
      double precision gi
      double precision one,zero
      parameter        (one=1.0d0,zero=0.0d0)

      sbgnrm = zero
      do 15 i = 1, n
        gi = g(i)
        if (nbd(i) .ne. 0) then
           if (gi .lt. zero) then
              if (nbd(i) .ge. 2) gi = max((x(i)-u(i)),gi)
           else
              if (nbd(i) .le. 2) gi = min((x(i)-l(i)),gi)
           endif
        endif
        sbgnrm = max(sbgnrm,abs(gi))
  15  continue

      return

      end

c======================= The end of projgr =============================

      subroutine subsm ( n, m, nsub, ind, l, u, nbd, x, d, xp, ws, wy,
     +                   theta, xx, gg,
     +                   col, head, iword, wv, wn, iprint, info )
      implicit none
      integer          n, m, nsub, col, head, iword, iprint, info, 
     +                 ind(nsub), nbd(n)
      double precision theta, 
     +                 l(n), u(n), x(n), d(n), xp(n), xx(n), gg(n),
     +                 ws(n, m), wy(n, m), 
     +                 wv(2*m), wn(2*m, 2*m)

c     **********************************************************************
c
c     This routine contains the major changes in the updated version.
c     The changes are described in the accompanying paper
c
c      Jose Luis Morales, Jorge Nocedal
c      "Remark On Algorithm 788: L-BFGS-B: Fortran Subroutines for Large-Scale
c       Bound Constrained Optimization". Decemmber 27, 2010.
c
c             J.L. Morales  Departamento de Matematicas, 
c                           Instituto Tecnologico Autonomo de Mexico
c                           Mexico D.F.
c
c             J, Nocedal    Department of Electrical Engineering and
c                           Computer Science.
c                           Northwestern University. Evanston, IL. USA
c
c                           January 17, 2011
c
c      **********************************************************************
c                           
c
c     Subroutine subsm
c
c     Given xcp, l, u, r, an index set that specifies
c       the active set at xcp, and an l-BFGS matrix B 
c       (in terms of WY, WS, SY, WT, head, col, and theta), 
c       this subroutine computes an approximate solution
c       of the subspace problem
c
c       (P)   min Q(x) = r'(x-xcp) + 1/2 (x-xcp)' B (x-xcp)
c
c             subject to l<=x<=u
c                       x_i=xcp_i for all i in A(xcp)
c                     
c       along the subspace unconstrained Newton direction 
c       
c          d = -(Z'BZ)^(-1) r.
c
c       The formula for the Newton direction, given the L-BFGS matrix
c       and the Sherman-Morrison formula, is
c
c          d = (1/theta)r + (1/theta*2) Z'WK^(-1)W'Z r.
c 
c       where
c                 K = [-D -Y'ZZ'Y/theta     L_a'-R_z'  ]
c                     [L_a -R_z           theta*S'AA'S ]
c
c     Note that this procedure for computing d differs 
c     from that described in [1]. One can show that the matrix K is
c     equal to the matrix M^[-1]N in that paper.
c
c     n is an integer variable.
c       On entry n is the dimension of the problem.
c       On exit n is unchanged.
c
c     m is an integer variable.
c       On entry m is the maximum number of variable metric corrections
c         used to define the limited memory matrix.
c       On exit m is unchanged.
c
c     nsub is an integer variable.
c       On entry nsub is the number of free variables.
c       On exit nsub is unchanged.
c
c     ind is an integer array of dimension nsub.
c       On entry ind specifies the coordinate indices of free variables.
c       On exit ind is unchanged.
c
c     l is a double precision array of dimension n.
c       On entry l is the lower bound of x.
c       On exit l is unchanged.
c
c     u is a double precision array of dimension n.
c       On entry u is the upper bound of x.
c       On exit u is unchanged.
c
c     nbd is a integer array of dimension n.
c       On entry nbd represents the type of bounds imposed on the
c         variables, and must be specified as follows:
c         nbd(i)=0 if x(i) is unbounded,
c                1 if x(i) has only a lower bound,
c                2 if x(i) has both lower and upper bounds, and
c                3 if x(i) has only an upper bound.
c       On exit nbd is unchanged.
c
c     x is a double precision array of dimension n.
c       On entry x specifies the Cauchy point xcp. 
c       On exit x(i) is the minimizer of Q over the subspace of
c                                                        free variables. 
c
c     d is a double precision array of dimension n.
c       On entry d is the reduced gradient of Q at xcp.
c       On exit d is the Newton direction of Q. 
c
c    xp is a double precision array of dimension n.
c       used to safeguard the projected Newton direction 
c
c    xx is a double precision array of dimension n
c       On entry it holds the current iterate
c       On output it is unchanged

c    gg is a double precision array of dimension n
c       On entry it holds the gradient at the current iterate
c       On output it is unchanged
c
c     ws and wy are double precision arrays;
c     theta is a double precision variable;
c     col is an integer variable;
c     head is an integer variable.
c       On entry they store the information defining the
c                                          limited memory BFGS matrix:
c         ws(n,m) stores S, a set of s-vectors;
c         wy(n,m) stores Y, a set of y-vectors;
c         theta is the scaling factor specifying B_0 = theta I;
c         col is the number of variable metric corrections stored;
c         head is the location of the 1st s- (or y-) vector in S (or Y).
c       On exit they are unchanged.
c
c     iword is an integer variable.
c       On entry iword is unspecified.
c       On exit iword specifies the status of the subspace solution.
c         iword = 0 if the solution is in the box,
c                 1 if some bound is encountered.
c
c     wv is a double precision working array of dimension 2m.
c
c     wn is a double precision array of dimension 2m x 2m.
c       On entry the upper triangle of wn stores the LEL^T factorization
c         of the indefinite matrix
c
c              K = [-D -Y'ZZ'Y/theta     L_a'-R_z'  ]
c                  [L_a -R_z           theta*S'AA'S ]
c                                                    where E = [-I  0]
c                                                              [ 0  I]
c       On exit wn is unchanged.
c
c     iprint is an INTEGER variable that must be set by the user.
c       It controls the frequency and type of output generated:
c        iprint<0    no output is generated;
c        iprint=0    print only one line at the last iteration;
c        0<iprint<99 print also f and |proj g| every iprint iterations;
c        iprint=99   print details of every iteration except n-vectors;
c        iprint=100  print also the changes of active set and final x;
c        iprint>100  print details of every iteration including x and g;
c       When iprint > 0, the file iterate.dat will be created to
c                        summarize the iteration.
c
c     info is an integer variable.
c       On entry info is unspecified.
c       On exit info = 0       for normal return,
c                    = nonzero for abnormal return 
c                                  when the matrix K is ill-conditioned.
c
c     Subprograms called:
c
c       Linpack dtrsl.
c
c
c     References:
c
c       [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited
c       memory algorithm for bound constrained optimization'',
c       SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.
c
c
c
c                           *  *  *
c
c     NEOS, November 1994. (Latest revision June 1996.)
c     Optimization Technology Center.
c     Argonne National Laboratory and Northwestern University.
c     Written by
c                        Ciyou Zhu
c     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
c
c
c     ************

      integer          pointr,m2,col2,ibd,jy,js,i,j,k
      double precision alpha, xk, dk, temp1, temp2 
      double precision one,zero
      parameter        (one=1.0d0,zero=0.0d0)
c
      double precision dd_p

      if (nsub .le. 0) return
      if (iprint .ge. 99) write (6,1001)

c     Compute wv = W'Zd.

      pointr = head 
      do 20 i = 1, col
         temp1 = zero
         temp2 = zero
         do 10 j = 1, nsub
            k = ind(j)
            temp1 = temp1 + wy(k,pointr)*d(j)
            temp2 = temp2 + ws(k,pointr)*d(j)
  10     continue
         wv(i) = temp1
         wv(col + i) = theta*temp2
         pointr = mod(pointr,m) + 1
  20  continue
 
c     Compute wv:=K^(-1)wv.

      m2 = 2*m
      col2 = 2*col
      call dtrsl(wn,m2,col2,wv,11,info)
      if (info .ne. 0) return
      do 25 i = 1, col
         wv(i) = -wv(i)
  25     continue
      call dtrsl(wn,m2,col2,wv,01,info)
      if (info .ne. 0) return
 
c     Compute d = (1/theta)d + (1/theta**2)Z'W wv.
 
      pointr = head
      do 40 jy = 1, col
         js = col + jy
         do 30 i = 1, nsub
            k = ind(i)
            d(i) = d(i) + wy(k,pointr)*wv(jy)/theta     
     +                  + ws(k,pointr)*wv(js)
  30     continue
         pointr = mod(pointr,m) + 1
  40  continue

      call dscal( nsub, one/theta, d, 1 )
c 
c-----------------------------------------------------------------
c     Let us try the projection, d is the Newton direction

      iword = 0

      call dcopy ( n, x, 1, xp, 1 )
c
      do 50 i=1, nsub
         k  = ind(i)
         dk = d(i)
         xk = x(k)
         if ( nbd(k) .ne. 0 ) then
c
            if ( nbd(k).eq.1 ) then          ! lower bounds only
               x(k) = max( l(k), xk + dk )
               if ( x(k).eq.l(k) ) iword = 1
            else 
c     
               if ( nbd(k).eq.2 ) then       ! upper and lower bounds
                  xk   = max( l(k), xk + dk ) 
                  x(k) = min( u(k), xk )
                  if ( x(k).eq.l(k) .or. x(k).eq.u(k) ) iword = 1
               else
c
                  if ( nbd(k).eq.3 ) then    ! upper bounds only
                     x(k) = min( u(k), xk + dk )
                     if ( x(k).eq.u(k) ) iword = 1
                  end if 
               end if
            end if
c            
         else                                ! free variables
            x(k) = xk + dk
         end if 
 50   continue
c
      if ( iword.eq.0 ) then
         go to 911
      end if
c
c     check sign of the directional derivative
c
      dd_p = zero
      do 55 i=1, n
         dd_p  = dd_p + (x(i) - xx(i))*gg(i)
 55   continue
      if ( dd_p .gt.zero ) then
         call dcopy( n, xp, 1, x, 1 )
         write(6,*) ' Positive dir derivative in projection '
         write(6,*) ' Using the backtracking step '
      else
         go to 911
      endif
c
c-----------------------------------------------------------------
c
      alpha = one
      temp1 = alpha
      ibd   = 0 
      do 60 i = 1, nsub
         k = ind(i)
         dk = d(i)
         if (nbd(k) .ne. 0) then
            if (dk .lt. zero .and. nbd(k) .le. 2) then
               temp2 = l(k) - x(k)
               if (temp2 .ge. zero) then
                  temp1 = zero
               else if (dk*alpha .lt. temp2) then
                  temp1 = temp2/dk
               endif
            else if (dk .gt. zero .and. nbd(k) .ge. 2) then
               temp2 = u(k) - x(k)
               if (temp2 .le. zero) then
                  temp1 = zero
               else if (dk*alpha .gt. temp2) then
                  temp1 = temp2/dk
               endif
            endif
            if (temp1 .lt. alpha) then
               alpha = temp1
               ibd = i
            endif
         endif
 60   continue
      
      if (alpha .lt. one) then
         dk = d(ibd)
         k = ind(ibd)
         if (dk .gt. zero) then
            x(k) = u(k)
            d(ibd) = zero
         else if (dk .lt. zero) then
            x(k) = l(k)
            d(ibd) = zero
         endif
      endif
      do 70 i = 1, nsub
         k    = ind(i)
         x(k) = x(k) + alpha*d(i)
 70   continue
cccccc
 911  continue

      if (iprint .ge. 99) write (6,1004)

 1001 format (/,'----------------SUBSM entered-----------------',/)
 1004 format (/,'----------------exit SUBSM --------------------',/)

      return

      end
c====================== The end of subsm ===============================

      subroutine dcsrch(f,g,stp,ftol,gtol,xtol,stpmin,stpmax,
     +                  task,isave,dsave)
      character*(*) task
      integer isave(2)
      double precision f,g,stp,ftol,gtol,xtol,stpmin,stpmax
      double precision dsave(13)
c     **********
c
c     Subroutine dcsrch
c
c     This subroutine finds a step that satisfies a sufficient
c     decrease condition and a curvature condition.
c
c     Each call of the subroutine updates an interval with 
c     endpoints stx and sty. The interval is initially chosen 
c     so that it contains a minimizer of the modified function
c
c           psi(stp) = f(stp) - f(0) - ftol*stp*f'(0).
c
c     If psi(stp) <= 0 and f'(stp) >= 0 for some step, then the
c     interval is chosen so that it contains a minimizer of f. 
c
c     The algorithm is designed to find a step that satisfies 
c     the sufficient decrease condition 
c
c           f(stp) <= f(0) + ftol*stp*f'(0),
c
c     and the curvature condition
c
c           abs(f'(stp)) <= gtol*abs(f'(0)).
c
c     If ftol is less than gtol and if, for example, the function
c     is bounded below, then there is always a step which satisfies
c     both conditions. 
c
c     If no step can be found that satisfies both conditions, then 
c     the algorithm stops with a warning. In this case stp only 
c     satisfies the sufficient decrease condition.
c
c     A typical invocation of dcsrch has the following outline:
c
c     task = 'START'
c  10 continue
c        call dcsrch( ... )
c        if (task .eq. 'FG') then
c           Evaluate the function and the gradient at stp 
c           goto 10
c           end if
c
c     NOTE: The user must no alter work arrays between calls.
c
c     The subroutine statement is
c
c        subroutine dcsrch(f,g,stp,ftol,gtol,xtol,stpmin,stpmax,
c                          task,isave,dsave)
c     where
c
c       f is a double precision variable.
c         On initial entry f is the value of the function at 0.
c            On subsequent entries f is the value of the 
c            function at stp.
c         On exit f is the value of the function at stp.
c
c       g is a double precision variable.
c         On initial entry g is the derivative of the function at 0.
c            On subsequent entries g is the derivative of the 
c            function at stp.
c         On exit g is the derivative of the function at stp.
c
c       stp is a double precision variable. 
c         On entry stp is the current estimate of a satisfactory 
c            step. On initial entry, a positive initial estimate 
c            must be provided. 
c         On exit stp is the current estimate of a satisfactory step
c            if task = 'FG'. If task = 'CONV' then stp satisfies
c            the sufficient decrease and curvature condition.
c
c       ftol is a double precision variable.
c         On entry ftol specifies a nonnegative tolerance for the 
c            sufficient decrease condition.
c         On exit ftol is unchanged.
c
c       gtol is a double precision variable.
c         On entry gtol specifies a nonnegative tolerance for the 
c            curvature condition. 
c         On exit gtol is unchanged.
c
c       xtol is a double precision variable.
c         On entry xtol specifies a nonnegative relative tolerance
c            for an acceptable step. The subroutine exits with a
c            warning if the relative difference between sty and stx
c            is less than xtol.
c         On exit xtol is unchanged.
c
c       stpmin is a double precision variable.
c         On entry stpmin is a nonnegative lower bound for the step.
c         On exit stpmin is unchanged.
c
c       stpmax is a double precision variable.
c         On entry stpmax is a nonnegative upper bound for the step.
c         On exit stpmax is unchanged.
c
c       task is a character variable of length at least 60.
c         On initial entry task must be set to 'START'.
c         On exit task indicates the required action:
c
c            If task(1:2) = 'FG' then evaluate the function and 
c            derivative at stp and call dcsrch again.
c
c            If task(1:4) = 'CONV' then the search is successful.
c
c            If task(1:4) = 'WARN' then the subroutine is not able
c            to satisfy the convergence conditions. The exit value of
c            stp contains the best point found during the search.
c
c            If task(1:5) = 'ERROR' then there is an error in the
c            input arguments.
c
c         On exit with convergence, a warning or an error, the
c            variable task contains additional information.
c
c       isave is an integer work array of dimension 2.
c         
c       dsave is a double precision work array of dimension 13.
c
c     Subprograms called
c
c       MINPACK-2 ... dcstep
c
c     MINPACK-1 Project. June 1983.
c     Argonne National Laboratory. 
c     Jorge J. More' and David J. Thuente.
c
c     MINPACK-2 Project. October 1993.
c     Argonne National Laboratory and University of Minnesota. 
c     Brett M. Averick, Richard G. Carter, and Jorge J. More'. 
c
c     **********
      double precision zero,p5,p66
      parameter(zero=0.0d0,p5=0.5d0,p66=0.66d0)
      double precision xtrapl,xtrapu
      parameter(xtrapl=1.1d0,xtrapu=4.0d0)

      logical brackt
      integer stage
      double precision finit,ftest,fm,fx,fxm,fy,fym,ginit,gtest,
     +       gm,gx,gxm,gy,gym,stx,sty,stmin,stmax,width,width1

c     Initialization block.

      if (task(1:5) .eq. 'START') then

c        Check the input arguments for errors.

         if (stp .lt. stpmin) task = 'ERROR: STP .LT. STPMIN'
         if (stp .gt. stpmax) task = 'ERROR: STP .GT. STPMAX'
         if (g .ge. zero) task = 'ERROR: INITIAL G .GE. ZERO'
         if (ftol .lt. zero) task = 'ERROR: FTOL .LT. ZERO'
         if (gtol .lt. zero) task = 'ERROR: GTOL .LT. ZERO'
         if (xtol .lt. zero) task = 'ERROR: XTOL .LT. ZERO'
         if (stpmin .lt. zero) task = 'ERROR: STPMIN .LT. ZERO'
         if (stpmax .lt. stpmin) task = 'ERROR: STPMAX .LT. STPMIN'

c        Exit if there are errors on input.

         if (task(1:5) .eq. 'ERROR') return

c        Initialize local variables.

         brackt = .false.
         stage = 1
         finit = f
         ginit = g
         gtest = ftol*ginit
         width = stpmax - stpmin
         width1 = width/p5

c        The variables stx, fx, gx contain the values of the step, 
c        function, and derivative at the best step. 
c        The variables sty, fy, gy contain the value of the step, 
c        function, and derivative at sty.
c        The variables stp, f, g contain the values of the step, 
c        function, and derivative at stp.

         stx = zero
         fx = finit
         gx = ginit
         sty = zero
         fy = finit
         gy = ginit
         stmin = zero
         stmax = stp + xtrapu*stp
         task = 'FG'

         goto 1000

      else

c        Restore local variables.

         if (isave(1) .eq. 1) then
            brackt = .true.
         else
            brackt = .false.
         endif
         stage = isave(2) 
         ginit = dsave(1) 
         gtest = dsave(2) 
         gx = dsave(3) 
         gy = dsave(4) 
         finit = dsave(5) 
         fx = dsave(6) 
         fy = dsave(7) 
         stx = dsave(8) 
         sty = dsave(9) 
         stmin = dsave(10) 
         stmax = dsave(11) 
         width = dsave(12) 
         width1 = dsave(13) 

      endif

c     If psi(stp) <= 0 and f'(stp) >= 0 for some step, then the
c     algorithm enters the second stage.

      ftest = finit + stp*gtest
      if (stage .eq. 1 .and. f .le. ftest .and. g .ge. zero) 
     +   stage = 2

c     Test for warnings.

      if (brackt .and. (stp .le. stmin .or. stp .ge. stmax))
     +   task = 'WARNING: ROUNDING ERRORS PREVENT PROGRESS'
      if (brackt .and. stmax - stmin .le. xtol*stmax) 
     +   task = 'WARNING: XTOL TEST SATISFIED'
      if (stp .eq. stpmax .and. f .le. ftest .and. g .le. gtest) 
     +   task = 'WARNING: STP = STPMAX'
      if (stp .eq. stpmin .and. (f .gt. ftest .or. g .ge. gtest)) 
     +   task = 'WARNING: STP = STPMIN'

c     Test for convergence.

      if (f .le. ftest .and. abs(g) .le. gtol*(-ginit)) 
     +   task = 'CONVERGENCE'

c     Test for termination.

      if (task(1:4) .eq. 'WARN' .or. task(1:4) .eq. 'CONV') goto 1000

c     A modified function is used to predict the step during the
c     first stage if a lower function value has been obtained but 
c     the decrease is not sufficient.

      if (stage .eq. 1 .and. f .le. fx .and. f .gt. ftest) then

c        Define the modified function and derivative values.

         fm = f - stp*gtest
         fxm = fx - stx*gtest
         fym = fy - sty*gtest
         gm = g - gtest
         gxm = gx - gtest
         gym = gy - gtest

c        Call dcstep to update stx, sty, and to compute the new step.

         call dcstep(stx,fxm,gxm,sty,fym,gym,stp,fm,gm,
     +               brackt,stmin,stmax)

c        Reset the function and derivative values for f.

         fx = fxm + stx*gtest
         fy = fym + sty*gtest
         gx = gxm + gtest
         gy = gym + gtest

      else

c       Call dcstep to update stx, sty, and to compute the new step.

        call dcstep(stx,fx,gx,sty,fy,gy,stp,f,g,
     +              brackt,stmin,stmax)

      endif

c     Decide if a bisection step is needed.

      if (brackt) then
         if (abs(sty-stx) .ge. p66*width1) stp = stx + p5*(sty - stx)
         width1 = width
         width = abs(sty-stx)
      endif

c     Set the minimum and maximum steps allowed for stp.

      if (brackt) then
         stmin = min(stx,sty)
         stmax = max(stx,sty)
      else
         stmin = stp + xtrapl*(stp - stx)
         stmax = stp + xtrapu*(stp - stx)
      endif
 
c     Force the step to be within the bounds stpmax and stpmin.
 
      stp = max(stp,stpmin)
      stp = min(stp,stpmax)

c     If further progress is not possible, let stp be the best
c     point obtained during the search.

      if (brackt .and. (stp .le. stmin .or. stp .ge. stmax)
     +   .or. (brackt .and. stmax-stmin .le. xtol*stmax)) stp = stx

c     Obtain another function and derivative.

      task = 'FG'

 1000 continue

c     Save local variables.

      if (brackt) then
         isave(1) = 1
      else
         isave(1) = 0
      endif
      isave(2) = stage
      dsave(1) =  ginit
      dsave(2) =  gtest
      dsave(3) =  gx
      dsave(4) =  gy
      dsave(5) =  finit
      dsave(6) =  fx
      dsave(7) =  fy
      dsave(8) =  stx
      dsave(9) =  sty
      dsave(10) = stmin
      dsave(11) = stmax
      dsave(12) = width
      dsave(13) = width1

      return
      end
      
c====================== The end of dcsrch ==============================

      subroutine dcstep(stx,fx,dx,sty,fy,dy,stp,fp,dp,brackt,
     +                  stpmin,stpmax)
      logical brackt
      double precision stx,fx,dx,sty,fy,dy,stp,fp,dp,stpmin,stpmax
c     **********
c
c     Subroutine dcstep
c
c     This subroutine computes a safeguarded step for a search
c     procedure and updates an interval that contains a step that
c     satisfies a sufficient decrease and a curvature condition.
c
c     The parameter stx contains the step with the least function
c     value. If brackt is set to .true. then a minimizer has
c     been bracketed in an interval with endpoints stx and sty.
c     The parameter stp contains the current step. 
c     The subroutine assumes that if brackt is set to .true. then
c
c           min(stx,sty) < stp < max(stx,sty),
c
c     and that the derivative at stx is negative in the direction 
c     of the step.
c
c     The subroutine statement is
c
c       subroutine dcstep(stx,fx,dx,sty,fy,dy,stp,fp,dp,brackt,
c                         stpmin,stpmax)
c
c     where
c
c       stx is a double precision variable.
c         On entry stx is the best step obtained so far and is an
c            endpoint of the interval that contains the minimizer. 
c         On exit stx is the updated best step.
c
c       fx is a double precision variable.
c         On entry fx is the function at stx.
c         On exit fx is the function at stx.
c
c       dx is a double precision variable.
c         On entry dx is the derivative of the function at 
c            stx. The derivative must be negative in the direction of 
c            the step, that is, dx and stp - stx must have opposite 
c            signs.
c         On exit dx is the derivative of the function at stx.
c
c       sty is a double precision variable.
c         On entry sty is the second endpoint of the interval that 
c            contains the minimizer.
c         On exit sty is the updated endpoint of the interval that 
c            contains the minimizer.
c
c       fy is a double precision variable.
c         On entry fy is the function at sty.
c         On exit fy is the function at sty.
c
c       dy is a double precision variable.
c         On entry dy is the derivative of the function at sty.
c         On exit dy is the derivative of the function at the exit sty.
c
c       stp is a double precision variable.
c         On entry stp is the current step. If brackt is set to .true.
c            then on input stp must be between stx and sty. 
c         On exit stp is a new trial step.
c
c       fp is a double precision variable.
c         On entry fp is the function at stp
c         On exit fp is unchanged.
c
c       dp is a double precision variable.
c         On entry dp is the the derivative of the function at stp.
c         On exit dp is unchanged.
c
c       brackt is an logical variable.
c         On entry brackt specifies if a minimizer has been bracketed.
c            Initially brackt must be set to .false.
c         On exit brackt specifies if a minimizer has been bracketed.
c            When a minimizer is bracketed brackt is set to .true.
c
c       stpmin is a double precision variable.
c         On entry stpmin is a lower bound for the step.
c         On exit stpmin is unchanged.
c
c       stpmax is a double precision variable.
c         On entry stpmax is an upper bound for the step.
c         On exit stpmax is unchanged.
c
c     MINPACK-1 Project. June 1983
c     Argonne National Laboratory. 
c     Jorge J. More' and David J. Thuente.
c
c     MINPACK-2 Project. October 1993.
c     Argonne National Laboratory and University of Minnesota. 
c     Brett M. Averick and Jorge J. More'.
c
c     **********
      double precision zero,p66,two,three
      parameter(zero=0.0d0,p66=0.66d0,two=2.0d0,three=3.0d0)
   
      double precision gamma,p,q,r,s,sgnd,stpc,stpf,stpq,theta

      sgnd = dp*(dx/abs(dx))

c     First case: A higher function value. The minimum is bracketed. 
c     If the cubic step is closer to stx than the quadratic step, the 
c     cubic step is taken, otherwise the average of the cubic and 
c     quadratic steps is taken.

      if (fp .gt. fx) then
         theta = three*(fx - fp)/(stp - stx) + dx + dp
         s = max(abs(theta),abs(dx),abs(dp))
         gamma = s*sqrt((theta/s)**2 - (dx/s)*(dp/s))
         if (stp .lt. stx) gamma = -gamma
         p = (gamma - dx) + theta
         q = ((gamma - dx) + gamma) + dp
         r = p/q
         stpc = stx + r*(stp - stx)
         stpq = stx + ((dx/((fx - fp)/(stp - stx) + dx))/two)*
     +                                                       (stp - stx)
         if (abs(stpc-stx) .lt. abs(stpq-stx)) then
            stpf = stpc
         else
            stpf = stpc + (stpq - stpc)/two
         endif
         brackt = .true.

c     Second case: A lower function value and derivatives of opposite 
c     sign. The minimum is bracketed. If the cubic step is farther from 
c     stp than the secant step, the cubic step is taken, otherwise the 
c     secant step is taken.

      else if (sgnd .lt. zero) then
         theta = three*(fx - fp)/(stp - stx) + dx + dp
         s = max(abs(theta),abs(dx),abs(dp))
         gamma = s*sqrt((theta/s)**2 - (dx/s)*(dp/s))
         if (stp .gt. stx) gamma = -gamma
         p = (gamma - dp) + theta
         q = ((gamma - dp) + gamma) + dx
         r = p/q
         stpc = stp + r*(stx - stp)
         stpq = stp + (dp/(dp - dx))*(stx - stp)
         if (abs(stpc-stp) .gt. abs(stpq-stp)) then
            stpf = stpc
         else
            stpf = stpq
         endif
         brackt = .true.

c     Third case: A lower function value, derivatives of the same sign,
c     and the magnitude of the derivative decreases.

      else if (abs(dp) .lt. abs(dx)) then

c        The cubic step is computed only if the cubic tends to infinity 
c        in the direction of the step or if the minimum of the cubic
c        is beyond stp. Otherwise the cubic step is defined to be the 
c        secant step.

         theta = three*(fx - fp)/(stp - stx) + dx + dp
         s = max(abs(theta),abs(dx),abs(dp))

c        The case gamma = 0 only arises if the cubic does not tend
c        to infinity in the direction of the step.

         gamma = s*sqrt(max(zero,(theta/s)**2-(dx/s)*(dp/s)))
         if (stp .gt. stx) gamma = -gamma
         p = (gamma - dp) + theta
         q = (gamma + (dx - dp)) + gamma
         r = p/q
         if (r .lt. zero .and. gamma .ne. zero) then
            stpc = stp + r*(stx - stp)
         else if (stp .gt. stx) then
            stpc = stpmax
         else
            stpc = stpmin
         endif
         stpq = stp + (dp/(dp - dx))*(stx - stp)

         if (brackt) then

c           A minimizer has been bracketed. If the cubic step is 
c           closer to stp than the secant step, the cubic step is 
c           taken, otherwise the secant step is taken.

            if (abs(stpc-stp) .lt. abs(stpq-stp)) then
               stpf = stpc
            else
               stpf = stpq
            endif
            if (stp .gt. stx) then
               stpf = min(stp+p66*(sty-stp),stpf)
            else
               stpf = max(stp+p66*(sty-stp),stpf)
            endif
         else

c           A minimizer has not been bracketed. If the cubic step is 
c           farther from stp than the secant step, the cubic step is 
c           taken, otherwise the secant step is taken.

            if (abs(stpc-stp) .gt. abs(stpq-stp)) then
               stpf = stpc
            else
               stpf = stpq
            endif
            stpf = min(stpmax,stpf)
            stpf = max(stpmin,stpf)
         endif

c     Fourth case: A lower function value, derivatives of the same sign, 
c     and the magnitude of the derivative does not decrease. If the 
c     minimum is not bracketed, the step is either stpmin or stpmax, 
c     otherwise the cubic step is taken.

      else
         if (brackt) then
            theta = three*(fp - fy)/(sty - stp) + dy + dp
            s = max(abs(theta),abs(dy),abs(dp))
            gamma = s*sqrt((theta/s)**2 - (dy/s)*(dp/s))
            if (stp .gt. sty) gamma = -gamma
            p = (gamma - dp) + theta
            q = ((gamma - dp) + gamma) + dy
            r = p/q
            stpc = stp + r*(sty - stp)
            stpf = stpc
         else if (stp .gt. stx) then
            stpf = stpmax
         else
            stpf = stpmin
         endif
      endif

c     Update the interval which contains a minimizer.

      if (fp .gt. fx) then
         sty = stp
         fy = fp
         dy = dp
      else
         if (sgnd .lt. zero) then
            sty = stx
            fy = fx
            dy = dx
         endif
         stx = stp
         fx = fp
         dx = dp
      endif

c     Compute the new step.

      stp = stpf

      return
      end
      
