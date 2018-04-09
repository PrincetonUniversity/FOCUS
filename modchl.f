C      ALGORITHM 695 , COLLECTED ALGORITHMS FROM ACM.
C      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
C      VOL. 17, NO. 3, SEPTEMBER, 1991, PP. 306-312.
c

C*********************************************************************
C
C       subroutine name: modchl
C
C       authors :  Elizabeth Eskow and Robert B. Schnabel
C
C       date    : December, 1988
C
C       purpose : perform a modified cholesky factorization
C                 of the form (Ptranspose)AP  + E = L(Ltranspose),
C       where L is stored in the lower triangle of the
C       original matrix A.
C       The factorization has 2 phases:
C        phase 1: Pivot on the maximum diagonal element.
C            Check that the normal cholesky update
C            would result in a positive diagonal
C            at the current iteration, and
C            if so, do the normal cholesky update,
C            otherwise switch to phase 2.
C        phase 2: Pivot on the minimum of the negatives
C            of the lower gerschgorin bound
C            estimates.
C            Compute the amount to add to the
C            pivot element and add this
C            to the pivot element.
C            Do the cholesky update.
C            Update the estimates of the
C            gerschgorin bounds.
C
C       input   : ndim    - largest dimension of matrix that
C                           will be used
C
C                 n       - dimension of matrix A
C
C                 A       - n*n symmetric matrix (only lower triangular
C            portion of A, including the main diagonal, is used)
C
C                 g       - n*1 work array
C
C                 mcheps - machine precision
C
C                tau1    - tolerance used for determining when to switch
C                          phase 2
C
C                tau2    - tolerance used for determining the maximum
C                          condition number of the final 2X2 submatrix.
C
C
C       output  : L     - stored in the matrix A (in lower triangular
C                           portion of A, including the main diagonal)
C
C                 P     - a record of how the rows and columns
C                         of the matrix were permuted while
C                         performing the decomposition
C
C                 E     - n*1 array, the ith element is the
C                         amount added to the diagonal of A
C                         at the ith iteration
C
C
C***********************************************************************
      subroutine modchl(ndim,n,A,g,mcheps,tau1,tau2,P,E)
*
      integer n,ndim
      double precision A(ndim,n),g(n),mcheps,tau1,tau2
      integer P(n)
      double precision E(n)
*
C
C  j        - current iteration number
C  iming    - index of the row with the min. of the
C           neg. lower Gersch. bounds
C  imaxd    - index of the row with the maximum diag.
C           element
C  i,itemp,jpl,k  - temporary integer variables
C  delta    - amount to add to Ajj at the jth iteration
C  gamma    - the maximum diagonal element of the original
C           matrix A.
C  normj    - the 1 norm of A(colj), rows j+1 --> n.
C  ming     - the minimum of the neg. lower Gersch. bounds
C  maxd     - the maximum diagonal element
C  taugam - tau1 * gamma
C  phase1      - logical, true if in phase1, otherwise false
C  delta1,temp,jdmin,tdmin,tempjj - temporary double precision vars.
C
*
      integer j,iming,i,imaxd,itemp,jp1,k
      double precision delta,gamma
      double precision normj, ming,maxd
      double precision delta1,temp,jdmin,tdmin,taugam,tempjj
      logical phase1
      intrinsic abs, max, sqrt, min
*
      call init(n, ndim, A, phase1, delta, P, g, E,
     *         ming,tau1,gamma,taugam)
C
C     check for n=1
C
      if (n.eq.1) then
         delta = (tau2 * abs(A(1,1))) - A(1,1)
         if (delta .gt. 0) E(1) = delta
         if (A(1,1) .eq. 0) E(1) = tau2
         A(1,1)=sqrt(A(1,1)+E(1))
      endif
C
      do 200 j = 1, n-1
C
C        PHASE 1
C
         if ( phase1 ) then
C
C           Find index of maximum diagonal element A(i,i) where i>=j
C
            maxd = A(j,j)
            imaxd = j
            do 20 i = j+1, n
               if (maxd .lt. A(i,i)) then
                  maxd = A(i,i)
                  imaxd = i
               end if
 20         continue
*
C
C           Pivot to the top the row and column with the max diag
C
            if (imaxd .ne. j) then
C
C              Swap row j with row of max diag
C
               do 30 i = 1, j-1
                  temp = A(j,i)
                  A(j,i) = A(imaxd,i)
                  A(imaxd,i) = temp
 30            continue
C
C              Swap colj and row maxdiag between j and maxdiag
C
               do 35 i = j+1,imaxd-1
                  temp = A(i,j)
                  A(i,j) = A(imaxd,i)
                  A(imaxd,i) = temp
 35            continue
C
C              Swap column j with column of max diag
C
               do 40 i = imaxd+1, n
                  temp = A(i,j)
                  A(i,j) = A(i,imaxd)
                  A(i,imaxd) = temp
 40            continue
C
C              Swap diag elements
C
               temp = A(j,j)
               A(j,j) = A(imaxd,imaxd)
               A(imaxd,imaxd) = temp
C
C              Swap elements of the permutation vector
C
               itemp = P(j)
               P(j) = P(imaxd)
               P(imaxd) = itemp
*
            end if
*
*
C           Check to see whether the normal cholesky update for this
C           iteration would result in a positive diagonal,
C           and if not then switch to phase 2.
*
            jp1 = j+1
            tempjj=A(j,j)
*
            if (tempjj.gt.0) then
*
               jdmin=A(jp1,jp1)
               do 60 i = jp1, n
                  temp = A(i,j) * A(i,j) / tempjj
                  tdmin = A(i,i) - temp
                  jdmin = min(jdmin, tdmin)
 60            continue
*
               if (jdmin .lt. taugam) phase1 = .false.
*
            else
*
               phase1 = .false.
*
            end if
*
            if (phase1) then
C
C              Do the normal cholesky update if still in phase 1
C
               A(j,j) = sqrt(A(j,j))
               tempjj = A(j,j)
               do 70 i = jp1, n
                  A(i,j) = A(i,j) / tempjj
 70            continue
               do 80 i=jp1,n
                  temp=A(i,j)
                  do 75 k = jp1, i
                     A(i,k) = A(i,k) - (temp * A(k,j))
 75               continue
 80            continue
*
               if (j .eq. n-1) A(n,n)=sqrt(A(n,n))
*
            else
*
C
C              Calculate the negatives of the lower gerschgorin bounds
C
               call gersch(ndim,n,A,j,g)
*
            end if
*
         end if
*
*
C
C        PHASE 2
C
         if (.not. phase1) then
*
            if (j .ne. n-1) then
C
C              Find the minimum negative gershgorin bound
C
*
               iming=j
               ming = g(j)
               do 90 i = j+1,n
                  if (ming .gt. g(i)) then
                     ming = g(i)
                     iming = i
                  end if
 90            continue
*
C
C               Pivot to the top the row and column with the
C               minimum negative gerschgorin bound
C
                if (iming .ne. j) then
C
C                  Swap row j with row of min gersch bound
C
                   do 100 i = 1, j-1
                      temp = A(j,i)
                       A(j,i) = A(iming,i)
                       A(iming,i) = temp
 100               continue
C
C                  Swap colj with row iming from j to iming
C
                   do 105 i = j+1,iming-1
                      temp = A(i,j)
                      A(i,j) = A(iming,i)
                      A(iming,i) = temp
 105              continue
C
C                 Swap column j with column of min gersch bound
C
                  do 110 i = iming+1, n
                     temp = A(i,j)
                     A(i,j) = A(i,iming)
                     A(i,iming) = temp
 110              continue
C
C                 Swap diagonal elements
C
                  temp = A(j,j)
                  A(j,j) = A(iming,iming)
                  A(iming,iming) = temp
C
C                 Swap elements of the permutation vector
C
                  itemp = P(j)
                  P(j) = P(iming)
                  P(iming) = itemp
C
C                 Swap elements of the negative gerschgorin bounds vecto
C
                  temp = g(j)
                  g(j) = g(iming)
                  g(iming) = temp
*
               end if
C
C              Calculate delta and add to the diagonal.
C              delta=max{0,-A(j,j) + max{normj,taugam},delta_previous}
C              where normj=sum of |A(i,j)|,for i=1,n,
C              delta_previous is the delta computed at the previous iter
C              and taugam is tau1*gamma.
C
*
               normj = 0.0
               do 140 i = j+1, n
                  normj = normj + abs(A(i,j))
 140           continue
*
               temp = max(normj,taugam)
               delta1 = temp - A(j,j)
               temp = 0.0
               delta1 = max(temp, delta1)
               delta = max(delta1,delta)
               E(j) =  delta
               A(j,j) = A(j,j) + E(j)
C
C              Update the gerschgorin bound estimates
C              (note: g(i) is the negative of the
C               Gerschgorin lower bound.)
C
               if (A(j,j) .ne. normj) then
                  temp = (normj/A(j,j)) - 1.0
*
                  do 150 i = j+1, n
                     g(i) = g(i) + abs(A(i,j)) * temp
 150              continue
*
               end if
C
C              Do the cholesky update
C
               A(j,j) = sqrt(A(j,j))
               tempjj = A(j,j)
               do 160 i = j+1, n
                  A(i,j) = A(i,j) / tempjj
 160           continue
               do 180 i = j+1, n
                  temp = A(i,j)
                  do 170 k = j+1, i
                     A(i,k) = A(i,k) - (temp * A(k,j))
 170              continue
 180           continue
*
            else
*
               call fin2x2(ndim, n, A, E, j, tau2, delta,gamma)
*
            end if
*
         end if
*
 200   continue
*
      return
      end
C***********************************************************************
C       subroutine name : init
C
C       purpose : set up for start of cholesky factorization
C
C       input : n, ndim, A, tau1
C
C       output : phase1    - boolean value set to true if in phase one,
C             otherwise false.
C      delta     - amount to add to Ajj at iteration j
C      P,g,E - described above in modchl
C      ming      - the minimum negative gerschgorin bound
C      gamma     - the maximum diagonal element of A
C      taugam  - tau1 * gamma
C
C***********************************************************************
      subroutine init(n,ndim,A,phase1,delta,P,g,E,ming,
     *                tau1,gamma,taugam)
*
      integer n,ndim
      double precision A(ndim,n)
      logical phase1
      double precision delta,g(n),E(n)
      integer P(n)
      double precision ming,tau1,gamma,taugam
      intrinsic abs, max
*
*
      phase1 = .true.
      delta = 0.0
      ming = 0.0
      do 10 i=1,n
         P(i)=i
         g(i)= 0.0
         E(i) = 0.0
 10   continue
*
C
C     Find the maximum magnitude of the diagonal elements.
C     If any diagonal element is negative, then phase1 is false.
C
      gamma = 0.0
      do 20 i=1,n
         gamma=max(gamma,abs(A(i,i)))
         if (A(i,i) .lt. 0.0) phase1 = .false.
 20   continue
*
      taugam = tau1 * gamma
*
C
C     If not in phase1, then calculate the initial gerschgorin bounds
C     needed for the start of phase2.
C
      if ( .not.phase1) call gersch(ndim,n,A,1,g)
*
      return
      end
C***********************************************************************
C
C       subroutine name : gersch
C
C       purpose : Calculate the negative of the gerschgorin bounds
C                 called once at the start of phase II.
C
C       input   : ndim, n, A, j
C
C       output  : g - an n vector containing the negatives of the
C           Gerschgorin bounds.
C
C***********************************************************************
      subroutine gersch(ndim, n, A, j, g)
*
      integer ndim, n, j
      double precision A(ndim,n), g(n)
*
      integer i, k
      double precision offrow
      intrinsic abs
*
      do 30 i = j, n
         offrow = 0.0
         do 10 k = j, i-1
            offrow = offrow + abs(A(i,k))
 10      continue
         do 20 k = i+1, n
            offrow = offrow + abs(A(k,i))
 20      continue
            g(i) = offrow - A(i,i)
 30   continue
*
      return
      end
C***********************************************************************
C
C  subroutine name : fin2x2
C
C  purpose : Handles final 2X2 submatrix in Phase II.
C            Finds eigenvalues of final 2 by 2 submatrix,
C            calculates the amount to add to the diagonal,
C            adds to the final 2 diagonal elements,
C            and does the final update.
C
C  input : ndim, n, A, E, j, tau2,
C          delta - amount added to the diagonal in the
C                  previous iteration
C
C  output : A - matrix with complete L factor in the lower triangle,
C           E - n*1 vector containing the amount added to the diagonal
C               at each iteration,
C           delta - amount added to diagonal elements n-1 and n.
C
C***********************************************************************
      subroutine fin2x2(ndim, n, A, E, j, tau2, delta,gamma)
*
      integer ndim, n, j
      double precision A(ndim,n), E(n), tau2, delta,gamma
*
      double precision t1, t2, t3,lmbd1,lmbd2,lmbdhi,lmbdlo
      double precision delta1, temp
      intrinsic sqrt, max, min
*
C
C     Find eigenvalues of final 2 by 2 submatrix
C
      t1 = A(n-1,n-1) + A(n,n)
      t2 = A(n-1,n-1) - A(n,n)
      t3 = sqrt(t2*t2 + 4.0*A(n,n-1)*A(n,n-1))
      lmbd1 = (t1 - t3)/2.
      lmbd2 = (t1 + t3)/2.
      lmbdhi = max(lmbd1,lmbd2)
      lmbdlo = min(lmbd1,lmbd2)
C
C     Find delta such that:
C     1.  the l2 condition number of the final
C     2X2 submatrix + delta*I <= tau2
C     2. delta >= previous delta,
C     3. lmbdlo + delta >= tau2 * gamma,
C     where lmbdlo is the smallest eigenvalue of the final
C     2X2 submatrix
C
*
      delta1=(lmbdhi-lmbdlo)/(1.0-tau2)
      delta1= max(delta1,gamma)
      delta1= tau2 * delta1 - lmbdlo
      temp = 0.0
      delta = max(delta, temp)
      delta = max(delta, delta1)
*
      if (delta .gt. 0.0) then
         A(n-1,n-1) = A(n-1,n-1) + delta
         A(n,n) = A(n,n) + delta
         E(n-1) = delta
         E(n) = delta
      end if
C
C     Final update
C
      A(n-1,n-1) = sqrt(A(n-1,n-1))
      A(n,n-1) = A(n,n-1)/A(n-1,n-1)
      A(n,n) = A(n,n) - (A(n,n-1)**2)
      A(n,n) = sqrt(A(n,n))
*
      return
      end

