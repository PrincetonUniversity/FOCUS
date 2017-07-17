C ----------------------------------------------------------
C SEGMENT 4:
C 4. TNPACK routines (new version)
c----------------------------------------------------------
c
C***********************************************************************
C***********************************************************************
C                                                                      *
C                  -    "TNPACK"   -                                   *
C                                                                      *
C          A  PACKAGE  OF  ROUTINES  FOR  A                            *
C       LARGE-SCALE TRUNCATED NEWTON ALGORITHM                         *
C                                                                      *
C                copyright   (c)   1990  by                            *
C                                                                      *
C            Tamar Schlick   &    Aaron Fogelson                       *
c                                                                      *
c                 Updated (November 1998) by                           *
c                                                                      *
c                 Dexuan Xie & Tamar Schlick                           *
C                                                                      *
C  The changes made involve the following components:                  *
C                                                                      *
C    (a) Modified Cholesky factorization                               *
C    (b) Negative curvature test                                       *
C    (c) Line search scheme                                            *
C                                                                      *
C  For details, see                                                    *
C                                                                      *
C    (1) Xie, D. and Schlick, T.  Remark on the Updated Truncated      *
C        Newton Minimization Package, Algorithm 702, ACM Trans. Math   *
C        Softw., in Press  (1998)                                      *
C                                                                      *
C    (2) Xie, D. and Schlick, T.  Efficient implementation of the      *
C        truncated-Newton algorithm for large-scale chemistry          *
C        applications, SIAM J. Opt., in Press (1998)                   *
C                                                                      *
C  The manuscripts are also available from the authors. Please send    *
C  requests by email to                                                *
C                                                                      *
C      dexuan@cims.nyu.edu    or   schlick@nyu.edu                     *
C                                                                      *
C***********************************************************************
C                                                                      *
C PLEASE NOTE:                                                         *
C                                                                      *
C (i)   Double Precision is used througout.                            *
C (ii)  Documentation for the program is given separately; here        *
C       only brief comments are provided.                              *
C (iii) It is essential that the user become familiar with all input   *
C       options and parameters before using this package. Entry to     *
C       TNPACK from the user's driver must first be done through a     *
C       call to routine SETLIS. SETLIS sets sample input options and   *
C       parameters, some of which the user may then modify to suit     *
C       the problem. Routine TNMIN should then be called.              *
C (iv)  For runs on Cray computers, the Yale Sparse Matrix Package     *
C       parameter RATIO must be changed from 2 to 1 - routine SDRVMD   *
C (v)   If running on VAX/VMS systems with built-in facilities         *
C       for BLAS, you may need to compile this program with something  *
C       like FOR/NOBLAS, since some standard BLAS are included here.   *
C                                                                      *
C***********************************************************************
C                                                                      *
C TNPACK COMPONENTS:                                                   *
C                                                                      *
C (A)  USER-SUPPLIED SUBROUTINES (generic names):                      *
C      -----------------------------------------                       *
C                                                                      *
C  *   CALFGH    - calculates F,G,H and M (preconditioner) at X        *
C  *   CALHDP    - calculates HD, a Hessian/vector product             *
C  *   CALPAT    - determines the sparsity pattern of M                *
C                                                                      *
C (B) TRUNCATED NEWTON SUBROUTINES:                                    *
C     ----------------------------                                     *
C                                                                      *
C  *   TNMIN              -    the TN driver routine                   *
C  *   OUTER              -    the outer Newton iteration              *
C  *   INNER              -    the inner PCG iteration                 *
C  *   OURHDP             -    computes Hessian/vector products        *
C  *   SETLIS             -    sets sample parameters for minimization *
C  *   CHKLIS             -    checks input options and parameters     *
C  *   DIVWRK             -    divides work space                      *
C  *   MLINES, NEWSTP     -    perform the line search                 *
C                                                                      *
C (C)  YSMP (YALE SPARSE MATRIX PACKAGE) SUBROUTINES:                  *
C      ---------------------------------------------                   *
C                                                                      *
C  *   SDRVMD             -    the driver for solving Mz=r             *
C  *   ODRV               -    the driver for finding M's reordering   *
C  *   MD,MDI,MDM,MDP,MDU -    the minimum-degree reordering routines  *
C  *   SRO                -    prepares for symmetric reordering of M  *
C  *   SSF                -    performs symbolic factorization of M    *
C  *   SNFMOD             -    performs numerical factorization of M   *
C  *   SNS                -    performs numerical solution of Mz=r     *
C                                                                      *
C (D)  OTHER SUBPROGRAMS:                                              *
C      ----------------                                                *
C                                                                      *
C  *   DCOPY, DAXPY         (blas subroutines)                         *
C  *   DNRM2, DDOT, DMACH   (blas functions)                           *
C                                                                      *
C***********************************************************************
C***********************************************************************
      SUBROUTINE TNMIN(N,X,F,G,OPLIST,PLIST,INFORM,NZ,W,LW,IW,LIW,
     +                 CALFGH,CALPAT,CALHDP)

C --------------------------------------------------------
C TNMIN:  Truncated Newton Interface Routine
C --------------------------------------------------------

C     .. Scalar Arguments ..
      DOUBLE PRECISION F
      INTEGER INFORM,LIW,LW,N,NZ
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION G(N),PLIST(20),W(LW),X(N)
      INTEGER IW(LIW),OPLIST(20)
C     ..
C     .. Subroutine Arguments ..
      EXTERNAL CALFGH,CALHDP,CALPAT
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION TAU
      INTEGER ICALLS,IXA,IXD,IXD1,IXHD,IXIA,IXIPER,IXISP,IXJA,IXP,IXPER,
     +        IXRSP,IXY,IXZ,MC,MP,NPROB,SRLS
C     ..
C     .. Local Scalars ..
      INTEGER I,MARK,NSP
C     ..
C     .. External Subroutines ..
      EXTERNAL CHKLIS,DIVWRK,OUTER
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MAX
C     ..
C     .. Common blocks ..
      COMMON NPROB
      COMMON /TAUINPUT/TAU,SRLS,MC
      COMMON /TN001/MP,ICALLS
      COMMON /TN004/IXP,IXY,IXD,IXHD,IXZ,IXA,IXRSP,IXD1
      COMMON /TN005/IXPER,IXIPER,IXIA,IXJA,IXISP
C     ..
C     .. Save statement ..

C --------------------------------------------------------
C check parameters and divide work space
C --------------------------------------------------------

c
c--- INFORM: indicates the status of TNPACK
c

      SAVE /TN001/,/TN004/,/TN005/
C     ..
      IF (INFORM.NE.-10) THEN
          WRITE (*,FMT=9000)
          INFORM = -1
          RETURN

      ELSE
          INFORM = 0
      END IF

      MARK = 0
c
c--- Check input options & parameters (OPLIST,PLIST)
c
c    INFORM = -2: input errors in N, OPLIST,PLIST
c
      CALL CHKLIS(N,OPLIST,PLIST,MARK)
      IF (MARK.LT.0) THEN
          INFORM = -2
          RETURN

      END IF
c
c--- Specifies the unit number for printing
c
      MP = OPLIST(11)

      IF (ICALLS.EQ.0) THEN

c--- xie: NPROB is used only for testing Alg.566 to print out
c         the input information only once for running all 18 problems.
C    To do so, uncomment "if (NPROB .eq. 1) then"  and "endif"   below
c    Then compile  with make test3''
c---------------------------------------------------------------------
c         IF (NPROB.EQ.1) THEN
              WRITE (MP,FMT=9010) N, (OPLIST(I),I=1,6)
              WRITE (MP,FMT=9020) (OPLIST(I),I=7,17)
              WRITE (MP,FMT=9030) (PLIST(I),I=1,8)
c         END IF

      ELSE
          WRITE (MP,FMT=9040) ICALLS + 1
      END IF

      ICALLS = ICALLS + 1
c
c--- Compute starting indices for work vectors W,IW
c
      MARK = 0
      CALL DIVWRK(N,NZ,LW,LIW,MARK)
      IF (MARK.LT.0) THEN
          INFORM = -3
          RETURN

      END IF
c
c--- Outer Newton iteration of the TN method
c
      NSP = 3*N + 4*MAX(N,NZ)

      CALL OUTER(N,X,F,G,OPLIST,PLIST,INFORM,NZ,NSP,W(IXP),W(IXY),
     +           W(IXD),W(IXHD),W(IXZ),W(IXA),W(IXRSP),IW(IXPER),
     +           IW(IXIPER),IW(IXIA),IW(IXJA),IW(IXISP),CALFGH,CALPAT,
     +           CALHDP,W(IXD1))

      RETURN

 9000 FORMAT (/,2X,'SETLIS MUST BE CALLED BEFORE TNMIN! ')
 9010 FORMAT (/,'  PARAMETER INFO, ENTERING TNPACK',/,8X,' N = ',I8,/,
     +       5X,'OPLIST:',/,8X,' 1) IPCG   = ',I8,8X,
     +       '(preconditioning option)',/,8X,' 2) IPRINT = ',I8,8X,
     +       '(printing option)',/,8X,' 3) IPFREQ = ',I8,8X,
     +       '(printing frequency)',/,8X,' 4) MXITN  = ',I8,8X,
     +       '(max Newton itns.)',/,8X,' 5) MXITCG = ',I8,8X,
     +       '(max PCG itns.)',/,8X,' 6) MAXNF  = ',I8,8X,
     +       '(max F&G evals.)')
 9020 FORMAT (8X,' 7) IORDER = ',I8,8X,
     +       '(M-reordering option - calc. per.)',/,8X,' 8) IPERPR = ',
     +       I8,8X,'(M-reordering printing option)',/,8X,
     +       ' 9) IPKNOW = ',I8,8X,'(M-reordering option -',
     +       ' per. known)',/,8X,'10) IHD    = ',I8,8X,
     +       '(Numeric HD option)',/,8X,'11) MP     = ',I8,8X,
     +       '(printing unit)',/,8X,'12) IEXIT  = ',I8,8X,
     +       '(descent option for neg.',' curvature)',/,8X,
     +       '13) ITT    = ',I8,8X,'(truncation test option)',/,8X,
     +       '14) MAXFEV = ',I8,8X,'(max F evals. in line search)',/,8X,
     +       '15) TA     = ',I8,8X,
     +       '(descent direction test option by PCG)',/,8X,
     +       '16) SRLS   = ',I8,8X,'(stopping rule for line search)',/,
     +       8X,'17) MC     = ',I8,8X,'(modified Cholesky)')
 9030 FORMAT (5X,'PLIST:',/,8X,' 1) EPSF   = ',1P,E10.3,6X,
     +       '(controls conv. test, F accuracy)',/,8X,' 2) EPSG   = ',
     +       1P,E10.3,6X,'(controls conv. test, G accuracy)',/,8X,
     +       ' 3) ETAR   = ',1P,E10.3,6X,
     +       '(controls res.-based truncation)',/,8X,' 4) ETAQ   = ',1P,
     +       E10.3,6X,'(controls quad.-model based ','truncation)',/,8X,
     +       ' 5) PRODK  = ',1P,E10.3,6X,'(controls PCG neg curv. test)'
     +       ,/,8X,' 6) FTOL   = ',1P,E10.3,6X,
     +       '(controls F test in line search)',/,8X,' 7) GTOL   = ',1P,
     +       E10.3,6X,'(controls G test in line ','search)',/,8X,
     +       ' 8) TAU    = ',1P,E10.3,6X,'(parameter in UMC)')
 9040 FORMAT (/,6X,'ICALLS =',I8)
      END
C***********************************************************************
C***********************************************************************
      SUBROUTINE OUTER(N,X,F,G,OPLIST,PLIST,INFORM,NZ,NSP,P,Y,D,HD,Z,A,
     +                 RSP,PER,IPER,IA,JA,ISP,CALFGH,CALPAT,CALHDP,D1)

C --------------------------------------------------------
C OUTER:  Outer Newton iteration of the TN method
C --------------------------------------------------------

C --------------------------------------------------------
C subroutines and functions called:
C        calfgh, calpat, calhdp (user-supplied)
C        inner, mlines, odrv, sdrvmd, dnrm2
C --------------------------------------------------------


C --------------------------------------------------------
C Initialize:
C
c           *  Set the parameters
c           *  Give an initial guess  X_0
c           *  Compute E(X_0) and G_0
C --------------------------------------------------------

C     .. Scalar Arguments ..
      DOUBLE PRECISION F
      INTEGER INFORM,N,NSP,NZ
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(N+NZ),D(N),D1(N),G(N),HD(N),P(N),PLIST(20),
     +                 RSP(NSP),X(N),Y(N),Z(N)
      INTEGER IA(N+1),IPER(N),ISP(NSP),JA(N+NZ),OPLIST(20),PER(N)
C     ..
C     .. Subroutine Arguments ..
      EXTERNAL CALFGH,CALHDP,CALPAT
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION EPSMCH,ETAQ,ETAR,FTOL,GTOL,PRODK,SQEPS2,SQRTN,TAU
      INTEGER ICALLS,IEXIT,IHD,ITT,MAXFEV,MC,MP,SRLS,TA
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION EPSF,EPSF2,EPSF3,EPSG,FOLD,GNORM,LAMBDA,ONE,ONEF,
     +                 PERCNT,RNORM,XNORM,XSNORM
      INTEGER ESP,I,ILINE,IORDER,IPCG,IPERPR,IPFREQ,IPKNOW,IPRINT,ITR,
     +        ITRCG,ITRMAJ,LENA,MAXNF,MODE,MXITCG,MXITN,NFEV,NFUN,NOUT,
     +        OFLAG,OPATH,SFLAG,SPATH
      LOGICAL T1,T123,T2,T3,T4
C     ..
C     .. External Functions ..
      DOUBLE PRECISION DNRM2
      EXTERNAL DNRM2
C     ..
C     .. External Subroutines ..
      EXTERNAL INNER,MLINES,ODRV,SDRVMD
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,DBLE,MOD,SQRT
C     ..
C     .. Common blocks ..
      COMMON /TAINPUT/TA
      COMMON /TAUINPUT/TAU,SRLS,MC
      COMMON /TN001/MP,ICALLS
      COMMON /TN002/PRODK,ETAR,ETAQ,IEXIT,IHD,ITT
      COMMON /TN003/FTOL,GTOL,MAXFEV
      COMMON /TN006/EPSMCH,SQEPS2,SQRTN
C     ..
C     .. Save statement ..
      SAVE /TN001/,/TN002/,/TN003/,/TN006/
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION RDUM(1),ZDUM(1)
C     ..
      IPCG = OPLIST(1)
      IPRINT = OPLIST(2)
      IPFREQ = OPLIST(3)
      MXITN = OPLIST(4)
      MXITCG = OPLIST(5)
      MAXNF = OPLIST(6)
      IORDER = OPLIST(7)
      IPERPR = OPLIST(8)
      IPKNOW = OPLIST(9)
      IHD = OPLIST(10)
      MP = OPLIST(11)
      IEXIT = OPLIST(12)
      ITT = OPLIST(13)
      MAXFEV = OPLIST(14)

      EPSF = PLIST(1)
      EPSG = PLIST(2)
      ETAR = PLIST(3)
      ETAQ = PLIST(4)
      PRODK = PLIST(5)
      FTOL = PLIST(6)
      GTOL = PLIST(7)

c-- xie: new parameters
      TA = OPLIST(15)
      SRLS = OPLIST(16)
      MC = OPLIST(17)
      TAU = PLIST(8)
c----------------------------

      ITRMAJ = 0
      ITRCG = 0
      NFUN = 0
      EPSF2 = SQRT(EPSF)
      EPSF3 = EPSF** (1./3.)
      ONE = 1.0D0

C --------------------------------------------------------
C determine the pattern of the preconditioner M
C --------------------------------------------------------

      IF (IPCG.EQ.1) THEN
          CALL CALPAT(N,X,A,IA,JA)
          LENA = IA(N+1) - 1
          PERCNT = (DBLE(LENA)/DBLE(N* (N+1)/2))*100.0
          IF (LENA.GT. (N+NZ)) THEN
              INFORM = -4
              WRITE (MP,FMT=9000) LENA,N + NZ
              GO TO 30

          END IF

          IF (ICALLS.LE.1) WRITE (MP,FMT=9010) LENA,PERCNT
      END IF

C --------------------------------------------------------
C compute f, g, H, and M at x0
C --------------------------------------------------------

      NOUT = 2
      IF (IPCG.EQ.1) NOUT = 3

      CALL CALFGH(N,X,F,G,A,IA,JA,NOUT)

      NFUN = NFUN + 1

      IF (IPRINT.GE.1) THEN
          WRITE (MP,FMT=9020)
          WRITE (MP,FMT=9030) (X(I),I=1,N)
          IF (IPRINT.GE.2) THEN
              WRITE (MP,FMT=9040)
              WRITE (MP,FMT=9030) (G(I),I=1,N)
          END IF

      END IF

C --------------------------------------------------------
C check for convergence at the first iteration
C --------------------------------------------------------

      XNORM = DNRM2(N,X,1)/SQRTN
      GNORM = DNRM2(N,G,1)/SQRTN
      IF (GNORM.LT. (EPSG* (ONE+ABS(F)))) THEN
          WRITE (MP,FMT=9050) ITRMAJ,F,GNORM
          INFORM = 0
          WRITE (MP,FMT=9060) INFORM
          ITRCG = 0
          GO TO 30

      END IF

C --------------------------------------------------------
C when reordering M, call ODRV to: i) compute the permutation
C arrays and then prepare for symmetric reordering of M, or
C ii) just prepare for the symmetric reordering of M when the
C permutation arrays are known (IPKNOW=1). Then call SDRVMD
C to factorize M symbolically.
C --------------------------------------------------------

      IF (IPCG.EQ.1) THEN

          IF (IPKNOW.EQ.1) THEN
              OPATH = 5

          ELSE
              OPATH = 4
              DO 10 I = 1,N
                  PER(I) = I
                  IPER(I) = I
   10         CONTINUE
          END IF

          OFLAG = 0
          IF (IORDER.EQ.1) THEN
              CALL ODRV(N,IA,JA,A,PER,IPER,NSP,ISP,OPATH,OFLAG)
              IF (IPERPR.EQ.1 .AND. IPKNOW.NE.1) THEN
                  WRITE (MP,FMT=9070)
                  WRITE (MP,FMT=9080) (PER(I),I=1,N)
              END IF

          END IF

          SPATH = 4
          IF (OFLAG.EQ.0) CALL SDRVMD(N,PER,IPER,IA,JA,A,RDUM,ZDUM,NSP,
     +                                ISP,RSP,ESP,SPATH,SFLAG)
          IF (OFLAG.NE.0 .OR. SFLAG.NE.0) THEN
              WRITE (MP,FMT=9090) OFLAG,SFLAG
              GO TO 30

          END IF

      END IF

      WRITE (MP,FMT=9050) ITRMAJ,F,GNORM

C --------------------------------------------------------
C  MAIN LOOP BEGINS
C --------------------------------------------------------

   20 CONTINUE
      MODE = 0
      ITRMAJ = ITRMAJ + 1

C --------------------------------------------------------
C print x and f if specified
C --------------------------------------------------------

      IF (IPFREQ.GT.0) THEN
          IF (MOD(ITRMAJ,IPFREQ).EQ.0) THEN
              WRITE (MP,FMT=9100) ITRMAJ
              WRITE (MP,FMT=9030) (X(I),I=1,N)
          END IF

      END IF

C --------------------------------------------------------
C check if either the number of maximum-allowed Newton
C iterations or function evaluations has been exceeded
C --------------------------------------------------------

      IF (ITRMAJ.GE.MXITN) THEN
          INFORM = -5
          WRITE (MP,FMT=9110) ITRMAJ,MXITN
          GO TO 30

      END IF

      IF (NFUN.GT.MAXNF) THEN
          INFORM = -6
          WRITE (MP,FMT=9120)
          GO TO 30

      END IF

C --------------------------------------------------------
C call INNER to solve for the search vector p
C --------------------------------------------------------

      CALL INNER(N,MODE,MXITCG,ITRMAJ,ITR,IPCG,NZ,NSP,X,G,P,Y,PER,IPER,
     +           IA,JA,A,ISP,RSP,D,HD,Z,XNORM,CALFGH,CALHDP,D1)

      RNORM = DNRM2(N,Y,1)

C --------------------------------------------------------
C save old value of f for convergence tests
C --------------------------------------------------------

      FOLD = F
      ITRCG = ITRCG + ITR

C --------------------------------------------------------
C call line search
C --------------------------------------------------------

      LAMBDA = ONE
      ILINE = 0
      CALL MLINES(N,X,F,G,P,LAMBDA,ILINE,Y,NFEV,CALFGH)
      NFUN = NFUN + NFEV

C --------------------------------------------------------
C summarize progress and prepare for convergence tests
C --------------------------------------------------------

      GNORM = DNRM2(N,G,1)/SQRTN
      XNORM = DNRM2(N,X,1)/SQRTN

      WRITE (MP,FMT=9130) ITRMAJ,F,GNORM,ITR,MODE,RNORM,LAMBDA,NFUN
      CALL OUTPUT

      XSNORM = LAMBDA* (DNRM2(N,P,1)/SQRTN)

C --------------------------------------------------------
C call a routine MONIT, as desired, to print further info
C --------------------------------------------------------
C
C     IF (ICALLS .EQ. 1) CALL MONIT(ITRMAJ,ICALLS,F,GNORM)

C --------------------------------------------------------
C four convergence tests are performed:
C tests 1 & 2 check for convergence of the f and x sequences
C tests 3 & 4 check that |g| is suff. small
C --------------------------------------------------------
c      EPSF2 = 1.E-8

      ONEF = ONE + ABS(F)
      T1 = (FOLD-F) .LT. (EPSF*ONEF)
      T2 = XSNORM .LT. (EPSF2* (ONE+XNORM))
      T3 = GNORM .LT. (EPSF3*ONEF)
      T4 = GNORM .LT. (EPSG*ONEF)
      T123 = T1 .AND. (T2 .AND. T3)

C -------------------------------------------------------------------
C to check conv. tests for difficult problems, uncomment next 4 lines
C
C     WRITE (MP,*) '       t1:',FOLD-F,EPSF*ONEF,' ',t1
C     WRITE (MP,*) '       t2:',XSNORM,EPSF2*(ONE + XNORM),' ',t2
C     WRITE (MP,*) '       t3:',GNORM,EPSF3*ONEF,' ',t3
C     WRITE (MP,*) '       t4:',GNORM,EPSG*ONEF,' ',t4
C -------------------------------------------------------------------

      IF (ILINE.EQ.1) THEN
          IF (T123 .OR. T4) INFORM = 1

      ELSE
          IF (T3 .OR. T4) THEN
              INFORM = 2

          ELSE
              INFORM = 3
          END IF

          IF (ILINE.EQ.0) THEN
              WRITE (MP,FMT=9140)

          ELSE IF (ILINE.EQ.2) THEN
              WRITE (MP,FMT=9150)

          ELSE IF (ILINE.EQ.3) THEN
              WRITE (MP,FMT=9160)

          ELSE IF (ILINE.EQ.4) THEN
              WRITE (MP,FMT=9170)

          ELSE IF (ILINE.EQ.5) THEN
              WRITE (MP,FMT=9180)

          ELSE IF (ILINE.EQ.6) THEN
              WRITE (MP,FMT=9190)

          ELSE IF (ILINE.EQ.7) THEN
              WRITE (MP,FMT=9200)
          END IF

      END IF


      IF (INFORM.EQ.1 .OR. ILINE.NE.1) THEN
          WRITE (MP,FMT=9210) ITRCG
          WRITE (MP,FMT=9220) INFORM
          IF (ILINE.NE.1 .AND. INFORM.NE.1) WRITE (MP,FMT=9230)
          IF (T1) WRITE (MP,FMT=9240)
          IF (T2) WRITE (MP,FMT=9250)
          IF (T3) WRITE (MP,FMT=9260)
          IF (T4) WRITE (MP,FMT=9270)
          WRITE (MP,FMT=9280)
          IF (IPRINT.GE.1) THEN
              WRITE (MP,FMT=9290)
              WRITE (MP,FMT=9030) (X(I),I=1,N)
              IF (IPRINT.GE.2) THEN
                  WRITE (MP,FMT=9300)
                  WRITE (MP,FMT=9030) (G(I),I=1,N)
              END IF

          END IF

          GO TO 30

      END IF

C --------------------------------------------------------
C line search OK and conv. criteria are not satisfied;
C compute H and M, prepare for symmetric reordering of M,
C and continue main loop
C --------------------------------------------------------

      IF (IPCG.EQ.0 .AND. IHD.EQ.1) GO TO 20
      IF (IPCG.EQ.0 .AND. IHD.EQ.0) THEN
          NOUT = 4

      ELSE IF (IPCG.EQ.1 .AND. IHD.EQ.0) THEN
          NOUT = 5

      ELSE
          NOUT = 6
      END IF

      CALL CALFGH(N,X,F,G,A,IA,JA,NOUT)

      IF (IPCG.EQ.1 .AND. IORDER.EQ.1) THEN
          OFLAG = 0
          OPATH = 5

c---xie: re-calculate IA and JA because they have been changed by
c        using reordering routine ODRV
          CALL CALPAT(N,X,A,IA,JA)
c---------------------------------------------------

          CALL ODRV(N,IA,JA,A,PER,IPER,NSP,ISP,OPATH,OFLAG)
          IF (OFLAG.NE.0) THEN
              WRITE (MP,FMT=9310) OFLAG
              GO TO 30

          END IF

      END IF

      GO TO 20
C
C --------------------------------------------------------
C  MAIN LOOP ENDS
C --------------------------------------------------------

   30 CONTINUE


      RETURN

 9000 FORMAT (/,1X,' INSUFFICIENT STORAGE FOR M. LENGTH of M  = ',I8,/,
     +       1X,' STORAGE AVAILABLE = N+NZ = ',I8,/)
 9010 FORMAT (/,1X,'*****   upper-M has ',I8,' nonzeros (',F6.2,
     +       ' % )  ***** ',/)
 9020 FORMAT (/,'  INITIAL X:  ',/)
 9030 FORMAT (6 (3X,F10.6))
 9040 FORMAT (/,'  INITIAL G:  ',/)
 9050 FORMAT (2X,'____________________________',/,'  ITN.',1X,
     +       '    F    ',5X,'   |G|   ',2X,' #CG ',3X,'mode',2X,
     +       '  |R|   ',2X,'  LAMBDA',3X,' NFG',/,I5,1P,E13.4,1P,E12.2)
 9060 FORMAT (' INFORM = ',I5,/,' (|G(X0)| sufficiently small)',/,1X,
     +       '____________________________',/)
 9070 FORMAT (/,1X,' PERMUTATION ARRAY FROM ODRV: ',/)
 9080 FORMAT (12 (1X,I5))
 9090 FORMAT (/,1X,' OFLAG =  ',I6,' SFLAG =  ',I6,/,1X,
     +       'INSUFF. STORAGE',
     +       ' FOR YSMP (INCREASE NZ) OR DUPLICATE ENTRY IN A',/)
 9100 FORMAT (/,'  X AT ITRMAJ',I5,':',/)
 9110 FORMAT (/,1X,' MXITN (Max Newton iterations) EXCEEDED!',2 (I8),/)
 9120 FORMAT (/,1X,' MAXNF (Max function evaluations) EXCEEDED!',/)
 9130 FORMAT (I5,1P,E13.4,1P,E12.2,1X,I5,4X,I2,1X,1P,E11.2,1P,E11.2,I5)
 9140 FORMAT (/,1X,' IMPROPER INPUT PARAMETERS DURING THE LINE SEARCH.',
     +       /)
 9150 FORMAT (/,1X,' RELATIVE WIDTH OF THE INTERVAL OF UNCERTAINTY',
     +       ' IN THE LINE SEARCH',/,'IS AT MOST XTOL.',/)
 9160 FORMAT (/,1X,' NUMBER OF CALLS TO FUNCTION IN THE',
     +       ' LINE SEARCH HAS REACHED MAXFEV.')
 9170 FORMAT (/,1X,' THE STEP IN THE LINE SEARCH IS TOO SMALL.',/)
 9180 FORMAT (/,1X,' THE STEP IN THE LINE SEARCH IS TOO LARGE.',/)
 9190 FORMAT (/,1X,' ROUNDING ERRORS PREVENT FURTHER PROGRESS IN ',/,
     +       '  THE LINE SEARCH.',/)
 9200 FORMAT (/,1X,' P IS NOT A DESCENT DIRECTION.',/)
 9210 FORMAT (/,1X,'____________________________',2X,'#CG,tot',/,29X,I7)
 9220 FORMAT (' INFORM = ',I5)
 9230 FORMAT (/,' CHECK LS PARAMETERS, TRY RELAXING CONV. CRITERIA.',/,
     +       ' AND TRY A DIFFERENT TRUNCATION PARAMETER).',/)
 9240 FORMAT (' CONV. TEST 1 SATISFIED ')
 9250 FORMAT (' CONV. TEST 2 SATISFIED ')
 9260 FORMAT (' CONV. TEST 3 SATISFIED ')
 9270 FORMAT (' CONV. TEST 4 SATISFIED ')
 9280 FORMAT (1X,'____________________________')
 9290 FORMAT (/,'  FINAL X:  ',/)
 9300 FORMAT (/,'  FINAL G:  ',/)
 9310 FORMAT (/,1X,' IN TNMIN, after ODRV for new M, OFLAG =  ',I6,/,1X,
     +       'INSUFFICIENT STORAGE FOR YSMP. INCREASE NZ.',/)
      END

C***********************************************************************
C***********************************************************************
      SUBROUTINE INNER(N,MODE,MXITCG,ITRMAJ,ITR,IPCG,NZ,NSP,X,G,P,RES,
     +                 PER,IPER,IA,JA,A,ISP,RSP,D,HD,Z,XNORM,CALFGH,
     +                 CALHDP,D1)

C --------------------------------------------------------
C INNER:  Inner PCG Iteration
C --------------------------------------------------------

C --------------------------------------------------------
C subroutines and functions called:
C        calhdp (user-supplied), ourhdp,
C        sdrvmd, dcopy, daxpy, dnrm2, ddot
C --------------------------------------------------------


C --------------------------------------------------------
C STEP 1 - Initialization
C
C   (a) Compute ||g||, set p=0, and res= -g
C   (b) If preconditioning, perform the numerical
C       factorization of M and solve Mz=res;
C       else (M=I) set z=res
C   (c) Set d=z and compute delta=d*d and z*res
C --------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION XNORM
      INTEGER IPCG,ITR,ITRMAJ,MODE,MXITCG,N,NSP,NZ
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(N+NZ),D(N),D1(N),G(N),HD(N),P(N),RES(N),
     +                 RSP(NSP),X(N),Z(N)
      INTEGER IA(N+1),IPER(N),ISP(NSP),JA(N+NZ),PER(N)
C     ..
C     .. Subroutine Arguments ..
      EXTERNAL CALFGH,CALHDP
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION ETAQ,ETAR,PRODK,TAU
      INTEGER ICALLS,IEXIT,IHD,ITT,MC,MP,SRLS,TA
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ALPHA,BETA,DELTA,ETA,GNORM,PRODCT,PTG,PTR,QNEW,
     +                 QOLD,RITR,RMAJ,RNORM,RTZ,RTZOLD,TEMP,TEMP_OLD,
     +                 ZETA
      INTEGER ESP,I,SFLAG,SPATH
C     ..
C     .. External Functions ..
      DOUBLE PRECISION DDOT,DNRM2
      EXTERNAL DDOT,DNRM2
C     ..
C     .. External Subroutines ..
      EXTERNAL DAXPY,DCOPY,OURHDP,SDRVMD
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MIN
C     ..
C     .. Common blocks ..
      COMMON /TAINPUT/TA
      COMMON /TAUINPUT/TAU,SRLS,MC
      COMMON /TN001/MP,ICALLS
      COMMON /TN002/PRODK,ETAR,ETAQ,IEXIT,IHD,ITT
C     ..
C     .. Save statement ..
      SAVE /TN001/,/TN002/
C     ..
      TEMP_OLD = 0.D0
      ITR = 0
      GNORM = DNRM2(N,G,1)
      DO I = 1,N
          P(I) = 0.0D0
          RES(I) = -G(I)
      END DO

      QOLD = 0.D0

      IF (IPCG.EQ.1) THEN
          SPATH = 2
          CALL SDRVMD(N,PER,IPER,IA,JA,A,RES,Z,NSP,ISP,RSP,ESP,SPATH,
     +                SFLAG)
          IF (SFLAG.GT.0) THEN
              WRITE (MP,FMT=9000) SFLAG
              GO TO 30

          END IF

      ELSE
          CALL DCOPY(N,RES,1,Z,1)
      END IF

      CALL DCOPY(N,Z,1,D,1)
      DELTA = DDOT(N,D,1,D,1)

      RTZ = DDOT(N,RES,1,Z,1)
      RNORM = DDOT(N,RES,1,RES,1)

      RMAJ = ITRMAJ

C --------------------------------------------------------
C STEP 2 -  MAIN LOOP BEGINS:
C          (a) compute the Hessian-vector product Hd
C          (b) if d*Hd is small, exit with a descent direction
C              chosen according to: On ITR 1, p=-g
C                                   on ITR>1, p= current p
c  Note: we use Test 1A' or Test 2A in (b).
C --------------------------------------------------------

   10 CONTINUE

c--- Test if ITR is larger than the maximum number of PCG iterations
      IF (ITR.GE.MXITCG) THEN
          MODE = 3
          GO TO 30

      END IF


C--- Compute a Hessian/vector product Hd
C IHD=1: Compute Hd by a finite-difference of gradient
C     0: Compute Hd directly
C
C      D  - input vector
C      HD - resultant product vector
C      X  - current coordinate
C      G  - gradient
C--------------------------------------
      IF (IHD.EQ.1) THEN
          CALL OURHDP(N,D,HD,X,G,RSP(1),ITR,XNORM,CALFGH)
      ELSE
          CALL CALHDP(N,D,HD,X,G)
      END IF

      ITR = ITR + 1
      RITR = ITR
      PRODCT = DDOT(N,D,1,HD,1)

      IF (MC.EQ.1) GO TO 20

C--- xie: The singularity test
      ZETA = 1.D-20
      IF (ABS(RTZ).LE.ZETA*RNORM*RNORM .OR. ABS(PRODCT).LE.ZETA) THEN
          IF (ITR.EQ.1) THEN
              MODE = 1
              DO I = 1,N
                  P(I) = -G(I)
              END DO

          ELSE
              MODE = 5
          END IF

          GO TO 30

      END IF
C---------------------------------------------

   20 CONTINUE

c--- xie: Determine descent search directions
c
c    TA = 1:  use Test 1A' (the modified negative curvature test).
c             A default choice
c    TA = 2:  use Test 2A (the strong negative curvature test)
c----------------------------------------------------------

c---xie: Test 1A':
      IF (TA.EQ.1) THEN
          IF (PRODCT.LE. (DELTA*PRODK)) THEN
              IF (ITR.EQ.1) THEN
                  MODE = 1
                  DO I = 1,N
                      P(I) = -G(I)
                  END DO

              ELSE
                  MODE = 2
c--- xie: 8/20/98  To test D as search direction
                  IF (IEXIT.EQ.2) THEN
                      CALL DCOPY(N,D,1,P,1)
                  END IF
c----------------------------------------------
              END IF

              GO TO 30

          END IF

      END IF
c--------------------- end of Test 1A'

      ALPHA = RTZ/PRODCT
      RTZOLD = RTZ

c-- xie: keep the current direction P in D1
      DO I = 1,N
          D1(I) = P(I)
      END DO

c--- Update P:  P = P + ALPHA * D
      CALL DAXPY(N,ALPHA,D,1,P,1)

c--- xie: Test 2A
c
c   MODE = 1: G*P > 0 at the first CG or PCG iteration
c   MODE = 4: G*P > 0 or G*P becomes increasing
c
c   PRODK = 10E-10 as the default value
c---------------------------------------------------------
      IF (TA.EQ.2) THEN
          TEMP = DDOT(N,P,1,G,1)
          IF (TEMP.GE.TEMP_OLD) THEN
              IF (ITR.EQ.1) THEN
                  MODE = 1
                  DO I = 1,N
                      P(I) = -G(I)
                  END DO

              ELSE
                  MODE = 4
c--- exit with the previous PCG iterate (i.e., D1) as the
C    search direction
                  DO I = 1,N
                      P(I) = D1(I)
                  END DO
              END IF

              GO TO 30

          ELSE
              TEMP_OLD = TEMP - ZETA
          END IF

      END IF
c-------------------------- End of Test 2A

c
c--- Update residual vector: RES = RES - ALPHA * HD
c
      CALL DAXPY(N,-ALPHA,HD,1,RES,1)
c
c--- STEP 3 - Truncation Test
c
      IF (ITT.EQ.1) THEN
          PTR = DDOT(N,P,1,RES,1)
          PTG = DDOT(N,P,1,G,1)
          QNEW = 0.5D0* (PTR+PTG)
          IF (RITR* (1.D0- (QOLD/QNEW)).LE.ETAQ) THEN
              MODE = 0
              GO TO 30

          END IF

          QOLD = QNEW

      ELSE
          RNORM = DNRM2(N,RES,1)
          ETA = MIN(ETAR/RMAJ,GNORM)
          IF (RNORM.LE. (GNORM*ETA)) THEN
              MODE = 0
              GO TO 30

          END IF

      END IF
c---------------------- End the test

C --------------------------------------------------------
C STEP 4 - Continuation of PCG.
C          (a)  if preconditioning, solve Mz=res by re-using
C               the factors of M; else (M=I) set z=res,
C          (b)  update d and delta
C          (c)  go on to the next PCG iteration.
C --------------------------------------------------------

      IF (IPCG.EQ.1) THEN
          SPATH = 3
c
c--- Solving Mz=r by the unusual modified Cholesky factorization
c
          CALL SDRVMD(N,PER,IPER,IA,JA,A,RES,Z,NSP,ISP,RSP,ESP,SPATH,
     +                SFLAG)

          IF (SFLAG.GT.0) THEN
              WRITE (MP,FMT=9000) SFLAG
              GO TO 30

          END IF

      ELSE
c
c---  COPIES A VECTOR, RES, TO A VECTOR, Z
c
          CALL DCOPY(N,RES,1,Z,1)
      END IF

      RTZ = DDOT(N,RES,1,Z,1)

      BETA = RTZ/RTZOLD

      DO I = 1,N
          D(I) = BETA*D(I) + Z(I)
      END DO

      DELTA = DDOT(N,D,1,D,1)

      GO TO 10

C --------------------------------------------------------
C  MAIN LOOP ENDS
C --------------------------------------------------------

   30 CONTINUE

      RETURN

 9000 FORMAT (/,1X,' IN TNPCG, SFLAG =  ',I6,/,1X,
     +'INSUFF. STORAGE          FOR YSMP (INCREASE NZ) OR DUPLICATE ENTR
     +Y IN A',/)
      END
C***********************************************************************
C***********************************************************************
      SUBROUTINE OURHDP(N,D,HD,X,G,Y,ITR,XNORM,CALFGH)

C ------------------------------------------------------
C OURHDP: Our numeric Hessian/vector multiplication
C ------------------------------------------------------



C --------------------------------------------------------
C DELH is set to:  2*SQRT(EPS)*(1+|X|) / |D|, with division
C precautions and upper/lower limits
C --------------------------------------------------------

C     .. Parameters ..
      DOUBLE PRECISION ONE,DMAX,DIMAX
      PARAMETER (ONE=1.D0,DMAX=0.1D0,DIMAX=ONE/DMAX)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION XNORM
      INTEGER ITR,N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION D(N),G(N),HD(N),X(N),Y(N)
C     ..
C     .. Subroutine Arguments ..
      EXTERNAL CALFGH
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION EPSMCH,SQEPS2,SQRTN
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION DDEN,DELH,DELHI,DNORM,DNUM,FNEW
      INTEGER I
C     ..
C     .. External Functions ..
      DOUBLE PRECISION DNRM2
      EXTERNAL DNRM2
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MAX
C     ..
C     .. Common blocks ..
      COMMON /TN006/EPSMCH,SQEPS2,SQRTN
C     ..
C     .. Save statement ..
      SAVE /TN006/,DNUM
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION ADUM(1)
      INTEGER IADUM(1),JADUM(1)
C     ..
      IF (ITR.EQ.0) DNUM = SQEPS2* (ONE+XNORM*SQRTN)

      DNORM = DNRM2(N,D,1)

      DDEN = MAX(DNUM*DIMAX,DNORM)
      DELH = MAX(DNUM/DDEN,DNUM*0.01)

      DELHI = ONE/DELH

      DO 10 I = 1,N
          Y(I) = X(I) + DELH*D(I)
   10 CONTINUE

      CALL CALFGH(N,Y,FNEW,HD,ADUM,IADUM,JADUM,1)

      DO 20 I = 1,N
          HD(I) = (HD(I)-G(I))*DELHI
   20 CONTINUE

      RETURN

      END
C***********************************************************************
C***********************************************************************
      SUBROUTINE SETLIS(N,OPLIST,PLIST,INFORM)

C --------------------------------------------------------
C SETLIS:  Set sample values for the TNMIN-call list of
C          options and parameters (OPLIST & PLIST).
C          The user should then change some of these items
C          - to suit the problem - before the TNMIN call
C --------------------------------------------------------


C --------------------------------------------------------
C function called:  dmach (blas)
C --------------------------------------------------------

C     .. Scalar Arguments ..
      INTEGER INFORM,N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION PLIST(20)
      INTEGER OPLIST(20)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION EPSMCH,SQEPS2,SQRTN
      INTEGER ICALLS,MP
C     ..
C     .. External Functions ..
      DOUBLE PRECISION DMACH
      EXTERNAL DMACH
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DBLE,SQRT
C     ..
C     .. Common blocks ..
      COMMON /TN001/MP,ICALLS
      COMMON /TN006/EPSMCH,SQEPS2,SQRTN
C     ..
C     .. Save statement ..
      SAVE /TN001/,/TN006/
C     ..
      ICALLS = 0

C --------------------------------------------------------
C Brief Description of parameters (see manual for details;
C M below refers to the preconditioner):
C --------------------------------------------------------
C INFORM     -   diagnostic/status parameter
C ICALLS     -   counter for number of TNPACK calls (multiple
C                calls may be made for a series of minimizations)
C OPLIST(1)  -   preconditioning option (0 - No, 1 - Yes)
C OPLIST(2)  -   general output controler (0,1,2)
C OPLIST(3)  -   output controler, frequency of printing X (0,1,...)
C OPLIST(4)  -   limit of total Newton iterations
C OPLIST(5)  -   limit of PCG iterations in each Newton iteration
C OPLIST(6)  -   limit of total function evaluations
C OPLIST(7)  -   M-reordering option (0 - Do not reorder, 1 - Reorder)
C OPLIST(8)  -   printing option of M's reordering (permutation array)
C                (0 - Do not print, 1 - Print)
C OPLIST(9)  -   option for specifying whether the permutation array
C                for reordering M is known (0 - Not known, 1 - Known).
C                (This option is useful for multiple TNPACK minimiza-
C                tions, involving reorderings of M, if the sparsity
C                structure of M does not change. The user can reset
C                OPLIST(9) from 0 to 1 before the second TNMIN call).
C OPLIST(10) -   Hessian/vector multiplication option (0 - use a user-
C                supplied routine, 1 - use our default finite-difference
C                design routine); IHD in COMMON/TN002
C OPLIST(11) -   unit number for printing; MP in COMMON/TN001
c
c--- xie: modify OPLIST(12)
C OPLIST(12) -   option selector for combination of exit directions in
C                case of detection of negative curvature for PCG itns.
c                Only 1: -G/P is used. The other selectors (such as
c                0: -G/-G,  2: -G/D, 3: -Y/-G, 4: -Y/P,
C                5: -Y/D, where Y=-[M**(-1)]G) are deleted.
c-----------------------------------------------------------------
c
C OPLIST(13) -   option selector for PCG truncation test
C                (0: residual-based, 1: quadratic-model based);
C                ITT in COMMON/TN002
C OPLIST(14) -   limit for number of function calls in the line search
C                (at each Newton iteration); MAXFEV in COMMON/TN003
c
C PLIST(1)   -   desired   F   accuracy in min. conv. test
C PLIST(2)   -   desired ||G|| accuracy in min. conv. test
C PLIST(3)   -   truncation parameter for residual test;
C                ETAR in COMMON/TN002
C PLIST(4)   -   truncation parameter for quadratic-model test;
C                ETAQ in COMMON/TN002
C PLIST(5)   -   tolerance for 'negative curvature' test;
C                PRODK in COMMON/TN002
C PLIST(6)   -   line-search conv. parameter, controlling function
C                decrease; FTOL in COMMON/TN003
C PLIST(7)   -   line-search conv. parameter, controlling reduction
C                of derivative magnitude; GTOL in COMMON/TN003
c
c--- xie: new parameters
c OPLIST(15) -   option selector for termination rule of PCG
C                1: the modified negative curvature test 1A'
C                2: the descent direction test 2A)
C                Default value is 1
C
C OPLIST(16) -   option selector for line search stopping rule
C                1: Criterion 1
C                2: Criterion 2 (more lenient stopping rule)
C                Default value is 1
C OPLIST(17) -   option selector for modified Cholesky factorization
C                1: the MC by Gill and Murray
C                2: the UMC by Xie and Schlick
C                Default value is 1
c PLIST(8)   -   unusual modified Cholesky factorization parameter TAU
c                Default value is 10.d0
C ---------------------------------------------------------------

      IF (N.LE.0) THEN
          WRITE (*,FMT=9000) N
          RETURN

      END IF

      INFORM = -10

      OPLIST(1) = 0
      OPLIST(2) = 0
      OPLIST(3) = 0
      OPLIST(4) = 100*N
      OPLIST(5) = 10*N
      OPLIST(6) = 100*N
      OPLIST(7) = 0
      OPLIST(8) = 0
      OPLIST(9) = 0
      OPLIST(10) = 1
      OPLIST(11) = 6
      OPLIST(12) = 1
      OPLIST(13) = 0
      OPLIST(14) = 20

      MP = OPLIST(11)
      EPSMCH = DMACH(1)
      SQEPS2 = 2.D0*SQRT(EPSMCH)
      SQRTN = SQRT(DBLE(N))

      PLIST(1) = 1.0D-10
      PLIST(2) = SQRT(EPSMCH)
      PLIST(3) = 0.5D0
      PLIST(4) = 0.5D0
      PLIST(5) = 1.0D-10
      PLIST(6) = 1.0D-04
      PLIST(7) = 0.9D0

c-- xie: defualt new parameters
      PLIST(8) = 1.0D1
      OPLIST(15) = 1
      OPLIST(16) = 1
      OPLIST(17) = 1
c--------------------------------

      RETURN

 9000 FORMAT (/,2X,' N <= 0, N = ',I10)
      END

C***********************************************************************
C***********************************************************************
      SUBROUTINE CHKLIS(N,OPLIST,PLIST,MARK)

C --------------------------------------------------------
C CHKLIS:  Check input options & parameters (OPLIST,PLIST)
C --------------------------------------------------------


C     .. Scalar Arguments ..
      INTEGER MARK,N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION PLIST(20)
      INTEGER OPLIST(20)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION EPSMCH,SQEPS2,SQRTN
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ONE,SQEPS,ZERO
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC SQRT
C     ..
C     .. Common blocks ..
      COMMON /TN006/EPSMCH,SQEPS2,SQRTN
C     ..
C     .. Save statement ..
      SAVE /TN006/
C     ..
      ZERO = 0.D0
      ONE = 1.D0
      SQEPS = SQRT(EPSMCH)

      IF (N.LE.0) THEN
          WRITE (*,FMT=9000) N
          MARK = -1
          RETURN

      END IF

      IF (OPLIST(1).NE.0 .AND. OPLIST(1).NE.1) THEN
          OPLIST(1) = 0
          WRITE (*,FMT=9010)
      END IF

      IF (OPLIST(2).LT.0 .OR. OPLIST(2).GT.2) THEN
          OPLIST(2) = 0
          WRITE (*,FMT=9020)
      END IF

      IF (OPLIST(3).LT.0) THEN
          OPLIST(3) = 0
          WRITE (*,FMT=9030)
      END IF

      IF (OPLIST(4).LE.0 .OR. OPLIST(5).LE.0 .OR. OPLIST(6).LE.0) THEN
          WRITE (*,FMT=9040)
          MARK = -4
          RETURN

      END IF

      IF (OPLIST(1).EQ.0) THEN
          OPLIST(7) = 0
          OPLIST(8) = 0
          OPLIST(9) = 0

      ELSE
          IF (OPLIST(7).NE.0 .AND. OPLIST(7).NE.1) THEN
              OPLIST(7) = 0
              WRITE (*,FMT=9050)
          END IF

          IF (OPLIST(7).EQ.0 .AND. (OPLIST(8).NE.0.OR.
     +        OPLIST(9).NE.0)) THEN
              OPLIST(8) = 0
              OPLIST(9) = 0
          END IF

      END IF

      IF (OPLIST(10).NE.0 .AND. OPLIST(10).NE.1) THEN
          OPLIST(10) = 1
          WRITE (*,FMT=9060)
      END IF

      IF (OPLIST(11).LE.0) THEN
          WRITE (*,FMT=9070)
          MARK = -11
          RETURN

      END IF

      IF (OPLIST(12).LT.0 .OR. OPLIST(12).GT.5) THEN
          OPLIST(12) = 1
          WRITE (*,FMT=9080) OPLIST(12)
      END IF

      IF (OPLIST(13).NE.0 .AND. OPLIST(13).NE.1) THEN
          OPLIST(13) = 0
          WRITE (*,FMT=9090)
      END IF

      IF (OPLIST(14).LE.0 .OR. OPLIST(14).GT.40) THEN
          OPLIST(14) = 20
          WRITE (*,FMT=9100)
      END IF

      IF (PLIST(1).LT.EPSMCH) THEN
          PLIST(1) = EPSMCH*100.D0
          WRITE (*,FMT=9110) PLIST(1)
      END IF

      IF (PLIST(2).LT.EPSMCH) THEN
          PLIST(2) = SQEPS*10.D0
          WRITE (*,FMT=9120) PLIST(2),SQEPS
      END IF

      IF (OPLIST(13).EQ.0) THEN
          IF (PLIST(3).LT.ZERO .OR. PLIST(3).GT.ONE) THEN
              PLIST(3) = 0.5D0
              WRITE (*,FMT=9130) PLIST(3)

          ELSE IF (PLIST(3).LT.1.D-4) THEN
              WRITE (*,FMT=9150)
          END IF

      END IF

      IF (OPLIST(13).EQ.1) THEN
          IF (PLIST(4).LT.ZERO .OR. PLIST(4).GT.ONE) THEN
              PLIST(4) = 0.5D0
              WRITE (*,FMT=9140) PLIST(4)

          ELSE IF (PLIST(4).LT.1.D-4) THEN
              WRITE (*,FMT=9160)
          END IF

      END IF

      IF (PLIST(5).LT.EPSMCH) THEN
          PLIST(5) = SQEPS*0.1D0
          WRITE (*,FMT=9170) PLIST(5)
      END IF

      IF (PLIST(6).LT.ZERO .OR. PLIST(7).LT.ZERO .OR.
     +    PLIST(6).GT.ONE .OR. PLIST(7).GT.ONE .OR.
     +    PLIST(6).GE.PLIST(7)) THEN
          PLIST(6) = 1.0D-04
          PLIST(7) = 0.9D0
          WRITE (*,FMT=9180)
      END IF

c-- xie: set default values of OPLIST(15), OPLIST(16), OPLIST(17),
c-- and PLIST(8)
      IF (OPLIST(15).LE.0 .OR. OPLIST(15).GT.2) OPLIST(15) = 1
      IF (OPLIST(16).LE.0 .OR. OPLIST(16).GT.2) OPLIST(16) = 1
      IF (OPLIST(17).LE.0 .OR. OPLIST(17).GT.2) OPLIST(17) = 1
      IF (PLIST(8).LT.ZERO) PLIST(8) = 1.0D1
c------------------------------------------------


      RETURN

 9000 FORMAT (/,2X,'N <= 0, N = ',I10)
 9010 FORMAT (/,2X,'OPLIST(1)  OUT-OF-RANGE, reset to 0 (no PCG)')
 9020 FORMAT (/,2X,'OPLIST(2)  OUT-OF-RANGE, reset to 0')
 9030 FORMAT (/,2X,'OPLIST(3)  OUT-OF-RANGE, reset to 0')
 9040 FORMAT (/,2X,'OPLIST (4,5, and/or 6) <= 0')
 9050 FORMAT (/,2X,'OPLIST(7)  OUT-OF-RANGE, reset to 0 (no reordering)'
     +       )
 9060 FORMAT (/,2X,
     +       'OPLIST(10) OUT-OF-RANGE, reset to 1 (our Hd routine)')
 9070 FORMAT (/,2X,'OPLIST(11), PRINTING UNIT, <= 0')
 9080 FORMAT (/,2X,'OPLIST(12)  OUT-OF-RANGE, reset to',I10)
 9090 FORMAT (/,2X,
     +       'OPLIST(13)  OUT-OF-RANGE, reset to 0 (residual test)')
 9100 FORMAT (/,2X,'OPLIST(14) <=0 or unreasonably high, reset to 20')
 9110 FORMAT (/,2X,'PLIST(1) < EPSMCH, reset to',1P,E15.5)
 9120 FORMAT (/,2X,'PLIST(2) < EPSMCH, reset to',1P,E15.5,/,2X,
     +       'However, around SQRT(EPSMCH) is recommended:',1P,E15.5)
 9130 FORMAT (/,2X,'PLIST(3)  OUT-OF-RANGE, reset to',1P,E15.5)
 9140 FORMAT (/,2X,'PLIST(4)  OUT-OF-RANGE, reset to',1P,E15.5)
 9150 FORMAT (/,2X,'NOTE: A SMALL INPUT VALUE FOR PLIST(3) MAY NOT',/,
     +       2X,'EXPLOIT THE BENEFIT OF TRUNCATION')
 9160 FORMAT (/,2X,'NOTE: A SMALL INPUT VALUE FOR PLIST(4) MAY NOT',/,
     +       2X,'EXPLOIT THE BENEFIT OF TRUNCATION')
 9170 FORMAT (/,2X,'PLIST(5) < EPSMCH, reset to',1P,E15.5)
 9180 FORMAT (/,2X,'PLIST(6 and/or 7)  OUT-OF-RANGE, both were reset',/,
     +       2X,'to default values')
      END

C***********************************************************************
C***********************************************************************
      SUBROUTINE DIVWRK(N,NZ,LW,LIW,MARK)

C --------------------------------------------------------
C DIVWRK:  Compute starting indices for work vectors W,IW
C --------------------------------------------------------


C --------------------------------------------------------
C W is partitioned to vectors:
C   P(n),Y(n),D(n),HD(n),Z(n),A(n+nz),RSP(nsp),
C IW is partitioned to vectors:
C   PER(n),IPER(n),IA(n+1),JA(n+nz),ISP(nsp),
C where   nz  =  # of nonzeros in strict upper-triang. of M,
C and     nsp =  3*n + 4*max(n,nz).
C --------------------------------------------------------

C     .. Scalar Arguments ..
      INTEGER LIW,LW,MARK,N,NZ
C     ..
C     .. Scalars in Common ..
      INTEGER ICALLS,IXA,IXD,IXD1,IXHD,IXIA,IXIPER,IXISP,IXJA,IXP,IXPER,
     +        IXRSP,IXY,IXZ,MP
C     ..
C     .. Local Scalars ..
      INTEGER NEEDIW,NEEDW,NSP
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MAX
C     ..
C     .. Common blocks ..
      COMMON /TN001/MP,ICALLS
      COMMON /TN004/IXP,IXY,IXD,IXHD,IXZ,IXA,IXRSP,IXD1
      COMMON /TN005/IXPER,IXIPER,IXIA,IXJA,IXISP
C     ..
C     .. Save statement ..
      SAVE /TN001/,/TN004/,/TN005/
C     ..
      NSP = 3*N + 4*MAX(N,NZ)
      NEEDW = 6*N + NZ + NSP
      NEEDIW = 4*N + NZ + NSP + 1
      IF (N.LE.0 .OR. NZ.LT.0 .OR. LW.LT.NEEDW .OR. LIW.LT.NEEDIW) THEN
          WRITE (*,FMT=9000) N,NZ,LW,NEEDW,LIW,NEEDIW
          MARK = -1
          RETURN

      END IF

      IXP = 1
      IXY = IXP + N
      IXD = IXY + N
      IXHD = IXD + N
      IXZ = IXHD + N
      IXA = IXZ + N
      IXRSP = IXA + N + NZ

      IXD1 = IXRSP + 3*N + 4*NZ

      IXPER = 1
      IXIPER = 1 + N
      IXIA = IXIPER + N
      IXJA = IXIA + N + 1
      IXISP = IXJA + N + NZ


      RETURN

 9000 FORMAT (/,2X,'ERROR in DIMENSIONS of N,NZ,LW,or LIW',/,2X,'N = ',
     +       I10,'NZ = ',I10,/,2X,'LW = ',I10,'(NEED ',I10,')','LIW = ',
     +       I10,'(NEED ',I10,')')
      END
C***********************************************************************
C***********************************************************************
      SUBROUTINE MLINES(N,X,F,G,S,STP,INFO,WA,NFEV,CALFGH)

C --------------------------------------------------------
C MLINES:  Line-search algorithm of More' and Thuente
C --------------------------------------------------------

C --------------------------------------------------------
C subroutines and functions called:
C       calfgh       (user-supplied)
C       newstp
C       daxpy, dcopy, ddot  (blas)
C --------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION F,STP
      INTEGER INFO,N,NFEV
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION G(N),S(N),WA(N),X(N)
C     ..
C     .. Subroutine Arguments ..
      EXTERNAL CALFGH
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION FTOL,GTOL,TAU
      INTEGER MAXFEV,MC,SRLS
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION DG,DGINIT,DGM,DGTEST,DGX,DGXM,DGY,DGYM,FINIT,FM,
     +                 FTEST1,FX,FXM,FY,FYM,P5,P66,STMAX,STMIN,STPMAX,
     +                 STPMIN,STX,STY,WIDTH,WIDTH1,XTOL,XTRAPF,ZERO
      INTEGER INFOC
      LOGICAL BRACKT,STAGE1
C     ..
C     .. External Functions ..
      DOUBLE PRECISION DDOT
      EXTERNAL DDOT
C     ..
C     .. External Subroutines ..
      EXTERNAL DAXPY,DCOPY,NEWSTP
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MAX,MIN
C     ..
C     .. Common blocks ..
      COMMON /TAUINPUT/TAU,SRLS,MC
      COMMON /TN003/FTOL,GTOL,MAXFEV
C     ..
C     .. Save statement ..
      SAVE
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION ADUM(1)
      INTEGER IDUM(1),JDUM(1)
C     ..
C     .. Data statements ..

      DATA P5,P66,XTRAPF,ZERO/0.5D0,0.66D0,4.0D0,0.0D0/
      DATA XTOL,STPMIN,STPMAX/1.0D-17,1.0D-20,1.0D+20/
C     ..
      INFOC = 1
C --------------------------------------------------------
C    CHECK THE INPUT PARAMETERS FOR ERRORS.
C --------------------------------------------------------
      IF (N.LE.0 .OR. STP.LE.ZERO .OR. FTOL.LT.ZERO .OR.
     +    GTOL.LT.ZERO .OR. XTOL.LT.ZERO .OR. STPMIN.LT.ZERO .OR.
     +    STPMAX.LT.STPMIN .OR. MAXFEV.LE.0) RETURN
C --------------------------------------------------------
C     COMPUTE THE INITIAL GRADIENT IN THE SEARCH DIRECTION
C     AND CHECK THAT S IS A DESCENT DIRECTION.
C --------------------------------------------------------
      DGINIT = DDOT(N,G,1,S,1)
      IF (DGINIT.GE.ZERO) THEN
          WRITE (*,FMT=*) '  GTP in MLINES = ',DGINIT
          INFO = 7
          RETURN

      END IF
C --------------------------------------------------------
C     INITIALIZE LOCAL VARIABLES.
C --------------------------------------------------------
      BRACKT = .FALSE.
      STAGE1 = .TRUE.
      NFEV = 0
      FINIT = F
      DGTEST = FTOL*DGINIT
      WIDTH = STPMAX - STPMIN
      WIDTH1 = WIDTH/P5
      CALL DCOPY(N,X,1,WA,1)
C --------------------------------------------------------
C     THE VARIABLES STX, FX, DGX CONTAIN THE VALUES OF THE STEP,
C     FUNCTION, AND DIRECTIONAL DERIVATIVE AT THE BEST STEP.
C     THE VARIABLES STY, FY, DGY CONTAIN THE VALUE OF THE STEP,
C     FUNCTION, AND DERIVATIVE AT THE OTHER ENDPOINT OF
C     THE INTERVAL OF UNCERTAINTY.
C     THE VARIABLES STP, F, DG CONTAIN THE VALUES OF THE STEP,
C     FUNCTION, AND DERIVATIVE AT THE CURRENT STEP.
C --------------------------------------------------------
      STX = ZERO
      FX = FINIT
      DGX = DGINIT
      STY = ZERO
      FY = FINIT
      DGY = DGINIT
C --------------------------------------------------------
C     START OF ITERATION.
C --------------------------------------------------------
   10 CONTINUE
C --------------------------------------------------------
C        SET THE MINIMUM AND MAXIMUM STEPS TO CORRESPOND
C        TO THE PRESENT INTERVAL OF UNCERTAINTY.
C --------------------------------------------------------
      IF (BRACKT) THEN
          STMIN = MIN(STX,STY)
          STMAX = MAX(STX,STY)

      ELSE
          STMIN = STX
          STMAX = STP + XTRAPF* (STP-STX)
      END IF
C --------------------------------------------------------
C        FORCE THE STEP TO BE WITHIN THE BOUNDS STPMAX AND STPMIN.
C --------------------------------------------------------
      STP = MAX(STP,STPMIN)
      STP = MIN(STP,STPMAX)
C --------------------------------------------------------
C        IF AN UNUSUAL TERMINATION IS TO OCCUR THEN LET
C        STP BE THE LOWEST POINT OBTAINED SO FAR.
C --------------------------------------------------------
      IF ((BRACKT.AND. (STP.LE.STMIN.OR.STP.GE.STMAX)) .OR.
     +    NFEV.GE.MAXFEV-1 .OR. INFOC.EQ.0 .OR.
     +    (BRACKT.AND.STMAX-STMIN.LE.XTOL*STMAX)) STP = STX
C --------------------------------------------------------
C        EVALUATE THE FUNCTION AND GRADIENT AT STP
C        AND COMPUTE THE DIRECTIONAL DERIVATIVE.
C --------------------------------------------------------
      CALL DCOPY(N,WA,1,X,1)
      CALL DAXPY(N,STP,S,1,X,1)
C
      CALL CALFGH(N,X,F,G,ADUM,IDUM,JDUM,1)
C
      INFO = 0
      NFEV = NFEV + 1
      DG = DDOT(N,G,1,S,1)
      FTEST1 = FINIT + STP*DGTEST
C --------------------------------------------------------
C        TEST FOR CONVERGENCE.
C --------------------------------------------------------
      IF ((BRACKT.AND. (STP.LE.STMIN.OR.STP.GE.STMAX)) .OR.
     +    INFOC.EQ.0) INFO = 6
      IF (STP.EQ.STPMAX .AND. F.LE.FTEST1 .AND. DG.LE.DGTEST) INFO = 5
      IF (STP.EQ.STPMIN .AND. (F.GT.FTEST1.OR.DG.GE.DGTEST)) INFO = 4
      IF (NFEV.GE.MAXFEV) INFO = 3
      IF (BRACKT .AND. STMAX-STMIN.LE.XTOL*STMAX) INFO = 2

c--- xie:  the stopping rule of the line search
c  The original one:
      IF (SRLS.EQ.1) THEN
          IF (F.LE.FTEST1 .AND. ABS(DG).LE.GTOL* (-DGINIT)) INFO = 1

      ELSE
c-- xie: the lenient stopping rule of the line search
          IF (F.LE.FTEST1 .AND. (DG.GE.GTOL*DGINIT.OR.
     +        DG.LE. (2.0D0-GTOL)*DGINIT)) INFO = 1
      END IF
C --------------------------------------------------------
C        CHECK FOR TERMINATION.
C --------------------------------------------------------
      IF (INFO.NE.0) RETURN
C --------------------------------------------------------
C        IN THE FIRST STAGE WE SEEK A STEP FOR WHICH THE MODIFIED
C        FUNCTION HAS A NONPOSITIVE VALUE AND NONNEGATIVE DERIVATIVE.
C --------------------------------------------------------
      IF (STAGE1 .AND. F.LE.FTEST1 .AND.
     +    DG.GE.MIN(FTOL,GTOL)*DGINIT) STAGE1 = .FALSE.
C --------------------------------------------------------
C        A MODIFIED FUNCTION IS USED TO PREDICT THE STEP ONLY IF
C        WE HAVE NOT OBTAINED A STEP FOR WHICH THE MODIFIED
C        FUNCTION HAS A NONPOSITIVE FUNCTION VALUE AND NONNEGATIVE
C        DERIVATIVE, AND IF A LOWER FUNCTION VALUE HAS BEEN
C        OBTAINED BUT THE DECREASE IS NOT SUFFICIENT.
C --------------------------------------------------------
      IF (STAGE1 .AND. F.LE.FX .AND. F.GT.FTEST1) THEN
C --------------------------------------------------------
C           DEFINE THE MODIFIED FUNCTION AND DERIVATIVE VALUES.
C --------------------------------------------------------
          FM = F - STP*DGTEST
          FXM = FX - STX*DGTEST
          FYM = FY - STY*DGTEST
          DGM = DG - DGTEST
          DGXM = DGX - DGTEST
          DGYM = DGY - DGTEST
C --------------------------------------------------------
C           CALL NEWSTP TO UPDATE THE INTERVAL OF UNCERTAINTY
C           AND TO COMPUTE THE NEW STEP.
C --------------------------------------------------------
          CALL NEWSTP(STX,FXM,DGXM,STY,FYM,DGYM,STP,FM,DGM,BRACKT,STMIN,
     +                STMAX,INFOC)
C --------------------------------------------------------
C           RESET THE FUNCTION AND GRADIENT VALUES FOR F.
C --------------------------------------------------------
          FX = FXM + STX*DGTEST
          FY = FYM + STY*DGTEST
          DGX = DGXM + DGTEST
          DGY = DGYM + DGTEST

      ELSE
C --------------------------------------------------------
C           CALL NEWSTP TO UPDATE THE INTERVAL OF UNCERTAINTY
C           AND TO COMPUTE THE NEW STEP.
C --------------------------------------------------------
          CALL NEWSTP(STX,FX,DGX,STY,FY,DGY,STP,F,DG,BRACKT,STMIN,STMAX,
     +                INFOC)
      END IF
C --------------------------------------------------------
C        FORCE A SUFFICIENT DECREASE IN THE SIZE OF THE
C        INTERVAL OF UNCERTAINTY.
C --------------------------------------------------------
      IF (BRACKT) THEN
          IF (ABS(STY-STX).GE.P66*WIDTH1) STP = STX + P5* (STY-STX)
          WIDTH1 = WIDTH
          WIDTH = ABS(STY-STX)
      END IF
C --------------------------------------------------------
C        END OF ITERATION.
C --------------------------------------------------------
      GO TO 10

      END
C
C***********************************************************************
C***********************************************************************
      SUBROUTINE NEWSTP(STX,FX,DX,STY,FY,DY,STP,FP,DP,BRACKT,STPMIN,
     +                  STPMAX,INFO)

C --------------------------------------------------------
C NEWSTP:  Compute safeguarded step for linesearch; update
C          interval of uncertainty
C --------------------------------------------------------

C argument-list variables:


C     .. Scalar Arguments ..
      DOUBLE PRECISION DP,DX,DY,FP,FX,FY,STP,STPMAX,STPMIN,STX,STY
      INTEGER INFO
      LOGICAL BRACKT
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION GAMMA,P,Q,R,S,SGND,STPC,STPF,STPQ,THETA,XSAFE
      LOGICAL BOUND
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MAX,MIN,SQRT
C     ..
C-- xie
C other variables:
C     .. Scalars in Common ..
      DOUBLE PRECISION TAU
      INTEGER MC,SRLS
C     ..
C     .. Common blocks ..
      COMMON /TAUINPUT/TAU,SRLS,MC
C     ..
      XSAFE = 1.0D-3
c-----------------------

      INFO = 0
C --------------------------------------------------------
C     CHECK THE INPUT PARAMETERS FOR ERRORS.
C --------------------------------------------------------
      IF ((BRACKT.AND. (STP.LE.MIN(STX,STY).OR.STP.GE.MAX(STX,
     +    STY))) .OR. DX* (STP-STX).GE.0.0 .OR. STPMAX.LT.STPMIN) RETURN
C --------------------------------------------------------
C     DETERMINE IF THE DERIVATIVES HAVE OPPOSITE SIGN.
C --------------------------------------------------------
      SGND = DP* (DX/ABS(DX))
C --------------------------------------------------------
C     FIRST CASE. A HIGHER FUNCTION VALUE.
C     THE MINIMUM IS BRACKETED. IF THE CUBIC STEP IS CLOSER
C     TO STX THAN THE QUADRATIC STEP, THE CUBIC STEP IS TAKEN,
C     ELSE THE AVERAGE OF THE CUBIC AND QUADRATIC STEPS IS TAKEN.
C --------------------------------------------------------
      IF (FP.GT.FX) THEN
          INFO = 1
          BOUND = .TRUE.
          THETA = 3* (FX-FP)/ (STP-STX) + DX + DP
          S = MAX(ABS(THETA),ABS(DX),ABS(DP))
          GAMMA = S*SQRT((THETA/S)**2- (DX/S)* (DP/S))
          IF (STP.LT.STX) GAMMA = -GAMMA
          P = (GAMMA-DX) + THETA
          Q = ((GAMMA-DX)+GAMMA) + DP
          R = P/Q
          STPC = STX + R* (STP-STX)
          STPQ = STX + ((DX/ ((FX-FP)/ (STP-STX)+DX))/2)* (STP-STX)
          IF (ABS(STPC-STX).LT.ABS(STPQ-STX)) THEN
              STPF = STPC

          ELSE
              STPF = STPC + (STPQ-STPC)/2
          END IF

c--- xie: to avoid a too small STPE
          IF (SRLS.EQ.1) GO TO 10

          IF (STP.GT.STX) THEN
              STPF = MAX(STX+XSAFE* (STP-STX),STPF)

          ELSE
              STPF = MIN(STX+XSAFE* (STP-STX),STPF)
          END IF
c---------------------------------------------------------------

   10     BRACKT = .TRUE.
C --------------------------------------------------------
C     SECOND CASE. A LOWER FUNCTION VALUE AND DERIVATIVES OF
C     OPPOSITE SIGN. THE MINIMUM IS BRACKETED. IF THE CUBIC
C     STEP IS CLOSER TO STX THAN THE QUADRATIC (SECANT) STEP,
C     THE CUBIC STEP IS TAKEN, ELSE THE QUADRATIC STEP IS TAKEN.
C --------------------------------------------------------
      ELSE IF (SGND.LT.0.0) THEN
          INFO = 2
          BOUND = .FALSE.
          THETA = 3* (FX-FP)/ (STP-STX) + DX + DP
          S = MAX(ABS(THETA),ABS(DX),ABS(DP))
          GAMMA = S*SQRT((THETA/S)**2- (DX/S)* (DP/S))
          IF (STP.GT.STX) GAMMA = -GAMMA
          P = (GAMMA-DP) + THETA
          Q = ((GAMMA-DP)+GAMMA) + DX
          R = P/Q
          STPC = STP + R* (STX-STP)
          STPQ = STP + (DP/ (DP-DX))* (STX-STP)
          IF (ABS(STPC-STP).GT.ABS(STPQ-STP)) THEN
              STPF = STPC

          ELSE
              STPF = STPQ
          END IF

          BRACKT = .TRUE.
C --------------------------------------------------------
C     THIRD CASE. A LOWER FUNCTION VALUE, DERIVATIVES OF THE
C     SAME SIGN, AND THE MAGNITUDE OF THE DERIVATIVE DECREASES.
C     THE CUBIC STEP IS ONLY USED IF THE CUBIC TENDS TO INFINITY
C     IN THE DIRECTION OF THE STEP OR IF THE MINIMUM OF THE CUBIC
C     IS BEYOND STP. OTHERWISE THE CUBIC STEP IS DEFINED TO BE
C     EITHER STPMIN OR STPMAX. THE QUADRATIC (SECANT) STEP IS ALSO
C     COMPUTED AND IF THE MINIMUM IS BRACKETED THEN THE THE STEP
C     CLOSEST TO STX IS TAKEN, ELSE THE STEP FARTHEST AWAY IS TAKEN.
C --------------------------------------------------------
      ELSE IF (ABS(DP).LT.ABS(DX)) THEN
          INFO = 3
          BOUND = .TRUE.
          THETA = 3* (FX-FP)/ (STP-STX) + DX + DP
          S = MAX(ABS(THETA),ABS(DX),ABS(DP))
C --------------------------------------------------------
C        THE CASE GAMMA = 0 ONLY ARISES IF THE CUBIC DOES NOT TEND
C        TO INFINITY IN THE DIRECTION OF THE STEP.
C --------------------------------------------------------
          GAMMA = S*SQRT(MAX(0.0D0, (THETA/S)**2- (DX/S)* (DP/S)))
          IF (STP.GT.STX) GAMMA = -GAMMA
          P = (GAMMA-DP) + THETA
          Q = (GAMMA+ (DX-DP)) + GAMMA
          R = P/Q
          IF (R.LT.0.0 .AND. GAMMA.NE.0.0) THEN
              STPC = STP + R* (STX-STP)

          ELSE IF (STP.GT.STX) THEN
              STPC = STPMAX

          ELSE
              STPC = STPMIN
          END IF

          STPQ = STP + (DP/ (DP-DX))* (STX-STP)
          IF (BRACKT) THEN
              IF (ABS(STP-STPC).LT.ABS(STP-STPQ)) THEN
                  STPF = STPC

              ELSE
                  STPF = STPQ
              END IF

          ELSE
              IF (ABS(STP-STPC).GT.ABS(STP-STPQ)) THEN
                  STPF = STPC

              ELSE
                  STPF = STPQ
              END IF

          END IF
C --------------------------------------------------------
C     FOURTH CASE. A LOWER FUNCTION VALUE, DERIVATIVES OF THE
C     SAME SIGN, AND THE MAGNITUDE OF THE DERIVATIVE DOES
C     NOT DECREASE. IF THE MINIMUM IS NOT BRACKETED, THE STEP
C     IS EITHER STPMIN OR STPMAX, ELSE THE CUBIC STEP IS TAKEN.
C --------------------------------------------------------

      ELSE
          INFO = 4
          BOUND = .FALSE.
          IF (BRACKT) THEN
              THETA = 3* (FP-FY)/ (STY-STP) + DY + DP
              S = MAX(ABS(THETA),ABS(DY),ABS(DP))
              GAMMA = S*SQRT((THETA/S)**2- (DY/S)* (DP/S))
              IF (STP.GT.STY) GAMMA = -GAMMA
              P = (GAMMA-DP) + THETA
              Q = ((GAMMA-DP)+GAMMA) + DY
              R = P/Q
              STPC = STP + R* (STY-STP)
              STPF = STPC

          ELSE IF (STP.GT.STX) THEN
              STPF = STPMAX

          ELSE
              STPF = STPMIN
          END IF

      END IF
C --------------------------------------------------------
C     UPDATE THE INTERVAL OF UNCERTAINTY. THIS UPDATE DOES NOT
C     DEPEND ON THE NEW STEP OR THE CASE ANALYSIS ABOVE.
C --------------------------------------------------------
c--- xie: according the scheme, SGND is defined by
c     SGND = DP*(STX - STP)
c    But in the original code it was missed, leading to
c     SGND = DP*(DX/ABS(DX))
c    So we add it as blow
      SGND = DP* (STX-STP)
c--------------------------------
      IF (FP.GT.FX) THEN
          STY = STP
          FY = FP
          DY = DP

      ELSE
          IF (SGND.LT.0.0) THEN
              STY = STX
              FY = FX
              DY = DX
          END IF

          STX = STP
          FX = FP
          DX = DP
      END IF
C --------------------------------------------------------
C     COMPUTE THE NEW STEP AND SAFEGUARD IT.
C --------------------------------------------------------
      STPF = MIN(STPMAX,STPF)
      STPF = MAX(STPMIN,STPF)
      STP = STPF
      IF (BRACKT .AND. BOUND) THEN
          IF (STY.GT.STX) THEN
              STP = MIN(STX+0.66* (STY-STX),STP)

          ELSE
              STP = MAX(STX+0.66* (STY-STX),STP)
          END IF

      END IF

      RETURN

      END
C***********************************************************************
C***********************************************************************
      DOUBLE PRECISION FUNCTION DMACH(JOB)
C --------------------------------------------------------
C  SMACH COMPUTES MACHINE PARAMETERS OF FLOATING POINT
C  ARITHMETIC FOR USE IN TESTING ONLY.  NOT REQUIRED BY
C  LINPACK PROPER.
C  IF TROUBLE WITH AUTOMATIC COMPUTATION OF THESE QUANTITIES,
C  THEY CAN BE SET BY DIRECT ASSIGNMENT STATEMENTS.
C  ASSUME THE COMPUTER HAS
C     B = BASE OF ARITHMETIC
C     T = NUMBER OF BASE  B  DIGITS
C     L = SMALLEST POSSIBLE EXPONENT
C     U = LARGEST POSSIBLE EXPONENT
C  THEN
C     EPS = B**(1-T)
C     TINY = 100.0*B**(-L+T)
C     HUGE = 0.01*B**(U-T)
C  DMACH SAME AS SMACH EXCEPT T, L, U APPLY TO
C  DOUBLE PRECISION.
C  CMACH SAME AS SMACH EXCEPT IF COMPLEX DIVISION
C  IS DONE BY
C     1/(X+I*Y) = (X-I*Y)/(X**2+Y**2)
C  THEN
C     TINY = SQRT(TINY)
C     HUGE = SQRT(HUGE)
C  JOB IS 1, 2 OR 3 FOR EPSILON, TINY AND HUGE, RESPECTIVELY.
C --------------------------------------------------------

C     .. Scalar Arguments ..
      INTEGER JOB
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION EPS,HUGE,S,TINY
C     ..
      EPS = 1.0D0
   10 EPS = EPS/2.0D0
      S = 1.0D0 + EPS
      IF (S.GT.1.0D0) GO TO 10
      EPS = 2.0D0*EPS

      S = 1.0D0
   20 TINY = S
      S = S/16.0D0
      IF (S*1.0.NE.0.0D0) GO TO 20
      TINY = (TINY/EPS)*100.0
      HUGE = 1.0D0/TINY

      IF (JOB.EQ.1) DMACH = EPS
      IF (JOB.EQ.2) DMACH = TINY
      IF (JOB.EQ.3) DMACH = HUGE
      RETURN

      END
C***********************************************************************
C                                                                1/15/81
C***********************************************************************
C  ODRV -- DRIVER FOR SPARSE MATRIX REORDERING ROUTINES
C***********************************************************************
      SUBROUTINE ODRV(N,IA,JA,A,P,IP,NSP,ISP,PATH,FLAG)
C
C  DESCRIPTION
C
C    ODRV FINDS A MINIMUM DEGREE ORDERING OF THE ROWS AND COLUMNS OF A
C    SYMMETRIC MATRIX M STORED IN (IA,JA,A) FORMAT (SEE BELOW).  FOR THE
C    REORDERED MATRIX, THE WORK AND STORAGE REQUIRED TO PERFORM GAUSSIAN
C    ELIMINATION IS (USUALLY) SIGNIFICANTLY LESS.
C
C    IF ONLY THE NONZERO ENTRIES IN THE UPPER TRIANGLE OF M ARE BEING
C    STORED, THEN ODRV SYMMETRICALLY REORDERS (IA,JA,A), (OPTIONALLY)
C    WITH THE DIAGONAL ENTRIES PLACED FIRST IN EACH ROW.  THIS IS TO
C    ENSURE THAT IF M(I,J) WILL BE IN THE UPPER TRIANGLE OF M WITH
C    RESPECT TO THE NEW ORDERING, THEN M(I,J) IS STORED IN ROW I (AND
C    THUS M(J,I) IS NOT STORED);  WHEREAS IF M(I,J) WILL BE IN THE
C    STRICT LOWER TRIANGLE OF M, THEN M(J,I) IS STORED IN ROW J (AND
C    THUS M(I,J) IS NOT STORED).
C
C
C  STORAGE OF SPARSE MATRICES
C
C    THE NONZERO ENTRIES OF THE MATRIX M ARE STORED ROW-BY-ROW IN THE
C    ARRAY A.  TO IDENTIFY THE INDIVIDUAL NONZERO ENTRIES IN EACH ROW,
C    WE NEED TO KNOW IN WHICH COLUMN EACH ENTRY LIES.  THESE COLUMN
C    INDICES ARE STORED IN THE ARRAY JA;  I.E., IF  A(K) = M(I,J),  THEN
C    JA(K) = J.  TO IDENTIFY THE INDIVIDUAL ROWS, WE NEED TO KNOW WHERE
C    EACH ROW STARTS.  THESE ROW POINTERS ARE STORED IN THE ARRAY IA;
C    I.E., IF M(I,J) IS THE FIRST NONZERO ENTRY (STORED) IN THE I-TH ROW
C    AND  A(K) = M(I,J),  THEN  IA(I) = K.  MOREOVER, IA(N+1) POINTS TO
C    THE FIRST LOCATION FOLLOWING THE LAST ELEMENT IN THE LAST ROW.
C    THUS, THE NUMBER OF ENTRIES IN THE I-TH ROW IS  IA(I+1) - IA(I),
C    THE NONZERO ENTRIES IN THE I-TH ROW ARE STORED CONSECUTIVELY IN
C
C            A(IA(I)),  A(IA(I)+1),  ..., A(IA(I+1)-1),
C
C    AND THE CORRESPONDING COLUMN INDICES ARE STORED CONSECUTIVELY IN
C
C            JA(IA(I)), JA(IA(I)+1), ..., JA(IA(I+1)-1).
C
C    SINCE THE COEFFICIENT MATRIX IS SYMMETRIC, ONLY THE NONZERO ENTRIES
C    IN THE UPPER TRIANGLE NEED BE STORED.  FOR EXAMPLE, THE MATRIX
C
C             ( 1  0  2  3  0 )
C             ( 0  4  0  0  0 )
C         M = ( 2  0  5  6  0 )
C             ( 3  0  6  7  8 )
C             ( 0  0  0  8  9 )
C
C    COULD BE STORED AS
C
C            \ 1  2  3  4  5  6  7  8  9 10 11 12 13
C         ---+--------------------------------------
C         IA \ 1  4  5  8 12 14
C         JA \ 1  3  4  2  1  3  4  1  3  4  5  4  5
C          A \ 1  2  3  4  2  5  6  3  6  7  8  8  9
C
C    OR (SYMMETRICALLY) AS
C
C            \ 1  2  3  4  5  6  7  8  9
C         ---+--------------------------
C         IA \ 1  4  5  7  9 10
C         JA \ 1  3  4  2  3  4  4  5  5
C          A \ 1  2  3  4  5  6  7  8  9          .
C
C
C  PARAMETERS
C
C    N    - ORDER OF THE MATRIX
C
C    IA   - INTEGER ONE-DIMENSIONAL ARRAY CONTAINING POINTERS TO DELIMIT
C           ROWS IN JA AND A;  DIMENSION = N+1
C
C    JA   - INTEGER ONE-DIMENSIONAL ARRAY CONTAINING THE COLUMN INDICES
C           CORRESPONDING TO THE ELEMENTS OF A;  DIMENSION = NUMBER OF
C           NONZERO ENTRIES IN (THE UPPER TRIANGLE OF) M
C
C    A    - REAL ONE-DIMENSIONAL ARRAY CONTAINING THE NONZERO ENTRIES IN
C           (THE UPPER TRIANGLE OF) M, STORED BY ROWS;  DIMENSION =
C           NUMBER OF NONZERO ENTRIES IN (THE UPPER TRIANGLE OF) M
C
C    P    - INTEGER ONE-DIMENSIONAL ARRAY USED TO RETURN THE PERMUTATION
C           OF THE ROWS AND COLUMNS OF M CORRESPONDING TO THE MINIMUM
C           DEGREE ORDERING;  DIMENSION = N
C
C    IP   - INTEGER ONE-DIMENSIONAL ARRAY USED TO RETURN THE INVERSE OF
C           THE PERMUTATION RETURNED IN P;  DIMENSION = N
C
C    NSP  - DECLARED DIMENSION OF THE ONE-DIMENSIONAL ARRAY ISP;  NSP
C           MUST BE AT LEAST  3N+4K,  WHERE K IS THE NUMBER OF NONZEROES
C           IN THE STRICT UPPER TRIANGLE OF M
C
C    ISP  - INTEGER ONE-DIMENSIONAL ARRAY USED FOR WORKING STORAGE;
C           DIMENSION = NSP
C
C    PATH - INTEGER PATH SPECIFICATION;  VALUES AND THEIR MEANINGS ARE -
C             1  FIND MINIMUM DEGREE ORDERING ONLY
C             2  FIND MINIMUM DEGREE ORDERING AND REORDER SYMMETRICALLY
C                  STORED MATRIX (USED WHEN ONLY THE NONZERO ENTRIES IN
C                  THE UPPER TRIANGLE OF M ARE BEING STORED)
C             3  REORDER SYMMETRICALLY STORED MATRIX AS SPECIFIED BY
C                  INPUT PERMUTATION (USED WHEN AN ORDERING HAS ALREADY
C                  BEEN DETERMINED AND ONLY THE NONZERO ENTRIES IN THE
C                  UPPER TRIANGLE OF M ARE BEING STORED)
C             4  SAME AS 2 BUT PUT DIAGONAL ENTRIES AT START OF EACH ROW
C             5  SAME AS 3 BUT PUT DIAGONAL ENTRIES AT START OF EACH ROW
C
C    FLAG - INTEGER ERROR FLAG;  VALUES AND THEIR MEANINGS ARE -
C               0    NO ERRORS DETECTED
C              9N+K  INSUFFICIENT STORAGE IN MD
C             10N+1  INSUFFICIENT STORAGE IN ODRV
C             11N+1  ILLEGAL PATH SPECIFICATION
C
C
C  CONVERSION FROM REAL TO DOUBLE PRECISION
C
C    CHANGE THE REAL DECLARATIONS IN ODRV AND SRO TO DOUBLE PRECISION
C    DECLARATIONS.
C
C-----------------------------------------------------------------------
C
C....   REAL  A(1)
C
C----INITIALIZE ERROR FLAG AND VALIDATE PATH SPECIFICATION
C     .. Scalar Arguments ..
      INTEGER FLAG,N,NSP,PATH
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(*)
      INTEGER IA(*),IP(*),ISP(*),JA(*),P(*)
C     ..
C     .. Local Scalars ..
      INTEGER HEAD,L,MAX,NEXT,Q,TMP,V
      LOGICAL DFLAG
C     ..
C     .. External Subroutines ..
      EXTERNAL MD,SRO
C     ..
      FLAG = 0
      IF (PATH.LT.1 .OR. 5.LT.PATH) GO TO 50
C
C----ALLOCATE STORAGE AND FIND MINIMUM DEGREE ORDERING
      IF ((PATH-1)* (PATH-2)* (PATH-4).NE.0) GO TO 10
      MAX = (NSP-N)/2
      V = 1
      L = V + MAX
      HEAD = L + MAX
      NEXT = HEAD + N
      IF (MAX.LT.N) GO TO 40
C
      CALL MD(N,IA,JA,MAX,ISP(V),ISP(L),ISP(HEAD),P,IP,ISP(V),FLAG)
      IF (FLAG.NE.0) GO TO 30
C
C----ALLOCATE STORAGE AND SYMMETRICALLY REORDER MATRIX
   10 IF ((PATH-2)* (PATH-3)* (PATH-4)* (PATH-5).NE.0) GO TO 20
      TMP = (NSP+1) - N
      Q = TMP - (IA(N+1)-1)
      IF (Q.LT.1) GO TO 40
C
      DFLAG = PATH .EQ. 4 .OR. PATH .EQ. 5
      CALL SRO(N,IP,IA,JA,A,ISP(TMP),ISP(Q),DFLAG)
C
   20 RETURN
C
C ** ERROR -- ERROR DETECTED IN MD
   30 RETURN
C ** ERROR -- INSUFFICIENT STORAGE
   40 FLAG = 10*N + 1
      RETURN
C ** ERROR -- ILLEGAL PATH SPECIFIED
   50 FLAG = 11*N + 1
      RETURN

      END
C***********************************************************************
C  SSF --  SYMBOLIC UT-D-U FACTORIZATION OF SPARSE SYMMETRIC MATRIX
C***********************************************************************
      SUBROUTINE SSF(N,P,IP,IA,JA,IJU,JU,IU,JUMAX,Q,MARK,JL,FLAG)
C
C  ADDITIONAL PARAMETERS
C
C    Q     - INTEGER ONE-DIMENSIONAL WORK ARRAY;  DIMENSION = N
C
C    MARK  - INTEGER ONE-DIMENSIONAL WORK ARRAY;  DIMENSION = N
C
C    JL    - INTEGER ONE-DIMENSIONAL WORK ARRAY;  DIMENSION = N
C
C
C  DEFINITIONS OF INTERNAL PARAMETERS (DURING K-TH STAGE OF ELIMINATION)
C
C    Q CONTAINS AN ORDERED LINKED LIST REPRESENTATION OF THE NONZERO
C      STRUCTURE OF THE K-TH ROW OF U --
C        Q(K) IS THE FIRST COLUMN WITH A NONZERO ENTRY
C        Q(I) IS THE NEXT COLUMN WITH A NONZERO ENTRY AFTER COLUMN I
C      IN EITHER CASE, Q(I) = N+1 INDICATES THE END OF THE LIST
C
C    JL CONTAINS LISTS OF ROWS TO BE MERGED INTO UNELIMINATED ROWS --
C        I GE K => JL(I) IS THE FIRST ROW TO BE MERGED INTO ROW I
C        I LT K => JL(I) IS THE ROW FOLLOWING ROW I IN SOME LIST OF ROWS
C      IN EITHER CASE, JL(I) = 0 INDICATES THE END OF A LIST
C
C    MARK(I) IS THE LAST ROW STORED IN JU FOR WHICH U(MARK(I),I) NE 0
C
C    JUMIN AND JUPTR ARE THE INDICES IN JU OF THE FIRST AND LAST
C      ELEMENTS IN THE LAST ROW SAVED IN JU
C
C    LUK IS THE NUMBER OF NONZERO ENTRIES IN THE K-TH ROW
C
C-----------------------------------------------------------------------
C
C
C----INITIALIZATION
C     .. Scalar Arguments ..
      INTEGER FLAG,JUMAX,N
C     ..
C     .. Array Arguments ..
      INTEGER IA(*),IJU(*),IP(*),IU(*),JA(*),JL(*),JU(*),MARK(*),P(*),
     +        Q(*)
C     ..
C     .. Local Scalars ..
      INTEGER I,J,JMAX,JMIN,JUMIN,JUPTR,K,LMAX,LUI,LUK,M,QM,TAG,VJ
      LOGICAL CLIQUE
C     ..
      JUMIN = 1
      JUPTR = 0
      IU(1) = 1
      DO 10 K = 1,N
          MARK(K) = 0
          JL(K) = 0
   10 CONTINUE
C
C----FOR EACH ROW K
      DO 190 K = 1,N
          LUK = 0
          Q(K) = N + 1
C
          TAG = MARK(K)
          CLIQUE = .FALSE.
          IF (JL(K).NE.0) CLIQUE = JL(JL(K)) .EQ. 0
C
C------INITIALIZE NONZERO STRUCTURE OF K-TH ROW TO ROW P(K) OF M
          JMIN = IA(P(K))
          JMAX = IA(P(K)+1) - 1
          IF (JMIN.GT.JMAX) GO TO 40
          DO 30 J = JMIN,JMAX
              VJ = IP(JA(J))
              IF (VJ.LE.K) GO TO 30
C
              QM = K
   20         M = QM
              QM = Q(M)
              IF (QM.LT.VJ) GO TO 20
              IF (QM.EQ.VJ) GO TO 200
              LUK = LUK + 1
              Q(M) = VJ
              Q(VJ) = QM
              IF (MARK(VJ).NE.TAG) CLIQUE = .FALSE.
C
   30     CONTINUE
C
C------IF EXACTLY ONE ROW IS TO BE MERGED INTO THE K-TH ROW AND THERE IS
C------A NONZERO ENTRY IN EVERY COLUMN IN THAT ROW IN WHICH THERE IS A
C------NONZERO ENTRY IN ROW P(K) OF M, THEN DO NOT COMPUTE FILL-IN, JUST
C------USE THE COLUMN INDICES FOR THE ROW WHICH WAS TO HAVE BEEN MERGED
   40     IF (.NOT.CLIQUE) GO TO 50
          IJU(K) = IJU(JL(K)) + 1
          LUK = IU(JL(K)+1) - (IU(JL(K))+1)
          GO TO 170
C
C------MODIFY NONZERO STRUCTURE OF K-TH ROW BY COMPUTING FILL-IN
C------FOR EACH ROW I TO BE MERGED IN
   50     LMAX = 0
          IJU(K) = JUPTR
C
          I = K
   60     I = JL(I)
          IF (I.EQ.0) GO TO 100
C
C--------MERGE ROW I INTO K-TH ROW
          LUI = IU(I+1) - (IU(I)+1)
          JMIN = IJU(I) + 1
          JMAX = IJU(I) + LUI
          QM = K
C
          DO 80 J = JMIN,JMAX
              VJ = JU(J)
   70         M = QM
              QM = Q(M)
              IF (QM.LT.VJ) GO TO 70
              IF (QM.EQ.VJ) GO TO 80
              LUK = LUK + 1
              Q(M) = VJ
              Q(VJ) = QM
              QM = VJ
   80     CONTINUE
C
C--------REMEMBER LENGTH AND POSITION IN JU OF LONGEST ROW MERGED
          IF (LUI.LE.LMAX) GO TO 90
          LMAX = LUI
          IJU(K) = JMIN
C
   90     GO TO 60
C
C------IF THE K-TH ROW IS THE SAME LENGTH AS THE LONGEST ROW MERGED,
C------THEN USE THE COLUMN INDICES FOR THAT ROW
  100     IF (LUK.EQ.LMAX) GO TO 170
C
C------IF THE TAIL OF THE LAST ROW SAVED IN JU IS THE SAME AS THE HEAD
C------OF THE K-TH ROW, THEN OVERLAP THE TWO SETS OF COLUMN INDICES --
C--------SEARCH LAST ROW SAVED FOR FIRST NONZERO ENTRY IN K-TH ROW ...
          I = Q(K)
          IF (JUMIN.GT.JUPTR) GO TO 120
          DO 110 JMIN = JUMIN,JUPTR
              IF (JU(JMIN)-I) 110,130,120
  110     CONTINUE
  120     GO TO 150
C
C--------... AND THEN TEST WHETHER TAIL MATCHES HEAD OF K-TH ROW
  130     IJU(K) = JMIN
          DO 140 J = JMIN,JUPTR
              IF (JU(J).NE.I) GO TO 150
              I = Q(I)
              IF (I.GT.N) GO TO 170
  140     CONTINUE
          JUPTR = JMIN - 1
C
C------SAVE NONZERO STRUCTURE OF K-TH ROW IN JU
  150     I = K
          JUMIN = JUPTR + 1
          JUPTR = JUPTR + LUK
          IF (JUPTR.GT.JUMAX) GO TO 210
          DO 160 J = JUMIN,JUPTR
              I = Q(I)
              JU(J) = I
              MARK(I) = K
  160     CONTINUE
          IJU(K) = JUMIN
C
C------ADD K TO ROW LIST FOR FIRST NONZERO ELEMENT IN K-TH ROW
  170     IF (LUK.LE.1) GO TO 180
          I = JU(IJU(K))
          JL(K) = JL(I)
          JL(I) = K
C
  180     IU(K+1) = IU(K) + LUK
  190 CONTINUE
C
      FLAG = 0
      RETURN
C
C ** ERROR -- DUPLICATE ENTRY IN A
  200 FLAG = 2*N + P(K)
      RETURN
C ** ERROR -- INSUFFICIENT STORAGE FOR JU
  210 FLAG = 6*N + K
      RETURN

      END
C
C***********************************************************************
C  SNS -- SOLUTION OF SPARSE SYMMETRIC POSITIVE DEFINITE SYSTEM OF
C         LINEAR EQUATIONS  MX = B  GIVEN UT-D-U FACTORIZATION OF M
C***********************************************************************
      SUBROUTINE SNS(N,P,D,IJU,JU,IU,U,Z,B,TMP)
C       REAL  D(*), U(*),  Z(*), B(*),  TMP(*),  TMPK, SUM
C
C  ADDITIONAL PARAMETERS
C
C    TMP   - REAL ONE-DIMENSIONAL WORK ARRAY;  DIMENSION = N
C
C-----------------------------------------------------------------------
C
C----SET TMP TO PERMUTED B
C     .. Scalar Arguments ..
      INTEGER N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION B(*),D(*),TMP(*),U(*),Z(*)
      INTEGER IJU(*),IU(*),JU(*),P(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION SUM,TMPK
      INTEGER I,J,JMAX,JMIN,K,MU
C     ..
      DO 10 K = 1,N
          TMP(K) = B(P(K))
   10 CONTINUE
C
C----SOLVE  UT D Y = B  BY FORWARD SUBSTITUTION
      DO 40 K = 1,N
          TMPK = TMP(K)
          JMIN = IU(K)
          JMAX = IU(K+1) - 1
          IF (JMIN.GT.JMAX) GO TO 30
          MU = IJU(K) - JMIN
          DO 20 J = JMIN,JMAX
              TMP(JU(MU+J)) = TMP(JU(MU+J)) + U(J)*TMPK
   20     CONTINUE
   30     TMP(K) = TMPK*D(K)
   40 CONTINUE
C
C----SOLVE  U X = Y  BY BACK SUBSTITUTION
      K = N
      DO 70 I = 1,N
          SUM = TMP(K)
          JMIN = IU(K)
          JMAX = IU(K+1) - 1
          IF (JMIN.GT.JMAX) GO TO 60
          MU = IJU(K) - JMIN
          DO 50 J = JMIN,JMAX
              SUM = SUM + U(J)*TMP(JU(MU+J))
   50     CONTINUE
   60     TMP(K) = SUM
          Z(P(K)) = SUM
          K = K - 1
   70 CONTINUE
C
      RETURN

      END
C***********************************************************************
C***********************************************************************
C***********************************************************************
C  SRO -- SYMMETRIC REORDERING OF SPARSE SYMMETRIC MATRIX
C***********************************************************************
      SUBROUTINE SRO(N,IP,IA,JA,A,Q,R,DFLAG)
C
C  DESCRIPTION
C
C    THE NONZERO ENTRIES OF THE MATRIX M ARE ASSUMED TO BE STORED
C    SYMMETRICALLY IN (IA,JA,A) FORMAT (I.E., NOT BOTH M(I,J) AND M(J,I)
C    ARE STORED IF I NE J).
C
C    SRO DOES NOT REARRANGE THE ORDER OF THE ROWS, BUT DOES MOVE
C    NONZEROES FROM ONE ROW TO ANOTHER TO ENSURE THAT IF M(I,J) WILL BE
C    IN THE UPPER TRIANGLE OF M WITH RESPECT TO THE NEW ORDERING, THEN
C    M(I,J) IS STORED IN ROW I (AND THUS M(J,I) IS NOT STORED);  WHEREAS
C    IF M(I,J) WILL BE IN THE STRICT LOWER TRIANGLE OF M, THEN M(J,I) IS
C    STORED IN ROW J (AND THUS M(I,J) IS NOT STORED).
C
C
C  ADDITIONAL PARAMETERS
C
C    Q     - INTEGER ONE-DIMENSIONAL WORK ARRAY;  DIMENSION = N
C
C    R     - INTEGER ONE-DIMENSIONAL WORK ARRAY;  DIMENSION = NUMBER OF
C            NONZERO ENTRIES IN THE UPPER TRIANGLE OF M
C
C    DFLAG - LOGICAL VARIABLE;  IF DFLAG = .TRUE., THEN STORE NONZERO
C            DIAGONAL ELEMENTS AT THE BEGINNING OF THE ROW
C
C-----------------------------------------------------------------------
C
C       REAL  A(*),  AK
C
C
C--PHASE 1 -- FIND ROW IN WHICH TO STORE EACH NONZERO
C----INITIALIZE COUNT OF NONZEROES TO BE STORED IN EACH ROW
C     .. Scalar Arguments ..
      INTEGER N
      LOGICAL DFLAG
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(*)
      INTEGER IA(*),IP(*),JA(*),Q(*),R(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION AK
      INTEGER I,ILAST,J,JAK,JDUMMY,JMAX,JMIN,K
C     ..
      DO 10 I = 1,N
          Q(I) = 0
   10 CONTINUE
C
C----FOR EACH NONZERO ELEMENT A(J)
      DO 30 I = 1,N
          JMIN = IA(I)
          JMAX = IA(I+1) - 1
          IF (JMIN.GT.JMAX) GO TO 30
          DO 20 J = JMIN,JMAX
C
C--------FIND ROW (=R(J)) AND COLUMN (=JA(J)) IN WHICH TO STORE A(J) ...
              K = JA(J)
              IF (IP(K).LT.IP(I)) JA(J) = I
              IF (IP(K).GE.IP(I)) K = I
              R(J) = K
C
C--------... AND INCREMENT COUNT OF NONZEROES (=Q(R(J)) IN THAT ROW
              Q(K) = Q(K) + 1
   20     CONTINUE
   30 CONTINUE
C
C
C--PHASE 2 -- FIND NEW IA AND PERMUTATION TO APPLY TO (JA,A)
C----DETERMINE POINTERS TO DELIMIT ROWS IN PERMUTED (JA,A)
      DO 40 I = 1,N
          IA(I+1) = IA(I) + Q(I)
          Q(I) = IA(I+1)
   40 CONTINUE
C
C----DETERMINE WHERE EACH (JA(J),A(J)) IS STORED IN PERMUTED (JA,A)
C----FOR EACH NONZERO ELEMENT (IN REVERSE ORDER)
      ILAST = 0
      JMIN = IA(1)
      JMAX = IA(N+1) - 1
      J = JMAX
      DO 70 JDUMMY = JMIN,JMAX
          I = R(J)
          IF (.NOT.DFLAG .OR. JA(J).NE.I .OR. I.EQ.ILAST) GO TO 50
C
C------IF DFLAG, THEN PUT DIAGONAL NONZERO AT BEGINNING OF ROW
          R(J) = IA(I)
          ILAST = I
          GO TO 60
C
C------PUT (OFF-DIAGONAL) NONZERO IN LAST UNUSED LOCATION IN ROW
   50     Q(I) = Q(I) - 1
          R(J) = Q(I)
C
   60     J = J - 1
   70 CONTINUE
C
C
C--PHASE 3 -- PERMUTE (JA,A) TO UPPER TRIANGULAR FORM (WRT NEW ORDERING)
      DO 90 J = JMIN,JMAX
   80     IF (R(J).EQ.J) GO TO 90
          K = R(J)
          R(J) = R(K)
          R(K) = K
          JAK = JA(K)
          JA(K) = JA(J)
          JA(J) = JAK
          AK = A(K)
          A(K) = A(J)
          A(J) = AK
          GO TO 80

   90 CONTINUE
C
      RETURN

      END
C
C***********************************************************************
C***********************************************************************
C  MD -- MINIMUM DEGREE ALGORITHM (BASED ON ELEMENT MODEL)
C***********************************************************************
      SUBROUTINE MD(N,IA,JA,MAX,V,L,HEAD,LAST,NEXT,MARK,FLAG)
C
C  DESCRIPTION
C
C    MD FINDS A MINIMUM DEGREE ORDERING OF THE ROWS AND COLUMNS OF A
C    SYMMETRIC MATRIX M STORED IN (IA,JA,A) FORMAT.
C
C
C  ADDITIONAL PARAMETERS
C
C    MAX  - DECLARED DIMENSION OF THE ONE-DIMENSIONAL ARRAYS V AND L;
C           MAX MUST BE AT LEAST  N+2K,  WHERE K IS THE NUMBER OF
C           NONZEROES IN THE STRICT UPPER TRIANGLE OF M
C
C    V    - INTEGER ONE-DIMENSIONAL WORK ARRAY;  DIMENSION = MAX
C
C    L    - INTEGER ONE-DIMENSIONAL WORK ARRAY;  DIMENSION = MAX
C
C    HEAD - INTEGER ONE-DIMENSIONAL WORK ARRAY;  DIMENSION = N
C
C    LAST - INTEGER ONE-DIMENSIONAL ARRAY USED TO RETURN THE PERMUTATION
C           OF THE ROWS AND COLUMNS OF M CORRESPONDING TO THE MINIMUM
C           DEGREE ORDERING;  DIMENSION = N
C
C    NEXT - INTEGER ONE-DIMENSIONAL ARRAY USED TO RETURN THE INVERSE OF
C           THE PERMUTATION RETURNED IN LAST;  DIMENSION = N
C
C    MARK - INTEGER ONE-DIMENSIONAL WORK ARRAY (MAY BE THE SAME AS V);
C           DIMENSION = N
C
C    FLAG - INTEGER ERROR FLAG;  VALUES AND THEIR MEANINGS ARE -
C             0      NO ERRORS DETECTED
C             11N+1  INSUFFICIENT STORAGE IN MD
C
C
C  DEFINITIONS OF INTERNAL PARAMETERS
C
C    ---------+---------------------------------------------------------
C    V(S)     \ VALUE FIELD OF LIST ENTRY
C    ---------+---------------------------------------------------------
C    L(S)     \ LINK FIELD OF LIST ENTRY  (0 => END OF LIST)
C    ---------+---------------------------------------------------------
C    L(VI)    \ POINTER TO ELEMENT LIST OF UNELIMINATED VERTEX VI
C    ---------+---------------------------------------------------------
C    L(EJ)    \ POINTER TO BOUNDARY LIST OF ACTIVE ELEMENT EJ
C    ---------+---------------------------------------------------------
C    HEAD(D)  \ VJ => VJ HEAD OF D-LIST D
C             \  0 => NO VERTEX IN D-LIST D
C
C
C             \                  VI UNELIMINATED VERTEX
C             \          VI IN EK           \       VI NOT IN EK
C    ---------+-----------------------------+---------------------------
C    NEXT(VI) \ UNDEFINED BUT NONNEGATIVE   \ VJ => VJ NEXT IN D-LIST
C             \                             \  0 => VI TAIL OF D-LIST
C    ---------+-----------------------------+---------------------------
C    LAST(VI) \ (NOT SET UNTIL MDP)         \ -D => VI HEAD OF D-LIST D
C             \-VK => COMPUTE DEGREE        \ VJ => VJ LAST IN D-LIST
C             \ EJ => VI PROTOTYPE OF EJ    \  0 => VI NOT IN ANY D-LIST
C             \  0 => DO NOT COMPUTE DEGREE C    ---------+-----------------------------+---------------------------
C    MARK(VI) \ MARK(VK)                    \ NONNEGATIVE TAG < MARK(VK)
C
C
C             \                   VI ELIMINATED VERTEX
C             \      EI ACTIVE ELEMENT      \           OTHERWISE
C    ---------+-----------------------------+---------------------------
C    NEXT(VI) \ -J => VI WAS J-TH VERTEX    \ -J => VI WAS J-TH VERTEX
C             \       TO BE ELIMINATED      \       TO BE ELIMINATED
C    ---------+-----------------------------+---------------------------
C    LAST(VI) \  M => SIZE OF EI = M        \ UNDEFINED
C    ---------+-----------------------------+---------------------------
C    MARK(VI) \ -M => OVERLAP COUNT OF EI   \ UNDEFINED
C             \       WITH EK = M           C             \ OTHERWISE NONNEGATIVE TAG   C             \       < MARK(VK)            C
C-----------------------------------------------------------------------
C
C
C----INITIALIZATION
C     .. Scalar Arguments ..
      INTEGER FLAG,MAX,N
C     ..
C     .. Array Arguments ..
      INTEGER HEAD(*),IA(*),JA(*),L(*),LAST(*),MARK(*),NEXT(*),V(*)
C     ..
C     .. Local Scalars ..
      INTEGER DMIN,EK,K,TAG,TAIL,VK
C     ..
C     .. External Subroutines ..
      EXTERNAL MDI,MDM,MDP,MDU
C     ..
C     .. Equivalences ..
      EQUIVALENCE (VK,EK)
C     ..
      TAG = 0
      CALL MDI(N,IA,JA,MAX,V,L,HEAD,LAST,NEXT,MARK,TAG,FLAG)
      IF (FLAG.NE.0) RETURN
C
      K = 0
      DMIN = 1
C
C----WHILE  K < N  DO
   10 IF (K.GE.N) GO TO 40
C
C------SEARCH FOR VERTEX OF MINIMUM DEGREE
   20 IF (HEAD(DMIN).GT.0) GO TO 30
      DMIN = DMIN + 1
      GO TO 20
C
C------REMOVE VERTEX VK OF MINIMUM DEGREE FROM DEGREE LIST
   30 VK = HEAD(DMIN)
      HEAD(DMIN) = NEXT(VK)
      IF (HEAD(DMIN).GT.0) LAST(HEAD(DMIN)) = -DMIN
C
C------NUMBER VERTEX VK, ADJUST TAG, AND TAG VK
      K = K + 1
      NEXT(VK) = -K
      LAST(EK) = DMIN - 1
      TAG = TAG + LAST(EK)
      MARK(VK) = TAG
C
C------FORM ELEMENT EK FROM UNELIMINATED NEIGHBORS OF VK
      CALL MDM(VK,TAIL,V,L,LAST,NEXT,MARK)
C
C------PURGE INACTIVE ELEMENTS AND DO MASS ELIMINATION
      CALL MDP(K,EK,TAIL,V,L,HEAD,LAST,NEXT,MARK)
C
C------UPDATE DEGREES OF UNELIMINATED VERTICES IN EK
      CALL MDU(EK,DMIN,V,L,HEAD,LAST,NEXT,MARK)
C
      GO TO 10
C
C----GENERATE INVERSE PERMUTATION FROM PERMUTATION
   40 DO 50 K = 1,N
          NEXT(K) = -NEXT(K)
          LAST(NEXT(K)) = K
   50 CONTINUE
C
      RETURN

      END
C
C***********************************************************************
C  MDI -- INITIALIZATION
C***********************************************************************
      SUBROUTINE MDI(N,IA,JA,MAX,V,L,HEAD,LAST,NEXT,MARK,TAG,FLAG)
C
C----INITIALIZE DEGREES, ELEMENT LISTS, AND DEGREE LISTS
C     .. Scalar Arguments ..
      INTEGER FLAG,MAX,N,TAG
C     ..
C     .. Array Arguments ..
      INTEGER HEAD(*),IA(*),JA(*),L(*),LAST(*),MARK(*),NEXT(*),V(*)
C     ..
C     .. Local Scalars ..
      INTEGER DVI,J,JMAX,JMIN,SFS,VI,VJ
C     ..
      DO 10 VI = 1,N
          MARK(VI) = 1
          L(VI) = 0
          HEAD(VI) = 0
   10 CONTINUE
      SFS = N + 1
C
C----CREATE NONZERO STRUCTURE
C----FOR EACH NONZERO ENTRY A(VI,VJ) IN STRICT UPPER TRIANGLE
      DO 30 VI = 1,N
          JMIN = IA(VI)
          JMAX = IA(VI+1) - 1
          IF (JMIN.GT.JMAX) GO TO 30
          DO 20 J = JMIN,JMAX
              VJ = JA(J)
              IF (VI.GE.VJ) GO TO 20
              IF (SFS.GE.MAX) GO TO 50
C
C------ENTER VJ IN ELEMENT LIST FOR VI
              MARK(VI) = MARK(VI) + 1
              V(SFS) = VJ
              L(SFS) = L(VI)
              L(VI) = SFS
              SFS = SFS + 1
C
C------ENTER VI IN ELEMENT LIST FOR VJ
              MARK(VJ) = MARK(VJ) + 1
              V(SFS) = VI
              L(SFS) = L(VJ)
              L(VJ) = SFS
              SFS = SFS + 1
   20     CONTINUE
   30 CONTINUE
C
C----CREATE DEGREE LISTS AND INITIALIZE MARK VECTOR
      DO 40 VI = 1,N
          DVI = MARK(VI)
          NEXT(VI) = HEAD(DVI)
          HEAD(DVI) = VI
          LAST(VI) = -DVI
          IF (NEXT(VI).GT.0) LAST(NEXT(VI)) = VI
          MARK(VI) = TAG
   40 CONTINUE
C
      RETURN
C
C ** ERROR -- INSUFFICIENT STORAGE
   50 FLAG = 9*N + VI
      RETURN

      END
C
C***********************************************************************
C  MDM -- FORM ELEMENT FROM UNELIMINATED NEIGHBORS OF VK
C***********************************************************************
      SUBROUTINE MDM(VK,TAIL,V,L,LAST,NEXT,MARK)
C
C----INITIALIZE TAG AND LIST OF UNELIMINATED NEIGHBORS
C     .. Scalar Arguments ..
      INTEGER TAIL,VK
C     ..
C     .. Array Arguments ..
      INTEGER L(*),LAST(*),MARK(*),NEXT(*),V(*)
C     ..
C     .. Local Scalars ..
      INTEGER B,BLP,BLPMAX,ES,LB,LS,S,TAG,VB,VS
C     ..
C     .. Equivalences ..
      EQUIVALENCE (VS,ES)
C     ..
      TAG = MARK(VK)
      TAIL = VK
C
C----FOR EACH VERTEX/ELEMENT VS/ES IN ELEMENT LIST OF VK
      LS = L(VK)
   10 S = LS
      IF (S.EQ.0) GO TO 50
      LS = L(S)
      VS = V(S)
      IF (NEXT(VS).LT.0) GO TO 20
C
C------IF VS IS UNELIMINATED VERTEX, THEN TAG AND APPEND TO LIST OF
C------UNELIMINATED NEIGHBORS
      MARK(VS) = TAG
      L(TAIL) = S
      TAIL = S
      GO TO 40
C
C------IF ES IS ACTIVE ELEMENT, THEN ...
C--------FOR EACH VERTEX VB IN BOUNDARY LIST OF ELEMENT ES
   20 LB = L(ES)
      BLPMAX = LAST(ES)
      DO 30 BLP = 1,BLPMAX
          B = LB
          LB = L(B)
          VB = V(B)
C
C----------IF VB IS UNTAGGED VERTEX, THEN TAG AND APPEND TO LIST OF
C----------UNELIMINATED NEIGHBORS
          IF (MARK(VB).GE.TAG) GO TO 30
          MARK(VB) = TAG
          L(TAIL) = B
          TAIL = B
   30 CONTINUE
C
C--------MARK ES INACTIVE
      MARK(ES) = TAG
C
   40 GO TO 10
C
C----TERMINATE LIST OF UNELIMINATED NEIGHBORS
   50 L(TAIL) = 0
C
      RETURN

      END
C
C***********************************************************************
C  MDP -- PURGE INACTIVE ELEMENTS AND DO MASS ELIMINATION
C***********************************************************************
      SUBROUTINE MDP(K,EK,TAIL,V,L,HEAD,LAST,NEXT,MARK)
C
C----INITIALIZE TAG
C     .. Scalar Arguments ..
      INTEGER EK,K,TAIL
C     ..
C     .. Array Arguments ..
      INTEGER HEAD(*),L(*),LAST(*),MARK(*),NEXT(*),V(*)
C     ..
C     .. Local Scalars ..
      INTEGER ES,EVI,FREE,I,ILP,ILPMAX,LI,LS,LVI,S,TAG,VI
C     ..
      TAG = MARK(EK)
C
C----FOR EACH VERTEX VI IN EK
      LI = EK
      ILPMAX = LAST(EK)
      IF (ILPMAX.LE.0) GO TO 120
      DO 110 ILP = 1,ILPMAX
          I = LI
          LI = L(I)
          VI = V(LI)
C
C------REMOVE VI FROM DEGREE LIST
          IF (LAST(VI).EQ.0) GO TO 30
          IF (LAST(VI).GT.0) GO TO 10
          HEAD(-LAST(VI)) = NEXT(VI)
          GO TO 20

   10     NEXT(LAST(VI)) = NEXT(VI)
   20     IF (NEXT(VI).GT.0) LAST(NEXT(VI)) = LAST(VI)
C
C------REMOVE INACTIVE ITEMS FROM ELEMENT LIST OF VI
   30     LS = VI
   40     S = LS
          LS = L(S)
          IF (LS.EQ.0) GO TO 60
          ES = V(LS)
          IF (MARK(ES).LT.TAG) GO TO 50
          FREE = LS
          L(S) = L(LS)
          LS = S
   50     GO TO 40
C
C------IF VI IS INTERIOR VERTEX, THEN REMOVE FROM LIST AND ELIMINATE
   60     LVI = L(VI)
          IF (LVI.NE.0) GO TO 70
          L(I) = L(LI)
          LI = I
C
          K = K + 1
          NEXT(VI) = -K
          LAST(EK) = LAST(EK) - 1
          GO TO 110
C
C------ELSE ...
C--------CLASSIFY VERTEX VI
   70     IF (L(LVI).NE.0) GO TO 90
          EVI = V(LVI)
          IF (NEXT(EVI).GE.0) GO TO 90
          IF (MARK(EVI).LT.0) GO TO 80
C
C----------IF VI IS PROTOTYPE VERTEX, THEN MARK AS SUCH, INITIALIZE
C----------OVERLAP COUNT FOR CORRESPONDING ELEMENT, AND MOVE VI TO END
C----------OF BOUNDARY LIST
          LAST(VI) = EVI
          MARK(EVI) = -1
          L(TAIL) = LI
          TAIL = LI
          L(I) = L(LI)
          LI = I
          GO TO 100
C
C----------ELSE IF VI IS DUPLICATE VERTEX, THEN MARK AS SUCH AND ADJUST
C----------OVERLAP COUNT FOR CORRESPONDING ELEMENT
   80     LAST(VI) = 0
          MARK(EVI) = MARK(EVI) - 1
          GO TO 100
C
C----------ELSE MARK VI TO COMPUTE DEGREE
   90     LAST(VI) = -EK
C
C--------INSERT EK IN ELEMENT LIST OF VI
  100     V(FREE) = EK
          L(FREE) = L(VI)
          L(VI) = FREE
  110 CONTINUE
C
C----TERMINATE BOUNDARY LIST
  120 L(TAIL) = 0
C
      RETURN

      END
C
C***********************************************************************
C  MDU -- UPDATE DEGREES OF UNELIMINATED VERTICES IN EK
C***********************************************************************
      SUBROUTINE MDU(EK,DMIN,V,L,HEAD,LAST,NEXT,MARK)
C
C----INITIALIZE TAG
C     .. Scalar Arguments ..
      INTEGER DMIN,EK
C     ..
C     .. Array Arguments ..
      INTEGER HEAD(*),L(*),LAST(*),MARK(*),NEXT(*),V(*)
C     ..
C     .. Local Scalars ..
      INTEGER B,BLP,BLPMAX,DVI,ES,EVI,I,ILP,ILPMAX,S,TAG,VB,VI,VS
C     ..
C     .. Equivalences ..
      EQUIVALENCE (VS,ES)
C     ..
      TAG = MARK(EK) - LAST(EK)
C
C----FOR EACH VERTEX VI IN EK
      I = EK
      ILPMAX = LAST(EK)
      IF (ILPMAX.LE.0) GO TO 110
      DO 100 ILP = 1,ILPMAX
          I = L(I)
          VI = V(I)
          IF (LAST(VI)) 10,100,80
C
C------IF VI NEITHER PROTOTYPE NOR DUPLICATE VERTEX, THEN MERGE ELEMENTS
C------TO COMPUTE DEGREE
   10     TAG = TAG + 1
          DVI = LAST(EK)
C
C--------FOR EACH VERTEX/ELEMENT VS/ES IN ELEMENT LIST OF VI
          S = L(VI)
   20     S = L(S)
          IF (S.EQ.0) GO TO 90
          VS = V(S)
          IF (NEXT(VS).LT.0) GO TO 30
C
C----------IF VS IS UNELIMINATED VERTEX, THEN TAG AND ADJUST DEGREE
          MARK(VS) = TAG
          DVI = DVI + 1
          GO TO 50
C
C----------IF ES IS ACTIVE ELEMENT, THEN EXPAND
C------------CHECK FOR OUTMATCHED VERTEX
   30     IF (MARK(ES).LT.0) GO TO 60
C
C------------FOR EACH VERTEX VB IN ES
          B = ES
          BLPMAX = LAST(ES)
          DO 40 BLP = 1,BLPMAX
              B = L(B)
              VB = V(B)
C
C--------------IF VB IS UNTAGGED, THEN TAG AND ADJUST DEGREE
              IF (MARK(VB).GE.TAG) GO TO 40
              MARK(VB) = TAG
              DVI = DVI + 1
   40     CONTINUE
C
   50     GO TO 20
C
C------ELSE IF VI IS OUTMATCHED VERTEX, THEN ADJUST OVERLAPS BUT DO NOT
C------COMPUTE DEGREE
   60     LAST(VI) = 0
          MARK(ES) = MARK(ES) - 1
   70     S = L(S)
          IF (S.EQ.0) GO TO 100
          ES = V(S)
          IF (MARK(ES).LT.0) MARK(ES) = MARK(ES) - 1
          GO TO 70
C
C------ELSE IF VI IS PROTOTYPE VERTEX, THEN CALCULATE DEGREE BY
C------INCLUSION/EXCLUSION AND RESET OVERLAP COUNT
   80     EVI = LAST(VI)
          DVI = LAST(EK) + LAST(EVI) + MARK(EVI)
          MARK(EVI) = 0
C
C------INSERT VI IN APPROPRIATE DEGREE LIST
   90     NEXT(VI) = HEAD(DVI)
          HEAD(DVI) = VI
          LAST(VI) = -DVI
          IF (NEXT(VI).GT.0) LAST(NEXT(VI)) = VI
          IF (DVI.LT.DMIN) DMIN = DVI
C
  100 CONTINUE
C
  110 RETURN

      END
C***********************************************************************
C***********************************************************************
C     Modified Routines  of the Yale Sparse Matrix Package (YSMP):
C
C         SDRVMD   - a modified form of the SDRV driver.
C         SNFMOD   - a modified form of the SNF routine.
C
C               copyright (c) 1990 by Tamar Schlick
C
C***********************************************************************
C        Yale's SDRV solves a linear system for (symmetric) positive-
C   definite matrices. It calls the following routines:
C                 SSF (for symbolic factorization)
C                 SNF (for numerical factorization) and
C                 SNS (for numerical solution).
C        Our goal is to solve large sparse symmetric linear systems for
C   matrices that are not necessarily pos-def. Thus, we
C   replace SNF by SNFMOD so that a modified-Cholesky (MCF), rather than
C   a Cholesky, factorization is performed. In SDRV, we replace
C   the statement "CALL SNF" by "CALL SNFMOD".
C        In Yale's SDRV, the diagnostic parameter FLAG is set to zero in
C   the non pos-def case (and control is returned to the main program).
C   Here, instead, we set FLAG in SNFMOD to a negative integer if the
C   matrix is not sufficiently pos-def.
C   Specifically, FLAG is set to minus
C   the position of the diagonal element in the original matrix whose
C   modification was of the largest magnitude. Recall that in MCF we
C   produce matrices E,D, and U so that
C                 M + E = UT-D-U  where E and D are
C   diagonal and U is unit upper-triangular. FLAG records the index k
C   for the largest modification in E: ( E(P(k)) = max over i {E(i)} ).
C
C   All modifications to the original YSMP code are indicated.
C                                                                1/15/81
C***********************************************************************
C  SDRV -- DRIVER FOR SPARSE SYMMETRIC POSITIVE DEFINITE MATRIX ROUTINES
C***********************************************************************

C ====================  change #1  (replacement) =====================1
C WAS:  SUBROUTINE SDRV
C =====================================================================
C ==================================================================end
      SUBROUTINE SDRVMD(N,P,IP,IA,JA,A,B,Z,NSP,ISP,RSP,ESP,PATH,FLAG)
C
C  DESCRIPTION
C
C ====================  change #2  (replacement) =====================2
C WAS: SDRV SOLVES SPARSE SYMMETRIC POSITIVE DEFINITE SYSTEMS OF LINEAR
C =====================================================================
C    SDRVMD SOLVES SPARSE SYMMETRIC SYSTEMS OF LINEAR
C ==================================================================end
C    EQUATIONS.  THE SOLUTION PROCESS IS DIVIDED INTO THREE STAGES --
C
C      SSF - THE COEFFICIENT MATRIX M IS FACTORED SYMBOLICALLY TO
C            DETERMINE WHERE FILLIN WILL OCCUR DURING THE NUMERIC
C            FACTORIZATION.
C
C ====================  change #3  (replacement) =====================3
C WAS: SNF - M IS FACTORED NUMERICALLY INTO THE PRODUCT UT-D-U, WHERE
C =====================================================================
C      SNFMOD - M+E IS FACTORED NUMERICALLY BY THE GILL/MURRAY/WRIGHT
C            MODIFIED CHOLESKY FACTORIZATION: M + E = UT-D-U, WHERE
C            E IS DIAGONAL,
C ==================================================================end
C            D IS DIAGONAL AND U IS UNIT UPPER TRIANGULAR.
C
C      SNS - THE LINEAR SYSTEM  MX = B  IS SOLVED USING THE UT-D-U
C            FACTORIZATION FROM SNF.
C
C    FOR SEVERAL SYSTEMS WITH THE SAME COEFFICIENT MATRIX, SSF AND SNF
C    NEED BE DONE ONLY ONCE (FOR THE FIRST SYSTEM);  THEN SNS IS DONE
C    ONCE FOR EACH ADDITIONAL RIGHT-HAND SIDE.  FOR SEVERAL SYSTEMS
C    WHOSE COEFFICIENT MATRICES HAVE THE SAME NONZERO STRUCTURE, SSF
C    NEED BE DONE ONLY ONCE (FOR THE FIRST SYSTEM);  THEN SNF AND SNS
C    ARE DONE ONCE FOR EACH ADDITIONAL SYSTEM.
C
C
C  STORAGE OF SPARSE MATRICES
C
C    THE NONZERO ENTRIES OF THE MATRIX M ARE STORED ROW-BY-ROW IN THE
C    ARRAY A.  TO IDENTIFY THE INDIVIDUAL NONZERO ENTRIES IN EACH ROW,
C    WE NEED TO KNOW IN WHICH COLUMN EACH ENTRY LIES.  THESE COLUMN
C    INDICES ARE STORED IN THE ARRAY JA;  I.E., IF  A(K) = M(I,J),  THEN
C    JA(K) = J.  TO IDENTIFY THE INDIVIDUAL ROWS, WE NEED TO KNOW WHERE
C    EACH ROW STARTS.  THESE ROW POINTERS ARE STORED IN THE ARRAY IA;
C    I.E., IF M(I,J) IS THE FIRST NONZERO ENTRY (STORED) IN THE I-TH ROW
C    AND  A(K) = M(I,J),  THEN  IA(I) = K.  MOREOVER, IA(N+1) POINTS TO
C    THE FIRST LOCATION FOLLOWING THE LAST ELEMENT IN THE LAST ROW.
C    THUS, THE NUMBER OF ENTRIES IN THE I-TH ROW IS  IA(I+1) - IA(I),
C    THE NONZERO ENTRIES IN THE I-TH ROW ARE STORED CONSECUTIVELY IN
C
C            A(IA(I)),  A(IA(I)+1),  ..., A(IA(I+1)-1),
C
C    AND THE CORRESPONDING COLUMN INDICES ARE STORED CONSECUTIVELY IN
C
C            JA(IA(I)), JA(IA(I)+1), ..., JA(IA(I+1)-1).
C
C    SINCE THE COEFFICIENT MATRIX IS SYMMETRIC, ONLY THE NONZERO ENTRIES
C    IN THE UPPER TRIANGLE NEED BE STORED, FOR EXAMPLE, THE MATRIX
C
C             ( 1  0  2  3  0 )
C             ( 0  4  0  0  0 )
C         M = ( 2  0  5  6  0 )
C             ( 3  0  6  7  8 )
C             ( 0  0  0  8  9 )
C
C    COULD BE STORED AS
C
C            \ 1  2  3  4  5  6  7  8  9 10 11 12 13
C         ---+--------------------------------------
C         IA \ 1  4  5  8 12 14
C         JA \ 1  3  4  2  1  3  4  1  3  4  5  4  5
C          A \ 1  2  3  4  2  5  6  3  6  7  8  8  9
C
C    OR (SYMMETRICALLY) AS
C
C            \ 1  2  3  4  5  6  7  8  9
C         ---+--------------------------
C         IA \ 1  4  5  7  9 10
C         JA \ 1  3  4  2  3  4  4  5  5
C          A \ 1  2  3  4  5  6  7  8  9          .
C
C
C  REORDERING THE ROWS AND COLUMNS OF M
C
C    A SYMMETRIC PERMUTATION OF THE ROWS AND COLUMNS OF THE COEFFICIENT
C    MATRIX M (E.G., WHICH REDUCES FILLIN OR ENHANCES NUMERICAL
C    STABILITY) MUST BE SPECIFIED.  THE SOLUTION Z IS RETURNED IN THE
C    ORIGINAL ORDER.
C
C    TO SPECIFY THE TRIVIAL ORDERING (I.E., THE IDENTITY PERMUTATION),
C    SET  P(I) = IP(I) = I,  I=1,...,N.  IN THIS CASE, P AND IP CAN BE
C    THE SAME ARRAY.
C
C    IF A NONTRIVIAL ORDERING (I.E., NOT THE IDENTITY PERMUTATION) IS
C    SPECIFIED AND M IS STORED SYMMETRICALLY (I.E., NOT BOTH M(I,J) AND
C    M(J,I) ARE STORED FOR I NE J), THEN ODRV SHOULD BE CALLED (WITH
C    PATH = 3 OR 5) TO SYMMETRICALLY REORDER (IA,JA,A) BEFORE CALLING
C    SDRV.  THIS IS TO ENSURE THAT IF M(I,J) WILL BE IN THE UPPER
C    TRIANGLE OF M WITH RESPECT TO THE NEW ORDERING, THEN M(I,J) IS
C    STORED IN ROW I (AND THUS M(J,I) IS NOT STORED);  WHEREAS IF M(I,J)
C    WILL BE IN THE STRICT LOWER TRIANGLE OF M, THEN M(J,I) IS STORED IN
C    ROW J (AND THUS M(I,J) IS NOT STORED).
C
C
C  PARAMETERS
C
C    N    - NUMBER OF VARIABLES/EQUATIONS
C
C    P    - INTEGER ONE-DIMENSIONAL ARRAY SPECIFYING A PERMUTATION OF
C           THE ROWS AND COLUMNS OF M;  DIMENSION = N
C
C    IP   - INTEGER ONE-DIMENSIONAL ARRAY CONTAINING THE INVERSE OF THE
C           PERMUTATION SPECIFIED IN P;  I.E., IP(P(I)) = I, I=1,...,N;
C           DIMENSION = N
C
C    IA   - INTEGER ONE-DIMENSIONAL ARRAY CONTAINING POINTERS TO DELIMIT
C           ROWS IN JA AND A;  DIMENSION = N+1
C
C    JA   - INTEGER ONE-DIMENSIONAL ARRAY CONTAINING THE COLUMN INDICES
C           CORRESPONDING TO THE ELEMENTS OF A;  DIMENSION = NUMBER OF
C           NONZERO ENTRIES IN M STORED
C
C    A    - REAL ONE-DIMENSIONAL ARRAY CONTAINING THE NONZERO ENTRIES IN
C           THE COEFFICIENT MATRIX M, STORED BY ROWS;  DIMENSION =
C           NUMBER OF NONZERO ENTRIES IN M STORED
C
C    B    - REAL ONE-DIMENSIONAL ARRAY CONTAINING THE RIGHT-HAND SIDE B;
C           B AND Z CAN BE THE SAME ARRAY;  DIMENSION = N
C
C    Z    - REAL ONE-DIMENSIONAL ARRAY CONTAINING THE SOLUTION X;  Z AND
C           B CAN BE THE SAME ARRAY;  DIMENSION = N
C
C    NSP  - DECLARED DIMENSION OF THE ONE-DIMENSIONAL ARRAYS ISP AND
C           RSP;  NSP MUST BE (SUBSTANTIALLY) LARGER THAN  3N+2K,  WHERE
C           K = NUMBER OF NONZERO ENTRIES IN THE UPPER TRIANGLE OF M
C
C    ISP  - INTEGER ONE-DIMENSIONAL ARRAY USED FOR WORKING STORAGE;  ISP
C           AND RSP SHOULD BE EQUIVALENCED;  DIMENSION = NSP
C
C    RSP  - REAL ONE-DIMENSIONAL ARRAY USED FOR WORKING STORAGE;  RSP
C           AND ISP SHOULD BE EQUIVALENCED;  DIMENSION = NSP
C
C    ESP  - INTEGER VARIABLE;  IF SUFFICIENT STORAGE WAS AVAILABLE TO
C           PERFORM THE SYMBOLIC FACTORIZATION (SSF), THEN ESP IS SET TO
C           THE AMOUNT OF EXCESS STORAGE PROVIDED (NEGATIVE IF
C           INSUFFICIENT STORAGE WAS AVAILABLE TO PERFORM THE NUMERIC
C           FACTORIZATION (SNF))
C
C    PATH - INTEGER PATH SPECIFICATION;  VALUES AND THEIR MEANINGS ARE -
C             1  PERFORM SSF, SNF, AND SNS
C             2  PERFORM SNF AND SNS (ISP/RSP IS ASSUMED TO HAVE BEEN
C                  SET UP IN AN EARLIER CALL TO SDRV (FOR SSF))
C             3  PERFORM SNS ONLY (ISP/RSP IS ASSUMED TO HAVE BEEN SET
C                  UP IN AN EARLIER CALL TO SDRV (FOR SSF AND SNF))
C             4  PERFORM SSF
C             5  PERFORM SSF AND SNF
C             6  PERFORM SNF ONLY (ISP/RSP IS ASSUMED TO HAVE BEEN SET
C                  UP IN AN EARLIER CALL TO SDRV (FOR SSF))
C
C    FLAG - INTEGER ERROR FLAG;  VALUES AND THEIR MEANINGS ARE -
C
C               0     NO ERRORS DETECTED
C              2N+K   DUPLICATE ENTRY IN A  --  ROW = K
C              6N+K   INSUFFICIENT STORAGE IN SSF  --  ROW = K
C              7N+1   INSUFFICIENT STORAGE IN SNF
C              8N+K   ZERO PIVOT  --  ROW = K
C             10N+1   INSUFFICIENT STORAGE IN SDRV
C             11N+1   ILLEGAL PATH SPECIFICATION
C
C ====================  change #4  (insertion) =======================4
C              <0     MATRIX NOT SUFF. POS-DEF (detected in SNFMOD)
C                     FLAG IS SET TO MINUS THE INDEX (IN THE ORIGINAL
C                     MATRIX) OF THE LARGEST ADDITION in E
C ==================================================================end
C
C
C  CONVERSION FROM REAL TO DOUBLE PRECISION
C
C    CHANGE THE REAL DECLARATIONS IN SDRV, SNF, AND SNS TO DOUBLE
C    PRECISION DECLARATIONS;  AND CHANGE THE VALUE IN THE DATA STATEMENT
C    FOR THE INTEGER VARIABLE RATIO (IN SDRV) FROM 1 TO 2.
C
C  NOTE: FOR CRAY, SET RATIO to 1!
C-----------------------------------------------------------------------
C
C       REAL  A(*),  B(*),  Z(*),  RSP(*)

C       DATA  RATIO/1/
C     .. Scalar Arguments ..
      INTEGER ESP,FLAG,N,NSP,PATH
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(*),B(*),RSP(*),Z(*)
      INTEGER IA(*),IP(*),ISP(*),JA(*),P(*)
C     ..
C     .. Local Scalars ..
      INTEGER D,IJU,IL,IU,JL,JU,JUMAX,MARK,Q,RATIO,TMP,U,UMAX
C     ..
C     .. External Subroutines ..
      EXTERNAL SNFMOD,SNS,SSF
C     ..
C     .. Data statements ..
      DATA RATIO/2/
C     ..
C
C----VALIDATE PATH SPECIFICATION
      IF (PATH.LT.1 .OR. 6.LT.PATH) GO TO 60
C
C----ALLOCATE STORAGE AND FACTOR M SYMBOLICALLY TO DETERMINE FILL-IN
      IJU = 1
      IU = IJU + N
      JL = IU + N + 1
      JU = JL + N
      Q = (NSP+1) - N
      MARK = Q - N
      JUMAX = MARK - JU
C
      IF ((PATH-1)* (PATH-4)* (PATH-5).NE.0) GO TO 10
      IF (JUMAX.LE.0) GO TO 50
      CALL SSF(N,P,IP,IA,JA,ISP(IJU),ISP(JU),ISP(IU),JUMAX,ISP(Q),
     +         ISP(MARK),ISP(JL),FLAG)
      IF (FLAG.NE.0) GO TO 40
C
C----ALLOCATE STORAGE AND FACTOR M NUMERICALLY
   10 IL = JU + ISP(IJU+ (N-1))
      TMP = ((IL-1)+ (RATIO-1))/RATIO + 1
      D = TMP + N
      U = D + N
      UMAX = (NSP+1) - U
      ESP = UMAX - (ISP(IU+N)-1)
C
      IF ((PATH-1)* (PATH-2)* (PATH-5)* (PATH-6).NE.0) GO TO 20
      IF (UMAX.LE.0) GO TO 50

C ====================  change #5  (replacement) =====================5
C  WAS:   CALL SNF
C ==================================================================end

      CALL SNFMOD(N,P,IP,IA,JA,A,RSP(D),ISP(IJU),ISP(JU),ISP(IU),RSP(U),
     +            UMAX,ISP(IL),ISP(JL),FLAG)

C ====================  change #6  (replacement) =====================6
C  WAS:         IF (FLAG.NE.0)  GO TO 100
C ==================================================================end

      IF (FLAG.GT.0) GO TO 40
C
C----SOLVE SYSTEM OF LINEAR EQUATIONS  MX = B
   20 IF ((PATH-1)* (PATH-2)* (PATH-3).NE.0) GO TO 30
      IF (UMAX.LE.0) GO TO 50
      CALL SNS(N,P,RSP(D),ISP(IJU),ISP(JU),ISP(IU),RSP(U),Z,B,RSP(TMP))
C
   30 RETURN
C
C ** ERROR -- ERROR DETECTED IN SSF, SNF, OR SNS
   40 RETURN
C ** ERROR -- INSUFFICIENT STORAGE
   50 FLAG = 10*N + 1
      RETURN
C ** ERROR -- ILLEGAL PATH SPECIFICATION
   60 FLAG = 11*N + 1
      RETURN

      END
C
C***********************************************************************
C***********************************************************************
C NUMERICAL FACTORIZATION OF SYMMETRIC MATRICES
C***********************************************************************
C
C ====================  change #1  (replacement) =====================1
C WAS:
C  SNF -- NUMERICAL UT-D-U FACTORIZATION OF SPARSE SYMMETRIC POSITIVE
C         DEFINITE MATRIX
C       SUBROUTINE  SNF
C =====================================================================
C
C  SNFMOD -- NUMERICAL FACTORIZATION OF SPARSE SYMMETRIC MATRICES M BY
C         THE GILL/MURRAY/WRIGHT MODIFIED CHOLESKY FACTORIZATION (GMW
C         MCF) WITHOUT PIVOTING.  THE FACTORIZATION PRODUCES U,D, AND
C         E SO THAT   M + E = UT-D-U,  WHERE  E AND D ARE DIAGONAL
C         MATRICES. THIS ROUTINE IS A MODIFICATION OF THE YSMP
C         routine SNF. ALL CHANGES ARE INDICATED.
C
C ==================================================================end
      SUBROUTINE SNFMOD(N,P,IP,IA,JA,A,D,IJU,JU,IU,U,UMAX,IL,JL,FLAG)
C
C  ADDITIONAL PARAMETERS
C
C    IL    - INTEGER ONE-DIMENSIONAL WORK ARRAY;  DIMENSION = N
C
C    JL    - INTEGER ONE-DIMENSIONAL WORK ARRAY;  DIMENSION = N
C
C
C  DEFINITIONS OF INTERNAL PARAMETERS (DURING K-TH STAGE OF ELIMINATION)
C
C    (D(I),I=K,N) CONTAINS THE K-TH ROW OF U (EXPANDED)
C
C    IL(I) POINTS TO THE FIRST NONZERO ELEMENT IN COLUMNS K,...,N OF
C      ROW I OF U
C
C    JL CONTAINS LISTS OF ROWS TO BE ADDED TO UNELIMINATED ROWS --
C      I GE K => JL(I) IS THE FIRST ROW TO BE ADDED TO ROW I
C      I LT K => JL(I) IS THE ROW FOLLOWING ROW I IN SOME LIST OF ROWS
C      IN EITHER CASE, JL(I) = 0 INDICATES THE END OF A LIST
C
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      INTEGER FLAG,N,UMAX
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(*),D(*),U(*)
      INTEGER IA(*),IJU(*),IL(*),IP(*),IU(*),JA(*),JL(*),JU(*),P(*)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION TAU
      INTEGER MC,SRLS
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION BOUND,DEL,DK,EK,ELT,ELT2,ELTNEW,EMAX,EPS,EPS1,
     +                 GAMMA,UKIDI,W,WW,XI,XIN,ZERO
      INTEGER I,ILI,IROW,J,JMAX,JMIN,JUMUJ,K,KK,KKMAX,KKMIN,MU,NEXTI,VJ
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MAX,MIN
C     ..
C     .. Common blocks ..
      COMMON /TAUINPUT/TAU,SRLS,MC
C     ..
      FLAG = 0
      ZERO = 0.0D0
      EMAX = ZERO
C ==================================================================end
C
C----CHECK FOR SUFFICIENT STORAGE FOR U
      IF (IU(N+1)-1.GT.UMAX) GO TO 140
C
C----INITIALIZATION
      DO 10 K = 1,N
          D(K) = 0
          JL(K) = 0
   10 CONTINUE
C ====================  change #3  (insertion) =======================3
C Calculate GAMMA and XI, the largest magnitudes of the diag. and  off-
C diag. elements, respectively. When the diag. elts. are stored first
C in A in each row (PATH = 4 or 5 in ODRV), GAMMA=max(GAMMA,A(IA(i))),
C i=1,...,IA(n+1)-1. We assume that this IS the case. (If this were
C later changed, then for each row I we would have to loop through KK
C = KKMIN,..., KKMAX  where KKMIN = IA(I), KKMAX = IA(I+1)-1, and test
C whether I = JA(KK), ie. row index = column index. If this equality
C holds, the element is a diagonal). Then calculate DEL and BOUND:
C DEL =  max ( max(XI,GAMMA)*EPS, EPS) where EPS is a small given
C number, and  BOUND = max ( XI/N, GAMMA, EPS).
C =====================================================================
      EPS = 1.0D-06
      GAMMA = ZERO
      XI = ZERO
      DO 30 IROW = 1,N
          GAMMA = MAX(GAMMA,ABS(A(IA(IROW))))
          KKMIN = IA(IROW) + 1
          KKMAX = IA(IROW+1) - 1
          IF (KKMIN.GT.KKMAX) GO TO 30
          DO 20 KK = KKMIN,KKMAX
              XI = MAX(XI,ABS(A(KK)))
   20     CONTINUE
   30 CONTINUE

      EPS1 = MAX(GAMMA,XI)*EPS
      DEL = MAX(EPS,EPS1)
      XIN = N
      XIN = XI/XIN
      BOUND = MAX(GAMMA,XIN,EPS)

C ==================================================================end
C
C----FOR EACH ROW K
      DO 130 K = 1,N
C
C------INITIALIZE K-TH ROW WITH ELEMENTS NONZERO IN ROW P(K) OF M
          JMIN = IA(P(K))
          JMAX = IA(P(K)+1) - 1
          IF (JMIN.GT.JMAX) GO TO 50
          DO 40 J = JMIN,JMAX
              VJ = IP(JA(J))
              IF (K.LE.VJ) D(VJ) = A(J)
   40     CONTINUE
C
C------MODIFY K-TH ROW BY ADDING IN THOSE ROWS I WITH U(I,K) NE 0
C------FOR EACH ROW I TO BE ADDED IN
   50     DK = D(K)
          I = JL(K)
   60     IF (I.EQ.0) GO TO 90
          NEXTI = JL(I)
C
C--------COMPUTE MULTIPLIER AND UPDATE DIAGONAL ELEMENT
          ILI = IL(I)
          UKIDI = -U(ILI)*D(I)
          DK = DK + UKIDI*U(ILI)
          U(ILI) = UKIDI
C
C--------ADD MULTIPLE OF ROW I TO K-TH ROW ...
          JMIN = ILI + 1
          JMAX = IU(I+1) - 1
          IF (JMIN.GT.JMAX) GO TO 80
          MU = IJU(I) - IU(I)
          DO 70 J = JMIN,JMAX
              D(JU(MU+J)) = D(JU(MU+J)) + UKIDI*U(J)
   70     CONTINUE
C
C--------... AND ADD I TO ROW LIST FOR NEXT NONZERO ENTRY
          IL(I) = JMIN
          J = JU(MU+JMIN)
          JL(I) = JL(J)
          JL(J) = I
C
   80     I = NEXTI
          GO TO 60
C
C ====================  change #4  (replacement) =====================4
C WAS:
C------CHECK FOR ZERO PIVOT
C  9      IF (DK.EQ.0)  GO TO 108
C =====================================================================
C STATEMENT 9 ABOVE WILL BE MODIFIED TO RESET Dk IN THE EVENT THE
C THE MATRIX IS NOT SUFF. POSITIVE-DEFINITE. NOTE THAT EVEN WHEN Dk>0,
C IT MAY BE MODIFIED IF THE MATRIX IS NOT POS. DEF!
C
C Dk is set as:  Dk = MAX ( ABS(Dk), DEL, (ELT**2)/BOUND), where
C ELT is the largest magnitude among the elements in the Kth row of U.
C This restriction guarantees that all elts. of D are strictly positive
C and that the elts. of the factors satisfy a uniform bound.
C [   Recall that we work with the auxiliary quantities  Vik = Uik * Dk.
C     The bound we want to impose on the elts. of U,
C          ( max(Uik)**2 )  * Dk  <=  BOUND, is equivalent to
C          ( max(Vik)**2 )  / Dk  <=  BOUND, or
C          Dk   >=    (max(Vik)**2) / BOUND.)
C     The value for ELT = max(Vik), max over i for fixed k, is found by
C     looping through the appropriate elements of U. These elements
C     are currently stored in part of D.  ]
C
C =====================================================================
C
   90     W = DK
C =====================================================================
          ELT = ZERO

          JMIN = IU(K)
          JMAX = IU(K+1) - 1
          MU = IJU(K) - JMIN
          IF (JMIN.GT.JMAX) GO TO 110
          DO 100 J = JMIN,JMAX
              ELTNEW = ABS(D(JU(MU+J)))
              ELT = MAX(ELT,ELTNEW)
  100     CONTINUE
  110     CONTINUE
C

          ELT2 = ELT*ELT

c--- xie: the modified Cholesky (MC) factorization
c         MC = 1: the original MC by Gill and Murray
c         MC = 2: the UMC by Xie & Schlick 6/12/96
c----------------------------------------------------
          IF (MC.EQ.1) THEN
              DK = MAX(ABS(DK),DEL,ELT2/BOUND)
              EK = DK - W

          ELSE
              WW = DK + TAU
              IF (WW.GT.DEL) THEN
                  DK = MAX(WW,ELT2/BOUND)

              ELSE IF (ABS(WW).LE.DEL) THEN
                  DK = DEL

              ELSE
                  DK = MIN(WW,-ELT2/BOUND)
              END IF

              EK = DK - W
          END IF
c-------------------------- End of UMC

          IF (EK.GT.EMAX) THEN
              EMAX = EK
              FLAG = -P(K)
          END IF

C------SAVE DIAGONAL ELEMENT
          D(K) = 1/DK

C------SAVE NONZERO ENTRIES IN K-TH ROW OF U ...
          JMIN = IU(K)
          JMAX = IU(K+1) - 1
          IF (JMIN.GT.JMAX) GO TO 130
          MU = IJU(K) - JMIN
          DO 120 J = JMIN,JMAX
              JUMUJ = JU(MU+J)
              U(J) = D(JUMUJ)
              D(JUMUJ) = 0
  120     CONTINUE

C------ AND ADD K TO ROW LIST FOR FIRST NONZERO ENTRY IN K-TH ROW
          IL(K) = JMIN
          I = JU(MU+JMIN)
          JL(K) = JL(I)
          JL(I) = K
  130 CONTINUE

C uncomment next line to check for the largest diagonal modification
C     IF (FLAG .LT. 0) WRITE (6,*) '      NMAX, EMAX', FLAG, EMAX
C
C ====================  change #5  (deletion) ========================5
C WAS:     FLAG = 0
C ==================================================================end
      RETURN
C
C ** ERROR -- INSUFFICIENT STORAGE FOR U
  140 FLAG = 7*N + 1
      RETURN

      END
!======================================================================
!
! BLAS SUBROUTINES

C***********************************************************************
C***********************************************************************
C     BLAS LEVEL 1, DOUBLE PRECISION
C     FROM NETLIB, FRI OCT 12 17:03:00 EDT 1990
C***********************************************************
      SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)
C --------------------------------------------------------
C  CONSTANT TIMES A VECTOR PLUS A VECTOR.
C  USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
C  JACK DONGARRA, LINPACK, 3/11/78.
C --------------------------------------------------------
      DOUBLE PRECISION DX(*),DY(*),DA
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N

      IF (N. LE. 0) RETURN
      IF (DA .EQ. 0.0D0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
C --------------------------------------------------------
C  CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C  NOT EQUAL TO 1
C --------------------------------------------------------
      IX = 1
      IY = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DY(IY) + DA*DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C --------------------------------------------------------
C  CODE FOR BOTH INCREMENTS EQUAL TO 1
C  CLEAN-UP LOOP
C --------------------------------------------------------
   20 M = MOD(N,4)
      IF (M .EQ. 0) GO TO 40
      DO 30 I = 1,M
        DY(I) = DY(I) + DA*DX(I)
   30 CONTINUE
      IF (N .LT. 4) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
        DY(I) = DY(I) + DA*DX(I)
        DY(I + 1) = DY(I + 1) + DA*DX(I + 1)
        DY(I + 2) = DY(I + 2) + DA*DX(I + 2)
        DY(I + 3) = DY(I + 3) + DA*DX(I + 3)
   50 CONTINUE
      RETURN
      END
C***********************************************************************
C***********************************************************************
      SUBROUTINE  DCOPY(N,DX,INCX,DY,INCY)
C --------------------------------------------------------
C  COPIES A VECTOR, X, TO A VECTOR, Y.
C  USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
C  JACK DONGARRA, LINPACK, 3/11/78.
C --------------------------------------------------------
      DOUBLE PRECISION DX(*),DY(*)
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N

      IF (N .LE. 0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
C --------------------------------------------------------
C  CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C  NOT EQUAL TO 1
C --------------------------------------------------------
      IX = 1
      IY = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C --------------------------------------------------------
C  CODE FOR BOTH INCREMENTS EQUAL TO 1
C  CLEAN-UP LOOP
C --------------------------------------------------------
   20 M = MOD(N,7)
      IF (M .EQ. 0) GO TO 40
      DO 30 I = 1,M
        DY(I) = DX(I)
   30 CONTINUE
      IF (N .LT. 7) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,7
        DY(I) = DX(I)
        DY(I + 1) = DX(I + 1)
        DY(I + 2) = DX(I + 2)
        DY(I + 3) = DX(I + 3)
        DY(I + 4) = DX(I + 4)
        DY(I + 5) = DX(I + 5)
        DY(I + 6) = DX(I + 6)
   50 CONTINUE
      RETURN
      END
C***********************************************************************
C***********************************************************************
      DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY)
C --------------------------------------------------------
C  FORMS THE DOT PRODUCT OF TWO VECTORS.
C  USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
C  JACK DONGARRA, LINPACK, 3/11/78.
C --------------------------------------------------------
      DOUBLE PRECISION DX(*),DY(*),DTEMP
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N

      DDOT = 0.0D0
      DTEMP = 0.0D0
      IF (N .LE. 0) RETURN
      IF (INCX.EQ.1. AND. INCY.EQ.1) GO TO 20
C --------------------------------------------------------
C  CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C  NOT EQUAL TO 1
C --------------------------------------------------------
      IX = 1
      IY = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DTEMP = DTEMP + DX(IX)*DY(IY)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      DDOT = DTEMP
      RETURN
C --------------------------------------------------------
C   CODE FOR BOTH INCREMENTS EQUAL TO 1
C   CLEAN-UP LOOP
C --------------------------------------------------------
   20 M = MOD(N,5)
      IF (M .EQ. 0) GO TO 40
      DO 30 I = 1,M
        DTEMP = DTEMP + DX(I)*DY(I)
   30 CONTINUE
      IF (N .LT. 5) GO TO 60
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        DTEMP = DTEMP + DX(I)*DY(I) + DX(I + 1)*DY(I + 1) +
     *   DX(I + 2)*DY(I + 2) + DX(I + 3)*DY(I + 3) + DX(I + 4)*DY(I + 4)
   50 CONTINUE
   60 DDOT = DTEMP
      RETURN
      END
C***********************************************************************
C***********************************************************************
      DOUBLE PRECISION FUNCTION DNRM2 (N,DX,INCX)
      INTEGER NEXT
      DOUBLE PRECISION DX(*),CUTLO,CUTHI,HITEST,SUM,XMAX,ZERO,ONE
      DATA ZERO,ONE /0.0D0,1.0D0/
C --------------------------------------------------------
C  EUCLIDEAN NORM OF THE N-VECTOR STORED IN DX() WITH STORAGE
C  INCREMENT INCX .
C  IF    N .LE. 0 RETURN WITH RESULT = 0.
C  IF N .GE. 1 THEN INCX MUST BE .GE. 1
C        C.L.LAWSON, 1978 JAN 08
C  FOUR PHASE METHOD     USING TWO BUILT-IN CONSTANTS THAT ARE
C  HOPEFULLY APPLICABLE TO ALL MACHINES.
C      CUTLO = MAXIMUM OF  DSQRT(U/EPS)  OVER ALL KNOWN MACHINES.
C      CUTHI = MINIMUM OF  DSQRT(V)      OVER ALL KNOWN MACHINES.
C  WHERE
C      EPS = SMALLEST NO. SUCH THAT EPS + 1. .GT. 1.
C      U   = SMALLEST POSITIVE NO.   (UNDERFLOW LIMIT)
C      V   = LARGEST  NO.            (OVERFLOW  LIMIT)
C
C  BRIEF OUTLINE OF ALGORITHM..
C  PHASE 1    SCANS ZERO COMPONENTS.
C  MOVE TO PHASE 2 WHEN A COMPONENT IS NONZERO AND .LE. CUTLO
C  MOVE TO PHASE 3 WHEN A COMPONENT IS .GT. CUTLO
C  MOVE TO PHASE 4 WHEN A COMPONENT IS .GE. CUTHI/M
C  WHERE M = N FOR X() REAL AND M = 2*N FOR COMPLEX.
C  VALUES FOR CUTLO AND CUTHI..
C  FROM THE ENVIRONMENTAL PARAMETERS LISTED IN THE IMSL CONVERTER
C  DOCUMENT THE LIMITING VALUES ARE AS FOLLOWS..
C  CUTLO, S.P.   U/EPS = 2**(-102) FOR  HONEYWELL.  CLOSE SECONDS ARE
C                UNIVAC AND DEC AT 2**(-103)
C                THUS CUTLO = 2**(-51) = 4.44089E-16
C  CUTHI, S.P.   V = 2**127 FOR UNIVAC, HONEYWELL, AND DEC.
C                THUS CUTHI = 2**(63.5) = 1.30438E19
C  CUTLO, D.P.   U/EPS = 2**(-67) FOR HONEYWELL AND DEC.
C                THUS CUTLO = 2**(-33.5) = 8.23181D-11
C  CUTHI, D.P.   SAME AS S.P.  CUTHI = 1.30438D19
C  DATA CUTLO, CUTHI / 8.232D-11, 1.304D19 /
C  DATA CUTLO, CUTHI / 4.441E-16, 1.304E19 /

      INTEGER N, NN, INCX, i, J
      DATA CUTLO, CUTHI / 8.232D-11, 1.304D19 /

      IF (N .GT. 0) GO TO 10
         DNRM2  = ZERO
         GO TO 300

   10 ASSIGN 30 TO NEXT
      SUM = ZERO
      NN = N * INCX
C ------------------
C BEGIN MAIN LOOP
C ------------------
      I = 1
   20    GO TO NEXT,(30, 50, 70, 110)
   30 IF (DABS(DX(I)) .GT. CUTLO) GO TO 85
      ASSIGN 50 TO NEXT
      XMAX = ZERO
C ------------------
C PHASE 1.  SUM IS ZERO
C ------------------
   50 IF (DX(I) .EQ. ZERO) GO TO 200
      IF (DABS(DX(I)) .GT. CUTLO) GO TO 85
C ------------------
C PREPARE FOR PHASE 2.
C ------------------
      ASSIGN 70 TO NEXT
      GO TO 105
C ------------------
C PREPARE FOR PHASE 4.
C ------------------
  100 I = J
      ASSIGN 110 TO NEXT
      SUM = (SUM / DX(I)) / DX(I)
  105 XMAX = DABS(DX(I))
      GO TO 115
C ------------------
C SUM IS SMALL.
C SCALE TO AVOID DESTRUCTIVE UNDERFLOW.
C ------------------
   70 IF (DABS(DX(I)) .GT. CUTLO) GO TO 75
C ------------------
C COMMON CODE FOR PHASES 2 AND 4.
C IN PHASE 4 SUM IS LARGE.  SCALE TO AVOID OVERFLOW.
C ------------------
  110 IF (DABS(DX(I)) .LE. XMAX) GO TO 115
         SUM = ONE + SUM * (XMAX / DX(I))**2
         XMAX = DABS(DX(I))
         GO TO 200
  115 SUM = SUM + (DX(I)/XMAX)**2
      GO TO 200
C ------------------
C PREPARE FOR PHASE 3.
C ------------------
   75 SUM = (SUM * XMAX) * XMAX
C ------------------
C FOR REAL OR D.P. SET HITEST = CUTHI/N
C FOR COMPLEX      SET HITEST = CUTHI/(2*N)
C ------------------
   85 HITEST = CUTHI/FLOAT(N)
C ------------------
C PHASE 3.  SUM IS MID-RANGE.  NO SCALING.
C ------------------
      DO 95 J =I,NN,INCX
      IF (DABS(DX(J)) .GE. HITEST) GO TO 100
   95    SUM = SUM + DX(J)**2
      DNRM2 = DSQRT(SUM)
      GO TO 300

  200 CONTINUE
      I = I + INCX
      IF (I .LE. NN) GO TO 20
C ------------------
C END OF MAIN LOOP.
C COMPUTE SQUARE ROOT AND ADJUST FOR SCALING.
C ------------------
      DNRM2 = XMAX * DSQRT(SUM)
  300 CONTINUE
      RETURN
      END
