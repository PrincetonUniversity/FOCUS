!title ( nlinear ) ! Nonlinear solver using NAG routine.
!latex \briefly{Using NAG nonlinear solver to optimize. Either use Powell hybrid method to find a solution for a system of nonlinear equations, or take 
!latex         modified-Newton algorithm for finding a minimum of cost functions.}
!latex 
!latex \calledby{\link{notopt}}
!latex 
!latex \tableofcontents
!latex 
!latex \subsection{Powell hybrid method}
!latex We are calling \nag{http://www.nag.co.uk/numeric/FL/manual/pdf/C05/c05pdf.pdf}{C05PDF} to find to find a solution for a system of nonlinear equations 
!latex which is composed with the first derivatives of cost functions. And also the first derivatives for the equations, which is the second derivatives of 
!latex cost funtions are provided. The maximum calling count is also decided by th \inputvar{Ntauout}.
!latex 
!latex \subsection{Modified-Newton method}
!latex not yet implemented. \nag{http://www.nag.co.uk/numeric/FL/manual/pdf/E04/e04lbf.pdf}{E04LBF}
!latex
!latex \subsection{Sigular-value decomposition}
!latex \nag{http://www.nag.co.uk/numeric/FL/manual/pdf/F08/f08kbf.pdf}{F08KBF}
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!---------------------------------------------------------------------------------------------------------
!  C05PDF:  Nonlinear equations solver with first derivatives           May/30/2016
!---------------------------------------------------------------------------------------------------------
SUBROUTINE Powell
  use kmodule, only: zero, sqrtmachprec, myid, ounit, ext, evolution, Ntauout, itau, &
                     totalenergy, n1E, n2E, NFcoil, Ncoils, Cdof, Ndof,   &
                     bnorm, tflux, ttlen, eqarc, xtol, Idisplay
  use hdf5
  implicit none
  include "mpif.h"


  INTEGER           :: N, IREVCM, MODE, LDFJAC, LR, IFAIL, idof, astat, ierr, mm, NN, icoil
  REAL, allocatable :: FVEC(:), DIAG(:), FJAC(:,:), R(:), QTF(:), W(:,:)
  REAL              :: FACTOR, X(Ndof)
  INTEGER           :: ii,jj,c1, n1, c2, n2, irestart
  
! the following are used by the macros HWRITEXX below; do not alter/remove;
  INTEGER            :: hdfier, rank
  integer(hid_t)     :: file_id, space_id, dset_id                 ! warning: the string "INTEGER" is a macro;
  integer(hsize_t)   :: onedims(1:1), twodims(1:2), threedims(1:3) ! warning: the string "INTEGER" is a macro;
  INTEGER            :: imn
  REAL               :: tvolume

  NN = NFcoil; N = Ndof

  LDFJAC = N
  LR     = N*(N+1)/2
  IFAIL  = Idisplay

  SALLOCATE( FVEC,(N       ), zero)
  SALLOCATE( DIAG,(N       ), zero)
  SALLOCATE( FJAC,(LDFJAC,N), zero)
  SALLOCATE( R   ,(LR      ), zero)
  SALLOCATE( QTF ,(N       ), zero)
  SALLOCATE( W   ,(N     ,4), zero)

  call ntpack(X(1:N))
  call costfun(2); call ntchdms(2)
  if (myid .eq. 0) write(ounit, '("nlinear :",12X,"; Initial Value " 11X,"; E = ",es23.15,"; D = ",es23.15,"; B = ",es23.15,"; F = ",es23.15,"; L = ",es23.15,"; A = ",es23.15)') &
                         totalenergy, sqrt(sum(n1E(1:Ncoils,0:Cdof)**2)/N), bnorm, tflux, ttlen, eqarc
  do jj = 1, N
     call ntconvert(jj,c2,n2)
     do ii = 1,N
        call ntconvert(ii,c1,n1)
        FJAC(ii,jj) = n2E(c1,n1,c2,n2)
     enddo
  enddo

!!$  !----------------------------------------------
!!$  if(myid .ne. 0) return
!!$  call h5open_f( hdfier ) ! initialize Fortran interface;
!!$  FATAL( SVD    , hdfier.ne.0, error calling h5open_f )
!!$
!!$  call h5fcreate_f( trim(ext)//".matrix.h5", H5F_ACC_TRUNC_F, file_id, hdfier ) ! create new file;
!!$  FATAL( SVD    , hdfier.ne.0, error calling h5fcreate_f )
!!$
!!$  HWRITERA( n, n       ,  Hessian  ,  FJAC(1:n,1:n))
!!$  
!!$  call h5fclose_f( file_id, hdfier ) ! terminate access;
!!$  FATAL( nlinear, hdfier.ne.0, error calling h5fclose_f )
!!$  
!!$  call h5close_f( hdfier ) ! close Fortran interface;
!!$  FATAL( nlinear, hdfier.ne.0, error calling h5close_f )
!!$
!!$  return
!!$  !----------------------------------------------
  
  itau = 0
  !XTOL   = sqrtmachprec
  MODE   = 1
  FACTOR = 100.0
 ! DIAG(1:N) = 1.0
  IREVCM = 0

  call  C05PDF(IREVCM,N,X,FVEC,FJAC,LDFJAC,XTOL,DIAG,MODE,FACTOR,R,LR,QTF,W,IFAIL)

  do while (IREVCM .ne. 0)

     if(itau .gt. Ntauout) exit

     if (IREVCM .eq. 1) then
        itau = itau + 1
        if (myid .eq. 0) write(ounit, '("nlinear :",12X,"; icount = ",i6,10X,"; E = ",es23.15,"; D = ",es23.15,"; B = ",es23.15,"; F = ",es23.15,"; L = ",es23.15,"; A = ",es23.15)') &
                         itau, totalenergy, sqrt(sum(FVEC(1:N)**2)/N), bnorm, tflux, ttlen, eqarc
        evolution(itau,0) = real(itau) / Ntauout
        evolution(itau,1) = totalenergy
        call restart(irestart)
     else if (IREVCM .eq. 2) then
        call ntunpack(X(1:N))
        call costfun(1); call ntchdms(1)
        do idof = 1,N
           call ntconvert(idof,c1,n1)
           FVEC(idof) = n1E(c1,n1)
        enddo
     else if (IREVCM .eq. 3) then
        call ntunpack(X(1:N))
        call costfun(2); call ntchdms(2)
        do jj = 1, N
           call ntconvert(jj,c2,n2)
           do ii = 1,N
              call ntconvert(ii,c1,n1)
              FJAC(ii,jj) = n2E(c1,n1,c2,n2)
           enddo
        enddo
     endif

     call  C05PDF(IREVCM,N,X,FVEC,FJAC,LDFJAC,XTOL,DIAG,MODE,FACTOR,R,LR,QTF,W,IFAIL)

  enddo


  if (ifail .eq. 0) then
     if (myid .eq. 0) write(ounit,*) "Find a good solution."
     if (myid .eq. 0) write(ounit, '("nlinear :",12X,"; icount = ",i6,10X,"; E = ",es23.15,"; D = ",es23.15)') &
                      itau, totalenergy, sqrt(sum(FVEC(1:N)**2)/N)
     call restart(irestart)
  else
     if (myid .eq. 0) write(ounit,*) "ifail not equal zero."
     if (myid .eq. 0) write(ounit,'("ifail = ",i2)') ifail
     if (myid .eq. 0) write(ounit, '("nlinear :",12X,"; icount = ",i6,10X,"; E = ",es23.15,"; D = ",es23.15)') &
                      itau, totalenergy, sqrt(sum(FVEC(1:N)**2)/N)
     call restart(irestart)
  endif

  DEALLOCATE( FVEC )
  DEALLOCATE( DIAG )
  DEALLOCATE( FJAC )
  DEALLOCATE( R    )
  DEALLOCATE( QTF  )
  DEALLOCATE( W    )

  return
  
END SUBROUTINE Powell

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!---------------------------------------------------------------------------------------------------------
!  NAG E04LBF: Minimizing subroutine using modified Newton algorithm with first and second derivatives
!---------------------------------------------------------------------------------------------------------
SUBROUTINE Newton

  use kmodule, only: zero, one, myid, ounit, Ndof, Ntauout, stepmx, xtol, eta, Idisplay, &
                     totalenergy, bnorm, tflux, ttlen, eqarc
  implicit none
  include "mpif.h"

  INTEGER              :: liw, lw, n, lh, ibound, ifail, iprint, MAXCAL, ierr, astat,  irestart
  INTEGER, allocatable :: istate(:), iw(:)
  REAL                 :: f
  REAL   , allocatable :: bl(:), bu(:), g(:), hesd(:), hesl(:), w(:), x(:)

  Intrinsic            :: count
  EXTERNAL             :: funct, H, monit, pack, unpak, costfun, e04lbf


  n = Ndof; iprint = -1 * Idisplay; MAXCAL = max(Ntauout, 50*n); ibound = 1 !no constraints
  lh = max(n*(n-1)/2,1); LIW = 2; LW = max(7*n+n*(n-1)/2,8); ifail = Idisplay; Ntauout = MAXCAL

  FATAL(Newton, n .lt. 1, n is wrong)
  FATAL(Newton, ifail .lt. -1 .or.  ifail .gt. 1, ifail unsupported)
  FATAL(Newton, ibound .lt. 0 .or. ibound .gt. 4, ibound value not supported)


  SALLOCATE(x ,     (n), zero)
  SALLOCATE(g ,     (n), zero)
  SALLOCATE(bl,     (n), zero)
  SALLOCATE(bu,     (n), zero)
  SALLOCATE(w ,    (lw), zero)
  SALLOCATE(iw,   (liw), zero)
  SALLOCATE(hesd  , (n), zero)
  SALLOCATE(hesl  ,(lh), zero)
  SALLOCATE(istate, (n), zero)

  call pack(X(1:n))
  call unpack(x(1:n))
  
  call costfun(0); if (myid .eq.0) write (ounit,9998) totalenergy, bnorm, tflux, ttlen, eqarc

  call e04lbf(n,funct,h,monit,iprint,maxcal,eta,xtol,stepmx,ibound,bl,bu,x,hesl,lh,hesd,istate,f,g,iw,liw,w,lw,ifail)

  call unpack(x(1:n))
  call restart(irestart)
  call costfun(0); if (myid .eq.0) write (ounit,9997) totalenergy, bnorm, tflux, ttlen, eqarc

  deallocate(x, g, bl, bu, w, iw, hesd, hesl, istate)

  if ( myid .ne. 0) return   ! diagnos
  select case (ifail) 
  case (:(-1))
     write (ounit,9999) ifail, 'Abnormal exit; Check IFLAG or other things.'
  case (0)
     write (ounit,9999) ifail, 'Success     0; E04LBF running successfully.'
  case (1)
     write (ounit,9999) ifail, 'Error       1; wrong input value.'
  case (2)
     write (ounit,9999) ifail, 'Warning     2; Maxcal reached.'
  case (3)
     write (ounit,9999) ifail, 'Warning     3; Not a minimum, just lower.'
  case (4)

  case (5)
     write (ounit,9999) ifail, "Warning     5; Not a minimum, can't keep going"   
  end select
 
9999 format ("E04LBF  : " 10X " : IFAIL = " I6, 5X, A )
9998 format ("E04LBF  : " 10X " : COST FUNCTIONS before entering E04LBF : E = "ES23.15" ; B = "ES23.15" ; F = "ES23.15" ; L = "ES23.15" ; A = "ES23.15)
9997 format ("E04LBF  : " 10X " : COST FUNCTIONS after  entering E04LBF : E = "ES23.15" ; B = "ES23.15" ; F = "ES23.15" ; L = "ES23.15" ; A = "ES23.15)

    return

end SUBROUTINE Newton

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

Subroutine funct(iflag,n,xc,fc,gc,iw,liw,w,lw)
  use kmodule, only: totalenergy, t1E
  Implicit None
  ! Routine to evaluate objective function and its 1st derivatives.

  Integer, Intent (Inout) :: iflag, iw(liw)
  Integer, Intent (In   ) :: liw, lw, n
  Real   , Intent (Inout) :: w(lw)
  Real   , Intent (In   ) :: xc(n)
  Real   , Intent (Out  ) :: fc, gc(n)

  INTEGER                 :: idof, c1, n1
  external                :: costfun
  ! .. Executable Statements ..
  call unpack(xc(1:n))
  call costfun(1)

  fc = totalenergy
  do idof = 1, n
     call DoFconvert(idof,c1,n1)
     gc(idof) = t1E(c1,n1)
  enddo

  Return
End Subroutine funct

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

Subroutine H(iflag,n,xc,fhesl,lh,fhesd,iw,liw,w,lw)
  use kmodule, only: t2E
  Implicit None

  ! Routine to evaluate objective function and its 1st derivatives.
  Integer, Intent (Inout) :: iflag, iw(liw)
  Integer, Intent (In   ) :: liw, lw, n, lh
  Real   , Intent (Out  ) :: fhesl(lh)
  Real   , Intent (Inout) :: w(lw), fhesd(n)
  Real   , Intent (In   ) :: xc(n)

  INTEGER                 :: idof, jdof, c1, n1, c2, n2,nn
  external                :: costfun
  ! .. Executable Statements ..

  call unpack(xc(1:n))
  call costfun(2)

  do idof = 2, n

     call DoFconvert(idof,c1,n1)
     do jdof = 1, idof-1
        call DoFconvert(jdof,c2,n2)
        nn = (idof-2)*(idof-1)/2 + jdof ! lower terms
        fhesl( nn ) = t2E(c1,n1,c2,n2)
     enddo
  enddo

  do idof = 1, n
     call DoFconvert(idof,c1,n1)
     fhesd(idof) = t2E(c1,n1,c1,n1)  !diagnal terms
  enddo

  Return
End Subroutine H

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

Subroutine monit(n,xc,fc,gc,istate,gpjnrm,cond,posdef,niter,nf,iw,liw,w,lw)

  use kmodule, only: myid, ounit, Idisplay, Ntauout, &
                     totalenergy, bnorm, tflux, ttlen, eqarc, evolution
  implicit none
  include "mpif.h"

  ! Monitoring routine

  INTEGER, INTENT ( inout ) :: iw(liw)
  INTEGER, INTENT ( in    ) :: liw, lw, n, nf, niter, istate(n)
  REAL   , INTENT ( inout ) :: w(lw)
  REAL   , INTENT ( in    ) :: cond, fc, gpjnrm, gc(n), xc(n)
  LOGICAL, INTENT ( in    ) :: posdef

  ! .. Local Scalars ..
  Integer :: isj, j, irestart
  ! .. Executable Statements ..

  evolution(niter, 0) = real(niter) / Ntauout
  evolution(niter, 1) = totalenergy
  call restart(irestart)

  if ( myid .ne. 0) return
  write (ounit,9999) niter, Ntauout, nf, fc, gpjnrm, bnorm, tflux, ttlen, eqarc

!!$  if ( Idisplay .eq. -1 ) then
!!$     write ( ounit, '("E04LBF  : "10X" : "3(A15," ; "))') 'J', 'X(J)', 'G(J)', 'STATUS'
!!$     do j = 1, n
!!$        isj = istate(j)
!!$        select case (isj)
!!$        case (1:)
!!$           write (ounit,9998) j, xc(j), gc(j), ' Free'
!!$        case (-1)
!!$           write (ounit,9998) j, xc(j), gc(j), ' Upper Bound'
!!$        case (-2)
!!$           write (ounit,9998) j, xc(j), gc(j), ' Lower Bound'
!!$        case (-3)
!!$           write (ounit,9998) j, xc(j), gc(j), ' Constant '
!!$        end select
!!$     enddo
!!$  endif

  If (cond .ne. 0.0) Then
     If (cond>1.0E6) Then
        Write (ounit,'("E04LBF  : " 10X " : "A)') 'Estimated condition number of projected Hessian is more than 1.0E+6'
     Else
        Write (ounit,'("E04LBF  : " 10X " : "A,ES23.15)') 'Estimated condition number of projected Hessian = ', cond
     End If
     If (.Not. posdef) Then
        Write (ounit,'("E04LBF  : " 10X " : "A)') 'Projected Hessian matrix is not positive definite'
     End If
  End If

9999 format ("E04LBF  : " I10 " / " I10  " : nfunt = "I6" ; E = "ES23.15" ; G = "ES23.15" ; B = "ES23.15" ; F = "ES23.15" ; L = "ES23.15" ; A = "ES23.15)
9998 format ("E04LBF  : " 10X " : "I15" ; "ES15.7" ; "ES15.7" ; "A15)


  return

End Subroutine monit

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

SUBROUTINE ntchdms(nderiv)   ! NewTon CHange DiMenSion
! change the number of DoF adding in lambda terms
  use kmodule, only: zero, sqrtmachprec, myid, ounit, Idisplay, iter, &
       Ncoils, NFcoil, Ndof, Cdof, &
       totalenergy, t1E, t2E, n1E, n2E, &
       eqarc      , t1A, t2A, weight_eqarc, dlc, dls


  implicit none
  include "mpif.h"

  INTEGER      :: nderiv

  INTEGER      :: astat, ierr
  INTEGER      :: lCdof, icoil, jcoil, ll, mm

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  ! begin change dimension

  FATAL(ntchdsm, nderiv .ne. 1 .and. nderiv .ne. 2, unsupported nderiv)
  FATAL(ntchdsm, .not. allocated(t1E) , should allocate t1E first)

  Ndof = (1+3+8*NFcoil)*Ncoils
  
  if(allocated(n1E)) deallocate(n1E)
  if(allocated(n2E)) deallocate(n2E)
  lCdof = 8*NFcoil+8
      
  SALLOCATE(n1E, (1:Ncoils, 0:lCdof), zero)
  do icoil = 1, Ncoils
     n1E(icoil, 0:Cdof) = t1E(icoil, 0:Cdof)
     n1E(icoil, Cdof+1:Cdof+NFcoil+1) = dlc(icoil,0:NFcoil,0)
     n1E(icoil, Cdof+NFcoil+2:Cdof+2*NFcoil+2) = dls(icoil,0:NFcoil,0)
  enddo
  
  if(nderiv .eq. 2) then

     SALLOCATE(n2E, (1:Ncoils, 0:lCdof, 1:Ncoils, 0:lCdof), zero)
     
     do icoil = 1, Ncoils
        
        do ll = 0, Cdof
           do jcoil = 1, Ncoils
              n2E(icoil, ll, jcoil, 0:Cdof) = t2E(icoil, ll, jcoil, 0:Cdof)
           enddo
        enddo
        
        do ll = 0, NFcoil

           n2E(icoil, ll+Cdof       +1, icoil, 0:Cdof) = dlc(icoil, ll, 0:Cdof)
           n2E(icoil, ll+Cdof+NFcoil+2, icoil, 0:Cdof) = dls(icoil, ll, 0:Cdof)

           do mm = 0, Cdof
              n2E(icoil, mm, icoil, ll+Cdof       +1) = dlc(icoil, ll, mm)
               n2E(icoil, mm, icoil, ll+Cdof+NFcoil+2) = dls(icoil, ll, mm)
           enddo

        enddo

        
     enddo    
  endif

  return
  
END SUBROUTINE ntchdms

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine SVD

  use kmodule, only: zero, machprec, myid, ncpu, ounit, lunit, ext,  Ndof, t2E
  use hdf5
  
  implicit none
  include "mpif.h"

  INTEGER           :: astat, ierr, n, info, lda, ldu, ldvt, lwork, isingular, &
                       i, j, c1, n1, c2, n2, ifail, neigen, ieigen, nstep
  REAL              :: dummy(1,1), step, start, finish
  REAL, allocatable :: a(:,:), ab(:,:), inver(:,:), b(:), s(:), u(:,:), vt(:,:), work(:), w(:)

! the following are used by the macros HWRITEXX below; do not alter/remove;
  INTEGER            :: hdfier, rank
  integer(hid_t)     :: file_id, space_id, dset_id                 ! warning: the string "INTEGER" is a macro;
  integer(hsize_t)   :: onedims(1:1), twodims(1:2), threedims(1:3) ! warning: the string "INTEGER" is a macro;
  INTEGER            :: imn
  REAL               :: tvolume

  if (myid .eq. 0) write(ounit,'("SVD     : "10X" : "A)') "Start to analyze the Hessian matrix using SVD (F08KBF)."

  if (myid .eq. 0) write(ounit,'("SVD     : "10X" : Dimension of Hessian matrix is "I6" X " I6)') Ndof, Ndof

  n = Ndof; lda = n; ldu = n; ldvt = n
  
  !Allocate arrays
  SALLOCATE(a, (1:lda, 1:n), zero)
  SALLOCATE(ab,(1:lda, 1:n), zero)
  SALLOCATE(b, (1:n       ), zero)
  SALLOCATE(s, (1:n       ), zero)
  SALLOCATE(w, (1:n       ), zero)
  SALLOCATE(u, (1:ldu, 1:n), zero)
  SALLOCATE(vt,(1:ldvt,1:n), zero)
  SALLOCATE(inver,(1:lda, 1:n), zero)

  !get the Hessian matrix
  start = MPI_WTIME()
  call costfun(2)
  finish = MPI_WTIME()
  if (myid .eq. 0) WRITE(ounit, '("SVD     : "10X" : costfun(2) finished in "ES12.5" seconds.")') finish - start

  do j=1,n
     call DoFconvert(j, c2, n2)
     do i=1, n
        call DoFconvert(i, c1, n1)
        a(i,j) = t2E(c1, n1, c2, n2)
     enddo
  enddo

  !if( ncpu .gt. 1) stop "SVD only works on single node."
  if(myid .ne. 0) return

  ab = a
  !write(ounit,'("SVD     : "10X" : "4ES23.15)') a(1,2), a(2,1), t2E(1,0,1,1), t2E(1,1,1,0)

  ! get optimal work dimension
  start = MPI_WTIME()
  lwork = -1
  call DGESVD('A','A', n, n, a, lda, s, u, ldu, vt, ldvt, dummy, lwork, info)
  lwork = max(5*n, nint(dummy(1,1)))
  SALLOCATE(work, (lwork), zero)

  call DGESVD('A','A', n, n, a, lda, s, u, ldu, vt, ldvt, work , lwork, info)

  finish = MPI_WTIME()  
  WRITE(ounit, '("SVD     : "10X" : SVD finished in "ES12.5" seconds.")') finish - start

  write(ounit, '("SVD     : "10X" : INFO = "I4" ; min(abs(s))="es23.15" ;")') info, minval(abs(s)) 

  isingular = 0
  do i = 1, n
     if (s(i) .ge. machprec) isingular = isingular + 1
  enddo

  write(ounit, '("SVD     : "10X" : Rank = " I6" ; Max = "ES23.15" ; Min = "ES23.15)') isingular, s(1), s(n)

!!$  a = ab
!!$  ! calculate eigenvalues and eigenvectors
!!$  lwork = -1; deallocate(work)
!!$  call F08FAF('V','U', n, a, lda, w,dummy(1,1),lwork,info)
!!$  lwork = max(3*n-1,nint(dummy(1,1)))
!!$  Allocate (work(lwork))
!!$  call F08FAF('V','U', n, a, lda, w,work,lwork,info)
!!$
!!$  write(ounit, '("Eigen   : "10X" : INFO = "I4" ; min(w)="es23.15" ; max(w)="es23.15" ;")') info, minval(w), maxval(w)
!!$  ieigen = 0 ; neigen = 0
!!$  do i = 1, n
!!$     if (w(i) .ge. machprec) then 
!!$        neigen = neigen + 1
!!$     else
!!$        ieigen = i
!!$        write(ounit, '("Eigen   : "10X" : Find a negtivie eigenvalues at " I6)') ieigen
!!$     endif        
!!$  enddo
!!$  write(ounit, '("Eigen   : "10X" : number of positive eigenvalues = "I6)') neigen
!!$
!!$  !test
!!$  !b(1:n) = matmul(ab(1:n,1:n), a(1:n,ieigen)) - w(ieigen) * a(1:n,ieigen)
!!$  !write(ounit, '("Eigen   : "10X" : The square summation of the vector is"ES23.15)') sum(b(1:n)**2)
!!$
!!$  write(ounit, '("Eigen   : "10X" : The 24-th eigenvalue is"ES23.15)') w(24)
!!$  write(ounit, '("Eigen   : "10X" : The  N-th eigenvalue is"ES23.15)') w(n)
!!$  step = 1.0E-4; nstep = 100  !evolution stepsize
!!$  call coilevl( a(1:n,n), step, nstep )  !array for the direction; step for the stepsize  
!!$  write(ounit, '("Eigen   : "10X" : Writing coils evolution finished")')
  
  ! call inverse procedure
  !call matrinv(ab(1:n,1:n), inver(1:n,1:n), n, ifail)
  !FATAL(SVD     , ifail .ne. 0, inversing error)

  call h5open_f( hdfier ) ! initialize Fortran interface;
  FATAL( SVD    , hdfier.ne.0, error calling h5open_f )

  call h5fcreate_f( trim(ext)//".matrix.h5", H5F_ACC_TRUNC_F, file_id, hdfier ) ! create new file;
  FATAL( SVD    , hdfier.ne.0, error calling h5fcreate_f )

  HWRITEIV( 1          ,  info     ,  info       )
  HWRITEIV( 1          ,  rank     ,  isingular  )
  HWRITERA( n, n       ,  Hessian  ,  ab(1:n,1:n))
  HWRITERA( n, n       ,  inverse  ,  inver(1:n,1:n))
  HWRITERA( n, n       ,  output   ,  a(1:n,1:n) )
  HWRITERA( n, n       ,  U        ,  u(1:n,1:n) )
  HWRITERA( n, n       ,  VT       ,  vt(1:n,1:n))
  HWRITERV( n          ,  S        ,  s(1:n     ))
  if(info .gt. 0) then
     HWRITERV( lwork      ,  Work     ,work(1:lwork))
  endif
  
  call h5fclose_f( file_id, hdfier ) ! terminate access;
  FATAL( nlinear, hdfier.ne.0, error calling h5fclose_f )
  
  call h5close_f( hdfier ) ! close Fortran interface;
  FATAL( nlinear, hdfier.ne.0, error calling h5close_f )

  DALLOCATE(a)
  DALLOCATE(ab)
  DALLOCATE(b)
  DALLOCATE(u)
  DALLOCATE(s)
  DALLOCATE(vt)
  DALLOCATE(work)


  write(ounit, '("SVD     : "10X" : SVD finished and please check ", A, ".matrix.h5")') trim(ext)

  return

end subroutine SVD

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine matrinv( a, ai, n, ifail)

  !-------------------------------------------------------
  ! Only suit for real symmetric positive definite matrix;
  ! a(1:n,1:n) is the input matrix keeping unchanged;
  ! ai(1:n, 1:n) is the inversed matrix for a;
  !--------------------------------------------------------

  implicit none

  real     :: a(1:n,1:n), ai(1:n,1:n)
  integer  :: n, ifail

  INTEGER  :: i, j, lda, ldb
  REAL     :: ab(1:n+1,1:n), z(1:n)

  lda = 1+n
  ldb = n

  ab(1:n,1:n) = a(1:n,1:n)
  
  Call F01ABF(ab,lda,n,ai,ldb,z,ifail)

  if( ifail .ne. 0) return
  
  do j = 1,n
     do i = 1,j
        ai(i, j) = ai(j,i)
     enddo
  enddo

  return
end subroutine matrinv

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine coilevl( dir, stepsize, nstep )

  use kmodule, only : Ndof, coil, totalenergy, evolution, bnorm, tflux, ttlen, eqarc, ccsep, itau, t1E
  implicit none

  INTEGER :: nstep
  REAL    :: dir(1:Ndof), stepsize

  INTEGER :: istep, irestart
  REAL    :: xdof(1:Ndof)

  call wrtdir(dir)
  call pack(xdof(1:Ndof))
  irestart = 1; call restart(irestart)
  do istep = 0, nstep
     xdof = xdof + stepsize*istep*dir
     call unpack(xdof(1:Ndof))
     call costfun(1)
     if(allocated(evolution)) then
        itau = istep
        evolution(itau,0) = itau
        evolution(itau,1) = totalenergy
        evolution(itau,2) = (sum(t1E**2))
        evolution(itau,3) = bnorm
        evolution(itau,4) = tflux
        evolution(itau,5) = ttlen
        evolution(itau,6) = eqarc
        evolution(itau,7) = ccsep
        evolution(itau,8) = coil(1)%I
        evolution(itau,9) = coil(1)%L
     endif
     call restart(irestart)
  enddo

end subroutine coilevl

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine wrtdir(dir)
  use kmodule, only: zero,Ndof, Ncoils, NFcoil, Cdof, shudson, myid, ext, lunit
  implicit none
  include "mpif.h"

  REAL    :: dir(1:Ndof)

  INTEGER :: icoil, idof, mm, i, inum
  REAL    :: coil(1:Ncoils, 0:Cdof) !local variable; different with coil in global

  if(myid .ne. 0) return
  
  idof = 0; coil = zero
  
  do icoil = 1, Ncoils
     do idof = 1, 1+3+6*NFcoil
        inum = shudson(idof)
        coil(icoil,inum) = dir(idof)
     enddo
  enddo

  open(lunit,file="direction."//trim(ext), status="unknown", form="formatted" )

  write(lunit,'( 2A23)', advance='no') 'coil', 'I'
  !write(lunit,'(    A)', advance='no') 'XC'
  write(lunit,'(10I23)', advance='no') (i, i=0,NFcoil)
 ! write(lunit,'(    A)', advance='no') 'XS'
  write(lunit,'(10I23)', advance='no') (i, i=0,NFcoil)
  !write(lunit,'(    A)', advance='no') 'YC'
  write(lunit,'(10I23)', advance='no') (i, i=0,NFcoil)
  !write(lunit,'(    A)', advance='no') 'YS'
  write(lunit,'(10I23)', advance='no') (i, i=0,NFcoil)
 ! write(lunit,'(    A)', advance='no') 'ZC'
  write(lunit,'(10I23)', advance='no') (i, i=0,NFcoil)
  !write(lunit,'(    A)', advance='no') 'ZS'
  write(lunit,'(10I23)', advance='no') (i, i=0,NFcoil)
  write(lunit, *) ''

  do icoil = 1, Ncoils     
     write(lunit,'(I23,ES23.15)', advance='no') icoil, coil(icoil,0)
     write(lunit,'(1000ES23.15)', advance='no') (coil(icoil,i), i=1,Cdof)
     write(lunit, *) ''
  enddo
  
  close(lunit)

end subroutine wrtdir
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

SUBROUTINE hybrid
  use kmodule, only: zero, myid, ncpu, ounit, sqrtmachprec, Ntauout, Ndof, Idisplay, coil
  implicit none
  include "mpif.h"

  INTEGER           :: n,ldfjac,maxfev,mode,nprint,info,nfev,njev,lr
  REAL, allocatable :: X(:), FVEC(:), DIAG(:), FJAC(:,:), R(:), QTF(:), W1(:), W2(:), W3(:), W4(:)
  REAL              :: FACTOR, xtol
  INTEGER           :: irestart = 1, astat, ierr

  EXTERNAL          :: fcn, pack


  N = Ndof
  maxfev = Ntauout
  xtol = sqrtmachprec
  LDFJAC = N
  LR     = N*(N+1)/2
  INFO   = Idisplay
  MODE   = 1
  FACTOR = 100.0
  nprint = 1

  SALLOCATE( X   ,(1:N     ), zero)
  SALLOCATE( FVEC,(N       ), zero)
  SALLOCATE( DIAG,(N       ), zero)
  SALLOCATE( FJAC,(LDFJAC,N), zero)
  SALLOCATE( R   ,(LR      ), zero)
  SALLOCATE( QTF ,(N       ), zero)
  SALLOCATE( W1  ,(N       ), zero)
  SALLOCATE( W2  ,(N       ), zero) 
  SALLOCATE( W3  ,(N       ), zero) 
  SALLOCATE( W4  ,(N       ), zero) 

  call pack(x(1:n))
  if (myid .eq. 0) write(ounit, '("hybrid  : "10X" : Begin using hybrid Newton method to solve nonlinear equations.")')
  call costfun(1); call output

  call hybrj(fcn,n,x,fvec,fjac,ldfjac,xtol,maxfev,diag,mode,factor,nprint,info,nfev,njev,r,lr,qtf,w1,w2,w3,w4)

  if (myid .eq. 0) then
     write(ounit, '("hybrid  : "10X" : Hybrid Newton method finished. INFO = "I3)') info
     select case ( info )
     case (0)
        write(ounit, '("hybrid  : "10X" : improper input parameters.")')
     case (1)
        write(ounit, '("hybrid  : "10X" : relative error between two consecutive iterates is at most xtol.")')
     case (2)
        write(ounit, '("hybrid  : "10X" : number of calls to fcn with iflag = 1 has reached maxfev.")')
     case (3)
        write(ounit, '("hybrid  : "10X" : xtol is too small. no further improvement in the approximate solution x is possible.")')
     case (4)
        write(ounit, '("hybrid  : "10X" : iteration is not making good progress, as measured by the improvement from the last five jacobian evaluations.")')
     case (5)
        write(ounit, '("hybrid  : "10X" : iteration is not making good progress, as measured by the improvement from the last ten iterations.")')
     end select
  endif

  deallocate(x, fvec, diag, fjac, r, qtf, w1, w2, w3, w4)

END SUBROUTINE hybrid

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

SUBROUTINE fcn(n, x, fvec, fjac, ldfjac, iflag)
  use kmodule, only: myid, ounit, ncpu, t2E
  implicit none
  include "mpif.h"

  INTEGER, INTENT( IN) :: n,ldfjac,iflag
  REAL,    INTENT( IN) :: x(n)
  REAL,    INTENT(OUT) :: fvec(n),fjac(ldfjac,n)

  INTEGER :: ii, jj, c1, c2, n1, n2, astat, ierr
  REAL    :: f0

  call unpack(x)
  if (iflag == 1) then
     call getdf(f0, fvec)
  else if (iflag == 2) then
     call costfun(2)
     do jj = 1, N
         call DoFconvert(jj,c2,n2)
         do ii = 1,N
            call DoFconvert(ii,c1,n1)
            FJAC(ii,jj) = t2E(c1,n1,c2,n2)
         enddo
     enddo
  else
     call output
  endif

END SUBROUTINE fcn
