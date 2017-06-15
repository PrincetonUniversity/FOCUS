SUBROUTINE congrad
  use kmodule, only: sqrtmachprec, myid, ounit, Ncoils, Ndof, t1E, itau, CG_Niter
  implicit none
  include "mpif.h"

  INTEGER                 :: idof, icoil, c1, n1, ierr, astat
  REAL                    :: alpha, beta, f
  REAL, dimension(1:Ndof) :: xdof, p, gradk, gradf

  if (myid .eq. 0) write(ounit, '("truncnt : "10X" : Begin using Nonlinear Conjugate Gradient to optimize.")')
  
  call getdf(f, gradk)

  call pack(xdof(1:Ndof)) ! initial xdof;
  p(1:Ndof) = -gradk ! initial step direction;
  alpha = 1.0 ! initial step size;
  call output

  do

     call wolfe(xdof, p, alpha) ! find a step size matching the Wolfe condiction;
     xdof = xdof + alpha*p(1:Ndof) ! next xdof
     call unpack(xdof)
     call discretecoil
     call getdf(f, gradf)     
     call output

     if ( sum(gradf**2) .lt. sqrtmachprec) exit  ! reach minimum 
     beta = sum(gradf**2) / sum( (gradf-gradk)*p )
     p = -gradf + beta*p ! direction for next step;
     gradk = gradf  !save for the current step;

     alpha = 1.0  ! reset alpha;

     if (itau .ge. CG_Niter) exit  ! reach maximum iterations;

  enddo

  if(myid .eq. 0) write(ounit, '("congrad : Computation using conjugate gradient finished.")')

  return
END SUBROUTINE congrad

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine output
  
  use kmodule, only : zero, ounit, myid, iter, &
                      NFcoil, Ndof, Ncoils, Cdof, t1E, &
                      itau, tauend, Ntauout, Savfreq, tautol, &
                      coil, icoil, Ncoils, totalenergy, evolution, bnorm, tflux, ttlen, eqarc, ccsep, tbn
  
  implicit none
  
  include "mpif.h"

  REAL                 :: f, sumdE, dE(1:Ndof)
  
  INTEGER              :: idof, irestart, astat, ierr

  irestart = 1

  !call getdf( f, dE)
  
  sumdE = sum(t1E**2)

  if( myid.eq.0 ) write(ounit,1000) itau, totalenergy, sumdE, bnorm, tflux, ttlen, eqarc, ccsep

  if(allocated(evolution)) then
     evolution(itau,0) = itau
     evolution(itau,1) = totalenergy
     evolution(itau,2) = sumdE
     evolution(itau,3) = bnorm
     evolution(itau,4) = tflux
     evolution(itau,5) = ttlen
     evolution(itau,6) = eqarc
     evolution(itau,7) = ccsep
     evolution(itau,8) = coil(1)%I
     evolution(itau,9) = coil(1)%L
  endif
  
  if(mod(itau,Savfreq) .eq. 0) call restart( irestart )

  itau = itau + 1

  return  

1000 format("progres :"i11" : E="es15.7" ; D="es15.7" ; B="es15.7" ; F="es15.7" ; L="es15.7" ; A="es15.7" ; C="es15.7" ;")

end subroutine output

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

SUBROUTINE wolfe( x0, p, alpha )

  use kmodule, only : zero, sqrtmachprec, ounit, myid, Ndof

  implicit none
  include "mpif.h"

  REAL, INTENT( in)       :: x0(1:Ndof), p(1:Ndof)
  REAL, INTENT(out)       :: alpha

  REAL                    :: zoom
  REAL, parameter         :: c1 = 1e-4, c2 = 0.1       ! c1 & c2
  INTEGER                 :: i
  REAL                    :: a0, am, ap, ac, f0, fc, fp, rd
  REAL, dimension(1:Ndof) :: xc, g0, gc, gp


  i = 0
  a0 = 0.0
  am = alpha
  if (am <= 0.0) STOP "alpha_max should be larger than alpha_0"

  call unpack(x0)
  call discretecoil
  call getdf(f0, g0) ! get f0 and df0;
  
  ap = a0            ! previous alpha;
  fp = f0            ! previous fnction;
  gp = g0            ! previous gradient;

!!$  call RANDOM_NUMBER(rd)
!!$  ac = am*rd   ! current alpha;
  ac = alpha

  do
     xc = x0 + ac*p  ! current xdof
     call unpack(xc)
     call discretecoil
     call getdf(fc, gc)

     if ( (fc > f0 + c1*ac*sum(p*g0)) .or. (fc >= fp .and. i > 1) ) then
        !TMPOUT("wolfe: case1")
        !if (myid .eq. 0) print *, ap, ac, am, rd
        alpha = zoom( x0, p, ap, ac )
        return
     endif

     if ( abs(sum(p*gc)) <= -c2*sum(p*g0) ) then
        !TMPOUT("wolfe: case2")
        !if (myid .eq. 0) print *, ap, ac, am, rd
        alpha = ac
        return
     endif

     if ( sum(p*gc) >= zero ) then
        !TMPOUT("wolfe: case3")
        !if (myid .eq. 0) print *, ap, ac
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

!!$     call RANDOM_NUMBER(rd)
!!$     ac = ac + (am-ac)*rd
  ac = 2.0 * ac

     i = i+1
     
  end do

  return

END SUBROUTINE wolfe

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

REAL FUNCTION zoom( x0, p, alo, ahi )

  use kmodule, only : zero, ounit, myid, Ndof

  implicit none
  include "mpif.h"

  REAL, INTENT(   in)     :: x0(1:Ndof), p(1:Ndof)
  REAL, INTENT(inout)     :: alo, ahi

  REAL, parameter         :: c1 = 1e-4, c2 = 0.1       ! c1 & c2
  REAL                    :: f0, fc, fl, fp, alpha
  REAL, dimension(1:Ndof) :: xc, xl, g0, gc, gp, gl

  call unpack(x0)
  call discretecoil
  call getdf(f0, g0) ! get f0 and df0;

  do
     alpha = 0.5*(alo + ahi)
     xc = x0 + alpha*p
     call unpack(xc)
     call discretecoil
     call getdf(fc, gc)

     xl = x0 + alo*p
     call unpack(xl)
     call discretecoil
     call getdf(fl, gl)

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

SUBROUTINE getdf(f ,g)
  use kmodule, only: myid, ounit, Ndof, Ncoils, coil, NFcoil, totalenergy, t1E
  implicit none
  include "mpif.h"

  REAL, INTENT(out) :: f, g(1:Ndof)

  INTEGER           :: mm, icoil, idof, ierr

  call costfun(1)
  f = totalenergy

  idof = 0
  do icoil = 1, Ncoils

     if(coil(icoil)%Ic .ne. 0) then
        idof = idof + 1 ; g(idof) = t1E(icoil,        0) ! coil(icoil)%I
     endif

     if(coil(icoil)%Lc .ne. 0) then
        idof = idof + 1 ; g(idof) = t1E(icoil,        1) ! coil(icoil)%xc( 0)
        idof = idof + 1 ; g(idof) = t1E(icoil,   2*NFcoil+3) ! coil(icoil)%yc( 0)
        idof = idof + 1 ; g(idof) = t1E(icoil,   4*NFcoil+5) ! coil(icoil)%zc( 0)
        do mm = 1, NFcoil 
           idof = idof + 1 ; g(idof) = t1E(icoil,mm+     1) ! coil(icoil)%xc(mm)
           idof = idof + 1 ; g(idof) = t1E(icoil,mm+2*NFcoil+3) ! coil(icoil)%yc(mm)
           idof = idof + 1 ; g(idof) = t1E(icoil,mm+4*NFcoil+5) ! coil(icoil)%zc(mm)
           idof = idof + 1 ; g(idof) = t1E(icoil,mm+  NFcoil+2) ! coil(icoil)%xs(mm)
           idof = idof + 1 ; g(idof) = t1E(icoil,mm+3*NFcoil+4) ! coil(icoil)%ys(mm)
           idof = idof + 1 ; g(idof) = t1E(icoil,mm+5*NFcoil+6) ! coil(icoil)%zs(mm)
        enddo ! end of do mm
     endif

  enddo ! end of do icoil
  FATAL( getdf, idof.ne.Ndof, counting error )

  return
END SUBROUTINE getdf
  
  
