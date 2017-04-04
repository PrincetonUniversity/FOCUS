!title (coilsep) ! Constraints for finding maximum coil-coil separation. 

!latex \briefly{In order to get large accsess to the plasma for diagnostics, it's better to have bigger separations between adjacent coils. 
!latex To find a continuous cost function, we assume two coils are charged electric lines and do a double integral to calculate their electric
!latex potential as the coil-coil separation cost function. We can also try to calculate magnetic forces between two coils which may be more
!latex pratical for engineering consideration.}

!latex \calledby{\link{denergy}}

!latex \tableofcontents

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!latex \subsection{Cost function (0 order)}
!latex For two adjacent coils, coil $i$ and coil $j$, the electric potential between them are written as,
!latex \begin{equation}
!latex C = \sum_{i,j = 1}^{Ncoils}\int_{c_i}\int_{C_j} \frac{dl_i \ dl_j}{\delta r^2}
!latex \end{equation}
!latex where,  $dl_i = \sqrt{\dot{x_i}^2+\dot{y_i}^2+\dot{z_i}^2}$, $dl_j = \sqrt{\dot{x_j}^2+\dot{y_j}^2+\dot{z_j}^2}$ and 
!latex $\delta r = \sqrt{(x_j - x_i)^2+(y_j - y_i)^2+(z_j - z_i)^2}$. And using the square of $\delta r$ actually makes the funtion dimensionless.\\
!latex \emph{Note: Since the differential arc length is included in the cost function $C$, minimizing the function $C$ would lead coils compressed to 
!latex point, as well as increasing their separation. Another idea is summing up the distance between two points of the same poloidal angel on the 
!latex two coils, rather than doing a double integral over the whole coils.}
!latex 
!latex \subsection{First derivatives}
!latex The first derivative of function $C$ on the $x_m$ term of coil $i$ can be written as,
!latex \begin{equation}
!latex \frac{\partial{C}}{\partial{x_m^i}} = \int_{c_{i-1}}\int_{C_i} \frac{dl_{i-1} \ dl_i}{\delta r^2} \; + \; \int_{c_i}\int_{C_{i+1}} \frac{dl_i \ dl_{i+1}}{\delta r^2}
!latex \end{equation}
!latex 
!latex 
!latex 
!latex
!latex 
!latex 
!latex 
!latex
!latex 
!latex 
!latex 
!latex


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine coilsep (nderiv)
  use kmodule, only: zero, mincc, Ncoils, Cdof, ccsep, t1C, t2C, myid, ncpu, ounit
  implicit none
  include "mpif.h"

  INTEGER          :: nderiv
  
  INTEGER          :: icoil, jcoil, ierr, astat,array2size, array4size
  REAL             :: dccsep, lccsep
  REAL, allocatable:: tmp(:), tmp2(:,:), l1C(:,:), l2C(:,:,:,:)

  if( .not. allocated(mincc) ) then
     SALLOCATE(mincc, (1:Ncoils, 1:Ncoils), zero)
  else
     mincc = zero
  endif
  
  select case (nderiv)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  case (0)

     dccsep = zero

     do icoil = 1, Ncoils
        
        if( myid.ne.modulo(icoil-1,ncpu) ) cycle ! parallelization loop; 

        jcoil = icoil+1                 ! assuming coil numbers are in order; 07/06/2016
        if (icoil .eq. Ncoils) jcoil = 1

        call epotent0(icoil, jcoil, lccsep)

        dccsep = dccsep + lccsep

     enddo
     
     call MPI_REDUCE( dccsep, ccsep, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
     RlBCAST(ccsep, 1, 0)
     
     ccsep = ccsep/Ncoils

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  case (1)

     dccsep = zero
     SALLOCATE(tmp, (         0:Cdof), zero)
     SALLOCATE(l1C, (1:Ncoils,0:Cdof), zero)
     if( .not. allocated(t1c)) allocate(t1C(1:Ncoils,0:Cdof))
     t1C = zero

     do icoil = 1, Ncoils

        if( myid.ne.modulo(icoil-1,ncpu) ) cycle ! parallelization loop;

        ! icoil & icoil+1
        jcoil = icoil+1; if (icoil .eq. Ncoils) jcoil = 1             ! assuming coil numbers are in order; 07/06/2016        
        call epotent1(icoil, jcoil, tmp(0:Cdof))
        l1c(icoil,1:Cdof) = l1c(icoil,1:Cdof) + tmp(1:Cdof)

        dccsep = dccsep + tmp(0) ! avoid double calculating

        ! icoil & icoil-1
        jcoil = icoil-1; if (icoil .eq. 1     ) jcoil = Ncoils        ! assuming coil numbers are in order; 07/06/2016        
        call epotent1(icoil, jcoil, tmp(0:Cdof))
        l1c(icoil,1:Cdof) = l1c(icoil,1:Cdof) + tmp(1:Cdof)

     enddo

     call MPI_REDUCE( dccsep, ccsep, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
     RlBCAST(ccsep, 1, 0)

     array2size = Ncoils * ( Cdof + 1 )
     call MPI_REDUCE( l1C, t1C, array2size, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
     RlBCAST( t1C, array2size, 0 )

     ccsep = ccsep/Ncoils
     t1C   = t1C  /Ncoils

     deallocate( tmp, l1C)
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!    

  case (2)

     dccsep = zero
     SALLOCATE(tmp, (         0:Cdof), zero)
     SALLOCATE(tmp2,(0:Cdof  ,0:Cdof), zero)
     SALLOCATE(l1C, (1:Ncoils,0:Cdof), zero)
     SALLOCATE(l2C, (1:Ncoils,0:Cdof,1:Ncoils,0:Cdof), zero)

     if( .not. allocated(t1c)) allocate(t1C(1:Ncoils,0:Cdof))
     if( .not. allocated(t2c)) allocate(t2C(1:Ncoils,0:Cdof,1:Ncoils,0:Cdof))
     t1C = zero; t2C = zero

     do icoil = 1, Ncoils

        if( myid.ne.modulo(icoil-1,ncpu) ) cycle ! parallelization loop;
        
        ! icoil & icoil+1 on icoil,icoil
        jcoil = icoil+1; if (icoil .eq. Ncoils) jcoil = 1             ! assuming coil numbers are in order; 07/07/2016        
        call epotent2(icoil, jcoil, icoil, tmp2(0:Cdof,0:Cdof))
        l1c(icoil,1:Cdof) = l1c(icoil,1:Cdof) + tmp2(1:Cdof,0)
        l2c(icoil,1:Cdof,icoil,1:Cdof) = l2c(icoil,1:Cdof,icoil,1:Cdof) + tmp2(1:Cdof,1:Cdof)

        dccsep = dccsep + tmp2(icoil,0) ! avoid double calculating

        ! icoil & icoil-1 on icoil,icoil
        jcoil = icoil-1; if (icoil .eq. 1     ) jcoil = Ncoils        ! assuming coil numbers are in order; 07/07/2016
        call epotent2(icoil, jcoil, icoil, tmp2(0:Cdof,0:Cdof))
        l1c(icoil,1:Cdof) = l1c(icoil,1:Cdof) + tmp2(1:Cdof,0)
        l2c(icoil,1:Cdof,icoil,1:Cdof) = l2c(icoil,1:Cdof,icoil,1:Cdof) + tmp2(1:Cdof,1:Cdof)    

        ! icoil & icoil+1 on icoil,icoil+1
        jcoil = icoil+1; if (icoil .eq. Ncoils) jcoil = 1             ! assuming coil numbers are in order; 07/07/2016        
        call epotent2(icoil, jcoil, jcoil, tmp2(0:Cdof,0:Cdof))
        l2c(icoil,1:Cdof,jcoil,1:Cdof) =                                  tmp2(1:Cdof,1:Cdof)

        ! icoil & icoil-1 on icoil,icoil-1
        jcoil = icoil-1; if (icoil .eq. 1     ) jcoil = Ncoils        ! assuming coil numbers are in order; 07/07/2016
        call epotent2(icoil, jcoil, jcoil, tmp2(0:Cdof,0:Cdof))
        l2c(icoil,1:Cdof,jcoil,1:Cdof) =                                  tmp2(1:Cdof,1:Cdof)

     enddo  

     call MPI_REDUCE( dccsep, ccsep, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
     RlBCAST(ccsep, 1, 0)

     array2size = Ncoils * ( Cdof + 1 )
     call MPI_REDUCE( l1C, t1C, array2size, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
     RlBCAST( t1C, array2size, 0 )
     
     array4size = array2size * array2size
     call MPI_REDUCE( l2C, t2C, array4size, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
     RlBCAST( t2C, array4size, 0 )
     
     ccsep = ccsep/Ncoils
     t1C   = t1C  /Ncoils 
     t2C   = t2C  /Ncoils

     deallocate( tmp, tmp2, l1c, l2c )
  end select

  return
end subroutine coilsep

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine epotent0(icoil, jcoil, ccpotent)
  use kmodule, only: zero, one, pi2, coil, Ncoils, NDcoil, mincc, ounit, myid
  implicit none
  include "mpif.h"

  INTEGER        :: icoil, jcoil
  REAL           :: ccpotent

  INTEGER        :: ierr, astat, k1, k2
  REAL           :: dli, dlj, dr
  REAL,parameter :: delta = 0.0                           ! delta for the case dr=0 avoiding dividing zero; 


  FATAL(epotent, icoil .lt. 1 .or. icoil .gt. Ncoils, illegal)
  FATAL(epotent, jcoil .lt. 1 .or. jcoil .gt. Ncoils, illegal)
  FATAL(epotent, icoil .eq. jcoil, illegal)
  FATAL(epotent, .not. allocated(mincc), illegal)

  ccpotent = zero;  mincc(icoil, jcoil) = one      ! mincc is set to a number larger than 0;

  do k1 = 0, NDcoil-1
     
     dli = sqrt(coil(icoil)%xt(k1)**2 + coil(icoil)%yt(k1)**2 + coil(icoil)%zt(k1)**2)
     dr  = zero

     do k2 = 0, NDcoil-1

        dlj = sqrt(coil(jcoil)%xt(k2)**2 + coil(jcoil)%yt(k2)**2 + coil(jcoil)%zt(k2)**2)
        dr  = sqrt( (coil(icoil)%xx(k1)-coil(jcoil)%xx(k2))**2 + (coil(icoil)%yy(k1)-coil(jcoil)%yy(k2))**2 + (coil(icoil)%zz(k1)-coil(jcoil)%zz(k2))**2 ) + delta

        mincc(icoil,jcoil) = min(dr, mincc(icoil,jcoil))
        ccpotent = ccpotent + dli*dlj/dr**2 ! if takes a flexible exponent, how to normalize?

     enddo
  enddo

  ccpotent = ccpotent * pi2/NDcoil * pi2/NDcoil

  return
end subroutine epotent0

!!$!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!!$
!!$subroutine epotent1(icoil, jcoil, ccpotent)   ! icoil-jcoil separation on icoil terms
!!$  use kmodule, only: zero, one, pi2, coil, Ncoils, NFcoil, NDcoil, Cdof, mincc, ounit, myid, cmt, smt
!!$  implicit none
!!$  include "mpif.h"
!!$
!!$  INTEGER        :: icoil, jcoil
!!$  REAL           :: ccpotent(0:Cdof)                !derivatives over icoil terms
!!$
!!$  INTEGER        :: ierr, astat, k1, k2, ll, NN
!!$  REAL           :: dli(0:Cdof), dlj(0:Cdof), dr, dri(0:Cdof)
!!$  REAL,parameter :: delta = 0.0                           ! delta for the case dr=0 avoiding dividing zero; 
!!$
!!$  FATAL(epotent, icoil .lt. 1 .or. icoil .gt. Ncoils, illegal)
!!$  FATAL(epotent, jcoil .lt. 1 .or. jcoil .gt. Ncoils, illegal)
!!$  FATAL(epotent, icoil .eq. jcoil, illegal)
!!$
!!$  NN = NFcoil
!!$  ccpotent = zero; mincc(icoil, jcoil) = one      !  mincc is set to a number larger than 0;
!!$
!!$  do k1 = 0, NDcoil-1
!!$     
!!$     dli(0) = sqrt(coil(icoil)%xt(k1)**2 + coil(icoil)%yt(k1)**2 + coil(icoil)%zt(k1)**2)
!!$
!!$     do ll = 0, NN
!!$        dli(ll        + 1) = -ll*smt(k1, ll)*coil(icoil)%xt(k1)/dli(0)
!!$        dli(ll +   NN + 2) =  ll*cmt(k1, ll)*coil(icoil)%xt(k1)/dli(0)
!!$        dli(ll + 2*NN + 3) = -ll*smt(k1, ll)*coil(icoil)%yt(k1)/dli(0)
!!$        dli(ll + 3*NN + 4) =  ll*cmt(k1, ll)*coil(icoil)%yt(k1)/dli(0)
!!$        dli(ll + 4*NN + 5) = -ll*smt(k1, ll)*coil(icoil)%zt(k1)/dli(0)
!!$        dli(ll + 5*NN + 6) =  ll*cmt(k1, ll)*coil(icoil)%zt(k1)/dli(0)
!!$     enddo ! end ll
!!$
!!$     do k2 = 0, NDcoil-1
!!$
!!$        dlj(0) = sqrt(coil(jcoil)%xt(k2)**2 + coil(jcoil)%yt(k2)**2 + coil(jcoil)%zt(k2)**2)
!!$        dr  = sqrt( (coil(icoil)%xx(k1)-coil(jcoil)%xx(k2))**2 + (coil(icoil)%yy(k1)-coil(jcoil)%yy(k2))**2 + (coil(icoil)%zz(k1)-coil(jcoil)%zz(k2))**2 ) + delta
!!$        mincc(icoil,jcoil) = min(dr, mincc(icoil,jcoil))
!!$
!!$        do ll = 0, NN
!!$           dri(ll        + 1) = cmt(k1,ll)*(coil(icoil)%xx(k1)-coil(jcoil)%xx(k2))/dr
!!$           dri(ll +   NN + 2) = smt(k1,ll)*(coil(icoil)%xx(k1)-coil(jcoil)%xx(k2))/dr
!!$           dri(ll + 2*NN + 3) = cmt(k1,ll)*(coil(icoil)%yy(k1)-coil(jcoil)%yy(k2))/dr
!!$           dri(ll + 3*NN + 4) = smt(k1,ll)*(coil(icoil)%yy(k1)-coil(jcoil)%yy(k2))/dr
!!$           dri(ll + 4*NN + 5) = cmt(k1,ll)*(coil(icoil)%zz(k1)-coil(jcoil)%zz(k2))/dr
!!$           dri(ll + 5*NN + 6) = smt(k1,ll)*(coil(icoil)%zz(k1)-coil(jcoil)%zz(k2))/dr
!!$        enddo ! end ll
!!$
!!$        ccpotent(0) = ccpotent(0) + dli(0)*dlj(0)/dr**2 ! if takes a flexible exponent, how to normalize?
!!$        do ll = 1, Cdof
!!$           ccpotent(ll) = ccpotent(ll) + dlj(0)*dli(ll)/dr**2 - 2*dli(0)*dlj(0)*dri(ll)/dr**3
!!$        enddo
!!$       
!!$     enddo
!!$  enddo
!!$
!!$  ccpotent = ccpotent * pi2/NDcoil * pi2/NDcoil
!!$
!!$  return
!!$end subroutine epotent1

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine epotent1(icoil, jcoil, ccpotent)   ! icoil-jcoil separation on icoil terms
  use kmodule, only: zero, one, pi2, coil, Ncoils, NFcoil, NDcoil, Cdof, cmt, smt, mincc, ounit, myid
  implicit none
  include "mpif.h"

  INTEGER        :: icoil, jcoil
  REAL           :: ccpotent(0:Cdof)                !derivatives over icoil terms

  INTEGER        :: ierr, astat, k1, k2, ll, NN
  REAL           :: tix, tiy, tiz, dix, diy, diz, tjx, tjy, tjz, djx, djy, djz, dr, dli, dlj
  REAL,parameter :: delta = 0.0                           ! delta for the case dr=0 avoiding dividing zero; 

  FATAL(epotent, icoil .lt. 1 .or. icoil .gt. Ncoils, illegal)
  FATAL(epotent, jcoil .lt. 1 .or. jcoil .gt. Ncoils, illegal)
  FATAL(epotent, icoil .eq. jcoil, illegal)
  FATAL(epotent, .not. allocated(mincc), illegal)

  NN = NFcoil
  ccpotent = zero; mincc(icoil, jcoil) = one      !  mincc is set to a number larger than 0;

  do k1 = 0, NDcoil-1

     tix = coil(icoil)%xt(k1); dix = coil(icoil)%xx(k1)
     tiy = coil(icoil)%yt(k1); diy = coil(icoil)%yy(k1)
     tiz = coil(icoil)%zt(k1); diz = coil(icoil)%zz(k1)

     dli = sqrt( tix**2 + tiy**2 + tiz**2 )
    
     do k2 = 0, NDcoil-1

        tjx = coil(jcoil)%xt(k2); djx = coil(jcoil)%xx(k2)
        tjy = coil(jcoil)%yt(k2); djy = coil(jcoil)%yy(k2)
        tjz = coil(jcoil)%zt(k2); djz = coil(jcoil)%zz(k2)

        dlj = sqrt( tjx**2 + tjy**2 + tjz**2 )
        dr  = sqrt( (djx-dix)**2 + (djy-diy)**2 + (djz-diz)**2 ) + delta
        mincc(icoil,jcoil) = min(dr, mincc(icoil,jcoil))

        ccpotent(0) = ccpotent(0) + dli*dlj/dr**2 ! if takes a flexible exponent, how to normalize?

        do ll = 0, NN
           ccpotent(ll     +1) = ccpotent(ll     +1) - ll*smt(k1,ll)*dlj*tix/(dr**2*dli) - 2*dli*dlj*(dix-djx)*cmt(k1,ll)/dr**4
           ccpotent(ll+  NN+2) = ccpotent(ll+  NN+2) + ll*cmt(k1,ll)*dlj*tix/(dr**2*dli) - 2*dli*dlj*(dix-djx)*smt(k1,ll)/dr**4

           ccpotent(ll+2*NN+3) = ccpotent(ll+2*NN+3) - ll*smt(k1,ll)*dlj*tiy/(dr**2*dli) - 2*dli*dlj*(diy-djy)*cmt(k1,ll)/dr**4
           ccpotent(ll+3*NN+4) = ccpotent(ll+3*NN+4) + ll*cmt(k1,ll)*dlj*tiy/(dr**2*dli) - 2*dli*dlj*(diy-djy)*smt(k1,ll)/dr**4

           ccpotent(ll+4*NN+5) = ccpotent(ll+4*NN+5) - ll*smt(k1,ll)*dlj*tiz/(dr**2*dli) - 2*dli*dlj*(diz-djz)*cmt(k1,ll)/dr**4
           ccpotent(ll+5*NN+6) = ccpotent(ll+5*NN+6) + ll*cmt(k1,ll)*dlj*tiz/(dr**2*dli) - 2*dli*dlj*(diz-djz)*smt(k1,ll)/dr**4

        enddo ! end ll
       
     enddo
  enddo

  ccpotent = ccpotent * pi2/NDcoil * pi2/NDcoil

  return
end subroutine epotent1

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine epotent2(icoil, jcoil, mcoil, ccpotent)   ! icoil-jcoil separation on icoil & mcoil terms
  use kmodule, only: zero, one, pi2, coil, Ncoils, NFcoil, NDcoil, Cdof, cmt, smt, mincc, ounit, myid
  implicit none
  include "mpif.h"

  INTEGER        :: icoil, jcoil, mcoil
  REAL           :: ccpotent(0:Cdof,0:Cdof)                !derivatives over icoil & mcoil terms

  INTEGER        :: ierr, astat, k1, k2, ll, mm, NN
  REAL           :: tix, tiy, tiz, dix, diy, diz, tjx, tjy, tjz, djx, djy, djz, dr, dli, dlj
  REAL,parameter :: delta = 0.0                           ! delta for the case dr=0 avoiding dividing zero; 


  FATAL(epotent, icoil .lt. 1 .or. icoil .gt. Ncoils, illegal)
  FATAL(epotent, jcoil .lt. 1 .or. jcoil .gt. Ncoils, illegal)
  FATAL(epotent, icoil .eq. jcoil, illegal)
  FATAL(epotent, .not. allocated(mincc), illegal)
  FATAL(epotent, mcoil .ne. icoil .and. mcoil .ne. jcoil, zero derivatives)

  NN = NFcoil
  ccpotent = zero;  mincc(icoil, jcoil) = one     ! mincc is set to a number larger than 0;

  do k1 = 0, NDcoil-1

     tix = coil(icoil)%xt(k1); dix = coil(icoil)%xx(k1)
     tiy = coil(icoil)%yt(k1); diy = coil(icoil)%yy(k1)
     tiz = coil(icoil)%zt(k1); diz = coil(icoil)%zz(k1)

     dli = sqrt( tix**2 + tiy**2 + tiz**2 )
    
     do k2 = 0, NDcoil-1

        tjx = coil(jcoil)%xt(k2); djx = coil(jcoil)%xx(k2)
        tjy = coil(jcoil)%yt(k2); djy = coil(jcoil)%yy(k2)
        tjz = coil(jcoil)%zt(k2); djz = coil(jcoil)%zz(k2)

        dlj = sqrt( tjx**2 + tjy**2 + tjz**2 )
        dr  = sqrt( (djx-dix)**2 + (djy-diy)**2 + (djz-diz)**2 ) + delta
        mincc(icoil,jcoil) = min(dr, mincc(icoil,jcoil))

        ccpotent(0,0) = ccpotent(0,0) + dli*dlj/dr**2 ! if takes a flexible exponent, how to normalize?

        do ll = 0, NN
           ccpotent(ll     +1,0) = ccpotent(ll     +1,0) - ll*smt(k1,ll)*dlj*tix/(dr**2*dli) - 2*dli*dlj*(dix-djx)*cmt(k1,ll)/dr**4
           ccpotent(ll+  NN+2,0) = ccpotent(ll+  NN+2,0) + ll*cmt(k1,ll)*dlj*tix/(dr**2*dli) - 2*dli*dlj*(dix-djx)*smt(k1,ll)/dr**4

           ccpotent(ll+2*NN+3,0) = ccpotent(ll+2*NN+3,0) - ll*smt(k1,ll)*dlj*tiy/(dr**2*dli) - 2*dli*dlj*(diy-djy)*cmt(k1,ll)/dr**4
           ccpotent(ll+3*NN+4,0) = ccpotent(ll+3*NN+4,0) + ll*cmt(k1,ll)*dlj*tiy/(dr**2*dli) - 2*dli*dlj*(diy-djy)*smt(k1,ll)/dr**4

           ccpotent(ll+4*NN+5,0) = ccpotent(ll+4*NN+5,0) - ll*smt(k1,ll)*dlj*tiz/(dr**2*dli) - 2*dli*dlj*(diz-djz)*cmt(k1,ll)/dr**4
           ccpotent(ll+5*NN+6,0) = ccpotent(ll+5*NN+6,0) + ll*cmt(k1,ll)*dlj*tiz/(dr**2*dli) - 2*dli*dlj*(diz-djz)*smt(k1,ll)/dr**4
        enddo ! end ll

        if ( mcoil .eq. jcoil) then  ! over jcoil term

           do mm = 0, NN
              do ll = 0, NN

                 !--------------------------------------------cc / xcj ----------------------------------------------------------------------------------------

                 ccpotent(ll     +1,mm     +1) = ccpotent(ll     +1,mm     +1) + (-ll*smt(k1,ll))*(-mm*smt(k2,mm))*tix*tjx/(dr**2*dli*dlj)   &      ! cc/xci/xcj
                                                                               - (-ll*smt(k1,ll))*cmt(k2,mm)*2*dlj*tix*(djx-dix)/(dr**4*dli) &
                                                                               - (-mm*smt(k2,mm))*cmt(k1,ll)*2*dli*tjx*(dix-djx)/(dr**4*dlj) &
                                                                               + (    cmt(k2,mm))*cmt(k1,ll)*2*dli*dlj          /(dr**4    ) &
                                                                               + (    cmt(k2,mm))*cmt(k1,ll)*8*dli*dlj*(dix-djx)*(djx-dix)/(dr**6 )

                 ccpotent(ll+  NN+2,mm     +1) = ccpotent(ll+  NN+2,mm     +1) + ( ll*cmt(k1,ll))*(-mm*smt(k2,mm))*tix*tjx/(dr**2*dli*dlj)   &      ! cc/xsi/xcj
                                                                               - ( ll*cmt(k1,ll))*cmt(k2,mm)*2*dlj*tix*(djx-dix)/(dr**4*dli) &
                                                                               - (-mm*smt(k2,mm))*smt(k1,ll)*2*dli*tjx*(dix-djx)/(dr**4*dlj) &
                                                                               + (    cmt(k2,mm))*smt(k1,ll)*2*dli*dlj          /(dr**4    ) &
                                                                               + (    cmt(k2,mm))*smt(k1,ll)*8*dli*dlj*(dix-djx)*(djx-dix)/(dr**6 )

                 ccpotent(ll+2*NN+3,mm     +1) = ccpotent(ll+2*NN+3,mm     +1) + (-ll*smt(k1,ll))*(-mm*smt(k2,mm))*tiy*tjx/(dr**2*dli*dlj)   &      ! cc/yci/xcj
                                                                               - (-ll*smt(k1,ll))*cmt(k2,mm)*2*dlj*tiy*(djx-dix)/(dr**4*dli) &
                                                                               - (-mm*smt(k2,mm))*cmt(k1,ll)*2*dli*tjx*(diy-djy)/(dr**4*dlj) &
                                                                              ! + (    cmt(k2,mm))*cmt(k1,ll)*2*dli*dlj          /(dr**4    ) &
                                                                               + (    cmt(k2,mm))*cmt(k1,ll)*8*dli*dlj*(diy-djy)*(djx-dix)/(dr**6 )

                 ccpotent(ll+3*NN+4,mm     +1) = ccpotent(ll+3*NN+4,mm     +1) + ( ll*cmt(k1,ll))*(-mm*smt(k2,mm))*tiy*tjx/(dr**2*dli*dlj)   &      ! cc/ysi/xcj
                                                                               - ( ll*cmt(k1,ll))*cmt(k2,mm)*2*dlj*tiy*(djx-dix)/(dr**4*dli) &
                                                                               - (-mm*smt(k2,mm))*smt(k1,ll)*2*dli*tjx*(diy-djy)/(dr**4*dlj) &
                                                                              ! + (    cmt(k2,mm))*smt(k1,ll)*2*dli*dlj          /(dr**4    ) &
                                                                               + (    cmt(k2,mm))*smt(k1,ll)*8*dli*dlj*(diy-djy)*(djx-dix)/(dr**6 )

                 ccpotent(ll+4*NN+5,mm     +1) = ccpotent(ll+4*NN+5,mm     +1) + (-ll*smt(k1,ll))*(-mm*smt(k2,mm))*tiz*tjx/(dr**2*dli*dlj)   &      ! cc/zci/xcj
                                                                               - (-ll*smt(k1,ll))*cmt(k2,mm)*2*dlj*tiz*(djx-dix)/(dr**4*dli) &
                                                                               - (-mm*smt(k2,mm))*cmt(k1,ll)*2*dli*tjx*(diz-djz)/(dr**4*dlj) &
                                                                              ! + (    cmt(k2,mm))*cmt(k1,ll)*2*dli*dlj          /(dr**4    ) &
                                                                               + (    cmt(k2,mm))*cmt(k1,ll)*8*dli*dlj*(diz-djz)*(djx-dix)/(dr**6 )

                 ccpotent(ll+5*NN+6,mm     +1) = ccpotent(ll+5*NN+6,mm     +1) + ( ll*cmt(k1,ll))*(-mm*smt(k2,mm))*tiz*tjx/(dr**2*dli*dlj)   &      ! cc/zsi/xcj
                                                                               - ( ll*cmt(k1,ll))*cmt(k2,mm)*2*dlj*tiz*(djx-dix)/(dr**4*dli) &
                                                                               - (-mm*smt(k2,mm))*smt(k1,ll)*2*dli*tjx*(diz-djz)/(dr**4*dlj) &
                                                                              ! + (    cmt(k2,mm))*smt(k1,ll)*2*dli*dlj          /(dr**4    ) &
                                                                               + (    cmt(k2,mm))*smt(k1,ll)*8*dli*dlj*(diz-djz)*(djx-dix)/(dr**6 )

                 !--------------------------------------------cc / xsj ----------------------------------------------------------------------------------------

                 ccpotent(ll     +1,mm+  NN+2) = ccpotent(ll     +1,mm+  NN+2) + (-ll*smt(k1,ll))*( mm*cmt(k2,mm))*tix*tjx/(dr**2*dli*dlj)   &      ! cc/xci/xsj
                                                                               - (-ll*smt(k1,ll))*smt(k2,mm)*2*dlj*tix*(djx-dix)/(dr**4*dli) &
                                                                               - ( mm*cmt(k2,mm))*cmt(k1,ll)*2*dli*tjx*(dix-djx)/(dr**4*dlj) &
                                                                               + (    smt(k2,mm))*cmt(k1,ll)*2*dli*dlj          /(dr**4    ) &
                                                                               + (    smt(k2,mm))*cmt(k1,ll)*8*dli*dlj*(dix-djx)*(djx-dix)/(dr**6 )

                 ccpotent(ll+  NN+2,mm+  NN+2) = ccpotent(ll+  NN+2,mm+  NN+2) + ( ll*cmt(k1,ll))*( mm*cmt(k2,mm))*tix*tjx/(dr**2*dli*dlj)   &      ! cc/xsi/xsj
                                                                               - ( ll*cmt(k1,ll))*smt(k2,mm)*2*dlj*tix*(djx-dix)/(dr**4*dli) &
                                                                               - ( mm*cmt(k2,mm))*smt(k1,ll)*2*dli*tjx*(dix-djx)/(dr**4*dlj) &
                                                                               + (    smt(k2,mm))*smt(k1,ll)*2*dli*dlj          /(dr**4    ) &
                                                                               + (    smt(k2,mm))*smt(k1,ll)*8*dli*dlj*(dix-djx)*(djx-dix)/(dr**6 )

                 ccpotent(ll+2*NN+3,mm+  NN+2) = ccpotent(ll+2*NN+3,mm+  NN+2) + (-ll*smt(k1,ll))*( mm*cmt(k2,mm))*tiy*tjx/(dr**2*dli*dlj)   &      ! cc/yci/xsj
                                                                               - (-ll*smt(k1,ll))*smt(k2,mm)*2*dlj*tiy*(djx-dix)/(dr**4*dli) &
                                                                               - ( mm*cmt(k2,mm))*cmt(k1,ll)*2*dli*tjx*(diy-djy)/(dr**4*dlj) &
                                                                               !+ (    smt(k2,mm))*cmt(k1,ll)*2*dli*dlj          /(dr**4    ) &
                                                                               + (    smt(k2,mm))*cmt(k1,ll)*8*dli*dlj*(diy-djy)*(djx-dix)/(dr**6 )

                 ccpotent(ll+3*NN+4,mm+  NN+2) = ccpotent(ll+3*NN+4,mm+  NN+2) + ( ll*cmt(k1,ll))*( mm*cmt(k2,mm))*tiy*tjx/(dr**2*dli*dlj)   &      ! cc/ysi/xsj
                                                                               - ( ll*cmt(k1,ll))*smt(k2,mm)*2*dlj*tiy*(djx-dix)/(dr**4*dli) &
                                                                               - ( mm*cmt(k2,mm))*smt(k1,ll)*2*dli*tjx*(diy-djy)/(dr**4*dlj) &
                                                                              ! + (    smt(k2,mm))*smt(k1,ll)*2*dli*dlj          /(dr**4    ) &
                                                                               + (    smt(k2,mm))*smt(k1,ll)*8*dli*dlj*(diy-djy)*(djx-dix)/(dr**6 )

                 ccpotent(ll+4*NN+5,mm+  NN+2) = ccpotent(ll+4*NN+5,mm+  NN+2) + (-ll*smt(k1,ll))*( mm*cmt(k2,mm))*tiz*tjx/(dr**2*dli*dlj)   &      ! cc/zci/xsj
                                                                               - (-ll*smt(k1,ll))*smt(k2,mm)*2*dlj*tiz*(djx-dix)/(dr**4*dli) &
                                                                               - ( mm*cmt(k2,mm))*cmt(k1,ll)*2*dli*tjx*(diz-djz)/(dr**4*dlj) &
                                                                               !+ (    smt(k2,mm))*cmt(k1,ll)*2*dli*dlj          /(dr**4    ) &
                                                                               + (    smt(k2,mm))*cmt(k1,ll)*8*dli*dlj*(diz-djz)*(djx-dix)/(dr**6 )

                 ccpotent(ll+5*NN+6,mm+  NN+2) = ccpotent(ll+5*NN+6,mm+  NN+2) + ( ll*cmt(k1,ll))*( mm*cmt(k2,mm))*tiz*tjx/(dr**2*dli*dlj)   &      ! cc/zsi/xsj
                                                                               - ( ll*cmt(k1,ll))*smt(k2,mm)*2*dlj*tiz*(djx-dix)/(dr**4*dli) &
                                                                               - ( mm*cmt(k2,mm))*smt(k1,ll)*2*dli*tjx*(diz-djz)/(dr**4*dlj) &
                                                                               !+ (    smt(k2,mm))*smt(k1,ll)*2*dli*dlj          /(dr**4    ) &
                                                                               + (    smt(k2,mm))*smt(k1,ll)*8*dli*dlj*(diz-djz)*(djx-dix)/(dr**6 )

                 !--------------------------------------------cc / ycj ----------------------------------------------------------------------------------------

                 ccpotent(ll     +1,mm+2*NN+3) = ccpotent(ll     +1,mm+2*NN+3) + (-ll*smt(k1,ll))*(-mm*smt(k2,mm))*tix*tjy/(dr**2*dli*dlj)   &      ! cc/xci/ycj
                                                                               - (-ll*smt(k1,ll))*cmt(k2,mm)*2*dlj*tix*(djy-diy)/(dr**4*dli) &
                                                                               - (-mm*smt(k2,mm))*cmt(k1,ll)*2*dli*tjy*(dix-djx)/(dr**4*dlj) &
                                                                              ! + (    cmt(k2,mm))*cmt(k1,ll)*2*dli*dlj          /(dr**4    ) &
                                                                               + (    cmt(k2,mm))*cmt(k1,ll)*8*dli*dlj*(dix-djx)*(djy-diy)/(dr**6 )

                 ccpotent(ll+  NN+2,mm+2*NN+3) = ccpotent(ll+  NN+2,mm+2*NN+3) + ( ll*cmt(k1,ll))*(-mm*smt(k2,mm))*tix*tjy/(dr**2*dli*dlj)   &      ! cc/xsi/ycj
                                                                               - ( ll*cmt(k1,ll))*cmt(k2,mm)*2*dlj*tix*(djy-diy)/(dr**4*dli) &
                                                                               - (-mm*smt(k2,mm))*smt(k1,ll)*2*dli*tjy*(dix-djx)/(dr**4*dlj) &
                                                                               !+ (    cmt(k2,mm))*smt(k1,ll)*2*dli*dlj          /(dr**4    ) &
                                                                               + (    cmt(k2,mm))*smt(k1,ll)*8*dli*dlj*(dix-djx)*(djy-diy)/(dr**6 )

                 ccpotent(ll+2*NN+3,mm+2*NN+3) = ccpotent(ll+2*NN+3,mm+2*NN+3) + (-ll*smt(k1,ll))*(-mm*smt(k2,mm))*tiy*tjy/(dr**2*dli*dlj)   &      ! cc/yci/ycj
                                                                               - (-ll*smt(k1,ll))*cmt(k2,mm)*2*dlj*tiy*(djy-diy)/(dr**4*dli) &
                                                                               - (-mm*smt(k2,mm))*cmt(k1,ll)*2*dli*tjy*(diy-djy)/(dr**4*dlj) &
                                                                               + (    cmt(k2,mm))*cmt(k1,ll)*2*dli*dlj          /(dr**4    ) &
                                                                               + (    cmt(k2,mm))*cmt(k1,ll)*8*dli*dlj*(diy-djy)*(djy-diy)/(dr**6 )

                 ccpotent(ll+3*NN+4,mm+2*NN+3) = ccpotent(ll+3*NN+4,mm+2*NN+3) + ( ll*cmt(k1,ll))*(-mm*smt(k2,mm))*tiy*tjy/(dr**2*dli*dlj)   &      ! cc/ysi/ycj
                                                                               - ( ll*cmt(k1,ll))*cmt(k2,mm)*2*dlj*tiy*(djy-diy)/(dr**4*dli) &
                                                                               - (-mm*smt(k2,mm))*smt(k1,ll)*2*dli*tjy*(diy-djy)/(dr**4*dlj) &
                                                                               + (    cmt(k2,mm))*smt(k1,ll)*2*dli*dlj          /(dr**4    ) &
                                                                               + (    cmt(k2,mm))*smt(k1,ll)*8*dli*dlj*(diy-djy)*(djy-diy)/(dr**6 )

                 ccpotent(ll+4*NN+5,mm+2*NN+3) = ccpotent(ll+4*NN+5,mm+2*NN+3) + (-ll*smt(k1,ll))*(-mm*smt(k2,mm))*tiz*tjy/(dr**2*dli*dlj)   &      ! cc/zci/ycj
                                                                               - (-ll*smt(k1,ll))*cmt(k2,mm)*2*dlj*tiz*(djy-diy)/(dr**4*dli) &
                                                                               - (-mm*smt(k2,mm))*cmt(k1,ll)*2*dli*tjy*(diz-djz)/(dr**4*dlj) &
                                                                              ! + (    cmt(k2,mm))*cmt(k1,ll)*2*dli*dlj          /(dr**4    ) &
                                                                               + (    cmt(k2,mm))*cmt(k1,ll)*8*dli*dlj*(diz-djz)*(djy-diy)/(dr**6 )

                 ccpotent(ll+5*NN+6,mm+2*NN+3) = ccpotent(ll+5*NN+6,mm+2*NN+3) + ( ll*cmt(k1,ll))*(-mm*smt(k2,mm))*tiz*tjy/(dr**2*dli*dlj)   &      ! cc/zsi/ycj
                                                                               - ( ll*cmt(k1,ll))*cmt(k2,mm)*2*dlj*tiz*(djy-diy)/(dr**4*dli) &
                                                                               - (-mm*smt(k2,mm))*smt(k1,ll)*2*dli*tjy*(diz-djz)/(dr**4*dlj) &
                                                                              ! + (    cmt(k2,mm))*smt(k1,ll)*2*dli*dlj          /(dr**4    ) &
                                                                               + (    cmt(k2,mm))*smt(k1,ll)*8*dli*dlj*(diz-djz)*(djy-diy)/(dr**6 )

                 !--------------------------------------------cc / ysj ----------------------------------------------------------------------------------------

                 ccpotent(ll     +1,mm+3*NN+4) = ccpotent(ll     +1,mm+3*NN+4) + (-ll*smt(k1,ll))*( mm*cmt(k2,mm))*tix*tjy/(dr**2*dli*dlj)   &      ! cc/xci/ysj
                                                                               - (-ll*smt(k1,ll))*smt(k2,mm)*2*dlj*tix*(djy-diy)/(dr**4*dli) &
                                                                               - ( mm*cmt(k2,mm))*cmt(k1,ll)*2*dli*tjy*(dix-djx)/(dr**4*dlj) &
                                                                              ! + (    smt(k2,mm))*cmt(k1,ll)*2*dli*dlj          /(dr**4    ) &
                                                                               + (    smt(k2,mm))*cmt(k1,ll)*8*dli*dlj*(dix-djx)*(djy-diy)/(dr**6 )

                 ccpotent(ll+  NN+2,mm+3*NN+4) = ccpotent(ll+  NN+2,mm+3*NN+4) + ( ll*cmt(k1,ll))*( mm*cmt(k2,mm))*tix*tjy/(dr**2*dli*dlj)   &      ! cc/xsi/ysj
                                                                               - ( ll*cmt(k1,ll))*smt(k2,mm)*2*dlj*tix*(djy-diy)/(dr**4*dli) &
                                                                               - ( mm*cmt(k2,mm))*smt(k1,ll)*2*dli*tjy*(dix-djx)/(dr**4*dlj) &
                                                                              ! + (    smt(k2,mm))*smt(k1,ll)*2*dli*dlj          /(dr**4    ) &
                                                                               + (    smt(k2,mm))*smt(k1,ll)*8*dli*dlj*(dix-djx)*(djy-diy)/(dr**6 )

                 ccpotent(ll+2*NN+3,mm+3*NN+4) = ccpotent(ll+2*NN+3,mm+3*NN+4) + (-ll*smt(k1,ll))*( mm*cmt(k2,mm))*tiy*tjy/(dr**2*dli*dlj)   &      ! cc/yci/ysj
                                                                               - (-ll*smt(k1,ll))*smt(k2,mm)*2*dlj*tiy*(djy-diy)/(dr**4*dli) &
                                                                               - ( mm*cmt(k2,mm))*cmt(k1,ll)*2*dli*tjy*(diy-djy)/(dr**4*dlj) &
                                                                               + (    smt(k2,mm))*cmt(k1,ll)*2*dli*dlj          /(dr**4    ) &
                                                                               + (    smt(k2,mm))*cmt(k1,ll)*8*dli*dlj*(diy-djy)*(djy-diy)/(dr**6 )

                 ccpotent(ll+3*NN+4,mm+3*NN+4) = ccpotent(ll+3*NN+4,mm+3*NN+4) + ( ll*cmt(k1,ll))*( mm*cmt(k2,mm))*tiy*tjy/(dr**2*dli*dlj)   &      ! cc/ysi/ysj
                                                                               - ( ll*cmt(k1,ll))*smt(k2,mm)*2*dlj*tiy*(djy-diy)/(dr**4*dli) &
                                                                               - ( mm*cmt(k2,mm))*smt(k1,ll)*2*dli*tjy*(diy-djy)/(dr**4*dlj) &
                                                                               + (    smt(k2,mm))*smt(k1,ll)*2*dli*dlj          /(dr**4    ) &
                                                                               + (    smt(k2,mm))*smt(k1,ll)*8*dli*dlj*(diy-djy)*(djy-diy)/(dr**6 )

                 ccpotent(ll+4*NN+5,mm+3*NN+4) = ccpotent(ll+4*NN+5,mm+3*NN+4) + (-ll*smt(k1,ll))*( mm*cmt(k2,mm))*tiz*tjy/(dr**2*dli*dlj)   &      ! cc/zci/ysj
                                                                               - (-ll*smt(k1,ll))*smt(k2,mm)*2*dlj*tiz*(djy-diy)/(dr**4*dli) &
                                                                               - ( mm*cmt(k2,mm))*cmt(k1,ll)*2*dli*tjy*(diz-djz)/(dr**4*dlj) &
                                                                              ! + (    smt(k2,mm))*cmt(k1,ll)*2*dli*dlj          /(dr**4    ) &
                                                                               + (    smt(k2,mm))*cmt(k1,ll)*8*dli*dlj*(diz-djz)*(djy-diy)/(dr**6 )

                 ccpotent(ll+5*NN+6,mm+3*NN+4) = ccpotent(ll+5*NN+6,mm+3*NN+4) + ( ll*cmt(k1,ll))*( mm*cmt(k2,mm))*tiz*tjy/(dr**2*dli*dlj)   &      ! cc/zsi/ysj
                                                                               - ( ll*cmt(k1,ll))*smt(k2,mm)*2*dlj*tiz*(djy-diy)/(dr**4*dli) &
                                                                               - ( mm*cmt(k2,mm))*smt(k1,ll)*2*dli*tjy*(diz-djz)/(dr**4*dlj) &
                                                                              ! + (    smt(k2,mm))*smt(k1,ll)*2*dli*dlj          /(dr**4    ) &
                                                                               + (    smt(k2,mm))*smt(k1,ll)*8*dli*dlj*(diz-djz)*(djy-diy)/(dr**6 )

                 !--------------------------------------------cc / zcj ----------------------------------------------------------------------------------------

                 ccpotent(ll     +1,mm+4*NN+5) = ccpotent(ll     +1,mm+4*NN+5) + (-ll*smt(k1,ll))*(-mm*smt(k2,mm))*tix*tjz/(dr**2*dli*dlj)   &      ! cc/xci/zcj
                                                                               - (-ll*smt(k1,ll))*cmt(k2,mm)*2*dlj*tix*(djz-diz)/(dr**4*dli) &
                                                                               - (-mm*smt(k2,mm))*cmt(k1,ll)*2*dli*tjz*(dix-djx)/(dr**4*dlj) &
                                                                              ! + (    cmt(k2,mm))*cmt(k1,ll)*2*dli*dlj          /(dr**4    ) &
                                                                               + (    cmt(k2,mm))*cmt(k1,ll)*8*dli*dlj*(dix-djx)*(djz-diz)/(dr**6 )

                 ccpotent(ll+  NN+2,mm+4*NN+5) = ccpotent(ll+  NN+2,mm+4*NN+5) + ( ll*cmt(k1,ll))*(-mm*smt(k2,mm))*tix*tjz/(dr**2*dli*dlj)   &      ! cc/xsi/zcj
                                                                               - ( ll*cmt(k1,ll))*cmt(k2,mm)*2*dlj*tix*(djz-diz)/(dr**4*dli) &
                                                                               - (-mm*smt(k2,mm))*smt(k1,ll)*2*dli*tjz*(dix-djx)/(dr**4*dlj) &
                                                                              ! + (    cmt(k2,mm))*smt(k1,ll)*2*dli*dlj          /(dr**4    ) &
                                                                               + (    cmt(k2,mm))*smt(k1,ll)*8*dli*dlj*(dix-djx)*(djz-diz)/(dr**6 )

                 ccpotent(ll+2*NN+3,mm+4*NN+5) = ccpotent(ll+2*NN+3,mm+4*NN+5) + (-ll*smt(k1,ll))*(-mm*smt(k2,mm))*tiy*tjz/(dr**2*dli*dlj)   &      ! cc/yci/zcj
                                                                               - (-ll*smt(k1,ll))*cmt(k2,mm)*2*dlj*tiy*(djz-diz)/(dr**4*dli) &
                                                                               - (-mm*smt(k2,mm))*cmt(k1,ll)*2*dli*tjz*(diy-djy)/(dr**4*dlj) &
                                                                              ! + (    cmt(k2,mm))*cmt(k1,ll)*2*dli*dlj          /(dr**4    ) &
                                                                               + (    cmt(k2,mm))*cmt(k1,ll)*8*dli*dlj*(diy-djy)*(djz-diz)/(dr**6 )

                 ccpotent(ll+3*NN+4,mm+4*NN+5) = ccpotent(ll+3*NN+4,mm+4*NN+5) + ( ll*cmt(k1,ll))*(-mm*smt(k2,mm))*tiy*tjz/(dr**2*dli*dlj)   &      ! cc/ysi/zcj
                                                                               - ( ll*cmt(k1,ll))*cmt(k2,mm)*2*dlj*tiy*(djz-diz)/(dr**4*dli) &
                                                                               - (-mm*smt(k2,mm))*smt(k1,ll)*2*dli*tjz*(diy-djy)/(dr**4*dlj) &
                                                                              ! + (    cmt(k2,mm))*smt(k1,ll)*2*dli*dlj          /(dr**4    ) &
                                                                               + (    cmt(k2,mm))*smt(k1,ll)*8*dli*dlj*(diy-djy)*(djz-diz)/(dr**6 )

                 ccpotent(ll+4*NN+5,mm+4*NN+5) = ccpotent(ll+4*NN+5,mm+4*NN+5) + (-ll*smt(k1,ll))*(-mm*smt(k2,mm))*tiz*tjz/(dr**2*dli*dlj)   &      ! cc/zci/zcj
                                                                               - (-ll*smt(k1,ll))*cmt(k2,mm)*2*dlj*tiz*(djz-diz)/(dr**4*dli) &
                                                                               - (-mm*smt(k2,mm))*cmt(k1,ll)*2*dli*tjz*(diz-djz)/(dr**4*dlj) &
                                                                               + (    cmt(k2,mm))*cmt(k1,ll)*2*dli*dlj          /(dr**4    ) &
                                                                               + (    cmt(k2,mm))*cmt(k1,ll)*8*dli*dlj*(diz-djz)*(djz-diz)/(dr**6 )

                 ccpotent(ll+5*NN+6,mm+4*NN+5) = ccpotent(ll+5*NN+6,mm+4*NN+5) + ( ll*cmt(k1,ll))*(-mm*smt(k2,mm))*tiz*tjz/(dr**2*dli*dlj)   &      ! cc/zsi/zcj
                                                                               - ( ll*cmt(k1,ll))*cmt(k2,mm)*2*dlj*tiz*(djz-diz)/(dr**4*dli) &
                                                                               - (-mm*smt(k2,mm))*smt(k1,ll)*2*dli*tjz*(diz-djz)/(dr**4*dlj) &
                                                                               + (    cmt(k2,mm))*smt(k1,ll)*2*dli*dlj          /(dr**4    ) &
                                                                               + (    cmt(k2,mm))*smt(k1,ll)*8*dli*dlj*(diz-djz)*(djz-diz)/(dr**6 )

                 !--------------------------------------------cc / zsj ----------------------------------------------------------------------------------------

                 ccpotent(ll     +1,mm+5*NN+6) = ccpotent(ll     +1,mm+5*NN+6) + (-ll*smt(k1,ll))*( mm*cmt(k2,mm))*tix*tjz/(dr**2*dli*dlj)   &      ! cc/xci/zsj
                                                                               - (-ll*smt(k1,ll))*smt(k2,mm)*2*dlj*tix*(djz-diz)/(dr**4*dli) &
                                                                               - ( mm*cmt(k2,mm))*cmt(k1,ll)*2*dli*tjz*(dix-djx)/(dr**4*dlj) &
                                                                              ! + (    smt(k2,mm))*cmt(k1,ll)*2*dli*dlj          /(dr**4    ) &
                                                                               + (    smt(k2,mm))*cmt(k1,ll)*8*dli*dlj*(dix-djx)*(djz-diz)/(dr**6 )

                 ccpotent(ll+  NN+2,mm+5*NN+6) = ccpotent(ll+  NN+2,mm+5*NN+6) + ( ll*cmt(k1,ll))*( mm*cmt(k2,mm))*tix*tjz/(dr**2*dli*dlj)   &      ! cc/xsi/zsj
                                                                               - ( ll*cmt(k1,ll))*smt(k2,mm)*2*dlj*tix*(djz-diz)/(dr**4*dli) &
                                                                               - ( mm*cmt(k2,mm))*smt(k1,ll)*2*dli*tjz*(dix-djx)/(dr**4*dlj) &
                                                                              ! + (    smt(k2,mm))*smt(k1,ll)*2*dli*dlj          /(dr**4    ) &
                                                                               + (    smt(k2,mm))*smt(k1,ll)*8*dli*dlj*(dix-djx)*(djz-diz)/(dr**6 )

                 ccpotent(ll+2*NN+3,mm+5*NN+6) = ccpotent(ll+2*NN+3,mm+5*NN+6) + (-ll*smt(k1,ll))*( mm*cmt(k2,mm))*tiy*tjz/(dr**2*dli*dlj)   &      ! cc/yci/zsj
                                                                               - (-ll*smt(k1,ll))*smt(k2,mm)*2*dlj*tiy*(djz-diz)/(dr**4*dli) &
                                                                               - ( mm*cmt(k2,mm))*cmt(k1,ll)*2*dli*tjz*(diy-djy)/(dr**4*dlj) &
                                                                              ! + (    smt(k2,mm))*cmt(k1,ll)*2*dli*dlj          /(dr**4    ) &
                                                                               + (    smt(k2,mm))*cmt(k1,ll)*8*dli*dlj*(diy-djy)*(djz-diz)/(dr**6 )

                 ccpotent(ll+3*NN+4,mm+5*NN+6) = ccpotent(ll+3*NN+4,mm+5*NN+6) + ( ll*cmt(k1,ll))*( mm*cmt(k2,mm))*tiy*tjz/(dr**2*dli*dlj)   &      ! cc/ysi/zsj
                                                                               - ( ll*cmt(k1,ll))*smt(k2,mm)*2*dlj*tiy*(djz-diz)/(dr**4*dli) &
                                                                               - ( mm*cmt(k2,mm))*smt(k1,ll)*2*dli*tjz*(diy-djy)/(dr**4*dlj) &
                                                                              ! + (    smt(k2,mm))*smt(k1,ll)*2*dli*dlj          /(dr**4    ) &
                                                                               + (    smt(k2,mm))*smt(k1,ll)*8*dli*dlj*(diy-djy)*(djz-diz)/(dr**6 )

                 ccpotent(ll+4*NN+5,mm+5*NN+6) = ccpotent(ll+4*NN+5,mm+5*NN+6) + (-ll*smt(k1,ll))*( mm*cmt(k2,mm))*tiz*tjz/(dr**2*dli*dlj)   &      ! cc/zci/zsj
                                                                               - (-ll*smt(k1,ll))*smt(k2,mm)*2*dlj*tiz*(djz-diz)/(dr**4*dli) &
                                                                               - ( mm*cmt(k2,mm))*cmt(k1,ll)*2*dli*tjz*(diz-djz)/(dr**4*dlj) &
                                                                               + (    smt(k2,mm))*cmt(k1,ll)*2*dli*dlj          /(dr**4    ) &
                                                                               + (    smt(k2,mm))*cmt(k1,ll)*8*dli*dlj*(diz-djz)*(djz-diz)/(dr**6 )

                 ccpotent(ll+5*NN+6,mm+5*NN+6) = ccpotent(ll+5*NN+6,mm+5*NN+6) + ( ll*cmt(k1,ll))*( mm*cmt(k2,mm))*tiz*tjz/(dr**2*dli*dlj)   &      ! cc/zsi/zsj
                                                                               - ( ll*cmt(k1,ll))*smt(k2,mm)*2*dlj*tiz*(djz-diz)/(dr**4*dli) &
                                                                               - ( mm*cmt(k2,mm))*smt(k1,ll)*2*dli*tjz*(diz-djz)/(dr**4*dlj) &
                                                                               + (    smt(k2,mm))*smt(k1,ll)*2*dli*dlj          /(dr**4    ) &
                                                                               + (    smt(k2,mm))*smt(k1,ll)*8*dli*dlj*(diz-djz)*(djz-diz)/(dr**6 )
              enddo
           enddo
           
        elseif ( mcoil .eq. icoil ) then

           do mm = 0, NN
              do ll = 0, NN

                 !--------------------------------------------cc / xci ----------------------------------------------------------------------------------------

                 ccpotent(ll     +1,mm     +1) = ccpotent(ll     +1,mm     +1) + (-ll*smt(k1,ll))*(-mm*smt(k1,mm))  *dlj                        /(dr**2*dli   )  &     ! cc/xci/xci
                                                                               - (-ll*smt(k1,ll))*(    cmt(k1,mm))*2*dlj*tix*(dix-djx)          /(dr**4*dli   )  &
                                                                               - (-ll*smt(k1,ll))*(-mm*smt(k1,mm))  *dlj*tix*tix                /(dr**2*dli**3)  &
                                                                               - (    cmt(k1,ll))*(-mm*smt(k1,mm))*2*dlj*tix*(dix-djx)          /(dr**4*dli   )  &
                                                                               - (    cmt(k1,ll))*(    cmt(k1,mm))*2*dli*dlj                    /(dr**4       )  &
                                                                               + (    cmt(k1,ll))*(    cmt(k1,mm))*8*dli*dlj*(dix-djx)*(dix-djx)/(dr**6       )

                 ccpotent(ll+  NN+2,mm     +1) = ccpotent(ll+  NN+2,mm     +1) + ( ll*cmt(k1,ll))*(-mm*smt(k1,mm))  *dlj                        /(dr**2*dli   )  &     ! cc/xsi/xci
                                                                               - ( ll*cmt(k1,ll))*(    cmt(k1,mm))*2*dlj*tix*(dix-djx)          /(dr**4*dli   )  &
                                                                               - ( ll*cmt(k1,ll))*(-mm*smt(k1,mm))  *dlj*tix*tix                /(dr**2*dli**3)  &
                                                                               - (    smt(k1,ll))*(-mm*smt(k1,mm))*2*dlj*tix*(dix-djx)          /(dr**4*dli   )  &
                                                                               - (    smt(k1,ll))*(    cmt(k1,mm))*2*dli*dlj                    /(dr**4       )  &
                                                                               + (    smt(k1,ll))*(    cmt(k1,mm))*8*dli*dlj*(dix-djx)*(dix-djx)/(dr**6       )

                 ccpotent(ll+2*NN+3,mm     +1) = ccpotent(ll+2*NN+3,mm     +1) + (-ll*smt(k1,ll))*(      0       )  *dlj                        /(dr**2*dli   )  &     ! cc/yci/xci
                                                                               - (-ll*smt(k1,ll))*(    cmt(k1,mm))*2*dlj*tiy*(dix-djx)          /(dr**4*dli   )  &
                                                                               - (-ll*smt(k1,ll))*(-mm*smt(k1,mm))  *dlj*tiy*tix                /(dr**2*dli**3)  &
                                                                               - (    cmt(k1,ll))*(-mm*smt(k1,mm))*2*dlj*tix*(diy-djy)          /(dr**4*dli   )  &
                                                                               - (    cmt(k1,ll))*(      0       )*2*dli*dlj                    /(dr**4       )  &
                                                                               + (    cmt(k1,ll))*(    cmt(k1,mm))*8*dli*dlj*(diy-djy)*(dix-djx)/(dr**6       )

                 ccpotent(ll+3*NN+4,mm     +1) = ccpotent(ll+3*NN+4,mm     +1) + ( ll*cmt(k1,ll))*(      0       )  *dlj                        /(dr**2*dli   )  &     ! cc/ysi/xci
                                                                               - ( ll*cmt(k1,ll))*(    cmt(k1,mm))*2*dlj*tiy*(dix-djx)          /(dr**4*dli   )  &
                                                                               - ( ll*cmt(k1,ll))*(-mm*smt(k1,mm))  *dlj*tiy*tix                /(dr**2*dli**3)  &
                                                                               - (    smt(k1,ll))*(-mm*smt(k1,mm))*2*dlj*tix*(diy-djy)          /(dr**4*dli   )  &
                                                                               - (    smt(k1,ll))*(      0       )*2*dli*dlj                    /(dr**4       )  &
                                                                               + (    smt(k1,ll))*(    cmt(k1,mm))*8*dli*dlj*(diy-djy)*(dix-djx)/(dr**6       ) 

                 ccpotent(ll+4*NN+5,mm     +1) = ccpotent(ll+4*NN+5,mm     +1) + (-ll*smt(k1,ll))*(      0       )  *dlj                        /(dr**2*dli   )  &     ! cc/zci/xci
                                                                               - (-ll*smt(k1,ll))*(    cmt(k1,mm))*2*dlj*tiz*(dix-djx)          /(dr**4*dli   )  &
                                                                               - (-ll*smt(k1,ll))*(-mm*smt(k1,mm))  *dlj*tiz*tix                /(dr**2*dli**3)  &
                                                                               - (    cmt(k1,ll))*(-mm*smt(k1,mm))*2*dlj*tix*(diz-djz)          /(dr**4*dli   )  &
                                                                               - (    cmt(k1,ll))*(      0       )*2*dli*dlj                    /(dr**4       )  &
                                                                               + (    cmt(k1,ll))*(    cmt(k1,mm))*8*dli*dlj*(diz-djz)*(dix-djx)/(dr**6       )

                 ccpotent(ll+5*NN+6,mm     +1) = ccpotent(ll+5*NN+6,mm     +1) + ( ll*cmt(k1,ll))*(      0       )  *dlj                        /(dr**2*dli   )  &     ! cc/zsi/xci
                                                                               - ( ll*cmt(k1,ll))*(    cmt(k1,mm))*2*dlj*tiz*(dix-djx)          /(dr**4*dli   )  &
                                                                               - ( ll*cmt(k1,ll))*(-mm*smt(k1,mm))  *dlj*tiz*tix                /(dr**2*dli**3)  &
                                                                               - (    smt(k1,ll))*(-mm*smt(k1,mm))*2*dlj*tix*(diz-djz)          /(dr**4*dli   )  &
                                                                               - (    smt(k1,ll))*(      0       )*2*dli*dlj                    /(dr**4       )  &
                                                                               + (    smt(k1,ll))*(    cmt(k1,mm))*8*dli*dlj*(diz-djz)*(dix-djx)/(dr**6       ) 

                 !--------------------------------------------cc / xsi ----------------------------------------------------------------------------------------

                 ccpotent(ll     +1,mm+  NN+2) = ccpotent(ll     +1,mm+  NN+2) + (-ll*smt(k1,ll))*( mm*cmt(k1,mm))  *dlj                        /(dr**2*dli   )  &     ! cc/xci/xsi
                                                                               - (-ll*smt(k1,ll))*(    smt(k1,mm))*2*dlj*tix*(dix-djx)          /(dr**4*dli   )  &
                                                                               - (-ll*smt(k1,ll))*( mm*cmt(k1,mm))  *dlj*tix*tix                /(dr**2*dli**3)  &
                                                                               - (    cmt(k1,ll))*( mm*cmt(k1,mm))*2*dlj*tix*(dix-djx)          /(dr**4*dli   )  &
                                                                               - (    cmt(k1,ll))*(    smt(k1,mm))*2*dli*dlj                    /(dr**4       )  &
                                                                               + (    cmt(k1,ll))*(    smt(k1,mm))*8*dli*dlj*(dix-djx)*(dix-djx)/(dr**6       )

                 ccpotent(ll+  NN+2,mm+  NN+2) = ccpotent(ll+  NN+2,mm+  NN+2) + ( ll*cmt(k1,ll))*( mm*cmt(k1,mm))  *dlj                        /(dr**2*dli   )  &     ! cc/xsi/xsi
                                                                               - ( ll*cmt(k1,ll))*(    smt(k1,mm))*2*dlj*tix*(dix-djx)          /(dr**4*dli   )  &
                                                                               - ( ll*cmt(k1,ll))*( mm*cmt(k1,mm))  *dlj*tix*tix                /(dr**2*dli**3)  &
                                                                               - (    smt(k1,ll))*( mm*cmt(k1,mm))*2*dlj*tix*(dix-djx)          /(dr**4*dli   )  &
                                                                               - (    smt(k1,ll))*(    smt(k1,mm))*2*dli*dlj                    /(dr**4       )  &
                                                                               + (    smt(k1,ll))*(    smt(k1,mm))*8*dli*dlj*(dix-djx)*(dix-djx)/(dr**6       )

                 ccpotent(ll+2*NN+3,mm+  NN+2) = ccpotent(ll+2*NN+3,mm+  NN+2) + (-ll*smt(k1,ll))*(      0       )  *dlj                        /(dr**2*dli   )  &     ! cc/yci/xsi
                                                                               - (-ll*smt(k1,ll))*(    cmt(k1,mm))*2*dlj*tiy*(dix-djx)          /(dr**4*dli   )  &
                                                                               - (-ll*smt(k1,ll))*(-mm*smt(k1,mm))  *dlj*tiy*tix                /(dr**2*dli**3)  &
                                                                               - (    cmt(k1,ll))*(-mm*smt(k1,mm))*2*dlj*tix*(diy-djy)          /(dr**4*dli   )  &
                                                                               - (    cmt(k1,ll))*(      0       )*2*dli*dlj                    /(dr**4       )  &
                                                                               + (    cmt(k1,ll))*(    cmt(k1,mm))*8*dli*dlj*(diy-djy)*(dix-djx)/(dr**6       )

                 ccpotent(ll+3*NN+4,mm+  NN+2) = ccpotent(ll+3*NN+4,mm+  NN+2) + ( ll*cmt(k1,ll))*(      0       )  *dlj                        /(dr**2*dli   )  &     ! cc/ysi/xsi
                                                                               - ( ll*cmt(k1,ll))*(    cmt(k1,mm))*2*dlj*tiy*(dix-djx)          /(dr**4*dli   )  &
                                                                               - ( ll*cmt(k1,ll))*(-mm*smt(k1,mm))  *dlj*tiy*tix                /(dr**2*dli**3)  &
                                                                               - (    smt(k1,ll))*(-mm*smt(k1,mm))*2*dlj*tix*(diy-djy)          /(dr**4*dli   )  &
                                                                               - (    smt(k1,ll))*(      0       )*2*dli*dlj                    /(dr**4       )  &
                                                                               + (    smt(k1,ll))*(    cmt(k1,mm))*8*dli*dlj*(diy-djy)*(dix-djx)/(dr**6       ) 

                 ccpotent(ll+4*NN+5,mm+  NN+2) = ccpotent(ll+4*NN+5,mm+  NN+2) + (-ll*smt(k1,ll))*(      0       )  *dlj                        /(dr**2*dli   )  &     ! cc/zci/xsi
                                                                               - (-ll*smt(k1,ll))*(    cmt(k1,mm))*2*dlj*tiz*(dix-djx)          /(dr**4*dli   )  &
                                                                               - (-ll*smt(k1,ll))*(-mm*smt(k1,mm))  *dlj*tiz*tix                /(dr**2*dli**3)  &
                                                                               - (    cmt(k1,ll))*(-mm*smt(k1,mm))*2*dlj*tix*(diz-djz)          /(dr**4*dli   )  &
                                                                               - (    cmt(k1,ll))*(      0       )*2*dli*dlj                    /(dr**4       )  &
                                                                               + (    cmt(k1,ll))*(    cmt(k1,mm))*8*dli*dlj*(diz-djz)*(dix-djx)/(dr**6       )

                 ccpotent(ll+5*NN+6,mm+  NN+2) = ccpotent(ll+5*NN+6,mm+  NN+2) + ( ll*cmt(k1,ll))*(      0       )  *dlj                        /(dr**2*dli   )  &     ! cc/zsi/xsi
                                                                               - ( ll*cmt(k1,ll))*(    cmt(k1,mm))*2*dlj*tiz*(dix-djx)          /(dr**4*dli   )  &
                                                                               - ( ll*cmt(k1,ll))*(-mm*smt(k1,mm))  *dlj*tiz*tix                /(dr**2*dli**3)  &
                                                                               - (    smt(k1,ll))*(-mm*smt(k1,mm))*2*dlj*tix*(diz-djz)          /(dr**4*dli   )  &
                                                                               - (    smt(k1,ll))*(      0       )*2*dli*dlj                    /(dr**4       )  &
                                                                               + (    smt(k1,ll))*(    cmt(k1,mm))*8*dli*dlj*(diz-djz)*(dix-djx)/(dr**6       ) 

                 !--------------------------------------------cc / yci ----------------------------------------------------------------------------------------

                 ccpotent(ll     +1,mm+2*NN+3) = ccpotent(ll     +1,mm+2*NN+3) + (-ll*smt(k1,ll))*(      0       )  *dlj                        /(dr**2*dli   )  &     ! cc/xci/yci
                                                                               - (-ll*smt(k1,ll))*(    cmt(k1,mm))*2*dlj*tix*(diy-djy)          /(dr**4*dli   )  &
                                                                               - (-ll*smt(k1,ll))*(-mm*smt(k1,mm))  *dlj*tix*tiy                /(dr**2*dli**3)  &
                                                                               - (    cmt(k1,ll))*(-mm*smt(k1,mm))*2*dlj*tiy*(dix-djx)          /(dr**4*dli   )  &
                                                                               - (    cmt(k1,ll))*(      0       )*2*dli*dlj                    /(dr**4       )  &
                                                                               + (    cmt(k1,ll))*(    cmt(k1,mm))*8*dli*dlj*(dix-djx)*(diy-djy)/(dr**6       )

                 ccpotent(ll+  NN+2,mm+2*NN+3) = ccpotent(ll+  NN+2,mm+2*NN+3) + ( ll*cmt(k1,ll))*(      0       )  *dlj                        /(dr**2*dli   )  &     ! cc/xsi/yci
                                                                               - ( ll*cmt(k1,ll))*(    cmt(k1,mm))*2*dlj*tix*(diy-djy)          /(dr**4*dli   )  &
                                                                               - ( ll*cmt(k1,ll))*(-mm*smt(k1,mm))  *dlj*tix*tiy                /(dr**2*dli**3)  &
                                                                               - (    smt(k1,ll))*(-mm*smt(k1,mm))*2*dlj*tiy*(dix-djx)          /(dr**4*dli   )  &
                                                                               - (    smt(k1,ll))*(      0       )*2*dli*dlj                    /(dr**4       )  &
                                                                               + (    smt(k1,ll))*(    cmt(k1,mm))*8*dli*dlj*(dix-djx)*(diy-djy)/(dr**6       )

                 ccpotent(ll+2*NN+3,mm+2*NN+3) = ccpotent(ll+2*NN+3,mm+2*NN+3) + (-ll*smt(k1,ll))*(-mm*smt(k1,mm))  *dlj                        /(dr**2*dli   )  &     ! cc/yci/yci
                                                                               - (-ll*smt(k1,ll))*(    cmt(k1,mm))*2*dlj*tiy*(diy-djy)          /(dr**4*dli   )  &
                                                                               - (-ll*smt(k1,ll))*(-mm*smt(k1,mm))  *dlj*tiy*tiy                /(dr**2*dli**3)  &
                                                                               - (    cmt(k1,ll))*(-mm*smt(k1,mm))*2*dlj*tiy*(diy-djy)          /(dr**4*dli   )  &
                                                                               - (    cmt(k1,ll))*(    cmt(k1,mm))*2*dli*dlj                    /(dr**4       )  &
                                                                               + (    cmt(k1,ll))*(    cmt(k1,mm))*8*dli*dlj*(diy-djy)*(diy-djy)/(dr**6       )

                 ccpotent(ll+3*NN+4,mm+2*NN+3) = ccpotent(ll+3*NN+4,mm+2*NN+3) + ( ll*cmt(k1,ll))*(-mm*smt(k1,mm))  *dlj                        /(dr**2*dli   )  &     ! cc/ysi/yci
                                                                               - ( ll*cmt(k1,ll))*(    cmt(k1,mm))*2*dlj*tiy*(diy-djy)          /(dr**4*dli   )  &
                                                                               - ( ll*cmt(k1,ll))*(-mm*smt(k1,mm))  *dlj*tiy*tiy                /(dr**2*dli**3)  &
                                                                               - (    smt(k1,ll))*(-mm*smt(k1,mm))*2*dlj*tiy*(diy-djy)          /(dr**4*dli   )  &
                                                                               - (    smt(k1,ll))*(    cmt(k1,mm))*2*dli*dlj                    /(dr**4       )  &
                                                                               + (    smt(k1,ll))*(    cmt(k1,mm))*8*dli*dlj*(diy-djy)*(diy-djy)/(dr**6       ) 

                 ccpotent(ll+4*NN+5,mm+2*NN+3) = ccpotent(ll+4*NN+5,mm+2*NN+3) + (-ll*smt(k1,ll))*(      0       )  *dlj                        /(dr**2*dli   )  &     ! cc/zci/yci
                                                                               - (-ll*smt(k1,ll))*(    cmt(k1,mm))*2*dlj*tiz*(diy-djy)          /(dr**4*dli   )  &
                                                                               - (-ll*smt(k1,ll))*(-mm*smt(k1,mm))  *dlj*tiz*tiy                /(dr**2*dli**3)  &
                                                                               - (    cmt(k1,ll))*(-mm*smt(k1,mm))*2*dlj*tiy*(diz-djz)          /(dr**4*dli   )  &
                                                                               - (    cmt(k1,ll))*(      0       )*2*dli*dlj                    /(dr**4       )  &
                                                                               + (    cmt(k1,ll))*(    cmt(k1,mm))*8*dli*dlj*(diz-djz)*(diy-djy)/(dr**6       )

                 ccpotent(ll+5*NN+6,mm+2*NN+3) = ccpotent(ll+5*NN+6,mm+2*NN+3) + ( ll*cmt(k1,ll))*(      0       )  *dlj                        /(dr**2*dli   )  &     ! cc/zsi/yci
                                                                               - ( ll*cmt(k1,ll))*(    cmt(k1,mm))*2*dlj*tiz*(diy-djy)          /(dr**4*dli   )  &
                                                                               - ( ll*cmt(k1,ll))*(-mm*smt(k1,mm))  *dlj*tiz*tiy                /(dr**2*dli**3)  &
                                                                               - (    smt(k1,ll))*(-mm*smt(k1,mm))*2*dlj*tiy*(diz-djz)          /(dr**4*dli   )  &
                                                                               - (    smt(k1,ll))*(      0       )*2*dli*dlj                    /(dr**4       )  &
                                                                               + (    smt(k1,ll))*(    cmt(k1,mm))*8*dli*dlj*(diz-djz)*(diy-djy)/(dr**6       ) 

                 !--------------------------------------------cc / ysi ----------------------------------------------------------------------------------------

                 ccpotent(ll     +1,mm+3*NN+4) = ccpotent(ll     +1,mm+3*NN+4) + (-ll*smt(k1,ll))*(      0       )  *dlj                        /(dr**2*dli   )  &     ! cc/xci/ysi
                                                                               - (-ll*smt(k1,ll))*(    smt(k1,mm))*2*dlj*tix*(diy-djy)          /(dr**4*dli   )  &
                                                                               - (-ll*smt(k1,ll))*( mm*cmt(k1,mm))  *dlj*tix*tiy                /(dr**2*dli**3)  &
                                                                               - (    cmt(k1,ll))*( mm*cmt(k1,mm))*2*dlj*tiy*(dix-djx)          /(dr**4*dli   )  &
                                                                               - (    cmt(k1,ll))*(      0       )*2*dli*dlj                    /(dr**4       )  &
                                                                               + (    cmt(k1,ll))*(    smt(k1,mm))*8*dli*dlj*(dix-djx)*(diy-djy)/(dr**6       )

                 ccpotent(ll+  NN+2,mm+3*NN+4) = ccpotent(ll+  NN+2,mm+3*NN+4) + ( ll*cmt(k1,ll))*(      0       )  *dlj                        /(dr**2*dli   )  &     ! cc/xsi/ysi
                                                                               - ( ll*cmt(k1,ll))*(    smt(k1,mm))*2*dlj*tix*(diy-djy)          /(dr**4*dli   )  &
                                                                               - ( ll*cmt(k1,ll))*( mm*cmt(k1,mm))  *dlj*tix*tiy                /(dr**2*dli**3)  &
                                                                               - (    smt(k1,ll))*( mm*cmt(k1,mm))*2*dlj*tiy*(dix-djx)          /(dr**4*dli   )  &
                                                                               - (    smt(k1,ll))*(      0       )*2*dli*dlj                    /(dr**4       )  &
                                                                               + (    smt(k1,ll))*(    smt(k1,mm))*8*dli*dlj*(dix-djx)*(diy-djy)/(dr**6       )

                 ccpotent(ll+2*NN+3,mm+3*NN+4) = ccpotent(ll+2*NN+3,mm+3*NN+4) + (-ll*smt(k1,ll))*( mm*cmt(k1,mm))  *dlj                        /(dr**2*dli   )  &     ! cc/yci/ysi
                                                                               - (-ll*smt(k1,ll))*(    cmt(k1,mm))*2*dlj*tiy*(diy-djy)          /(dr**4*dli   )  &
                                                                               - (-ll*smt(k1,ll))*(-mm*smt(k1,mm))  *dlj*tiy*tiy                /(dr**2*dli**3)  &
                                                                               - (    cmt(k1,ll))*(-mm*smt(k1,mm))*2*dlj*tiy*(diy-djy)          /(dr**4*dli   )  &
                                                                               - (    cmt(k1,ll))*(    smt(k1,mm))*2*dli*dlj                    /(dr**4       )  &
                                                                               + (    cmt(k1,ll))*(    cmt(k1,mm))*8*dli*dlj*(diy-djy)*(diy-djy)/(dr**6       )

                 ccpotent(ll+3*NN+4,mm+3*NN+4) = ccpotent(ll+3*NN+4,mm+3*NN+4) + ( ll*cmt(k1,ll))*( mm*cmt(k1,mm))  *dlj                        /(dr**2*dli   )  &     ! cc/ysi/ysi
                                                                               - ( ll*cmt(k1,ll))*(    cmt(k1,mm))*2*dlj*tiy*(dix-djx)          /(dr**4*dli   )  &
                                                                               - ( ll*cmt(k1,ll))*(-mm*smt(k1,mm))  *dlj*tiy*tix                /(dr**2*dli**3)  &
                                                                               - (    smt(k1,ll))*(-mm*smt(k1,mm))*2*dlj*tix*(diy-djy)          /(dr**4*dli   )  &
                                                                               - (    smt(k1,ll))*(    smt(k1,mm))*2*dli*dlj                    /(dr**4       )  &
                                                                               + (    smt(k1,ll))*(    cmt(k1,mm))*8*dli*dlj*(diy-djy)*(dix-djx)/(dr**6       ) 

                 ccpotent(ll+4*NN+5,mm+3*NN+4) = ccpotent(ll+4*NN+5,mm+3*NN+4) + (-ll*smt(k1,ll))*(      0       )  *dlj                        /(dr**2*dli   )  &     ! cc/zci/ysi
                                                                               - (-ll*smt(k1,ll))*(    cmt(k1,mm))*2*dlj*tiz*(diy-djy)          /(dr**4*dli   )  &
                                                                               - (-ll*smt(k1,ll))*(-mm*smt(k1,mm))  *dlj*tiz*tiy                /(dr**2*dli**3)  &
                                                                               - (    cmt(k1,ll))*(-mm*smt(k1,mm))*2*dlj*tiy*(diz-djz)          /(dr**4*dli   )  &
                                                                               - (    cmt(k1,ll))*(      0       )*2*dli*dlj                    /(dr**4       )  &
                                                                               + (    cmt(k1,ll))*(    cmt(k1,mm))*8*dli*dlj*(diz-djz)*(diy-djy)/(dr**6       )

                 ccpotent(ll+5*NN+6,mm+3*NN+4) = ccpotent(ll+5*NN+6,mm+3*NN+4) + ( ll*cmt(k1,ll))*(      0       )  *dlj                        /(dr**2*dli   )  &     ! cc/zsi/ysi
                                                                               - ( ll*cmt(k1,ll))*(    cmt(k1,mm))*2*dlj*tiz*(diy-djy)          /(dr**4*dli   )  &
                                                                               - ( ll*cmt(k1,ll))*(-mm*smt(k1,mm))  *dlj*tiz*tiy                /(dr**2*dli**3)  &
                                                                               - (    smt(k1,ll))*(-mm*smt(k1,mm))*2*dlj*tiy*(diz-djz)          /(dr**4*dli   )  &
                                                                               - (    smt(k1,ll))*(      0       )*2*dli*dlj                    /(dr**4       )  &
                                                                               + (    smt(k1,ll))*(    cmt(k1,mm))*8*dli*dlj*(diz-djz)*(diy-djy)/(dr**6       ) 

                 !--------------------------------------------cc / zci ----------------------------------------------------------------------------------------

                 ccpotent(ll     +1,mm+4*NN+5) = ccpotent(ll     +1,mm+4*NN+5) + (-ll*smt(k1,ll))*(      0       )  *dlj                        /(dr**2*dli   )  &     ! cc/xci/zci
                                                                               - (-ll*smt(k1,ll))*(    cmt(k1,mm))*2*dlj*tix*(diz-djz)          /(dr**4*dli   )  &
                                                                               - (-ll*smt(k1,ll))*(-mm*smt(k1,mm))  *dlj*tix*tiz                /(dr**2*dli**3)  &
                                                                               - (    cmt(k1,ll))*(-mm*smt(k1,mm))*2*dlj*tiz*(dix-djx)          /(dr**4*dli   )  &
                                                                               - (    cmt(k1,ll))*(      0       )*2*dli*dlj                    /(dr**4       )  &
                                                                               + (    cmt(k1,ll))*(    cmt(k1,mm))*8*dli*dlj*(dix-djx)*(diz-djz)/(dr**6       )

                 ccpotent(ll+  NN+2,mm+4*NN+5) = ccpotent(ll+  NN+2,mm+4*NN+5) + ( ll*cmt(k1,ll))*(      0       )  *dlj                        /(dr**2*dli   )  &     ! cc/xsi/zci
                                                                               - ( ll*cmt(k1,ll))*(    cmt(k1,mm))*2*dlj*tix*(diz-djz)          /(dr**4*dli   )  &
                                                                               - ( ll*cmt(k1,ll))*(-mm*smt(k1,mm))  *dlj*tix*tiz                /(dr**2*dli**3)  &
                                                                               - (    smt(k1,ll))*(-mm*smt(k1,mm))*2*dlj*tiz*(dix-djx)          /(dr**4*dli   )  &
                                                                               - (    smt(k1,ll))*(      0       )*2*dli*dlj                    /(dr**4       )  &
                                                                               + (    smt(k1,ll))*(    cmt(k1,mm))*8*dli*dlj*(dix-djx)*(diz-djz)/(dr**6       )

                 ccpotent(ll+2*NN+3,mm+4*NN+5) = ccpotent(ll+2*NN+3,mm+4*NN+5) + (-ll*smt(k1,ll))*(      0       )  *dlj                        /(dr**2*dli   )  &     ! cc/yci/zci
                                                                               - (-ll*smt(k1,ll))*(    cmt(k1,mm))*2*dlj*tiy*(diz-djz)          /(dr**4*dli   )  &
                                                                               - (-ll*smt(k1,ll))*(-mm*smt(k1,mm))  *dlj*tiy*tiz                /(dr**2*dli**3)  &
                                                                               - (    cmt(k1,ll))*(-mm*smt(k1,mm))*2*dlj*tiz*(diy-djy)          /(dr**4*dli   )  &
                                                                               - (    cmt(k1,ll))*(      0       )*2*dli*dlj                    /(dr**4       )  &
                                                                               + (    cmt(k1,ll))*(    cmt(k1,mm))*8*dli*dlj*(diy-djy)*(diz-djz)/(dr**6       )

                 ccpotent(ll+3*NN+4,mm+4*NN+5) = ccpotent(ll+3*NN+4,mm+4*NN+5) + ( ll*cmt(k1,ll))*(      0       )  *dlj                        /(dr**2*dli   )  &     ! cc/ysi/zci
                                                                               - ( ll*cmt(k1,ll))*(    cmt(k1,mm))*2*dlj*tiy*(diz-djz)          /(dr**4*dli   )  &
                                                                               - ( ll*cmt(k1,ll))*(-mm*smt(k1,mm))  *dlj*tiy*tiz                /(dr**2*dli**3)  &
                                                                               - (    smt(k1,ll))*(-mm*smt(k1,mm))*2*dlj*tiz*(diy-djy)          /(dr**4*dli   )  &
                                                                               - (    smt(k1,ll))*(      0       )*2*dli*dlj                    /(dr**4       )  &
                                                                               + (    smt(k1,ll))*(    cmt(k1,mm))*8*dli*dlj*(diy-djy)*(diz-djz)/(dr**6       ) 

                 ccpotent(ll+4*NN+5,mm+4*NN+5) = ccpotent(ll+4*NN+5,mm+4*NN+5) + (-ll*smt(k1,ll))*(-mm*smt(k1,mm))  *dlj                        /(dr**2*dli   )  &     ! cc/zci/zci
                                                                               - (-ll*smt(k1,ll))*(    cmt(k1,mm))*2*dlj*tiz*(diz-djz)          /(dr**4*dli   )  &
                                                                               - (-ll*smt(k1,ll))*(-mm*smt(k1,mm))  *dlj*tiz*tiz                /(dr**2*dli**3)  &
                                                                               - (    cmt(k1,ll))*(-mm*smt(k1,mm))*2*dlj*tiz*(diz-djz)          /(dr**4*dli   )  &
                                                                               - (    cmt(k1,ll))*(    cmt(k1,mm))*2*dli*dlj                    /(dr**4       )  &
                                                                               + (    cmt(k1,ll))*(    cmt(k1,mm))*8*dli*dlj*(diz-djz)*(diz-djz)/(dr**6       )

                 ccpotent(ll+5*NN+6,mm+4*NN+5) = ccpotent(ll+5*NN+6,mm+4*NN+5) + ( ll*cmt(k1,ll))*(-mm*smt(k1,mm))  *dlj                        /(dr**2*dli   )  &     ! cc/zsi/zci
                                                                               - ( ll*cmt(k1,ll))*(    cmt(k1,mm))*2*dlj*tiz*(diz-djz)          /(dr**4*dli   )  &
                                                                               - ( ll*cmt(k1,ll))*(-mm*smt(k1,mm))  *dlj*tiz*tiz                /(dr**2*dli**3)  &
                                                                               - (    smt(k1,ll))*(-mm*smt(k1,mm))*2*dlj*tiz*(diz-djz)          /(dr**4*dli   )  &
                                                                               - (    smt(k1,ll))*(    cmt(k1,mm))*2*dli*dlj                    /(dr**4       )  &
                                                                               + (    smt(k1,ll))*(    cmt(k1,mm))*8*dli*dlj*(diz-djz)*(diz-djz)/(dr**6       ) 

                 !--------------------------------------------cc / zsi ----------------------------------------------------------------------------------------

                 ccpotent(ll     +1,mm+5*NN+6) = ccpotent(ll     +1,mm+5*NN+6) + (-ll*smt(k1,ll))*(      0       )  *dlj                        /(dr**2*dli   )  &     ! cc/xci/zsi
                                                                               - (-ll*smt(k1,ll))*(    smt(k1,mm))*2*dlj*tix*(diz-djz)          /(dr**4*dli   )  &
                                                                               - (-ll*smt(k1,ll))*( mm*cmt(k1,mm))  *dlj*tix*tiz                /(dr**2*dli**3)  &
                                                                               - (    cmt(k1,ll))*( mm*cmt(k1,mm))*2*dlj*tiz*(dix-djx)          /(dr**4*dli   )  &
                                                                               - (    cmt(k1,ll))*(      0       )*2*dli*dlj                    /(dr**4       )  &
                                                                               + (    cmt(k1,ll))*(    smt(k1,mm))*8*dli*dlj*(dix-djx)*(diz-djz)/(dr**6       )

                 ccpotent(ll+  NN+2,mm+5*NN+6) = ccpotent(ll+  NN+2,mm+5*NN+6) + ( ll*cmt(k1,ll))*(      0       )  *dlj                        /(dr**2*dli   )  &     ! cc/xsi/zsi
                                                                               - ( ll*cmt(k1,ll))*(    smt(k1,mm))*2*dlj*tix*(diz-djz)          /(dr**4*dli   )  &
                                                                               - ( ll*cmt(k1,ll))*( mm*cmt(k1,mm))  *dlj*tix*tiz                /(dr**2*dli**3)  &
                                                                               - (    smt(k1,ll))*( mm*cmt(k1,mm))*2*dlj*tiz*(dix-djx)          /(dr**4*dli   )  &
                                                                               - (    smt(k1,ll))*(      0       )*2*dli*dlj                    /(dr**4       )  &
                                                                               + (    smt(k1,ll))*(    smt(k1,mm))*8*dli*dlj*(dix-djx)*(diz-djz)/(dr**6       )

                 ccpotent(ll+2*NN+3,mm+5*NN+6) = ccpotent(ll+2*NN+3,mm+5*NN+6) + (-ll*smt(k1,ll))*(      0       )  *dlj                        /(dr**2*dli   )  &     ! cc/yci/zsi
                                                                               - (-ll*smt(k1,ll))*(    cmt(k1,mm))*2*dlj*tiy*(diz-djz)          /(dr**4*dli   )  &
                                                                               - (-ll*smt(k1,ll))*(-mm*smt(k1,mm))  *dlj*tiy*tiz                /(dr**2*dli**3)  &
                                                                               - (    cmt(k1,ll))*(-mm*smt(k1,mm))*2*dlj*tiz*(diy-djy)          /(dr**4*dli   )  &
                                                                               - (    cmt(k1,ll))*(      0       )*2*dli*dlj                    /(dr**4       )  &
                                                                               + (    cmt(k1,ll))*(    cmt(k1,mm))*8*dli*dlj*(diy-djy)*(diz-djz)/(dr**6       )

                 ccpotent(ll+3*NN+4,mm+5*NN+6) = ccpotent(ll+3*NN+4,mm+5*NN+6) + ( ll*cmt(k1,ll))*(      0       )  *dlj                        /(dr**2*dli   )  &     ! cc/ysi/zsi
                                                                               - ( ll*cmt(k1,ll))*(    cmt(k1,mm))*2*dlj*tiy*(diz-djz)          /(dr**4*dli   )  &
                                                                               - ( ll*cmt(k1,ll))*(-mm*smt(k1,mm))  *dlj*tiy*tiz                /(dr**2*dli**3)  &
                                                                               - (    smt(k1,ll))*(-mm*smt(k1,mm))*2*dlj*tiz*(diy-djy)          /(dr**4*dli   )  &
                                                                               - (    smt(k1,ll))*(      0       )*2*dli*dlj                    /(dr**4       )  &
                                                                               + (    smt(k1,ll))*(    cmt(k1,mm))*8*dli*dlj*(diy-djy)*(diz-djz)/(dr**6       ) 

                 ccpotent(ll+4*NN+5,mm+5*NN+6) = ccpotent(ll+4*NN+5,mm+5*NN+6) + (-ll*smt(k1,ll))*( mm*cmt(k1,mm))  *dlj                        /(dr**2*dli   )  &     ! cc/zci/zsi
                                                                               - (-ll*smt(k1,ll))*(    cmt(k1,mm))*2*dlj*tiz*(diz-djz)          /(dr**4*dli   )  &
                                                                               - (-ll*smt(k1,ll))*(-mm*smt(k1,mm))  *dlj*tiz*tiz                /(dr**2*dli**3)  &
                                                                               - (    cmt(k1,ll))*(-mm*smt(k1,mm))*2*dlj*tiz*(diz-djz)          /(dr**4*dli   )  &
                                                                               - (    cmt(k1,ll))*(    smt(k1,mm))*2*dli*dlj                    /(dr**4       )  &
                                                                               + (    cmt(k1,ll))*(    cmt(k1,mm))*8*dli*dlj*(diz-djz)*(diz-djz)/(dr**6       )

                 ccpotent(ll+5*NN+6,mm+5*NN+6) = ccpotent(ll+5*NN+6,mm+5*NN+6) + ( ll*cmt(k1,ll))*( mm*cmt(k1,mm))  *dlj                        /(dr**2*dli   )  &     ! cc/zsi/zsi
                                                                               - ( ll*cmt(k1,ll))*(    cmt(k1,mm))*2*dlj*tiz*(diz-djz)          /(dr**4*dli   )  &
                                                                               - ( ll*cmt(k1,ll))*(-mm*smt(k1,mm))  *dlj*tiz*tiz                /(dr**2*dli**3)  &
                                                                               - (    smt(k1,ll))*(-mm*smt(k1,mm))*2*dlj*tiz*(diz-djz)          /(dr**4*dli   )  &
                                                                               - (    smt(k1,ll))*(    smt(k1,mm))*2*dli*dlj                    /(dr**4       )  &
                                                                               + (    smt(k1,ll))*(    cmt(k1,mm))*8*dli*dlj*(diz-djz)*(diz-djz)/(dr**6       ) 

              enddo
           enddo
        else
           if(myid .eq. 0) write(ounit,'("epotent : "10X" : mcoil .ne. icoil or jcoil ; mcoil = "i4" ; icoil = "i4" ; jcoil = "i4)') mcoil, icoil, jcoil
           call MPI_ABORT( MPI_COMM_WORLD, 1, ierr )
        endif

     enddo
  enddo
                
  ccpotent = ccpotent * pi2/NDcoil * pi2/NDcoil

  return
end subroutine epotent2

!!$  do k1 = 0, NDcoil-1
!!$     
!!$     dli(0) = sqrt(coil(icoil)%xt(k1)**2 + coil(icoil)%yt(k1)**2 + coil(icoil)%zt(k1)**2)
!!$
!!$     do ll = 0, NN
!!$        dli(ll        + 1) = -ll*smt(k1, ll)*coil(icoil)%xt(k1)/dli(0)
!!$        dli(ll +   NN + 2) =  ll*cmt(k1, ll)*coil(icoil)%xt(k1)/dli(0)
!!$        dli(ll + 2*NN + 3) = -ll*smt(k1, ll)*coil(icoil)%yt(k1)/dli(0)
!!$        dli(ll + 3*NN + 4) =  ll*cmt(k1, ll)*coil(icoil)%yt(k1)/dli(0)
!!$        dli(ll + 4*NN + 5) = -ll*smt(k1, ll)*coil(icoil)%zt(k1)/dli(0)
!!$        dli(ll + 5*NN + 6) =  ll*cmt(k1, ll)*coil(icoil)%zt(k1)/dli(0)
!!$     enddo ! end ll
!!$
!!$     do k2 = 0, NDcoil-1
!!$
!!$        dlj(0) = sqrt(coil(jcoil)%xt(k2)**2 + coil(jcoil)%yt(k2)**2 + coil(jcoil)%zt(k2)**2)
!!$        dr  = sqrt( (coil(icoil)%xx(k1)-coil(jcoil)%xx(k2))**2 + (coil(icoil)%yy(k1)-coil(jcoil)%yy(k2))**2 + (coil(icoil)%zz(k1)-coil(jcoil)%zz(k2))**2 ) + delta
!!$        mincc(icoil,jcoil) = min(dr, mincc(icoil,jcoil))
!!$        
!!$        do ll = 0, NN
!!$           dlj(ll        + 1) = -ll*smt(k2, ll)*coil(jcoil)%xt(k2)/dlj(0)
!!$           dlj(ll +   NN + 2) =  ll*cmt(k2, ll)*coil(jcoil)%xt(k2)/dlj(0)
!!$           dlj(ll + 2*NN + 3) = -ll*smt(k2, ll)*coil(jcoil)%yt(k2)/dlj(0)
!!$           dlj(ll + 3*NN + 4) =  ll*cmt(k2, ll)*coil(jcoil)%yt(k2)/dlj(0)
!!$           dlj(ll + 4*NN + 5) = -ll*smt(k2, ll)*coil(jcoil)%zt(k2)/dlj(0)
!!$           dlj(ll + 5*NN + 6) =  ll*cmt(k2, ll)*coil(jcoil)%zt(k2)/dlj(0)
!!$        enddo ! end ll
!!$
!!$        do ll = 0, NN
!!$           dri(ll        + 1,0) = cmt(k1,ll)*(coil(icoil)%xx(k1)-coil(jcoil)%xx(k2))/dr; drj(ll        + 1,0) = cmt(k2,ll)*(coil(jcoil)%xx(k2)-coil(icoil)%xx(k1))/dr
!!$           dri(ll +   NN + 2,0) = smt(k1,ll)*(coil(icoil)%xx(k1)-coil(jcoil)%xx(k2))/dr; drj(ll +   NN + 2,0) = smt(k2,ll)*(coil(jcoil)%xx(k2)-coil(icoil)%xx(k1))/dr
!!$           dri(ll + 2*NN + 3,0) = cmt(k1,ll)*(coil(icoil)%yy(k1)-coil(jcoil)%yy(k2))/dr; drj(ll + 2*NN + 3,0) = cmt(k2,ll)*(coil(jcoil)%yy(k2)-coil(icoil)%yy(k1))/dr
!!$           dri(ll + 3*NN + 4,0) = smt(k1,ll)*(coil(icoil)%yy(k1)-coil(jcoil)%yy(k2))/dr; drj(ll + 3*NN + 4,0) = smt(k2,ll)*(coil(jcoil)%yy(k2)-coil(icoil)%yy(k1))/dr
!!$           dri(ll + 4*NN + 5,0) = cmt(k1,ll)*(coil(icoil)%zz(k1)-coil(jcoil)%zz(k2))/dr; drj(ll + 4*NN + 5,0) = cmt(k2,ll)*(coil(jcoil)%zz(k2)-coil(icoil)%zz(k1))/dr
!!$           dri(ll + 5*NN + 6,0) = smt(k1,ll)*(coil(icoil)%zz(k1)-coil(jcoil)%zz(k2))/dr; drj(ll + 5*NN + 6,0) = smt(k2,ll)*(coil(jcoil)%zz(k2)-coil(icoil)%zz(k1))/dr
!!$        enddo ! end ll
!!$
!!$
!!$        ccpotent(0,0) = ccpotent(0,0) + dli(0)*dlj(0)/dr**2 ! if takes a flexible exponent, how to normalize?
!!$        do ll = 1, Cdof
!!$           ccpotent(ll,0) = ccpotent(ll,0) + dlj(0)*dli(ll)/dr**2 - 2*dli(0)*dlj(0)*dri(ll)/dr**3
!!$        enddo
!!$       
!!$     enddo
!!$  enddo
!!$
!!$
!!$
!!$  do k1 = 0, NDcoil-1
!!$     
!!$     dli(0) = sqrt(coil(icoil)%xt(k1)**2 + coil(icoil)%yt(k1)**2 + coil(icoil)%zt(k1)**2)
!!$     dr  = zero
!!$     
!!$     tix(0) = coil(icoil)%xt(k1); dix(0) = coil(icoil)%xx(k1)
!!$     tiy(0) = coil(icoil)%yt(k1); diy(0) = coil(icoil)%yy(k1)
!!$     tiz(0) = coil(icoil)%zt(k1); diz(0) = coil(icoil)%zz(k1)
!!$
!!$     do ll = 0, NN
!!$        tix(ll        + 1) = -ll*smt(k1, ll); dix(ll        + 1) = cmt(k1, ll)
!!$        tix(ll +   NN + 2) =  ll*cmt(k1, ll); dix(ll +   NN + 2) = smt(k1, ll)
!!$        tiy(ll + 2*NN + 3) = -ll*smt(k1, ll); diy(ll + 2*NN + 3) = cmt(k1, ll)
!!$        tiy(ll + 3*NN + 4) =  ll*cmt(k1, ll); diy(ll + 3*NN + 4) = smt(k1, ll)
!!$        tiz(ll + 4*NN + 5) = -ll*smt(k1, ll); diz(ll + 4*NN + 5) = cmt(k1, ll)
!!$        tiz(ll + 5*NN + 6) =  ll*cmt(k1, ll); diz(ll + 5*NN + 6) = smt(k1, ll)
!!$     enddo ! end ll
!!$
!!$     do k2 = 0, NDcoil-1
!!$
!!$        dlj(0) = sqrt(coil(jcoil)%xt(k2)**2 + coil(jcoil)%yt(k2)**2 + coil(jcoil)%zt(k2)**2)
!!$
!!$        tjx(0) = coil(jcoil)%xt(k2); djx(0) = coil(jcoil)%xx(k2)
!!$        tjy(0) = coil(jcoil)%yt(k2); djy(0) = coil(jcoil)%yy(k2)
!!$        tjz(0) = coil(jcoil)%zt(k2); djz(0) = coil(jcoil)%zz(k2)
!!$
!!$        do ll = 0, NN
!!$           tjx(ll        + 1) = -ll*smt(k1, ll); djx(ll        + 1) = cmt(k1, ll)
!!$           tjx(ll +   NN + 2) =  ll*cmt(k1, ll); djx(ll +   NN + 2) = smt(k1, ll)
!!$           tjy(ll + 2*NN + 3) = -ll*smt(k1, ll); djy(ll + 2*NN + 3) = cmt(k1, ll)
!!$           tjy(ll + 3*NN + 4) =  ll*cmt(k1, ll); djy(ll + 3*NN + 4) = smt(k1, ll)
!!$           tjz(ll + 4*NN + 5) = -ll*smt(k1, ll); djz(ll + 4*NN + 5) = cmt(k1, ll)
!!$           tjz(ll + 5*NN + 6) =  ll*cmt(k1, ll); djz(ll + 5*NN + 6) = smt(k1, ll)
!!$        enddo ! end ll
!!$
!!$        dr  = sqrt( (dix(0)-djx(0))**2 + (diy(0)-djy(0))**2 + (diz(0)-djz(0))**2 ) + delta
!!$        mincc(icoil,jcoil) = min(dr, mincc(icoil,jcoil))
!!$
!!$        ccpotent(0,0) = ccpotent(0) + dli(0)*dlj(0)/dr**2 ! if takes a flexible exponent, how to normalize?
!!$
!!$        do ll = 1, 2*NN+2
!!$           ccpotent(ll       ,0) = ccpotent(ll       ,0) + dlj(0)*tix(0)*tix(ll       )/(dr**2*dli(0)) + 2*dli(0)*dlj(0)*(djx(0)-dix(0))*dix(ll       )/dr**4
!!$           ccpotent(ll+2*NN+2,0) = ccpotent(ll+2*NN+2,0) + dlj(0)*tiy(0)*tiy(ll+2*NN+2)/(dr**2*dli(0)) + 2*dli(0)*dlj(0)*(djy(0)-diy(0))*diy(ll+2*NN+2)/dr**4
!!$           ccpotent(ll+4*NN+4,0) = ccpotent(ll+4*NN+4,0) + dlj(0)*tiz(0)*tiz(ll+4*NN+4)/(dr**2*dli(0)) + 2*dli(0)*dlj(0)*(djz(0)-diz(0))*diz(ll+4*NN+4)/dr**4
!!$        enddo
!!$
!!$        if(mcoil .eq. icoil) then
!!$           do ll = 1, 2*NN+2
!!$              do mm = 1, 2*NN+2
!!$              ccpotent(ll       ,mm       ) = ccpotent(ll       ,mm       ) + dlj(0)*tix(ll       )*tix(mm       )/(dr**2*dli(0)) + &
!!$                                              dlj(0)*tix(0)*tix(ll       )*(2*dli(0)*(dix(0)-djx(0))*djx(mm       ) + dr**2*tix(0)*tix(mm       ))/(dr**4*dli(0)**2) - &
!!$                                              2*dli(0)*dlj(0)*dix(ll       )*dix(mm       )/dr**4 + 2*dlj(0)*(djx(0)-dix(0))*dix(ll       )*tix(0)*tix(mm       )/dr**4 -&
!!$                                              8*dli(0)*dlj(0)*(djx(0)-dix(0))**2*dix(ll       )*dix(mm       )/dr**5
!!$
!!$
!!$2*dli(0)*dlj(0)*(djx(0)-dix(0))*dix(ll       )/dr**4
!!$              ccpotent(ll+2*NN+2,mm       ) = ccpotent(ll+2*NN+2,mm       ) + dlj(0)*tiy(0)*tiy(ll+2*NN+2)/(dr**2*dli(0)) + 2*dli(0)*dlj(0)*(djy(0)-diy(0))*diy(ll+2*NN+2)/dr**4
!!$              ccpotent(ll+4*NN+4,mm       ) = ccpotent(ll+4*NN+4,mm       ) + dlj(0)*tiz(0)*tiz(ll+4*NN+4)/(dr**2*dli(0)) + 2*dli(0)*dlj(0)*(djz(0)-diz(0))*diz(ll+4*NN+4)/dr**4
!!$           enddo
!!$     enddo
!!$  enddo
!!$
!!$  ccpotent = ccpotent * pi2/NDcoil * pi2/NDcoil
!!$
!!$  return

