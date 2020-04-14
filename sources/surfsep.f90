!title (coil-surface separation) ! Prevent coils too close to an arbitrary surface.

!latex \briefly{In this subroutine, we calculate the potential between the coil set and an arbitrary surface, which could be the plasma boundary. 
!latex By minimizing the potential, the distance between the coils and the ``prevent'' surface will be increased. 
!latex Due to the singularity of the formula, the coils would not be able to have intersections with the prevent surface.}

!latex \calledby{\link{solvers}}

!latex \tableofcontents

!latex \section{Interface}
!latex (reserved for future documentation)
!latex \section{Prevent surface}
!latex By default, the ``prevent" surface is the target plasma boundary. 
!latex In following updates, it will be modified to enable an arbitrary toroidal surface, which could be a surface away from the plasma with certain distance.
!latex 
!latex \section{Potential calculation}
!latex For each coil, the potential energy between coil and the prevent surface is calculated as
!latex \be \label{eq:surfsep}
!latex \ds E({\cal C}, {\cal S}) \equiv \oint_{\cal S} \int_{\cal C} \frac{\dd{s} \  \dd{l}}{|\vect{x}_c - \vect{x}_s|^q} \ ,
!latex \ee
!latex where ${\cal C}$ denotes the coil,  ${\cal S}$ denotes the prevent surface and $q$ is an user-specified exponential factor.
!latex The total energy for the entire coil set is then a summation over the coils,
!latex \be
!latex \ds f_S(\vect{X}) = \sum_{i = 1}^{Ncoils} E ({\cal C}_i, {\cal S}) \ .
!latex \ee
!latex 
!latex By minimizing $f_S$, the coil set would be driven away from the prevent surface.
!latex The exponential factor $q$ can accelerate/decelerate the ``preventing'' effect.
!latex Because of the singularity, the coils cannot cross the prevent surface.
!latex In such way, the coils will be well-separated from the plasma boundary.
!latex Of course, this requires the initial guesses have no intersections with the surface.
!latex 
!latex The discretized expression for \Eqn{surfsep} is
!latex \be
!latex \ds E({\cal C}, {\cal S}) \equiv \sum_{icoil=1}^{Ncoils} \sum_{jzeta=0}^{Nzeta-1} \sum_{iteta=0}^{Nteta-1} 
!latex \frac{\sqrt{g} \Delta \t \Delta \z \sqrt{{x'_c}^2 + {y'_c}^2 + {z'_c}^2} \Delta t}{\left ( \sqrt{(x_c - x_s)^2 + (y_c - y_s)^2 + (z_c - z_s)^2} \right )^q} \ .
!latex \ee
!latex 
!latex \section{First derivatives}
!latex The first order functional derivatives are 
!latex \begin{align}
!latex \ds \frac{\delta {E}}{\delta x_c} & = & \frac{1}{P} \left( -  \frac{q (x_c-{x_s}) \left({x'_c}^2+{y'_c}^2+{z'_c}^2\right)^2}{({x_s}-x_c)^2+({y_s}-y_c)^2+({z_s}-z_c)^2} + \frac{q x'_c \left({x'_c}^2+{y'_c}^2+{z'_c}^2\right) \left((x_c-{x_s}) x'_c+(y_c-{y_s}) y'_c+(z_c-{z_s}) z'_c\right)}{({x_s}-x_c)^2+({y_s}-y_c)^2+({z_s}-z_c)^2} \right. \nonumber \\
!latex & & \left. -   x''_c \left({x'_c}^2+{y'_c}^2+{z'_c}^2\right) + x'_c \left(x'_c x''_c+y'_c y''_c+z'_c z''_c\right)\right) \ , \\
!latex \ds \frac{\delta {E}}{\delta y_c} & = & \frac{1}{P} \left( -  \frac{q (y_c-{y_s}) \left({x'_c}^2+{y'_c}^2+{z'_c}^2\right)^2}{({x_s}-x_c)^2+({y_s}-y_c)^2+({z_s}-z_c)^2} +\frac{q y'_c \left({x'_c}^2+{y'_c}^2+{z'_c}^2\right) \left((x_c-{x_s}) x'_c+(y_c-{y_s}) y'_c+(z_c-{z_s}) z'_c\right)}{({x_s}-x_c)^2+({y_s}-y_c)^2+({z_s}-z_c)^2} \right. \nonumber \\
!latex & & \left. -  y''_c \left({x'_c}^2+{y'_c}^2+{z'_c}^2\right) + y'_c \left(x'_c x''_c+y'_c y''_c+z'_c z''_c\right)\right)  \ , \\
!latex \frac{\delta {E}}{\delta z_c} & = & \ds \frac{1}{P} \left( -  \frac{q (z_c-{z_s}) \left({x'_c}^2+{y'_c}^2+{z'_c}^2\right)^2}{({x_s}-x_c)^2+({y_s}-y_c)^2+({z_s}-z_c)^2} + \frac{q z'_c \left({x'_c}^2+{y'_c}^2+{z'_c}^2\right) \left((x_c-{x_s}) x'_c+(y_c-{y_s}) y'_c+(z_c-{z_s}) z'_c\right)}{({x_s}-x_c)^2+({y_s}-y_c)^2+({z_s}-z_c)^2} \right. \nonumber \\
!latex & & \left. -  z''_c \left({x'_c}^2+{y'_c}^2+{z'_c}^2\right) + z'_c \left(x'_c x''_c+y'_c y''_c+z'_c z''_c\right)\right) \ ,
!latex \end{align}
!latex where $P = [({x_s}-x_c)^2 + ({y_s}-y_c)^2 + ({z_s}-z_c)^2]^{q/2} [{x'_c}^2+{y'_c}^2+{z'_c}^2]^{3/2}$.
!latex In above equations, only the integrads are written out.
!latex \section{Second derivatives}
!latex (reserved for future development)
!latex \section{Comments and notes}
!latex \bi
!latex \item[1.] In the numerator of \Eqn{surfsep}, the coil length $\dd{l}$ is presented. 
!latex Thus, minimizing $f_S$ could lead to reducing coil lengths.
!latex Thinking about if the coil becomes one single point and $\dd{l} = 0$, $f_S = 0$.
!latex Fortunately, minimizing $f_B$ will normally drive coils to be longer and $f_L$ quadrtic term can keep the coil length to be around a constant value.
!latex In addition, the exponential factor $q$ (\inputvar{cssep\_factor}) could determine the piority between pushing coils away or condensate coils.
!latex \red{It's recommend to keep $q>=2$. }
!latex 
!latex \item[2.] A brief document about the function of cssep can be seen at {\href{https://docs.google.com/document/d/1vY5YNwe7T2OTKQsBIjotOe8QTtd3Z88qOV8FiJ3lztg/edit?usp=sharing}{google doc}}.
!latex \ei


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

SUBROUTINE surfsep(ideriv)
!------------------------------------------------------------------------------------------------------ 
! DATE: 04/05/2018
! calculate the potential energy (and derivatives) between coils and the "prevent" surface
!------------------------------------------------------------------------------------------------------  
  use globals, only: dp, zero, half, pi2, machprec, ncpu, myid, ounit, &
       coil, DoF, Ncoils, Nfixgeo, Ndof, cssep, t1S, t2S, psurf, surf, &
       icssep, mcssep, LM_fvec, LM_fjac, weight_cssep, MPI_COMM_FOCUS

  implicit none
  include "mpif.h"
  INTEGER, INTENT(in) :: ideriv
  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  INTEGER             :: astat, ierr, icoil, iteta, jzeta, NumGrid, Nteta, Nzeta, idof, ND, ivec
  REAL                :: dcssep, discretefactor, d1S(1:Ndof), L1S(1:Ndof), coilsum
  REAL                :: lcssep(Ncoils), jac(Ncoils, Ndof)

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  cssep = zero ; lcssep = zero;
  Nteta = surf(psurf)%Nteta ! please note psurf could be other surfaces than the plasma boundary.
  Nzeta = surf(psurf)%Nzeta ! Nteta & Nzeta are local variables here,
  NumGrid = Nteta*Nzeta     ! specifying the resolution of the prevent surface.
  discretefactor = (pi2/Nteta) * (pi2/Nzeta)
  !num_free = Ncoils - Nfixgeo ! number of free coils
  lcssep = zero
  
  if( ideriv >= 0 ) then
     ivec = 1
     do icoil = 1, Ncoils
        if (coil(icoil)%type /= 1) cycle ! skip for other coils
        coilsum = zero
        if ( coil(icoil)%Lc /= 0 ) then 
           do jzeta = 0, Nzeta - 1
              do iteta = 0, Nteta - 1           
                 if( myid.ne.modulo(jzeta*Nteta+iteta,ncpu) ) cycle ! parallelization loop;                     
                 call CSPotential0(icoil, iteta, jzeta, dcssep)
                 coilsum = coilsum + dcssep*surf(psurf)%ds(iteta, jzeta)  ! local cssep
              enddo ! end do iteta
           enddo ! end do jzeta
           call MPI_BARRIER( MPI_COMM_FOCUS, ierr )
           call MPI_REDUCE( coilsum, lcssep(icoil), 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_FOCUS, ierr )
           RlBCAST(lcssep(icoil), 1, 0 )
           ! L-M format of targets
           if (mcssep > 0) LM_fvec(ivec) = weight_cssep * lcssep(icoil)
           ivec = ivec + 1
        endif       
     enddo ! end do icoil

     cssep = sum(lcssep(1:Ncoils)) * discretefactor /  (Ncoils - Nfixgeo + machprec) ! average value

     ! L-M format of targets
     if (mcssep > 0) then
        FATAL( surfsep, ivec == (Ncoils-Nfixgeo), Errors in counting ivec )
     endif

  endif

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if ( ideriv >= 1 ) then

     t1S = zero ; jac = zero
     idof = 0 ; ivec = 1

     do icoil = 1, Ncoils
        ND = DoF(icoil)%ND
        d1S = zero
        l1S = zero
        if ( coil(icoil)%Ic /= 0 ) then ! if current is free;
           idof = idof +1
        endif

        if ( coil(icoil)%Lc /= 0 ) then ! if geometry is free;
           if (coil(icoil)%type /= 1) then  ! skip for other coils
              do jzeta = 0, Nzeta - 1
                 do iteta = 0, Nteta - 1
                    if( myid.ne.modulo(jzeta*Nteta+iteta,ncpu) ) cycle ! parallelization loop;
                    call CSPotential1(icoil, iteta, jzeta, d1S(idof+1:idof+ND), ND)
                    l1S(idof+1:idof+ND) = l1S(idof+1:idof+ND) + d1S(idof+1:idof+ND) * surf(psurf)%ds(iteta, jzeta)
                 enddo ! end do iteta
              enddo ! end do jzeta
              call MPI_BARRIER( MPI_COMM_FOCUS, ierr )
              call MPI_REDUCE( l1S, jac(icoil, 1:Ndof), Ndof, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_FOCUS, ierr )
              RlBCAST( jac(icoil, 1:Ndof), Ndof, 0 )
              ! L-M format of targets
              if (mcssep > 0)  LM_fjac(ivec, 1:Ndof) = weight_cssep * jac(icoil, 1:Ndof)           
              ivec = ivec + 1  
           endif
           idof = idof + ND ! ND should be zero if Lc==0
        endif        
                   
     enddo ! end do icoil
     FATAL( surfsep , idof .ne. Ndof, counting error in packing )

     ! L-M format of targets
     if (mcssep > 0) then
        FATAL( surfsep, ivec == (Ncoils-Nfixgeo), Errors in counting ivec )
     endif
     
     do idof = 1, Ndof       
        t1S(idof) = sum(jac(1:Ncoils, idof)) * discretefactor /  (Ncoils - Nfixgeo + machprec) 
     enddo
     
  endif

  return

END SUBROUTINE surfsep

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

SUBROUTINE CSPotential0(icoil, iteta, jzeta, dcssep)
!------------------------------------------------------------------------------------------------------ 
! DATE: 04/05/2018
! calculate the potential energy between the i-th coil and the (iteta, jzeta) point on the surface
!------------------------------------------------------------------------------------------------------  
  use globals, only: dp, zero, coil, myid, ounit, Ncoils, surf, psurf, cssep_factor, MPI_COMM_FOCUS
  implicit none
  include "mpif.h"

  INTEGER, intent(in)  :: icoil, iteta, jzeta
  REAL   , intent(out) :: dcssep
  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  INTEGER              :: kseg, astat, ierr
  REAL                 :: dl, xt, yt, zt, xc, yc, zc, xs, ys, zs, dx, dy, dz, lr
  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  FATAL( CSPotential0, icoil .lt. 1 .or. icoil .gt. Ncoils, icoil not in right range )
!  FATAL( CSPotential0, iteta .lt. 0 .or. iteta .gt. Nteta , iteta not in right range )
!  FATAL( CSPotential0, jzeta .lt. 0 .or. jzeta .gt. Nzeta , jzeta not in right range )
  
  dcssep = zero 

  xs = surf(psurf)%xx(iteta, jzeta) ; ys = surf(psurf)%yy(iteta, jzeta) ; zs = surf(psurf)%zz(iteta, jzeta)
  
  do kseg = 0, coil(icoil)%NS-1

     ! easy convention
     xc = coil(icoil)%xx(kseg) ; yc = coil(icoil)%yy(kseg) ; zc = coil(icoil)%zz(kseg)
     xt = coil(icoil)%xt(kseg) ; yt = coil(icoil)%yt(kseg) ; zt = coil(icoil)%zt(kseg)
     
     dl = xt*xt + yt*yt + zt*zt ! elment length (square)
     dx = xc - xs ; dy = yc - ys; dz = zc - zs ! r_vector, distance
     lr = dx*dx + dy*dy + dz*dz ! length**2 of r_vector

     dcssep = dcssep + sqrt(dl) / (lr**(cssep_factor/2.0)) * coil(icoil)%dd(kseg)

  enddo ! end kseg

  return
END SUBROUTINE CSPotential0

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

SUBROUTINE CSPotential1(icoil, iteta, jzeta, d1S, ND)
!------------------------------------------------------------------------------------------------------ 
! DATE: 04/05/2018
! calculate the derivatives of the potential energy 
! between the i-th coil and the (iteta, jzeta) point on the surface
!------------------------------------------------------------------------------------------------------  
  use globals, only: dp, zero, coil, myid, ounit, Ncoils, surf, psurf, cssep_factor, DoF, MPI_COMM_FOCUS
  implicit none
  include "mpif.h"

  INTEGER, intent(in)  :: icoil, iteta, jzeta, ND
  REAL   , intent(out) :: d1S(1:1, 1:ND)
  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  INTEGER              :: kseg, astat, ierr
  REAL                 :: q, xt, yt, zt, xa, ya, za, xc, yc, zc, xs, ys, zs
  REAL                 :: dl, dx, dy, dz, lr, pm
  REAL, dimension(1:1, 0:coil(icoil)%NS-1) :: dSx, dSy, dSz
  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  FATAL( CSPotential0, icoil .lt. 1 .or. icoil .gt. Ncoils, icoil not in right range )
  d1S = zero
  q = cssep_factor ! easy convention
  xs = surf(psurf)%xx(iteta, jzeta) ; ys = surf(psurf)%yy(iteta, jzeta) ; zs = surf(psurf)%zz(iteta, jzeta)

  do kseg = 0, coil(icoil)%NS-1

     ! easy convention
     xc = coil(icoil)%xx(kseg) ; yc = coil(icoil)%yy(kseg) ; zc = coil(icoil)%zz(kseg)
     xt = coil(icoil)%xt(kseg) ; yt = coil(icoil)%yt(kseg) ; zt = coil(icoil)%zt(kseg) 
     xa = coil(icoil)%xa(kseg) ; ya = coil(icoil)%ya(kseg) ; za = coil(icoil)%za(kseg)

     dl = xt*xt + yt*yt + zt*zt ! elment length (square)
     dx = xc - xs ; dy = yc - ys; dz = zc - zs ! r_vector, distance
     lr = dx*dx + dy*dy + dz*dz ! length**2 of r_vector
     pm = 1.0 / ( (lr**(q/2.0)) * (dl**(1.5)) ) * coil(icoil)%dd(kseg) ! 1/P * Delta t in the documentation

     dSx(1, kseg) = pm * ( -q*dx*dl*dl/lr + q*xt*dl*(dx*xt+dy*yt+dz*zt)/lr - xa*dl + xt*(xt*xa+yt*ya+zt*za) )
     dSy(1, kseg) = pm * ( -q*dy*dl*dl/lr + q*yt*dl*(dx*xt+dy*yt+dz*zt)/lr - ya*dl + yt*(xt*xa+yt*ya+zt*za) )
     dSz(1, kseg) = pm * ( -q*dz*dl*dl/lr + q*zt*dl*(dx*xt+dy*yt+dz*zt)/lr - za*dl + zt*(xt*xa+yt*ya+zt*za) )
     
  enddo ! end kseg

  d1S(1:1, 1:ND) = matmul(dSx, DoF(icoil)%xof) + matmul(dSy, DoF(icoil)%yof) + matmul(dSz, DoF(icoil)%zof)

  return

END SUBROUTINE CSPotential1

