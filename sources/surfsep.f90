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
       icssep, mcssep, LM_fvec, LM_fjac, weight_cssep, Nteta, Nzeta, Nfp, MPI_COMM_FOCUS

  use mpi
  implicit none
  INTEGER, INTENT(in) :: ideriv
  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  INTEGER             :: astat, ierr, icoil, iteta, jzeta, NumGrid, idof, ND, ivec, totalcoil
  REAL                :: dcssep, discretefactor, d1S(1:Ndof), L1S(1:Ndof), coilsum
  REAL                :: lcssep(Ncoils), jac(Ncoils, Ndof)

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  cssep = zero ; lcssep = zero;
  NumGrid = Nteta*Nzeta     ! specifying the resolution of the prevent surface.
  discretefactor = (pi2/Nteta) * (pi2/Nzeta)
  lcssep = zero
  totalcoil = zero
  
  if( ideriv >= 0 ) then
     ivec = 1
     do icoil = 1, Ncoils
        if (coil(icoil)%type /= 1) cycle ! skip for other coils
        coilsum = zero
        if ( coil(icoil)%Lc /= 0 ) then
           if ( coil(icoil)%symm == 0 ) then
              totalcoil = totalcoil + 1
           elseif ( coil(icoil)%symm == 1 ) then
              totalcoil = totalcoil + Nfp
           elseif ( coil(icoil)%symm == 2 ) then 
              totalcoil = totalcoil + Nfp*2
           else
              FATAL( diagnos, .true. , Errors in coil symmetry )
           endif
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

     cssep = sum(lcssep(1:Ncoils)) * discretefactor /  (totalcoil + machprec) ! average value

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
           if (coil(icoil)%type == 1) then  ! skip for other coils
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
        t1S(idof) = sum(jac(1:Ncoils, idof)) * discretefactor /  (totalcoil + machprec)
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
  use globals, only: dp, zero, coil, myid, ounit, Ncoils, surf, psurf, cssep_factor, case_cssep, &
                     mincssep, cssep_alpha, cssep_beta, cssep_gamma, cssep_sigma, Nfp, pi2, MPI_COMM_FOCUS
  use mpi
  implicit none

  INTEGER, intent(in)  :: icoil, iteta, jzeta
  REAL   , intent(out) :: dcssep
  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  INTEGER              :: kseg, astat, ierr, j0, per0, l0, ss0, NS, cond
  REAL                 :: dl, xt, yt, zt, xc, yc, zc, xs, ys, zs, dx, dy, dz, lr
  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  FATAL( CSPotential0, icoil .lt. 1 .or. icoil .gt. Ncoils, icoil not in right range )

  FATAL( CSPotential0, cssep_alpha .lt. 0.0, coil to surface alpha cannot be negative )
  FATAL( CSPotential0, cssep_beta .lt. 2.0, coil to surface beta cannot be less than 2 )
  FATAL( CSPotential0, cssep_sigma .lt. 0.0, coil to surface sigma cannot be negative )
  FATAL( CSPotential0, cssep_gamma .lt. 1.0, coil to surface gamma cannot be less than 1 )
  FATAL( CSPotential0, mincssep .lt. 0.0, coil to surface separation cannot be negative )

  if( case_cssep .eq. 1 ) then
     cssep_alpha = 0.0
     cssep_sigma = 1.0
     cssep_gamma = cssep_factor
  endif
 
  dcssep = zero
  cond = 0

  xs = surf(psurf)%xx(iteta, jzeta); ys = surf(psurf)%yy(iteta, jzeta); zs = surf(psurf)%zz(iteta, jzeta)

  if ( coil(icoil)%symm == 0 ) then
     per0 = 1
     ss0 = 0
  elseif ( coil(icoil)%symm == 1 ) then
     per0 = Nfp
     ss0 = 0
  elseif ( coil(icoil)%symm == 2 ) then
     per0 = Nfp
     ss0 = 1
  else
     FATAL( diagnos, .true. , Errors in coil symmetry )
  endif
  do j0 = 1, per0
  do l0 = 0, ss0

  do kseg = 0, coil(icoil)%NS-1
     
     xc = (coil(icoil)%xx(kseg))*cos(pi2*(j0-1)/Nfp) - (coil(icoil)%yy(kseg))*sin(pi2*(j0-1)/Nfp)
     yc = ((-1.0)**l0)*((coil(icoil)%yy(kseg))*cos(pi2*(j0-1)/Nfp) + (coil(icoil)%xx(kseg))*sin(pi2*(j0-1)/Nfp))
     zc = (coil(icoil)%zz(kseg))*((-1.0)**l0)
     
     xt = (coil(icoil)%xt(kseg))*cos(pi2*(j0-1)/Nfp) - (coil(icoil)%yt(kseg))*sin(pi2*(j0-1)/Nfp)
     yt = ((-1.0)**l0)*((coil(icoil)%yt(kseg))*cos(pi2*(j0-1)/Nfp) + (coil(icoil)%xt(kseg))*sin(pi2*(j0-1)/Nfp))
     zt = (coil(icoil)%zt(kseg))*((-1.0)**l0)

     dl = xt*xt + yt*yt + zt*zt ! elment length (square)
     dx = xc - xs ; dy = yc - ys; dz = zc - zs ! r_vector, distance
     lr = dx*dx + dy*dy + dz*dz ! length**2 of r_vector

     if ( sqrt(lr) .lt. mincssep ) then
        if( log10(cssep_alpha*(mincssep-sqrt(lr))) .gt. 308.0 / cssep_beta ) then
           dcssep = HUGE(dcssep)
           return
        endif
        dcssep = dcssep + ( cssep_alpha*(mincssep-sqrt(lr)) )**cssep_beta*coil(icoil)%dd(kseg)*sqrt(dl)
     endif
     dcssep = dcssep + cssep_sigma*lr**(-0.5*cssep_gamma)*sqrt(dl)*coil(icoil)%dd(kseg)

  enddo ! end kseg

  enddo
  enddo

  call lenDeriv0( icoil, coil(icoil)%L )
  dcssep = dcssep / coil(icoil)%L
  
  return
END SUBROUTINE CSPotential0

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

SUBROUTINE CSPotential1(icoil, iteta, jzeta, d1S, ND)
!------------------------------------------------------------------------------------------------------ 
! DATE: 04/05/2018
! calculate the derivatives of the potential energy 
! between the i-th coil and the (iteta, jzeta) point on the surface
!------------------------------------------------------------------------------------------------------  
  use globals, only: dp, zero, coil, myid, ounit, Ncoils, surf, psurf, cssep_factor, DoF, case_cssep, &
                     mincssep, cssep_alpha, cssep_beta, cssep_gamma, cssep_sigma, Nfp, pi2, MPI_COMM_FOCUS
  use mpi
  implicit none

  INTEGER, intent(in)  :: icoil, iteta, jzeta, ND
  REAL   , intent(out) :: d1S(1:1, 1:ND)
  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  INTEGER              :: kseg, astat, ierr, j0, per0, l0, ss0
  REAL                 :: q, xt, yt, zt, xa, ya, za, xc, yc, zc, xs, ys, zs
  REAL                 :: dl, dx, dy, dz, lr, pm, holdd
  REAL, dimension(1:1, 0:coil(icoil)%NS-1) :: dSx, dSy, dSz
  REAL, dimension(0:coil(icoil)%NS-1, 1:ND) :: DoFx, DoFy, DoFz
  REAL, dimension(1:1,1:ND) :: d1L
  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  FATAL( CSPotential0, icoil .lt. 1 .or. icoil .gt. Ncoils, icoil not in right range )

  if( case_cssep .eq. 1 ) then
     cssep_alpha = 0.0
     cssep_sigma = 1.0
     cssep_gamma = cssep_factor
  endif

  d1S = zero
  q = cssep_factor ! easy convention
  xs = surf(psurf)%xx(iteta, jzeta) ; ys = surf(psurf)%yy(iteta, jzeta) ; zs = surf(psurf)%zz(iteta, jzeta)

  d1L = zero
  call lenDeriv0( icoil, coil(icoil)%L )
  call lenDeriv1( icoil, d1L(1:1,1:ND), ND )

  if ( coil(icoil)%symm == 0 ) then
     per0 = 1
     ss0 = 0
  elseif ( coil(icoil)%symm == 1 ) then
     per0 = Nfp
     ss0 = 0
  elseif ( coil(icoil)%symm == 2 ) then
     per0 = Nfp
     ss0 = 1
  else
     FATAL( diagnos, .true. , Errors in coil symmetry )
  endif
  
  do j0 = 1, per0
  do l0 = 0, ss0
 
  dSx = zero
  dSy = zero
  dSz = zero
  holdd = zero
 
  do kseg = 0, coil(icoil)%NS-1
     
     xc = (coil(icoil)%xx(kseg))*cos(pi2*(j0-1)/Nfp) - (coil(icoil)%yy(kseg))*sin(pi2*(j0-1)/Nfp)
     yc = ((-1.0)**l0)*((coil(icoil)%yy(kseg))*cos(pi2*(j0-1)/Nfp) + (coil(icoil)%xx(kseg))*sin(pi2*(j0-1)/Nfp))
     zc = (coil(icoil)%zz(kseg))*((-1.0)**l0)

     xt = (coil(icoil)%xt(kseg))*cos(pi2*(j0-1)/Nfp) - (coil(icoil)%yt(kseg))*sin(pi2*(j0-1)/Nfp)
     yt = ((-1.0)**l0)*((coil(icoil)%yt(kseg))*cos(pi2*(j0-1)/Nfp) + (coil(icoil)%xt(kseg))*sin(pi2*(j0-1)/Nfp))
     zt = (coil(icoil)%zt(kseg))*((-1.0)**l0)
     
     xa = (coil(icoil)%xa(kseg))*cos(pi2*(j0-1)/Nfp) - (coil(icoil)%ya(kseg))*sin(pi2*(j0-1)/Nfp)
     ya = ((-1.0)**l0)*((coil(icoil)%ya(kseg))*cos(pi2*(j0-1)/Nfp) + (coil(icoil)%xa(kseg))*sin(pi2*(j0-1)/Nfp))
     za = (coil(icoil)%za(kseg))*((-1.0)**l0)

     dl = xt*xt + yt*yt + zt*zt ! elment length (square)
     dx = xc - xs ; dy = yc - ys; dz = zc - zs ! r_vector, distance
     lr = dx*dx + dy*dy + dz*dz ! length**2 of r_vector
     
     if ( sqrt(lr) .lt. mincssep ) then
        dSx(1,kseg)=dSx(1,kseg)-1.0*cssep_alpha*cssep_beta*(cssep_alpha*(mincssep-sqrt(lr)))**(cssep_beta-1.0)*lr**(-0.5)*dx*sqrt(dl)
        dSy(1,kseg)=dSy(1,kseg)-1.0*cssep_alpha*cssep_beta*(cssep_alpha*(mincssep-sqrt(lr)))**(cssep_beta-1.0)*lr**(-0.5)*dy*sqrt(dl)
        dSz(1,kseg)=dSz(1,kseg)-1.0*cssep_alpha*cssep_beta*(cssep_alpha*(mincssep-sqrt(lr)))**(cssep_beta-1.0)*lr**(-0.5)*dz*sqrt(dl)
           
        dSx(1,kseg)=dSx(1,kseg)+cssep_alpha*cssep_beta*(cssep_alpha*(mincssep-sqrt(lr)))**(cssep_beta-1.0)*lr**(-0.5)*(dx*xt+dy*yt+dz*zt)*dl**(-0.5)*xt
        dSy(1,kseg)=dSy(1,kseg)+cssep_alpha*cssep_beta*(cssep_alpha*(mincssep-sqrt(lr)))**(cssep_beta-1.0)*lr**(-0.5)*(dx*xt+dy*yt+dz*zt)*dl**(-0.5)*yt
        dSz(1,kseg)=dSz(1,kseg)+cssep_alpha*cssep_beta*(cssep_alpha*(mincssep-sqrt(lr)))**(cssep_beta-1.0)*lr**(-0.5)*(dx*xt+dy*yt+dz*zt)*dl**(-0.5)*zt

        dSx(1,kseg)=dSx(1,kseg)-(cssep_alpha*(mincssep-sqrt(lr)))**cssep_beta*(-1.0*dl**(-1.5)*(xt*xa+yt*ya+zt*za)*xt+dl**(-0.5)*xa)
        dSy(1,kseg)=dSy(1,kseg)-(cssep_alpha*(mincssep-sqrt(lr)))**cssep_beta*(-1.0*dl**(-1.5)*(xt*xa+yt*ya+zt*za)*yt+dl**(-0.5)*ya)
        dSz(1,kseg)=dSz(1,kseg)-(cssep_alpha*(mincssep-sqrt(lr)))**cssep_beta*(-1.0*dl**(-1.5)*(xt*xa+yt*ya+zt*za)*zt+dl**(-0.5)*za)

        holdd = holdd + (cssep_alpha*(mincssep-sqrt(lr)))**cssep_beta*sqrt(dl)*coil(icoil)%dd(kseg)
     endif
     dSx(1,kseg)=dSx(1,kseg)-1.0*cssep_sigma*cssep_gamma*lr**(-0.5*(cssep_gamma+2.0))*dx*sqrt(dl)
     dSy(1,kseg)=dSy(1,kseg)-1.0*cssep_sigma*cssep_gamma*lr**(-0.5*(cssep_gamma+2.0))*dy*sqrt(dl)
     dSz(1,kseg)=dSz(1,kseg)-1.0*cssep_sigma*cssep_gamma*lr**(-0.5*(cssep_gamma+2.0))*dz*sqrt(dl)
        
     dSx(1,kseg)=dSx(1,kseg)+cssep_sigma*cssep_gamma*lr**(-0.5*(cssep_gamma+2.0))*(dx*xt+dy*yt+dz*zt)*dl**(-0.5)*xt
     dSy(1,kseg)=dSy(1,kseg)+cssep_sigma*cssep_gamma*lr**(-0.5*(cssep_gamma+2.0))*(dx*xt+dy*yt+dz*zt)*dl**(-0.5)*yt
     dSz(1,kseg)=dSz(1,kseg)+cssep_sigma*cssep_gamma*lr**(-0.5*(cssep_gamma+2.0))*(dx*xt+dy*yt+dz*zt)*dl**(-0.5)*zt
        
     dSx(1,kseg)=dSx(1,kseg)-cssep_sigma*lr**(-0.5*cssep_gamma)*(-1.0*dl**(-1.5)*(xt*xa+yt*ya+zt*za)*xt+dl**(-0.5)*xa)
     dSy(1,kseg)=dSy(1,kseg)-cssep_sigma*lr**(-0.5*cssep_gamma)*(-1.0*dl**(-1.5)*(xt*xa+yt*ya+zt*za)*yt+dl**(-0.5)*ya)
     dSz(1,kseg)=dSz(1,kseg)-cssep_sigma*lr**(-0.5*cssep_gamma)*(-1.0*dl**(-1.5)*(xt*xa+yt*ya+zt*za)*zt+dl**(-0.5)*za)
        
     dSx(1,kseg) = dSx(1,kseg)*coil(icoil)%dd(kseg)
     dSy(1,kseg) = dSy(1,kseg)*coil(icoil)%dd(kseg)
     dSz(1,kseg) = dSz(1,kseg)*coil(icoil)%dd(kseg)
        
     holdd = holdd + cssep_sigma*lr**(-0.5*cssep_gamma)*sqrt(dl)*coil(icoil)%dd(kseg)
     
  enddo ! end kseg

  DoFx = (DoF(icoil)%xof)*cos(pi2*(j0-1)/Nfp) - (DoF(icoil)%yof)*sin(pi2*(j0-1)/Nfp)
  DoFy = ((-1.0)**l0)*((DoF(icoil)%yof)*cos(pi2*(j0-1)/Nfp) + (DoF(icoil)%xof)*sin(pi2*(j0-1)/Nfp))
  DoFz = (DoF(icoil)%zof)*((-1.0)**l0)
  
  d1S(1:1, 1:ND) = d1S(1:1, 1:ND) + (matmul(dSx,DoFx)+matmul(dSy,DoFy)+matmul(dSz,DoFz))/coil(icoil)%L - holdd*d1L(1:1,1:ND)/(coil(icoil)%L)**2.0

  enddo 
  enddo 
  
  return

END SUBROUTINE CSPotential1
