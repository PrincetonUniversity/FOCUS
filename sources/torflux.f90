
!title (toroidal flux) ! Construct toroidal flux constraints.

!latex \briefly{The toroidal magnetic flux of any cross-sections on the magnetic surface is constant. 
!latex          This can be used to constrain coil currents to avoid trival solutions.}

!latex \calledby{\link{solvers}}

!latex \section{General}
!latex  To avoid trivial solutions, like when $I_i \rightarrow 0$ $\forall i$, $f_B \rightarrow 0$, 
!latex  it is sufficient to constrain the enclosed toroidal flux.
!latex  If ${\bf B}\cdot{\bf n}=0$ on the boundary, then the toroidal flux through any poloidal 
!latex  cross-sectional surfaces is constant.
!latex  We include an objective function defined as
!latex  \be
!latex  \ds f_\Psi(\vect{X}) \equiv \frac{1}{2\pi} \int_0^{2\pi} \frac{1}{2} \left( \frac{\Psi_\z \; 
!latex   - \; \Psi_o}{\Psi_o} \right)^2 \dd{\z},
!latex  \ee
!latex  where the flux through a poloidal surface, ${\cal T}$, produced by cutting the boundary with plane 
!latex  $\z=const.$ is computed using Stokes' theorem,
!latex  \begin{align}
!latex  \ds \Psi_\z({\bf X}) & \equiv \int_{\cal T} \! \vect{B} \cdot \dd{\vect{S}} = 
!latex  \oint_{\partial {\cal T}} \!\! \vect{A} \cdot \dd{\vect{l}}.
!latex  \end{align}
!latex  Here $\dd \vect{l}$ is on the boundary curve of the poloidal surface and the total magnetic vector 
!latex  potential $\vect{A}$ is 
!latex  \begin{align} \vect{A}(\vect{X})  = \frac{\mu_0}{4\pi}\sum_{i=1}^{N_C} I_i \ \int_{C_i} 
!latex  \frac{\dd{\vect{l}_i}}{r}.
!latex  \end{align}
!latex  The variation of $f_\Psi$ resulting from $\delta \vect{x_i}$ is 
!latex  \begin{align}\label{eq:Otderiv}
!latex  \ds \delta f_\Psi(\vect{X}) &= \frac{1}{2\pi} \int_0^{2\pi} \left (\frac{\Psi_\z \; - \; 
!latex  \Psi_o}{\Psi_o} \right ) \ \frac{\delta \Psi_\z}{\Psi_o}\dd{\z} ,
!latex  \end{align}
!latex  where
!latex  $\ds \delta \Psi_\z = \int_{\partial {\cal T}} \delta \vect{A} \cdot \dd{\vect{l}}$ and 
!latex  \begin{align} \label{eq:var_psi}
!latex  \ds \delta \vect{A}(\vect{x}) &= \frac{\mu_0}{4\pi} I_i \ \int_0^{2\pi} \left [ -\frac{\vect{r} 
!latex  \cdot \vect{x}_i'}{r^3} \ \delta \vect{x}_i \, + \, \frac{\vect{r} \cdot \delta \vect{x}_i}{r^3} \ 
!latex  \vect{x}_i' \right ] \dd{t} .
!latex  \end{align}

!latex  \section{First derivatives}
!latex  We can write \Eqn{var_psi} into $x,y,z$ components, (subscript $i$ is omitting here)
!latex  \begin{align}
!latex  \ds \delta A_x & = \frac{\mu_0}{4\pi} I_i \ \int_0^{2\pi} \left [ \frac{\Delta x \dd{x} + 
!latex  \Delta y \dd{y} + \Delta z \dd{z}}{r^3} \ x' -  \frac{\Delta x x' + \Delta y y' + \Delta z z'}{r^3} 
!latex  \ \dd{x}  \right ] \dd{t} ; \\
!latex  \ds \delta A_y & = \frac{\mu_0}{4\pi} I_i \ \int_0^{2\pi} \left [ \frac{\Delta x \dd{x} + 
!latex  \Delta y \dd{y} + \Delta z \dd{z}}{r^3} \ y' - \frac{\Delta x x' + \Delta y y' + \Delta z z'}{r^3} 
!latex  \ \dd{y}  \right ] \dd{t} ; \\
!latex  \ds \delta A_z & = \frac{\mu_0}{4\pi} I_i \ \int_0^{2\pi} \left [ \frac{\Delta x \dd{x} + 
!latex  \Delta y \dd{y} + \Delta z \dd{z}}{r^3} \ z' - \frac{\Delta x x' + \Delta y y' + \Delta z z'}{r^3} 
!latex  \ \dd{z}  \right ] \dd{t} . 
!latex  \end{align}
!latex  
!latex  Here, we are applying $\vect{r} = \Delta x \ \vect{e_x} $ and $\vect{x}' = x' \vect{e_x}$. More
!latex  specifically, $\Delta x = x_{surf} - x_{coil} $ and $x' = \dd{x} / \dd{t}$.
!latex 
!latex  The first derivatives can be calculated as
!latex  \begin{align}
!latex  \ds \pdv{A_x}{x} & = \frac{\mu_0}{4\pi} I_i \ \int_0^{2\pi} 
!latex      \left [ -  \frac{\Delta y y' + \Delta z z'}{r^3}  \right ] \dd{t} ; \\
!latex  \ds \pdv{A_x}{y} & = \frac{\mu_0}{4\pi} I_i \ \int_0^{2\pi} 
!latex      \left [ \frac{\Delta y x'}{r^3} \right ] \dd{t} ; \\
!latex  \ds \pdv{A_x}{z} & = \frac{\mu_0}{4\pi} I_i \ \int_0^{2\pi} 
!latex      \left [ \frac{\Delta z x'}{r^3} \right ] \dd{t} .
!latex  \end{align}

!latex  \begin{align}
!latex  \ds \pdv{A_y}{x} & = \frac{\mu_0}{4\pi} I_i \ \int_0^{2\pi} 
!latex      \left [ \frac{\Delta x y'}{r^3}  \right ] \dd{t} ; \\
!latex  \ds \pdv{A_y}{y} & = \frac{\mu_0}{4\pi} I_i \ \int_0^{2\pi} 
!latex      \left [ - \frac{\Delta x x' + \Delta z z'}{r^3} \right ] \dd{t} ; \\
!latex  \ds \pdv{A_y}{z} & = \frac{\mu_0}{4\pi} I_i \ \int_0^{2\pi} 
!latex      \left [ \frac{\Delta z y'}{r^3} \right ] \dd{t} .
!latex  \end{align}

!latex  \begin{align}
!latex  \ds \pdv{A_z}{x} & = \frac{\mu_0}{4\pi} I_i \ \int_0^{2\pi} 
!latex      \left [ \frac{\Delta x z'}{r^3}  \right ] \dd{t} ; \\
!latex  \ds \pdv{A_z}{y} & = \frac{\mu_0}{4\pi} I_i \ \int_0^{2\pi} 
!latex      \left [ \frac{\Delta y z'}{r^3} \right ] \dd{t} ; \\
!latex  \ds \pdv{A_z}{z} & = \frac{\mu_0}{4\pi} I_i \ \int_0^{2\pi} 
!latex      \left [ \frac{\Delta x x' + \Delta y y'}{r^3} \right ] \dd{t} .
!latex  \end{align}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
subroutine torflux( ideriv )
!------------------------------------------------------------------------------------------------------ 
! DATE:  06/15/2017;
! Calculate the toroidal flux integal over several toroidal cross-sections on the flux surface;
! Nomarlization term target_tflux**2 will be divided outside in subroutine costfun;
! ideriv = 0 -> only calculate the toroidal flux constraint;
! ideriv = 1 -> calculate the toroidal flux constraint and its first derivatives;
! ideriv = 2 -> calculate the toroidal flux constraint and its first & second derivatives;
!------------------------------------------------------------------------------------------------------   
  use globals, only: dp, zero, half, one, pi2, sqrtmachprec, bsconstant, ncpu, myid, ounit, &
       coil, DoF, surf, Ncoils, Nteta, Nzeta, discretefactor, Cdof, &
       tflux, t1F, t2F, Ndof, psi_avg, target_tflux, tflux_sign, &
       itflux, mtflux, LM_fvec, LM_fjac, weight_tflux, plasma, MPI_COMM_FOCUS
  use mpi
  implicit none

  INTEGER, INTENT(in)                   :: ideriv
  !--------------------------------------------------------------------------------------------
  INTEGER                               :: astat, ierr
  INTEGER                               :: icoil, iteta, jzeta, idof, ND, isurf
  REAL                                  :: dflux, lflux, lsum
  REAL                                  :: lax, lay, laz          ! local Ax, Ay and Az
  REAL, dimension(0:Cdof, 0:Cdof)       :: dAx, dAy, dAz          ! dA of each coil;
  REAL, dimension(1:Ndof, 0:Nzeta-1)    :: ldF, dF
  REAL, dimension(0:Nzeta-1)            :: ldiff, psi_diff
  !--------------------------initialize and allocate arrays------------------------------------- 
  isurf = plasma
  tflux = zero ; lsum = zero ; psi_avg = zero ; dflux = zero ; psi_diff = zero
  ldiff = zero ; lax = zero; lay = zero; laz = zero      !already allocted; reset to zero;

  !-------------------------------calculate Bn-------------------------------------------------- 
  if( ideriv >= 0 ) then
     
     do jzeta = 0, Nzeta - 1
        if( myid.ne.modulo(jzeta,ncpu) ) cycle ! parallelization loop; 
         
        lflux = zero
        do iteta = 0, Nteta - 1
           lax = zero; lay = zero; laz = zero
           do icoil = 1, Ncoils
              call bpotential0(icoil, iteta, jzeta, dAx(0,0), dAy(0,0), dAz(0,0))
              lax = lax + dAx( 0, 0) * coil(icoil)%I * bsconstant
              lay = lay + dAy( 0, 0) * coil(icoil)%I * bsconstant
              laz = laz + dAz( 0, 0) * coil(icoil)%I * bsconstant
           enddo ! end do icoil
           lflux = lflux + lax * surf(isurf)%xt(iteta,jzeta) + &    ! local flux;
                           lay * surf(isurf)%yt(iteta,jzeta) + &
                           laz * surf(isurf)%zt(iteta,jzeta)
        enddo ! end do iteta
        lflux = lflux * pi2/Nteta * tflux_sign ! discretization factor;
        lsum  = lsum + lflux
        ldiff(jzeta) = lflux - target_tflux
        dflux = dflux + ldiff(jzeta)**2
     enddo ! end do jzeta

     call MPI_BARRIER( MPI_COMM_FOCUS, ierr )
     call MPI_REDUCE( dflux, tflux  , 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_FOCUS, ierr )
     call MPI_REDUCE( lsum , psi_avg, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_FOCUS, ierr )
     call MPI_REDUCE( ldiff, psi_diff, Nzeta, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_FOCUS, ierr )
               
     RlBCAST( psi_avg, 1, 0)
     RlBCAST( tflux, 1, 0)
     RlBCAST( psi_diff, Nzeta, 0)
     
     psi_avg = psi_avg / Nzeta
     tflux = half * tflux / Nzeta 

     ! Another type of target functions
     if (mtflux > 0) then
        LM_fvec(itflux+1:itflux+mtflux) = weight_tflux * psi_diff(0:Nzeta-1)
     endif
  
  endif

  !-------------------------------calculate tflux/x------------------------------------------------
  if ( ideriv >= 1 ) then

     ldF = zero; dF = zero; t1F = zero

     do jzeta = 0, Nzeta - 1
        if( myid.ne.modulo(jzeta,ncpu) ) cycle ! parallelization loop;
        do iteta = 0, Nteta - 1             
           idof = 0
           do icoil = 1, Ncoils
              ND = DoF(icoil)%ND
              if ( coil(icoil)%Ic /= 0 ) then !if current is free;
                 call bpotential0(icoil, iteta, jzeta, &
                      & dAx(0,0), dAy(0,0), dAz(0,0))

                 ldF(idof+1, jzeta) = ldF(idof+1, jzeta) &
                      & + bsconstant * ( dAx(0,0)*surf(isurf)%xt(iteta,jzeta)   &
                      &                + dAy(0,0)*surf(isurf)%yt(iteta,jzeta)   &
                      &                + dAz(0,0)*surf(isurf)%zt(iteta,jzeta) )
                 idof = idof +1
              endif

              if ( coil(icoil)%Lc /= 0 ) then !if geometry is free;
                 call bpotential1(icoil, iteta, jzeta, &
                      &       dAx(1:ND,0), dAy(1:ND,0), dAz(1:ND,0), ND)

                 ldF(idof+1:idof+ND, jzeta) = ldF(idof+1:idof+ND, jzeta) &
                      & + bsconstant * coil(icoil)%I * ( dAx(1:ND,0)*surf(isurf)%xt(iteta,jzeta)   &
                      &                                + dAy(1:ND,0)*surf(isurf)%yt(iteta,jzeta)   &
                      &                                + dAz(1:ND,0)*surf(isurf)%zt(iteta,jzeta) )

                 idof = idof + ND
              endif

           enddo !end icoil;
           FATAL( torflux , idof .ne. Ndof, counting error in packing )

        enddo !end iteta;
     enddo !end jzeta

     ldF = ldF * pi2/Nteta * tflux_sign

     call MPI_BARRIER( MPI_COMM_FOCUS, ierr )
     call MPI_REDUCE(ldF, dF, Ndof*Nzeta, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_FOCUS, ierr )
     RlBCAST( dF, Ndof*Nzeta, 0 )

     do idof = 1, Ndof
        t1F(idof) = sum( psi_diff(0:Nzeta-1) * dF(idof, 0:Nzeta-1) ) / Nzeta
     enddo

     ! Another type of target functions
     if (mtflux > 0) then
        do idof = 1, Ndof
           LM_fjac(itflux+1:itflux+mtflux, idof) = dF(idof, 0:Nzeta-1)
        enddo
     endif
     
  endif

  !--------------------------------------------------------------------------------------------

  call MPI_barrier( MPI_COMM_FOCUS, ierr )

  return
end subroutine torflux

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine bpotential0(icoil, iteta, jzeta, tAx, tAy, tAz)
!------------------------------------------------------------------------------------------------------ 
! DATE:  06/15/2017; 01/20/2020
! calculate the magnetic potential from coil(icoil) at the evaluation point (iteta, jzeta);
! Biot-Savart constant and currents are not included for later simplication.
! Discretizing factor is includeed; coil(icoil)%dd(kseg) 
!------------------------------------------------------------------------------------------------------   
  use globals, only: dp, coil, surf, Ncoils, Nteta, Nzeta, MPI_COMM_FOCUS,  &
                     zero, myid, ounit, plasma, Nfp, cosnfp, sinnfp, two, bsconstant,coil_type_spline
  use mpi
  implicit none

  INTEGER, intent(in ) :: icoil, iteta, jzeta
  REAL   , intent(out) :: tAx, tAy, tAz

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  INTEGER              :: ierr, astat, kseg, isurf, ip, is, cs, Npc
  REAL                 :: dlx, dly, dlz, rm, ltx, lty, ltz, &
       &                  Ax, Ay, Az, xx, yy, zz, rr

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  isurf = plasma
  FATAL( bpotential0, icoil .lt. 1 .or. icoil .gt. Ncoils, icoil not in right range )
  FATAL( bpotential0, iteta .lt. 0 .or. iteta .gt. Nteta , iteta not in right range )
  FATAL( bpotential0, jzeta .lt. 0 .or. jzeta .gt. Nzeta , jzeta not in right range )
  ! initialization
  Npc = 1 ; cs = 0  
  tAx = zero ; tAy = zero ; tAz = zero
  dlx = zero ; dly = zero ; dlz = zero
  ltx = zero ; lty = zero ; ltz = zero
  ! check if the coil is stellarator symmetric
  select case (coil(icoil)%symm) 
  case ( 0 )
     cs  = 0
     Npc = 1
  case ( 1 )
     cs  = 0
     Npc = Nfp
  case ( 2) 
     cs  = 1
     Npc = Nfp
  end select

  ! periodicity and stellarator symmetry
  do ip = 1, Npc
     do is = 0, cs
        ! find the point on plasma by rotating in reverse direction. + symmetric
        xx = ( surf(isurf)%xx(iteta,jzeta)*cosnfp(ip) + surf(isurf)%yy(iteta,jzeta)*sinnfp(ip) )
        yy = (-surf(isurf)%xx(iteta,jzeta)*sinnfp(ip) + surf(isurf)%yy(iteta,jzeta)*cosnfp(ip) ) * (-1)**is
        zz =   surf(isurf)%zz(iteta,jzeta) * (-1)**is
        Ax = zero; Ay = zero; Az = zero
        select case (coil(icoil)%type)        
        case(1,coil_type_spline)        
           ! Fourier coils
           do kseg = 0, coil(icoil)%NS-1
              dlx = xx - coil(icoil)%xx(kseg)
              dly = yy - coil(icoil)%yy(kseg)
              dlz = zz - coil(icoil)%zz(kseg)
              rm  = 1.0 / sqrt(dlx**2 + dly**2 + dlz**2)
              ltx = coil(icoil)%xt(kseg)
              lty = coil(icoil)%yt(kseg)
              ltz = coil(icoil)%zt(kseg)
              Ax = Ax + ltx * rm * coil(icoil)%dd(kseg)
              Ay = Ay + lty * rm * coil(icoil)%dd(kseg)
              Az = Az + ltz * rm * coil(icoil)%dd(kseg)
           enddo    ! enddo kseg
        case(3)
           ! central current and vertical field (zero contribution)
           rr = sqrt( xx**2 + yy**2 )
           ! \vec A=-\frac{\mu_0I}{2\pi} \ln(r) \hat e_z
           Az = -two*bsconstant*log(rr) 
        case default
           FATAL(bpotential0, .true., not supported coil types)
        end select 
        ! sum all the contributions
        tAx = tAx + (Ax*cosnfp(ip) - Ay*sinnfp(ip))*(-1)**is
        tAy = tAy + (Ay*cosnfp(ip) + Ax*sinnfp(ip))
        tAz = tAz +  Az
     enddo
  enddo

  return

end subroutine bpotential0

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine bpotential1(icoil, iteta, jzeta, tAx, tAy, tAz, ND)
!------------------------------------------------------------------------------------------------------ 
! DATE:  06/15/2017
! calculate the magnetic potential and its 1st derivatives from coil(icoil) at the evaluation point;
! Biot-Savart constant and currents are not included for later simplication.
! Discretizing factor is includeed; coil(icoil)%dd(kseg) 
!------------------------------------------------------------------------------------------------------    
  use globals, only: dp, coil, DoF, surf, NFcoil, Ncoils, Nteta, Nzeta, &
                     zero, myid, ounit, plasma, Nfp, cosnfp, sinnfp, MPI_COMM_FOCUS,coil_type_spline
  use mpi
  implicit none

  INTEGER, intent(in ) :: icoil, iteta, jzeta, ND
  REAL, dimension(1:1, 1:ND), intent(inout) :: tAx, tAy, tAz

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  INTEGER              :: ierr, astat, kseg, NS, isurf, ip, is, cs, Npc
  REAL                 :: dlx, dly, dlz, r, rm3, ltx, lty, ltz, xx, yy, zz
  REAL, dimension(1:1, 1:ND) :: Ax, Ay, Az
  REAL, dimension(1:1, 0:coil(icoil)%NS-1)   :: dAxx, dAxy, dAxz, dAyx, dAyy, dAyz, dAzx, dAzy, dAzz

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  isurf = plasma 
  FATAL( bpotential1, icoil .lt. 1 .or. icoil .gt. Ncoils, icoil not in right range )
  FATAL( bpotential1, iteta .lt. 0 .or. iteta .gt. Nteta , iteta not in right range )
  FATAL( bpotential1, jzeta .lt. 0 .or. jzeta .gt. Nzeta , jzeta not in right range )
  FATAL( bpotential1, ND <= 0, wrong inout dimension of ND )

  ! initialization
  Npc = 1 ; cs = 0  
  tAx = zero ; tAy = zero ; tAz = zero
  dlx = zero ; dly = zero ; dlz = zero
  ltx = zero ; lty = zero ; ltz = zero
  ! check if the coil is stellarator symmetric
  select case (coil(icoil)%symm) 
  case ( 0 )
     cs  = 0
     Npc = 1
  case ( 1 )
     cs  = 0
     Npc = Nfp
  case ( 2) 
     cs  = 1
     Npc = Nfp
  end select
  
  NS = coil(icoil)%NS
  ! periodicity and stellarator symmetry
  do ip = 1, Npc
     do is = 0, cs
        ! find the point on plasma by rotating in reverse direction. + symmetric
        xx = ( surf(isurf)%xx(iteta,jzeta)*cosnfp(ip) + surf(isurf)%yy(iteta,jzeta)*sinnfp(ip) )
        yy = (-surf(isurf)%xx(iteta,jzeta)*sinnfp(ip) + surf(isurf)%yy(iteta,jzeta)*cosnfp(ip) ) * (-1)**is
        zz =   surf(isurf)%zz(iteta,jzeta) * (-1)**is
        Ax = zero; Ay = zero; Az = zero
        select case (coil(icoil)%type)        
        case(1,coil_type_spline)        
           ! Fourier coils
           do kseg = 0, NS-1
              dlx = xx - coil(icoil)%xx(kseg)
              dly = yy - coil(icoil)%yy(kseg)
              dlz = zz - coil(icoil)%zz(kseg)     
              r = sqrt(dlx**2 + dly**2 + dlz**2); rm3 = r**(-3)
              ltx = coil(icoil)%xt(kseg)
              lty = coil(icoil)%yt(kseg)
              ltz = coil(icoil)%zt(kseg)
              dAxx(1,kseg) = - (dly*lty + dlz*ltz) * rm3 * coil(icoil)%dd(kseg) !Ax/x
              dAxy(1,kseg) =    dly*ltx            * rm3 * coil(icoil)%dd(kseg) !Ax/y
              dAxz(1,kseg) =    dlz*ltx            * rm3 * coil(icoil)%dd(kseg) !Ax/z
              dAyx(1,kseg) =    dlx*lty            * rm3 * coil(icoil)%dd(kseg) !Ay/x
              dAyy(1,kseg) = - (dlx*ltx + dlz*ltz) * rm3 * coil(icoil)%dd(kseg) !Ay/y
              dAyz(1,kseg) =    dlz*lty            * rm3 * coil(icoil)%dd(kseg) !Ay/z
              dAzx(1,kseg) =    dlx*ltz            * rm3 * coil(icoil)%dd(kseg) !Az/x
              dAzy(1,kseg) =    dly*ltz            * rm3 * coil(icoil)%dd(kseg) !Az/y
              dAzz(1,kseg) = - (dlx*ltx + dly*lty) * rm3 * coil(icoil)%dd(kseg) !Az/z
           enddo    ! enddo kseg
           Ax(1:1, 1:ND) = matmul(dAxx, DoF(icoil)%xof) + matmul(dAxy, DoF(icoil)%yof) + matmul(dAxz, DoF(icoil)%zof)
           Ay(1:1, 1:ND) = matmul(dAyx, DoF(icoil)%xof) + matmul(dAyy, DoF(icoil)%yof) + matmul(dAyz, DoF(icoil)%zof)
           Az(1:1, 1:ND) = matmul(dAzx, DoF(icoil)%xof) + matmul(dAzy, DoF(icoil)%yof) + matmul(dAzz, DoF(icoil)%zof)
        case(3) 
           continue
        case default
           FATAL(bpotential1, .true., not supported coil types)
        end select
        ! sum all the contributions
        tAx = tAx + (Ax*cosnfp(ip) - Ay*sinnfp(ip))*(-1)**is
        tAy = tAy + (Ay*cosnfp(ip) + Ax*sinnfp(ip))
        tAz = tAz +  Az
     enddo
  enddo
  return

end subroutine bpotential1

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
