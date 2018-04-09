
!title (toroidal flux) ! Construct toroidal flux constraints.

!latex \briefly{Based on the theory of MHD, the toroidal flux on each plasma flux surface should be identical. Thus, we can put a constraint on the total 
!latex          toroidal flux, rather than fixing the currents in coils.}

!latex \calledby{\link{denergy}}

!latex \tableofcontents

!latex \subsection{Toroidal flux}
!latex \bi
!latex \item[1.] Toroidal flux at a poloidal surface which is produced by cutting the plasma surface with a $\zeta$ = constant surface equals the line integral of magnetic  
!latex           vector potential over the boundary, based on the Stokes' theorem. That is,
!latex \be
!latex  \Phi_i & = & \int_{S_i} \vec{B} \cdot d\vec{s} \, = \, \int_{S_i} \nabla \times \vec{A} \cdot d\vec{s} \, = \, \int_{l_i} \vec{A} \cdot d\vec{l} \\
!latex \vec{A} & = & \sum_{j=1}^{Ncoils} I_j \int_{coil\_j} \frac{d\vec{l}}{r}
!latex \ee
!latex Here $i$ is denoted to the poloidal surface label.
!latex \item[2.] The total toroidal flux constraints then can be represented as,
!latex \be
!latex tflux & \equiv & \frac{1}{nzeta} \sum_{i=1}^{nzeta} \frac{1}{2} (\Phi_i - \Phi_o)^2
!latex \ee
!latex \ei
!latex \subsection{First derivatives}
!latex The first derivatives of toroidal flux cost function $tflux$ are derived as,
!latex \be
!latex \ds \frac{\partial{tflux}}{\partial{I^j}} & = & \frac{1}{nzeta} \sum_{i=1}^{nzeta} (\Phi_i - \Phi_o) \frac{\partial{\Phi_i}}{\partial{I^j}} \\ \nonumber
!latex \ds                                       & = & \frac{1}{nzeta} \sum_{i=1}^{nzeta} (\Phi_i - \Phi_o) \int_{l_i}  \int_{coil\_j} \frac{d\vec{l}}{r} \cdot d\vec{l}\\
!latex \ds \frac{\partial{tflux}}{\partial{x^j_n}} & = & \frac{1}{nzeta} \sum_{i=1}^{nzeta} (\Phi_i - \Phi_o) \frac{\partial{\Phi_i}}{\partial{x^j_n}} \\ \nonumber
!latex \ds                                         & = & \frac{1}{nzeta} \sum_{i=1}^{nzeta} (\Phi_i - \Phi_o) \int_{l_i} I^j \int_{coil\_j} 
!latex                                                                                           \frac{\partial{\frac{d\vec{l}}{r}}}{\partial{x^j_n}} \cdot d\vec{l}
!latex \ee
!latex Here, $j$ means the argument is about the $jth$ coil.
!latex \subsection{Second derivatives}
!latex Similarly, the second derivatives of $tflux$ can be written as,
!latex \be
!latex \ds \frac{\partial^2{tflux}}{\partial{X^j}\partial{X^k}} & = & \frac{1}{nzeta} \sum_{i=1}^{nzeta} \frac{\partial{\Phi_i}}{\partial{X^j}} \frac{\partial{\Phi_i}}{\partial{X^k}} 
!latex                                                               + \delta_j^k (\Phi_i - \Phi_o) \frac{\partial^2{\Phi_i}}{\partial{X^j}\partial{X^k}}
!latex \ee
!latex Here, $X$ represents all the DoFs, both the currents and geometry parameters.
!latex 
!latex \subsection{Normalization}
!latex Since $\Phi_o$ is a user sepcified constant and identical at each cross-section, the normalization for toroidal flux cost function can be implemented by dividing all the functions and derivatives with ${\Phi_o}^2$.
!latex In order to calculate the flux value first ( in which case, target\_flux would be reset to zero and can not be divied. ), this normalization is finished in costfun subroutine in \link{denergy}.
 
subroutine torflux(nderiv)
  use kmodule, only: zero, half, one, pi2, sqrtmachprec, bsconstant, &
                     coil, surf, cmt, smt, NFcoil, NDcoil, Ncoils, Nteta, Nzeta, discretefactor, Cdof, &
                     tflux, t1F, t2F, target_tflux, isign, &
                     ncpu, myid, ounit
  implicit none
  include "mpif.h"

  INTEGER           :: nderiv


  INTEGER           :: astat, ierr
  INTEGER           :: icoil, mm, NN, iteta, jzeta, kseg, ii, jj, kk, ll, c1, c2, n1, n2   ! icoil, iteta, jzeta, kseg are local
  REAL              :: r, rm2, rm3, rm4, bx, by, bz, ax, ay, az, lm, c12, lbnorm, lflux, dflux, flux_error, lsumflux, sumflux
  REAL, allocatable :: l1F( :, :), l2F( :, :, :, :), d1F( :, :), d2F( :, :, :, :)
  INTEGER           :: array2size, array4size

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  
  NN    = NFcoil
  Cdof  = 6 * NN + 6
  tflux = zero; lsumflux = zero; sumflux = zero ! sum of flux to decide the sign of isign
  
  do icoil = 1, Ncoils

   if ( .not. allocated(coil(icoil)%Ax) ) then
    allocate(coil(icoil)%Ax(0:Cdof, 0:Cdof), stat=astat)
    allocate(coil(icoil)%Ay(0:Cdof, 0:Cdof), stat=astat)
    allocate(coil(icoil)%Az(0:Cdof, 0:Cdof), stat=astat)
   endif

   coil(icoil)%Ax = zero
   coil(icoil)%Ay = zero
   coil(icoil)%Az = zero

  enddo

  select case ( nderiv )
   !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  case(0)  

   dflux = zero

   do jzeta = 0, Nzeta - 1

    if( myid.ne.modulo(jzeta,ncpu) ) cycle ! parallelization loop; 
    lflux = zero 

    do iteta = 0, Nteta - 1

     ax = zero; ay = zero; az = zero

     do icoil = 1, Ncoils

      call bpotential0(icoil,iteta,jzeta,coil(icoil)%Ax(0,0),coil(icoil)%Ay(0,0),coil(icoil)%Az(0,0))

      ax = ax + coil(icoil)%Ax(0, 0) * pi2 / NDcoil * coil(icoil)%I * bsconstant
      ay = ay + coil(icoil)%Ay(0, 0) * pi2 / NDcoil * coil(icoil)%I * bsconstant
      az = az + coil(icoil)%Az(0, 0) * pi2 / NDcoil * coil(icoil)%I * bsconstant

     enddo ! end do icoil

     lflux = lflux + ax * surf(1)%xt(iteta,jzeta) + ay * surf(1)%yt(iteta,jzeta) + az * surf(1)%zt(iteta,jzeta)

    enddo ! end do iteta

    lflux = lflux * pi2/Nteta
    lsumflux = lsumflux + lflux

    flux_error = lflux - target_tflux
    dflux = dflux + flux_error**2

   enddo ! end do jzeta
 
   call MPI_REDUCE( dflux, tflux, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
   RlBCAST(tflux, 1, 0)

   tflux = tflux * half / Nzeta

   !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  case(1)

   SALLOCATE( l1F, (1:Ncoils, 0:Cdof), zero )
   SALLOCATE( d1F, (1:Ncoils, 0:Cdof), zero )

   if ( .not. allocated(t1F) ) allocate(t1F(1:Ncoils, 0:Cdof), stat=astat)
   t1F = zero; dflux = zero

   do jzeta = 0, Nzeta - 1

    if( myid.ne.modulo(jzeta,ncpu) ) cycle ! parallelization loop; 
    lflux = zero; l1F   = zero 

    do iteta = 0, Nteta - 1

     ax = zero; ay = zero; az = zero

     do icoil = 1, Ncoils

      call bpotential1(icoil,iteta,jzeta,coil(icoil)%Ax(0,0:Cdof),coil(icoil)%Ay(0,0:Cdof),coil(icoil)%Az(0,0:Cdof))

      coil(icoil)%Ax = coil(icoil)%Ax * pi2 / NDcoil * bsconstant
      coil(icoil)%Ay = coil(icoil)%Ay * pi2 / NDcoil * bsconstant
      coil(icoil)%Az = coil(icoil)%Az * pi2 / NDcoil * bsconstant

      ax = ax + coil(icoil)%Ax(0, 0) * coil(icoil)%I
      ay = ay + coil(icoil)%Ay(0, 0) * coil(icoil)%I
      az = az + coil(icoil)%Az(0, 0) * coil(icoil)%I

     enddo ! end do icoil

     lflux = lflux + ax * surf(1)%xt(iteta,jzeta) + ay * surf(1)%yt(iteta,jzeta) + az * surf(1)%zt(iteta,jzeta)

     ! first derivatives of local flux
     do c1 = 1, Ncoils
      l1F(c1,0) = l1F(c1,0) + coil(c1)%Ax(0,0) * surf(1)%xt(iteta,jzeta) + coil(c1)%Ay(0,0) * surf(1)%yt(iteta,jzeta) + coil(c1)%Az(0,0) * surf(1)%zt(iteta,jzeta)
     enddo

     do n1 = 1, Cdof
      do c1 = 1, Ncoils       
       l1F(c1,n1) = l1F(c1,n1) + ( coil(c1)%Ax(0,n1) * surf(1)%xt(iteta,jzeta) + coil(c1)%Ay(0,n1) * surf(1)%yt(iteta,jzeta) + coil(c1)%Az(0,n1) * surf(1)%zt(iteta,jzeta) ) * coil(c1)%I 
      enddo
     enddo

    enddo ! end do iteta

    lflux = lflux * pi2/Nteta
    l1F   = l1F   * pi2/Nteta
    lsumflux = lsumflux + lflux

    flux_error = lflux - target_tflux
    dflux = dflux + flux_error**2

    do c1 = 1, Ncoils
     do n1 = 0, Cdof
      d1F(c1,n1) = d1F(c1,n1) + flux_error*l1F(c1,n1)
     enddo
    enddo
    
 !   d1F = d1F + flux_error * l1F

   enddo ! end do jzeta

   call MPI_REDUCE( dflux, tflux, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
   RlBCAST(tflux, 1, 0)

   array2size = Ncoils * ( Cdof + 1 )
   call MPI_REDUCE( d1F, t1F, array2size, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
   RlBCAST( t1F, array2size, 0 )

   tflux = tflux * half / Nzeta
   t1F   = t1F          / Nzeta

   deallocate( l1F )
   deallocate( d1F )

   !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  case(2)

   SALLOCATE( l1F, (1:Ncoils, 0:Cdof), zero )
   SALLOCATE( d1F, (1:Ncoils, 0:Cdof), zero )

   SALLOCATE( l2F, (1:Ncoils, 0:Cdof, 1:Ncoils, 0:Cdof), zero )
   SALLOCATE( d2F, (1:Ncoils, 0:Cdof, 1:Ncoils, 0:Cdof), zero )

   if ( .not. allocated(t1F) ) allocate(t1F(1:Ncoils, 0:Cdof), stat=astat)
   if ( .not. allocated(t2F) ) allocate(t2F(1:Ncoils, 0:Cdof, 1:Ncoils, 0:Cdof), stat=astat)
   t1F = zero; t2F = zero; dflux = zero

   do jzeta = 0, Nzeta - 1

    if( myid.ne.modulo(jzeta,ncpu) ) cycle ! parallelization loop; 
    lflux = zero; l1F = zero; l2F = zero

    do iteta = 0, Nteta - 1

     ax = zero; ay = zero; az = zero

     do icoil = 1, Ncoils

      call bpotential2(icoil,iteta,jzeta,coil(icoil)%Ax(0:Cdof,0:Cdof),coil(icoil)%Ay(0:Cdof,0:Cdof),coil(icoil)%Az(0:Cdof,0:Cdof))
      
      coil(icoil)%Ax = coil(icoil)%Ax * pi2 / NDcoil * bsconstant
      coil(icoil)%Ay = coil(icoil)%Ay * pi2 / NDcoil * bsconstant
      coil(icoil)%Az = coil(icoil)%Az * pi2 / NDcoil * bsconstant

      ax = ax + coil(icoil)%Ax(0, 0) * coil(icoil)%I
      ay = ay + coil(icoil)%Ay(0, 0) * coil(icoil)%I
      az = az + coil(icoil)%Az(0, 0) * coil(icoil)%I

     enddo ! end do icoil

     lflux = lflux + ax * surf(1)%xt(iteta,jzeta) + ay * surf(1)%yt(iteta,jzeta) + az * surf(1)%zt(iteta,jzeta)

     ! first derivatives
     do c1 = 1, Ncoils
      l1F(c1,0) = l1F(c1,0)               +   coil(c1)%Ax( 0, 0) * surf(1)%xt(iteta,jzeta) + coil(c1)%Ay( 0, 0) * surf(1)%yt(iteta,jzeta) + coil(c1)%Az( 0, 0) * surf(1)%zt(iteta,jzeta)
     enddo
     do n1 = 1, Cdof
      do c1 = 1, Ncoils
       l1F(c1,n1) = l1F(c1,n1)            + ( coil(c1)%Ax(n1, 0) * surf(1)%xt(iteta,jzeta) + coil(c1)%Ay(n1, 0) * surf(1)%yt(iteta,jzeta) + coil(c1)%Az(n1, 0) * surf(1)%zt(iteta,jzeta) ) * coil(c1)%I
      enddo
     enddo
     ! second derivatives
     do n1 = 1, Cdof
      do c1 = 1, Ncoils
       l2F(c1, 0,c1,n1) = l2F(c1, 0,c1,n1) +  coil(c1)%Ax(n1, 0) * surf(1)%xt(iteta,jzeta) + coil(c1)%Ay(n1, 0) * surf(1)%yt(iteta,jzeta) + coil(c1)%Az(n1, 0) * surf(1)%zt(iteta,jzeta)
       l2F(c1,n1,c1, 0) = l2F(c1,n1,c1, 0) +  coil(c1)%Ax(n1, 0) * surf(1)%xt(iteta,jzeta) + coil(c1)%Ay(n1, 0) * surf(1)%yt(iteta,jzeta) + coil(c1)%Az(n1, 0) * surf(1)%zt(iteta,jzeta)
      enddo
     enddo

     do n2 = 1, Cdof
      do n1 = 1, Cdof
       do c1 = 1, Ncoils
        l2F(c1,n1,c1,n2) = l2F(c1,n1,c1,n2) +(coil(c1)%Ax(n1,n2) * surf(1)%xt(iteta,jzeta) + coil(c1)%Ay(n1,n2) * surf(1)%yt(iteta,jzeta) + coil(c1)%Az(n1,n2) * surf(1)%zt(iteta,jzeta) ) * coil(c1)%I
       enddo
      enddo
     enddo

    enddo ! end do iteta

    lflux = lflux * pi2/Nteta
    l1F   = l1F   * pi2/Nteta
    l2F   = l2F   * pi2/Nteta
    lsumflux = lsumflux + lflux

    flux_error = lflux - target_tflux
    dflux = dflux + flux_error**2
    d1F = d1F + flux_error * l1F
           
    do c1 = 1, Ncoils
     do n1 = 0, Cdof
      do c2 = 1, Ncoils
       do n2 = 0, Cdof
        d2F(c1,n1,c2,n2) = d2F(c1,n1,c2,n2) + l1F(c2,n2)*l1F(c1,n1) + flux_error * l2F(c1,n1,c2,n2)
       enddo
      enddo
     enddo
    enddo


   enddo !enddo jzeta
    


   call MPI_REDUCE( dflux, tflux, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
   RlBCAST(tflux, 1, 0)

   array2size = Ncoils * ( Cdof + 1 )
   call MPI_REDUCE( d1F, t1F, array2size, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
   RlBCAST( t1F, array2size, 0 )

   array4size = array2size * array2size
   call MPI_REDUCE( d2F, t2F, array4size, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
   RlBCAST( t2F, array4size, 0 )


   tflux = tflux * half / Nzeta
   t1F   = t1F          / Nzeta
   t2F   = t2F          / Nzeta


   deallocate( l1F )
   deallocate( d1F )
   deallocate( l2F )
   deallocate( d2F )

  end select

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  do icoil = 1, Ncoils
   if (allocated(coil(icoil)%Ax) ) deallocate( coil(icoil)%Ax, coil(icoil)%Ay, coil(icoil)%Az )
  enddo

  call MPI_REDUCE( lsumflux, sumflux, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
  RlBCAST(sumflux, 1, 0)
  if(sumflux .lt. zero) isign = -1     ! target_tflux sign is decided by the sum of lflux

  return
end subroutine torflux

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine bpotential0(icoil, iteta, jzeta, Ax, Ay, Az)
  
  use kmodule, only: coil, surf, NDcoil, Ncoils, Nteta, Nzeta, &
                     zero, myid, ounit
  implicit none
  include "mpif.h"

  INTEGER, intent(in ) :: icoil, iteta, jzeta
  REAL   , intent(out) :: Ax, Ay, Az

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  INTEGER              :: ierr, astat, kseg
  REAL                 :: dlx, dly, dlz, r, ltx, lty, ltz

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  FATAL( bpotential0, icoil .lt. 1 .or. icoil .gt. Ncoils, icoil not in right range )
  FATAL( bpotential0, iteta .lt. 0 .or. iteta .gt. Nteta , iteta not in right range )
  FATAL( bpotential0, jzeta .lt. 0 .or. jzeta .gt. Nzeta , jzeta not in right range )
  
  dlx = zero; ltx = zero; Ax = zero
  dly = zero; lty = zero; Ay = zero
  dlz = zero; ltz = zero; Az = zero

  do kseg = 0, NDcoil - 1       ! NDcoil can be changed to coil(icoil)%D; 06/15/2016
        
   dlx = coil(icoil)%xx(kseg) - surf(1)%xx(iteta,jzeta)
   dly = coil(icoil)%yy(kseg) - surf(1)%yy(iteta,jzeta)
   dlz = coil(icoil)%zz(kseg) - surf(1)%zz(iteta,jzeta)

   r = sqrt(dlx**2 + dly**2 + dlz**2)

   ltx = coil(icoil)%xt(kseg)
   lty = coil(icoil)%yt(kseg)
   ltz = coil(icoil)%zt(kseg)
   
   Ax = Ax + ltx / r
   Ay = Ay + lty / r
   Az = Az + ltz / r

  enddo    ! enddo kseg

  return

end subroutine bpotential0

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine bpotential1(icoil, iteta, jzeta, Ax, Ay, Az)
  
  use kmodule, only: coil, surf, NFcoil, NDcoil, Ncoils, Cdof, Nteta, Nzeta, cmt, smt, &
                     zero, myid, ounit
  implicit none
  include "mpif.h"

  INTEGER, intent(in ) :: icoil, iteta, jzeta
  REAL   , intent(out) :: Ax(0:Cdof), Ay(0:Cdof), Az(0:Cdof)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  INTEGER              :: ierr, astat, kseg, NN, ll
  REAL                 :: dlx, dly, dlz, r, rm3, ltx, lty, ltz

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  FATAL( bpotential1, icoil .lt. 1 .or. icoil .gt. Ncoils, icoil not in right range )
  FATAL( bpotential1, iteta .lt. 0 .or. iteta .gt. Nteta , iteta not in right range )
  FATAL( bpotential1, jzeta .lt. 0 .or. jzeta .gt. Nzeta , jzeta not in right range )
  
  NN = NFcoil
  dlx = zero; ltx = zero; Ax = zero
  dly = zero; lty = zero; Ay = zero
  dlz = zero; ltz = zero; Az = zero

  do kseg = 0, NDcoil - 1       ! NDcoil can be changed to coil(icoil)%D; 06/15/2016
        
   dlx = coil(icoil)%xx(kseg) - surf(1)%xx(iteta,jzeta)
   dly = coil(icoil)%yy(kseg) - surf(1)%yy(iteta,jzeta)
   dlz = coil(icoil)%zz(kseg) - surf(1)%zz(iteta,jzeta)
   r = sqrt(dlx**2 + dly**2 + dlz**2); rm3 = r**(-3)

   ltx = coil(icoil)%xt(kseg)
   lty = coil(icoil)%yt(kseg)
   ltz = coil(icoil)%zt(kseg)

   Ax(0) = Ax(0) + ltx / r
   Ay(0) = Ay(0) + lty / r
   Az(0) = Az(0) + ltz / r

   do ll = 0, NN

    Ax(ll     +1) = Ax(ll     +1) - ll*smt(kseg,ll)/r - ltx*dlx*cmt(kseg,ll)*rm3 !Ax/xc
    Ax(ll+  NN+2) = Ax(ll+  NN+2) + ll*cmt(kseg,ll)/r - ltx*dlx*smt(kseg,ll)*rm3 !Ax/xs
    Ax(ll+2*NN+3) = Ax(ll+2*NN+3)                     - ltx*dly*cmt(kseg,ll)*rm3 !Ax/yc
    Ax(ll+3*NN+4) = Ax(ll+3*NN+4)                     - ltx*dly*smt(kseg,ll)*rm3 !Ax/ys
    Ax(ll+4*NN+5) = Ax(ll+4*NN+5)                     - ltx*dlz*cmt(kseg,ll)*rm3 !Ax/zc
    Ax(ll+5*NN+6) = Ax(ll+5*NN+6)                     - ltx*dlz*smt(kseg,ll)*rm3 !Ax/zs

    Ay(ll     +1) = Ay(ll     +1)                     - lty*dlx*cmt(kseg,ll)*rm3 !Ay/xc
    Ay(ll+  NN+2) = Ay(ll+  NN+2)                     - lty*dlx*smt(kseg,ll)*rm3 !Ay/xs
    Ay(ll+2*NN+3) = Ay(ll+2*NN+3) - ll*smt(kseg,ll)/r - lty*dly*cmt(kseg,ll)*rm3 !Ay/yc
    Ay(ll+3*NN+4) = Ay(ll+3*NN+4) + ll*cmt(kseg,ll)/r - lty*dly*smt(kseg,ll)*rm3 !Ay/ys
    Ay(ll+4*NN+5) = Ay(ll+4*NN+5)                     - lty*dlz*cmt(kseg,ll)*rm3 !Ay/zc
    Ay(ll+5*NN+6) = Ay(ll+5*NN+6)                     - lty*dlz*smt(kseg,ll)*rm3 !Ay/zs

    Az(ll     +1) = Az(ll     +1)                     - ltz*dlx*cmt(kseg,ll)*rm3 !Az/xc
    Az(ll+  NN+2) = Az(ll+  NN+2)                     - ltz*dlx*smt(kseg,ll)*rm3 !Az/xs
    Az(ll+2*NN+3) = Az(ll+2*NN+3)                     - ltz*dly*cmt(kseg,ll)*rm3 !Az/yc
    Az(ll+3*NN+4) = Az(ll+3*NN+4)                     - ltz*dly*smt(kseg,ll)*rm3 !Az/ys
    Az(ll+4*NN+5) = Az(ll+4*NN+5) - ll*smt(kseg,ll)/r - ltz*dlz*cmt(kseg,ll)*rm3 !Az/zc
    Az(ll+5*NN+6) = Az(ll+5*NN+6) + ll*cmt(kseg,ll)/r - ltz*dlz*smt(kseg,ll)*rm3 !Az/zs

   enddo ! enddo ll
   
  enddo    ! enddo kseg

  return

end subroutine bpotential1

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine bpotential2(icoil, iteta, jzeta, Ax, Ay, Az)

  use kmodule, only: coil, surf, NFcoil, NDcoil, Ncoils, Cdof, Nteta, Nzeta, cmt, smt, &
       zero, myid, ounit
  implicit none
  include "mpif.h"

  INTEGER, intent(in ) :: icoil, iteta, jzeta
  REAL   , intent(out) :: Ax(0:Cdof,0:Cdof), Ay(0:Cdof,0:Cdof), Az(0:Cdof,0:Cdof)

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  INTEGER              :: ierr, astat, kseg, NN, ll, mm
  REAL                 :: dlx(0:Cdof), dly(0:Cdof), dlz(0:Cdof), r, rm2, rm3, ltx(0:Cdof), lty(0:Cdof), ltz(0:Cdof), rp(0:Cdof,0:Cdof)

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  FATAL( bpotential1, icoil .lt. 1 .or. icoil .gt. Ncoils, icoil not in right range )
  FATAL( bpotential1, iteta .lt. 0 .or. iteta .gt. Nteta , iteta not in right range )
  FATAL( bpotential1, jzeta .lt. 0 .or. jzeta .gt. Nzeta , jzeta not in right range )

  NN = NFcoil
  dlx = zero; ltx = zero; Ax = zero
  dly = zero; lty = zero; Ay = zero
  dlz = zero; ltz = zero; Az = zero

  do kseg = 0, NDcoil - 1       ! NDcoil can be changed to coil(icoil)%D; 06/15/2016

   dlx(0) = coil(icoil)%xx(kseg) - surf(1)%xx(iteta,jzeta)
   dly(0) = coil(icoil)%yy(kseg) - surf(1)%yy(iteta,jzeta)
   dlz(0) = coil(icoil)%zz(kseg) - surf(1)%zz(iteta,jzeta)

   r = sqrt(dlx(0)**2 + dly(0)**2 + dlz(0)**2); rm2 = r**(-2); rm3 = r**(-3)

   ltx(0) = coil(icoil)%xt(kseg)
   lty(0) = coil(icoil)%yt(kseg)
   ltz(0) = coil(icoil)%zt(kseg)

   do ll = 0, NN
    dlx( ll        + 1 ) = cmt(kseg,ll)
    dlx( ll +   NN + 2 ) = smt(kseg,ll)
    dly( ll + 2*NN + 3 ) = cmt(kseg,ll)
    dly( ll + 3*NN + 4 ) = smt(kseg,ll)
    dlz( ll + 4*NN + 5 ) = cmt(kseg,ll)
    dlz( ll + 5*NN + 6 ) = smt(kseg,ll)

    ltx(ll        + 1) = -ll*smt(kseg, ll)
    ltx(ll +   NN + 2) =  ll*cmt(kseg, ll)
    lty(ll + 2*NN + 3) = -ll*smt(kseg, ll)
    lty(ll + 3*NN + 4) =  ll*cmt(kseg, ll)
    ltz(ll + 4*NN + 5) = -ll*smt(kseg, ll)
    ltz(ll + 5*NN + 6) =  ll*cmt(kseg, ll)
   enddo

   ! The derivatives of r 
   do ll = 1, 2*NN+2
    rp(ll       , 0) = dlx(0) * dlx(ll       ) / r
    rp(ll+2*NN+2, 0) = dly(0) * dly(ll+2*NN+2) / r
    rp(ll+4*NN+4, 0) = dlz(0) * dlz(ll+4*NN+4) / r
   enddo

   do mm = 1, 2*NN+2
    do ll = 1, 2*NN+2
     rp(ll       ,mm       ) = dlx(ll       ) * dlx(mm      ) / r  - dlx(0) * dlx(ll       ) * rp(mm,       0) * rm2
     rp(ll       ,mm+2*NN+2) =                                     - dlx(0) * dlx(ll       ) * rp(mm+2*Nn+2,0) * rm2
     rp(ll       ,mm+4*NN+4) =                                     - dlx(0) * dlx(ll       ) * rp(mm+4*Nn+4,0) * rm2

     rp(ll+2*NN+2,mm       ) =                                     - dly(0) * dly(ll+2*NN+2) * rp(mm       ,0) * rm2
     rp(ll+2*NN+2,mm+2*NN+2) = dly(ll+2*NN+2) * dly(mm+2*NN+2) / r - dly(0) * dly(ll+2*NN+2) * rp(mm+2*Nn+2,0) * rm2
     rp(ll+2*NN+2,mm+4*nn+4) =                                     - dly(0) * dly(ll+2*NN+2) * rp(mm+4*nn+4,0) * rm2

     rp(ll+4*NN+4,mm       ) =                                     - dlz(0) * dlz(ll+4*NN+4) * rp(mm       ,0) * rm2
     rp(ll+4*NN+4,mm+2*NN+2) =                                     - dlz(0) * dlz(ll+4*NN+4) * rp(mm+2*Nn+2,0) * rm2
     rp(ll+4*NN+4,mm+4*nn+4) = dlz(ll+4*NN+4) * dlz(mm+4*NN+4) / r - dlz(0) * dlz(ll+4*NN+4) * rp(mm+4*nn+4,0) * rm2
    enddo
   enddo

   ! derivatives of A
   Ax(0,0) = Ax(0,0) + ltx(0) / r
   Ay(0,0) = Ay(0,0) + lty(0) / r
   Az(0,0) = Az(0,0) + ltz(0) / r

   do ll = 1, Cdof
    Ax(ll, 0) = Ax(ll,0) + ltx(ll) / r - ltx(0) * rp(ll,0) * rm2
    Ay(ll, 0) = Ay(ll,0) + lty(ll) / r - lty(0) * rp(ll,0) * rm2
    Az(ll, 0) = Az(ll,0) + ltz(ll) / r - ltz(0) * rp(ll,0) * rm2
   enddo

   do mm = 1, Cdof
    do ll = 1, mm

     Ax(ll, mm) = Ax(ll,mm) - ( ltx(ll)*rp(mm, 0) + ltx(mm)*rp(ll,0) + ltx( 0)*rp(ll,mm) ) * rm2 + 2*ltx( 0)*rp(mm,0)*rp(ll,0)*rm3

     Ay(ll, mm) = Ay(ll,mm) - ( lty(ll)*rp(mm, 0) + lty(mm)*rp(ll,0) + lty( 0)*rp(ll,mm) ) * rm2 + 2*lty( 0)*rp(mm,0)*rp(ll,0)*rm3

     Az(ll, mm) = Az(ll,mm) - ( ltz(ll)*rp(mm, 0) + ltz(mm)*rp(ll,0) + ltz( 0)*rp(ll,mm) ) * rm2 + 2*ltz( 0)*rp(mm,0)*rp(ll,0)*rm3

    enddo
   enddo
  
   ! symmetric matrix
   do mm = 1, Cdof-1
    do ll = mm+1, Cdof

     Ax(ll, mm) = Ax(mm, ll)
     Ay(ll, mm) = Ay(mm, ll)
     Az(ll, mm) = Az(mm, ll)

    enddo
   enddo

 enddo    ! enddo kseg

 return

end subroutine bpotential2

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
