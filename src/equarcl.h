!title (equarcl) ! EQUAl ARC Length constraints. 

!latex \briefly{Force each discretized coil arc has the same length to eliminate the non-uniqueness of Fourier harmonics.}

!latex \calledby{\link{denergy}}

!latex \tableofcontents

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsection{Introductiom}
!latex The Fourier harmonics depends on the selection of Fourier angel, which is not unique. Thus, it makes the solution for coil optimization not unique and optimizors applying
!latex  Newton method would fail to find the minimum of cost functions. To eliminate the non-uniqueness of Fourier harmonics, one general way is to add an additional 
!latex spectral constraints. Details can be found in \href{http://w3.pppl.gov/~shudson/Spec/bb00aa.pdf}{\blu{Spec:bb00aa}}. Here, what we used is just simply forcing each discretized 
!latex coil arc has the same length.
!latex 
!latex \subsection{Equal arc constraint}
!latex \bi
!latex \item[1.] The arc length for a discretized coil can be represented as, 
!latex \begin{equation} dl = \sqrt{ \dot{x}^2 + \dot{y}^2 + \dot{z}^2} \end{equation}
!latex \item[2.] Assuming that all arcs in a coil have the same length. Then, the coil length is $L = 2\pi * dl$, where $2\pi$ is the discretization number.
!latex \item[3.] Then our equal-arc constraint can be written in,
!latex \begin{equation}
!latex eqarc = \int_{coils} \int_{icoil} \frac{1}{2}( \sqrt{ \dot{x}^2 + \dot{y}^2 + \dot{z}^2} - \frac{L_i}{2\pi})^2 dt
!latex \end{equation}
!latex After normalized,
!latex \begin{equation}
!latex eqarc = \int_{coils} \int_{icoil} \frac{1}{2}( \frac{2\pi dl}{L_i} - 1 )^2 dt
!latex \end{equation}
!latex \ei
!latex 
!latex \subsection{First derivatives}
!latex \bi
!latex \item[1.]
!latex The first derivatives for equal-arc constraint are,
!latex \be
!latex \frac{\partial{eqarc}}{\partial{x_m^i}} & = & \int_{icoil} ( \frac{2\pi dl}{L_i} - 1 ) (\frac{2\pi}{L_i} \frac{\partial{dl}}{\partial{x_m^i}} \, - \, 
!latex  \frac{2\pi dl}{{L_i}^2} \frac{\partial{L_i}}{\partial{x_m^i}} ) \\
!latex \frac{\partial{dl}}{\partial{x_m^i}} & = &  \frac{\dot{x}}{dl} \frac{\partial{\dot{x}}}{\partial{x_m^i}} \\
!latex \frac{\partial{L_i}}{\partial{x_m^i}} & = & \int \frac{\partial{dl}}{\partial{x_m^i}} dt \\
!latex \ee
!latex \ei
!latex 
!latex \subsection{Second derivatives}
!latex \bi
!latex \item[1.] The second derivatives for equal-arc constraint are,
!latex \be 
!latex \ds \frac{\partial^2{eqarc}}{\partial{x_m^i} \partial{x_n^i}}   =   \int_{icoil} (\frac{2\pi}{L_i} \frac{\partial{dl}}{\partial{x_m^i}} \, - \, 
!latex  \frac{2\pi dl}{{L_i}^2} \frac{\partial{L_i}}{\partial{x_m^i}} ) (\frac{2\pi}{L_i} \frac{\partial{dl}}{\partial{x_n^i}} \, - \, 
!latex  \frac{2\pi dl}{{L_i}^2} \frac{\partial{L_i}}{\partial{x_n^i}} ) \; + \; \nonumber \\ 
!latex ( \frac{2\pi dl}{L_i} - 1 ) \left ( 
!latex  \frac{2\pi}{L_i} \frac{\partial^2{dl}}{\partial{x_m^i}\partial{x_n^i}} - \frac{2\pi}{{L_i}^2} \frac{\partial{dl}}{\partial{x_m^i}} \frac{\partial{L_i}}{\partial{x_n^i}}
!latex - \frac{2\pi}{{L_i}^2} \frac{\partial{dl}}{\partial{x_n^i}} \frac{\partial{L_i}}{\partial{x_m^i}} - \frac{2\pi dl}{{L_i}^2} \frac{\partial^2{L_i}}{\partial{x_m^i}\partial{x_n^i}}
!latex + \frac{4\pi dl}{{L_i}^3}\frac{\partial{L_i}}{\partial{x_m^i}} \frac{\partial{L_i}}{\partial{x_n^i}} \right )
!latex \ee
!latex \ei
!latex \subsection{Notes}
!latex \emph{For unfixed bugs with normalized version and also for the consideration of changing minimum condition, I only constructed the version without normalization. The normalization
!latex may be added sometime later.}
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!!$
!!$SUBROUTINE equarcl( nderiv )
!!$  use kmodule, only : zero, half, one, pi2, myid, ncpu, ounit, &
!!$                      Cdof, Ncoils, NDcoil, NFcoil, coil, smt, cmt, &
!!$                      eqarc, t1A, t2A
!!$  implicit none
!!$  include "mpif.h"
!!$
!!$
!!$  INTEGER           :: nderiv
!!$
!!$!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!!$  INTEGER           :: astat, ierr, icoil, kseg, ll, mm, NN
!!$  REAL              :: arclength, dlength
!!$  REAL, allocatable :: d1A(:,:), d1L(:), ddl1(:), ddl2(:,:), d2A(:,:,:,:), d2L(:,:)
!!$
!!$  NN      = NFcoil
!!$  
!!$  ! calculate each coil's length first
!!$  do icoil = 1, Ncoils
!!$   call LenDeriv0( icoil, coil(icoil)%L )
!!$  enddo
!!$  
!!$  eqarc = zero
!!$
!!$  select case ( nderiv )
!!$
!!$!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!!$
!!$  case(0)
!!$
!!$   do icoil = 1, Ncoils
!!$
!!$    arclength = coil(icoil)%L / pi2
!!$
!!$    do kseg = 0, NDcoil-1
!!$     dlength = sqrt( coil(icoil)%xt(kseg)**2 + coil(icoil)%yt(kseg)**2 + coil(icoil)%zt(kseg)**2 )
!!$     eqarc = eqarc + ( (dlength - arclength) / arclength )**2 !normalized
!!$    enddo ! end kseg
!!$
!!$   enddo ! end icoil
!!$
!!$   eqarc = eqarc * half * pi2 / NDcoil / Ncoils
!!$
!!$   !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!!$
!!$  case(1)
!!$
!!$   SALLOCATE(d1A, (1:Ncoils, 0:Cdof), zero)
!!$   SALLOCATE(d1L, (          0:Cdof), zero)
!!$   SALLOCATE(ddl1,(          0:Cdof), zero)
!!$
!!$   if ( .not. allocated(t1A) ) then
!!$    SALLOCATE(t1A, (1:Ncoils, 0:Cdof), zero)
!!$   else
!!$    t1A = zero
!!$   endif
!!$
!!$   do icoil = 1, Ncoils
!!$
!!$    arclength = coil(icoil)%L/pi2
!!$    call LenDeriv1( icoil, d1L(0:Cdof) )
!!$
!!$    do kseg = 0, NDcoil-1
!!$
!!$     dlength = sqrt( coil(icoil)%xt(kseg)**2 + coil(icoil)%yt(kseg)**2 + coil(icoil)%zt(kseg)**2 )
!!$     eqarc = eqarc + ( (dlength - arclength) / arclength )**2
!!$
!!$     do ll = 0, NN 
!!$
!!$      ddl1(ll        + 1) = coil(icoil)%xt(kseg)*(-ll*smt(kseg,ll))/dlength
!!$      ddl1(ll +   NN + 2) = coil(icoil)%xt(kseg)*( ll*cmt(kseg,ll))/dlength
!!$      ddl1(ll + 2*NN + 3) = coil(icoil)%yt(kseg)*(-ll*smt(kseg,ll))/dlength
!!$      ddl1(ll + 3*NN + 4) = coil(icoil)%yt(kseg)*( ll*cmt(kseg,ll))/dlength
!!$      ddl1(ll + 4*NN + 5) = coil(icoil)%zt(kseg)*(-ll*smt(kseg,ll))/dlength
!!$      ddl1(ll + 5*NN + 6) = coil(icoil)%zt(kseg)*( ll*cmt(kseg,ll))/dlength
!!$
!!$     enddo
!!$     
!!$     do ll = 1, Cdof
!!$      d1A(icoil, ll) = d1A(icoil, ll) + ( dlength/arclength-one ) * ( ddl1(ll)/arclength - pi2*dlength*d1L(ll)/(coil(icoil)%L)**2 )
!!$     enddo
!!$
!!$    enddo ! end kseg
!!$
!!$   enddo ! end icoil
!!$
!!$   eqarc = eqarc * half * pi2 / NDcoil / Ncoils
!!$   d1A   = d1A          * pi2 / Ndcoil / Ncoils
!!$
!!$   t1A = d1A
!!$
!!$   DALLOCATE( d1A  )
!!$   DALLOCATE( d1L  )
!!$   DALLOCATE( ddl1 )
!!$   !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!!$
!!$  case(2)  !bugs exist; 0629/2016
!!$
!!$   SALLOCATE(d1A, (1:Ncoils, 0:Cdof), zero)
!!$   SALLOCATE(ddl1,(          0:Cdof), zero)
!!$
!!$   SALLOCATE(d2L, (0:Cdof  , 0:Cdof), zero)
!!$   SALLOCATE(d2A, (1:Ncoils, 0:Cdof, 1:Ncoils, 0:Cdof), zero)
!!$   SALLOCATE(ddl2,(0:Cdof  , 0:Cdof), zero)
!!$
!!$   if ( .not. allocated(t1A) ) then
!!$    SALLOCATE(t1A, (1:Ncoils, 0:Cdof), zero)
!!$    SALLOCATE(t2A, (1:Ncoils, 0:Cdof, 1:Ncoils, 0:Cdof), zero)
!!$   else
!!$    t1A = zero
!!$    t2A = zero
!!$   endif
!!$
!!$   do icoil = 1, Ncoils
!!$
!!$    arclength = coil(icoil)%L/pi2
!!$    call LenDeriv2( icoil, d2L(0:Cdof,0:Cdof) )
!!$
!!$    do kseg = 0, NDcoil-1
!!$
!!$     dlength = sqrt( coil(icoil)%xt(kseg)**2 + coil(icoil)%yt(kseg)**2 + coil(icoil)%zt(kseg)**2 )
!!$     eqarc = eqarc + ( (dlength - arclength) / arclength )**2
!!$
!!$     do ll = 0, NN 
!!$
!!$      ddl1(ll        + 1) = coil(icoil)%xt(kseg)*(-ll*smt(kseg,ll))/dlength
!!$      ddl1(ll +   NN + 2) = coil(icoil)%xt(kseg)*( ll*cmt(kseg,ll))/dlength
!!$      ddl1(ll + 2*NN + 3) = coil(icoil)%yt(kseg)*(-ll*smt(kseg,ll))/dlength
!!$      ddl1(ll + 3*NN + 4) = coil(icoil)%yt(kseg)*( ll*cmt(kseg,ll))/dlength
!!$      ddl1(ll + 4*NN + 5) = coil(icoil)%zt(kseg)*(-ll*smt(kseg,ll))/dlength
!!$      ddl1(ll + 5*NN + 6) = coil(icoil)%zt(kseg)*( ll*cmt(kseg,ll))/dlength
!!$
!!$     enddo
!!$
!!$     do ll = 1, Cdof
!!$      do mm = 1, Cdof
!!$       ddl2(ll,mm) = -ddl1(ll) * ddl1(mm) / dlength
!!$      enddo
!!$     enddo
!!$
!!$     do ll = 0, NN
!!$      do mm = 0, NN
!!$
!!$       ddl2(ll     +1,mm     +1) = ddl2(ll     +1,mm     +1) + (-mm*smt(kseg,mm))*(-ll*smt(kseg,ll))/dlength
!!$       ddl2(ll+  NN+2,mm     +1) = ddl2(ll+  NN+2,mm     +1) + (-mm*smt(kseg,mm))*( ll*cmt(kseg,ll))/dlength
!!$       ddl2(ll     +1,mm+  NN+2) = ddl2(ll     +1,mm+  NN+2) + ( mm*cmt(kseg,mm))*(-ll*smt(kseg,ll))/dlength
!!$       ddl2(ll+  NN+2,mm+  NN+2) = ddl2(ll+  NN+2,mm+  NN+2) + ( mm*cmt(kseg,mm))*( ll*cmt(kseg,ll))/dlength
!!$
!!$       ddl2(ll+2*NN+3,mm+2*NN+3) = ddl2(ll+2*NN+3,mm+2*NN+3) + (-mm*smt(kseg,mm))*(-ll*smt(kseg,ll))/dlength
!!$       ddl2(ll+3*NN+4,mm+2*NN+3) = ddl2(ll+3*NN+4,mm+2*NN+3) + (-mm*smt(kseg,mm))*( ll*cmt(kseg,ll))/dlength
!!$       ddl2(ll+2*NN+3,mm+3*NN+4) = ddl2(ll+2*NN+3,mm+3*NN+4) + ( mm*cmt(kseg,mm))*(-ll*smt(kseg,ll))/dlength
!!$       ddl2(ll+3*NN+4,mm+3*NN+4) = ddl2(ll+3*NN+4,mm+3*NN+4) + ( mm*cmt(kseg,mm))*( ll*cmt(kseg,ll))/dlength
!!$
!!$       ddl2(ll+4*NN+5,mm+4*NN+5) = ddl2(ll+4*NN+5,mm+4*NN+5) + (-mm*smt(kseg,mm))*(-ll*smt(kseg,ll))/dlength
!!$       ddl2(ll+5*NN+6,mm+4*NN+5) = ddl2(ll+5*NN+6,mm+4*NN+5) + (-mm*smt(kseg,mm))*( ll*cmt(kseg,ll))/dlength
!!$       ddl2(ll+4*NN+5,mm+5*NN+6) = ddl2(ll+4*NN+5,mm+5*NN+6) + ( mm*cmt(kseg,mm))*(-ll*smt(kseg,ll))/dlength
!!$       ddl2(ll+5*NN+6,mm+5*NN+6) = ddl2(ll+5*NN+6,mm+5*NN+6) + ( mm*cmt(kseg,mm))*( ll*cmt(kseg,ll))/dlength
!!$
!!$      enddo
!!$     enddo
!!$          
!!$     do ll = 1, Cdof
!!$      d1A(icoil, ll) = d1A(icoil, ll) + ( dlength/arclength-one ) * ( ddl1(ll)/arclength - pi2*dlength*d2L(ll,0)/(coil(icoil)%L)**2 )
!!$     enddo
!!$    
!!$     do ll = 1, Cdof
!!$      do mm = 1, Cdof
!!$
!!$       d2A(icoil,ll,icoil,mm) = d2A(icoil,ll,icoil,mm) + ( ddl1(ll)/arclength - pi2*dlength*d2L(ll,0)/(coil(icoil)%L)**2 ) * &
!!$                                                         ( ddl1(mm)/arclength - pi2*dlength*d2L(mm,0)/(coil(icoil)%L)**2 ) + &
!!$                    ( dlength/arclength-one ) * ( ddl2(ll,mm)/arclength - pi2*ddl1(ll)*d2L(mm,0)/(coil(icoil)%L)**2 - pi2*ddl1(mm)*d2L(ll,0)/(coil(icoil)%L)**2 - &
!!$                                                - pi2*dlength*d2L(ll,mm)/(coil(icoil)%L)**2 + 2*pi2*dlength*d2L(ll,0)*d2L(mm,0)/(coil(icoil)%L)**3 )
!!$      enddo
!!$     enddo
!!$
!!$
!!$
!!$
!!$
!!$
!!$    enddo ! end kseg
!!$
!!$   enddo ! end icoil
!!$
!!$   eqarc = eqarc * half * pi2 / NDcoil / Ncoils
!!$   d1A   = d1A          * pi2 / NDcoil / Ncoils
!!$   d2A   = d2A          * pi2 / NDcoil / Ncoils
!!$
!!$   t1A = d1A
!!$   t2A = d2A
!!$
!!$   DALLOCATE( d1A  )
!!$   DALLOCATE( ddl1 )
!!$   DALLOCATE( ddl2 )
!!$   DALLOCATE( d2A  )
!!$   DALLOCATE( d2L  )
!!$
!!$  end select
!!$end SUBROUTINE equarcl

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

SUBROUTINE equarcl( nderiv )
  use kmodule, only : zero, half, one, pi2, myid, ncpu, ounit, &
                      Cdof, Ncoils, NDcoil, NFcoil, coil, smt, cmt, &
                      eqarc, t1A, t2A
  implicit none
  include "mpif.h"


  INTEGER           :: nderiv

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  INTEGER           :: astat, ierr, icoil, kseg, ll, mm, NN, array2size, array4size
  REAL              :: arclength, dlength, leqarc
  REAL, allocatable :: d1A(:,:), d1L(:), ddl1(:), ddl2(:,:), d2A(:,:,:,:), d2L(:,:)

  NN      = NFcoil
  
  ! calculate each coil's length first
  do icoil = 1, Ncoils
   call LenDeriv0( icoil, coil(icoil)%L )
  enddo
  
  eqarc = zero; leqarc = zero

  select case ( nderiv )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  case(0)

   do icoil = 1, Ncoils

    if( myid.ne.modulo(icoil-1,ncpu) ) cycle ! parallelization loop; 

    arclength = coil(icoil)%L / pi2 ! L/pi2 * pi2/NDcoil

    do kseg = 0, NDcoil-1
     dlength = sqrt( coil(icoil)%xt(kseg)**2 + coil(icoil)%yt(kseg)**2 + coil(icoil)%zt(kseg)**2 )
     leqarc = leqarc + ( dlength - arclength )**2 
    enddo ! end kseg

   enddo ! end icoil

   call MPI_REDUCE( leqarc, eqarc, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
   RlBCAST(eqarc, 1, 0)

   eqarc = eqarc * half * pi2 / NDcoil / Ncoils

   !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  case(1)

   SALLOCATE(d1A, (1:Ncoils, 0:Cdof), zero)
   SALLOCATE(d1L, (          0:Cdof), zero)
   SALLOCATE(ddl1,(          0:Cdof), zero)

   if ( .not. allocated(t1A) ) then
    SALLOCATE(t1A, (1:Ncoils, 0:Cdof), zero)
   else
    t1A = zero
   endif

   do icoil = 1, Ncoils

    if( myid.ne.modulo(icoil-1,ncpu) ) cycle ! parallelization loop;

    arclength = coil(icoil)%L / pi2
    call LenDeriv1( icoil, d1L(0:Cdof) )

    do kseg = 0, NDcoil-1

     dlength = sqrt( coil(icoil)%xt(kseg)**2 + coil(icoil)%yt(kseg)**2 + coil(icoil)%zt(kseg)**2 )
     leqarc = leqarc + ( dlength - arclength )**2

     do ll = 0, NN 

      ddl1(ll        + 1) = coil(icoil)%xt(kseg)*(-ll*smt(kseg,ll))/dlength
      ddl1(ll +   NN + 2) = coil(icoil)%xt(kseg)*( ll*cmt(kseg,ll))/dlength
      ddl1(ll + 2*NN + 3) = coil(icoil)%yt(kseg)*(-ll*smt(kseg,ll))/dlength
      ddl1(ll + 3*NN + 4) = coil(icoil)%yt(kseg)*( ll*cmt(kseg,ll))/dlength
      ddl1(ll + 4*NN + 5) = coil(icoil)%zt(kseg)*(-ll*smt(kseg,ll))/dlength
      ddl1(ll + 5*NN + 6) = coil(icoil)%zt(kseg)*( ll*cmt(kseg,ll))/dlength

     enddo
     
     do ll = 1, Cdof
      d1A(icoil, ll) = d1A(icoil, ll) + ( dlength-arclength ) * ( ddl1(ll)- d1L(ll)/pi2 )
     enddo

    enddo ! end kseg

   enddo ! end icoil

   call MPI_REDUCE( leqarc, eqarc, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
   RlBCAST(eqarc, 1, 0)

   array2size = Ncoils * ( Cdof + 1 )
   call MPI_REDUCE( d1A, t1A, array2size, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
   RlBCAST( t1A, array2size, 0 )

   eqarc = eqarc * half * pi2 / NDcoil / Ncoils
   t1A   = t1A          * pi2 / Ndcoil / Ncoils


   DALLOCATE( d1A  )
   DALLOCATE( d1L  )
   DALLOCATE( ddl1 )
   !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  case(2)

   SALLOCATE(d1A, (1:Ncoils, 0:Cdof), zero)
   SALLOCATE(ddl1,(          0:Cdof), zero)

   SALLOCATE(d2L, (0:Cdof  , 0:Cdof), zero)
   SALLOCATE(d2A, (1:Ncoils, 0:Cdof, 1:Ncoils, 0:Cdof), zero)
   SALLOCATE(ddl2,(0:Cdof  , 0:Cdof), zero)

   if ( .not. allocated(t1A) ) then
    SALLOCATE(t1A, (1:Ncoils, 0:Cdof), zero)
    SALLOCATE(t2A, (1:Ncoils, 0:Cdof, 1:Ncoils, 0:Cdof), zero)
   else
    t1A = zero
    t2A = zero
   endif

   do icoil = 1, Ncoils

    if( myid.ne.modulo(icoil-1,ncpu) ) cycle ! parallelization loop;

    arclength = coil(icoil)%L / pi2
    call LenDeriv2( icoil, d2L(0:Cdof,0:Cdof) )

    do kseg = 0, NDcoil-1

     dlength = sqrt( coil(icoil)%xt(kseg)**2 + coil(icoil)%yt(kseg)**2 + coil(icoil)%zt(kseg)**2 )
     leqarc = leqarc + ( (dlength - arclength) / arclength )**2

     do ll = 0, NN 

      ddl1(ll        + 1) = coil(icoil)%xt(kseg)*(-ll*smt(kseg,ll))/dlength
      ddl1(ll +   NN + 2) = coil(icoil)%xt(kseg)*( ll*cmt(kseg,ll))/dlength
      ddl1(ll + 2*NN + 3) = coil(icoil)%yt(kseg)*(-ll*smt(kseg,ll))/dlength
      ddl1(ll + 3*NN + 4) = coil(icoil)%yt(kseg)*( ll*cmt(kseg,ll))/dlength
      ddl1(ll + 4*NN + 5) = coil(icoil)%zt(kseg)*(-ll*smt(kseg,ll))/dlength
      ddl1(ll + 5*NN + 6) = coil(icoil)%zt(kseg)*( ll*cmt(kseg,ll))/dlength

     enddo

     do ll = 1, Cdof
      do mm = 1, Cdof
       ddl2(ll,mm) = -ddl1(ll) * ddl1(mm) / dlength
      enddo
     enddo

     do ll = 0, NN
      do mm = 0, NN

       ddl2(ll     +1,mm     +1) = ddl2(ll     +1,mm     +1) + (-mm*smt(kseg,mm))*(-ll*smt(kseg,ll))/dlength
       ddl2(ll+  NN+2,mm     +1) = ddl2(ll+  NN+2,mm     +1) + (-mm*smt(kseg,mm))*( ll*cmt(kseg,ll))/dlength
       ddl2(ll     +1,mm+  NN+2) = ddl2(ll     +1,mm+  NN+2) + ( mm*cmt(kseg,mm))*(-ll*smt(kseg,ll))/dlength
       ddl2(ll+  NN+2,mm+  NN+2) = ddl2(ll+  NN+2,mm+  NN+2) + ( mm*cmt(kseg,mm))*( ll*cmt(kseg,ll))/dlength

       ddl2(ll+2*NN+3,mm+2*NN+3) = ddl2(ll+2*NN+3,mm+2*NN+3) + (-mm*smt(kseg,mm))*(-ll*smt(kseg,ll))/dlength
       ddl2(ll+3*NN+4,mm+2*NN+3) = ddl2(ll+3*NN+4,mm+2*NN+3) + (-mm*smt(kseg,mm))*( ll*cmt(kseg,ll))/dlength
       ddl2(ll+2*NN+3,mm+3*NN+4) = ddl2(ll+2*NN+3,mm+3*NN+4) + ( mm*cmt(kseg,mm))*(-ll*smt(kseg,ll))/dlength
       ddl2(ll+3*NN+4,mm+3*NN+4) = ddl2(ll+3*NN+4,mm+3*NN+4) + ( mm*cmt(kseg,mm))*( ll*cmt(kseg,ll))/dlength

       ddl2(ll+4*NN+5,mm+4*NN+5) = ddl2(ll+4*NN+5,mm+4*NN+5) + (-mm*smt(kseg,mm))*(-ll*smt(kseg,ll))/dlength
       ddl2(ll+5*NN+6,mm+4*NN+5) = ddl2(ll+5*NN+6,mm+4*NN+5) + (-mm*smt(kseg,mm))*( ll*cmt(kseg,ll))/dlength
       ddl2(ll+4*NN+5,mm+5*NN+6) = ddl2(ll+4*NN+5,mm+5*NN+6) + ( mm*cmt(kseg,mm))*(-ll*smt(kseg,ll))/dlength
       ddl2(ll+5*NN+6,mm+5*NN+6) = ddl2(ll+5*NN+6,mm+5*NN+6) + ( mm*cmt(kseg,mm))*( ll*cmt(kseg,ll))/dlength

      enddo
     enddo
          
     do ll = 1, Cdof
      d1A(icoil, ll) = d1A(icoil, ll) + ( dlength-arclength ) * ( ddl1(ll)- d2L(ll,0)/pi2 )
     enddo
    
     do ll = 1, Cdof
      do mm = 1, Cdof

       d2A(icoil,ll,icoil,mm) = d2A(icoil,ll,icoil,mm) + ( ddl1(ll)- d2L(ll,0)/pi2 ) * ( ddl1(mm)- d2L(mm,0)/pi2 ) + &
                                ( dlength-arclength ) * ( ddl2(ll,mm) - d2L(ll,mm)/pi2 )
      enddo
     enddo

    enddo ! end kseg

   enddo ! end icoil

   call MPI_REDUCE( leqarc, eqarc, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
   RlBCAST(eqarc, 1, 0)

   array2size = Ncoils * ( Cdof + 1 )
   call MPI_REDUCE( d1A, t1A, array2size, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
   RlBCAST( t1A, array2size, 0 )

   array4size = array2size * array2size
   call MPI_REDUCE( d2A, t2A, array4size, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
   RlBCAST( t2A, array4size, 0 )

   eqarc = eqarc * half * pi2 / NDcoil / Ncoils
   t1A   = t1A          * pi2 / NDcoil / Ncoils
   t2A   = t2A          * pi2 / NDcoil / Ncoils

   DALLOCATE( d1A  )
   DALLOCATE( ddl1 )
   DALLOCATE( ddl2 )
   DALLOCATE( d2A  )
   DALLOCATE( d2L  )

  end select
end SUBROUTINE equarcl
