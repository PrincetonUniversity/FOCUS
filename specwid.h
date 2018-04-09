!title (specwid) ! Spectral width condensation.

!latex \briefly{Condense the Fourier spectra of X, Y and Z to eliminate the nonuniqueness of the poloidal angel 
!latex          in the parametric representation of space curves. The spectral width can be represented as 
!latex          $M = \frac{\sum_{m=1,M}m^{p+q}({X_{c,m}}^2)+{X_{s,m}}^2)+{Y_{c,m}}^2)+{Y_{s,m}}^2)+{Z_{c,m}}^2)+{Z_{s,m}}^2)}
!latex          {\sum_{m=1,M}m^p({X_{c,m}}^2)+{X_{s,m}}^2)+{Y_{c,m}}^2)+{Y_{s,m}}^2)+{Z_{c,m}}^2)+{Z_{s,m}}^2)}$. 
!latex          We have mupltiple approaches for deriving the minimum of $M$.}

!latex \calledby{\link{denergy}}
!latex \tableofcontents

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \subsection{Introduction}
!latex \bi
!latex \item[1.] For $p=0, q=2$ and ignoring the normalizing denominator, the spectral width can be written as,
!latex $ M = \sum_{m=1,M}m^2({X_{c,m}}^2 +{X_{s,m}}^2 +{Y_{c,m}}^2 +{Y_{s,m}}^2 +{Z_{c,m}}^2 +{Z_{s,m}}^2) $ 
!latex (\emph{\color{red}Here, I'm not quite sure if this simplification is appropriate.}) And we do not want 
!latex this term changes our physics properties. So it's implemented in tangential variation.
!latex \be
!latex \delta x & = & \dot{x} \delta u \\
!latex \delta y & = & \dot{y} \delta u \\
!latex \delta z & = & \dot{z} \delta u \\
!latex \delta X_{c,m} & = & \int_0^{2\pi}\delta x \cos(m\theta) \; \mathrm{d} \theta \\
!latex \delta M & = & \sum_{m=1,M}2m^2\left ( X_{c,m} \delta X_{c,m} + X_{s,m} \delta X_{s,m} + Y_{c,m} \delta Y_{c,m}
!latex  + Y_{s,m} \delta Y_{s,m} + X_{c,m} \delta X_{c,m} + X_{s,m} \delta X_{s,m} \right ) \nonumber \\
!latex          & = & \sum_{m=1,M}2m^2\left ( X_{c,m} \int_0^{2\pi} \dot{x} \delta u\cos(m\theta) \; \mathrm{d} \theta + \ldots 
!latex                + Z_{s,m} \int_0^{2\pi} \dot{x} \delta u \cos(m\theta) \; \mathrm{d} \theta \right ) \nonumber \\
!latex          & = & 2 \int_0^{2\pi} \left ( \sum_{m=1,M} m^2 X_{c,m} \cos(m\theta) \dot{x} + \ldots + 
!latex                 \sum_{m=1,M} m^2 Z_{c,m} \cos(m\theta) \dot{z} \right ) \delta u \; \mathrm{d} \theta \nonumber \\
!latex   & = & -2 \int_0^{2\pi} \left ( \ddot{x} \dot{x} + \ldots + \ddot{z} \dot{z} \right ) \delta u \; \mathrm{d} \theta
!latex \ee
!latex \item[2.] Thus, if we want to minimize $ M $, make $\delta M = 0$, which equals making 
!latex $\ddot{x} \dot{x} + \ddot{y} \dot{y} + \ddot{z} \dot{z} = 0$ or $\dot{x}^2 + \dot{y}^2 + \dot{z}^2 = const.$ .
!latex That's so-called equal-arclength.
!latex \item[3.] In later constructions, we will use the equal-arclength constraint, rather than minimize spectral width directly. 
!latex \ei
!latex
!latex \subsection{Differential flow A}
!latex \bi
!latex \item[1.] For differential flow, we can just include the equal-arclength condition and to avoid the length term, we use 
!latex $\ddot{x} \dot{x} + \ddot{y} \dot{y} + \ddot{z} \dot{z} = 0$. So the equal-arclength constraint can be written
!latex as,
!latex \begin{equation}
!latex A = \frac{1}{2} \sum_i \int_0^{2\pi} ( \ddot{x} \dot{x} + \ddot{y} \dot{y} + \ddot{z} \dot{z} )^2 \; \mathrm{d} \theta
!latex \end{equation}
!latex \item[2.] The equal-arclength contraint's first derivatives can be written as,
!latex \be
!latex \ds \frac{\partial{A}}{\partial{x^i_m}} & = & \int_0^{2\pi} \left ( \frac{\partial{\ddot{x}}}{\partial{x^i_m}} \dot{x} + 
!latex                                                         \frac{\partial{\dot{x}}}{\partial{x^i_m}} \ddot{x} \right ) 
!latex                         \left ( \ddot{x} \dot{x} + \ddot{y} \dot{y} + \ddot{z} \dot{z} \right )\; \mathrm{d} \theta\\
!latex       x  & = & \sum_n      X_{c,n}\cos(n \theta) +     X_{s,n} \sin{n \theta} \\
!latex \dot{x}  & = & \sum_n -n   X_{c,n}\sin(n \theta) + n   X_{s,n} \cos{n \theta} \\
!latex \ddot{x} & = & \sum_n -n^2 X_{c,n}\cos(n \theta) - n^2 X_{s,n} \sin{n \theta}
!latex \ee
!latex \item[3.] And the second derivatives can be written as, (\emph{Actually, there is no need constructing the second derivatives
!latex            for differential flow method.})
!latex \be\begin{array}{cccccccc}
!latex \ds \frac{\partial^2{A}}{\partial{x^i_m} \partial{x^i_n}} & = & \int_0^{2\pi} & \left ( \frac{\partial{\ddot{x}}}{\partial{x^i_m}}
!latex     \frac{\partial{\dot{x}}}{\partial{x^i_n}} + \frac{\partial{\dot{x}}}{\partial{x^i_m}} \frac{\partial{\ddot{x}}}
!latex          {\partial{x^i_n}} \right ) \left ( \ddot{x} \dot{x} + \ddot{y} \dot{y} + \ddot{z} \dot{z} \right ) \, & + & \,   
!latex    \left  ( \frac{\partial{\ddot{x}}}{\partial{x^i_m}} \dot{x} + \frac{\partial{\dot{x}}}{\partial{x^i_m}} \ddot{x} \right ) 
!latex    \left  ( \frac{\partial{\ddot{x}}}{\partial{x^i_n}} \dot{x} + \frac{\partial{\dot{x}}}{\partial{x^i_n}} \ddot{x} \right ) 
!latex    \; \mathrm{d} \theta \\
!latex 
!latex \ds \frac{\partial^2{A}}{\partial{x^i_m} \partial{y^i_n}} & = & \int_0^{2\pi} &   &&
!latex    \left  ( \frac{\partial{\ddot{x}}}{\partial{x^i_m}} \dot{x} + \frac{\partial{\dot{x}}}{\partial{x^i_m}} \ddot{x} \right ) 
!latex    \left  ( \frac{\partial{\ddot{y}}}{\partial{y^i_n}} \dot{y} + \frac{\partial{\dot{y}}}{\partial{y^i_n}} \ddot{y} \right ) 
!latex    \; \mathrm{d} \theta 
!latex \end{array} \ee
!latex \ei
!latex 
!latex \subsection{Differential flow B}
!latex \bi
!latex \item[1.]In the above method, we are trying to minimize the line integral of equal arclength constraint, which will be zero  
!latex only when the integrand equals zero in the whole domain. That is hard to accomplish. Thus, similar as what we did with the  
!latex Langrange multipiler, we discretize the distribution function of the integrand using Fourier transformation and try to minimize  
!latex the square sum of the Fourier coefficients to zero. 
!latex \be
!latex \ds  f(\theta) & = & \ddot{x} \dot{x} + \ddot{y} \dot{y} + \ddot{z} \dot{z} \nonumber \\
!latex                & = & \sum_{n=0,N} f_{c,n} \cos(n \theta) + f_{s,n} \sin(n \theta)  \\
!latex \ds f_{c,n} & = & \frac{1}{\pi} \int_0^{2\pi} f(\theta) \cos(n \theta) \; \mathrm{d} \theta  \\
!latex \ds f_{s,n} & = & \frac{1}{\pi} \int_0^{2\pi} f(\theta) \sin(n \theta) \; \mathrm{d} \theta  \\
!latex A & = & \frac{1}{2} \sum_i \sum_n ({f_{c,n}}^2 + {f_{s,n}}^2 )
!latex \ee
!latex \item[2.] If the equal arclength constraint is written as above equation, then its derivatives can be derived as,
!latex \be
!latex \ds \frac{\partial{A}}{\partial{x^i_m}} & = & \sum_n \left( f_{c,n}^i \frac{\partial{f_{c,n}}}{\partial{x^i_m}} + 
!latex                                             f_{s,n}^i \frac{\partial{f_{s,n}}}{\partial{x^i_m}} \right ) \\
!latex \ds \frac{\partial^2{A}}{\partial{x^i_p} \partial{x^i_q}} & = & \sum_n \left( \frac{\partial{f_{c,n}}}{\partial{x^i_p}}
!latex     \frac{\partial{f_{c,n}}}{\partial{x^i_q}} + f_{c,n}^i \frac{\partial^2{f_{c,n}}}{\partial{x^i_p}\partial{x^i_q}} + 
!latex     \frac{\partial{f_{s,n}}}{\partial{x^i_p}} \frac{\partial{f_{s,n}}}{\partial{x^i_q}} + 
!latex     f_{s,n}^i \frac{\partial^2{f_{s,n}}}{\partial{x^i_p}\partial{x^i_q}} \right )
!latex \ee
!latex
!latex \ei
!latex 
!latex \subsection{Newton method}
!latex \bi
!latex \item[1.] For minimizing problems subjecting to some conditions, the most effective method is Langrange multipiler. 
!latex In this case, we can include :
!latex \be 
!latex \ds E_{total} & =  & E_{physics} + \sum_i \lambda_i(\theta) ( \ddot{x} \dot{x} + \ddot{y} \dot{y} + \ddot{z} \dot{z} )
!latex \ee
!latex Here, $\lambda_i(\theta)$ is a function of $\theta$. In the code, we do a FFT on the constraint function as
!latex \be
!latex \ds  f(\theta) & = & \ddot{x} \dot{x} + \ddot{y} \dot{y} + \ddot{z} \dot{z} \nonumber \\
!latex             & = & \sum_{n=0,N} f_{c,n} \cos(n \theta) + f_{s,n} \sin(n \theta)  \\
!latex \ds f_{c,n} & = & \frac{1}{\pi} \int_0^{2\pi} f(\theta) \cos(n \theta) \; \mathrm{d} \theta  \\
!latex \ds f_{s,n} & = & \frac{1}{\pi} \int_0^{2\pi} f(\theta) \sin(n \theta) \; \mathrm{d} \theta 
!latex \ee
!latex And the new Langrange multipiler is,
!latex \be
!latex \ds E_{total} & =  & E_{physics} + \sum_i \sum_n ( \lambda_{i,c} f_{c,n} + \lambda_{i,s} f_{s,n} )
!latex \ee
!latex Here,$i$ is denoted as coil label and $n$ for fourier modes number.
!latex \item[2.] First derivatives of the total ``energy" function are,
!latex \be
!latex \ds \frac{\partial{E_{total}}}{\partial{x^i_m}} & = & \frac{\partial{E_{physics}}}{\partial{x^i_m}} +  \sum_n \left (
!latex     \lambda_{i,c} \frac{\partial{f^i_{c,n}}}{\partial{x^i_m}} + \lambda_{i,s} \frac{\partial{f^i_{s,n}}}{\partial{x^i_m}} \right ) \\
!latex \ds \frac{\partial{E_{total}}}{\partial{\lambda_{i,c}}} & = &  f^i_{c,n}  \\
!latex \ds \frac{\partial{E_{total}}}{\partial{\lambda_{i,s}}} & = &  f^i_{s,n}  
!latex \ee
!latex And we have
!latex \be
!latex \ds \frac{\partial{f_{c,n}}}{\partial{x^i_m}} & = & \frac{1}{\pi} \int_0^{2\pi} \frac{\partial{f(\theta)}}{\partial{x^i_m}} 
!latex     \cos(n \theta) \; \mathrm{d} \theta  \\
!latex \ds \frac{\partial{f_{s,n}}}{\partial{x^i_m}} & = & \frac{1}{\pi} \int_0^{2\pi} \frac{\partial{f(\theta)}}{\partial{x^i_m}} 
!latex     \sin(n \theta) \; \mathrm{d} \theta  \\
!latex \frac{\partial{f(\theta)}}{\partial{x^i_m}} & = & \frac{\partial{\ddot{x}}}{\partial{x^i_m}} \dot{x} + 
!latex                                                   \frac{\partial{\dot{x}}}{\partial{x^i_m}} \ddot{x}
!latex \ee
!latex \item[3.] The second derivatives are written as,
!latex \be
!latex \ds \frac{\partial^2{E_{total}}}{\partial{\lambda_{i,c}} \partial{x^i_m}} & = & \frac{\partial{f^i_{c,n}}}{\partial{x^i_m}} \\
!latex \ds \frac{\partial^2{E_{total}}}{\partial{x^i_m} \partial{x^i_n}} & = & \frac{\partial^2{E_{phsyics}}}{ \partial{x^i_m} 
!latex     \partial{x^i_n}} +  \sum_n \left ( \lambda_{i,c} \frac{\partial^2{f^i_{c,n}}}{\partial{x^i_m} \partial{x^i_n}} + 
!latex     \lambda_{i,s} \frac{\partial^2{f^i_{s,n}}}{\partial{x^i_m} \partial{x^i_n}} \right )
!latex \ee
!latex And we have
!latex \be
!latex \ds \frac{\partial^2{f_{c,n}}}{\partial{x^i_m}\partial{x^i_n}} & = & \frac{1}{\pi} \int_0^{2\pi} \frac{\partial^2{f(\theta)}}
!latex  {\partial{x^i_m} \partial{x^i_n}} \cos(n \theta) \; \mathrm{d} \theta  \\
!latex \ds \frac{\partial^2{f_{s,n}}}{\partial{x^i_m}\partial{x^i_n}} & = & \frac{1}{\pi} \int_0^{2\pi} \frac{\partial^2{f(\theta)}}
!latex  {\partial{x^i_m} \partial{x^i_n}} \sin(n \theta) \; \mathrm{d} \theta  \\
!latex \ds \frac{\partial^2{f(\theta)}}{\partial{x^i_m}\partial{x^i_n}} & = & \frac{\partial{\ddot{x}}}{\partial{x^i_m}}
!latex   \frac{\partial{\dot{x}}}{\partial{x^i_n}} + \frac{\partial{\dot{x}}}{\partial{x^i_m}} \frac{\partial{\ddot{x}}} {\partial{x^i_n}} 
!latex \ee
!latex \ei

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
SUBROUTINE specwid(nderiv)
  use kmodule, only: zero, half, myid, ncpu, ounit, &
                     coil, smt, cmt, Ncoils, Cdof, NDcoil, NFcoil, &
                     eqarc, t1A, t2A

  implicit none
  include "mpif.h"

  INTEGER           :: nderiv
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  INTEGER           :: astat, ierr, icoil, kseg, ll, mm, nd, nf, array2size, array4size
  REAL              :: leqarc, normalfactor
  REAL, allocatable :: fc(:), fs(:), fA(:), f1A(:,:), f2A(:,:,:), l1A(:,:), l2A(:,:,:,:), ffc(:,:), ffs(:,:), &
                       fcxx(:,:,:), fsxx(:,:,:), df1(:), df2(:), lfc(:), lfs(:)

  !initialization
  !allocation
  
  array2size = Ncoils * ( Cdof + 1 )
  array4size = array2size * array2size

  nd = NDcoil-1; nf = NFcoil
  eqarc = zero; leqarc = zero
! normalfactor = 2.0*Ncoils*(nf+1)
  normalfactor = 1.0

  SALLOCATE( fA, (0:nd), zero )
  SALLOCATE( fc, (0:nf), zero )
  SALLOCATE( fs, (0:nf), zero )       

  select case(nderiv)
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  case(0)

     do icoil = 1, Ncoils

        if( myid.ne.modulo(icoil-1,ncpu) ) cycle ! parallelization loop;
        
        do kseg = 0, nd

           fA(kseg) = coil(icoil)%xa(kseg)*coil(icoil)%xt(kseg) + coil(icoil)%ya(kseg)*coil(icoil)%yt(kseg) + &
                      coil(icoil)%za(kseg)*coil(icoil)%zt(kseg)
        enddo ! end kseg

        call coildft( fA, fc, fs, nd, nf)

        leqarc = leqarc + sum(fc(0:nf)**2) + sum(fs(0:nf)**2) !+ sum(coil(icoil)%zc(0:NFcoil))**2

     enddo !end icoil
     
     call MPI_REDUCE( leqarc, eqarc, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
     RlBCAST(eqarc, 1, 0)
     eqarc = eqarc * half / normalfactor

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  case(1)

     SALLOCATE(f1A, (0:nd,0:Cdof), zero)
     SALLOCATE(l1A, (1:Ncoils, 0: Cdof), zero)
     SALLOCATE(lfc, (0: nf), zero)
     SALLOCATE(lfs, (0: nf), zero)

     if ( .not. allocated(t1A) ) then
        SALLOCATE(t1A, (1:Ncoils, 0:Cdof), zero)
     else
        t1A = zero
     endif
          
     do icoil = 1, Ncoils

        if( myid.ne.modulo(icoil-1,ncpu) ) cycle ! parallelization loop;

        do kseg = 0, nd

           fA(kseg) = coil(icoil)%xa(kseg)*coil(icoil)%xt(kseg) + coil(icoil)%ya(kseg)*coil(icoil)%yt(kseg) + &
                      coil(icoil)%za(kseg)*coil(icoil)%zt(kseg)
        enddo ! end kseg

        call coildft( fA, fc, fs, nd, nf)

        lfc = fc; lfs = fs

        leqarc = leqarc + sum(lfc(0:nf)**2) + sum(lfs(0:nf)**2) !+ sum(coil(icoil)%zc(0:NFcoil))**2

        do ll = 0, nf

           do kseg = 0, nd

              f1A(kseg,ll     +1) =  -ll**2*cmt(kseg,ll)*coil(icoil)%xt(kseg) - ll*smt(kseg,ll)*coil(icoil)%xa(kseg)
              f1A(kseg,ll+  nf+2) =  -ll**2*smt(kseg,ll)*coil(icoil)%xt(kseg) + ll*cmt(kseg,ll)*coil(icoil)%xa(kseg)
              f1A(kseg,ll+2*nf+3) =  -ll**2*cmt(kseg,ll)*coil(icoil)%yt(kseg) - ll*smt(kseg,ll)*coil(icoil)%ya(kseg)
              f1A(kseg,ll+3*nf+4) =  -ll**2*smt(kseg,ll)*coil(icoil)%yt(kseg) + ll*cmt(kseg,ll)*coil(icoil)%ya(kseg)
              f1A(kseg,ll+4*nf+5) =  -ll**2*cmt(kseg,ll)*coil(icoil)%zt(kseg) - ll*smt(kseg,ll)*coil(icoil)%za(kseg)
              f1A(kseg,ll+5*nf+6) =  -ll**2*smt(kseg,ll)*coil(icoil)%zt(kseg) + ll*cmt(kseg,ll)*coil(icoil)%za(kseg)

           enddo ! end kseg

        enddo ! end ll

        do ll = 1, Cdof
          
           call coildft( f1A(0:Nd,ll), fc, fs, nd, nf)

           l1A(icoil, ll) = dot_product(lfc,fc) + dot_product(lfs,fs)

        enddo ! end ll

!!$        do ll = 4*NFcoil+4, 5*NFcoil+5
!!$           l1A(icoil, ll) = l1A(icoil,ll) + sum(coil(icoil)%zc(0:NFcoil))   !z(0) = 0 terms
!!$        enddo

     enddo ! end icoil
     
     call MPI_REDUCE( leqarc, eqarc, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
     RlBCAST(eqarc, 1, 0)
     eqarc = eqarc * half / normalfactor

     call MPI_REDUCE( l1A, t1A, array2size, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
     RlBCAST( t1A, array2size, 0 )
     t1A = t1A  / normalfactor

     DALLOCATE( f1A )
     DALLOCATE( l1A )
     DALLOCATE( lfc )
     DALLOCATE( lfs )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  case(2)

     SALLOCATE(f1A, (0:nd,0:Cdof), zero)
     SALLOCATE(l1A, (1:Ncoils, 0: Cdof), zero)
     SALLOCATE(ffc, (0:nf, 0:Cdof), zero)
     SALLOCATE(ffs, (0:nf, 0:Cdof), zero)
     SALLOCATE(df1, (0:Cdof), zero)  ! df1(n) = \partial(\dot{x}) / \partial{x_n}
     SALLOCATE(df2, (0:Cdof), zero)  ! df2(n) = \partial(\ddot{x}) / \partial{x_n}
     SALLOCATE(fcxx, (0:nf, 0:Cdof, 0:Cdof), zero )
     SALLOCATE(fsxx, (0:nf, 0:Cdof, 0:Cdof), zero )
     SALLOCATE(f2A, (0:nd, 0:Cdof, 0:Cdof), zero )
     SALLOCATE(l2A, (1:Ncoils, 0:Cdof, 1:Ncoils, 0:Cdof), zero)

     if ( .not. allocated(t1A) ) then
        SALLOCATE(t1A, (1:Ncoils, 0:Cdof), zero)
     else
        t1A = zero
     endif
     
     if ( .not. allocated(t2A) ) then
        SALLOCATE(t2A, (1:Ncoils, 0:Cdof, 1:Ncoils, 0:Cdof), zero)
     else
        t2A = zero
     endif

     do kseg = 0, nd
        
        do ll = 0, nf
           df1(ll     +1) = -ll*smt(kseg,ll); df2(ll     +1) = -ll**2*cmt(kseg,ll)
           df1(ll+  nf+2) =  ll*cmt(kseg,ll); df2(ll+  nf+2) = -ll**2*smt(kseg,ll)
           df1(ll+2*nf+3) = -ll*smt(kseg,ll); df2(ll+2*nf+3) = -ll**2*cmt(kseg,ll)
           df1(ll+3*nf+4) =  ll*cmt(kseg,ll); df2(ll+3*nf+4) = -ll**2*smt(kseg,ll)
           df1(ll+4*nf+5) = -ll*smt(kseg,ll); df2(ll+4*nf+5) = -ll**2*cmt(kseg,ll)
           df1(ll+5*nf+6) =  ll*cmt(kseg,ll); df2(ll+5*nf+6) = -ll**2*smt(kseg,ll)
        enddo !end ll

        do ll = 1, 2*nf+2
           do mm = 1, 2*nf+2
              f2A(kseg,ll       ,mm       ) = df1(ll       )*df2(mm       ) + df1(mm       )*df2(ll       )
              f2A(kseg,ll+2*nf+2,mm+2*nf+2) = df1(ll+2*nf+2)*df2(mm+2*nf+2) + df1(mm+2*nf+2)*df2(ll+2*nf+2)
              f2A(kseg,ll+4*nf+4,mm+4*nf+4) = df1(ll+4*nf+4)*df2(mm+4*nf+4) + df1(mm+4*nf+4)*df2(ll+4*nf+4)
           enddo !end mm
        enddo !end ll

     enddo !end kseg

     do ll = 1, Cdof
        do mm = 1, Cdof
           call coildft(f2A(0:nd, ll, mm), fcxx(0:nf,ll,mm), fsxx(0:nf,ll,mm), nd, nf)
        enddo
     enddo
          
     do icoil = 1, Ncoils

        if( myid.ne.modulo(icoil-1,ncpu) ) cycle ! parallelization loop;

        do kseg = 0, nd

           fA(kseg) = coil(icoil)%xa(kseg)*coil(icoil)%xt(kseg) + coil(icoil)%ya(kseg)*coil(icoil)%yt(kseg) + &
                      coil(icoil)%za(kseg)*coil(icoil)%zt(kseg)
        enddo ! end kseg

        call coildft( fA, fc, fs, nd, nf)

        leqarc = leqarc +  sum(fc(0:nf)**2) + sum(fs(0:nf)**2) !+ sum(coil(icoil)%zc(0:NFcoil))**2

        ffc(0:nf,0) = fc(0:nf)
        ffs(0:nf,0) = fs(0:nf)

        do ll = 0, nf

           do kseg = 0, nd

              f1A(kseg,ll     +1) =  -ll**2*cmt(kseg,ll)*coil(icoil)%xt(kseg) - ll*smt(kseg,ll)*coil(icoil)%xa(kseg)
              f1A(kseg,ll+  nf+2) =  -ll**2*smt(kseg,ll)*coil(icoil)%xt(kseg) + ll*cmt(kseg,ll)*coil(icoil)%xa(kseg)
              f1A(kseg,ll+2*nf+3) =  -ll**2*cmt(kseg,ll)*coil(icoil)%yt(kseg) - ll*smt(kseg,ll)*coil(icoil)%ya(kseg)
              f1A(kseg,ll+3*nf+4) =  -ll**2*smt(kseg,ll)*coil(icoil)%yt(kseg) + ll*cmt(kseg,ll)*coil(icoil)%ya(kseg)
              f1A(kseg,ll+4*nf+5) =  -ll**2*cmt(kseg,ll)*coil(icoil)%zt(kseg) - ll*smt(kseg,ll)*coil(icoil)%za(kseg)
              f1A(kseg,ll+5*nf+6) =  -ll**2*smt(kseg,ll)*coil(icoil)%zt(kseg) + ll*cmt(kseg,ll)*coil(icoil)%za(kseg)

           enddo ! end kseg

        enddo ! end ll

        do ll = 1, Cdof
          
           call coildft( f1A(0:Nd,ll), fc, fs, nd, nf)

           ffc(0:nf,ll) = fc(0:nf)
           ffs(0:nf,ll) = fs(0:nf)

           l1A(icoil, ll) =  dot_product(ffc(0:nf,0),fc) + dot_product(ffs(0:nf,0),fs)

        enddo ! end ll

!!$        do ll = 4*NFcoil+4, 5*NFcoil+5
!!$           l1A(icoil, ll) = l1A(icoil,ll) + sum(coil(icoil)%zc(0:NFcoil))   !z(0) = 0 terms
!!$        enddo

        do ll = 1, Cdof          
           do mm = 1, Cdof
              l2A(icoil, ll, icoil, mm) = dot_product(ffc(0:nf,ll), ffc(0:nf,mm)) + dot_product(ffc(0:nf,0), fcxx(0:nf,ll,mm)) &
                                        + dot_product(ffs(0:nf,ll), ffs(0:nf,mm)) + dot_product(ffs(0:nf,0), fsxx(0:nf,ll,mm))
           enddo  !end mm
        enddo !end ll

!!$        do ll = 4*NFcoil+5, 5*NFcoil+5
!!$           do mm = 4*NFcoil+5, 5*NFcoil+5
!!$              l2A(icoil,ll, icoil, mm) = l2A(icoil,ll, icoil, mm) + 1
!!$           enddo
!!$        enddo
                                   
     enddo ! end icoil

     call MPI_REDUCE( leqarc, eqarc, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
     RlBCAST(eqarc, 1, 0)
     eqarc = eqarc * half / normalfactor

     call MPI_REDUCE( l1A, t1A, array2size, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
     RlBCAST( t1A, array2size, 0 )
     t1A = t1A / normalfactor

     call MPI_REDUCE( l2A, t2A, array4size, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
     RlBCAST( t2A, array4size, 0 )
     t2A = t2A / normalfactor

     DALLOCATE( f1A )
     DALLOCATE( l1A )
     DALLOCATE( ffc )
     DALLOCATE( ffs )
     DALLOCATE( df1 )
     DALLOCATE( df2 )
     DALLOCATE( fcxx)
     DALLOCATE( fsxx)
     DALLOCATE( f2A )
     DALLOCATE( l2A )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  end select

  DALLOCATE( fA )
  DALLOCATE( fc )
  DALLOCATE( fs )

  return

END SUBROUTINE specwid

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine specwid_df(nderiv)
! spectral width for differential flow method; 07/19/2016
  use kmodule, only: zero, half, myid, ounit, ncpu, &
                     Cdof, Ncoils, NFcoil, coil, eqarc, t1A
  implicit none
  include "mpif.h"

  INTEGER, intent(in) :: nderiv

  INTEGER             :: astat, ierr, icoil, NN, array2size, array4size
  REAL                :: deqarc, leqarc
  REAl, allocatable   :: l1A(:,:)
  !--------------------------------------------------------

  eqarc = zero; leqarc = zero; NN = NFcoil
  array2size = Ncoils * ( Cdof + 1 ); array4size = array2size * array2size

  select case (nderiv)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  case(0)

     do icoil = 1, Ncoils

        if( myid.ne.modulo(icoil-1,ncpu) ) cycle ! parallelization loop;

        call dfeqarc0(icoil, deqarc)

        leqarc = leqarc + deqarc

     enddo !end icoil

     call MPI_REDUCE( leqarc, eqarc, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
     RlBCAST(eqarc, 1, 0)

     eqarc = eqarc * half / Ncoils

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  case(1)

     if ( .not. allocated(t1A) ) then
        SALLOCATE(t1A, (1:Ncoils, 0:Cdof), zero)
     else
        t1A = zero
     endif

     SALLOCATE(l1A, (1:Ncoils, 0:Cdof), zero)

     do icoil = 1, Ncoils

        if( myid.ne.modulo(icoil-1,ncpu) ) cycle ! parallelization loop;

        call dfeqarc1( icoil, l1A(icoil, 0:Cdof) )

        leqarc = leqarc + l1A(icoil,0)

     enddo !end icoil
     
     call MPI_REDUCE( leqarc, eqarc, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
     RlBCAST(eqarc, 1, 0)

     call MPI_REDUCE( l1A, t1A, array2size, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
     RlBCAST( t1A, array2size, 0 )

     eqarc = eqarc * half  / Ncoils
     t1A   = t1A           / Ncoils

     DALLOCATE( l1A  )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  case default

     FATAL( specwid, .true., selected nderiv is not supported )

  end select
     
  return

end subroutine specwid_df

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!!$subroutine specwid_nt(nderiv)
!!$! spectral width for newton method; 07/19/2016
!!$
!!$end subroutine specwid_nt

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine dfeqarc0(icoil, derivs)

!calculate the line integral of differential arc on icoil
  use kmodule, only: zero, pi2, coil, NDcoil, Ncoils, &
                     myid, ncpu, ounit
  implicit none
  include "mpif.h"

  INTEGER, intent(in )  :: icoil
  REAL   , intent(out)  :: derivs

  INTEGER :: astat, ierr, kseg

  !-----------------------------------------------------------
  FATAL( dfeqarc0, icoil .lt. 1 .or. icoil .gt. Ncoils, icoil not in right range )
  
  derivs = zero

  do kseg = 0, NDcoil-1

     derivs = derivs + ( coil(icoil)%xa(kseg)*coil(icoil)%xt(kseg) + coil(icoil)%ya(kseg)*coil(icoil)%yt(kseg) + &
                         coil(icoil)%za(kseg)*coil(icoil)%zt(kseg) )**2
     
  enddo
  
  derivs = derivs * pi2/NDcoil

  return

end subroutine dfeqarc0

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine dfeqarc1(icoil, derivs)

!calculate the line integral of differential arc and its first dirivatives on icoil
  use kmodule, only: zero, pi2, Cdof, coil, smt, cmt, NFcoil, NDcoil, Ncoils, &
                     myid, ncpu, ounit
  implicit none
  include "mpif.h"

  INTEGER, intent(in )  :: icoil
  REAL   , intent(out)  :: derivs(0:Cdof)

  INTEGER               :: astat, ierr, kseg, ll, NN
  REAl                  :: dlength

  !-----------------------------------------------------------
  FATAL( dfeqarc1, icoil .lt. 1 .or. icoil .gt. Ncoils, icoil not in right range )
  
  derivs = zero; NN = NFcoil

  do kseg = 0, NDcoil-1

     dlength = coil(icoil)%xa(kseg)*coil(icoil)%xt(kseg) + coil(icoil)%ya(kseg)*coil(icoil)%yt(kseg) + &
               coil(icoil)%za(kseg)*coil(icoil)%zt(kseg)

     derivs(0) = derivs(0) + dlength**2

     do ll = 0, NN
        derivs(ll        + 1) = derivs(ll        + 1) + ( -ll**2*cmt(kseg,ll)*coil(icoil)%xt(kseg) - ll*smt(kseg,ll)*coil(icoil)%xa(kseg) ) * dlength
        derivs(ll +   NN + 2) = derivs(ll +   NN + 2) + ( -ll**2*smt(kseg,ll)*coil(icoil)%xt(kseg) + ll*cmt(kseg,ll)*coil(icoil)%xa(kseg) ) * dlength
        derivs(ll + 2*NN + 3) = derivs(ll + 2*NN + 3) + ( -ll**2*cmt(kseg,ll)*coil(icoil)%yt(kseg) - ll*smt(kseg,ll)*coil(icoil)%ya(kseg) ) * dlength
        derivs(ll + 3*NN + 4) = derivs(ll + 3*NN + 4) + ( -ll**2*smt(kseg,ll)*coil(icoil)%yt(kseg) + ll*cmt(kseg,ll)*coil(icoil)%ya(kseg) ) * dlength
        derivs(ll + 4*NN + 5) = derivs(ll + 4*NN + 5) + ( -ll**2*cmt(kseg,ll)*coil(icoil)%zt(kseg) - ll*smt(kseg,ll)*coil(icoil)%za(kseg) ) * dlength
        derivs(ll + 5*NN + 6) = derivs(ll + 5*NN + 6) + ( -ll**2*smt(kseg,ll)*coil(icoil)%zt(kseg) + ll*cmt(kseg,ll)*coil(icoil)%za(kseg) ) * dlength
     enddo ! end ll  
     
  enddo !end kseg
  
  derivs = derivs * pi2/NDcoil

  return

end subroutine dfeqarc1

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine langrange(nderiv)
  use kmodule, only: zero, half, myid, ncpu, ounit, &
                     coil, smt, cmt, Ncoils, Cdof, NDcoil, NFcoil, &
                     eqarc, t1A, t2A, dlc, dls

  implicit none
  include "mpif.h"

  INTEGER           :: nderiv
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  INTEGER           :: astat, ierr, icoil, kseg, ll, mm, nd, nf, array2size, array4size, lmdsize
  REAL              :: leqarc
  REAL, allocatable :: fc(:), fs(:), fA(:), f1A(:,:), f2A(:,:,:), l1A(:,:), l2A(:,:,:,:), llc(:,:,:), lls(:,:,:), &
                       fcxx(:,:,:), fsxx(:,:,:), df1(:), df2(:)

  !initialization
  !allocation
  
  array2size = Ncoils * ( Cdof + 1 )
  array4size = array2size * array2size
  lmdsize = Ncoils*(NFcoil+1)*(Cdof+1)
  nd = NDcoil-1; nf = NFcoil
  eqarc = zero; leqarc = zero

  SALLOCATE( fA, (0:nd), zero )
  SALLOCATE( fc, (0:nf), zero )
  SALLOCATE( fs, (0:nf), zero )       

  select case(nderiv)
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  case(0)

     do icoil = 1, Ncoils

        if( myid.ne.modulo(icoil-1,ncpu) ) cycle ! parallelization loop;
        
        do kseg = 0, nd

           fA(kseg) = coil(icoil)%xa(kseg)*coil(icoil)%xt(kseg) + coil(icoil)%ya(kseg)*coil(icoil)%yt(kseg) + &
                      coil(icoil)%za(kseg)*coil(icoil)%zt(kseg)
        enddo ! end kseg

        call coildft( fA, fc, fs, nd, nf)

        leqarc = leqarc + dot_product(coil(icoil)%lmdc,fc) + dot_product(coil(icoil)%lmds,fs)

     enddo !end icoil
     
     call MPI_REDUCE( leqarc, eqarc, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
     RlBCAST(eqarc, 1, 0)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  case(1)

     SALLOCATE(f1A, (0:nd,0:Cdof), zero)
     SALLOCATE(l1A, (1:Ncoils, 0: Cdof), zero)
     SALLOCATE(llc, (1:Ncoils, 0: nf, 0:Cdof), zero)
     SALLOCATE(lls, (1:Ncoils, 0: nf, 0:Cdof), zero)

     if ( .not. allocated(t1A) ) then
        SALLOCATE(t1A, (1:Ncoils, 0:Cdof), zero)
     else
        t1A = zero
     endif
     
     if ( .not. allocated(dlc) ) then
        SALLOCATE(dlc, (1:Ncoils, 0:nf, 0:Cdof), zero)
        SALLOCATE(dls, (1:Ncoils, 0:nf, 0:Cdof), zero)
     else
        dlc = zero
        dls = zero
     endif
     
     do icoil = 1, Ncoils

        if( myid.ne.modulo(icoil-1,ncpu) ) cycle ! parallelization loop;

        do kseg = 0, nd

           fA(kseg) = coil(icoil)%xa(kseg)*coil(icoil)%xt(kseg) + coil(icoil)%ya(kseg)*coil(icoil)%yt(kseg) + &
                      coil(icoil)%za(kseg)*coil(icoil)%zt(kseg)
        enddo ! end kseg

        call coildft( fA, fc, fs, nd, nf)

        llc(icoil, 0:nf, 0) = fc(0:nf)
        lls(icoil, 0:nf, 0) = fs(0:nf)

        leqarc = leqarc + dot_product(coil(icoil)%lmdc,fc) + dot_product(coil(icoil)%lmds,fs)

        do ll = 0, nf

           do kseg = 0, nd

              f1A(kseg,ll     +1) =  -ll**2*cmt(kseg,ll)*coil(icoil)%xt(kseg) - ll*smt(kseg,ll)*coil(icoil)%xa(kseg)
              f1A(kseg,ll+  nf+2) =  -ll**2*smt(kseg,ll)*coil(icoil)%xt(kseg) + ll*cmt(kseg,ll)*coil(icoil)%xa(kseg)
              f1A(kseg,ll+2*nf+3) =  -ll**2*cmt(kseg,ll)*coil(icoil)%yt(kseg) - ll*smt(kseg,ll)*coil(icoil)%ya(kseg)
              f1A(kseg,ll+3*nf+4) =  -ll**2*smt(kseg,ll)*coil(icoil)%yt(kseg) + ll*cmt(kseg,ll)*coil(icoil)%ya(kseg)
              f1A(kseg,ll+4*nf+5) =  -ll**2*cmt(kseg,ll)*coil(icoil)%zt(kseg) - ll*smt(kseg,ll)*coil(icoil)%za(kseg)
              f1A(kseg,ll+5*nf+6) =  -ll**2*smt(kseg,ll)*coil(icoil)%zt(kseg) + ll*cmt(kseg,ll)*coil(icoil)%za(kseg)

           enddo ! end kseg

        enddo ! end ll

        do ll = 1, Cdof
          
           call coildft( f1A(0:Nd,ll), fc, fs, nd, nf)

           l1A(icoil, ll) = dot_product(coil(icoil)%lmdc,fc) + dot_product(coil(icoil)%lmds,fs)

        enddo ! end ll

     enddo ! end icoil
     
     call MPI_REDUCE( leqarc, eqarc, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
     RlBCAST(eqarc, 1, 0)

     call MPI_REDUCE( l1A, t1A, array2size, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
     RlBCAST( t1A, array2size, 0 )

     call MPI_REDUCE( llc, dlc, lmdsize, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr ) 
     RlBCAST( dlc, lmdsize, 0 )

     call MPI_REDUCE( lls, dls, lmdsize, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr ) 
     RlBCAST( dls, lmdsize, 0 )
    
     DALLOCATE( f1A )
     DALLOCATE( l1A )
     DALLOCATE( llc )
     DALLOCATE( lls )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  case(2)

     SALLOCATE(f1A, (0:nd,0:Cdof), zero)
     SALLOCATE(l1A, (1:Ncoils, 0: Cdof), zero)
     SALLOCATE(llc, (1:Ncoils, 0: nf, 0:Cdof), zero)
     SALLOCATE(lls, (1:Ncoils, 0: nf, 0:Cdof), zero)

     SALLOCATE(df1, (0:Cdof), zero)  ! df1(n) = \partial(\dot{x}) / \partial{x_n}
     SALLOCATE(df2, (0:Cdof), zero)  ! df2(n) = \partial(\ddot{x}) / \partial{x_n}
     SALLOCATE(f2A, (0:nd, 0:Cdof, 0:Cdof), zero )
     SALLOCATE(fcxx, (0:nf, 0:Cdof, 0:Cdof), zero )
     SALLOCATE(fsxx, (0:nf, 0:Cdof, 0:Cdof), zero )
     SALLOCATE(l2A, (1:Ncoils, 0:Cdof, 1:Ncoils, 0:Cdof), zero)

     if ( .not. allocated(t1A) ) then
        SALLOCATE(t1A, (1:Ncoils, 0:Cdof), zero)
     else
        t1A = zero
     endif
     
     if ( .not. allocated(dlc) ) then
        SALLOCATE(dlc, (1:Ncoils, 0:nf, 0:Cdof), zero)
        SALLOCATE(dls, (1:Ncoils, 0:nf, 0:Cdof), zero)
     else
        dlc = zero
        dls = zero
     endif

     if ( .not. allocated(t2A) ) then
        SALLOCATE(t2A, (1:Ncoils, 0:Cdof, 1:Ncoils, 0:Cdof), zero)
     else
        t2A = zero
     endif

     do kseg = 0, nd
        
        do ll = 0, nf
           df1(ll     +1) = -ll*smt(kseg,ll); df2(ll     +1) = -ll**2*cmt(kseg,ll)
           df1(ll+  nf+2) =  ll*cmt(kseg,ll); df2(ll+  nf+2) = -ll**2*smt(kseg,ll)
           df1(ll+2*nf+3) = -ll*smt(kseg,ll); df2(ll+2*nf+3) = -ll**2*cmt(kseg,ll)
           df1(ll+3*nf+4) =  ll*cmt(kseg,ll); df2(ll+3*nf+4) = -ll**2*smt(kseg,ll)
           df1(ll+4*nf+5) = -ll*smt(kseg,ll); df2(ll+4*nf+5) = -ll**2*cmt(kseg,ll)
           df1(ll+5*nf+6) =  ll*cmt(kseg,ll); df2(ll+5*nf+6) = -ll**2*smt(kseg,ll)
        enddo !end ll

        do ll = 1, 2*nf+2
           do mm = 1, 2*nf+2
              f2A(kseg,ll       ,mm       ) = df1(ll       )*df2(mm       ) + df1(mm       )*df2(ll       )
              f2A(kseg,ll+2*nf+2,mm+2*nf+2) = df1(ll+2*nf+2)*df2(mm+2*nf+2) + df1(mm+2*nf+2)*df2(ll+2*nf+2)
              f2A(kseg,ll+4*nf+4,mm+4*nf+4) = df1(ll+4*nf+4)*df2(mm+4*nf+4) + df1(mm+4*nf+4)*df2(ll+4*nf+4)
           enddo !end mm
        enddo !end ll

     enddo !end kseg

     do ll = 1, Cdof
        do mm = 1, Cdof
           call coildft(f2A(0:nd, ll, mm), fcxx(0:nf,ll,mm), fsxx(0:nf,ll,mm), nd, nf)
        enddo
     enddo
          
     do icoil = 1, Ncoils

        if( myid.ne.modulo(icoil-1,ncpu) ) cycle ! parallelization loop;

        do kseg = 0, nd

           fA(kseg) = coil(icoil)%xa(kseg)*coil(icoil)%xt(kseg) + coil(icoil)%ya(kseg)*coil(icoil)%yt(kseg) + &
                      coil(icoil)%za(kseg)*coil(icoil)%zt(kseg)
        enddo ! end kseg

        call coildft( fA, fc, fs, nd, nf)

        llc(icoil, 0:nf, 0) = fc(0:nf)
        lls(icoil, 0:nf, 0) = fs(0:nf)

        leqarc = leqarc + dot_product(coil(icoil)%lmdc,fc) + dot_product(coil(icoil)%lmds,fs)

        do ll = 0, nf

           do kseg = 0, nd

              f1A(kseg,ll     +1) =  -ll**2*cmt(kseg,ll)*coil(icoil)%xt(kseg) - ll*smt(kseg,ll)*coil(icoil)%xa(kseg)
              f1A(kseg,ll+  nf+2) =  -ll**2*smt(kseg,ll)*coil(icoil)%xt(kseg) + ll*cmt(kseg,ll)*coil(icoil)%xa(kseg)
              f1A(kseg,ll+2*nf+3) =  -ll**2*cmt(kseg,ll)*coil(icoil)%yt(kseg) - ll*smt(kseg,ll)*coil(icoil)%ya(kseg)
              f1A(kseg,ll+3*nf+4) =  -ll**2*smt(kseg,ll)*coil(icoil)%yt(kseg) + ll*cmt(kseg,ll)*coil(icoil)%ya(kseg)
              f1A(kseg,ll+4*nf+5) =  -ll**2*cmt(kseg,ll)*coil(icoil)%zt(kseg) - ll*smt(kseg,ll)*coil(icoil)%za(kseg)
              f1A(kseg,ll+5*nf+6) =  -ll**2*smt(kseg,ll)*coil(icoil)%zt(kseg) + ll*cmt(kseg,ll)*coil(icoil)%za(kseg)

           enddo ! end kseg

        enddo ! end ll

        do ll = 1, Cdof
          
           call coildft( f1A(0:Nd,ll), fc, fs, nd, nf)

           l1A(icoil, ll) =  dot_product(coil(icoil)%lmdc,fc) + dot_product(coil(icoil)%lmds,fs)

           llc(icoil, 0:nf, ll) = fc(0:nf)
           lls(icoil, 0:nf, ll) = fs(0:nf)

        enddo ! end ll

        do ll = 1, Cdof          
           do mm = 1, Cdof
              l2A(icoil, ll, icoil, mm) = dot_product(coil(icoil)%lmdc, fcxx(0:nf,ll,mm)) + dot_product(coil(icoil)%lmdc, fsxx(0:nf,ll,mm))
           enddo  !end mm
        enddo !end ll
                                   
     enddo ! end icoil

     call MPI_REDUCE( leqarc, eqarc, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
     RlBCAST(eqarc, 1, 0)

     call MPI_REDUCE( l1A, t1A, array2size, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
     RlBCAST( t1A, array2size, 0 )

     call MPI_REDUCE( llc, dlc, lmdsize, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr ) 
     RlBCAST( dlc, lmdsize, 0 )

     call MPI_REDUCE( lls, dls, lmdsize, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr ) 
     RlBCAST( dls, lmdsize, 0 )

     call MPI_REDUCE( l2A, t2A, array4size, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
     RlBCAST( t2A, array4size, 0 )

     DALLOCATE( f1A )
     DALLOCATE( l1A )
     DALLOCATE( llc )
     DALLOCATE( lls )
     DALLOCATE( df1 )
     DALLOCATE( df2 )
     DALLOCATE( fcxx)
     DALLOCATE( fsxx)
     DALLOCATE( f2A )
     DALLOCATE( l2A )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  end select

  DALLOCATE( fA )
  DALLOCATE( fc )
  DALLOCATE( fs )

  return

end subroutine langrange
     
     
        
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine coildft(data, fcos, fsin, nd, nf)

!----------------------------------------------
! f(x) = data(0:N)
! fsin(n) = 1/pi \int_0^{2pi} f(x) sin(nx) dx
! fcos(n) = 1/pi \int_0^{2pi} f(x) cos(nx) dx ! zero term need additional half factor
!----------------------------------------------

  use kmodule, only: zero, pi2, smt, cmt, myid, ounit
  implicit none
  include "mpif.h"

  REAL,    intent(in ) :: data(0:nd)
  REAL,    intent(out) :: fsin(0:nf), fcos(0:nf)
  INTEGER, intent(in ) :: nd, nf

  INTEGER              ::  astat, ierr, ff, dd


  fsin = zero; fcos = zero

  do ff = 0, nf
     do dd = 0, nd

        fcos(ff) = fcos(ff) + data(dd)*cmt(dd,ff)
        fsin(ff) = fsin(ff) + data(dd)*smt(dd,ff)

     enddo
  enddo

  fcos = fcos / (nd + 1)
  fsin = fsin / (nd + 1)

  fcos(1:nf) = fcos(1:nf)*2 ! nonzero terms
  fsin(1:nf) = fsin(1:nf)*2 ! nonzero terms

  return

end subroutine coildft

! test coildft
!!$subroutine test
!!$
!!$  use kmodule, only: ounit, myid, coil, NDcoil, NFcoil
!!$  implicit none
!!$  include "mpif.h"
!!$
!!$  INTEGER  :: astat, ierr, nn
!!$  REAL     :: fsin1(0:NFcoil), fcos1(0:NFcoil), fsin2(0:NFcoil), fcos2(0:NFcoil), tmp(0:NDcoil-1)
!!$
!!$  if(myid .eq. 0) write(ounit, *) "total segments of COIL 1 is ", size(coil(1)%xx)
!!$
!!$  call FFT( coil(1)%xx, fcos1, fsin1, NDcoil, NFcoil )
!!$
!!$  if(myid .eq. 0) write(ounit,'("nd = "I4," nf = ",I4)') size(tmp), size(fcos2)
!!$  tmp(0:NDcoil-1) = coil(1)%xx(0:NDcoil-1)
!!$  call coildft( tmp, fcos2, fsin2, NDcoil-1, NFcoil )
!!$
!!$  if( myid .eq. 0) then
!!$     write(ounit,*) "fcos1     , fcos2     ,fsin1     , fsin2     ,"
!!$     do nn = 0, NFcoil
!!$        write(ounit,'(4(F10.5,", "))') fcos1(nn), fcos2(nn), fsin1(nn), fsin2(nn)
!!$     enddo
!!$  endif
!!$
!!$  return
!!$
!!$end subroutine test



  
