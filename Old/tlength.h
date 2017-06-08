!title (tlength) ! Calculate total length cost functon and its derivatives. (czhu)

!latex \briefly{Calling tlength(nderiv) to calculate the length cost function w/o its derivatives on all degrees of freedom. When \emph{nderiv = 0}, it only 
!latex          calculates the $0-order$ length cost function, which is $/tlength = \sum_{i=1,Ncoils}\frac{1}{2}(L_i - L_o^i)^2/$. When \emph{nderiv = 0}, it will 
!latex          calculate both the $0-order$ and $1^{st}-order$ derivatives in array \emph{t1L(1:Ncoils,0:Codf)}. When \emph{nderiv = 0}, it will calculate all the 
!latex          $0-order$, $1^{st}-order$ and $2^{nd}-order$ derivatives $t2L(1:Ncoils,0:Codf,1:Ncoils,0:Codf)$.}

!latex \calledby{\link{denergy}}

!latex \tableofcontents

!latex \newcommand{\cmt}{\cos(mt)}
!latex \newcommand{\smt}{\sin(mt)}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!latex \subsection{Total length cost function}
!latex \bi
!latex \item[1.] The length of each coils can be calculated through,
!latex \be
!latex L_i = \int_{icoil} \sqrt{\dot{x}(t)^2 + \dot{y}(t)^2 + \dot{z}(t)^2} dt \approx 
!latex      \sum_{kseg=0,Ndcoil-1} \sqrt{xt(kseg)^2+yt(kseg)^2+zt(kseg)^2} \frac{2\pi}{NDcoil}
!latex \ee
!latex The results of user discretized method and \nag{http://www.nag.co.uk/numeric/FL/manual19/pdf/D01/d01eaf_fl19.pdf}{D01EAF} adaptive routine 
!latex can be compared through subroutine \link{descent}.\\
!latex \item[2.] And then the cost funtion on length is
!latex \be  
!latex ttlen = \sum_{i=1,Ncoils}\frac{1}{2}Lw_i (L_i - L_o^i)^2
!latex \ee
!latex where $Lw_i$ is the weight of length of each coil, while there is another weight on the length cost function $weight\_ttlen$
!latex \ei
!latex \subsection{First dirivatives}
!latex Since the length of coils has no relationships with the current in the coil, the first derivatives of $ttlen$ on currents are all zero. 
!latex An the derivatives on the geometry variables can be represented as,
!latex \be
!latex \frac{\partial{ttlen}}{\partial{x_n^i}} \equiv (L_i - L_o^i) \, \frac{\partial{L_i}}{\partial{x_n^i}}
!latex \ee
!latex Here $i$ is denoted to the $i^{th}$ coil and $x_n$ means the $n^{th}$ DoF. And
!latex \be
!latex \frac{\partial{L_i}}{\partial{x_n^i}} \equiv \int_{icoil} \frac{\dot{x}}{\sqrt{\dot{x}(t)^2 + \dot{y}(t)^2 + \dot{z}(t)^2}} 
!latex                                                          \frac{\partial{\dot{x}}}{\partial{x_n^i}} \, dt
!latex \ee
!latex \subsection{Second derivatives}
!latex The second derivatives of length cost function are only related to geometry variables of the same coil. That is,
!latex \be
!latex \frac{\partial^2{ttlen}}{\partial{x_n^i}\partial{x_m^i}} \equiv \frac{\partial{L_i}}{\partial{x_m^i}} \, \frac{\partial{L_i}}{\partial{x_n^i}} + 
!latex                                                                  (L_i - L_o^i)\frac{\partial^2{L_i}}{\partial{x_n^i}\partial{x_m^i}}
!latex                                                           \equiv \int_{icoil} \frac{1}{\sqrt{\dot{x}(t)^2 + \dot{y}(t)^2 + \dot{z}(t)^2}} 
!latex                                                                               \frac{\partial{\dot{x}}}{\partial{x_n^i}}  \frac{\partial{\dot{x}}}{\partial{x_m^i}} \,-\, 
!latex                                                                               \frac{\dot{x} \, \dot{x}}{(\sqrt{\dot{x}(t)^2 + \dot{y}(t)^2 + \dot{z}(t)^2})^3} \, dt
!latex \ee
!latex 
!latex \subsection{Normalization}
!latex While normalizing length cost function, we should divide it with the target length of each coils. That is,
!latex \begin{equation}
!latex ttlen = \sum_{i=1,Ncoils}\frac{1}{2}Lw_i \frac{(L_i - L_o^i)^2}{{L_o^i}^2}
!latex \end{equation}
!latex Since the target lengths are user specified constant, the normalization of length cost function derivatives can be implemented by dividing them with each coil's target length.
!latex 



! seems paralization works slower than series. 06/24
subroutine tlength(nderiv)
  use kmodule, only : zero, half, pi2, &
                      coil, Ncoils, NFcoil, NDcoil, smt, cmt, Cdof, &
                      ttlen, t1L, t2L, &
                      ncpu, myid, ounit
  implicit none
  include "mpif.h"
  INTEGER           :: nderiv
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  INTEGER           :: astat, ierr
  INTEGER           :: NN, icoil, kseg, ll, c1, n1, c2, n2, array2size , array4size
  REAL              :: dlength, llength
  REAL, allocatable :: d1L(:), d2L(:,:), l1L(:,:), l2L(:,:,:,:)
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  ttlen = zero; llength = zero; NN = NFcoil

  array2size = Ncoils * ( Cdof + 1 ); array4size = array2size * array2size
  
  select case ( nderiv )
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  case (0)

     do icoil = 1, Ncoils

        if( myid.ne.modulo(icoil-1,ncpu) ) cycle ! parallelization loop;
        call LenDeriv0(icoil, coil(icoil)%L)
        llength = llength + coil(icoil)%Lw * ( coil(icoil)%L - coil(icoil)%Lo )**2 / coil(icoil)%Lo**2

     enddo

     do icoil = 1, Ncoils
      RlBCAST( coil(icoil)%L, 1, modulo(icoil-1,ncpu) ) !broadcast each coil's length
     enddo

     call MPI_REDUCE( llength, ttlen, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
     RlBCAST( ttlen, 1, 0)
     ttlen = half * ttlen / Ncoils ! Averaged on each coils; need to be normalized;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  case (1)
     if ( .not. allocated(t1L) ) then
        SALLOCATE(t1L, (1:Ncoils, 0:Cdof), zero)
     else
        t1L = zero
     endif

     SALLOCATE(d1L, (0:Cdof), zero)
     SALLOCATE(l1L, (1:Ncoils, 0:Cdof), zero)

     do icoil = 1, Ncoils

        if( myid.ne.modulo(icoil-1,ncpu) ) cycle ! parallelization loop;

        call LenDeriv1( icoil, d1L(0:Cdof) )
        coil(icoil)%L  = d1L(0)
  
        llength = llength + coil(icoil)%Lw * ( coil(icoil)%L - coil(icoil)%Lo )**2 / coil(icoil)%Lo**2

        do ll = 1, Cdof
           l1L(icoil,ll) = coil(icoil)%Lw * ( coil(icoil)%L - coil(icoil)%Lo ) * d1L(ll) / coil(icoil)%Lo**2
        enddo ! end ll

     enddo ! end icoil

     do icoil = 1, Ncoils
      RlBCAST( coil(icoil)%L, 1, modulo(icoil-1,ncpu) )
     enddo

     call MPI_REDUCE( llength, ttlen, 1         , MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
     RlBCAST( ttlen,     1, 0)

     call MPI_REDUCE( l1L    , t1L  , array2size, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
     RlBCAST( t1L, array2size, 0)
         
     ttlen = half * ttlen / Ncoils ! Averaged on each coils; need to be normalized;
     t1L   =          t1L / Ncoils

     DALLOCATE(d1L)
     DALLOCATE(l1L)
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  case (2)
     if ( .not. allocated(t1L) ) then
        SALLOCATE(t1L, (1:Ncoils, 0:Cdof), zero)
     else
        t1L = zero
     endif
     if ( .not. allocated(t2L) ) then
        SALLOCATE(t2L, (1:Ncoils, 0:Cdof, 1:Ncoils, 0:Cdof), zero)
     else
        t2L = zero
     endif

     SALLOCATE(l1L, (1:Ncoils, 0:Cdof), zero)
     SALLOCATE(l2L, (1:Ncoils, 0:Cdof, 1:Ncoils, 0:Cdof), zero)
     SALLOCATE(d2L, (0:Cdof  , 0:Cdof), zero)


     do icoil = 1, Ncoils

        if( myid.ne.modulo(icoil-1,ncpu) ) cycle ! parallelization loop;

        call LenDeriv2( icoil, d2L(0:Cdof,0:Cdof) )

        coil(icoil)%L = d2L(0,0)
             
        llength = llength + coil(icoil)%Lw * ( coil(icoil)%L - coil(icoil)%Lo )**2 / coil(icoil)%Lo**2
        !print *, "myid = ", myid, " ; llength = ", llength,  coil(icoil)%Lw, coil(icoil)%L, coil(icoil)%Lo
        do ll = 1, Cdof
           l1L(icoil,ll) =  ( coil(icoil)%L - coil(icoil)%Lo ) * d2L(ll,0) / coil(icoil)%Lo**2
        enddo ! end ll

        do n1 = 1, Cdof
         do n2 = 1, Cdof
          l2L(icoil,n1,icoil,n2) = ( d2L(n1,0) * d2L(n2,0) + (coil(icoil)%L - coil(icoil)%Lo) * d2L(n1,n2) ) / coil(icoil)%Lo**2
         enddo
        enddo
        
     enddo ! end icoil

     call MPI_REDUCE( llength, ttlen, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
     RlBCAST( ttlen, 1, 0)

     call MPI_REDUCE( l1L    , t1L  , array2size, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
     RlBCAST( t1L, array2size, 0)

     do icoil = 1, Ncoils
      RlBCAST( coil(icoil)%L, 1, modulo(icoil-1,ncpu) )
     enddo
     
     call MPI_REDUCE( l2L    , t2L  , array4size, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
     RlBCAST( t2L, array4size, 0 )
  

     ttlen = half * ttlen / Ncoils ! Averaged on each coils; need to be normalized;
     t1L     =        t1L / Ncoils
     t2L     =        t2L / Ncoils

     
     DALLOCATE(l1L)
     DALLOCATE(d2L)
     DALLOCATE(l2L)
     
  end select

  return
end subroutine tlength

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! calculate coil length

subroutine LenDeriv0(icoil, length)

  use kmodule, only: zero, pi2, coil, Ncoils, NFcoil, NDcoil, myid, ounit
  implicit none
  include "mpif.h"

  INTEGER, intent(in)  :: icoil
  REAL   , intent(out) :: length
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  INTEGER              :: kseg, astat, ierr
  REAL                 :: dlength

  FATAL( LenDeriv0, icoil .lt. 1 .or. icoil .gt. Ncoils, icoil not in right range )
  
  length = zero
  
  do kseg = 0, NDcoil-1

     dlength = sqrt(coil(icoil)%xt(kseg)**2 + coil(icoil)%yt(kseg)**2 + coil(icoil)%zt(kseg)**2)
     length  = length + dlength  

  enddo ! end kseg

  length = length * pi2 / NDcoil

  return

end subroutine LenDeriv0

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!calculate coil length in derivs(0) and first derivatives in derivs(1:Cdof)

subroutine LenDeriv1(icoil, derivs)

  use kmodule, only: zero, pi2, coil, Ncoils, NFcoil, NDcoil, smt, cmt, Cdof, myid, ounit
  implicit none
  include "mpif.h"

  INTEGER, intent(in)  :: icoil
  REAL   , intent(out) :: derivs(0:Cdof)
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  INTEGER              :: kseg, ll, NN, astat, ierr
  REAL                 :: dlength

  FATAL( LenDeriv1, icoil .lt. 1 .or. icoil .gt. Ncoils, icoil not in right range )
  
  NN = NFcoil; derivs = zero
  
  do kseg = 0, NDcoil-1

     dlength = sqrt(coil(icoil)%xt(kseg)**2 + coil(icoil)%yt(kseg)**2 + coil(icoil)%zt(kseg)**2)
     derivs(0) = derivs(0) + dlength

     do ll = 0, NN
        derivs(ll        + 1) = derivs(ll        + 1) + coil(icoil)%xt(kseg) * ( -ll*smt(kseg, ll) ) / dlength
        derivs(ll +   NN + 2) = derivs(ll +   NN + 2) + coil(icoil)%xt(kseg) * (  ll*cmt(kseg, ll) ) / dlength
        derivs(ll + 2*NN + 3) = derivs(ll + 2*NN + 3) + coil(icoil)%yt(kseg) * ( -ll*smt(kseg, ll) ) / dlength
        derivs(ll + 3*NN + 4) = derivs(ll + 3*NN + 4) + coil(icoil)%yt(kseg) * (  ll*cmt(kseg, ll) ) / dlength
        derivs(ll + 4*NN + 5) = derivs(ll + 4*NN + 5) + coil(icoil)%zt(kseg) * ( -ll*smt(kseg, ll) ) / dlength
        derivs(ll + 5*NN + 6) = derivs(ll + 5*NN + 6) + coil(icoil)%zt(kseg) * (  ll*cmt(kseg, ll) ) / dlength
     enddo ! end ll  

  enddo ! end kseg
  derivs = derivs * pi2 / NDcoil
  return

end subroutine LenDeriv1

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! calculate coil length in derivs(0,0), first derivatives in derivs(1:Cdof,0) and second derivatives in derivs(1:Cdof, 1:Cdof)

subroutine LenDeriv2( icoil, derivs )

  use kmodule, only: zero, pi2, coil, Ncoils, NFcoil, NDcoil, smt, cmt, Cdof, myid, ounit
  implicit none
  include "mpif.h"

  INTEGER, intent(in)  :: icoil
  REAL   , intent(out) :: derivs(0:Cdof, 0:Cdof)
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  INTEGER              :: kseg, ll, NN, n1, n2, astat, ierr
  REAL                 :: dlength, ltx(0:Cdof), lty(0:Cdof), ltz(0:Cdof)

  FATAL( LenDeriv2, icoil .lt. 1 .or. icoil .gt. Ncoils, icoil not in right range )

  NN = NFcoil
  ltx = zero; lty = zero; ltz = zero; derivs = zero

  do kseg = 0, NDcoil-1

   ltx(0) = coil(icoil)%xt(kseg)
   lty(0) = coil(icoil)%yt(kseg)
   ltz(0) = coil(icoil)%zt(kseg)

   do ll = 0, NN
    ltx(ll        + 1) = -ll*smt(kseg, ll)
    ltx(ll +   NN + 2) =  ll*cmt(kseg, ll)
    lty(ll + 2*NN + 3) = -ll*smt(kseg, ll)
    lty(ll + 3*NN + 4) =  ll*cmt(kseg, ll)
    ltz(ll + 4*NN + 5) = -ll*smt(kseg, ll)
    ltz(ll + 5*NN + 6) =  ll*cmt(kseg, ll)
   enddo ! end ll

   dlength = sqrt( ltx(0)**2 + lty(0)**2 + ltz(0)**2 )
   derivs(0,0)  = derivs(0,0) + dlength

   do ll = 1, 2*NN+2
    derivs(ll           ,0) = derivs(ll           ,0) + ltx(0) * ltx(ll           ) / dlength
    derivs(ll + 2*NN + 2,0) = derivs(ll + 2*NN + 2,0) + lty(0) * lty(ll + 2*NN + 2) / dlength
    derivs(ll + 4*NN + 4,0) = derivs(ll + 4*NN + 4,0) + ltz(0) * ltz(ll + 4*NN + 4) / dlength
   enddo ! end ll        

   do n1 = 1, 2*NN+2
    do n2 = 1, 2*NN+2

     derivs(n1       ,n2       ) = derivs(n1       ,n2       ) + ltx(n2       )*ltx(n1       )/dlength  - ltx(n2       )*ltx(n1       )*ltx(0)*ltx(0)/dlength**3 !l/x/x
     derivs(n1       ,n2+2*NN+2) = derivs(n1       ,n2+2*NN+2) +                                     0  - lty(n2+2*NN+2)*ltx(n1       )*lty(0)*ltx(0)/dlength**3 !l/x/y
     derivs(n1       ,n2+4*NN+4) = derivs(n1       ,n2+4*NN+4) +                                     0  - ltz(n2+4*NN+4)*ltx(n1       )*ltz(0)*ltx(0)/dlength**3 !l/x/z

     derivs(n1+2*NN+2,n2       ) = derivs(n1+2*NN+2,n2       ) +                                     0  - ltx(n2       )*lty(n1+2*NN+2)*ltx(0)*lty(0)/dlength**3 !l/y/x
     derivs(n1+2*NN+2,n2+2*NN+2) = derivs(n1+2*NN+2,n2+2*NN+2) + lty(n2+2*Nn+2)*lty(n1+2*NN+2)/dlength  - lty(n2+2*NN+2)*lty(n1+2*NN+2)*lty(0)*lty(0)/dlength**3 !l/y/y
     derivs(n1+2*NN+2,n2+4*NN+4) = derivs(n1+2*NN+2,n2+4*NN+4) +                                     0  - ltz(n2+4*NN+4)*lty(n1+2*NN+2)*ltz(0)*lty(0)/dlength**3 !l/y/z

     derivs(n1+4*NN+4,n2       ) = derivs(n1+4*NN+4,n2       ) +                                     0  - ltx(n2       )*ltz(n1+4*NN+4)*ltx(0)*ltz(0)/dlength**3 !l/z/x
     derivs(n1+4*NN+4,n2+2*NN+2) = derivs(n1+4*NN+4,n2+2*NN+2) +                                     0  - lty(n2+2*NN+2)*ltz(n1+4*NN+4)*lty(0)*ltz(0)/dlength**3 !l/z/y
     derivs(n1+4*NN+4,n2+4*NN+4) = derivs(n1+4*NN+4,n2+4*NN+4) + ltz(n2+4*NN+4)*ltz(n1+4*NN+4)/dlength  - ltz(n2+4*NN+4)*ltz(n1       )*ltz(0)*ltz(0)/dlength**3 !l/z/z

    enddo
   enddo

  enddo ! end kseg

   derivs = derivs * pi2 / NDcoil

  return
end subroutine LenDeriv2


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine tlengthExp(nderiv)

  use kmodule, only : zero, half, pi2, &
                      coil, Ncoils, NFcoil, NDcoil, smt, cmt, Cdof, &
                      ttlen, t1L, t2L, Lo, &
                      ncpu, myid, ounit
  implicit none
  include "mpif.h"
  INTEGER           :: nderiv
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  INTEGER           :: astat, ierr
  INTEGER           :: NN, icoil, kseg, ll, c1, n1, c2, n2, array2size , array4size
  REAL              :: dlength, llength, normalize
  REAL, allocatable :: d1L(:), d2L(:,:), l1L(:,:), l2L(:,:,:,:)
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  ttlen = zero; llength = zero; NN = NFcoil; normalize = zero

  array2size = Ncoils * ( Cdof + 1 ); array4size = array2size * array2size

  do icoil = 1, Ncoils
     normalize = normalize + exp(coil(icoil)%Lw*Lo)
  enddo
   
  select case ( nderiv )
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  case (0)

     do icoil = 1, Ncoils

        if( myid.ne.modulo(icoil-1,ncpu) ) cycle ! parallelization loop;
        call LenDeriv0(icoil, coil(icoil)%L)
        llength = llength + exp(coil(icoil)%Lw*coil(icoil)%L)

     enddo

     do icoil = 1, Ncoils
      RlBCAST( coil(icoil)%L, 1, modulo(icoil-1,ncpu) ) !broadcast each coil's length
     enddo

     call MPI_REDUCE( llength, ttlen, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
     RlBCAST( ttlen, 1, 0)
     ttlen = ttlen / normalize! Averaged on each coils and normalized with target Length;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  case (1)
     if ( .not. allocated(t1L) ) then
        SALLOCATE(t1L, (1:Ncoils, 0:Cdof), zero)
     else
        t1L = zero
     endif

     SALLOCATE(d1L, (0:Cdof), zero)
     SALLOCATE(l1L, (1:Ncoils, 0:Cdof), zero)

     do icoil = 1, Ncoils

        if( myid.ne.modulo(icoil-1,ncpu) ) cycle ! parallelization loop;

        call LenDeriv1( icoil, d1L(0:Cdof) )
        coil(icoil)%L  = d1L(0)
  
        llength = llength + exp(coil(icoil)%Lw*coil(icoil)%L)

        do ll = 1, Cdof
           l1L(icoil,ll) = exp(coil(icoil)%Lw*coil(icoil)%L) * coil(icoil)%Lw * d1L(ll)
        enddo ! end ll

     enddo ! end icoil

     do icoil = 1, Ncoils
      RlBCAST( coil(icoil)%L, 1, modulo(icoil-1,ncpu) )
     enddo

     call MPI_REDUCE( llength, ttlen, 1         , MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
     RlBCAST( ttlen,     1, 0)

     call MPI_REDUCE( l1L    , t1L  , array2size, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
     RlBCAST( t1L, array2size, 0)
         
     ttlen =  ttlen / normalize ! Averaged on each coils; need to be normalized;
     t1L   =    t1L / normalize

     DALLOCATE(d1L)
     DALLOCATE(l1L)
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  case (2)
     if ( .not. allocated(t1L) ) then
        SALLOCATE(t1L, (1:Ncoils, 0:Cdof), zero)
     else
        t1L = zero
     endif
     if ( .not. allocated(t2L) ) then
        SALLOCATE(t2L, (1:Ncoils, 0:Cdof, 1:Ncoils, 0:Cdof), zero)
     else
        t2L = zero
     endif

     SALLOCATE(l1L, (1:Ncoils, 0:Cdof), zero)
     SALLOCATE(l2L, (1:Ncoils, 0:Cdof, 1:Ncoils, 0:Cdof), zero)
     SALLOCATE(d2L, (0:Cdof  , 0:Cdof), zero)


     do icoil = 1, Ncoils

        if( myid.ne.modulo(icoil-1,ncpu) ) cycle ! parallelization loop;

        call LenDeriv2( icoil, d2L(0:Cdof,0:Cdof) )

        coil(icoil)%L = d2L(0,0)
             
        llength = llength + exp(coil(icoil)%Lw*coil(icoil)%L)
        do ll = 1, Cdof
           l1L(icoil,ll) =  exp(coil(icoil)%Lw*coil(icoil)%L) * coil(icoil)%Lw * d2L(ll,0)
        enddo ! end ll

        do n1 = 1, Cdof
         do n2 = 1, Cdof
          l2L(icoil,n1,icoil,n2) = coil(icoil)%Lw * exp(coil(icoil)%Lw*coil(icoil)%L) * ( d2L(n1,0)*d2L(n2,0) + d2L(n1,n2))
         enddo
        enddo
        
     enddo ! end icoil

     call MPI_REDUCE( llength, ttlen, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
     RlBCAST( ttlen, 1, 0)

     call MPI_REDUCE( l1L    , t1L  , array2size, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
     RlBCAST( t1L, array2size, 0)

     do icoil = 1, Ncoils
      RlBCAST( coil(icoil)%L, 1, modulo(icoil-1,ncpu) )
     enddo
     
     call MPI_REDUCE( l2L    , t2L  , array4size, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
     RlBCAST( t2L, array4size, 0 )
  

     ttlen = ttlen / normalize ! Averaged on each coils and be normalized;
     t1L   =   t1L / normalize
     t2L   =   t2L / normalize

     
     DALLOCATE(l1L)
     DALLOCATE(d2L)
     DALLOCATE(l2L)
     
  end select

  return
end subroutine tlengthExp
