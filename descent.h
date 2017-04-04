
!title (optimize) ! Minimize &ldquo;energy&rdquo; via a differential flow.

!latex \briefly{Minimize ``energy'' via a differential flow.}

!latex \calledby{\link{notopt}}
!latex \calls{\link{knotxx}}

!latex \tableofcontents

!latex \subsection{overview}
!latex \bi
!latex \item[1.] The geometry of the coils is ``evolved'' under the ``energy gradient flow'', where the energy function, $E({\bf x})$, defined below,
!latex           is considered to be a function of the coil geometry, 
!latex           ${\bf x} = \{ x^{i}_{c,m}, x^{i}_{s,m}, y^{i}_{c,m}, y^{i}_{s,m}, z^{i}_{c,m}, z^{i}_{s,m} \}$, where
!latex           the $x^{i}_{c,m}$ etc. are the Fourier harmonics of the $i$-th coil.
!latex           The Fourier representation of the coils is described in \link{iccoil}.
!latex \item[2.] The evolution is described mathematically as a system of coupled, first-order equations:
!latex           \be \frac{\partial {\bf x}}{\partial \tau} = - \frac{\partial E}{\partial {\bf x}},
!latex           \ee
!latex           where $\tau$ is an artifical time.
!latex \item[3.] The integration is performed using \nag{http://www.nag.co.uk/numeric/FL/manual19/pdf/D02/d02bjf_fl19.pdf}{D02BJF}, 
!latex           and is controlled by \inputvar{tauend}, \inputvar{tautol} and \inputvar{Ntauout}.
!latex \ei

!latex \subsection{derivarives test}
!latex \bi
!latex \item[1.] When \inputvar{Loptimize = -1}, the code will test the first derivatives comparing the results of finite difference, shudson`s and czhu`s method.
!latex Here, minus sign in shudson`s method was temporarily eliminated.
!latex \item[2.] When \inputvar{Loptimize = -2}, the code will test the second derivatives comparing the results of finite difference and czhu`s method. The results of
!latex finite difference come from using the first derivatives of shudson`s method differentiating on the intervals. That is,
!latex \be  \frac{\partial^2{E}}{\partial{x_1} \partial{x_2}} \equiv \frac{\partial{F_1}}{\partial{x_2}} 
!latex      \equiv \frac{F_1(x_2 + \frac{1}{2}\Delta) - F_1(x_2 - \frac{1}{2}\Delta)}{\Delta}
!latex \ee
!latex Right now, the results of $\frac{\partial^2{E}}{\partial{I_i}^2}$ don`t match. But the results of czhu`s method and finite difference result from 0-order 
!latex functions, which is $\frac{\partial^2{E}}{\partial{I_i}^2} \approx \frac{E(I_i-\frac{1}{2}\Delta)-2E(I_i)+E(I_i-\frac{1}{2}\Delta)}{\Delta ^2}$, do match.\\
!latex \item[3.] \emph{ The second derivatives of energy function $E$ on which variable were decided by the value jcoil and jdof, representing the coil number and 
!latex the jdof-th DoF in each coil. Be careful with the differenced of DoF array arrangements between czhu`s method and 
!latex shudson`s method.}
!latex \item[3.] When \inputvar{Loptimize = 1}, the code will perform differential flow method using 
!latex                \nag{http://www.nag.co.uk/numeric/FL/manual19/pdf/D02/d02bjf_fl19.pdf}{D02BJF}, with shudson's method to get the first derivatives. It's
!latex                really fast and you can turn on/off Io, Lo to control the constrains of currents and length.
!latex \item[4.] When \inputvar{Loptimize = 2}, the code will perform differential flow method using 
!latex                \nag{http://www.nag.co.uk/numeric/FL/manual19/pdf/D02/d02bjf_fl19.pdf}{D02BJF}, with czhu's method to get the first derivatives. It's
!latex                not so fast, but it can calculate the second derivatives. Like shudson's method, you can turn on/off the constrains of bnormal, toroidal flux 
!latex                and length (and more inthe later) through changing the value of \inputvar{weight\_bnorm}, \inputvar{weight\_tflux} and \inputvar{weight\_ttlen}.
!latex \ei


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine descent
  
  use kmodule, only : zero, half, myid, ncpu, ounit, &
                      NFcoil, &
                      Loptimize, Ndof, &
                      tauend, tautol, itau, &
                      icoil, Ncoils, coil
  
  implicit none
  
  include "mpif.h"
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  INTEGER              :: astat, ierr

  INTEGER              :: NN, idof, irelab, id02bjf, Lwork, mm
  REAL                 :: tstart, tend, tol, tau
  REAL   , allocatable :: xdof(:), work(:)
  CHARACTER            :: relabs
  
  external             :: denergy, progres, denergy2, progres2, D02BJW
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  NN = NFcoil  ! this assumes the Fourier resolution of each coil is the same; not general; 14 Apr 16;
  
  SALLOCATE( xdof, (1:Ndof), zero )
  
  call pack( xdof )

  if(myid .eq. 0) write(ounit,'("descent : " 10x " : Begin using differential flow with simple energy function.")')   
   
   Lwork = Ndof * 20
   SALLOCATE( work, (1:Lwork), zero )
   
   itau = 0 ; tstart = zero ; tend = tauend ; tol = tautol ; relabs = 'D'
   
   if( myid.eq.0 ) write(ounit,1000) "calling", Ndof
   
   id02bjf = 1 ; call D02BJF( tstart, tend, Ndof, xdof(1:Ndof), denergy, tol, relabs, progres, D02BJW, work(1:Lwork), id02bjf )
   
   if( myid.eq.0 ) write(ounit,1000) "called ", Ndof, id02bjf
   
1000 format("descent : " 10x " : "a7" D02BJF ; Ndof =",i7," ;":" id02bjf ="i3" ;")
   
   DALLOCATE( work )
  
   call unpack( xdof )

  return
  
end subroutine descent
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine descent2
  
  use kmodule, only : zero, half, myid, ncpu, ounit, &
                      NFcoil, &
                      Loptimize, Ndof, &
                      tauend, tautol, itau, &
                      icoil, Ncoils, coil
  
  implicit none
  
  include "mpif.h"
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  INTEGER              :: comm, astat, ierr

  INTEGER              :: NN, idof, irelab, id02bjf, Lwork, mm
  REAL                 :: tstart, tend, tol, tau
  REAL   , allocatable :: xdof(:), work(:)
  CHARACTER            :: relabs
  
  external             :: denergy2, progres2, D02BJW
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  NN = NFcoil  ! this assumes the Fourier resolution of each coil is the same; not general; 14 Apr 16;
  
  SALLOCATE( xdof, (1:Ndof), zero )
  
  call pack( xdof )

   if(myid .eq. 0) write(ounit,'("descent : " 10x " : Begin using differential flow with integrated energy function.")')
   
   Lwork = Ndof * 20
   SALLOCATE( work, (1:Lwork), zero )
   
   itau = 0 ; tstart = zero ; tend = tauend ; tol = tautol ; relabs = 'D' 
   
   if( myid.eq.0 ) write(ounit,1001) "calling", Ndof
   
   id02bjf = 1 ; call D02BJF( tstart, tend, Ndof, xdof(1:Ndof), denergy2, tol, relabs, progres2, D02BJW, work(1:Lwork), id02bjf )
   
   if( myid.eq.0 ) write(ounit,1001) "called ", Ndof, id02bjf
   
1001 format("descent : " 10x " : "a7" D02BJF ; Ndof =",i7," ;":" id02bjf ="i3" ;")
 
   DALLOCATE( work )
  
   call unpack( xdof )

   call mpi_barrier(MPI_COMM_WORLD, ierr)

  return
  
end subroutine descent2

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine progres( tau, xdof )
  
  use kmodule, only : zero, ounit, myid, &
                      NFcoil, &
                      itau, tauend, Ntauout, &
                      coil, icoil, Ncoils, totalenergy, evolution
  
  implicit none
  
  include "mpif.h"

  REAL                 :: tau, xdof(*)
  
  INTEGER              :: irestart, NN, Ndof, astat
  REAL   , allocatable :: dE(:)

  NN = NFcoil ; Ndof = Ncoils * ( 1 + 3 * ( 1 + NN + NN ) ) ; irestart = 1

  SALLOCATE( dE, (1:Ndof), zero )
  
  call denergy( tau, xdof(1:Ndof), dE(1:Ndof) )

  if( myid.eq.0 ) write(ounit,1000) tau, totalenergy, sqrt(sum(dE(1:Ndof)**2)/Ndof)

  DALLOCATE( dE )

  evolution(itau,0) = tau
  evolution(itau,1) = totalenergy

  irestart = 1 ; call restart( irestart )

  itau = itau + 1 ; tau = itau * tauend / Ntauout

  return  

1000 format("progres :"es11.4" : E="es23.15" ; D="es23.15" ;")

end subroutine progres


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine progres2( tau, xdof )
  
  use kmodule, only : zero, ounit, myid, iter, &
                      NFcoil, &
                      itau, tauend, Ntauout, Savfreq, tautol, &
                      coil, icoil, Ncoils, totalenergy, evolution, bnorm, tflux, ttlen, eqarc, ccsep, tbn
  
  implicit none
  
  include "mpif.h"

  REAL                 :: tau, xdof(*), sumdE
  
  INTEGER              :: irestart, NN, Ndof, astat, ierr
  REAL   , allocatable :: dE(:)

  NN = NFcoil ; Ndof = Ncoils * ( 1 + 3 * ( 1 + NN + NN ) ) ; irestart = 1

  SALLOCATE( dE, (1:Ndof), zero )
  
  call denergy2( tau, xdof(1:Ndof), dE(1:Ndof) )

  sumdE = sqrt(sum(dE(1:Ndof)**2)/Ndof)

  if( myid.eq.0 ) write(ounit,1000) tau, totalenergy, sumdE, bnorm, tflux, ttlen, eqarc, ccsep

  DALLOCATE( dE )

  FATAL(progres2, itau.gt.Ntauout, exceeds allocation)

  call unpack(xdof(1:Ndof))

  if(allocated(evolution)) then
     evolution(itau,0) = tau
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

  if( itau .gt. 1 ) then 
     if ( abs(totalenergy-evolution(itau-1,1))/totalenergy .lt. 1E-8 ) then
     tauend = itau
     call restart( irestart )
     tau = (itau-1) * tauend / Ntauout ! go backwards, which will cause an ifail=5 error;
     iter = 0 ; return
    !call MPI_ABORT( MPI_COMM_WORLD, 1, ierr )
    !STOP "Adaptive precision (tautol) reached!"
  endif  
  endif
  
  if(mod(itau,Savfreq) .eq. 0) call restart( irestart )

  itau = itau + 1 ; tau = itau * tauend / Ntauout

  iter = 0 !reset iter

  return  

1000 format("progres :"es11.4" : E="es15.7" ; D="es15.7" ; B="es15.7" ; F="es15.7" ; L="es15.7" ; A="es15.7" ; C="es15.7" ;")

end subroutine progres2
