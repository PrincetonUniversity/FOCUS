
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


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine descent
  !---------------------------------------------------------------------------------------------
  ! Using the differential flow to optimize the coils;
  ! DATE: 2017/04/05
  !---------------------------------------------------------------------------------------------    
  use globals, only : zero, half, myid, ncpu, ounit, &
        Ndof, itau, SD_tauend, SD_tausta, SD_tautol, SD_Nout, itau

  implicit none  
  include "mpif.h"

  !---------------------------------------------------------------------------------------------     
  INTEGER              :: comm, astat, ierr

  INTEGER              :: irelab, id02bjf, Lwork, mm
  REAL                 :: tstart, tend, tol, lxdof(1:Ndof)
  REAL   , allocatable :: work(:)
  CHARACTER            :: relabs

  external             :: denergy, progres, D02BJW
  !---------------------------------------------------------------------------------------------    

  if(myid == 0) write(ounit, *) "-----------Differential Flow Optimizer-----------------------"
  FATAL( descent, Ndof < 1, INVALID Ndof value )
  if(myid .eq. 0) write(ounit,'("descent : Begin optimizations using differential flow with " &
       &                         I7 " degrees of freedom.")') Ndof

  call packdof(lxdof) !copy to local;

  Lwork = Ndof * 20
  SALLOCATE( work, (1:Lwork), zero )

  itau = 0 ; tstart = SD_tausta ; tend = SD_tauend ; tol = SD_tautol ; relabs = 'D' 

  id02bjf = 1
  call D02BJF( tstart, tend, Ndof, lxdof(1:Ndof), denergy, tol, relabs, progres, D02BJW, work(1:Lwork), id02bjf )

  if(myid .eq. 0) write(ounit,'("descent : Differential flow finished. id02bjf = " I3)') id02bjf

  DALLOCATE( work )

  call unpacking(lxdof)

  call mpi_barrier(MPI_COMM_WORLD, ierr)

  return

end subroutine descent

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine denergy( tau, lxdof, dE )
  
  use globals, only :  Ndof, myid, ncpu, ounit, t1E, dofnorm
  implicit none
  include "mpif.h"
  !---------------------------------------------------------------------------------------------      
  REAL                 :: tau, lxdof(*), dE(*)
  
  INTEGER              :: ierr, astat, iorder
  external             :: unpacking, costfun
  !---------------------------------------------------------------------------------------------    

  call unpacking(lxdof(1:Ndof))

  iorder = 1
  call costfun(iorder)
 
  dE(1:Ndof) = - (t1E(1:Ndof) * dofnorm(1:Ndof) )

  return
end subroutine denergy


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine progres( tau, lxdof )
  
  use globals, only : zero, sqrtmachprec, ounit, myid, iter, Ndof, Ncoils, Tdof, coil, FouCoil, coilspace, &
                      itau, SD_tausta, SD_tauend, SD_Nout, SD_SaveFreq, SD_tautol, &
                      totalenergy, evolution, bnorm, tflux, ttlen, specw, ccsep, bharm, t1E, dofnorm
  
  implicit none  
  include "mpif.h"
  !---------------------------------------------------------------------------------------------    
  REAL    :: tau, lxdof(*)
  
  INTEGER :: irestart, iorder, astat, ierr, idof, NF, icoil
  REAL    :: sumdE, dE(1:Ndof)
  !---------------------------------------------------------------------------------------------    

  irestart = 1 ;  iorder = 1

  call unpacking(lxdof(1:Ndof))
  call costfun(iorder)
  dE = - (t1E * dofnorm)
  sumdE = sqrt(sum(dE(1:Ndof)**2)/Ndof)
  if( myid==0 ) write(ounit,1000) tau, totalenergy, sumdE, bnorm, tflux, ttlen, specw, ccsep, bharm

  FATAL(progres, itau > SD_Nout, exceeds allocation)

  !save evolution;
  if(allocated(evolution)) then
     evolution(itau,0) = tau
     evolution(itau,1) = totalenergy
     evolution(itau,2) = sumdE
     evolution(itau,3) = bnorm
     evolution(itau,4) = tflux
     evolution(itau,5) = ttlen
     evolution(itau,6) = specw
     evolution(itau,7) = ccsep
     evolution(itau,8) = bharm
  endif

  !save all the coil parameters;
  idof = 0
  do icoil = 1, Ncoils
     coilspace(itau, idof+1          ) = coil(icoil)%I          ;  idof = idof + 1

     select case (coil(icoil)%itype)
     case (1)
        NF = FouCoil(icoil)%NF

        coilspace(itau, idof+1:idof+NF+1) = FouCoil(icoil)%xc(0:NF) ; idof = idof + NF +1
        coilspace(itau, idof+1:idof+NF  ) = FouCoil(icoil)%xs(1:NF) ; idof = idof + NF
        coilspace(itau, idof+1:idof+NF+1) = FouCoil(icoil)%yc(0:NF) ; idof = idof + NF +1
        coilspace(itau, idof+1:idof+NF  ) = FouCoil(icoil)%ys(1:NF) ; idof = idof + NF
        coilspace(itau, idof+1:idof+NF+1) = FouCoil(icoil)%zc(0:NF) ; idof = idof + NF +1
        coilspace(itau, idof+1:idof+NF  ) = FouCoil(icoil)%zs(1:NF) ; idof = idof + NF
     case default
        FATAL(descent, .true., not supported coil types)
     end select
  enddo
  FATAL( descent , idof .ne. Tdof, counting error in restart )

  if ( itau > 1 ) then
     if (abs(totalenergy-evolution(itau-1,1))/totalenergy < sqrtmachprec ) then
        SD_tauend = itau
        call restart( irestart )
        tau = (itau-1) * (SD_tauend - SD_tausta) / SD_Nout ! go backwards, which will cause an ifail=5 error;
        iter = 0 ; return
     endif
  endif
  
  if(mod(itau,SD_SaveFreq) == 0) call restart( irestart )

  itau = itau + 1; tau = itau * (SD_tauend - SD_tausta) / SD_Nout

  iter = 0 !reset iter

  return  

1000 format("progres :"es11.4" : E="es15.7" ; D="es15.7" ; B="es15.7" ; F="es15.7" ; L="es15.7 &
                             " ; A="es15.7" ; C="es15.7" ; H="es15.7)

end subroutine progres
