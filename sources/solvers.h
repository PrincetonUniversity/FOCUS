!title (solvers) ! Integrated interface for all solvers

!latex \briefly{After initializing the surface and coils, the core part is calling minization algorithms 
!latex to optimize coil parameters.}

!latex \calledby{\link{focus}}
!latex \calls{\link{bnormal}, \link{bmnharm}, \link{torflux}, \link{length}, \link{coilsep}, 
!latex         \link{descent}, \link{congrad}, \link{hybridnt}, \link{truncnt}}

!latex  \section{Cost functions}
!latex  The chi-squared optimization method is used here. The single target function is composed of 
!latex  several chosen object functions with user-supplied weights. The general formula is
!latex  \be
!latex  \ds \chi^2(\vect{X}) = \sum_j w_j \left( \frac{f_j({\vect X}) - f_{j,o}}{f_{j,o}} \right)^2 .
!latex  \ee
!latex  
!latex  Currently, we have implemented constraints on Bnormal, Bmn harmonics, toroidal fulx, coil length, 
!latex  coil-coil separation. For details, please view the docuentation of each constraint.
!latex 
!latex  \section{Normalization}
!latex  Besides the normalization terms in each constraint, like $|{\bf B}|$ in Bnormal, there is also an 
!latex  option to normalize the object function values to its initial value.
!latex  
!latex  When \inputvar{IsNormWeight = 1}, all the nonzero weights will be divided by the current object
!latex  function values. For example, in the beginning, the Bnormal error is $f_{B_0} = 0.1$ and 
!latex  input $w_B = 1.0$. Then the updated $w'_B = w_B/f_{B_0} = 10.0$, such that at every step
!latex  \be
!latex  \ds w'_B f_B = w_B \frac{f_B}{f_{B_0}} \ .
!latex  \ee
!latex  
!latex  \emph{* Please note that when writing the output file, the original weights (as same as input) 
!latex          and {\bf IsNormWeight=1} are stored. So when you restart, the updated weights could be
!latex          different.}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine solvers
  use globals, only : ierr, iout, myid, ounit, IsQuiet, IsNormWeight, Ndof, Nouts, xdof, &
       & case_optimize, DF_maxiter, CG_maxiter, HN_maxiter, TN_maxiter, coil, DoF
  implicit none
  include "mpif.h"

  REAL :: start, finish, ii, lxdof(1:Ndof), f, g(1:Ndof)


  if (myid == 0) write(ounit, *) "-----------OPTIMIZATIONS-------------------------------------"

  if (myid == 0 .and. IsQuiet < 1) write(ounit, '("solvers : total number of DOF is " I6)') Ndof

  if (abs(case_optimize) >= 1) call AllocData(1)
  if (abs(case_optimize) >= 2) call AllocData(2)

  if (case_optimize < 0) then          ! finite difference checking derivatives;
     call fdcheck(case_optimize)
     return
  endif

  if (IsNormWeight /= 0) call normweight

  if (myid == 0 .and. IsQuiet < 0) write(ounit, *) "---------------------------------------------------"
  if (myid == 0 .and. IsQuiet < 0) write(ounit, *) " Initial status"
  if (myid == 0) write(ounit, '("output  : "A6" : "8(A12," ; "))') "iout", "mark", "chi", "dE_norm", &
       "Bnormal", "Bmn harmonics", "tor. flux", "coil length", "c-s sep." 
  call costfun(1)
  call saveBmn    ! in bmnharm.h;
  iout = 0 ! reset output counter;
  call output(0.0)
  
  !--------------------------------DF--------------------------------------------------------------------
  if (DF_maxiter > 0)  then
     if (myid == 0 .and. IsQuiet < 0) write(ounit, *) "---------------------------------------------------"
     if (myid == 0 .and. IsQuiet < 0) write(ounit, *) " Optimizing with Differential Flow (DF)"
     call MPI_BARRIER( MPI_COMM_WORLD, ierr ) ! wait all cpus;
     start = MPI_Wtime()
     call unpacking(xdof)
     call descent
     call packdof(xdof)
     finish = MPI_Wtime()
     if (myid  ==  0) write(ounit,'("solvers : DF takes ", es23.15," seconds;")') finish - start
  endif

  !--------------------------------CG--------------------------------------------------------------------
  if (CG_maxiter > 0)  then
     if (myid == 0 .and. IsQuiet < 0) write(ounit, *) "---------------------------------------------------"
     if (myid == 0 .and. IsQuiet < 0) write(ounit, *) " Optimizing with Nonlinear Conjugate Gradient (CG)"
     call MPI_BARRIER( MPI_COMM_WORLD, ierr ) ! wait all cpus;
     start = MPI_Wtime()
     call unpacking(xdof)
     call congrad
     call packdof(xdof)
     finish = MPI_Wtime()
     if (myid  ==  0) write(ounit,'("solvers : CG takes ", es23.15," seconds;")') finish - start
  endif
  
  !--------------------------------HN--------------------------------------------------------------------
  if (HN_maxiter > 0)  then
     if (myid == 0 .and. IsQuiet < 0) write(ounit, *) "---------------------------------------------------"
     if (myid == 0 .and. IsQuiet < 0) write(ounit, *) " Optimizing with Hybrid Newton Method (HN) "
     call MPI_BARRIER( MPI_COMM_WORLD, ierr ) ! wait all cpus;
     start = MPI_Wtime()
     call unpacking(xdof)
     !call hybridnt
     finish = MPI_Wtime()
     if (myid  ==  0) write(ounit,'("solvers : HN takes ", es23.15," seconds;")') finish - start
  endif
  
  !--------------------------------TN--------------------------------------------------------------------
  if (TN_maxiter > 0)  then
     if (myid == 0 .and. IsQuiet < 0) write(ounit, *) "---------------------------------------------------"
     if (myid == 0 .and. IsQuiet < 0) write(ounit, *) " Optimizing with Truncated Newton Method (TN) "
     call MPI_BARRIER( MPI_COMM_WORLD, ierr ) ! wait all cpus;
     start = MPI_Wtime()
     call unpacking(xdof)
     !call truncnt
     call packdof(xdof)
     finish = MPI_Wtime()
     if (myid  ==  0) write(ounit,'("solvers : TN takes ", es23.15," seconds;")') finish - start
  endif
  
  !------------------------------------------------------------------------------------------------------
  
  return
end subroutine solvers

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine costfun(ideriv)
  use globals, only: zero, one, sqrtmachprec, myid, ounit, astat, ierr, IsQuiet, &
       Ncoils, deriv, Ndof, xdof, dofnorm, coil, &
       chi, t1E, t2E, &
       bnorm      , t1B, t2B, weight_bnorm,  &
       bharm      , t1H, t2H, weight_bharm,  &
       tflux      , t1F, t2F, weight_tflux, target_tflux, psi_avg, &
       ttlen      , t1L, t2L, weight_ttlen, case_length, &
       cssep      , t1S, t2S, weight_cssep, &
       specw      , t1P, t2P, weight_specw, &
       ccsep      , t1C, t2C, weight_ccsep

  implicit none
  include "mpif.h"

  INTEGER, INTENT(in) :: ideriv

  REAL                :: start, finish
  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  call mpi_barrier(MPI_COMM_WORLD, ierr)

  if (IsQuiet <= -2) then

     call bnormal(0)
     if (weight_bharm > sqrtmachprec) call bmnharm(0)

     if ( abs(target_tflux) < sqrtmachprec ) then
        call torflux(0)
        target_tflux = psi_avg
        if (myid == 0) write(ounit,'("costfun : Reset target toroidal flux to "ES23.15)') target_tflux
     endif

     call torflux(0)

     if ( (case_length == 1) .and. (sum(coil(1:Ncoils)%Lo) < sqrtmachprec) ) then
        coil(1:Ncoils)%Lo = one
        call length(0)
        coil(1:Ncoils)%Lo = coil(1:Ncoils)%L
        if(myid .eq. 0) write(ounit,'("solvers : reset target coil length to the current actual length. ")')
     endif

     call length(0)

  endif

  chi = zero


  if ( ideriv == 1 ) then
     t1E = zero
  elseif ( ideriv == 2 ) then
     t1E = zero; t2E = zero
  endif

  !call unpacking(xdof)
  ! Bnormal surface intergration;
  if (weight_bnorm > sqrtmachprec) then
 
     call bnormal(ideriv)
     chi = chi + weight_bnorm * bnorm
     if     ( ideriv == 1 ) then
        t1E = t1E +  weight_bnorm * t1B
     elseif ( ideriv == 2 ) then
        t1E = t1E +  weight_bnorm * t1B
        t2E = t2E +  weight_bnorm * t2B
     endif
  endif

  ! individual Bn harmonics;
  if (weight_bharm > sqrtmachprec) then

     call bmnharm(ideriv)
     chi = chi + weight_bharm * bharm
     if     ( ideriv == 1 ) then
        t1E = t1E +  weight_bharm * t1H
     elseif ( ideriv == 2 ) then
        t1E = t1E +  weight_bharm * t1H
        t2E = t2E +  weight_bharm * t2H
     endif
  endif

  ! toroidal flux;
  if (weight_tflux > sqrtmachprec) then
     if ( abs(target_tflux) < sqrtmachprec ) then
        call torflux(0)
        target_tflux = psi_avg        
        if (myid==0) write(ounit,'("solvers : Reset target toroidal flux to "ES23.15)') target_tflux
     endif

     call torflux(ideriv)
     chi = chi + weight_tflux * tflux / target_tflux**2 ! normalization
     if     ( ideriv == 1 ) then
        t1E = t1E +  weight_tflux * t1F / target_tflux**2
     elseif ( ideriv == 2 ) then
        t1E = t1E +  weight_tflux * t1F / target_tflux**2
        t2E = t2E +  weight_tflux * t2F / target_tflux**2
     endif
  endif

  ! coil length;
  if (weight_ttlen > sqrtmachprec) then

     if ( (case_length == 1) .and. (sum(coil(1:Ncoils)%Lo) < sqrtmachprec) ) then
        coil(1:Ncoils)%Lo = one
        call length(0)
        coil(1:Ncoils)%Lo = coil(1:Ncoils)%L
        if(myid .eq. 0) write(ounit,'("solvers : reset target coil length to the current actual length. ")')
     endif

     call length(ideriv)

     chi = chi + weight_ttlen * ttlen
     if     ( ideriv == 1 ) then
        t1E = t1E +  weight_ttlen * t1L
     elseif ( ideriv == 2 ) then
        t1E = t1E +  weight_ttlen * t1L
        t2E = t2E +  weight_ttlen * t2L
     endif

  endif

  ! coil surface separation;
  if (weight_cssep > sqrtmachprec) then
 
     call surfsep(ideriv)
     chi = chi + weight_cssep * cssep
     if     ( ideriv == 1 ) then
        t1E = t1E +  weight_cssep * t1S
     elseif ( ideriv == 2 ) then
        t1E = t1E +  weight_cssep * t1S
        t2E = t2E +  weight_cssep * t2S
     endif
  endif

!!$
!!$  ! if (myid == 0) write(ounit,'("calling tlength used",f10.5,"seconds.")') finish-start
!!$
!!$  if (weight_eqarc .ge. sqrtmachprec) then
!!$
!!$  ! call equarcl(ideriv)
!!$  ! call specwid_df(ideriv)
!!$   if ( Loptimize == 3) then
!!$    call langrange(ideriv)
!!$   else
!!$    !call specwid(ideriv)
!!$    call specwid_df(ideriv)
!!$   endif
!!$   
!!$   chi = chi + weight_eqarc * eqarc
!!$   if     ( ideriv == 1 ) then
!!$    t1E = t1E + weight_eqarc * t1A
!!$   elseif ( ideriv == 2 ) then
!!$    t1E = t1E + weight_eqarc * t1A
!!$    t2E = t2E + weight_eqarc * t2A
!!$   endif
!!$
!!$  endif
!!$
!!$  if (weight_ccsep .ge. sqrtmachprec) then
!!$
!!$   call coilsep(ideriv)
!!$   chi = chi + weight_ccsep * ccsep
!!$   if     ( ideriv == 1 ) then
!!$    t1E = t1E + weight_ccsep * t1C
!!$   elseif ( ideriv == 2 ) then
!!$    t1E = t1E + weight_ccsep * t1C
!!$    t2E = t2E + weight_ccsep * t2C
!!$   endif
!!$
!!$  endif
!!$
!!$  if (allocated(deriv)) then
!!$     deriv = zero
!!$     do ii = 1, Ndof
!!$        call DoFconvert(ii, icl, inf)
!!$        if (allocated(t1E)) deriv(ii, 0) = t1E(icl, inf)
!!$        if (allocated(t1B)) deriv(ii, 1) = t1B(icl, inf)
!!$        if (allocated(t1F)) deriv(ii, 2) = t1F(icl, inf)
!!$        if (allocated(t1L)) deriv(ii, 3) = t1L(icl, inf)
!!$        if (allocated(t1A)) deriv(ii, 4) = t1A(icl, inf)
!!$        if (allocated(t1C)) deriv(ii, 5) = t1C(icl, inf)
!!$     enddo
!!$  endif

  if ( ideriv == 1 ) then                ! multiply t1E & t2E with normalized terms; 06/09/2017
     t1E = t1E * dofnorm
!!$   elseif ( ideriv == 2 ) then
!!$      t1E = t1E * dofnorm
!!$      t2E = t2E * hesnorm
  endif

  call mpi_barrier(MPI_COMM_WORLD, ierr)

  return
end subroutine costfun

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine normweight
  use globals, only : zero, one, machprec, ounit, myid, xdof, bnorm, bharm, tflux, ttlen, cssep, specw, ccsep, &
       weight_bnorm, weight_bharm, weight_tflux, weight_ttlen, weight_cssep, weight_specw, weight_ccsep, &
       target_tflux, psi_avg, coil, Ncoils, case_length, Bmnc, Bmns, tBmnc, tBmns

  implicit none  
  include "mpif.h"

  INTEGER    :: ierr, icoil
  REAL       :: tmp, cur_tflux, modBn, modtBn

  !----------------------------------------------------------------------------------------------------

  call unpacking(xdof)

  !-!-!-!-!-!-!-!-!-!-tflux-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( weight_tflux .ge. machprec ) then

     call torflux(0)
     if ( abs(target_tflux) < machprec ) then
        target_tflux = psi_avg
        if(myid .eq. 0) write(ounit,'("solvers : Reset target toroidal flux to "ES12.5)') target_tflux
     else if (sum(abs(coil(1:Ncoils)%Ic)) == Ncoils) then !only valid when all currents are free;
           do icoil = 1, Ncoils
              coil(icoil)%I = coil(icoil)%I * target_tflux / psi_avg
           enddo
           if(myid .eq. 0) write(ounit,'("solvers : rescale coil currents with a factor of "ES12.5)') &
                target_tflux / psi_avg
     endif

     call torflux(0)
     if (abs(tflux) > machprec) weight_tflux = weight_tflux / tflux * target_tflux**2
     if( myid .eq. 0 ) write(ounit, 1000) "weight_tflux", weight_tflux

  endif

  !-!-!-!-!-!-!-!-!-!-bharm-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( weight_bharm >= machprec ) then

     call bmnharm(0)
     modBn = sqrt(sum(Bmnc**2 + Bmns**2))
     modtBn = sqrt(sum(tBmnc**2 + tBmns**2))
     do icoil = 1, Ncoils
        coil(icoil)%I = coil(icoil)%I * modtBn / modBn
     enddo
     if(myid .eq. 0) write(ounit,'("solvers : rescale coil currents with a factor of "ES12.5)') &
          modtBn / modBn
     call bmnharm(0)
     if (abs(bharm) > machprec) weight_bharm = weight_bharm / bharm
     if( myid == 0 ) write(ounit, 1000) "weight_bharm", weight_bharm

  endif

  !-!-!-!-!-!-!-!-!-!-bnorm-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( weight_bnorm >= machprec ) then

     call bnormal(0)   
     if (abs(bnorm) > machprec) weight_bnorm = weight_bnorm / bnorm
     if( myid == 0 ) write(ounit, 1000) "weight_bnorm", weight_bnorm

  endif

  !-!-!-!-!-!-!-!-!-!-ttlen-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( weight_ttlen .ge. machprec ) then

     if ( (case_length == 1) .and. (sum(coil(1:Ncoils)%Lo) < machprec) ) then
        coil(1:Ncoils)%Lo = one
        call length(0)
        coil(1:Ncoils)%Lo = coil(1:Ncoils)%L
        if(myid .eq. 0) write(ounit,'("solvers : reset target coil length to the current actual length. ")')
     endif

     call length(0)

     if (abs(ttlen) .gt. machprec) weight_ttlen = weight_ttlen / ttlen
     if( myid .eq. 0 ) write(ounit, 1000) "weight_ttlen", weight_ttlen

  endif

  !-!-!-!-!-!-!-!-!-!-cssep-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( weight_cssep >= machprec ) then

     call surfsep(0)   
     if (abs(cssep) > machprec) weight_cssep = weight_cssep / cssep
     if( myid == 0 ) write(ounit, 1000) "weight_cssep", weight_cssep

  endif

!!$!-!-!-!-!-!-!-!-!-!-eqarc-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!!$
!!$  if( weight_eqarc .ge. machprec ) then
!!$
!!$   call specwid(0)
!!$   if (abs(eqarc) .gt. machprec) weight_eqarc = weight_eqarc / eqarc
!!$   if( myid .eq. 0 ) write(ounit, 1000) "weight_eqarc", weight_eqarc
!!$   
!!$  endif 
!!$
!!$!-!-!-!-!-!-!-!-!-!-ccsep-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!!$
!!$  if( weight_ccsep .ge. machprec ) then
!!$
!!$   call coilsep(0)
!!$   if (abs(ccsep) .gt. machprec) weight_ccsep = weight_ccsep / ccsep
!!$   if( myid .eq. 0 ) write(ounit, 1000) "weight_ccsep", weight_ccsep
!!$   
!!$  endif

1000 format("solvers : "A12" is normalized to " ES23.15)

  call packdof(xdof)

  return

end subroutine normweight

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine output (mark)

  use globals, only : zero, ounit, myid, ierr, astat, iout, Nouts, Ncoils, save_freq, Tdof, &
       coil, coilspace, FouCoil, chi, t1E, bnorm, bharm, tflux, ttlen, cssep, specw, ccsep, &
       evolution, xdof, DoF, exit_tol, exit_signal

  implicit none  
  include "mpif.h"

  REAL, INTENT( IN ) :: mark

  INTEGER            :: idof, NF, icoil
  REAL               :: sumdE


  iout = iout + 1
  
  FATAL( output , iout > Nouts+2, maximum iteration reached )

  sumdE = sqrt(sum(t1E**2)) ! Eucliean norm 2; 

  if (myid == 0) write(ounit, '("output  : "I6" : "8(ES12.5," ; "))') iout, mark, chi, sumdE, bnorm, bharm, &
       tflux, ttlen, cssep

  ! save evolution data;
  if (allocated(evolution)) then
     evolution(iout,0) = mark
     evolution(iout,1) = chi
     evolution(iout,2) = sumdE
     evolution(iout,3) = bnorm
     evolution(iout,4) = bharm
     evolution(iout,5) = tflux
     evolution(iout,6) = ttlen
     evolution(iout,7) = cssep
     !evolution(iout,8) = 0.0
     !evolution(iout,8) = ccsep
  endif

  ! exit the optimization if no obvious changes in past 5 outputs; 07/20/2017
  if (iout>5) then
     if ( abs(evolution(iout,1) - evolution(iout-5, 1)) / evolution(iout,1) < exit_tol ) exit_signal = .True.
  end if
  
  !save all the coil parameters;
  if (allocated(coilspace)) then
     idof = 0
     do icoil = 1, Ncoils
        coilspace(iout, idof+1 ) = coil(icoil)%I ;  idof = idof + 1

        select case (coil(icoil)%itype)
        case (1)
           NF = FouCoil(icoil)%NF
           coilspace(iout, idof+1:idof+NF+1) = FouCoil(icoil)%xc(0:NF) ; idof = idof + NF +1
           coilspace(iout, idof+1:idof+NF  ) = FouCoil(icoil)%xs(1:NF) ; idof = idof + NF
           coilspace(iout, idof+1:idof+NF+1) = FouCoil(icoil)%yc(0:NF) ; idof = idof + NF +1
           coilspace(iout, idof+1:idof+NF  ) = FouCoil(icoil)%ys(1:NF) ; idof = idof + NF
           coilspace(iout, idof+1:idof+NF+1) = FouCoil(icoil)%zc(0:NF) ; idof = idof + NF +1
           coilspace(iout, idof+1:idof+NF  ) = FouCoil(icoil)%zs(1:NF) ; idof = idof + NF
        case default
           FATAL(descent, .true., not supported coil types)
        end select
     enddo
     FATAL( output , idof .ne. Tdof, counting error in restart )
  endif

  if(mod(iout,save_freq) .eq. 0) call saving

  call MPI_BARRIER( MPI_COMM_WORLD, ierr ) ! wait all cpus;

  return  

end subroutine output

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
