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
       & case_optimize, DF_maxiter, CG_maxiter, HN_maxiter, TN_maxiter
  implicit none
  include "mpif.h"

  REAL :: start, finish


  iout = 0 ! reset output counter;

  if (myid == 0) write(ounit, *) "-----------OPTIMIZATIONS-------------------------------------"

  if (myid == 0 .and. IsQuiet < 1) write(ounit, '("solvers : total number of DOF is " I6)') Ndof

  if (abs(case_optimize) >= 1) call AllocData(1)
  if (abs(case_optimize) >= 2) call AllocData(2)

  if (case_optimize < 0) then          ! finite difference checking derivatives;
     call fdcheck(case_optimize)
     return
  endif

  if (IsNormWeight /= 0) call normweight
  
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
  use globals, only: zero, sqrtmachprec, myid, ounit, astat, ierr, IsQuiet, &
       Ncoils, deriv, Ndof, xdof, dofnorm, &
       chi, t1E, t2E, &
       bnorm      , t1B, t2B, weight_bnorm,  &
       bharm      , t1H, t2H, weight_bharm,  &
       tflux      , t1F, t2F, weight_tflux, target_tflux, isign, &
       ttlen      , t1L, t2L, weight_ttlen, &
       specw      , t1S, t2S, weight_specw, &
       ccsep      , t1C, t2C, weight_ccsep

  implicit none
  include "mpif.h"

  INTEGER, INTENT(in) :: ideriv

  REAL                :: start, finish
  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if (IsQuiet <= -1) then
     call bnormal(0)
     if (weight_bharm > sqrtmachprec) call bmnharm(0)
!!$
!!$   if ( target_tflux == 0.0 ) then
!!$    call torflux(0)
!!$    target_tflux = isign*sqrt(2.0*tflux)             
!!$    if (myid == 0) write(ounit,'("Costfun :"11X" : Reset target toroidal flux to"ES23.15)') target_tflux
!!$   endif
!!$
!!$   call torflux(0)
!!$
!!$   if (lc == 1) then
!!$    call tlength(0)
!!$   else
!!$    call tlengthExp(0)
!!$   endif
!!$   
!!$   call specwid(0)
!!$   call coilsep(0)
  endif
  
  chi = zero
  
  !if ( ideriv .ge. 1 .and. .not. allocated(t1E) ) ALLOCATE(t1E(1:Ndof))
  !if ( ideriv == 2 .and. .not. allocated(t2E) ) ALLOCATE(t2E(1:Ncoils, 0:Cdof, 1:Ncoils, 0:Cdof))

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

!!$  ! if (myid == 0) write(ounit,'("calling bnormal used",f10.5,"seconds.")') finish-start
!!$
!!$  if (weight_tflux > sqrtmachprec) then
!!$
!!$   if ( target_tflux == 0.0 ) then
!!$    call torflux(0)
!!$    target_tflux = isign*sqrt(2.0*tflux)             
!!$    if (myid == 0) write(ounit,'("Costfun :"11X" : Reset target toroidal flux to"ES23.15)') target_tflux
!!$   endif
!!$
!!$   call torflux(ideriv)
!!$   chi = chi + weight_tflux * tflux / target_tflux**2 ! normalization
!!$   if     ( ideriv == 1 ) then
!!$    t1E = t1E +  weight_tflux * t1F / target_tflux**2
!!$   elseif ( ideriv == 2 ) then
!!$    t1E = t1E +  weight_tflux * t1F / target_tflux**2
!!$    t2E = t2E +  weight_tflux * t2F / target_tflux**2
!!$   endif
!!$
!!$  endif
!!$
!!$  ! if (myid == 0) write(ounit,'("calling torflux used",f10.5,"seconds.")') start-finish
!!$
!!$  if (weight_ttlen > sqrtmachprec) then
!!$
!!$   if ( lc == 1 ) then 
!!$    call tlength(ideriv)
!!$   elseif (lc == 2) then
!!$    call tlengthExp(ideriv)
!!$   else
!!$    FATAL( dnergy, .true., Conflicts between lc and weight_ttlen )
!!$   !call MPI_ABORT( MPI_COMM_WORLD, 1, ierr )
!!$   !stop "COSTFUN: Conflicts between lc and weight_ttlen"
!!$   endif
!!$   
!!$   chi = chi + weight_ttlen * ttlen
!!$   if     ( ideriv == 1 ) then
!!$    t1E = t1E +  weight_ttlen * t1L
!!$   elseif ( ideriv == 2 ) then
!!$    t1E = t1E +  weight_ttlen * t1L
!!$    t2E = t2E +  weight_ttlen * t2L
!!$   endif
!!$
!!$  endif
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

  return
end subroutine costfun

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine normweight
  use globals, only : zero, sqrtmachprec, ounit, myid, xdof, bnorm, bharm, tflux, ttlen, specw, ccsep, &
       weight_bnorm, weight_bharm, weight_tflux, weight_ttlen, weight_specw, weight_ccsep

  implicit none  
  include "mpif.h"

  INTEGER    :: ierr, icoil
  REAL       :: tmp, cur_tflux

  !----------------------------------------------------------------------------------------------------

  call unpacking(xdof)

!-!-!-!-!-!-!-!-!-!-bnorm-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( weight_bnorm >= sqrtmachprec ) then

   call bnormal(0)   
   if (abs(bnorm) > sqrtmachprec) weight_bnorm = weight_bnorm / bnorm
   if( myid == 0 ) write(ounit, 1000) "weight_bnorm", weight_bnorm
   
  endif
!!$
!!$!-!-!-!-!-!-!-!-!-!-tflux-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!!$
!!$  if( weight_tflux .ge. sqrtmachprec ) then
!!$
!!$   if ( target_tflux .eq. 0.0 ) then
!!$    call torflux(0)
!!$    target_tflux = isign*sqrt(2.0*tflux)             
!!$    if(myid .eq. 0) write(ounit,'("costfun :"11X" : Reset target toroidal flux to"ES23.15)') target_tflux
!!$   else
!!$      tmp = target_tflux
!!$      target_tflux = zero
!!$      call torflux(0)
!!$      target_tflux = tmp
!!$      cur_tflux = isign*sqrt(2.0*tflux)
!!$      Io = Io * target_tflux / cur_tflux
!!$      do icoil = 1, Ncoils
!!$         coil(icoil)%I = Io
!!$         coil(icoil)%Io = Io
!!$      enddo
!!$      if(myid .eq. 0) write(ounit,'("costfun :"11X" : rescale coil currents with a factor of"ES23.15)') target_tflux / cur_tflux
!!$   endif
!!$
!!$   call torflux(0)
!!$   if (abs(tflux) .gt. sqrtmachprec) weight_tflux = weight_tflux / tflux * target_tflux**2
!!$   if( myid .eq. 0 ) write(ounit, 1000) "weight_tflux", weight_tflux
!!$   
!!$  endif  
!!$
!!$!-!-!-!-!-!-!-!-!-!-ttlen-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!!$
!!$  if( weight_ttlen .ge. sqrtmachprec ) then
!!$
!!$   if( lc .eq. 1 ) then 
!!$    call tlength(0)
!!$   elseif(lc .eq. 2) then
!!$    call tlengthExp(0)
!!$   else
!!$    FATAL( dnergy, .true., Conflicts between lc and weight_ttlen )
!!$   !call MPI_ABORT( MPI_COMM_WORLD, 1, ierr )
!!$   !stop "Weights_normalize: Conflicts between lc and weight_ttlen"
!!$   endif
!!$
!!$   if (abs(ttlen) .gt. sqrtmachprec) weight_ttlen = weight_ttlen / ttlen
!!$   if( myid .eq. 0 ) write(ounit, 1000) "weight_ttlen", weight_ttlen
!!$   
!!$  endif 
!!$
!!$!-!-!-!-!-!-!-!-!-!-eqarc-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!!$
!!$  if( weight_eqarc .ge. sqrtmachprec ) then
!!$
!!$   call specwid(0)
!!$   if (abs(eqarc) .gt. sqrtmachprec) weight_eqarc = weight_eqarc / eqarc
!!$   if( myid .eq. 0 ) write(ounit, 1000) "weight_eqarc", weight_eqarc
!!$   
!!$  endif 
!!$
!!$!-!-!-!-!-!-!-!-!-!-ccsep-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!!$
!!$  if( weight_ccsep .ge. sqrtmachprec ) then
!!$
!!$   call coilsep(0)
!!$   if (abs(ccsep) .gt. sqrtmachprec) weight_ccsep = weight_ccsep / ccsep
!!$   if( myid .eq. 0 ) write(ounit, 1000) "weight_ccsep", weight_ccsep
!!$   
!!$  endif

1000 format("solvers : "A12" is normalized to" ES23.15)

  return

end subroutine normweight

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine output (mark)

  use globals, only : zero, ounit, myid, ierr, astat, iout, Nouts, Ncoils, save_freq, Tdof, &
       coil, coilspace, FouCoil, chi, t1E, bnorm, bharm, tflux, ttlen, specw, ccsep, evolution, xdof, DoF

  implicit none  
  include "mpif.h"

  REAL, INTENT( IN ) :: mark

  INTEGER            :: idof, NF, icoil
  REAL               :: sumdE

  iout = iout + 1

  FATAL( output , iout > Nouts+1, maximum iteration reached )

  sumdE = sqrt(sum(t1E**2)) ! Eucliean norm 2; 

  if (myid == 0) write(ounit, '("output  : "I6" : "9(ES12.5," ; "))') iout, mark, chi, sumdE, bnorm, bharm, &
       tflux, ttlen, specw, ccsep

  ! save evolution data;
  if (allocated(evolution)) then
     evolution(iout,0) = mark
     evolution(iout,1) = chi
     evolution(iout,2) = sumdE
     evolution(iout,3) = bnorm
     evolution(iout,4) = bharm
     evolution(iout,5) = tflux
     evolution(iout,6) = ttlen
     evolution(iout,7) = specw
     evolution(iout,8) = ccsep
  endif

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


  return  

end subroutine output

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
