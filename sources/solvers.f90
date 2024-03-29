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
  use globals
  use mpi
  implicit none

  REAL :: start, finish


  if (myid == 0) write(ounit, *) "-----------OPTIMIZATIONS-------------------------------------"
  if (myid == 0) write(ounit, '("solvers : Total number of DOF is " I6)') Ndof
  if (myid == 0 .and. IsQuiet < 1) then
     write(ounit, '(8X,": Initial weights are: ")')
     ! Bnormal
     if (weight_bnorm > machprec) then
        write(ounit, '(8X,": weight_bnorm = ", ES12.5, "; case_bnormal = ", I1)') weight_bnorm, case_bnormal
     endif 
     ! Bmn harmonics
     if (weight_bharm > machprec) then 
        write(ounit, '(8X,": weight_bharm = ", ES12.5, "; bharm_jsurf = ", I1)') weight_bharm, bharm_jsurf
     endif       
     ! toroidal flux
     if (weight_tflux > machprec) then 
        write(ounit, '(8X,": weight_tflux = ", ES12.5, "; target_tflux = ", ES12.5)') weight_tflux, target_tflux
     endif      
     ! total current
     if (weight_isum > machprec) then 
        write(ounit, '(8X,": weight_isum  = ", ES12.5, "; target_isum  = ", ES12.5, " Ampere")') weight_isum, target_isum
     endif          
     ! coil length
     if (weight_ttlen > machprec) then 
        write(ounit, '(8X,": weight_ttlen    = ", ES12.5, "; case_length  = ", I1)') weight_ttlen, case_length
        write(ounit, '(8X,":   target_length = ", ES12.5, "; length_delta = ", ES12.5)') target_length, length_delta
     endif
     ! curvature 
     if (weight_curv > machprec) then 
        write(ounit, '(8X,": weight_curv  = ", ES12.5, ";  case_curv = ", I1, &
        & "; penfun_curv = ", I1)') weight_curv, case_curv, penfun_curv
        write(ounit, '(8X,":   curv_alpha = ", ES12.5, "; curv_beta  = ", ES12.5)') curv_alpha, curv_beta
        write(ounit, '(8X,":   curv_gamma = ", ES12.5, "; curv_sigma = ", ES12.5)')  curv_gamma, curv_sigma
        write(ounit, '(8X,":   curv_k0    = ", ES12.5, "; curv_k1    = ", ES12.5, &
        & "; curv_k1len = ", I1)') curv_k0, curv_k1, curv_k1len
     endif 
     ! torsion
     if (weight_tors > machprec) then 
        write(ounit, '(8X,": weight_tors  = ", ES12.5, "; penfun_tors = ", I1)') weight_tors, penfun_tors
        write(ounit, '(8X,":   tors_alpha = ", ES12.5, "; tors_beta   = ", ES12.5)') tors_alpha, tors_beta
        write(ounit, '(8X,":   tors_gamma = ", ES12.5, "; tors0       = ", ES12.5)') tors_gamma, tors0
     endif 
     ! Nissin complexicity 
     if (weight_nissin > 0.0_dp) then 
        write(ounit, '(8X,": weight_nissin  = ", ES12.5, ";  penfun_nissin = ", I1, &
        & "; nissin0 = ", ES12.5)') weight_nissin, penfun_nissin, nissin0
        write(ounit, '(8X,":   nissin_alpha = ", ES12.5, "; nissin_beta    = ", ES12.5)') nissin_alpha, nissin_beta
        write(ounit, '(8X,":   nissin_gamma = ", ES12.5, "; nissin_sigma   = ", ES12.5)')  nissin_gamma, nissin_sigma
     endif 
     ! coil-surface separation
     if (weight_cssep > machprec) then 
        write(ounit, '(8X,": weight_cssep = ", ES12.5, "; cssep_factor = ", ES12.5)') weight_cssep, cssep_factor
     endif 
     ! coil-coil separation
     if (weight_ccsep > machprec) then 
        write(ounit, '(8X,": weight_ccsep  = ", ES12.5, "; penfun_ccsep = ", I1 &
        & "ccsep_skip = ", I1)') weight_ccsep, penfun_ccsep, ccsep_skip
        write(ounit, '(8X,":   ccsep_alpha = ", ES12.5, "; ccsep_beta   = ", ES12.5, &
        & "; r_delta = ", ES12.5)') ccsep_alpha, ccsep_beta, r_delta
     endif
     ! stochastic bnorm
     if (weight_sbnorm > machprec) then
        write(ounit, '(8X,": weight_sbnorm = ", ES12.5, "; Npert = ", I10)') weight_sbnorm, Npert
     endif 
     ! straight-out
     if (weight_straight > machprec) then 
        write(ounit, '(8X,": weight_straight  = ", ES12.5, ";  case_straight = ", I1)') weight_straight, case_straight
        write(ounit, '(8X,":   str_alpha = ", ES12.5)') str_alpha
        write(ounit, '(8X,":   str_k0    = ", ES12.5)') str_k0
     endif
  endif

  ! Call stochastic construction if neccessary
  if ( weight_sbnorm .gt. 0.0 .or. Npert .gt. 0 ) call perturbation(0)

  if (abs(case_optimize) >= 1) call AllocData(1)
  if (abs(case_optimize) >= 2) call AllocData(2)

  ! evaluate the initial coils, in case coils intersect with plasma
  if (myid==0) write(ounit, *) "------------- Initial status ------------------------"
  call diagnos

  if (IsNormWeight /= 0) call normweight

  if (case_optimize < 0) then          ! finite difference checking derivatives;
     call fdcheck(case_optimize)
     return
  endif

  if (myid ==0) then
     write(ounit, *) "------------- Iterations ------------------------"
  endif 
  
  if (lim_out .eq. 0) then
     if (myid == 0) then
        if (IsQuiet < 0) then
           write(ounit, '("output  : "A6" : "14(A12," ; "))') "iout", "mark", "chi", "dE_norm", &
                "Bnormal", "Bmn harmonics", "tor. flux", "coil length", "c-s sep.", "curvature", &
                "c-c sep.", "torsion", "nissin", "sBnormal", "straight"
        else
           write(ounit, '("output  : "A6" : "4(A15," ; "))') "iout", "mark", "total function", "||gradient||", "Bnormal"
        endif
     endif
  else
     if (myid .eq. 0) then
        write(ounit, '("output  : "A6" : "3(A12," ; "))',advance="no") "iout", "mark", "chi", "dE_norm"
        if ( weight_bnorm .ne. 0.0 ) then
           write(ounit, '(A12," ; ")',advance="no") "Bnormal"
        endif
        if ( weight_bharm .ne. 0.0 ) then
           write(ounit, '(A12," ; ")',advance="no") "Bmn harmonics"
        endif
        if ( weight_tflux .ne. 0.0 ) then
           write(ounit, '(A12," ; ")',advance="no") "tor. flux"
        endif
        if ( weight_isum .ne. 0.0 ) then
           write(ounit, '(A12," ; ")',advance="no") "total current"
        endif        
        if ( weight_ttlen .ne. 0.0 ) then
           write(ounit, '(A12," ; ")',advance="no") "coil length"
        endif
        if ( weight_cssep .ne. 0.0 ) then
           write(ounit, '(A12," ; ")',advance="no") "c-s sep."
        endif
        if ( weight_curv .ne. 0.0 ) then
           write(ounit, '(A12," ; ")',advance="no") "curvature"
        endif
        if ( weight_ccsep .ne. 0.0 ) then
           write(ounit, '(A12," ; ")',advance="no") "c-c sep."
        endif
        if ( weight_tors .ne. 0.0 ) then
           write(ounit, '(A12," ; ")',advance="no") "torsion"
        endif
        if ( weight_nissin .ne. 0.0 ) then
           write(ounit, '(A12," ; ")',advance="no") "nissin"
        endif
        if ( weight_sbnorm .ne. 0.0 ) then
           write(ounit, '(A12," ; ")',advance="no") "sBnormal"
        endif
        if ( weight_straight .ne. 0.0 ) then
            write(ounit, '(A12," ; ")',advance="no") "straight"
        endif
        write(*,*)
     endif
  endif
  
  call costfun(1)
  call saveBmn    ! in bmnharm.h;
  iout = 0 ! reset output counter;
  call output(zero)
  
  !--------------------------------DF--------------------------------------------------------------------
  if (DF_maxiter > 0)  then
     if (myid == 0 .and. IsQuiet < 0) write(ounit, *) "------------- Differential Flow (DF) ----------------"
     call MPI_BARRIER( MPI_COMM_FOCUS, ierr ) ! wait all cpus;
     start = MPI_Wtime()
     call unpacking(xdof)
     call descent
     call packdof(xdof)
     finish = MPI_Wtime()
     if (myid  ==  0) write(ounit,'(8X,": DF takes ", es23.15," seconds;")') finish - start
  endif

  !--------------------------------CG--------------------------------------------------------------------
  if (CG_maxiter > 0)  then
     if (myid == 0 .and. IsQuiet < 0) write(ounit, *) "------------- Nonlinear Conjugate Gradient (CG) -----"
     call MPI_BARRIER( MPI_COMM_FOCUS, ierr ) ! wait all cpus;
     start = MPI_Wtime()
     call unpacking(xdof)
     call congrad
     call packdof(xdof)
     finish = MPI_Wtime()
     if (myid  ==  0) write(ounit,'(8X,": CG takes ", es23.15," seconds;")') finish - start
  endif

  !--------------------------------LM--------------------------------------------------------------------
  if (LM_maxiter > 0)  then
     if (myid == 0 .and. IsQuiet < 0) write(ounit, *) "----------- Levenberg-Marquardt algorithm (L-M) -----"
     call MPI_BARRIER( MPI_COMM_FOCUS, ierr ) ! wait all cpus;
     start = MPI_Wtime()
     call unpacking(xdof)
     call lmalg
     call packdof(xdof)
     finish = MPI_Wtime()
     if (myid  ==  0) write(ounit,'(8X,": LM takes ", es23.15," seconds;")') finish - start
  endif
  
  !--------------------------------HN--------------------------------------------------------------------
  if (HN_maxiter > 0)  then
     if (myid == 0 .and. IsQuiet < 0) write(ounit, *) "---------------------------------------------------"
     if (myid == 0 .and. IsQuiet < 0) write(ounit, *) " Optimizing with Hybrid Newton Method (HN) "
     call MPI_BARRIER( MPI_COMM_FOCUS, ierr ) ! wait all cpus;
     start = MPI_Wtime()
     call unpacking(xdof)
     !call hybridnt
     finish = MPI_Wtime()
     if (myid  ==  0) write(ounit,'(8X,": HN takes ", es23.15," seconds;")') finish - start
  endif
  
  !--------------------------------TN--------------------------------------------------------------------
  if (TN_maxiter > 0)  then
     if (myid == 0 .and. IsQuiet < 0) write(ounit, *) "---------------------------------------------------"
     if (myid == 0 .and. IsQuiet < 0) write(ounit, *) " Optimizing with Truncated Newton Method (TN) "
     call MPI_BARRIER( MPI_COMM_FOCUS, ierr ) ! wait all cpus;
     start = MPI_Wtime()
     call unpacking(xdof)
     !call truncnt
     call packdof(xdof)
     finish = MPI_Wtime()
     if (myid  ==  0) write(ounit,'(8X,": TN takes ", es23.15," seconds;")') finish - start
  endif
  
  !------------------------------------------------------------------------------------------------------
  
  return
end subroutine solvers

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine costfun(ideriv)
  use globals, only: dp, zero, one, machprec, myid, ounit, astat, ierr, IsQuiet, &
       Ncoils, deriv, Ndof, xdof, dofnorm, coil, &
       chi, t1E, t2E, LM_maxiter, LM_fjac, LM_mfvec, sumdE, LM_output, LM_fvec, curv_k0, &
       bnorm      , t1B, t2B, weight_bnorm,  &
       bharm      , t1H, t2H, weight_bharm,  &
       tflux      , t1F, t2F, weight_tflux, target_tflux, psi_avg, &
       isum       , t1I, target_isum, weight_isum, &
       ttlen      , t1L, t2L, weight_ttlen, case_length, &
       cssep      , t1S, t2S, weight_cssep, &
       specw      , t1P, t2P, weight_specw, &
       ccsep      , t1C, t2C, weight_ccsep, &
       tors       , t1T, t2T, weight_tors, &
       nissin     , t1N, t2N, weight_nissin, &
       curv       , t1K, t2K, weight_curv, &
       str        , t1Str,t2Str,weight_straight, &
       bnormavg   , t1Bavg, weight_sbnorm, & 
       MPI_COMM_FOCUS

  use mpi
  implicit none

  INTEGER, INTENT(in) :: ideriv
  
  INTEGER             :: ivec
  REAL                :: start, finish
  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  call mpi_barrier(MPI_COMM_FOCUS, ierr)
  
  bnorm = zero 
  bharm = zero
  tflux = zero
  ttlen = zero
  cssep = zero
  specw = zero
  ccsep = zero
  curv  = zero
  str   = zero
  tors  = zero
  nissin   = zero
  bnormavg = zero

  if (IsQuiet <= -2) then

     call bnormal(0)

     if ( abs(target_tflux) < machprec ) then
        call torflux(0)
        target_tflux = psi_avg
        if (myid == 0) write(ounit,'("costfun : Reset target toroidal flux to "ES23.15)') target_tflux
     endif

     call torflux(0)

     if ( (case_length == 1) .and. (sum(coil(1:Ncoils)%Lo) < machprec) ) then
        coil(1:Ncoils)%Lo = one
        call length(0)
        coil(1:Ncoils)%Lo = coil(1:Ncoils)%L
        if(myid .eq. 0) write(ounit,'(8X,": reset target coil length to the current actual length. ")')
     endif

     call length(0)
     call surfsep(0)
     call curvature(0)
     call straight(0)   
     call coilsep(0)
     call torsion(0)
     call nisscom(0)
     call sbnormal(0)

  endif

  chi = zero
  if ( ideriv == 1 ) then
     t1E = zero
  elseif ( ideriv == 2 ) then
     t1E = zero; t2E = zero
  endif

!if(myid==0) write(ounit, '("-------i=1, chi = "ES23.15)') chi

  ! Bn cost functions
  if (weight_bnorm > machprec .or. weight_bharm > machprec) then
 
     call bnormal(ideriv)
     ! Bnormal surface intergration;
     if (weight_bnorm > machprec) then
        chi = chi + weight_bnorm * bnorm
        if     ( ideriv == 1 ) then
           t1E = t1E +  weight_bnorm * t1B
        elseif ( ideriv == 2 ) then
           t1E = t1E +  weight_bnorm * t1B
           t2E = t2E +  weight_bnorm * t2B
        endif
     endif

     ! individual Bn harmonics;
     if (weight_bharm > machprec) then

        chi = chi + weight_bharm * bharm
        if     ( ideriv == 1 ) then
           t1E = t1E +  weight_bharm * t1H
        elseif ( ideriv == 2 ) then
           t1E = t1E +  weight_bharm * t1H
           t2E = t2E +  weight_bharm * t2H
        endif
     endif
  endif

!if(myid==0) write(ounit, '("-------i=2, chi = "ES23.15)') chi

  ! toroidal flux;
  if (weight_tflux > machprec) then
     if ( abs(target_tflux) < machprec ) then
        call torflux(0)
        target_tflux = psi_avg        
        if (myid==0) write(ounit,'(8X,": Reset target toroidal flux to "ES23.15)') target_tflux
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

  ! total current;
  if (weight_isum > machprec) then

     call current(ideriv)
     chi = chi + weight_isum * isum ! normalization
     if     ( ideriv == 1 ) then
        t1E = t1E +  weight_isum * t1I 
     elseif ( ideriv == 2 ) then
        t1E = t1E +  weight_isum * t1I 
      ! t2E = t2E +  weight_isum * t2I 
     endif
  endif

!if(myid==0) write(ounit, '("-------i=3, chi = "ES23.15)') chi

  ! coil length;
  if (weight_ttlen > machprec) then

     if ( (case_length == 1) .and. (sum(coil(1:Ncoils)%Lo) < machprec) ) then
        coil(1:Ncoils)%Lo = one
        call length(0)
        coil(1:Ncoils)%Lo = coil(1:Ncoils)%L
        if(myid .eq. 0) write(ounit,'(8X,": reset target coil length to the current actual length. ")')
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
  !if(myid==0) write(ounit, '("-------i=4, chi = "ES23.15)') chi
  ! coil curvature
  if (weight_curv > machprec) then

     call curvature(ideriv)
     chi = chi + weight_curv * curv
     if     ( ideriv == 1 ) then
        t1E = t1E +  weight_curv * t1K
     elseif ( ideriv == 2 ) then
        t1E = t1E +  weight_curv * t1K
        t2E = t2E +  weight_curv * t2K
     endif
  endif
  !if(myid==0) write(ounit, '("-------i=5, chi = "ES23.15)') chi
  ! straight-out coil
  if (weight_straight > machprec) then

     call straight(ideriv)
     chi = chi + weight_straight * str
     if     ( ideriv == 1 ) then
        t1E = t1E +  weight_straight * t1Str
     elseif ( ideriv == 2 ) then
        t1E = t1E +  weight_straight * t1Str
        t2E = t2E +  weight_straight * t2Str
     endif
   endif  
   !if(myid==0) write(ounit, '("-------i=6, chi = "ES23.15)') chi
  ! coil surface separation;
  if (weight_cssep > machprec) then
 
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
!!$  if (weight_eqarc .ge. machprec) then
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
  !if(myid==0) write(ounit, '("-------i=7, chi = "ES23.15)') chi
  ! coil to coil separation 
  if (weight_ccsep > machprec) then

     call coilsep(ideriv)
     chi = chi + weight_ccsep * ccsep
     if     ( ideriv == 1 ) then
        t1E = t1E + weight_ccsep * t1C
     elseif ( ideriv == 2 ) then
        t1E = t1E + weight_ccsep * t1C
        t2E = t2E + weight_ccsep * t2C
     endif
  endif
  !if(myid==0) write(ounit, '("-------i=8, chi = "ES23.15)') chi
  ! torsion
  if (weight_tors > machprec) then

     call torsion(ideriv)
     chi = chi + weight_tors * tors
     if     ( ideriv == 1 ) then
        t1E = t1E + weight_tors * t1T
     elseif ( ideriv == 2 ) then
        t1E = t1E + weight_tors * t1T
        t2E = t2E + weight_tors * t2T
     endif
  endif
  !if(myid==0) write(ounit, '("-------i=9, chi = "ES23.15)') chi
  ! nissin
  if (weight_nissin > 0.0_dp) then
     call nisscom(ideriv)
     chi = chi + weight_nissin * nissin
     if     ( ideriv == 1 ) then
        t1E = t1E + weight_nissin * t1N
     elseif ( ideriv == 2 ) then
        t1E = t1E + weight_nissin * t1N
        t2E = t2E + weight_nissin * t2N
     endif
  endif

  ! stochastic bnorm
  if (weight_sbnorm > 0.0_dp) then
  
     call sbnormal(ideriv) 
     chi = chi + weight_sbnorm * bnormavg
     if     ( ideriv == 1 ) then
        t1E = t1E + weight_sbnorm * t1Bavg
     !elseif ( ideriv == 2 ) then
     !   t1E = t1E + weight_nissin * t1N
     !   t2E = t2E + weight_nissin * t2N
     endif
  endif

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

     ! L-M format
     if (LM_mfvec > 0) then
        do ivec = 1, LM_mfvec
           LM_fjac(ivec, 1:Ndof) = LM_fjac(ivec, 1:Ndof) * dofnorm(1:Ndof)
        enddo
     endif
!!$   elseif ( ideriv == 2 ) then
!!$      t1E = t1E * dofnorm
!!$      t2E = t2E * hesnorm
  endif
  
  if (.not. LM_output) then ! regular output
     if ( ideriv == 1 )  sumdE = sqrt(sum(t1E**2)) ! Eucliean norm 2;
  else                      ! L-M output
     chi = sum(LM_fvec**2)
     if ( ideriv == 1 ) sumDE = sqrt(sum(LM_fjac**2))
  endif  
  
  call mpi_barrier(MPI_COMM_FOCUS, ierr)

  return
end subroutine costfun

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine normweight
  use globals, only: dp, zero, one, machprec, ounit, myid, xdof, bnorm, bharm, tflux, ttlen, cssep, specw, ccsep, &
       weight_bnorm, weight_bharm, weight_tflux, weight_ttlen, weight_cssep, weight_specw, weight_ccsep, &
       target_tflux, psi_avg, coil, Ncoils, case_length, Bmnc, Bmns, tBmnc, tBmns, weight_curv, curv, &
       weight_tors, tors, weight_nissin, nissin, str, weight_straight, weight_sbnorm, bnormavg, & 
       isum, weight_isum, MPI_COMM_FOCUS

  use mpi
  implicit none

  INTEGER    :: ierr, icoil
  REAL       :: tmp, cur_tflux, modBn, modtBn

  !----------------------------------------------------------------------------------------------------

  call unpacking(xdof)

  !-!-!-!-!-!-!-!-!-!-tflux-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( weight_tflux .ge. machprec ) then

     call torflux(0)
     if ( abs(target_tflux) < machprec ) then
        target_tflux = psi_avg
        if(myid .eq. 0) write(ounit,'(8X,": Reset target toroidal flux to "ES12.5)') target_tflux
     else if (sum(abs(coil(1:Ncoils)%Ic)) == Ncoils) then !only valid when all currents are free;
           do icoil = 1, Ncoils
              coil(icoil)%I = coil(icoil)%I * target_tflux / psi_avg
           enddo
           if(myid .eq. 0) write(ounit,'(8X,": rescale coil currents with a factor of "ES12.5)') &
                target_tflux / psi_avg
     endif

     call torflux(0)
     if (abs(tflux) > machprec) weight_tflux = weight_tflux / tflux * target_tflux**2
     if( myid .eq. 0 ) write(ounit, 1000) "weight_tflux", weight_tflux
     if( myid .eq. 0 .and. weight_tflux < machprec) write(ounit, '("warning : weight_tflux < machine_precision, tflux will not be used.")')

  endif

  !-!-!-!-!-!-!-!-!-!-bharm or bnorm!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( weight_bharm >= machprec .or. weight_bnorm >= machprec ) then

     call bnormal(0)

     if ( weight_bharm >= machprec ) then 
        modBn = sqrt(sum(Bmnc**2 + Bmns**2))
        modtBn = sqrt(sum(tBmnc**2 + tBmns**2))
!!$        do icoil = 1, Ncoils
!!$           coil(icoil)%I = coil(icoil)%I * modtBn / modBn
!!$        enddo
        if(myid .eq. 0) write(ounit,'(8X,": Please rescale coil currents with a factor of "ES12.5)') &
             modtBn / modBn
        call bnormal(0)
        if (abs(bharm) > machprec) weight_bharm = weight_bharm / bharm
        if( myid == 0 ) write(ounit, 1000) "weight_bharm", weight_bharm
        if( myid .eq. 0 .and. weight_bharm < machprec) write(ounit, '("warning : weight_bharm < machine_precision, bharm will not be used.")')
     endif

     if ( weight_bnorm >= machprec ) then
        if (abs(bnorm) > machprec) weight_bnorm = weight_bnorm / bnorm
        if( myid == 0 ) write(ounit, 1000) "weight_bnorm", weight_bnorm
        if( myid .eq. 0 .and. weight_bnorm < machprec) write(ounit, '("warning : weight_bnorm < machine_precision, bnorm will not be used.")')
     endif

  endif

  !-!-!-!-!-!-!-!-!-!-total current-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( weight_isum >= machprec ) then

     call current(0) 
     if (abs(isum) > machprec) weight_isum = weight_isum / isum
     if( myid == 0 ) write(ounit, 1000) "weight_isum", weight_isum
     if( myid .eq. 0 .and. weight_isum < machprec) write(ounit, '("warning : weight_isum < machine_precision, total current constranit will not be used.")')

  endif

  !-!-!-!-!-!-!-!-!-!-ttlen-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( weight_ttlen .ge. machprec ) then

     if ( (case_length == 1) .and. (sum(coil(1:Ncoils)%Lo) < machprec) ) then
        coil(1:Ncoils)%Lo = one
        call length(0)
        coil(1:Ncoils)%Lo = coil(1:Ncoils)%L
        if(myid .eq. 0) write(ounit,'(8X,": reset target coil length to the current actual length. ")')
     endif

     call length(0)

     if (abs(ttlen) .gt. machprec) weight_ttlen = weight_ttlen / ttlen
     if( myid .eq. 0 ) write(ounit, 1000) "weight_ttlen", weight_ttlen
     if( myid .eq. 0 .and. weight_ttlen < machprec) write(ounit, '("warning : weight_ttlen < machine_precision, ttlen will not be used.")')

  endif

  !-!-!-!-!-!-!-!-!-!-curv-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( weight_curv >= machprec ) then

     call curvature(0) 
     if (abs(curv) > machprec) weight_curv = weight_curv / curv
     if( myid == 0 ) write(ounit, 1000) "weight_curv", weight_curv
     if( myid .eq. 0 .and. weight_curv < machprec) write(ounit, '("warning : weight_curv < machine_precision, curvature will not be used.")')

  endif
  !-!-!-!-!-!-!-!-!-!-curv-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( weight_straight >= machprec ) then

     call straight(0) 
     if (abs(str) > machprec) weight_straight = weight_straight / str
     if( myid == 0 ) write(ounit, 1000) "weight_straight", weight_straight
     if( myid .eq. 0 .and. weight_straight < machprec) write(ounit, '("warning : weight_straight < machine_precision, straight-out coils&
									 will not be used.")')
  endif

  !-!-!-!-!-!-!-!-!-!-cssep-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( weight_cssep >= machprec ) then

     call surfsep(0)   
     if (abs(cssep) > machprec) weight_cssep = weight_cssep / cssep
     if( myid == 0 ) write(ounit, 1000) "weight_cssep", weight_cssep
     if( myid .eq. 0 .and. weight_cssep < machprec) write(ounit, '("warning : weight_cssep < machine_precision, cssep will not be used.")')

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

  if( weight_ccsep >= machprec ) then

     call coilsep(0)
     if (abs(ccsep) > machprec) weight_ccsep = weight_ccsep / ccsep
     if( myid .eq. 0 ) write(ounit, 1000) "weight_ccsep", weight_ccsep
     if( myid .eq. 0 .and. weight_curv < machprec) write(ounit, '("warning : weight_ccsep < machine_precision, ccsep will not be used.")')
   
  endif

  !-!-!-!-!-!-!-!-!-!-tors-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( weight_tors >= machprec ) then

     call torsion(0)
     if (abs(tors) > machprec) weight_tors = weight_tors / tors
     if( myid == 0 ) write(ounit, 1000) "weight_tors", weight_tors
     if( myid .eq. 0 .and. weight_tors < machprec) write(ounit, '("warning : weight_tors < machine_precision, tors will not be used.")')

  endif

  !-!-!-!-!-!-!-!-!-!-nissin-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( weight_nissin > 0) then

     call nisscom(0)
     if (abs(nissin) > machprec) weight_nissin = weight_nissin / nissin
     if( myid == 0 ) write(ounit, 1000) "weight_nissin", weight_nissin
     if( myid .eq. 0 .and. weight_nissin < machprec) write(ounit, '("warning : weight_nissin < machine_precision, nissin will not be used.")')

  endif

  !-!-!-!-!-!-!-!-!-!-sbnorm-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( weight_sbnorm > 0.0_dp ) then

     call sbnormal(0)
     if (abs(bnormavg) > machprec) weight_sbnorm = weight_sbnorm / bnormavg
     if( myid == 0 ) write(ounit, 1000) "weight_sbnorm", weight_sbnorm
     if( myid .eq. 0 .and. weight_sbnorm < machprec) write(ounit, '("warning : weight_sbnorm < machine_precision, bnormavg will not be used.")')
  
  endif

1000 format(8X,": "A15" is normalized to " ES12.5)

  call packdof(xdof)

  return

end subroutine normweight

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine output (mark)

  use globals, only: dp, zero, ounit, myid, ierr, astat, iout, Nouts, Ncoils, save_freq, Tdof, &
       coil, coilspace, FouCoil, chi, t1E, bnorm, bnormavg, bharm, tflux, ttlen, cssep, specw, ccsep, &
       evolution, xdof, DoF, exit_tol, exit_signal, sumDE, curv, tors, nissin, MPI_COMM_FOCUS, IsQuiet, &
       weight_bnorm, weight_bharm, weight_tflux, weight_ttlen, weight_cssep, weight_curv, weight_ccsep, &
       weight_tors, weight_nissin, weight_sbnorm, lim_out, weight_straight, str, coil_type_spline, Splines, &
       isum, weight_isum 

  use mpi
  implicit none
  REAL, INTENT( IN ) :: mark
  INTEGER            :: idof, NF, icoil, NCP

  iout = iout + 1
  FATAL( output , iout > Nouts+2, maximum iteration reached )
  if (lim_out .eq. 0) then
     if (myid == 0) then
        if (IsQuiet < 0) then
           write(ounit, '("output  : "I6" : "15(ES12.5," ; "))') iout, mark, chi, sumdE, bnorm, bharm, &
                   tflux, isum, ttlen, cssep, curv, ccsep, tors, nissin, bnormavg, str
        else
           write(ounit, '("output  : "I6" : "4(ES15.8," ; "))') iout, mark, chi, sumdE, bnorm
        endif
     endif
  else
     if (myid == 0) then
        write(ounit, '("output  : "I6" : "3(ES12.5," ; "))',advance="no") iout, mark, chi, sumdE
        if ( weight_bnorm .ne. 0.0 ) then
           write(ounit, '(ES12.5," ; ")',advance="no") bnorm
        endif
        if ( weight_bharm .ne. 0.0 ) then
           write(ounit, '(ES12.5," ; ")',advance="no") bharm
        endif
        if ( weight_tflux .ne. 0.0 ) then
           write(ounit, '(ES12.5," ; ")',advance="no") tflux
        endif
        if ( weight_isum .ne. 0.0 ) then
           write(ounit, '(ES12.5," ; ")',advance="no") isum
        endif               
        if ( weight_ttlen .ne. 0.0 ) then
           write(ounit, '(ES12.5," ; ")',advance="no") ttlen
        endif 
        if ( weight_cssep .ne. 0.0 ) then
           write(ounit, '(ES12.5," ; ")',advance="no") cssep
        endif
        if ( weight_curv .ne. 0.0 ) then
           write(ounit, '(ES12.5," ; ")',advance="no") curv
        endif
        if ( weight_ccsep .ne. 0.0 ) then
           write(ounit, '(ES12.5," ; ")',advance="no") ccsep
        endif
        if ( weight_tors .ne. 0.0 ) then
           write(ounit, '(ES12.5," ; ")',advance="no") tors
        endif
        if ( weight_nissin .ne. 0.0 ) then
           write(ounit, '(ES12.5," ; ")',advance="no") nissin
        endif
        if ( weight_sbnorm .ne. 0.0 ) then
           write(ounit, '(ES12.5," ; ")',advance="no") bnormavg
        endif
        if ( weight_straight .ne. 0.0 ) then
            write(ounit, '(ES12.5," ; ")',advance="no") str
        endif
        write(*,*)
     endif
  endif

  ! save evolution data;
  if (allocated(evolution)) then
     evolution(iout,0) = mark
     evolution(iout,1) = chi
     evolution(iout,2) = sumdE
     evolution(iout,3) = bnorm
     evolution(iout,4) = bharm
     evolution(iout,5) = tflux
     evolution(iout,6) = isum
     evolution(iout,7) = ttlen
     evolution(iout,8) = cssep
     evolution(iout,9) = curv
     evolution(iout,10) = ccsep
     evolution(iout,11)= tors
     evolution(iout,12)= nissin
     evolution(iout,13)= bnormavg
     evolution(iout,14) = str
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

        select case (coil(icoil)%type)
        case (1)
           NF = FouCoil(icoil)%NF
           coilspace(iout, idof+1:idof+NF+1) = FouCoil(icoil)%xc(0:NF) ; idof = idof + NF +1
           coilspace(iout, idof+1:idof+NF  ) = FouCoil(icoil)%xs(1:NF) ; idof = idof + NF
           coilspace(iout, idof+1:idof+NF+1) = FouCoil(icoil)%yc(0:NF) ; idof = idof + NF +1
           coilspace(iout, idof+1:idof+NF  ) = FouCoil(icoil)%ys(1:NF) ; idof = idof + NF
           coilspace(iout, idof+1:idof+NF+1) = FouCoil(icoil)%zc(0:NF) ; idof = idof + NF +1
           coilspace(iout, idof+1:idof+NF  ) = FouCoil(icoil)%zs(1:NF) ; idof = idof + NF
!!$        case default
!!$           FATAL(output, .true., not supported coil types)
        case (coil_type_spline)
           NCP = Splines(icoil)%NCP
           coilspace(iout, idof+1:idof+3*NCP) = Splines(icoil)%Cpoints(0:3*NCP-1) ; idof = idof + 3*NCP
        end select
     enddo
     FATAL( output , idof .ne. Tdof, counting error in restart )
  endif

  if(mod(iout,save_freq) .eq. 0) call saving

  call MPI_BARRIER( MPI_COMM_FOCUS, ierr ) ! wait all cpus;

  return  

end subroutine output

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
