
! subroutine list
!
! initial
! check_input
! read_namelist
! write_namelist
! mute
!
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
subroutine initial
   use globals
   use mpi
   implicit none
 
   INTEGER :: index_dot
 
   !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
 
   myid = 0 ; ncpu = 1
 
   ! MPI initialize
   call MPI_init( ierr )
   MPI_COMM_FAMUS = MPI_COMM_WORLD
   call MPI_COMM_RANK( MPI_COMM_FAMUS, myid, ierr )
   call MPI_COMM_SIZE( MPI_COMM_FAMUS, ncpu, ierr )
 
   if(myid == 0) write(ounit, *) "---------------------  FAMUS ", version, "------------------------------"
   if(myid == 0) write(ounit,'("famus   : Begin execution with ncpu =",i5)') ncpu
 
   !-------------read input namelist----------------------------------------------------------------------
   if(myid == 0) then ! only the master node reads the input; 25 Mar 15;
       call getarg(1,ext) ! get argument from command line
       select case(trim(ext))
       case ( '-h', '--help' )
           write(ounit,*)'-------HELP INFORMATION--------------------------'
           write(ounit,*)' Usage: xfocus <options> input_file'
           write(ounit,*)'    <options>'
           write(ounit,*)'     --init / -i  :  Write an example input file'
           write(ounit,*)'     --help / -h  :  Output help message'
           write(ounit,*)'-------------------------------------------------'
           call MPI_ABORT( MPI_COMM_FAMUS, 0, ierr )
       case ( '-i', '--init' )
           call write_namelist ! in initial.h
       case default
           index_dot = INDEX(ext,'.input')
           IF (index_dot .gt. 0)  ext = ext(1:index_dot-1)
#ifdef DEBUG
           write(ounit, '("DEBUG info: extension from command line is "A)') trim(ext)
#endif
       end select
   endif
 
   ClBCAST( ext,  100,  0 )
   inputfile = trim(ext)//".input"
   
   call read_namelist(inputfile)
 
   return
 
end subroutine initial

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine read_namelist(filename)
   use globals, only : myid, ncpu, focusin, ounit, runit, MPI_COMM_FAMUS
   use ncsx_ports_mod, only: ncsx_ports, ncsx_ports_on
   use mpi
   implicit none
 
   CHARACTER(100), INTENT(IN) :: filename
 
   LOGICAL :: exist
   INTEGER :: icpu, ierr, ports_nml_stat
    
   if( myid == 0 ) then
    inquire(file=trim(filename), EXIST=exist) ! inquire if inputfile existed;
    FATAL( initial, .not.exist, input file ext.input not provided )
#ifdef DEBUG
    write(ounit, '("        : read namelist from ", A)') trim(filename)
#endif
   endif
 
   do icpu = 1, ncpu
      call MPI_BARRIER( MPI_COMM_FAMUS, ierr )
      if (myid == icpu-1) then                              ! each cpu read the namelist in turn;
        open(runit, file=trim(filename), status="old", action='read')
        read(runit, focusin)

        ! turn on ncsx_ports functions only if namelist is found:
        read(unit=runit, nml=ncsx_ports, iostat=ports_nml_stat)
        if (ports_nml_stat == 0) ncsx_ports_on = .true.   

        close(runit)
      endif ! end of if( myid == 0 )
   enddo
 
   return
end subroutine read_namelist

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine check_input
  use globals
  use mpi
  implicit none

  LOGICAL :: exist

  !-------------machine constants -----------------------------------------------------------------------
  machprec = epsilon(pi)         ! get the machine precision
  sqrtmachprec = sqrt(machprec)  ! sqrt of machine precision
  vsmall = ten * machprec        ! very small number
  small = thousand * machprec    ! small number

  if (myid == master) then
   write(ounit, '("initial : machine_prec   = ", ES12.5, " ; sqrtmachprec   = ", ES12.5)') & 
      machprec, sqrtmachprec
  endif 
  
  !-------------output files name ---------------------------------------------------------------------------

  hdf5file   = "focus_"//trim(ext)//".h5"
  out_famus  = trim(ext)//".focus"
  out_coils  = trim(ext)//".coils"
  out_harm   = trim(ext)//".harmonics"
  out_plasma = trim(ext)//".plasma"

  !-------------show the namelist for checking----------------------------------------------------------

  if (myid == 0) then ! Not quiet to output more informations;

     write(ounit, *) "-----------INPUT NAMELIST------------------------------------"
     write(ounit, '("initial : Read namelist focusin from : ", A)') trim(inputfile)
     write(ounit, '("        : Read plasma boundary  from : ", A)') trim(input_surf)
     if (weight_bharm > machprec) then
        write(ounit, '("        : Read Bmn harmonics    from : ", A)') trim(input_harm)
     endif

     select case( case_init )
     case(-1 )
        if (trim(input_coils) == 'none') input_coils = "coils."//trim(ext)
        inquire( file=trim(input_coils), exist=exist )
        FATAL( initial, .not.exist, coils file coils.ext not provided )
        FATAL( initial, NFcoil <= 0    , no enough harmonics )
        FATAL( initial, Nseg   <= 0    , no enough segments  )
        FATAL( initial, target_length  < zero, illegal )        
        write(ounit, '("        : Read initial coils    from : ", A, A)') trim(input_coils), '(MAKEGRID format)'
     case( 0 )
        if (trim(input_coils) == 'none') input_coils = trim(ext)//".focus"
        inquire( file=trim(input_coils), exist=exist )
        FATAL( initial, .not.exist, FOCUS coil file ext.focus not provided )
        ! FATAL( initial, mod(momentq, 2) .ne. 0, momentq can only be even number )        
        write(ounit, '("        : Read initial dipole   from : ", A, A)') trim(input_coils), '(Parameters only)'
        inquire( file=trim(fixed_coils), exist=exist )
        !FATAL( initial, .not.exist, fixed coil file ext.focus not provided )
        if(exist) write(ounit, '("        : Read fixed coils      from : ", A, A)') trim(fixed_coils), '(FOCUS format)'
     case( 1 )
        FATAL( initial, Ncoils < 1, should provide the No. of coils)
        FATAL( initial, init_current == zero, invalid coil current)
        FATAL( initial, init_radius < zero, invalid coil radius)
        FATAL( initial, NFcoil <= 0    , no enough harmonics )
        FATAL( initial, Nseg   <= 0    , no enough segments  )
        FATAL( initial, target_length  < zero, illegal )
        if (IsQuiet < 1) write(ounit, 1000) 'case_init', case_init, 'Initialize circular coils.'
     case( 2 )
        FATAL( initial, Ncoils < 1, should provide the No. of coils)
        FATAL( initial, init_current == zero, invalid coil current)
        FATAL( initial, init_radius < zero, invalid coil radius)
        FATAL( initial, target_length  < zero, illegal )
        if (IsQuiet < 1) write(ounit, 1000) 'case_init', case_init, 'Initialize magnetic dipoles.'
     case default
        FATAL( initial, .true., selected case_init is not supported )
     end select

     select case (IsQuiet)
     case (:-2)
        write(ounit, 1000) 'IsQuiet', IsQuiet, 'Output all information.'
     case (-1 )
        write(ounit, 1000) 'IsQuiet', IsQuiet, 'Output detailed information.'
     case ( 0 )
        write(ounit, 1000) 'IsQuiet', IsQuiet, 'Output essential information.'
     case ( 1:)
        write(ounit, 1000) 'IsQuiet', IsQuiet, 'Output concise information.'
     case default
        FATAL( initial, .true., IsQuiet /= integer unspported option)
     end select
     
1000 format(8X, ": ", A15, " = ", I6, " ; ", A)        

     select case (IsSymmetric)
     case (0)
        if (IsQuiet < 0) write(ounit, 1000) 'IsSymmetric', IsSymmetric, & 
             &  'No stellarator symmetry or periodicity enforced.'
     case (1)
        if (IsQuiet < 0) write(ounit, 1000) 'IsSymmetric', IsSymmetric, &
             &  'Plasma boundary periodicity is enforced.'
     case (2)
        if (IsQuiet < 0) write(ounit, 1000) 'IsSymmetric', IsSymmetric, &
             &  'Stellarator symmetry for magnetic dipoles are enforced.' 
     case default
        FATAL( initial, .true., IsSymmetric /= 0 or 2 unspported option)
     end select

     FATAL( initial, Nteta   <=   0, illegal surface resolution )
     FATAL( initial, Nzeta   <=   0, illegal surface resolution )     

     FATAL( initial, case_coils /= 1, only fourier representation is valid )
     if (IsQuiet < 0) write(ounit, 1000) 'case_coils', case_coils, 'Using Fourier series as the basic representation.'

     select case (case_optimize)
     case ( -2 )
        write(ounit, 1000) 'case_optimize', case_optimize, 'Test the 2nd derivatives.'
        FATAL( initial, .true., 2nd deriavtives are not ready.)
     case ( -1 )
        write(ounit, 1000) 'case_optimize', case_optimize, 'Test the 1st derivatives.'
     case (  0 )
        write(ounit, 1000) 'case_optimize', case_optimize, 'No optimization will be performed.'
     case (  1 )
        write(ounit, 1000) 'case_optimize', case_optimize, 'Gradient-based optimizations will be performed.'

        if (DF_maxiter > 0) then
           if (IsQuiet < 1) write(ounit, '(26X,": Differential flow will be used, maxiter=", I6)') DF_maxiter
           if (IsQuiet < 0) write(ounit, '(26X,": DF_tausta = "ES12.5" ; DF_tauend = "ES12.5 &
                &  " ; DF_xtol = "ES12.5)') DF_tausta, DF_tauend, DF_xtol
        endif

        if (CG_maxiter > 0) then        
           FATAL( Initial, CG_wolfe_c1 <= zero .or. CG_wolfe_c1 >= one, should be 0<c1<1 )
           FATAL( Initial, CG_wolfe_c2 <= zero .or. CG_wolfe_c2 >= one, should be 0<c2<1 )
           FATAL( Initial, CG_wolfe_c1 >= CG_wolfe_c2, should be c1<c2)
           if (IsQuiet < 1) write(ounit, '(26X,": Nonlinear Conjugate Gradient method will be used, maxiter="&
                &  , I6)') CG_maxiter
           if (IsQuiet < 0) write(ounit, '(26X,": CG_wolfe_c1 = "ES12.5" ; CG_wolfe_c2 = "ES12.5 &
                &  " ; CG_xtol = " ES12.5)') CG_wolfe_c1, CG_wolfe_c2, CG_xtol
        endif

     case (  2 )
        write(ounit, 1000) 'case_optimize', case_optimize, 'Hessian-based optimizations will be performed.'
        FATAL( initial, .true., Hessian-based optimizations are not ready.)

        if (DF_maxiter > 0) then
           if (IsQuiet < 1) write(ounit, '(26X,": Differential flow will be used, maxiter=", I6)') DF_maxiter
           if (IsQuiet < 0) write(ounit,'(26X,": DF_tausta = "ES12.5" ; DF_tauend = "ES12.5 &
                &  " ; DF_xtol = "ES12.5)') DF_tausta, DF_tauend, DF_xtol
        endif

        if (CG_maxiter > 0) then        
           FATAL( Initial, CG_wolfe_c1 <= zero .or. CG_wolfe_c1 >= one, should be 0<c1<1 )
           FATAL( Initial, CG_wolfe_c2 <= zero .or. CG_wolfe_c2 >= one, should be 0<c2<1 )
           FATAL( Initial, CG_wolfe_c1 >= CG_wolfe_c2, should be c1<c2)
           if (IsQuiet < 1) write(ounit, '(26X,": Nonlinear Conjugate Gradient method will be used, maxiter="&
                &  , I6)') CG_maxiter
           if (IsQuiet < 0) write(ounit,'(26X,": CG_wolfe_c1 = "ES12.5" ; CG_wolfe_c2 = "ES12.5 &
                &  " ; CG_xtol = "ES12.5)') CG_wolfe_c1, CG_wolfe_c2, CG_xtol
        endif

     case default
        FATAL( initial, .true., selected case_optimize is not supported )
     end select

!!$     if (allow_inverse .and. mod(momentq, 2)==0) then
!!$        write(ounit, '("warnning: You are using an even number for MomentQ and will not be able to have negative directions")')
!!$     endif
     
     if (case_optimize > 0) then
        FATAL( initial, DF_maxiter < 0, must be non-negative )
        FATAL( initial, CG_maxiter < 0, must be non-negative )
     endif

     select case (IsNormalize)
     case (0)
        write(ounit, 1000) 'IsNormalize', IsNormalize, 'No normalization on coil parameters.'
     case default
        write(ounit, 1000) 'IsNormalize', IsNormalize, 'Normalization on coil parameters will be performed.'
     end select

     select case (IsNormWeight)
     case (0)
        write(ounit, 1000) 'IsNormWeight', IsNormWeight, 'No normalization on weights for constraints.'
     case default
        write(ounit, 1000) 'IsNormWeight', IsNormWeight, 'Normalization on weights for constraints will be performed.'
     end select

     select case ( case_bnormal )
     case ( 0 )
        if (IsQuiet < 1) write(ounit, 1000) 'case_bnormal', case_bnormal, 'No normalization on Bn.'
     case ( 1 )
        if (IsQuiet < 1) write(ounit, 1000) 'case_bnormal', case_bnormal, 'Bn normalized to |B|.'
     case default
        FATAL( initial, .true., selected case_bnormal is not supported )
     end select

     select case ( case_length )
     case ( 1 )
        if (IsQuiet < 1) write(ounit, 1000) 'case_length', case_length, 'Quadratic format of length penalty.'
     case ( 2 )
        if (IsQuiet < 1) write(ounit, 1000) 'case_length', case_length, 'Exponential format of length penalty.'
     case default
        FATAL( initial, .true., selected case_length is not supported )
     end select

     FATAL( initial, weight_bnorm  < zero, illegal )
     FATAL( initial, weight_bharm  < zero, illegal )
     FATAL( initial, weight_tflux  < zero, illegal )
     FATAL( initial, weight_dpbin  < zero, illegal )
     FATAL( initial, weight_pmvol  < zero, illegal )
     FATAL( initial, weight_pmsum  < zero, illegal )
     FATAL( initial, weight_resbn  < zero, illegal )

     if (weight_resbn > machprec) then
        write(ounit, '("Res. Bn : resbn_m = ", I2," , resbn_n = ", I2)') resbn_m, resbn_n
        write(ounit, '(8X, ": target_resbn = ", ES12.5)') target_resbn
     endif

     select case ( case_postproc )
     case ( 0 )
        if (IsQuiet < 1) write(ounit, 1000) 'case_postproc', case_postproc, 'No extra post-processings.'
     case ( 1 )
        if (IsQuiet < 1) write(ounit, 1000) 'case_postproc', case_postproc, 'Coil evaluations will be performed.'
     case ( 2 )
        if (IsQuiet < 1) write(ounit, 1000) 'case_postproc', case_postproc, & 
             &  'Coil evaluations and writing SPEC input will be performed.'
     case ( 3 )
        if (IsQuiet < 1) write(ounit, 1000) 'case_postproc', case_postproc, & 
             &  'Coil evaluations and field-line tracing will be performed.'
     case ( 4 )
        if (IsQuiet < 1) write(ounit, 1000) 'case_postproc', case_postproc, & 
             &  'Vacuum Boozer coordinates decompostion will be performed.'
     case ( 5 )
        if (IsQuiet < 1) write(ounit, 1000) 'case_postproc', case_postproc, & 
             &  'A binary mgrid file will be saved.'
     case ( 7 )
        if (IsQuiet < 1) write(ounit, 1000) 'case_postproc', case_postproc, & 
             &  'Diagnostics, Write binary mgrid file, Trace poincare plot.'
     case default
        FATAL( initial, .true., selected case_postproc is not supported )
     end select

     FATAL( initial, save_freq <= 0, should not be negative )
     write(ounit, '("outputs : HDF5 outputs           are saved in : ", A)') trim(hdf5file)
     if (save_coils /= 0) then
        write(ounit, '("outputs : Optimizated coils      are saved in : ", A, " ; ", A)') &
             trim(out_famus), trim(out_coils)
     endif
     if (weight_bharm > machprec) then
        write(ounit, '("outputs : Realized Bn harmonics  are saved in : ", A)') trim(out_harm)
     endif
     if (update_plasma/=0) then
        write(ounit, '("outputs : Updated plasma boundary is saved in : ", A)') trim(out_plasma)
     endif

  endif

  ClBCAST( input_coils, 100, 0 )

  ! initialize iteration and total iterations;
  iout = 1 ; Nouts = 1
  if (case_optimize >0) then
     Nouts = DF_maxiter + CG_maxiter + +QN_maxiter + LM_maxiter + SA_maxiter
     Nouts = max(1, HY_maxiter) * Nouts
  endif

  !save weights before normalized
  tmpw_bnorm = weight_bnorm
  tmpw_bharm = weight_bharm
  tmpw_tflux = weight_tflux
  tmpt_tflux = target_tflux
  tmpw_ttlen = weight_ttlen
 !tmpw_specw = weight_specw
  tmpw_cssep = weight_cssep

  call MPI_BARRIER( MPI_COMM_FAMUS, ierr )

  return

end subroutine check_input

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

SUBROUTINE write_namelist
  use globals
  use ncsx_ports_mod, only: ncsx_ports
  use mgrid_focus
  use poincare_mod
  use mpi
  implicit none

  LOGICAL :: exist
  CHARACTER(LEN=100) :: example = 'example.input'

  inquire(file=trim(example), EXIST=exist) ! inquire if inputfile existed;
  FATAL( initial, exist, example input file example.input already existed )
  write(ounit, *) 'Writing an template input file in ', trim(example)
  open(wunit, file=trim(example), status='unknown', action='write')
  write(wunit, focusin)
  write(wunit, mgrid)
  write(wunit, poincare)
  write(wunit, ncsx_ports)
  close(wunit)

  call MPI_ABORT( MPI_COMM_FAMUS, 0, ierr )

  return
END SUBROUTINE write_namelist

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine mute(action)
   use globals, only: ounit
   implicit none 
 
   INTEGER,intent(in) :: action
   INTEGER, parameter :: iopen = 1, iclose = 0
   INTEGER            :: ios
   CHARACTER(LEN=*), parameter  :: tmp_out = 'tmp.famus_output' 
 
   ! open a tmp file for screen output
   if (action == iopen) then
     ounit = 37 ! something not used
     open(ounit, file=tmp_out, status="unknown", action="write", iostat=ios) ! create a scratch file?
     if (ios.ne.0) print *, "something wrong with open a tmp file in FAMUS.mute. IOSTAT=", ios 
   else
     close(ounit)
     ounit = 6 ! recover to screen output
   endif
   return
 end subroutine mute
