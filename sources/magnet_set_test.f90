!-------------------------------------------------------------------------------
! magnet_set_test.f90
!
! Test client for the magnet_set modules
!-------------------------------------------------------------------------------
program magnet_set_test

    use magnet_set_globals, only: magnet_set, plates_initialized, &
                                  nMagnets_total, pol_err_count, &
                                  construction
    use magnet_set_io,      only: read_namelist, read_vessel_coeffs, &
                                  read_lcfs_coeffs, print_magnets_to_file
    use magnet_set_build,   only: count_stacks, count_magnets, magnet_properties
    use magnet_set_trapz,   only: trapezoid_properties, subdivide_trapezoids

    implicit none

    INTEGER :: nl_unit
    INTEGER :: cmd_arg_stat, open_stat, nl_stat
    CHARACTER(len=200) :: nl, out_fname, out_type, usage_str, line

    ! Obtain the command arguments
    usage_str = 'Usage: ./magnet_set_test [namelist] [output_filename] ' &
                // '[output_type]'   

    call get_command_argument(1, nl, status=cmd_arg_stat)
    if (cmd_arg_stat /= 0) stop trim(usage_str)

    call get_command_argument(2, out_fname, status=cmd_arg_stat)
    if (cmd_arg_stat /= 0) stop trim(usage_str)

    call get_command_argument(3, out_type, status=cmd_arg_stat)
    if (cmd_arg_stat /= 0) stop trim(usage_str)

    ! Read the namelist
    call read_namelist(nl)
    write(*,*) 'Namelist read'

    call read_vessel_coeffs()
    write(*,*) 'Vessel coefficients read'

    call read_lcfs_coeffs()
    write(*,*) 'LCFS coefficients read'

    if (trim(construction) == 'planar') then

        call count_stacks()
        write(*,*) 'Stacks counted'
    
        call count_magnets()
        write(*,*) 'Magnets counted'
    
        call magnet_properties()
        write(*,*) 'Magnet properties obtained'
        write(*,'(I0,a)') pol_err_count, ' magnets with erroneous polarizations'

    else if (trim(construction) == 'trapezoidal') then

        call trapezoid_properties()
        write(*,*) 'Trapezoids defined'

        call subdivide_trapezoids()
        write(*,*) 'Trapezoids subdivided and magnet properties determined'

    end if

    call print_magnets_to_file(out_fname, out_type)
    write(*,'(I0,A)') nMagnets_total, ' magnets printed to file'

end program magnet_set_test

