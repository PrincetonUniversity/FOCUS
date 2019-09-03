!-------------------------------------------------------------------------------
! magnet_set_io.f90
!
! Subroutines input and output relating to the magnet set
!-------------------------------------------------------------------------------
module magnet_set_io

implicit none

contains

subroutine arrays_zero()

    use magnet_set_globals, only: np, nv, nw, nl, vert_r, vert_z, &
                                  avail_mag_wd, avail_mag_lg, &
                                  arrays_zeroed

    implicit none

    INTEGER :: i, j

    do i = 1, np
        do j = 1, nv
            vert_r(i,j) = 0
            vert_z(i,j) = 0
        end do
    end do

    do i = 1, nw
        avail_mag_wd(i) = 0
    end do

    do i = 1, nl
        avail_mag_lg(i) = 0
    end do

    arrays_zeroed = .true.

end subroutine arrays_zero

!-------------------------------------------------------------------------------
! read_namelist(nl)
! 
! Reads a namelist file specified by the string (nl)
!
! TODO: perform some checks on the values after the namelist has been read
!-------------------------------------------------------------------------------
subroutine read_namelist(nl)

    use magnet_set_globals, only: magnet_set, plates_initialized

    implicit none

    CHARACTER(len=200), intent(IN) :: nl
    INTEGER :: nl_unit = 12, open_stat, nl_stat
    INTEGER :: i, j, imax, jmax
    CHARACTER(len=200) :: line

    ! zero out the arrays
    call arrays_zero()

    ! Read the namelist
    open(unit=nl_unit, file=nl, status='old', action='read', iostat=open_stat)
    if (open_stat /= 0) stop 'Unable to open namelist file'
    read(unit=nl_unit, nml=magnet_set, iostat=nl_stat)
    if (nl_stat /= 0) then
        backspace(nl_unit)
        read(nl_unit, fmt='(A)') line
        write(*, *) 'Invalid line in namelist : ' // trim(line)
        stop
    end if
    close(nl_unit)

    plates_initialized = .true.

end subroutine read_namelist


!-------------------------------------------------------------------------------
! read_vessel_coeffs()
!
! Reads in Fourier coefficients parametrizing the vacuum vessel surface from
! a file. The file must have one row containing the integer number of modes N,
! followed by N rows with 6 columns:
!     1. poloidal (m) mode number
!     2. toroidal (n) mode number
!     3. cosine coefficient, radial coordinate
!     4. sine coefficient, z coordinate
!     5. sine coefficient, r coordinate
!     6. cosine coefficient, z coordinate
! All other rows must begin with the commenting character (!)
!
! Input parameters
!     CHARACTER :: vfile  -> Path to the file with the table of coefficients.
!-------------------------------------------------------------------------------
subroutine read_vessel_coeffs()

    use magnet_set_globals, only: vfile, nModes, vrc, vzs, vrs, vzc, vm, vn, &
                                  vessel_loaded, ves_r00

    implicit none

    INTEGER :: vunit = 13, i, open_stat, read_stat, vm_in, vn_in
    REAL    :: vrc_in, vzs_in, vrs_in, vzc_in
    logical :: nModes_read = .false.

    if (vfile == 'none') then
        write(*, *) 'Warning: no vessel file given'
        write(*, *) 'Vessel-dependent operations will not be available.'
        return
    end if

    ! Open the file
    open(unit=vunit, file=vfile, status='old', action='read', iostat=open_stat)
    if (open_stat /= 0) then
        write(*, *) 'Warning: unable to open vessel file ' // trim(vfile)
        write(*, *) 'Vessel-dependent operations will not be available.'
        return
    end if

    ! Read the number of modes
    do while (.not. nModes_read)
        read(unit=vunit, fmt=*, iostat=read_stat) nModes
        if (read_stat > 0) cycle
        if (read_stat < 0) exit
        nModes_read = .true.
    end do
    if (.not. nModes_read) then
        write(*, *) 'Warning: unable to read number of modes from file ' &
                    // trim(vfile)
        write(*, *) 'Vessel-dependent operations will not be available.'
        return
    end if

    ! Allocate the mode arrays
    allocate( vrc(nModes), vzs(nModes), vrs(nModes), vzc(nModes), &
              vm(nModes),  vn(nModes) )

    ! Read in the coefficients
    i = 0
    do while (i < nModes)
        read(unit=vunit, fmt=*, iostat=read_stat) &
            vm_in, vn_in, vrc_in, vzs_in, vrs_in, vzc_in
        if (read_stat > 0) cycle
        if (read_stat < 0) then
            write(*, *) 'Warning: unable to obtain expected number of ' &
                        // 'coefficients from file ' // trim(vfile)
            write(*, *) '    Number expected: ', nModes
            write(*, *) '    Number obtained: ', i
            write(*, *) 'Vessel-dependent operations will not be available.'
            return
        end if
        i = i + 1
        vm(i) = vm_in
        vn(i) = vn_in
        vrc(i) = vrc_in
        vzs(i) = vzs_in
        vrs(i) = vrs_in
        vzc(i) = vzc_in
        if (vm_in == 0 .and. vn_in == 0) ves_r00 = vrc_in
    end do

    vessel_loaded = .true.
        
end subroutine read_vessel_coeffs

!-------------------------------------------------------------------------------
! print_magnets_to_file(filename, output_type)
! 
! Prints magnet parameters to file for use by other programs.
!
! Input parameters:
!     CHARACTER :: filename    -> name of the file to be saved
!     CHARACTER :: output_type -> output file specifier:
!                                     'focus'     : for input to FOCUS
!                                     'geometric' : location, orientation, size
!-------------------------------------------------------------------------------
subroutine print_magnets_to_file(filename, output_type)

    use magnet_set_globals, only: magnets, nMagnets_total, magnets_initialized

    implicit none

    CHARACTER(len=200), intent(IN) :: filename
    CHARACTER(len=20),  intent(IN) :: output_type
    CHARACTER(len=13) :: mag_name
    INTEGER :: i, file_unit = 12, open_status, write_status

    ! Verify that the magnets have been initialized
    if (.not. magnets_initialized) then
        stop 'print_magnets_to_file: magnet set has not been initialized'
    end if

    ! Open the file
    open(unit=file_unit, file=trim(filename), action='write', &
         iostat=open_status)
    if (open_status /= 0) then
        write(*, *) 'print_magnets_to_file: unable to open file ' // filename
        stop
    end if

    !---------------------------------------------------------------------------
    ! File for input to FOCUS
    !---------------------------------------------------------------------------
    if (trim(output_type) == 'focus') then

        write(unit=file_unit, fmt=*) '# Total number of dipoles '
        write(unit=file_unit, fmt='(I0)', iostat=write_status) nMagnets_total
        if (write_status /= 0) then
            write(*, *) 'print_magnets_to_file: unable to write ' &
                        // 'first line of ' // filename
            stop
        end if

        write(unit=file_unit, fmt=*) &
             '# coiltype, symmetry,  coilname,  ox,  oy,  oz,  Ic,  M_0,  ' &
             // 'pho,  Lc,  mp,  mt '
        
        do i = 1, nMagnets_total
    
            write(mag_name, fmt='(A, I10.10)') 'pm_', i
            write(unit=file_unit, &
                  fmt='(I0, A, I0, 3A, E15.8, A, E15.8, A, E15.8, A, ' // &
                                      'E15.8, A, E15.8, A, E15.8, A, E15.8)', &
                  iostat=write_status) &
                      magnets(i)%coiltype,       ', ', &
                      magnets(i)%symm,           ', ', &
                      trim(mag_name),            ', ', &
                      magnets(i)%ox,             ', ', &
                      magnets(i)%oy,             ', ', &
                      magnets(i)%oz,          ', 1, ', &
                      magnets(i)%moment,         ', ', &
                      1.0,                    ', 0, ', &
                      magnets(i)%n_phi,          ', ', &
                      magnets(i)%n_theta
            if (write_status /= 0) then
                write(*, fmt='(A, I0, A)') &
                    'print_magnets_to_file: unable to write line ', i, &
                    'to ' // trim(filename)
                stop
            end if
        end do

    !---------------------------------------------------------------------------
    ! File with geometric data
    !---------------------------------------------------------------------------
    else if (trim(output_type) == 'geometric') then
       
        do i = 1, nMagnets_total
    
            write(mag_name, fmt='(A, I10.10)') 'pm_', i
            write(unit=file_unit, &
                  fmt='(E15.8, X, E15.8, X, E15.8, X, E15.8, X, ' // &
                      ' E15.8, X, E15.8, X, E15.8, X, E15.8, X, ' // &
                      ' E15.8, X, E15.8, X, E15.8, X, E15.8      )', &
                  iostat=write_status) &
                      magnets(i)%ox, magnets(i)%oy, magnets(i)%oz, &
                      magnets(i)%nx, magnets(i)%ny, magnets(i)%nz, &
                      magnets(i)%lx, magnets(i)%ly, magnets(i)%lz, &
                      magnets(i)%lg, magnets(i)%wd, magnets(i)%ht
            if (write_status /= 0) then
                write(*, fmt='(A, I0, A)') &
                    'print_magnets_to_file: unable to write line ', i, &
                    'to ' // trim(filename)
                stop
            end if
        end do

    else

        write(*, *) 'print_magnets_to_file: unrecognized output_type'
        stop

    end if

    close(file_unit)

end subroutine print_magnets_to_file
 

end module magnet_set_io

