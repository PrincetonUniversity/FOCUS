!-------------------------------------------------------------------------------
! magnet_set_io.f90
!
! Subroutines input and output relating to the magnet set
!
! Author:       K. C. Hammond
! Contact:      khammond@pppl.gov
! Last updated: 2019-10-18
!-------------------------------------------------------------------------------
module magnet_set_io

implicit none

contains

subroutine arrays_zero()

    use magnet_set_globals, only: np, nv, nw, nl, vert_r, vert_z, &
                                  vert_theta, vert_sep, vert_tlcfs, &
                                  vessel_r, vessel_z, &
                                  avail_mag_wd, avail_mag_lg, &
                                  arrays_zeroed

    implicit none

    INTEGER :: i, j

    do i = 1, np
        do j = 1, nv
            vert_r(i,j) = 0
            vert_z(i,j) = 0
            vert_theta(i,j) = 0
            vert_sep(i,j) = 0
            vert_tlcfs(i,j) = 0
            vessel_r(i,j) = 0
            vessel_z(i,j) = 0
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
                                  vessel_loaded, ves_r00, ves_r10

    implicit none

    INTEGER :: vunit = 13, i, open_stat, read_stat, vm_in, vn_in
    REAL(8)    :: vrc_in, vzs_in, vrs_in, vzc_in
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
        if (vm_in == 1 .and. vn_in == 0) ves_r10 = vrc_in
    end do

    vessel_loaded = .true.
        
end subroutine read_vessel_coeffs

!-------------------------------------------------------------------------------
! read_lcfs_coeffs()
!
! Reads in Fourier coefficients parametrizing the plasma lcfs from 
! a file. The file must have one row containing the integer number of modes N,
! followed by N rows with 6 columns:
!     1. toroidal (n) mode number
!     2. poloidal (m) mode number -- NEGATIVE helicity is assumed
!     3. cosine coefficient, radial coordinate
!     4. sine coefficient, z coordinate
!     5. sine coefficient, r coordinate
!     6. cosine coefficient, z coordinate
! All other rows must begin with the commenting character (!)
!
! Input parameters
!     CHARACTER :: pfile  -> Path to the file with the table of coefficients.
!-------------------------------------------------------------------------------
subroutine read_lcfs_coeffs()

    use magnet_set_globals, only: pfile, nModesPl, prc, pzs, prs, pzc, pm, pn, &
                                  plasma_loaded, pla_r00

    implicit none

    INTEGER :: punit = 13, i, open_stat, read_stat, pm_in, pn_in
    REAL(8)    :: prc_in, pzs_in, prs_in, pzc_in
    logical :: nModes_read = .false.

    if (pfile == 'none') then
        write(*, *) 'Warning: no LCFS file given'
        write(*, *) 'LCFS-dependent operations will not be available.'
        return
    end if

    ! Open the file
    open(unit=punit, file=pfile, status='old', action='read', iostat=open_stat)
    if (open_stat /= 0) then
        write(*, *) 'Warning: unable to open vessel file ' // trim(pfile)
        write(*, *) 'Vessel-dependent operations will not be available.'
        return
    end if

    ! Read the number of modes
    do while (.not. nModes_read)
        read(unit=punit, fmt=*, iostat=read_stat) nModesPl
        if (read_stat > 0) cycle
        if (read_stat < 0) exit
        nModes_read = .true.
    end do
    if (.not. nModes_read) then
        write(*, *) 'Warning: unable to read number of modes from file ' &
                    // trim(pfile)
        write(*, *) 'Vessel-dependent operations will not be available.'
        return
    end if

    ! Allocate the mode arrays
    allocate( prc(nModesPl), pzs(nModesPl), prs(nModesPl), pzc(nModesPl), &
              pm(nModesPl),  pn(nModesPl) )

    ! Read in the coefficients
    i = 0
    do while (i < nModesPl)
        read(unit=punit, fmt=*, iostat=read_stat) &
            pn_in, pm_in, prc_in, prs_in, pzc_in, pzs_in
        if (read_stat > 0) cycle
        if (read_stat < 0) then
            write(*, *) 'Warning: unable to obtain expected number of ' &
                        // 'coefficients from file ' // trim(pfile)
            write(*, *) '    Number expected: ', nModesPl
            write(*, *) '    Number obtained: ', i
            write(*, *) 'Vessel-dependent operations will not be available.'
            return
        end if
        i = i + 1
        pm(i) = pm_in
        pn(i) = pn_in
        prc(i) = prc_in
        pzs(i) = pzs_in
        prs(i) = prs_in
        pzc(i) = pzc_in
        if (pm_in == 0 .and. pn_in == 0) pla_r00 = prc_in
    end do

    plasma_loaded = .true.
        
end subroutine read_lcfs_coeffs

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

    use magnet_set_globals, only: magnets, nMagnets_total, &
                                  subtraps, nSubtraps, &
                                  magnets_initialized, subtraps_initialized

    implicit none

    CHARACTER(len=200), intent(IN) :: filename
    CHARACTER(len=20),  intent(IN) :: output_type
    CHARACTER(len=13) :: mag_name
    INTEGER :: i, file_unit = 12, open_status, write_status

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

        ! Verify that the magnets have been initialized
        if (.not. magnets_initialized) then
            stop 'print_magnets_to_file: magnet set has not been initialized'
        end if

        write(unit=file_unit, fmt='(A)') '# Total number of dipoles '
        write(unit=file_unit, fmt='(I0)', iostat=write_status) nMagnets_total
        if (write_status /= 0) then
            write(*, *) 'print_magnets_to_file: unable to write ' &
                        // 'first line of ' // filename
            stop
        end if

        write(unit=file_unit, fmt='(A)') &
             '# coiltype, symmetry,  coilname,  ox,  oy,  oz,  Ic,  M_0,  ' &
             // 'pho,  Lc,  mp,  mt '
        
        do i = 1, nMagnets_total
    
            write(mag_name, fmt='(A, I10.10)') 'pm_', i
            write(unit=file_unit, &
                  fmt='(I0, A, I0, 3A, '                         // &
                        'ES15.8, A, ES15.8, A, ES15.8, A, '      // &
                        'ES15.8, A, ES15.8, A, ES15.8, A, ES15.8)', &
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
       
        ! Verify that the magnets have been initialized
        if (.not. magnets_initialized) then
            stop 'print_magnets_to_file: magnet set has not been initialized'
        end if

        do i = 1, nMagnets_total
    
            write(unit=file_unit, &
                  fmt='(ES15.8, X, ES15.8, X, ES15.8, X, ES15.8, X, ' // &
                      ' ES15.8, X, ES15.8, X, ES15.8, X, ES15.8, X, ' // &
                      ' ES15.8, X, ES15.8, X, ES15.8, X, ES15.8, X, ' // &
                      ' ES15.8, X, ES15.8 )', &
                  iostat=write_status) &
                      magnets(i)%ox, magnets(i)%oy, magnets(i)%oz, &
                      magnets(i)%nx, magnets(i)%ny, magnets(i)%nz, &
                      magnets(i)%lx, magnets(i)%ly, magnets(i)%lz, &
                      magnets(i)%lg, magnets(i)%wd, magnets(i)%ht, &
                      magnets(i)%n_phi, magnets(i)%n_theta
            if (write_status /= 0) then
                write(*, fmt='(A, I0, A)') &
                    'print_magnets_to_file: unable to write line ', i, &
                    'to ' // trim(filename)
                stop
            end if
        end do

    !---------------------------------------------------------------------------
    ! File with geometric data (trapezoid corners)
    !---------------------------------------------------------------------------
    else if (trim(output_type) == 'trapezoids') then

        ! Verify that the trapezoids have been initialized
        if (.not. subtraps_initialized) then
            stop 'print_magnets_to_file: trapezoids not initialized'
        end if

        do i = 1, nSubtraps

            write(unit=file_unit, &
                 fmt='(ES15.8, X, ES15.8, X, ES15.8, X, ES15.8, X, ' // &
                     ' ES15.8, X, ES15.8, X, ES15.8, X, ES15.8, X, ' // &
                     ' ES15.8, X, ES15.8, X, ES15.8, X, ES15.8, X, ' // &
                     ' ES15.8, X, ES15.8, X, ES15.8, X, ES15.8, X, ' // &
                     ' ES15.8, X, ES15.8, X, ES15.8, X, ES15.8, X, ' // &
                     !' ES15.8, X )', &
                     ' ES15.8, X, ES15.8, X, ES15.8, X, ES15.8, X ) ', &
                 iostat=write_status) &
                     !traps(i)%otbx, traps(i)%otby, traps(i)%otbz, &
                     !traps(i)%otfx, traps(i)%otfy, traps(i)%otfz, &
                     !traps(i)%opbx, traps(i)%opby, traps(i)%opbz, &
                     !traps(i)%opfx, traps(i)%opfy, traps(i)%opfz, &
                     !traps(i)%vx,   traps(i)%vy,   traps(i)%vz,   &
                     !traps(i)%otx,  traps(i)%oty,  traps(i)%otz,  &
                     !traps(i)%obx,  traps(i)%oby,  traps(i)%obz
                     subtraps(i)%ctx1, subtraps(i)%cty1, subtraps(i)%ctz1, &
                     subtraps(i)%ctx2, subtraps(i)%cty2, subtraps(i)%ctz2, &
                     subtraps(i)%ctx3, subtraps(i)%cty3, subtraps(i)%ctz3, &
                     subtraps(i)%ctx4, subtraps(i)%cty4, subtraps(i)%ctz4, &
                     subtraps(i)%cbx1, subtraps(i)%cby1, subtraps(i)%cbz1, &
                     subtraps(i)%cbx2, subtraps(i)%cby2, subtraps(i)%cbz2, &
                     subtraps(i)%cbx3, subtraps(i)%cby3, subtraps(i)%cbz3, &
                     subtraps(i)%cbx4, subtraps(i)%cby4, subtraps(i)%cbz4 
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

