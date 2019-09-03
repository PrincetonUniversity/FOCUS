!-------------------------------------------------------------------------------
! magnet_set_build.f90
!
! Routines for defining the locations, orientations, etc. of the magnets based
! on input geometric data
!-------------------------------------------------------------------------------
module magnet_set_build

implicit none

contains

!-------------------------------------------------------------------------------
! vertex_theta_to_rz
!
! Calculates the r and z coordinates of plate vertices based on poloidal angles
! and separation distances from the vessel.
!-------------------------------------------------------------------------------
subroutine vertex_theta_to_rz()

    use magnet_set_globals, only: pi, nPlates, nVertices, plate_phi, &
                                  vert_theta, vert_sep, vert_r, vert_z, &
                                  vessel_r, vessel_z
    use magnet_set_calc,    only: ves_r, ves_z, ves_drdt, ves_dzdt, cross_prod

    implicit none

    INTEGER :: i, j
    REAL :: drdt, dzdt, alpha, nr, nz

    do i = 1, nPlates

        do j = 1, nVertices

            ! Assume vertices with (0,0) are meant to not exist
            if (vert_theta(i,j) == 0 .and. vert_sep(i,j) == 0) then

                vert_r(i,j) = 0
                vert_z(i,j) = 0
                vessel_r(i,j) = 0
                vessel_z(i,j) = 0

            else

                ! Coordinates and derivatives on the vessel
                vessel_r(i,j) = ves_r( vert_theta(i,j), plate_phi(i) )
                vessel_z(i,j) = ves_z( vert_theta(i,j), plate_phi(i) )
                drdt = ves_drdt( vert_theta(i,j), plate_phi(i) )
                dzdt = ves_dzdt( vert_theta(i,j), plate_phi(i) )

                ! Normal vector pointing away from the vessel
                alpha = atan2( dzdt, drdt )  - 0.5*pi
                nr = cos(alpha)
                nz = sin(alpha)

                ! r and z coordinates of the vertex
                vert_r(i,j) = vessel_r(i,j) + vert_sep(i,j)*nr
                vert_z(i,j) = vessel_z(i,j) + vert_sep(i,j)*nz

            end if

        end do

    end do

end subroutine vertex_theta_to_rz

!-------------------------------------------------------------------------------
! vertex_rz_to_theta_sep
!
! Calculates the length (separation distance) and poloidal angle of the line
! that (1) connects the vessel to each plate vertex given by r and z coordinates
! and (2) is perpendicular to the vessel in the phi plane of the plate vertex.
!-------------------------------------------------------------------------------
subroutine vertex_rz_to_theta_sep()

    use magnet_set_globals, only: pi, ves_r00, nPlates, nVertices, plate_phi, &
                                  vert_theta, vert_sep, vert_r, vert_z, &
                                  vessel_r, vessel_z
    use magnet_set_calc,    only: ves_perp_intersect_2d

    implicit none

    INTEGER :: i, j
    REAL :: l0, theta0, chi2

    do i = 1, nPlates

        do j = 1, nVertices

            ! Assume vertices with (0,0) are meant to not exist
            if (vert_r(i,j) == 0 .and. vert_z(i,j) == 0) then

                vert_sep(i,j) = 0
                vert_theta(i,j) = 0
                vessel_r(i,j) = 0
                vessel_z(i,j) = 0

            else

                l0 = 0.0
                theta0 = atan2(vert_z(i,j), vert_r(i,j)-ves_r00)
                call ves_perp_intersect_2d(plate_phi(i), &
                         vert_r(i,j), vert_z(i,j), l0, theta0, &
                         vert_sep(i,j), vert_theta(i,j), &
                         vessel_r(i,j), vessel_z(i,j), chi2)

            end if

        end do

    end do

end subroutine vertex_rz_to_theta_sep

!-------------------------------------------------------------------------------
! count_stacks()
! 
! Determines the number of available places for magnets (or stacks of magnets)
! based the number of plates,
! the length of each plate, and the length of each magnet/stack
!-------------------------------------------------------------------------------
subroutine count_stacks()
 
    use magnet_set_globals, only: pi, vertex_mode, nPlates, nVertices, &
                                  plate_phi, plate_dphi, &
                                  vert_r, vert_z, vert_theta, vert_sep, &
                                  vessel_r, vessel_z, &
                                  segs_on_plate, nSegments, segments, &
                                  gap_tor, gap_pol, gap_end, &
                                  stack_ht_max, avail_mag_lg, &
                                  nStacks_total, &
                                  plates_initialized, segments_initialized

    implicit none

    INTEGER :: i, j, n
    REAL :: vert_dr, vert_dz
    REAL, allocatable :: seg_lg(:,:), alpha(:,:), beta(:,:), avail_dims(:)

    ! Initial status checks
    if (.not. plates_initialized) then
        stop 'count_magnets: plates not yet initialized'
    end if
    if (segments_initialized) then
        stop 'count_magnets: plate segments have already been initialized'
    end if

    ! Calculate vertex r and z coords if not explicitly supplied
    if (trim(vertex_mode) == 'theta_sep') then
        call vertex_theta_to_rz
    else if (trim(vertex_mode) == 'rz') then
        call vertex_rz_to_theta_sep
    else
        write(*,*) 'count_stacks: unrecognized vertex_mode'
        stop
    end if

    ! Determine the lengths and angles of each vertex
    nSegments = 0
    allocate( seg_lg(nPlates, nVertices-1),  alpha(nPlates, nVertices-1) , &
              beta(nPlates, nVertices-1) )
    do i = 1, nPlates
        segs_on_plate(i) = 0
        do j = 1, nVertices-1

            ! If no more nonzero vertices, move on to the next plate
            if (vert_r(i,j+1) == 0 .and. vert_z(i,j+1) == 0) exit

            ! Segment length and angle relative to r-hat
            vert_dr = vert_r(i,j+1) - vert_r(i,j)
            vert_dz = vert_z(i,j+1) - vert_z(i,j)
            seg_lg(i,j) = sqrt( vert_dr*vert_dr + vert_dz*vert_dz )
            alpha(i,j) = atan2( vert_dz, vert_dr ) - 0.5 * pi

            ! Angle of the line to the segment from the vessel reference point
            beta(i,j) = atan2(vert_z(i,j) - vessel_z(i,j), &
                              vert_r(i,j) - vessel_r(i,j))

            ! Running total of the segments on the current support plate
            nSegments = nSegments + 1
            segs_on_plate(i) = segs_on_plate(i) + 1

        end do
    end do

    ! Determine the gaps required by the concavity at the ends of each segment
    n = 0
    allocate(segments(nSegments))
    do i = 1, nPlates
        do j = 1, segs_on_plate(i)

            n = n + 1
            segments(n)%lg = seg_lg(i,j)
            segments(n)%alpha = alpha(i,j)
            segments(n)%r1 = vert_r(i,j)
            segments(n)%r2 = vert_r(i,j+1)
            segments(n)%z1 = vert_z(i,j)
            segments(n)%z2 = vert_z(i,j+1)
            segments(n)%phi = plate_phi(i)
            segments(n)%dphi = plate_dphi(i)

            if (j == 1) then
                segments(n)%gap_con_1 = 0
            else
                segments(n)%gap_con_1 = &
                    concavity_gap(alpha(i,j-1), alpha(i,j), stack_ht_max, &
                                  1, vert_sep(i,j), beta(i,j))
            end if

            if (j == segs_on_plate(i)) then
                segments(n)%gap_con_2 = 0
            else
                segments(n)%gap_con_2 = &
                    concavity_gap(alpha(i,j), alpha(i,j+1), stack_ht_max, &
                                  2, vert_sep(i,j+1), beta(i,j+1))
           end if

        end do
    end do

    ! Determine the magnet length and number of magnets for each segment
    nStacks_total = 0
    allocate(avail_dims(size(avail_mag_lg)))
    do i = 1, size(avail_mag_lg)
        avail_dims(i) = avail_mag_lg(i)
    end do
    do i = 1, nSegments
        call best_magnet_dim(avail_dims, gap_pol, &
                             gap_end + segments(i)%gap_con_1, &
                             gap_end + segments(i)%gap_con_2, &
                             segments(i)%lg, &
                             segments(i)%mag_lg_ind, segments(i)%nStacks)
        nStacks_total = nStacks_total + segments(i)%nStacks
    end do   

    deallocate(seg_lg, alpha, beta, avail_dims)

    segments_initialized = .true.

end subroutine count_stacks

!-------------------------------------------------------------------------------
! count_magnets()
! 
! Determines the defining properties of each individual magnet in the 
! configuration
!-------------------------------------------------------------------------------
subroutine count_magnets()

    use magnet_set_globals, only: pi, ves_r00, ves_tol, &
                                  nPlates, nVertices, plate_phi, plate_dphi, &
                                  nSegments, segments, &
                                  gap_rad, gap_tor, gap_pol, gap_end, &
                                  gap_seg, gap_vac, &
                                  stacks, nStacks_total, stack_ht_max, &
                                  fixed_stack_ht, &
                                  avail_mag_lg, avail_mag_ht, nMagnets_total, &
                                  segments_initialized, stacks_initialized
    use magnet_set_calc,    only: ves_dist_2d

    implicit none

    INTEGER :: i, j, k, n
    REAL :: rmin, width, length, seg_dist, height, stack_dist
    REAL :: row_lg, avail_lg, offset
    REAL, allocatable :: avail_dims(:)
    REAL :: x1, y1, z1                     ! coordinates of segment start point
    REAL :: x2, y2, z2                     ! coordinates of segment end point
    REAL :: ux, uy, uz                     ! unit vector along segment
    REAL :: nx, ny, nz, nr                 ! normal vector to segment
    REAL :: r, z, theta0, vdist, theta, vr, vz, chi2

    ! Initialization checks
    if (.not. segments_initialized) then
        stop 'count_magnets: plate segments not initialized'
    end if
    if (stacks_initialized) then
        stop 'count_magnets: stacks already initialized'
    end if

    ! Array of available heights
    allocate(avail_dims(size(avail_mag_ht)))
    do i = 1, size(avail_mag_lg)
        avail_dims(i) = avail_mag_ht(i)
    end do

    ! Allocate the array of magnet stacks
    allocate(stacks(nStacks_total))

    n = 0
    nMagnets_total = 0
    do i = 1, nSegments

        rmin = min(segments(i)%r1, segments(i)%r2)
        width = rmin * segments(i)%dphi - gap_tor
        length = avail_mag_lg( segments(i)%mag_lg_ind )

        row_lg = (length + gap_pol)*segments(i)%nStacks - gap_pol + 2*gap_end
        avail_lg = segments(i)%lg - &
                       (segments(i)%gap_con_1 + segments(i)%gap_con_2)
        offset = 0.5*(avail_lg - row_lg)

        x1 = segments(i)%r1 * cos( segments(i)%phi )
        y1 = segments(i)%r1 * sin( segments(i)%phi )
        z1 = segments(i)%z1

        x2 = segments(i)%r2 * cos( segments(i)%phi )
        y2 = segments(i)%r2 * sin( segments(i)%phi )
        z2 = segments(i)%z2

        ux = (x2 - x1) / segments(i)%lg
        uy = (y2 - y1) / segments(i)%lg
        uz = (z2 - z1) / segments(i)%lg

        nx = sin( 0.5*pi - segments(i)%alpha ) * cos( segments(i)%phi )  
        ny = sin( 0.5*pi - segments(i)%alpha ) * sin( segments(i)%phi )
        nz = cos( 0.5*pi - segments(i)%alpha )
        nr = nx * cos(segments(i)%phi) + ny * sin(segments(i)%phi)

        do j = 1, segments(i)%nStacks
 
            n = n + 1

            stacks(n)%nx = nx
            stacks(n)%ny = ny
            stacks(n)%nz = nz

            stacks(n)%lx = ux
            stacks(n)%ly = uy
            stacks(n)%lz = uz

            seg_dist = segments(i)%gap_con_1 + offset + gap_end + 0.5*length &
                           + (j - 1)*(length + gap_pol)

            stacks(n)%ox = x1 + seg_dist * ux
            stacks(n)%oy = y1 + seg_dist * uy
            stacks(n)%oz = z1 + seg_dist * uz

            stacks(n)%wd = width
            stacks(n)%lg = length

            ! Available height within the stack in which to fit magnets
            if (fixed_stack_ht) then

                stacks(n)%ht = stack_ht_max

            else

                ! If height is not fixed calculate clearance to vessel for stack
                r = sqrt(stacks(n)%ox**2 + stacks(n)%oy**2)
                z = stacks(n)%oz
                nz = stacks(n)%nz
                nr = stacks(n)%nx * cos(segments(i)%phi) &
                     + stacks(n)%ny * sin(segments(i)%phi)
                theta0 = atan2(z, r-ves_r00)
                call ves_dist_2d(segments(i)%phi, r, z, nr, nz, stack_ht_max, &
                                 theta0, vdist, theta, vr, vz, chi2)
                if (vdist > 0) then
                    stacks(n)%ht = 0
                elseif (vdist < -stack_ht_max .or. chi2 > ves_tol**2) then
                    stacks(n)%ht = stack_ht_max
                else
                    stacks(n)%ht = -vdist
                endif
            end if

            ! Number of magnets per stack based on available ht.
            call best_magnet_dim(avail_dims, gap_rad, gap_seg, gap_vac, &
                                 stacks(n)%ht, stacks(n)%mag_ht_ind,  &
                                 stacks(n)%nMagnets)
            nMagnets_total = nMagnets_total + stacks(n)%nMagnets

        end do

    end do

    deallocate(avail_dims)

    stacks_initialized = .true.

end subroutine count_magnets 

!-------------------------------------------------------------------------------
! magnet_properties
!
! Populates the array of individual magnets with properties relating to their
! locations and dimensions
!-------------------------------------------------------------------------------
subroutine magnet_properties

    use magnet_set_globals, only: m_bulk, symm, coiltype, &
                                  stacks, nStacks_total, &
                                  magnets, nMagnets_total, avail_mag_ht, &
                                  gap_seg, gap_rad, gap_vac, &
                                  stacks_initialized, magnets_initialized

    implicit none

    INTEGER :: i, j, n
    REAL :: height, stack_dist

    if (.not. stacks_initialized) then
        stop 'count_magnets: plate segments not initialized'
    end if
    if (magnets_initialized) then
        stop 'count_magnets: magnets already initialized'
    end if

    allocate( magnets(nMagnets_total) )

    n = 0
    do i = 1, nStacks_total

        height = avail_mag_ht( stacks(i)%mag_ht_ind )

        if (stacks(i)%nMagnets == 0) cycle

        do j = 1, stacks(i)%nMagnets

            n = n + 1

            if (n > nMagnets_total) stop '    MAGNET OVERLOAD!!!'
    
            stack_dist = gap_seg + 0.5*height + (j - 1)*(height + gap_rad)
    
            !height = (stacks(i)%ht - gap_seg &
            !                 - gap_vac + gap_rad)/stacks(i)%nMagnets - gap_rad
    
            magnets(n)%ox = stacks(i)%ox - stack_dist * stacks(i)%nx 
            magnets(n)%oy = stacks(i)%oy - stack_dist * stacks(i)%ny
            magnets(n)%oz = stacks(i)%oz - stack_dist * stacks(i)%nz
    
            magnets(n)%nx = stacks(i)%nx
            magnets(n)%ny = stacks(i)%ny
            magnets(n)%nz = stacks(i)%nz
    
            magnets(n)%lx = stacks(i)%lx
            magnets(n)%ly = stacks(i)%ly
            magnets(n)%lz = stacks(i)%lz
    
            magnets(n)%lg  = stacks(i)%lg
            magnets(n)%wd  = stacks(i)%wd
            magnets(n)%ht  = height
            magnets(n)%vol = magnets(n)%lg * magnets(n)%wd * magnets(n)%ht
    
            magnets(n)%moment = m_bulk * magnets(n)%vol
            magnets(n)%coiltype = coiltype
            magnets(n)%symm = symm

            magnets(n)%n_phi = atan2(magnets(n)%ny, magnets(n)%nx)
            magnets(n)%n_theta = &
                atan2(sqrt(magnets(n)%nx**2 + magnets(n)%ny**2), magnets(n)%nz)

        end do

    end do

    magnets_initialized = .true.

end subroutine magnet_properties


!-------------------------------------------------------------------------------
! concavity_gap(alpha1, alpha2, h, vert_sep, beta)
!
! Determines the gap spacing required between a rectangular magnet block and the
! vertex of a plate segment that makes a concave angle with an adjoining
! plate segment. If the the angle between the segments is straight or convex,
! no gap is required ( = 0 ).
!
! Note: normal vectors to plate segments are assumed to point away from the
! vessel, whereas magnets are assumed to be mounted between the plate segments
! and the vessel.
!
! Input parameters:
!     REAL :: alpha1, alpha2  -> angles (in radians) that the normal vectors
!                                of the two segments make relative to the
!                                major-radius unit vector.
!     REAL :: h               -> height of the magnet block above the plate
!                                segment
!     INTEGER :: which_alpha  -> 1 or 2, depending on which value of alpha
!                                corresponds to the segment under consideration
!     REAL :: vert_sep        -> separation distance between the vertex and
!                                the point on the vessel along a vessel normal
!     REAL :: beta            -> angle of the normal relative to the radial
!                                unit vector
!-------------------------------------------------------------------------------
REAL function concavity_gap(alpha1, alpha2, h, which_alpha, vert_sep, beta)

    use magnet_set_globals, only: pi

    implicit none

    ! Input parameters are as defined in the documentation string
    REAL, intent(IN) :: alpha1, alpha2, h, vert_sep, beta
    INTEGER, intent(IN) :: which_alpha
    REAL :: d_alpha, angle, allowable_ht

    ! Angle btw. vessel normal and perp. bisector of segment from pt. on vessel
    if (which_alpha == 1) then
        angle = beta - alpha1
    else if (which_alpha == 2) then
        angle = alpha2 - beta
    else
        write(*,*) 'concavity_gap: bad input for which_alpha parameter'
        stop
    end if
    
    ! Determine the allowable height from the segment to avoid vessel collision
    allowable_ht = vert_sep * cos(angle)

    if (allowable_ht < h) then

        concavity_gap = max(0.0, vert_sep * sin(angle))

    else

        d_alpha = alpha2 - alpha1
    
        ! Adjust d_alpha to account for phase jumps
        if (d_alpha < -pi) d_alpha = d_alpha + 2 * pi
        if (d_alpha >  pi) d_alpha = d_alpha - 2 * pi
    
        if (d_alpha <= 0) then
            ! Case of a convex or straight angle
            concavity_gap = 0.0
        else
            concavity_gap = h / tan( 0.5 * (pi - d_alpha) )
        end if

    end if

end function concavity_gap


!-------------------------------------------------------------------------------
! best_magnet_length
!
! Given a choice of available magnet lengths, selects the option that results
! in the optimal total magnet length on a particular support plate segment.
!
! Input parameters:
!     REAL :: avail_dims -> array of possible dimensions
!     REAL :: gap_btwn    -> required gap between magnets 
!     REAL :: gap_end1  -> required gap, begin. of plate segment
!     REAL :: gap_end2  -> required gap, end of plate segment
!     REAL :: l_seg      -> length of the plate segment
!
! Output parameters:
!     INTEGER :: ind_best -> index of the optimal length in avail_mag_lg
!     INTEGER :: nMagnets -> number of magnets that can fit on the segment
!-------------------------------------------------------------------------------
subroutine best_magnet_dim(avail_dims, gap_btwn, gap_end1, gap_end2, &
                           l_seg, ind_best, nMagnets)

    use magnet_set_globals, only: avail_mag_lg

    implicit none

    REAL, allocatable, intent(IN) :: avail_dims(:)
    REAL, intent(IN) :: gap_btwn, gap_end1, gap_end2, l_seg
    REAL :: reserve_per_seg, reserve_per_mag, total_length, opt_total_length 
    INTEGER :: i, nAvail, nMagnets_i
    INTEGER, intent(OUT) :: ind_best, nMagnets

    nAvail = size(avail_dims)
    opt_total_length = 0
    ind_best = 1
    nMagnets = 0

    do i = 1, nAvail

        ! Determine the number of magnets of the given length that will fit
        reserve_per_mag = avail_dims(i) + gap_btwn
        reserve_per_seg = gap_end1 + gap_end2 - gap_btwn
        nMagnets_i = floor((l_seg - reserve_per_seg)/reserve_per_mag)
        if (nMagnets_i < 0) cycle
        total_length = avail_mag_lg(i) * nMagnets_i

        ! Update the optimum magnet size & number if applicable
        if (total_length >= opt_total_length) then
            opt_total_length = total_length
            ind_best = i
            nMagnets = nMagnets_i
        end if

    end do

end subroutine best_magnet_dim

end module magnet_set_build
