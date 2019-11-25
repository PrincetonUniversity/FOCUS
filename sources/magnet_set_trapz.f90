!-------------------------------------------------------------------------------
! magnet_set_trapz.f90
! 
! Subroutines for generating a magnet layout in which magnets are oriented to
! face the vessel and grouped in trapezoidal volumes whose bases (the side
! far from the vessel) align to fixed toroidal angles.
!
! Author: K. C. Hammond
! Contact: khammond@pppl.gov
! Last updated 2019-10-29
!-------------------------------------------------------------------------------
module magnet_set_trapz

implicit none

contains

!-------------------------------------------------------------------------------
! trapezoid_properties()
!
! Determines the bounding planes, stack-normal vectors, and associated 
! parameters according to the base vertex specifications
!-------------------------------------------------------------------------------
subroutine trapezoid_properties()

    use magnet_set_globals, only: pi, vertex_mode, nPlates, nVertices, &
                                  plate_phi, plate_dphi, plates_initialized, &
                                  vert_r, vert_z, vert_theta, vert_sep, &
                                  vessel_r, vessel_z, &
                                  traps, nTrapezoids, trapezoids_initialized, &
                                  trap_err_count, gap_vac, gap_trp
    use magnet_set_build, only: vertex_theta_to_rz, lcfs_theta_to_ves_theta, &
                                vertex_rz_to_theta_sep
    use magnet_set_calc, only: ves_perp_intersect_3d, ves_unorm, &
                               ves_r, ves_z, &
                               ves_drdt, ves_drdp, ves_dzdt, ves_dzdp, &
                               plane_4point, plane_elev, plane_intersect

    implicit none

    INTEGER :: i, j, n, np
    REAL(8) :: vert_phi, vert_l0, vert_l, chi2 
    REAL(8) :: drdt, drdp, dzdt, dzdp, vr, vnx, vny, vnz
    REAL(8) :: elev1t, elev2t, elev3t, elev4t, elev5t, elev_top
    REAL(8) :: elev1b, elev2b, elev3b, elev4b, elev_bot
    REAL(8), allocatable :: vert_x(:,:), vert_y(:,:)
    REAL(8), allocatable :: vert_nx(:,:), vert_ny(:,:), vert_nz(:,:)
    REAL(8), allocatable :: vx(:,:), vy(:,:), vz(:,:)
    REAL(8), allocatable :: ves_nx(:,:), ves_ny(:,:), ves_nz(:,:)
    REAL(8), allocatable :: ves_theta(:,:), ves_phi(:,:)

    ! Initial status checks
    if (.not. plates_initialized) then
        stop 'trapezoid_properties: plates not yet initialized'
    end if
    if (trapezoids_initialized) then
        stop 'trapezoid_properties: trapezoids have already been initialized'
    end if

    ! Calculate vertex r and z coords if not explicitly supplied
    if (trim(vertex_mode) == 'theta_sep') then
        call vertex_theta_to_rz
    else if (trim(vertex_mode) == 'tlcfs_sep') then
        call lcfs_theta_to_ves_theta
        call vertex_theta_to_rz
    else if (trim(vertex_mode) == 'rz') then
        call vertex_rz_to_theta_sep
    else
        write(*,*) 'trapezoid_properties: unrecognized vertex_mode'
        stop
    end if

    ! Allocate variables according to the number of trapezoids
    allocate( vert_x(nPlates+1, nVertices),  vert_y(nPlates+1, nVertices),    &
              vert_nx(nPlates+1, nVertices), vert_ny(nPlates+1, nVertices),   &
              vert_nz(nPlates+1, nVertices), vx(nPlates+1, nVertices),        &
              vy(nPlates+1, nVertices),      vz(nPlates+1, nVertices),        &
              ves_nx(nPlates+1, nVertices),  ves_ny(nPlates+1, nVertices),    &
              ves_nz(nPlates+1, nVertices),  ves_theta(nPlates+1, nVertices), &
              ves_phi(nPlates+1, nVertices) )
    nTrapezoids = nPlates * (nVertices-1)
    allocate( traps(nTrapezoids) )

    write(*,*) '    Arrays allocated'

    ! Determine orthogonal projections of each vertex onto the vessel
    do i = 1, nPlates+1
        if (i == 1) then
            vert_phi = plate_phi(i) - 0.5 * plate_dphi(i)
        else
            vert_phi = plate_phi(i-1) + 0.5 * plate_dphi(i-1)
        end if
        do j = 1, nVertices
            vert_x(i,j) = vert_r(i,j)*cos(vert_phi)
            vert_y(i,j) = vert_r(i,j)*sin(vert_phi)
            vert_l0 = 0.0
            call ves_perp_intersect_3d(vert_x(i,j), vert_y(i,j), vert_z(i,j), &
                     vert_l0, vert_theta(i,j), vert_phi, &
                     vert_l, ves_theta(i,j), ves_phi(i,j), &
                     vx(i,j), vy(i,j), vz(i,j), &
                     ves_nx(i,j), ves_ny(i,j), ves_nz(i,j), chi2)
        end do
    end do

    write(*,*) '    Intersections determined for trapezoid bounds'

    ! Determine the trapezoid properties
    n = 0
    trap_err_count = 0
    do i = 1, nPlates
        do j = 1, nVertices-1

            n = n + 1

            ! Reference point on vessel (from which normal vector originates)
            traps(n)%ves_theta = 0.5 * ( ves_theta(i,j) + ves_theta(i,j+1) )
            traps(n)%ves_phi   = 0.5 * ( ves_phi(i,j)   + ves_phi(i+1,j)   )
            vr = ves_r(traps(n)%ves_theta, traps(n)%ves_phi)
            traps(n)%vx = vr * cos(traps(n)%ves_phi)
            traps(n)%vy = vr * sin(traps(n)%ves_phi)
            traps(n)%vz = ves_z(traps(n)%ves_theta, traps(n)%ves_phi)

            ! Normal vector for the top/bottom bounding surfs (& magnet stack)
            drdt = ves_drdt(traps(n)%ves_theta, traps(n)%ves_phi)
            dzdt = ves_dzdt(traps(n)%ves_theta, traps(n)%ves_phi)
            drdp = ves_drdp(traps(n)%ves_theta, traps(n)%ves_phi)
            dzdp = ves_dzdp(traps(n)%ves_theta, traps(n)%ves_phi)
            call ves_unorm(traps(n)%ves_phi, vr, drdt, drdp, dzdt, dzdp, &
                           vnx, vny, vnz)
            traps(n)%snx = -vnx
            traps(n)%sny = -vny
            traps(n)%snz = -vnz

            ! Determine the lateral bounding planes for the trapezoid
            !call bounding_planes_avg_4pt(i, j, vert_x, vert_y, vert_z, &
            !                             vx, vy, vz, traps(n))
            call bounding_planes_lateral(i, j, vert_x, vert_y, vert_z, &
                                         vx, vy, vz, traps(n))

            ! Reference point for top plane (based on highest vessel point)
            elev1t = plane_elev(vx(i,j),     vy(i,j),     vz(i,j),            &
                                traps(n)%vx, traps(n)%vy, traps(n)%vz,        &
                                vnx, vny, vnz)
            elev2t = plane_elev(vx(i,j+1),     vy(i,j+1),     vz(i,j+1),      &
                                traps(n)%vx,   traps(n)%vy,   traps(n)%vz,    &
                                vnx, vny, vnz)
            elev3t = plane_elev(vx(i+1,j),     vy(i+1,j),     vz(i+1,j),      &
                                traps(n)%vx,   traps(n)%vy,   traps(n)%vz,    &
                                vnx, vny, vnz)
            elev4t = plane_elev(vx(i+1,j+1),   vy(i+1,j+1),   vz(i+1,j+1),    &
                                traps(n)%vx,    traps(n)%vy,    traps(n)%vz,  &
                                vnx, vny, vnz)
            elev5t = 0.0   ! The elevation of the reference point on the vessel
            elev_top = max(elev1t, elev2t, elev3t, elev4t, elev5t) + gap_vac
            traps(n)%otx = traps(n)%vx + vnx * elev_top
            traps(n)%oty = traps(n)%vy + vny * elev_top
            traps(n)%otz = traps(n)%vz + vnz * elev_top

            ! Reference point for base plane (from lowest vertex point)
            elev1b = plane_elev(vert_x(i,j), vert_y(i,j), vert_z(i,j),  &
                                traps(n)%vx, traps(n)%vy, traps(n)%vz,  &
                                vnx, vny, vnz)
            elev2b = plane_elev(vert_x(i,j+1), vert_y(i,j+1), vert_z(i,j+1), &
                                traps(n)%vx,   traps(n)%vy,   traps(n)%vz,   &
                                vnx, vny, vnz)
            elev3b = plane_elev(vert_x(i+1,j), vert_y(i+1,j), vert_z(i+1,j), &
                                traps(n)%vx,   traps(n)%vy,   traps(n)%vz,   &
                                vnx, vny, vnz)
            elev4b = plane_elev( &
                            vert_x(i+1,j+1), vert_y(i+1,j+1), vert_z(i+1,j+1), &
                            traps(n)%vx,    traps(n)%vy,    traps(n)%vz,       &
                            vnx, vny, vnz)
            elev_bot = min(elev1b, elev2b, elev3b, elev4b)
            traps(n)%obx = traps(n)%vx + vnx * elev_bot
            traps(n)%oby = traps(n)%vy + vny * elev_bot
            traps(n)%obz = traps(n)%vz + vnz * elev_bot

            ! Determine the eight corners of the trapezoid
            call trapezoid_corners(traps(n))

            ! With the trapezoid edge lines defined, re-compute allowable hgts
            call reset_trap_height(traps(n))

            ! Calculate the trapezoid volume and centroid
            call trap_volume_ctr(traps(n))

            ! Check for erroneous or undesirable properties
            call check_trapezoid(traps(n))
            if (traps(n)%err) trap_err_count = trap_err_count + 1

            ! Output for debugging
            if (n >= 0 .and. n <= 0) then
                write(*,fmt='(ES11.4, X, ES11.4, X, ES11.4, X, ' // &
                             'ES11.4, X, ES11.4, X, ES11.4, X, ' // &
                             'ES11.4, X, ES11.4, X, ES11.4, X, ' // &
                             'ES11.4, X, ES11.4, X, ES11.4, X, ' // &
                             'ES11.4, X, ES11.4, X, ES11.4, X, ' // &
                             'ES11.4, X, ES11.4, X, ES11.4, X, ' // &
                             'ES11.4, X, ES11.4, X, ES11.4, X, ' // &
                             'ES11.4, X, ES11.4, X, ES11.4 )') &
                      traps(n)%ctx1, traps(n)%cty1, traps(n)%ctz1, & 
                      traps(n)%ctx2, traps(n)%cty2, traps(n)%ctz2, &
                      traps(n)%ctx3, traps(n)%cty3, traps(n)%ctz3, &
                      traps(n)%ctx4, traps(n)%cty4, traps(n)%ctz4, &
                      traps(n)%cbx1, traps(n)%cby1, traps(n)%cbz1, & 
                      traps(n)%cbx2, traps(n)%cby2, traps(n)%cbz2, &
                      traps(n)%cbx3, traps(n)%cby3, traps(n)%cbz3, &
                      traps(n)%cbx4, traps(n)%cby4, traps(n)%cbz4
            end if
            if (n >= 0 .and. n <= 0) then
                write(*, fmt='(ES11.4, X, ES11.4, X, ES11.4, A)') &
                      traps(n)%ox, traps(n)%oy, traps(n)%oz, '; ...'
            end if

        end do
    end do

    deallocate(vert_x, vert_y, vert_nx, vert_ny, vert_nz, &
               vx, vy, vz, ves_nx, ves_ny, ves_nz, &
               ves_theta, ves_phi)

    trapezoids_initialized = .true.

end subroutine trapezoid_properties

!-------------------------------------------------------------------------------
! bounding_planes_avg_4pt(i, j, vert_x, vert_y, vert_z, ves_x, ves_y, ves_z, 
!                         trap)
!
! Calculates the parameters of the lateral planes of a trapezoidal box at a 
! given grid location according to the 4-point average technique (i.e., not
! accounting for the properties of adjacent trapezoids).
!
! The inputs i and j are assumed not to exceed the dimensions of the subsequent
! arrays in the argument list.
!-------------------------------------------------------------------------------
subroutine bounding_planes_avg_4pt(i, j, vert_x, vert_y, vert_z, &
                                   ves_x, ves_y, ves_z, trap)

    use magnet_set_globals, only: trapezoid, gap_trp

    implicit none

    integer, intent(IN) :: i, j
    REAL(8), intent(IN) :: vert_x(:,:), vert_y(:,:), vert_z(:,:)
    REAL(8), intent(IN) :: ves_x(:,:), ves_y(:,:), ves_z(:,:)
    type(trapezoid) :: trap

    ! Reference point & normal vector, toroidal back plane 
    call avg_plane_4pt(vert_x(i,j),  vert_y(i,j),   vert_z(i,j),   &
                      ves_x(i,j),    ves_y(i,j),    ves_z(i,j),       &
                      ves_x(i,j+1),  ves_y(i,j+1),  ves_z(i,j+1),     &
                      vert_x(i,j+1), vert_y(i,j+1), vert_z(i,j+1), &
                      trap%ntbx, trap%ntby, trap%ntbz, &
                      trap%otbx, trap%otby, trap%otbz)
    trap%otbx = trap%otbx + 0.5 * gap_trp * trap%ntbx
    trap%otby = trap%otby + 0.5 * gap_trp * trap%ntby
    trap%otbz = trap%otbz + 0.5 * gap_trp * trap%ntbz

    ! Reference point & normal vector, toroidal front plane 
    call avg_plane_4pt( &
                 vert_x(i+1,j),   vert_y(i+1,j),   vert_z(i+1,j),   &
                 vert_x(i+1,j+1), vert_y(i+1,j+1), vert_z(i+1,j+1), &
                 ves_x(i+1,j+1),  ves_y(i+1,j+1),  ves_z(i+1,j+1),     &
                 ves_x(i+1,j),    ves_y(i+1,j),    ves_z(i+1,j),       &
                 trap%ntfx, trap%ntfy, trap%ntfz, &
                 trap%otfx, trap%otfy, trap%otfz)
    trap%otfx = trap%otfx + 0.5 * gap_trp * trap%ntfx
    trap%otfy = trap%otfy + 0.5 * gap_trp * trap%ntfy
    trap%otfz = trap%otfz + 0.5 * gap_trp * trap%ntfz

    ! Reference point & normal vector, poloidal back plane 
    call avg_plane_4pt(vert_x(i,j),   vert_y(i,j),   vert_z(i,j),   &
                      vert_x(i+1,j), vert_y(i+1,j), vert_z(i+1,j), &
                      ves_x(i+1,j),  ves_y(i+1,j),  ves_z(i+1,j),     &
                      ves_x(i,j),    ves_y(i,j),    ves_z(i,j),       &
                      trap%npbx, trap%npby, trap%npbz, &
                      trap%opbx, trap%opby, trap%opbz)
    trap%opbx = trap%opbx + 0.5 * gap_trp * trap%npbx
    trap%opby = trap%opby + 0.5 * gap_trp * trap%npby
    trap%opbz = trap%opbz + 0.5 * gap_trp * trap%npbz

    ! Reference point & normal vector, poloidal front plane 
    call avg_plane_4pt( &
                vert_x(i,j+1),   vert_y(i,j+1),   vert_z(i,j+1),   &
                ves_x(i,j+1),       ves_y(i,j+1),       ves_z(i,j+1),       &
                ves_x(i+1,j+1),     ves_y(i+1,j+1),     ves_z(i+1,j+1),     &
                vert_x(i+1,j+1), vert_y(i+1,j+1), vert_z(i+1,j+1), &
                trap%npfx, trap%npfy, trap%npfz, &
                trap%opfx, trap%opfy, trap%opfz)
    trap%opfx = trap%opfx + 0.5 * gap_trp * trap%npfx
    trap%opfy = trap%opfy + 0.5 * gap_trp * trap%npfy
    trap%opfz = trap%opfz + 0.5 * gap_trp * trap%npfz

end subroutine bounding_planes_avg_4pt

!-------------------------------------------------------------------------------
! bounding_planes_lateral()
! 
! Determines parameters of the bounding planes for the toroidal and poloidal 
! sides of each trapezoidal box through calls to the plane_boxcar subroutine.
!
! Most of the work in this wrapper subroutine concerns selecting subsets of
! points in the vert_{x,y,z} and ves_{x,y,z} arrays to be used in the boxcar
! averaging, handling cases in which the input trapezoid is near the edges
! of the grid of planes and poloidal vertices.
!-------------------------------------------------------------------------------
subroutine bounding_planes_lateral(i_plate, j_vertex, vert_x, vert_y, vert_z, &
                                   ves_x, ves_y, ves_z, trap)

    use magnet_set_globals, only: trapezoid, gap_trp, nPlates, nVertices,  &
                                  pi, nfp, n_box_pol, n_box_tor, &
                                  traps_poloidally_periodic, &
                                  traps_toroidally_periodic

    implicit none

    integer, intent(IN) :: i_plate, j_vertex
    REAL(8), intent(IN) :: vert_x(:,:), vert_y(:,:), vert_z(:,:)
    REAL(8), intent(IN) :: ves_x(:,:), ves_y(:,:), ves_z(:,:)
    type(trapezoid) :: trap
    integer :: n_arr_pol_max, n_arr_pol, n_arr_tor_max, n_arr_tor
    integer :: ind_min_pol, ind_max_pol, ind_min_tor, ind_max_tor
    integer :: m, n, tor_index, disp, disp_min
    integer :: pol_index, pol_ind_offs, vert_index
    logical :: reflect
    REAL(8) :: phiref
    REAL(8), allocatable :: arr_base_pb_x(:), arr_base_pb_y(:), arr_base_pb_z(:)
    REAL(8), allocatable :: arr_base_pf_x(:), arr_base_pf_y(:), arr_base_pf_z(:)
    REAL(8), allocatable :: arr_base_tb_x(:), arr_base_tb_y(:), arr_base_tb_z(:)
    REAL(8), allocatable :: arr_base_tf_x(:), arr_base_tf_y(:), arr_base_tf_z(:)
    REAL(8), allocatable :: arr_top_pb_x(:),  arr_top_pb_y(:),  arr_top_pb_z(:)
    REAL(8), allocatable :: arr_top_pf_x(:),  arr_top_pf_y(:),  arr_top_pf_z(:)
    REAL(8), allocatable :: arr_top_tb_x(:),  arr_top_tb_y(:),  arr_top_tb_z(:)
    REAL(8), allocatable :: arr_top_tf_x(:),  arr_top_tf_y(:),  arr_top_tf_z(:)

    !---------------------------------------------------------------------------
    ! Part 0: Allocate arrays to store vertex points for plane boxcar averaging.
    !         Use the largest possible size for each array.
    !---------------------------------------------------------------------------

    n_arr_pol_max = 2*n_box_pol + 2
    n_arr_tor_max = 2*n_box_tor + 2

    allocate(arr_base_tb_x(n_arr_pol_max), arr_base_tb_y(n_arr_pol_max), &
             arr_base_tb_z(n_arr_pol_max), arr_base_tf_x(n_arr_pol_max), &
             arr_base_tf_y(n_arr_pol_max), arr_base_tf_z(n_arr_pol_max), &
             arr_top_tb_x(n_arr_pol_max), arr_top_tb_y(n_arr_pol_max),   &
             arr_top_tb_z(n_arr_pol_max), arr_top_tf_x(n_arr_pol_max),   &
             arr_top_tf_y(n_arr_pol_max), arr_top_tf_z(n_arr_pol_max),   &
             arr_base_pb_x(n_arr_tor_max), arr_base_pb_y(n_arr_tor_max), &
             arr_base_pb_z(n_arr_tor_max), arr_base_pf_x(n_arr_tor_max), &
             arr_base_pf_y(n_arr_tor_max), arr_base_pf_z(n_arr_tor_max), &
             arr_top_pb_x(n_arr_tor_max), arr_top_pb_y(n_arr_tor_max),   &
             arr_top_pb_z(n_arr_tor_max), arr_top_pf_x(n_arr_tor_max),   &
             arr_top_pf_y(n_arr_tor_max), arr_top_pf_z(n_arr_tor_max)     )

    !---------------------------------------------------------------------------
    ! Part 1: For the input trapezoid with a location on the grid defined
    !         by i_plane and j_vertex, determine the number of points to 
    !         average over in the poloidal dimension (for the toroidal 
    !         back/front planes) and in the toroidal dimension (for the 
    !         poloidal back/front planes)
    !
    !         In cases where periodicity is assumed, the number of points will
    !         always be the same. Otherwise, the number of points is reduced
    !         for trapezoids near the edge of the grid (depending on the
    !         n_boxcar variables).
    !---------------------------------------------------------------------------
    if (traps_poloidally_periodic) then
        n_arr_pol = n_arr_pol_max
    else
        ind_min_pol = min(0, j_vertex-n_box_pol)
        ind_max_pol = max(nVertices, j_vertex+1+n_box_pol)
        n_arr_pol = ind_max_pol - ind_min_pol + 1
    end if

    if (traps_toroidally_periodic) then
        n_arr_tor = n_arr_tor_max
    else
        ind_min_tor = min(0, i_plate-n_box_tor)
        ind_max_tor = max(nPlates+1, i_plate+1+n_box_tor)
        n_arr_tor = ind_max_tor - ind_min_tor + 1
    end if

    !---------------------------------------------------------------------------
    ! Part 2: Populate the arrays of points to be used for boxcar averaging
    !         for the toroidal back/front and poloidal back/front planes.
    !
    !         In cases where periodicity is assumed, boxcar points that lie
    !         "over the edge of the grid" are taken from stellarator-symmetric
    !         locations in the grid. Otherwise, the arrays of points are 
    !         truncated at the edge of the grid.
    !---------------------------------------------------------------------------
    do m = 1, n_arr_pol

        ! Determine indices of the points to be taken from the vert/top arrays
        if (traps_poloidally_periodic) then

            disp = (m-1) - n_box_pol
            if (j_vertex + disp <= 0) then
                pol_index = (nVertices-1) + j_vertex + disp
            else if (j_vertex + disp > nVertices) then
                pol_index = j_vertex + disp - (nVertices-1)
            else
                pol_index = j_vertex + disp
            end if

            if (pol_index > nVertices .or. pol_index <= 0) then
                stop 'bounding_planes_lateral: n_box_pol must be < nVertices'
            end if

        else

            pol_index = ind_min_pol + m - 1

        end if

        !if (i_plate == 1 .and. j_vertex < 10) then
        !    write(*,fmt='(A, I2, A, I2, A, I2)') &
        !        '    m = ', m, '; i_plate = ', i_plate, '; pol_index = ', &
        !        pol_index
        !end if

        ! Populate array parameters for the toroidal front and back planes
        arr_base_tb_x(m) = vert_x(i_plate, pol_index)
        arr_base_tb_y(m) = vert_y(i_plate, pol_index)
        arr_base_tb_z(m) = vert_z(i_plate, pol_index)
        arr_base_tf_x(m) = vert_x(i_plate+1, pol_index)
        arr_base_tf_y(m) = vert_y(i_plate+1, pol_index)
        arr_base_tf_z(m) = vert_z(i_plate+1, pol_index)
        arr_top_tb_x(m) = ves_x(i_plate, pol_index)
        arr_top_tb_y(m) = ves_y(i_plate, pol_index)
        arr_top_tb_z(m) = ves_z(i_plate, pol_index)
        arr_top_tf_x(m) = ves_x(i_plate+1, pol_index)
        arr_top_tf_y(m) = ves_y(i_plate+1, pol_index)
        arr_top_tf_z(m) = ves_z(i_plate+1, pol_index)

    end do

    do n = 1, n_arr_tor

        ! Determine indices of the points to be taken from the vert/top arrays
        if (traps_toroidally_periodic) then

            disp = (n-1) - n_box_tor
            if (i_plate + disp <= 0) then
                tor_index = -disp - (i_plate-2)
                pol_index = nVertices - (j_vertex-1)
                pol_ind_offs  = -1
                reflect = .true.
                phiref = 0.0
            else if (i_plate + disp > nPlates+1) then
                tor_index = nPlates+1 - (i_plate + disp - (nPlates+1))
                pol_index = nVertices - (j_vertex-1)
                pol_ind_offs = -1
                reflect = .true.
                phiref = pi/nfp
            else
                tor_index = i_plate + disp
                pol_index = j_vertex
                pol_ind_offs = 1
                reflect = .false.
                phiref = 0.0
            end if

            if (tor_index > (nPlates+1) .or. tor_index <= 0) then
                stop 'bounding_planes_lateral: n_box_tor must be < nPlates+1'
            end if

        else

            vert_index = ind_min_tor + m - 1

        end if

        if (i_plate == 1 .and. j_vertex < 10) then
            write(*,fmt='(A, I2, A, I2, A, I2, A, I2)') &
                '    n = ', n, '; tor_index= ', tor_index, '; pol_index = ', &
                pol_index, '; pol_ind_offs = ', pol_ind_offs
        end if

        ! Populate array parameters for the toroidal front and back planes
        call stell_reflect(reflect, phiref, tor_index, pol_index, &
                           vert_x, vert_y, vert_z, &
                           arr_base_pb_x(n), arr_base_pb_y(n), arr_base_pb_z(n))
        call stell_reflect(reflect, phiref, tor_index, pol_index+pol_ind_offs, &
                           vert_x, vert_y, vert_z, &
                           arr_base_pf_x(n), arr_base_pf_y(n), arr_base_pf_z(n))
        call stell_reflect(reflect, phiref, tor_index, pol_index, &
                           ves_x, ves_y, ves_z, &
                           arr_top_pb_x(n), arr_top_pb_y(n), arr_top_pb_z(n))
        call stell_reflect(reflect, phiref, tor_index, pol_index+pol_ind_offs, &
                           ves_x, ves_y, ves_z, &
                           arr_top_pf_x(n), arr_top_pf_y(n), arr_top_pf_z(n))
        !arr_base_pb_x(n) = vert_x(tor_index, pol_index)
        !arr_base_pb_y(n) = vert_y(tor_index, pol_index)
        !arr_base_pb_z(n) = sgn_z * vert_z(tor_index, pol_index)
        !arr_base_pf_x(n) = vert_x(tor_index, pol_index+pol_ind_offs)
        !arr_base_pf_y(n) = vert_y(tor_index, pol_index+pol_ind_offs)
        !arr_base_pf_z(n) = sgn_z * vert_z(tor_index, pol_index+pol_ind_offs)
        !arr_top_pb_x(n) = ves_x(tor_index, pol_index)
        !arr_top_pb_y(n) = ves_y(tor_index, pol_index)
        !arr_top_pb_z(n) = sgn_z * ves_z(tor_index, pol_index)
        !arr_top_pf_x(n) = ves_x(tor_index, pol_index+pol_ind_offs)
        !arr_top_pf_y(n) = ves_y(tor_index, pol_index+pol_ind_offs)
        !arr_top_pf_z(n) = sgn_z * ves_z(tor_index, pol_index+pol_ind_offs)

    end do

    !---------------------------------------------------------------------------
    ! Part 3: Calculate the plane parameters using the plane_boxcar subroutine,
    !         given the input points determined in the previous steps.
    !---------------------------------------------------------------------------

    ! Toroidal back plane
    call plane_boxcar(n_arr_pol,                                   &
                      arr_base_tb_x, arr_base_tb_y, arr_base_tb_z, &
                      arr_top_tb_x, arr_top_tb_y, arr_top_tb_z,    &
                      trap%ntbx, trap%ntby, trap%ntbz,             &
                      trap%otbx, trap%otby, trap%otbz               )
    trap%otbx = trap%otbx + 0.5 * gap_trp * trap%ntbx
    trap%otby = trap%otby + 0.5 * gap_trp * trap%ntby
    trap%otbz = trap%otbz + 0.5 * gap_trp * trap%ntbz

    ! Toroidal front plane
    call plane_boxcar(n_arr_pol,                                   &
                      arr_base_tf_x, arr_base_tf_y, arr_base_tf_z, &
                      arr_top_tf_x, arr_top_tf_y, arr_top_tf_z,    &
                      trap%ntfx, trap%ntfy, trap%ntfz,             &
                      trap%otfx, trap%otfy, trap%otfz               )
    trap%ntfx = -trap%ntfx
    trap%ntfy = -trap%ntfy
    trap%ntfz = -trap%ntfz
    trap%otfx = trap%otfx + 0.5 * gap_trp * trap%ntfx
    trap%otfy = trap%otfy + 0.5 * gap_trp * trap%ntfy
    trap%otfz = trap%otfz + 0.5 * gap_trp * trap%ntfz

    ! Poloidal back plane
    call plane_boxcar(n_arr_tor,                                   &
                      arr_base_pb_x, arr_base_pb_y, arr_base_pb_z, &
                      arr_top_pb_x, arr_top_pb_y, arr_top_pb_z,    &
                      trap%npbx, trap%npby, trap%npbz,             &
                      trap%opbx, trap%opby, trap%opbz               )
    trap%npbx = -trap%npbx
    trap%npby = -trap%npby
    trap%npbz = -trap%npbz
    trap%opbx = trap%opbx + 0.5 * gap_trp * trap%npbx
    trap%opby = trap%opby + 0.5 * gap_trp * trap%npby
    trap%opbz = trap%opbz + 0.5 * gap_trp * trap%npbz
         
    ! Poloidal front plane
    call plane_boxcar(n_arr_tor,                                   &
                      arr_base_pf_x, arr_base_pf_y, arr_base_pf_z, &
                      arr_top_pf_x, arr_top_pf_y, arr_top_pf_z,    &
                      trap%npfx, trap%npfy, trap%npfz,             &
                      trap%opfx, trap%opfy, trap%opfz               )
    trap%opfx = trap%opfx + 0.5 * gap_trp * trap%npfx
    trap%opfy = trap%opfy + 0.5 * gap_trp * trap%npfy
    trap%opfz = trap%opfz + 0.5 * gap_trp * trap%npfz
 
    deallocate(arr_base_tb_x, arr_base_tb_y, arr_base_tb_z, &
               arr_base_tf_x, arr_base_tf_y, arr_base_tf_z, &
               arr_top_tb_x, arr_top_tb_y, arr_top_tb_z,    &
               arr_top_tf_x, arr_top_tf_y, arr_top_tf_z,    &
               arr_base_pb_x, arr_base_pb_y, arr_base_pb_z, &
               arr_base_pf_x, arr_base_pf_y, arr_base_pf_z, &
               arr_top_pb_x, arr_top_pb_y, arr_top_pb_z,    &
               arr_top_pf_x, arr_top_pf_y, arr_top_pf_z      )

    if (i_plate == 1 .and. j_vertex < 10) then
        write(*,fmt='(A, F7.4, A, F7.4, A, F7.4, A, F7.4, A, F7.4, A, F7.4)') &
            '        otbx = ', trap%otbx, '; otby = ', trap%otby, &
            '; otbz = ', trap%otbz, '; ntbx = ', trap%ntbx, &
            '; ntby = ', trap%ntby, '; ntbz = ', trap%ntbz
    end if

end subroutine bounding_planes_lateral

!-------------------------------------------------------------------------------
! stell_reflect(reflect, i, j, x, y, z, xr, yr, zr)
!
! Wrapper function that selects elements from input coordinate arrays, 
! performing a stellarator-symmetric reflection if desired
!-------------------------------------------------------------------------------
subroutine stell_reflect(reflect, phiref, i, j, x, y, z, xout, yout, zout)

    implicit none

    logical, intent(IN) :: reflect
    integer, intent(IN) :: i, j
    REAL(8), intent(IN) :: phiref, x(:,:), y(:,:), z(:,:)
    REAL(8), intent(OUT) :: xout, yout, zout
    REAL(8) :: phi, r
    
    if (reflect) then

        r = sqrt(x(i,j)**2 + y(i,j)**2)
        phi = atan2(y(i,j), x(i,j))
        xout = r * cos(phiref + (phiref-phi))
        yout = r * sin(phiref + (phiref-phi))
        zout = -z(i,j)

    else

        xout = x(i,j)
        yout = y(i,j)
        zout = z(i,j)
 
    end if

end subroutine stell_reflect 
        
        
!-------------------------------------------------------------------------------
! plane_boxcar(n, base_x, base_y, base_z, top_x, top_y, top_z, 
!              nx, ny, nz, ox, oy, oz)
!
! Determines the normal vector and reference point for a trapezoid bounding
! plane by averaging preselected arrays of input points
!-------------------------------------------------------------------------------
subroutine plane_boxcar(n, base_x, base_y, base_z, top_x, top_y, top_z, &
                        nx, ny, nz, ox, oy, oz)

    use magnet_set_calc, only: cross_prod

    implicit none

    integer, intent(IN) :: n
    REAL(8), allocatable, intent(IN) :: base_x(:), base_y(:), base_z(:)
    REAL(8), allocatable, intent(IN) :: top_x(:), top_y(:), top_z(:)
    REAL(8), intent(OUT) :: nx, ny, nz, ox, oy, oz
    integer :: i
    REAL(8) :: radial_dx, radial_dy, radial_dz, norm_r
    REAL(8) :: angular_dx, angular_dy, angular_dz, norm_a, adx, ady, adz
    REAL(8) :: perp_x, perp_y, perp_z, norm_perp
    REAL(8), allocatable :: rdx(:), rdy(:), rdz(:)

    allocate(rdx(n), rdy(n), rdz(n))
    radial_dx  = 0.0; radial_dy  = 0.0; radial_dz  = 0.0
    angular_dx = 0.0; angular_dy = 0.0; angular_dz = 0.0
    ox         = 0.0; oy         = 0.0; oz         = 0.0

    do i = 1, n

        ! Construct the reference point as the average of all input points
        ox = ox + base_x(i) + top_x(i)
        oy = oy + base_y(i) + top_y(i)
        oz = oz + base_z(i) + top_z(i)

        ! Construct the radial vector as the avg. of top-to-base unit vectors
        rdx(i) = base_x(i) - top_x(i)
        rdy(i) = base_y(i) - top_y(i)
        rdz(i) = base_z(i) - top_z(i)

        norm_r = sqrt(rdx(i)**2 + rdy(i)**2 + rdz(i)**2)

        radial_dx = radial_dx + (rdx(i)/norm_r)
        radial_dy = radial_dy + (rdy(i)/norm_r)
        radial_dz = radial_dz + (rdz(i)/norm_r)

        if (i > 1) then

            ! Construct the angular vector as the avg. of angular unit vectors
            adx = (top_x(i) + 0.5*rdx(i)) - (top_x(i-1) + 0.5*rdx(i-1))
            ady = (top_y(i) + 0.5*rdy(i)) - (top_y(i-1) + 0.5*rdy(i-1))
            adz = (top_z(i) + 0.5*rdz(i)) - (top_z(i-1) + 0.5*rdz(i-1))

            norm_a = sqrt(adx**2 + ady**2 + adz**2)

            angular_dx = angular_dx + adx/norm_a
            angular_dy = angular_dy + ady/norm_a
            angular_dz = angular_dz + adz/norm_a

        end if

    end do

    ! Division for the average calculations
    ox = ox/(2*n)
    oy = oy/(2*n)
    oz = oz/(2*n)
    radial_dx = radial_dx/n
    radial_dy = radial_dy/n
    radial_dz = radial_dz/n
    angular_dx = angular_dx/(n-1)
    angular_dy = angular_dy/(n-1)
    angular_dz = angular_dz/(n-1)

    ! Plane-normal vectors
    call cross_prod(angular_dx, angular_dy, angular_dz, &
                    radial_dx,  radial_dy,  radial_dz, perp_x, perp_y, perp_z)
    norm_perp = sqrt(perp_x**2 + perp_y**2 + perp_z**2)
    nx = perp_x/norm_perp
    ny = perp_y/norm_perp
    nz = perp_z/norm_perp

    deallocate(rdx, rdy, rdz)

end subroutine plane_boxcar

!-------------------------------------------------------------------------------
! subdivide_trapezoids()
!
! Given an array of trapezoidal boxes, subdivides them perpendicular to the
! stack normal to create stacks of trapezoids that fit within each original 
! trapezoid
!-------------------------------------------------------------------------------
subroutine subdivide_trapezoids()

    use magnet_set_globals, only: trapezoids_initialized, nTrapezoids, traps, &
                                  trap_err_count, &
                                  subtraps_initialized, nTraps_per_stack, &
                                  nSubtraps, subtraps, &
                                  nMagnets_total, magnets, symm, coiltype, &
                                  magnets_initialized, m_bulk, gap_rad
    
    implicit none

    integer :: i, j, n
    REAL(8) :: height, sub_height, h_bot, h_top

    ! Initial status checks
    if (.not. trapezoids_initialized) then
        stop 'subdivide_trapezoids: trapezoids not yet initialized'
    end if
    if (subtraps_initialized) then
        stop 'subdivide_trapezoids: sub-trapezoids are already initialized'
    end if
    if (magnets_initialized) then
        stop 'subdivide_trapezoids: magnets are already initialized'
    end if
   
    ! Allocate the array of sub-trapezoids (only for non-erroneous trapezoids)
    nSubtraps = (nTrapezoids - trap_err_count) * nTraps_per_stack 
    nMagnets_total = nSubtraps
    allocate( subtraps(nSubtraps), magnets(nSubtraps) )

    n = 0
    do i = 1, nTrapezoids

        height = ( traps(i)%otx - traps(i)%obx ) * traps(i)%snx + &
                 ( traps(i)%oty - traps(i)%oby ) * traps(i)%sny + &
                 ( traps(i)%otz - traps(i)%obz ) * traps(i)%snz

        sub_height = height / nTraps_per_stack

        if (.not. traps(i)%err) then

        do j = 1, nTraps_per_stack
    
                n = n + 1
                subtraps(n) = traps(i)
    
                ! Referernce points for top/base planes of each sub-trapezoid
                h_bot = sub_height * (j-1) + 0.5 * gap_rad
                h_top = sub_height * j     - 0.5 * gap_rad
                subtraps(n)%obx = traps(i)%obx + traps(i)%snx * h_bot
                subtraps(n)%oby = traps(i)%oby + traps(i)%sny * h_bot
                subtraps(n)%obz = traps(i)%obz + traps(i)%snz * h_bot
                subtraps(n)%otx = traps(i)%obx + traps(i)%snx * h_top
                subtraps(n)%oty = traps(i)%oby + traps(i)%sny * h_top
                subtraps(n)%otz = traps(i)%obz + traps(i)%snz * h_top
    
                ! Update corner points, volume, and centr. for the sub-trapezoid
                call trapezoid_corners(subtraps(n))
                call trap_volume_ctr(subtraps(n))
                call check_trapezoid(subtraps(n))
    
                ! Populate magnet structure with relevant data
                magnets(n)%ht = sub_height - 2*gap_rad
                magnets(n)%wd = 0.0  ! not well-defined for trapezoidal magnets
                magnets(n)%lg = 0.0  ! not well-defined for trapezoidal magnets
                magnets(n)%vol = subtraps(n)%vol
                magnets(n)%ox = subtraps(n)%ox
                magnets(n)%oy = subtraps(n)%oy
                magnets(n)%oz = subtraps(n)%oz
                magnets(n)%nx = subtraps(n)%snx
                magnets(n)%ny = subtraps(n)%sny
                magnets(n)%nz = subtraps(n)%snz
                magnets(n)%lx = 0.0  ! not well-defined for trapezoidal magnets
                magnets(n)%ly = 0.0 
                magnets(n)%lz = 0.0
                magnets(n)%n_phi = atan2(magnets(n)%ny, magnets(n)%nx)
                magnets(n)%n_theta = &
                    atan2(sqrt(magnets(n)%nx**2 + magnets(n)%ny**2), &
                          magnets(n)%nz)
                magnets(n)%moment = m_bulk * magnets(n)%vol
                magnets(n)%coiltype = coiltype
                magnets(n)%symm = symm
    
        end do

    end if

    end do

    subtraps_initialized = .true.
    magnets_initialized = .true.

end subroutine subdivide_trapezoids

!-------------------------------------------------------------------------------
! trapezoid_corners(trap)
!
! Populates the coordinates of the eight corners of a trapezoid with defined
! normal vectors and reference points for each of its six faces.
!-------------------------------------------------------------------------------
subroutine trapezoid_corners(tp)

    use magnet_set_globals, only: trapezoid
    use magnet_set_calc,    only: plane_intersect

    implicit none

    type(trapezoid) :: tp
    
    ! Top face, corner 1: toroidal back face and poloidal back face
    call plane_intersect(tp%otx,  tp%oty,  tp%otz,  tp%snx,  tp%sny,  tp%snz,  &
                         tp%otbx, tp%otby, tp%otbz, tp%ntbx, tp%ntby, tp%ntbz, &
                         tp%opbx, tp%opby, tp%opbz, tp%npbx, tp%npby, tp%npbz, &
                         tp%ctx1, tp%cty1, tp%ctz1)

    ! Top face, corner 2: toroidal front face and poloidal back face
    call plane_intersect(tp%otx,  tp%oty,  tp%otz,  tp%snx,  tp%sny,  tp%snz,  &
                         tp%otfx, tp%otfy, tp%otfz, tp%ntfx, tp%ntfy, tp%ntfz, &
                         tp%opbx, tp%opby, tp%opbz, tp%npbx, tp%npby, tp%npbz, &
                         tp%ctx2, tp%cty2, tp%ctz2)

    ! Top face, corner 3: toroidal front face and poloidal front face
    call plane_intersect(tp%otx,  tp%oty,  tp%otz,  tp%snx,  tp%sny,  tp%snz,  &
                         tp%otfx, tp%otfy, tp%otfz, tp%ntfx, tp%ntfy, tp%ntfz, &
                         tp%opfx, tp%opfy, tp%opfz, tp%npfx, tp%npfy, tp%npfz, &
                         tp%ctx3, tp%cty3, tp%ctz3)

    ! Top face, corner 4: toroidal back face and poloidal front face
    call plane_intersect(tp%otx,  tp%oty,  tp%otz,  tp%snx,  tp%sny,  tp%snz,  &
                         tp%otbx, tp%otby, tp%otbz, tp%ntbx, tp%ntby, tp%ntbz, &
                         tp%opfx, tp%opfy, tp%opfz, tp%npfx, tp%npfy, tp%npfz, &
                         tp%ctx4, tp%cty4, tp%ctz4)

    ! Bottom face, corner 1: toroidal back face and poloidal back face
    call plane_intersect(tp%obx,  tp%oby,  tp%obz,  tp%snx,  tp%sny,  tp%snz,  &
                         tp%otbx, tp%otby, tp%otbz, tp%ntbx, tp%ntby, tp%ntbz, &
                         tp%opbx, tp%opby, tp%opbz, tp%npbx, tp%npby, tp%npbz, &
                         tp%cbx1, tp%cby1, tp%cbz1)

    ! Top face, corner 2: toroidal front face and poloidal back face
    call plane_intersect(tp%obx,  tp%oby,  tp%obz,  tp%snx,  tp%sny,  tp%snz,  &
                         tp%otfx, tp%otfy, tp%otfz, tp%ntfx, tp%ntfy, tp%ntfz, &
                         tp%opbx, tp%opby, tp%opbz, tp%npbx, tp%npby, tp%npbz, &
                         tp%cbx2, tp%cby2, tp%cbz2)

    ! Top face, corner 3: toroidal front face and poloidal front face
    call plane_intersect(tp%obx,  tp%oby,  tp%obz,  tp%snx,  tp%sny,  tp%snz,  &
                         tp%otfx, tp%otfy, tp%otfz, tp%ntfx, tp%ntfy, tp%ntfz, &
                         tp%opfx, tp%opfy, tp%opfz, tp%npfx, tp%npfy, tp%npfz, &
                         tp%cbx3, tp%cby3, tp%cbz3)

    ! Top face, corner 4: toroidal back face and poloidal front face
    call plane_intersect(tp%obx,  tp%oby,  tp%obz,  tp%snx,  tp%sny,  tp%snz,  &
                         tp%otbx, tp%otby, tp%otbz, tp%ntbx, tp%ntby, tp%ntbz, &
                         tp%opfx, tp%opfy, tp%opfz, tp%npfx, tp%npfy, tp%npfz, &
                         tp%cbx4, tp%cby4, tp%cbz4)

end subroutine trapezoid_corners

!-------------------------------------------------------------------------------
! query_points(np, d, x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, 
!              xq, yq, zq)
!
! Generates a set of query points for fitting a plane to a set of four line
! segments (defined by four "corner" points). 
!
! Input parameters:
!     integer :: np         -> number of query points per segment 
!     REAL(8) :: d          -> distance between query points along each segment
!     REAL(8) :: x1, y1, z1 -> x, y, and z coordinates of the first corner
!     REAL(8) :: x2, y2, z2 -> x, y, and z coordinates of the second corner
!     REAL(8) :: x3, y3, z3 -> x, y, and z coordinates of the third corner
!     REAL(8) :: x4, y4, z4 -> x, y, and z coordinates of the fourth corner
!
! Output parameters:
!     integer :: nq                       -> total number of query points (4*np)
!     REAL(8), allocatable :: xq, yq, zq  -> arrays of x, y, and z coords of
!                                            the query points
!-------------------------------------------------------------------------------

subroutine query_points(np, d, x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, &
                        nq, xq, yq, zq)

    implicit none

    integer, intent(IN) :: np
    REAL(8), intent(IN) :: d
    REAL(8), intent(IN) :: x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4
    REAL(8), allocatable :: xq(:), yq(:), zq(:)
    integer, intent(OUT) :: nq
    integer :: i
    REAL(8) :: dx1, dy1, dz1, dx2, dy2, dz2, dx3, dy3, dz3, dx4, dy4, dz4
    REAL(8) :: nx1, ny1, nz1, nx2, ny2, nz2, nx3, ny3, nz3, nx4, ny4, nz4
    REAL(8) :: norm_1, norm_2, norm_3, norm_4

    ! Compute normal vectors for each segment
    dx1 = x2 - x1
    dy1 = y2 - y1
    dz1 = z2 - z1
    norm_1 = sqrt(dx1**2 + dy1**2 + dz1**2)
    nx1 = dx1/norm_1
    ny1 = dy1/norm_1
    nz1 = dz1/norm_1
    
    dx2 = x3 - x2
    dy2 = y3 - y2
    dz2 = z3 - z2
    norm_2 = sqrt(dx2**2 + dy2**2 + dz2**2)
    nx2 = dx2/norm_2
    ny2 = dy2/norm_2
    nz2 = dz2/norm_2

    dx3 = x4 - x3
    dy3 = y4 - y3
    dz3 = z4 - z3
    norm_3 = sqrt(dx3**2 + dy3**2 + dz3**2)
    nx3 = dx3/norm_3
    ny3 = dy3/norm_3
    nz3 = dz3/norm_3

    dx4 = x1 - x4
    dy4 = y1 - y4
    dz4 = z1 - z4
    norm_4 = sqrt(dz4**2 + dy4**2 + dz4**2)
    nx4 = dx4/norm_4
    ny4 = dy4/norm_4
    nz4 = dz4/norm_4

    ! Generate query points along each segment
    do i = 1, np
        xq(i       ) = x1 + (i-1)*d*nx1
        xq(i +   np) = x2 + (i-1)*d*nx2
        xq(i + 2*np) = x3 + (i-1)*d*nx3
        xq(i + 3*np) = x4 + (i-1)*d*nx4
        yq(i       ) = y1 + (i-1)*d*ny1
        yq(i +   np) = y2 + (i-1)*d*ny2
        yq(i + 2*np) = y3 + (i-1)*d*ny3
        yq(i + 3*np) = y4 + (i-1)*d*ny4
        zq(i       ) = z1 + (i-1)*d*nz1
        zq(i +   np) = z2 + (i-1)*d*nz2
        zq(i + 2*np) = z3 + (i-1)*d*nz3
        zq(i + 3*np) = z4 + (i-1)*d*nz4
    end do

    nq = np * 4

end subroutine query_points


!-------------------------------------------------------------------------------
! avg_plane_4pt(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, 
!               nx, ny, nz, ox, oy, oz)
!
! Computes an "average plane" betwen four points based on the four unique
! planes corresponding to the four unique combinations of three points.
! The normal vector of the "average plane" is computed from the arithmetic means
! of the components of the normal vectors of the four combination planes. The 
! reference point is the centroid of the four input points.
!
! The input points are assumed to be given in perimeter order; this determines
! the direction of the normal vector.
!
! Input parameters
!     REAL(8) :: x1, y1, z1  -> x, y, and z coordinates of the first point
!     REAL(8) :: x2, y2, z2  -> x, y, and z coordinates of the second point
!     REAL(8) :: x3, y3, z3  -> x, y, and z coordinates of the third point
!     REAL(8) :: x4, y4, z4  -> x, y, and z coordinates of the fourth point
!
! Output parameters
!     REAL(8) :: nx, ny, nz  -> x, y, and z components of the avg plane normal
!     REAL(8) :: ox, oy, oz  -> x, y, and z coordinates of the reference point
!-------------------------------------------------------------------------------
subroutine avg_plane_4pt(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, &
                         nx, ny, nz, ox, oy, oz)

    use magnet_set_calc, only: plane_3point

    implicit none

    REAL(8), intent(IN)  :: x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4
    REAL(8), intent(OUT) :: nx, ny, nz, ox, oy, oz
    REAL(8) :: nx1, ny1, nz1, nx2, ny2, nz2, nx3, ny3, nz3, nx4, ny4, nz4, norm
    
    ! Normal vectors of the four subsets of three points
    call plane_3point(x3, y3, z3, x2, y2, z2, x1, y1, z1, nx1, ny1, nz1)
    call plane_3point(x4, y4, z4, x3, y3, z3, x1, y1, z1, nx2, ny2, nz2)
    call plane_3point(x1, y1, z1, x4, y4, z4, x2, y2, z2, nx3, ny3, nz3)
    call plane_3point(x4, y4, z4, x3, y3, z3, x2, y2, z2, nx4, ny4, nz4)

    ! Mean normal vector
    nx = 0.25*(nx1 + nx2 + nx3 + nx4)
    ny = 0.25*(ny1 + ny2 + ny3 + ny4)
    nz = 0.25*(nz1 + nz2 + nz3 + nz4)
    norm = sqrt(nx**2 + ny**2 + nz**2)
    nx = nx / norm
    ny = ny / norm
    nz = nz / norm

    if (abs(nx**2 + ny**2 + nz**2 - 1.0) > 1e-10) then
        write(*,*) 'Normal vector not unity'
    end if

    ! Centroid as reference point
    ox = 0.25*(x1 + x2 + x3 + x4)
    oy = 0.25*(y1 + y2 + y3 + y4)
    oz = 0.25*(z1 + z2 + z3 + z4)

end subroutine avg_plane_4pt

!-------------------------------------------------------------------------------
! reset_trap_height(tp)
!
! Determines the allowable height for a trapezoidal box to avoid collision with
! the vessel, based on its edge lines
!-------------------------------------------------------------------------------
subroutine reset_trap_height(tp)

    use magnet_set_globals, only: ves_r00, gap_vac, trapezoid
    use magnet_set_calc,    only: ves_dist_3d, plane_elev

    implicit none

    type(trapezoid) :: tp
    REAL(8) :: dx1, dy1, dz1, dx2, dy2, dz2, dx3, dy3, dz3, dx4, dy4, dz4
    REAL(8) :: ax1, ay1, az1, ax2, ay2, az2, ax3, ay3, az3, ax4, ay4, az4
    REAL(8) :: norm_1, norm_2, norm_3, norm_4
    REAL(8) :: length_1, length_2, length_3, length_4
    REAL(8) :: l0, theta0, phi0, theta, phi, chi2
    REAL(8) :: xv1, yv1, zv1, xv2, yv2, zv2, xv3, yv3, zv3, xv4, yv4, zv4
    REAL(8) :: elev1t, elev2t, elev3t, elev4t, elev5t, elev_top

    ! Directional unit vectors for edge 1
    dx1 = tp%cbx1 - tp%ctx1
    dy1 = tp%cby1 - tp%cty1
    dz1 = tp%cbz1 - tp%ctz1
    norm_1 = sqrt(dx1**2 + dy1**2 + dz1**2)
    ax1 = dx1 / norm_1
    ay1 = dy1 / norm_1
    az1 = dz1 / norm_1
    
    ! Directional unit vectors for edge 2
    dx2 = tp%cbx2 - tp%ctx2
    dy2 = tp%cby2 - tp%cty2
    dz2 = tp%cbz2 - tp%ctz2
    norm_2 = sqrt(dx2**2 + dy2**2 + dz2**2)
    ax2 = dx2 / norm_2
    ay2 = dy2 / norm_2
    az2 = dz2 / norm_2
    
    ! Directional unit vectors for edge 3
    dx3 = tp%cbx3 - tp%ctx3
    dy3 = tp%cby3 - tp%cty3
    dz3 = tp%cbz3 - tp%ctz3
    norm_3 = sqrt(dx3**2 + dy3**2 + dz3**2)
    ax3 = dx3 / norm_3
    ay3 = dy3 / norm_3
    az3 = dz3 / norm_3
    
    ! Directional unit vectors for edge 4
    dx4 = tp%cbx4 - tp%ctx4
    dy4 = tp%cby4 - tp%cty4
    dz4 = tp%cbz4 - tp%ctz4
    norm_4 = sqrt(dx4**2 + dy4**2 + dz4**2)
    ax4 = dx4 / norm_4
    ay4 = dy4 / norm_4
    az4 = dz4 / norm_4
    
    ! Determine the intersection points on the vessel for each edge
    l0 = 0.0
    phi0 = atan2(tp%cty1, tp%ctx1)
    theta0 = atan2(tp%ctz1, sqrt(tp%ctx1**2 + tp%cty1**2) - ves_r00)
    call ves_dist_3d(tp%ctx1, tp%cty1, tp%ctz1, ax1, ay1, az1, l0, theta0, &
                     phi0, length_1, theta, phi, xv1, yv1, zv1, chi2)
    l0 = 0.0
    phi0 = atan2(tp%cty2, tp%ctx2)
    theta0 = atan2(tp%ctz2, sqrt(tp%ctx2**2 + tp%cty2**2) - ves_r00)
    call ves_dist_3d(tp%ctx2, tp%cty2, tp%ctz2, ax2, ay2, az2, l0, theta0, &
                     phi0, length_2, theta, phi, xv2, yv2, zv2, chi2)
    l0 = 0.0
    phi0 = atan2(tp%cty3, tp%ctx3)
    theta0 = atan2(tp%ctz3, sqrt(tp%ctx3**2 + tp%cty3**2) - ves_r00)
    call ves_dist_3d(tp%ctx3, tp%cty3, tp%ctz3, ax3, ay3, az3, l0, theta0, &
                     phi0, length_3, theta, phi, xv3, yv3, zv3, chi2)
    l0 = 0.0
    phi0 = atan2(tp%cty4, tp%ctx4)
    theta0 = atan2(tp%ctz4, sqrt(tp%ctx4**2 + tp%cty4**2) - ves_r00)
    call ves_dist_3d(tp%ctx4, tp%cty4, tp%ctz4, ax4, ay4, az4, l0, theta0, &
                     phi0, length_4, theta, phi, xv4, yv4, zv4, chi2)

    ! Evaluate the elevations of each vessel point above the vessel
    elev1t = plane_elev(xv1, yv1, zv1, tp%vx, tp%vy, tp%vz,  &
                        -tp%snx, -tp%sny, -tp%snz)
    elev2t = plane_elev(xv2, yv2, zv2, tp%vx, tp%vy, tp%vz,  &
                        -tp%snx, -tp%sny, -tp%snz)
    elev3t = plane_elev(xv3, yv3, zv3, tp%vx, tp%vy, tp%vz,  &
                        -tp%snx, -tp%sny, -tp%snz)
    elev4t = plane_elev(xv4, yv4, zv4, tp%vx, tp%vy, tp%vz,  &
                        -tp%snx, -tp%sny, -tp%snz)
    elev5t = 0.0   ! The elevation of the reference point on the vessel
    elev_top = max(elev1t, elev2t, elev3t, elev4t, elev5t) + gap_vac

    ! Reset the reference point for the top of the trapezoid/stack
    tp%otx = tp%vx - tp%snx * elev_top
    tp%oty = tp%vy - tp%sny * elev_top
    tp%otz = tp%vz - tp%snz * elev_top

    ! Re-calculate the corners according to the adjusted top plane height
    call trapezoid_corners(tp)

end subroutine reset_trap_height

!-------------------------------------------------------------------------------
! trap_volume_ctr(tp)
!
! Calculates the volume enclosed in a trapezoid, as well as the centroid.
! 
! The trapezoid structure tp is taken as an input. All 8 corner coordinates are
! assumed to be defined within the structure prior to the call to this 
! subroutine. During execution, the following parameters of tp will be updated:
!     tp%vol: the enclosed volume
!     tp%ox:  the x coordinate of the centroid
!     tp%oy:  the y coordinate of the centroid
!     tp%oz:  the z coordinate of the centroid
! 
! Based on the formulas for volume and centroid in lecture notes 
! by R. Nuernberg, Imperial College of London, accessed 2019-11-03 from
! wwwf.imperial.ac.uk/~rn/centroid.pdf
!-------------------------------------------------------------------------------
subroutine trap_volume_ctr(tp)

    use magnet_set_globals, only: trapezoid
    use magnet_set_calc,    only: cross_prod

    implicit none

    type(trapezoid) :: tp
    integer :: i
    REAL(8) :: ntx1, nty1, ntz1, ntx2, nty2, ntz2
    REAL(8) :: nbx1, nby1, nbz1, nbx2, nby2, nbz2
    REAL(8) :: ntbx1, ntby1, ntbz1, ntbx2, ntby2, ntbz2
    REAL(8) :: ntfx1, ntfy1, ntfz1, ntfx2, ntfy2, ntfz2
    REAL(8) :: npbx1, npby1, npbz1, npbx2, npby2, npbz2
    REAL(8) :: npfx1, npfy1, npfz1, npfx2, npfy2, npfz2
    REAL(8) :: trip_prod(12), quadsum_x(12), quadsum_y(12), quadsum_z(12)
    
    ! Normal vectors (not unit length) for the two triangles composing each face
    !
    !            Poloidal
    !             front
    !          4----------3 
    !          |   Base   |
    ! Toroidal |  (view   | Toroidal    
    !   back   |  toward  |  front     
    !          |  vessel) |
    !          1----------2
    !            Poloidal
    !              back

    ! Base face, triangle 1: base2-base1 x base4-base1
    call cross_prod( tp%cbx2 - tp%cbx1, tp%cby2 - tp%cby1, tp%cbz2 - tp%cbz1, &
                     tp%cbx4 - tp%cbx1, tp%cby4 - tp%cby1, tp%cbz4 - tp%cbz1, &
                     nbx1, nby1, nbz1 )
    trip_prod(1) = tp%cbx1 *  nbx1 + tp%cby1 *  nby1 + tp%cbz1 *  nbz1
    quadsum_x(1) = nbx1 * ((tp%cbx1 + tp%cbx2)**2 + (tp%cbx1 + tp%cbx4)**2 &
                                                  + (tp%cbx2 + tp%cbx4)**2  )
    quadsum_y(1) = nby1 * ((tp%cby1 + tp%cby2)**2 + (tp%cby1 + tp%cby4)**2 &
                                                  + (tp%cby2 + tp%cby4)**2  )
    quadsum_z(1) = nbz1 * ((tp%cbz1 + tp%cbz2)**2 + (tp%cbz1 + tp%cbz4)**2 &
                                                  + (tp%cbz2 + tp%cbz4)**2  )

    ! Base face, triangle 2: base4-base3 x base2-base3
    call cross_prod( tp%cbx4 - tp%cbx3, tp%cby4 - tp%cby3, tp%cbz4 - tp%cbz3, &
                     tp%cbx2 - tp%cbx3, tp%cby2 - tp%cby3, tp%cbz2 - tp%cbz3, &
                     nbx2, nby2, nbz2 )
    trip_prod(2) = tp%cbx3 *  nbx2 + tp%cby3 *  nby2 + tp%cbz3 *  nbz2
    quadsum_x(2) = nbx2 * ((tp%cbx2 + tp%cbx3)**2 + (tp%cbx2 + tp%cbx4)**2 &
                                                  + (tp%cbx3 + tp%cbx4)**2  )
    quadsum_y(2) = nby2 * ((tp%cby2 + tp%cby3)**2 + (tp%cby2 + tp%cby4)**2 &
                                                  + (tp%cby3 + tp%cby4)**2  )
    quadsum_z(2) = nbz2 * ((tp%cbz2 + tp%cbz3)**2 + (tp%cbz2 + tp%cbz4)**2 &
                                                  + (tp%cbz3 + tp%cbz4)**2  )

    ! Poloidal back face, triangle 1: top1-base1 x base2-base1
    call cross_prod( tp%ctx1 - tp%cbx1, tp%cty1 - tp%cby1, tp%ctz1 - tp%cbz1, &
                     tp%cbx2 - tp%cbx1, tp%cby2 - tp%cby1, tp%cbz2 - tp%cbz1, &
                     npbx1, npby1, npbz1 )
    trip_prod(3) = tp%cbx1 * npbx1 + tp%cby1 * npby1 + tp%cbz1 * npbz1
    quadsum_x(3) = npbx1 * ((tp%ctx1 + tp%cbx2)**2 + (tp%ctx1 + tp%cbx1)**2 &
                                                   + (tp%cbx2 + tp%cbx1)**2  )
    quadsum_y(3) = npby1 * ((tp%cty1 + tp%cby2)**2 + (tp%cty1 + tp%cby1)**2 &
                                                   + (tp%cby2 + tp%cby1)**2  )
    quadsum_z(3) = npbz1 * ((tp%ctz1 + tp%cbz2)**2 + (tp%ctz1 + tp%cbz1)**2 &
                                                   + (tp%cbz2 + tp%cbz1)**2  )

    ! Poloidal back face, triangle 2: base2-top2 x top1-top2
    call cross_prod( tp%cbx2 - tp%ctx2, tp%cby2 - tp%cty2, tp%cbz2 - tp%ctz2, &
                     tp%ctx1 - tp%ctx2, tp%cty1 - tp%cty2, tp%ctz1 - tp%ctz2, &
                     npbx2, npby2, npbz2 )
    trip_prod(4) = tp%ctx2 * npbx2 + tp%cty2 * npby2 + tp%ctz2 * npbz2
    quadsum_x(4) = npbx2 * ((tp%ctx1 + tp%cbx2)**2 + (tp%ctx1 + tp%ctx2)**2 &
                                                   + (tp%cbx2 + tp%ctx2)**2  )
    quadsum_y(4) = npby2 * ((tp%cty1 + tp%cby2)**2 + (tp%cty1 + tp%cty2)**2 &
                                                   + (tp%cby2 + tp%cty2)**2  )
    quadsum_z(4) = npbz2 * ((tp%ctz1 + tp%cbz2)**2 + (tp%ctz1 + tp%ctz2)**2 &
                                                   + (tp%cbz2 + tp%ctz2)**2  )

    ! Poloidal front face, triangle 1: top3-base3 x base4-base3
    call cross_prod( tp%ctx3 - tp%cbx3, tp%cty3 - tp%cby3, tp%ctz3 - tp%cbz3, &
                     tp%cbx4 - tp%cbx3, tp%cby4 - tp%cby3, tp%cbz4 - tp%cbz3, &
                     npfx1, npfy1, npfz1 )
    trip_prod(5) = tp%cbx3 * npfx1 + tp%cby3 * npfy1 + tp%cbz3 * npfz1
    quadsum_x(5) = npfx1 * ((tp%ctx3 + tp%cbx3)**2 + (tp%ctx3 + tp%cbx4)**2 &
                                                   + (tp%cbx3 + tp%cbx4)**2  )
    quadsum_y(5) = npfy1 * ((tp%cty3 + tp%cby3)**2 + (tp%cty3 + tp%cby4)**2 &
                                                   + (tp%cby3 + tp%cby4)**2  )
    quadsum_z(5) = npfz1 * ((tp%ctz3 + tp%cbz3)**2 + (tp%ctz3 + tp%cbz4)**2 &
                                                   + (tp%cbz3 + tp%cbz4)**2  )

    ! Poloidal front face, triangle 2: base4-top4 x top3-top4
    call cross_prod( tp%cbx4 - tp%ctx4, tp%cby4 - tp%cty4, tp%cbz4 - tp%ctz4, &
                     tp%ctx3 - tp%ctx4, tp%cty3 - tp%cty4, tp%ctz3 - tp%ctz4, &
                     npfx2, npfy2, npfz2 )
    trip_prod(6) = tp%ctx4 * npfx2 + tp%cty4 * npfy2 + tp%ctz4 * npfz2
    quadsum_x(6) = npfx2 * ((tp%ctx3 + tp%ctx4)**2 + (tp%ctx3 + tp%cbx4)**2 &
                                                   + (tp%ctx4 + tp%cbx4)**2  )
    quadsum_y(6) = npfy2 * ((tp%cty3 + tp%cty4)**2 + (tp%cty3 + tp%cby4)**2 &
                                                   + (tp%cty4 + tp%cby4)**2  )
    quadsum_z(6) = npfz2 * ((tp%ctz3 + tp%ctz4)**2 + (tp%ctz3 + tp%cbz4)**2 &
                                                   + (tp%ctz4 + tp%cbz4)**2  )

    ! Toroidal back face, triangle 1: base4-base1 x top1-base1
    call cross_prod( tp%cbx4 - tp%cbx1, tp%cby4 - tp%cby1, tp%cbz4 - tp%cbz1, &
                     tp%ctx1 - tp%cbx1, tp%cty1 - tp%cby1, tp%ctz1 - tp%cbz1, &
                     ntbx1, ntby1, ntbz1 )
    trip_prod(7) = tp%cbx1 * ntbx1 + tp%cby1 * ntby1 + tp%cbz1 * ntbz1
    quadsum_x(7) = ntbx1 * ((tp%cbx1 + tp%ctx1)**2 + (tp%cbx1 + tp%cbx4)**2 &
                                                   + (tp%ctx1 + tp%cbx4)**2  )
    quadsum_y(7) = ntby1 * ((tp%cby1 + tp%cty1)**2 + (tp%cby1 + tp%cby4)**2 &
                                                   + (tp%cty1 + tp%cby4)**2  )
    quadsum_z(7) = ntbz1 * ((tp%cbz1 + tp%ctz1)**2 + (tp%cbz1 + tp%cbz4)**2 &
                                                   + (tp%ctz1 + tp%cbz4)**2  )

    ! Toroidal back face, triangle 2: top1-top4 x base4-top4 
    call cross_prod( tp%ctx1 - tp%ctx4, tp%cty1 - tp%cty4, tp%ctz1 - tp%ctz4, &
                     tp%cbx4 - tp%ctx4, tp%cby4 - tp%cty4, tp%cbz4 - tp%ctz4, &
                     ntbx2, ntby2, ntbz2 )
    trip_prod(8) = tp%ctx4 * ntbx2 + tp%cty4 * ntby2 + tp%ctz4 * ntbz2
    quadsum_x(8) = ntbx2 * ((tp%ctx1 + tp%ctx4)**2 + (tp%ctx1 + tp%cbx4)**2 &
                                                   + (tp%ctx4 + tp%cbx4)**2  )
    quadsum_y(8) = ntby2 * ((tp%cty1 + tp%cty4)**2 + (tp%cty1 + tp%cby4)**2 &
                                                   + (tp%cty4 + tp%cby4)**2  )
    quadsum_z(8) = ntbz2 * ((tp%ctz1 + tp%ctz4)**2 + (tp%ctz1 + tp%cbz4)**2 &
                                                   + (tp%ctz4 + tp%cbz4)**2  )

    ! Toroidal front face, triangle 1: top2-base2 x base3-base2
    call cross_prod( tp%ctx2 - tp%cbx2, tp%cty2 - tp%cby2, tp%ctz2 - tp%cbz2, &
                     tp%cbx3 - tp%cbx2, tp%cby3 - tp%cby2, tp%cbz3 - tp%cbz2, &
                     ntfx1, ntfy1, ntfz1 )
    trip_prod(9) = tp%cbx2 * ntfx1 + tp%cby2 * ntfy1 + tp%cbz2 * ntfz1
    quadsum_x(9) = ntfx1 * ((tp%ctx2 + tp%cbx2)**2 + (tp%ctx2 + tp%cbx3)**2 &
                                                   + (tp%cbx2 + tp%cbx3)**2  )
    quadsum_y(9) = ntfy1 * ((tp%cty2 + tp%cby2)**2 + (tp%cty2 + tp%cby3)**2 &
                                                   + (tp%cby2 + tp%cby3)**2  )
    quadsum_z(9) = ntfz1 * ((tp%ctz2 + tp%cbz2)**2 + (tp%ctz2 + tp%cbz3)**2 &
                                                   + (tp%cbz2 + tp%cbz3)**2  )

    ! Toroidal front face, triangle 2: base3-top3 x top2-top3
    call cross_prod( tp%cbx3 - tp%ctx3, tp%cby3 - tp%cty3, tp%cbz3 - tp%ctz3, &
                     tp%ctx2 - tp%ctx3, tp%cty2 - tp%cty3, tp%ctz2 - tp%ctz3, &
                     ntfx2, ntfy2, ntfz2 )
    trip_prod(10) = tp%ctx3 * ntfx2 + tp%cty3 * ntfy2 + tp%ctz3 * ntfz2
    quadsum_x(10) = ntfx2 * ((tp%cbx3 + tp%ctx2)**2 + (tp%cbx3 + tp%ctx3)**2 &
                                                    + (tp%ctx2 + tp%ctx3)**2  )
    quadsum_y(10) = ntfy2 * ((tp%cby3 + tp%cty2)**2 + (tp%cby3 + tp%cty3)**2 &
                                                    + (tp%cty2 + tp%cty3)**2  )
    quadsum_z(10) = ntfz2 * ((tp%cbz3 + tp%ctz2)**2 + (tp%cbz3 + tp%ctz3)**2 &
                                                    + (tp%ctz2 + tp%ctz3)**2  )

    ! Top face, triangle 1: top4-top1 x top2-top1
    call cross_prod( tp%ctx4 - tp%ctx1, tp%cty4 - tp%cty1, tp%ctz4 - tp%ctz1, &
                     tp%ctx2 - tp%ctx1, tp%cty2 - tp%cty1, tp%ctz2 - tp%ctz1, &
                     ntx1, nty1, ntz1 )
    trip_prod(11) = tp%ctx1 *  ntx1 + tp%cty1 *  nty1 + tp%ctz1 *  ntz1
    quadsum_x(11) = ntx1 * ((tp%ctx1 + tp%ctx2)**2 + (tp%ctx1 + tp%ctx4)**2 &
                                                   + (tp%ctx2 + tp%ctx4)**2  )
    quadsum_y(11) = nty1 * ((tp%cty1 + tp%cty2)**2 + (tp%cty1 + tp%cty4)**2 &
                                                   + (tp%cty2 + tp%cty4)**2  )
    quadsum_z(11) = ntz1 * ((tp%ctz1 + tp%ctz2)**2 + (tp%ctz1 + tp%ctz4)**2 &
                                                   + (tp%ctz2 + tp%ctz4)**2  )

    ! Top face, triangle 2: top2-top3 x top4-top3
    call cross_prod( tp%ctx2 - tp%ctx3, tp%cty2 - tp%cty3, tp%ctz2 - tp%ctz3, &
                     tp%ctx4 - tp%ctx3, tp%cty4 - tp%cty3, tp%ctz4 - tp%ctz3, &
                     ntx2, nty2, ntz2 )
    trip_prod(12) = tp%ctx3 *  ntx2 + tp%cty3 *  nty2 + tp%ctz3 *  ntz2
    quadsum_x(12) = ntx2 * ((tp%ctx4 + tp%ctx2)**2 + (tp%ctx4 + tp%ctx3)**2 &
                                                   + (tp%ctx2 + tp%ctx3)**2  )
    quadsum_y(12) = nty2 * ((tp%cty4 + tp%cty2)**2 + (tp%cty4 + tp%cty3)**2 &
                                                   + (tp%cty2 + tp%cty3)**2  )
    quadsum_z(12) = ntz2 * ((tp%ctz4 + tp%ctz2)**2 + (tp%ctz4 + tp%ctz3)**2 &
                                                   + (tp%ctz2 + tp%ctz3)**2  )

    ! For the volume, add the contributions from each cross product, 
    ! dotted with a vector from the origin to one of the vertices
    tp%vol = 0.0 
    do i = 1, 12
        !write(*, fmt='(A, I2, A, ES11.4)') &
        !      'triple_product ', i, ': ', trip_prod(i)
        tp%vol = tp%vol + trip_prod(i)
    end do
    tp%vol = tp%vol / 6.0

    ! Calculate the centroid based on the quadrature sum contributions
    tp%ox = 0.0
    tp%oy = 0.0
    tp%oz = 0.0
    do i = 1, 12
        !write(*, fmt='(A, I2, A, ES11.4)') &
        !      'quadsum ', i, ': ', quadsum_x(i)
        tp%ox = tp%ox + quadsum_x(i)
        tp%oy = tp%oy + quadsum_y(i)
        tp%oz = tp%oz + quadsum_z(i)
    end do
    tp%ox = tp%ox / (tp%vol * 48.0)
    tp%oy = tp%oy / (tp%vol * 48.0)
    tp%oz = tp%oz / (tp%vol * 48.0)

end subroutine trap_volume_ctr

!-------------------------------------------------------------------------------
! check_trapezoid(tp)
! 
! Checks a trapezoid for erroneous or otherwise undesirable properties. If one
! or more such properties are found, the error flag of the input trapezoid is
! set to true.
!
! List of erroneous properties:
!     - The height of the top face above the base face is negative
!     - Any pair of edges of the top face intersect one another
!     - Any pair of edges of the bottom face intersect one another
!-------------------------------------------------------------------------------
subroutine check_trapezoid(tp)

    use magnet_set_globals, only: trapezoid
    use magnet_set_calc,    only: cross_prod

    implicit none

    type(trapezoid) :: tp
    REAL(8) :: height
    REAL(8) :: ut1, ut2, ut3, ut4, vt1, vt2, vt3, vt4
    REAL(8) :: ub1, ub2, ub3, ub4, vb1, vb2, vb3, vb4
    REAL(8) :: unx, uny, unz, vnx, vny, vnz, norm

    tp%err = .false.

    !---------------------------------------------------------------------------
    ! Ensure that the top plane's height along the stack normal is positive
    !---------------------------------------------------------------------------

    height = (tp%otx-tp%obx)*tp%snx + (tp%oty-tp%oby)*tp%sny + &
             (tp%otz-tp%obz)*tp%snz
    if (height < 0.0) tp%err = .true.

    !---------------------------------------------------------------------------
    ! Check for crossed segments after transform to 2D coords in top/base planes
    !---------------------------------------------------------------------------

    ! Orthogonal basis vectors in the base (and also top) planes
    norm = sqrt( (tp%ctx2-tp%ctx1)**2 + (tp%cty2-tp%cty1)**2 + &
                 (tp%ctz2-tp%ctz1)**2 )
    unx = (tp%ctx2-tp%ctx1) / norm
    uny = (tp%cty2-tp%cty1) / norm
    unz = (tp%ctz2-tp%ctz1) / norm
    call cross_prod( tp%snx, tp%sny, tp%snz, unx, uny, unz, vnx, vny, vnz )

    ! Corner points in the 2D orthogonal (u,v) basis
    ut1 = (tp%ctx1-tp%otx)*unx + (tp%cty1-tp%oty)*uny + (tp%ctz1-tp%otz)*unz
    ut2 = (tp%ctx2-tp%otx)*unx + (tp%cty2-tp%oty)*uny + (tp%ctz2-tp%otz)*unz
    ut3 = (tp%ctx3-tp%otx)*unx + (tp%cty3-tp%oty)*uny + (tp%ctz3-tp%otz)*unz
    ut4 = (tp%ctx4-tp%otx)*unx + (tp%cty4-tp%oty)*uny + (tp%ctz4-tp%otz)*unz
    ub1 = (tp%cbx1-tp%obx)*unx + (tp%cby1-tp%oby)*uny + (tp%cbz1-tp%obz)*unz
    ub2 = (tp%cbx2-tp%obx)*unx + (tp%cby2-tp%oby)*uny + (tp%cbz2-tp%obz)*unz
    ub3 = (tp%cbx3-tp%obx)*unx + (tp%cby3-tp%oby)*uny + (tp%cbz3-tp%obz)*unz
    ub4 = (tp%cbx4-tp%obx)*unx + (tp%cby4-tp%oby)*uny + (tp%cbz4-tp%obz)*unz
    vt1 = (tp%ctx1-tp%otx)*vnx + (tp%cty1-tp%oty)*vny + (tp%ctz1-tp%otz)*vnz
    vt2 = (tp%ctx2-tp%otx)*vnx + (tp%cty2-tp%oty)*vny + (tp%ctz2-tp%otz)*vnz
    vt3 = (tp%ctx3-tp%otx)*vnx + (tp%cty3-tp%oty)*vny + (tp%ctz3-tp%otz)*vnz
    vt4 = (tp%ctx4-tp%otx)*vnx + (tp%cty4-tp%oty)*vny + (tp%ctz4-tp%otz)*vnz
    vb1 = (tp%cbx1-tp%obx)*vnx + (tp%cby1-tp%oby)*vny + (tp%cbz1-tp%obz)*vnz
    vb2 = (tp%cbx2-tp%obx)*vnx + (tp%cby2-tp%oby)*vny + (tp%cbz2-tp%obz)*vnz
    vb3 = (tp%cbx3-tp%obx)*vnx + (tp%cby3-tp%oby)*vny + (tp%cbz3-tp%obz)*vnz
    vb4 = (tp%cbx4-tp%obx)*vnx + (tp%cby4-tp%oby)*vny + (tp%cbz4-tp%obz)*vnz

    ! In the orthogonal basis, check if any opposite edge segments cross
    if (segments_crossed(ut1, vt1, ut3, vt3, ut2, vt2, ut4, vt4) .or.    &
        segments_crossed(ut2, vt2, ut4, vt4, ut3, vt3, ut1, vt1) .or.    &
        segments_crossed(ub1, vb1, ub3, vb3, ub2, vb2, ub4, vb4) .or.    &
        segments_crossed(ub2, vb2, ub4, vb4, ub3, vb3, ub1, vb1)       ) &
            tp%err = .true.

end subroutine check_trapezoid

!-------------------------------------------------------------------------------
! segments_crossed(x1s, y1s, x2s, y2s, x1e, y1e, x2e, y2e)
!
! Determines whether two line segments in a 2d plane intersect one another.
!
! Input parameters:
!     REAL(8) x1s, y1s   -> x, y coordinates, start point, segment 1
!     REAL(8) x2s, y2s   -> x, y coordinates, start point, segment 2
!     REAL(8) x1e, y1e   -> x, y coordinates, end point, segment 1
!     REAL(8) x2e, y2e   -> x, y coordinates, end point, segment 2
!-------------------------------------------------------------------------------
logical function segments_crossed(x1s, y1s, x2s, y2s, x1e, y1e, x2e, y2e)

    use magnet_set_calc, only: inverse_2x2

    implicit none

    REAL(8), intent(IN) :: x1s, y1s, x2s, y2s, x1e, y1e, x2e, y2e
    REAL(8) :: ax1, ay1, ax2, ay2, length1, length2
    REAL(8) :: inv11, inv12, inv21, inv22
    REAL(8) :: s1, s2

    ! Lengths of each segment
    length1 = sqrt( (x1e-x1s)**2 + (y1e-y1s)**2 )
    length2 = sqrt( (x2e-x2s)**2 + (y2e-y2s)**2 )

    ! Unit axis vectors for each segment
    ax1 = (x1e-x1s) / length1
    ay1 = (y1e-y1s) / length1
    ax2 = (x2e-x2s) / length2
    ay2 = (y2e-y2s) / length2

    ! Inverse of the matrix defined by [ ax1  -ax2;  ay1  -ay2  ]
    call inverse_2x2(ax1, -ax2, ay1, -ay2, inv11, inv12, inv21, inv22)

    ! Distances of the intersection from start points of each line segment
    s1 = inv11*(x2s-x1s) + inv12*(y2s-y1s)
    s2 = inv21*(x2s-x1s) + inv22*(y2s-y1s)

    ! Determine whether the line intersection point is within both segments
    if (s1 >= 0.0 .and. s1 <= length1 .and. s2 >= 0.0 .and. s2 <= length2) then
        segments_crossed = .true.
    else
        segments_crossed = .false.
    end if

    ! Uncomment for debugging
    !write(*,fmt='(A, F8.4, A, F8.4)') &
    !    'Intersection: ', x1s+s1*ax1, ', ', y1s+s1*ay1

end function segments_crossed

end module magnet_set_trapz

