!-------------------------------------------------------------------------------
! magnet_set_globals.f90
!
! Global variables pertaining to the spatial layout of the set of magnets for
! the SAS device.
!
! Author:       K. C. Hammond
! Contact:      khammond@pppl.gov
! Last updated: 2019-10-18
!-------------------------------------------------------------------------------
module magnet_set_globals

    implicit none

    REAL(8), parameter :: in2m = 0.0254
    REAL(8), parameter :: pi = 3.141592653589793
    REAL(8), parameter :: eps = 1.0e-7
    REAL(8), parameter :: deg2rad = pi/180.

    ! Upper limits on some array bounds
    INTEGER, parameter :: np = 200
    INTEGER, parameter :: nv = 200
    INTEGER, parameter :: ns = nv-1
    INTEGER, parameter :: nw = 10
    INTEGER, parameter :: nl = 10
    INTEGER, parameter :: nh = 10

    LOGICAL :: arrays_zeroed = .false.

    ! General parameters
    INTEGER :: coiltype = 2 ! FOCUS setting (should be 2 to indicate dipole)
    INTEGER :: symm = 2     ! whether FOCUS should copy magnets to symmetric
                            ! locations
    CHARACTER(len=15) :: construction = 'planar'

    ! Vacuum vessel geometric information
    logical :: vessel_loaded = .false.
    CHARACTER(len=500) :: vfile = 'none' ! path to the file with the coeffs
    REAL(8) :: ves_tol = 1.0e-5        ! tolerance for vessel-related calculations
    INTEGER :: maxIter = 10         ! max iterations for vessel calculations
    INTEGER :: nModes               ! total number of Fourier modes
    INTEGER :: nfp = 3              ! number of field periods
    REAL(8) :: ves_r00                 ! R_00 coefficient
    REAL(8) :: ves_r10                 ! R_01 coefficient
    REAL(8), allocatable :: vrc(:)     ! cosine coefficients, vessel surface r
    REAL(8), allocatable :: vzs(:)     ! sine coefficients, vessel surface z
    REAL(8), allocatable :: vrs(:)     ! sine coefficients, vessel surface r
    REAL(8), allocatable :: vzc(:)     ! cosine coefficients, vessel surface z
    INTEGER, allocatable :: vm(:)   ! poloidal mode numbers for the above coeffs
    INTEGER, allocatable :: vn(:)   ! toroidal mode numbers for the above coeffs

    ! Plasma lcfs geometric information
    logical :: plasma_loaded = .false.
    CHARACTER(len=500) :: pfile = 'none' ! path to the file with the coeffs
    INTEGER :: nModesPl             ! total number of Fourier modes
    REAL(8) :: pla_r00                 ! R_00 coefficient
    REAL(8), allocatable :: prc(:)     ! cosine coefficients, lcfs r
    REAL(8), allocatable :: pzs(:)     ! sine coefficients, lcfs z
    REAL(8), allocatable :: prs(:)     ! sine coefficients, lcfs r
    REAL(8), allocatable :: pzc(:)     ! cosine coefficients, lcfs z
    INTEGER, allocatable :: pm(:)   ! poloidal mode numbers for the above coeffs
    INTEGER, allocatable :: pn(:)   ! toroidal mode numbers for the above coeffs

    ! Details about the support plates 
    logical :: plates_initialized = .false.
    logical :: fill_wedge_gaps = .true.  ! if wedge-shaped gaps between segments
                                         ! are to be filled with magnets
    CHARACTER(len=10) :: vertex_mode = 'rz'
    INTEGER :: nPlates                ! number of support plates per half-module
    INTEGER :: nVertices              ! (maximum) number of vertices per plate
    REAL(8), dimension(np) :: plate_phi  ! toroidal angle (central) of each plate
    REAL(8), dimension(np) :: plate_dphi ! toroidal angle subtended by the plate
    INTEGER, dimension(np) :: segs_on_plate ! number of finite-length segmts
    REAL(8), dimension(np,nv) :: vert_r        ! r coordinates of each vertex
    REAL(8), dimension(np,nv) :: vert_z        ! z coordinates of each vertex
    REAL(8), dimension(np,nv) :: vessel_r      ! r coordinate on vessel near vertex
    REAL(8), dimension(np,nv) :: vessel_z      ! z coordinate on vessel near vertex
    REAL(8), dimension(np,nv) :: vert_theta    ! theta coord of each vertex (rad)
    REAL(8), dimension(np,nv) :: vert_sep      ! sep. btwn vertex and vessel (m)
    REAL(8), dimension(np,nv) :: vert_tlcfs    ! lcfs theta coord of vertex (rad)

    ! Details about each plate segment (populated in call to count_magnets)
    logical :: segments_initialized = .false.
    INTEGER :: nSegments
    type segment
        REAL(8) :: lg             ! length of the segment between vertices
        REAL(8) :: wd             ! width (toroidal dimension) of the segment
        REAL(8) :: r1, r2         ! r coordinates of each end
        REAL(8) :: z1, z2         ! z coordinates of each end
        REAL(8) :: phi            ! toroidal angle of the segment (center)
        REAL(8) :: dphi           ! toroidal angle subtended by the segment
        REAL(8) :: alpha          ! angle between plate normal and r-hat
        REAL(8) :: gap_con_1      ! concavity gap, begin. of plate segmt
        REAL(8) :: gap_con_2      ! concavity gap, end of plate segmt
        REAL(8) :: divider_elev_1 ! elev angle of div line above segmt, beginning
        REAL(8) :: divider_elev_2 ! elev angle of div line above segmt, end
        INTEGER :: nStacks     ! num. of mag. stacks on each segmt
        INTEGER :: mag_lg_ind  ! index of magnet length for segmt
    end type segment
    type(segment), allocatable :: segments(:)

    ! Details about each trapezoidal box
    logical :: trapezoids_initialized = .false.
    logical :: subtraps_initialized = .false.
    integer :: nTrapezoids
    integer :: trap_err_count
    integer :: nTraps_per_stack = 1
    integer :: nSubtraps
    type trapezoid
        REAL(8) :: ntbx, ntby, ntbz ! unit normal vector of toroidal back plane
        REAL(8) :: ntfx, ntfy, ntfz ! unit normal vector of toroidal front plane
        REAL(8) :: npbx, npby, npbz ! unit normal vector of poloidal back plane
        REAL(8) :: npfx, npfy, npfz ! unit normal vector of poloidal front plane
        REAL(8) :: snx,  sny,  snz  ! stack normal vector (toward vessel)
        REAL(8) :: otbx, otby, otbz ! reference point on toroidal back plane
        REAL(8) :: otfx, otfy, otfz ! reference point on toroidal front plane
        REAL(8) :: opbx, opby, opbz ! reference point on poloidal back plane
        REAL(8) :: opfx, opfy, opfz ! reference point on poloidal front plane
        REAL(8) :: otx,  oty,  otz  ! reference point on top plane
        REAL(8) :: obx,  oby,  obz  ! reference point ("origin") on bottom plane
        REAL(8) :: vx, vy, vz       ! intersection of stack normal with vessel
        REAL(8) :: ves_theta        ! theta, intersect., stack norm with vessel
        REAL(8) :: ves_phi          ! phi, intersection of norm with vessel
        REAL(8) :: ctx1, cty1, ctz1 ! coords, corner 1, top face (b tor b pol)
        REAL(8) :: ctx2, cty2, ctz2 ! coords, corner 2, top face (f tor b pol)
        REAL(8) :: ctx3, cty3, ctz3 ! coords, corner 3, top face (f tor f pol)
        REAL(8) :: ctx4, cty4, ctz4 ! coords, corner 4, top face (b tor f pol)
        REAL(8) :: cbx1, cby1, cbz1 ! coords, corner 1, bot face (b tor b pol)
        REAL(8) :: cbx2, cby2, cbz2 ! coords, corner 2, bot face (f tor b pol)
        REAL(8) :: cbx3, cby3, cbz3 ! coords, corner 3, bot face (f tor f pol)
        REAL(8) :: cbx4, cby4, cbz4 ! coords, corner 4, bot face (b tor f pol)
        REAL(8) :: vol              ! volume contained within the trapezoid
        REAL(8) :: ox, oy, oz       ! coords, centroid (center of mass)
        logical :: err              ! true if geometry is erroneous/undesirable
    end type trapezoid
    type(trapezoid), allocatable :: traps(:)
    type(trapezoid), allocatable :: subtraps(:)

    ! Specifications for spacing between magnets
    REAL(8) :: gap_rad = 0  ! gap between magnets within a stack (minor radial)
    REAL(8) :: gap_tor = 0  ! gap between magnets on adjacent plates (toroidal)
    REAL(8) :: gap_pol = 0  ! gap between magnets along a plate segment (poloidal)
    REAL(8) :: gap_end = 0  ! gap on each end of a plate edge (not incl. concavity)
    REAL(8) :: gap_seg = 0  ! gap between plate segment and bottom of magnet stack
    REAL(8) :: gap_vac = 0  ! gap between top of magnet stack and vacuum vessel
    REAL(8) :: gap_top = 0  ! clearance for top magnet (additional to gap_vac)
    REAL(8) :: gap_trp = 0  ! gap between adjacent trapezoidal boxes

    ! Allowable parameters, assuming a discrete selection is provided
    REAL(8), dimension(nw) :: avail_mag_wd  ! allowable magnet widths
    REAL(8), dimension(nh) :: avail_mag_ht  ! allowable magnet heights
    REAL(8), dimension(nl) :: avail_mag_lg  ! allowable magnet lengths

    ! Parameters for stacks of magnets (within each position on a segment)
    logical :: stacks_initialized = .false.
    INTEGER :: nStacks_total          ! total number of magnet stacks
    REAL(8)    :: stack_ht_max           ! max height of each stack of magnets
    logical :: fixed_stack_ht = .false. ! true if stack heights are all fixed
    type stack
        REAL(8) :: ht
        REAL(8) :: wd
        REAL(8) :: lg
        REAL(8) :: ox
        REAL(8) :: oy
        REAL(8) :: oz
        REAL(8) :: nx
        REAL(8) :: ny
        REAL(8) :: nz
        REAL(8) :: lx
        REAL(8) :: ly
        REAL(8) :: lz
        INTEGER :: mag_ht_ind
        INTEGER :: nMagnets
    end type stack    
    type(stack), allocatable :: stacks(:)

    ! Individual magnet parameters
    logical :: magnets_initialized = .false.
    CHARACTER(len=20) :: polarization_mode = 'stack'
    INTEGER :: nMagnets_total         ! total number of magnets
    REAL(8)    :: m_bulk                 ! magnetization per unit volume
    INTEGER :: pol_err_count = 0      ! number of erroneous magnet vectors
    type magnet
        REAL(8) :: ht     ! height of the magnet (minor radial)
        REAL(8) :: wd     ! width of the magnet (toroidal)
        REAL(8) :: lg     ! length of the magnet (poloidal)
        REAL(8) :: vol    ! volume of the magnet
        REAL(8) :: ox     ! x coordinate of magnet center
        REAL(8) :: oy     ! y coordinate of magnet center
        REAL(8) :: oz     ! z coordinate of magnet center
        REAL(8) :: nx     ! x component of unit normal vector
        REAL(8) :: ny     ! y component of unit normal vector
        REAL(8) :: nz     ! z component of unit normal vector
        REAL(8) :: lx     ! x comp. poloidal unit vector (for rect. orientation)
        REAL(8) :: ly     ! y component of poloidal unit vector
        REAL(8) :: lz     ! z component of poloidal unit vector
        REAL(8) :: n_phi  ! azimuthal angle of the unit polarization vector
        REAL(8) :: n_theta ! polar angle of the unit polarization vector
        REAL(8) :: moment  ! magnetic moment
        INTEGER :: coiltype ! coil type number for input to FOCUS
        INTEGER :: symm     ! symmetry number for input to FOCUS
        logical :: err
    end type magnet
    type(magnet), allocatable :: magnets(:)

    NAMELIST /magnet_set/ vfile,             &
                          pfile,             &
                          nfp,               &
                          ves_tol,           &
                          maxIter,           &
                          nPlates,           &
                          construction,      &
                          nVertices,         &
                          plate_phi,         &
                          plate_dphi,        &
                          nTraps_per_stack,  &
                          coiltype,          &
                          symm,              &
                          vertex_mode,       &
                          vert_r,            &
                          vert_z,            &
                          vert_theta,        &
                          vert_tlcfs,        &
                          vert_sep,          &
                          gap_rad,           &
                          gap_tor,           &
                          gap_pol,           &
                          gap_end,           &
                          gap_seg,           &
                          gap_vac,           &
                          gap_trp,           &
                          fill_wedge_gaps,   &
                          stack_ht_max,      &
                          fixed_stack_ht,    &
                          avail_mag_wd,      &
                          avail_mag_lg,      &
                          avail_mag_ht,      &
                          m_bulk,            &
                          polarization_mode       
    
end module magnet_set_globals

