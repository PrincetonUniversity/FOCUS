!-------------------------------------------------------------------------------
! magnet_set_globals.f90
!
! Global variables pertaining to the spatial layout of the set of magnets for
! the SAS device.
!-------------------------------------------------------------------------------
module magnet_set_globals

    implicit none

    REAL, parameter :: in2m = 0.0254
    REAL, parameter :: pi = 3.141592653589793
    REAL, parameter :: eps = 1.0e-7
    REAL, parameter :: deg2rad = pi/180.

    ! Upper limits on some array bounds
    INTEGER, parameter :: np = 50
    INTEGER, parameter :: nv = 50
    INTEGER, parameter :: ns = nv-1
    INTEGER, parameter :: nw = 10
    INTEGER, parameter :: nl = 10
    INTEGER, parameter :: nh = 10

    LOGICAL :: arrays_zeroed = .false.

    ! General parameters
    INTEGER :: coiltype = 2 ! FOCUS setting (should be 2 to indicate dipole)
    INTEGER :: symm = 2     ! whether FOCUS should copy magnets to symmetric
                            ! locations

    ! Vacuum vessel geometric information
    logical :: vessel_loaded = .false.
    CHARACTER(len=500) :: vfile = 'none' ! path to the file with the coeffs
    REAL :: ves_tol = 1.0e-5        ! tolerance for vessel-related calculations
    INTEGER :: maxIter = 10         ! max iterations for vessel calculations
    INTEGER :: nModes               ! total number of Fourier modes
    INTEGER :: nfp = 3              ! number of field periods
    REAL :: ves_r00                 ! R_00 coefficient
    REAL, allocatable :: vrc(:)     ! cosine coefficients, vessel surface r
    REAL, allocatable :: vzs(:)     ! sine coefficients, vessel surface z
    REAL, allocatable :: vrs(:)     ! sine coefficients, vessel surface r
    REAL, allocatable :: vzc(:)     ! cosine coefficients, vessel surface z
    INTEGER, allocatable :: vm(:)   ! poloidal mode numbers for the above coeffs
    INTEGER, allocatable :: vn(:)   ! toroidal mode numbers for the above coeffs

    ! Details about the support plates 
    logical :: plates_initialized = .false.
    CHARACTER(len=10) :: vertex_mode = 'rz'
    INTEGER :: nPlates                ! number of support plates per half-module
    INTEGER :: nVertices              ! (maximum) number of vertices per plate
    REAL, dimension(np) :: plate_phi  ! toroidal angle (central) of each plate
    REAL, dimension(np) :: plate_dphi ! toroidal angle subtended by the plate
    INTEGER, dimension(np) :: segs_on_plate ! number of finite-length segmts
    REAL, dimension(np,nv) :: vert_r        ! r coordinates of each vertex
    REAL, dimension(np,nv) :: vert_z        ! z coordinates of each vertex
    REAL, dimension(np,nv) :: vessel_r      ! r coordinate on vessel near vertex
    REAL, dimension(np,nv) :: vessel_z      ! z coordinate on vessel near vertex
    REAL, dimension(np,nv) :: vert_theta    ! theta coord of each vertex (rad)
    REAL, dimension(np,nv) :: vert_sep      ! sep. btwn vertex and vessel (m)

    ! Details about each plate segment (populated in call to count_magnets)
    logical :: segments_initialized = .false.
    INTEGER :: nSegments
    type segment
        REAL :: lg            ! length of the segment between vertices
        REAL :: wd            ! width (toroidal dimension) of the segment
        REAL :: r1, r2        ! r coordinates of each end
        REAL :: z1, z2        ! z coordinates of each end
        REAL :: phi           ! toroidal angle of the segment (center)
        REAL :: dphi          ! toroidal angle subtended by the segment
        REAL :: alpha         ! angle between plate normal and r-hat
        REAL :: gap_con_1     ! concavity gap, begin. of plate segmt
        REAL :: gap_con_2     ! concavity gap, end of plate segmt
        INTEGER :: nStacks    ! num. of mag. stacks on each segmt
        INTEGER :: mag_lg_ind ! index of magnet length for segmt
    end type segment
    type(segment), allocatable :: segments(:)

    ! Specifications for spacing between magnets
    REAL :: gap_rad   ! gap between magnets within a stack (minor radial)
    REAL :: gap_tor   ! gap between magnets on adjacent plates (toroidal)
    REAL :: gap_pol   ! gap between magnets along a plate segment (poloidal)
    REAL :: gap_end   ! gap on each end of a plate edge (not incl. concavity)
    REAL :: gap_seg   ! gap between plate segment and bottom of magnet stack
    REAL :: gap_vac   ! gap between top of magnet stack and vacuum vessel

    ! Allowable parameters, assuming a discrete selection is provided
    REAL, dimension(nw) :: avail_mag_wd  ! allowable magnet widths
    REAL, dimension(nh) :: avail_mag_ht  ! allowable magnet heights
    REAL, dimension(nl) :: avail_mag_lg  ! allowable magnet lengths

    ! Parameters for stacks of magnets (within each position on a segment)
    logical :: stacks_initialized = .false.
    INTEGER :: nStacks_total          ! total number of magnet stacks
    REAL    :: stack_ht_max           ! max height of each stack of magnets
    logical :: fixed_stack_ht = .false. ! true if stack heights are all fixed
    type stack
        REAL :: ht
        REAL :: wd
        REAL :: lg
        REAL :: ox
        REAL :: oy
        REAL :: oz
        REAL :: nx
        REAL :: ny
        REAL :: nz
        REAL :: lx
        REAL :: ly
        REAL :: lz
        INTEGER :: mag_ht_ind
        INTEGER :: nMagnets
    end type stack    
    type(stack), allocatable :: stacks(:)

    ! Individual magnet parameters
    logical :: magnets_initialized = .false.
    INTEGER :: nMagnets_total         ! total number of magnets
    REAL    :: m_bulk                 ! magnetization per unit volume
    type magnet
        REAL :: ht     ! height of the magnet (minor radial)
        REAL :: wd     ! width of the magnet (toroidal)
        REAL :: lg     ! length of the magnet (poloidal)
        REAL :: vol    ! volume of the magnet
        REAL :: ox     ! x coordinate of magnet center
        REAL :: oy     ! y coordinate of magnet center
        REAL :: oz     ! z coordinate of magnet center
        REAL :: nx     ! x component of unit normal vector
        REAL :: ny     ! y component of unit normal vector
        REAL :: nz     ! z component of unit normal vector
        REAL :: lx     ! x comp. poloidal unit vector (for rect. orientation)
        REAL :: ly     ! y component of poloidal unit vector
        REAL :: lz     ! z component of poloidal unit vector
        REAL :: n_phi  ! azimuthal angle of normal vector
        REAL :: n_theta ! polar angle of normal vector
        REAL :: moment  ! magnetic moment
        INTEGER :: coiltype ! coil type number for input to FOCUS
        INTEGER :: symm     ! symmetry number for input to FOCUS
    end type magnet
    type(magnet), allocatable :: magnets(:)

    NAMELIST /magnet_set/ vfile,          &
                          nfp,            &
                          ves_tol,        &
                          maxIter,        &
                          nPlates,        &
                          nVertices,      &
                          plate_phi,      &
                          plate_dphi,     &
                          coiltype,       &
                          symm,           &
                          vertex_mode,    &
                          vert_r,         &
                          vert_z,         &
                          vert_theta,     &
                          vert_sep,       &
                          gap_rad,        &
                          gap_tor,        &
                          gap_pol,        &
                          gap_end,        &
                          gap_seg,        &
                          gap_vac,        &
                          stack_ht_max,   &
                          fixed_stack_ht, &
                          avail_mag_wd,   &
                          avail_mag_lg,   &
                          avail_mag_ht,   &
                          m_bulk       
    
end module magnet_set_globals

