!-------------------------------------------------------------------------------
! ncsx_ports_mod.f90
! 
! Global variables used by the subroutines that determine whether points
! are inside NCSX ports
!-------------------------------------------------------------------------------
module ncsx_ports_mod

    use focus_globals, only: dp

    implicit none

    REAL, parameter :: in2m = 0.0254  ! number of meters in one inch
    REAL, parameter :: pi = 3.141592653589793
    REAL, parameter :: deg2rad = pi/180. ! degrees in one radian

    ! If false, skip all ncsx_ports subroutines (e.g., if no namelist found)
    LOGICAL :: ncsx_ports_on = .false.

    ! Checks if the values in this module have been initialized
    LOGICAL :: ncsx_ports_mod_initialized = .false.

    ! Required clearance between magnet locations and port walls (meters)
    REAL :: port_gap = 0.

    ! Number of field periods in NCSX
    INTEGER :: nfp_ncsx = 3

    ! Number of ports included in the model
    INTEGER :: nPortsNCSX = 15

    ! Locations of files with port parameter data
    CHARACTER(len=100) :: dir_ncsx_param = "ncsx_ports_geometry/"
    CHARACTER(len=100) :: file_param_cir = "Circular_Port_Parameters.csv"
    CHARACTER(len=100) :: file_param_nbp = "NB_Port_Parameters.csv"
    CHARACTER(len=100) :: file_param_p12 = "Port_12_Parameters.csv"
    CHARACTER(len=100) :: file_param_p04 = "Port_4_Parameters.csv"
    CHARACTER(len=100) :: file_param_dom = "Dome_Cylindrical_Parameters.csv"
    INTEGER, parameter :: port_unit = 18   ! unit for accessing port data files

    ! Circular port parameters
    INTEGER, parameter          :: nCirPrt = 11        ! number of circ. ports
    INTEGER, dimension(nCirPrt) :: id_cir              ! ID number of port
    REAL, dimension(nCirPrt) :: xo_cir, yo_cir, zo_cir ! origin of port axis
    REAL, dimension(nCirPrt) :: ax_cir, ay_cir, az_cir ! axis dir. from origin
    REAL, dimension(nCirPrt) :: ir_cir                 ! inner radius
    REAL, dimension(nCirPrt) :: thick_cir       ! wall thickness incl. clearance
    REAL, dimension(nCirPrt) :: l0_cir, l1_cir  ! port endpts. along axis

    ! Dome parameters (presently modeled as an additional circular port)
    REAL :: xo_dom, yo_dom, zo_dom   ! origin of dome axis
    REAL :: ax_dom, ay_dom, az_dom   ! axis dir. from origin
    REAL :: ir_dom                   ! inner radius
    REAL :: thick_dom                ! wall thickness
    REAL :: l0_dom, l1_dom           ! endpoints of dome along axis

    ! NB port parameters
    ! Note: there are 6 "ref. circles" @ rounded edges, but only 1-3, 6 are used
    REAL :: yo1_pnb, zo1_pnb  ! y, z coords. of 1st ref. circle (at midplane)
    REAL :: yo2_pnb, zo2_pnb  ! y, z coords. of 2nd ref. circle (above midplane)
    REAL :: yo3_pnb, zo3_pnb  ! y, z coords. of 3rd ref. circle (above midplane)
    REAL :: yo6_pnb, zo6_pnb  ! y, z coords. of 6th ref. circle (below midplane)
    REAL :: or1sq_pnb          ! square of outer radius of circle 1 incl. clear.
    REAL :: or2sq_pnb          ! square of outer radius of circle 2 incl. clear.
    REAL :: or3sq_pnb          ! square of outer radius of circle 3 incl. clear.
    REAL :: or6sq_pnb          ! square of outer radius of circle 6 incl. clear.
    REAL :: l0_pnb, l1_pnb    ! distances of inner & outer endpoints from z-axis
    REAL :: zMax_pnb           ! maximum extent in z of the port
    REAL :: yMax_pnb           ! maximum extent in y of the port
    REAL :: yTop_pnb, zTop_pnb ! intersect btwn top circle and upper line seg
    REAL ::  yUp_pnb,  zUp_pnb ! intersect btwn mid circle and upper line seg
    REAL ::  yLo_pnb,  zLo_pnb ! intersect btwn mid circle and lower line seg
    REAL :: yBot_pnb, zBot_pnb ! intersect btwn bot circle and lower line seg
    REAL :: slopeUp_pnb        ! slope, line seg connecting top to mid circle
    REAL :: slopeLo_pnb        ! slope, line seg connecting mid to bot circle

    ! Port 4 parameters 
    ! Note: "a" refers to the CW-facing side; "b" refers to the CCW-facing side
    REAL :: phi_p04        ! azimuthal angle of the port's reference axis
    REAL :: theta_a_p04    ! angle btwn. a-side inner wall and ref. axis
    REAL :: theta_b_p04    ! neg. angle btwn. b-side inner wall and ref. axis
    REAL :: w1a_p04        ! distance btwn. axis and a-side wall @ narrow xsect 
    REAL :: w1b_p04        ! distance btwn. axis and b-side wall @ narrow xsect
    REAL :: w2a_p04        ! distance btwn. axis and a-side wall @ port entrance
    REAL :: w2b_p04        ! distance btwn. axis and b-side wall @ port entrance
    REAL :: slope_a_p04    ! slope of outer wall segment on a-side
    REAL :: slope_b_p04    ! slope of outer wall segment on b-side
    REAL :: h_p04          ! height (z-dimension) of port
    REAL :: zMax_p04       ! maximum value of z attained (0.5*h + thickness)
    REAL :: z_circ_p04     ! height of centers of circ. corners above midplane
    REAL :: r_circ_p04     ! radius of rounded corners (in plane perp. to axis)
    REAL :: orc_p04        ! outer radius of rounded corners
    REAL :: orc2_p04       ! outer radius^2, rounded corners 
    REAL :: l0_p04, l1_p04 ! distances of inner and outer endpoints from z-axis
    REAL :: l_narrow_p04   ! distance of narrow cross-section from z-axis
    REAL :: thick_p04      ! wall thickness including clearance
    REAL :: thick_ai_p04   ! wall thickness perp. to port axis, a side, inner
    REAL :: thick_ao_p04   ! wall thickness perp. to port axis, a side, outer
    REAL :: thick_bi_p04   ! wall thickness perp. to port axis, b side, inner
    REAL :: thick_bo_p04   ! wall thickness perp. to port axis, b side, outer
    REAL :: rhat_x_p04     ! unit vector along the port axis, x component
    REAL :: rhat_y_p04     ! unit vector along the port axis, y component
    REAL :: phat_x_p04     ! unit vector perpendicular to port axis, x component
    REAL :: phat_y_p04     ! unit vector perpendicular to port axis, y component

    ! Port 12 parameters
    REAL :: xMin_p12       ! Minimum possible value of the x coordinate
    REAL :: xMax_p12       ! Maximum possible value of the x coordinate
    REAL :: yMax_p12       ! Maximum possible value of the y coordinate
    REAL :: xo1_p12        ! x coord. of the center of the inboard circle
    REAL :: xo2_p12        ! x coord. of the center of the outboard circle
    REAL :: xInb_p12       ! x coord., intersect. of line seg w/ inboard circ.
    REAL :: yInb_p12       ! y coord., intersect. of line seg w/ inboard circ.
    REAL :: xOut_p12       ! x coord., intersect. of line seg w/ outboard circ.
    REAL :: yOut_p12       ! y coord., intersect. of line seg w/ outboard circ.
    REAL :: or1_p12        ! outer radius of the inboard circle incl. clearance
    REAL :: or2_p12        ! outer radius of the outboard circle incl. clearance
    REAL :: or1sq_p12      ! or1_p12**2
    REAL :: or2sq_p12      ! or2_p12**2
    REAL :: l0_p12, l1_p12 ! port endpoints (as distances from the xy plane)
    REAL :: slope_p12      ! slope of line segment

    ! Ports to be included (none by default)
    LOGICAL :: incl_port_nb = .false.
    LOGICAL :: incl_port_02 = .false.
    LOGICAL :: incl_port_04 = .false.
    LOGICAL :: incl_port_05 = .false.
    LOGICAL :: incl_port_06 = .false.
    LOGICAL :: incl_port_07 = .false.
    LOGICAL :: incl_port_08 = .false.
    LOGICAL :: incl_port_09 = .false.
    LOGICAL :: incl_port_10 = .false.
    LOGICAL :: incl_port_11 = .false.
    LOGICAL :: incl_port_12 = .false.
    LOGICAL :: incl_port_15 = .false.
    LOGICAL :: incl_port_17 = .false.
    LOGICAL :: incl_port_18 = .false.
    LOGICAL :: incl_dome    = .false.   ! Note: in practice, exclusion of 
                                          ! the dome would imply the exclusion 
                                          ! of ports 17 and 18 as well.
    LOGICAL, dimension(nCirPrt) :: incl_circular_ports

    NAMELIST /ncsx_ports/ dir_ncsx_param, &
                          file_param_cir, &
                          file_param_nbp, &
                          file_param_p04, &
                          file_param_p12, &
                          file_param_dom, &
                          port_gap,       &
                          incl_port_nb,   &
                          incl_port_02,   &
                          incl_port_04,   &
                          incl_port_05,   &
                          incl_port_06,   &
                          incl_port_07,   &
                          incl_port_08,   &
                          incl_port_09,   &
                          incl_port_10,   &
                          incl_port_11,   &
                          incl_port_12,   &
                          incl_port_15,   &
                          incl_port_17,   &
                          incl_port_18,   &
                          incl_dome


end module ncsx_ports_mod

