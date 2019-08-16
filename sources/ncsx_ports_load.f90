!-------------------------------------------------------------------------------
! ncsx_ports_load.f90
!
! Subroutines for loading the geometric data on the NCSX ports from files.
! In the normal usage case, these should be called only once at runtime.
!-------------------------------------------------------------------------------
module ncsx_ports_load

implicit none

contains

!-------------------------------------------------------------------------------
! load_ports
!
! Loads data on all NCSX ports by calling the other subroutines in this file
!
! Subroutines:
!     load_circular_ports()
!     load_dome_cyl()
!     load_nb_port()
!     load_port_4()
!     load_port_12()
!     populate_circular_port_incl_array()
!-------------------------------------------------------------------------------
subroutine load_ports()

    implicit none

    call load_circular_ports()
    call load_dome_cyl()
    call load_nb_port()
    call load_port_4()
    call load_port_12()
    call populate_circular_port_incl_array()

end subroutine load_ports


!-------------------------------------------------------------------------------
! load_circular_ports
!
! Loads the data on the circular ports from a file and determines the 
! relevant parameters for each port
!
! Modules:
!     ncsx_ports_mod
!-------------------------------------------------------------------------------
subroutine load_circular_ports()

    use ncsx_ports_mod, only: in2m, deg2rad, port_gap, &
                              dir_ncsx_param, file_param_cir, port_unit, &
                              nCirPrt, id_cir, &
                              xo_cir, yo_cir, zo_cir, ax_cir, ay_cir, az_cir, &
                              ir_cir, thick_cir, l0_cir, l1_cir

    implicit none

    !---------------------------------------------------------------------------
    ! Quantities stored in the data files
    !---------------------------------------------------------------------------

    ! y-intercept (inches) of the vertical plane containing the port axis:
    REAL                     :: plane_y0_in
    REAL, dimension(nCirPrt) :: plane_y0 

    ! azimuthal angle (degrees) of the vertical plane (i.e., perpendicular to
    ! the xy plane) in which the axis of the port resides:
    REAL                     :: plane_phi_deg
    REAL, dimension(nCirPrt) :: plane_phi

    ! radial coordinate (inches) of the reference point at which the port axis
    ! intersects the xy plane:
    REAL                     :: r_intersect_in
    REAL, dimension(nCirPrt) :: r_intersect

    ! indicates which root of the algebraic solution for the x, y coordinates
    ! of the intersection of the port axis with the xy plane to use.
    ! (This should be +1 for all circular ports in the first quadrant except
    ! for 17, which should be negative.)
    INTEGER, dimension(nCirPrt) :: root_cir

    ! angle (degrees) of the port axis above the xy plane
    REAL                     :: elev_theta_deg
    REAL, dimension(nCirPrt) :: elev_theta

    ! placeholders in inches/degrees for quantities in the module
    REAL :: ir_cir_in, thick_cir_in, l0_cir_in, l1_cir_in

    ! variables for file i/o
    CHARACTER(len=150) :: filename_cir
    INTEGER :: open_status, read_status
    INTEGER :: i
    REAL    :: ignore
    CHARACTER(len=16) :: errpfx = "load_circular_ports: "

    !---------------------------------------------------------------------------
    ! Loading of parameters from file
    !---------------------------------------------------------------------------

    filename_cir = trim(dir_ncsx_param) // file_param_cir
    open(unit=port_unit, file=filename_cir, status="old", action="read", &
         iostat=open_status)
    if (open_status /= 0) &
        stop trim(errpfx) // " Unable to open file " // filename_cir

    ! Skip over the first row (titles)
    read(unit=port_unit, fmt='(A200)', iostat=read_status) 
    if (read_status /= 0) &
        stop trim(errpfx) // " Error reading title row of " // filename_cir

    ! Read data from each row
    do i = 1, nCirPrt
        read(unit=port_unit, fmt=*, iostat=read_status) &
             id_cir(i), plane_y0_in, plane_phi_deg, r_intersect_in, &
             root_cir(i), elev_theta_deg, ignore, ignore, ignore, &
             ir_cir_in, thick_cir_in, l0_cir_in, l1_cir_in
        if (read_status /= 0) &
            stop trim(errpfx) // " Error reading file " // filename_cir 

        ! Convert to metric units where necessary
        plane_y0(i)    = in2m    * plane_y0_in
        plane_phi(i)   = deg2rad * plane_phi_deg
        r_intersect(i) = in2m    * r_intersect_in
        elev_theta(i)  = deg2rad * elev_theta_deg
        ir_cir(i)      = in2m    * ir_cir_in
        thick_cir(i)   = in2m    * thick_cir_in
        l0_cir(i)      = in2m    * l0_cir_in 
        l1_cir(i)      = in2m    * l1_cir_in

        ! Incorporate port clearance requirement in port thickness
        thick_cir(i) = thick_cir(i) + port_gap
    end do

    close(port_unit)

    !---------------------------------------------------------------------------
    ! Conversion to global parameters
    !---------------------------------------------------------------------------

    do i = 1, nCirPrt

        ! x coordinate of axis reference point
        xo_cir(i) = &
            cos(plane_phi(i)) &
                * ( -plane_y0(i) * sin(plane_phi(i)) &
                   + root_cir(i) * sqrt( (plane_y0(i)*sin(plane_phi(i)))**2  &
                                         + r_intersect(i)**2 - plane_y0(i)**2 ))

        ! y coordinate of axis reference point
        yo_cir(i) = xo_cir(i) * tan(plane_phi(i)) + plane_y0(i)

        ! z coordinate of axis reference point
        zo_cir(i) = 0.

        ! Axis vector
        ax_cir(i) = cos(elev_theta(i)) * cos(plane_phi(i))
        ay_cir(i) = cos(elev_theta(i)) * sin(plane_phi(i))
        az_cir(i) = sin(elev_theta(i))

    end do
    
end subroutine load_circular_ports


!-------------------------------------------------------------------------------
! load_dome_cyl
!
! Loads the data on the dome (modeled as a cylinder, equivalent to a circular
! port) from a file and determines the relevant parameters 
!
! Modules:
!     ncsx_ports_mod
!-------------------------------------------------------------------------------
subroutine load_dome_cyl

    use ncsx_ports_mod, only: in2m, dir_ncsx_param, file_param_dom, port_unit, &
                              port_gap, &
                              xo_dom, yo_dom, zo_dom, ax_dom, ay_dom, az_dom, &
                              ir_dom, thick_dom, l0_dom, l1_dom

    implicit none

    ! Placeholders for quantities in inches
    REAL :: xo_dom_in, yo_dom_in, zo_dom_in
    REAL :: ir_dom_in, thick_dom_in, l0_dom_in, l1_dom_in

    ! Dome axis parameters prior to normalization
    REAL :: ax_dom_pre, ay_dom_pre, az_dom_pre, length_a_dom_pre

    ! File i/o and helper variables
    CHARACTER(len=150) :: filename
    INTEGER :: open_status, read_status
    CHARACTER(len=16) :: errpfx = "load_dome_cyl: "

    ! Load data from file
    filename = trim(dir_ncsx_param) // file_param_dom
    open(unit=port_unit, file=filename, status="old", action="read", &
         iostat=open_status)
    if (open_status /= 0) &
        stop trim(errpfx) // " Unable to open file " // filename

    ! Skip the first line (headers)
    read(unit=port_unit, fmt=*, iostat=read_status)
    if (read_status /= 0) &
        stop trim(errpfx) // " Error reading header of " // filename

    ! Extract data from second line
    read(unit=port_unit, fmt=*, iostat=read_status) &
        xo_dom_in, yo_dom_in, zo_dom_in, ax_dom_pre, ay_dom_pre, az_dom_pre, &
        ir_dom_in, l0_dom_in, l1_dom_in, thick_dom_in
    if (read_status /= 0) &
        stop trim(errpfx) // " Error extracting data from " // filename

    close(port_unit)

    ! Convert to metric units
    xo_dom = in2m * xo_dom_in
    yo_dom = in2m * yo_dom_in
    zo_dom = in2m * zo_dom_in
    ir_dom = in2m * ir_dom_in
    thick_dom = in2m * thick_dom_in
    l0_dom = in2m * l0_dom_in
    l1_dom = in2m * l1_dom_in

    ! Incorporate clearance requirement into thickness
    thick_dom = thick_dom + port_gap

    ! Ensure normalization of the unit vector
    length_a_dom_pre = sqrt(ax_dom_pre**2 + ay_dom_pre**2 + az_dom_pre**2)
    ax_dom = ax_dom_pre / length_a_dom_pre
    ay_dom = ay_dom_pre / length_a_dom_pre
    az_dom = az_dom_pre / length_a_dom_pre

end subroutine load_dome_cyl

!-------------------------------------------------------------------------------
! load_nb_port
!
! Loads the data on the nb port
!
! Modules:
!     ncsx_ports_mod
!-------------------------------------------------------------------------------
subroutine load_nb_port()

    use ncsx_ports_mod, only: in2m, port_gap, &
                              dir_ncsx_param, file_param_nbp, port_unit, &
                              yo1_pnb, zo1_pnb, yo2_pnb, zo2_pnb, &
                              yo3_pnb, zo3_pnb, yo6_pnb, zo6_pnb, &
                              or1sq_pnb, or2sq_pnb, or3sq_pnb, or6sq_pnb, &
                              l0_pnb, l1_pnb, zMax_pnb, yMax_pnb, &
                              yTop_pnb, zTop_pnb, yBot_pnb, zBot_pnb, &
                              yUp_pnb,  zUp_pnb,  yLo_pnb,  zLo_pnb, &
                              slopeUp_pnb, slopeLo_pnb

    implicit none

    ! placeholder variables for quantities in inches
    REAL :: yo1_pnb_in, zo1_pnb_in, yo2_pnb_in, zo2_pnb_in, &
            yo3_pnb_in, zo3_pnb_in, yo6_pnb_in, zo6_pnb_in, &
            ir1_pnb_in, ir2_pnb_in, ir3_pnb_in, ir6_pnb_in, &
            l0_pnb_in, l1_pnb_in, thick_pnb_in

    ! Helper variables
    REAL :: or1, or2, or3, or6 ! outer radii of rounded edges
    REAL :: dy12, dz12, dist12, theta12, dTheta12
    REAL :: dy61, dz61, dist61, theta61, dTheta61

    ! helpers for file i/o
    CHARACTER(len=150) :: filename
    INTEGER :: open_status, read_status
    REAL    :: ignore
    CHARACTER(len=20) :: errpfx = "load_nb_port: "

    ! Open the relevant file
    filename = trim(dir_ncsx_param) // file_param_nbp
    open(unit=port_unit, file=filename, status="old", action="read", &
         iostat=open_status)
    if (open_status /= 0) &
        stop trim(errpfx) // " Unable to open file " // filename

    ! Ignore the header line
    read(unit=port_unit, fmt=*, iostat=read_status)
    if (read_status /= 0) &
        stop trim(errpfx) // " Error reading header line of " // filename

    ! Read data from second line
    read(unit=port_unit, fmt=*, iostat=read_status) &
        yo1_pnb_in, zo1_pnb_in, ir1_pnb_in, &
        yo2_pnb_in, zo2_pnb_in, ir2_pnb_in, &
        yo3_pnb_in, zo3_pnb_in, ir3_pnb_in, &
        ignore, ignore, ignore, ignore, ignore, ignore, &
        yo6_pnb_in, zo6_pnb_in, ir6_pnb_in, &
        ignore, thick_pnb_in, l0_pnb_in, l1_pnb_in
    if (read_status /= 0) &
        stop trim(errpfx) // " Error extracting data from " // filename

    close(port_unit)

    ! Convert quantities to meters
    yo1_pnb = in2m * yo1_pnb_in
    zo1_pnb = in2m * zo1_pnb_in
    yo2_pnb = in2m * yo2_pnb_in
    zo2_pnb = in2m * zo2_pnb_in
    yo3_pnb = in2m * yo3_pnb_in
    zo3_pnb = in2m * zo3_pnb_in
    yo6_pnb = in2m * yo6_pnb_in
    zo6_pnb = in2m * zo6_pnb_in
    l0_pnb  = in2m * l0_pnb_in
    l1_pnb  = in2m * l1_pnb_in

    ! Outer radii (& squares) of rounded edges, incorporating clearance req'mt
    or1 = in2m * (ir1_pnb_in + thick_pnb_in) + port_gap
    or2 = in2m * (ir2_pnb_in + thick_pnb_in) + port_gap
    or3 = in2m * (ir3_pnb_in + thick_pnb_in) + port_gap
    or6 = in2m * (ir6_pnb_in + thick_pnb_in) + port_gap
    or1sq_pnb = or1*or1
    or2sq_pnb = or2*or2
    or3sq_pnb = or3*or3
    or6sq_pnb = or6*or6

    ! Intersection points of upper segment with top (#2) and middle (#1) circles
    dy12 = yo2_pnb - yo1_pnb
    dz12 = zo2_pnb - zo1_pnb
    dist12 = sqrt(dy12*dy12 + dz12*dz12)
    dTheta12 = atan(dy12, dz12)
    call tangent_connector_angles(or1, or2, dist12, theta12, ignore, ignore)
    yTop_pnb = yo2_pnb + or2*sin(dTheta12 + theta12)
    zTop_pnb = zo2_pnb + or2*cos(dTheta12 + theta12)
    yUp_pnb  = yo1_pnb + or1*sin(dTheta12 + theta12)
    zUp_pnb  = zo1_pnb + or1*cos(dTheta12 + theta12)

    ! Intersection point of lower segment with mid (#1) and bot (#6) circles
    dy61 = yo1_pnb - yo6_pnb
    dz61 = zo1_pnb - zo6_pnb
    dist61 = sqrt(dy61*dy61 + dz61*dz61)
    dTheta61 = atan(dy61, dz61)
    call tangent_connector_angles(or6, or1, dist61, theta61, ignore, ignore)
    yBot_pnb = yo6_pnb + or6*sin(dTheta61 + theta61)
    zBot_pnb = zo6_pnb + or6*cos(dTheta61 + theta61)
    yLo_pnb  = yo1_pnb + or1*sin(dTheta61 + theta61)
    zLo_pnb  = zo1_pnb + or1*cos(dTheta61 + theta61)

    ! Slopes of the upper and lower segments
    slopeUp_pnb = ( yUp_pnb - yTop_pnb ) / ( zUp_pnb - zTop_pnb )
    slopeLo_pnb = ( yBot_pnb - yLo_pnb ) / ( zBot_pnb - zLo_pnb )

    ! Maximal extents of y and z
    yMax_pnb = yo1_pnb + or1
    zMax_pnb = zo2_pnb + or2

end subroutine load_nb_port


!-------------------------------------------------------------------------------
! load_port_4
!
! Loads the data on Port 4 and converts to global parameters
!
! Modules:
!     ncsx_ports_mod
!-------------------------------------------------------------------------------
subroutine load_port_4()

    use ncsx_ports_mod, only: in2m, deg2rad, port_gap, &
                              dir_ncsx_param, file_param_p04, port_unit, &
                              phi_p04, theta_a_p04, theta_b_p04, &
                              w1a_p04, w1b_p04, w2a_p04, w2b_p04, &
                              slope_a_p04, slope_b_p04, h_p04, zMax_p04, &
                              z_circ_p04, r_circ_p04, orc_p04, orc2_p04, &
                              l0_p04, l1_p04, l_narrow_p04, thick_p04, &
                              thick_ai_p04, thick_bi_p04, &
                              thick_ao_p04, thick_bo_p04, &
                              rhat_x_p04, rhat_y_p04, phat_x_p04, phat_y_p04

    implicit none

    ! Placeholder variables for quantities in inches
    REAL :: w1a_p04_in, w1b_p04_in, w2a_p04_in, w2b_p04_in, h_p04_in, &
            r_circ_p04_in, l0_p04_in, l1_p04_in, l_narrow_p04_in, thick_p04_in

    ! Placeholder variables for quantities in degrees
    REAL :: phi_p04_deg, theta_a_p04_deg, theta_b_p04_deg

    ! Additional data values in file
    REAL :: d_in   ! Distance between outboard end of port and narrow point 
    REAL :: d0_in  ! Distance between narrow point and inboard end of port

    ! Helpers for file i/o
    CHARACTER(len=150) :: filename
    INTEGER :: open_status, read_status
    CHARACTER(len=25) :: errpfx = "load_port_4: "

    ! Open the file
    filename = trim(dir_ncsx_param) // file_param_p04
    open(unit=port_unit, file=filename, status="old", action="read", &
         iostat=open_status)
    if (open_status /= 0) &
        stop trim(errpfx) // " Unable to open file " // filename

    ! Skip the header line
    read(unit=port_unit, fmt=*, iostat=read_status)
    if (read_status /= 0) &
        stop trim(errpfx) // " Error reading header line of " // filename

    ! Extract the data
    read(unit=port_unit, fmt=*, iostat=read_status) &
        phi_p04_deg, theta_a_p04_deg, theta_b_p04_deg, &
        w1a_p04_in, w1b_p04_in, w2a_p04_in, w2b_p04_in, h_p04_in, &
        r_circ_p04_in, l1_p04_in, d_in, d0_in, thick_p04_in
    if (read_status /= 0) &
        stop trim(errpfx) // " Error loading data from " // filename

    close(port_unit)

    l0_p04_in = l1_p04_in - (d_in + d0_in)
    l_narrow_p04_in = l1_p04_in - d_in

    ! Quantities in meters
    phi_p04      = deg2rad * phi_p04_deg
    theta_a_p04  = deg2rad * theta_a_p04_deg
    theta_b_p04  = deg2rad * theta_b_p04_deg
    w1a_p04      = in2m    * w1a_p04_in
    w1b_p04      = in2m    * w1b_p04_in
    w2a_p04      = in2m    * w2a_p04_in
    w2b_p04      = in2m    * w2b_p04_in
    h_p04        = in2m    * h_p04_in
    r_circ_p04   = in2m    * r_circ_p04_in
    l0_p04       = in2m    * l0_p04_in
    l1_p04       = in2m    * l1_p04_in
    l_narrow_p04 = in2m    * l_narrow_p04_in
    thick_p04    = in2m    * thick_p04_in

    ! Port clearance requirement
    thick_p04 = thick_p04 + port_gap

    ! Derived quantities
    slope_a_p04  = (w2a_p04 - w1a_p04) / (l1_p04 - l_narrow_p04)
    slope_b_p04  = (w2b_p04 - w1b_p04) / (l1_p04 - l_narrow_p04)
    orc_p04      = r_circ_p04 + thick_p04
    orc2_p04     = orc_p04**2
    zMax_p04     = 0.5*h_p04 + thick_p04
    z_circ_p04   = zMax_p04 - orc_p04

    ! Thickness in the direction perpendicular to the port axis
    thick_ai_p04 = thick_p04 / cos(theta_a_p04)
    thick_bi_p04 = thick_p04 / cos(theta_b_p04)
    thick_ao_p04 = thick_p04 / cos(atan(slope_a_p04))
    thick_bo_p04 = thick_p04 / cos(atan(slope_b_p04))

    ! Unit vectors parallel (r) and perpendicular (p) to the port axis
    rhat_x_p04   = cos(phi_p04)
    rhat_y_p04   = sin(phi_p04)
    phat_x_p04   = -sin(phi_p04)
    phat_y_p04   = cos(phi_p04)

end subroutine load_port_4

!-------------------------------------------------------------------------------
! load_port_12
!
! Loads the data on port 12 and converts to global parameters
!
! Modules: 
!     ncsx_ports_mod
!-------------------------------------------------------------------------------
subroutine load_port_12

    use ncsx_ports_mod, only: in2m, port_gap, &
                              dir_ncsx_param, file_param_p12, port_unit, &
                              xMin_p12, xMax_p12, yMax_p12, xo1_p12, xo2_p12, &
                              xInb_p12, yInb_p12, xOut_p12, yOut_p12, &
                              or1_p12, or2_p12, or1sq_p12, or2sq_p12, &
                              l0_p12, l1_p12, slope_p12

    implicit none

    ! Local variables for storing quantities in inches
    REAL :: xo1_p12_in, xo2_p12_in, ir1_p12_in, ir2_p12_in, &
            l0_p12_in, l1_p12_in, thick_p12_in

    ! Helper variables
    REAL :: theta, ignore

    ! Helpers for file i/o
    CHARACTER(len=150) filename
    INTEGER :: open_status, read_status
    CHARACTER(len=25) :: errpfx = "load_port_12: "

    ! Open the source file
    filename = trim(dir_ncsx_param) // file_param_p12
    open(unit=port_unit, file=filename, status="old", action="read", &
         iostat=open_status)
    if (open_status /= 0) &
        stop trim(errpfx) // " Unable to open file " // filename

    ! Ignore the header row
    read(unit=port_unit, fmt=*, iostat=read_status)
    if (read_status /= 0) &
        stop trim(errpfx) // " Error reading header of " // filename

    ! Read data from file
    read(unit=port_unit, fmt=*, iostat=read_status) &
         xo1_p12_in, ir1_p12_in, xo2_p12_in, ir2_p12_in, &
         l0_p12_in, l1_p12_in, thick_p12_in
    if (read_status /= 0) &
        stop trim(errpfx) // " Error importing data from " // filename

    close(port_unit)

    ! Convert inches to metric
    xo1_p12   = in2m * xo1_p12_in
    xo2_p12   = in2m * xo2_p12_in
    or1_p12   = in2m * (ir1_p12_in + thick_p12_in)
    or2_p12   = in2m * (ir2_p12_in + thick_p12_in)
    l0_p12    = in2m * l0_p12_in
    l1_p12    = in2m * l1_p12_in

    ! Incorporate port clearance requirement
    or1_p12 = or1_p12 + port_gap
    or2_p12 = or2_p12 + port_gap

    ! Derived quantities
    xMin_p12  = xo1_p12 - or1_p12
    xMax_p12  = xo2_p12 + or2_p12
    yMax_p12  = max(or1_p12, or2_p12)
    or1sq_p12 = or1_p12*or1_p12
    or2sq_p12 = or2_p12*or2_p12
    call tangent_connector_angles(or1_p12, or2_p12, (xo2_p12 - xo1_p12), &
                                  theta, ignore, ignore)
    xInb_p12 = xo1_p12 + or1_p12 * cos(theta)
    yInb_p12 = or1_p12 * sin(theta)
    xOut_p12 = xo2_p12 + or2_p12 * cos(theta)
    yOut_p12 = or2_p12 * sin(theta)
    slope_p12 = (yOut_p12 - yInb_p12) / (xOut_p12 - xInb_p12)

end subroutine load_port_12

!-------------------------------------------------------------------------------
! populate_circular_port_incl_array
!
! Populates an array of logicals in which each element dictates whether a
! given circular port should be included, based on previously-defined logical
! variables applying to each individual circular port.
!-------------------------------------------------------------------------------
subroutine populate_circular_port_incl_array()

    use ncsx_ports_mod, only: incl_port_02, incl_port_05, incl_port_06, &
                              incl_port_07, incl_port_08, incl_port_09, &
                              incl_port_10, incl_port_11, incl_port_15, &
                              incl_port_17, incl_port_18, &
                              incl_circular_ports

    implicit none

    incl_circular_ports = &
        (/ incl_port_02, incl_port_05, incl_port_06, incl_port_07, &
           incl_port_08, incl_port_09, incl_port_10, incl_port_11, &
           incl_port_15, incl_port_17, incl_port_18                   /)

end subroutine populate_circular_port_incl_array

!-------------------------------------------------------------------------------
! tangent_connector_angles(r1, r2, d, thetaO, thetaI1, thetaI2)
!
! Calculates the angular locations of the intersection points of a line 
! that is tangent to two circles. The unit vector in the direction from the 
! center of circle 1 to the center of circle 2 is defined to have an angle of 
! zero for both circles.
!
! Input parameters:
!     REAL :: r1, r2      -> radii of circles 1 and 2, respectively
!     REAL :: d           -> distance between circles 1 and 2
!
! Output parameters:
!     REAL :: thetaO      -> angular coordinate (same for both circles) of the 
!                            tangent intersection point for an "outer" tangent
!                            line (i.e., does not pass between the circle 
!                            centers). Note that -thetaO is also a solution.
!     REAL :: thetaI1,    -> angular coordinates of the tangent intersection
!             thetaI2        points for an "inner" tangent line (i.e., passes
!                            between the circle centers). Note that 
!                            (-thetaI1, -thetaI2) is also a solution.
!-------------------------------------------------------------------------------
subroutine tangent_connector_angles(r1, r2, d, thetaO, thetaI1, thetaI2)

    implicit none

    REAL, intent(IN)  :: r1, r2, d
    REAL, intent(OUT) :: thetaO, thetaI1, thetaI2
    REAL              :: cosThetaO, cosThetaI1, cosThetaI2

    cosThetaO  =  (r1 - r2) / d
    cosThetaI1 =  (r1 + r2) / d
    cosThetaI2 = -cosThetaI1

    thetaO  = acos(cosThetaO)
    thetaI1 = acos(cosThetaI1)
    thetaI2 = acos(cosThetaI2)

    return

end subroutine tangent_connector_angles

end module ncsx_ports_load
