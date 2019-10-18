!-------------------------------------------------------------------------------
! magnet_set_calc.f90
!
! Contains subroutines for calculations relating the magnet structure to the
! vessel geometry
!
! Author:       K. C. Hammond
! Contact:      khammond@pppl.gov
! Last updated: 2019-10-18
!-------------------------------------------------------------------------------
module magnet_set_calc

implicit none

contains

!-------------------------------------------------------------------------------
! ves_dist_3d(ox, oy, oz, ax, ay, az, l0, theta0, phi0, l, theta, phi, x, y, z)
!
! Computes the distance between a point and the vessel along a particular
! line in 3D space, using the Newton-Raphson method.
!
! Input parameters:
!     REAL :: ox, oy, oz  -> x, y, and z coordinates of the reference point
!     REAL :: ax, ay, az  -> x, y, and z components of a unit vector defining
!                            the direction from the test point to the vessel
!     REAL :: l0, theta0, 
!             phi0        -> initial guesses of the distance l, as well as
!                            the theta and phi coordinates of the intersection
!                            point on the vessel
!
! Return parameters:
!     REAL :: l           -> distance from the point to the vessel
!     REAL :: theta, phi  -> theta and phi coordinates of the intersection point
!     REAL :: x, y, z     -> x, y, and z coordinates of the intersection point
!-------------------------------------------------------------------------------
subroutine ves_dist_3d(ox, oy, oz, ax, ay, az, l0, theta0, phi0, &
                       l, theta, phi, x, y, z, chi2)

    use magnet_set_globals, only: ves_tol, maxIter

    implicit none

    REAL, intent(IN)  :: ox, oy, oz, ax, ay, az, l0, theta0, phi0
    REAL, intent(OUT) :: l, theta, phi, x, y, z, chi2
    INTEGER :: i
    REAL :: vx, vy, vz, vr, fx, fy, fz, sl, st, sp
    REAL :: drdt, drdp, dzdt, dzdp
    REAL ::  jac_xl,  jac_xt,  jac_xp, &
             jac_yl,  jac_yt,  jac_yp, &
             jac_zl,  jac_zt,  jac_zp
    REAL :: ijac_lx, ijac_ly, ijac_lz, &
            ijac_px, ijac_py, ijac_pz, &
            ijac_tx, ijac_ty, ijac_tz

    ! constant components of the Jacobian (first column)
    jac_xl = -ax
    jac_yl = -ay
    jac_zl = -az

    ! Initialize values for which to solve
    l     = l0
    theta = theta0
    phi   = phi0

    do i = 1, maxIter

        vr = ves_r(theta, phi)
        vx = vr * cos(phi)
        vy = vr * sin(phi)
        vz = ves_z(theta, phi)

        ! components that should ideally be zero
        fx = vx - (ox + l*ax)
        fy = vy - (oy + l*ay)
        fz = vz - (oz + l*az)

        ! terminate if solution has sufficiently converged
        chi2 = fx*fx + fy*fy + fz*fz
        if (chi2 < ves_tol*ves_tol) exit
        
        ! compute the Jacobian, second column
        drdt = ves_drdt(theta, phi)
        dzdt = ves_dzdt(theta, phi)
        jac_xt = drdt * cos(phi)
        jac_yt = drdt * sin(phi)
        jac_zt = dzdt

        ! compute the Jacobian, third column
        drdp = ves_drdp(theta, phi)
        dzdp = ves_dzdp(theta, phi)
        jac_xp = drdp * cos(phi) - vr * sin(phi)
        jac_yp = drdp * sin(phi) + vr * cos(phi)
        jac_zp = dzdp

        ! invert the Jacobian
        call inverse_3x3(jac_xl, jac_xt, jac_xp, jac_yl, jac_yt, jac_yp, &
                         jac_zl, jac_zt, jac_zp, ijac_lx, ijac_ly, ijac_lz, &
                         ijac_tx, ijac_ty, ijac_tz, ijac_px, ijac_py, ijac_pz)

        ! The Newton step
        sl = -ijac_lx * fx - ijac_ly * fy - ijac_lz * fz
        st = -ijac_tx * fx - ijac_ty * fy - ijac_tz * fz
        sp = -ijac_px * fx - ijac_py * fy - ijac_pz * fz

        ! Update the parameters
        l     = l     + sl
        theta = theta + st
        phi   = phi   + sp

    end do

end subroutine ves_dist_3d

!-------------------------------------------------------------------------------
! ves_dist_2d(phi, or, oz, ar, az, l0, theta0, l, theta, r, z, chi2)
!
! Computes the distance between a point and the vessel along a particular
! line in a 2D (poloidal, constant-phi) cross-section
!
! Input parameters:
!     REAL :: phi         -> the toroidal angle (radians) of the cross-section
!     REAL :: or, oz      -> r and z coordinates (meters) of the ref. point
!     REAL :: ar, az      -> r and z components of a unit vector defining
!                            the direction from the test point to the vessel
!     REAL :: l0, theta0  -> initial guesses of the distance l, as well as
!                            the theta coordinate of the intersection point
!                            on the vessel
!
! Return parameters:
!     REAL :: l           -> distance from the point to the vessel
!     REAL :: theta       -> theta coordinate of the intersection point
!     REAL :: vr, vz      -> r and z coordinates of the intersection point
!-------------------------------------------------------------------------------
subroutine ves_dist_2d(phi, or, oz, ar, az, l0, theta0, l, theta, vr, vz, chi2)

    use magnet_set_globals, only: ves_tol, maxIter

    implicit none

    REAL, intent(IN)  :: phi, or, oz, ar, az, l0, theta0
    REAL, intent(OUT) :: l, theta, vr, vz, chi2
    INTEGER :: i
    REAL :: fr, fz, sl, st
    REAL :: drdt, dzdt
    REAL ::  jac_rl,  jac_rt,  jac_zl,  jac_zt
    REAL :: ijac_lr, ijac_lz, ijac_tr, ijac_tz

    ! constant components of the Jacobian (first column)
    jac_rl = -ar
    jac_zl = -az

    ! Initialize values for which to solve
    l     = l0
    theta = theta0

    do i = 1, maxIter

        vr = ves_r(theta, phi)
        vz = ves_z(theta, phi)

        ! components that should ideally be zero
        fr = vr - (or + l*ar)
        fz = vz - (oz + l*az)

        ! terminate if solution has sufficiently converged
        chi2 = fr*fr + fz*fz
        if (chi2 < ves_tol*ves_tol) exit
        
        ! compute the Jacobian, second column
        drdt = ves_drdt(theta, phi)
        dzdt = ves_dzdt(theta, phi)
        jac_rt = drdt 
        jac_zt = dzdt

        ! invert the Jacobian
        call inverse_2x2( jac_rl,  jac_rt,  jac_zl,  jac_zt, &
                         ijac_lr, ijac_lz, ijac_tr, ijac_tz )

        ! The Newton step
        sl = -ijac_lr * fr - ijac_lz * fz
        st = -ijac_tr * fr - ijac_tz * fz

        ! Update the parameters
        l     = l     + sl
        theta = theta + st

    end do

end subroutine ves_dist_2d

!-------------------------------------------------------------------------------
! ves_perp_intersect_3d(ox, oy, oz, l0, theta0, phi0, 
!                       l, theta, phi, x, y, z, nx, ny, nz)
!
! Determines the parameters of a line segment perpendicular to the vessel that
! intersects a given query point.
!
! Input parameters:
!     REAL :: ox, oy, oz  -> x, y, and z coordinates of the reference point
!     REAL :: l0, theta0, 
!             phi0        -> initial guesses of the distance l, as well as
!                            the theta and phi coordinates of the intersection
!                            point on the vessel
!
! Return parameters:
!     REAL :: l           -> distance from the point to the vessel
!     REAL :: theta, phi  -> theta and phi coordinates of the intersection point
!     REAL :: x, y, z     -> x, y, and z coordinates of the intersection point
!     REAL :: ux, uy, uz  -> x, y, and z components of the unit normal vector
!                            at the vessel intersection point
!-------------------------------------------------------------------------------
subroutine ves_perp_intersect_3d(ox, oy, oz, l0, theta0, phi0, &
                                 l, theta, phi, x, y, z, ux, uy, uz, chi2)

    use magnet_set_globals, only: ves_tol, maxIter

    implicit none

    REAL, intent(IN)  :: ox, oy, oz, l0, theta0, phi0
    REAL, intent(OUT) :: l, theta, phi, x, y, z, ux, uy, uz, chi2
    INTEGER :: i
    REAL :: fx, fy, fz, sl, st, sp
    REAL :: vx, vy, vz, vr, nx, ny, nz, ln 
    REAL :: drdt, drdp, dzdt, dzdp, dxdt, dxdp, dydt, dydp
    REAL :: d2rdt2, d2rdp2, d2zdt2, d2zdp2, d2rdtdp, d2zdtdp
    REAL :: dnxdt, dnydt, dnzdt, dnxdp, dnydp, dnzdp, dlndt, dlndp
    REAL :: duxdt, duydt, duzdt, duxdp, duydp, duzdp
    REAL :: grad_l, grad_t, grad_p, slope
    REAL ::  jac_xl,  jac_xt,  jac_xp, &
             jac_yl,  jac_yt,  jac_yp, &
             jac_zl,  jac_zt,  jac_zp
    REAL :: ijac_lx, ijac_ly, ijac_lz, &
            ijac_px, ijac_py, ijac_pz, &
            ijac_tx, ijac_ty, ijac_tz

    ! Initialize values for which to solve
    l     = l0
    theta = theta0
    phi   = phi0

    call calcF_ves_perp_intersect_3d(l, theta, phi, ox, oy, oz, &
             fx, fy, fz, vr, vx, vy, vz, drdt, dzdt, drdp, dzdp, &
             nx, ny, nz, ln, ux, uy, uz)

    do i = 1, maxIter

        ! terminate if solution has sufficiently converged
        chi2 = fx*fx + fy*fy + fz*fz
        if (chi2 < ves_tol*ves_tol) exit

        ! second derivatives and cross-derivatives for the Jacobian
        d2rdt2  = ves_d2rdt2(theta, phi)
        d2zdt2  = ves_d2zdt2(theta, phi)
        d2rdp2  = ves_d2rdp2(theta, phi)
        d2zdp2  = ves_d2zdp2(theta, phi)
        d2rdtdp = ves_d2rdtdp(theta, phi)
        d2zdtdp = ves_d2zdtdp(theta, phi)

        ! derivatives of x and y coordinates on the vessel
        dxdt = drdt * cos(phi)
        dydt = drdt * sin(phi)
        dxdp = drdp * cos(phi) - vr * sin(phi)
        dydp = drdp * sin(phi) + vr * cos(phi)
        
        ! compute the Jacobian, first column
        jac_xl = ux
        jac_yl = uy
        jac_zl = uz

        ! compute the Jacobian, second column
        dlndt = ( (dzdt*drdp - dzdp*drdt)               &
                      *(  d2zdt2*drdp + dzdt*d2rdtdp &
                        - d2zdtdp*drdt - dzdp*d2rdt2  ) &
                 + vr*drdt * (dzdt**2 + drdt**2)        &
                 + vr**2   * (dzdt*d2zdt2 + drdt*d2rdt2)  ) / ln
        dnxdt = (d2zdt2*drdp + dzdt*d2rdtdp - d2zdtdp*drdt - dzdp*d2rdt2) &
                    * sin(phi) + (d2zdt2*vr + dzdt*drdt) * cos(phi)
        dnydt = (-d2zdt2*drdp - dzdt*d2rdtdp + d2zdtdp*drdt + dzdp*d2rdt2) &
                    * cos(phi) + (d2zdt2*vr + dzdt*drdt) * sin(phi)
        dnzdt = -d2rdt2*vr - drdt*drdt
        duxdt = (ln*dnxdt - nx*dlndt) / ln**2
        duydt = (ln*dnydt - ny*dlndt) / ln**2
        duzdt = (ln*dnzdt - nz*dlndt) / ln**2
        jac_xt = dxdt + l*duxdt
        jac_yt = dydt + l*duydt
        jac_zt = dzdt + l*duzdt

        ! compute the Jacobian, third column
        dlndp = ( (dzdt*drdp - dzdp*drdt)               &
                      *(  d2zdtdp*drdp + dzdt*d2rdp2 &
                        - d2zdp2*drdt - dzdp*d2rdtdp  ) &
                 + vr*drdp * (dzdt**2 + drdt**2)        &
                 + vr**2   * (dzdt*d2zdtdp + drdt*d2rdtdp) ) / ln
        dnxdp = (d2zdtdp*drdp + dzdt*d2rdp2 - dzdt*vr + dzdp*drdt) * sin(phi) &
                 + (d2zdtdp*vr + 2*dzdt*drdp - d2zdp2*drdt - dzdp*d2rdtdp) &
                     * cos(phi)
        dnydp = (d2zdtdp*vr + 2*dzdt*drdp - dzdp*drdt) * sin(phi) &
                 + (-d2zdtdp*drdp - dzdt*d2rdp2 + dzdt*vr + d2zdp2*drdt &
                    + dzdp*d2rdtdp) * cos(phi)
        dnzdp = -d2rdtdp*vr - drdt*drdp
        duxdp = (ln*dnxdp - nx*dlndp) / ln**2
        duydp = (ln*dnydp - ny*dlndp) / ln**2
        duzdp = (ln*dnzdp - nz*dlndp) / ln**2
        jac_xp = dxdp + l*duxdp
        jac_yp = dydp + l*duydp
        jac_zp = dzdp + l*duzdp

        ! invert the Jacobian
        call inverse_3x3(jac_xl, jac_xt, jac_xp, jac_yl, jac_yt, jac_yp, &
                         jac_zl, jac_zt, jac_zp, ijac_lx, ijac_ly, ijac_lz, &
                         ijac_tx, ijac_ty, ijac_tz, ijac_px, ijac_py, ijac_pz)

        ! The Newton step
        sl = -ijac_lx * fx - ijac_ly * fy - ijac_lz * fz
        st = -ijac_tx * fx - ijac_ty * fy - ijac_tz * fz
        sp = -ijac_px * fx - ijac_py * fy - ijac_pz * fz

        grad_l = fx*jac_xl + fy*jac_yl + fz*jac_zl
        grad_t = fx*jac_xt + fy*jac_yt + fz*jac_zt
        grad_p = fx*jac_xp + fy*jac_yp + fz*jac_zp
        slope = grad_l*sl + grad_t*st + grad_p*sp 
        call updateStep_ves_perp_intersect_3d(sl, st, sp, l, theta, phi, &
                 slope, ox, oy, oz, fx, fy, fz, vr, vx, vy, vz, &
                 drdt, dzdt, drdp, dzdp, nx, ny, nz, ln, ux, uy, uz)

    end do

end subroutine ves_perp_intersect_3d

subroutine calcF_ves_perp_intersect_3d(l, theta, phi, ox, oy, oz, &
               fx, fy, fz, vr, vx, vy, vz, drdt, dzdt, drdp, dzdp, &
               nx, ny, nz, ln, ux, uy, uz)

    implicit none

    REAL, intent(IN)  :: l, theta, phi, ox, oy, oz
    REAL, intent(OUT) :: fx, fy, fz, vr, vx, vy, vz, drdt, dzdt, drdp, dzdp, &
                         nx, ny, nz, ln, ux, uy, uz

    vr = ves_r(theta, phi)
    vx = vr * cos(phi)
    vy = vr * sin(phi)
    vz = ves_z(theta, phi)

    ! Derivatives for computing the normal vector
    drdt = ves_drdt(theta, phi)
    dzdt = ves_dzdt(theta, phi)
    drdp = ves_drdp(theta, phi)
    dzdp = ves_dzdp(theta, phi)
    
    ! normal vectors (not unit vectors)
    nx =   dzdt*(drdp*sin(phi) + vr*cos(phi)) - dzdp*drdt*sin(phi)
    ny =  -dzdt*(drdp*cos(phi) - vr*sin(phi)) + dzdp*drdt*cos(phi)
    nz =   drdt*(drdp*cos(phi) - vr*sin(phi))*sin(phi) &
         - drdt*(drdp*sin(phi) + vr*cos(phi))*cos(phi)

    ! unit normal vectors
    ln = sqrt( (dzdt*drdp - dzdp*drdt)**2 + vr**2*(dzdt**2 + drdt**2) )
    ux = nx / ln
    uy = ny / ln
    uz = nz / ln

    ! components that should ideally be zero
    fx = vx + l*ux - ox
    fy = vy + l*uy - oy
    fz = vz + l*uz - oz

end subroutine calcF_ves_perp_intersect_3d

subroutine updateStep_ves_perp_intersect_3d(sl, st, sp, l, theta, phi, &
               slope, ox, oy, oz, fx, fy, fz, vr, vx, vy, vz, &
               drdt, dzdt, drdp, dzdp, nx, ny, nz, ln, ux, uy, uz)

    use magnet_set_globals, only: ves_tol, ves_r00, ves_r10

    implicit none

    REAL, intent(IN)  :: sl, st, sp, slope, ox, oy, oz
    REAL, intent(OUT) :: l, theta, phi, fx, fy, fz, vr, vx, vy, vz, &
                         drdt, dzdt, drdp, dzdp, nx, ny, nz, ln, ux, uy, uz
    REAL :: lambda, fsq, fsq_prev, l_prev, theta_prev, phi_prev, newlambda, &
            lambda2, rhs1, rhs2, a, b, disc, fsq2, test
    REAL :: alpha = 1.0e-4, tol = 1.0e-7
    INTEGER :: j

    ! Adjust the step according to the line search algorithm if necessary
    lambda = 1.0
    j = 0
    fsq_prev = 0.5 * (fx**2 + fy**2 + fz**2)
    l_prev = l
    theta_prev = theta
    phi_prev = phi

    do while (.true.)
        j = j + 1
        l     =     l_prev + lambda * sl
        theta = theta_prev + lambda * st
        phi   =   phi_prev + lambda * sp

        call calcF_ves_perp_intersect_3d(l, theta, phi, ox, oy, oz, &
                 fx, fy, fz, vr, vx, vy, vz, drdt, dzdt, drdp, dzdp, &
                 nx, ny, nz, ln, ux, uy, uz)

        fsq = 0.5 * (fx**2 + fy**2 + fz**2)

        !if (.true.) exit
        if (fsq < fsq_prev + alpha*lambda*slope) exit

        if (lambda*max(abs(sl), abs(st), abs(sp)) < tol) exit

        if (j == 1) then
            newlambda = -slope/(2 * (fsq - fsq_prev - slope))
        else
            rhs1 = fsq - fsq_prev - lambda*slope
            rhs2 = fsq2 - fsq_prev - lambda2*slope
            a = (rhs1/lambda**2 - rhs2/lambda2**2)/(lambda - lambda2)
            b = (-lambda2*rhs1/lambda**2 + lambda*rhs2/lambda2**2) &
                    /(lambda - lambda2)
            if (a == 0.0) then
                newlambda = -slope/(2.0 * b) 
            else
                disc = b**2 - 3.0*a*slope
                if (disc < 0.0) then
                    stop 'vessel_perp_intersect_3d: roundoff problem'
                else
                    newlambda = (-b + sqrt(disc))/(3.0*a)
                end if
            end if
            if (newlambda > 0.5 * lambda) newlambda = 0.5 * lambda
        end if

        lambda2 = lambda
        lambda = max(0.1*lambda, newlambda)
        fsq2 = fsq
        j = j + 1
    end do

end subroutine updateStep_ves_perp_intersect_3d

!-------------------------------------------------------------------------------
! ves_perp_intersect_2d(phi, or, oz, l0, theta0, l, theta, r, z)
!
! Computes the poloidal angle at which a line extending from a given point
! (or, oz) in a poloidal plane interscts the vessel wall at a perpendicular
! angle. Also computes the distance from the input point and the vessel 
! intersection point.
!
! Input parameters:
!     REAL :: phi         -> the toroidal angle (radians) of the cross-section
!     REAL :: or, oz      -> r and z coordinates (meters) of the ref. point
!     REAL :: l0, theta0  -> initial guesses of the distance l, as well as
!                            the theta coordinate of the intersection point
!                            on the vessel
!
! Return parameters:
!     REAL :: l           -> distance from the point to the vessel
!     REAL :: theta       -> theta coordinate of the intersection point
!     REAL :: r, z        -> r and z coordinates of the intersection point
!-------------------------------------------------------------------------------
subroutine ves_perp_intersect_2d(phi, or, oz, l0, theta0, l, theta, r, z, chi2)

    use magnet_set_globals, only: ves_tol, maxIter

    implicit none

    REAL, intent(IN)  :: phi, or, oz, l0, theta0
    REAL, intent(OUT) :: l, theta, r, z, chi2
    INTEGER :: i
    REAL :: vr, vz, nr, nz, fr, fz, sl, st
    REAL :: drdt, dzdt, d2rdt2, d2zdt2, norm, dnorm_dt, dnr_dt, dnz_dt
    REAL ::  jac_rl,  jac_rt,  jac_zl,  jac_zt
    REAL :: ijac_lr, ijac_lz, ijac_tr, ijac_tz

    ! Initialize values for which to solve
    l     = l0
    theta = theta0

    do i = 1, maxIter

        ! Compute vessel coordinates and derivatives
        vr = ves_r(theta, phi)
        vz = ves_z(theta, phi)
        drdt = ves_drdt(theta, phi)
        dzdt = ves_dzdt(theta, phi)
        d2rdt2 = ves_d2rdt2(theta, phi)
        d2zdt2 = ves_d2zdt2(theta, phi)

        ! Unit vector normal to the vessel
        norm = sqrt(drdt**2 + dzdt**2)
        nr = dzdt/norm
        nz = -drdt/norm

        ! components that should ideally be zero
        fr = vr + l*nr - or
        fz = vz + l*nz - oz

        ! terminate if solution has sufficiently converged
        chi2 = fr*fr + fz*fz
        if (chi2 < ves_tol*ves_tol) exit

        ! Compute the Jacobian, first column
        jac_rl = nr
        jac_zl = nz
        
        ! compute the Jacobian, second column
        dnorm_dt = (drdt*d2rdt2 + dzdt*d2zdt2) / norm
        dnr_dt =  d2zdt2/norm - dzdt*dnorm_dt/norm**2
        dnz_dt = -d2rdt2/norm + drdt*dnorm_dt/norm**2
        jac_rt = drdt + l*dnr_dt
        jac_zt = dzdt + l*dnz_dt

        ! invert the Jacobian
        call inverse_2x2( jac_rl,  jac_rt,  jac_zl,  jac_zt, &
                         ijac_lr, ijac_lz, ijac_tr, ijac_tz )

        ! The Newton step
        sl = -ijac_lr * fr - ijac_lz * fz
        st = -ijac_tr * fr - ijac_tz * fz

        ! Update the parameters
        l     = l     + sl
        theta = theta + st

    end do

end subroutine ves_perp_intersect_2d

!-------------------------------------------------------------------------------
! ves_unorm(phi, r, drdt, drdp, dzdt, dzdp, ux, uy, uz)
!
! Computes the unit normal vector at a given point characterized by theta and
! phi. 
!
! Note: theta is not explicitly provided as an input, other inputs depend 
! implicitly on theta.
!
! Input parameters
!     REAL :: phi        -> toroidal angle (radians) where the vector is to be 
!                           computed
!     REAL :: r          -> radial coordinate (meters) where is the vector is 
!                           computed, output by ves_r(theta, phi)
!     REAL :: drdt, drdp -> derivatives of the radial coordinate with respect
!                           to theta and phi, output by ves_drdt(theta, phi)
!                           and ves_drdp(theta, phi)
!     REAL :: dzdt, dzdp -> derivatives of the z coordinate with respect 
!                           to theta and phi, output by ves_dzdt(theta, phi)
!                           and ves_dzdp(theta, phi)
!
! Output parameters:
!     REAL :: ux, uy, uz -> x, y, and z components of the unit normal vector
!                           at poloidal angle theta, and toroidal angle phi
!-------------------------------------------------------------------------------
subroutine ves_unorm(phi, r, drdt, drdp, dzdt, dzdp, ux, uy, uz)

    implicit none

    REAL, intent(IN) :: phi, r, drdt, drdp, dzdt, dzdp
    REAL, intent(OUT) :: ux, uy, uz
    REAL :: dxdt, dxdp, dydp, dydt, vx, vy, vz, length

    ! x = r*cos(phi)
    ! y = r*sin(phi)
    dxdt = drdt*cos(phi)
    dxdp = drdp*cos(phi) - r*sin(phi)

    dydt = drdt*sin(phi)
    dydp = drdp*sin(phi) + r*cos(phi)

    call cross_prod(dxdp, dydp, dzdp, dxdt, dydt, dzdt, vx, vy, vz)

    length = sqrt(vx*vx + vy*vy + vz*vz)

    ux = vx / length
    uy = vy / length
    uz = vz / length

end subroutine ves_unorm

!-------------------------------------------------------------------------------
! ves_r(theta, phi)
!
! Computes the r coordinate on the vessel at poloidal angle theta and toroidal
! angle phi.
!-------------------------------------------------------------------------------
REAL function ves_r(theta, phi)

    use magnet_set_globals, only: nModes, nfp, vrc, vrs, vm, vn

    implicit none

    REAL, intent(IN) :: theta, phi
    INTEGER :: i

    ves_r = 0
    do i = 1, nModes
        ves_r = ves_r + vrc(i)*cos(vm(i)*theta + vn(i)*nfp*phi) &
                      + vrs(i)*sin(vm(i)*theta + vn(i)*nfp*phi)
    end do

end function ves_r

!-------------------------------------------------------------------------------
! ves_z(theta, phi)
!
! Computes the r coordinate on the vessel at poloidal angle theta and toroidal
! angle phi.
!-------------------------------------------------------------------------------
REAL function ves_z(theta, phi)

    use magnet_set_globals, only: nModes, nfp, vzs, vzc, vm, vn

    implicit none

    REAL, intent(IN) :: theta, phi
    INTEGER :: i

    ves_z = 0
    do i = 1, nModes
        ves_z = ves_z + vzs(i)*sin(vm(i)*theta + vn(i)*nfp*phi) &
                      + vzc(i)*cos(vm(i)*theta + vn(i)*nfp*phi)
    end do

end function ves_z

!-------------------------------------------------------------------------------
! ves_drdt(theta, phi)
!
! Computes the derivative of the r coordinate of the vessel boundary with
! respect to the poloidal angle at the location given by theta and phi.
!-------------------------------------------------------------------------------
REAL function ves_drdt(theta, phi)

    use magnet_set_globals, only: nModes, nfp, vrc, vrs, vm, vn

    implicit none

    REAL, intent(IN) :: theta, phi
    INTEGER :: i

    ves_drdt = 0
    do i = 1, nModes
        ves_drdt = ves_drdt - vm(i)*vrc(i)*sin(vm(i)*theta + vn(i)*nfp*phi) &
                            + vm(i)*vrs(i)*cos(vm(i)*theta + vn(i)*nfp*phi)
    end do

end function ves_drdt

!-------------------------------------------------------------------------------
! ves_drdp(theta, phi)
!
! Computes the derivative of the r coordinate of the vessel boundary with
! respect to the toroidal angle at the location given by theta and phi.
!-------------------------------------------------------------------------------
REAL function ves_drdp(theta, phi)

    use magnet_set_globals, only: nModes, nfp, vrc, vrs, vm, vn

    implicit none

    REAL, intent(IN) :: theta, phi
    INTEGER :: i

    ves_drdp = 0
    do i = 1, nModes
        ves_drdp = ves_drdp - vn(i)*vrc(i)*sin(vm(i)*theta + vn(i)*nfp*phi) &
                            + vn(i)*vrs(i)*cos(vm(i)*theta + vn(i)*nfp*phi)
    end do

    ves_drdp = ves_drdp * nfp

end function ves_drdp

!-------------------------------------------------------------------------------
! ves_dzdt(theta, phi)
!
! Computes the derivative of the z coordinate of the vessel boundary with
! respect to the poloidal angle at the location given by theta and phi.
!-------------------------------------------------------------------------------
REAL function ves_dzdt(theta, phi)

    use magnet_set_globals, only: nModes, nfp, vzs, vzc, vm, vn

    implicit none

    REAL, intent(IN) :: theta, phi
    INTEGER :: i

    ves_dzdt = 0
    do i = 1, nModes
        ves_dzdt = ves_dzdt + vm(i)*vzs(i)*cos(vm(i)*theta + vn(i)*nfp*phi) &
                            - vm(i)*vzc(i)*sin(vm(i)*theta + vn(i)*nfp*phi)
    end do

end function ves_dzdt

!-------------------------------------------------------------------------------
! ves_dzdp(theta, phi)
!
! Computes the derivative of the z coordinate of the vessel boundary with
! respect to the toroidal angle at the location given by theta and phi.
!-------------------------------------------------------------------------------
REAL function ves_dzdp(theta, phi)

    use magnet_set_globals, only: nModes, nfp, vzs, vzc, vm, vn

    implicit none

    REAL, intent(IN) :: theta, phi
    INTEGER :: i

    ves_dzdp = 0
    do i = 1, nModes
        ves_dzdp = ves_dzdp + vn(i)*vzs(i)*cos(vm(i)*theta + vn(i)*nfp*phi) &
                            - vn(i)*vzc(i)*sin(vm(i)*theta + vn(i)*nfp*phi)
    end do

    ves_dzdp = ves_dzdp * nfp

end function ves_dzdp

!-------------------------------------------------------------------------------
! ves_d2rdt2(theta, phi)
!
! Computes the second derivative of the r coordinate of the vessel boundary with
! respect to the poloidal angle at the location given by theta and phi.
!-------------------------------------------------------------------------------
REAL function ves_d2rdt2(theta, phi)

    use magnet_set_globals, only: nModes, nfp, vrc, vrs, vm, vn

    implicit none

    REAL, intent(IN) :: theta, phi
    INTEGER :: i

    ves_d2rdt2 = 0
    do i = 1, nModes
        ves_d2rdt2 = ves_d2rdt2 &
                         - vm(i)**2*vrc(i)*cos(vm(i)*theta + vn(i)*nfp*phi) &
                         - vm(i)**2*vrs(i)*sin(vm(i)*theta + vn(i)*nfp*phi)
    end do


end function ves_d2rdt2

!-------------------------------------------------------------------------------
! ves_d2rdp2(theta, phi)
!
! Computes the second derivative of the r coordinate of the vessel boundary with
! respect to the toroidal angle at the location given by theta and phi.
!-------------------------------------------------------------------------------
REAL function ves_d2rdp2(theta, phi)

    use magnet_set_globals, only: nModes, nfp, vrc, vrs, vm, vn

    implicit none

    REAL, intent(IN) :: theta, phi
    INTEGER :: i

    ves_d2rdp2 = 0
    do i = 1, nModes
        ves_d2rdp2 = ves_d2rdp2 &
                         - vn(i)**2*vrc(i)*cos(vm(i)*theta + vn(i)*nfp*phi) &
                         - vn(i)**2*vrs(i)*sin(vm(i)*theta + vn(i)*nfp*phi)
    end do

    ves_d2rdp2 = ves_d2rdp2 * nfp**2

end function ves_d2rdp2

!-------------------------------------------------------------------------------
! ves_d2zdt2(theta, phi)
!
! Computes the second derivative of the z coordinate of the vessel boundary with
! respect to the toroidal angle at the location given by theta and phi.
!-------------------------------------------------------------------------------
REAL function ves_d2zdt2(theta, phi)

    use magnet_set_globals, only: nModes, nfp, vzs, vzc, vm, vn

    implicit none

    REAL, intent(IN) :: theta, phi
    INTEGER :: i

    ves_d2zdt2 = 0
    do i = 1, nModes
        ves_d2zdt2 = ves_d2zdt2 &
                         - vm(i)**2*vzs(i)*sin(vm(i)*theta + vn(i)*nfp*phi) &
                         - vm(i)**2*vzc(i)*cos(vm(i)*theta + vn(i)*nfp*phi)
    end do

end function ves_d2zdt2

!-------------------------------------------------------------------------------
! ves_d2zdp2(theta, phi)
!
! Computes the second derivative of the z coordinate of the vessel boundary with
! respect to the toroidal angle at the location given by theta and phi.
!-------------------------------------------------------------------------------
REAL function ves_d2zdp2(theta, phi)

    use magnet_set_globals, only: nModes, nfp, vzs, vzc, vm, vn

    implicit none

    REAL, intent(IN) :: theta, phi
    INTEGER :: i

    ves_d2zdp2 = 0
    do i = 1, nModes
        ves_d2zdp2 = ves_d2zdp2 &
                         - vn(i)**2*vzs(i)*sin(vm(i)*theta + vn(i)*nfp*phi) &
                         - vn(i)**2*vzc(i)*cos(vm(i)*theta + vn(i)*nfp*phi)
    end do

    ves_d2zdp2 = ves_d2zdp2 * nfp**2

end function ves_d2zdp2

!-------------------------------------------------------------------------------
! ves_d2rdtdp(theta, phi)
!
! Computes the cross-derivative of the r coordinate of the vessel boundary with
! respect to the poloidal and toroidal angles at the location given by theta 
! and phi.
!-------------------------------------------------------------------------------
REAL function ves_d2rdtdp(theta, phi)

    use magnet_set_globals, only: nModes, nfp, vrc, vrs, vm, vn

    implicit none

    REAL, intent(IN) :: theta, phi
    INTEGER :: i

    ves_d2rdtdp = 0
    do i = 1, nModes
        ves_d2rdtdp = ves_d2rdtdp &
            - nfp*vn(i)*vm(i)*vrc(i)*cos(vm(i)*theta + vn(i)*nfp*phi) &
            - nfp*vn(i)*vm(i)*vrs(i)*sin(vm(i)*theta + vn(i)*nfp*phi)
    end do

end function ves_d2rdtdp

!-------------------------------------------------------------------------------
! ves_d2zdtdp(theta, phi)
!
! Computes the cross-derivative of the z coordinate of the vessel boundary with
! respect to the poloidal and toroidal angles at the location given by theta 
! and phi.
!-------------------------------------------------------------------------------
REAL function ves_d2zdtdp(theta, phi)

    use magnet_set_globals, only: nModes, nfp, vzs, vzc, vm, vn

    implicit none

    REAL, intent(IN) :: theta, phi
    INTEGER :: i

    ves_d2zdtdp = 0
    do i = 1, nModes
        ves_d2zdtdp = ves_d2zdtdp &
            - nfp*vn(i)*vm(i)*vzs(i)*sin(vm(i)*theta + vn(i)*nfp*phi) &
            - nfp*vn(i)*vm(i)*vzc(i)*cos(vm(i)*theta + vn(i)*nfp*phi)
    end do

end function ves_d2zdtdp

!-------------------------------------------------------------------------------
! plas_r(theta, phi)
!
! Computes the r coordinate on the plasma lcfs at poloidal angle theta and 
! toroidal angle phi.
!-------------------------------------------------------------------------------
REAL function plas_r(theta, phi)

    use magnet_set_globals, only: nModesPl, nfp, prc, prs, pm, pn

    implicit none

    REAL, intent(IN) :: theta, phi
    INTEGER :: i

    plas_r = 0
    do i = 1, nModesPl
        plas_r = plas_r + prc(i)*cos(pm(i)*theta - pn(i)*nfp*phi) &
                        + prs(i)*sin(pm(i)*theta - pn(i)*nfp*phi)
    end do

end function plas_r

!-------------------------------------------------------------------------------
! plas_z(theta, phi)
!
! Computes the z coordinate on the plasma lcfs at poloidal angle theta and 
! toroidal angle phi.
!-------------------------------------------------------------------------------
REAL function plas_z(theta, phi)

    use magnet_set_globals, only: nModesPl, nfp, pzs, pzc, pm, pn

    implicit none

    REAL, intent(IN) :: theta, phi
    INTEGER :: i

    plas_z = 0
    do i = 1, nModesPl
        plas_z = plas_z + pzs(i)*sin(pm(i)*theta - pn(i)*nfp*phi) &
                        + pzc(i)*cos(pm(i)*theta - pn(i)*nfp*phi)
    end do

end function plas_z


!-------------------------------------------------------------------------------
! plas_drdt(theta, phi)
!
! Computes the derivative of the r coordinate on the plasma lcfs with respect
! to poloigal angle at poloidal angle theta and toroidal angle phi.
!-------------------------------------------------------------------------------
REAL function plas_drdt(theta, phi)

    use magnet_set_globals, only: nModesPl, nfp, prc, prs, pm, pn

    implicit none

    REAL, intent(IN) :: theta, phi
    INTEGER :: i

    plas_drdt = 0
    do i = 1, nModesPl
        plas_drdt = plas_drdt - pm(i)*prc(i)*sin(pm(i)*theta - pn(i)*nfp*phi) &
                              + pm(i)*prs(i)*cos(pm(i)*theta - pn(i)*nfp*phi)
    end do

end function plas_drdt

!-------------------------------------------------------------------------------
! plas_dzdt(theta, phi)
!
! Computes the derivative of the z coordinate on the plasma lcfs with respect
! to poloidal angle at poloidal angle theta and toroidal angle phi.
!-------------------------------------------------------------------------------
REAL function plas_dzdt(theta, phi)

    use magnet_set_globals, only: nModesPl, nfp, pzs, pzc, pm, pn

    implicit none

    REAL, intent(IN) :: theta, phi
    INTEGER :: i

    plas_dzdt = 0
    do i = 1, nModesPl
        plas_dzdt = plas_dzdt + pm(i)*pzs(i)*cos(pm(i)*theta - pn(i)*nfp*phi) &
                              - pm(i)*pzc(i)*sin(pm(i)*theta - pn(i)*nfp*phi)
    end do

end function plas_dzdt



!-------------------------------------------------------------------------------
! cross_prod(ax, ay, az, bx, by, bz, cx, cy, cz)
! 
! Computes the cross product of two Cartesian 3-dimensional vectors
!
! Input parameters
!     REAL :: ax, ay, az  -> x, y, and z components of the first vector
!     REAL :: bx, by, bz  -> x, y, and z components of the second vector
!
! Output parameters
!     REAL :: cx, cy, cz  -> x, y, and z components of a cross b
!-------------------------------------------------------------------------------
subroutine cross_prod(ax, ay, az, bx, by, bz, cx, cy, cz)

    implicit none

    REAL, intent(IN) :: ax, ay, az, bx, by, bz
    REAL, intent(OUT) :: cx, cy, cz

    cx =  ay*bz - az*by
    cy = -ax*bz + az*bx
    cz =  ax*by - ay*bx

end subroutine cross_prod

!-------------------------------------------------------------------------------
! inverse_3x3(m11, m12, m13, m21, m22, m23, m31, m32, m33, 
!             i11, i12, i13, i21, i22, i23, i31, i32, i33)
!
! Computes the inverse of a 3x3 matrix.
! Input variables are the components m_ij of the matrix to invert.
! Output variables are the components i_ij of the inverse.
!
! Source: Wolfram MathWorld
!-------------------------------------------------------------------------------
subroutine inverse_3x3(m11, m12, m13, m21, m22, m23, m31, m32, m33, &
                       i11, i12, i13, i21, i22, i23, i31, i32, i33)

    implicit none

    REAL, intent(IN)  :: m11, m12, m13, m21, m22, m23, m31, m32, m33
    REAL, intent(OUT) :: i11, i12, i13, i21, i22, i23, i31, i32, i33
    REAL :: det

    det = m11*(m22*m33-m23*m32) - m12*(m21*m33-m23*m31) + m13*(m21*m32-m22*m31)

    i11 = (m22*m33 - m23*m32) / det
    i12 = (m13*m32 - m12*m33) / det
    i13 = (m12*m23 - m13*m22) / det
    i21 = (m23*m31 - m21*m33) / det
    i22 = (m11*m33 - m13*m31) / det
    i23 = (m13*m21 - m11*m23) / det
    i31 = (m21*m32 - m22*m31) / det
    i32 = (m12*m31 - m11*m32) / det
    i33 = (m11*m22 - m12*m21) / det

end subroutine inverse_3x3

!-------------------------------------------------------------------------------
! inverse_2x2(m11, m12, m21, m22, i11, i12, i21, i22)
!
! Computes the inverse of a 2x2 matrix.
! Input variables are the components m_ij of the matrix to invert.
! Output variables are the components i_ij of the inverse.
!
! Source: Wolfram MathWorld
!-------------------------------------------------------------------------------
subroutine inverse_2x2(m11, m12, m21, m22, i11, i12, i21, i22)

    implicit none

    REAL, intent(IN)  :: m11, m12, m21, m22
    REAL, intent(OUT) :: i11, i12, i21, i22
    REAL :: det

    det = m11*m22 - m12*m21

    i11 =  m22 / det
    i12 = -m12 / det
    i21 = -m21 / det
    i22 =  m11 / det

end subroutine inverse_2x2

!-------------------------------------------------------------------------------
! product_3x3(a11, a12, a13, a21, a22, a23, a31, a32, a33, 
!             b11, b12, b13, b21, b22, b23, b31, b32, b33,
!             p11, p12, p13, p21, p22, p23, p31, p32, p33)
!
! Computes the product p of two 3x3 matrices a and b
!-------------------------------------------------------------------------------
subroutine product_3x3(a11, a12, a13, a21, a22, a23, a31, a32, a33, &
                       b11, b12, b13, b21, b22, b23, b31, b32, b33, &
                       p11, p12, p13, p21, p22, p23, p31, p32, p33)

    implicit none

    REAL, intent(IN)  :: a11, a12, a13, a21, a22, a23, a31, a32, a33
    REAL, intent(IN)  :: b11, b12, b13, b21, b22, b23, b31, b32, b33
    REAL, intent(OUT) :: p11, p12, p13, p21, p22, p23, p31, p32, p33

    p11 = a11*b11 + a12*b21 + a13*b31
    p12 = a11*b12 + a12*b22 + a13*b32
    p13 = a11*b13 + a12*b23 + a13*b33
    p21 = a21*b11 + a22*b21 + a23*b31
    p22 = a21*b12 + a22*b22 + a23*b32
    p23 = a21*b13 + a22*b23 + a23*b33
    p31 = a31*b11 + a32*b21 + a33*b31
    p32 = a31*b12 + a32*b22 + a33*b32
    p33 = a31*b13 + a32*b23 + a33*b33

end subroutine product_3x3

end module magnet_set_calc

