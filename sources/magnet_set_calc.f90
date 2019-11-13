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
!     REAL(8) :: ox, oy, oz  -> x, y, and z coordinates of the reference point
!     REAL(8) :: ax, ay, az  -> x, y, and z components of a unit vector defining
!                            the direction from the test point to the vessel
!     REAL(8) :: l0, theta0, 
!             phi0        -> initial guesses of the distance l, as well as
!                            the theta and phi coordinates of the intersection
!                            point on the vessel
!
! Return parameters:
!     REAL(8) :: l           -> distance from the point to the vessel
!     REAL(8) :: theta, phi  -> theta and phi coordinates of the intersection point
!     REAL(8) :: vx, vy, vz  -> x, y, and z coordinates of the intersection point
!-------------------------------------------------------------------------------
subroutine ves_dist_3d(ox, oy, oz, ax, ay, az, l0, theta0, phi0, &
                       l, theta, phi, vx, vy, vz, chi2)

    use magnet_set_globals, only: ves_tol, maxIter

    implicit none

    REAL(8), intent(IN)  :: ox, oy, oz, ax, ay, az, l0, theta0, phi0
    REAL(8), intent(OUT) :: l, theta, phi, vx, vy, vz, chi2
    INTEGER :: i
    REAL(8) :: vr, fx, fy, fz, sl, st, sp
    REAL(8) :: drdt, drdp, dzdt, dzdp
    REAL(8) ::  jac_xl,  jac_xt,  jac_xp, &
             jac_yl,  jac_yt,  jac_yp, &
             jac_zl,  jac_zt,  jac_zp
    REAL(8) :: ijac_lx, ijac_ly, ijac_lz, &
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
!     REAL(8) :: phi         -> the toroidal angle (radians) of the cross-section
!     REAL(8) :: or, oz      -> r and z coordinates (meters) of the ref. point
!     REAL(8) :: ar, az      -> r and z components of a unit vector defining
!                            the direction from the test point to the vessel
!     REAL(8) :: l0, theta0  -> initial guesses of the distance l, as well as
!                            the theta coordinate of the intersection point
!                            on the vessel
!
! Return parameters:
!     REAL(8) :: l           -> distance from the point to the vessel
!     REAL(8) :: theta       -> theta coordinate of the intersection point
!     REAL(8) :: vr, vz      -> r and z coordinates of the intersection point
!-------------------------------------------------------------------------------
subroutine ves_dist_2d(phi, or, oz, ar, az, l0, theta0, l, theta, vr, vz, chi2)

    use magnet_set_globals, only: ves_tol, maxIter

    implicit none

    REAL(8), intent(IN)  :: phi, or, oz, ar, az, l0, theta0
    REAL(8), intent(OUT) :: l, theta, vr, vz, chi2
    INTEGER :: i
    REAL(8) :: fr, fz, sl, st
    REAL(8) :: drdt, dzdt
    REAL(8) ::  jac_rl,  jac_rt,  jac_zl,  jac_zt
    REAL(8) :: ijac_lr, ijac_lz, ijac_tr, ijac_tz

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
!     REAL(8) :: ox, oy, oz  -> x, y, and z coordinates of the reference point
!     REAL(8) :: l0, theta0, 
!             phi0        -> initial guesses of the distance l, as well as
!                            the theta and phi coordinates of the intersection
!                            point on the vessel
!
! Return parameters:
!     REAL(8) :: l           -> distance from the point to the vessel
!     REAL(8) :: theta, phi  -> theta and phi coordinates of the intersection point
!     REAL(8) :: vx, vy, vz  -> x, y, and z coordinates of the intersection point
!     REAL(8) :: ux, uy, uz  -> x, y, and z components of the unit normal vector
!                            at the vessel intersection point
!-------------------------------------------------------------------------------
subroutine ves_perp_intersect_3d(ox, oy, oz, l0, theta0, phi0, &
                                 l, theta, phi, vx, vy, vz, ux, uy, uz, chi2)

    use magnet_set_globals, only: ves_tol, maxIter

    implicit none

    REAL(8), intent(IN)  :: ox, oy, oz, l0, theta0, phi0
    REAL(8), intent(OUT) :: l, theta, phi, vx, vy, vz, ux, uy, uz, chi2
    INTEGER :: i
    REAL(8) :: fx, fy, fz, sl, st, sp
    REAL(8) :: vr, nx, ny, nz, ln 
    REAL(8) :: drdt, drdp, dzdt, dzdp, dxdt, dxdp, dydt, dydp
    REAL(8) :: d2rdt2, d2rdp2, d2zdt2, d2zdp2, d2rdtdp, d2zdtdp
    REAL(8) :: dnxdt, dnydt, dnzdt, dnxdp, dnydp, dnzdp, dlndt, dlndp
    REAL(8) :: duxdt, duydt, duzdt, duxdp, duydp, duzdp
    REAL(8) :: grad_l, grad_t, grad_p, slope
    REAL(8) ::  jac_xl,  jac_xt,  jac_xp, &
             jac_yl,  jac_yt,  jac_yp, &
             jac_zl,  jac_zt,  jac_zp
    REAL(8) :: ijac_lx, ijac_ly, ijac_lz, &
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

    REAL(8), intent(IN)  :: l, theta, phi, ox, oy, oz
    REAL(8), intent(OUT) :: fx, fy, fz, vr, vx, vy, vz, drdt, dzdt, drdp, dzdp, &
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

    REAL(8), intent(IN)  :: sl, st, sp, slope, ox, oy, oz
    REAL(8), intent(OUT) :: l, theta, phi, fx, fy, fz, vr, vx, vy, vz, &
                         drdt, dzdt, drdp, dzdp, nx, ny, nz, ln, ux, uy, uz
    REAL(8) :: lambda, fsq, fsq_prev, l_prev, theta_prev, phi_prev, newlambda, &
            lambda2, rhs1, rhs2, a, b, disc, fsq2, test
    REAL(8) :: alpha = 1.0e-4, tol = 1.0e-7
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
!     REAL(8) :: phi         -> the toroidal angle (radians) of the cross-section
!     REAL(8) :: or, oz      -> r and z coordinates (meters) of the ref. point
!     REAL(8) :: l0, theta0  -> initial guesses of the distance l, as well as
!                            the theta coordinate of the intersection point
!                            on the vessel
!
! Return parameters:
!     REAL(8) :: l           -> distance from the point to the vessel
!     REAL(8) :: theta       -> theta coordinate of the intersection point
!     REAL(8) :: r, z        -> r and z coordinates of the intersection point
!-------------------------------------------------------------------------------
subroutine ves_perp_intersect_2d(phi, or, oz, l0, theta0, l, theta, r, z, chi2)

    use magnet_set_globals, only: ves_tol, maxIter

    implicit none

    REAL(8), intent(IN)  :: phi, or, oz, l0, theta0
    REAL(8), intent(OUT) :: l, theta, r, z, chi2
    INTEGER :: i
    REAL(8) :: vr, vz, nr, nz, fr, fz, sl, st
    REAL(8) :: drdt, dzdt, d2rdt2, d2zdt2, norm, dnorm_dt, dnr_dt, dnz_dt
    REAL(8) ::  jac_rl,  jac_rt,  jac_zl,  jac_zt
    REAL(8) :: ijac_lr, ijac_lz, ijac_tr, ijac_tz

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
!     REAL(8) :: phi        -> toroidal angle (radians) where the vector is to be 
!                           computed
!     REAL(8) :: r          -> radial coordinate (meters) where is the vector is 
!                           computed, output by ves_r(theta, phi)
!     REAL(8) :: drdt, drdp -> derivatives of the radial coordinate with respect
!                           to theta and phi, output by ves_drdt(theta, phi)
!                           and ves_drdp(theta, phi)
!     REAL(8) :: dzdt, dzdp -> derivatives of the z coordinate with respect 
!                           to theta and phi, output by ves_dzdt(theta, phi)
!                           and ves_dzdp(theta, phi)
!
! Output parameters:
!     REAL(8) :: ux, uy, uz -> x, y, and z components of the unit normal vector
!                           at poloidal angle theta, and toroidal angle phi
!-------------------------------------------------------------------------------
subroutine ves_unorm(phi, r, drdt, drdp, dzdt, dzdp, ux, uy, uz)

    implicit none

    REAL(8), intent(IN) :: phi, r, drdt, drdp, dzdt, dzdp
    REAL(8), intent(OUT) :: ux, uy, uz
    REAL(8) :: dxdt, dxdp, dydp, dydt, vx, vy, vz, length

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
REAL(8) function ves_r(theta, phi)

    use magnet_set_globals, only: nModes, nfp, vrc, vrs, vm, vn

    implicit none

    REAL(8), intent(IN) :: theta, phi
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
REAL(8) function ves_z(theta, phi)

    use magnet_set_globals, only: nModes, nfp, vzs, vzc, vm, vn

    implicit none

    REAL(8), intent(IN) :: theta, phi
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
REAL(8) function ves_drdt(theta, phi)

    use magnet_set_globals, only: nModes, nfp, vrc, vrs, vm, vn

    implicit none

    REAL(8), intent(IN) :: theta, phi
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
REAL(8) function ves_drdp(theta, phi)

    use magnet_set_globals, only: nModes, nfp, vrc, vrs, vm, vn

    implicit none

    REAL(8), intent(IN) :: theta, phi
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
REAL(8) function ves_dzdt(theta, phi)

    use magnet_set_globals, only: nModes, nfp, vzs, vzc, vm, vn

    implicit none

    REAL(8), intent(IN) :: theta, phi
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
REAL(8) function ves_dzdp(theta, phi)

    use magnet_set_globals, only: nModes, nfp, vzs, vzc, vm, vn

    implicit none

    REAL(8), intent(IN) :: theta, phi
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
REAL(8) function ves_d2rdt2(theta, phi)

    use magnet_set_globals, only: nModes, nfp, vrc, vrs, vm, vn

    implicit none

    REAL(8), intent(IN) :: theta, phi
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
REAL(8) function ves_d2rdp2(theta, phi)

    use magnet_set_globals, only: nModes, nfp, vrc, vrs, vm, vn

    implicit none

    REAL(8), intent(IN) :: theta, phi
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
REAL(8) function ves_d2zdt2(theta, phi)

    use magnet_set_globals, only: nModes, nfp, vzs, vzc, vm, vn

    implicit none

    REAL(8), intent(IN) :: theta, phi
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
REAL(8) function ves_d2zdp2(theta, phi)

    use magnet_set_globals, only: nModes, nfp, vzs, vzc, vm, vn

    implicit none

    REAL(8), intent(IN) :: theta, phi
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
REAL(8) function ves_d2rdtdp(theta, phi)

    use magnet_set_globals, only: nModes, nfp, vrc, vrs, vm, vn

    implicit none

    REAL(8), intent(IN) :: theta, phi
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
REAL(8) function ves_d2zdtdp(theta, phi)

    use magnet_set_globals, only: nModes, nfp, vzs, vzc, vm, vn

    implicit none

    REAL(8), intent(IN) :: theta, phi
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
REAL(8) function plas_r(theta, phi)

    use magnet_set_globals, only: nModesPl, nfp, prc, prs, pm, pn

    implicit none

    REAL(8), intent(IN) :: theta, phi
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
REAL(8) function plas_z(theta, phi)

    use magnet_set_globals, only: nModesPl, nfp, pzs, pzc, pm, pn

    implicit none

    REAL(8), intent(IN) :: theta, phi
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
REAL(8) function plas_drdt(theta, phi)

    use magnet_set_globals, only: nModesPl, nfp, prc, prs, pm, pn

    implicit none

    REAL(8), intent(IN) :: theta, phi
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
REAL(8) function plas_dzdt(theta, phi)

    use magnet_set_globals, only: nModesPl, nfp, pzs, pzc, pm, pn

    implicit none

    REAL(8), intent(IN) :: theta, phi
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
!     REAL(8) :: ax, ay, az  -> x, y, and z components of the first vector
!     REAL(8) :: bx, by, bz  -> x, y, and z components of the second vector
!
! Output parameters
!     REAL(8) :: cx, cy, cz  -> x, y, and z components of a cross b
!-------------------------------------------------------------------------------
subroutine cross_prod(ax, ay, az, bx, by, bz, cx, cy, cz)

    implicit none

    REAL(8), intent(IN) :: ax, ay, az, bx, by, bz
    REAL(8), intent(OUT) :: cx, cy, cz

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

    REAL(8), intent(IN)  :: m11, m12, m13, m21, m22, m23, m31, m32, m33
    REAL(8), intent(OUT) :: i11, i12, i13, i21, i22, i23, i31, i32, i33
    REAL(8) :: det

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

    REAL(8), intent(IN)  :: m11, m12, m21, m22
    REAL(8), intent(OUT) :: i11, i12, i21, i22
    REAL(8) :: det

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

    REAL(8), intent(IN)  :: a11, a12, a13, a21, a22, a23, a31, a32, a33
    REAL(8), intent(IN)  :: b11, b12, b13, b21, b22, b23, b31, b32, b33
    REAL(8), intent(OUT) :: p11, p12, p13, p21, p22, p23, p31, p32, p33

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

!-------------------------------------------------------------------------------
! plane_3point(x1, y1, z1, x2, y2, z2, x3, y3, z3, nx, ny, nz)
!
! Calculates the normal vector of a plane that includes three input points.
! The vector is oriented as the cross-product of the vectors from the second
! point to the first point, and from the second point to the third point,
! respectively.
! 
! Input parameters:
!     REAL(8) :: x1, x2, x3   -> x coordinates of the three points
!     REAL(8) :: y1, y2, y3   -> y coordinates of the three points
!     REAL(8) :: z1, z2, z3   -> z coordinates of the three points
!
! Return parameters:
!     REAL(8) :: nx, ny, nz   -> x, y, and z components of the unit normal
!-------------------------------------------------------------------------------
subroutine plane_3point(x1, y1, z1, x2, y2, z2, x3, y3, z3, nx, ny, nz)

    implicit none

    REAL(8), intent(IN)  :: x1, y1, z1, x2, y2, z2, x3, y3, z3
    REAL(8), intent(OUT) :: nx, ny, nz
    REAL(8) :: v1_x, v1_y, v1_z, v2_x, v2_y, v2_z, cprod_x, cprod_y, cprod_z, norm

    v1_x = x1 - x2
    v1_y = y1 - y2
    v1_z = z1 - z2
    v2_x = x3 - x2
    v2_y = y3 - y2
    v2_z = z3 - z2

    cprod_x =  v1_y * v2_z - v1_z * v2_y
    cprod_y = -v1_x * v2_z + v1_z * v2_x
    cprod_z =  v1_x * v2_y - v1_y * v2_x

    norm = sqrt(cprod_x**2 + cprod_y**2 + cprod_z**2)

    nx = cprod_x/norm
    ny = cprod_y/norm
    nz = cprod_z/norm   

end subroutine plane_3point

!-------------------------------------------------------------------------------
! plane_4point(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, 
!              nx, ny, nz, ox, oy, oz)
!
! Calculates the normal vector and a reference point for a plane that fits
! best (in a least-squares sense) to a set of four input points. The sign of
! the normal vector will be determined by the handedness of the first three
! points supplied.
!
! Input parameters:
!     REAL(8) :: x1, y1, z1   -> x, y, and z coordinates of the first point
!     REAL(8) :: x2, y2, z2   -> x, y, and z coordinates of the second point
!     REAL(8) :: x3, y3, z3   -> x, y, and z coordinates of the third point
!     REAL(8) :: x4, y4, z4   -> x, y, and z coordinates of the fourth point
!
! Output parameters:
!     REAL(8) :: nx, ny, nz   -> x, y, and z components of the unit normal
!     REAL(8) :: ox, oy, oz   -> x, y, and z coordinates of a reference point
!-------------------------------------------------------------------------------
subroutine plane_4point(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, &
                        nx, ny, nz, ox, oy, oz)

    implicit none

    REAL(8), intent(IN)  :: x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4
    REAL(8), intent(OUT) :: nx, ny, nz, ox, oy, oz
    REAL(8) :: maxdx, maxdy, maxdz
    REAL(8) :: x1m, y1m, z1m, x2m, y2m, z2m, x3m, y3m, z3m, x4m, y4m, z4m
    REAL(8) :: M11, M12, M13, M21, M22, M23, M31, M32, M33, M41, M42, M43
    REAL(8) :: MT11, MT12, MT13, MT14, MT21, MT22, MT23, MT24, &
            MT31, MT32, MT33, MT34
    REAL(8) :: MTM11, MTM12, MTM13, MTM21, MTM22, MTM23, MTM31, MTM32, MTM33
    REAL(8) :: MTMI11, MTMI12, MTMI13, MTMI21, MTMI22, MTMI23, &
            MTMI31, MTMI32, MTMI33
    REAL(8) :: S11, S12, S13, S14, S21, S22, S23, S24, S31, S32, S33, S34
    REAL(8) :: norm, a, b, c, nxm, nym, nzm, oxm, oym, ozm, cprodx, cprody, cprodz
    integer :: permute

    ! Permute the coordinates to minimize the expected slopes
    maxdx = abs(max(x1, x2, x3, x4) - min(x1, x2, x3, x4))
    maxdy = abs(max(y1, y2, y3, y4) - min(y1, y2, y3, y4))
    maxdz = abs(max(z1, z2, z3, z4) - min(z1, z2, z3, z4))
    if (maxdz <= maxdx .and. maxdz <= maxdy) then
        permute = 0
    else if (maxdy < maxdx .and. maxdy < maxdz) then
        permute = 1
    else
        permute = 2
    end if

    ! cancel permutation
    permute = 0

    if (permute == 0) then
        x1m = x1; x2m = x2; x3m = x3; x4m = x4
        y1m = y1; y2m = y2; y3m = y3; y4m = y4
        z1m = z1; z2m = z2; z3m = z3; z4m = z4
    else if (permute == 1) then
        x1m = y1; x2m = y2; x3m = y3; x4m = y4
        y1m = z1; y2m = z2; y3m = z3; y4m = z4
        z1m = x1; z2m = x2; z3m = x3; z4m = x4
    else if (permute == 2) then
        x1m = z1; x2m = z2; x3m = z3; x4m = z4
        y1m = x1; y2m = x2; y3m = x3; y4m = x4
        z1m = y1; z2m = y2; z3m = y3; z4m = y4
    end if
    
    ! Populate the "M" matrix for the least-squares projection
    M11 = x1m; M12 = y1m; M13 = 1.0
    M21 = x2m; M22 = y2m; M23 = 1.0
    M31 = x3m; M32 = y3m; M33 = 1.0
    M41 = x4m; M42 = y4m; M43 = 1.0

    ! M-transpose
    MT11 = M11; MT12 = M21; MT13 = M31; MT14 = M41
    MT21 = M12; MT22 = M22; MT23 = M32; MT24 = M42
    MT31 = M13; MT32 = M23; MT33 = M33; MT34 = M43

    ! Product of M-transpose with M
    MTM11 = MT11*M11 + MT12*M21 + MT13*M31 + MT14*M41
    MTM12 = MT11*M12 + MT12*M22 + MT13*M32 + MT14*M42
    MTM13 = MT11*M13 + MT12*M23 + MT13*M33 + MT14*M43
    MTM21 = MT21*M11 + MT22*M21 + MT23*M31 + MT24*M41
    MTM22 = MT21*M12 + MT22*M22 + MT23*M32 + MT24*M42
    MTM23 = MT21*M13 + MT22*M23 + MT23*M33 + MT24*M43
    MTM31 = MT31*M11 + MT32*M21 + MT33*M31 + MT34*M41
    MTM32 = MT31*M12 + MT32*M22 + MT33*M32 + MT34*M42
    MTM33 = MT31*M13 + MT32*M23 + MT33*M33 + MT34*M43

    ! Inverse of MTM
    call inverse_3x3(MTM11, MTM12, MTM13, MTM21, MTM22, MTM23, &
                     MTM31, MTM32, MTM33, MTMI11, MTMI12, MTMI13, &
                     MTMI21, MTMI22, MTMI23, MTMI31, MTMI32, MTMI33)

    ! "Solution" matrix: inverse(M-transpose * M) * M-transpose
    S11 = MTMI11*MT11 + MTMI12*MT21 + MTMI13*MT31
    S12 = MTMI11*MT12 + MTMI12*MT22 + MTMI13*MT32
    S13 = MTMI11*MT13 + MTMI12*MT23 + MTMI13*MT33
    S14 = MTMI11*MT14 + MTMI12*MT24 + MTMI13*MT34
    S21 = MTMI21*MT11 + MTMI22*MT21 + MTMI23*MT31
    S22 = MTMI21*MT12 + MTMI22*MT22 + MTMI23*MT32
    S23 = MTMI21*MT13 + MTMI22*MT23 + MTMI23*MT33
    S24 = MTMI21*MT14 + MTMI22*MT24 + MTMI23*MT34
    S31 = MTMI31*MT11 + MTMI32*MT21 + MTMI33*MT31
    S32 = MTMI31*MT12 + MTMI32*MT22 + MTMI33*MT32
    S33 = MTMI31*MT13 + MTMI32*MT23 + MTMI33*MT33
    S34 = MTMI31*MT14 + MTMI32*MT24 + MTMI33*MT34

    ! Best-fit coefficients: zm = a*xm + b*ym + c
    a = S11*z1m + S12*z2m + S13*z3m + S14*z4m
    b = S21*z1m + S22*z2m + S23*z3m + S24*z4m
    c = S31*z1m + S32*z2m + S33*z3m + S34*z4m

    ! Normal vector (in permuted coordinates)
    norm = sqrt(a**2 + b**2 + 1.0)
    nxm = -a/norm
    nym = -b/norm
    nzm = 1.0/norm

    ! Reference point (in permuted coordinates)
    oxm = x1m
    oym = y1m
    ozm = a*oxm + b*oym + c

    ! De-permute if necessary
    if (permute == 0) then
        nx = nxm; ox = oxm
        ny = nym; oy = oym
        nz = nzm; oz = ozm
    else if (permute == 1) then
        nx = nzm; ox = ozm
        ny = nxm; oy = oxm
        nz = oym; oz = oym
    else if (permute == 2) then
        nx = nym; ox = oym
        ny = nzm; oy = ozm
        nz = nxm; oz = oxm
    else
    end if

    ! Reverse the normal vector depending on the handedness of the input points
    call cross_prod((x2-x1), (y2-y1), (z2-z1), (x3-x1), (y3-y1), (z3-z1), &
                    cprodx, cprody, cprodz)
    if (cprodx*nx + cprody*ny + cprodz*nz < 0) then 
        nx = -nx
        ny = -ny
        nz = -nz
    end if

end subroutine plane_4point
!-------------------------------------------------------------------------------
! plane_6point(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, 
!              nx, ny, nz, ox, oy, oz)
!
! Calculates the normal vector and a reference point for a plane that fits
! best (in a least-squares sense) to a set of four input points. 
!
! Input parameters:
!     ingeger :: np                                -> number of query points
!     REAL(8), allocatable :: xq(:), yq(:), zq(:)  -> x, y, and z coordinates of
!                                                     the query points
!
! Output parameters:
!     REAL(8) :: nx, ny, nz   -> x, y, and z components of the unit normal
!     REAL(8) :: ox, oy, oz   -> x, y, and z coordinates of a reference point
!-------------------------------------------------------------------------------
subroutine plane_Npoint(np, xq, yq, zq, nx, ny, nz, ox, oy, oz, a, b, c)

    implicit none

    integer, intent(IN)  :: np
    REAL(8), allocatable, intent(IN)  :: xq(:), yq(:), zq(:)
    REAL(8), intent(OUT) :: nx, ny, nz, ox, oy, oz
    REAL(8) :: maxdx, maxdy, maxdz
    REAL(8), allocatable :: xm(:), ym(:), zm(:)
    REAL(8), allocatable :: M(:,:), MT(:,:)
    REAL(8) :: MTM11, MTM12, MTM13, MTM21, MTM22, MTM23, MTM31, MTM32, MTM33
    REAL(8) :: MTMI11, MTMI12, MTMI13, MTMI21, MTMI22, MTMI23, &
            MTMI31, MTMI32, MTMI33
    REAL(8), allocatable :: S(:,:)
    REAL(8) :: norm, a, b, c, nxm, nym, nzm, oxm, oym, ozm, cprodx, cprody, cprodz
    integer :: permute, i

    allocate( xm(np), ym(np), zm(np), M(np, 3), MT(3, np), S(3, np) )

    ! Permute the coordinates to minimize the expected slopes
    maxdx = abs(maxval(xq) - minval(xq))
    maxdy = abs(maxval(yq) - minval(yq))
    maxdz = abs(maxval(zq) - minval(zq))
    if (maxdz <= maxdx .and. maxdz <= maxdy) then
        permute = 0
    else if (maxdy < maxdx .and. maxdy < maxdz) then
        permute = 1
    else
        permute = 2
    end if

    ! uncomment to cancel permutation
    permute = 0

    do i = 1, np
        if (permute == 0) then
            xm(i) = xq(i)
            ym(i) = yq(i)
            zm(i) = zq(i)
        else if (permute == 1) then
            xm(i) = yq(i)
            ym(i) = zq(i)
            zm(i) = xq(i)
        else if (permute == 2) then
            xm(i) = zq(i)
            ym(i) = xq(i)
            zm(i) = yq(i)
        end if
    end do
    
    ! Populate the "M" matrix and its transpose for the least-squares projection
    do i = 1, np
        M(i,1) = xq(i); MT(1,i) = xq(i)
        M(i,2) = yq(i); MT(2,i) = yq(i)
        M(i,3) = 1.0;   MT(3,i) = 1.0
    end do

    ! Product of M-transpose with M
    MTM11 = 0.0; MTM12 = 0.0; MTM13 = 0.0
    MTM21 = 0.0; MTM22 = 0.0; MTM23 = 0.0
    MTM31 = 0.0; MTM32 = 0.0; MTM33 = 0.0
    do i = 1, np
        MTM11 = MTM11 + MT(1,i)*M(i,1)
        MTM12 = MTM12 + MT(1,i)*M(i,2)
        MTM13 = MTM13 + MT(1,i)*M(i,3)
        MTM21 = MTM21 + MT(2,i)*M(i,1)
        MTM22 = MTM22 + MT(2,i)*M(i,2)
        MTM23 = MTM23 + MT(2,i)*M(i,3)
        MTM31 = MTM31 + MT(3,i)*M(i,1)
        MTM32 = MTM32 + MT(3,i)*M(i,2)
        MTM33 = MTM33 + MT(3,i)*M(i,3)
    end do

    ! Inverse of MTM
    call inverse_3x3(MTM11, MTM12, MTM13, MTM21, MTM22, MTM23, &
                     MTM31, MTM32, MTM33, MTMI11, MTMI12, MTMI13, &
                     MTMI21, MTMI22, MTMI23, MTMI31, MTMI32, MTMI33)

    ! "Solution" matrix: inverse(M-transpose * M) * M-transpose
    do i = 1, np
        S(1,i) = MTMI11*MT(1,i) + MTMI12*MT(2,i) + MTMI13*MT(3,i)
        S(2,i) = MTMI21*MT(1,i) + MTMI22*MT(2,i) + MTMI23*MT(3,i)
        S(3,i) = MTMI31*MT(1,i) + MTMI32*MT(2,i) + MTMI33*MT(3,i)
    end do

    ! Best-fit coefficients: zm = a*xm + b*ym + c
    a = 0.0
    b = 0.0
    c = 0.0
    do i = 1, np
        a = a + S(1,i)*zm(i) 
        b = b + S(2,i)*zm(i) 
        c = c + S(3,i)*zm(i)
    end do

    ! Normal vector (in permuted coordinates)
    norm = sqrt(a**2 + b**2 + 1.0)
    nxm = -a/norm
    nym = -b/norm
    nzm = 1.0/norm

    ! Reference point (in permuted coordinates)
    oxm = xm(1)
    oym = ym(1)
    ozm = a*oxm + b*oym + c

    ! De-permute if necessary
    if (permute == 0) then
        nx = nxm; ox = oxm
        ny = nym; oy = oym
        nz = nzm; oz = ozm
    else if (permute == 1) then
        nx = nzm; ox = ozm
        ny = nxm; oy = oxm
        nz = oym; oz = oym
    else if (permute == 2) then
        nx = nym; ox = oym
        ny = nzm; oy = ozm
        nz = nxm; oz = oxm
    else
    end if

    deallocate( xm, ym, zm, M, MT, S )

end subroutine plane_Npoint


!-------------------------------------------------------------------------------
! plane_elev(ox, oy, oz, nx, ny, nz)
!
! Calculates the elevation of a point above a plane.
!
! Input parameters:
!     REAL(8) :: x, y, z    -> x, y, z coords of the query/input point
!     REAL(8) :: ox, oy, oz -> x, y, z coords of a reference point on the plane
!     REAL(8) :: nx, ny, nz -> x, y, z components of the plane normal vector
!
! Returns:
!     REAL(8) :: plane_elev -> elevation of the input point above the plane
!-------------------------------------------------------------------------------
REAL(8) function plane_elev(x, y, z, ox, oy, oz, nx, ny, nz)

    implicit none

    REAL(8), intent(IN)  :: x, y, z, ox, oy, oz, nx, ny, nz
    
    plane_elev = (x - ox) * nx + (y - oy) * ny + (z - oz) * nz
    
end function plane_elev

!-------------------------------------------------------------------------------
! plane_intersect(ox1, oy1, oz1, nx1, ny1, nz1, &
!                 ox2, oy2, oz2, nx2, ny2, nz2, &
!                 ox3, oy3, oz3, nx3, ny3, nz3, &
!                 xi, yi, zi)
!
! Calculates the point of intersection between three planes. The approach is
! to solve a matrix equation of the form Ax = b, where:
!     A is a matrix in which each row is a normal vector
!     x contains the components of the intersection point
!     b vector in which each component is the dot product of a normal vector
!       with its corresponding intersection point.
!
! Equivalently, solves the system of 3 equations of the form 
!     (I - On) dot Nn = 0,
! where I is the intersection point, On is the origin point of plane n, and 
! Nn is the normal vector of plane n.
!
! Note: normal vectors are assumed to be linearly independent from one another;
! i.e., a solution is assumed to exist.
!
! Input parameters:
!     REAL(8) :: ox1, oy1, oz1 -> x, y, and z coords, ref point on plane 1
!     REAL(8) :: nx1, ny1, nz1 -> x, y, and z components, norm vector to plane 1
!     REAL(8) :: ox2, oy2, oz2 -> x, y, and z coords, ref point on plane 2
!     REAL(8) :: nx2, ny2, nz2 -> x, y, and z components, norm vector to plane 2
!     REAL(8) :: ox3, oy3, oz3 -> x, y, and z coords, ref point on plane 3
!     REAL(8) :: nx3, ny3, nz3 -> x, y, and z components, norm vector to plane 3
!
! Output parameters:
!     REAL(8) :: xi, yi, zi      -> x, y, and z coords, intersection point
!-------------------------------------------------------------------------------
subroutine plane_intersect(ox1, oy1, oz1, nx1, ny1, nz1, &
                           ox2, oy2, oz2, nx2, ny2, nz2, &
                           ox3, oy3, oz3, nx3, ny3, nz3, &
                           xi, yi, zi)

    implicit none

    REAL(8), intent(IN)  :: ox1, oy1, oz1, nx1, ny1, nz1
    REAL(8), intent(IN)  :: ox2, oy2, oz2, nx2, ny2, nz2
    REAL(8), intent(IN)  :: ox3, oy3, oz3, nx3, ny3, nz3
    REAL(8), intent(OUT) :: xi, yi, zi
    REAL(8) :: ni11, ni12, ni13, ni21, ni22, ni23, ni31, ni32, ni33
    REAL(8) :: b1, b2, b3

    ! The right-hand side of the matrix equation (see doc string)
    b1 = nx1*ox1 + ny1*oy1 + nz1*oz1
    b2 = nx2*ox2 + ny2*oy2 + nz2*oz2
    b3 = nx3*ox3 + ny3*oy3 + nz3*oz3

    ! The inverse of the matrix of normal vectors
    call inverse_3x3(nx1,  ny1,  nz1,  nx2,  ny2,  nz2,  nx3,  ny3,  nz3, &
                     ni11, ni12, ni13, ni21, ni22, ni23, ni31, ni32, ni33)

    ! Solution
    xi = ni11*b1 + ni12*b2 + ni13*b3
    yi = ni21*b1 + ni22*b2 + ni23*b3
    zi = ni31*b1 + ni32*b2 + ni33*b3

end subroutine plane_intersect
    

end module magnet_set_calc

