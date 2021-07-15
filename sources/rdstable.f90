!title (boundary) ! The plasma O and X points are read from file 
!latex \briefly{A Fourier representation for the plasma O and X points are read from file }

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine rdstable(filename, index)
  use globals, only : dp, zero, half, pi2, myid, ounit, runit, IsQuiet, period_Nseg, resbn_m, &
          gsurf, MPI_COMM_FOCUS
  use mpi
  implicit none

  CHARACTER*100, INTENT(IN) :: filename
  INTEGER, INTENT(IN) :: index
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  INTEGER :: iosta, astat, ierr, NF_stable, NF_axis, i, Nseg_stable

  ! read the header
  if( myid == 0 ) then
     open(runit, file=trim(filename), status='old', action='read')
     read(runit,*) !empty line
     read(runit,*) gsurf(index)%NF_stable
  endif

  !Broadcast the values
  IlBCAST( gsurf(index)%NF_stable, 1, 0 )
  FATAL( stable, gsurf(index)%NF_stable <= 0, invalid )
  NF_stable = gsurf(index)%NF_stable

  !allocate arrays
  SALLOCATE( gsurf(index)%on,       (1:NF_stable), 0   )
  SALLOCATE( gsurf(index)%osnc,     (1:NF_stable), 0.0 )
  SALLOCATE( gsurf(index)%osns,     (1:NF_stable), 0.0 )
  SALLOCATE( gsurf(index)%othetanc, (1:NF_stable), 0.0 )
  SALLOCATE( gsurf(index)%othetans, (1:NF_stable), 0.0 )
  SALLOCATE( gsurf(index)%xn,       (1:NF_stable), 0   )
  SALLOCATE( gsurf(index)%xsnc,     (1:NF_stable), 0.0 )
  SALLOCATE( gsurf(index)%xsns,     (1:NF_stable), 0.0 )
  SALLOCATE( gsurf(index)%xthetanc, (1:NF_stable), 0.0 )
  SALLOCATE( gsurf(index)%xthetans, (1:NF_stable), 0.0 )

  if( myid == 0 ) then
     read(runit,*) !empty line
     read(runit,*) !empty line
     do i = 1, NF_stable
        read(runit,*) gsurf(index)%on(i), gsurf(index)%osnc(i), gsurf(index)%osns(i), &
                gsurf(index)%othetanc(i), gsurf(index)%othetans(i)
     enddo
  endif  
  IlBCAST( gsurf(index)%on(1:NF_stable), NF_stable, 0 )
  RlBCAST( gsurf(index)%osnc(1:NF_stable), NF_stable, 0 )
  RlBCAST( gsurf(index)%osns(1:NF_stable), NF_stable, 0 )
  RlBCAST( gsurf(index)%othetanc(1:NF_stable), NF_stable, 0 )
  RlBCAST( gsurf(index)%othetans(1:NF_stable), NF_stable, 0 )

  if( myid == 0 ) then
     read(runit,*) !empty line
     read(runit,*) !empty line
     do i = 1, NF_stable
        read(runit,*) gsurf(index)%xn(i), gsurf(index)%xsnc(i), gsurf(index)%xsns(i), &
                gsurf(index)%xthetanc(i), gsurf(index)%xthetans(i)
     enddo
  endif
  IlBCAST( gsurf(index)%xn(1:NF_stable), NF_stable, 0 )
  RlBCAST( gsurf(index)%xsnc(1:NF_stable), NF_stable, 0 )
  RlBCAST( gsurf(index)%xsns(1:NF_stable), NF_stable, 0 )
  RlBCAST( gsurf(index)%xthetanc(1:NF_stable), NF_stable, 0 )
  RlBCAST( gsurf(index)%xthetans(1:NF_stable), NF_stable, 0 )

  if( myid == 0 ) then
     read(runit,*) ! empty line
     read(runit,*) gsurf(index)%NF_axis
  endif
  IlBCAST( gsurf(index)%NF_axis, 1, 0 )
  FATAL( stable, gsurf(index)%NF_axis <= 0, invalid ) 

  NF_axis = gsurf(index)%NF_axis
  SALLOCATE( gsurf(index)%axisn,   (1:NF_axis), 0   )
  SALLOCATE( gsurf(index)%axisrnc, (1:NF_axis), 0.0 )
  SALLOCATE( gsurf(index)%axiszns, (1:NF_axis), 0.0 )

  if( myid == 0 ) then
     read(runit,*) !empty line
     read(runit,*) !empty line
     do i = 1, NF_axis
        read(runit,*) gsurf(index)%axisn(i), gsurf(index)%axisrnc(i), gsurf(index)%axiszns(i)
     enddo
  endif
  IlBCAST( gsurf(index)%axisn(1:NF_axis), NF_axis, 0 )
  RlBCAST( gsurf(index)%axisrnc(1:NF_axis), NF_axis, 0 )
  RlBCAST( gsurf(index)%axiszns(1:NF_axis), NF_axis, 0 )

  if( myid == 0 ) close(runit,iostat=iosta)
  IlBCAST( iosta, 1, 0 )
  FATAL( stable, iosta.ne.0, error closing the stablefile )

  !-------------output for check-------------------------------------------------------------------------
  if( myid == 0 .and. IsQuiet <= 0) then
     write(ounit, *) "-----------Reading surface-----------------------------------"
     write(ounit, '("ghost : Stable field lines read from ", A)') trim(filename)
     write(ounit, '(8X": NF_stable = " I06 " ; NF_axis = " I06 " ;" )') gsurf(index)%NF_stable, gsurf(index)%NF_axis
  endif

  if( myid == 0 .and. IsQuiet <= -2) then ! very detailed output;
     write(ounit,'("        : " 10x " : on ="10i13   )')       gsurf(index)%on(1:NF_stable)
     write(ounit,'("        : " 10x " : osnc ="10es13.5)')     gsurf(index)%osnc(1:NF_stable)
     write(ounit,'("        : " 10x " : osns ="10es13.5)')     gsurf(index)%osns(1:NF_stable)
     write(ounit,'("        : " 10x " : othetanc ="10es13.5)') gsurf(index)%othetanc(1:NF_stable)
     write(ounit,'("        : " 10x " : othetans ="10es13.5)') gsurf(index)%othetans(1:NF_stable)
     write(ounit,'("        : " 10x " : xn ="10i13   )')       gsurf(index)%xn(1:NF_stable)
     write(ounit,'("        : " 10x " : xsnc ="10es13.5)')     gsurf(index)%xsnc(1:NF_stable)
     write(ounit,'("        : " 10x " : xsns ="10es13.5)')     gsurf(index)%xsns(1:NF_stable)
     write(ounit,'("        : " 10x " : xthetanc ="10es13.5)') gsurf(index)%xthetanc(1:NF_stable)
     write(ounit,'("        : " 10x " : xthetans ="10es13.5)') gsurf(index)%xthetans(1:NF_stable)
     write(ounit,'("        : " 10x " : axisn ="10i13   )')    gsurf(index)%axisn(1:NF_axis)
     write(ounit,'("        : " 10x " : axisrnc ="10es13.5)')  gsurf(index)%axisrnc(1:NF_axis)
     write(ounit,'("        : " 10x " : axiszns ="10es13.5)')  gsurf(index)%axiszns(1:NF_axis)
  endif

  ! Probably delete everything below
  ! Probably allocate all the gsurf variables here

  gsurf(index)%Nseg_stable = period_Nseg*resbn_m + 1 ! At some point resbn_m needs to be in input file 
  Nseg_stable = gsurf(index)%Nseg_stable

  SALLOCATE( gsurf(index)%zeta, (1:Nseg_stable), 0.0 )
  SALLOCATE( gsurf(index)%Ra, (1:Nseg_stable), 0.0 )
  SALLOCATE( gsurf(index)%Za, (1:Nseg_stable), 0.0 )
  SALLOCATE( gsurf(index)%os, (1:Nseg_stable), 0.0 )
  SALLOCATE( gsurf(index)%xs, (1:Nseg_stable), 0.0 )
  SALLOCATE( gsurf(index)%otheta, (1:Nseg_stable), 0.0 )
  SALLOCATE( gsurf(index)%xtheta, (1:Nseg_stable), 0.0 )
  SALLOCATE( gsurf(index)%osdot, (1:Nseg_stable), 0.0 )
  SALLOCATE( gsurf(index)%xsdot, (1:Nseg_stable), 0.0 )
  SALLOCATE( gsurf(index)%othetadot, (1:Nseg_stable), 0.0 )
  SALLOCATE( gsurf(index)%xthetadot, (1:Nseg_stable), 0.0 )
  SALLOCATE( gsurf(index)%ox, (1:Nseg_stable), 0.0 )
  SALLOCATE( gsurf(index)%oy, (1:Nseg_stable), 0.0 )
  SALLOCATE( gsurf(index)%oz, (1:Nseg_stable), 0.0 )
  SALLOCATE( gsurf(index)%xx, (1:Nseg_stable), 0.0 )
  SALLOCATE( gsurf(index)%xy, (1:Nseg_stable), 0.0 )
  SALLOCATE( gsurf(index)%xz, (1:Nseg_stable), 0.0 )
  SALLOCATE( gsurf(index)%of, (1:Nseg_stable), 0.0 )
  SALLOCATE( gsurf(index)%og, (1:Nseg_stable), 0.0 )
  SALLOCATE( gsurf(index)%xf, (1:Nseg_stable), 0.0 )
  SALLOCATE( gsurf(index)%xg, (1:Nseg_stable), 0.0 )
  SALLOCATE( gsurf(index)%oxdot, (1:Nseg_stable), 0.0 )
  SALLOCATE( gsurf(index)%oydot, (1:Nseg_stable), 0.0 )
  SALLOCATE( gsurf(index)%ozdot, (1:Nseg_stable), 0.0 )
  SALLOCATE( gsurf(index)%xxdot, (1:Nseg_stable), 0.0 )
  SALLOCATE( gsurf(index)%xydot, (1:Nseg_stable), 0.0 )
  SALLOCATE( gsurf(index)%xzdot, (1:Nseg_stable), 0.0 )
  SALLOCATE( gsurf(index)%obsups, (1:Nseg_stable), 0.0 )
  SALLOCATE( gsurf(index)%obsupzeta, (1:Nseg_stable), 0.0 )
  SALLOCATE( gsurf(index)%obsuptheta, (1:Nseg_stable), 0.0 )
  SALLOCATE( gsurf(index)%xbsups, (1:Nseg_stable), 0.0 )
  SALLOCATE( gsurf(index)%xbsupzeta, (1:Nseg_stable), 0.0 )
  SALLOCATE( gsurf(index)%xbsuptheta, (1:Nseg_stable), 0.0 )


  ! do this once for now
  !if (ghost_use .eq. 1) then
  !   call ghost(1)
  !endif
  gsurf(index)%donee = 0

  do i = 1,gsurf(index)%Nseg_stable
     gsurf(index)%zeta(i) = (i-1)*pi2*resbn_m/(Nseg_stable-1)
  enddo

  gsurf(index)%Ndof_stable = 8*gsurf(index)%NF_stable
  SALLOCATE( gsurf(index)%xdof_stable, (1:gsurf(index)%Ndof_stable), 0.0 )
  SALLOCATE( gsurf(index)%dFdxdof_stable, (1:gsurf(index)%Ndof_stable), 0.0 )
  
  gsurf(index)%xdof_stable(            1:  NF_stable) = gsurf(index)%osnc(1:NF_stable)
  gsurf(index)%xdof_stable(  NF_stable+1:2*NF_stable) = gsurf(index)%osns(1:NF_stable)
  gsurf(index)%xdof_stable(2*NF_stable+1:3*NF_stable) = gsurf(index)%othetanc(1:NF_stable)
  gsurf(index)%xdof_stable(3*NF_stable+1:4*NF_stable) = gsurf(index)%othetans(1:NF_stable)
  gsurf(index)%xdof_stable(4*NF_stable+1:5*NF_stable) = gsurf(index)%xsnc(1:NF_stable)
  gsurf(index)%xdof_stable(5*NF_stable+1:6*NF_stable) = gsurf(index)%xsns(1:NF_stable)
  gsurf(index)%xdof_stable(6*NF_stable+1:7*NF_stable) = gsurf(index)%xthetanc(1:NF_stable)
  gsurf(index)%xdof_stable(7*NF_stable+1:8*NF_stable) = gsurf(index)%xthetans(1:NF_stable)

  gsurf(index)%iter_track = 1
  SALLOCATE( gsurf(index)%xdof_stable_hold, (1:gsurf(index)%Ndof_stable), 0.0 )
  gsurf(index)%xdof_stable_hold(1:gsurf(index)%Ndof_stable) = gsurf(index)%xdof_stable(1:gsurf(index)%Ndof_stable)






  !surf(index)%Nteta = Nteta 
  !surf(index)%Nzeta = Nzeta

  !if (index == plasma) then
  !   Nfp = surf(plasma)%Nfp 
  !   surf_Nfp = Nfp ! local surface Nfp
  !   select case (IsSymmetric)
  !   case ( 0 )
  !      surf_Nfp = 1                    ! reset Nfp to 1;
  !      symmetry = 0
  !   case ( 1 )                    ! plasma and coil periodicity enabled;
  !      symmetry = 0
  !   case ( 2 )                    ! stellarator symmetry enforced;
  !      symmetry = 1     
  !   end select

  !   SALLOCATE( cosnfp, (1:Nfp), zero )
  !   SALLOCATE( sinnfp, (1:Nfp), zero )  
  !   do ip = 1, Nfp
  !      cosnfp(ip) = cos((ip-1)*pi2/Nfp)
  !      sinnfp(ip) = sin((ip-1)*pi2/Nfp)
  !   enddo
  !   ! discretefactor = discretefactor/Nfp
  !   surf(index)%Nzeta = Nzeta * surf_Nfp * 2**symmetry ! the total number from [0, 2pi]
  !   discretefactor = (pi2/surf(plasma)%Nteta) * (pi2/surf(plasma)%Nzeta)
  !endif
  
  !SALLOCATE( surf(index)%xx, (0:Nteta-1,0:Nzeta-1), zero ) !x coordinates;
  !SALLOCATE( surf(index)%yy, (0:Nteta-1,0:Nzeta-1), zero ) !y coordinates
  !SALLOCATE( surf(index)%zz, (0:Nteta-1,0:Nzeta-1), zero ) !z coordinates
  !SALLOCATE( surf(index)%nx, (0:Nteta-1,0:Nzeta-1), zero ) !unit nx;
  !SALLOCATE( surf(index)%ny, (0:Nteta-1,0:Nzeta-1), zero ) !unit ny;
  !SALLOCATE( surf(index)%nz, (0:Nteta-1,0:Nzeta-1), zero ) !unit nz;
  !SALLOCATE( surf(index)%ds, (0:Nteta-1,0:Nzeta-1), zero ) !jacobian;
  !SALLOCATE( surf(index)%xt, (0:Nteta-1,0:Nzeta-1), zero ) !dx/dtheta;
  !SALLOCATE( surf(index)%yt, (0:Nteta-1,0:Nzeta-1), zero ) !dy/dtheta;
  !SALLOCATE( surf(index)%zt, (0:Nteta-1,0:Nzeta-1), zero ) !dz/dtheta;
  !SALLOCATE( surf(index)%pb, (0:Nteta-1,0:Nzeta-1), zero ) !target Bn;
  !SALLOCATE( surf(index)%xp, (0:Nteta-1,0:Nzeta-1), zero ) !dx/dzeta;
  !SALLOCATE( surf(index)%yp, (0:Nteta-1,0:Nzeta-1), zero ) !dy/dzeta;
  !SALLOCATE( surf(index)%zp, (0:Nteta-1,0:Nzeta-1), zero ) !dz/dzeta;
  
  !surf(index)%vol = zero  ! volume enclosed by plasma boundary
  !surf(index)%area = zero ! surface area
 
  ! The center point value was used to discretize grid;
  !do ii = 0, Nteta-1
  !   teta = ( ii + half ) * pi2 / surf(index)%Nteta
  !   do jj = 0, Nzeta-1
  !      zeta = ( jj + half ) * pi2 / surf(index)%Nzeta
  !      RR(0:2) = zero ; ZZ(0:2) = zero
  !      do imn = 1, surf(index)%Nfou
  !         arg = surf(index)%bim(imn) * teta - surf(index)%bin(imn) * zeta
  !         RR(0) =  RR(0) +     surf(index)%Rbc(imn) * cos(arg) + surf(index)%Rbs(imn) * sin(arg)
  !         ZZ(0) =  ZZ(0) +     surf(index)%Zbc(imn) * cos(arg) + surf(index)%Zbs(imn) * sin(arg)
  !         RR(1) =  RR(1) + ( - surf(index)%Rbc(imn) * sin(arg) + surf(index)%Rbs(imn) * cos(arg) ) * surf(index)%bim(imn)
  !         ZZ(1) =  ZZ(1) + ( - surf(index)%Zbc(imn) * sin(arg) + surf(index)%Zbs(imn) * cos(arg) ) * surf(index)%bim(imn)
  !         RR(2) =  RR(2) - ( - surf(index)%Rbc(imn) * sin(arg) + surf(index)%Rbs(imn) * cos(arg) ) * surf(index)%bin(imn)
  !         ZZ(2) =  ZZ(2) - ( - surf(index)%Zbc(imn) * sin(arg) + surf(index)%Zbs(imn) * cos(arg) ) * surf(index)%bin(imn)
  !      enddo ! end of do imn; 30 Oct 15;
  !      szeta = sin(zeta)
  !      czeta = cos(zeta)
  !      xx(1:3) = (/   RR(0) * czeta,   RR(0) * szeta, ZZ(0) /)
  !      xt(1:3) = (/   RR(1) * czeta,   RR(1) * szeta, ZZ(1) /)
  !      xz(1:3) = (/   RR(2) * czeta,   RR(2) * szeta, ZZ(2) /) &
  !              + (/ - RR(0) * szeta,   RR(0) * czeta, zero  /)
  !      ! minus sign for theta counterclockwise direction;
  !      ds(1:3) = -(/ xt(2) * xz(3) - xt(3) * xz(2), & 
  !                    xt(3) * xz(1) - xt(1) * xz(3), &
  !                    xt(1) * xz(2) - xt(2) * xz(1) /)
  !      dd = sqrt( sum( ds(1:3)*ds(1:3) ) )
  !      ! x, y, z coordinates for the surface;
  !      surf(index)%xx(ii,jj) = xx(1)
  !      surf(index)%yy(ii,jj) = xx(2)
  !      surf(index)%zz(ii,jj) = xx(3)
  !      ! dx/dt, dy/dt, dz/dt (dt for d theta)
  !      surf(index)%xt(ii,jj) = xt(1)
  !      surf(index)%yt(ii,jj) = xt(2)
  !      surf(index)%zt(ii,jj) = xt(3)
  !      ! dx/dp, dy/dp, dz/dp (dp for d zeta(phi))
  !      surf(index)%xp(ii,jj) = xz(1)
  !      surf(index)%yp(ii,jj) = xz(2)
  !      surf(index)%zp(ii,jj) = xz(3)
  !      ! surface normal vectors and ds for the jacobian;
  !      surf(index)%nx(ii,jj) = ds(1) / dd
  !      surf(index)%ny(ii,jj) = ds(2) / dd
  !      surf(index)%nz(ii,jj) = ds(3) / dd
  !      surf(index)%ds(ii,jj) =         dd
  !      ! using Gauss theorom; V = \int_S x \cdot n dt dz
  !      surf(index)%vol = surf(index)%vol + surf(index)%xx(ii,jj) * ds(1) &
  !           & + surf(index)%yy(ii,jj) * ds(2) + surf(index)%zz(ii,jj) * ds(3)
  !      ! surface area 
  !      surf(index)%area = surf(index)%area + surf(index)%ds(ii,jj)
  !   enddo ! end of do jj; 14 Apr 16;
  !enddo ! end of do ii; 14 Apr 16;

  !! print volume and area
  !surf(index)%vol  = abs(surf(index)%vol)/3 * (pi2/surf(index)%Nteta) * (pi2/surf(index)%Nzeta)
  !surf(index)%area = abs(surf(index)%area)  * (pi2/surf(index)%Nteta) * (pi2/surf(index)%Nzeta)
  !if (index == plasma) then
  !   surf(index)%vol  = surf(index)%vol  * surf_Nfp * 2**symmetry
  !   surf(index)%area = surf(index)%area * surf_Nfp * 2**symmetry
  !endif    
     
  !if( myid == 0 .and. IsQuiet <= 0) then
  !   write(ounit, '(8X": Enclosed total surface volume ="ES12.5" m^3 ; area ="ES12.5" m^2." )') &
  !        surf(index)%vol, surf(index)%area
  !endif

  ! check theta direction for the plasma surface and determine the toroidal flux sign
  !if (index == plasma) then
  !   dz = surf(plasma)%zz(1,0) - surf(plasma)%zz(0,0)
  !   if (dz > 0) then
  !      ! counter-clockwise
  !      if( myid == 0) write(ounit, '(8X": The theta angle used is counter-clockwise.")')
  !      tflux_sign = -1
  !   else
  !      ! clockwise
  !      if( myid == 0) write(ounit, '(8X": The theta angle used is clockwise.")')
  !      tflux_sign =  1 
  !   endif
  !endif

  !calculate target Bn with input harmonics; 05 Jan 17;
  !if(surf(index)%NBnf >  0) then
  !   do jj = 0, Nzeta-1
  !      zeta = ( jj + half ) * pi2 / surf(index)%Nzeta
  !      do ii = 0, Nteta-1
  !         teta = ( ii + half ) * pi2 / surf(index)%Nteta
  !         do imn = 1, surf(index)%NBnf
  !            arg = surf(index)%Bnim(imn) * teta - surf(index)%Bnin(imn) * zeta
  !            surf(index)%pb(ii,jj) = surf(index)%pb(ii,jj) + surf(index)%Bnc(imn)*cos(arg) + surf(index)%Bns(imn)*sin(arg)
  !         enddo
  !      enddo
  !   enddo
  !endif

  return
  
end subroutine rdstable
