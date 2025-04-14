subroutine rdhdf5(filename, index)
   ! read surface data from the HDF5 file
   ! data includes nfp, issymmetric, xsurf, ysurf, zsurf 
   !      nx, ny, nz, nn, plasma_Bn
   use globals, only: dp, zero, half, pi2, myid, ounit, runit, IsQuiet, IsSymmetric, &
                      Nteta, Nzeta, surf, discretefactor, Nfp, plasma, symmetry, &
                      tflux_sign, cosnfp, sinnfp, MPI_COMM_FOCUS, surf_Nfp
   use mpi
   use hdf5
   implicit none

   CHARACTER*100, INTENT(IN) :: filename
   INTEGER, INTENT(IN) :: index

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

   INTEGER :: iosta, astat, ierr, hdfier, dims(1:2), ip
   integer(hid_t) :: file_id, dset_id, dspace_id, mem_space_id
   integer(hsize_t) :: onedims(1:1), twodims(1:2), threedims(1:3), maxdims(1:1)
   REAL :: dz

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  !-------------output for check-------------------------------------------------------------------------
   if (myid == 0 .and. IsQuiet <= 0) then
      write (ounit, *) "-----------Reading surface-----------------------------------"
      write (ounit, '("surface : Surface data will be read from ", A)') trim(filename)
   endif

   ! read the data
   if (myid == 0) then
        call h5open_f(ierr) ! initialize
        call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, ierr) ! open file
        ! read and assign the data
        HREADIS(Nfp, surf(index)%Nfp) ! number of toroidal periodicity
        HREADIS(IsSymmetric, IsSymmetric) ! stellarator symmetry option
        HREADRS(disfac, discretefactor) ! read discret_factor for surface integration
        ! the array should be (Ntheta, Nzeta)
        HREADRA(xsurf, surf(index)%xx) ! x coordinates of the surface elements
        HREADRA(ysurf, surf(index)%yy) ! y coordinates of the surface elements
        HREADRA(zsurf, surf(index)%zz) ! z coordinates of the surface elements
        HREADRA(nx, surf(index)%nx) ! x coordinates of the surface unit normal
        HREADRA(ny, surf(index)%ny) ! y coordinates of the surface unit normal
        HREADRA(nz, surf(index)%nz) ! z coordinates of the surface unit normal
        HREADRA(nn, surf(index)%ds) ! the surface jacobian (normal vector magnitude)
        HREADRA(plas_Bn, surf(index)%pb) ! plasma Bn information
        ! check the dimensions
        dims = SHAPE(surf(index)%xx)
        Nteta = dims(1)
        Nzeta = dims(2)
        if (IsQuiet <= 0) then
            write (ounit, '("surface : Nfp = " I06 " ; IsSymmetric = " I06)') surf(index)%Nfp, IsSymmetric
            write (ounit, '("surface : Surface resolution: Nteta = "I6", Nzeta = "I6" .")')  Nteta, Nzeta
            write (ounit, '("surface : discretization factor disfac = "ES13.5" .")')  discretefactor            
        endif 
    endif

    ! other CPUs allocate data first
    IlBCAST(Nteta, 1, 0)
    IlBCAST(Nzeta, 1, 0)   
    if (myid /= 0) then 
        SALLOCATE(surf(index)%xx, (0:Nteta - 1, 0:Nzeta - 1), zero) !x coordinates;
        SALLOCATE(surf(index)%yy, (0:Nteta - 1, 0:Nzeta - 1), zero) !y coordinates
        SALLOCATE(surf(index)%zz, (0:Nteta - 1, 0:Nzeta - 1), zero) !z coordinates
        SALLOCATE(surf(index)%nx, (0:Nteta - 1, 0:Nzeta - 1), zero) !unit nx;
        SALLOCATE(surf(index)%ny, (0:Nteta - 1, 0:Nzeta - 1), zero) !unit ny;
        SALLOCATE(surf(index)%nz, (0:Nteta - 1, 0:Nzeta - 1), zero) !unit nz;
        SALLOCATE(surf(index)%ds, (0:Nteta - 1, 0:Nzeta - 1), zero) !jacobian;
        SALLOCATE(surf(index)%pb, (0:Nteta - 1, 0:Nzeta - 1), zero) !target Bn;
    endif 

    ! broadcast everything
    IlBCAST(surf(index)%Nfp, 1, 0)
    IlBCAST(IsSymmetric, 1, 0)
    RlBCAST(discretefactor, 1, 0)
    RlBCAST(surf(index)%xx, Nteta*Nzeta, 0)
    RlBCAST(surf(index)%yy, Nteta*Nzeta, 0)
    RlBCAST(surf(index)%zz, Nteta*Nzeta, 0)
    RlBCAST(surf(index)%nx, Nteta*Nzeta, 0)
    RlBCAST(surf(index)%ny, Nteta*Nzeta, 0)
    RlBCAST(surf(index)%nz, Nteta*Nzeta, 0)
    RlBCAST(surf(index)%ds, Nteta*Nzeta, 0)
    RlBCAST(surf(index)%pb, Nteta*Nzeta, 0)    
    ! dealing with Nfp
    Nfp = surf(plasma)%Nfp
    surf_Nfp = Nfp ! local surface Nfp
    select case (IsSymmetric)
    case (0)
        surf_Nfp = 1             ! reset Nfp to 1;
        symmetry = 0
    case (1)                    ! plasma and coil periodicity enabled;
        symmetry = 0
    case (2)                    ! stellarator symmetry enforced;
        symmetry = 1
    end select
    surf(index)%Nteta = Nteta
    surf(index)%Nzeta = Nzeta*surf_Nfp*2**symmetry ! the total number from [0, 2pi]
    ! calculate the area and volumn
    surf(index)%area = SUM(surf(index)%ds) * discretefactor * surf_Nfp * 2**symmetry
    surf(index)%vol = SUM((surf(index)%xx*surf(index)%nx+ surf(index)%yy*surf(index)%ny &
                         + surf(index)%zz*surf(index)%nz)*surf(index)%ds) / 3 &
                         * discretefactor * surf_Nfp * 2**symmetry
   if (myid == 0 .and. IsQuiet <= 0) then
      write (ounit, '(8X": Enclosed total surface volume ="ES12.5" m^3 ; area ="ES12.5" m^2." )') &
         surf(index)%vol, surf(index)%area
   endif
   ! allocate cosnfp and sinnfp
   if (index == plasma) then
      SALLOCATE(cosnfp, (1:Nfp), zero)
      SALLOCATE(sinnfp, (1:Nfp), zero)
      do ip = 1, Nfp
         cosnfp(ip) = cos((ip - 1)*pi2/Nfp)
         sinnfp(ip) = sin((ip - 1)*pi2/Nfp)
      enddo
    endif
   ! check theta direction for the plasma surface and determine the toroidal flux sign
   if (index == plasma) then
      dz = surf(plasma)%zz(1, 0) - surf(plasma)%zz(0, 0)
      if (dz > 0) then
         ! counter-clockwise
         if (myid == 0) write (ounit, '(8X": The theta angle used is counter-clockwise.")')
         tflux_sign = -1
      else
         ! clockwise
         if (myid == 0) write (ounit, '(8X": The theta angle used is clockwise.")')
         tflux_sign = 1
      endif
   endif
   ! some additional quantities; not sure if need read from files
   SALLOCATE(surf(index)%xt, (0:Nteta - 1, 0:Nzeta - 1), zero) !dx/dtheta;
   SALLOCATE(surf(index)%yt, (0:Nteta - 1, 0:Nzeta - 1), zero) !dy/dtheta;
   SALLOCATE(surf(index)%zt, (0:Nteta - 1, 0:Nzeta - 1), zero) !dz/dtheta;
   SALLOCATE(surf(index)%xp, (0:Nteta - 1, 0:Nzeta - 1), zero) !dx/dzeta;
   SALLOCATE(surf(index)%yp, (0:Nteta - 1, 0:Nzeta - 1), zero) !dy/dzeta;
   SALLOCATE(surf(index)%zp, (0:Nteta - 1, 0:Nzeta - 1), zero) !dz/dzeta;
end subroutine rdhdf5
