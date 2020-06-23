
!title (bfield) ! Computes magnetic field.

!latex \briefly{Computes magnetic field given coil geometry.}

!latex \calledby{\link{bnormal}}
!latex \calls{}

!latex \tableofcontents

!latex \subsection{magnetic field}
!latex \bi
!latex \item The magnetic field of filamentary coils is calculated bt Biot-Savart Law, involving a line integral.
!latex J. Hanson and S. Hirshman had a better representation for straight segments to avoid unnecessary sigularities
!latex and improve numerical error at points neary the coil.
!latex \item But currently, we use the normal expression of Biot-Savart Law and derivatives of B with repsect to 
!latex x, y, z is also calculated.
!latex \item Later, error analysis and comparison to Hanson's method should be carried out. 
!latex \ei

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
subroutine bfield(icoil, x, y, z, tBx, tBy, tBz)
   !------------------------------------------------------------------------------------------------------ 
   ! DATE:  06/15/2016; 03/26/2017
   ! calculate the magnetic field of icoil using manually discretized coils. 
   ! Biot-Savart constant and currents are not included for later simplication. 
   ! Be careful if coils have different resolutions.
   !------------------------------------------------------------------------------------------------------   
     use globals, only: dp, coil, surf, Ncoils, Nteta, Nzeta, cosnfp, sinnfp, &
                        zero, myid, ounit, Nfp, pi2, half, two, one, bsconstant, MPI_COMM_FOCUS
     use mpi
     implicit none
   
     INTEGER, intent(in ) :: icoil
     REAL,    intent(in ) :: x, y, z
     REAL   , intent(out) :: tBx, tBy, tBz
   
   !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
     INTEGER              :: ierr, astat, kseg, ns, i
     REAL                 :: Rix, Riy, Riz, Ri, Rfx, Rfy, Rfz, Rf, lx, ly, lz, ll, Rfac

     ns = coil(icoil)%NS
     tBx = zero; tBy = zero; tBz = zero
     do i = 0, ns-1
        Rix = x - coil(icoil)%xx(i); Rfx = x - coil(icoil)%xx(i+1); lx = coil(icoil)%xx(i+1) - coil(icoil)%xx(i)
        Riy = y - coil(icoil)%yy(i); Rfy = y - coil(icoil)%yy(i+1); ly = coil(icoil)%yy(i+1) - coil(icoil)%yy(i)
        Riz = z - coil(icoil)%zz(i); Rfz = z - coil(icoil)%zz(i+1); lz = coil(icoil)%zz(i+1) - coil(icoil)%zz(i)
        Ri = sqrt(Rix*Rix+Riy*Riy+Riz*Riz)
        Rf = sqrt(Rfx*Rfx+Rfy*Rfy+Rfz*Rfz)
        ll = sqrt(lx*lx+ly*ly+lz*lz)
        Rfac = 2*(Ri+Rf)/(Ri*Rf)/((Ri+Rf)**2 - ll**2)
        tBx = tBx + Rfac*(ly*Riz - lz*Riy)
        tBy = tBy + Rfac*(lz*Rix - lx*Riz)
        tBz = tBz + Rfac*(lx*Riy - ly*Rix)
     enddo
     tBx = tBx*bsconstant*coil(icoil)%I
     tBy = tBy*bsconstant*coil(icoil)%I
     tBz = tBz*bsconstant*coil(icoil)%I

     return
end subroutine bfield
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine coils_bfield(B,x,y,z)
  use globals, only: dp, coil, surf, Ncoils, Nteta, Nzeta, &
       zero, myid, ounit, Nfp, bsconstant, one, two, ncpu, &
       master, nworker, myworkid, MPI_COMM_MASTERS, MPI_COMM_MYWORLD, MPI_COMM_WORKERS, MPI_COMM_FOCUS
  use mpi
  implicit none

  REAL  , intent( in)     :: x, y, z
  REAL  , intent(inout)   :: B(3)
  !INTEGER, INTENT(in)     :: comm ! MPI communicator

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  INTEGER              :: ierr, astat
  REAL                 :: Bx, By, Bz
  INTEGER              :: icoil

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  call MPI_BARRIER(MPI_COMM_MYWORLD, ierr ) ! wait all cpus;

  B = zero
  do icoil = 1, Ncoils
     if ( myworkid /= modulo(icoil-1, nworker) ) cycle ! MPI
     ! Bx = zero; By = zero; Bz = zero
     call bfield( icoil, x, y, z, Bx, By, Bz )
     B(1) = B(1) + Bx
     B(2) = B(2) + By
     B(3) = B(3) + Bz
  enddo

  call MPI_ALLREDUCE(MPI_IN_PLACE, B, 3, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_MYWORLD, ierr )

  return

end subroutine coils_bfield

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
