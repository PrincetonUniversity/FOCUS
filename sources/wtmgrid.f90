! write binary mgrid file
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine wtmgrid
  use globals, only : dp, zero, half, pi2, ext, ncpu, myid, ounit, wunit, runit, surf, plasma, &
       sqrtmachprec, master, nmaster, nworker, masterid, color, myworkid, &
       MPI_COMM_MASTERS, MPI_COMM_MYWORLD, MPI_COMM_WORKERS, MPI_COMM_FOCUS
  use mgrid_mod
  implicit none
  include "mpif.h"

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  LOGICAL              :: exist
  INTEGER              :: ierr, astat, iostat, ip, iz, ir, nextcur, icpu
  REAL                 :: RpZ(1:3), R, P, Z, B(1:3), pressure, gap, &
       czeta, szeta, xx, yy, zz, dx, dy, dz, dBx, dBy, dBz
  REAL, allocatable    :: BRZp(:,:,:,:), BRpZ(:,:,:,:)
  CHARACTER(LEN=100)   :: mgrid_name
  CHARACTER(LEN=30)    :: curlabel(1:1)

  do icpu = 1, ncpu
     call MPI_BARRIER( MPI_COMM_FOCUS, ierr )
     if (myid == icpu-1) then                              ! each cpu read the namelist in turn;
        open(runit, file=trim(trim(ext)//".input"), status="old", action='read')
        read(runit, mgrid)
        close(runit)
     endif ! end of if( myid == 0 )
  enddo
  
  mgrid_name = "mgrid.focus_"//trim(ext) ! filename, could be user input
  if (Mfp <= 0) Mfp = surf(plasma)%Nfp   ! overrid to nfp_raw if not specified
  B = zero  ; dx = 1E-4 ; dy = 1E-4 ; dz = 1E-4

  FATAL( wrmgrid, abs(Rmin)+abs(Rmax)<sqrtmachprec, please provide effective dimensions in mgrid nml )

  SALLOCATE( BRZp, (1:3,1:Nr,1:Nz,1:Np), zero )
#ifdef DIV_CHECK
  SALLOCATE( BRpZ, (1:2,1:Nr,1:Nz,1:Np), zero )
#endif

  if( myid.eq.0 ) write( ounit,'("wtmgrid : writing mgrid file at grid of [ "4(ES12.5,2X)" ]",3i6)') Rmin, Rmax, Zmin, Zmax, np, nr, nz

  ! split cores
  color = modulo(myid, np)
  CALL MPI_COMM_SPLIT(MPI_COMM_FOCUS, color, myid, MPI_COMM_MYWORLD, ierr)
  CALL MPI_COMM_RANK(MPI_COMM_MYWORLD, myworkid, ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_MYWORLD, nworker, ierr)

  if (myworkid /= 0) then
     color = MPI_UNDEFINED
     masterid = -1
  else 
     color = 0
  endif
  CALL MPI_COMM_SPLIT(MPI_COMM_FOCUS, color, myid, MPI_COMM_MASTERS, ierr)
  if ( myworkid == 0 ) then
     CALL MPI_COMM_RANK(MPI_COMM_MASTERS, masterid, ierr)
     CALL MPI_COMM_SIZE(MPI_COMM_MASTERS, nmaster, ierr)
  endif
  IlBCAST( nmaster, 1, master )

  do ip = 1, np 
     RpZ(2) = Pmin + ( Pmax - Pmin ) * ( ip - 1 ) / ( np - 0 ) / Mfp
     if ( modulo(myid, np) .ne. modulo(ip-1,nmaster) ) cycle
     do iz = 1, nz ; RpZ(3) = Zmin + ( Zmax - Zmin ) * ( iz - 1 ) / ( nz - 1 )
        do ir = 1, nr ; RpZ(1) = Rmin + ( Rmax - Rmin ) * ( ir - 1 ) / ( nr - 1 )
           czeta = cos(RpZ(2))
           szeta = sin(RpZ(2))
           xx = RpZ(1) * czeta
           yy = RpZ(1) * szeta
           zz = RpZ(3)
           call coils_bfield(B,xx,yy,zz)
           BRZp(1,ir,iz,ip) = (   B(1) * czeta + B(2) * szeta )
           BRZp(3,ir,iz,ip) = ( - B(1) * szeta + B(2) * czeta ) 
           BRZp(2,ir,iz,ip) =     B(3)
#ifdef DIV_CHECK
           dBx = B(1) ; dBy = B(2) ; dBz = B(3) 
           BRpZ(2,ir,iz,ip) = sqrt( B(1)*B(1) + B(2)*B(2) + B(3)*B(3) )

           xx = xx + dx
           call coils_bfield(B,xx,yy,zz)
           dBx = ( B(1) - dBx ) / dx
           xx = xx - dx

           yy = yy + dy
           call coils_bfield(B,xx,yy,zz)
           dBy = ( B(2) - dBy ) / dy
           yy = yy - dy

           zz = zz + dz
           call coils_bfield(B,xx,yy,zz)
           dBz = ( B(3) - dBz ) / dz 
           zz = zz - dz

           BRpZ(1,ir,iz,ip) = dBx + dBy + dBz
           BRpZ(2,ir,iz,ip) = BRpZ(1,ir,iz,ip) / BRpZ(2,ir,iz,ip)
#endif
        enddo
     enddo
  enddo

  if (masterid >=0) then
     call MPI_ALLREDUCE( MPI_IN_PLACE, BRZp, 3*nr*nz*np, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_MASTERS, ierr)
#ifdef DIV_CHECK
     call MPI_ALLREDUCE( MPI_IN_PLACE, BRpZ, 2*nr*nz*np, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_MASTERS, ierr)
#endif
     CALL MPI_COMM_FREE(MPI_COMM_MASTERS, ierr)
  endif

  CALL MPI_COMM_FREE(MPI_COMM_MYWORLD, ierr)

  if( myid.eq.0 ) then
#ifdef DIV_CHECK
     write(ounit, '("wtmgrid : max. div B = "ES23.15 " ; max. div B / |B| = "ES23.15 )') maxval(BRpZ(1,1:Nr,1:Nz,1:Np)),  maxval(BRpZ(2,1:Nr,1:Nz,1:Np))
#endif
     nextcur = 1 ; curlabel(1) = "focus-coils"

     write( ounit,'("wtmgrid : writing ",A," ; Mfp="i3" ;")')  trim(mgrid_name), Mfp

     !open( wunit, file=trim(ext)//".fo.mgrid", status="unknown", form="unformatted", iostat=iostat )
     open( wunit, file=trim(mgrid_name), status="unknown", form="unformatted", iostat=iostat )
     FATAL( wtmgrid, iostat.ne.0, error opening ext.fo.mgrid )
     write(wunit) Nr, Nz, Np, Mfp, nextcur
     write(wunit) Rmin, Zmin, Rmax, Zmax
     write(wunit) curlabel(1:nextcur)
     write(wunit) BRZp(1:3,1:Nr,1:Nz,1:Np)
     close(wunit)

  endif

  DEALLOCATE(BRZp)
#ifdef DIV_CHECK
  DEALLOCATE(BRpZ)
#endif

  return

end subroutine wtmgrid

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
