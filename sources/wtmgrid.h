! write binary mgrid file

subroutine wtmgrid
  use globals, only : dp, zero, half, pi2, ext, ncpu, myid, ounit, wunit, nfp_raw,  &
       pp_raxis, pp_zaxis, pp_rmax, pp_zmax, mgrid_rmin, mgrid_rmax, mgrid_zmin, mgrid_zmax, &
       master, nmaster, nworker, masterid, color, myworkid, MPI_COMM_MASTERS, MPI_COMM_MYWORLD, MPI_COMM_WORKERS
  implicit none
  include "mpif.h"

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  LOGICAL              :: exist
  INTEGER              :: ierr, astat, iostat, ip, iz, ir, np, nz, nr, Mfp, nextcur, icommand
  REAL                 :: RpZ(1:3), R, P, Z, Pmin, Pmax, Zmin, Zmax, Rmin, Rmax, B(1:3), pressure, gap, &
       czeta, szeta, xx, yy, zz, dx, dy, dz, dBx, dBy, dBz
  REAL, allocatable    :: BRZp(:,:,:,:), BRpZ(:,:,:,:)
  CHARACTER(LEN=100)   :: mgrid_name
  CHARACTER(LEN=30)    :: curlabel(1:1)

  icommand = 0
  mgrid_name = "mgrid.focus_"//trim(ext) ! filename, could be user input
  np =  72  ; nz = 121  ; nr = 121 ; Mfp = nfp_raw ! SHOULD BE USER INPUT; 04 Aug 16;  
  !np =  12  ; nz = 11  ; nr = 11 ; Mfp = 2 ! SHOULD BE USER INPUT; 04 Aug 16;
  B = zero  ; dx = 1E-4 ; dy = 1E-4 ; dz = 1E-4

  SALLOCATE( BRZp, (1:3,1:Nr,1:Nz,1:Np), zero )
#ifdef DIV_CHECK
  SALLOCATE( BRpZ, (1:2,1:Nr,1:Nz,1:Np), zero )
#endif

  Pmin = zero ; Pmax = pi2 ! DO NOT CHANGE; 04 Aug 16;

 ! Rmin =  4.8 ; Rmax =  6.8
 ! Zmin = -1.5 ; Zmax =  1.5
  
  Rmin = mgrid_rmin ; Zmin = mgrid_zmin
  Rmax = mgrid_rmax ; Zmax = mgrid_zmax

  if( myid.eq.0 ) write( ounit,'("wtmgrid : writing mgrid file at grid of [ "4(ES12.5,2X)" ]",3i6)') Rmin, Rmax, Zmin, Zmax, np, nr, nz

  do ip = 1, np 
     RpZ(2) = Pmin + ( Pmax - Pmin ) * ( ip - 1 ) / ( np - 0 ) / Mfp
     ! if ( modulo(myid, np) .ne. modulo(ip-1,nmaster) ) cycle
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

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
