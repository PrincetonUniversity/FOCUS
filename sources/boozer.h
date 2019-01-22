subroutine last_surface
  USE globals, only : dp, myid, ncpu, zero, half, pi, pi2, ounit, pi, total_num, pp_maxiter, XYZB
  USE mpi
  IMPLICIT NONE

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  INTEGER              :: ierr, astat, iflag
  INTEGER              :: tor_num
  REAL                 :: theta, zeta, r, x, y, z

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  ! call wtmgrid
  call MPI_BARRIER( MPI_COMM_WORLD, ierr ) ! wait all cpus;

  tor_num = 360    ! toroidal planes number
  total_num = pp_maxiter * tor_num

  SALLOCATE( XYZB, (1:total_num, 1:4), zero)

  ! starting point
  theta = zero ; zeta = zero
  call surfcoord( theta, zeta, r, z)
  x = r*cos(zeta)
  y = r*sin(zeta)
  
  if ( myid /= 0 ) return

  write(ounit, '("poincare: starting filed line tracing at x="F5.2, ", y="F5.2, ", z="F5.2)') x, y, z

  ! filedline tracing
  call fieldline_tracing(x,y,z,total_num,pp_maxiter,XYZB)

  write(ounit, '("poincare: Fieldline tracing finished")')
  
  return
end subroutine last_surface

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine fieldline_tracing(x,y,z,imax,n2,H)
  implicit none
  integer*4 ::n2,imax,j,i
  real*8 :: x,y,z,dphi,pi,dt,B,Bx,By,Bz,x0,y0,z0,g,iota
  real*8 :: s(4), k1x,k2x,k3x,k4x,k5x,k6x,k7x,k8x,k9x,k10x
  real*8 :: k1y,k2y,k3y,k4y,k5y,k6y,k7y,k8y,k9y,k10y,xr
  real*8 :: k1z,k2z,k3z,k4z,k5z,k6z,k7z,k8z,k9z,k10z
  real*8,dimension(imax,4):: H

  real*8,dimension(imax+1,4):: f
  real*8,dimension(2*n2,3):: f2

  pi=3.141592653589793239
  dphi=2*pi/(float(imax)/n2)

  do j=1,imax
     H(j,1)=x ; H(j,2)=y ; H(j,3)=z
     call coils_bfield(s,x,y,z)
     Bx=s(1)  ; By=s(2)  ; Bz=s(3)  ; B=s(4)
     H(j,4)=B
     dt=(y-x*tan(j*dphi))/(tan(j*dphi)*Bx/sqrt(Bx**2+By**2)-By/sqrt(Bx**2+By**2))*sqrt(B**2/(Bx**2+By**2))

     f(j,1)=Bx/B
     f(j,2)=By/B
     f(j,3)=Bz/B
     x0=x  ; y0=y  ; z0=z

     if(j<8)then
        k1x=Bx/B ; k1y=By/B ; k1z=Bz/B

        call coils_bfield(s,x+dt*4/27*k1x,y+dt*4/27*k1y,z+dt*4/27*k1z)
        Bx=s(1) ; By=s(2) ; Bz=s(3) ; B=s(4)
        k2x=Bx/B ; k2y=By/B ; k2z=Bz/B

        call coils_bfield(s,x+dt/18*(k1x+3*k2x),y+dt/18*(k1y+3*k2y),z+dt/18*(k1z+3*k2z))
        Bx=s(1) ; By=s(2) ; Bz=s(3) ; B=s(4)
        k3x=Bx/B ; k3y=By/B ; k3z=Bz/B

        call coils_bfield(s,x+dt/12*(k1x+3*k3x),y+dt/12*(k1y+3*k3y),z+dt/12*(k1z+3*k3z))
        Bx=s(1) ; By=s(2) ; Bz=s(3) ; B=s(4)
        k4x=Bx/B ; k4y=By/B ; k4z=Bz/B

        call coils_bfield(s,x+dt/8*(k1x+3*k4x),y+dt/8*(k1y+3*k4y),z+dt/8*(k1z+3*k4z))
        Bx=s(1) ; By=s(2) ; Bz=s(3) ; B=s(4)
        k5x=Bx/B ; k5y=By/B ; k5z=Bz/B

        call coils_bfield(s,x+dt/54*(13*k1x-27*k3x+42*k4x+8*k5x),y+dt/54*(13*k1y-27*k3y+&
             42*k4y+8*k5y),z+dt/54*(13*k1z-27*k3z+42*k4z+8*k5z))
        Bx=s(1) ; By=s(2) ; Bz=s(3) ; B=s(4)
        k6x=Bx/B ; k6y=By/B ; k6z=Bz/B

        call coils_bfield(s,x+dt/4320*(389*k1x-54*k3x+966*k4x-824*k5x+243*k6x),y+dt/4320*(389*k1y-&
             54*k3y+966*k4y-824*k5y+243*k6y),z+dt/4320*(389*k1z-54*k3z+966*k4z-824*k5z+243*k6z))
        Bx=s(1) ; By=s(2) ; Bz=s(3) ; B=s(4)
        k7x=Bx/B ; k7y=By/B ; k7z=Bz/B

        call coils_bfield(s,x+dt/20*(-234*k1x+81*k3x-1164*k4x+656*k5x-122*k6x+800*k7x),y+dt/20*(-234*k1y+81*k3y-&
             1164*k4y+656*k5y-122*k6y+800*k7y),z+dt/20*(-234*k1z+81*k3z-1164*k4z+656*k5z-122*k6z+800*k7z))
        Bx=s(1) ; By=s(2) ; Bz=s(3) ; B=s(4)
        k8x=Bx/B ; k8y=By/B ; k8z=Bz/B

        call coils_bfield(s,x+dt/288*(-127*k1x+18*k3x-678*k4x+456*k5x-9*k6x+576*k7x+4*k8x),y+&
             dt/288*(-127*k1y+18*k3y-678*k4y+456*k5y-9*k6y+576*k7y+4*k8y),z+dt/288*(-127*k1z+&
             18*k3z-678*k4z+456*k5z-9*k6z+576*k7z+4*k8z))
        Bx=s(1) ; By=s(2) ; Bz=s(3) ; B=s(4)
        k9x=Bx/B ; k9y=By/B ; k9z=Bz/B

        call coils_bfield(s,x+dt/820*(1481*k1x-81*k3x+7104*k4x-3376*k5x+&
             72*k6x-5040*k7x-60*k8x+720*k9x),y+dt/820*(1481*k1y-81*k3y+&
             7104*k4y-3376*k5y+72*k6y-5040*k7y-60*k8y+720*k9y),z+dt/820*(1481*k1z-&
             81*k3z+7104*k4z-3376*k5z+72*k6z-5040*k7z-60*k8z+720*k9z))
        Bx=s(1) ; By=s(2) ; Bz=s(3) ; B=s(4)
        k10x=Bx/B ; k10y=By/B ; k10z=Bz/B

        x=x+dt/840*(41*k1x+27*k4x+272*k5x+27*k6x+216*k7x+216*k9x+41*k10x)
        y=y+dt/840*(41*k1y+27*k4y+272*k5y+27*k6y+216*k7y+216*k9y+41*k10y)
        z=z+dt/840*(41*k1z+27*k4z+272*k5z+27*k6z+216*k7z+216*k9z+41*k10z)

     else
        x=x+dt/120960*(-36799.0*f(j-7,1)+295767.0*f(j-6,1)-1041723.0*f(j-5,1)&
             +2102243.0*f(j-4,1)-2664477.0*f(j-3,1)+2183877.0*f(j-2,1)-1152169.0*f(j-1,1)+434241.0*f(j,1))
        y=y+dt/120960*(-36799.0*f(j-7,2)+295767.0*f(j-6,2)-1041723.0*f(j-5,2)&
             +2102243.0*f(j-4,2)-2664477.0*f(j-3,2)+2183877.0*f(j-2,2)-1152169.0*f(j-1,2)+434241.0*f(j,2))
        z=z+dt/120960*(-36799.0*f(j-7,3)+295767.0*f(j-6,3)-1041723.0*f(j-5,3)&
             +2102243.0*f(j-4,3)-2664477.0*f(j-3,3)+2183877.0*f(j-2,3)-1152169.0*f(j-1,3)+434241.0*f(j,3))
     end if

     call coils_bfield(s,x,y,z)
     Bx=s(1) ; By=s(2) ; Bz=s(3) ; B=s(4)
     f(j+1,1)=Bx/B
     f(j+1,2)=By/B
     f(j+1,3)=Bz/B

     if (j>7) then
        x=x0+dt/120960*(1375.0*f(j-6,1)-11351.0*f(j-5,1)+41499.0*f(j-4,1)-88547.0*f(j-3,1)&
             +123133.0*f(j-2,1)-121797.0*f(j-1,1)+139849.0*f(j,1)+36799.0*f(j+1,1))
        y=y0+dt/120960*(1375.0*f(j-6,2)-11351.0*f(j-5,2)+41499.0*f(j-4,2)-88547.0*f(j-3,2)&
             +123133.0*f(j-2,2)-121797.0*f(j-1,2)+139849.0*f(j,2)+36799.0*f(j+1,2))
        z=z0+dt/120960*(1375.0*f(j-6,3)-11351.0*f(j-5,3)+41499.0*f(j-4,3)-88547.0*f(j-3,3)&
             +123133.0*f(j-2,3)-121797.0*f(j-1,3)+139849.0*f(j,3)+36799.0*f(j+1,3))
     end if

  end do
  return

end subroutine fieldline_tracing


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine wtmgrid
  use globals, only : dp, zero, half, pi2, ext, ncpu, myid, ounit, wunit
  implicit none

  include "mpif.h"


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  LOGICAL              :: exist
  INTEGER              :: ierr, astat, iostat, imn, ibn00aa, ip, iz, ir, np, nz, nr, itangent, ibfield, Mfp, nextcur
  REAL                 :: RpZ(1:3), R, P, Z, Pmin, Pmax, Zmin, Zmax, Rmin, Rmax, B(1:4), pressure, gap, &
                          czeta, szeta, xx, yy, zz, dx, dy, dz, dBx, dBy, dBz
  REAL, allocatable    :: BRZp(:,:,:,:), dBRZp(:,:,:,:), BRpZ(:,:,:,:), dBRpZ(:,:,:,:)
  CHARACTER*13         :: suffix
  CHARACTER*30         :: curlabel(1:1)

  np =  72  ; nz = 121  ; nr = 121 ; Mfp = 2 ! SHOULD BE USER INPUT; 04 Aug 16;  
  !np =  12  ; nz = 11  ; nr = 11 ; Mfp = 2 ! SHOULD BE USER INPUT; 04 Aug 16;
  B = zero  ; dx = 1E-4 ; dy = 1E-4 ; dz = 1E-4
  
  SALLOCATE( BRZp, (1:3,1:Nr,1:Nz,1:Np), zero )
  SALLOCATE(dBRZp, (1:3,1:Nr,1:Nz,1:Np), zero )
  SALLOCATE( BRpZ, (1:2,1:Nr,1:Nz,1:Np), zero )
  SALLOCATE(dBRpZ, (1:2,1:Nr,1:Nz,1:Np), zero )

  Pmin = zero ; Pmax = pi2 ! DO NOT CHANGE; 04 Aug 16;

 !call plasdim(Rmin, Rmax, Zmin, Zmax)  !calculate plasma surface boundary  ;09/11/2016
 !call coildim(Rmin, Rmax, Zmin, Zmax)

 !gap = 0.3
 !Rmin = Rmin -gap; Rmax = Rmax + gap
 !Zmin = Zmin -gap; Zmax = Zmax + gap

  Rmin =  2.8 ; Rmax =  3.2
  Zmin = -0.2 ; Zmax =  0.2

  if( myid.eq.0 ) write( ounit,'("wtmgrid : writing mgrid file at grid of [ "4(ES12.5,2X)" ]",3i6)') Rmin, Rmax, Zmin, Zmax, np, nr, nz

  do ip = 1, np ; RpZ(2) = Pmin + ( Pmax - Pmin ) * ( ip - 1 ) / ( np - 0 ) / Mfp

     if ( myid.ne.modulo(ip,ncpu) ) cycle

     do iz = 1, nz ; RpZ(3) = Zmin + ( Zmax - Zmin ) * ( iz - 1 ) / ( nz - 1 )

        do ir = 1, nr ; RpZ(1) = Rmin + ( Rmax - Rmin ) * ( ir - 1 ) / ( nr - 1 )
                      
           czeta = cos(RpZ(2))
           szeta = sin(RpZ(2))

           xx = RpZ(1) * czeta
           yy = RpZ(1) * szeta
           zz = RpZ(3)

           call coils_bfield(B,xx,yy,zz) 

           dBRZp(1,ir,iz,ip) = (   B(1) * czeta + B(2) * szeta )
           dBRZp(3,ir,iz,ip) = ( - B(1) * szeta + B(2) * czeta ) 
           dBRZp(2,ir,iz,ip) =     B(3)

           dBx = B(1) ; dBy = B(2) ; dBz = B(3) 
           dBRpZ(2,ir,iz,ip) = B(4)

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
           
           ! write(ounit, '("(x, y, z) = "3ES12.5" ; div B = " ES12.5)') xx, yy, zz, dBx + dBy + dBz
           dBRpZ(1,ir,iz,ip) = dBx + dBy + dBz
           dBRpZ(2,ir,iz,ip) = dBRpZ(1,ir,iz,ip) / dBRpZ(2,ir,iz,ip)

        enddo

     enddo

  enddo

  call MPI_Reduce(dBRZp, BRZp, 3*nr*nz*np, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_Reduce(dBRpZ, BRpZ, 2*nr*nz*np, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

  if( myid.eq.0 ) then

     write(ounit, '("wtmgrid : max. div B = "ES23.15 " ; max. div B / |B| = "ES23.15 )') maxval(BRpZ(1,1:Nr,1:Nz,1:Np)),  maxval(BRpZ(2,1:Nr,1:Nz,1:Np))

     nextcur = 1 ; curlabel(1) = "focus-space-curves"

     write(suffix,'(i3.3,".",i4.4,".",i4.4)') Np, Nr, Nz

     write( ounit,'("wtmgrid : writing mgrid.ext.f."i3.3"."i4.4"."i4.4" ; Mfp="i3" ;")') np, nr, nz, Mfp

    !open( wunit, file=trim(ext)//".fo.mgrid", status="unknown", form="unformatted", iostat=iostat )
     open( wunit, file="mgrid."//trim(ext)//".f."//suffix, status="unknown", form="unformatted", iostat=iostat )
     FATAL( wtmgrid, iostat.ne.0, error opening ext.fo.mgrid )
     write(wunit) Nr, Nz, Np, Mfp, nextcur
     write(wunit) Rmin, Zmin, Rmax, Zmax
     write(wunit) curlabel(1:nextcur)
     write(wunit) BRZp(1:3,1:Nr,1:Nz,1:Np)
     close(wunit)

  endif

  DEALLOCATE(dBRZp)
  DEALLOCATE( BRZp)
  DEALLOCATE(dBRpZ)
  DEALLOCATE( BRpZ)


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  return

end subroutine wtmgrid
