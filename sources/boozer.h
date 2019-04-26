subroutine boozmn
  USE globals, only : dp, myid, ncpu, zero, ounit, total_num, pp_maxiter, pp_ns, &
       XYZB, lboozmn, bmin, bmim, booz_mnc, booz_mns, booz_mpol, booz_ntor, booz_mn, nfp_raw
  USE mpi
  IMPLICIT NONE

  ! allocate data for following calculations

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  INTEGER              :: ierr, astat, iflag
  INTEGER              :: tor_num, in, im, imn
  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  call MPI_BARRIER( MPI_COMM_WORLD, ierr ) ! wait all cpus;

  FATAL( boozmn_01, booz_mpol < 0, invalid poloidal mode resolution )
  FATAL( boozmn_02, booz_ntor < 0, invalid toroidal mode resolution )

  lboozmn = .true. ! turn on boozmn
  tor_num = 360    ! toroidal planes number
  total_num = pp_maxiter * tor_num ! total data points per line

  booz_mpol = 16 ; booz_ntor = 32
  booz_mn = booz_mpol*(2*booz_ntor+1) + (booz_ntor+1) ! (1:M, -N:N) + (0, 0:N)
  SALLOCATE( bmin, (1:booz_mn), 0)
  SALLOCATE( bmim, (1:booz_mn), 0)
  SALLOCATE( booz_mnc, (1:booz_mn, 1:pp_ns), zero)
  SALLOCATE( booz_mns, (1:booz_mn, 1:pp_ns), zero)
  SALLOCATE( XYZB, (1:total_num, 1:4, 1:pp_ns), zero)

  ! prepare bmin & bmim
  imn = 0
  do im = 0, booz_mpol
     do in = -booz_ntor, booz_ntor
        if ( im==0 .and. in<0 ) cycle
        imn = imn + 1
        bmim(imn) = im
        bmin(imn) = in  !*Nfp_raw
     enddo
  enddo
  
  FATAL( boozmn_03, imn .ne. booz_mn, packing error )

end subroutine boozmn

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine boozsurf(XYZB, x, y, z, iota, isurf)
  USE globals, only : dp, myid, ncpu, zero, half, two, pi, pi2, ounit, total_num, pp_maxiter, &
                      bmin, bmim, booz_mnc, booz_mns, booz_mn, machprec, &
                      masterid
  USE mpi
  IMPLICIT NONE

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  REAL, dimension(total_num, 4) :: XYZB ! XYZB on one surface
  REAL   , intent( in) :: x, y, z ! starting point
  REAL   , intent( in) :: iota ! calculated iota
  INTEGER, intent( in) :: isurf ! fieldline ordering

  INTEGER              :: ierr, astat, iflag
  INTEGER              :: i, imn, tor_num, pol_num, iteta, jzeta
  REAL,dimension(total_num) :: chi, zeta, teta
  REAL                 :: Gpol, ang, dteta, dzeta
  INTEGER, allocatable :: weight(:,:)
  REAL, allocatable    :: Btz(:,:)

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  ! write(ounit, '("poincare: starting filed line tracing at x="F5.2, ", y="F5.2, ", z="F5.2)') x, y, z

  ! filedline tracing
  call fieldline_tracing(x, y, z, total_num, pp_maxiter, XYZB(1:total_num, 1:4))
 
  ! calculate chi = \int B dl
  chi = zero
  do i = 2, total_num 
     chi(i) = chi(i-1) + (XYZB(i, 4) + XYZB(i-1, 4))*half &
                 * sqrt( (XYZB(i, 1) - XYZB(i-1, 1))**2   &
                       + (XYZB(i, 2) - XYZB(i-1, 2))**2   & 
                       + (XYZB(i, 3) - XYZB(i-1, 3))**2   )
  enddo

  ! calculate poloidal current Gpol = \int \vetc{B} \cdot d \zeta
  ! Gpol = 2.0E-7 * total_current
  Gpol = chi(total_num) / (pi2*pp_maxiter)
  FATAL(booz_04 , abs(Gpol) < machprec, zero external poloidal currents)

!!$  ! Fourier decomposition
!!$  do imn = 1, booz_mn     
!!$     
!!$     do i = 1, total_num 
!!$        ang = (bmin(imn) - bmim(imn)*abs(iota))/Gpol * chi(i)
!!$        booz_mnc(imn, isurf) = booz_mnc(imn, isurf) + XYZB(i, 4) * cos(ang)
!!$        booz_mns(imn, isurf) = booz_mns(imn, isurf) + XYZB(i, 4) * sin(ang)
!!$     enddo
!!$
!!$     if ( bmim(imn) == 0 .and. bmin(imn) == 0 ) then
!!$        booz_mnc(imn, isurf) = booz_mnc(imn, isurf) * half
!!$        booz_mns(imn, isurf) = booz_mns(imn, isurf) * half
!!$     endif
!!$  enddo
!!$
!!$  booz_mnc(1:booz_mn, isurf) = booz_mnc(1:booz_mn, isurf) * two / total_num
!!$  booz_mns(1:booz_mn, isurf) = booz_mns(1:booz_mn, isurf) * two / total_num


  ! Boozer angles
  zeta = mod(chi/Gpol          , pi2)
  teta = mod(chi/Gpol*abs(iota), pi2)

  ! map back to two dimensional grid
  ! tor_num = total_num/pp_maxiter
  ! pol_num = pp_maxiter
  tor_num = 256
  pol_num = 128
  SALLOCATE( Btz   , (0:pol_num, 0:tor_num), zero )
  SALLOCATE( weight, (0:pol_num, 0:tor_num), 0    )
  dzeta = pi2/tor_num
  dteta = pi2/pol_num
  do i = 1, total_num
     iteta = int(teta(i)/dteta)
     jzeta = int(zeta(i)/dzeta)
     Btz(iteta, jzeta) =  Btz(iteta, jzeta) + XYZB(i, 4)
     weight(iteta, jzeta) = weight(iteta, jzeta) + 1
  enddo
 
  ! Fourier decomposition
  do jzeta = 0, tor_num-1
     do iteta = 0, pol_num-1
        ! weight(iteta, jzeta) = max(weight(iteta, jzeta), 1) ! avoida dividing zero
        if ( weight(iteta,jzeta) /= 0)  Btz(iteta,jzeta) = Btz(iteta,jzeta) / weight(iteta,jzeta)
        do imn = 1, booz_mn
           ang = bmim(imn) * iteta*dteta - bmin(imn) * jzeta*dzeta
           booz_mnc(imn, isurf) = booz_mnc(imn, isurf) + Btz(iteta,jzeta) * cos(ang)
           booz_mns(imn, isurf) = booz_mns(imn, isurf) + Btz(iteta,jzeta) * sin(ang)
        enddo
     enddo
  enddo

  booz_mnc(1:booz_mn, isurf) = booz_mnc(1:booz_mn, isurf) * two / (tor_num*pol_num)
  booz_mns(1:booz_mn, isurf) = booz_mns(1:booz_mn, isurf) * two / (tor_num*pol_num)

  imn = 1
  FATAL( boozer_05, bmim(imn) /= 0 .or. bmin(imn) /= 0, wrong mn initialization )
  booz_mnc(imn, isurf) = booz_mnc(imn, isurf) * half
  booz_mns(imn, isurf) = booz_mns(imn, isurf) * half
     
  DALLOCATE( Btz )
  DALLOCATE( weight )
!!$
!!$  ! Fourier decomposition
!!$  do imn = 1, booz_mn
!!$
!!$     booz_mnc(imn, isurf) = zero ; booz_mns(imn, isurf) = zero
!!$
!!$     do i = 1, total_num
!!$        ang = bmim(imn) * teta(i) - bmin(imn) * zeta(i)
!!$        booz_mnc(imn, isurf) = booz_mnc(imn, isurf) + XYZB(i, 4) * cos(ang)
!!$        booz_mns(imn, isurf) = booz_mns(imn, isurf) + XYZB(i, 4) * sin(ang)
!!$     enddo
!!$
!!$     if ( bmim(imn) == 0 .and. bmin(imn) == 0 ) then
!!$        booz_mnc(imn, isurf) = booz_mnc(imn, isurf) * half
!!$        booz_mns(imn, isurf) = booz_mns(imn, isurf) * half
!!$     endif
!!$  enddo
!!$
!!$  booz_mnc(1:booz_mn, isurf) = booz_mnc(1:booz_mn, isurf) * two / total_num
!!$  booz_mns(1:booz_mn, isurf) = booz_mns(1:booz_mn, isurf) * two / total_num

  ! finish decomposition

  write(ounit, '("boozmn  : myid="I6" ; Gpol="ES12.5" ; iota="ES12.5" ; Booz_mnc(1)="ES12.5 &
       " ; Booz_mns(1)="ES12.5)') masterid, Gpol, iota, booz_mnc(1, isurf), booz_mns(1, isurf)
  
  return
end subroutine boozsurf

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE Bmn_clt(phi_B,theta_B,B0,Bf)
implicit none

real*8,dimension(300000,1):: B0,phi_B,theta_B
real*8,dimension(30,30)::Bf
real*8,dimension(60,60)::Bmn
real*8,dimension(30,60)::theta,phi,ang,Bf1
real*8,dimension(61,61)::B1,n_B,B2
real*8,dimension(30,1)::phi1
real*8,dimension(1,60)::theta1
integer*4 m,n,i,j,k,l
real*8 pi,s1,s2,s3,s4,B_mn,B_c,B_s, dteta, dzeta

pi=3.141592653589793239
m=60
n=60
dzeta = 2*pi/n
dteta = 2*pi/m

do k=2,300000
    j=floor(theta_B(k,1)/(dzeta))
    i=floor(phi_B(k,1)/(dteta))

    
    s1=sqrt((theta_B(k,1)-dzeta* j   )**2+(phi_B(k,1)-dteta* i   )**2)
    s2=sqrt((theta_B(k,1)-dzeta* j   )**2+(phi_B(k,1)-dteta*(i+1))**2)
    s3=sqrt((theta_B(k,1)-dzeta*(j+1))**2+(phi_B(k,1)-dteta*(i+1))**2)
    s4=sqrt((theta_B(k,1)-dzeta*(j+1))**2+(phi_B(k,1)-dteta* i   )**2)
    i=i+1
    j=j+1

    B1(i,j)=B1(i,j)+B0(k,1)/s1
    n_B(i,j)=n_B(i,j)+1/s1
    B1(i+1,j)=B1(i+1,j)+B0(k,1)/s2
    n_B(i+1,j)=n_B(i+1,j)+1/s2
    B1(i+1,j+1)=B1(i+1,j+1)+B0(k,1)/s3
    n_B(i+1,j+1)=n_B(i+1,j+1)+1/s3
    B1(i,j+1)=B1(i,j+1)+B0(k,1)/s4
    n_B(i,j+1)=n_B(i,j+1)+1/s4
enddo
B1(1:m,1)=B1(1:m,1)+B1(1:m,n+1)
n_B(1:m,1)=n_B(1:m,1)+n_B(1:m,n+1)
B1(1,1:n)=B1(1,1:n)+B1(m+1,1:n)
n_B(1,1:n)=n_B(1,1:n)+n_B(m+1,1:n)
do i=1,m+1
do j=1,n+1
 B2(i,j)=B1(i,j)/n_B(i,j)
enddo
enddo

Bmn=B2(1:m,1:n)


do i=1,m/2
phi1(i,1)=dteta*2*(i-1)
enddo

do i=1,n
theta1(1,i)=dzeta*(i-1)
enddo

do i=1,n
phi(:,i)=phi1(:,1)
enddo

do i=1,m/2
theta(i,:)=theta1(1,:)
enddo


do i=0,m/2-1
    do j=0,n-1
        ang=-i*Phi+j*Theta
        B_mn=0
        B_c=0
        B_s=0
        do k=1,m/2
        do  l=1,n 
            B_c=B_c+Bmn(k,l)*cos(ang(k,l))
            B_s=B_s+Bmn(k,l)*sin(ang(k,l))
            B_mn=sqrt(B_c**2+B_s**2)/float(n*m/2)
        enddo
        enddo
        Bf1(i+1,j+1)=B_mn
 !       write(*,*)B_mn
    enddo
enddo

Bf=Bf1(:,1:n/2)*2.0
Bf(1,1)=Bf1(1,1)
return
end SUBROUTINE Bmn_clt

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
