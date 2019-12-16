SUBROUTINE specinp
  !-------------------------------------------------------------------------------!
  ! prepare inputs for SPEC
  ! This subroutine contains three functions:
  ! 1. Decompose the Bn on the computational boundary (read from plasma.boundary)
  ! 2. Calculate the poloidal and toroidal closed currents (Itor and Gpol)
  ! 3. Write down a xxx.Vns file with all the information for SPEC
  !-------------------------------------------------------------------------------!
  use globals, only: dp, zero, half, two, pi2, mu0,  myid, wunit, ounit,  surf, bn, ext, &
                     Nfp_raw, Nteta, Nzeta, plasma
  implicit none
  include "mpif.h"
  !-------------------------------------------------------------------------------
  INTEGER             :: mf, nf  ! Fourier modes size
  INTEGER             :: imn=0, ii, jj, im, in, astat, ierr, Nbf, iteta, jzeta, isurf
  REAL                :: teta, zeta, arg, tol, tmpc, tmps, curtor, curpol
  INTEGER, allocatable:: bnim(:), bnin(:)
  REAL   , allocatable:: bnc(:), bns(:)

  ! use the plasma for now; could be the limiter surface; 2019/12/15
  isurf = plasma
  ! default Fourier resolution; could be customized
  mf = 24 ; nf = 12
  ! compute Bn
  call bnormal(0) ! calculate Bn

  if (myid .ne. 0) return

  Nbf = (mf+1)*(2*nf+1) ! (0:mf)*(-nf:nf)

  SALLOCATE( bnim, (1:Nbf), 0    )
  SALLOCATE( bnin, (1:Nbf), 0    )
  SALLOCATE( bnc , (1:Nbf), zero )
  SALLOCATE( bns , (1:Nbf), zero )  

  do im = 0, mf
     do in = -nf, nf
        tmpc = zero ; tmps = zero
        do ii = 0, Nteta-1 
           teta = ( ii + half ) * pi2 / Nteta
           do jj = 0, Nzeta-1
              zeta = ( jj + half ) * pi2 / Nzeta
              arg = im*teta - in*Nfp_raw*zeta
              tmpc = tmpc + (-bn(ii, jj)*surf(isurf)%ds(ii,jj))*cos(arg)  ! minus sign is required because
              tmps = tmps + (-bn(ii, jj)*surf(isurf)%ds(ii,jj))*sin(arg)  ! the normal vector in SPEC is e_t x e_z
           enddo ! end jj
        enddo ! end ii

      ! if ( (abs(tmpc) + abs(tmps)) .lt. tol ) cycle
        if ( im .eq. 0 .and. in .lt. 0) cycle ! neglect m=0, n<0 terms

        imn = imn + 1
       ! bnin(imn) = in * bNfp ; bnim(imn) = im
        bnin(imn) = in         ; bnim(imn) = im

      ! if (im .eq. 0  ) then
      !    tmpc = tmpc*half
      !    tmps = tmps*half
      ! endif

        bnc(imn) = tmpc
        bns(imn) = tmps

     enddo ! end im
  enddo ! end in

  Nbf = imn

  bnc = bnc * two / (Nteta*Nzeta)
  bns = bns * two / (Nteta*Nzeta)

  ! compute curpol and curtor
  curtor = zero ; curpol = zero

  jzeta = 0
  do iteta = 0, Nteta-1
     curtor = curtor + surf(isurf)%Bx(iteta,jzeta)*surf(isurf)%xt(iteta,jzeta) &
                   & + surf(isurf)%By(iteta,jzeta)*surf(isurf)%yt(iteta,jzeta) &
                   & + surf(isurf)%Bz(iteta,jzeta)*surf(isurf)%zt(iteta,jzeta) 
  enddo
  curtor = curtor * pi2/Nteta ! / mu0  ! SPEC currents are normalized with mu0

  iteta = 0
  do jzeta = 0, Nzeta-1
     curpol = curpol + surf(isurf)%Bx(iteta,jzeta)*surf(isurf)%xp(iteta,jzeta) &
                   & + surf(isurf)%By(iteta,jzeta)*surf(isurf)%yp(iteta,jzeta) &
                   & + surf(isurf)%Bz(iteta,jzeta)*surf(isurf)%zp(iteta,jzeta) 
  enddo
  curpol = curpol * pi2/Nzeta 

  ! write SPEC input
  write(ounit, '("postproc: preparing input information for SPEC. Please view "A)'), 'focus_'//trim(ext)//'.sp'
  open(wunit, file='focus_'//trim(ext)//'.sp', status='unknown', action='write')

!  write(wunit, '(2I)') Nfou, Nfp_raw
!  do imn = 1, Nfou
!     write(wunit, '(2I, 4ES23.15)') bim(imn), bin(imn), Rbc(imn), Rbs(imn), Zbc(imn), Zbs(imn)
!  end do
  
!  write(wunit, '(2I, 2ES23.15)') Nbf, Nfp_raw, curtor, curpol
!  do imn = 1, Nbf
!     write(wunit, '(2I, 2ES23.15)') bnim(imn), bnin(imn), bnc(imn), bns(imn)
!  end do

  write(wunit,'(" curtor      = ",es23.15   )') curtor
  write(wunit,'(" curpol      = ",es23.15   )') curpol
  write(wunit,'(" Nfp         = ",i9        )') Nfp_raw
  do imn = 1, surf(isurf)%Nfou
     write(wunit,1010) surf(isurf)%bin(imn)/Nfp_raw, surf(isurf)%bim(imn), surf(isurf)%Rbc(imn), &
                       surf(isurf)%bin(imn)/Nfp_raw, surf(isurf)%bim(imn), surf(isurf)%Zbs(imn), &
                       surf(isurf)%bin(imn/Nfp_raw), surf(isurf)%bim(imn), surf(isurf)%Rbs(imn), &
                       surf(isurf)%bin(imn)/Nfp_raw, surf(isurf)%bim(imn), surf(isurf)%Zbc(imn)  ! wall is read as plasma boundary
  enddo
  do imn = 1, Nbf
     write(wunit,1020) bnin(imn), bnim(imn), bns(imn), bnin(imn), bnim(imn), zero, &
                       bnin(imn), bnim(imn), bnc(imn), bnin(imn), bnim(imn), zero  ! only vacuum data is calculated
  enddo  

1010 format("Rwc(",i3,",",i3,")",2x,"=",es23.15," Zws(",i3,",",i3,")",2x,"=",es23.15," Rws(",i3,",",i3,")",2x,"=",es23.15," Zwc(",i3,",",i3,")",2x,"=",es23.15)
1020 format("Vns(",i3,",",i3,")",2x,"=",es23.15," Bns(",i3,",",i3,")",2x,"=",es23.15," Vnc(",i3,",",i3,")",2x,"=",es23.15," Bnc(",i3,",",i3,")",2x,"=",es23.15)

  close(wunit)

  return
END SUBROUTINE specinp
