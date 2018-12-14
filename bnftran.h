!title(bnftran) ! Transforming the Bnormal distribution into Fourier harmonics.

!latex \briefly{Transforming former saved Bnormal data on predefined surface grid into Fourier harmonics.}
!latex \calledby{\link{identfy}}
!latex 
!latex \subsection{overview}
!latex \bi
!latex \item[1.] Discrete Bnormal distribution on the surface grid can be represented by a two-dimensional
!latex           Forurier series.
!latex    \be
!latex    \ds B_n(\t,\z) \equiv \sum_{m,n} B_{mn}^c cos(m\t-n\z) \ + \ B_{mn}^s sin(m\t-n\z)
!latex    \ee
!latex \item[2.] So, the Fourier harmonics $ B_{mn}^c$ and  $ B_{mn}^s$ can be computed as,
!latex    \be
!latex    \ds  B_{mn}^c & = & \frac{1}{2\pi^2} \int_0^{2\pi} \int_0^{2\pi} B_n(\t,\z) cos(m\t-n\z)
!l tex                         d\t d\z \quad \\ %\text{(for (0,0) term should be $\frac{1}{4\pi^2}$)} \\
!latex                         d\t d\z \quad     \mbox{\rm (for $(0,0)$ term should be $\frac{1}{4\pi^2}$)} \\
!latex                  &\approx & \frac{2}{Nteta \times Nzeta} \sum_0^{Nteta-1} \sum_0^{Nzeta-1} 
!latex                           B_n(iteta, jzeta) cos(m\t-n\z) \\
!latex    \ds  B_{mn}^s & = &\frac{1}{2\pi^2} \int_0^{2\pi} \int_0^{2\pi} B_n(\t,\z) sin(m\t-n\z)
!l tex                         d\t d\z quad \\ %\text{(for (0,0) term should be $\frac{1}{4\pi^2}$)} \\
!latex                         d\t d\z quad     \mbox{(for $(0,0)$ term should be $\frac{1}{4\pi^2}$)} \\
!latex                  &\approx & \frac{2}{Nteta \times Nzeta} \sum_0^{Nteta-1} \sum_0^{Nzeta-1} 
!latex                           B_n(iteta, jzeta) sin(m\t-n\z)
!latex    \ee
!latex  \item[3.] If there exist conjugated terms (e.g. m,n and  -m,-n), the computed terms will be 
!latex            timed with a factor of 2 or 0;
!latex  \item[4.] What should be noted is that the maximum modes in Fouries harmonics should not be 
!latex            less/equal than the half of surface grid resolutions.
!latex \ei

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

SUBROUTINE bnftran
  !-------------------------------------------------------------------------------!
  ! Calculate Fourier harmonics bnc,bns for bn(0:Nteta, 0:Nzeta);                 !
  ! CZHU; first version: 2017/01/11; last revised: 2017/01/11                     !
  !-------------------------------------------------------------------------------!
  use kmodule, only: zero, half, two, pi2, myid, ounit, &
       Nteta, Nzeta, surf, tbn, Nbf, bnim, bnin, bnc, bns, bNfp, Mpol, Ntor
  implicit none
  include "mpif.h"
  !-------------------------------------------------------------------------------
  INTEGER             :: mf, nf  ! predefined Fourier modes size
  INTEGER             :: imn=0, ii, jj, im, in, astat, ierr, maxN, maxM
  REAL                :: teta, zeta, arg
  !------------------------------------------------------------------------------- 

  mf = Mpol ;  nf = Ntor
  FATAL(bnftran, mf .le. 0 .and. nf .le. 0, INVALID size for Fourier harmonics)

#ifdef RMP
  !----------------------------------------------
  FATAL(bnftran, nbf .le. 0, RMP option requires nonzero input Bn)

  maxN = maxval(abs(bnin))
  maxM = maxval(abs(bnim))

  FATAL(bnftran, maxN .ge. Nzeta/2, toroidal grid resolution not enough)
  FATAL(bnftran, maxM .ge. Nteta/2, poloidal grid resolution not enough)

  bnc = zero ; bns = zero

  do imn = 1, Nbf

     in = bnin(imn)
     im = bnim(imn)

     do ii = 0, Nteta-1 
        teta = ( ii + half ) * pi2 / Nteta
        do jj = 0, Nzeta-1
           zeta = ( jj + half ) * pi2 / Nzeta
           arg = im*teta - in*zeta
           bnc(imn) = bnc(imn) + tbn(ii,jj)*cos(arg)
           bns(imn) = bns(imn) + tbn(ii,jj)*sin(arg)
        enddo ! end jj
     enddo ! end ii

     if (im .eq. 0 .and. in .eq. 0) then
        bnc(imn) = bnc(imn)*half
        bns(imn) = bns(imn)*half
     endif

  enddo

  bnc = bnc * two / (Nteta*Nzeta)
  bns = bns * two / (Nteta*Nzeta)

  !----------------------------------------------
#else
  !----------------------------------------------
  if(nbf .gt. 0) then  ! if there is input Bn target
     DALLOCATE(bnim)
     DALLOCATE(bnin)
     DALLOCATE(bnc )
     DALLOCATE(bns )
  endif

  Nbf = (mf+1)*(2*nf+1) ! (0:mf)*(-nf:nf)

  SALLOCATE( bnim, (1:Nbf), 0    )
  SALLOCATE( bnin, (1:Nbf), 0    )
  SALLOCATE( bnc , (1:Nbf), zero )
  SALLOCATE( bns , (1:Nbf), zero )  

  do in = -nf, nf
     do im = 0, mf

        imn = imn + 1
        bnin(imn) = in * bNfp ; bnim(imn) = im

        do ii = 0, Nteta-1 
           teta = ( ii + half ) * pi2 / Nteta
           do jj = 0, Nzeta-1
              zeta = ( jj + half ) * pi2 / Nzeta
              arg = im*teta - in*bNfp*zeta
              !bnc(imn) = bnc(imn) + (tbn(ii, jj)-surf(1)%bnt(ii,jj))*cos(arg)
              !bns(imn) = bns(imn) + (tbn(ii, jj)-surf(1)%bnt(ii,jj))*sin(arg)

              !bnc(imn) = bnc(imn) + (surf(1)%bnt(ii,jj))*cos(arg)
              !bns(imn) = bns(imn) + (surf(1)%bnt(ii,jj))*sin(arg)

              bnc(imn) = bnc(imn) + (tbn(ii, jj))*cos(arg)
              bns(imn) = bns(imn) + (tbn(ii, jj))*sin(arg)
           enddo ! end jj
        enddo ! end ii

        if (im .eq. 0 ) then
           bnc(imn) = bnc(imn)*half
           bns(imn) = bns(imn)*half
        endif

     enddo ! end im
  enddo ! end in

  bnc = bnc * two / (Nteta*Nzeta)
  bns = bns * two / (Nteta*Nzeta)
  !----------------------------------------------
#endif

  call write_plasma
END SUBROUTINE bnftran


SUBROUTINE Bmodule
  !-------------------------------------------------------------------------------!
  ! Calculate Fourier harmonics Bmc, Bms for bm(0:Nteta, 0:Nzeta);                !
  !-------------------------------------------------------------------------------!
  use kmodule, only: zero, half, two, pi2, myid, ounit, ext, lunit, &
       Nteta, Nzeta, surf, SaveBx, SaveBy, SaveBz, Mpol, Ntor, Bmod_n, &
       bmn, bNfp, bin, bim, rbc, rbs, zbc, zbs
  implicit none
  include "mpif.h"
  !-------------------------------------------------------------------------------
  INTEGER             :: mf, nf  ! Fourier modes size
  INTEGER             :: imn=0, ii, jj, im, in, mn, astat, ierr
  REAL                :: teta, zeta, arg
  INTEGER, allocatable:: bmin(:)
  REAL, allocatable   :: bm(:,:), bm_zeta(:)
  !------------------------------------------------------------------------------- 

  if (myid .ne. 0) return

  nf = Ntor

  FATAL(Bmodule, mf .le. 0 .and. nf .le. 0, INVALID size for Fourier harmonics)

  FATAL(Bmodule, nf .ge. Nzeta/2, toroidal grid resolution not enough)

  SALLOCATE( bm, (0:Nteta, 0:Nzeta), zero )

  bm = sqrt(SaveBx**2 + SaveBy**2 + SaveBz**2)

  SALLOCATE( bmin, (0:nf), 0    )
  SALLOCATE( bm_zeta , (0:Nzeta-1), zero ) 
  SALLOCATE( Bmod_n , (0:nf), zero ) 

  do jj = 0, Nzeta-1
     do ii = 0, Nteta-1
        bm_zeta(jj) = bm_zeta(jj) + bm(ii, jj)
     enddo
  enddo

  bm_zeta = bm_zeta / Nteta
  
  do in = 0, nf
     bmin(in) = in
     do jj = 0, Nzeta-1
        zeta = ( jj + half ) * pi2 / Nzeta
        Bmod_n(in) = Bmod_n(in) + bm_zeta(jj)*cos(bmin(in)*zeta)
     enddo
  enddo
  
  Bmod_n = Bmod_n * two / Nzeta
  Bmod_n(0) = Bmod_n(0) * half

!!$  do im = 0, mf
!!$     do in = -nf, nf
!!$
!!$        imn = imn + 1
!!$        bmin(imn) = in ; bmim(imn) = im
!!$
!!$        do ii = 0, Nteta-1 
!!$           teta = ( ii + half ) * pi2 / Nteta
!!$           do jj = 0, Nzeta-1
!!$              zeta = ( jj + half ) * pi2 / Nzeta
!!$              arg = im*teta - in*zeta
!!$
!!$              bmc(imn) = bmc(imn) + (bm(ii, jj))*cos(arg)
!!$              bms(imn) = bms(imn) + (bm(ii, jj))*sin(arg)
!!$           enddo ! end jj
!!$        enddo ! end ii
!!$
!!$        if (im .eq. 0 ) then
!!$           bmc(imn) = bmc(imn)*half
!!$           bms(imn) = bms(imn)*half
!!$        endif
!!$
!!$     enddo ! end in
!!$  enddo ! end im
!!$
!!$  bmc = bmc * two / (Nteta*Nzeta)
!!$  bms = bms * two / (Nteta*Nzeta)
!!$
!!$  imn = 0
!!$  do im = 0, mf
!!$     do in = -nf, nf
!!$        imn = imn + 1
!!$        Bmod_n(in) = Bmod_n(in) + sqrt(bmc(imn)**2 + bms(imn)**2)
!!$     enddo
!!$  enddo

!!$  open(lunit, file=trim(ext)//".Bmod", status='unknown', action='write')
!!$
!!$
!!$  !write(lunit,*      ) "#n   Bm_n"
!!$  !do in = -nf, nf
!!$     write(lunit, '(I4, ES23.15)')  in, Bmod_n(in)
!!$  !enddo
!!$
!!$  !------------------------------------
!!$
!!$  write(lunit,*      ) "#bmn bNfp  mn"
!!$  write(lunit,'(3I)' ) bmn, bNfp,  mn
!!$
!!$  write(lunit,*      ) "#------- plasma boundary------"
!!$  write(lunit,*      ) "#  n   m   Rbc   Rbs    Zbc   Zbs"
!!$  do imn = 1, bmn
!!$     write(lunit,'(2I, 4ES15.6)') bin(imn)/bNfp, bim(imn), Rbc(imn), Rbs(imn), Zbc(imn), Zbs(imn)
!!$  enddo
!!$
!!$  write(lunit,*      ) "#-------|B| harmonics----------"
!!$  write(lunit,*      ) "#   n       m        bmc        bms   "
!!$  if (mn .gt. 0) then
!!$  do imn = 1, mn
!!$     write(lunit,'(2I, 2ES15.6)') bmin(imn), bmim(imn), bmc(imn), bms(imn)
!!$  enddo
!!$  else
!!$     write(lunit,'(2I, 2ES15.6)') 0, 0, 0.0, 0.0
!!$  endif
!!$
!!$  close(lunit)

  deallocate(bm, bmin, bm_zeta)

  return

END SUBROUTINE Bmodule
  

SUBROUTINE Bmodule2
  !-------------------------------------------------------------------------------!
  ! Calculate Fourier harmonics Bmc, Bms for bm(0:Nteta, 0:Nzeta);                !
  !-------------------------------------------------------------------------------!
  use kmodule, only: zero, half, two, pi2, myid, ounit, ext, lunit, &
       Nteta, Nzeta, surf, SaveBx, SaveBy, SaveBz, Mpol, Ntor, Bmod_n, &
       bmn, bNfp, bin, bim, rbc, rbs, zbc, zbs
  implicit none
  include "mpif.h"
  !-------------------------------------------------------------------------------
  INTEGER             :: mf, nf  ! Fourier modes size
  INTEGER             :: imn=0, ii, jj, im, in, mn, astat, ierr
  REAL                :: teta, zeta, arg
  INTEGER, allocatable:: bmin(:), bmim(:)
  REAL, allocatable   :: bm(:,:), bmc(:), bms(:)
  !------------------------------------------------------------------------------- 

  if (myid .ne. 0) return

  nf = Ntor ; mf = Mpol

  mn = (mf+1)*(2*nf+1) ! (0:mf)*(-nf:nf)

  FATAL(Bmodule, mf .le. 0 .and. nf .le. 0, INVALID size for Fourier harmonics)

  SALLOCATE( bm, (0:Nteta, 0:Nzeta), zero )

  bm = sqrt(SaveBx**2 + SaveBy**2 + SaveBz**2)

  SALLOCATE( bmim, (1:mn), 0    )
  SALLOCATE( bmin, (1:mn), 0    )
  SALLOCATE( bmc , (1:mn), zero ) 
  SALLOCATE( bms , (1:mn), zero ) 


  do im = 0, mf
     do in = -nf, nf
        
        if (im==0 .and. in<0) cycle  ! skip m=0, n<0 terms
        
        imn = imn + 1
        bmin(imn) = in ; bmim(imn) = im

        do ii = 0, Nteta-1 
           teta = ( ii + half ) * pi2 / Nteta
           do jj = 0, Nzeta-1
              zeta = ( jj + half ) * pi2 / Nzeta
              arg = im*teta - in*zeta

              bmc(imn) = bmc(imn) + (bm(ii, jj))*cos(arg)
              bms(imn) = bms(imn) + (bm(ii, jj))*sin(arg)
           enddo ! end jj
        enddo ! end ii

        if (im .eq. 0 .and. in .eq. 0 ) then
           bmc(imn) = bmc(imn)*half
           bms(imn) = bms(imn)*half
        endif

     enddo ! end in
  enddo ! end im

  bmc = bmc * two / (Nteta*Nzeta)
  bms = bms * two / (Nteta*Nzeta)

  !------------------------------------
  open(lunit, file=trim(ext)//".Bmod", status='unknown', action='write')

  write(lunit,*      ) "#bmn bNfp  mn"
  write(lunit,'(3I)' ) bmn, bNfp,  mn

  write(lunit,*      ) "#------- plasma boundary------"
  write(lunit,*      ) "#  n   m   Rbc   Rbs    Zbc   Zbs"
  do imn = 1, bmn
     write(lunit,'(2I, 4ES15.6)') bin(imn)/bNfp, bim(imn), Rbc(imn), Rbs(imn), Zbc(imn), Zbs(imn)
  enddo

  write(lunit,*      ) "#-------|B| harmonics----------"
  write(lunit,*      ) "#   n       m        bmc        bms   "
  if (mn .gt. 0) then
     do imn = 1, mn
        write(lunit,'(2I, 2ES15.6)') bmin(imn), bmim(imn), bmc(imn), bms(imn)
     enddo
  else
     write(lunit,'(2I, 2ES15.6)') 0, 0, 0.0, 0.0
  endif

  close(lunit)

  deallocate(bm, bmin, bmim,  bmc, bms)

  return

END SUBROUTINE Bmodule2
