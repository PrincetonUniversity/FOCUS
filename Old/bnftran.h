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

SUBROUTINE bnftran(mf, nf)
  !-------------------------------------------------------------------------------!
  ! Calculate Fourier harmonics bnc,bns for bn(0:Nteta, 0:Nzeta);                 !
  ! CZHU; first version: 2017/01/11; last revised: 2017/01/11                     !
  !-------------------------------------------------------------------------------!
  use kmodule, only: zero, half, two, pi2, myid, ounit, &
       Nteta, Nzeta, surf, tbn, Nbf, bnim, bnin, bnc, bns, surf, bNfp
  implicit none
  include "mpif.h"
  !-------------------------------------------------------------------------------
  INTEGER, INTENT(in) :: mf, nf  ! predefined Fourier modes size
  INTEGER             :: imn=0, ii, jj, im, in, astat, ierr, maxN, maxM
  REAL                :: teta, zeta, arg
  !------------------------------------------------------------------------------- 

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
