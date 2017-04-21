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

SUBROUTINE BnFTran
  !-------------------------------------------------------------------------------!
  ! Calculate Fourier harmonics bnc,bns for bn(0:Nteta, 0:Nzeta);                 !
  ! CZHU; first version: 2017/01/11; last revised: 2017/04/06                     !
  !-------------------------------------------------------------------------------!
  use globals, only: zero, half, two, pi2, myid, ounit, &
       Nteta, Nzeta, surf, bn, NBnf, Bnim, Bnin, Cur_Bnc, Cur_Bns
  implicit none
  include "mpif.h"
  !-------------------------------------------------------------------------------
  INTEGER,parameter   :: mf=20, nf = 20  ! predefined Fourier modes size
  INTEGER             :: imn=0, ii, jj, im, in, astat, ierr, maxN, maxM
  REAL                :: teta, zeta, arg
  !------------------------------------------------------------------------------- 

  FATAL(bnftran, mf .le. 0 .and. nf .le. 0, INVALID size for Fourier harmonics)

  if( allocated( Cur_Bnc ) ) then
     DALLOCATE( Cur_Bnc )
     DALLOCATE( Cur_Bns )
  endif

  if ( NBnf >= 1 ) then

     maxN = maxval(abs(bnin))
     maxM = maxval(abs(bnim))

     FATAL(bnftran, maxN .ge. Nzeta/2, toroidal grid resolution not enough)
     FATAL(bnftran, maxM .ge. Nteta/2, poloidal grid resolution not enough)

     SALLOCATE( Cur_Bnc, (1:NBnf), zero )
     SALLOCATE( Cur_Bns, (1:NBnf), zero )

     do imn = 1, NBnf

        in = bnin(imn)
        im = bnim(imn)

        do ii = 0, Nteta-1 
           teta = ( ii + half ) * pi2 / Nteta
           do jj = 0, Nzeta-1
              zeta = ( jj + half ) * pi2 / Nzeta
              arg = im*teta - in*zeta
              Cur_Bnc(imn) = Cur_Bnc(imn) + surf(1)%bn(ii,jj)*cos(arg)
              Cur_Bns(imn) = Cur_Bns(imn) + surf(1)%bn(ii,jj)*sin(arg)
           enddo ! end jj
        enddo ! end ii

        if (im .eq. 0 .and. in .eq. 0) then
           Cur_Bnc(imn) = Cur_Bnc(imn)*half
           Cur_Bns(imn) = Cur_Bns(imn)*half
        endif

     enddo

     Cur_Bnc = Cur_Bnc * two / (Nteta*Nzeta)
     Cur_Bns = Cur_Bns * two / (Nteta*Nzeta)

  else

     NBnf = (mf+1)*(2*nf+1) ! (0:mf)*(-nf:nf)

     SALLOCATE( bnim, (1:NBnf), 0    )
     SALLOCATE( bnin, (1:NBnf), 0    )
     SALLOCATE( Cur_Bnc , (1:NBnf), zero )
     SALLOCATE( Cur_Bns , (1:NBnf), zero )  

     do in = -nf, nf
        do im = 0, mf

           imn = imn + 1
           bnin(imn) = in ; bnim(imn) = im

           do ii = 0, Nteta-1 
              teta = ( ii + half ) * pi2 / Nteta
              do jj = 0, Nzeta-1
                 zeta = ( jj + half ) * pi2 / Nzeta
                 arg = im*teta - in*zeta

                 Cur_Bnc(imn) = Cur_Bnc(imn) + (surf(1)%bn(ii, jj))*cos(arg)
                 Cur_Bns(imn) = Cur_Bns(imn) + (surf(1)%bn(ii, jj))*sin(arg)
              enddo ! end jj
           enddo ! end ii

           if (im .eq. 0 ) then
              Cur_Bnc(imn) = Cur_Bnc(imn)*half
              Cur_Bns(imn) = Cur_Bns(imn)*half
           endif

        enddo ! end im
     enddo ! end in

     Cur_Bnc = Cur_Bnc * two / (Nteta*Nzeta)
     Cur_Bns = Cur_Bns * two / (Nteta*Nzeta)
  endif

  call write_plasma
  
  return
END SUBROUTINE bnftran

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

SUBROUTINE write_plasma
!-------------------------------------------------------------------------------!
! write down the unpdated plasma boundary information;                          !
! CZHU; first version: 2017/01/11; last revised: 2017/01/11                     !
!-------------------------------------------------------------------------------!
  use globals, only : myid, wunit, ext, &
                      Nfou, Nfp, NBnf, bim, bin, Bnim, Bnin, Rbc, Rbs, Zbc, Zbs, Cur_Bnc, Cur_Bns
  
  implicit none  
  include "mpif.h"
!-------------------------------------------------------------------------------
  INTEGER :: imn
!-------------------------------------------------------------------------------
  if(myid .ne. 0) return

  open(wunit, file=trim(ext)//".boundary", status='unknown', action='write')

  write(wunit,*      ) "#Nfou Nfp NBnf"
  write(wunit,'(3I)' ) Nfou, Nfp, NBnf

  write(wunit,*      ) "#------- plasma boundary------"
  write(wunit,*      ) "#  n   m   Rbc   Rbs    Zbc   Zbs"
  do imn = 1, Nfou
     write(wunit,'(2I, 4ES15.6)') bin(imn), bim(imn), Rbc(imn), Rbs(imn), Zbc(imn), Zbs(imn)
  enddo

  write(wunit,*      ) "#-------Bn harmonics----------"
  write(wunit,*      ) "#  n  m  bnc   bns"
  if (NBnf .gt. 0) then
  do imn = 1, NBnf
     write(wunit,'(2I, 2ES15.6)') Bnin(imn), Bnim(imn), Cur_Bnc(imn), Cur_Bns(imn)
  enddo
  else
     write(wunit,'(2I, 2ES15.6)') 0, 0, 0.0, 0.0
  endif

  close(wunit)
END SUBROUTINE write_plasma

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
