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
  REAL                :: teta, zeta, arg, t1, t2
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

     Cur_Bnc = zero; Cur_Bns = zero
     t1 = MPI_Wtime()
     call twodft(surf(1)%tn, Cur_Bnc, Cur_Bns, bnim, bnin, NBnf)
     t2 = MPI_Wtime()
     write(*, *) "twodft takes ", t2-t1

     Cur_Bnc = zero; Cur_Bns = zero
     t1 = MPI_Wtime()
     call tfft(surf(1)%tn, Cur_Bnc, Cur_Bns, bnim, bnin, NBnf)
     t2 = MPI_Wtime()
     write(*, *) "tfft takes ", t2-t1

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

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item \verb+tfft+

subroutine tfft(func, hc, hs, im, in, mn)

 !use constants
 !use inputlist, only : Nfp
 !use allglobal, only : pi2nfp

  use globals, only : zero, one, pi2, ounit, myid, Nteta, Nzeta

  implicit none
  
  REAL   , INTENT(in ) :: func(0:Nteta, 0:Nzeta)
  REAL   , INTENT(out) :: hc(1:mn), hs(1:mn)
  INTEGER, INTENT(in ) :: mn, im(1:mn), in(1:mn)

  INTEGER   :: Nt, Nz, Ntz, imn, ifail, mm, nn, NFP
  REAL      :: ijreal(1:Nteta*Nzeta), ijimag(1:Nteta*Nzeta), trigm(2*Nteta), trign(2*Nzeta), &
               trigwk(2*Nteta*Nzeta), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn)
  CHARACTER :: isr
  
  LOGICAL   :: check=.false.
  INTEGER   :: jj, kk, ii
  REAL      :: jireal(1:Nteta*Nzeta), jiimag(1:Nteta*Nzeta), arg, ca, sa

  NFP = 1 ; Nt = Nteta; Nz = Nzeta

  if( check ) then ; jireal=ijreal ; jiimag=ijimag
  endif

  Ntz=Nteta*Nzeta
  ii = 0
  do jj = 0, Nteta-1
     do kk = 0, Nzeta-1
        ii = ii+1
        ijreal(ii) = func(jj, kk)
     enddo
  enddo

  TMPOUT(ijreal(1:10)) 
  TMPOUT(ijimag(1:10))
  ijimag = zero

  isr = 'I' ; ifail = 0

  call C06FUF( Nt , Nz , ijreal(1:Ntz) , ijimag(1:Ntz) , isr , trigm , trign , trigwk , ifail )

  isr = 'S' ; ifail = 0

 !call C06FUF( Nt , Nz , ijreal(1:Ntz) , ijimag(1:Ntz) , isr , trigm , trign , trigwk , ifail )
  
  !ijreal(:) = ijreal(:)/sqrt(one*Ntz) ; ijimag(:) = ijimag(:)/sqrt(one*Ntz)
  TMPOUT(ijreal(1:10))
  TMPOUT(ijimag(1:10))
  
  cfmn=zero ; sfmn=zero ; efmn=zero ; ofmn=zero
  
  do imn = 1,mn ; mm = im(imn) ; nn = in(imn) / Nfp
   
   if    ( mm.gt.0 .and. nn.gt.0 ) then
    
    efmn(imn) =   ijreal(1+(Nt-mm)+(   nn)*Nt) + ijreal(1+(   mm)+(Nz-nn)*Nt)
    ofmn(imn) =   ijimag(1+(Nt-mm)+(   nn)*Nt) - ijimag(1+(   mm)+(Nz-nn)*Nt)
    cfmn(imn) =   ijimag(1+(Nt-mm)+(   nn)*Nt) + ijimag(1+(   mm)+(Nz-nn)*Nt)
    sfmn(imn) = - ijreal(1+(Nt-mm)+(   nn)*Nt) + ijreal(1+(   mm)+(Nz-nn)*Nt) 
    
   elseif( mm.gt.0 .and. nn.eq.0 ) then
    
    efmn(imn) =   ijreal(1+(Nt-mm)           ) + ijreal(1+(   mm)           )
    ofmn(imn) =   ijimag(1+(Nt-mm)           ) - ijimag(1+(   mm)           )
    cfmn(imn) =   ijimag(1+(Nt-mm)           ) + ijimag(1+(   mm)           )
    sfmn(imn) = - ijreal(1+(Nt-mm)           ) + ijreal(1+(   mm)           )
    
   elseif( mm.gt.0 .and. nn.lt.0 ) then
    
    efmn(imn) =   ijreal(1+(Nt-mm)+(Nz+nn)*Nt) + ijreal(1+(   mm)+(  -nn)*Nt)
    ofmn(imn) =   ijimag(1+(Nt-mm)+(Nz+nn)*Nt) - ijimag(1+(   mm)+(  -nn)*Nt)
    cfmn(imn) =   ijimag(1+(Nt-mm)+(Nz+nn)*Nt) + ijimag(1+(   mm)+(  -nn)*Nt)
    sfmn(imn) = - ijreal(1+(Nt-mm)+(Nz+nn)*Nt) + ijreal(1+(   mm)+(  -nn)*Nt) 
    
   elseif( mm.eq.0 .and. nn.gt.0 ) then
    
    efmn(imn) =   ijreal(1+        (Nz-nn)*Nt) + ijreal(1+        (   nn)*Nt)
    ofmn(imn) = - ijimag(1+        (Nz-nn)*Nt) + ijimag(1+        (   nn)*Nt)
    cfmn(imn) =   ijimag(1+        (Nz-nn)*Nt) + ijimag(1+        (   nn)*Nt)
    sfmn(imn) =   ijreal(1+        (Nz-nn)*Nt) - ijreal(1+        (   nn)*Nt)
    
   elseif( mm.eq.0 .and. nn.eq.0 ) then
    
    efmn(imn) =   ijreal(1)
    ofmn(imn) =     zero
    cfmn(imn) =   ijimag(1)
    sfmn(imn) =     zero
    
   endif
   
  enddo

  hc = efmn
  hs = ofmn

  call C06GCF(ijimag, Ntz, ifail)
  call C06FUF( Nt , Nz , ijreal(1:Ntz) , ijimag(1:Ntz) , isr , trigm , trign , trigwk , ifail )
  call C06GCF(ijimag, Ntz, ifail)
  TMPOUT(ijreal(1:10)) 
  TMPOUT(ijimag(1:10))

  if( .not.check ) return
  
  ijreal(1:Ntz)=zero ; ijimag(1:Ntz)=zero
  
  do jj=0,Nt-1
   do kk=0,Nz-1
    do imn=1,mn ; arg=im(imn)*jj*pi2/Nt-in(imn)*kk*(pi2/nfp)/Nz ; ca=cos(arg) ; sa=sin(arg)
     
     ijreal(1+jj+kk*Nt) = ijreal(1+jj+kk*Nt) + efmn(imn)*ca + ofmn(imn)*sa
     ijimag(1+jj+kk*Nt) = ijimag(1+jj+kk*Nt) + cfmn(imn)*ca + sfmn(imn)*sa
     
    enddo
   enddo
  enddo
  
  write(ounit,'("tfft : Fourier reconstruction error = "2es15.5)')sqrt(sum((ijreal-jireal)**2)/Ntz),sqrt(sum((ijimag-jiimag)**2)/Ntz)

  return

end subroutine tfft

