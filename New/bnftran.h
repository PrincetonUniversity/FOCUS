!title(bnftran) ! Transform the Bnormal distribution into Fourier harmonics.

!latex \briefly{There are severl subroutines in this file. The main purpose is to calculate the Fourier 
!latex harmonics of Bn on the plasma surface. Right now, the Fourier decompositions are carried out using  
!latex normal Fourier Transformation in polar coordinates. In the future, FFT and flux coordinates capability 
!latex will be included.}
!latex \calledby{\link{costfun}}

!latex \subsection{Background}
!latex The width of magnetic island is proportional to the square root of resonant perturbation amplitude. 
!latex Directly optimizing Bn harmonics (with well-chose sensitivity matrix) rather than minimizing the 
!latex surface integration would be a better idea. Especially for designing RMP coils in tokamaks.
!latex
!latex \subsection{0-order cost function}
!latex In the flux coordinate ($\psi, \theta, \pi$), the normal magnetic field perturbation can be written
!latex as,
!latex \be
!latex \ds B_n(\theta, \phi) = \vec{B} \cdot \nabla{\psi} = \sum_{m,n} \Delta_{mn} exp^{-i(m\theta - n\phi)}\ .
!latex \ee
!latex If we define the cost function {\bf H = bharm} as,
!latex \be
!latex \ds H = \frac{1}{2} \sum_{m,n} w_{mn} \ ( \Delta_{mn} - \Delta_{mn}^o )^2
!latex \ee
!latex Using trigonometric functions,
!latex \be
!latex \Delta_{mn} & = & \sqrt{ {\Delta_{mn}^c}^2 + {\Delta_{mn}^s}^2 } \\
!latex \Delta_{mn}^c &=& \frac{1}{2\pi^2} \int_0^{2\pi} \int_0^{2\pi} Bn \ \cos(m\theta-n\phi) d\theta d\phi \\
!latex \Delta_{mn}^s & = & \frac{1}{2\pi^2} \int_0^{2\pi} \int_0^{2\pi} Bn \ \sin(m\theta-n\phi) d\theta d\phi
!latex \ee
!latex 
!latex \subsection{1st-order derivatives}
!latex The derivatives of $H$ with respect to coil parameters can be calculated as,
!latex \be
!latex \pdv{H}{x} & = & \sum_{m,n} w_mn  \ ( \Delta_{mn} - \Delta_{mn}^o ) \pdv{\Delta_{mn}}{x} \\
!latex \pdv{\Delta_{mn}}{x} & = & \frac{\Delta_{mn}^c \pdv{\Delta_{mn}^c}{x} + 
!latex                                \Delta_{mn}^c \pdv{\Delta_{mn}^c}{x}}{\Delta_{mn}} \\
!latex \pdv{\Delta_{mn}^c}{x} & = & \frac{1}{2\pi^2} \int_0^{2\pi} \int_0^{2\pi} 
!latex         [\pdv{B_x}{x} n_x + [\pdv{B_y}{y} n_y + [\pdv{B_z}{z} n_z] \ \cos(m\theta-n\phi) d\theta d\phi \\
!latex \pdv{\Delta_{mn}^s}{x} & = &  \frac{1}{2\pi^2} \int_0^{2\pi} \int_0^{2\pi} 
!latex         [\pdv{B_x}{x} n_x + [\pdv{B_y}{y} n_y + [\pdv{B_z}{z} n_z] \ \sin(m\theta-n\phi) d\theta d\phi
!latex \ee

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

SUBROUTINE bmnharm( ideriv )
  !----------------------------------------------------------------------------------------
  ! calculate the bharm cost function
  !----------------------------------------------------------------------------------------
  use globals, only: zero, half, myid, Ndof, Nteta, Nzeta, surf, &
                     dBx, bharm, t1H, Bmn, wBmn, tBmn, Bmnim, Bmnin, NBmn
  implicit none
  include "mpif.h"

  INTEGER, INTENT(in) :: ideriv
  !----------------------------------------------------------------------------------------
  
  INTEGER             :: idof, ierr, astat, imn
  REAL, allocatable   :: Bmnc(:), Bmns(:), dBxc(:), dBxs(:) 

  !--------------------------initialize and allocate arrays-------------------------------- 
  SALLOCATE( Bmnc, (1:NBmn), zero )  ! 0-order Bmn_cos
  SALLOCATE( Bmns, (1:NBmn), zero )  ! 0-order Bmn_sin
  SALLOCATE( dBxc, (1:NBmn), zero )  ! temporary dBx_mn_cos
  SALLOCATE( dBxs, (1:NBmn), zero )  ! temporary dBx_mn_sin

  call bnormal(ideriv) ! calculate Bn info;

  !-------------------------------calculate H--------------------------------------------- 
  if( ideriv >= 0 ) then

     call twodft( surf(1)%bn, Bmns, Bmnc, Bmnim, Bmnin, NBmn )
!!$     if (myid == 0) then
!!$        do imn = 1, NBmn
!!$           write(*, '("n="I3,"m="I3, "Bmnc="ES12.5, "Bmns="ES12.5)') &
!!$                         Bmnin(imn), Bmnim(imn), Bmnc(imn), Bmns(imn)
!!$        enddo
!!$     endif
                 
     Bmn = sqrt( Bmns**2 + Bmnc**2 )
     bharm = half * sum( wBmn * (Bmn - tBmn)**2)
        
  endif

  !-------------------------------calculate H/x------------------------------------------------
  if ( ideriv >= 1 ) then
     ! can parallelize in Ndof direction;
     do idof = 1, Ndof
        call twodft( dBx(idof,  0:Nteta-1, 0:Nzeta-1), dBxs, dBxc, Bmnim, Bmnin, NBmn )
        t1H(idof) = sum( wBmn * (Bmn - tBmn) * (Bmnc*dBxc + Bmns*dBxs)/Bmn )
     enddo
     
  endif

  !--------------------------------------------------------------------------------------------
  DALLOCATE( Bmnc )
  DALLOCATE( Bmns )
  DALLOCATE( dBxc )
  DALLOCATE( dBxs )
  
  return

END SUBROUTINE bmnharm

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

SUBROUTINE readBmn
  !----------------------------------------------------------------------------------------
  ! read Bmn harmonics related arrays;
  ! allocate trig functions;
  !----------------------------------------------------------------------------------------
  use globals, only: zero, half, pi2, myid, ounit, runit, ext, IsQuiet, Nteta, Nzeta,  &
                     NBnf, Bnin, Bnim, NBmn, Bmnin, Bmnim, wBmn, tBmn, carg, sarg
  implicit none
  include "mpif.h"

  INTEGER  :: ii, jj, ij, imn, ierr, astat
  REAL     :: teta, zeta, arg
  LOGICAL  :: exist

  !----------------------------------------------------------------------------------------
  inquire( file=trim(ext)//".harmonics", exist=exist)  
  FATAL( readBmn, .not.exist, ext.harmonics does not exist ) 

  if (myid == 0) then
     open(runit, file=trim(ext)//".harmonics", status='old', action='read')
     read(runit,*) ! comment line;
     read(runit,*) NBmn !read dimensions
  endif
  IlBCAST( NBmn , 1, 0 )
  FATAL( readBmn, NBmn <= 0, at leasy one harmonis)
   
  SALLOCATE( Bmnim, (1:NBmn), 0 )     ! n values for Bmn;
  SALLOCATE( Bmnin, (1:NBmn), 0 )     ! m values for Bmn;
  SALLOCATE( wBmn , (1:NBmn), zero )  ! weight for different Bmn;
  SALLOCATE( tBmn , (1:NBmn), zero )  ! target Bmn values;

  if (myid == 0) then
     read(runit,*) ! comment line;
     do imn = 1, NBmn
        read(runit,*) Bmnin(imn), Bmnim(imn), wBmn(imn), tBmn(imn)
     enddo
  endif

  IlBCAST( Bmnin, NBmn, 0 )
  IlBCAST( Bmnim, NBmn, 0 )
  RlBCAST( wBmn , NBmn, 0 )
  RlBCAST( tBmn , NBmn, 0 )

  if( myid == 0 .and. IsQuiet <= 0) then
     write(ounit, *) "-----------Reading Bn harmonics weights and targets---------------------------"
     write(ounit, '("readbmn : NBmn = "I6)') NBmn
     if (IsQuiet <= -2) then
        do imn=1, NBmn
            write(ounit, '("readbmn : n = " I3 " ; m = " I3 " ; w_mn = "ES12.5" ; tBmn = "ES12.5)'), &
                 Bmnin(imn), Bmnim(imn), wBmn(imn), tBmn(imn)
         enddo
      endif
   endif

  !-------------------------store trig functions-------------------------------------------
  SALLOCATE( carg,  (1:Nteta*Nzeta, 1:NBmn), zero )
  SALLOCATE( sarg,  (1:Nteta*Nzeta, 1:NBmn), zero )
  ij = 0
  do jj = 0, Nzeta-1 ; zeta = ( jj + half ) * pi2 / Nzeta
     do ii = 0, Nteta-1 ; teta = ( ii + half ) * pi2 / Nteta
        ij = ij + 1
        do imn = 1, NBmn
           arg = Bmnim(imn) * teta - Bmnin(imn) * zeta
           carg(ij, imn) = cos(arg)
           sarg(ij, imn) = sin(arg)
        enddo
     enddo
  enddo

  return
END SUBROUTINE readBmn

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

SUBROUTINE twodft(func, hs, hc, im, in, mn)
  !-------------------------------------------------------------------------------!
  ! 2D discrete Fourier transformation;
  ! Calculate Fourier harmonics hc,hs for func(0:Nteta-1, 0:Nzeta-1);
  ! im(1:mn), in(1:mn) stores the predefined m, n values;
  ! Assuming there are no conjugate terms in im & in;
  ! carg and sarg stored the trig functions.
  ! Right now, it's using normal Fourier transforming, later FFT will be enabled.
  !-------------------------------------------------------------------------------!
  use globals, only: zero, half, two, pi2, myid, ounit, Nteta, Nzeta, carg, sarg
  implicit none
  include "mpif.h"
  !-------------------------------------------------------------------------------
  REAL   , INTENT(in ) :: func(1:Nteta*Nzeta) ! 2D array into 1D array;
  REAL   , INTENT(out) :: hc(1:mn), hs(1:mn)
  INTEGER, INTENT(in ) :: mn, im(1:mn), in(1:mn)

  INTEGER              :: m, n, imn, maxN, maxM, astat, ierr
  !------------------------------------------------------------------------------- 

  FATAL(twodft, mn < 1, invalid size for 2D Fourier transformation)

  maxN = maxval(abs(in))
  maxM = maxval(abs(im))
  FATAL(twodft, maxN >= Nzeta/2, toroidal grid resolution not enough)
  FATAL(twodft, maxM >= Nteta/2, poloidal grid resolution not enough)

  do imn = 1, mn
     m = im(imn); n = in(imn)

     hc(imn) = sum(func(1:Nteta*Nzeta) * carg(1:Nteta*Nzeta, imn))
     hs(imn) = sum(func(1:Nteta*Nzeta) * sarg(1:Nteta*Nzeta, imn))

     if (m==0 .and. n==0) then  ! for (0,0) term, times a half factor;
        hc(imn) = hc(imn)*half
        hs(imn) = hs(imn)*half
     endif

  enddo

  hc = hc * two/(Nteta*Nzeta)  ! Discretizing factor;
  hs = hs * two/(Nteta*Nzeta)  ! Discretizing factor;

  return
END SUBROUTINE twodft

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

SUBROUTINE twoift(func, hs, hc, im, in, mn)
  !-------------------------------------------------------------------------------!
  ! 2D INVERSE discrete Fourier transformation;
  ! Calculate Fourier harmonics func for given hc,hs harmonics;
  ! im(1:mn), in(1:mn) stores the predefined m, n values;
  ! carg and sarg stored the trig functions.
  ! Right now, it's using normal Fourier transforming, later FFT will be enabled.
  !-------------------------------------------------------------------------------!
  use globals, only: zero, half, two, pi2, myid, ounit, Nteta, Nzeta, carg, sarg
  implicit none
  include "mpif.h"
  !-------------------------------------------------------------------------------
  REAL   , INTENT(out) :: func(1:Nteta*Nzeta) ! 2D array into 1D array;
  REAL   , INTENT(in ) :: hc(1:mn), hs(1:mn)
  INTEGER, INTENT(in ) :: mn, im(1:mn), in(1:mn)

  INTEGER              :: itz, astat, ierr
  !------------------------------------------------------------------------------- 

  FATAL(twodft, mn < 1, invalid size for 2D Fourier transformation)

  do itz = 1, Nteta*Nzeta
     func(itz) = sum(hc(1:mn)*carg(itz, 1:mn)) + sum(hs(1:mn)*sarg(itz, 1:mn))
  enddo

  return
END SUBROUTINE twoift

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
