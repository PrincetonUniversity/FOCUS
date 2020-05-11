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
!latex In the flux coordinate ($\psi, \theta, \phi$), the normal magnetic field perturbation can be written
!latex as,
!latex \be
!latex \ds B_n(\theta,\phi) = \vec{B} \cdot \nabla{\psi} = \sum_{m,n} \Delta_{mn} e^{-i(m\theta - n\phi)}\ .
!latex \ee
!latex If we define the cost function {\bf H = bharm} as,
!latex \be
!latex \ds H & = & \frac{1}{2} \sum_{m,n} w_{mn} \ ( \Delta_{mn} - \Delta_{mn}^o )^2 \nonumber \\
!latex       & = & \frac{1}{2} \sum_{m,n} w_{mn} \left [ ( \Delta^c_{mn} - \Delta^{c,o}_{mn} )^2 
!latex                                                 + ( \Delta^s_{mn} - \Delta^{s,o}_{mn} )^2 \right ]
!latex \ee
!latex Using trigonometric functions,
!latex \be
!      \Delta_{mn} & = & \sqrt{ {\Delta_{mn}^c}^2 + {\Delta_{mn}^s}^2 } \\
!latex \Delta_{mn}^c &=& \frac{1}{2\pi^2} \int_0^{2\pi} \int_0^{2\pi} Bn \ \cos(m\theta-n\phi) d\t d\phi \\
!latex \Delta_{mn}^s & = & \frac{1}{2\pi^2} \int_0^{2\pi} \int_0^{2\pi} Bn \ \sin(m\theta-n\phi) d\t d\phi
!latex \ee
!latex 
!latex \subsection{1st-order derivatives}
!latex The derivatives of $H$ with respect to coil parameters can be calculated as,
!latex \be
!latex \pdv{H}{x} & = & \sum_{m,n} w_{mn} \left [ ( \Delta^c_{mn}-\Delta^{c,o}_{mn} ) \pdv{\Delta^c_{mn}}{x}
!latex                \  + \ ( \Delta^s_{mn}-\Delta^{s,o}_{mn} ) \pdv{\Delta^s_{mn}}{x} \right ] \\

!latex \pdv{\Delta_{mn}^c}{x} & = & \frac{1}{2\pi^2} \int_0^{2\pi} \int_0^{2\pi} 
!latex         \left [\pdv{B_x}{x} n_x + \pdv{B_y}{y} n_y + \pdv{B_z}{z} n_z \right ] 
!latex         \ \cos(m\theta-n\phi) d\theta d\phi \\

!latex \pdv{\Delta_{mn}^s}{x} & = &  \frac{1}{2\pi^2} \int_0^{2\pi} \int_0^{2\pi} 
!latex         \left [\pdv{B_x}{x} n_x + \pdv{B_y}{y} n_y + \pdv{B_z}{z} n_z \right ] 
!latex         \ \sin(m\theta-n\phi) d\theta d\phi
!latex \ee

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
!!$
!!$SUBROUTINE bmnharm( ideriv )
!!$  !----------------------------------------------------------------------------------------
!!$  ! calculate the bharm cost function
!!$  !----------------------------------------------------------------------------------------
!!$  use globals, only: zero, half, myid, Ndof, Nteta, Nzeta, surf, &
!!$                     bn, dB, bharm, t1H, Bmnc, Bmns, wBmn, tBmnc, tBmns, Bmnim, Bmnin, NBmn
!!$  implicit none
!!$  include "mpif.h"
!!$
!!$  INTEGER, INTENT(in) :: ideriv
!!$  !----------------------------------------------------------------------------------------
!!$  
!!$  INTEGER             :: idof, ierr, astat, imn
!!$  REAL, allocatable   :: dBc(:), dBs(:) 
!!$
!!$  !--------------------------initialize and allocate arrays-------------------------------- 
!!$
!!$  SALLOCATE( dBc, (1:NBmn), zero )  ! temporary dB_mn_cos
!!$  SALLOCATE( dBs, (1:NBmn), zero )  ! temporary dB_mn_sin
!!$
!!$  call bnormal(ideriv) ! calculate Bn info;
!!$
!!$  !-------------------------------calculate H--------------------------------------------- 
!!$  if( ideriv >= 0 ) then
!!$
!!$    !call twodft( surf(1)%bn, Bmns, Bmnc, Bmnim, Bmnin, NBmn ) ! total Bn
!!$     call twodft(         bn, Bmns, Bmnc, Bmnim, Bmnin, NBmn ) ! Bn from coils
!!$!!!$     if (myid == 0) then
!!$!!!$        do imn = 1, NBmn
!!$!!!$           write(*, '("n="I3,"m="I3, "Bmnc="ES12.5, "Bmns="ES12.5)') &
!!$!!!$                         Bmnin(imn), Bmnim(imn), Bmnc(imn), Bmns(imn)
!!$!!!$        enddo
!!$!!!$     endif
!!$                 
!!$     bharm = half * sum( wBmn * ((Bmnc - tBmnc)**2 + (Bmns - tBmns)**2) )
!!$        
!!$  endif
!!$
!!$  !-------------------------------calculate H/x------------------------------------------------
!!$  if ( ideriv >= 1 ) then
!!$     ! can parallelize in Ndof direction;
!!$     do idof = 1, Ndof
!!$        call twodft( dB(idof,  0:Nteta-1, 0:Nzeta-1), dBs, dBc, Bmnim, Bmnin, NBmn )
!!$        t1H(idof) = sum( wBmn * ( (Bmnc - tBmnc)*dBc + (Bmns - tBmns)*dBs ) )
!!$     enddo
!!$     
!!$  endif
!!$
!!$  !--------------------------------------------------------------------------------------------
!!$
!!$  DALLOCATE( dBc )
!!$  DALLOCATE( dBs )
!!$  
!!$  return
!!$
!!$END SUBROUTINE bmnharm

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
module bharm_mod
  ! contains some common variables used in subroutine bnormal
  ! allocating once and re-using them will save allocation time
  use globals, only : dp
  implicit none

  ! 0-order
  ! none for now; in future, others should be moved to here. 03/30/2019
  ! 1st-order
  REAL, allocatable :: dBc(:), dBs(:)

end module bharm_mod
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

SUBROUTINE readBmn
  !----------------------------------------------------------------------------------------
  ! read Bmn harmonics related arrays;
  ! allocate trig functions;
  !----------------------------------------------------------------------------------------
  use globals, only: dp, zero, half, pi2, myid, ounit, runit, ext, IsQuiet, Nteta, Nzeta, Nfp, &
                     NBmn, Bmnin, Bmnim, wBmn, tBmnc, tBmns, carg, sarg, case_bnormal, &
                     input_harm, bharm_jsurf, surf, plasma, MPI_COMM_FOCUS, bharm_factor, &
                     two, weight_bnorm, machprec
  use bharm_mod
  implicit none
  include "mpif.h"

  INTEGER  :: ii, jj, ij, imn, ierr, astat, isurf
  REAL     :: teta, zeta, arg
  LOGICAL  :: exist

  !----------------------------------------------------------------------------------------
  isurf = plasma
  inquire( file=trim(input_harm), exist=exist)  
  FATAL( readBmn, .not.exist, ext.harmonics does not exist ) 

  if (myid == 0) then
     open(runit, file=trim(input_harm), status='old', action='read')
     read(runit,*) ! comment line;
     read(runit,*) NBmn !read dimensions
  endif
  IlBCAST( NBmn , 1, 0 )
  FATAL( readBmn, NBmn <= 0, at leasy one harmonis)
   
  SALLOCATE( Bmnim, (1:NBmn), 0 )     ! n values for Bmn;
  SALLOCATE( Bmnin, (1:NBmn), 0 )     ! m values for Bmn;
  SALLOCATE( wBmn , (1:NBmn), zero )  ! weight for different Bmn;
  SALLOCATE( tBmnc, (1:NBmn), zero )  ! target Bmn cos values;
  SALLOCATE( tBmns, (1:NBmn), zero )  ! target Bmn sin values;

  if (myid == 0) then
     read(runit,*) ! comment line;
     do imn = 1, NBmn
        read(runit,*) Bmnin(imn), Bmnim(imn),tBmnc(imn), tBmns(imn), wBmn(imn)
     enddo
  endif

  IlBCAST( Bmnin, NBmn, 0 )
  IlBCAST( Bmnim, NBmn, 0 )
  RlBCAST( wBmn , NBmn, 0 )
  RlBCAST( tBmnc, NBmn, 0 )
  RlBCAST( tBmns, NBmn, 0 )

  if( myid == 0 .and. IsQuiet <= 0) then
     write(ounit, *) "-----------Reading Bn harmonics weights and targets---------------------------"
     write(ounit, '("readbmn : NBmn = "I6)') NBmn
     if (IsQuiet <= -2) then
        do imn=1, NBmn
            write(ounit, '("readbmn : n = " I3 " ; m = " I3 " ; tBmnc = "ES12.5 " ; tBmns = "ES12.5 &
              & " ; w_mn = "ES12.5)') Bmnin(imn), Bmnim(imn), tBmnc(imn),tBmns(imn), wBmn(imn)
         enddo
      endif
      close(runit)
      write(ounit, '("******* : case_bnormal has been reset to 0, since Bn harmonics is turned on.")')
   endif
   case_bnormal = 0

   !-------------------------store trig functions-------------------------------------------
   SALLOCATE( carg,  (1:Nteta*Nzeta, 1:NBmn), zero )
   SALLOCATE( sarg,  (1:Nteta*Nzeta, 1:NBmn), zero )

   Bmnin(1:NBmn) = Bmnin(1:NBmn) * surf(isurf)%Nfp

   ij = 0
   ! the same as in rdsurf.h
   do jj = 0, Nzeta-1
      zeta = ( jj + half ) * pi2 / surf(isurf)%Nzeta
      do ii = 0, Nteta-1
         teta = ( ii + half ) * pi2 / surf(isurf)%Nteta
         ij = ij + 1
         do imn = 1, NBmn
            arg = Bmnim(imn) * teta - Bmnin(imn) * zeta
            carg(ij, imn) = cos(arg)
            sarg(ij, imn) = sin(arg)
         enddo
         ! Additional weighting
         if (bharm_jsurf == 0) then
            continue
         else if (bharm_jsurf == 1) then ! Bn * dA
            carg(ij, 1:NBmn) = carg(ij, 1:NBmn) * (surf(isurf)%ds(ii, jj))
            sarg(ij, 1:NBmn) = sarg(ij, 1:NBmn) * (surf(isurf)%ds(ii, jj))
         else if ( bharm_jsurf == 2) then ! Bn * sqrt(dA)
            carg(ij, 1:NBmn) = carg(ij, 1:NBmn) * sqrt(surf(isurf)%ds(ii, jj))
            sarg(ij, 1:NBmn) = sarg(ij, 1:NBmn) * sqrt(surf(isurf)%ds(ii, jj))
         end if
      enddo
   enddo

   SALLOCATE( dBc, (1:NBmn), zero )  ! dB_mn_cos
   SALLOCATE( dBs, (1:NBmn), zero )  ! dB_mn_sin

   ! calculate the additional normalized factor
   bharm_factor = two
   if (bharm_jsurf == 0) then
      continue
   else if (bharm_jsurf == 1) then 
      bharm_factor = two * pi2**2 / surf(isurf)%area 
   else if (bharm_jsurf == 2) then
      bharm_factor = two * pi2 / sqrt(surf(isurf)%area )
   end if

   ! assign plasma Bn if bnorm is optimized
   if (weight_bnorm > machprec) then
      if (myid==0) write(ounit, '(8X,": Plasma Bn is replaced by the target harmonics.")') 
      call twoift(surf(plasma)%pb(0:Nteta-1, 0:Nzeta-1), tBmns, tBmnc, Bmnim, Bmnin, NBmn )
   end if 
  return
END SUBROUTINE readBmn

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
  use globals, only: dp, zero, half, two, pi2, myid, ounit, bharm_factor, &
       Nteta, Nzeta, carg, sarg, bharm_jsurf, surf, plasma, MPI_COMM_FOCUS 
  implicit none
  include "mpif.h"
  !-------------------------------------------------------------------------------
  REAL   , INTENT(in ) :: func(1:Nteta*Nzeta) ! 2D array into 1D array;
  REAL   , INTENT(out) :: hc(1:mn), hs(1:mn)
  INTEGER, INTENT(in ) :: mn, im(1:mn), in(1:mn)

  INTEGER              :: m, n, imn, maxN, maxM, astat, ierr, isurf
  !------------------------------------------------------------------------------- 

  FATAL(twodft, mn < 1, invalid size for 2D Fourier transformation)

  isurf = plasma
  maxN = maxval(abs(in))
  maxM = maxval(abs(im))
  FATAL(twodft, maxN >= Nzeta/2, toroidal grid resolution not enough)
  FATAL(twodft, maxM >= Nteta/2, poloidal grid resolution not enough)

  do imn = 1, mn
     m = im(imn); n = in(imn)

     hc(imn) = sum(func(1:Nteta*Nzeta) * carg(1:Nteta*Nzeta, imn))
     hs(imn) = sum(func(1:Nteta*Nzeta) * sarg(1:Nteta*Nzeta, imn))

     if (m==0 .and. n==0) then  ! for (0,0) term, times a half factor;
     ! if (m==0) then  ! for (0,0) term, times a half factor;
        hc(imn) = hc(imn)*half
        hs(imn) = hs(imn)*half
     endif

  enddo

  hc = hc * two/(Nteta*Nzeta)  ! Discretizing factor;
  hs = hs * two/(Nteta*Nzeta)  ! Discretizing factor;

  ! Additional weighting
  hc = hc * bharm_factor
  hs = hs * bharm_factor

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
  use globals, only: dp, zero, half, two, pi2, myid, ounit, Nteta, Nzeta, &
   carg, sarg, MPI_COMM_FOCUS, bharm_factor, bharm_jsurf, plasma, surf
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

  ! divide by the jacobians
  if (bharm_jsurf==0) then
      continue
  else if (bharm_jsurf==1) then
      func = func / reshape(surf(plasma)%ds*surf(plasma)%ds, (/Nteta*Nzeta/))
  else if (bharm_jsurf==2) then
      func = func / reshape(surf(plasma)%ds, (/Nteta*Nzeta/))
  end if 

  ! Additional weighting
  func = func / bharm_factor

  return
END SUBROUTINE twoift

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

SUBROUTINE saveBmn
  !----------------------------------------------------------------------------------------
  ! save the present Bmn harmonics in iBmnc and iBmns;
  !----------------------------------------------------------------------------------------
  use globals, only: dp, zero, ierr, astat, myid, machprec, MPI_COMM_FOCUS, &
   & weight_Bharm, NBmn, Bmnc, Bmns, iBmnc, iBmns
  implicit none
  include "mpif.h"

  if (weight_bharm > machprec) then

     FATAL( saveBmn, .not. allocated( Bmnc), you should allocate  Bmnc first. )
     FATAL( saveBmn, .not. allocated(iBmnc), you should allocate iBmnc first. )

     iBmnc = Bmnc
     iBmns = Bmns

  endif

  return

END SUBROUTINE saveBmn
