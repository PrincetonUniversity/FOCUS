
!title (boundary) ! The plasma boundary is read from file 

!latex \briefly{A Fourier representation for the plasma boundary is read from file }

!latex \calledby{\link{xfocus}}
!latex \calls{\link{}}

!latex \section{General representation (stellarator)}
!latex \subsection{overview}
!latex The general representation for plasma boundary is in \subroutine{generic}. The basic fomulation
!latex is
!latex \be
!latex \ds R &= \sum R_{mn}^c \, \cos(m\t - n\z) + R_{mn}^s \, \sin(m\t - n\z) \nonumber \\
!latex \ds Z &= \sum Z_{mn}^c \, \cos(m\t - n\z) + Z_{mn}^s \, \sin(m\t - n\z) \nonumber
!latex \ee
!latex Usually, if we imply stellarator symmetry, then $R_{mn}^s$ and $Z_{mn}^c$ would be zero.
!latex 
!latex The positive driection for poloidal angle $\t$ is \red{counterclockwise} and for toroidal angle is also
!latex \red{counterclockwise} from the top view. The positive surface normal should be pointed outwards.
!latex \subsection{Variables}
!latex The Fourier harmonics of the plasma boundary are  reqired in \verb+plasma.boundary+, 
!latex and the format of this file is as follows: 
!latex \begin{raw}
!latex Nfou       ! integer: number of Fourier harmonics for the plasma boundary;
!latex Nfp        ! integer: number of field periodicity;
!latex NBnf       ! integer: number of Fourier harmonics for Bn;
!latex ---------------------------------------------------------
!latex bim(1:bmn) ! integer: poloidal mode identification;
!latex bin(1:bmn) ! integer: toroidal mode identification;
!latex Bnim(1:bmn)! integer: poloidal mode identification, for Bn;
!latex Bnin(1:bmn)! integer: toroidal mode identification, for Bn;
!latex ---------------------------------------------------------
!latex Rbc(1:bmn) ! real   : cylindrical R cosine harmonics;
!latex Rbs(1:bmn) ! real   : cylindrical R   sine harmonics;
!latex Zbc(1:bmn) ! real   : cylindrical Z cosine harmonics;
!latex Zbs(1:bmn) ! real   : cylindrical Z   sine harmonics;
!latex Bns(1:nbf) ! real   : B normal sin harmonics;
!latex Bnc(1:nbf) ! real   : B normal cos harmonics; \end{raw}
!latex Note that immediately after reading (and broadcasting) 
!latex \verb+bin+, the field periodicity factor is included, i.e. \verb+bin = bin * Nfp+.
!latex \subsection{Sample file}
!latex Example of the plasma.boundary file:
!latex  { \begin{raw}
!latex #Nfou Nfp NBnf
!latex 4 2 1
!latex #plasma boundary
!latex # n m Rbc Rbs Zbc Zbs
!latex 0 0  3.00 0.0 0.0  0.00
!latex 0 1  0.30 0.0 0.0 -0.30
!latex 1 0  0.00 0.0 0.0 -0.06
!latex 1 1 -0.06 0.0 0.0 -0.06
!latex #Bn harmonics
!latex # n m bnc bns
!latex 0 0 0.0 0.0
!latex \end{raw}
!latex }
!latex \section{Knotran}
!latex The input surface file for knotrans is descriped in \code{knotxx}.
!latex \section{Tokamak}
!latex This part is reserved for later development of the interface for tokamaks.

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine fousurf
  
  use globals, only : dp, zero, half, pi2, myid, ounit, runit, input_surf, IsQuiet, IsSymmetric, &
                      Nfou, Nfp, NBnf, bim, bin, Bnim, Bnin, Rbc, Rbs, Zbc, Zbs, Bnc, Bns,  &
                      Nteta, Nzeta, surf, Npc, discretefactor, Nfp_raw, cosnfp, sinnfp, &
                      half_shift
  
  implicit none
  
  include "mpif.h"
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOGICAL :: exist
  INTEGER :: iosta, astat, ierr, ii, jj, imn, ip, symmetry
  REAL    :: RR(0:2), ZZ(0:2), szeta, czeta, xx(1:3), xt(1:3), xz(1:3), ds(1:3), &
             teta, zeta, arg, dd, shift
  
  !-------------read plasma.boundary---------------------------------------------------------------------  
  inquire( file=trim(input_surf), exist=exist)  
  FATAL( surface, .not.exist, plasma.boundary does not exist ) 
  if( myid == 0 ) then
     open(runit, file=trim(input_surf), status='old', action='read')
     read(runit,*) !empty line
     read(runit,*) Nfou, Nfp, NBnf !read dimensions
  endif
  
  !Broadcast the values
  IlBCAST( Nfou , 1, 0 )
  IlBCAST( Nfp  , 1, 0 )
  IlBCAST( NBnf , 1, 0 )  
  FATAL( surface, Nfou <= 0, invalid )
  FATAL( surface, Nfp  <= 0, invalid )
  FATAL( surface, NBnf <  0, invalid )

  !allocate arrays
  SALLOCATE( bim, (1:Nfou), 0 )
  SALLOCATE( bin, (1:Nfou), 0 )  
  SALLOCATE( Rbc, (1:Nfou), zero )
  SALLOCATE( Rbs, (1:Nfou), zero )
  SALLOCATE( Zbc, (1:Nfou), zero )
  SALLOCATE( Zbs, (1:Nfou), zero )

  if( myid == 0 ) then
   read(runit,*) !empty line
   read(runit,*) !empty line
   do imn = 1, Nfou
      read(runit,*) bin(imn), bim(imn), Rbc(imn), Rbs(imn), Zbc(imn), Zbs(imn)
   enddo
  endif  

  IlBCAST( bim(1:Nfou), Nfou, 0 )
  IlBCAST( bin(1:Nfou), Nfou, 0 )
 
  bin(1:Nfou) = bin(1:Nfou) * Nfp  !The full plasma;
     
  RlBCAST( Rbc(1:Nfou), Nfou, 0 )
  RlBCAST( Rbs(1:Nfou), Nfou, 0 )
  RlBCAST( Zbc(1:Nfou), Nfou, 0 )
  RlBCAST( Zbs(1:Nfou), Nfou, 0 )

  !read Bnormal ditributions
  if( NBnf  >   0) then
     SALLOCATE( Bnim, (1:NBnf), 0    )
     SALLOCATE( Bnin, (1:NBnf), 0    )
     SALLOCATE( Bnc , (1:NBnf), zero )
     SALLOCATE( Bns , (1:NBnf), zero )

     if( myid == 0 ) then
        read(runit,*) !empty line
        read(runit,*) !empty line
        do imn = 1, NBnf 
           read(runit,*) Bnin(imn), Bnim(imn), Bnc(imn), Bns(imn)
        enddo
     endif

     IlBCAST( Bnim(1:NBnf), NBnf, 0 )
     IlBCAST( Bnin(1:NBnf), NBnf, 0 )

     !if (IsSymmetric  ==  0)
     Bnin(1:NBnf) = Bnin(1:NBnf) * Nfp ! Disarde periodicity;
     ! This should be consistent with bnftran; Before fully constructed the stellarator symmetry,
     ! it's turned off;
     
     RlBCAST( Bnc(1:NBnf) , NBnf, 0 )
     RlBCAST( Bns(1:NBnf) , NBnf, 0 )
  endif

  if( myid == 0 ) close(runit,iostat=iosta)
  
  IlBCAST( iosta, 1, 0 )
  
  FATAL( surface, iosta.ne.0, error closing plasma.boundary )
  
  !-------------output for check-------------------------------------------------------------------------
  if( myid == 0 .and. IsQuiet <= 0) then
     write(ounit, *) "-----------Reading surface-----------------------------------"
     write(ounit, '("surface : Plasma boundary will be discretized in "I6" X "I6" elements.")') Nteta, Nzeta
     write(ounit, '(8X": Nfou = " I06 " ; Nfp = " I06 " ; NBnf = " I06 " ;" )') Nfou, Nfp, NBnf
  endif

  if( myid == 0 .and. IsQuiet <= -2) then !very detailed output;
     write(ounit,'("        : " 10x " : bim ="10i13   )') bim(1:Nfou)
     write(ounit,'("        : " 10x " : bin ="10i13   )') bin(1:Nfou)
     write(ounit,'("        : " 10x " : Rbc ="10es13.5)') Rbc(1:Nfou)
     write(ounit,'("        : " 10x " : Rbs ="10es13.5)') Rbs(1:Nfou)
     write(ounit,'("        : " 10x " : Zbc ="10es13.5)') Zbc(1:Nfou)
     write(ounit,'("        : " 10x " : Zbs ="10es13.5)') Zbs(1:Nfou)
     if(Nbnf > 0) then
        write(ounit,'("        : " 10x " : Bnim ="10i13  )') Bnim(1:NBnf)
        write(ounit,'("        : " 10x " : Bnin ="10i13  )') Bnin(1:NBnf)
        write(ounit,'("        : " 10x " : Bnc ="10es13.5)') Bnc (1:NBnf)
        write(ounit,'("        : " 10x " : Bns ="10es13.5)') Bns (1:NBnf)
     endif
  endif

  !-------------discretize surface data------------------------------------------------------------------  
  
  Nfp_raw = Nfp ! save the raw value of Nfp
  select case (IsSymmetric)
  case ( 0 )
     Nfp = 1                          !reset Nfp to 1;
     Npc = Nfp_raw                    !number of coils periodicity
     symmetry = 0
  case ( 1 )                          !plasma and coil periodicity enabled;
     Npc = Nfp
     symmetry = 0
  case ( 2 )                          ! stellarator symmetry enforced;
     Npc = Nfp
     symmetry = 1     
  end select
  ! discretefactor = discretefactor/Nfp

  SALLOCATE( cosnfp, (1:Npc), zero )
  SALLOCATE( sinnfp, (1:Npc), zero )  
  
  do ip = 1, Npc
     cosnfp(ip) = cos((ip-1)*pi2/Npc)
     sinnfp(ip) = sin((ip-1)*pi2/Npc)
  enddo

  allocate( surf(1:1) ) ! can allow for myltiple plasma boundaries 
                        ! if multiple currents are allowed; 14 Apr 16;
  
  surf(1)%Nteta = Nteta ! not used yet; used for multiple surfaces; 20170307;
  surf(1)%Nzeta = Nzeta * Nfp * 2**symmetry ! the total number from [0, 2pi]
  
  SALLOCATE( surf(1)%xx, (0:surf(1)%Nteta-1,0:surf(1)%Nzeta-1), zero ) !x coordinates;
  SALLOCATE( surf(1)%yy, (0:surf(1)%Nteta-1,0:surf(1)%Nzeta-1), zero ) !y coordinates
  SALLOCATE( surf(1)%zz, (0:surf(1)%Nteta-1,0:surf(1)%Nzeta-1), zero ) !z coordinates
  SALLOCATE( surf(1)%nx, (0:surf(1)%Nteta-1,0:surf(1)%Nzeta-1), zero ) !unit nx;
  SALLOCATE( surf(1)%ny, (0:surf(1)%Nteta-1,0:surf(1)%Nzeta-1), zero ) !unit ny;
  SALLOCATE( surf(1)%nz, (0:surf(1)%Nteta-1,0:surf(1)%Nzeta-1), zero ) !unit nz;
  SALLOCATE( surf(1)%ds, (0:surf(1)%Nteta-1,0:surf(1)%Nzeta-1), zero ) !jacobian;
  SALLOCATE( surf(1)%xt, (0:surf(1)%Nteta-1,0:surf(1)%Nzeta-1), zero ) !dx/dtheta;
  SALLOCATE( surf(1)%yt, (0:surf(1)%Nteta-1,0:surf(1)%Nzeta-1), zero ) !dy/dtheta;
  SALLOCATE( surf(1)%zt, (0:surf(1)%Nteta-1,0:surf(1)%Nzeta-1), zero ) !dz/dtheta;
  SALLOCATE( surf(1)%pb, (0:        Nteta-1,0:        Nzeta-1), zero ) !target Bn;
  SALLOCATE( surf(1)%xp, (0:surf(1)%Nteta-1,0:surf(1)%Nzeta-1), zero ) !dx/dzeta;
  SALLOCATE( surf(1)%yp, (0:surf(1)%Nteta-1,0:surf(1)%Nzeta-1), zero ) !dy/dzeta;
  SALLOCATE( surf(1)%zp, (0:surf(1)%Nteta-1,0:surf(1)%Nzeta-1), zero ) !dz/dzeta;
  
  surf(1)%vol = zero  ! volume enclosed by plasma boundary
 
  discretefactor = (pi2/surf(1)%Nteta) * (pi2/surf(1)%Nzeta)

  if (half_shift) then
     shift = half
  else 
     shift = zero
     if(myid.eq.0) write(ounit, '(8X": half-shift in surface evaluation is turned off." )')
  endif
    
! The center point value was used to discretize grid;
  do ii = 0, surf(1)%Nteta-1
     teta = ( ii + shift ) * pi2 / surf(1)%Nteta
     do jj = 0, surf(1)%Nzeta-1
        zeta = ( jj + shift ) * pi2 / surf(1)%Nzeta
        RR(0:2) = zero ; ZZ(0:2) = zero

        do imn = 1, Nfou ; arg = bim(imn) * teta - bin(imn) * zeta

           RR(0) =  RR(0) +     Rbc(imn) * cos(arg) + Rbs(imn) * sin(arg)
           ZZ(0) =  ZZ(0) +     Zbc(imn) * cos(arg) + Zbs(imn) * sin(arg)

           RR(1) =  RR(1) + ( - Rbc(imn) * sin(arg) + Rbs(imn) * cos(arg) ) * bim(imn)
           ZZ(1) =  ZZ(1) + ( - Zbc(imn) * sin(arg) + Zbs(imn) * cos(arg) ) * bim(imn)

           RR(2) =  RR(2) - ( - Rbc(imn) * sin(arg) + Rbs(imn) * cos(arg) ) * bin(imn)
           ZZ(2) =  ZZ(2) - ( - Zbc(imn) * sin(arg) + Zbs(imn) * cos(arg) ) * bin(imn)

        enddo ! end of do imn; 30 Oct 15;

        szeta = sin(zeta)
        czeta = cos(zeta)

        xx(1:3) = (/   RR(0) * czeta,   RR(0) * szeta, ZZ(0) /)
        xt(1:3) = (/   RR(1) * czeta,   RR(1) * szeta, ZZ(1) /)
        xz(1:3) = (/   RR(2) * czeta,   RR(2) * szeta, ZZ(2) /) + (/ - RR(0) * szeta,   RR(0) * czeta, zero  /)

        ds(1:3) = -(/ xt(2) * xz(3) - xt(3) * xz(2), & ! minus sign for theta counterclockwise direction;
             xt(3) * xz(1) - xt(1) * xz(3), &
             xt(1) * xz(2) - xt(2) * xz(1) /)

        dd = sqrt( sum( ds(1:3)*ds(1:3) ) )

        ! x, y, z coordinates for the surface;
        surf(1)%xx(ii,jj) = xx(1)
        surf(1)%yy(ii,jj) = xx(2)
        surf(1)%zz(ii,jj) = xx(3)

        ! dx/dt, dy/dt, dz/dt (dt for d theta)
        surf(1)%xt(ii,jj) = xt(1)
        surf(1)%yt(ii,jj) = xt(2)
        surf(1)%zt(ii,jj) = xt(3)

        ! dx/dp, dy/dp, dz/dp (dp for d zeta(phi))
        surf(1)%xp(ii,jj) = xz(1)
        surf(1)%yp(ii,jj) = xz(2)
        surf(1)%zp(ii,jj) = xz(3)

        ! surface normal vectors and ds for the jacobian;
        surf(1)%nx(ii,jj) = ds(1) / dd
        surf(1)%ny(ii,jj) = ds(2) / dd
        surf(1)%nz(ii,jj) = ds(3) / dd
        surf(1)%ds(ii,jj) =         dd

        ! using Gauss theorom; V = \int_S x \cdot n dt dz
        surf(1)%vol = surf(1)%vol + surf(1)%xx(ii,jj) * ds(1)

     enddo ! end of do jj; 14 Apr 16;
  enddo ! end of do ii; 14 Apr 16;

  surf(1)%vol = abs(surf(1)%vol) * discretefactor
  if( myid == 0 .and. IsQuiet <= 0) write(ounit, '(8X": Enclosed volume ="ES12.5" m^3 ;" )') surf(1)%vol

  !calculate target Bn with input harmonics; 05 Jan 17;
  if(NBnf >  0) then

     do jj = 0, Nzeta-1 ; zeta = ( jj + shift ) * pi2 / (Nzeta*Nfp)
        do ii = 0, Nteta-1 ; teta = ( ii + shift ) * pi2 / Nteta
           do imn = 1, NBnf
              arg = Bnim(imn) * teta - Bnin(imn) * zeta
              surf(1)%pb(ii,jj) = surf(1)%pb(ii,jj) + Bnc(imn)*cos(arg) + Bns(imn)*sin(arg)
           enddo
        enddo
     enddo

  endif
  
  
  return
  
end subroutine fousurf

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine surfcoord( theta, zeta, r, z)
  use globals, only: dp, zero, Nfou, bim, bin, Rbc, Rbs, Zbc, Zbs
  implicit none
  include "mpif.h"

  REAL, INTENT(in ) :: theta, zeta
  REAL, INTENT(out) :: r, z

  INTEGER           :: imn
  REAL              :: arg
  !-------------calculate r, z coodinates for theta, zeta------------------------------------------------  
  if( .not. allocated(bim) ) STOP  "please allocate surface data first!"

  r = zero; z = zero  
  do imn = 1, Nfou
     arg = bim(imn) * theta - bin(imn) * zeta
     R =  R + Rbc(imn) * cos(arg) + Rbs(imn) * sin(arg)
     Z =  Z + Zbc(imn) * cos(arg) + Zbs(imn) * sin(arg)
  enddo

  return
end subroutine surfcoord

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
