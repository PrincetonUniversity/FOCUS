
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

subroutine fousurf(filename, index)
  use globals, only : dp, zero, half, pi2, myid, ounit, runit, IsQuiet, IsSymmetric,  &
                      Nteta, Nzeta, surf, discretefactor, Nfp, plasma, symmetry,      &
                      tflux_sign, cosnfp, sinnfp
  use mpi
  implicit none

  CHARACTER*100, INTENT(IN) :: filename
  INTEGER, INTENT(IN) :: index
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  INTEGER :: iosta, astat, ierr, ii, jj, imn, Nfou, Nbnf, ip
  REAL    :: RR(0:2), ZZ(0:2), szeta, czeta, xx(1:3), xt(1:3), xz(1:3), ds(1:3), &
             teta, zeta, arg, dd, dz
  
  ! read the header
  if( myid == 0 ) then
     open(runit, file=trim(filename), status='old', action='read')
     read(runit,*) !empty line
     read(runit,*) surf(index)%Nfou, surf(index)%Nfp, surf(index)%NBnf !read dimensions
  endif
  
  !Broadcast the values
  IlBCAST( surf(index)%Nfou , 1, 0 )
  IlBCAST( surf(index)%Nfp  , 1, 0 )
  IlBCAST( surf(index)%NBnf , 1, 0 )  
  FATAL( rdsurf, surf(index)%Nfou <= 0, invalid )
  FATAL( rdsurf, surf(index)%Nfp  <= 0, invalid )
  FATAL( rdsurf, surf(index)%NBnf <  0, invalid )
  Nfou = surf(index)%Nfou
  NBnf = surf(index)%NBnf

  !allocate arrays
  SALLOCATE( surf(index)%bim, (1:Nfou), 0 )
  SALLOCATE( surf(index)%bin, (1:Nfou), 0 )  
  SALLOCATE( surf(index)%Rbc, (1:Nfou), zero )
  SALLOCATE( surf(index)%Rbs, (1:Nfou), zero )
  SALLOCATE( surf(index)%Zbc, (1:Nfou), zero )
  SALLOCATE( surf(index)%Zbs, (1:Nfou), zero )

  if( myid == 0 ) then
   read(runit,*) !empty line
   read(runit,*) !empty line
   do imn = 1, surf(index)%Nfou
      read(runit,*) surf(index)%bin(imn), surf(index)%bim(imn), surf(index)%Rbc(imn), &
           & surf(index)%Rbs(imn), surf(index)%Zbc(imn), surf(index)%Zbs(imn)
   enddo
  endif  

  IlBCAST( surf(index)%bim(1:Nfou), surf(index)%Nfou, 0 )
  IlBCAST( surf(index)%bin(1:Nfou), surf(index)%Nfou, 0 )
 
  surf(index)%bin(1:Nfou) = surf(index)%bin(1:Nfou) * surf(index)%Nfp  !The full plasma;
     
  RlBCAST( surf(index)%Rbc(1:Nfou), surf(index)%Nfou, 0 )
  RlBCAST( surf(index)%Rbs(1:Nfou), surf(index)%Nfou, 0 )
  RlBCAST( surf(index)%Zbc(1:Nfou), surf(index)%Nfou, 0 )
  RlBCAST( surf(index)%Zbs(1:Nfou), surf(index)%Nfou, 0 )

  !read Bnormal ditributions
  if( surf(index)%NBnf > 0) then
     SALLOCATE( surf(index)%Bnim, (1:NBnf), 0    )
     SALLOCATE( surf(index)%Bnin, (1:NBnf), 0    )
     SALLOCATE( surf(index)%Bnc , (1:NBnf), zero )
     SALLOCATE( surf(index)%Bns , (1:NBnf), zero )

     if( myid == 0 ) then
        read(runit,*) !empty line
        read(runit,*) !empty line
        do imn = 1, surf(index)%NBnf 
           read(runit,*) surf(index)%Bnin(imn), surf(index)%Bnim(imn), surf(index)%Bnc(imn), surf(index)%Bns(imn)
        enddo
     endif

     IlBCAST( surf(index)%Bnim(1:NBnf), surf(index)%NBnf, 0 )
     IlBCAST( surf(index)%Bnin(1:NBnf), surf(index)%NBnf, 0 )

     !if (IsSymmetric  ==  0)
     surf(index)%Bnin(1:NBnf) = surf(index)%Bnin(1:NBnf) * surf(index)%Nfp ! periodicity;
     ! This should be consistent with bnftran; Before fully constructed the stellarator symmetry,
     ! it's turned off;
     
     RlBCAST( surf(index)%Bnc(1:NBnf) , surf(index)%NBnf, 0 )
     RlBCAST( surf(index)%Bns(1:NBnf) , surf(index)%NBnf, 0 )
  endif

  if( myid == 0 ) close(runit,iostat=iosta)
  
  IlBCAST( iosta, 1, 0 )
  
  FATAL( surface, iosta.ne.0, error closing the surface )
  
  !-------------output for check-------------------------------------------------------------------------
  if( myid == 0 .and. IsQuiet <= 0) then
     write(ounit, *) "-----------Reading surface-----------------------------------"
     write(ounit, '("surface : The surface ", A," will be discretized in "I6" X "I6" elements.")') trim(filename), Nteta, Nzeta
     write(ounit, '(8X": Nfou = " I06 " ; Nfp = " I06 " ; NBnf = " I06 " ;" )') surf(index)%Nfou, surf(index)%Nfp, surf(index)%NBnf
  endif

  if( myid == 0 .and. IsQuiet <= -2) then ! very detailed output;
     write(ounit,'("        : " 10x " : bim ="10i13   )') surf(index)%bim(1:Nfou)
     write(ounit,'("        : " 10x " : bin ="10i13   )') surf(index)%bin(1:Nfou)
     write(ounit,'("        : " 10x " : Rbc ="10es13.5)') surf(index)%Rbc(1:Nfou)
     write(ounit,'("        : " 10x " : Rbs ="10es13.5)') surf(index)%Rbs(1:Nfou)
     write(ounit,'("        : " 10x " : Zbc ="10es13.5)') surf(index)%Zbc(1:Nfou)
     write(ounit,'("        : " 10x " : Zbs ="10es13.5)') surf(index)%Zbs(1:Nfou)
     if(Nbnf > 0) then
        write(ounit,'("        : " 10x " : Bnim ="10i13  )') surf(index)%Bnim(1:NBnf)
        write(ounit,'("        : " 10x " : Bnin ="10i13  )') surf(index)%Bnin(1:NBnf)
        write(ounit,'("        : " 10x " : Bnc ="10es13.5)') surf(index)%Bnc (1:NBnf)
        write(ounit,'("        : " 10x " : Bns ="10es13.5)') surf(index)%Bns (1:NBnf)
     endif
  endif

  surf(index)%Nteta = Nteta 
  surf(index)%Nzeta = Nzeta

  if (index == plasma) then    
     select case (IsSymmetric)
     case ( 0 )
        Nfp = 1                    ! reset Nfp to 1;
        symmetry = 0
     case ( 1 )                    ! plasma and coil periodicity enabled;
        Nfp = surf(plasma)%Nfp     ! use the raw Nfp
        symmetry = 0
     case ( 2 )                    ! stellarator symmetry enforced;
        Nfp = surf(plasma)%Nfp     ! use the raw Nfp
        symmetry = 1     
     end select

     SALLOCATE( cosnfp, (1:Nfp), zero )
     SALLOCATE( sinnfp, (1:Nfp), zero )  
     do ip = 1, Nfp
        cosnfp(ip) = cos((ip-1)*pi2/Nfp)
        sinnfp(ip) = sin((ip-1)*pi2/Nfp)
     enddo
     ! discretefactor = discretefactor/Nfp
     surf(index)%Nzeta = Nzeta * Nfp * 2**symmetry ! the total number from [0, 2pi]
     discretefactor = (pi2/surf(plasma)%Nteta) * (pi2/surf(plasma)%Nzeta)
  endif
  
  SALLOCATE( surf(index)%xx, (0:Nteta-1,0:Nzeta-1), zero ) !x coordinates;
  SALLOCATE( surf(index)%yy, (0:Nteta-1,0:Nzeta-1), zero ) !y coordinates
  SALLOCATE( surf(index)%zz, (0:Nteta-1,0:Nzeta-1), zero ) !z coordinates
  SALLOCATE( surf(index)%nx, (0:Nteta-1,0:Nzeta-1), zero ) !unit nx;
  SALLOCATE( surf(index)%ny, (0:Nteta-1,0:Nzeta-1), zero ) !unit ny;
  SALLOCATE( surf(index)%nz, (0:Nteta-1,0:Nzeta-1), zero ) !unit nz;
  SALLOCATE( surf(index)%ds, (0:Nteta-1,0:Nzeta-1), zero ) !jacobian;
  SALLOCATE( surf(index)%xt, (0:Nteta-1,0:Nzeta-1), zero ) !dx/dtheta;
  SALLOCATE( surf(index)%yt, (0:Nteta-1,0:Nzeta-1), zero ) !dy/dtheta;
  SALLOCATE( surf(index)%zt, (0:Nteta-1,0:Nzeta-1), zero ) !dz/dtheta;
  SALLOCATE( surf(index)%pb, (0:Nteta-1,0:Nzeta-1), zero ) !target Bn;
  SALLOCATE( surf(index)%xp, (0:Nteta-1,0:Nzeta-1), zero ) !dx/dzeta;
  SALLOCATE( surf(index)%yp, (0:Nteta-1,0:Nzeta-1), zero ) !dy/dzeta;
  SALLOCATE( surf(index)%zp, (0:Nteta-1,0:Nzeta-1), zero ) !dz/dzeta;
  
  surf(index)%vol = zero  ! volume enclosed by plasma boundary
  surf(index)%area = zero ! surface area
 
  ! The center point value was used to discretize grid;
  do ii = 0, Nteta-1
     teta = ( ii + half ) * pi2 / surf(index)%Nteta
     do jj = 0, Nzeta-1
        zeta = ( jj + half ) * pi2 / surf(index)%Nzeta
        RR(0:2) = zero ; ZZ(0:2) = zero
        do imn = 1, surf(index)%Nfou
           arg = surf(index)%bim(imn) * teta - surf(index)%bin(imn) * zeta
           RR(0) =  RR(0) +     surf(index)%Rbc(imn) * cos(arg) + surf(index)%Rbs(imn) * sin(arg)
           ZZ(0) =  ZZ(0) +     surf(index)%Zbc(imn) * cos(arg) + surf(index)%Zbs(imn) * sin(arg)
           RR(1) =  RR(1) + ( - surf(index)%Rbc(imn) * sin(arg) + surf(index)%Rbs(imn) * cos(arg) ) * surf(index)%bim(imn)
           ZZ(1) =  ZZ(1) + ( - surf(index)%Zbc(imn) * sin(arg) + surf(index)%Zbs(imn) * cos(arg) ) * surf(index)%bim(imn)
           RR(2) =  RR(2) - ( - surf(index)%Rbc(imn) * sin(arg) + surf(index)%Rbs(imn) * cos(arg) ) * surf(index)%bin(imn)
           ZZ(2) =  ZZ(2) - ( - surf(index)%Zbc(imn) * sin(arg) + surf(index)%Zbs(imn) * cos(arg) ) * surf(index)%bin(imn)
        enddo ! end of do imn; 30 Oct 15;
        szeta = sin(zeta)
        czeta = cos(zeta)
        xx(1:3) = (/   RR(0) * czeta,   RR(0) * szeta, ZZ(0) /)
        xt(1:3) = (/   RR(1) * czeta,   RR(1) * szeta, ZZ(1) /)
        xz(1:3) = (/   RR(2) * czeta,   RR(2) * szeta, ZZ(2) /) &
                + (/ - RR(0) * szeta,   RR(0) * czeta, zero  /)
        ! minus sign for theta counterclockwise direction;
        ds(1:3) = -(/ xt(2) * xz(3) - xt(3) * xz(2), & 
                      xt(3) * xz(1) - xt(1) * xz(3), &
                      xt(1) * xz(2) - xt(2) * xz(1) /)
        dd = sqrt( sum( ds(1:3)*ds(1:3) ) )
        ! x, y, z coordinates for the surface;
        surf(index)%xx(ii,jj) = xx(1)
        surf(index)%yy(ii,jj) = xx(2)
        surf(index)%zz(ii,jj) = xx(3)
        ! dx/dt, dy/dt, dz/dt (dt for d theta)
        surf(index)%xt(ii,jj) = xt(1)
        surf(index)%yt(ii,jj) = xt(2)
        surf(index)%zt(ii,jj) = xt(3)
        ! dx/dp, dy/dp, dz/dp (dp for d zeta(phi))
        surf(index)%xp(ii,jj) = xz(1)
        surf(index)%yp(ii,jj) = xz(2)
        surf(index)%zp(ii,jj) = xz(3)
        ! surface normal vectors and ds for the jacobian;
        surf(index)%nx(ii,jj) = ds(1) / dd
        surf(index)%ny(ii,jj) = ds(2) / dd
        surf(index)%nz(ii,jj) = ds(3) / dd
        surf(index)%ds(ii,jj) =         dd
        ! using Gauss theorom; V = \int_S x \cdot n dt dz
        surf(index)%vol = surf(index)%vol + surf(index)%xx(ii,jj) * ds(1) &
             & + surf(index)%yy(ii,jj) * ds(2) + surf(index)%zz(ii,jj) * ds(3)
        ! surface area 
        surf(index)%area = surf(index)%area + surf(index)%ds(ii,jj)
     enddo ! end of do jj; 14 Apr 16;
  enddo ! end of do ii; 14 Apr 16;

  ! print volume and area
  surf(index)%vol  = abs(surf(index)%vol)/3 * (pi2/surf(index)%Nteta) * (pi2/surf(index)%Nzeta)
  surf(index)%area = abs(surf(index)%area)  * (pi2/surf(index)%Nteta) * (pi2/surf(index)%Nzeta)
  if (index == plasma) then
     surf(index)%vol  = surf(index)%vol  * Nfp * 2**symmetry
     surf(index)%area = surf(index)%area * Nfp * 2**symmetry
  endif    
     
  if( myid == 0 .and. IsQuiet <= 0) then
     write(ounit, '(8X": Enclosed total surface volume ="ES12.5" m^3 ; area ="ES12.5" m^2." )') &
          surf(index)%vol, surf(index)%area
  endif

  ! check theta direction for the plasma surface and determine the toroidal flux sign
  if (index == plasma) then
     dz = surf(plasma)%zz(1,0) - surf(plasma)%zz(0,0)
     if (dz > 0) then
        ! counter-clockwise
        if( myid == 0) write(ounit, '(8X": The theta angle used is counter-clockwise.")')
        tflux_sign = -1
     else
        ! clockwise
        if( myid == 0) write(ounit, '(8X": The theta angle used is clockwise.")')
        tflux_sign =  1 
     endif
  endif

  !calculate target Bn with input harmonics; 05 Jan 17;
  if(surf(index)%NBnf >  0) then
     do jj = 0, Nzeta-1
        zeta = ( jj + half ) * pi2 / surf(index)%Nzeta
        do ii = 0, Nteta-1
           teta = ( ii + half ) * pi2 / surf(index)%Nteta
           do imn = 1, surf(index)%NBnf
              arg = surf(index)%Bnim(imn) * teta - surf(index)%Bnin(imn) * zeta
              surf(index)%pb(ii,jj) = surf(index)%pb(ii,jj) + surf(index)%Bnc(imn)*cos(arg) + surf(index)%Bns(imn)*sin(arg)
           enddo
        enddo
     enddo
  endif

  return
  
end subroutine fousurf

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine surfcoord( index, theta, zeta, r, z)
  use globals, only: dp, zero, surf
  use mpi
  implicit none

  INTEGER, INTENT(in) :: index
  REAL, INTENT(in ) :: theta, zeta
  REAL, INTENT(out) :: r, z

  INTEGER           :: imn
  REAL              :: arg
  !-------------calculate r, z coodinates for theta, zeta------------------------------------------------  
  if( .not. allocated(surf(index)%bim) ) STOP  "please allocate surface data first!"

  r = zero; z = zero  
  do imn = 1, surf(index)%Nfou
     arg = surf(index)%bim(imn) * theta - surf(index)%bin(imn) * zeta
     R =  R + surf(index)%Rbc(imn) * cos(arg) + surf(index)%Rbs(imn) * sin(arg)
     Z =  Z + surf(index)%Zbc(imn) * cos(arg) + surf(index)%Zbs(imn) * sin(arg)
  enddo

  return
end subroutine surfcoord

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
