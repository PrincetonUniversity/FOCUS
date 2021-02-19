!title (boundary) ! Plasma boundary in Boozer coordinates

!latex \briefly{A Fourier representation for the plasma boundary in Boozer coordinates.}

!latex \calledby{\link{xfocus}}
!l tex \calls{\link{}}

!latex \tableofcontents

!latex \subsection{Overview}
!latex \bi
!latex \item[1.] The transformation between the Boozer and VMEC angles is calculated by the code
!latex \href{https://princetonuniversity.github.io/STELLOPT/BOOZ_XFORM}{\blu{BOOZ\_XFORM}}.
!latex  By reading the output of BOOZ\_XFORM, we can easily get the information we need: $Rmnc\_B$,
!latex  $Zmns\_B$ and $Pmns\_B$. The three quantities are all decomposed in Boozer coordinates $(\t_B, \z_B)$.
!latex  (Please note: here and the following we assume stellarator symmetry, $lasym = F$.)
!latex
!latex  \item[2.] To inverse the transformation and map to cylindrical coordinates, we need to do
!latex  \be
!latex  R(\t_B, \z_B) & = & \sum_{m,n} Rmnc\_B_{mn} \cos(m\t_B - n\z_B) \nonumber \\
!latex  \phi & = & \z_B + p(\t_B, \z_B) \\
!latex  Z(\t_B, \z_B) & = & \sum_{m,n} Zmns\_B_{mn} \sin(m\t_B - n\z_B) \nonumber
!latex  \ee
!latex  Here, $\phi$ is the cylindrical angle and $p(\t_B, \z_B) = \sum_{mn} Pmns\_B \sin(m\t_B - n\z_B)$.
!latex
!latex  \item[3.] Now, we can deriv the transformation to Cartesian coordinates, as
!latex  \be
!latex  x(\t_B, \z_B) & = & R(\t_B, \z_B) * \cos(\phi) \nonumber \\
!latex  y(\t_B, \z_B) & = & R(\t_B, \z_B) * \sin(\phi)  \\
!latex  z(\t_B, \z_B) & = & Z(\t_B, \z_B) \nonumber
!latex  \ee
!latex
!latex  \ei
!latex
!latex  \subsection{More details}
!latex  For numerical applications in FOCUS, more quantities are computed.
!latex  \bi
!latex  \item[1.] Tangential derivatives (in the following, subscripts $B$ in $\t_B$, $\z_B$ and $Rmnc\_B$ etc. are neglected.)
!latex  \be
!latex  \pdv{x}{\t} & = & \pdv{R}{\t} \cos(\phi) - R \sin(\phi) \pdv{\phi}{\t} \nonumber \\
!latex              & = & \sum -m Rmnc \sin(m\t-n\z) \cos(\phi) - R \sin(\phi) [\sum \ m Pmns \cos(m\t-n\z)] \\
!latex  \pdv{x}{\z} & = & \pdv{R}{\z} \cos(\phi) - R \sin(\phi) \pdv{\phi}{\z} \nonumber \\
!latex              & = & \sum \, n Rmnc \sin(m\t-n\z) \cos(\phi) - R \sin(\phi) [1 + \sum -n Pmns \cos(m\t-n\z)]
!latex  \ee
!latex
!latex  \be
!latex  \pdv{y}{\t} & = & \pdv{R}{\t} \sin(\phi) + R \cos(\phi) \pdv{\phi}{\t} \nonumber \\
!latex              & = & \sum -m Rmnc \sin(m\t-n\z) \sin(\phi) + R \cos(\phi) [\sum \ m Pmns \cos(m\t-n\z)] \\
!latex  \pdv{y}{\z} & = & \pdv{R}{\z} \sin(\phi) + R \cos(\phi) \pdv{\phi}{\z} \nonumber \\
!latex              & = & \sum \, n Rmnc \sin(m\t-n\z) \sin(\phi) + R \cos(\phi) [1 + \sum -n Pmns \cos(m\t-n\z)]
!latex  \ee
!latex
!latex  \be
!latex  \pdv{z}{\t} & = & \sum \ m Zmns \cos(m\t-n\z) \\
!latex  \pdv{z}{\z} & = & \sum -n Zmns \cos(m\t-n\z)
!latex  \ee
!latex
!latex  \item[2.] The \emph{unit} normal vector of surface is defined as
!latex  \be
!latex  \vect{n} = \frac{\vect{x_{\z}} \times \vect{x_{\t}}}{|\vect{x_{\z}} \times \vect{x_{\t}}|} \ .
!latex  \ee
!latex  If poloidal and toroidal angle are both in counter-clockwise direction, the normal vector will point outwards.
!latex
!latex  \item[3.] The Jacobian is $J = |\vect{x_{\z}} \times \vect{x_{\t}}|$.
!latex
!latex  \ei

!latex  \subsection{Numerical implemention}
!latex  When \inputvar{Itopoloy = 2}, \FOCUS will read \emph{plasma.boundary}.
!latex  An example is showing as,
!latex  { \begin{verbatim}
!latex #bmn bNfp nbf
!latex 4 2 0
!latex #plasma boundary
!latex # n m   Rbc Rbs Zbc  Zbs  Pmnc Pmns
!latex   0 0  3.00 0.0 0.0  0.00 0.0  0.1
!latex   0 1  0.30 0.0 0.0 -0.30 0.0  0.1
!latex   1 0  0.00 0.0 0.0 -0.06 0.0  0.1
!latex   1 1 -0.06 0.0 0.0 -0.06 0.0  0.1
!latex #Bn harmonics
!latex # n m bnc bns
!latex 0 0 0.0 0.0
!latex \end{verbatim}
!latex }

!latex  $Rbc$, $Zbs$ and $Pmnc$ are all zero in stellarator symmetry.
!latex  Number of $B_n$ coefficients, $nbf$, should also be zero, except for RMP coils.
!latex

!latex \section{Boozer coordinates integration in FOCUS} \label{booz_xform}
!latex BOOZ\_XFORM transforms the VMEC coordinates into the straight field line coordinates introduced by Boozer.
!latex Although Boozer coordinates are normally used for representing magnetic fields, we are going to use Boozer angles ($\t_B, \z_B$) parameterizing the flux surface, on which magnetic field lines are lying.
!latex The surface in cylindrical coordinates is transformed into Boozer coordinates by
!latex \begin{equation}
!latex \left \{
!latex \begin{array}{lll}
!latex \ds R(\t_B, \z_B)   =  \sum_{m,n} R^B_{mn} \cos(m\t_B - n\z_B) ; \\
!latex \ds \phi(\t_B, \z_B)  =  \z_B +\sum_{m,n} P^B_{mn} \sin(m\t_B - n\z_B) ; \\
!latex \ds Z(\t_B, \z_B)  =  \sum_{m,n} Z^B_{mn} \sin(m\t_B - n\z_B) \ .
!latex \end{array}
!latex \right .
!latex \end{equation}
!latex Here, stellarator symmetry is assumed.
!latex Fourier coefficients, $R^B_{mn}$, $Z^B_{mn}$ and $P^B_{mn}$, are provided in BOOZ\_XFORM outputs.
!latex
!latex FOCUS uses Cartesian coordinates $\vect{x}(\t_B, \z_B)$, therefore, the flux surface is parameterized as
!latex \begin{equation}
!latex \left \{
!latex \begin{array}{lll}
!latex x(\t_B, \z_B)  =  R(\t_B, \z_B) \cos(\phi) ;  \\
!latex y(\t_B, \z_B)  =  R(\t_B, \z_B) \sin(\phi) ;  \\
!latex z(\t_B, \z_B)  =  Z(\t_B, \z_B) \ .
!latex \end{array}
!latex \right .
!latex \end{equation}
!latex Then the tangential derivatives are
!latex \begin{align}
!latex \pdv{x}{\t_B} = & \pdv{R}{\t_B} \cos(\phi) - R \sin(\phi) \pdv{\phi}{\t_B} ;\\
!latex \pdv{x}{\z_B} = & \pdv{R}{\z_B} \cos(\phi) - R \sin(\phi) \pdv{\phi}{\z_B} .
!latex \end{align}
!latex Likewise, $\partial y / \partial {\t_B}$, $\partial y / \partial {\z_B}$, $\partial z / \partial {\t_B}$ and $\partial z / \partial {\z_B}$ are computed.
!latex The surface normal is $\vect{n} = \vect{x}_{\z} \times \vect{x}_{\t}$ and so is the Jacobian $J = \frac{\partial (x, y, z)}{\partial (\t_B, \z_B)} = |\vect{x}_{\t} \times \vect{x}_{\z}|$.

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine rdbooz

    use globals, only : dp, zero, half, one, pi2, myid, ounit, runit, input_surf, IsQuiet, IsSymmetric, &
    Nfou, Nfp, NBnf, bim, bin, Bnim, Bnin, Rbc, Rbs, Zbc, Zbs, Bnc, Bns, Pmnc, Pmns, &
    Nteta, Nzeta, surf, discretefactor, Nfp_raw, cosnfp, sinnfp, &
    half_shift, shift, MPI_COMM_FAMUS, symm_factor
 
    use mpi
 
    implicit none

 
 !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
 
    INTEGER :: iosta, astat, ierr, ii, jj, imn, ip, index, symmetry
    REAL    :: RR(0:2), ZZ(0:2), Phi(0:2), szeta, czeta, &
               xx(1:3), xt(1:3), xz(1:3), ds(1:3), teta, zeta, arg, carg, sarg, dd, tmp, dz
 
 !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
 
    index = 1

    if (myid .eq. 0) then
       open (runit, file=trim(input_surf), status='old', action='read')
       read (runit, *) !empty line
       read (runit, *) Nfou, Nfp, NBnf !read dimensions
    endif
 
    !Broadcast the values
    IlBCAST(Nfou, 1, 0)
    IlBCAST(Nfp, 1, 0)
    IlBCAST(NBnf, 1, 0)
 
    FATAL(rdbooz, Nfou .le. 0, invalid)
    FATAL(rdbooz, Nfp .le. 0, invalid)
    FATAL(rdbooz, NBnf < 0, invalid)
 
    ! Allocate arrays
    SALLOCATE(bim, (1:Nfou), 0)
    SALLOCATE(bin, (1:Nfou), 0)
    SALLOCATE(Rbc, (1:Nfou), zero)
    SALLOCATE(Rbs, (1:Nfou), zero)
    SALLOCATE(Zbc, (1:Nfou), zero)
    SALLOCATE(Zbs, (1:Nfou), zero)
    SALLOCATE(Pmnc, (1:Nfou), zero)
    SALLOCATE(Pmns, (1:Nfou), zero)
 
    if (myid .eq. 0) then
       read (runit, *) !empty line
       read (runit, *) !empty line
       do imn = 1, Nfou
          read (runit, *) bin(imn), bim(imn), Rbc(imn), Rbs(imn), &
             Zbc(imn), Zbs(imn), Pmnc(imn), Pmns(imn)
       enddo
    endif
 
    IlBCAST(bim(1:Nfou), Nfou, 0)
    IlBCAST(bin(1:Nfou), Nfou, 0)
 
    bin(1:Nfou) = bin(1:Nfou)*Nfp ! Disable periodicity
 
    RlBCAST(Rbc(1:Nfou), Nfou, 0)
    RlBCAST(Rbs(1:Nfou), Nfou, 0)
    RlBCAST(Zbc(1:Nfou), Nfou, 0)
    RlBCAST(Zbs(1:Nfou), Nfou, 0)
    RlBCAST(Pmnc(1:Nfou), Nfou, 0)
    RlBCAST(Pmns(1:Nfou), Nfou, 0)
 
    if (NBnf .gt. 0) then  !read Bn terms
       SALLOCATE(Bnim, (1:NBnf), 0)
       SALLOCATE(Bnin, (1:NBnf), 0)
       SALLOCATE(Bnc, (1:NBnf), zero)
       SALLOCATE(Bns, (1:NBnf), zero)
 
       if (myid .eq. 0) then
          read (runit, *) !empty line
          read (runit, *) !empty line
          do imn = 1, NBnf
             read (runit, *) Bnin(imn), Bnim(imn), Bnc(imn), Bns(imn)
          enddo
       endif
 
       IlBCAST(Bnim(1:NBnf), NBnf, 0)
       IlBCAST(Bnin(1:NBnf), NBnf, 0)
 
       Bnin(1:NBnf) = Bnin(1:NBnf)*Nfp ! Disable periodicity
 
       RlBCAST(Bnc(1:NBnf), NBnf, 0)
       RlBCAST(Bns(1:NBnf), NBnf, 0)
    endif
 
    if (myid .eq. 0) close (runit, iostat=iosta)
 
    IlBCAST(iosta, 1, 0)
 
    FATAL(rdbooz, iosta .ne. 0, error closing surface)
 
 !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
    if (myid == 0 .and. IsQuiet <= 0) then
       write (ounit, *) "-----------Reading surface in Boozer coordinates--------------"
       write (ounit, '("surface : The surface ", A," will be discretized in "I6" X "I6" elements.")') trim(input_surf), Nteta, Nzeta
       write (ounit, '(8X": Nfou = " I06 " ; Nfp = " I06 " ; NBnf = " I06 " ;" )') Nfou, Nfp, NBnf
    endif
 
    if (myid == 0 .and. IsQuiet <= -2) then ! very detailed output;
       write (ounit, '("        : " 10x " : bim ="10i13   )') bim(1:Nfou)
       write (ounit, '("        : " 10x " : bin ="10i13   )') bin(1:Nfou)
       write (ounit, '("        : " 10x " : Rbc ="10es13.5)') Rbc(1:Nfou)
       write (ounit, '("        : " 10x " : Rbs ="10es13.5)') Rbs(1:Nfou)
       write (ounit, '("        : " 10x " : Zbc ="10es13.5)') Zbc(1:Nfou)
       write (ounit, '("        : " 10x " : Zbs ="10es13.5)') Zbs(1:Nfou)
       write (ounit, '("surface : " 10x " : Pmnc ="10es13.5)') Pmnc(1:Nfou)
       write (ounit, '("surface : " 10x " : Pmns ="10es13.5)') Pmns(1:Nfou)
 
       if (NBnf > 0) then
          write (ounit, '("        : " 10x " : Bnim ="10i13  )') Bnim(1:NBnf)
          write (ounit, '("        : " 10x " : Bnin ="10i13  )') Bnin(1:NBnf)
          write (ounit, '("        : " 10x " : Bnc ="10es13.5)') Bnc(1:NBnf)
          write (ounit, '("        : " 10x " : Bns ="10es13.5)') Bns(1:NBnf)
       endif
    endif
 
  !-------------discretize surface data------------------------------------------------------------------  
  
    Nfp_raw = Nfp ! save the raw value of Nfp
    select case (IsSymmetric)
    case ( 0 )
       Nfp = 1                          !reset Nfp to 1;
       symmetry = 0
    case ( 1 )                          !plasma and coil periodicity enabled;
       symmetry = 0
    case ( 2 )                          ! stellarator symmetry enforced;
       symmetry = 1     
    end select
    ! discretefactor = discretefactor/Nfp
    symm_factor = nfp * 2**symmetry
  
    SALLOCATE( cosnfp, (1:Nfp_raw), zero )
    SALLOCATE( sinnfp, (1:Nfp_raw), zero )  
    
    do ip = 1, Nfp_raw
       cosnfp(ip) = cos((ip-1)*pi2/Nfp_raw)
       sinnfp(ip) = sin((ip-1)*pi2/Nfp_raw)
    enddo
  
    allocate( surf(1:1) ) ! can allow for myltiple plasma boundaries 
                          ! if multiple currents are allowed; 14 Apr 16;
    
    surf(1)%Nteta = Nteta ! not used yet; used for multiple surfaces; 20170307;
    surf(1)%Nzeta = Nzeta * Nfp * 2**symmetry ! the total number from [0, 2pi]
 
    SALLOCATE(surf(index)%xx, (0:Nteta - 1, 0:Nzeta - 1), zero) !x coordinates;
    SALLOCATE(surf(index)%yy, (0:Nteta - 1, 0:Nzeta - 1), zero) !y coordinates
    SALLOCATE(surf(index)%zz, (0:Nteta - 1, 0:Nzeta - 1), zero) !z coordinates
    SALLOCATE(surf(index)%nx, (0:Nteta - 1, 0:Nzeta - 1), zero) !unit nx;
    SALLOCATE(surf(index)%ny, (0:Nteta - 1, 0:Nzeta - 1), zero) !unit ny;
    SALLOCATE(surf(index)%nz, (0:Nteta - 1, 0:Nzeta - 1), zero) !unit nz;
    SALLOCATE(surf(index)%ds, (0:Nteta - 1, 0:Nzeta - 1), zero) !jacobian;
    SALLOCATE(surf(index)%xt, (0:Nteta - 1, 0:Nzeta - 1), zero) !dx/dtheta;
    SALLOCATE(surf(index)%yt, (0:Nteta - 1, 0:Nzeta - 1), zero) !dy/dtheta;
    SALLOCATE(surf(index)%zt, (0:Nteta - 1, 0:Nzeta - 1), zero) !dz/dtheta;
    SALLOCATE(surf(index)%pb, (0:Nteta - 1, 0:Nzeta - 1), zero) !target Bn;
    SALLOCATE(surf(index)%xp, (0:Nteta - 1, 0:Nzeta - 1), zero) !dx/dzeta;
    SALLOCATE(surf(index)%yp, (0:Nteta - 1, 0:Nzeta - 1), zero) !dy/dzeta;
    SALLOCATE(surf(index)%zp, (0:Nteta - 1, 0:Nzeta - 1), zero) !dz/dzeta;
 
    surf(index)%vol = zero  ! volume enclosed by plasma boundary
    surf(index)%area = zero ! surface area
 
    discretefactor = (pi2/surf(1)%Nteta) * (pi2/surf(1)%Nzeta)

    if (half_shift) then
        shift = half
     else 
        shift = zero
        if(myid.eq.0) write(ounit, '(8X": half-shift in surface evaluation is turned off." )')
     endif
       
   ! The center point value was used to discretize grid;
     do ii = 0, Nteta-1
        teta = ( ii + shift ) * pi2 / surf(1)%Nteta
        do jj = 0, Nzeta-1
           zeta = ( jj + shift ) * pi2 / surf(1)%Nzeta
          RR(0:2) = zero; ZZ(0:2) = zero; Phi(0:2) = zero
          do imn = 1, Nfou
             arg = bim(imn)*teta - bin(imn)*zeta
             carg = cos(arg); sarg = sin(arg)
 
             RR(0) = RR(0) + Rbc(imn)*carg + Rbs(imn)*sarg
             ZZ(0) = ZZ(0) + Zbc(imn)*carg + Zbs(imn)*sarg
             Phi(0) = Phi(0) + Pmnc(imn)*carg + Pmns(imn)*sarg
 
             RR(1) = RR(1) + (-Rbc(imn)*sarg + Rbs(imn)*carg)*bim(imn)
             ZZ(1) = ZZ(1) + (-Zbc(imn)*sarg + Zbs(imn)*carg)*bim(imn)
             Phi(1) = Phi(1) + (-Pmnc(imn)*sarg + Pmns(imn)*carg)*bim(imn)
 
             RR(2) = RR(2) - (-Rbc(imn)*sarg + Rbs(imn)*carg)*bin(imn)
             ZZ(2) = ZZ(2) - (-Zbc(imn)*sarg + Zbs(imn)*carg)*bin(imn)
             Phi(2) = Phi(2) - (-Pmnc(imn)*sarg + Pmns(imn)*carg)*bin(imn)
          enddo ! end of do imn; 12/12/2018;
 
          Phi(0) = Phi(0) + zeta
          !Phi(1) = Phi(1)
          Phi(2) = Phi(2) + one
          szeta = sin(Phi(0))
          czeta = cos(Phi(0))
          xx(1:3) = (/RR(0)*czeta, RR(0)*szeta, ZZ(0)/)
          xt(1:3) = (/RR(1)*czeta - RR(0)*szeta*Phi(1), &
                      RR(1)*szeta + RR(0)*czeta*Phi(1), &
                      ZZ(1)/)
          xz(1:3) = (/RR(2)*czeta - RR(0)*szeta*Phi(2), &
                      RR(2)*szeta + RR(0)*czeta*Phi(2), &
                      ZZ(2)/)
          ds(1:3) = -(/xt(2)*xz(3) - xt(3)*xz(2), &
                       xt(3)*xz(1) - xt(1)*xz(3), &
                       xt(1)*xz(2) - xt(2)*xz(1)/) !careful with the negative sign; means counterclockwise;
          dd = sqrt(sum(ds(1:3)*ds(1:3)))
          ! x, y, z coordinates for the surface;
          surf(index)%xx(ii, jj) = xx(1)
          surf(index)%yy(ii, jj) = xx(2)
          surf(index)%zz(ii, jj) = xx(3)
          ! dx/dt, dy/dt, dz/dt (dt for d theta)
          surf(index)%xt(ii, jj) = xt(1)
          surf(index)%yt(ii, jj) = xt(2)
          surf(index)%zt(ii, jj) = xt(3)
          ! dx/dp, dy/dp, dz/dp (dp for d zeta(phi))
          surf(index)%xp(ii, jj) = xz(1)
          surf(index)%yp(ii, jj) = xz(2)
          surf(index)%zp(ii, jj) = xz(3)
          ! surface normal vectors and ds for the jacobian;
          surf(index)%nx(ii, jj) = ds(1)/dd
          surf(index)%ny(ii, jj) = ds(2)/dd
          surf(index)%nz(ii, jj) = ds(3)/dd
          surf(index)%ds(ii, jj) = dd
          ! using Gauss theorom; V = \int_S x \cdot n dt dz
          surf(index)%vol = surf(index)%vol + surf(index)%xx(ii, jj)*ds(1) &
                                          & + surf(index)%yy(ii, jj)*ds(2) &
                                          & + surf(index)%zz(ii, jj)*ds(3)
          ! surface area
          surf(index)%area = surf(index)%area + surf(index)%ds(ii, jj)
       enddo ! end of do jj; 14 Apr 16;
    enddo ! end of do ii; 14 Apr 16;
 
    ! print volume and area
    surf(index)%vol = abs(surf(index)%vol)/3 * discretefactor * Nfp * 2**symmetry
    surf(index)%area = abs(surf(index)%area) * discretefactor * Nfp * 2**symmetry
 
    if (myid == 0 .and. IsQuiet <= 0) then
       write (ounit, '(8X": Enclosed total surface volume ="ES12.5" m^3 ; area ="ES12.5" m^2." )') &
          surf(index)%vol, surf(index)%area
    endif
 
    ! ! check theta direction for the plasma surface and determine the toroidal flux sign
    ! if (index == plasma) then
    !    dz = surf(plasma)%zz(1, 0) - surf(plasma)%zz(0, 0)
    !    if (dz > 0) then
    !       ! counter-clockwise
    !       if (myid == 0) write (ounit, '(8X": The theta angle used is counter-clockwise.")')
    !       tflux_sign = -1
    !    else
    !       ! clockwise
    !       if (myid == 0) write (ounit, '(8X": The theta angle used is clockwise.")')
    !       tflux_sign = 1
    !    endif
    ! endif
 
  !calculate target Bn with input harmonics; 05 Jan 17;
    if(NBnf >  0) then
        do jj = 0, Nzeta-1
           zeta = ( jj + shift ) * pi2 / surf(1)%Nzeta
           do ii = 0, Nteta-1
              teta = ( ii + shift ) * pi2 / surf(1)%Nteta
              do imn = 1, NBnf
                 arg = Bnim(imn) * teta - Bnin(imn) * zeta
                 surf(1)%pb(ii,jj) = surf(1)%pb(ii,jj) + Bnc(imn)*cos(arg) + Bns(imn)*sin(arg)
              enddo
           enddo
        enddo
     endif
 
 !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
 
    return
 
 end subroutine rdbooz
 