!title (boundary) ! Plasma boundary in Boozer coordinates

!latex \briefly{A Fourier representation for the plasma boundary in Boozer coordinates.}

!latex \calledby{\link{xfocus}}
!l tex \calls{\link{}}

!latex \tableofcontents

!latex \subsection{Overview}
!latex \bi
!latex \item[1.] The transformation between the Boozer and VMEC angles is calculated by the code 
!latex \href{https://bitbucket.org/lazerson_princeton/stellopt/wiki/BOOZ_XFORM}{\blu{BOOZ\_XFORM}}.
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

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine boozsurf
  
  use kmodule, only : zero, half, one, pi2, myid, ncpu, ounit, lunit, ext, Isymmetric, &
                      bmn, bNfp, nbf, bim, bin, bnim, bnin, Rbc, Rbs, Zbc, Zbs, Pmnc, Pmns, bnc, bns, &
                      Nteta, Nzeta, surf
  
  implicit none
  
  include "mpif.h"
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOGICAL :: exist
  INTEGER :: iosta, iostb, iostc, iostd, ioste, iostf, astat, ierr, ii, jj, imn
  REAL    :: RR(0:2), ZZ(0:2), Phi(0:2), szeta, czeta, &
             xx(1:3), xt(1:3), xz(1:3), ds(1:3), teta, zeta, arg, carg, sarg, dd, tmp
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  nbf = 0 !initialize
  
  inquire( file="plasma.boundary", exist=exist)
  
  FATAL( surface, .not.exist, plasma.boundary does not exist )
  
  if( myid.eq.0 ) open(lunit, file="plasma.boundary", status='old', action='read')
  
  if( myid.eq.0 ) read(lunit,*) !empty line
  if( myid.eq.0 ) read(lunit,*) bmn, bNfp, nbf !read dimensions
  
  !Broadcast the values
  IlBCAST( bmn , 1, 0 )
  IlBCAST( bNfp, 1, 0 )
  IlBCAST( nbf , 1, 0 )
  
  FATAL( surface, bmn .le.0, invalid )
  FATAL( surface, bNfp.le.0, invalid )
  
  SALLOCATE( bim, (1:bmn), 0 )
  SALLOCATE( bin, (1:bmn), 0 )
  
  SALLOCATE( Rbc, (1:bmn), zero )
  SALLOCATE( Rbs, (1:bmn), zero )
  SALLOCATE( Zbc, (1:bmn), zero )
  SALLOCATE( Zbs, (1:bmn), zero )
  SALLOCATE(Pmnc, (1:bmn), zero )
  SALLOCATE(Pmns, (1:bmn), zero )

  if( myid.eq.0 ) then
   read(lunit,*) !empty line
   read(lunit,*) !empty line
   do imn = 1, bmn
      read(lunit,*) bin(imn), bim(imn), Rbc(imn), Rbs(imn), Zbc(imn), Zbs(imn), Pmnc(imn), Pmns(imn)
   enddo
  endif  

  IlBCAST( bim(1:bmn), bmn, 0 )
  IlBCAST( bin(1:bmn), bmn, 0 )
  
  if (Isymmetric .eq. 0) bin(1:bmn) = bin(1:bmn) * bNfp  ! Disable periodicity
  
  RlBCAST( Rbc(1:bmn), bmn, 0 )
  RlBCAST( Rbs(1:bmn), bmn, 0 )
  RlBCAST( Zbc(1:bmn), bmn, 0 )
  RlBCAST( Zbs(1:bmn), bmn, 0 )
  RlBCAST(Pmnc(1:bmn), bmn, 0 )
  RlBCAST(Pmns(1:bmn), bmn, 0 )

  if( nbf .gt. 0) then  !read Bn terms
     SALLOCATE( bnim, (1:nbf), 0    )
     SALLOCATE( bnin, (1:nbf), 0    )
     SALLOCATE( bnc , (1:nbf), zero )
     SALLOCATE( bns , (1:nbf), zero )

     if( myid.eq.0 ) then
        read(lunit,*) !empty line
        read(lunit,*) !empty line
        do imn = 1, nbf ; read(lunit,*) bnin(imn), bnim(imn), bnc(imn), bns(imn)
        enddo
     endif

     IlBCAST( bnim(1:nbf), nbf, 0 )
     IlBCAST( bnin(1:nbf), nbf, 0 )

     if (Isymmetric .eq. 0) bnin(1:nbf) = bnin(1:nbf) * bNfp ! Disable periodicity
     
     RlBCAST( bnc(1:nbf) , nbf, 0 )
     RlBCAST( bns(1:nbf) , nbf, 0 )
  endif

  if( myid.eq.0 ) close(lunit,iostat=iosta)
  
  IlBCAST( iosta, 1, 0 )
  
  FATAL( surface, iosta.ne.0, error closing plasma.boundary )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( myid.eq.0 ) then
   write(ounit,'("surface : " 10x " : bmn ="i06 " ; bNfp = "i06 " ; nbf ="i06)') bmn, bNfp, nbf
#ifdef DEBUG
   write(ounit,'("surface : " 10x " : bim ="10i13   )') bim(1:bmn)
   write(ounit,'("surface : " 10x " : bin ="10i13   )') bin(1:bmn)
   write(ounit,'("surface : " 10x " : Rbc ="10es13.5)') Rbc(1:bmn)
   write(ounit,'("surface : " 10x " : Rbs ="10es13.5)') Rbs(1:bmn)
   write(ounit,'("surface : " 10x " : Zbc ="10es13.5)') Zbc(1:bmn)
   write(ounit,'("surface : " 10x " : Zbs ="10es13.5)') Zbs(1:bmn)
   write(ounit,'("surface : " 10x " : Pmnc ="10es13.5)') Pmnc(1:bmn)
   write(ounit,'("surface : " 10x " : Pmns ="10es13.5)') Pmns(1:bmn)
  if (nbf .gt. 0) then
   write(ounit,'("surface : " 10x " : bnim ="10i13  )') bnim(1:nbf)
   write(ounit,'("surface : " 10x " : bnin ="10i13  )') bnin(1:nbf)
   write(ounit,'("surface : " 10x " : bnc ="10es13.5)') bnc(1:nbf)
   write(ounit,'("surface : " 10x " : bns ="10es13.5)') bns(1:nbf)
   endif
#endif
  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  allocate( surf(1:1) ) ! can allow for multiple target boundaries if multiple currents are allowed; 14 Apr 16;
  
  surf(1)%Nteta = Nteta
  surf(1)%Nzeta = Nzeta
  
  SALLOCATE( surf(1)%xx, (0:Nteta,0:Nzeta), zero )
  SALLOCATE( surf(1)%yy, (0:Nteta,0:Nzeta), zero )
  SALLOCATE( surf(1)%zz, (0:Nteta,0:Nzeta), zero )
  SALLOCATE( surf(1)%nx, (0:Nteta,0:Nzeta), zero )
  SALLOCATE( surf(1)%ny, (0:Nteta,0:Nzeta), zero )
  SALLOCATE( surf(1)%nz, (0:Nteta,0:Nzeta), zero )
  SALLOCATE( surf(1)%ds, (0:Nteta,0:Nzeta), zero )
  SALLOCATE( surf(1)%xt, (0:Nteta,0:Nzeta), zero )
  SALLOCATE( surf(1)%yt, (0:Nteta,0:Nzeta), zero )
  SALLOCATE( surf(1)%zt, (0:Nteta,0:Nzeta), zero )
  SALLOCATE( surf(1)%bnt,(0:Nteta,0:Nzeta), zero )
 
! The center point value was used to discretize grid 
  do ii = 0, Nteta ; teta = ( ii + half ) * pi2 / Nteta
   do jj = 0, Nzeta ; zeta = ( jj + half ) * pi2 / Nzeta
    
    RR(0:2) = zero ; ZZ(0:2) = zero ; Phi(0:2) = zero
    
    do imn = 1, bmn 
       
     arg = bim(imn) * teta - bin(imn) * zeta 
     carg = cos(arg) ; sarg = sin(arg)
     
      RR(0) =  RR(0) +     Rbc(imn) * carg +  Rbs(imn) * sarg
      ZZ(0) =  ZZ(0) +     Zbc(imn) * carg +  Zbs(imn) * sarg
     Phi(0) = Phi(0) +    Pmnc(imn) * carg + Pmns(imn) * sarg
     
      RR(1) =  RR(1) + ( - Rbc(imn) * sarg +  Rbs(imn) * carg ) * bim(imn)
      ZZ(1) =  ZZ(1) + ( - Zbc(imn) * sarg +  Zbs(imn) * carg ) * bim(imn)
     Phi(1) = Phi(1) + ( -Pmnc(imn) * sarg + Pmns(imn) * carg ) * bim(imn)

      RR(2) =  RR(2) - ( - Rbc(imn) * sarg +  Rbs(imn) * carg ) * bin(imn)
      ZZ(2) =  ZZ(2) - ( - Zbc(imn) * sarg +  Zbs(imn) * carg ) * bin(imn)
     Phi(2) = Phi(2) - ( -Pmnc(imn) * sarg + Pmns(imn) * carg ) * bin(imn)
     
    enddo ! end of do imn; 12/12/2018;
  
    Phi(0) = Phi(0) + zeta
   !Phi(1) = Phi(1)
    Phi(2) = Phi(2) + one
  
    szeta = sin(Phi(0))
    czeta = cos(Phi(0))
    
    xx(1:3) = (/  RR(0) * czeta,  RR(0) * szeta, ZZ(0)     /)

    xt(1:3) = (/  RR(1) * czeta - RR(0) * szeta * Phi(1),  &
                  RR(1) * szeta + RR(0) * czeta * Phi(1),  & 
                  ZZ(1) /)
    
    xz(1:3) = (/  RR(2) * czeta - RR(0) * szeta * Phi(2),  &
                  RR(2) * szeta + RR(0) * czeta * Phi(2),  &
                  ZZ(2) /) 

    ds(1:3) = - (/ xt(2) * xz(3) - xt(3) * xz(2), &
                   xt(3) * xz(1) - xt(1) * xz(3), &
                   xt(1) * xz(2) - xt(2) * xz(1) /) !careful with the negative sign; means counterclockwise;

    dd = sqrt( sum( ds(1:3)*ds(1:3) ) )
   
    surf(1)%xx(ii,jj) = xx(1)
    surf(1)%yy(ii,jj) = xx(2)
    surf(1)%zz(ii,jj) = xx(3)

    surf(1)%xt(ii,jj) = xt(1)
    surf(1)%yt(ii,jj) = xt(2)
    surf(1)%zt(ii,jj) = xt(3)

    surf(1)%nx(ii,jj) = ds(1) / dd
    surf(1)%ny(ii,jj) = ds(2) / dd
    surf(1)%nz(ii,jj) = ds(3) / dd
    surf(1)%ds(ii,jj) =         dd

   enddo ! end of do jj; 14 Apr 16;
  enddo ! end of do ii; 14 Apr 16;

  !calculate target Bn with input harmonics; 05 Jan 17;
  if(nbf .gt. 0) then
     
     do ii = 0, Nteta ; teta = ( ii + half ) * pi2 / Nteta
        do jj = 0, Nzeta ; zeta = ( jj + half ) * pi2 / Nzeta

           tmp = 0.0
           do imn = 1, nbf

           !if( ii.eq.0 .and. jj.eq.0 ) write(ounit,'(2i6,2es23.15)') bnim(imn), bnin(imn), bnc(imn), bns(imn)
           !if( ii.eq.8 .and. jj.eq.8 ) write(ounit,'(2i6,2es23.15)') bnim(imn), bnin(imn), bnc(imn), bns(imn)

              arg = bnim(imn) * teta - bnin(imn) * zeta ; carg = cos(arg) ; sarg = sin(arg)
              tmp = tmp + bnc(imn) * carg + bns(imn) * sarg
              !if (ii.eq.0 .and. jj .eq. 0) write(*,*) tmp
              !surf(1)%bnt(ii, jj) = surf(1)%bnt(ii, jj) + bnc(imn) * carg + bns(imn) * sarg !there is bug for this; 05 Jan 17;
           enddo
           surf(1)%bnt(ii, jj) = tmp

          !if( ii.eq.0 .and. jj.eq.0 ) write(ounit,*) surf(1)%bnt(ii, jj)
          !if( ii.eq.8 .and. jj.eq.8 ) write(ounit,*) teta, zeta, surf(1)%bnt(ii, jj)

        enddo
     enddo

  endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  return
  
end subroutine boozsurf
