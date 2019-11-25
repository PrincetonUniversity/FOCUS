
!title (energy) ! Construct &ldquo;energy&rdquo; (cost function). 

!latex \briefly{Construct ``energy'' function, and constraints.}

!latex \calledby{\link{descent}}
!latex \calls{\link{bnormal}, \link{torflux}, \link{tlength}, \link{specwid}, \link{ccsep}}

!latex \tableofcontents

!latex \newcommand{\cmt}{\cos(mt)}
!latex \newcommand{\smt}{\sin(mt)}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
   
!latex \subsection{energy functional}
!latex The total energy can be represented as,
!latex \be E & = & \sum_i \omega_i \, ( \chi_i - O_i )^2 \nonumber \\
!latex        & = & \omega_{bnorm} \int_s \frac{1}{2} \frac{(\vec{B} \cdot \vec{n})^2}{|B|^2} ds \, + \, \omega_{tflux} \sum_{i=1}^{Nzeta} \frac{1}{2}
!latex              \frac{( \Psi_i-\Psi_o )^2}{{\Psi_o}^2} \, + \, \omega_{length} \frac{\sum_{i=1}^{Ncoils} \exp(L_w L_i)}{N_c \exp(L_w L_o)} \\
!latex %  \omega_{length} \sum_{i=1}^{Ncoils} \frac{1}{2} \frac{( L_i-L_o )^2}{{L_o}^2} \, \nonumber \\
!latex       & & \, + \, \omega_{eqarc} \sum_{i=1}^{Ncoils} \sum_{n=0}^{NFcoil} ({\lambda^i_n}^2) \, + \,
!latex            \omega_{ccsep} \sum_{i,j = 1}^{Ncoils}\int_{c_i}\int_{C_j} \frac{dl_i \ dl_j}{\Delta r^2} \, \ldots \nonumber \\
!latex \ee
!latex Right now, the constraints of bnormal, toroidal flux, length, equal-arc angle and coil-coil separation are constructed. Later, more constraints, like
!latex coil-plasma separation, coil curvature etc.\ , will be added.\\
!latex \emph{ The subroutines denergy, dBnxyz, dlength were originally written by Dr. S. Hudson, which has a relly fast speed calculating the energy functional 
!latex        and the first derivatives. It's temporarily turned off.}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine denergy( tau, xdof, dE )
  
  use kmodule, only : zero, one, half, two, ten, pi2, sqrtmachprec, ncpu, myid, ounit, &
                      NFcoil, NDcoil, &
                     !xkc, xks, ykc, yks, zkc, zks, &
                     !alpha, kpotent, kspring, kcenter, Npoints, &
                      coil, icoil, Ncoils, cmt, smt, &
                      Nteta, Nzeta, absreq, relreq, itime, iteta, jzeta, surf, &
                      totalenergy, discretefactor, ttlen
  
  implicit none
  
  include "mpif.h"
  
  REAL                 :: tau, xdof(*), dE(*)

  INTEGER, parameter   :: Ndimensions = 1
  
  INTEGER              :: ierr, astat, Ndof, idof, NN, ii, mm, Ndim, Nfun, id01eaf, lmincalls, lmaxcalls, RR, Lwk, kk, llmodnp, ifail, itt
  REAL                 :: tt(1:Ndimensions), xx(0:1), yy(0:1), zz(0:1)
! REAL                 :: angle, teta, zeta, at(1:3), az(1:3), xx(1:3), xt(1:3), xz(1:3), dx(1:3), dxo, ddx, ca, sa, dca, dsa, cxx, cyy, czz 

  REAL                 :: AA(1:Ndimensions), BB(1:Ndimensions), labsreq, lrelreq, localenergy, Ifactor, Lfactor
  
  REAL   , allocatable :: FF(:,:), DD(:), dB(:,:,:,:,:), Bn(:), wk(:), lE(:), dL(:,:,:,:), Fi(:)
  
  external             :: dBndxyz, dlength
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  NN = NFcoil ! shorthand; 14 Apr 16;
  
  Ndof = 0
  do icoil = 1, Ncoils
   ;             ; Ndof = Ndof + 1 ; coil(icoil)%I      = xdof(Ndof)
   ;  mm = 0     ; Ndof = Ndof + 1 ; coil(icoil)%xc(mm) = xdof(Ndof)
   ;             ; Ndof = Ndof + 1 ; coil(icoil)%yc(mm) = xdof(Ndof)
   ;             ; Ndof = Ndof + 1 ; coil(icoil)%zc(mm) = xdof(Ndof)
   do mm = 1, NN ; Ndof = Ndof + 1 ; coil(icoil)%xc(mm) = xdof(Ndof)
    ;            ; Ndof = Ndof + 1 ; coil(icoil)%yc(mm) = xdof(Ndof)
    ;            ; Ndof = Ndof + 1 ; coil(icoil)%zc(mm) = xdof(Ndof)
    ;            ; Ndof = Ndof + 1 ; coil(icoil)%xs(mm) = xdof(Ndof)
    ;            ; Ndof = Ndof + 1 ; coil(icoil)%ys(mm) = xdof(Ndof)
    ;            ; Ndof = Ndof + 1 ; coil(icoil)%zs(mm) = xdof(Ndof)
   enddo

  enddo

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  do icoil = 1, Ncoils
   
   if( myid.ne.modulo(icoil-1,ncpu) ) cycle ! parallelization loop; 14 Apr 16;
   
   ;  mm = 0
   ;coil(icoil)%xx(0:NDcoil) =                            cmt(0:NDcoil,mm) * coil(icoil)%xc(mm)
   ;coil(icoil)%yy(0:NDcoil) =                            cmt(0:NDcoil,mm) * coil(icoil)%yc(mm)
   ;coil(icoil)%zz(0:NDcoil) =                            cmt(0:NDcoil,mm) * coil(icoil)%zc(mm)
   ;coil(icoil)%xt(0:NDcoil) =                                               zero
   ;coil(icoil)%yt(0:NDcoil) =                                               zero
   ;coil(icoil)%zt(0:NDcoil) =                                               zero
   do mm = 1, NFcoil
    coil(icoil)%xx(0:NDcoil) = coil(icoil)%xx(0:NDcoil) + (   cmt(0:NDcoil,mm) * coil(icoil)%xc(mm) + smt(0:NDcoil,mm) * coil(icoil)%xs(mm) )
    coil(icoil)%yy(0:NDcoil) = coil(icoil)%yy(0:NDcoil) + (   cmt(0:NDcoil,mm) * coil(icoil)%yc(mm) + smt(0:NDcoil,mm) * coil(icoil)%ys(mm) )
    coil(icoil)%zz(0:NDcoil) = coil(icoil)%zz(0:NDcoil) + (   cmt(0:NDcoil,mm) * coil(icoil)%zc(mm) + smt(0:NDcoil,mm) * coil(icoil)%zs(mm) )
    coil(icoil)%xt(0:NDcoil) = coil(icoil)%xt(0:NDcoil) + ( - smt(0:NDcoil,mm) * coil(icoil)%xc(mm) + cmt(0:NDcoil,mm) * coil(icoil)%xs(mm) ) * mm
    coil(icoil)%yt(0:NDcoil) = coil(icoil)%yt(0:NDcoil) + ( - smt(0:NDcoil,mm) * coil(icoil)%yc(mm) + cmt(0:NDcoil,mm) * coil(icoil)%ys(mm) ) * mm
    coil(icoil)%zt(0:NDcoil) = coil(icoil)%zt(0:NDcoil) + ( - smt(0:NDcoil,mm) * coil(icoil)%zc(mm) + cmt(0:NDcoil,mm) * coil(icoil)%zs(mm) ) * mm
   enddo ! end of do ii; 14 Apr 16;
   
  enddo ! end of do icoil; 14 Apr 16;
  
  do icoil = 1, Ncoils ; llmodnp = modulo(icoil-1,ncpu)
   RlBCAST( coil(icoil)%xx(0:NDcoil), NDcoil+1, llmodnp )
   RlBCAST( coil(icoil)%yy(0:NDcoil), NDcoil+1, llmodnp )
   RlBCAST( coil(icoil)%zz(0:NDcoil), NDcoil+1, llmodnp )
   RlBCAST( coil(icoil)%xt(0:NDcoil), NDcoil+1, llmodnp )
   RlBCAST( coil(icoil)%yt(0:NDcoil), NDcoil+1, llmodnp )
   RlBCAST( coil(icoil)%zt(0:NDcoil), NDcoil+1, llmodnp )
  enddo
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  Ndim = Ndimensions 
  
  Nfun = 3 * ( 1 + 3 * ( 1 + NN + NN ) )

  SALLOCATE( FF, (1:Nfun,1:Ncoils), zero ) ! in this section, FF is not required to have additional subscript indicating coil, but this is required below;
  SALLOCATE( Fi, (1:Nfun         ), zero )

#ifdef CONTINUOUS
  AA(1:Ndim) = zero ; BB(1:Ndim) = pi2 ; labsreq = absreq ; lrelreq = relreq
  lmincalls = 2**4 ; RR = 2**Ndim + 2*Ndim**2 + 2*Ndim + 1 ; lmaxcalls = 8 * max( lmincalls, RR ) ! what is suitable value for lmincalls; 14 Apr 16;
  Lwk = 6 * Ndim + 9 * Nfun + ( Ndim + Nfun + 2 ) * ( 1 + lmaxcalls / RR )
  SALLOCATE( wk, (1:Lwk          ), zero )
  SALLOCATE( DD, (1:Nfun         ), zero )
#endif

  SALLOCATE( dB, (1:Ncoils,1:3,0:3,0:1,0:NN), zero )
  
  SALLOCATE( Bn, (0:Ncoils), zero ) ! extra workspace: Bn(0) = sum( Bn(1:Ncoils) ) ; 14 Apr 16;

  SALLOCATE( lE, (1:Ndof), zero )

  localenergy = zero


  do iteta = 0, Nteta-1 ! iteta is global;
   do jzeta = 0, Nzeta-1 ! jzeta is global;
    
    if( myid.ne.modulo(iteta*Nzeta+jzeta,ncpu) ) cycle ! parallelization loop; 14 Apr 16;
    
    Bn(1:Ncoils) = zero ! normal magnetic field due to each coil; 14 Apr 16;
    
    do icoil = 1, Ncoils ! icoil is global;

#ifdef CONTINUOUS
     
      lmincalls = 2**4 ; RR = 2**Ndim + 2*Ndim**2 + 2*Ndim + 1 ; lmaxcalls = 8 * max( lmincalls, RR ) ! what is suitable value for lmincalls; 14 Apr 16;
      
      id01eaf = 1
      call D01EAF( Ndim, AA(1:Ndim), BB(1:Ndim), lmincalls, lmaxcalls, Nfun, dBndxyz, labsreq, lrelreq, Lwk, wk(1:Lwk), FF(1:Nfun,icoil), DD(1:Nfun), id01eaf )
      
      select case( id01eaf )
      case( 0 ) ! write(ounit,'("denergy : " 10x " : myid="i3" ; icoil="i4" ; id01daf="i2" ; success          ;")') myid, icoil, id01eaf
      case( 1 ) ; write(ounit,'("denergy : " 10x " : myid="i3" ; icoil="i4" ; id01daf="i2" ; maxcls too small ;")') myid, icoil, id01eaf
      case( 2 ) ; write(ounit,'("denergy : " 10x " : myid="i3" ; icoil="i4" ; id01daf="i2" ; Lwk    too small ;")') myid, icoil, id01eaf
      case( 3 ) ; write(ounit,'("denergy : " 10x " : myid="i3" ; icoil="i4" ; id01daf="i2" ; maxcls too small ;")') myid, icoil, id01eaf
      case( 4 ) ; write(ounit,'("denergy : " 10x " : myid="i3" ; icoil="i4" ; id01daf="i2" ; input error ;     ")') myid, icoil, id01eaf
      end select
 
#else

     FF(1:Nfun,icoil) = zero
     
     do itime = 0, NDcoil-1 ! itime is global; 14 Apr 16;

      tt(1) = itime * pi2 / NDcoil ! this can be made redundant; 14 Apr 16;

      call dBndxyz( Ndim, tt(1:Ndim), Nfun, Fi(1:Nfun) ) ! can change argument list; 14 Apr 16;

      FF(1:Nfun,icoil) = FF(1:Nfun,icoil) + Fi(1:Nfun)

     enddo

     FF(1:Nfun,icoil) = FF(1:Nfun,icoil) * (pi2/NDcoil) ! can make this factor global; 14 Apr 16;
     
#endif

     ;kk =      0
     ;kk = kk + 1 ; dB(icoil,1,0,0, 0) = FF(kk,icoil) ! this can be removed if dBndxyz arguments are changed; 14 Apr 16;
     do mm = 0, NN 
      kk = kk + 1 ; dB(icoil,1,1,0,mm) = FF(kk,icoil)
      kk = kk + 1 ; dB(icoil,1,2,0,mm) = FF(kk,icoil)
      kk = kk + 1 ; dB(icoil,1,3,0,mm) = FF(kk,icoil)
     enddo
     do mm = 1, NN 
      kk = kk + 1 ; dB(icoil,1,1,1,mm) = FF(kk,icoil)
      kk = kk + 1 ; dB(icoil,1,2,1,mm) = FF(kk,icoil)
      kk = kk + 1 ; dB(icoil,1,3,1,mm) = FF(kk,icoil)
     enddo
     ;kk = kk + 1 ; dB(icoil,2,0,0, 0) = FF(kk,icoil)
     do mm = 0, NN
      kk = kk + 1 ; dB(icoil,2,1,0,mm) = FF(kk,icoil)
      kk = kk + 1 ; dB(icoil,2,2,0,mm) = FF(kk,icoil)
      kk = kk + 1 ; dB(icoil,2,3,0,mm) = FF(kk,icoil)
     enddo
     do mm = 1, NN
      kk = kk + 1 ; dB(icoil,2,1,1,mm) = FF(kk,icoil)
      kk = kk + 1 ; dB(icoil,2,2,1,mm) = FF(kk,icoil)
      kk = kk + 1 ; dB(icoil,2,3,1,mm) = FF(kk,icoil)
     enddo
     ;kk = kk + 1 ; dB(icoil,3,0,0, 0) = FF(kk,icoil)
     do mm = 0, NN
      kk = kk + 1 ; dB(icoil,3,1,0,mm) = FF(kk,icoil)
      kk = kk + 1 ; dB(icoil,3,2,0,mm) = FF(kk,icoil)
      kk = kk + 1 ; dB(icoil,3,3,0,mm) = FF(kk,icoil)
     enddo
     do mm = 1, NN
      kk = kk + 1 ; dB(icoil,3,1,1,mm) = FF(kk,icoil)
      kk = kk + 1 ; dB(icoil,3,2,1,mm) = FF(kk,icoil)
      kk = kk + 1 ; dB(icoil,3,3,1,mm) = FF(kk,icoil)
     enddo
     FATAL( denergy, kk.ne.Nfun, counting error )

     Bn(icoil) = coil(icoil)%I * ( surf(1)%nx(iteta,jzeta) * dB(icoil,1,0,0,0) &
                                 + surf(1)%ny(iteta,jzeta) * dB(icoil,2,0,0,0) &
                                 + surf(1)%nz(iteta,jzeta) * dB(icoil,3,0,0,0) )
     
    enddo ! end of do icoil; 14 Apr 16;
    
    Bn(0) = sum( Bn(1:Ncoils) )

    localenergy = localenergy + Bn(0)**2 * surf(1)%ds(iteta,jzeta) ! factor of half included below; 14 Apr 16;

    Bn(0) = Bn(0) * surf(1)%ds(iteta,jzeta) ! re-defined Bn(0) to include Jacobian factor;

    idof = 0
    do icoil = 1, Ncoils
     ;             ; idof = idof + 1 ; lE(idof) = lE(idof) + Bn(0) * Bn(icoil)                                                        ! coil(icoil)%I
     ;  mm = 0     ; idof = idof + 1 ; lE(idof) = lE(idof) + Bn(0) * coil(icoil)%I * ( surf(1)%nx(iteta,jzeta) * dB(icoil,1,1,0,mm) &
                                                                                     + surf(1)%ny(iteta,jzeta) * dB(icoil,2,1,0,mm) &
                                                                                     + surf(1)%nz(iteta,jzeta) * dB(icoil,3,1,0,mm) ) ! coil(icoil)%xc(mm)
     ;             ; idof = idof + 1 ; lE(idof) = lE(idof) + Bn(0) * coil(icoil)%I * ( surf(1)%nx(iteta,jzeta) * dB(icoil,1,2,0,mm) &
                                                                                     + surf(1)%ny(iteta,jzeta) * dB(icoil,2,2,0,mm) &
                                                                                     + surf(1)%nz(iteta,jzeta) * dB(icoil,3,2,0,mm) ) ! coil(icoil)%yc(mm)
     ;             ; idof = idof + 1 ; lE(idof) = lE(idof) + Bn(0) * coil(icoil)%I * ( surf(1)%nx(iteta,jzeta) * dB(icoil,1,3,0,mm) &
                                                                                     + surf(1)%ny(iteta,jzeta) * dB(icoil,2,3,0,mm) &
                                                                                     + surf(1)%nz(iteta,jzeta) * dB(icoil,3,3,0,mm) ) ! coil(icoil)%zc(mm)
     do mm = 1, NN ; idof = idof + 1 ; lE(idof) = lE(idof) + Bn(0) * coil(icoil)%I * ( surf(1)%nx(iteta,jzeta) * dB(icoil,1,1,0,mm) &
                                                                                     + surf(1)%ny(iteta,jzeta) * dB(icoil,2,1,0,mm) &
                                                                                     + surf(1)%nz(iteta,jzeta) * dB(icoil,3,1,0,mm) ) ! coil(icoil)%xc(mm)
      ;            ; idof = idof + 1 ; lE(idof) = lE(idof) + Bn(0) * coil(icoil)%I * ( surf(1)%nx(iteta,jzeta) * dB(icoil,1,2,0,mm) &
                                                                                     + surf(1)%ny(iteta,jzeta) * dB(icoil,2,2,0,mm) &
                                                                                     + surf(1)%nz(iteta,jzeta) * dB(icoil,3,2,0,mm) ) ! coil(icoil)%yc(mm)
      ;            ; idof = idof + 1 ; lE(idof) = lE(idof) + Bn(0) * coil(icoil)%I * ( surf(1)%nx(iteta,jzeta) * dB(icoil,1,3,0,mm) &
                                                                                     + surf(1)%ny(iteta,jzeta) * dB(icoil,2,3,0,mm) &
                                                                                     + surf(1)%nz(iteta,jzeta) * dB(icoil,3,3,0,mm) ) ! coil(icoil)%zc(mm)
      ;            ; idof = idof + 1 ; lE(idof) = lE(idof) + Bn(0) * coil(icoil)%I * ( surf(1)%nx(iteta,jzeta) * dB(icoil,1,1,1,mm) &
                                                                                     + surf(1)%ny(iteta,jzeta) * dB(icoil,2,1,1,mm) &
                                                                                     + surf(1)%nz(iteta,jzeta) * dB(icoil,3,1,1,mm) ) ! coil(icoil)%xs(mm)
      ;            ; idof = idof + 1 ; lE(idof) = lE(idof) + Bn(0) * coil(icoil)%I * ( surf(1)%nx(iteta,jzeta) * dB(icoil,1,2,1,mm) &
                                                                                     + surf(1)%ny(iteta,jzeta) * dB(icoil,2,2,1,mm) &
                                                                                     + surf(1)%nz(iteta,jzeta) * dB(icoil,3,2,1,mm) ) ! coil(icoil)%ys(mm)
      ;            ; idof = idof + 1 ; lE(idof) = lE(idof) + Bn(0) * coil(icoil)%I * ( surf(1)%nx(iteta,jzeta) * dB(icoil,1,3,1,mm) &
                                                                                     + surf(1)%ny(iteta,jzeta) * dB(icoil,2,3,1,mm) &
                                                                                     + surf(1)%nz(iteta,jzeta) * dB(icoil,3,3,1,mm) ) ! coil(icoil)%zs(mm)
     enddo ! end of do mm; 14 Apr 16;
    enddo ! end of do icoil; 14 Apr 16;
    FATAL( denergy, idof.ne.Ndof, counting error )

   enddo ! end of do jzeta; 14 Apr 16;
  enddo ! end of do iteta; 14 Apr 16;

  DALLOCATE( Bn )
  DALLOCATE( dB )

  DALLOCATE( FF )
  DALLOCATE( Fi )

#ifdef CONTINUOUS
  DALLOCATE( wk )
  DALLOCATE( DD )
#endif
  call MPI_REDUCE( lE(1:Ndof) , dE(1:Ndof) , Ndof, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
  call MPI_REDUCE( localenergy, totalenergy,    1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
  DALLOCATE( lE )
  RlBCAST( dE(1:Ndof) , Ndof, 0 )
  RlBCAST( totalenergy,    1, 0 )
  totalenergy = half * totalenergy * discretefactor
  dE(1:Ndof)  =        dE(1:Ndof)  * discretefactor ! factor of - one included below; 14 Apr 16;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! enforce current contraints; 14 Apr 16;

  idof = 0
  do icoil = 1, Ncoils ; Ifactor = coil(icoil)%Ic * coil(icoil)%Iw * ( coil(icoil)%I - coil(icoil)%Io ) * coil(icoil)%I
   ;             ; idof = idof + 1 ; dE(idof) = dE(idof) + Ifactor                                                               ! coil(icoil)%I
   ;  mm = 0     ; idof = idof + 1 ! dE(idof) = dE(idof)                                                                         ! coil(icoil)%xc(mm)
   ;             ; idof = idof + 1 ! dE(idof) = dE(idof)                                                                         ! coil(icoil)%yc(mm)
   ;             ; idof = idof + 1 ! dE(idof) = dE(idof)                                                                         ! coil(icoil)%zc(mm)
   do mm = 1, NN ; idof = idof + 1 ! dE(idof) = dE(idof)                                                                         ! coil(icoil)%xc(mm)
    ;            ; idof = idof + 1 ! dE(idof) = dE(idof)                                                                         ! coil(icoil)%yc(mm)
    ;            ; idof = idof + 1 ! dE(idof) = dE(idof)                                                                         ! coil(icoil)%zc(mm)
    ;            ; idof = idof + 1 ! dE(idof) = dE(idof)                                                                         ! coil(icoil)%xs(mm)
    ;            ; idof = idof + 1 ! dE(idof) = dE(idof)                                                                         ! coil(icoil)%ys(mm)
    ;            ; idof = idof + 1 ! dE(idof) = dE(idof)                                                                         ! coil(icoil)%zs(mm)
   enddo
  enddo

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
! enforce length constraints 14 Apr 16;

  Ndim = Ndimensions ; AA(1:Ndim) = zero ; BB(1:Ndim) = pi2 ; labsreq = absreq ; lrelreq = relreq

  lmincalls = 2**4 ; RR = 2**Ndim + 2*Ndim**2 + 2*Ndim + 1 ; lmaxcalls = 8 * max( lmincalls, RR ) ! what is suitable value for lmincalls; 14 Apr 16;
  
  Nfun = ( 1 + 3 * ( 1 + NN + NN ) )

  Lwk = 6 * Ndim + 9 * Nfun + ( Ndim + Nfun + 2 ) * ( 1 + lmaxcalls / RR )
  
  SALLOCATE( wk, (1:Lwk          ), zero )
  SALLOCATE( FF, (1:Nfun,1:Ncoils), zero )
  SALLOCATE( DD, (1:Nfun         ), zero )

  do icoil = 1, Ncoils ! global;
   
   if( myid.ne.modulo(icoil-1,ncpu) ) cycle ! parallelization loop; 14 Apr 16;
 
   lmincalls = 2**4 ; RR = 2**Ndim + 2*Ndim**2 + 2*Ndim + 1 ; lmaxcalls = 8 * max( lmincalls, RR ) ! what is suitable value for lmincalls; 14 Apr 16;

   id01eaf = 1
   call D01EAF( Ndim, AA(1:Ndim), BB(1:Ndim), lmincalls, lmaxcalls, Nfun, dlength, labsreq, lrelreq, Lwk, wk(1:Lwk), FF(1:Nfun,icoil), DD(1:Nfun), id01eaf )
 
   select case( id01eaf )
   case( 0 ) ! write(ounit,'("denergy : " 10x " : myid="i3" ; icoil="i4" ; id01daf="i2" ; success          ;")') myid, icoil, id01eaf
   case( 1 ) ; write(ounit,'("denergy : " 10x " : myid="i3" ; icoil="i4" ; id01daf="i2" ; maxcls too small ;")') myid, icoil, id01eaf
   case( 2 ) ; write(ounit,'("denergy : " 10x " : myid="i3" ; icoil="i4" ; id01daf="i2" ; Lwk    too small ;")') myid, icoil, id01eaf
   case( 3 ) ; write(ounit,'("denergy : " 10x " : myid="i3" ; icoil="i4" ; id01daf="i2" ; maxcls too small ;")') myid, icoil, id01eaf
   case( 4 ) ; write(ounit,'("denergy : " 10x " : myid="i3" ; icoil="i4" ; id01daf="i2" ; input error ;     ")') myid, icoil, id01eaf
   end select

  enddo ! end of loop over icoil; 14 Apr 16;

  do icoil = 1, Ncoils ; llmodnp = modulo(icoil-1,ncpu)
   RlBCAST( FF(1:Nfun,icoil), Nfun, llmodnp )
  enddo

  SALLOCATE( dL, (1:Ncoils,0:3,0:1,0:NN), zero )
  
  do icoil = 1, Ncoils
   ;kk =      0
   ;kk = kk + 1 ; dL(icoil,0,0, 0) = FF(kk,icoil) ; coil(icoil)%L = dL(icoil,0,0, 0)
   do mm = 0, NN ! derivative wrt xc, yc, zc; 14 Apr 16;
    kk = kk + 1 ; dL(icoil,1,0,mm) = FF(kk,icoil)
    kk = kk + 1 ; dL(icoil,2,0,mm) = FF(kk,icoil)
    kk = kk + 1 ; dL(icoil,3,0,mm) = FF(kk,icoil)
   enddo
   do mm = 1, NN ! derivative wrt xs, ys, zs; 14 Apr 16;
    kk = kk + 1 ; dL(icoil,1,1,mm) = FF(kk,icoil)
    kk = kk + 1 ; dL(icoil,2,1,mm) = FF(kk,icoil)
    kk = kk + 1 ; dL(icoil,3,1,mm) = FF(kk,icoil)
   enddo
   FATAL( denergy, kk.ne.Nfun, counting error )
  enddo ! end of do icoil; 14 Apr 16;

  ttlen = zero

  idof = 0
  do icoil = 1, Ncoils ; Lfactor = coil(icoil)%Lc * coil(icoil)%Lw * ( coil(icoil)%L - coil(icoil)%Lo ) * coil(icoil)%L
   ;             ; idof = idof + 1 ! dE(idof) = dE(idof) +                                                                       ! coil(icoil)%I
   ;  mm = 0     ; idof = idof + 1 ; dE(idof) = dE(idof) + Lfactor * dL(icoil,1,0,mm)                                            ! coil(icoil)%xc(mm) 
   ;             ; idof = idof + 1 ; dE(idof) = dE(idof) + Lfactor * dL(icoil,2,0,mm)                                            ! coil(icoil)%yc(mm) 
   ;             ; idof = idof + 1 ; dE(idof) = dE(idof) + Lfactor * dL(icoil,3,0,mm)                                            ! coil(icoil)%zc(mm) 
   do mm = 1, NN ; idof = idof + 1 ; dE(idof) = dE(idof) + Lfactor * dL(icoil,1,0,mm)                                            ! coil(icoil)%xc(mm) 
    ;            ; idof = idof + 1 ; dE(idof) = dE(idof) + Lfactor * dL(icoil,2,0,mm)                                            ! coil(icoil)%yc(mm) 
    ;            ; idof = idof + 1 ; dE(idof) = dE(idof) + Lfactor * dL(icoil,3,0,mm)                                            ! coil(icoil)%zc(mm) 
    ;            ; idof = idof + 1 ; dE(idof) = dE(idof) + Lfactor * dL(icoil,1,1,mm)                                            ! coil(icoil)%xs(mm) 
    ;            ; idof = idof + 1 ; dE(idof) = dE(idof) + Lfactor * dL(icoil,2,1,mm)                                            ! coil(icoil)%ys(mm) 
    ;            ; idof = idof + 1 ; dE(idof) = dE(idof) + Lfactor * dL(icoil,3,1,mm)                                            ! coil(icoil)%zs(mm) 
   enddo
  ttlen = ttlen + half * coil(icoil)%Lw * ( coil(icoil)%L - coil(icoil)%Lo )**2
  enddo
  ttlen = ttlen / Ncoils
  DALLOCATE( dL )

  DALLOCATE( wk )
  DALLOCATE( FF )
  DALLOCATE( DD )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  dE(1:Ndof) = - dE(1:Ndof) ! energy descent; 14 Apr 16;

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

 !write(ounit,'("denergy : " 10x " : tau ="es12.05" ; EE="es12.05" ; |dE|="es12.05" ;")') tau, totalenergy, sqrt(sum(dE(1:Ndof)**2)/Ndof)

!  knotenergy       = zero
!  FF(1:6*NFcoil+3) = zero
!
!  SALLOCATE( ax, (1:3,0:Npoints), zero )
!  
!  do ii = 0, Npoints
!   aa = - one ; angle = zero ; teta = ii * pi2 / Npoints ; call knotxx( aa, angle, teta, ax(1:3,ii), at(1:3), az(1:3), xx(1:3), xt(1:3), xz(1:3) )
!  enddo
!
!  do ii = 0, Npoints-1 ; teta = ii * pi2 / Npoints
!   do jj = 0, Npoints-1 ; zeta = jj * pi2 / Npoints
!    
!    if( ii.eq.jj ) cycle
!
!    dx(1:3) = ax(1:3,ii) - ax(1:3,jj) ; dxo = sqrt( dx(1)**2 + dx(2)**2 + dx(3)**2 ) ; ddx = - alpha / dxo**(alpha+two)
!    
!    ;knotenergy        = knotenergy        + kpotent * one / dxo**alpha
!    do mm = 0, NFcoil ; dca = cos(mm*teta) - cos(mm*zeta) ; dsa = sin(mm*teta) - sin(mm*zeta)
!     FF(         1+mm) = FF(         1+mm) + kpotent * ddx * dx(1) * dca
!     FF(2*NFcoil+2+mm) = FF(2*NFcoil+2+mm) + kpotent * ddx * dx(2) * dca
!     FF(4*NFcoil+3+mm) = FF(4*NFcoil+3+mm) + kpotent * ddx * dx(3) * dca
!     if( mm.gt.0 ) then
!     FF(  NFcoil+1+mm) = FF(  NFcoil+1+mm) + kpotent * ddx * dx(1) * dsa
!     FF(3*NFcoil+2+mm) = FF(3*NFcoil+2+mm) + kpotent * ddx * dx(2) * dsa
!     FF(5*NFcoil+3+mm) = FF(5*NFcoil+3+mm) + kpotent * ddx * dx(3) * dsa
!     endif
!    enddo ! end of do mm; 14 Apr 16;
!    
!   enddo ! end of do jj; 14 Apr 16;
!  enddo ! end of do ii; 14 Apr 16;
  
!  cxx = zero ; cyy = zero ; czz = zero
!
!  do ii = 0   , Npoints-1 ; teta = ii * pi2 / Npoints
!      jj = ii+1              ; zeta = jj * pi2 / Npoints
!   
!    dx(1:3) = ax(1:3,ii) - ax(1:3,jj) ; dxo =       dx(1)**2 + dx(2)**2 + dx(3)**2
!
!    cxx = cxx + ax(1,ii) * sqrt(dxo) ; cyy = cyy + ax(2,ii) * sqrt(dxo) ; czz = czz + ax(3,ii) * sqrt(dxo) ! knot center; 14 Apr 16;
!    
!    ;knotenergy        = knotenergy        + kspring * dxo * half
!    do mm = 0, NFcoil ; dca = cos(mm*teta) - cos(mm*zeta) ; dsa = sin(mm*teta) - sin(mm*zeta)
!     FF(         1+mm) = FF(         1+mm) + kspring * dx(1) * dca
!     FF(2*NFcoil+2+mm) = FF(2*NFcoil+2+mm) + kspring * dx(2) * dca
!     FF(4*NFcoil+3+mm) = FF(4*NFcoil+3+mm) + kspring * dx(3) * dca
!     if( mm.gt.0 ) then
!     FF(  NFcoil+1+mm) = FF(  NFcoil+1+mm) + kspring * dx(1) * dsa
!     FF(3*NFcoil+2+mm) = FF(3*NFcoil+2+mm) + kspring * dx(2) * dsa
!     FF(5*NFcoil+3+mm) = FF(5*NFcoil+3+mm) + kspring * dx(3) * dsa
!     endif
!    enddo ! end of do mm; 14 Apr 16;
!        
!  enddo ! end of do ii; 14 Apr 16;
  
!  do ii = 0   , Npoints-1 ; teta = ii * pi2 / Npoints
!      jj = ii+1              ; zeta = jj * pi2 / Npoints
!   
!    dx(1:3) = ax(1:3,ii) - ax(1:3,jj) ; dxo = sqrt( dx(1)**2 + dx(2)**2 + dx(3)**2 )
!    
!    do mm = 0, NFcoil ;  ca = cos(mm*teta)                ;  sa = sin(mm*teta)
!     ;                ; dca = cos(mm*teta) - cos(mm*zeta) ; dsa = sin(mm*teta) - sin(mm*zeta)
!     FF(         1+mm) = FF(         1+mm) + kcenter * ( cxx * ( ca * dxo + ax(1,ii) * dx(1) * dca / dxo ) + &
!                                                         cyy * (            ax(2,ii) * dx(1) * dca / dxo ) + &
!                                                         czz * (            ax(3,ii) * dx(1) * dca / dxo ) )
!     FF(2*NFcoil+2+mm) = FF(2*NFcoil+2+mm) + kcenter * ( cxx * (            ax(1,ii) * dx(2) * dca / dxo ) + &
!                                                         cyy * ( ca * dxo + ax(2,ii) * dx(2) * dca / dxo ) + &
!                                                         czz * (            ax(3,ii) * dx(2) * dca / dxo ) )
!     FF(4*NFcoil+3+mm) = FF(4*NFcoil+3+mm) + kcenter * ( cxx * (            ax(1,ii) * dx(3) * dca / dxo ) + &
!                                                         cyy * (            ax(2,ii) * dx(3) * dca / dxo ) + &
!                                                         czz * ( ca * dxo + ax(3,ii) * dx(3) * dca / dxo ) )
!     if( mm.gt.0 ) then
!     FF(  NFcoil+1+mm) = FF(  NFcoil+1+mm) + kcenter * ( cxx * ( sa * dxo + ax(1,ii) * dx(1) * dsa / dxo ) + &
!                                                         cyy * (            ax(2,ii) * dx(1) * dsa / dxo ) + &
!                                                         czz * (            ax(3,ii) * dx(1) * dsa / dxo ) )
!     FF(3*NFcoil+2+mm) = FF(3*NFcoil+2+mm) + kcenter * ( cxx * (            ax(1,ii) * dx(2) * dsa / dxo ) + &
!                                                         cyy * ( sa * dxo + ax(2,ii) * dx(2) * dsa / dxo ) + &
!                                                         czz * (            ax(3,ii) * dx(2) * dsa / dxo ) )
!     FF(5*NFcoil+3+mm) = FF(5*NFcoil+3+mm) + kcenter * ( cxx * (            ax(1,ii) * dx(3) * dsa / dxo ) + &
!                                                         cyy * (            ax(2,ii) * dx(3) * dsa / dxo ) + &
!                                                         czz * ( sa * dxo + ax(3,ii) * dx(3) * dsa / dxo ) )
!     endif
!    enddo ! end of do mm; 14 Apr 16;
!        
!  enddo ! end of do ii; 14 Apr 16;
!
!  DALLOCATE(ax)
!  
!  FF(1:6*NFcoil+3) = - FF(1:6*NFcoil+3)
  return
  
end subroutine denergy

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine dBndxyz( Ndim, tt, Nfun, dFi )
  
  use kmodule, only : zero, one, small, myid, ounit, surf, coil, cmt, smt, icoil, itime, iteta, jzeta, NFcoil
  
  implicit none
  
  include "mpif.h"

  INTEGER              :: Ndim, Nfun
  REAL                 :: tt(1:Ndim), dFi(1:Nfun)

  INTEGER              :: ierr, astat, ifail, ii, NN, kk, mm
  REAL                 :: xx(0:1), yy(0:1), zz(0:1), dx, dy, dz, r1, r2, r3, Bx, By, Bz

  REAL   , allocatable :: ct(:), st(:)

  dFi(1:Nfun) = zero ; NN = NFcoil

  SALLOCATE( ct, (0:NN), zero )
  SALLOCATE( st, (0:NN), zero )
  
#ifdef CONTINUOUS
  call iccoil( tt(1), xx(0:1), yy(0:1), zz(0:1), ifail ) ! icoil is set in denergy; 14 Apr 16;  
  do mm = 0, NN ; ct(mm) = cos(mm*tt(1)) ; st(mm) = sin(mm*tt(1)) ! pre-compute trigonometric factors; 14 Apr 16;
  enddo
  dx = xx(0) - surf(1)%xx(iteta,jzeta)
  dy = yy(0) - surf(1)%yy(iteta,jzeta)
  dz = zz(0) - surf(1)%zz(iteta,jzeta) ; r2 = dx**2 + dy**2 + dz**2 ; r1 = sqrt(r2) ; r3 = r2 * r1  
#else
  do mm = 0, NN ; ct(mm) = cmt(itime,mm) ; st(mm) = smt(itime,mm)
  enddo
  xx(0:1) = (/ coil(icoil)%xx(itime), coil(icoil)%xt(itime) /)
  yy(0:1) = (/ coil(icoil)%yy(itime), coil(icoil)%yt(itime) /)
  zz(0:1) = (/ coil(icoil)%zz(itime), coil(icoil)%zt(itime) /)
#endif
  
  dx = xx(0) - surf(1)%xx(iteta,jzeta)
  dy = yy(0) - surf(1)%yy(iteta,jzeta)
  dz = zz(0) - surf(1)%zz(iteta,jzeta) ; r2 = dx**2 + dy**2 + dz**2 ; r1 = sqrt(r2) ; r3 = r2 * r1
  
  FATAL( dBndxyz, r1.lt.small, too close )
  
  ;kk =      0

  Bx                     = (   dy     * zz(1)     - dz     * yy(1)     ) / r3
  ;kk = kk + 1 ; dFi(kk) =                                                      Bx                       
  do mm = 0, NN ! derivative wrt xc, yc, zc; 14 Apr 16;
   kk = kk + 1 ; dFi(kk) =                                                    - Bx * 3 * dx * ct(mm) / r2
   kk = kk + 1 ; dFi(kk) = (   ct(mm) * zz(1)     + dz     * st(mm)*mm ) / r3 - Bx * 3 * dy * ct(mm) / r2
   kk = kk + 1 ; dFi(kk) = ( - dy     * st(mm)*mm - ct(mm) * yy(1)     ) / r3 - Bx * 3 * dz * ct(mm) / r2
  enddo
  do mm = 1, NN ! derivative wrt xs, ys, zs; 14 Apr 16;
   kk = kk + 1 ; dFi(kk) =                                                    - Bx * 3 * dx * st(mm) / r2
   kk = kk + 1 ; dFi(kk) = (   st(mm) * zz(1)     - dz     * ct(mm)*mm ) / r3 - Bx * 3 * dy * st(mm) / r2
   kk = kk + 1 ; dFi(kk) = (   dy     * ct(mm)*mm - st(mm) * yy(1)     ) / r3 - Bx * 3 * dz * st(mm) / r2
  enddo

  By                     = (   dz     * xx(1)     - dx     * zz(1)     ) / r3
  ;kk = kk + 1 ; dFi(kk) =                                                      By                       
  do mm = 0, NN
   kk = kk + 1 ; dFi(kk) = ( - dz     * st(mm)*mm - ct(mm) * zz(1)     ) / r3 - By * 3 * dx * ct(mm) / r2
   kk = kk + 1 ; dFi(kk) =                                                    - By * 3 * dy * ct(mm) / r2
   kk = kk + 1 ; dFi(kk) = (   ct(mm) * xx(1)     + dx     * st(mm)*mm ) / r3 - By * 3 * dz * ct(mm) / r2
  enddo
  do mm = 1, NN
   kk = kk + 1 ; dFi(kk) = (   dz     * ct(mm)*mm - st(mm) * zz(1)     ) / r3 - By * 3 * dx * st(mm) / r2
   kk = kk + 1 ; dFi(kk) =                                                    - By * 3 * dy * st(mm) / r2
   kk = kk + 1 ; dFi(kk) = (   st(mm) * xx(1)     - dx     * ct(mm)*mm ) / r3 - By * 3 * dz * st(mm) / r2
  enddo

  Bz                     = (   dx     * yy(1)     - dy     * xx(1)     ) / r3
  ;kk = kk + 1 ; dFi(kk) =                                                      Bz                       
  do mm = 0, NN
   kk = kk + 1 ; dFi(kk) = (   ct(mm) * yy(1)     + dy     * st(mm)*mm ) / r3 - Bz * 3 * dx * ct(mm) / r2
   kk = kk + 1 ; dFi(kk) = ( - dx     * st(mm)*mm - ct(mm) * xx(1)     ) / r3 - Bz * 3 * dy * ct(mm) / r2
   kk = kk + 1 ; dFi(kk) =                                                    - Bz * 3 * dz * ct(mm) / r2
  enddo
  do mm = 1, NN
   kk = kk + 1 ; dFi(kk) = (   st(mm) * yy(1)     - dy     * ct(mm)*mm ) / r3 - Bz * 3 * dx * st(mm) / r2
   kk = kk + 1 ; dFi(kk) = (   dx     * ct(mm)*mm - st(mm) * xx(1)     ) / r3 - Bz * 3 * dy * st(mm) / r2
   kk = kk + 1 ; dFi(kk) =                                                    - Bz * 3 * dz * st(mm) / r2
  enddo

  FATAL( dBndxyz, kk.ne.Nfun, counting error )
 
!  Bx                     = (   dy     * zz(1)     - dz     * yy(1)     ) / r3
!  By                     = (   dz     * xx(1)     - dx     * zz(1)     ) / r3
!  Bz                     = (   dx     * yy(1)     - dy     * xx(1)     ) / r3
!
!  do ii = 1, 3
!
!   if    ( ii.eq.1 ) then ; jj = 2 ; kk = 3
!   elseif( ii.eq.2 ) then ; jj = 3 ; kk = 1
!   elseif( ii.eq.3 ) then ; jj = 1 ; kk = 2
!   endif
!
!   dB(ii, 0) = dl(jj, 0) * Dx(kk, 0) - dl(kk, 0) * Dx(jj, 0) ! numerator
!   
!   do dd = 1, 3
!   
!   dB(ii,dd) = dl(jj,dd) * Dx(kk ,0) - dl(kk,dd) * Dx(jj, 0) &
!             + dl(jj, 0) * Dx(kk,dd) - dl(kk, 0) * Dx(jj,dd)
!
!  enddo



  DALLOCATE( ct )
  DALLOCATE( st )

  return
  
end subroutine dBndxyz

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine dlength( Ndim, tt, Nfun, dFi )
  
  use kmodule, only : zero, myid, ounit, NFcoil
  
  implicit none
  
  include "mpif.h"

  INTEGER              :: Ndim, Nfun
  REAL                 :: tt(1:Ndim), dFi(1:Nfun)

  INTEGER              :: ierr, astat, ifail, NN, kk, mm
  REAL                 :: xx(0:1), yy(0:1), zz(0:1), d1, d2
  REAL   , allocatable :: ct(:), st(:)

  dFi(1:Nfun) = zero ; NN = NFcoil
  
  call iccoil( tt(1), xx(0:1), yy(0:1), zz(0:1), ifail ) ! icoil is set in denergy; 14 Apr 16;
  
  d2 = xx(1)**2 + yy(1)**2 + zz(1)**2 ; d1 = sqrt( d2 )

  SALLOCATE( ct, (0:NN), zero )
  SALLOCATE( st, (0:NN), zero )
  do mm = 0, NN ; ct(mm) = cos(mm*tt(1)) ; st(mm) = sin(mm*tt(1)) ! pre-compute trigonometric factors; probably no need for this; 14 Apr 16;
  enddo

  ;kk =      0
  ;kk = kk + 1 ; dFi(kk) =                           d1
  do mm = 0, NN ! derivative wrt xc, yc, zc; 14 Apr 16;
   kk = kk + 1 ; dFi(kk) = xx(1) * ( - st(mm)*mm ) / d1
   kk = kk + 1 ; dFi(kk) = yy(1) * ( - st(mm)*mm ) / d1
   kk = kk + 1 ; dFi(kk) = yy(1) * ( - st(mm)*mm ) / d1
  enddo
  do mm = 1, NN ! derivative wrt xs, ys, zs; 14 Apr 16;
   kk = kk + 1 ; dFi(kk) = xx(1) * ( + ct(mm)*mm ) / d1
   kk = kk + 1 ; dFi(kk) = yy(1) * ( + ct(mm)*mm ) / d1
   kk = kk + 1 ; dFi(kk) = zz(1) * ( + ct(mm)*mm ) / d1
  enddo
  FATAL( dlength, kk.ne.Nfun, counting error )
  
  DALLOCATE( ct )
  DALLOCATE( st )

  return
  
end subroutine dlength

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!


subroutine costfun(nderiv)
  use kmodule, only: zero, sqrtmachprec, myid, ounit, Idisplay, iter, Loptimize, Lnormalize, &
       Ncoils, NFcoil, Cdof, deriv, Ndof, lc, Inorm, Gnorm, &
       totalenergy, t1E, t2E, &
       bnorm      , t1B, t2B, weight_bnorm,  &
       tflux      , t1F, t2F, weight_tflux, target_tflux, isign, &
       ttlen      , t1L, t2L, weight_ttlen, &
       eqarc      , t1A, t2A, weight_eqarc, &
       ccsep      , t1C, t2C, weight_ccsep, &
       qasym      , t1S, t2S, weight_qasym, &
       resbn      , t1R, t2R, weight_resbn

  implicit none
  include "mpif.h"

  INTEGER      :: nderiv

  INTEGER      :: astat, ierr, ii, icl, inf, icoil, jcoil

  REAL         :: start, finish
  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  iter = iter + 1
  if(Idisplay .eq. -2 .and. myid .eq. 0) write(ounit,'("Costfun :"12X": begin the "I6"th run.")') iter
  if(Idisplay .eq. -1) then

#ifdef NORM
    call bnormal(0)
#else
    call bnormal2(0)
#endif

   if ( target_tflux .eq. 0.0 ) then
    call torflux(0)
    target_tflux = isign*sqrt(2.0*tflux)             
    if(myid .eq. 0) write(ounit,'("Costfun :"11X" : Reset target toroidal flux to"ES23.15)') target_tflux
   endif

   call torflux(0)

   if (lc .eq. 1) then
    call tlength(0)
   else
    call tlengthExp(0)
   endif
   
   call specwid(0)
   call coilsep(0)
  endif
  
  totalenergy = zero
  
  if ( nderiv .ge. 1 .and. .not. allocated(t1E) ) ALLOCATE(t1E(1:Ncoils, 0:Cdof))
  if ( nderiv .eq. 2 .and. .not. allocated(t2E) ) ALLOCATE(t2E(1:Ncoils, 0:Cdof, 1:Ncoils, 0:Cdof))

  if     ( nderiv .eq. 1 ) then
   t1E = zero
  elseif ( nderiv .eq. 2 ) then
   t1E = zero; t2E = zero
  endif

  ! call unpack(xdof(1:Ndof))
  call discretecoil


  if (weight_bnorm .gt. sqrtmachprec .or. weight_qasym .gt. sqrtmachprec .or. weight_resbn .gt. sqrtmachprec) then

#ifdef NORM
    call bnormal(nderiv)
#else
    call bnormal2(nderiv)
#endif
    
  endif

  if (weight_bnorm .gt. sqrtmachprec) then
   totalenergy = totalenergy + weight_bnorm * bnorm
   if     ( nderiv .eq. 1 ) then
    t1E = t1E +  weight_bnorm * t1B
   elseif ( nderiv .eq. 2 ) then
    t1E = t1E +  weight_bnorm * t1B
    t2E = t2E +  weight_bnorm * t2B
   endif
  endif

  if (weight_qasym .gt. sqrtmachprec) then
   totalenergy = totalenergy + weight_qasym * qasym
   if     ( nderiv .eq. 1 ) then
    t1E = t1E +  weight_qasym * t1S
   elseif ( nderiv .eq. 2 ) then
    t1E = t1E +  weight_qasym * t1S
    t2E = t2E +  weight_qasym * t2S
   endif
  endif

  if (weight_resbn .gt. sqrtmachprec) then
   totalenergy = totalenergy + weight_resbn * resbn
   if     ( nderiv .eq. 1 ) then
    t1E = t1E +  weight_resbn * t1R
   elseif ( nderiv .eq. 2 ) then
    t1E = t1E +  weight_resbn * t1R
    t2E = t2E +  weight_resbn * t2R
   endif
  endif

  ! if(myid .eq. 0) write(ounit,'("calling bnormal used",f10.5,"seconds.")') finish-start

  if (weight_tflux .gt. sqrtmachprec) then

   if ( target_tflux .eq. 0.0 ) then
    call torflux(0)
    target_tflux = isign*sqrt(2.0*tflux)             
    if(myid .eq. 0) write(ounit,'("Costfun :"11X" : Reset target toroidal flux to"ES23.15)') target_tflux
   endif

   call torflux(nderiv)
   totalenergy = totalenergy + weight_tflux * tflux / target_tflux**2 ! normalization
   if     ( nderiv .eq. 1 ) then
    t1E = t1E +  weight_tflux * t1F / target_tflux**2
   elseif ( nderiv .eq. 2 ) then
    t1E = t1E +  weight_tflux * t1F / target_tflux**2
    t2E = t2E +  weight_tflux * t2F / target_tflux**2
   endif

  endif

  ! if(myid .eq. 0) write(ounit,'("calling torflux used",f10.5,"seconds.")') start-finish

  if (weight_ttlen .gt. sqrtmachprec) then

   if( lc .eq. 1 ) then 
    call tlength(nderiv)
   elseif (lc .eq. 2) then
    call tlengthExp(nderiv)
   else
    FATAL( dnergy, .true., Conflicts between lc and weight_ttlen )
   !call MPI_ABORT( MPI_COMM_WORLD, 1, ierr )
   !stop "COSTFUN: Conflicts between lc and weight_ttlen"
   endif
   
   totalenergy = totalenergy + weight_ttlen * ttlen
   if     ( nderiv .eq. 1 ) then
    t1E = t1E +  weight_ttlen * t1L
   elseif ( nderiv .eq. 2 ) then
    t1E = t1E +  weight_ttlen * t1L
    t2E = t2E +  weight_ttlen * t2L
   endif

  endif

  ! if(myid .eq. 0) write(ounit,'("calling tlength used",f10.5,"seconds.")') finish-start

  if (weight_eqarc .ge. sqrtmachprec) then

  ! call equarcl(nderiv)
  ! call specwid_df(nderiv)
   !if ( Loptimize .eq. 3 ) then
    !call langrange(nderiv)
   !else
    call specwid(nderiv)
    !call specwid_df(nderiv)
  ! endif
   
   totalenergy = totalenergy + weight_eqarc * eqarc
   if     ( nderiv .eq. 1 ) then
    t1E = t1E + weight_eqarc * t1A
   elseif ( nderiv .eq. 2 ) then
    t1E = t1E + weight_eqarc * t1A
    t2E = t2E + weight_eqarc * t2A
   endif

  endif

  if (weight_ccsep .ge. sqrtmachprec) then

   call coilsep(nderiv)
   totalenergy = totalenergy + weight_ccsep * ccsep
   if     ( nderiv .eq. 1 ) then
    t1E = t1E + weight_ccsep * t1C
   elseif ( nderiv .eq. 2 ) then
    t1E = t1E + weight_ccsep * t1C
    t2E = t2E + weight_ccsep * t2C
   endif

  endif
  
  ! normalization
  if (nderiv .ge. 1) then
     do icoil = 1, Ncoils
        t1E(icoil, 0) = t1E(icoil, 0)*Inorm
        t1E(icoil, 1:Cdof) = t1E(icoil, 1:Cdof)*Gnorm
     enddo
  endif
  if (nderiv .ge. 2) then
     do icoil = 1, Ncoils
        do jcoil = 1, Ncoils
           t2E(icoil,      0, jcoil,      0) = t2E(icoil,      0, jcoil,      0)*Inorm*Inorm
           t2E(icoil,      0, jcoil, 1:Cdof) = t2E(icoil,      0, jcoil, 1:Cdof)*Inorm*Gnorm
           t2E(icoil, 1:Cdof, jcoil,      0) = t2E(icoil, 1:Cdof, jcoil,      0)*Gnorm*Inorm
           t2E(icoil, 1:Cdof, jcoil, 1:Cdof) = t2E(icoil, 1:Cdof, jcoil, 1:Cdof)*Gnorm*Gnorm
        enddo
     enddo
  endif
  

  if(allocated(deriv)) then
     deriv = zero
     do ii = 1, Ndof
        !if (Loptimize .eq. 3) then
        !   call ntconvert(ii, icl, inf)
        !else
           call DoFconvert(ii, icl, inf)
        !endif
        if(allocated(t1E)) deriv(ii, 0) = t1E(icl, inf)
        if(allocated(t1B)) deriv(ii, 1) = t1B(icl, inf)
        if(allocated(t1F)) deriv(ii, 2) = t1F(icl, inf)
        if(allocated(t1L)) deriv(ii, 3) = t1L(icl, inf)
        if(allocated(t1A)) deriv(ii, 4) = t1A(icl, inf)
        if(allocated(t1C)) deriv(ii, 5) = t1C(icl, inf)
     enddo
  endif

  if( allocated(t1B) ) deallocate(t1B)
  if( allocated(t1F) ) deallocate(t1F)
  if( allocated(t1L) ) deallocate(t1L)
  if( allocated(t1A) ) deallocate(t1A)
  if( allocated(t1C) ) deallocate(t1C)
  if( allocated(t1S) ) deallocate(t1S)
  if( allocated(t1R) ) deallocate(t1R)

  if( allocated(t2B) ) deallocate(t2B)
  if( allocated(t2F) ) deallocate(t2F)
  if( allocated(t2L) ) deallocate(t2L)
  if( allocated(t2A) ) deallocate(t2A)
  if( allocated(t2C) ) deallocate(t2C)
  if( allocated(t2S) ) deallocate(t2S)
  if( allocated(t2R) ) deallocate(t2R)

  FATAL( denergy, iter.gt.100000, too many iterations )
  !if(iter .ge. 1E5) call MPI_ABORT( MPI_COMM_WORLD, 10, ierr )

  return
end subroutine costfun


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine denergy2( tau, xdof, dE )
  
  use kmodule, only : NFcoil, coil, icoil, Ncoils, Ndof, myid, ncpu, ounit, t1E

  
  implicit none
  
  include "mpif.h"
  
  REAL                 :: tau, xdof(*), dE(*)
  
  INTEGER              :: ierr, astat, idof, NN, mm
  
  external             :: costfun
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  NN = NFcoil ! shorthand; 

  call unpack(xdof(1:Ndof))

  call costfun(1)
     idof = 0
     do icoil = 1, Ncoils

        if(coil(icoil)%Ic .ne. 0) then
           idof = idof + 1 ; dE(idof) = t1E(icoil,        0) ! coil(icoil)%I
        endif
 
        if(coil(icoil)%Lc .ne. 0) then
           idof = idof + 1 ; dE(idof) = t1E(icoil,        1) ! coil(icoil)%xc( 0)
           idof = idof + 1 ; dE(idof) = t1E(icoil,   2*NN+3) ! coil(icoil)%yc( 0)
           idof = idof + 1 ; dE(idof) = t1E(icoil,   4*NN+5) ! coil(icoil)%zc( 0)
          do mm = 1, NN 
           idof = idof + 1 ; dE(idof) = t1E(icoil,mm+     1) ! coil(icoil)%xc(mm)
           idof = idof + 1 ; dE(idof) = t1E(icoil,mm+2*NN+3) ! coil(icoil)%yc(mm)
           idof = idof + 1 ; dE(idof) = t1E(icoil,mm+4*NN+5) ! coil(icoil)%zc(mm)
           idof = idof + 1 ; dE(idof) = t1E(icoil,mm+  NN+2) ! coil(icoil)%xs(mm)
           idof = idof + 1 ; dE(idof) = t1E(icoil,mm+3*NN+4) ! coil(icoil)%ys(mm)
           idof = idof + 1 ; dE(idof) = t1E(icoil,mm+5*NN+6) ! coil(icoil)%zs(mm)
          enddo ! end of do mm
        endif

    enddo ! end of do icoil
    FATAL( denergy, idof.ne.Ndof, counting error )

    DALLOCATE( t1E )
    dE(1:Ndof) = - dE(1:Ndof)
    return
  end subroutine denergy2

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine weightnormalize
  use kmodule, only:  zero, Io,  coil, Ncoils,  &
                      weight_bnorm, weight_tflux, weight_ttlen, weight_eqarc, weight_ccsep, weight_qasym, weight_resbn, &
                             bnorm,        tflux,        ttlen,        eqarc,        ccsep,        qasym,        resbn, &
                     target_tflux, isign, Lnormalize, sqrtmachprec, myid, ounit, lc          
  implicit none
  include "mpif.h"


  INTEGER    :: ierr, icoil
  REAL       :: tmp, cur_tflux
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  call discretecoil


!-!-!-!-!-!-!-!-!-!-tflux-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( weight_tflux .ge. sqrtmachprec ) then

   if ( target_tflux .eq. 0.0 ) then
    call torflux(0)
    target_tflux = isign*sqrt(2.0*tflux)             
    if(myid .eq. 0) write(ounit,'("costfun :"11X" : Reset target toroidal flux to"ES23.15)') target_tflux
   else
      tmp = target_tflux
      target_tflux = zero
      call torflux(0)
      target_tflux = tmp
      cur_tflux = isign*sqrt(2.0*tflux)
      Io = Io * target_tflux / cur_tflux
      do icoil = 1, Ncoils
         coil(icoil)%I = Io
         coil(icoil)%Io = Io
      enddo
      if(myid .eq. 0) write(ounit,'("costfun :"11X" : rescale coil currents with a factor of"ES23.15)') target_tflux / cur_tflux
   endif

   call torflux(0)
   if (abs(tflux) .gt. sqrtmachprec) weight_tflux = weight_tflux / tflux * target_tflux**2
   if( myid .eq. 0 ) write(ounit, 1000) "weight_tflux", weight_tflux
   
  endif  


!-!-!-!-!-!-!-!-!-!-bnorm-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if (weight_bnorm .gt. sqrtmachprec .or. weight_qasym .gt. sqrtmachprec .or. weight_resbn .gt. sqrtmachprec) then
  ! if( weight_bnorm .ge. sqrtmachprec ) then

#ifdef NORM
     call bnormal(0)
#else
     call bnormal2(0)
#endif
  
     if( weight_bnorm .ge. sqrtmachprec ) then
        if (abs(bnorm) .gt. sqrtmachprec) weight_bnorm = weight_bnorm / bnorm
        if( myid .eq. 0 ) write(ounit, 1000) "weight_bnorm", weight_bnorm
     endif

     if( weight_qasym .ge. sqrtmachprec ) then
        if (abs(qasym) .gt. sqrtmachprec) weight_qasym = weight_qasym / qasym
        if( myid .eq. 0 ) write(ounit, 1000) "weight_qasym", weight_qasym
     endif
   
     if( weight_resbn .ge. sqrtmachprec ) then
        if (abs(resbn) .gt. sqrtmachprec) weight_resbn = weight_resbn / resbn
        if( myid .eq. 0 ) write(ounit, 1000) "weight_resbn", weight_resbn
     endif

  endif

!-!-!-!-!-!-!-!-!-!-ttlen-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( weight_ttlen .ge. sqrtmachprec ) then

   if( lc .eq. 1 ) then 
    call tlength(0)
   elseif(lc .eq. 2) then
    call tlengthExp(0)
   else
    FATAL( dnergy, .true., Conflicts between lc and weight_ttlen )
   !call MPI_ABORT( MPI_COMM_WORLD, 1, ierr )
   !stop "Weights_normalize: Conflicts between lc and weight_ttlen"
   endif

   if (abs(ttlen) .gt. sqrtmachprec) weight_ttlen = weight_ttlen / ttlen
   if( myid .eq. 0 ) write(ounit, 1000) "weight_ttlen", weight_ttlen
   
  endif 

!-!-!-!-!-!-!-!-!-!-eqarc-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( weight_eqarc .ge. sqrtmachprec ) then

   call specwid(0)
   if (abs(eqarc) .gt. sqrtmachprec) weight_eqarc = weight_eqarc / eqarc
   if( myid .eq. 0 ) write(ounit, 1000) "weight_eqarc", weight_eqarc
   
  endif 

!-!-!-!-!-!-!-!-!-!-ccsep-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if( weight_ccsep .ge. sqrtmachprec ) then

   call coilsep(0)
   if (abs(ccsep) .gt. sqrtmachprec) weight_ccsep = weight_ccsep / ccsep
   if( myid .eq. 0 ) write(ounit, 1000) "weight_ccsep", weight_ccsep
   
  endif 
1000 format("wnormal : " 10x " : "a12" is normalized to" ES23.15)


  !Lnormalize = 0          ! turn back;
  return

end subroutine weightnormalize
