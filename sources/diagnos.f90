!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

SUBROUTINE diagnos
!------------------------------------------------------------------------------------------------------ 
! DATE: 07/13/2017
! diagonose the coil performance
!------------------------------------------------------------------------------------------------------   

  use globals
  use mpi
  implicit none

  INTEGER           :: icoil, itmp, NF, idof, i, j, isurf, cs, ip, is, Npc, coilInd0, coilInd1, j0, per0, l0, ss0, iseg, NS, NCP

  LOGICAL           :: lwbnorm, l_raw
  REAL              :: MaxCurv, AvgLength, MinCCdist, MinCPdist, tmp_dist, ReDot, ImDot, dum, AvgCurv, AvgTors, &
                       MinLambda, MaxS, torsRet, Bxhold, Byhold, Bzhold, xc, yc, zc, C
  REAL, parameter   :: infmax = 1.0E6
  REAL, allocatable :: Atmp(:,:), Btmp(:,:), absrp(:), deltax(:), deltay(:), deltaz(:), norm(:)
  REAL, allocatable :: nx(:), ny(:), nz(:), npx(:), npy(:), npz(:), nhatpx(:), nhatpy(:), nhatpz(:), normt(:), normn(:), zeta(:), alpha(:), thatx(:), &
     thaty(:), thatz(:), nhatx(:), nhaty(:), nhatz(:), bhatx(:), bhaty(:), bhatz(:), thatpx(:), thatpy(:), thatpz(:), &
     nhatax(:), nhatay(:), nhataz(:), bhatax(:), bhatay(:), bhataz(:), alphac(:), alphas(:)

  isurf = plasma
  itmp = 0
  lwbnorm = .True. 
  l_raw = .False. ! if use raw coils data
  if (myid == 0 .and. IsQuiet < 0) write(ounit, *) "-----------COIL DIAGNOSTICS----------------------------------"

  !--------------------------------cost functions-------------------------------------------------------  
  if (case_optimize == 0) call AllocData(0) ! if not allocate data;
  call costfun(0)

  if (myid == 0) then
      write(ounit, '("diagnos : Individual cost function values: ")') 
      write(ounit, 100) "Bnormal", bnorm
      write(ounit, 100) "B_mn harmonics", bharm
      write(ounit, 100) "Toroidal flux", tflux
      write(ounit, 100) "Total current penalty", tflux
      write(ounit, 100) "Coil length", ttlen
      write(ounit, 100) "Coil curvature", curv
      write(ounit, 100) "Coil torsion", tors
      write(ounit, 100) "Nissin complexity", nissin
      write(ounit, 100) "Coil-coil separation", ccsep
      write(ounit, 100) "Coil-surf separation", cssep
      write(ounit, 100) "Stochastic Bnormal", bnormavg
      write(ounit, 100) "Straight-out", str
  endif 

100  format (8X, ": ", A20, " = ", ES15.8) 

  !save all the coil parameters;
  if (allocated(coilspace)) then
     idof = 0
     do icoil = 1, Ncoils
        coilspace(iout, idof+1 ) = coil(icoil)%I ;  idof = idof + 1

        select case (coil(icoil)%type)
        case (1)
           NF = FouCoil(icoil)%NF
           coilspace(iout, idof+1:idof+NF+1) = FouCoil(icoil)%xc(0:NF) ; idof = idof + NF +1
           coilspace(iout, idof+1:idof+NF  ) = FouCoil(icoil)%xs(1:NF) ; idof = idof + NF
           coilspace(iout, idof+1:idof+NF+1) = FouCoil(icoil)%yc(0:NF) ; idof = idof + NF +1
           coilspace(iout, idof+1:idof+NF  ) = FouCoil(icoil)%ys(1:NF) ; idof = idof + NF
           coilspace(iout, idof+1:idof+NF+1) = FouCoil(icoil)%zc(0:NF) ; idof = idof + NF +1
           coilspace(iout, idof+1:idof+NF  ) = FouCoil(icoil)%zs(1:NF) ; idof = idof + NF
	     case (coil_type_spline)
           NCP = Splines(icoil)%NCP
           coilspace(iout, idof+1:idof+NCP*3) = Splines(icoil)%Cpoints(0:3*NCP-1) ; idof = idof + 3*NCP 
!!$     case default
!!$           FATAL(descent, .true., not supported coil types)
        end select
     enddo
     FATAL( output , idof .ne. Tdof, counting error in restart )
  endif
  !-------------------------------average coil length-------------------------------------------------------  
  AvgLength = zero
  if ( (case_length == 1) .and. (sum(coil(1:Ncoils)%Lo) < sqrtmachprec) ) coil(1:Ncoils)%Lo = one
  call length(0)
  do icoil = 1, Ncoils
     if(coil(icoil)%type .ne. 1 .AND. (coil(icoil)%type .ne. coil_type_spline)) cycle ! only for Fourier
     AvgLength = AvgLength + coil(icoil)%L
  enddo
  AvgLength = AvgLength / Ncoils
  if(myid .eq. 0) write(ounit, '(8X": Average length of the coils is"8X" :" ES23.15)') AvgLength

  !-------------------------------coil maximum curvature----------------------------------------------------  
  MaxCurv = zero
  do icoil = 1, Ncoils
     if(coil(icoil)%type .ne. 1 .AND. (coil(icoil)%type .ne. coil_type_spline)) cycle ! only for Fourier
     call CurvDeriv0(icoil,dum) !dummy return
     if (coil(icoil)%maxcurv .ge. MaxCurv) then
        MaxCurv = coil(icoil)%maxcurv
        itmp = icoil !record the number
     endif
#ifdef DEBUG
     if(myid .eq. 0) write(ounit, '(8X": Maximum curvature of "I3 "-th coil is : " ES23.15)') &
          icoil, coil(icoil)%maxcurv
#endif
  enddo
  if(myid .eq. 0) write(ounit, '(8X": Maximum curvature of all the coils is  :" ES23.15 &
       " ; at coil " I3)') MaxCurv, itmp

  !-------------------------------average coil curvature-------------------------------------------------------
  AvgCurv = zero
  do icoil = 1, Ncoils
     if(coil(icoil)%type .ne. 1 .AND. (coil(icoil)%type .ne. coil_type_spline)) exit ! only for Fourier
     call avgcurvature(icoil)
     AvgCurv = AvgCurv + coil(icoil)%avgcurv
#ifdef DEBUG
     if(myid .eq. 0) write(ounit, '(8X": Average curvature of "I3 "-th coil is : " ES23.15)') &
          icoil, coil(icoil)%avgcurv
#endif
  enddo
  AvgCurv = AvgCurv / Ncoils
  if(myid .eq. 0) write(ounit, '(8X": Average curvature of the coils is"5X" :" ES23.15)') AvgCurv

  !-------------------------------average coil torsion-----------------------------------------------------
  AvgTors = zero
  do icoil = 1, Ncoils
     if(coil(icoil)%type .ne. 1) exit ! only for Fourier 
     call TorsDeriv0(icoil,torsRet,dum) !dummy return
     AvgTors = AvgTors + abs(torsRet) 
#ifdef DEBUG
     if(myid .eq. 0) write(ounit, '(8X": Average torsion of "I3 "-th coil is   : " ES23.15)') &
        icoil, abs(torsRet)
#endif
  enddo
  AvgTors = AvgTors / Ncoils
  if(myid .eq. 0) write(ounit, '(8X": Average torsion of the coils is"5X"   :" ES23.15)') AvgTors

  !-------------------------------maximum coil c------------------------------------------------------------
  if( weight_nissin > 0.0_dp ) then
  call nissinDeriv0(1,dum)
  MaxS = coil(1)%maxs
  do icoil = 1, Ncoils
     if(coil(icoil)%type .ne. 1) exit ! only for Fourier
     call nissinDeriv0(icoil,dum) !dummy return
     if( coil(icoil)%maxs .ge. MaxS) then
        MaxS = coil(icoil)%maxs
     endif
#ifdef DEBUG
     if(myid .eq. 0) write(ounit, '(8X": Maximum c of "I3 "-th coil is         : " ES23.15)') &
        icoil, coil(icoil)%maxs
#endif
  enddo
  if(myid .eq. 0) write(ounit, '(8X": Maximum c of the coils is"5X"         :" ES23.15)') MaxS
  endif


  !-----------------------------minimum coil coil separation------------------------------------  
  ! coils are supposed to be placed in order
  minCCdist = infmax
  coilInd0 = 0
  coilInd1 = 1
  do icoil = 1, Ncoils
     if(coil(icoil)%type .ne. 1 .AND. (coil(icoil)%type .ne. coil_type_spline)) cycle
     if(Ncoils .eq. 1) exit ! if only one coil
     ! Data for the first coil
     SALLOCATE(Atmp, (1:3,0:coil(icoil)%NS-1), zero)
     Atmp(1, 0:coil(icoil)%NS-1) = coil(icoil)%xx(0:coil(icoil)%NS-1)
     Atmp(2, 0:coil(icoil)%NS-1) = coil(icoil)%yy(0:coil(icoil)%NS-1)
     Atmp(3, 0:coil(icoil)%NS-1) = coil(icoil)%zz(0:coil(icoil)%NS-1)
     ! Check distances for coil and its symmetric coils 
     if ( coil(icoil)%symm /= 0 ) then
        SALLOCATE(Btmp, (1:3,0:coil(icoil)%NS-1), zero)
        if ( coil(icoil)%symm .eq. 1 ) then
           cs = 0
           Npc = Nfp
        elseif ( coil(icoil)%symm .eq. 2 ) then
           cs = 1
           Npc = Nfp
        endif
        do ip = 1, Npc
           do is = 0, cs
              if ( ip .eq. 1 .and. is .eq. 0 ) cycle
              Btmp(1, 0:coil(icoil)%NS-1) = (coil(icoil)%xx(0:coil(icoil)%NS-1)*cosnfp(ip) &
                                        & - coil(icoil)%yy(0:coil(icoil)%NS-1)*sinnfp(ip) ) 
              Btmp(2, 0:coil(icoil)%NS-1) = (-1)**is * (coil(icoil)%xx(0:coil(icoil)%NS-1)*sinnfp(ip) &
                                                   & + coil(icoil)%yy(0:coil(icoil)%NS-1)*cosnfp(ip) )
              Btmp(3, 0:coil(icoil)%NS-1) = (-1)**is * (coil(icoil)%zz(0:coil(icoil)%NS-1))
              call mindist(Atmp, coil(icoil)%NS, Btmp, coil(icoil)%NS, tmp_dist)
#ifdef DEBUG
              if(myid .eq. 0) write(ounit, '(8X": distance between  "I3.3"-th and "I3.3"-th coil (ip="I2.2 &
                   ", is="I1") is : " ES23.15)') icoil, icoil, ip, is, tmp_dist
#endif
              if (minCCdist .ge. tmp_dist) then
                 minCCdist = tmp_dist
                 coilInd0 = icoil
                 coilInd1 = icoil
              endif
           enddo
        enddo
        DALLOCATE(Btmp)
     endif
     do itmp = 1, Ncoils
        ! skip self and non-Fourier coils
        if (itmp == icoil .or. ((coil(icoil)%type .ne. 1) .AND. (coil(icoil)%type .ne. coil_type_spline))) cycle
        SALLOCATE(Btmp, (1:3,0:coil(itmp)%NS-1), zero)
        ! check if the coil is stellarator symmetric
        !select case (coil(icoil)%symm)
        select case (coil(itmp)%symm) 
        case ( 0 )
           cs  = 0
           Npc = 1
        case ( 1 )
           cs  = 0
           Npc = Nfp
        case ( 2 ) 
           cs  = 1
           Npc = Nfp
        end select
        ! periodicity and stellarator symmetry
        do ip = 1, Npc
           do is = 0, cs
              Btmp(1, 0:coil(itmp)%NS-1) =            (coil(itmp)%xx(0:coil(itmp)%NS-1)*cosnfp(ip) &
                                                   & - coil(itmp)%yy(0:coil(itmp)%NS-1)*sinnfp(ip) )
              Btmp(2, 0:coil(itmp)%NS-1) = (-1)**is * (coil(itmp)%xx(0:coil(itmp)%NS-1)*sinnfp(ip) &
                                                   & + coil(itmp)%yy(0:coil(itmp)%NS-1)*cosnfp(ip) )
              Btmp(3, 0:coil(itmp)%NS-1) = (-1)**is * (coil(itmp)%zz(0:coil(itmp)%NS-1))
              call mindist(Atmp, coil(icoil)%NS, Btmp, coil(itmp)%NS, tmp_dist)
#ifdef DEBUG
              if(myid .eq. 0) write(ounit, '(8X": distance between  "I3.3"-th and "I3.3"-th coil (ip="I2.2 & 
                   ", is="I1") is : " ES23.15)') icoil, itmp, ip, is, tmp_dist
#endif
              if (minCCdist .ge. tmp_dist) then 
                 minCCdist = tmp_dist
                 coilInd0 = icoil
                 coilInd1 = itmp
              endif
           enddo
        enddo
        DALLOCATE(Btmp)
     enddo
     DALLOCATE(Atmp)
  enddo

  if(myid .eq. 0) write(ounit, '(8X": The minimum coil-coil distance is "4X" :" ES23.15" ; at coils"I3" &
          ,"I3)') minCCdist, coilInd0, coilInd1

  !--------------------------------minimum coil plasma separation-------------------------------  
  minCPdist = infmax
  do icoil = 1, Ncoils

     if(coil(icoil)%type .ne. 1  .AND. (coil(icoil)%type .ne. coil_type_spline)) cycle

     SALLOCATE(Atmp, (1:3,0:coil(icoil)%NS-1), zero)
     SALLOCATE(Btmp, (1:3,1:(Nteta*Nzeta)), zero)

     if ( coil(icoil)%symm == 0 ) then
        per0 = 1
        ss0 = 0
     elseif ( coil(icoil)%symm == 1 ) then
        per0 = Nfp
        ss0 = 0
     elseif ( coil(icoil)%symm == 2 ) then
        per0 = Nfp
        ss0 = 1
     else
        FATAL( diagnos, .true. , Errors in coil symmetry )
     endif
     do j0 = 1, per0
        do l0 = 0, ss0
           Atmp(1, 0:coil(icoil)%NS-1) = (coil(icoil)%xx(0:coil(icoil)%NS-1))*cos(pi2*(j0-1)/Nfp) - (coil(icoil)%yy(0:coil(icoil)%NS-1))*sin(pi2*(j0-1)/Nfp)
           Atmp(2, 0:coil(icoil)%NS-1) = ((-1.0)**l0)*((coil(icoil)%yy(0:coil(icoil)%NS-1))*cos(pi2*(j0-1)/Nfp) + (coil(icoil)%xx(0:coil(icoil)%NS-1))*sin(pi2*(j0-1)/Nfp)) 
           Atmp(3, 0:coil(icoil)%NS-1) = (coil(icoil)%zz(0:coil(icoil)%NS-1))*((-1.0)**l0)

           Btmp(1, 1:(Nteta*Nzeta)) = reshape(surf(isurf)%xx(0:Nteta-1, 0:Nzeta-1), (/Nteta*Nzeta/))
           Btmp(2, 1:(Nteta*Nzeta)) = reshape(surf(isurf)%yy(0:Nteta-1, 0:Nzeta-1), (/Nteta*Nzeta/))
           Btmp(3, 1:(Nteta*Nzeta)) = reshape(surf(isurf)%zz(0:Nteta-1, 0:Nzeta-1), (/Nteta*Nzeta/))

           call mindist(Atmp, coil(icoil)%NS, Btmp, Nteta*Nzeta, tmp_dist)

           if (minCPdist .ge. tmp_dist) then 
              minCPdist=tmp_dist
              itmp = icoil
           endif
        enddo
     enddo

     DALLOCATE(Atmp)
     DALLOCATE(Btmp)

  enddo

  if(myid .eq. 0) then
     write(ounit, '(8X": The minimum coil-plasma distance is    :" ES23.15 &
          " ; at coil " I3)') minCPdist, itmp
     if (minCPdist < 0.1) then
        write(ounit, '(8X":-----------WARNING!!!-------------------------")')
        write(ounit, '(8X": Initial coils might intersect with the plasma")')
        write(ounit, '(8X":----------------------------------------------")') 
     endif
  endif
  
  !--------------------------------overlap of Bn_mn harmonics-----------------------------------  
  ReDot = zero ; ImDot = zero
  if (weight_bharm > sqrtmachprec) then  ! do not care n,m; use cos - i sin
     ReDot = sum(Bmnc*tBmnc) + sum(Bmns*tBmns)
     ImDot = sum(Bmnc*tBmns) - sum(Bmns*tBmnc)
     overlap = sqrt( (ReDot*ReDot + ImDot*ImDot) &
                 / ( (sum(Bmnc*Bmnc)+sum(Bmns*Bmns)) * (sum(tBmnc*tBmnc)+sum(tBmns*tBmns)) ) )
     if(myid .eq. 0) write(ounit, '(8X": Overlap of the realized Bn harmonics is:" F8.3 " %")') 100*overlap
  endif

  !--------------------------------calculate the average Bn error-------------------------------
  if (allocated(surf(isurf)%bn)) then
     ! \sum{ |Bn| / |B| }/ (Nt*Nz)
     if(myid .eq. 0) then
        write(ounit, '(8X": Ave. relative absolute Bn error |Bn|/B : " ES12.5"; max(|Bn|)="ES12.5)') &
             sum(abs(surf(plasma)%bn/sqrt(surf(plasma)%Bx**2+surf(plasma)%By**2+surf(plasma)%Bz**2))) &
             / (Nteta*Nzeta), maxval(abs(surf(plasma)%bn))
        write(ounit, '(8X": Surface area normalized Bn error int(|Bn|/B*ds)/A : "ES23.15)') &
             sum(abs(surf(plasma)%bn)/sqrt(surf(plasma)%Bx**2+surf(plasma)%By**2+surf(plasma)%Bz**2) &
             *surf(plasma)%ds)*discretefactor/(surf(plasma)%area/(surf_Nfp*2**symmetry))
     endif
  endif

  !--------------------------------calculate the stochastic Bn error----------------------------
  if ( Npert .ge. 1 .and. allocated(surf(isurf)%bn) ) then
     
     if (case_optimize .eq. 0) call perturbation(0)
     call sbnormal( 0 )

     if(myid .eq. 0) write(ounit, '(8X": Maximum field error after perturbations: "ES23.15)') bnormmax
     if(myid .eq. 0) write(ounit, '(8X": Average field error after perturbations: "ES23.15)') bnormavg

  endif

  !--------------------------------calculate filamentary body force-----------------------------
  
  if ( filforce .eq. 1 ) then
     ! Calculation can be parallelized 
     do icoil = 1, Ncoils
        ! Skip type .ne. 1 due to datalloc, same with finite-build below
        if ( coil(icoil)%type .ne. 1 ) cycle
        do iseg = 0, coil(icoil)%NS
           do i = 1, Ncoils
              call bfield0(i, coil(icoil)%xx(iseg), coil(icoil)%yy(iseg), coil(icoil)%zz(iseg), Bxhold, Byhold, Bzhold)
              coil(icoil)%Bxx(iseg) = coil(icoil)%Bxx(iseg) + Bxhold
              coil(icoil)%Byy(iseg) = coil(icoil)%Byy(iseg) + Byhold
              coil(icoil)%Bzz(iseg) = coil(icoil)%Bzz(iseg) + Bzhold
           enddo
        enddo
        coil(icoil)%Fx = coil(icoil)%I*(coil(icoil)%yt*coil(icoil)%Bzz-coil(icoil)%zt*coil(icoil)%Byy)
        coil(icoil)%Fy = coil(icoil)%I*(coil(icoil)%zt*coil(icoil)%Bxx-coil(icoil)%xt*coil(icoil)%Bzz)
        coil(icoil)%Fz = coil(icoil)%I*(coil(icoil)%xt*coil(icoil)%Byy-coil(icoil)%yt*coil(icoil)%Bxx)
        coil(icoil)%Fx = coil(icoil)%Fx * coil(icoil)%dd
        coil(icoil)%Fy = coil(icoil)%Fy * coil(icoil)%dd
        coil(icoil)%Fz = coil(icoil)%Fz * coil(icoil)%dd
     enddo
     ! Include loop that calculates max force per unit length 
     !if(myid .eq. 0) write(ounit, '(8X": Maximum Force per Unit Length"8X" :"ES23.15)') mfpul
  endif

  !--------------------------------Calculate Finite-Build Frame---------------------------------
  if ( calcfb .eq. 1 ) then
     do icoil = 1, Ncoils
        NS = coil(icoil)%NS
        SALLOCATE( absrp, (0:NS-1), 0.0 )
        SALLOCATE( deltax, (0:NS), 0.0 )
        SALLOCATE( deltay, (0:NS), 0.0 )
        SALLOCATE( deltaz, (0:NS), 0.0 )
        SALLOCATE( norm, (0:NS), 0.0 )
        SALLOCATE( nx, (0:NS), 0.0 )
        SALLOCATE( ny, (0:NS), 0.0 )
        SALLOCATE( nz, (0:NS), 0.0 )
        SALLOCATE( npx, (0:NS), 0.0 )
        SALLOCATE( npy, (0:NS), 0.0 )
        SALLOCATE( npz, (0:NS), 0.0 )
        SALLOCATE( nhatpx, (0:NS), 0.0 )
        SALLOCATE( nhatpy, (0:NS), 0.0 )
        SALLOCATE( nhatpz, (0:NS), 0.0 )
        SALLOCATE( normt, (0:NS), 0.0 )
        SALLOCATE( normn, (0:NS), 0.0 )
        SALLOCATE( zeta, (0:NS), 0.0 )
        SALLOCATE( thatx, (0:NS), 0.0 )
        SALLOCATE( thaty, (0:NS), 0.0 )
        SALLOCATE( thatz, (0:NS), 0.0 )
        SALLOCATE( nhatx, (0:NS), 0.0 )
        SALLOCATE( nhaty, (0:NS), 0.0 )
        SALLOCATE( nhatz, (0:NS), 0.0 )
        SALLOCATE( bhatx, (0:NS), 0.0 )
        SALLOCATE( bhaty, (0:NS), 0.0 )
        SALLOCATE( bhatz, (0:NS), 0.0 )
        SALLOCATE( thatpx, (0:NS), 0.0 )
        SALLOCATE( thatpy, (0:NS), 0.0 )
        SALLOCATE( thatpz, (0:NS), 0.0 )
        SALLOCATE( nhatax, (0:NS), 0.0 )
        SALLOCATE( nhatay, (0:NS), 0.0 )
        SALLOCATE( nhataz, (0:NS), 0.0 )
        SALLOCATE( bhatax, (0:NS), 0.0 )
        SALLOCATE( bhatay, (0:NS), 0.0 )
        SALLOCATE( bhataz, (0:NS), 0.0 )
        
        absrp(0:NS-1) = sqrt( (coil(icoil)%xt(0:NS-1))**2 + (coil(icoil)%yt(0:NS-1))**2 + (coil(icoil)%zt(0:NS-1))**2 )
        xc = sum( coil(icoil)%xx(0:NS-1)*absrp(0:NS-1) ) / sum( absrp(0:NS-1) )
        yc = sum( coil(icoil)%yy(0:NS-1)*absrp(0:NS-1) ) / sum( absrp(0:NS-1) )
        zc = sum( coil(icoil)%zz(0:NS-1)*absrp(0:NS-1) ) / sum( absrp(0:NS-1) )
        
        deltax(0:NS) = coil(icoil)%xx(0:NS) - xc
        deltay(0:NS) = coil(icoil)%yy(0:NS) - yc
        deltaz(0:NS) = coil(icoil)%zz(0:NS) - zc
        
        coil(icoil)%nfbx(0:NS) = deltax(0:NS) - (deltax(0:NS)*coil(icoil)%xt(0:NS) + deltay(0:NS)*coil(icoil)%yt(0:NS) + deltaz(0:NS)*coil(icoil)%zt(0:NS)) &
             * coil(icoil)%xt(0:NS) / ( coil(icoil)%xt(0:NS)**2.0 + coil(icoil)%yt(0:NS)**2.0 + coil(icoil)%zt(0:NS)**2.0 )
        coil(icoil)%nfby(0:NS) = deltay(0:NS) - (deltax(0:NS)*coil(icoil)%xt(0:NS) + deltay(0:NS)*coil(icoil)%yt(0:NS) + deltaz(0:NS)*coil(icoil)%zt(0:NS)) &
             * coil(icoil)%yt(0:NS) / ( coil(icoil)%xt(0:NS)**2.0 + coil(icoil)%yt(0:NS)**2.0 + coil(icoil)%zt(0:NS)**2.0 )
        coil(icoil)%nfbz(0:NS) = deltaz(0:NS) - (deltax(0:NS)*coil(icoil)%xt(0:NS) + deltay(0:NS)*coil(icoil)%yt(0:NS) + deltaz(0:NS)*coil(icoil)%zt(0:NS)) &
             * coil(icoil)%zt(0:NS) / ( coil(icoil)%xt(0:NS)**2.0 + coil(icoil)%yt(0:NS)**2.0 + coil(icoil)%zt(0:NS)**2.0 )

        nx(0:NS) = coil(icoil)%nfbx(0:NS)
        ny(0:NS) = coil(icoil)%nfby(0:NS)
        nz(0:NS) = coil(icoil)%nfbz(0:NS)
        normn(0:NS) = sqrt( nx(0:NS)**2.0 + ny(0:NS)**2.0 + nz(0:NS)**2.0 )

        norm(0:NS) = sqrt( coil(icoil)%nfbx(0:NS)**2.0 + coil(icoil)%nfby(0:NS)**2.0 + coil(icoil)%nfbz(0:NS)**2.0 )
        coil(icoil)%nfbx(0:NS) = coil(icoil)%nfbx(0:NS) / norm(0:NS)
        coil(icoil)%nfby(0:NS) = coil(icoil)%nfby(0:NS) / norm(0:NS)
        coil(icoil)%nfbz(0:NS) = coil(icoil)%nfbz(0:NS) / norm(0:NS)

        nhatx(0:NS) = coil(icoil)%nfbx(0:NS)
        nhaty(0:NS) = coil(icoil)%nfby(0:NS)
        nhatz(0:NS) = coil(icoil)%nfbz(0:NS)

        coil(icoil)%bfbx(0:NS) = coil(icoil)%yt(0:NS)*coil(icoil)%nfbz(0:NS) - coil(icoil)%zt(0:NS)*coil(icoil)%nfby(0:NS)
        coil(icoil)%bfby(0:NS) = coil(icoil)%zt(0:NS)*coil(icoil)%nfbx(0:NS) - coil(icoil)%xt(0:NS)*coil(icoil)%nfbz(0:NS)
        coil(icoil)%bfbz(0:NS) = coil(icoil)%xt(0:NS)*coil(icoil)%nfby(0:NS) - coil(icoil)%yt(0:NS)*coil(icoil)%nfbx(0:NS)
        
        norm(0:NS) = sqrt( coil(icoil)%bfbx(0:NS)**2.0 + coil(icoil)%bfby(0:NS)**2.0 + coil(icoil)%bfbz(0:NS)**2.0 )
        coil(icoil)%bfbx(0:NS) = coil(icoil)%bfbx(0:NS) / norm(0:NS)
        coil(icoil)%bfby(0:NS) = coil(icoil)%bfby(0:NS) / norm(0:NS)
        coil(icoil)%bfbz(0:NS) = coil(icoil)%bfbz(0:NS) / norm(0:NS)
        
        bhatx(0:NS) = coil(icoil)%bfbx(0:NS)
        bhaty(0:NS) = coil(icoil)%bfby(0:NS)
        bhatz(0:NS) = coil(icoil)%bfbz(0:NS)

        normt(0:NS) = sqrt( coil(icoil)%xt(0:NS)**2.0 + coil(icoil)%yt(0:NS)**2.0 + coil(icoil)%zt(0:NS)**2.0 )
        
        thatx(0:NS) = coil(icoil)%xt(0:NS) / normt(0:NS)
        thaty(0:NS) = coil(icoil)%yt(0:NS) / normt(0:NS)
        thatz(0:NS) = coil(icoil)%zt(0:NS) / normt(0:NS)
        
        thatpx(0:NS) = coil(icoil)%xa(0:NS)/normt(0:NS) - ( coil(icoil)%xt(0:NS)*coil(icoil)%xa(0:NS) + coil(icoil)%yt(0:NS)*coil(icoil)%ya(0:NS) &
                     + coil(icoil)%zt(0:NS)*coil(icoil)%za(0:NS) )*coil(icoil)%xt(0:NS)/(normt(0:NS)**3.0)
        
        thatpy(0:NS) = coil(icoil)%ya(0:NS)/normt(0:NS) - ( coil(icoil)%xt(0:NS)*coil(icoil)%xa(0:NS) + coil(icoil)%yt(0:NS)*coil(icoil)%ya(0:NS) &
                     + coil(icoil)%zt(0:NS)*coil(icoil)%za(0:NS) )*coil(icoil)%yt(0:NS)/(normt(0:NS)**3.0)
        
        thatpz(0:NS) = coil(icoil)%za(0:NS)/normt(0:NS) - ( coil(icoil)%xt(0:NS)*coil(icoil)%xa(0:NS) + coil(icoil)%yt(0:NS)*coil(icoil)%ya(0:NS) &
                     + coil(icoil)%zt(0:NS)*coil(icoil)%za(0:NS) )*coil(icoil)%zt(0:NS)/(normt(0:NS)**3.0)
        
        npx(0:NS) = coil(icoil)%xt(0:NS) - ( coil(icoil)%xt(0:NS)*thatx(0:NS) + coil(icoil)%yt(0:NS)*thaty(0:NS) + &
                    coil(icoil)%zt(0:NS)*thatz(0:NS) )*thatx(0:NS) - ( deltax(0:NS)*thatpx(0:NS) + deltay(0:NS)*thatpy(0:NS) + &
                    deltaz(0:NS)*thatpz(0:NS) )*thatx(0:NS) - ( deltax(0:NS)*thatx(0:NS) + deltay(0:NS)*thaty(0:NS) + deltaz(0:NS)*thatz(0:NS) )*thatpx(0:NS)
        npy(0:NS) = coil(icoil)%yt(0:NS) - ( coil(icoil)%xt(0:NS)*thatx(0:NS) + coil(icoil)%yt(0:NS)*thaty(0:NS) + &
                    coil(icoil)%zt(0:NS)*thatz(0:NS) )*thaty(0:NS) - ( deltax(0:NS)*thatpx(0:NS) + deltay(0:NS)*thatpy(0:NS) + &
                    deltaz(0:NS)*thatpz(0:NS) )*thaty(0:NS) - ( deltax(0:NS)*thatx(0:NS) + deltay(0:NS)*thaty(0:NS) + deltaz(0:NS)*thatz(0:NS) )*thatpy(0:NS)
        npz(0:NS) = coil(icoil)%zt(0:NS) - ( coil(icoil)%xt(0:NS)*thatx(0:NS) + coil(icoil)%yt(0:NS)*thaty(0:NS) + &
                    coil(icoil)%zt(0:NS)*thatz(0:NS) )*thatz(0:NS) - ( deltax(0:NS)*thatpx(0:NS) + deltay(0:NS)*thatpy(0:NS) + &
                    deltaz(0:NS)*thatpz(0:NS) )*thatz(0:NS) - ( deltax(0:NS)*thatx(0:NS) + deltay(0:NS)*thaty(0:NS) + deltaz(0:NS)*thatz(0:NS) )*thatpz(0:NS)
        
        nhatpx(0:NS) = npx(0:NS)/normn(0:NS) - ( nx(0:NS)*npx(0:NS) + ny(0:NS)*npy(0:NS) + nz(0:NS)*npz(0:NS) )*nx(0:NS)/(normn(0:NS)**3.0)
        nhatpy(0:NS) = npy(0:NS)/normn(0:NS) - ( nx(0:NS)*npx(0:NS) + ny(0:NS)*npy(0:NS) + nz(0:NS)*npz(0:NS) )*ny(0:NS)/(normn(0:NS)**3.0)
        nhatpz(0:NS) = npz(0:NS)/normn(0:NS) - ( nx(0:NS)*npx(0:NS) + ny(0:NS)*npy(0:NS) + nz(0:NS)*npz(0:NS) )*nz(0:NS)/(normn(0:NS)**3.0)
        
        C = sum( nhatpx(0:NS-1)*bhatx(0:NS-1) + nhatpy(0:NS-1)*bhaty(0:NS-1) + nhatpz(0:NS-1)*bhatz(0:NS-1) ) / sum( normt(0:NS-1) )
        
        SALLOCATE( alphac, (0:Nalpha), 0.0 )
        SALLOCATE( alphas, (0:Nalpha), 0.0 )
        
        do i = 0, NS
           zeta(i) = real(i)*pi2/real(NS)
        enddo
        do i = 1, Nalpha
           alphac(i) = -2.0*sum( sin(real(i)*zeta(0:NS-1))*( nhatpx(0:NS-1)*bhatx(0:NS-1) + nhatpy(0:NS-1)*bhaty(0:NS-1) + nhatpz(0:NS-1)*bhatz(0:NS-1) - &
                       C*normt(0:NS-1) ) )/(real(NS)*real(i))
           alphas(i) =  2.0*sum( cos(real(i)*zeta(0:NS-1))*( nhatpx(0:NS-1)*bhatx(0:NS-1) + nhatpy(0:NS-1)*bhaty(0:NS-1) + nhatpz(0:NS-1)*bhatz(0:NS-1) - &
                       C*normt(0:NS-1) ) )/(real(NS)*real(i))
        enddo
        
        SALLOCATE( alpha, (0:NS), 0.0 )
        do i = 1, Nalpha
           alpha(0:NS) = alpha(0:NS) + alphac(i)*cos(real(i)*zeta) + alphas(i)*sin(real(i)*zeta)
        enddo
       
        nhatax(0:NS) = -1.0*sin(alpha(0:NS))*bhatx(0:NS) + cos(alpha(0:NS))*nhatx(0:NS)
        nhatay(0:NS) = -1.0*sin(alpha(0:NS))*bhaty(0:NS) + cos(alpha(0:NS))*nhaty(0:NS)
        nhataz(0:NS) = -1.0*sin(alpha(0:NS))*bhatz(0:NS) + cos(alpha(0:NS))*nhatz(0:NS)

        bhatax(0:NS) =      sin(alpha(0:NS))*nhatx(0:NS) + cos(alpha(0:NS))*bhatx(0:NS)
        bhatay(0:NS) =      sin(alpha(0:NS))*nhaty(0:NS) + cos(alpha(0:NS))*bhaty(0:NS)
        bhataz(0:NS) =      sin(alpha(0:NS))*nhatz(0:NS) + cos(alpha(0:NS))*bhatz(0:NS)
        
        coil(icoil)%nfbx(0:NS) = nhatax(0:NS)
        coil(icoil)%nfby(0:NS) = nhatay(0:NS)
        coil(icoil)%nfbz(0:NS) = nhataz(0:NS)
        
        coil(icoil)%bfbx(0:NS) = bhatax(0:NS)
        coil(icoil)%bfby(0:NS) = bhatay(0:NS)
        coil(icoil)%bfbz(0:NS) = bhataz(0:NS)

        DALLOCATE(absrp)
        DALLOCATE(deltax)
        DALLOCATE(deltay)
        DALLOCATE(deltaz)
        DALLOCATE(norm)
        DALLOCATE(nx)
        DALLOCATE(ny)
        DALLOCATE(nz)
        DALLOCATE(npx)
        DALLOCATE(npy)
        DALLOCATE(npz)
        DALLOCATE(nhatpx)
        DALLOCATE(nhatpy)
        DALLOCATE(nhatpz)
        DALLOCATE(normt)
        DALLOCATE(normn)
        DALLOCATE(zeta)
        DALLOCATE(thatx)
        DALLOCATE(thaty)
        DALLOCATE(thatz)
        DALLOCATE(nhatx)
        DALLOCATE(nhaty)
        DALLOCATE(nhatz)
        DALLOCATE(bhatx)
        DALLOCATE(bhaty)
        DALLOCATE(bhatz)
        DALLOCATE(thatpx)
        DALLOCATE(thatpy)
        DALLOCATE(thatpz)
        DALLOCATE(nhatax)
        DALLOCATE(nhatay)
        DALLOCATE(nhataz)
        DALLOCATE(bhatax)
        DALLOCATE(bhatay)
        DALLOCATE(bhataz)
        DALLOCATE(alphac)
        DALLOCATE(alphas)
        DALLOCATE(alpha)
     enddo ! Loop over coils
  endif

  return ! Should below be deleted?

  !--------------------------------calculate coil importance------------------------------------  
  if (.not. allocated(coil_importance)) then
     SALLOCATE( coil_importance, (1:Ncoils), zero )
  endif

  if (weight_bnorm > sqrtmachprec .or. weight_bharm > sqrtmachprec) then  ! make sure data_allocated
     do icoil = 1, Ncoils
        call importance(icoil)
     enddo
  
     if(myid .eq. 0) write(ounit, '(8X": The most and least important coils are :  " & 
          ES12.5" at coil" I4 " ; " ES12.5" at coil "I4)')      &
      maxval(coil_importance), maxloc(coil_importance), &
      minval(coil_importance), minloc(coil_importance)

  endif
  !--------------------------------------------------------------------------------------------- 

  if (myid == 0 .and. IsQuiet < 0) write(ounit, *) "--------------------------------------------"
 
  return

END SUBROUTINE diagnos

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine avgcurvature(icoil)


  use globals, only: dp, zero, pi2, ncpu, astat, ierr, myid, ounit, coil, NFcoil, Nseg, Ncoils,Splines,coil_type_spline
  use mpi
  implicit none

  INTEGER, INTENT(in) :: icoil

  REAL,allocatable    :: davgcurv(:)
  call length(0)
  SALLOCATE(davgcurv, (0:coil(icoil)%NS-1), zero)

  davgcurv = sqrt( (coil(icoil)%za*coil(icoil)%yt-coil(icoil)%zt*coil(icoil)%ya)**2  &
             + (coil(icoil)%xa*coil(icoil)%zt-coil(icoil)%xt*coil(icoil)%za)**2  & 
             + (coil(icoil)%ya*coil(icoil)%xt-coil(icoil)%yt*coil(icoil)%xa)**2 )& 
             / ((coil(icoil)%xt)**2+(coil(icoil)%yt)**2+(coil(icoil)%zt)**2)**(1.5)
  davgcurv = davgcurv*sqrt(coil(icoil)%xt**2+coil(icoil)%yt**2+coil(icoil)%zt**2)

  if (coil(icoil)%type == coil_type_spline) coil(icoil)%avgcurv = 1.0*sum(davgcurv)/size(davgcurv)
  if (coil(icoil)%type == 1)coil(icoil)%avgcurv = pi2*sum(davgcurv)/size(davgcurv)

  coil(icoil)%avgcurv = coil(icoil)%avgcurv/coil(icoil)%L

  return
end subroutine avgcurvature

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine mindist(array_A, dim_A, array_B, dim_B, minimum)
 
  use globals, only: dp
  implicit none

  INTEGER, INTENT(IN ) :: dim_A, dim_B
  REAL   , INTENT(IN ) :: array_A(1:3,1:dim_A), array_B(1:3,1:dim_B)
  REAL   , INTENT(OUT) :: minimum

  INTEGER :: ipoint, jpoint, itmp, jtmp
  REAL, parameter :: infmax = 1.0E6
  REAL    :: distance

  minimum = infmax
  do ipoint = 1, dim_A
     do jpoint = 1, dim_B

        distance = ( array_A(1, ipoint) - array_B(1, jpoint) )**2 &
                  +( array_A(2, ipoint) - array_B(2, jpoint) )**2 &
                  +( array_A(3, ipoint) - array_B(3, jpoint) )**2
        if (distance .le. minimum) minimum = distance
        
     enddo
  enddo

  minimum = sqrt(minimum)
  
  return

end subroutine mindist

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine importance(icoil)
  use globals, only: dp,  zero, pi2, ncpu, astat, ierr, myid, ounit, coil, NFcoil, Nseg, Ncoils, &
                     surf, Nteta, Nzeta, bsconstant, coil_importance, plasma, MPI_COMM_FOCUS
  use mpi
  implicit none

  INTEGER, INTENT(in) :: icoil  

  INTEGER                               :: iteta, jzeta, NumGrid, isurf
  REAL                                  :: dBx, dBy, dBz
  REAL, dimension(0:Nteta-1, 0:Nzeta-1) :: tbx, tby, tbz        ! summed Bx, By and Bz

  !--------------------------initialize and allocate arrays------------------------------------- 

  isurf = plasma
  NumGrid = Nteta*Nzeta
  tbx = zero; tby = zero; tbz = zero        !already allocted; reset to zero;

  do jzeta = 0, Nzeta - 1
     do iteta = 0, Nteta - 1

        if( myid.ne.modulo(jzeta*Nteta+iteta,ncpu) ) cycle ! parallelization loop;
        call bfield0(icoil, surf(isurf)%xx(iteta, jzeta), surf(isurf)%yy(iteta, jzeta), &
                   & surf(isurf)%zz(iteta, jzeta), tbx(iteta, jzeta), tby(iteta, jzeta), tbz(iteta, jzeta))

     enddo ! end do iteta
  enddo ! end do jzeta

  call MPI_BARRIER( MPI_COMM_FOCUS, ierr )     
  call MPI_ALLREDUCE( MPI_IN_PLACE, tbx, NumGrid, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FOCUS, ierr )
  call MPI_ALLREDUCE( MPI_IN_PLACE, tby, NumGrid, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FOCUS, ierr )
  call MPI_ALLREDUCE( MPI_IN_PLACE, tbz, NumGrid, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FOCUS, ierr )

  coil_importance(icoil) = sum( (tbx*surf(isurf)%Bx + tby*surf(isurf)%By + tbz*surf(isurf)%Bz) / &
                                (surf(isurf)%Bx**2 + surf(isurf)%By**2 + surf(isurf)%Bz**2) ) / NumGrid

  return

end subroutine importance
