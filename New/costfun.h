!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine costfun(ideriv)
  use globals, only: zero, sqrtmachprec, myid, ounit, IsQuiet, iter, &
       Ncoils, deriv, Ndof, xdof, &
       totalenergy, t1E, t2E, &
       bnorm      , t1B, t2B, weight_bnorm,  &
       tflux      , t1F, t2F, weight_tflux, target_tflux, isign, &
       ttlen      , t1L, t2L, weight_ttlen, &
       specw      , t1S, t2S, weight_specw, &
       ccsep      , t1C, t2C, weight_ccsep

  implicit none
  include "mpif.h"

  INTEGER, INTENT(in) :: ideriv

  INTEGER      :: astat, ierr

  REAL         :: start, finish
  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  iter = iter + 1
  ! if(IsQuiet .eq. -1 .and. myid .eq. 0) write(ounit,'("Costfun :"12X": begin the "I6"th run.")') iter
  if(IsQuiet <= -1) then
   call bnormal(0)
!!$
!!$   if ( target_tflux .eq. 0.0 ) then
!!$    call torflux(0)
!!$    target_tflux = isign*sqrt(2.0*tflux)             
!!$    if(myid .eq. 0) write(ounit,'("Costfun :"11X" : Reset target toroidal flux to"ES23.15)') target_tflux
!!$   endif
!!$
!!$   call torflux(0)
!!$
!!$   if (lc .eq. 1) then
!!$    call tlength(0)
!!$   else
!!$    call tlengthExp(0)
!!$   endif
!!$   
!!$   call specwid(0)
!!$   call coilsep(0)
  endif
  
  totalenergy = zero
  
  !if ( ideriv .ge. 1 .and. .not. allocated(t1E) ) ALLOCATE(t1E(1:Ndof))
  !if ( ideriv .eq. 2 .and. .not. allocated(t2E) ) ALLOCATE(t2E(1:Ncoils, 0:Cdof, 1:Ncoils, 0:Cdof))

  if     ( ideriv .eq. 1 ) then
   t1E = zero
  elseif ( ideriv .eq. 2 ) then
   t1E = zero; t2E = zero
  endif

  !call unpacking(xdof)

  if (weight_bnorm .gt. sqrtmachprec) then

   call bnormal(ideriv)
   totalenergy = totalenergy + weight_bnorm * bnorm
   if     ( ideriv .eq. 1 ) then
    t1E = t1E +  weight_bnorm * t1B
   elseif ( ideriv .eq. 2 ) then
    t1E = t1E +  weight_bnorm * t1B
    t2E = t2E +  weight_bnorm * t2B
   endif

  endif

!!$  ! if(myid .eq. 0) write(ounit,'("calling bnormal used",f10.5,"seconds.")') finish-start
!!$
!!$  if (weight_tflux .gt. sqrtmachprec) then
!!$
!!$   if ( target_tflux .eq. 0.0 ) then
!!$    call torflux(0)
!!$    target_tflux = isign*sqrt(2.0*tflux)             
!!$    if(myid .eq. 0) write(ounit,'("Costfun :"11X" : Reset target toroidal flux to"ES23.15)') target_tflux
!!$   endif
!!$
!!$   call torflux(ideriv)
!!$   totalenergy = totalenergy + weight_tflux * tflux / target_tflux**2 ! normalization
!!$   if     ( ideriv .eq. 1 ) then
!!$    t1E = t1E +  weight_tflux * t1F / target_tflux**2
!!$   elseif ( ideriv .eq. 2 ) then
!!$    t1E = t1E +  weight_tflux * t1F / target_tflux**2
!!$    t2E = t2E +  weight_tflux * t2F / target_tflux**2
!!$   endif
!!$
!!$  endif
!!$
!!$  ! if(myid .eq. 0) write(ounit,'("calling torflux used",f10.5,"seconds.")') start-finish
!!$
!!$  if (weight_ttlen .gt. sqrtmachprec) then
!!$
!!$   if( lc .eq. 1 ) then 
!!$    call tlength(ideriv)
!!$   elseif (lc .eq. 2) then
!!$    call tlengthExp(ideriv)
!!$   else
!!$    FATAL( dnergy, .true., Conflicts between lc and weight_ttlen )
!!$   !call MPI_ABORT( MPI_COMM_WORLD, 1, ierr )
!!$   !stop "COSTFUN: Conflicts between lc and weight_ttlen"
!!$   endif
!!$   
!!$   totalenergy = totalenergy + weight_ttlen * ttlen
!!$   if     ( ideriv .eq. 1 ) then
!!$    t1E = t1E +  weight_ttlen * t1L
!!$   elseif ( ideriv .eq. 2 ) then
!!$    t1E = t1E +  weight_ttlen * t1L
!!$    t2E = t2E +  weight_ttlen * t2L
!!$   endif
!!$
!!$  endif
!!$
!!$  ! if(myid .eq. 0) write(ounit,'("calling tlength used",f10.5,"seconds.")') finish-start
!!$
!!$  if (weight_eqarc .ge. sqrtmachprec) then
!!$
!!$  ! call equarcl(ideriv)
!!$  ! call specwid_df(ideriv)
!!$   if ( Loptimize .eq. 3) then
!!$    call langrange(ideriv)
!!$   else
!!$    !call specwid(ideriv)
!!$    call specwid_df(ideriv)
!!$   endif
!!$   
!!$   totalenergy = totalenergy + weight_eqarc * eqarc
!!$   if     ( ideriv .eq. 1 ) then
!!$    t1E = t1E + weight_eqarc * t1A
!!$   elseif ( ideriv .eq. 2 ) then
!!$    t1E = t1E + weight_eqarc * t1A
!!$    t2E = t2E + weight_eqarc * t2A
!!$   endif
!!$
!!$  endif
!!$
!!$  if (weight_ccsep .ge. sqrtmachprec) then
!!$
!!$   call coilsep(ideriv)
!!$   totalenergy = totalenergy + weight_ccsep * ccsep
!!$   if     ( ideriv .eq. 1 ) then
!!$    t1E = t1E + weight_ccsep * t1C
!!$   elseif ( ideriv .eq. 2 ) then
!!$    t1E = t1E + weight_ccsep * t1C
!!$    t2E = t2E + weight_ccsep * t2C
!!$   endif
!!$
!!$  endif
!!$
!!$  if(allocated(deriv)) then
!!$     deriv = zero
!!$     do ii = 1, Ndof
!!$        call DoFconvert(ii, icl, inf)
!!$        if(allocated(t1E)) deriv(ii, 0) = t1E(icl, inf)
!!$        if(allocated(t1B)) deriv(ii, 1) = t1B(icl, inf)
!!$        if(allocated(t1F)) deriv(ii, 2) = t1F(icl, inf)
!!$        if(allocated(t1L)) deriv(ii, 3) = t1L(icl, inf)
!!$        if(allocated(t1A)) deriv(ii, 4) = t1A(icl, inf)
!!$        if(allocated(t1C)) deriv(ii, 5) = t1C(icl, inf)
!!$     enddo
!!$  endif


  FATAL( costfun, iter > 100000, too many iterations for one single call)
 !if(iter .ge. 1E5) call MPI_ABORT( MPI_COMM_WORLD, 10, ierr )

  return
end subroutine costfun

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine AllocData(part)
!------------------------------------------------------------------------------------------------------ 
! DATE:  04/05/2017
! Allocate data before using them, especially for those used several times;
! part can be : 'dof', 'costfun0', 'costfun1'
!------------------------------------------------------------------------------------------------------
  use globals, only :  zero, sqrtmachprec, myid, ounit, ncpu, case_optimizer, &
       &               coil, DoF, FouCoil, surf, Ncoils, Nteta, Nzeta, Ndof, xdof, Inorm, Gnorm, dofnorm, &
       &               Tdof, evolution, deriv, coilspace, SD_Nout,  &
       &               totalenergy, t1E, t2E, evolution, coilspace, deriv, &
       &               bnorm      , t1B, t2B, weight_bnorm, bn, &
       &               tflux      , t1F, t2F, weight_tflux, target_tflux, isign, &
       &               ttlen      , t1L, t2L, weight_ttlen, &
       &               specw      , t1S, t2S, weight_specw, &
       &               ccsep      , t1C, t2C, weight_ccsep
  implicit none
  include "mpif.h"

  CHARACTER(*) , intent(in) :: part

  INTEGER                   :: ierr, astat, icoil, idof, NS, ND, NF

  !-------------------------------------------------------------------------------------------
  
  select case (part)

  !--------------------------------------------------------------------------------------------- 
  case ('dof')

     Ndof = 0; Tdof = 0

     do icoil = 1, Ncoils
        select case (coil(icoil)%itype)
        
        case(1)
           ! get number of DoF for each coil and allocate arrays;
           NF = FouCoil(icoil)%NF
           ND = (6*NF + 3) ! total variables for geometry
           DoF(icoil)%ND = coil(icoil)%Lc * ND !# of DoF for icoil;
           SALLOCATE(DoF(icoil)%xdof, (1:DoF(icoil)%ND), zero)
           SALLOCATE(DoF(icoil)%xof , (1:coil(icoil)%NS, 1:ND), zero)
           SALLOCATE(DoF(icoil)%yof , (1:coil(icoil)%NS, 1:ND), zero)
           SALLOCATE(DoF(icoil)%zof , (1:coil(icoil)%NS, 1:ND), zero)
        case default
           FATAL(AllocData, .true., not supported coil types)
        end select

        Ndof = Ndof + coil(icoil)%Ic + DoF(icoil)%ND
        Tdof = Tdof + 1              + ND

     enddo

     if(Ndof == 0) then ! no DOF;
        case_optimizer = 0
        if(myid==0) write(ounit, *) "AllocData : No free variables; no optimization performed."
     endif

     SALLOCATE(    xdof, (1:Ndof), zero ) ! dof vector;
     SALLOCATE( dofnorm, (1:Ndof), zero ) ! dof normalized value vector;
     SALLOCATE( evolution, (0:SD_Nout, 0:7), zero ) !evolution array;
     SALLOCATE( coilspace, (0:SD_Nout, 1:Tdof), zero ) ! all the coil parameters;
     
     idof = 0
     do icoil = 1, Ncoils

        if(coil(icoil)%Ic /= 0) then
           dofnorm(idof+1) = Inorm
           idof = idof + 1
        endif

        ND = DoF(icoil)%ND
        if(coil(icoil)%Lc /= 0) then
           dofnorm(idof+1:idof+ND) = Gnorm
           idof = idof + ND
        endif

     enddo !end do icoil;
     FATAL( AllocData , idof .ne. Ndof, counting error in unpacking )
  !--------------------------------------------------------------------------------------------- 
  case ('costfun0')

     if (weight_bnorm > sqrtmachprec) then
        SALLOCATE(         bn, (0:Nteta,0:Nzeta), zero ) !Bn from coils;        

        SALLOCATE( surf(1)%bn, (0:Nteta,0:Nzeta), zero ) !total Bn;
        SALLOCATE( surf(1)%Bx, (0:Nteta,0:Nzeta), zero ) !Bx on the surface;
        SALLOCATE( surf(1)%By, (0:Nteta,0:Nzeta), zero ) !By on the surface;
        SALLOCATE( surf(1)%Bz, (0:Nteta,0:Nzeta), zero ) !Bz on the surface;

        do icoil = 1, Ncoils
           ND = DoF(icoil)%ND
           SALLOCATE( coil(icoil)%Bx, (0:ND, 0:ND), zero )
           SALLOCATE( coil(icoil)%By, (0:ND, 0:ND), zero )
           SALLOCATE( coil(icoil)%Bz, (0:ND, 0:ND), zero )
        enddo
     endif

  !--------------------------------------------------------------------------------------------- 
  case ('costfun1')
     
     FATAL( AllocData, Ndof < 1, INVALID Ndof value )
     SALLOCATE( t1E, (1:Ndof), zero )
     SALLOCATE( deriv, (1:Ndof, 0:5), zero )

     if (weight_bnorm .gt. sqrtmachprec) then
        SALLOCATE(         bn, (0:Nteta,0:Nzeta), zero ) !Bn from coils;        

        SALLOCATE( surf(1)%bn, (0:Nteta,0:Nzeta), zero ) !total Bn;
        SALLOCATE( surf(1)%Bx, (0:Nteta,0:Nzeta), zero ) !Bx on the surface;
        SALLOCATE( surf(1)%By, (0:Nteta,0:Nzeta), zero ) !By on the surface;
        SALLOCATE( surf(1)%Bz, (0:Nteta,0:Nzeta), zero ) !Bz on the surface;

        do icoil = 1, Ncoils
           ND = DoF(icoil)%ND
           SALLOCATE( coil(icoil)%Bx, (0:ND, 0:ND), zero )
           SALLOCATE( coil(icoil)%By, (0:ND, 0:ND), zero )
           SALLOCATE( coil(icoil)%Bz, (0:ND, 0:ND), zero )
        enddo

        SALLOCATE( t1B, (1:Ndof), zero )
     endif

  !--------------------------------------------------------------------------------------------- 
  case default
     FATAL( AllocData, .true., unsupported allocate type )

  end select


  return
end subroutine AllocData

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
