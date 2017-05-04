
!title (restart) ! Writes output and restart files.

!latex \briefly{Writes output and restart files.}

!latex \calledby{\link{knotopt} and \link{descent} }
!l tex \calls{\link{iccoil},\link{knotxx}}

!latex \tableofcontents

!latex \subsection{h5 format}
!latex \bi
!latex \item[1.] All the restart information etc. is written to file.
!latex \ei

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine restart( irestart )

  use globals

  use hdf5

  implicit none

  include "mpif.h"

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  INTEGER, intent(in):: irestart

  INTEGER            :: astat, ierr, ii, jj, ifail, llmodnp, idof, mm, NF

  REAL               :: tt, aa, teta, zeta, ax(1:3), at(1:3), az(1:3), xx(1:3), xt(1:3), xz(1:3)
  REAL               :: x(0:1), y(0:1), z(0:1)
  ! REAL, allocatable  :: perA(:)

  ! the following are used by the macros HWRITEXX below; do not alter/remove;
  INTEGER            :: hdfier, rank
  integer(hid_t)     :: file_id, space_id, dset_id                 ! warning: the string "INTEGER" is a macro;
  integer(hsize_t)   :: onedims(1:1), twodims(1:2), threedims(1:3) ! warning: the string "INTEGER" is a macro;
  INTEGER            :: imn
  REAL               :: tvolume

  CHARACTER          :: suffix*3, srestart*6


  if(itau > SD_Nout) itau = SD_Nout ! what's this for?

  if( myid.ne.0 ) return

  !SALLOCATE(perA, (1:Tdof), zero)

  !store the derivatives;
  if (allocated(deriv)) then
     deriv = zero
     if(allocated(t1E)) deriv(1:Ndof,0) = t1E(1:Ndof)
     if(allocated(t1B)) deriv(1:Ndof,1) = t1B(1:Ndof)
     if(allocated(t1F)) deriv(1:Ndof,2) = t1F(1:Ndof)
     if(allocated(t1L)) deriv(1:Ndof,3) = t1L(1:Ndof)
     if(allocated(t1S)) deriv(1:Ndof,4) = t1S(1:Ndof)
     if(allocated(t1C)) deriv(1:Ndof,5) = t1C(1:Ndof)
  endif
  !calculate the new Bn
  if (allocated(surf(1)%bn)) call BnFTran

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  call h5open_f( hdfier ) ! initialize Fortran interface;
  FATAL( restart, hdfier.ne.0, error calling h5open_f )

  call h5fcreate_f( "focus_"//trim(ext)//".h5", H5F_ACC_TRUNC_F, file_id, hdfier ) ! create new file;
  FATAL( restart, hdfier.ne.0, error calling h5fcreate_f )

  !INPUT namelist;
  HWRITEIV( 1                ,   IsQuiet       ,   IsQuiet                       )
  HWRITEIV( 1                ,   IsSymmetric   ,   IsSymmetric                   )
  HWRITEIV( 1                ,   case_surface  ,   case_surface                  )
  HWRITERV( 1                ,   knotsurf      ,   knotsurf                      )
  HWRITEIV( 1                ,   Nteta         ,   Nteta                         )
  HWRITEIV( 1                ,   Nzeta         ,   Nzeta                         )
  HWRITEIV( 1                ,   case_init     ,   case_init                     )
  HWRITEIV( 1                ,   case_coils    ,   case_coils                    )
  HWRITEIV( 1                ,   Ncoils        ,   Ncoils                        )
  HWRITERV( 1                ,   init_current  ,   init_current                  )
  HWRITERV( 1                ,   init_radius   ,   init_radius                   )
  HWRITEIV( 1                ,   IsVaryCurrent ,   IsVaryCurrent                 )
  HWRITEIV( 1                ,   IsVaryGeometry,   IsvaryGeometry                )
  HWRITERV( 1                ,   target_length ,   target_length                 )
  HWRITEIV( 1                ,   NFcoil        ,   NFcoil                        )
  HWRITEIV( 1                ,   Nseg          ,   Nseg                          )
  HWRITEIV( 1                ,   case_optimizer,   case_optimizer                )
  HWRITEIV( 1                ,   IsNormalize   ,   IsNormalize                   )
  HWRITEIV( 1                ,   ISNormWeight  ,   IsNormWeight                  )
  HWRITERV( 1                ,   weight_bnorm  ,   weight_bnorm                  )
  HWRITERV( 1                ,   weight_tflux  ,   weight_tflux                  )
  HWRITERV( 1                ,   target_tflux  ,   target_tflux                  )
  HWRITERV( 1                ,   weight_ttlen  ,   weight_ttlen                  )
  HWRITERV( 1                ,   weight_specw  ,   weight_specw                  )
  HWRITERV( 1                ,   weight_ccsep  ,   weight_ccsep                  )
  HWRITERV( 1                ,   SD_tausta     ,   SD_tausta                     )
  HWRITERV( 1                ,   SD_tauend     ,   SD_tauend                     )
  HWRITEIV( 1                ,   SD_Nout       ,   SD_Nout                       )
  HWRITERV( 1                ,   SD_tautol     ,   SD_tautol                     )
  HWRITEIV( 1                ,   SD_SaveFreq   ,   SD_SaveFreq                   )
  HWRITERV( 1                ,   NT_xtol       ,   NT_xtol                       )
  HWRITERV( 1                ,   NT_eta        ,   NT_eta                        )
  HWRITERV( 1                ,   NT_stepmx     ,   NT_stepmx                     )
  HWRITEIV( 1                ,   case_postproc ,   case_postproc                 )
  HWRITERV( 1                ,   PP_odetol     ,   PP_odetol                     )
  HWRITEIV( 1                ,   PP_Ppts       ,   PP_Ppts                       )
  HWRITEIV( 1                ,   PP_Ptrj       ,   PP_Ptrj                       )
  HWRITERV( 1                ,   PP_phi        ,   PP_phi                        )
  HWRITERV( 1                ,   PP_Rmin       ,   PP_Rmin                       )
  HWRITERV( 1                ,   PP_Rmax       ,   PP_Rmax                       )
  HWRITERV( 1                ,   PP_Zmin       ,   PP_Zmin                       )
  HWRITERV( 1                ,   PP_Zmax       ,   PP_Zmax                       )
  HWRITERV( 1                ,   PP_bstol      ,   PP_bstol                      )
  HWRITEIV( 1                ,   PP_bsnlimit   ,   PP_bsnlimit                   )

  HWRITERA( Nteta,Nzeta      ,   xsurf         ,   surf(1)%xx(0:Nteta-1,0:Nzeta-1) )
  HWRITERA( Nteta,Nzeta      ,   ysurf         ,   surf(1)%yy(0:Nteta-1,0:Nzeta-1) )
  HWRITERA( Nteta,Nzeta      ,   zsurf         ,   surf(1)%zz(0:Nteta-1,0:Nzeta-1) )
  HWRITERA( Nteta,Nzeta      ,   nx            ,   surf(1)%nx(0:Nteta-1,0:Nzeta-1) )
  HWRITERA( Nteta,Nzeta      ,   ny            ,   surf(1)%ny(0:Nteta-1,0:Nzeta-1) )
  HWRITERA( Nteta,Nzeta      ,   nz            ,   surf(1)%nz(0:Nteta-1,0:Nzeta-1) )

  if (allocated(bn)) then
     HWRITERA( Nteta,Nzeta      ,   tgtBn         ,   surf(1)%tn(0:Nteta-1,0:Nzeta-1) )
     HWRITERA( Nteta,Nzeta      ,   curBn         ,           bn(0:Nteta-1,0:Nzeta-1) )
     HWRITERA( Nteta,Nzeta      ,   Bx            ,   surf(1)%Bx(0:Nteta-1,0:Nzeta-1) )
     HWRITERA( Nteta,Nzeta      ,   By            ,   surf(1)%By(0:Nteta-1,0:Nzeta-1) )
     HWRITERA( Nteta,Nzeta      ,   Bz            ,   surf(1)%Bz(0:Nteta-1,0:Nzeta-1) )
  endif

  HWRITEIV( 1                ,   itau          ,   itau                          )
  HWRITERA( itau+1, 8        ,   evolution     ,   evolution(0:itau, 0:7)        )
  HWRITERA( itau+1, Tdof     ,   coilspace     ,   coilspace(0:itau, 1:Tdof)     )
  if (allocated(deriv)) then
     HWRITERA( Ndof, 6       ,   deriv         ,   deriv(1:Ndof, 0:5)            )
  endif

  if (allocated(Cur_Bns)) then
     HWRITEIV( NBnf          ,   Bnin          ,   Bnin                          )
     HWRITEIV( NBnf          ,   Bnim          ,   Bnim                          )
     HWRITERV( NBnf          ,   Cur_Bnc       ,   Cur_Bnc                       )
     HWRITERV( NBnf          ,   Cur_Bns       ,   Cur_Bns                       )
  endif

  HWRITERV( 1                ,  time_initialize,   time_initialize               )
  HWRITERV( 1                ,  time_optimize  ,   time_optimize                 )
  HWRITERV( 1                ,  time_postproc  ,   time_postproc                 )

  call h5fclose_f( file_id, hdfier ) ! terminate access;
  FATAL( restart, hdfier.ne.0, error calling h5fclose_f )

  call h5close_f( hdfier ) ! close Fortran interface;
  FATAL( restart, hdfier.ne.0, error calling h5close_f )


  !--------------------------write focus coil file-----------------------------------------
  open( wunit, file=trim(ext)//".focus", status="unknown" )
  write(wunit, *), "# Total number of coils"
  write(wunit, '(I6)'), Ncoils

  do icoil = 1, Ncoils

     write(wunit, *), "#--------------------------------------------" 
     write(wunit, *), "# coil_type     coil_name"
     write(wunit,'(I3,4X, A10)'), coil(icoil)%itype, coil(icoil)%name
     write(wunit, *), "# Nseg   current  I_flag  Length L_flag target_length"
     write(wunit,'(I4, ES23.15, I3, ES23.15, I3, ES23.15)'), &
          coil(icoil)%NS, coil(icoil)%I, coil(icoil)%Ic, coil(icoil)%L, coil(icoil)%Lc, coil(icoil)%Lo
     select case (coil(icoil)%itype)
     case (1)
        NF = FouCoil(icoil)%NF ! shorthand;
        write(wunit, *) "# NFcoil"
        write(wunit, '(I3)') NF
        write(wunit, *) "# Fourier harmonics for coils ( xc; xs; yc; ys; zc; zs) "
        write(wunit, 1000) FouCoil(icoil)%xc(0:NF)
        write(wunit, 1000) FouCoil(icoil)%xs(0:NF)
        write(wunit, 1000) FouCoil(icoil)%yc(0:NF)
        write(wunit, 1000) FouCoil(icoil)%ys(0:NF)
        write(wunit, 1000) FouCoil(icoil)%zc(0:NF)
        write(wunit, 1000) FouCoil(icoil)%zs(0:NF)
     case default
        FATAL(restart, .true., not supported coil types)
     end select
  enddo
  close(wunit)
1000 format(9999ES23.15)

  !--------------------------write coils.ext file-----------------------------------------------  
  if( irestart == 1 ) then

     open(funit,file=trim(ext)//".coils", status="unknown", form="formatted" )
     write(funit,'("periods 1")')
     write(funit,'("begin filament")')
     write(funit,'("mirror NIL")')
     do icoil = 1, Ncoils
        do ii = 0, coil(icoil)%NS-1
           write(funit,1010) coil(icoil)%xx(ii), coil(icoil)%yy(ii), coil(icoil)%zz(ii), coil(icoil)%I
        enddo
        ii =  coil(icoil)%NS
        write(funit,1010) coil(icoil)%xx(ii), coil(icoil)%yy(ii), coil(icoil)%zz(ii), &
             zero, icoil, coil(icoil)%name
     enddo ! end of do icoil; 14 Apr 16;
     write(funit,'("end")')
     close(funit)

  endif

1010 format(4es23.15,:,i9,"  ",a10)


  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  return

end subroutine restart

