
!title (restart) ! Writes output and restart files.

!latex \briefly{Writes output and restart files.}

!latex \calledby{\link{solvers}}
!l tex \calls{}

!latex \tableofcontents

!latex \section{h5 format}
!latex \bi
!latex \item[1.] All the restart information etc. is written to file.
!latex \ei

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine saving

  use globals

  use hdf5

  implicit none

  include "mpif.h"

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!


  INTEGER            :: ii, jj, icoil, NF

  ! the following are used by the macros HWRITEXX below; do not alter/remove;
  INTEGER            :: hdfier, rank
  integer(hid_t)     :: file_id, space_id, dset_id                 ! warning: the string "INTEGER" is a macro;
  integer(hsize_t)   :: onedims(1:1), twodims(1:2), threedims(1:3) ! warning: the string "INTEGER" is a macro;
  INTEGER            :: imn
  REAL               :: tvolume

  CHARACTER          :: srestart*6


  if( myid.ne.0 ) return

  !SALLOCATE(perA, (1:Tdof), zero)

  !store the current derivatives;
  if (allocated(deriv)) then
     deriv = zero
     if(allocated(t1E)) deriv(1:Ndof,0) = t1E(1:Ndof)
     if(allocated(t1B)) deriv(1:Ndof,1) = t1B(1:Ndof)
     if(allocated(t1F)) deriv(1:Ndof,2) = t1F(1:Ndof)
     if(allocated(t1L)) deriv(1:Ndof,3) = t1L(1:Ndof)
     if(allocated(t1S)) deriv(1:Ndof,4) = t1S(1:Ndof)
     if(allocated(t1C)) deriv(1:Ndof,5) = t1C(1:Ndof)
     if(allocated(t1H)) deriv(1:Ndof,6) = t1H(1:Ndof)
  endif

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  call h5open_f( hdfier ) ! initialize Fortran interface;
  FATAL( restart, hdfier.ne.0, error calling h5open_f )

  call h5fcreate_f( trim(hdf5file), H5F_ACC_TRUNC_F, file_id, hdfier ) ! create new file;
  FATAL( restart, hdfier.ne.0, error calling h5fcreate_f )

  HWRITECH( 10               ,   version       ,   version                       )

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
  HWRITEIV( 1                ,   NFcoil        ,   NFcoil                        )
  HWRITEIV( 1                ,   Nseg          ,   Nseg                          )
  HWRITEIV( 1                ,   case_optimize ,   case_optimize                 )
  HWRITERV( 1                ,   exit_xtol     ,   exit_tol                      )
  HWRITEIV( 1                ,   IsNormalize   ,   IsNormalize                   )
  HWRITEIV( 1                ,   ISNormWeight  ,   IsNormWeight                  )
  HWRITEIV( 1                ,   case_bnormal  ,   case_bnormal                  )
  HWRITEIV( 1                ,   case_length   ,   case_length                   )
  HWRITERV( 1                ,   weight_bnorm  ,   weight_bnorm                  )
  HWRITERV( 1                ,   weight_bharm  ,   weight_bharm                  )
  HWRITERV( 1                ,   weight_tflux  ,   weight_tflux                  )
  HWRITERV( 1                ,   target_tflux  ,   target_tflux                  )
  HWRITERV( 1                ,   weight_ttlen  ,   weight_ttlen                  )
  HWRITERV( 1                ,   target_length ,   target_length                 )
  HWRITERV( 1                ,   weight_specw  ,   weight_specw                  )
  HWRITERV( 1                ,   weight_ccsep  ,   weight_ccsep                  )
  HWRITERV( 1                ,   weight_gnorm  ,   weight_gnorm                  )
  HWRITERV( 1                ,   weight_inorm  ,   weight_inorm                  )
  HWRITERV( 1                ,   DF_tausta     ,   DF_tausta                     )
  HWRITERV( 1                ,   DF_tauend     ,   DF_tauend                     )
  HWRITERV( 1                ,   DF_xtol       ,   DF_xtol                       )
  HWRITEIV( 1                ,   DF_maxiter    ,   DF_maxiter                    )
  HWRITEIV( 1                ,   CG_maxiter    ,   CG_maxiter                    )
  HWRITERV( 1                ,   CG_xtol       ,   CG_xtol                       )
  HWRITERV( 1                ,   CG_wolfe_c1   ,   CG_wolfe_c1                   )
  HWRITERV( 1                ,   CG_wolfe_c2   ,   CG_wolfe_c2                   )
  HWRITEIV( 1                ,   HN_maxiter    ,   HN_maxiter                    )
  HWRITERV( 1                ,   HN_xtol       ,   HN_xtol                       )
  HWRITERV( 1                ,   HN_factor     ,   HN_factor                     )
  HWRITEIV( 1                ,   TN_maxiter    ,   TN_maxiter                    )
  HWRITEIV( 1                ,   TN_reorder    ,   TN_reorder                    )
  HWRITERV( 1                ,   TN_xtol       ,   TN_xtol                       )
  HWRITERV( 1                ,   TN_cr         ,   TN_cr                         )
  HWRITEIV( 1                ,   case_postproc ,   case_postproc                 )
  HWRITEIV( 1                ,   save_freq     ,   save_freq                     )
  HWRITEIV( 1                ,   save_coils    ,   save_coils                    )
  HWRITEIV( 1                ,   save_harmonics,   save_harmonics                )
  HWRITEIV( 1                ,   save_filaments,   save_filaments                )

  HWRITEIV( 1                ,   Nfp           ,   Nfp_raw                         )
  HWRITERV( 1                ,   surf_vol      ,   surf(1)%vol                     )
  HWRITERA( Nteta,Nzeta      ,   xsurf         ,   surf(1)%xx(0:Nteta-1,0:Nzeta-1) )
  HWRITERA( Nteta,Nzeta      ,   ysurf         ,   surf(1)%yy(0:Nteta-1,0:Nzeta-1) )
  HWRITERA( Nteta,Nzeta      ,   zsurf         ,   surf(1)%zz(0:Nteta-1,0:Nzeta-1) )
  HWRITERA( Nteta,Nzeta      ,   nx            ,   surf(1)%nx(0:Nteta-1,0:Nzeta-1) )
  HWRITERA( Nteta,Nzeta      ,   ny            ,   surf(1)%ny(0:Nteta-1,0:Nzeta-1) )
  HWRITERA( Nteta,Nzeta      ,   nz            ,   surf(1)%nz(0:Nteta-1,0:Nzeta-1) )

  if (allocated(bn)) then
     HWRITERA( Nteta,Nzeta      ,   plas_Bn       ,   surf(1)%pb(0:Nteta-1,0:Nzeta-1) )
     HWRITERA( Nteta,Nzeta      ,        Bn       ,   surf(1)%bn(0:Nteta-1,0:Nzeta-1) )
     HWRITERA( Nteta,Nzeta      ,   Bx            ,   surf(1)%Bx(0:Nteta-1,0:Nzeta-1) )
     HWRITERA( Nteta,Nzeta      ,   By            ,   surf(1)%By(0:Nteta-1,0:Nzeta-1) )
     HWRITERA( Nteta,Nzeta      ,   Bz            ,   surf(1)%Bz(0:Nteta-1,0:Nzeta-1) )
  endif

  HWRITEIV( 1                ,   iout          ,   iout                          )
  HWRITERV( 1                ,   Inorm         ,   Inorm                         )
  HWRITERV( 1                ,   Gnorm         ,   Gnorm                         )
  HWRITERV( 1                ,   overlap       ,   overlap                       )
  HWRITERA( iout, 8          ,   evolution     ,   evolution(1:iout, 0:7)        )
  HWRITERA( iout, Tdof       ,   coilspace     ,   coilspace(1:iout, 1:Tdof)     )

  if (allocated(deriv)) then
     HWRITERA( Ndof, 6       ,   deriv         ,   deriv(1:Ndof, 0:6)            )
  endif

  if (allocated(Bmnc)) then
     HWRITEIV( NBmn          ,   Bmnin         ,   Bmnin                         )
     HWRITEIV( NBmn          ,   Bmnim         ,   Bmnim                         )
     HWRITERV( NBmn          ,   initial_Bmnc  ,   iBmnc                         )
     HWRITERV( NBmn          ,   initial_Bmns  ,   iBmns                         )
     HWRITERV( NBmn          ,   target_Bmnc   ,   tBmnc                         )
     HWRITERV( NBmn          ,   target_Bmns   ,   tBmns                         )
     HWRITERV( NBmn          ,          Bmnc   ,    Bmnc                         )
     HWRITERV( NBmn          ,          Bmns   ,    Bmns                         )
  endif

  if (allocated(coil_importance)) then
     HWRITERV( Ncoils*Npc    , coil_importance ,  coil_importance                )
  endif
  
  if (allocated(LM_fvec)) then
     HWRITEIV( 1                ,   ibnorm        ,   ibnorm                     )
     HWRITEIV( 1                ,   mbnorm        ,   mbnorm                     )     
     HWRITEIV( 1                ,   ibharm        ,   ibharm                     )
     HWRITEIV( 1                ,   mbharm        ,   mbharm                     )     
     HWRITEIV( 1                ,   itflux        ,   itflux                     )
     HWRITEIV( 1                ,   mtflux        ,   mtflux                     )     
     HWRITEIV( 1                ,   ittlen        ,   ittlen                     )
     HWRITEIV( 1                ,   mttlen        ,   mttlen                     )     
     HWRITEIV( 1                ,   icssep        ,   icssep                     )
     HWRITEIV( 1                ,   mcssep        ,   mcssep                     )
     HWRITERV( LM_mfvec         ,   LM_fvec       ,   LM_fvec                    )
     HWRITERA( LM_mfvec, Ndof   ,   LM_fjac       ,   LM_fjac                    )     
  endif

  if (allocated(ppr)) then
     HWRITERA( pp_ns, pp_maxiter+1,   ppr         ,  ppr(1:pp_ns, 0:pp_maxiter) )
     HWRITERA( pp_ns, pp_maxiter+1,   ppz         ,  ppz(1:pp_ns, 0:pp_maxiter) )
     HWRITEIV( 1                  ,   pp_ns       ,  pp_ns                      )
     HWRITERV( pp_ns              ,   iota        ,  iota(1:pp_ns)              )
  endif

  if (allocated(XYZB)) then
     HWRITERA( total_num,4      ,   XYZB          ,   XYZB(1:total_num, 1:4)     )
  endif

  HWRITERV( 1                ,  time_initialize,   time_initialize               )
  HWRITERV( 1                ,  time_optimize  ,   time_optimize                 )
  HWRITERV( 1                ,  time_postproc  ,   time_postproc                 )


  call h5fclose_f( file_id, hdfier ) ! terminate access;
  FATAL( restart, hdfier.ne.0, error calling h5fclose_f )

  call h5close_f( hdfier ) ! close Fortran interface;
  FATAL( restart, hdfier.ne.0, error calling h5close_f )


  !--------------------------write focus coil file-----------------------------------------
  if( save_coils == 1 ) then
     open( wunit, file=trim(coilfile), status="unknown", form="formatted")
     write(wunit, *) "# Total number of coils"
     write(wunit, '(8X,I6)') Ncoils

     do icoil = 1, Ncoils

        write(wunit, *) "#-----------------", icoil, "---------------------------" 
        write(wunit, *) "#coil_type     coil_name"
        write(wunit,'(3X,I3,4X, A10)') coil(icoil)%itype, coil(icoil)%name
        write(wunit, '(3(A6, A15, 8X))') " #Nseg", "current",  "Ifree", "Length", "Lfree", "target_length"
        write(wunit,'(2X, I4, ES23.15, 3X, I3, ES23.15, 3X, I3, ES23.15)') &
             coil(icoil)%NS, coil(icoil)%I, coil(icoil)%Ic, coil(icoil)%L, coil(icoil)%Lc, coil(icoil)%Lo
        select case (coil(icoil)%itype)
        case (1)
           NF = FouCoil(icoil)%NF ! shorthand;
           write(wunit, *) "#NFcoil"
           write(wunit, '(I3)') NF
           write(wunit, *) "#Fourier harmonics for coils ( xc; xs; yc; ys; zc; zs) "
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
  endif

  !--------------------------write coils.ext file-----------------------------------------------  

  if( save_coils == 1 ) then

     open(funit,file=trim(outcoils), status="unknown", form="formatted" )
     write(funit,'("periods "I3)') Nfp_raw
     write(funit,'("begin filament")')
     write(funit,'("mirror NIL")')
     do icoil = 1, Ncoils*Npc
        do ii = 0, coil(icoil)%NS-1
           write(funit,1010) coil(icoil)%xx(ii), coil(icoil)%yy(ii), coil(icoil)%zz(ii), coil(icoil)%I
        enddo
        ii =  coil(icoil)%NS
        write(funit,1010) coil(icoil)%xx(ii), coil(icoil)%yy(ii), coil(icoil)%zz(ii), &
             zero, icoil, coil(icoil)%name
     enddo
     write(funit,'("end")')
     close(funit)

  endif

1010 format(4es23.15,:,i9,"  ",a10)

  !--------------------------write .ext.fo.filaments.xxx file-----------------------------------

  write(srestart,'(i6.6)') iout

  if( save_filaments == 1 ) then

     open( funit, file="."//trim(ext)//".filaments."//srestart, status="unknown", form="unformatted" )  
     write(funit) Ncoils, Nseg
     do icoil = 1, Ncoils*Npc
        write(funit) coil(icoil)%xx(0:coil(icoil)%NS)
        write(funit) coil(icoil)%yy(0:coil(icoil)%NS)
        write(funit) coil(icoil)%zz(0:coil(icoil)%NS)
     enddo
     close( funit )

  endif


  !--------------------------write ext.harmonics file-----------------------------------

  if (save_harmonics == 1 .and. allocated(Bmnc)) then

     open(wunit, file=trim(harmfile), status='unknown', action='write')
     write(wunit,'("#NBmn")')                     ! comment line;
     write(wunit,'(I6)') NBmn                     ! write dimensions
     write(wunit,'("# n  m   Bmnc  Bmns  wBmn")') ! comment line;
     do imn = 1, NBmn
        write(wunit,'(2(I3, 4X), 3(ES23.15,4X))') Bmnin(imn)/Nfp_raw, Bmnim(imn), Bmnc(imn), Bmns(imn), wBmn(imn)
     enddo
     close(wunit)

  endif

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  return

end subroutine saving

