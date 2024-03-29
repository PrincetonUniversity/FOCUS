
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
  use mpi
  use hdf5

  implicit none

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!


  INTEGER            :: ii, jj, icoil, NF,NCP, ip, is, cs, Npc

  ! the following are used by the macros HWRITEXX below; do not alter/remove;
  INTEGER            :: hdfier, rank
  integer(hid_t)     :: file_id, space_id, dset_id                 ! warning: the string "INTEGER" is a macro;
  integer(hsize_t)   :: onedims(1:1), twodims(1:2), threedims(1:3) ! warning: the string "INTEGER" is a macro;
  INTEGER            :: imn, NSmax
  REAL               :: tvolume
  CHARACTER          :: srestart*6
  REAL, allocatable  :: tempvar(:,:)


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
     if(allocated(t1H)) deriv(1:Ndof,5) = t1H(1:Ndof)
     if(allocated(t1K)) deriv(1:Ndof,6) = t1K(1:Ndof)
     if(allocated(t1C)) deriv(1:Ndof,7) = t1C(1:Ndof)
     if(allocated(t1T)) deriv(1:Ndof,8) = t1T(1:Ndof)
     if(allocated(t1N)) deriv(1:Ndof,9) = t1N(1:Ndof)
     if(allocated(t1Bavg)) deriv(1:Ndof,10) = t1Bavg(1:Ndof)
     if(allocated(t1Str)) deriv(1:Ndof,11)=t1Str(1:Ndof)
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
  HWRITECH( 100              ,   input_surf    ,   input_surf                    )
  HWRITECH( 100              ,   input_coils   ,   input_coils                   )
  HWRITECH( 100              ,   input_harm    ,   input_harm                    )
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
  HWRITEIV( 1                ,   case_curv     ,   case_curv                     )
  HWRITEIV( 1                ,   case_straight ,   case_straight                 )
  HWRITEIV( 1                ,   curv_alpha    ,   curv_alpha                    )
  HWRITERV( 1                ,   weight_bnorm  ,   weight_bnorm                  )
  HWRITERV( 1                ,   weight_sbnorm ,   weight_sbnorm                 )
  HWRITERV( 1                ,   weight_bharm  ,   weight_bharm                  )
  HWRITERV( 1                ,   weight_tflux  ,   weight_tflux                  )
  HWRITERV( 1                ,   target_tflux  ,   target_tflux                  )
  HWRITERV( 1                ,   weight_isum   ,   weight_isum                   )
  HWRITERV( 1                ,   target_isum   ,   target_isum                   )
  HWRITERV( 1                ,   weight_ttlen  ,   weight_ttlen                  )
  HWRITERV( 1                ,   target_length ,   target_length                 )
  HWRITERV( 1                ,   curv_k0       ,   curv_k0                       ) 
  HWRITERV( 1                ,   weight_specw  ,   weight_specw                  )
  HWRITERV( 1                ,   weight_cssep  ,   weight_cssep                  )
  HWRITERV( 1                ,   weight_gnorm  ,   weight_gnorm                  )
  HWRITERV( 1                ,   weight_inorm  ,   weight_inorm                  )
  HWRITERV( 1                ,   weight_mnorm  ,   weight_mnorm                  )
  HWRITERV( 1                ,   weight_curv   ,   weight_curv                   )
  HWRITERV( 1                ,   weight_straight,   weight_straight              )
  HWRITERV( 1                ,   weight_ccsep  ,   weight_ccsep                  )
  HWRITERV( 1                ,   weight_tors   ,   weight_tors                   )
  HWRITERV( 1                ,   ccsep_alpha   ,   ccsep_alpha                   )
  HWRITERV( 1                ,   ccsep_beta    ,   ccsep_beta                    )
  HWRITERV( 1                ,   ccsep_skip    ,   ccsep_skip                    )
  HWRITERV( 1                ,   tors_alpha    ,   tors_alpha                    )
  HWRITERV( 1                ,   case_tors     ,   case_tors                     )
  HWRITERV( 1                ,   tors0         ,   tors0                         )
  HWRITERV( 1                ,   nissin_alpha  ,   nissin_alpha                  )
  HWRITERV( 1                ,   nissin_beta   ,   nissin_beta                   )
  HWRITERV( 1                ,   penfun_nissin ,   penfun_nissin                 )
  HWRITERV( 1                ,   nissin0       ,   nissin0                       )
  HWRITERV( 1                ,   nissin_sigma  ,   nissin_sigma                  )
  HWRITERV( 1                ,   nissin_gamma  ,   nissin_gamma                  )
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
  HWRITEIV( 1                ,   update_plasma ,   update_plasma                 )
  HWRITERV( 1                ,   pp_phi        ,   pp_phi                        )
  HWRITERV( 1                ,   pp_raxis      ,   pp_raxis                      )
  HWRITERV( 1                ,   pp_zaxis      ,   pp_zaxis                      )
  HWRITERV( 1                ,   pp_rmax       ,   pp_rmax                       )
  HWRITERV( 1                ,   pp_zmax       ,   pp_zmax                       )
  HWRITEIV( 1                ,   pp_ns         ,   pp_ns                         )
  HWRITEIV( 1                ,   pp_maxiter    ,   pp_maxiter                    )
  HWRITERV( 1                ,   pp_xtol       ,   pp_xtol                       )
  HWRITEIV( 1                ,   axis_npoints  ,   axis_npoints                  )

  HWRITEIV( 1                ,   Nfp           ,   surf(plasma)%Nfp                     )
  HWRITERV( 1                ,   surf_vol      ,   surf(plasma)%vol                     )
  HWRITERA( Nteta,Nzeta      ,   xsurf         ,   surf(plasma)%xx(0:Nteta-1,0:Nzeta-1) )
  HWRITERA( Nteta,Nzeta      ,   ysurf         ,   surf(plasma)%yy(0:Nteta-1,0:Nzeta-1) )
  HWRITERA( Nteta,Nzeta      ,   zsurf         ,   surf(plasma)%zz(0:Nteta-1,0:Nzeta-1) )
  HWRITERA( Nteta,Nzeta      ,   nx            ,   surf(plasma)%nx(0:Nteta-1,0:Nzeta-1) )
  HWRITERA( Nteta,Nzeta      ,   ny            ,   surf(plasma)%ny(0:Nteta-1,0:Nzeta-1) )
  HWRITERA( Nteta,Nzeta      ,   nz            ,   surf(plasma)%nz(0:Nteta-1,0:Nzeta-1) )
  HWRITERA( Nteta,Nzeta      ,   nn            ,   surf(plasma)%ds(0:Nteta-1,0:Nzeta-1) )

  if (allocated(bn)) then
     HWRITERA( Nteta,Nzeta      ,   plas_Bn       ,   surf(plasma)%pb(0:Nteta-1,0:Nzeta-1) )
     HWRITERA( Nteta,Nzeta      ,        Bn       ,   surf(plasma)%Bn(0:Nteta-1,0:Nzeta-1) )
     HWRITERA( Nteta,Nzeta      ,   Bx            ,   surf(plasma)%Bx(0:Nteta-1,0:Nzeta-1) )
     HWRITERA( Nteta,Nzeta      ,   By            ,   surf(plasma)%By(0:Nteta-1,0:Nzeta-1) )
     HWRITERA( Nteta,Nzeta      ,   Bz            ,   surf(plasma)%Bz(0:Nteta-1,0:Nzeta-1) )
  endif

  if (weight_sbnorm > 0) then
      NSmax = 0
      do icoil = 1, Ncoils
         if(coil(icoil)%NS .gt. NSmax) NSmax = coil(icoil)%NS
      enddo
      SALLOCATE(tempvar, (1:Ncoils, 0:NSmax) , 0.0)

      ! Save coil data
      do icoil = 1, Ncoils
         tempvar(icoil,0:coil(icoil)%NS) = coil(icoil)%xx(0:coil(icoil)%NS)
      enddo
      HWRITERA( Ncoils, NSmax+1, xx  , tempvar(1:Ncoils,0:NSmax) )
      tempvar(1:Ncoils,0:NSmax) = 0.0
      do icoil = 1, Ncoils
         tempvar(icoil,0:coil(icoil)%NS) = coil(icoil)%yy(0:coil(icoil)%NS)
      enddo
      HWRITERA( Ncoils, NSmax+1, yy  , tempvar(1:Ncoils,0:NSmax) )
      tempvar(1:Ncoils,0:NSmax) = 0.0
      do icoil = 1, Ncoils
         tempvar(icoil,0:coil(icoil)%NS) = coil(icoil)%zz(0:coil(icoil)%NS)
      enddo
      HWRITERA( Ncoils, NSmax+1, zz  , tempvar(1:Ncoils,0:NSmax) )
      tempvar(1:Ncoils,0:NSmax) = 0.0
      do icoil = 1, Ncoils
         tempvar(icoil,0:coil(icoil)%NS) = coil(icoil)%xt(0:coil(icoil)%NS)
      enddo
      HWRITERA( Ncoils, NSmax+1, xt  , tempvar(1:Ncoils,0:NSmax) )
      tempvar(1:Ncoils,0:NSmax) = 0.0
      do icoil = 1, Ncoils
         tempvar(icoil,0:coil(icoil)%NS) = coil(icoil)%yt(0:coil(icoil)%NS)
      enddo
      HWRITERA( Ncoils, NSmax+1, yt  , tempvar(1:Ncoils,0:NSmax) )
      tempvar(1:Ncoils,0:NSmax) = 0.0
      do icoil = 1, Ncoils
         tempvar(icoil,0:coil(icoil)%NS) = coil(icoil)%zt(0:coil(icoil)%NS)
      enddo
      HWRITERA( Ncoils, NSmax+1, zt  , tempvar(1:Ncoils,0:NSmax) )
      tempvar(1:Ncoils,0:NSmax) = 0.0
  endif 

  ! Save filamentary body force loads
  if (filforce .eq. 1) then
     do icoil = 1, Ncoils
        tempvar(icoil,0:coil(icoil)%NS) = coil(icoil)%Bxx(0:coil(icoil)%NS)
     enddo
     HWRITERA( Ncoils, NSmax+1, Bxx  , tempvar(1:Ncoils,0:NSmax) )
     tempvar(1:Ncoils,0:NSmax) = 0.0
     do icoil = 1, Ncoils
        tempvar(icoil,0:coil(icoil)%NS) = coil(icoil)%Byy(0:coil(icoil)%NS)
     enddo
     HWRITERA( Ncoils, NSmax+1, Byy  , tempvar(1:Ncoils,0:NSmax) )
     tempvar(1:Ncoils,0:NSmax) = 0.0
     do icoil = 1, Ncoils
        tempvar(icoil,0:coil(icoil)%NS) = coil(icoil)%Bzz(0:coil(icoil)%NS)
     enddo
     HWRITERA( Ncoils, NSmax+1, Bzz  , tempvar(1:Ncoils,0:NSmax) )
     tempvar(1:Ncoils,0:NSmax) = 0.0
     do icoil = 1, Ncoils
        tempvar(icoil,0:coil(icoil)%NS) = coil(icoil)%Fx(0:coil(icoil)%NS)
     enddo
     HWRITERA( Ncoils, NSmax+1, Fx  , tempvar(1:Ncoils,0:NSmax) )
     tempvar(1:Ncoils,0:NSmax) = 0.0
     do icoil = 1, Ncoils
        tempvar(icoil,0:coil(icoil)%NS) = coil(icoil)%Fy(0:coil(icoil)%NS)
     enddo
     HWRITERA( Ncoils, NSmax+1, Fy  , tempvar(1:Ncoils,0:NSmax) )
     tempvar(1:Ncoils,0:NSmax) = 0.0
     do icoil = 1, Ncoils
        tempvar(icoil,0:coil(icoil)%NS) = coil(icoil)%Fz(0:coil(icoil)%NS)
     enddo
     HWRITERA( Ncoils, NSmax+1, Fz  , tempvar(1:Ncoils,0:NSmax) )
     tempvar(1:Ncoils,0:NSmax) = 0.0
  endif
  
  ! Save finite-build coil frame
  if (calcfb .eq. 1) then
    do icoil = 1, Ncoils
       tempvar(icoil,0:coil(icoil)%NS) = coil(icoil)%nfbx(0:coil(icoil)%NS)
    enddo
    HWRITERA( Ncoils, NSmax+1, nfbx  , tempvar(1:Ncoils,0:NSmax) )
    tempvar(1:Ncoils,0:NSmax) = 0.0
    do icoil = 1, Ncoils
       tempvar(icoil,0:coil(icoil)%NS) = coil(icoil)%nfby(0:coil(icoil)%NS)
    enddo
    HWRITERA( Ncoils, NSmax+1, nfby  , tempvar(1:Ncoils,0:NSmax) )
    tempvar(1:Ncoils,0:NSmax) = 0.0
    do icoil = 1, Ncoils
       tempvar(icoil,0:coil(icoil)%NS) = coil(icoil)%nfbz(0:coil(icoil)%NS)
    enddo
    HWRITERA( Ncoils, NSmax+1, nfbz  , tempvar(1:Ncoils,0:NSmax) )
    tempvar(1:Ncoils,0:NSmax) = 0.0
    do icoil = 1, Ncoils
       tempvar(icoil,0:coil(icoil)%NS) = coil(icoil)%bfbx(0:coil(icoil)%NS)
    enddo
    HWRITERA( Ncoils, NSmax+1, bfbx  , tempvar(1:Ncoils,0:NSmax) )
    tempvar(1:Ncoils,0:NSmax) = 0.0
    do icoil = 1, Ncoils
       tempvar(icoil,0:coil(icoil)%NS) = coil(icoil)%bfby(0:coil(icoil)%NS)
    enddo
    HWRITERA( Ncoils, NSmax+1, bfby  , tempvar(1:Ncoils,0:NSmax) )
    tempvar(1:Ncoils,0:NSmax) = 0.0
    do icoil = 1, Ncoils
       tempvar(icoil,0:coil(icoil)%NS) = coil(icoil)%bfbz(0:coil(icoil)%NS)
    enddo
    HWRITERA( Ncoils, NSmax+1, bfbz  , tempvar(1:Ncoils,0:NSmax) )
    tempvar(1:Ncoils,0:NSmax) = 0.0

    DALLOCATE(tempvar)
  endif

  HWRITEIV( 1                ,   iout          ,   iout                          )
  HWRITERV( 1                ,   Inorm         ,   Inorm                         )
  HWRITERV( 1                ,   Gnorm         ,   Gnorm                         )
  HWRITERV( 1                ,   Mnorm         ,   Mnorm                         )
  HWRITERV( 1                ,   overlap       ,   overlap                       )
  HWRITERA( iout, 15         ,   evolution     ,   evolution(1:iout, 0:13)       )
  if (allocated(coilspace)) then
     HWRITERA( iout, Tdof       ,   coilspace     ,   coilspace(1:iout, 1:Tdof)     )
  endif

  if (allocated(deriv)) then
     HWRITERA( Ndof, 12      ,   deriv        ,    deriv(1:Ndof, 0:11)          )
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
     HWRITERV( Ncoils        , coil_importance ,  coil_importance                )
  endif
  
  if (allocated(LM_fvec)) then
     HWRITEIV( 1                ,   ibnorm        ,   ibnorm                     )
     HWRITEIV( 1                ,   mbnorm        ,   mbnorm                     )     
     HWRITEIV( 1                ,   ibharm        ,   ibharm                     )
     HWRITEIV( 1                ,   mbharm        ,   mbharm                     )     
     HWRITEIV( 1                ,   itflux        ,   itflux                     )
     HWRITEIV( 1                ,   mtflux        ,   mtflux                     ) 
     HWRITEIV( 1                ,   iisum         ,   iisum                      )
     HWRITEIV( 1                ,   misum         ,   misum                      )          
     HWRITEIV( 1                ,   ittlen        ,   ittlen                     )
     HWRITEIV( 1                ,   mttlen        ,   mttlen                     )     
     HWRITEIV( 1                ,   icssep        ,   icssep                     )
     HWRITEIV( 1                ,   mcssep        ,   mcssep                     )
     HWRITEIV( 1                ,   icurv         ,   icurv                      )
     HWRITEIV( 1                ,   mcurv         ,   mcurv                      )
     HWRITEIV( 1                ,   istr          ,   istr                       )
     HWRITEIV( 1                ,   mstr          ,   mstr                       )
     HWRITEIV( 1                ,   iccsep        ,   iccsep                     )
     HWRITEIV( 1                ,   mccsep        ,   mccsep                     )
     HWRITEIV( 1                ,   itors         ,   itors                      )
     HWRITEIV( 1                ,   mtors         ,   mtors                      )
     HWRITEIV( 1                ,   inissin       ,   inissin                    )
     HWRITEIV( 1                ,   mnissin       ,   mnissin                    )
     HWRITERV( LM_mfvec         ,   LM_fvec       ,   LM_fvec                    )
     HWRITERA( LM_mfvec, Ndof   ,   LM_fjac       ,   LM_fjac                    )     
  endif

  if (allocated(ppr)) then
     HWRITERA( pp_ns, pp_maxiter+1,   ppr         ,  ppr(1:pp_ns, 0:pp_maxiter) )
     HWRITERA( pp_ns, pp_maxiter+1,   ppz         ,  ppz(1:pp_ns, 0:pp_maxiter) )
     HWRITERV( pp_ns              ,   iota        ,  iota(1:pp_ns)              )
     HWRITERV( axis_npoints       ,   axis_phi    ,  axis_phi(1:axis_npoints)   )
     HWRITERV( axis_npoints       ,   axis_r      ,  axis_r(1:axis_npoints)     )
     HWRITERV( axis_npoints       ,   axis_z      ,  axis_z(1:axis_npoints)     )
  endif

  if (allocated(XYZB)) then
     HWRITERC( total_num,4, pp_ns ,   XYZB        ,   XYZB(1:total_num, 1:4, 1:pp_ns) )
     HWRITERA( booz_mn,      pp_ns,  booz_mnc     ,   booz_mnc(1:booz_mn, 1:pp_ns)    )    
     HWRITERA( booz_mn,      pp_ns,  booz_mns     ,   booz_mns(1:booz_mn, 1:pp_ns)    )    
     HWRITEIV( booz_mn,              bmim         ,   bmim(1:booz_mn)                 )  
     HWRITEIV( booz_mn,              bmin         ,   bmin(1:booz_mn)                 )  
  endif

  HWRITERV( 1                ,  time_initialize,   time_initialize               )
  HWRITERV( 1                ,  time_optimize  ,   time_optimize                 )
  HWRITERV( 1                ,  time_postproc  ,   time_postproc                 )

  ! do icoil = 1,Ncoils
  !   select case(icoil)
  !     case(1)
  !       HWRITERV( coil(icoil)%NS            ,  curvature_1,   coil(icoil)%curvature(1:coil(icoil)%NS )              )
  !       HWRITERV( coil(icoil)%NS            ,  straight_1,   coil(icoil)%straight (1:coil(icoil)%NS )              ) 
  !     case(2)
  !       HWRITERV( coil(icoil)%NS            ,  curvature_2,   coil(icoil)%curvature(1:coil(icoil)%NS )              )
  !       HWRITERV( coil(icoil)%NS            ,  straight_2,   coil(icoil)%straight (1:coil(icoil)%NS )              ) 
  !     case(3)
  !       HWRITERV( coil(icoil)%NS            ,  curvature_3,   coil(icoil)%curvature(1:coil(icoil)%NS )              )
  !       HWRITERV( coil(icoil)%NS            ,  straight_3,   coil(icoil)%straight (1:coil(icoil)%NS )              ) 
  !     case(4)
  !       HWRITERV( coil(icoil)%NS            ,  curvature_4,   coil(icoil)%curvature(1:coil(icoil)%NS )              )
  !       HWRITERV( coil(icoil)%NS            ,  straight_4,   coil(icoil)%straight (1:coil(icoil)%NS )              ) 
  !     case(5)
  !       HWRITERV( coil(icoil)%NS            ,  curvature_5,   coil(icoil)%curvature(1:coil(icoil)%NS )              )
  !       HWRITERV( coil(icoil)%NS            ,  straight_5,   coil(icoil)%straight (1:coil(icoil)%NS )              ) 
  !     case(6)
  !       HWRITERV( coil(icoil)%NS            ,  curvature_6,   coil(icoil)%curvature(1:coil(icoil)%NS )              )
  !       HWRITERV( coil(icoil)%NS            ,  straight_6,   coil(icoil)%straight (1:coil(icoil)%NS )              ) 
  !     case(7)
  !       HWRITERV( coil(icoil)%NS            ,  curvature_7,   coil(icoil)%curvature(1:coil(icoil)%NS )              )
  !       HWRITERV( coil(icoil)%NS            ,  straight_7,   coil(icoil)%straight (1:coil(icoil)%NS )              ) 
  !     case(8)
  !       HWRITERV( coil(icoil)%NS            ,  curvature_8,   coil(icoil)%curvature(1:coil(icoil)%NS )              )
  !       HWRITERV( coil(icoil)%NS            ,  straight_8,   coil(icoil)%straight (1:coil(icoil)%NS )              ) 
  !     case(9)
  !       HWRITERV( coil(icoil)%NS            ,  curvature_9,   coil(icoil)%curvature(1:coil(icoil)%NS )              )
  !       HWRITERV( coil(icoil)%NS            ,  straight_9,   coil(icoil)%straight (1:coil(icoil)%NS )              ) 
  !     end select
  ! end do



  call h5fclose_f( file_id, hdfier ) ! terminate access;
  FATAL( restart, hdfier.ne.0, error calling h5fclose_f )

  call h5close_f( hdfier ) ! close Fortran interface;
  FATAL( restart, hdfier.ne.0, error calling h5close_f )


  !--------------------------write focus coil file-----------------------------------------
  if( save_coils == 1 ) then
     open( wunit, file=trim(out_focus), status="unknown", form="formatted")
     write(wunit, *) "# Total number of coils"
     write(wunit, '(8X,I6)') Ncoils

     do icoil = 1, Ncoils

        write(wunit, *) "#-----------------", icoil, "---------------------------" 
        write(wunit, *) "#coil_type   coil_symm  coil_name"
        write(wunit,'(3X,I3,4X,I3,4X,A10)') coil(icoil)%type, coil(icoil)%symm, coil(icoil)%name

        select case (coil(icoil)%type)
        case (1)
           write(wunit, '(4(A6, A15, 8X))') " #Nseg", "current",  "Ifree", "Length", "Lfree", "target_length"!, "curv_k0"
           !write(wunit,'(2X, I4, ES23.15, 3X, I3, ES23.15, 3X, I3, ES23.15, ES23.15)') &
                !coil(icoil)%NS, coil(icoil)%I, coil(icoil)%Ic, coil(icoil)%L, coil(icoil)%Lc, coil(icoil)%Lo, coil(icoil)%curv_k0
           write(wunit,'(2X, I4, ES23.15, 3X, I3, ES23.15, 3X, I3, ES23.15)') &
                coil(icoil)%NS, coil(icoil)%I, coil(icoil)%Ic, coil(icoil)%L, coil(icoil)%Lc, coil(icoil)%Lo !,coil(icoil)%curv_k0  
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
        case (2) 
           write(wunit, *) "#  Lc  ox   oy   oz  Ic  I  mt  mp"
           write(wunit,'(2(I3, 3ES23.15))') coil(icoil)%Lc, coil(icoil)%ox, coil(icoil)%oy, coil(icoil)%oz, &
                                         coil(icoil)%Ic, coil(icoil)%I , coil(icoil)%mt, coil(icoil)%mp  
        case (3)
           write(wunit, *) "# Ic     I    Lc  Bz  (Ic control I; Lc control Bz)"
           write(wunit,'(I3, ES23.15, I3, ES23.15)') coil(icoil)%Ic, coil(icoil)%I, &
                                                     coil(icoil)%Lc, coil(icoil)%Bz
	case (coil_type_spline)
           write(wunit, '(4(A6, A15, 8X))') " #Nseg", "current",  "Ifree", "Length", "Lfree", "target_length"!, "k0"
           !write(wunit,'(2X, I4, ES23.15, 3X, I3, ES23.15, 3X, I3, ES23.15, ES23.15)') &
                !coil(icoil)%NS, coil(icoil)%I, coil(icoil)%Ic, coil(icoil)%L, coil(icoil)%Lc, coil(icoil)%Lo, coil(icoil)%k0
           write(wunit,'(2X, I4, ES23.15, 3X, I3, ES23.15, 3X, I3, ES23.15)') &
                coil(icoil)%NS, coil(icoil)%I, coil(icoil)%Ic, coil(icoil)%L, coil(icoil)%Lc, coil(icoil)%Lo !,coil(icoil)%k0  
           NCP = Splines(icoil)%NCP ! shorthand;
           write(wunit, *) "#NCP "
           write(wunit, '(I3,I3)') NCP 
           write(wunit, *) "knot vector"
           write(wunit, 1000) Splines(icoil)%vect
           write(wunit, *) "#Control points for coils ( x;y;z) "
           write(wunit, 1000) Splines(icoil)%Cpoints(0:NCP-1)
           write(wunit, 1000) Splines(icoil)%Cpoints(NCP:2*NCP-1)
           write(wunit, 1000) Splines(icoil)%Cpoints(NCP*2:3*NCP-1)
        case default
           FATAL(restart, .true., not supported coil types)
        end select
     enddo
     close(wunit)
1000 format(9999ES23.15)
  endif

  !--------------------------write coils.ext file-----------------------------------------------  

  if( save_coils == 1 ) then

     open(funit,file=trim(out_coils), status="unknown", form="formatted" )
     write(funit,'("periods "I3)') surf(plasma)%Nfp
     write(funit,'("begin filament")')
     write(funit,'("mirror NIL")')
     do icoil = 1, Ncoils
        ! will only write x,y,z in cartesian coordinates
        if ((coil(icoil)%type /= 1) .AND. (coil(icoil)%type /= coil_type_spline)) cycle
        ! check if the coil is stellarator symmetric
        select case (coil(icoil)%symm) 
        case ( 0 )
           cs  = 0
           Npc = 1
        case ( 1 )
           cs  = 0
           Npc = Nfp
        case ( 2) 
           cs  = 1
           Npc = Nfp
        end select
        ! periodicity and stellarator symmetry
        do ip = 1, Npc
           do is = 0, cs
              do ii = 0, coil(icoil)%NS-1
                 write(funit,1010) coil(icoil)%xx(ii)*cosnfp(ip)-coil(icoil)%yy(ii)*sinnfp(ip), &
                      &  (-1)**is*(coil(icoil)%xx(ii)*sinnfp(ip)+coil(icoil)%yy(ii)*cosnfp(ip)), &
                      &  (-1)**is*coil(icoil)%zz(ii), (-1)**is*coil(icoil)%I
              enddo
              ii =  0
              write(funit,1010) coil(icoil)%xx(ii)*cosnfp(ip)-coil(icoil)%yy(ii)*sinnfp(ip), &
                   &  (-1)**is*(coil(icoil)%xx(ii)*sinnfp(ip)+coil(icoil)%yy(ii)*cosnfp(ip)),& 
                   &  (-1)**is*coil(icoil)%zz(ii), zero, icoil, coil(icoil)%name
           enddo
        enddo
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
     do icoil = 1, Ncoils
        write(funit) coil(icoil)%xx(0:coil(icoil)%NS)
        write(funit) coil(icoil)%yy(0:coil(icoil)%NS)
        write(funit) coil(icoil)%zz(0:coil(icoil)%NS)
     enddo
     close( funit )

  endif


  !--------------------------write ext.harmonics file-----------------------------------

  if (save_harmonics == 1 .and. allocated(Bmnc)) then

     open(wunit, file=trim(out_harm), status='unknown', action='write')
     write(wunit,'("#NBmn")')                     ! comment line;
     write(wunit,'(I6)') NBmn                     ! write dimensions
     write(wunit,'("# n  m   Bmnc  Bmns  wBmn")') ! comment line;
     do imn = 1, NBmn
        write(wunit,'(2(I3, 4X), 3(ES23.15,4X))') Bmnin(imn)/surf(plasma)%Nfp, & 
             Bmnim(imn), Bmnc(imn), Bmns(imn), wBmn(imn)
     enddo
     close(wunit)

  endif

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if (update_plasma == 1 ) call write_plasma

  return

end subroutine saving


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

SUBROUTINE write_plasma
!-------------------------------------------------------------------------------!
! write down the unpdated plasma boundary information;                          !
! CZHU; first version: 2017/01/11; last revised: 2017/01/11                     !
!-------------------------------------------------------------------------------!
  use globals, only : dp, zero, half, two, pi2, myid, ncpu, ounit, wunit, ext, &
                      plasma, MPI_COMM_FOCUS, IsSymmetric, &
                      Nteta, Nzeta, surf, bnorm, sqrtmachprec, out_plasma
  use mpi
  implicit none  

  !-------------------------------------------------------------------------------
  INTEGER             :: mf, nf  ! predefined Fourier modes size
  INTEGER             :: imn=0, ii, jj, im, in, astat, ierr, maxN, maxM, isurf
  REAL                :: teta, zeta, arg, tol, tmpc, tmps
  !------------------------------------------------------------------------------- 

  ! use plasma as default
  isurf = plasma
  mf = 24 ;  nf = 24
  FATAL(bnftran, mf .le. 0 .and. nf .le. 0, INVALID size for Fourier harmonics)

  tmpc = zero ; tmps = zero

  if (bnorm .gt. sqrtmachprec ) then
     tol = 1.0E-8 * bnorm
  else
     tol = 1.0E-8
  endif

  if(myid .ne. 0) return
  FATAL( write_plasma, IsSymmetric==2, option not supported for now)

  if(surf(isurf)%Nbnf .gt. 0) then  ! if there is input Bn target
     DALLOCATE(surf(isurf)%bnim)
     DALLOCATE(surf(isurf)%bnin)
     DALLOCATE(surf(isurf)%bnc )
     DALLOCATE(surf(isurf)%bns )
  endif

  surf(isurf)%Nbnf = (mf+1)*(2*nf+1) ! (0:mf)*(-nf:nf)

  SALLOCATE( surf(isurf)%bnim, (1:surf(isurf)%Nbnf), 0    )
  SALLOCATE( surf(isurf)%bnin, (1:surf(isurf)%Nbnf), 0    )
  SALLOCATE( surf(isurf)%bnc , (1:surf(isurf)%Nbnf), zero )
  SALLOCATE( surf(isurf)%bns , (1:surf(isurf)%Nbnf), zero )  
  
  imn = 0
  do in = -nf, nf
     do im = 0, mf
        tmpc = zero ; tmps = zero
        do jj = 0, Nzeta-1
           zeta = ( jj + half ) * pi2 / surf(isurf)%Nzeta
           do ii = 0, Nteta-1
              teta = ( ii + half ) * pi2 / surf(isurf)%Nteta
              arg = im*teta - in*surf(isurf)%Nfp*zeta
              tmpc = tmpc + surf(isurf)%bn(ii,jj)*cos(arg)
              tmps = tmps + surf(isurf)%bn(ii,jj)*sin(arg)

           enddo ! end jj
        enddo ! end ii

        if ( (abs(tmpc) + abs(tmps)) .lt. tol ) cycle

        imn = imn + 1
        surf(isurf)%bnin(imn) = in * surf(isurf)%Nfp
        surf(isurf)%bnim(imn) = im

        if (im .eq. 0  ) then
           tmpc = tmpc*half
           tmps = tmps*half
        endif
        surf(isurf)%bnc(imn) = tmpc
        surf(isurf)%bns(imn) = tmps

     enddo ! end im
  enddo ! end in

  surf(isurf)%Nbnf = imn

  surf(isurf)%bnc = surf(isurf)%bnc * two / (Nteta*Nzeta)
  surf(isurf)%bns = surf(isurf)%bns * two / (Nteta*Nzeta)
  !----------------------------------------------


  open(wunit, file=trim(out_plasma), status='unknown', action='write')

  write(wunit,*      ) "#Nfou Nfp  Nbnf"
  write(wunit,'(3I6)' ) surf(isurf)%Nfou, surf(isurf)%Nfp, surf(isurf)%Nbnf

  write(wunit,*      ) "#------- plasma boundary------"
  write(wunit,*      ) "#  n   m   Rbc   Rbs    Zbc   Zbs"
  do imn = 1, surf(isurf)%Nfou
     write(wunit,'(2I6, 4ES15.6)') surf(isurf)%bin(imn)/surf(isurf)%Nfp, surf(isurf)%bim(imn), & 
          surf(isurf)%Rbc(imn), surf(isurf)%Rbs(imn), surf(isurf)%Zbc(imn), surf(isurf)%Zbs(imn)
  enddo

  write(wunit,*      ) "#-------Bn harmonics----------"
  write(wunit,*      ) "#  n  m  bnc   bns"
  if (surf(isurf)%Nbnf .gt. 0) then
     do imn = 1, surf(isurf)%Nbnf
        write(wunit,'(2I6, 2ES15.6)') surf(isurf)%bnin(imn)/surf(isurf)%Nfp, surf(isurf)%bnim(imn), &
             surf(isurf)%bnc(imn), surf(isurf)%bns(imn)
     enddo
  else
     write(wunit,'(2I6, 2ES15.6)') 0, 0, 0.0, 0.0
  endif

  close(wunit)
END SUBROUTINE write_plasma

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
